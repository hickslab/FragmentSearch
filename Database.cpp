#define _CRT_SECURE_NO_WARNINGS
#include "Database.h"

#include <sstream>
#include <stdexcept>
#include <chrono>
#include <thread>

#include <sys/stat.h>

struct result_struct {
    std::vector<Protein> load_results;
};

struct thread_args_struct {
    thread_args_struct(result_struct &results, const std::string &source_path, const char *file_buffer, size_t file_size, size_t start_position, size_t n_bytes) : 
        results(results), source_path(source_path), file_buffer(file_buffer), file_size(file_size), start_position(start_position), n_bytes(n_bytes)  { }

    result_struct &results;
    const std::string &source_path;
    const char *file_buffer;
    size_t file_size;
    size_t start_position;
    size_t n_bytes;
};

void LoadThreadProc(thread_args_struct args)
{
    const std::string &source_path = args.source_path;
    const char *file_buffer = args.file_buffer;
    size_t file_size = args.file_size;
    size_t start_position = args.start_position;
    size_t n_bytes = args.n_bytes;

    std::vector<Protein> results;
    size_t estimated_result_count = (n_bytes / 500);
    results.reserve(estimated_result_count);

    // Allocate a container for the sequences
    Protein protein;
    std::string description;
    bool in_description = false;
    bool found_sequence = false;

    // Run through the file and process protein sequences
    for (uint64_t i = start_position; i < file_size; i++) {

        // We shouldn't have any NULs, but just in case...
        if (file_buffer[i] == '\0') continue;

        // Check if this is the beginning of a new sequence
        if (!in_description && file_buffer[i] == '>') {
            // Process the existing sequence
            if (!protein.sequence.empty()) {
                protein.description = std::move(description);
                results.push_back(std::move(protein));
            }
            // Check if we're at the end of what we should be processing
            if (i >= start_position + n_bytes) break;

            in_description = true;
            found_sequence = true;
            protein.source_database = source_path;
            protein.sequence.clear();
            description.clear();
            // Average sequence length is 360 (but max is around 30k). Reserve a healthy amount of sequence memory
            protein.sequence.reserve(360 * 2);
        } else if (!found_sequence) {
            // We're at the end of the prior thread's sequence, nothing to do here
        } else {
            // In an existing sequence. Check if we're still populating the description
            if (in_description) {
                // Check if we've reached the end yet, signified by a newline character
                if (file_buffer[i] == '\r' || file_buffer[i] == '\n') {
                    in_description = false;
                } else {
                    description += file_buffer[i];
                }
            } else {
                // Otherwise we're adding to the sequence.
                // Skip any newline or whitespace characters.
                if (isalpha(file_buffer[i])) protein.sequence.push_back(file_buffer[i]);
            }
        }
    }


    // Process the last sequence
    if (!protein.sequence.empty()) {
        protein.description = std::move(description);
        results.push_back(protein);
    }
    args.results.load_results = std::move(results);
}

Database::Database(std::string path)
{
    this->source_path = path;

    FILE *database_fp = fopen(path.c_str(), "rb");
    if (!database_fp) {
        int err = errno;
        std::stringstream message;
        message << "Unable to open " << path << " for reading: errno " << err << ".\n";
        throw std::invalid_argument(message.str());
    }

    // Read the whole file into ram
    struct stat file_stat;
    stat(path.c_str(), &file_stat);
    uint64_t file_size = file_stat.st_size;
    char *file_buffer = new char[file_size];

    auto start_read = std::chrono::high_resolution_clock::now();
    if (file_size != fread(file_buffer, 1, file_size, database_fp)) {
        int err = errno;
        fclose(database_fp);
        std::stringstream message;
        message << "Unable to read " << file_size << " bytes from " << path << ": errno " << err << ".\n";
        delete[] file_buffer;
        throw std::invalid_argument(message.str());
    }
    auto finish_read = std::chrono::high_resolution_clock::now();
    fclose(database_fp);

    // Pre-allocate space for our sequences
    // The SwissProt database has about 1 sequence / 500 bytes. Use that as a ballpark estimate
    // to reserve some memory.
    size_t n_estimated_sequences = file_size / 500;
    this->sequences.reserve(n_estimated_sequences);

    // Allocate a container for the sequences
    Protein protein;
    std::string description;
    bool in_description = false;

    auto start_processing = std::chrono::high_resolution_clock::now();

    int n_threads = std::thread::hardware_concurrency();
    if (n_threads == 0) n_threads = 1;

    if (n_threads > 1) {
        std::vector<std::thread> threads;
        std::vector<result_struct> results;
        results.resize(n_threads);
        threads.reserve(n_threads);
        auto bytes_per_thread = file_size / n_threads;
        auto remainder = file_size - (bytes_per_thread * n_threads);
        for (int i = 0; i < n_threads; i++) {
            auto start_position = bytes_per_thread * i;
            auto bytes_to_process = bytes_per_thread;
            if (i == n_threads) bytes_to_process += remainder;
            thread_args_struct args(results[i], source_path, file_buffer, file_size, start_position, bytes_to_process);
            auto t = std::thread(LoadThreadProc, args);
            threads.push_back(std::move(t));
        }

        for (int i = 0; i < n_threads; i++) {
            threads[i].join();
            for (auto sequence : results[i].load_results) {
                this->sequences.push_back(std::move(sequence));
            }
        }
        
    } else {
        // Run through the file and process protein sequences
        for (uint64_t i = 0; i < file_size; i++) {

            // We shouldn't have any NULs, but just in case...
            if (file_buffer[i] == '\0') continue;

            // Check if this is the beginning of a new sequence
            if (!in_description && file_buffer[i] == '>') {
                // Process the existing sequence
                if (!protein.sequence.empty()) {
                    protein.description = std::move(description);
                    this->sequences.push_back(std::move(protein));
                }
                in_description = true;
                protein.source_database = this->source_path;
                protein.sequence.clear();
                description.clear();
            } else {
                // In an existing sequence. Check if we're still populating the description
                if (in_description) {
                    // Check if we've reached the end yet, signified by a newline character
                    if (file_buffer[i] == '\r' || file_buffer[i] == '\n') {
                        in_description = false;
                    } else {
                        description += file_buffer[i];
                    }
                } else {
                    // Otherwise we're adding to the sequence.
                    // Skip any newline or whitespace characters.
                    if (isalpha(file_buffer[i])) protein.sequence.push_back(file_buffer[i]);
                }
            }
        }

        // Process the last sequence
        if (!protein.sequence.empty()) {
            protein.description = std::move(description);
            this->sequences.push_back(protein);
        }
    }
    auto finish_processing = std::chrono::high_resolution_clock::now();

    delete[] file_buffer;

    long long read_time = std::chrono::duration_cast<std::chrono::milliseconds>(finish_read - start_read).count();
    long long process_time = std::chrono::duration_cast<std::chrono::milliseconds>(finish_processing - start_processing).count();
    fprintf(stderr, "Imported %zd sequences. Reading time = %lld ms, processing time = %lld ms.\n", this->sequences.size(), read_time, process_time);
}
