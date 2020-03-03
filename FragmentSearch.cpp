#include <thread>
#include <stdexcept>
#include <algorithm>
#include <chrono>

#include "FragmentSearch.h"

// Implementation file for the main program logic


// Amino acid list
double get_amino_acid_mass(char aa)
{
    switch (aa) {
    case 'A':
        return 71.03711;
    case 'R':
        return 156.10111;
    case 'N':
        return 114.04293;
    case 'D':
        return 115.02694;
    case 'C':
        return 103.00919 + 57.0214; // Assume reduction-alkylation
    case 'E':
        return 129.04259;
    case 'Q':
        return 128.05858;
    case 'G':
        return 57.02146;
    case 'H':
        return 137.05891;
    case 'I':
        return 113.08406;
    case 'L':
        return 113.08406;
    case 'K':
        return 128.09496;
    case 'M':
        return 131.04049;
    case 'F':
        return 147.06841;
    case 'P':
        return 97.05276;
    case 'S':
        return 87.03203;
    case 'T':
        return 101.04768;
    case 'W':
        return 186.07931;
    case 'Y':
        return 163.06333;
    case 'V':
        return 99.06841;
    default:
        throw std::invalid_argument("Unknown amino acid");
    }

}

std::vector<std::vector<char>> run_gluc_digest(const std::vector<char> &sequence)
{
    // Split on 'E' if present
    std::vector<std::vector<char>> results;
    auto current = sequence.begin();
    auto it = current;
    while ((it = std::find(current, sequence.end(), 'E')) != sequence.end()) {
        std::vector<char> subsequence(current, it + 1);
        current = it + 1;
        results.push_back(std::move(subsequence));
    }
    results.push_back(std::move(std::vector<char>(current, sequence.end())));
    return results;
}

void fragment_sequence(const std::vector<char> &sequence, std::vector<double> &fragment_list)
{
    // Run forwards; get B ions
    double mass = 0;
    for (int i = 0; i < sequence.size(); i++) {
        mass += get_amino_acid_mass(sequence[i]);
        // handle terminus; not necessary on N terminus in neutral state?
        fragment_list.push_back(mass);
    }
    // Run backwards; get Y ions
    mass = 0;
    for (int i = (int)sequence.size() - 1; i >= 0; i--) {
        mass += get_amino_acid_mass(sequence[i]);
        // handle C terminus
        mass += 18.01088;
        fragment_list.push_back(mass);
    }
}

int search_sequence(const std::vector<char> &sequence, const std::vector<double> &mass_list, double tolerance, bool gluc_digest, std::vector<std::vector<char>> &searched_sequences)
{
    std::vector<std::vector<char>> sequences;
    int match_count = 0;

    // Check if we want to do a GluC digest first
    if (gluc_digest) {
        sequences = run_gluc_digest(sequence);
    } else {
        sequences.push_back(sequence);
    }
    searched_sequences = sequences;

    // Catch exceptions for unknown amino acids. FASTA format supports X for unknown, B/Z for ambiguous, etc. Just skip those...
    for (const auto &seq : sequences) {
        try {
            std::vector<double> fragments;
            fragment_sequence(seq, fragments);

            // Return true if all of the target masses are in the fragment list (false otherwise)
            bool all_found = true;
            for (auto mass : mass_list) {
                bool mass_found = false;
                for (auto f : fragments) {
                    if (std::abs(mass - f) < tolerance) {
                        mass_found = true;
                        break;
                    }
                }
                if (!mass_found) {
                    all_found = false;
                    break;
                }
            }
            if (all_found) ++match_count;
        } catch (const std::invalid_argument) {
            // We could print a warning, but there are too many to make that practical.
        }
    }

    return match_count;
}

// Argument struct for each thread
struct helper_thread_args_struct {
    helper_thread_args_struct(std::vector<Protein> &sequences, const std::vector<double> &mass_list, bool gluc_digest) : sequences(sequences), mass_list(mass_list), gluc_digest(gluc_digest)
    {
        start_index = -1;
        stop_index = -1;
        tolerance = 0;
    }

    std::vector<Protein> &sequences;
    const std::vector<double> &mass_list;
    bool gluc_digest;
    double tolerance;
    int start_index;
    int stop_index;
    Results result;
};


// Helper thread to run calculations
int SearchThreadProc(struct helper_thread_args_struct *p_args)
{
    int n_matches = 0;
    int result_count;
    int n_sequences = p_args->stop_index - p_args->start_index + 1;
    p_args->result.matches.reserve(n_sequences);
    p_args->result.searched_sequence_lengths.reserve(n_sequences);
    p_args->result.searched_digest_lengths.reserve(n_sequences);
    p_args->result.digests_per_sequence.reserve(n_sequences);
    for (int i = p_args->start_index; i <= p_args->stop_index; i++) {
        if (p_args->sequences[i].sequence_valid()) {
            p_args->result.n_searched_sequences++;
            if (result_count = search_sequence(p_args->sequences[i].sequence, p_args->mass_list, p_args->tolerance, p_args->gluc_digest, p_args->sequences[i].digested_sequences)) {
                p_args->result.matches.push_back(p_args->sequences[i]);
                p_args->result.n_matched_sequences += result_count;
            }
        } else {
            p_args->result.n_skipped_sequences++;
        }
        p_args->result.n_digest_sequences += (int)p_args->sequences[i].digested_sequences.size();
        p_args->result.searched_sequence_lengths.push_back((int)p_args->sequences[i].sequence.size());
        p_args->result.digests_per_sequence.push_back((int)p_args->sequences[i].digested_sequences.size());
        std::vector<int> digest_sizes;
        digest_sizes.reserve(p_args->sequences[i].digested_sequences.size());
        for (const auto &seq : p_args->sequences[i].digested_sequences) {
            digest_sizes.push_back((int)seq.size());
        }
        p_args->result.searched_digest_lengths.push_back(std::move(digest_sizes));
    }

    return 0;
}

std::vector<Database> read_databases(const Configuration &config)
{
    std::vector<Database> databases;
    for (auto path : config.databases) {
        Database database = Database(path);
        databases.push_back(std::move(database));
    }
    return databases;
}

Results search_fragments(const Configuration &config, std::vector<Database> &databases)
{
    int n_threads = config.num_search_threads;
    // If unspecified, query how many hardware threads we can use at once
    if (n_threads == 0) {
        n_threads = std::thread::hardware_concurrency();
    }
    // In case hardware_concurrency() isn't supported, default to single-threaded
    if (n_threads == 0) n_threads = 1;

    // Allocate memory to store the input information for each thread + its results
    std::vector<struct helper_thread_args_struct> thread_args;
    std::vector<Results> database_results;
    database_results.reserve(databases.size());
    for (auto &db : databases) {
        int sequences_seen = (int)db.sequences.size();
        // Setup how many sequences each thread needs to search
        int n_sequences_per_thread = sequences_seen / n_threads;
        int excess = sequences_seen - (n_sequences_per_thread * n_threads);

        for (int i = 0; i < n_threads; i++) {
            // Setup input arguments for each thread
            struct helper_thread_args_struct args(db.sequences, config.target_masses, config.gluc_digest);
            args.tolerance = config.mass_tolerance;
            args.start_index = i * n_sequences_per_thread;
            args.stop_index = (i + 1) * n_sequences_per_thread - 1;
            if (i + 1 == n_threads) args.stop_index += excess;
            thread_args.push_back(args);
        }

        // Create a separate thread for each hardware thread
        std::vector<std::thread> threads;
        threads.resize(n_threads);
        for (int i = 0; i < n_threads; i++) {
            threads[i] = std::thread(SearchThreadProc, &thread_args[i]);
        }

        // Wait for everything to finish
        for (int i = 0; i < n_threads; i++) {
            threads[i].join();
        }

        // All the threads are done now, add up the results
        std::vector<Results> thread_results;
        thread_results.reserve(n_threads);
        for (int i = 0; i < n_threads; i++) {
            thread_results.push_back(std::move(thread_args[i].result));
        }
        database_results.push_back(Results::Combine(thread_results));
    }
    return Results::Combine(database_results);
}

void write_results(const Configuration &config, const Results &results, FILE *output_file)
{
    int current_seq = 0;
    for (const auto &protein : results.matches) {
        ++current_seq;
        fprintf(output_file, "%d: %s [%s]\n", current_seq, protein.description.c_str(), protein.source_database.c_str());
        fputc('\t', output_file);
        for (const auto c : protein.sequence) {
            fputc(c, output_file);
        }
        fputc('\n', output_file);
        if (config.gluc_digest) {
            fprintf(output_file, "\t\tGluC fragments:\n");
            for (const auto &fragment : protein.digested_sequences) {
                fputs("\t\t", output_file);
                for (const auto c : fragment) {
                    fputc(c, output_file);
                }
                fputc('\n', output_file);
            }
        }
        fputc('\n', output_file);
    }
}

Results run_fragment_search(const Configuration &config, FILE *output_file)
{
    // Open up the database files
    auto start_database_reading = std::chrono::high_resolution_clock::now();
    auto databases = read_databases(config);
    auto finish_database_reading = std::chrono::high_resolution_clock::now();

    // Run the database search
    auto start_fragment_search = std::chrono::high_resolution_clock::now();
    auto results = search_fragments(config, databases);
    auto finish_fragment_search = std::chrono::high_resolution_clock::now();

    // Write the results to disk
    auto start_writing_results = std::chrono::high_resolution_clock::now();
    write_results(config, results, output_file);
    auto finish_writing_results = std::chrono::high_resolution_clock::now();

    auto file_reading_time = finish_database_reading - start_database_reading;
    auto searching_time = finish_fragment_search - start_fragment_search;
    auto writing_time = finish_writing_results - start_writing_results;
    auto total_time = finish_writing_results - start_database_reading;

    long long file_millis = std::chrono::duration_cast<std::chrono::milliseconds>(file_reading_time).count();
    long long searching_millis = std::chrono::duration_cast<std::chrono::milliseconds>(searching_time).count();
    long long writing_millis = std::chrono::duration_cast<std::chrono::milliseconds>(writing_time).count();
    long long total_millis = std::chrono::duration_cast<std::chrono::milliseconds>(total_time).count();
    fprintf(stderr, "Elapsed time: %lld ms reading, %lld ms searching, %lld ms writing (%lld ms total).\n", file_millis, searching_millis, writing_millis, total_millis);

    return std::move(results);
}
