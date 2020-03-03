
#define _CRT_SECURE_NO_WARNINGS // Silence unnecessary warnings

#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <algorithm>
#include <numeric>

#include "FragmentSearch.h"
#include "Configuration.h"
#include "Database.h"
#include "Results.h"
#include "Protein.h"


// Print out correct command-line usage
void print_usage(const char *app_path)
{
    fprintf(stderr, "Usage: %s input_file output_file\n", app_path);
}

// Entry point for the application
int main(int argc, char *argv[])
{
    // Check command line is correct
    if (argc < 3) {
        print_usage(argv[0]);
        return 1;
    }

    try {
        // Read in configuration
        Configuration configuration = Configuration(argv[1]);
        if (configuration.databases.size() == 0) {
            fprintf(stderr, "No databases found.\n");
            return 1;
        }

        // Open up output file
        FILE *output_fp = fopen(argv[2], "w");
        if (!output_fp) {
            int err = errno;
            fprintf(stderr, "Failed to open output file '%s' for writing: errno %d.\n", argv[2], err);
            return 1;
        }

        // Run the database search
        auto final_results = run_fragment_search(configuration, output_fp);

        int matching_sequences = final_results.n_matched_sequences;
        int searched_sequences = final_results.n_searched_sequences;
        int skipped_sequences = final_results.n_skipped_sequences;
        int digest_sequences = final_results.n_digest_sequences;
        int total_sequences = searched_sequences + skipped_sequences;

        fputs("Search complete for mass list:", stdout);
        for (const auto mass : configuration.target_masses) {
            fprintf(stdout, " %.2lf", mass);
        }
        fputc('\n', stdout);

        fprintf(stdout, "%d matches, %d searched sequences, %d digests, %d skipped (%d total) (%lf%%).\n",
            matching_sequences, searched_sequences, digest_sequences, skipped_sequences, total_sequences, (double)matching_sequences / digest_sequences * 100);

        // Stats on sequence lengths
        auto minmax_sequence_lengths = std::minmax_element(final_results.searched_sequence_lengths.begin(), final_results.searched_sequence_lengths.end());
        auto min_sequence_length = *minmax_sequence_lengths.first;
        auto max_sequence_length = *minmax_sequence_lengths.second;
        auto average_sequence_length = std::accumulate(final_results.searched_sequence_lengths.begin(), final_results.searched_sequence_lengths.end(), 0.0) / final_results.searched_sequence_lengths.size();

        auto minmax_num_digests = std::minmax_element(final_results.digests_per_sequence.begin(), final_results.digests_per_sequence.end());
        auto min_num_digests = *minmax_num_digests.first;
        auto max_num_digests = *minmax_num_digests.second;
        auto average_num_digests = std::accumulate(final_results.digests_per_sequence.begin(), final_results.digests_per_sequence.end(), 0.0) / final_results.digests_per_sequence.size();

        fprintf(stdout, "Min sequence length = %d, max sequence length = %d, average sequence length = %lf\n", min_sequence_length, max_sequence_length, average_sequence_length);
        fprintf(stdout, "Min num digests = %d, max num digests = %d, average num digests = %lf\n", min_num_digests, max_num_digests, average_num_digests);


    } catch (const std::exception &ex) {
        fprintf(stderr, "Program execution terminated with exception: %s\n", ex.what());
        return 1;
    }

    return 0;
}
