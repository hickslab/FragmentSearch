#include "Results.h"

Results::Results()
{
    this->n_searched_sequences = 0;
    this->n_matched_sequences = 0;
    this->n_skipped_sequences = 0;
    this->n_digest_sequences = 0;
}

Results Results::Combine(std::vector<Results> results)
{
    Results final_result;
    for (const auto &result : results) {
        final_result.n_matched_sequences += result.n_matched_sequences;
        final_result.n_searched_sequences += result.n_searched_sequences;
        final_result.n_skipped_sequences += result.n_skipped_sequences;
        final_result.n_digest_sequences += result.n_digest_sequences;
        final_result.matches.insert(final_result.matches.end(), result.matches.begin(), result.matches.end());

        final_result.searched_sequence_lengths.insert(final_result.searched_sequence_lengths.end(), result.searched_sequence_lengths.begin(), result.searched_sequence_lengths.end());
        final_result.searched_digest_lengths.insert(final_result.searched_digest_lengths.end(), result.searched_digest_lengths.begin(), result.searched_digest_lengths.end());
        final_result.digests_per_sequence.insert(final_result.digests_per_sequence.end(), result.digests_per_sequence.begin(), result.digests_per_sequence.end());
    }
    return final_result;
}
