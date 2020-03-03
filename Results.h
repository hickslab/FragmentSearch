#ifndef RESULTS_H
#define RESULTS_H

#include <vector>

#include "Protein.h"

class Results
{
public:
    int n_searched_sequences;
    int n_skipped_sequences;
    int n_matched_sequences;
    int n_digest_sequences;
    std::vector<Protein> matches;
    std::vector<int> searched_sequence_lengths;
    std::vector<std::vector<int>> searched_digest_lengths;
    std::vector<int> digests_per_sequence;

public:
    Results();

    static Results Combine(std::vector<Results> results);

};

#endif // RESULTS_H
