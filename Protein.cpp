#include "Protein.h"

#include <cstring>

bool Protein::sequence_valid() const
{
    const char valid_amino_acids[] = "ARNDCcEQGHILKMFPSTWYV";

    for (auto c : this->sequence) {
        if (!strchr(valid_amino_acids, c)) return false;
    }
    return true;
}
