#ifndef PROTEIN_H
#define PROTEIN_H

#include <string>
#include <vector>

class Protein
{

public:
    std::string source_database;
    std::string description;
    std::vector<char> sequence;
    std::vector<std::vector<char>> digested_sequences;

public:
    bool sequence_valid() const;
};

#endif // PROTEIN_H
