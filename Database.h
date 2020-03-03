#ifndef DATABASE_H
#define DATABASE_H

#include <string>
#include <vector>

#include "Protein.h"

class Database
{
public:
    std::string source_path;
    std::vector<Protein> sequences;

public:
    Database(std::string path);

};

#endif // DATABASE_H
