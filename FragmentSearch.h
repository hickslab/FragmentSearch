#ifndef FRAGMENT_SEARCH_H
#define FRAGMENT_SEARCH_H

#include <string>
#include <vector>
#include "Configuration.h"
#include "Database.h"
#include "Results.h"

// Core header file with API definitions
// Main execution function
Results run_fragment_search(const Configuration &config, FILE *output_file);

std::vector<Database> read_databases(const Configuration &config);
Results search_fragments(const Configuration &config, std::vector<Database> &databases);
void write_results(const Configuration &config, const Results &results, FILE *output_file);

#endif
