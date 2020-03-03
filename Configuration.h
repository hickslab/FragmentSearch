#ifndef CONFIGURATION_H
#define CONFIGURATION_H

#include <vector>
#include <string>

class Configuration
{
    // Configuration variables
public:
    std::vector<std::string> databases;
    std::vector<double> target_masses;
    double mass_tolerance;
    bool gluc_digest;

    int num_search_threads;
    bool read_database_multithreaded;

public:
    // Constructor: from file
    Configuration(const char *filename);

private:
    // Constructor: default configuration
    Configuration();
};

#endif //CONFIGURATION_H
