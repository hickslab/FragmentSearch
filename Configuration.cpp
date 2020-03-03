#include "Configuration.h"

#include <cstdlib>
#include <cstring>
#include <fstream>
#include <algorithm>

#ifdef _MSC_VER
#define strcmpi(a, b) _strcmpi(a, b)
#define strncmpi(a, b, n) _strnicmp(a, b, n)
#else
#define strcmpi(a, b) strcasecmp(a, b)
#define strncmpi(a, b, n) strncasecmp(a, b, n)
#endif

// Default configuration
Configuration::Configuration()
{
    this->mass_tolerance = 0;
    this->gluc_digest = true;
    this->num_search_threads = 0;
    this->read_database_multithreaded = false;
}

Configuration::Configuration(const char *filename)
{
    Configuration();

    std::ifstream input = std::ifstream(filename);
    if (!input) {
        fprintf(stderr, "Failed to open '%s' for reading.\n", filename);
        return;
    }
    std::string buffer;
    while (std::getline(input, buffer)) {

        // Filter out any '\r' characters (important if a Linux system is reading a Windows-created input file)
        buffer.erase(std::remove(buffer.begin(), buffer.end(), '\r'), buffer.end());

        // Trim off whitespace
        const char *ptr = buffer.c_str();
        while (*ptr && isblank(*ptr)) ++ptr;

        // Skip empty lines
        if (*ptr == '\0') continue;

        // Skip comment lines
        if (*ptr == '#') continue;

        // Check if the '=' is present (warning otherwise)
        const char *equals = strchr(ptr, '=');
        if (!equals) {
            fprintf(stderr, "Invalid configuration line: '%s'\n", buffer.c_str());
            continue;
        }

        // Split into key, value
        std::string key = std::string(ptr, equals - 1);
        const char *value = equals + 1;

        // Trim whitespace from value
        while (*value && isblank(*value)) ++value;

        // Check against possible configuration (case insensitive)
        if (strcmpi(key.c_str(), "database") == 0) {
            // Trim off ending whitespace
            size_t len = strlen(value);
            const char *end_value = value + len - 1;
            while (*end_value && isblank(*end_value)) --end_value;
            // Get rid of quotes
            if (*value == '"') ++value;
            if (*end_value == '"') --end_value;
            if (value >= end_value) {
                fprintf(stderr, "Invalid database value: '%s'\n", value - 1);
                continue;
            }
            std::string value_str = std::string(value, end_value + 1);
            this->databases.push_back(value_str);
        }
        else if (strcmpi(key.c_str(), "mass_tolerance") == 0) {
            char *endptr;
            double value_double = strtod(value, &endptr);
            if (endptr == value) {
                fprintf(stderr, "Invalid double value for mass tolerance: '%s'\n", value);
                continue;
            }
            this->mass_tolerance = value_double;
        }
        else if (strcmpi(key.c_str(), "target_masses") == 0) {
            // Comma-separated list of doubles
            char *endptr;
            double value_double;
            do {
                value_double = strtod(value, &endptr);
                if (endptr == value) {
                    fprintf(stderr, "Invalid double value for target mass list: '%s'\n", value);
                    break;
                }
                value = endptr;
                this->target_masses.push_back(value_double);
            } while (*value && *value != '\n');
        }
        else if (strcmpi(key.c_str(), "num_search_threads") == 0) {
            char *endptr;
            int value_int = strtol(value, &endptr, 0);
            if (endptr == value) {
                fprintf(stderr, "Invalid int value for num_search_threads: '%s'\n", value);
                continue;
            }
            this->num_search_threads = value_int;
        }
        else if (strcmpi(key.c_str(), "read_database_multithreaded") == 0) {
            bool value_bool;
            if (strncmpi(value, "true", sizeof("true") - 1) == 0) {
                value_bool = true;
            }
            else if (strncmpi(value, "false", sizeof("false") - 1) == 0) {
                value_bool = false;
            }
            else {
                fprintf(stderr, "Invalid bool value for read_database_multithreaded: '%s'\n", value);
                continue;
            }
            this->read_database_multithreaded = value_bool;
        }
        else if (strcmpi(key.c_str(), "gluc_digest") == 0) {
            bool value_bool;
            if (strncmpi(value, "true", sizeof("true") - 1) == 0) {
                value_bool = true;
            } else if (strncmpi(value, "false", sizeof("false") - 1) == 0) {
                value_bool = false;
            }
            else {
                fprintf(stderr, "Invalid bool value for gluc_digest: '%s'\n", value);
                continue;
            }
            this->gluc_digest = value_bool;
        }
        else {
            fprintf(stderr, "Unrecognized option '%s'\n", key.c_str());
            continue;
        }
    }
}
