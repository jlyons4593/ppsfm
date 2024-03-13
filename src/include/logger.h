

#include <iostream>
#include <string>

#pragma once
namespace Logger {

     // Logs a section header with more spacing to signify a major section
    static void logSection(const std::string& sectionName) {
        std::cout << "\n\n========= " << sectionName << " =========\n" << std::endl;
    }

    // Logs a subsection header with less spacing, indicating it's part of a larger section
    static void logSubsection(const std::string& subsectionName) {
        std::cout << "\n--- " << subsectionName << " ---\n" << std::endl;
    }

    // Logs an error message
    static void logError(const std::string& errorMessage) {
        std::cout << "ERROR: " << errorMessage << std::endl;
    }
};

