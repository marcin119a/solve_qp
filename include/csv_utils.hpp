#pragma once
#include <string>
#include <vector>

enum class FileFormat {
    CSV,
    TSV,
    MutatedCSV,
    MutatedTSV,
    Unknown
};

std::pair<FileFormat, char> detect_format(const std::string& content);

std::pair<std::vector<std::vector<double>>, std::vector<std::string>>
load_samples_csv(const std::string& path);

std::pair<std::vector<std::vector<double>>, std::vector<std::string>>
load_signatures_csv(const std::string& path);