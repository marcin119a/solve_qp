#include "csv_utils.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <regex>
#include <stdexcept>

FileFormat detect_file_format(const std::string& content) {
    if (content.find('\t') != std::string::npos) {
        if (content.find("Mutated") != std::string::npos) return FileFormat::MutatedTSV;
        return FileFormat::TSV;
    }
    if (content.find(',') != std::string::npos) {
        if (content.find("Mutated") != std::string::npos) return FileFormat::MutatedCSV;
        return FileFormat::CSV;
    }
    return FileFormat::Unknown;
}

std::pair<FileFormat, char> detect_format(const std::string& content) {
    FileFormat fmt = detect_file_format(content);
    if (fmt == FileFormat::CSV || fmt == FileFormat::MutatedCSV) return {fmt, ','};
    if (fmt == FileFormat::TSV || fmt == FileFormat::MutatedTSV) return {fmt, '\t'};
    return {FileFormat::Unknown, ','};
}

std::vector<std::string> split_line(const std::string& line, char sep) {
    std::stringstream ss(line);
    std::string item;
    std::vector<std::string> result;
    while (std::getline(ss, item, sep)) {
        result.push_back(item);
    }
    return result;
}

std::vector<double> parse_line_to_doubles(const std::string& line, char sep, int skip_columns = 0) {
    std::vector<std::string> parts = split_line(line, sep);
    std::vector<double> result;

    for (size_t i = skip_columns; i < parts.size(); ++i) {
        try {
            result.push_back(std::stod(parts[i]));
        } catch (const std::invalid_argument& e) {
            std::cerr << "Nie można sparsować wartości '" << parts[i] << "' jako liczby w linii: " << line << std::endl;
            result.push_back(0.0);  // lub możesz rzucić wyjątek ponownie
        }
    }

    return result;
}

std::pair<std::vector<std::vector<double>>, std::vector<std::string>>
load_samples_csv(const std::string& path) {
    std::ifstream file(path);
    if (!file.is_open()) throw std::runtime_error("Cannot open samples file: " + path);

    std::stringstream buffer;
    buffer << file.rdbuf();
    std::string content = buffer.str();

    auto [format, sep] = detect_format(content);
    if (format == FileFormat::Unknown)
        throw std::runtime_error("Unknown Format in samples file");

    file.clear();
    file.seekg(0);
    std::string header_line;
    std::getline(file, header_line);
    std::vector<std::string> patient_names = split_line(header_line, sep);

    int skip_columns = 0;
    if (format == FileFormat::TSV || format == FileFormat::MutatedTSV)
        skip_columns = 1;
    else if (format == FileFormat::CSV || format == FileFormat::MutatedCSV)
        skip_columns = 2;

    std::vector<std::vector<double>> data;
    std::string line;
    while (std::getline(file, line)) {
        std::vector<double> row = parse_line_to_doubles(line, sep, skip_columns);
        if (!row.empty()) data.push_back(row);
    }

    if (format == FileFormat::TSV || format == FileFormat::MutatedTSV) {
        patient_names.erase(patient_names.begin());  // usuń pierwszy nagłówek (np. "Samples")
    } else if (format == FileFormat::CSV || format == FileFormat::MutatedCSV) {
        patient_names.erase(patient_names.begin(), patient_names.begin() + 2);
    }

    return {data, patient_names};
}

std::pair<std::vector<std::vector<double>>, std::vector<std::string>>
load_signatures_csv(const std::string& path) {
    std::ifstream file(path);
    if (!file.is_open()) throw std::runtime_error("Cannot open signatures file: " + path);

    std::stringstream buffer;
    buffer << file.rdbuf();
    std::string content = buffer.str();

    auto [format, sep] = detect_format(content);
    if (format == FileFormat::Unknown)
        throw std::runtime_error("Unknown Format in signatures file");

    file.clear();
    file.seekg(0);
    std::string header_line;
    std::getline(file, header_line);
    std::vector<std::string> sig_names = split_line(header_line, sep);
    sig_names.erase(sig_names.begin());  // usuń pierwszy nagłówek (np. "Type")
    sig_names.insert(sig_names.begin(), "Samples");

    std::vector<std::vector<double>> data;
    std::string line;
    while (std::getline(file, line)) {
        std::vector<double> row = parse_line_to_doubles(line, sep, 1);  // pomiń pierwszą kolumnę
        if (!row.empty()) data.push_back(row);
    }

    return {data, sig_names};
}
