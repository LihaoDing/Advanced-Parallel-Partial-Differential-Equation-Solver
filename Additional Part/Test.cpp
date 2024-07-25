#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>

void read_and_print(const std::string& filename, int imax, int jmax, const std::string& filename2) {
    std::ifstream file(filename, std::ios::binary); // Open in binary mode
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }

    std::vector<double> data(imax * jmax);
    file.read(reinterpret_cast<char*>(data.data()), sizeof(double) * data.size());
    file.close();

    std::ifstream file2(filename2);  // Open in text mode
    if (!file2.is_open()) {
        std::cerr << "Failed to open file: " << filename2 << std::endl;
        return;
    }

    std::vector<double> data2;
    std::string line;
    while (std::getline(file2, line)) {
        std::stringstream ss(line);
        double value;
        while (ss >> value) {
            data2.push_back(value);
            if (ss.peek() == ' ') {
                ss.ignore();
            }
        }
    }
    file2.close();

    int temp = 1;
    for (int i = 0; i < imax; i++) {
        for (int j = 0; j < jmax; j++) {
            if (i * jmax + j < data2.size()) {
                if (std::abs(data[i * jmax + j] - data2[i * jmax + j]) < 1e-6) {
                    int a = 1;
                } else {
                    std::cout << "Found a difference at (" << i << ", " << j << "): " << data[i * jmax + j] << " != " << data2[i * jmax + j] << std::endl;
                    temp = 0;
                }
            }
        }
    }
    if (temp == 1) {
        std::cout << "All values are equal" << std::endl;
    }
}

int main() {
    std::cout << std::endl;
    int imax = 301;
    int jmax = 301;

    for (int i = 0; i <= 300; ++i) {
        std::cout << "\nComparing output_C2_" << i << ".dat" << std::endl;
        std::string file1 = "../serial_per/output_C2_" + std::to_string(i) + ".dat";
        std::string file2 = "../result_per/output_C2_" + std::to_string(i) + ".dat";
        read_and_print(file2, imax, jmax, file1);
    }
    std::cout << std::endl;
    return 0;
}
