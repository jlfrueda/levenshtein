// g++ -Ofast -ffast-math -ftree-vectorize op3.cc

#include <algorithm>
#include <cassert>
#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>

int levenshtein(int m, const char* x, int n, const char* y)
{
    if (m < n) {
        std::swap(m, n);
        std::swap(x, y);
    }
    if (0 == n) {
        return m;
    } else {
        std::vector<int> buffer(3 * n);
        int* zp = buffer.data();
        int* z = zp + n;
        int* w = z + n;

        z[0] = 1;

        for (int i = 0; i < m + n - 1; ++i) {
            w[0] = std::min(1 + z[0], i + (x[i] == y[0] ? 0 : 1));
            int j0 = std::max(1, i - m + 1);
            int j1 = std::min(n, i + 1);
            for (int j = j0; j < j1; ++j) {
                int c1 = 1 + z[j - 1];
                int c2 = 1 + z[j];
                int c3 = zp[j - 1] + (x[i - j] == y[j] ? 0 : 1);
                w[j] = std::min(std::min(c1, c2), c3);
            }
            if (i + 1 < n) {
                w[i + 1] = i + 2;
            }
            int* tmp = zp;
            zp = z;
            z = w;
            w = tmp;
        }

        return z[n - 1];
    }
}

const std::string readFasta(const char* const path)
{
    std::ifstream stream(path, std::ios_base::in | std::ios_base::binary);
    std::string sequence;
    std::getline(stream, sequence);
    std::getline(stream, sequence);
    return sequence;
}

int main(int argc, char* argv[])
{
    assert(3 == argc);
    const auto sequence1 = argv[1];
    const auto sequence2 = argv[2];

    double cost = levenshtein(strlen(sequence1), sequence1, strlen(sequence2), sequence2);
    std::cout << "Cost = " << cost << std::endl;

    return EXIT_SUCCESS;
}