// g++ -Ofast -ffast-math -ftree-vectorize kernelized.cc

#include <algorithm>
#include <cassert>
#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>

void kernel(int i, int j, int m, const char* const x, int n, const char* const y, int* w, int* z, int* zp)
{
    int j0 = std::max(0, i - m + 1);
    int j1 = std::min(n - 1, i); // included, [j0, j1]

    if (j0 <= j && j <= j1) {
        int v = (0 == j || i == j ? i : zp[j - 1]) + (x[i - j] == y[j] ? 0 : 1);
        v = 0 == j ? v : std::min(v, 1 + z[j - 1]);
        v = i == j ? v : std::min(v, 1 + z[j]);
        w[j] = v;
    }
}

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
            int j0 = std::max(0, i - m + 1);
            int j1 = std::min(n - 1, i); // included, [j0, j1]
            for (int j = j0; j <= j1; ++j) {
                kernel(i, j, m, x, n, y, w, z, zp);
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