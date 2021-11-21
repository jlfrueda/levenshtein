// 2021 JLFR
// experiment with levenshtein distance on FASTA sequences on GPU
// result: not fast enough
// next: use more than one SM
// almost 75% slower on V100 wrt GTX1080!

#include <algorithm>
#include <cassert>
#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>

//==============================================================================
// pre: m >= n, len(buffer) >= 3 * n
//==============================================================================

__global__ void levenshteinKernel(int m, const uint8_t* const x, int n, const uint8_t* const y, int32_t* w, int32_t* z, int32_t* zp)
{
    // so, each thread should be responsible of ceil(n / numThreads) elements
    int numThreads = blockDim.x;
    int numSteps = n / numThreads + (n % numThreads != 0);

    int32_t* result = z;
    result[0] = 1.0f;

    for (int i = 0; i < m + n - 1; ++i) {

        int j0 = max(0, i - m + 1);
        int j1 = min(n - 1, i); // included, [j0, j1]

        for (int step = 0; step < numSteps; ++step) {
            int j = step * blockDim.x + threadIdx.x;
            if (j0 <= j && j <= j1) {
                int32_t v = (0 == j || i == j ? i : zp[j - 1]) + (x[i - j] == y[j] ? 0 : 1);
                v = 0 == j ? v : min(v, 1 + z[j - 1]);
                v = i == j ? v : min(v, 1 + z[j]);
                w[j] = v;
            }
        }
        int32_t* tmp = zp;
        zp = z;
        z = w;
        w = tmp;
        __syncthreads();
    }

    if (0 == threadIdx.x) {
        result[0] = z[n - 1];
    }
}

/*
const std::string readFasta(const char* const path)
{
    std::ifstream stream(path, std::ios_base::in | std::ios_base::binary);
    std::string sequence;
    std::getline(stream, sequence);
    std::getline(stream, sequence);
    return sequence;
}
*/

void runExperiment(const uint8_t* pX, const uint8_t* pY)
{
    int m = strlen(reinterpret_cast<const char*>(pX));
    int n = strlen(reinterpret_cast<const char*>(pY));

    if (m < n) {
        std::swap(m, n);
        std::swap(pX, pY);
    }

    size_t xPitch = 0;
    uint8_t* x = nullptr;
    cudaMallocPitch(&x, &xPitch, m, 2);

    uint8_t* y = x + xPitch;
    cudaMemcpy(x, pX, m * sizeof(uint8_t), cudaMemcpyHostToDevice);
    cudaMemcpy(y, pY, n * sizeof(uint8_t), cudaMemcpyHostToDevice);

    size_t bufferPitch = 0;
    int32_t* buffer = nullptr;
    cudaMallocPitch(&buffer, &bufferPitch, n * sizeof(int32_t), 3);

    int32_t* const w = buffer;
    int32_t* const z = reinterpret_cast<int32_t*>(reinterpret_cast<uint8_t*>(w) + bufferPitch);
    int32_t* const zp = reinterpret_cast<int32_t*>(reinterpret_cast<uint8_t*>(z) + bufferPitch);

    levenshteinKernel<<<1, 1024>>>(m, x, n, y, w, z, zp);

    cudaDeviceSynchronize();

    int32_t result = 0;
    cudaMemcpy(&result, z, sizeof(int32_t), cudaMemcpyDeviceToHost);

    cudaFree(buffer);
    cudaFree(x);

    std::cout << "Result: " << result << std::endl;
}

int main(int argc, char* argv[])
{
    assert(3 == argc);
    const auto sequence1 = reinterpret_cast<const uint8_t*>(argv[1]);
    const auto sequence2 = reinterpret_cast<const uint8_t*>(argv[2]);

    runExperiment(sequence1, sequence2);

    return EXIT_SUCCESS;
}
