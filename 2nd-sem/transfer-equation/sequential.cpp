
#include <vector>
#include <iostream>
#include <chrono>

#include "matrix.hpp"
#include "helper.hpp"

int main(int argc, char** argv) {

    std::size_t K, M;
    gtype tau, h; 

    K = std::atoi(argv[1]);
    M = std::atoi(argv[2]);
    tau = std::atof(argv[3]);
    h = std::atof(argv[4]);

#ifdef TIME
    auto start = std::chrono::steady_clock::now();
#endif
    matrix::matrix_t<gtype> grid(K, M);

    initGrid(grid, tau, h);
    initLayer(grid, tau, h);
        
    for(std::size_t i = 1; i < K - 1; ++i) {
        for(std::size_t j = 1; j < M - 1; ++j) {
            krest(grid, K - 1 - i, j, tau, h);
        }
        angle(grid, K - 1 - i, M - 1, tau, h);
    }

#ifdef TIME
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";
#endif

#ifdef DUMP    
    for(std::size_t i = K - 2; i > 0; --i) {
        for(std::size_t j = 1; j < M - 1; ++j) {
            std::cout << std::setw(10) << grid[i][j];
        }
        std::cout << std::endl;
    }
#endif
}
