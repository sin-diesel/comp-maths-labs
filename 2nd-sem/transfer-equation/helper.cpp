#include "helper.hpp"

gtype phi(gtype x) {
    return x;
}

gtype psi(gtype t) {
    return t;
}

gtype f(gtype t, gtype x) {
    return t + x;
}

void angle(matrix::matrix_t<gtype>& grid, std::size_t k, std::size_t m, gtype tau, gtype h) {
    std::size_t K = grid.get_rows_number();
    grid[k-1][m] = f((K - 1 - k) * tau, m * h) * tau - (grid[k][m] - grid[k][m-1]) / h * tau + grid[k][m];
}

void krest(matrix::matrix_t<gtype>& grid, std::size_t k, std::size_t m, gtype tau, gtype h) {
    std::size_t K = grid.get_rows_number();
    grid[k-1][m] = f((K - 1 - k) * tau, m * h) * 2 * tau - (grid[k][m+1] - grid[k][m-1]) * tau / h + grid[k+1][m];
}

void initGrid(matrix::matrix_t<gtype>& grid, gtype tau, gtype h) {
    std::size_t K = grid.get_rows_number();
    std::size_t M = grid.get_cols_number();

    for(std::size_t i = 0; i < M; ++i) {
        grid[K-1][i] = phi(i * h);
    }

    for(std::size_t i = 0; i < K; ++i) {
        grid[K - 1 - i][0] = psi(i * tau);
    }
}

void initLayer(matrix::matrix_t<gtype>& grid, gtype tau, gtype h) {

    std::size_t K = grid.get_rows_number();
    std::size_t M = grid.get_cols_number();
    
    for(std::size_t i = 1; i < M; ++i) {
        angle(grid, K - 1, i, tau, h);
    }

}


