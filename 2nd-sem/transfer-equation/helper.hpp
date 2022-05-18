#pragma once

#include "matrix.hpp"

using gtype = double;

gtype phi(gtype x);
gtype psi(gtype t);
gtype f(gtype t, gtype x);

void angle(matrix::matrix_t<gtype>& grid, std::size_t k, std::size_t m, gtype tau, gtype h);
void krest(matrix::matrix_t<gtype>& grid, std::size_t k, std::size_t m, gtype tau, gtype h);

void initGrid(matrix::matrix_t<gtype>& grid, gtype tau, gtype h);
void initLayer(matrix::matrix_t<gtype>& grid, gtype tau, gtype h);
