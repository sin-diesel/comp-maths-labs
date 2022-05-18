#pragma once

#include <utility>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <optional>
#include <iomanip>

#include "matrix_buffer.hpp"

namespace matrix {

const long double tolerance = 1e-5;

template<typename T, typename U>
bool equal(const T& lhs, const U& rhs) {
    return (std::abs(lhs - rhs) < tolerance); 
}

template<typename T>
class matrix_t final : private matrix_buff_t<T> {

    using matrix_buff_t<T>::swap_buffers;
    using matrix_buff_t<T>::construct_at;
    using matrix_buff_t<T>::at;
    using matrix_buff_t<T>::get_ptr;
    using matrix_buff_t<T>::resize_buffer;

public:

    using matrix_buff_t<T>::get_rows_number;
    using matrix_buff_t<T>::get_cols_number;
    using matrix_buff_t<T>::get_elements_number;

    matrix_t() : matrix_buff_t<T>(0u, 0u) {}
    matrix_t(const matrix_t& rhs);
    matrix_t(std::size_t rows, std::size_t cols, T val = T{});
    matrix_t(const std::initializer_list<std::initializer_list<T>>& init);
    matrix_t& operator=(const matrix_t& rhs);
    ~matrix_t() = default;

    void resize(std::size_t rows, std::size_t cols);

    /*  
        insert colomn on idx position

         function contract:

            1) *this : p * q
               m     : p * 1
            2) 0 <= idx <= m.get_cols_number()
            
    */
    void insert_col(std::size_t idx, const matrix_t<T>& m);

    void swap_rows(std::size_t lhs_idx, std::size_t rhs_idx);
    void swap_cols(std::size_t lhs_idx, std::size_t rhs_idx);
    std::size_t max_abs_col_elem(std::size_t idx, std::size_t start, std::size_t end) const;
    T det() const;

    class proxy_row_t {
        public:
            proxy_row_t(T* pointer) : row_(pointer) {}
            T& operator[](std::size_t idx) {return row_[idx];}
            const T& operator[](std::size_t idx) const {return row_[idx];}
        private:
            T* row_;
    };
        
    bool operator==(const matrix_t<T>& rhs) const;
    bool operator!=(const matrix_t<T>& rhs) const;

    proxy_row_t operator[](std::size_t idx);
    const proxy_row_t operator[](std::size_t idx) const;
};

template<typename T>
std::ostream& operator<<(std::ostream& out,  const matrix::matrix_t<T>& m);

template<typename T>
std::istream& operator>>(std::istream& in,  matrix::matrix_t<T>& m);

/*
    ***multiplication of 2 matrices***

    function contract:

        1)  lhs matrix: m * r
            rhs matrix: r * n
      
*/
template<typename T>
matrix_t<T> multiplication(const matrix_t<T>& lhs, const matrix_t<T>& rhs);

/* 
    ***solve linear system***

    fuction contract:

        1) right matrix must be m * 1
           left matrix must be m * n

        2) return value: partial solution + fundamental matrix 
           partial solution is a matrix m * 1

*/
template<typename T>
std::pair<matrix_t<T>, matrix_t<T>> solve_linear_system(const matrix_t<T>& left, const matrix_t<T>& right);


template<typename T>
void gauss_straight(matrix_t<T>& m);

template<typename T>
void gauss_reverse(matrix_t<T>& m);

/* ------------------------------------------------------------------
                        IMPLEMENTATION
---------------------------------------------------------------------*/

template<typename T>
matrix_t<T>::matrix_t(const matrix_t& rhs) : matrix_buff_t<T>(rhs.get_rows_number(), rhs.get_cols_number()) {
    for(std::size_t i = 0, maxi = get_rows_number(); i < maxi; ++i) {
        for(std::size_t j = 0, maxj = get_cols_number(); j < maxj; ++j) {
            construct_at(i, j, rhs[i][j]);
        }
    }
}

template<typename T>
matrix_t<T>::matrix_t(std::size_t rows, std::size_t cols, T val /* = T() */) : matrix_buff_t<T>(rows, cols) {
    for(std::size_t i = 0, maxi = rows; i < maxi; ++i) {
        for(std::size_t j = 0, maxj = cols; j < maxj; ++j) {
            construct_at(i, j, val);
        }
    }
}

template<typename T>
matrix_t<T>::matrix_t(const std::initializer_list<std::initializer_list<T>>& init) : matrix_buff_t<T>(0u, 0u) {
    std::size_t m = init.size();
    std::size_t n = 0u;

    for(const auto& row : init) {
        n = std::max(n, row.size());
    }

    std::size_t current_n = 0u;
    std::size_t current_m = 0u;
    resize_buffer(m, n);
    for(const auto& row : init) {
        for(const auto& elem : row) {
            construct_at(current_m, current_n, elem);
            ++current_n;
        }

        for(std::size_t i = current_n; i < n; ++i) {
            construct_at(current_m, i, T());
        }

        ++current_m; current_n = 0;
    }
}

template<typename T>
void matrix_t<T>::resize(std::size_t rows, std::size_t cols) {
    matrix_t<T> tmp;
    tmp.resize_buffer(rows, cols);

    for(std::size_t i = 0, maxi = rows; i < maxi; ++i) {
        for(std::size_t j = 0, maxj = cols; j < maxj; ++j) {
            T val = T();
            if((i < get_rows_number()) && (j < get_cols_number())) {
               val = at(i, j); 
            }

            tmp.construct_at(i, j, val);
        }
    }

    swap_buffers(tmp);
}    

template<typename T>
void matrix_t<T>::insert_col(std::size_t idx, const matrix_t<T>& m) {
    matrix_t<T> tmp;
    tmp.resize_buffer(get_rows_number(), get_cols_number() + 1);

    for(std::size_t i = 0, maxi = tmp.get_rows_number(); i < maxi; ++i) {
        for(std::size_t j = 0, maxj = tmp.get_cols_number(); j < maxj; ++j) {
            if(j < idx) {
                tmp.construct_at(i, j, at(i, j));
            } else if(j == idx) {
                tmp.construct_at(i, j, m[i][0]);
            } else {
                tmp.construct_at(i, j, m[i][j-1]);
            }
        }
    }

    swap_buffers(tmp);
}

template<typename T>
matrix_t<T>& matrix_t<T>::operator=(const matrix_t& rhs) {
    matrix_t tmp(rhs);
    swap_buffers(tmp);
    return *this;
}

template<typename T>
void matrix::matrix_t<T>::swap_rows(std::size_t lhs_idx, std::size_t rhs_idx) {
    if(lhs_idx == rhs_idx) {
        return;
    }

    for(std::size_t i = 0, max = get_cols_number(); i < max; ++i) {
        T tmp = at(lhs_idx, i);
        at(lhs_idx, i) = at(rhs_idx, i);
        at(rhs_idx, i) = tmp;
    }

    return;
}

template<typename T>
void matrix::matrix_t<T>::swap_cols(std::size_t lhs_idx, std::size_t rhs_idx) {
    if(lhs_idx == rhs_idx) {
        return;
    }

    std::size_t max_col = get_cols_number();

    for(std::size_t i = 0, max = get_rows_number(); i < max; ++i) {
        T tmp = at(i, lhs_idx);
        at(i, lhs_idx) = at(i, rhs_idx);
        at(i, rhs_idx) = tmp;
    }

    return;
}

template<typename T>
std::size_t matrix::matrix_t<T>::max_abs_col_elem(std::size_t idx, std::size_t start, std::size_t end) const {
    std::size_t ret = start;

    std::size_t max_col = get_cols_number();
    for(std::size_t i = start + 1; i < end; ++i) {
        if(std::abs(at(i, idx)) > std::abs(at(ret, idx))) {
            ret = i;
        }
    }

    return ret;
}

template<typename T>
T matrix::matrix_t<T>::det() const {
    long double ret = 1.0;

    matrix_t<double> copy = *this;

    std::size_t current_m = 0;
    std::size_t current_n = 0;

    std::size_t max_m = copy.get_rows_number();
    std::size_t max_n = copy.get_cols_number();

    while((current_m < max_m) && (current_n < max_n)) {
        std::size_t max_elem = copy.max_abs_col_elem(current_n, current_m, max_m);
        if(copy[max_elem][current_n] == 0) {
            ++current_n;
            continue;
        }

        if(max_elem != current_m) {
            ret *= -1;
            copy.swap_rows(max_elem, current_m);
        }

        for(std::size_t i = current_m + 1; i < max_m; ++i) {
            double f = copy[i][current_n] / copy[current_m][current_n];
            copy[i][current_n] = 0.0;
            for (std::size_t j = current_n + 1; j < max_n; ++j) {
                copy[i][j] = copy[i][j] - f * copy[current_m][j];
            }
        }

        ++current_m;
        ++current_n;
    }

    for(std::size_t i = 0lu, max = copy.get_cols_number(); i < max; ++i) {
        if(std::abs(ret) < tolerance) {
            return 0.0;
        }

        ret *= copy[i][i];
    }

    return ret;
}

template<typename T>
const typename matrix::matrix_t<T>::proxy_row_t matrix::matrix_t<T>::operator[](std::size_t idx) const {
    return const_cast<matrix_t<T>&>(*this)[idx];
}

template<typename T>
typename matrix::matrix_t<T>::proxy_row_t matrix::matrix_t<T>::operator[](std::size_t idx) {
    proxy_row_t ret(get_ptr(idx, 0));
    return ret;
}

template<typename T>
bool matrix::matrix_t<T>::operator==(const matrix_t<T>& rhs) const {
    if((get_cols_number() != rhs.get_cols_number()) || (get_rows_number() != rhs.get_rows_number())) {
        return false;
    }

    for(std::size_t i = 0, maxi = get_rows_number(); i < maxi; ++i) {
        for(std::size_t j = 0, maxj = get_cols_number(); j < maxj; ++j) {
            if(!equal(at(i, j), rhs[i][j])) {
                return false;
            }
        }
    }

    return true;
}

template<typename T>
bool matrix::matrix_t<T>::operator!=(const matrix_t<T>& rhs) const {
    return !(*this == rhs);
}

template<typename T>
std::ostream& matrix::operator<<(std::ostream& out,  const matrix::matrix_t<T>& m) {
    if(m.get_elements_number() == 0) {
        return out;
    }

    for(std::size_t i = 0; i < m.get_rows_number(); ++i) {
        for(std::size_t j = 0; j < m.get_cols_number(); ++j) {
            out << std::setw(8) << m[i][j] << ' ';
        }

        out << std::endl;
    }

    return out;
}

template <typename T>
std::istream& matrix::operator>>(std::istream& in,  matrix::matrix_t<T>& m) {
    if(m.get_elements_number() == 0) {
        return in;
    }

    for(std::size_t i = 0, maxi = m.get_rows_number(); i < maxi; ++i) {
        for(std::size_t j = 0, maxj = m.get_cols_number(); j < maxj; ++j) {
            in >> m[i][j];
        }
    }

    return in;
}

template<typename T>
matrix::matrix_t<T> matrix::multiplication(const matrix_t<T>& lhs, const matrix_t<T>& rhs) {
    std::size_t ret_m = lhs.get_rows_number();
    std::size_t ret_n = rhs.get_cols_number();

    matrix_t<T> ret{ret_m, ret_n};
    for(std::size_t i = 0, maxi = ret_m; i < maxi; ++i) {
        for(std::size_t j = 0, maxj = ret_n; j < maxj; ++j) {
            T elem{};
            for(std::size_t k = 0, maxk = lhs.get_cols_number(); k < maxk; ++k) {
                elem += (lhs[i][k] * rhs[k][j]);
            }

            ret[i][j] = elem;
        }

    }

    return ret;
}

template<typename T>
std::pair<matrix_t<T>, matrix_t<T>> solve_linear_system(const matrix_t<T>& left, const matrix_t<T>& right) {
    matrix_t<T> tmp(left);
    tmp.insert_col(left.get_cols_number(), right);

    gauss_straight(tmp);
    gauss_reverse(tmp);

    /* check for joint system */
    for(std::size_t i = 0, maxi = tmp.get_rows_number(); i < maxi; ++i) {
        bool flag = true; /* if all elements equal to 0 - true */
        for(std::size_t j = 0, maxj = tmp.get_cols_number() - 1; j < maxj; ++j) {
            if(!equal(tmp[i][j], 0.0)) {
                flag = false;
                break;
            }
        }

        if((flag && !equal(tmp[i][tmp.get_cols_number() - 1], 0.0))) {
            return std::make_pair(matrix_t<T>(0u, 0u), matrix_t<T>(0u, 0u));
        }
    }

    std::size_t current_m = 0;
    std::size_t current_n = 0;

    std::size_t max_m = tmp.get_rows_number();
    std::size_t max_n = tmp.get_cols_number();

    std::size_t rank = 0u;
    while((current_m < max_m) && (current_n < max_n - 1)) {
        if(matrix::equal(tmp[current_m][current_n], 0.0)) {
            ++current_n;
            continue;
        }

        /* normalize */
        T div = tmp[current_m][current_n];
        for(std::size_t i = current_n; i < max_n; ++i) {
            if(equal(tmp[current_m][i], 0.0)) {
                continue;
            }

            tmp[current_m][i] /= div;
        }

        tmp.swap_cols(rank, current_n);
        ++rank;
        ++current_n; ++current_m;
    }

    matrix_t<T> partial_solution(tmp.get_cols_number() - 1, 1);
    for(std::size_t i = 0, maxi = partial_solution.get_rows_number(); i < maxi; ++i) {
        partial_solution[i][0] = tmp[i][tmp.get_cols_number() - 1]; 
    }

    if((tmp.get_cols_number() == 1) || (tmp.get_cols_number() == (rank + 1))) {
        return {partial_solution, matrix_t<T>()};
    }

    matrix_t<T> fundamental_matrix(tmp.get_cols_number() - 1, tmp.get_cols_number() - rank - 1);
    for(std::size_t i = 0, maxi = rank; i < maxi; ++i) {
        for(std::size_t j = 0, maxj = fundamental_matrix.get_cols_number(); j < maxj; ++j) {
            fundamental_matrix[i][j] = -tmp[i][rank + j];
        }
    }

    current_m = rank;
    current_n = 0;

    max_m = fundamental_matrix.get_rows_number();
    max_n = fundamental_matrix.get_cols_number();
    while((current_m < max_m) && (current_n < max_n)) {
        fundamental_matrix[current_m][current_n] = 1.0;
        ++current_m; ++current_n;
    }

    return {partial_solution, fundamental_matrix};
}

template<typename T>
void gauss_straight(matrix_t<T>& m) {
    std::size_t current_m = 0;
    std::size_t current_n = 0;

    std::size_t max_m = m.get_rows_number();
    std::size_t max_n = m.get_cols_number();

    while((current_m < max_m) && (current_n < max_n)) {
        std::size_t max_elem = m.max_abs_col_elem(current_n, current_m, max_m);
        if(equal(m[max_elem][current_n], 0.0)) {
            ++current_n;
            continue;
        }

        if(max_elem != current_m) {
            m.swap_rows(max_elem, current_m);
        }

        for(std::size_t i = current_m + 1; i < max_m; ++i) {
            long double f = m[i][current_n] / m[current_m][current_n];
            m[i][current_n] = 0.0;
            for (std::size_t j = current_n + 1; j < max_n; ++j) {
                m[i][j] = m[i][j] - f * m[current_m][j];

                /* for accuracy of calculations */
                if(equal(m[i][j], 0.0)) {
                    m[i][j] = 0.0;
                }
            }
        }

        ++current_m;
        ++current_n;
    }
}

template<typename T>
void gauss_reverse(matrix_t<T>& m) {
    std::size_t current_m = 0;
    std::size_t current_n = 0;

    std::size_t max_m = m.get_rows_number();
    std::size_t max_n = m.get_cols_number();

    while((current_m < max_m) && (current_n < max_n)) {
        if(equal(m[current_m][current_m], 0.0)) {
            ++current_n;
            continue;
        }

        for(std::size_t i = 0; i < current_m; ++i) {
            if(equal(m[i][current_n], 0.0)) {
                continue;
            }

            long double f = m[i][current_n] / m[current_m][current_n];
            for (std::size_t j = current_n + 1; j < max_n; ++j) {
                m[i][j] = m[i][j] - f * m[current_m][j];

                /* for accuracy of calculations */
                if(equal(m[i][j], 0.0)) {
                    m[i][j] = 0.0;
                }
            }

            m[i][current_n] = 0.0;
        }

        ++current_m;
        ++current_n;
    }
}

} /* namespace matrix */
