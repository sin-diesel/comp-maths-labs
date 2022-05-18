#pragma once

namespace matrix {

template<typename T>
class matrix_buff_t {

protected:

    /* we need only size constructor, because for exeptions safe programming elements are constructed from matrix class */
    matrix_buff_t(std::size_t rows, std::size_t cols);
    matrix_buff_t(const matrix_buff_t& rhs) = delete;
    matrix_buff_t& operator=(const matrix_buff_t&) = delete; 
    ~matrix_buff_t();

    void swap_buffers(matrix_buff_t& rhs) noexcept;
    void resize_buffer(std::size_t rows, std::size_t cols);
    void construct_at(std::size_t row, std::size_t col, const T& elem);

    /* this method need only for operator[] overload, bad idia, but i don't know any other solution */
    T* get_ptr(std::size_t row, std::size_t col) const {return array_ + col + matrix_size_.cols_ * row;}

    T& at(std::size_t row, std::size_t col) {return array_[col + matrix_size_.cols_ * row];}
    const T& at(std::size_t row, std::size_t col) const {return array_[col + matrix_size_.cols_ * row];}

    std::size_t get_rows_number() const {return matrix_size_.rows_;}
    std::size_t get_cols_number() const {return matrix_size_.cols_;}
    std::size_t get_elements_number() const {return size_;}

private:

    void copy_construct(T* p, const T& val) { new (p) T(val); }
    void destroy(T* p) { p->~T(); }

    template <typename FwdIter>
    void destroy(FwdIter first, FwdIter last);

    T* array_;
    std::size_t size_; /* number of constructed elements */
    std::size_t capacity_;
    struct matrix_size_t {
        std::size_t rows_ = 0;
        std::size_t cols_ = 0;
    } matrix_size_;

}; /* matrix_buff_t */

/* ------------------------------------------------------------------
                        IMPLEMENTATION
---------------------------------------------------------------------*/

template<typename T>
matrix_buff_t<T>::matrix_buff_t(std::size_t rows, std::size_t cols) :
        array_((0 == (cols * rows)) ? nullptr : static_cast<T*>(::operator new((rows * cols) * sizeof(T)))),
        size_{0u},
        capacity_{rows * cols},
        matrix_size_{rows, cols} {}

template<typename T>
matrix_buff_t<T>::~matrix_buff_t() {
    if(array_) {
        destroy(array_, array_ + size_);
        ::operator delete(array_);
    }    
}

template<typename T>
void matrix_buff_t<T>::construct_at(std::size_t row, std::size_t col, const T& elem) {
    copy_construct(array_ + col + matrix_size_.cols_ * row, elem);
    ++size_;
}

template<typename T>
void matrix_buff_t<T>::swap_buffers(matrix_buff_t& rhs) noexcept {
    std::swap(array_, rhs.array_);
    std::swap(size_, rhs.size_);
    std::swap(capacity_, rhs.capacity_);
    std::swap(matrix_size_, rhs.matrix_size_);
}

template<typename T>
void matrix_buff_t<T>::resize_buffer(std::size_t rows, std::size_t cols) {
    matrix_buff_t<T> tmp(rows, cols);
    swap_buffers(tmp);
}

template<typename T>
template <typename FwdIter>
void matrix_buff_t<T>::destroy(FwdIter first, FwdIter last) {  
    while (first != last) { 
        destroy(&*first++);
    }
}

} /* namespace matrix */