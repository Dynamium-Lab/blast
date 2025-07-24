#pragma once
#include <blast>
#include <vector>

namespace blast {

/// @brief A simple matrix class that stores objects of type T. Data are stored in a column-major order.
template <typename T> struct ObjMatrix {
    std::vector<T> data;
    int cols = 0;
    int rows = 0;
    int size = 0;

    inline ObjMatrix() = default;

    inline ObjMatrix(int rows, int cols);
    inline ObjMatrix(const ObjMatrix& other);
    inline ObjMatrix(ObjMatrix&& other) noexcept;
    inline ObjMatrix& operator=(const ObjMatrix& other);
    inline ObjMatrix& operator=(ObjMatrix&& other) noexcept;

    inline T& operator()(int r, int c);
    inline const T& operator()(int r, int c) const;

    inline void resize(int rows, int cols);

    inline std::vector<T> col_copy(int c) const;

    inline ~ObjMatrix() = default;
};




/// @brief Constructs an ObjMatrix with the specified number of columns and rows.
/// @param rows The number of rows.
/// @param cols The number of columns.
template <typename T>
inline ObjMatrix<T>::ObjMatrix(int rows, int cols) : cols(cols), rows(rows), size(rows*cols) {
    Assert(r >= 0 && c >= 0);
    data.resize(cols * rows);
}

/// @brief Copy constructor. Does a DEEP copy of the data.
/// @param other The matrix to copy.
template <typename T>
inline ObjMatrix<T>::ObjMatrix(const ObjMatrix& other) : cols(other.cols), rows(other.rows), size(other.size), data(other.data) {}

/// @brief Move constructor. Does NOT copy the data.
/// @param other The matrix to move.
template <typename T>
inline ObjMatrix<T>::ObjMatrix(ObjMatrix&& other) noexcept : cols(other.cols), rows(other.rows), size(other.size), data(std::move(other.data)) {
    other.cols = 0;
    other.rows = 0;
    other.size = 0;
}

/// @brief Copy assignment operator. Does a DEEP copy of the data.
/// @param other The matrix to copy.
/// @return A reference to the copied matrix.
template <typename T>
inline ObjMatrix<T>& ObjMatrix<T>::operator=(const ObjMatrix& other) {
    if (this != &other) {
        cols = other.cols;
        rows = other.rows;
        size = other.size;

        // note (andre): am I stupid? the following line gives warning, but I don't know why!
        //       data = other.data;
        data.resize(other.size);
        std::copy(other.data.begin(), other.data.end(), data.begin());
    }
    return *this;
}

/// @brief Move assignment operator. Does NOT copy the data.
/// @param other The matrix to move.
/// @return A reference to the moved matrix.
template <typename T>
inline ObjMatrix<T>& ObjMatrix<T>::operator=(ObjMatrix&& other) noexcept {
    if (this != &other) {
        cols = other.cols;
        rows = other.rows;
        size = other.size;
        data = std::move(other.data);
        other.cols = 0;
        other.rows = 0;
        other.size = 0;
    }
    return *this;
}

/// @brief Accesses the element at the specified row and column.
/// @param r The row index.
/// @param c The column index.
/// @return A reference to the element at the specified row and column.
template <typename T>
inline T& ObjMatrix<T>::operator()(int r, int c) {
    Assert(r >= 0 && c >= 0);
    Assert(r < rows && c < cols);
    return data[c * rows + r];
}

/// @brief Accesses the element at the specified row and column (const version).
/// @param r The row index.
/// @param c The column index.
/// @return A const reference to the element at the specified row and column.
template <typename T>
inline const T& ObjMatrix<T>::operator()(int r, int c) const {
    Assert(r >= 0 && c >= 0);
    Assert(r < rows && c < cols);
    return data[c * rows + r];
}

/// @brief Resizes the matrix to the specified number of columns and rows.
/// @param rows Number of rows.
/// @param cols Number of columns.
template<typename T>
inline void ObjMatrix<T>::resize(int r, int c) {
    Assert(r >= 0 && c >= 0);
    if (size == 0) {
        this->rows = r;
        this->cols = c;
        this->size = r * c;
        data.resize(c * r);
    }
    else {
        auto old_matrix = *this;

        this->rows = r;
        this->cols = c;
        this->size = r * c;
        data.clear();
        data.resize(c * r);

        // Copy the old data into the new matrix but only up to the minimum of the old and new dimensions
        for (int i = 0; i < std::min(old_matrix.rows, r); i++) {
            for (int j = 0; j < std::min(old_matrix.cols, c); j++) {
                (*this)(j, i) = old_matrix(j, i);
            }
        }
    }
}


/// @brief Copy an entire column and return in a vector.
/// @param c Column to copy.
template<typename T>
inline blast_fn std::vector<T> ObjMatrix<T>::col_copy(int c) const {
    Assert(c < cols);
    Assert(c >= 0);
    std::vector<T> result(rows);
    int index = c * rows;
    for (int r = 0; r < rows; r++)
        result[r] = data[index++];
    return result;
}



} // namespace blast