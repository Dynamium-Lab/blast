#pragma once
#include <blast>

namespace blast {




inline blast_fn Matrix::Matrix(u32 r, u32 c) {
    size = r * c;
    rows = r;
    cols = c;
    if (size) {
        data = (real*)Malloc(ALIGN, size * sizeof(real));
        memset(data, 0, size * sizeof(real));
    }
}

inline blast_fn Matrix::Matrix(const Matrix& m) : size(m.size), cols(m.cols), rows(m.rows), is_alias(false) {
    if (size) {
        data = (real*)Malloc(ALIGN, size * sizeof(real));
        memcpy(data, m.data, size*sizeof(real));
    }
}

inline blast_fn Matrix::Matrix(Matrix&& m) noexcept : data(m.data), size(m.size), cols(m.cols), rows(m.rows), is_alias(m.is_alias) {
    m.data = nullptr;
    m.size = 0;
    m.rows = 0;
    m.cols = 0;
}

inline blast_fn Matrix::Matrix(real* d, u32 r, u32 c) {
    data = d;
    rows = r;
    cols = c;
    size = r*c;
    is_alias = true;
}

inline blast_fn Matrix::Matrix(const Array& v) : size(v.size), cols(1), rows(v.size), is_alias(false) {
    if (size) {
        data = (real*)Malloc(ALIGN, size * sizeof(real));
        memcpy(data, v.data, size*sizeof(real));
    }
}

inline blast_fn Matrix::~Matrix() {
    if (!is_alias && data)
        Free(data);
}

inline blast_fn void Matrix::resize(u32 r, u32 c) {
    Assert(!is_alias);

    if (data == nullptr) {
        data = (real*)Malloc(ALIGN, r*c * sizeof(real));
        memset(data, 0, r*c * sizeof(real));
    }
    else if (size < r*c) {
        real* tmp = (real*)Malloc(ALIGN, r*c * sizeof(real));
        memset(tmp, 0, r*c * sizeof(real));
        memcpy(tmp, data, size*sizeof(real));
        Free(data);
        data = tmp;
    }
    size = r*c;
    rows = r;
    cols = c;
    Assert(data);
}

inline blast_fn Matrix& Matrix::operator=(const Matrix& m) {
    if (this != &m) {
        if (m.size == 0) {
            size = 0;
            rows = 0;
            cols = 0;
            if (data)
                Free(data);
            data = nullptr;
            return *this;
        }
        else if (data == nullptr || is_alias) {
            data = (real*)Malloc(ALIGN, m.size * sizeof(real));
        }
        else {
            Free(data);
            data = (real*)Malloc(ALIGN, m.size * sizeof(real));
        }
        Assert(data);
        size = m.size;
        cols = m.cols;
        rows = m.rows;
        is_alias = false;
        memcpy(data, m.data, size*sizeof(real));
    }
    return *this;
}

inline blast_fn Matrix& Matrix::operator=(Matrix&& m) noexcept {
    if (this != &m) {
        if (data && !is_alias)
            Free(data);
        data = m.data;
        size = m.size;
        rows = m.rows;
        cols = m.cols;
        is_alias = m.is_alias;
        m.is_alias = true; // m is about to be destroid but we don't want it to free the memory
    }
    return *this;
}

inline blast_fn Matrix Matrix::operator-() {
    Matrix result(rows, cols);
    for (u32 i = 0; i < size; i++)
        result.data[i] = -data[i];
    return result;
}

inline blast_fn bool Matrix::operator==(const Matrix& m) const {
    Assert(cols == m.cols && rows == m.rows);
    return is_close(*this, m);
}

inline blast_fn bool Matrix::operator!=(const Matrix& m) const {
    Assert(cols == m.cols && rows == m.rows);
    auto result = is_close(*this, m);
    return !result;
}

inline blast_fn Matrix& Matrix::alias(Array& a) {
    Assert(a.data);
    if (data && !is_alias)
        Free(data);
    size = a.size;
    rows = size;
    cols = 1;
    data = a.data;
    is_alias = true;
    return *this;
}

inline blast_fn Matrix& Matrix::alias(std::vector<real>& v) {
    Assert(!v.empty());
    if (data && !is_alias)
        Free(data);
    size = (u32)v.size();
    rows = size;
    cols = 1;
    data = v.data();
    is_alias = true;
    return *this;
}

inline blast_fn Matrix& Matrix::alias(real* p, u32 r, u32 c) {
    Assert(p);
    if (data && !is_alias)
        Free(data);
    size = r*c;
    rows = r;
    cols = c;
    data = p;
    is_alias = true;
    return *this;
}

inline blast_fn real& Matrix::operator()(u32 row, u32 col) {
    Assert(row < rows && col < cols);
    return data[row + rows*col];
}

inline blast_fn real Matrix::operator()(u32 row, u32 col) const {
    Assert(row < rows && col < cols);
    return data[row + rows*col];
}

inline blast_fn real* Matrix::address(u32 row, u32 col) const {
    Assert(row < rows && col < cols);
    return &data[row + rows*col];
}

inline blast_fn Array Matrix::col(u32 c) const {
    Assert(c < this->cols);
    Array result;
    result.data = data + rows*c;
    result.size = rows;
    result.is_alias = true;
    return result;
}


// Misc operators

inline blast_fn Matrix operator+(const Matrix& m1, const Matrix& m2) {
    Assert(m1.cols == m2.cols && m1.rows == m2.rows);
    Matrix result = m1;
    for (u32 i = 0; i < m1.size; i++)
        result.data[i] += m2.data[i];

    return result;
}

inline blast_fn Matrix operator-(const Matrix& m1, const Matrix& m2) {
    Assert(m1.cols == m2.cols && m1.rows == m2.rows);
    Matrix result = m1;
    for (u32 i = 0; i < m1.size; i++)
        result.data[i] -= m2.data[i];
    return result;
}

inline blast_fn Matrix operator*(const Matrix& lhs, const Matrix& rhs) {
    Assert(lhs.cols == rhs.rows);
    Matrix result(lhs.rows, rhs.cols);
    for (u32 i = 0; i < lhs.rows; i++) {
        for (u32 j = 0; j < rhs.cols; j++) {
            real sum = 0;
            for (u32 k = 0; k < lhs.cols; k++) {
                sum += lhs(i, k) * rhs(k, j);
            }
            result(i, j) = sum;
        }
    }
    return result;
}

inline blast_fn Array operator*(const Matrix& m, const Array& v) {
    Assert(m.cols == v.size);
    Array result(m.rows);
    for (u32 i = 0; i < m.rows; i++) {
        real sum = 0;
        for (u32 j = 0; j < m.cols; j++) {
            sum += m(i, j) * v[j];
        }
        result[i] = sum;
    }
    return result;
}

inline blast_fn Matrix operator*(real& r, const Matrix& m) {
    Matrix result = m;
    for (u32 i = 0; i < m.size; i++)
        result.data[i] *= r;
    return result;
}

inline blast_fn Matrix operator/(const Matrix& m, real& r) {
    Matrix result =  m;
    for (u32 i = 0; i < m.size; i++)
        result.data[i] /= r;
    return result;
}


// Fundamental functions

inline blast_fn bool is_close(const Matrix& m1, const Matrix& m2, real eps) {
    Assert(m1.cols == m2.cols && m1.rows == m2.rows);
    for (u32 i =0; i < m1.size; i++)
        if(m1.data[i] - m2.data[i] > eps || m1.data[i] - m2.data[i] < -eps)
            return false;
    return true;
}

inline blast_fn bool is_small(const Matrix& m, real eps){
    for (u32 i = 0; i < m.size; i++)
        if (std::abs(m.data[i]) > eps)
            return false;
    return true;
}

inline blast_fn Matrix& zero(Matrix& m) {
    for (u32 i = 0; i < m.size; i++)
        m.data[i] = (real)0.0;
    return m;
}

inline blast_fn Matrix& constant(Matrix& m, real val) {
    for (u32 i = 0; i < m.size; i++)
        m.data[i] = val;
    return m;
}

inline blast_fn Matrix& fill_random(Matrix& m, real A) {
#if defined(__CUDA_ARCH__)
    curandState state;
    unsigned long long seed = clock64()+blockIdx.x*blockDim.x+threadIdx.x;
    curand_init(seed, 0, 0, &state);
    for (int i = 0; i < (int)m.size; i++)
        m.data[i] = A * curand_uniform(&state);
#else
    for (int i = 0; i < (int)m.size; i++)
        m.data[i] = A * get_random();
#endif
    return m;
}

// Linear algebra functions

inline blast_fn Matrix pw_mult(const Matrix& m1, const Matrix& m2) {
    Assert(m1.cols == m2.cols && m1.rows == m2.rows);
    auto result = m1;
    for (u32 i = 0; i < m1.size; i++)
        result.data[i] *= m2.data[i];

    return result;
}

inline blast_fn Matrix eye(int s) {
    Matrix result(s, s);
    for (u32 i = 0; i < s; ++i) {
        for (u32 j = 0; j < s; ++j) {
            result(j, i) = i == j ? 1 : 0;
        }
    }
    return result;
}

inline blast_fn Matrix transpose(const Matrix& m) {
    Matrix result(m.cols, m.rows);
    for (u32 i = 0; i < m.rows; ++i) {
        for (u32 j = 0; j < m.cols; ++j) {
            result(j, i) = m(i, j);
        }
    }
    return result;
}

inline blast_fn real determinant(const Matrix & m) {
    int n = m.rows;
    Matrix LU = m;

    for (int k = 0; k < n; k++) {
        for (int i = k + 1; i < n; i++) {
            LU(i, k) /= LU(k, k);
            for (int j = k + 1; j < n; j++) {
                LU(i, j) -= LU(i, k) * LU(k, j);
            }
        }
    }

    real det = 1.0;
    for (int i = 0; i < n; i++) {
        det *= LU(i, i);
    }

    return det;
}

inline blast_fn Matrix LU_decomp(const Matrix& m) {
    int n = m.rows;
    Matrix LU = m;

    for (int k = 0; k < n; k++) {
        for (int i = k + 1; i < n; i++) {
            LU(i, k) /= LU(k, k);
            for (int j = k + 1; j < n; j++) {
                LU(i, j) -= LU(i, k) * LU(k, j);
            }
        }
    }
    return LU;
}

inline blast_fn Array solveLU(const Matrix& LU, const Array& b) {
    int n = LU.rows;
    Array x(n, 0.0);

    for (int i = 0; i < n; i++) {
        real sum = 0.0;
        for (int j = 0; j < i; j++)
            sum += LU(i, j) * x[j];

        x[i] = b[i] - sum;
    }

    // Solve Ux = y (backsubstitution)
    for (int i = n - 1; i >= 0; i--) {
        real sum = 0.0;
        for (int j = i + 1; j < n; j++)
            sum += LU(i, j) * x[j];

        x[i] = (x[i] - sum) / LU(i, i);
    }

    return x;
}

inline blast_fn Matrix inverse(const Matrix& m) {
    Assert(m.rows == m.cols); // square matrix

    int n = m.rows;
    Matrix LU = LU_decomp(m);
    Matrix inv(n, n);

    for (int i = 0; i < n; i++) {
        // unit vector in one dimension
        Array unit(n);
        unit[i] = 1.0;

        // Solve for the i-th column of the inverse
        Array column = solveLU(LU, unit);

        // insert the column into the inverse matrix
        inv.col(i) = column;
    }
    return inv;
}

inline blast_fn Matrix pinv(const Matrix& m) {
    Matrix res(m.cols, m.rows);
    if (m.cols >= m.rows) {
        res = transpose(m) * inverse(m * transpose(m));
    }
    else {
        res = inverse(transpose(m) * m) * transpose(m);
    }
    return res;
}

inline blast_fn void QR_decomp(const Matrix& A, Matrix& Q, Matrix& R) {
    Matrix u(A.rows, A.cols);
    Array e(A.rows);
    Array sum(A.cols);

    for(u32 i = 0; i < A.cols; i++) {
        auto a = A.col(0);
        u(i, 0) = a[i];
    }
    const auto norm_u = norm(u.col(0));
    const auto ucol = u.col(0);
    e = (ucol / norm_u);
    for(u32 i = 0; i < A.rows; i++) {
        Q(i, 0) = e[i];
    }

    for (int j = 1; j < (int)A.cols; j ++) {
        const auto a = A.col(j);

        zero(sum);
        for (int i = 0; i < j; i++) {

            auto proj_u_a = dot(u.col(i), a) / dot(u.col(i), u.col(i)) * u.col(i);
            sum = sum + proj_u_a;
        }
        for (u32 i = 0; i < A.rows; i++)
            u(i, j) = a[i] - sum[i];

        if (norm(u.col(j)) != 0) {
            e = u.col(j) / norm(u.col(j));
            for(u32 i = 0; i < A.rows; i++) {
                Q(i, j) = e[i];
            }
        }
    }
    real oper = -1;
    Q = oper*Q;
    auto Q_t = transpose(Q);
    R = Q_t*A;
}

inline blast_fn Array eigen_values(const Matrix& A, Matrix& Q, Matrix& R, Matrix& eig_vec) {
    Array eig_val(A.rows);
    Array sing_val(A.rows);

    real diff = INF_REAL;
    real tol = (real)1e-8;
    int k = 0;
    int maxiter = 1000;

    auto A_new = A;
    // find eigenvalues from QR decomposition
    while (diff > tol && k < maxiter) {
        auto A_old = A_new;
        QR_decomp(A_old, Q, R);

        A_new = R*Q;
        if (k == 0)
            eig_vec = Q;
        else
            eig_vec = eig_vec*Q;

        auto abs_max = -INF_REAL;
        for (u32 i = 0; i < A.size; i++) {
            auto abs_tmp = std::abs(A_new.data[i] - A_old.data[i]);
            if (abs_tmp > abs_max) {
                diff = abs_tmp;
                abs_max = abs_tmp;
            }
        }
        k++;
    }
    for (u32 i = 0; i < A.rows; i++)
        eig_val[i] = A_new(i, i);

    // todo: Reorder eigenvectors and eigenvalues

    return eig_val;
}

inline blast_fn void SVD_decomp(const Matrix& A, Matrix& U, Matrix& S, Matrix& V) {
    Matrix Q1(A.rows, A.cols);
    Matrix Q2(A.rows, A.cols);
    Matrix R1(A.rows, A.cols);
    Matrix R2(A.rows, A.cols);

    auto A_t = transpose(A);
    auto AAT = A*A_t; // for U
    auto ATA = A_t*A; // for V

    auto eig_val = eigen_values(AAT, Q1, R1, U);
    eigen_values(ATA, Q2, R2, V);

    for(u32 i = 0; i < A.rows; i++) {
        if (eig_val[i] > 0)
            S(i, i) = sqrt(eig_val[i]);
    }
}

inline blast_fn Matrix inverse_svd(const Matrix& m, Matrix& U, Matrix& S, Matrix& V) {
    SVD_decomp(m, U, S, V);
    auto V_t = transpose(V);
    auto S_inv = S;
    for(u32 i = 0; i < S.cols; i++) {
        if(S_inv(i, i) != 0)
            S_inv(i, i) = 1/S(i, i);
    }
    S_inv = transpose(S_inv);
    auto W_new = U*S*V_t;
    auto W_inv = U*S_inv*V_t;
    W_inv = transpose(W_inv);

    return W_inv;
}

inline blast_fn Matrix pinv_svd(const Matrix& m) {
    Matrix res(m.cols, m.rows);
    if (m.cols >= m.rows) {
        auto W = m * transpose(m);
        Matrix S(W.rows, W.cols);
        Matrix U(W.rows, W.cols);
        Matrix V(W.rows, W.cols);
        auto W_inv = inverse_svd(W, U, S, V);
        res = transpose(m) * inverse_svd(W, U, S, V);
    }
    else {
        auto W = transpose(m) * m;
        Matrix S(W.rows, W.cols);
        Matrix U(W.rows, W.cols);
        Matrix V(W.rows, W.cols);
        res = inverse_svd(W, U, S, V) * transpose(m);
    }

    return res;
}


}// namespace blast