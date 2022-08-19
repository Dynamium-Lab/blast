#pragma once

#include "cuda/blast_cuda.cuh"

namespace blast {


struct cuVec3 {
    double x, y, z;
};
struct cuVec3Alias {
    double* data;
};
__device__ cuVec3 operator+(cuVec3, cuVec3);
__device__ cuVec3 operator-(cuVec3, cuVec3);
__device__ cuVec3 operator*(cuVec3, double);
__device__ cuVec3 operator*(double, cuVec3);
__device__ cuVec3 cross(cuVec3 v1, cuVec3 v2);
__device__ double dot(cuVec3 a, cuVec3 b);
__device__ cuVec3 operator+(cuVec3, cuVec3);
__device__ cuVec3 operator-(cuVec3, cuVec3);
__device__ cuVec3 operator*(cuVec3, double);
__device__ cuVec3 operator*(double, cuVec3);


struct cuMat3 {
    double data[9];
    __device__ double& operator[](unsigned i);
    __device__ double& operator()(unsigned i, unsigned j);
    __device__ void zero();
};
struct cuMat3Alias {
    double* data;
    __device__ double& operator[](unsigned i);
    __device__ double& operator()(unsigned i, unsigned j);
    __device__ void zero();
};
__device__ cuMat3 transpose(cuMat3&);
__device__ cuMat3 transpose(cuMat3Alias&);
__device__ cuVec3 operator*(cuMat3&, cuVec3&);
__device__ cuVec3 operator*(cuMat3&, cuVec3Alias&);


struct cuArray {

};

struct cuMatrix {
    double* data;
};



__device__ __host__
inline double dot(double* v1, double* v2, unsigned n) {
    double r = 0;
    for (int i = 0; i < n; i++)
        r += v1[i]*v2[i];
    return r;
}




//--- cuVec3 ---
inline __device__ cuVec3 cross(cuVec3 v1, cuVec3 v2) {
    cuVec3 r;
    r.x = v1.y*v2.z - v1.z*v2.y;
    r.y = v1.z*v2.x - v1.x*v2.z;
    r.z = v1.x*v2.y - v1.y*v2.x;
    return r;
}

inline __device__ double dot(cuVec3 a, cuVec3 b) {
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

inline __device__ cuVec3 operator+(cuVec3 a, cuVec3 b) {
    return {
        a.x + b.x,
        a.y + b.y,
        a.z + b.z
    };
}

inline __device__ cuVec3 operator-(cuVec3 a, cuVec3 b) {
    return {
        a.x - b.x,
        a.y - b.y,
        a.z - b.z
    };
}

inline __device__ cuVec3 operator*(cuVec3 a, double b) {
    return {
        a.x*b,
        a.y*b,
        a.z*b
    };
}

inline __device__ cuVec3 operator*(double a, cuVec3 b) {
    return {
        b.x*a,
        b.y*a,
        b.z*a
    };
}



//--- cuMat3 ---
inline __device__ cuMat3 transpose(cuMat3& a) {
    cuMat3 r;
    r.data[0] = a.data[0];
    r.data[1] = a.data[3];
    r.data[2] = a.data[6];
    r.data[3] = a.data[1];
    r.data[4] = a.data[4];
    r.data[5] = a.data[7];
    r.data[6] = a.data[2];
    r.data[7] = a.data[5];
    r.data[8] = a.data[8];
    return r;
}

inline __device__ cuMat3 transpose(cuMat3Alias& a) {
    cuMat3 r;
    r.data[0] = a.data[0];
    r.data[1] = a.data[3];
    r.data[2] = a.data[6];
    r.data[3] = a.data[1];
    r.data[4] = a.data[4];
    r.data[5] = a.data[7];
    r.data[6] = a.data[2];
    r.data[7] = a.data[5];
    r.data[8] = a.data[8];
    return r;
}

inline __device__ void cuMat3::zero() {
    for (int i = 0; i < 9; i++)
        data[i] = 0.0f;
}

inline __device__ double& cuMat3::operator[](unsigned i) {
    return data[i];
}

inline __device__ double& cuMat3::operator()(unsigned i, unsigned j) {
    return data[3*j + i];
}

inline __device__ cuVec3 operator*(cuMat3& m, cuVec3& v) {
    cuVec3 r;
    r.x = m[0]*v.x + m[3]*v.y + m[6]*v.z;
    r.y = m[1]*v.x + m[4]*v.y + m[7]*v.z;
    r.z = m[2]*v.x + m[5]*v.y + m[8]*v.z;
    return r;
}

} // namespace blast
