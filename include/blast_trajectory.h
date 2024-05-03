#pragma once

#include "blast.h"

namespace blast {

struct Trajectory;
struct Bspline;

struct Trajectory {
    Matrix pos; // njoints x npoints
    Matrix vel; // njoints x npoints
    Matrix acc; // njoints x npoints
    Array t;    // npoints

    Trajectory(u32 npoints, u32 njoints) :
        pos(njoints, npoints),
        vel(njoints, npoints),
        acc(njoints, npoints),
        t(npoints) {}

    Trajectory() = default;
};
Trajectory compute_5order_trajectory(real T, Matrix &task);

struct Bspline {
    Trajectory traj;
    Matrix control; // nctrl x njoints
    Matrix basis_p; // nctrl x npoints
    Matrix basis_v; // nctrl x npoints
    Matrix basis_a; // nctrl x npoints
    u32 joints;
    u32 points;
    u32 nctrl;
    u32 p;

    Bspline(u32 ncontrol, u32 npoints, u32 p, u32 njoints);
    Bspline() = delete;

    // Compute a trajectory from the given optimization vector
    //  - note: fastest when 'ncontrol' is a multiple of 4 (SIMD)
    void compute_trajectory(const Array &x, const Matrix &task);
    u32 xlen(const Matrix &task);

    void compute_basis();
    void compute_basis_open();
    void compute_control(const Array &x, const Matrix &task);
    void compute_control(const Array &x, const Matrix &task, real *dst);
};


}