#pragma once

#include "blast.h"
#include "blast_error.h"
#include <functional>
#include <cstring>

namespace blast {

template <typename T_manip>
struct Optimization;

template <typename T_manip>
struct Optimization {
    T_manip*    manip;
    Matrix*     task    = nullptr;
    Bspline*    bspline = nullptr;
    World*      world   = nullptr;
    u32         n_collision_constraints = 5;
    u32         n_collision_skip = 2;
};

double obj_time(unsigned n, const double* x, double* grad, void*);

template <typename T_manip>
real penalty_obj_time(Array x, Optimization<T_manip> optim) {
    const auto points = optim.bspline->points;

    Array cstr(optim.manip->ncon(points));
    auto f = obj_time(x.size, x.data, nullptr, nullptr);
    cstr_manip(optim.manip->ncon(points), cstr.data, x.size, x.data, nullptr, &optim); // todo: Fix for Optimization<T_manip> ??

    for (int i = 0; i < (int)cstr.size; i++)
        f += cstr[i] > 0 ? ((i+1) * cstr[i]) : 0; // todo: explore alternatives

    return f;
}

template <typename T_manip>
void internal_cstr_manip_single(double* result, unsigned n, const double* x, Optimization<T_manip>* opt) {
    Array xv;
    xv.alias(x, n);
    opt->bspline->compute_trajectory(xv, *opt->task);
    opt->manip->internal_constraints(opt->bspline->traj, result);
}

// Simple manipulator defined constraints function.
//  The given manipulator's constraints function is called and no additionnal constraints are added.
//
//  m      = number of constraints
//  result = constraints are inserted here
//  n      = number of elements in optimization vector
//  x      = optimization vector
//  grad   = if not NULL, gradient is inserted here
//  data   = is cast to Optimization struct
//
// Returns nothing
template <typename T_manip>
void cstr_manip(unsigned m, double *result, unsigned xlen, const double* x, double* grad, void* f_data) {
    Optimization<T_manip>* opt = (Optimization<T_manip>*)f_data;

    internal_cstr_manip_single(result, xlen, x, opt);

    if (grad) {
        const real eps = 1e-5;
        Array x_plus(xlen);
        Array r_plus(m);
        // todo: parallel?
        for (u32 j = 0; j < xlen; j++) {
            memcpy(x_plus.data, x, xlen * sizeof(real));
            x_plus[j] += eps;
            internal_cstr_manip_single(r_plus.data, xlen, x_plus.data, opt);
            for (u32 i = 0; i < m; i++)
                grad[i*xlen + j] = (r_plus[i]-result[i])/eps;
        }
    }
}

// Constraints function that applies internal constraints from T_manip and collision constraints with the world
template <typename T_manip>
void cstr_world(unsigned m, double *result, unsigned xlen, const double* x, double* grad, void* f_data) {
    Optimization<T_manip>* opt = (Optimization<T_manip>*)f_data;
    auto manip = opt->manip;
    const int points = opt->bspline->points;
    const auto ncon = manip->ncon(points);

    internal_cstr_manip_single(result, xlen, x, opt);

    std::vector<Capsule> capsules = manip->robot_capsules(opt->bspline->traj.pos, opt->n_collision_skip);
    double* r = &result[ncon];
    auto collisions = test_collision(capsules, opt->world, opt->n_collision_constraints);
    for (u32 i = 0; i < opt->n_collision_constraints; i ++) {
        *r = -collisions[i];
        r++;
    }

    if (grad) {
        const real eps = 1e-5;
        Array x_plus(xlen);
        Array r_plus(m);
        for (u32 j = 0; j < xlen; j++) {
            memcpy(x_plus.data, x, xlen * sizeof(real));
            x_plus[j] += eps;
            internal_cstr_manip_single(r_plus.data, xlen, x_plus.data, opt);
            capsules = manip->robot_capsules(opt->bspline->traj.pos, opt->n_collision_skip);

            collisions = test_collision(capsules, opt->world, opt->n_collision_constraints);
            for (u32 i = 0; i < opt->n_collision_constraints; i ++) {
                r_plus[ncon + i] = -collisions[i];
            }
            for (u32 i = 0; i < m; i++)
                grad[i*xlen + j] = (r_plus[i]-result[i])/eps;
        }
    }
}

//--- Initial guess generators -------------------------------------------------------------------------

// Random initial guess
Array guess_random(Bspline& bspline, Matrix& task);

// Best of 'nshotgun' random guesses. Best means the lowest max constraints violation.
template <typename T_manip>
Array guess_shot_max(Optimization<T_manip>& opt, int nshotgun) {
    auto manip = opt.manip;
    auto bspline = opt.bspline;
    auto task = opt.task;
    auto world = opt.world;

    Array best_x(bspline->xlen(task));
    real best_val = INF_REAL;
    for (int i = 0; i < nshotgun; i++) {
        auto x = guess_random(*bspline, *task);
        bspline->compute_trajectory(x, *task);
        Array c1(manip->ncon(bspline->traj.t.size));
        manip->internal_constraints(bspline->traj, c1.data);

        auto capsules = manip->robot_capsules(bspline->traj.pos, opt.n_collision_skip);

        Array c2(opt.n_collision_constraints);
        std::vector<real> collisions = test_collision(capsules, world, opt.n_collision_constraints);
        for (int i = 0; i < opt.n_collision_constraints; i ++) {
            c2[i] = -collisions[i];
        }
        auto r = std::max(max(c1), max(c2));
        if (r < best_val) {
            best_x = x;
            best_val = r;
        }
    }
    return best_x;
}

// Best of 'nshotgun' random guesses. Best means the lowest average constraints violation.
template <typename T_manip>
Array guess_shot_mean(Optimization<T_manip>& opt, int nshotgun) {
    auto manip = opt.manip;
    auto bspline = opt.bspline;
    auto task = opt.task;
    auto world = opt.world;

    Array best_x(bspline->xlen(*task));
    real best_val = INF_REAL;
    for (int idx_nshot = 0; idx_nshot < nshotgun; idx_nshot++) {
        auto x = guess_random(*bspline, *task);
        bspline->compute_trajectory(x, *task);
        Array c1(manip->ncon(bspline->traj.t.size)); // todo: double check it is the right size
        manip->internal_constraints(bspline->traj, c1.data);

        std::vector<Capsule> capsules = manip->robot_capsules(bspline->traj.pos, opt.n_collision_skip);

        Array c2(opt.n_collision_constraints);
        auto collisions = test_collision(capsules, world, opt.n_collision_constraints);
        for (u32 i = 0; i < opt.n_collision_constraints; i ++) {
            c2[i] = -collisions[i];
        }
        real r = 0;
        for (u32 i = 0; i < c1.size; i++)
            r += std::max({c1[i], 0.0});
        for (u32 i = 0; i < c2.size; i++)
            r += std::max({c2[i], 0.0});
        r *= x.data[x.size - 1]; // todo: Add the right format to penalize for longer trajectories (avoid local minimums)
        Assert( ! isnan(r));
        if (r < best_val) {
            best_x = x;
            best_val = r;
        }
    }
    return best_x;
}

}

