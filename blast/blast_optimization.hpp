#pragma once
#include "blast_math.hpp"
#include "collisions/world.h"
#include "blast_manipulator.hpp"
#include "blast_trajectory.hpp"

namespace blast {

struct Optimisation {
    Manipulator* manip   = nullptr;
    Matrix*      task    = nullptr;
    Bspline*     bspline = nullptr;
    objlist*     world   = nullptr;
    int          n_collision_constraints = 5;
    int          n_collision_skip = 2;
};

//--------- OBJECTIVES AND CONSTRAINTS ----------------------------------------------------------

host_fn double obj_time(unsigned n, const double* x, double* grad, void*);
host_fn void cstr_manip(unsigned m, double *result, unsigned n, const double* x, double* grad, void* data);

//--------- INITIAL GUESS GENERATION ----------------------------------------------------------

// Fill an optimization vector with random values
// Returns the filled array
host_fn Array guess_random(Gen3_7DOF& manip, Bspline& bspline, Matrix& task);

// Fill an optimization vector by trying 'nshotgun' random vectors
//  The best optimzation vector is determined by the max (worst) value of the manip.contraints
//
//    manip     = manipulator
//    bspline   = trajectory generator
//    task      = initial and final conditions of the task
//    nshotgun  = number of tries to determine the best
//
// Return the best vector
host_fn Array guess_shot_max(Gen3_7DOF& manip, Bspline& bspline, Matrix& task, int nshotgun);

// Fill an optimization vector by trying 'nshotgun' random vectors.
//  The best optimzation vector is determined by the sum of the manip.contraints (only the violated, i.e., positive, constraints are considered in the sum).
//
//    manip     = manipulator
//    bspline   = trajectory generator
//    task      = initial and final conditions of the task
//    nshotgun  = number of tries to determine the best
//
// Return the best vector
host_fn Array guess_shot_mean(Gen3_7DOF& manip, Bspline& bspline, Matrix& task, int nshotgun);

//--------- OBJECTIVES AND CONSTRAINTS ----------------------------------------------------------

// Simple objective function of total trajectory time.
//
//  n      = number of elements in optimization vector
//  x      = optimization vector
//  grad   = if not NULL, gradient is inserted here
//
// Returns the objective value
host_fn double obj_time(unsigned n, const double* x, double* grad, void*) {
    using namespace blast;
    double result = x[n-1];
    if (grad) {
        for(u32 i = 0; i < n-1; i++)
            grad[i] = 0;
        grad[n-1] = 1;
    }
    return result;
}

host_fn real penalty_obj_time(Array x, Optimisation optim) {
    const auto points = optim.bspline->points;

    Array cstr(optim.manip->ncon(points));
    auto f = obj_time(x.size, x.data, nullptr, nullptr);
    cstr_manip(optim.manip->ncon(points), cstr.data, x.size, x.data, nullptr, &optim);

    for (int i = 0; i < (int)cstr.size; i++)
        // f += cstr[i] > 0 ? (cstr[i]* cstr[i]) : 0; // todo: explore alternatives
        // f += cstr[i] > 0 ? (abs(cstr[i])) : 0; // todo: explore alternatives
        // f += cstr[i] > 0 ? ((i+1) * cstr[i] * cstr[i]) : 0; // todo: explore alternatives
        // logarithmic penalty of abs(cstr[i])
        // exponential
        f += cstr[i] > 0 ? ((i+1) * cstr[i]) : 0; // todo: explore alternatives

    return f;
}

host_fn void internal_cstr_manip_single(unsigned m, double* result, unsigned n, const double* x, blast::Optimisation* opt) {
    Array xv;
    xv.alias(x, n);
    opt->bspline->compute_trajectory(xv, *opt->task);
    opt->manip->constraints(opt->bspline->traj, result);
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
host_fn void cstr_manip(unsigned m, double *result, unsigned xlen, const double* x, double* grad, void* f_data) {
    Optimisation* opt = (Optimisation*)f_data;

    internal_cstr_manip_single(m, result, xlen, x, opt);

    if (grad) {
        const real eps = 1e-5;
        Array x_plus(xlen);
        Array r_plus(m);
        // todo: parallel?
        for (u32 j = 0; j < xlen; j++) {
            memcpy(x_plus.data, x, xlen * sizeof(real));
            x_plus[j] += eps;
            internal_cstr_manip_single(m, r_plus.data, xlen, x_plus.data, opt);
            for (u32 i = 0; i < m; i++)
                grad[i*xlen + j] = (r_plus[i]-result[i])/eps;
        }
    }
}

host_fn void cstr_world_gen3(unsigned m, double *result, unsigned xlen, const double* x, double* grad, void* f_data) {
    Optimisation* opt = (Optimisation*)f_data;
    const int points = opt->bspline->points;
    const auto ncon = opt->manip->ncon(points);

    internal_cstr_manip_single(m, result, xlen, x, opt);

    Gen3_7DOF* manip = (Gen3_7DOF*) opt->manip;
    Matrix caps_matrix = manip->robot_capsules(opt->bspline->traj.pos, opt->n_collision_skip);
    const int caps_size = caps_matrix.cols;

    capslist capsules;
    capsules.caps.resize(caps_size * 3); // 3 capsules for each point along the trajectory
    real radius = 0.055; // Hard coded radius of all robot capsules
    for (int i = 0; i < caps_size; i++) {
        auto caps_tmp = caps_matrix.col(i);

        Vec3 p_j2 = {caps_tmp[0], caps_tmp[1], caps_tmp[2]};
        Vec3 p_j4 = {caps_tmp[3], caps_tmp[4], caps_tmp[5]};
        Vec3 p_j6 = {caps_tmp[6], caps_tmp[7], caps_tmp[8]};
        Vec3 p_ee = {caps_tmp[9], caps_tmp[10], caps_tmp[11]};

        capsules.caps[i*3].p1 = p_j2;
        capsules.caps[i*3].p2 = p_j4;
        capsules.caps[i*3].r = radius;
        capsules.caps[i*3 + 1].p1 = p_j4;
        capsules.caps[i*3 + 1].p2 = p_j6;
        capsules.caps[i*3 + 1].r = radius;
        capsules.caps[i*3 + 2].p1 = p_j6;
        capsules.caps[i*3 + 2].p2 = p_ee;
        capsules.caps[i*3 + 2].r = radius;
    }
    double* r = &result[ncon];
    std::vector<real> collisions = test_collision(&capsules, opt->world, opt->n_collision_constraints);
    for (int i = 0; i < opt->n_collision_constraints; i ++) {
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
            internal_cstr_manip_single(m, r_plus.data, xlen, x_plus.data, opt);
            caps_matrix = manip->robot_capsules(opt->bspline->traj.pos, opt->n_collision_skip);
            for (int i = 0; i < caps_size; i++) {
                auto caps_tmp = caps_matrix.col(i);

                Vec3 p_j2 = {caps_tmp[0], caps_tmp[1], caps_tmp[2]};
                Vec3 p_j4 = {caps_tmp[3], caps_tmp[4], caps_tmp[5]};
                Vec3 p_j6 = {caps_tmp[6], caps_tmp[7], caps_tmp[8]};
                Vec3 p_ee = {caps_tmp[9], caps_tmp[10], caps_tmp[11]};

                capsules.caps[i*3].p1 = p_j2;
                capsules.caps[i*3].p2 = p_j4;
                capsules.caps[i*3].r = radius;
                capsules.caps[i*3 + 1].p1 = p_j4;
                capsules.caps[i*3 + 1].p2 = p_j6;
                capsules.caps[i*3 + 1].r = radius;
                capsules.caps[i*3 + 2].p1 = p_j6;
                capsules.caps[i*3 + 2].p2 = p_ee;
                capsules.caps[i*3 + 2].r = radius;
            }

            collisions = test_collision(&capsules, opt->world, opt->n_collision_constraints);
            for (int i = 0; i < opt->n_collision_constraints; i ++) {
                r_plus[ncon + i] = -collisions[i];
            }
            for (u32 i = 0; i < m; i++)
                grad[i*xlen + j] = (r_plus[i]-result[i])/eps;
        }
    }
}

host_fn void cstr_world_link6(unsigned m, double *result, unsigned xlen, const double* x, double* grad, void* f_data) {
    Optimisation* opt = (Optimisation*)f_data;
    const int points = opt->bspline->points;
    internal_cstr_manip_single(m, result, xlen, x, opt);

    capslist capsules;
    capsules.caps.resize(points * 3); // 3 capsules for each point along the trajectory
    for (int i = 0; i < points; i++) {
        // todo: compute the link capsule positions
        // capsules.caps[i*3].p1 = ;
        // capsules.caps[i*3].p2 = ;
        // capsules.caps[i*3].r = ;
        // capsules.caps[i*3 + 1].p1 = ;
        // capsules.caps[i*3 + 1].p2 = ;
        // capsules.caps[i*3 + 1].r = ;
        // capsules.caps[i*3 + 2].p1 = ;
        // capsules.caps[i*3 + 2].p2 = ;
        // capsules.caps[i*3 + 2].r = ;
    }

    double* r = &result[opt->manip->ncon(points)];
    std::vector<real> collisions = test_collision(&capsules, opt->world, opt->n_collision_constraints);
    for (int i = 0; i < opt->n_collision_constraints; i ++) {
        *r = collisions[i];
        r++;
    }
}

//--------- INITIAL GUESS GENERATION ----------------------------------------------------------
host_fn Array guess_random(Gen3_7DOF& manip, Bspline& bspline, Matrix& task) {
    Array x(bspline.xlen(task));
    fill_random(x, 1);
    x.back() = abs(x.back()) * 5 + 0.1;
    return x;
}

host_fn Array guess_shot_max(Gen3_7DOF& manip, Bspline& bspline, Matrix& task, int nshotgun) {
    Array best_x(bspline.xlen(task));
    real best_val = INF_REAL;
    for (int i = 0; i < nshotgun; i++) {
        auto x = guess_random(manip, bspline, task);
        bspline.compute_trajectory(x, task);
        auto c = manip.constraints(bspline.traj);
        auto r = array_max(c);
        if (r < best_val) {
            best_x = x;
            best_val = r;
        }
    }
    return best_x;
}

host_fn Array guess_shot_mean(Gen3_7DOF& manip, Bspline& bspline, Matrix& task, int nshotgun) {
    Array best_x(bspline.xlen(task));
    real best_val = INF_REAL;
    for (int i = 0; i < nshotgun; i++) {
        auto x = guess_random(manip, bspline, task);
        bspline.compute_trajectory(x, task);
        auto c = manip.constraints(bspline.traj);
        real r = 0;
        for (u32 i = 0; i < c.size; i++)
            r += c[i] > 0 ? c[i] : 0;
        Assert( ! isnan(r));
        if (r < best_val) {
            best_x = x;
            best_val = r;
        }
    }
    return best_x;
}

} // namespace blast

#ifdef BLAST_ENABLE_TESTS
#include "utilities/utilities.hpp"
#ifdef __NVCC__
#endif // nvcc
#endif // tests
