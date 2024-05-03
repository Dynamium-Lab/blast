#pragma once
#include "blast.hpp"
#include "blast_math.hpp"
#include "collisions/world.h"
#include "manipulators/kinova_gen3.h"
#include "manipulators/kinova_link6.h"
#include "blast_manipulator.hpp"
#include "collisions/Collision_pso.hpp"
#include "collisions/Collision_gwo.hpp"

namespace blast {

/* ----------- NEW OPTIMIZATION ---------------------------------------------------------------

// ----------- blast_optimization.h -----------------------------------------------------------
template <typename T_manip>
struct Optimization {
    T_manip*    manip;
    Matrix*     task    = nullptr;
    Bspline*    bspline = nullptr;
    objlist*    world   = nullptr;
    int         n_collision_constraints = 5;
    int         n_collision_skip = 2;
};

host_fn double obj_time(unsigned n, const double* x, double* grad, void*);

template <typename T_manip>
host_fn real penalty_obj_time(Array x, Optimization<T_manip> optim);

template <typename T_manip>
host_fn void internal_cstr_manip_single(unsigned m, double* result, unsigned n, const double* x, Optimization<T_manip>* opt);

template <typename T_manip>
host_fn void cstr_manip(unsigned m, double *result, unsigned n, const double* x, double* grad, void* data);

template <typename T_manip>
host_fn void cstr_world(unsigned m, double *result, unsigned xlen, const double* x, double* grad, void* f_data);

host_fn Array guess_random(Bspline& bspline, Matrix& task);
template <typename T_manip>
host_fn Array guess_shot_max(Optimization<T_manip>& opt, int nshotgun);
template <typename T_manip>
host_fn Array guess_shot_mean_collisions_link6(Optimization<T_manip>& opt, int nshotgun);

// ----------- blast_optimization.cpp ---------------------------------------------------------
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

template <typename T_manip>
host_fn real penalty_obj_time(Array x, Optimization<T_manip> optim) {
    const auto points = optim.bspline->points;

    Array cstr(optim.manip->ncon(points));
    auto f = obj_time(x.size, x.data, nullptr, nullptr);
    cstr_manip(optim.manip->ncon(points), cstr.data, x.size, x.data, nullptr, &optim);

    for (int i = 0; i < (int)cstr.size; i++)
        f += cstr[i] > 0 ? ((i+1) * cstr[i]) : 0; // todo: explore alternatives

    return f;
}

template <typename T_manip>
host_fn void internal_cstr_manip_single(unsigned m, double* result, unsigned n, const double* x, Optimization<T_manip>* opt) {
    Array xv;
    xv.alias(x, n);
    opt->bspline->compute_trajectory(xv, *opt->task);
    opt->manip->internal_constraints(opt->bspline->traj, result);
}

template <typename T_manip>
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
    Optimization<T_manip>* opt = (Optimization<T_manip>*)f_data;

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

template <typename T_manip>
host_fn void cstr_world(unsigned m, double *result, unsigned xlen, const double* x, double* grad, void* f_data) {
    Optimization<T_manip>* opt = (Optimization<T_manip>*)f_data;
    const int points = opt->bspline->points;
    const auto ncon = opt->manip->ncon(points);

    internal_cstr_manip_single(m, result, xlen, x, opt);

    auto manip = opt->manip;
    capslist capsules = manip->robot_capsules(opt->bspline->traj.pos, opt->n_collision_skip);

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
            capsules = manip->robot_capsules(opt->bspline->traj.pos, opt->n_collision_skip);

            collisions = test_collision(&capsules, opt->world, opt->n_collision_constraints);
            for (int i = 0; i < opt->n_collision_constraints; i ++) {
                r_plus[ncon + i] = -collisions[i];
            }
            for (u32 i = 0; i < m; i++)
                grad[i*xlen + j] = (r_plus[i]-result[i])/eps;
        }
    }
}


host_fn Array guess_random(Bspline& bspline, Matrix& task) {
    Array x(bspline.xlen(task));
    fill_random(x, 1);
    x.back() = abs(x.back()) * 5 + 0.1;
    return x;
}

template <typename T_manip>
host_fn Array guess_shot_max(Optimization<T_manip>& opt, int nshotgun) {
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
        std::vector<real> collisions = test_collision(&capsules, world, opt.n_collision_constraints);
        for (int i = 0; i < opt.n_collision_constraints; i ++) {
            c2[i] = -collisions[i];
        }
        auto r = max(array_max(c1), array_max(c2));
        if (r < best_val) {
            best_x = x;
            best_val = r;
        }
    }
    return best_x;
}

template <typename T_manip>
host_fn Array guess_shot_mean_collisions_link6(Optimization<T_manip>& opt, int nshotgun) {
    auto manip = opt.manip;
    auto bspline = opt.bspline;
    auto task = opt.task;
    auto world = opt.world;

    Array best_x(bspline->xlen(*task));
    real best_val = INF_REAL;
    for (int idx_nshot = 0; idx_nshot < nshotgun; idx_nshot++) {
        auto x = guess_random(*manip, *bspline, *task);
        bspline->compute_trajectory(x, *task);
        Array c1(manip->ncon(bspline->traj.t.size)); // todo: double check it is the right size
        manip->internal_constraints(bspline->traj, c1.data);

        auto capsules = manip->robot_capsules(bspline->traj.pos, opt.n_collision_skip);

        Array c2(opt.n_collision_constraints);
        std::vector<real> collisions = test_collision(&capsules, world, opt.n_collision_constraints);
        for (int i = 0; i < opt.n_collision_constraints; i ++) {
            c2[i] = -collisions[i];
        }
        real r = 0;
        for (u32 i = 0; i < c1.size; i++)
            r += c1[i] > 0 ? c1[i] : 0;
        for (u32 i = 0; i < c2.size; i++)
            r += c2[i] > 0 ? c2[i] : 0;
        Assert( ! isnan(r));
        if (r < best_val) {
            best_x = x;
            best_val = r;
        }
    }
    return best_x;
}


// -------- END OF NEW OPTIMIZATION -------------------------------------------------------------*/



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
host_fn Array guess_random(Gen3& manip, Bspline& bspline, Matrix& task);

// Fill an optimization vector by trying 'nshotgun' random vectors
//  The best optimzation vector is determined by the max (worst) value of the manip.contraints
//
//    manip     = manipulator
//    bspline   = trajectory generator
//    task      = initial and final conditions of the task
//    nshotgun  = number of tries to determine the best
//
// Return the best vector
host_fn Array guess_shot_max(Gen3& manip, Bspline& bspline, Matrix& task, int nshotgun);

// Fill an optimization vector by trying 'nshotgun' random vectors.
//  The best optimzation vector is determined by the sum of the manip.contraints (only the violated, i.e., positive, constraints are considered in the sum).
//
//    manip     = manipulator
//    bspline   = trajectory generator
//    task      = initial and final conditions of the task
//    nshotgun  = number of tries to determine the best
//
// Return the best vector
host_fn Array guess_shot_mean(Gen3& manip, Bspline& bspline, Matrix& task, int nshotgun);

// note: CUDA stuff, only enabled if compiling for Nvidia GPUs
#ifdef __NVCC__
// compute the trajectory and the constraints (slower, but access to trajectory)
// __global__ void pva_constraints_kernel(cuBspline pva);

// compute the constraints, but don't store the trajectory (faster, but no access to trajectory)
// __global__ void constraints_no_pva_kernel(cuBspline pva);
#endif

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

host_fn void cstr_world_gen3_pso(unsigned m, double *result, unsigned xlen, const double* x, double* grad, void* f_data) {
    Optimisation* opt = (Optimisation*)f_data;
    const int points = opt->bspline->points;
    const auto ncon = opt->manip->ncon(points);

    internal_cstr_manip_single(m, result, xlen, x, opt);

    Gen3_7DOF* manip = (Gen3_7DOF*) opt->manip;
    Matrix caps_matrix = manip->robot_capsules(opt->bspline->traj.pos, opt->n_collision_skip);
    const int caps_size = caps_matrix.cols;
    Matrix test_caps (6, caps_size); 

    real radius = 0.055; // Hard coded radius of all robot capsules
    Array min_dist(opt->n_collision_constraints);
    double* r = &result[ncon];
    // const int ncol = 3*size(opt->world->OBBlist);
    const int ncol = ((caps_matrix.rows/3)-1)*size(opt->world->OBBlist);

    std::vector<real> temp(ncol); 
    int idx = 0;
    for (int i = 0; i < caps_matrix.rows/3 - 1; i++) {
        for (int k = 0; k < caps_size; k++) {
            auto caps_tmp = caps_matrix.col(k);
            Vec3 elem1 = {caps_tmp[3*i], caps_tmp[3*i+1], caps_tmp[3*i+2]};
            Vec3 elem2 = {caps_tmp[3*i + 3], caps_tmp[3*i+4], caps_tmp[3*i+5]};
            
            test_caps(0, k) = elem1.x;
            test_caps(1, k) = elem1.y;
            test_caps(2, k) = elem1.z;
            test_caps(3, k) = elem2.x;
            test_caps(4, k) = elem2.y;
            test_caps(5, k) = elem2.z;
    
        }
        for (auto OBB:opt->world->OBBlist) {
            // temp[idx++] = collision_pso(test_caps, OBB) - radius;
            // temp[idx++] = collision_gwo(test_caps, OBB) - radius;

        }
    }
    std::sort(temp.begin(), temp.end());
    for (int i = 0; i < opt->n_collision_constraints; i ++) {
        *r = -temp[i];
        r++;
    }
    // if (grad) {
    //     const real eps = 1e-5;
    //     Array x_plus(xlen);
    //     Array r_plus(m);
    //     for (u32 j = 0; j < xlen; j++) {
    //         memcpy(x_plus.data, x, xlen * sizeof(real));
    //         x_plus[j] += eps;
    //         internal_cstr_manip_single(m, r_plus.data, xlen, x_plus.data, opt);
    //         caps_matrix = manip->robot_capsules(opt->bspline->traj.pos, opt->n_collision_skip);
    //         const int ncol = 3*size(opt->world->OBBlist);
    //         std::vector<real> temp(ncol); //not working 
    //         int idx = 0;
    //         for (int i = 0; i < caps_matrix.rows - 1; i++) {
    //             for (int k = 0; k < caps_size; k++) {
    //                 auto caps_tmp = caps_matrix.col(k);
    //                 Vec3 elem1 = {caps_tmp[i], caps_tmp[i+1], caps_tmp[i+2]};
    //                 Vec3 elem2 = {caps_tmp[3*i], caps_tmp[3*i+1], caps_tmp[3*i+2]};
                    
    //                 test_caps(0, k) = elem1.x;
    //                 test_caps(1, k) = elem1.y;
    //                 test_caps(2, k) = elem1.z;
    //                 test_caps(3, k) = elem2.x;
    //                 test_caps(4, k) = elem2.y;
    //                 test_caps(5, k) = elem2.z;
            
    //             }
    //             for (auto OBB:opt->world->OBBlist) {
    //                 temp[idx++] = collision_pso(test_caps, OBB);
    //             }
    //         }
    //         std::sort(temp.begin(), temp.end());
    //         for (int i = 0; i < opt->n_collision_constraints; i ++) {
    //             r_plus[ncon + i] = -temp[i];
    //         }
    //         for (u32 i = 0; i < m; i++)
    //         grad[i*xlen + j] = (r_plus[i]-result[i])/eps;
    //     }
    // }
}
host_fn void cstr_world_gen3(unsigned m, double *result, unsigned xlen, const double* x, double* grad, void* f_data) {
    Optimisation* opt = (Optimisation*)f_data;
    const int points = opt->bspline->points;
    const auto ncon = opt->manip->ncon(points);

    internal_cstr_manip_single(m, result, xlen, x, opt);

    Gen3* manip = (Gen3*) opt->manip;
    auto capsules = manip->robot_capsules(opt->bspline->traj.pos, opt->n_collision_skip);

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
            capsules = manip->robot_capsules(opt->bspline->traj.pos, opt->n_collision_skip);


            collisions = test_collision(&capsules, opt->world, opt->n_collision_constraints);
            for (int i = 0; i < opt->n_collision_constraints; i ++) {
                r_plus[ncon + i] = -collisions[i];
            }
            for (u32 i = 0; i < m; i++)
                grad[i*xlen + j] = (r_plus[i]-result[i])/eps;
        }
    }
}

// host_fn void cstr_world_link6(unsigned m, double *result, unsigned xlen, const double* x, double* grad, void* f_data) {
//     Optimisation* opt = (Optimisation*)f_data;
//     const int points = opt->bspline->points;
//     const auto ncon = opt->manip->ncon(points);

//     internal_cstr_manip_single(m, result, xlen, x, opt);

//     Link6_dev* manip = (Link6_dev*) opt->manip;
//     const int n_collision_tests = opt->bspline->traj.pos.cols / opt->n_collision_skip;
//     Matrix tmp_matrix(7, 5);
//     Matrix caps_matrix(7, n_collision_tests);
//     for (int i = 0; i < n_collision_tests; i++) {
//         tmp_matrix = manip->robot_capsules(opt->bspline->traj.pos.col(i*opt->n_collision_skip));
//         for (int j = 0; j < tmp_matrix.cols; j++) {
//             caps_matrix.col(i*tmp_matrix.cols + j) = tmp_matrix.col(j);
//         }
//     }
//     const int caps_size = caps_matrix.cols;

//     capslist capsules;
//     capsules.caps.resize(caps_size * 5); // 5 capsules for each point along the trajectory
//     for (int i = 0; i < caps_size; i++) {
//         auto caps_tmp = caps_matrix.col(i);

//         capsules.caps[i*5].p1 = { caps_tmp[0, 0], caps_tmp[0, 1], caps_tmp[0, 2] };
//         capsules.caps[i*5].p2 = { caps_tmp[0, 3], caps_tmp[0, 4], caps_tmp[0, 5] };
//         capsules.caps[i*5].r = caps_tmp[0, 6];
//         capsules.caps[i*5 + 1].p1 = { caps_tmp[1, 0], caps_tmp[1, 1], caps_tmp[1, 2] };
//         capsules.caps[i*5 + 1].p2 = { caps_tmp[1, 3], caps_tmp[1, 4], caps_tmp[1, 5] };
//         capsules.caps[i*5 + 1].r = caps_tmp[1, 6];
//         capsules.caps[i*5 + 2].p1 = { caps_tmp[2, 0], caps_tmp[2, 1], caps_tmp[2, 2] };
//         capsules.caps[i*5 + 2].p2 = { caps_tmp[2, 3], caps_tmp[2, 4], caps_tmp[2, 5] };
//         capsules.caps[i*5 + 2].r = caps_tmp[2, 6];
//         capsules.caps[i*5 + 3].p1 = { caps_tmp[3, 0], caps_tmp[3, 1], caps_tmp[3, 2] };
//         capsules.caps[i*5 + 3].p2 = { caps_tmp[3, 3], caps_tmp[3, 4], caps_tmp[3, 5] };
//         capsules.caps[i*5 + 3].r = caps_tmp[3, 6];
//         capsules.caps[i*5 + 4].p1 = { caps_tmp[4, 0], caps_tmp[4, 1], caps_tmp[4, 2] };
//         capsules.caps[i*5 + 4].p2 = { caps_tmp[4, 3], caps_tmp[4, 4], caps_tmp[4, 5] };
//         capsules.caps[i*5 + 4].r = caps_tmp[4, 6];
//     }
//     double* r = &result[ncon];
//     std::vector<real> collisions = test_collision(&capsules, opt->world, opt->n_collision_constraints);
//     for (int i = 0; i < opt->n_collision_constraints; i ++) {
//         *r = -collisions[i];
//         r++;
//     }

//     if (grad) {
//         const real eps = 1e-5;
//         Array x_plus(xlen);
//         Array r_plus(m);
//         for (u32 j = 0; j < xlen; j++) {
//             memcpy(x_plus.data, x, xlen * sizeof(real));
//             x_plus[j] += eps;
//             internal_cstr_manip_single(m, r_plus.data, xlen, x_plus.data, opt);
//             for (int i = 0; i < n_collision_tests; i++) {
//                 tmp_matrix = manip->robot_capsules(opt->bspline->traj.pos.col(i*opt->n_collision_skip));
//                 for (int j = 0; j < tmp_matrix.cols; j++) {
//                     caps_matrix.col(i*tmp_matrix.cols + j) = tmp_matrix.col(j);
//                 }
//             }
//             const int caps_size = caps_matrix.cols;

//             capslist capsules;
//             capsules.caps.resize(caps_size * 5); // 5 capsules for each point along the trajectory
//             for (int i = 0; i < caps_size; i++) {
//                 auto caps_tmp = caps_matrix.col(i);

//                 capsules.caps[i*5].p1 = { caps_tmp[0, 0], caps_tmp[0, 1], caps_tmp[0, 2] };
//                 capsules.caps[i*5].p2 = { caps_tmp[0, 3], caps_tmp[0, 4], caps_tmp[0, 5] };
//                 capsules.caps[i*5].r = caps_tmp[0, 6];
//                 capsules.caps[i*5 + 1].p1 = { caps_tmp[1, 0], caps_tmp[1, 1], caps_tmp[1, 2] };
//                 capsules.caps[i*5 + 1].p2 = { caps_tmp[1, 3], caps_tmp[1, 4], caps_tmp[1, 5] };
//                 capsules.caps[i*5 + 1].r = caps_tmp[1, 6];
//                 capsules.caps[i*5 + 2].p1 = { caps_tmp[2, 0], caps_tmp[2, 1], caps_tmp[2, 2] };
//                 capsules.caps[i*5 + 2].p2 = { caps_tmp[2, 3], caps_tmp[2, 4], caps_tmp[2, 5] };
//                 capsules.caps[i*5 + 2].r = caps_tmp[2, 6];
//                 capsules.caps[i*5 + 3].p1 = { caps_tmp[3, 0], caps_tmp[3, 1], caps_tmp[3, 2] };
//                 capsules.caps[i*5 + 3].p2 = { caps_tmp[3, 3], caps_tmp[3, 4], caps_tmp[3, 5] };
//                 capsules.caps[i*5 + 3].r = caps_tmp[3, 6];
//                 capsules.caps[i*5 + 4].p1 = { caps_tmp[4, 0], caps_tmp[4, 1], caps_tmp[4, 2] };
//                 capsules.caps[i*5 + 4].p2 = { caps_tmp[4, 3], caps_tmp[4, 4], caps_tmp[4, 5] };
//                 capsules.caps[i*5 + 4].r = caps_tmp[4, 6];
//             }

//             collisions = test_collision(&capsules, opt->world, opt->n_collision_constraints);
//             for (int i = 0; i < opt->n_collision_constraints; i ++) {
//                 r_plus[ncon + i] = -collisions[i];
//             }
//             for (u32 i = 0; i < m; i++)
//                 grad[i*xlen + j] = (r_plus[i]-result[i])/eps;
//         }
//     }
// }

//--------- INITIAL GUESS GENERATION ----------------------------------------------------------
host_fn Array guess_random(Gen3& manip, Bspline& bspline, Matrix& task) {
    Array x(bspline.xlen(task));
    fill_random(x, 1);
    x.back() = abs(x.back()) * 5 + 0.1;
    return x;
}

host_fn Array guess_shot_max(Gen3& manip, Bspline& bspline, Matrix& task, int nshotgun) {
    Array best_x(bspline.xlen(task));
    real best_val = INF_REAL;
    for (int i = 0; i < nshotgun; i++) {
        auto x = guess_random(manip, bspline, task);
        bspline.compute_trajectory(x, task);
        Array c(manip.ncon(bspline.traj.t.size));
        manip.internal_constraints(bspline.traj, c.data);
        auto r = array_max(c);
        if (r < best_val) {
            best_x = x;
            best_val = r;
        }
    }
    return best_x;
}

host_fn Array guess_shot_mean(Gen3& manip, Bspline& bspline, Matrix& task, int nshotgun) {
    Array best_x(bspline.xlen(task));
    real best_val = INF_REAL;
    for (int i = 0; i < nshotgun; i++) {
        auto x = guess_random(manip, bspline, task);
        bspline.compute_trajectory(x, task);
        Array c(manip.ncon(bspline.traj.t.size));
        manip.internal_constraints(bspline.traj, c.data);
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

// note: CUDA stuff, only enabled if compiling for Nvidia GPUs
// #ifdef __NVCC__
// __global__ void pva_constraints_kernel(cuBspline pva) {
//     const u32 point = blockIdx.x * blockDim.x + threadIdx.x;
//     pva.compute_trajectory(point);
//     blast::cuGen3MultiTraj* manip = (blast::cuGen3MultiTraj*)blast::manip_broadcast_arena;
//     const u32 pva_offset = point * 7;
//     const u32 constraints_offset = point * 21;
//     manip->compute_constraints(
//         pva.device_pos+pva_offset,
//         pva.device_vel+pva_offset,
//         pva.device_acc+pva_offset,
//         manip->dev_constraints + constraints_offset
//     );
// }

// __global__ void constraints_no_pva_kernel(cuBspline pva) {
//     const u32 point = blockIdx.x * blockDim.x + threadIdx.x;
//     const u32 nctrl = pva.ncontrol;
//     const real one_over_T = pva.one_over_T;
//     const real one_over_T2 = pva.one_over_T2;
//     if (point == 0)
//         pva.device_t[point] = pva.dt * point;
//     __syncthreads();
//     real pos[7];
//     real vel[7];
//     real acc[7];
//     const auto bp = &pva.device_basis_p[point*nctrl];
//     const auto bv = &pva.device_basis_v[point*nctrl];
//     const auto ba = &pva.device_basis_a[point*nctrl];
//     for (int joint = 0; joint < 7; joint++) {
//         const auto c = pva.device_control + joint*nctrl;
//         real position = 0;
//         real velocity = 0;
//         real acceleration = 0;
//         for (int i = 0; i < nctrl; i++) {
//             position     += c[i] * bp[i];
//             velocity     += c[i] * bv[i];
//             acceleration += c[i] * ba[i];
//         }
//         pos[joint] = position;
//         vel[joint] = velocity     * one_over_T;
//         acc[joint] = acceleration * one_over_T2;
//     }
//     blast::cuGen3MultiTraj* manip = (blast::cuGen3MultiTraj*)blast::manip_broadcast_arena;
//     const u32 constraints_offset = point * 21;
//     manip->compute_constraints(pos, vel, acc, manip->dev_constraints + constraints_offset);
// }

// // kernel that uses shared memory to speed up constraint computation
// __global__ void constraints_shared_kernel(cuBspline pva) {
//     const u32 point = blockIdx.x * blockDim.x + threadIdx.x;
//     const u32 nctrl = pva.ncontrol;

//     // __shared__ real bp[];
//     // __shared__ real bv[];
//     // __shared__ real ba[];
//     // __shared__ real control[12*7];
//     // if (threadIdx.x == 0)
//     //     cudaMemcpyAsync(control, pva.device_control, sizeof(real) * 12 * 7, cudaMemcpyDeviceToDevice);

//     // __syncthreads();

//     // cudaDeviceSynchronize();

//     real pos[7];
//     real vel[7];
//     real acc[7];
//     const auto bp = &pva.device_basis_p[point*nctrl];
//     const auto bv = &pva.device_basis_v[point*nctrl];
//     const auto ba = &pva.device_basis_a[point*nctrl];
//     for (int joint = 0; joint < 7; joint++) {
//         const auto c = pva.device_control + joint*nctrl;
//         real position = 0;
//         real velocity = 0;
//         real acceleration = 0;
//         for (int i = 0; i < nctrl; i++) {
//             position     += c[i] * bp[i];
//             velocity     += c[i] * bv[i];
//             acceleration += c[i] * ba[i];
//         }
//         pos[joint] = position;
//         vel[joint] = velocity     * pva.one_over_T;
//         acc[joint] = acceleration * pva.one_over_T2;
//     }

//     blast::cuGen3MultiTraj* manip = (blast::cuGen3MultiTraj*)blast::manip_broadcast_arena;
//     const u32 constraints_offset = point * 21;
//     manip->compute_constraints(pos, vel, acc, manip->dev_constraints + constraints_offset);
// }
// #endif

} // namespace blast

#ifdef BLAST_ENABLE_TESTS
#include "utilities/utilities.hpp"
#ifdef __NVCC__
// TEST_CASE("GpuCpuCorrectness", "[Manipulator]") {
//     using namespace blast;

//     printf("Testing on GPU with the following properties:\n");
//     print_device_properties();

//     const u32 points = 256;
//     const u32 joints = 7;
//     const u32 p = 5;
//     const u32 ncontrol = 24;
//     const u32 nblocks = 4;
//     static_assert(points % nblocks == 0);

//     cuBspline device_pva;
//     cuGen3MultiTraj device_manip;
//     Bspline bspline(ncontrol, points, p, joints);
//     Gen3 host_manip;

//     // random task and input vector
//     Matrix task(joints, 6);
//     real amp = 10;
//     for (u32 i = 0; i < task.rows; i++)
//         for (u32 j = 0; j < task.cols; j++)
//             task(i, j) = amp * get_random();

//     Array x(joints*(ncontrol-6) + 1);
//     for (u32 i = 0; i < x.size; i++)
//         x[i] = amp * get_random();
//     x.back() = abs(x.back());

//     // compute constraints on host
//     bspline.compute_trajectory(x, task);
//     auto host_con = host_manip.constraints(bspline.traj);

//     // compute constraints on GPU
//     device_pva.init(points, joints, p, ncontrol);
//     device_pva.compute_control_and_send(x, task);
//     device_manip.init(0, points);
//     pva_constraints_kernel<<< 1, points >>>(device_pva);
//     device_pva.fetch_pva();
//     device_manip.fetch_constraints(points);

//     // Test correctness of the trajectory
//     for (int i = 0; i < (int)points; i++) {
//         for (int j = 0; j < joints; j++) {
//             CHECK((float)bspline.traj.pos(j, i) == (float)device_pva.host->traj.pos(j, i));
//             CHECK((float)bspline.traj.vel(j, i) == (float)device_pva.host->traj.vel(j, i));
//             CHECK((float)bspline.traj.acc(j, i) == (float)device_pva.host->traj.acc(j, i));
//         }
//     }

//     // host_con should be the same (ish) as device_manip.host_constraints
//     for (int i = 0; i < (int)host_con.size; i++)
//         CHECK(abs(host_con[i] - device_manip.host_constraints[i]) < 0.0006);

//     BENCHMARK("Objective function and constraints - CPU only") {
//         // random optimization vector
//         auto x = blast::random_array(device_pva.host->xlen(task), amp);
//         x.back() = std::abs(x.back());
//         // compute trajectory
//         bspline.compute_trajectory(x, task);
//         host_manip.constraints(bspline.traj);
//     };

//     // BENCHMARK("Objective function and constraints - GPU contraints and trajectory") {
//     //     // random optimization vector
//     //     auto x = blast::random_array(device_pva.host->xlen(task), amp);
//     //     x.back() = std::abs(x.back());
//     //     // compute trajectory
//     //     device_pva.compute_control_and_send(x, task);
//     //     pva_constraints_kernel<<< nblocks, points/nblocks >>>(device_pva);
//     //     cuda_check_kernel;
//     //     device_pva.fetch_pva();
//     //     device_manip.fetch_constraints(points);
//     //     cudaDeviceSynchronize();
//     // };

//     // BENCHMARK("Objective function and constraints - GPU contraints only") {
//     //     // random optimization vector
//     //     auto x = blast::random_array(device_pva.host->xlen(task), amp);
//     //     x.back() = std::abs(x.back());
//     //     // compute trajectory
//     //     device_pva.compute_control_and_send(x, task);
//     //     constraints_no_pva_kernel<<< nblocks, points/nblocks >>>(device_pva);
//     //     cuda_check_kernel;
//     //     device_manip.fetch_constraints(points);
//     //     cudaDeviceSynchronize();
//     // };

//     double in_41[41] {};
//     double out[7] {};
//     // BENCHMARK("Manipulator dynamics with MDA") {
//     //     for (int i = 0; i < (int)points; i++) {
//     //         in_41[0] = bspline.traj.pos(0, i);
//     //         in_41[1] = bspline.traj.vel(0, i);
//     //         in_41[2] = bspline.traj.acc(0, i);
//     //         in_41[3] = bspline.traj.pos(1, i);
//     //         in_41[4] = bspline.traj.vel(1, i);
//     //         in_41[5] = bspline.traj.acc(1, i);
//     //         in_41[6] = bspline.traj.pos(2, i);
//     //         in_41[7] = bspline.traj.vel(2, i);
//     //         in_41[8] = bspline.traj.acc(2, i);
//     //         in_41[9] = bspline.traj.pos(3, i);
//     //         in_41[10] = bspline.traj.vel(3, i);
//     //         in_41[11] = bspline.traj.acc(3, i);
//     //         in_41[12] = bspline.traj.pos(4, i);
//     //         in_41[13] = bspline.traj.vel(4, i);
//     //         in_41[14] = bspline.traj.acc(4, i);
//     //         in_41[15] = bspline.traj.pos(5, i);
//     //         in_41[16] = bspline.traj.vel(5, i);
//     //         in_41[17] = bspline.traj.acc(5, i);
//     //         in_41[18] = bspline.traj.pos(6, i);
//     //         in_41[19] = bspline.traj.vel(6, i);
//     //         in_41[20] = bspline.traj.acc(6, i);

//     //         dynamics_mda(in_41, out);
//     //     }
//     // };

//     // BENCHMARK("Manipulator dynamics with MDA2") {
//     //     for (int i = 0; i < (int)points; i++) {
//     //         in_41[0] = bspline.traj.pos(0, i);
//     //         in_41[1] = bspline.traj.vel(0, i);
//     //         in_41[2] = bspline.traj.acc(0, i);
//     //         in_41[3] = bspline.traj.pos(1, i);
//     //         in_41[4] = bspline.traj.vel(1, i);
//     //         in_41[5] = bspline.traj.acc(1, i);
//     //         in_41[6] = bspline.traj.pos(2, i);
//     //         in_41[7] = bspline.traj.vel(2, i);
//     //         in_41[8] = bspline.traj.acc(2, i);
//     //         in_41[9] = bspline.traj.pos(3, i);
//     //         in_41[10] = bspline.traj.vel(3, i);
//     //         in_41[11] = bspline.traj.acc(3, i);
//     //         in_41[12] = bspline.traj.pos(4, i);
//     //         in_41[13] = bspline.traj.vel(4, i);
//     //         in_41[14] = bspline.traj.acc(4, i);
//     //         in_41[15] = bspline.traj.pos(5, i);
//     //         in_41[16] = bspline.traj.vel(5, i);
//     //         in_41[17] = bspline.traj.acc(5, i);
//     //         in_41[18] = bspline.traj.pos(6, i);
//     //         in_41[19] = bspline.traj.vel(6, i);
//     //         in_41[20] = bspline.traj.acc(6, i);

//     //         dynamics_mda2(in_41, out);
//     //     }
//     //     return out;
//     // };

//     // BENCHMARK("Manipulator dynamics with MDA3") {
//     //     return dynamics_mda3(bspline.traj, in_41, out);
//     // };

//     Array in_21 = random_array(21, 1.0);
//     double out_7[7] {};
//     BENCHMARK("Manipulator dynamics with MDA reduced NoSimp Opt1") {
//         for (int i = 0; i < (int)points; i++)
//             dynamics_mda_reduct_nosimp_opt1(in_21.data, out_7);
//         return out_7;
//     };

//     // BENCHMARK("Manipulator dynamics with MDA NoSimp Opt1") {
//     //     for (int i = 0; i < (int)points; i++) {
//     //         in_41[0] = bspline.traj.pos(0, i);
//     //         in_41[1] = bspline.traj.vel(0, i);
//     //         in_41[2] = bspline.traj.acc(0, i);
//     //         in_41[3] = bspline.traj.pos(1, i);
//     //         in_41[4] = bspline.traj.vel(1, i);
//     //         in_41[5] = bspline.traj.acc(1, i);
//     //         in_41[6] = bspline.traj.pos(2, i);
//     //         in_41[7] = bspline.traj.vel(2, i);
//     //         in_41[8] = bspline.traj.acc(2, i);
//     //         in_41[9] = bspline.traj.pos(3, i);
//     //         in_41[10] = bspline.traj.vel(3, i);
//     //         in_41[11] = bspline.traj.acc(3, i);
//     //         in_41[12] = bspline.traj.pos(4, i);
//     //         in_41[13] = bspline.traj.vel(4, i);
//     //         in_41[14] = bspline.traj.acc(4, i);
//     //         in_41[15] = bspline.traj.pos(5, i);
//     //         in_41[16] = bspline.traj.vel(5, i);
//     //         in_41[17] = bspline.traj.acc(5, i);
//     //         in_41[18] = bspline.traj.pos(6, i);
//     //         in_41[19] = bspline.traj.vel(6, i);
//     //         in_41[20] = bspline.traj.acc(6, i);

//     //         dynamics_mda_nosimp_opt1(in_41, out);
//     //     }
//     //     return out;
//     // };

//     BENCHMARK("Manipulator dynamics with MDA NoSimp Opt2") {
//         for (int i = 0; i < (int)points; i++) {
//             in_41[0] = bspline.traj.pos(0, i);
//             in_41[1] = bspline.traj.vel(0, i);
//             in_41[2] = bspline.traj.acc(0, i);
//             in_41[3] = bspline.traj.pos(1, i);
//             in_41[4] = bspline.traj.vel(1, i);
//             in_41[5] = bspline.traj.acc(1, i);
//             in_41[6] = bspline.traj.pos(2, i);
//             in_41[7] = bspline.traj.vel(2, i);
//             in_41[8] = bspline.traj.acc(2, i);
//             in_41[9] = bspline.traj.pos(3, i);
//             in_41[10] = bspline.traj.vel(3, i);
//             in_41[11] = bspline.traj.acc(3, i);
//             in_41[12] = bspline.traj.pos(4, i);
//             in_41[13] = bspline.traj.vel(4, i);
//             in_41[14] = bspline.traj.acc(4, i);
//             in_41[15] = bspline.traj.pos(5, i);
//             in_41[16] = bspline.traj.vel(5, i);
//             in_41[17] = bspline.traj.acc(5, i);
//             in_41[18] = bspline.traj.pos(6, i);
//             in_41[19] = bspline.traj.vel(6, i);
//             in_41[20] = bspline.traj.acc(6, i);

//             dynamics_mda_nosimp_opt2(in_41, out);
//         }
//         return out;
//     };

//     BENCHMARK("Manipulator dynamics with MDA NoSimp Opt2 with simd trig functions") {
//         dynamics_mda_nosimp_opt2_custom(bspline.traj, in_41, out);
//         return out;
//     };

// }
#endif // nvcc
#endif // tests
