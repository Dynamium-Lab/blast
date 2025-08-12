#pragma once

#include <blast>
#include "nlopt.h"

namespace blast {

// Create n pseudo-random arrays of d dimensions form : d-by-n
// (https://extremelearning.com.au/unreasonable-effectiveness-of-quasirandom-sequences/#GeneralizingGoldenRatio)
inline blast::Matrix lds_golden_ratio(const u32 n, const u32 d) {
    blast::Matrix result(d, n);

    blast::real g = 1.6180339887498948482; // Recommanded by wikipedia article

    // precompute factors a
    std::vector<blast::real> a;
    a.reserve(d);
    a.push_back(1.0 / g);
    for (u32 i = 1; i < d; i++)
        a.push_back(a.at(i-1) * a.at(0));

    for (u32 i = 0; i < n; i++) {
        for (u32 j = 0; j < d; j++) {
            result(j, i) = fmod(0.5 + a[j]*(i+1), 1);
        }
    }

    return result;
}

Array guess_random(Bspline& bspline, Matrix& task) {
    Array x(bspline.x_len(task));
    fill_random(x, 1);
    x.back() = std::abs(x.back()) * 5 + 0.1;
    return x;
}

inline Matrix get_random_guesses(const u32 n, const u32 d) {
    Array x(d);
    Matrix result(d, n);
    for (int i = 0; i < n; i++) {
        fill_random(x, 1);
        x.back() = std::abs(x.back()) * 5 + 0.1;
        for (int j = 0; j < d; j++) {
            result(j, i) = x[j];
        }
    }
    return result;
}

// todo: Remove collisions
inline Array guess_shot_mean(Optimization* opt) { // no collisions
    Array best_x(opt->bspline.x_len(opt->task));
    real best_val = INF_REAL;
    for (int idx_nshot = 0; idx_nshot < opt->guess.n_shot; idx_nshot++) {
        auto x = guess_random((opt->bspline), opt->task);
        opt->bspline.compute_trajectory(x, opt->task);
        Array c1(opt->constraints.n_constraints);
        compute_constraints(c1.data, x, opt);
        real r = 0;
        for (u32 i = 0; i < c1.size; i++)
            r += std::max({c1[i], 0.0});
        Assert( ! isnan(r));
        r = r*x.back(); // todo: Evaluate time estimate impact on trajectory
        if (r < best_val) {
            best_x = x;
            best_val = r;
        }
    }
    return best_x;
}

inline Array get_best_x(Optimization* opt) {
    real best_f = INF_REAL;
    Array best_x(opt->guess.candidates.rows);
    for (u32 i = 0; i < opt->guess.candidates.cols; i++) {
        auto current_x = opt->guess.candidates.col(i);
        auto current_f = compute_objective(current_x, opt);

        best_x = current_f < best_f ? current_x : best_x;
        best_f = current_f < best_f ? current_f : best_f;
    }

    return best_x;
}

inline Array init_guess(Optimization* opt) {
    Array x(opt->bspline.x_len(opt->task));
    switch(opt->guess.type) {
    case Guess::rrt_connect: {
        // auto x = guess_rrt<og::RRTConnect>(opt); todo: create and use opt.guess.parameter for range
        std::cout << "RRTConnect is not yet supported" << std::endl;
        break;
    }
    case Guess::random: {
        x = guess_shot_mean(opt);
        break;
    }
    case Guess::custom: {
        x = opt->guess.x0;
        break;
    }
    case Guess::from_list: {
        x = get_best_x(opt);
        break;
    }
    default:
        Assert(false);
    }
    return x;
}

} // namespace blast