#include "blast.h"

namespace blast {

Trajectory compute_5order_trajectory(real T, Matrix &task) {
    const u32 joints = task.rows;
    const u32 points = (u32)ceil(T * 1000 + 1);

    Trajectory result(points, joints);

    const real delta_T = 0.001f;
    T = (points - 1) * delta_T;

    Matrix A(6, joints);

    for (u32 j = 0; j < joints; j++) {
        // todo: Fix for trajectory that has NaN values
        const auto p0 = task(j, 0);
        const auto v0 = task(j, 1) * T;
        const auto a0 = task(j, 2) * T * T;
        const auto pf = task(j, 3);
        const auto vf = task(j, 4) * T;
        const auto af = task(j, 5) * T * T;

        A(0, j) = p0;
        A(1, j) = v0;
        A(2, j) = a0 * 0.5f;
        A(3, j) = 0.5f * af - 1.5f * a0 - 10 * (p0 - pf) - 6 * v0 - 4 * vf;
        A(4, j) = 1.5f * a0 - af + 15 * (p0 - pf) + 8 * v0 + 7 * vf;
        A(5, j) = 0.5f * (af - a0) - 6 * (p0 - pf) - 3 * (v0 + vf);
    }

    for (u32 i = 0; i < points; i++) {
        const auto s = (real)i / (real)(points - 1);
        const auto s2 = s * s;
        const auto s3 = s2 * s;
        const auto s4 = s3 * s;
        const auto s5 = s4 * s;
        result.t[i] = i * delta_T;
        for (u32 j = 0; j < joints; j++) {
            const auto a = &A(0, j);

            const auto q = a[0] + a[1] * s + a[2] * s2 + a[3] * s3 + a[4] * s4 + a[5] * s5;
            const auto qd = a[1] + 2 * a[2] * s + 3 * a[3] * s2 + 4 * a[4] * s3 + 5 * a[5] * s4;
            const auto qdd = 2 * a[2] + 6 * a[3] * s + 12 * a[4] * s2 + 20 * a[5] * s3;

            result.pos(j, i) = q;
            result.vel(j, i) = qd / T;
            result.acc(j, i) = qdd / (T * T);
        }
    }

    return result;
}

}