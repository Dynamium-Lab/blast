#pragma once
#include <blast>

namespace blast {

struct R2 {
  const int joints = 2;

  Vec3 p_base;
  Vec3 dv[2]; // vector to next joint
  Vec3 ev[2]; // direction vectors of joint
  real m[2];  // link masses
  Mat3 I[2];  // Inertial tensors
  Vec3 av[2]; // centers of mass
  Vec3 sv[2]; // centers of mass from next joint

  R2();
  Matrix dynamics(const Trajectory& traj);
  void   dynamics(const Trajectory& traj, Matrix& result);
};

R2::R2() {
  p_base = {0, 0, 0};
  // link mass
  m[0] = 1;
  m[1] = 2;

  // inertial tensors
  // I[0] = {0.01f, 0, 0, 0, 0.01f, 0, 0, 0, 0.01f};
  // I[1] = {0.01f, 0, 0, 0, 0.01f, 0, 0, 0, 0.01f};

  I[0] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
  I[1] = {0, 0, 0, 0, 0, 0, 0, 0, 0};

  // center of mass
  av[0] = {0.1, 0, 0};
  av[1] = {0.1, 0, 0};

  // vector to next joint
  dv[0] = {0.2, 0, 0};
  dv[1] = {0.2, 0, 0};

  // unit joint direction
  ev[0] = {0, 0, 1};
  ev[1] = {0, 0, 1};

  // center of mass from next joint
  sv[0] = -dv[0] + av[0];
  sv[1] = -dv[1] + av[1];
}

Matrix R2::dynamics(const Trajectory& traj) {
  using std::cout;
  using std::endl;
  const auto points = traj.pos.cols;
  const auto joints = traj.pos.rows;

  Matrix result(joints * 6, points);

  Mat3 Q1, Q2;
  Mat3 Q1t, Q2t;
  Vec3 w1, w2;
  Vec3 wd1, wd2;
  Vec3 cdd0 = {0, 0, (real) 9.81};
  Vec3 cdd1, cdd2;
  Vec3 f1, f2;
  Vec3 n1, n2;

  // loop all points
  for (u32 i = 0; i < points; i++) {
    auto p = traj.pos.col(i);
    auto v = traj.vel.col(i);
    auto a = traj.acc.col(i);

    Array s(2);
    Array c(2);
    blast::sincos(p, s, c);

    // note: these are stored column-wise
    Q1  = {c[0], s[0], 0, -s[0], c[0], 0, 0, 0, 1};
    Q2  = {c[1], 0, s[1], -s[1], 0, c[1], 0, -1, 0};
    Q1t = transpose(Q1);
    Q2t = transpose(Q2);

    // note: This is the Newton algorithm in 'Element de robotique' course notes.
    //       Careful because some variables are named differently and uses a slightly different conventions.
    //       For example, the ith coordinate frame turns with the ith joint, where in the course notes, the
    //       joint turns with respect to the coordinate frame.
    //-- kinematics
    w1 = v[0] * ev[0];
    w2 = Q2t * w1 + v[1] * ev[1];

    wd1  = a[0] * ev[0];
    cdd1 = Q1t * cdd0 + cross(wd1, av[0]) + cross(w1, cross(w1, av[0]));

    wd2  = Q2t * wd1 + a[1] * ev[1] + v[1] * cross(Q2t * w1, ev[1]);
    cdd2 = Q2t * cdd1 + cross(wd2, av[1]) + cross(w2, cross(w2, av[1])) - Q2t * cross(wd1, sv[0]) - Q2t * cross(w1, cross(w1, sv[0]));

    //-- dynamics
    f2 = m[1] * cdd2;
    n2 = I[1] * wd2 + cross(w2, I[1] * w2) + cross(av[1], f2);

    f1 = m[0] * cdd1 + Q2 * f2;
    n1 = I[0] * wd1 + cross(w1, I[0] * w1) + Q2 * n2 + cross(av[0], f1) - cross(sv[0], (Q2 * f2));
#if 0
        cout << "w1 "; print(w1);
        cout << "wd1 "; print(wd1);
        cout << "cdd1 "; print(cdd1);
        cout << "f1 "; print(f1);
        cout << "n1 "; print(n1);

        cout << "w2 "; print(w2);
        cout << "wd2 "; print(wd2);
        cout << "cdd2 "; print(cdd2);
        cout << "f2 "; print(f2);
        cout << "n2 "; print(n2);
#endif
    result(0, i) = f1.x;
    result(1, i) = f1.y;
    result(2, i) = f1.z;
    result(3, i) = n1.x;
    result(4, i) = n1.y;
    result(5, i) = n1.z;

    result(6, i)  = f2.x;
    result(7, i)  = f2.y;
    result(8, i)  = f2.z;
    result(9, i)  = n2.x;
    result(10, i) = n2.y;
    result(11, i) = n2.z;
  }
  return result;
}

void R2::dynamics(const Trajectory& traj, Matrix& result) {
  using std::cout;
  using std::endl;
  const auto points = traj.pos.cols;
  const auto joints = traj.pos.rows;

  Mat3 Q1, Q2;
  Mat3 Q1t, Q2t;
  Vec3 w1, w2;
  Vec3 wd1, wd2;
  Vec3 cdd0 = {0, 0, (real) 9.81};
  Vec3 cdd1, cdd2;
  Vec3 f1, f2;
  Vec3 n1, n2;

  // loop all points
  for (u32 i = 0; i < points; i++) {
    auto p = traj.pos.col(i);
    auto v = traj.vel.col(i);
    auto a = traj.acc.col(i);

#if 0
        // note: this is absurdly slow on linux
        // todo: check on windows
        Array s(2);
        Array c(2);
        blast::sincos(p, s, c);
#elif 1
    real s[] = {std::sin(p[0]), std::sin(p[1])};
    real c[] = {std::cos(p[0]), std::cos(p[1])};
#else
    __m128d s_tmp;
    __m128d c_tmp;
    __m128d angle = _mm_loadu_pd(p.data);
    s_tmp         = _mm_sincos_pd(&c_tmp, angle);
    real s[2];
    real c[2];
    _mm_store_pd(s, s_tmp);
    _mm_store_pd(c, c_tmp);
#endif

    // note: these are stored column-wise
    Q1  = {c[0], s[0], 0, -s[0], c[0], 0, 0, 0, 1};
    Q2  = {c[1], 0, s[1], -s[1], 0, c[1], 0, -1, 0};
    Q1t = transpose(Q1);
    Q2t = transpose(Q2);

    // note: This is the Newton algorithm in 'Element de robotique' course notes.
    //       Careful because some variables are named differently and uses a slightly different conventions.
    //       For example, the ith coordinate frame turns with the ith joint, where in the course notes, the
    //       joint turns with respect to the coordinate frame.
    //-- kinematics
    w1 = v[0] * ev[0];
    w2 = Q2t * w1 + v[1] * ev[1];

    wd1  = a[0] * ev[0];
    cdd1 = Q1t * cdd0 + cross(wd1, av[0]) + cross(w1, cross(w1, av[0]));

    wd2  = Q2t * wd1 + a[1] * ev[1] + v[1] * cross(Q2t * w1, ev[1]);
    cdd2 = Q2t * cdd1 + cross(wd2, av[1]) + cross(w2, cross(w2, av[1])) - Q2t * cross(wd1, sv[0]) - Q2t * cross(w1, cross(w1, sv[0]));

    //-- dynamics
    f2 = m[1] * cdd2;
    n2 = I[1] * wd2 + cross(w2, I[1] * w2) + cross(av[1], f2);

    f1 = m[0] * cdd1 + Q2 * f2;
    n1 = I[0] * wd1 + cross(w1, I[0] * w1) + Q2 * n2 + cross(av[0], f1) - cross(sv[0], (Q2 * f2));
#if 0
        cout << "w1 "; print(w1);
        cout << "wd1 "; print(wd1);
        cout << "cdd1 "; print(cdd1);
        cout << "f1 "; print(f1);
        cout << "n1 "; print(n1);

        cout << "w2 "; print(w2);
        cout << "wd2 "; print(wd2);
        cout << "cdd2 "; print(cdd2);
        cout << "f2 "; print(f2);
        cout << "n2 "; print(n2);
#endif
    result(0, i) = f1.x;
    result(1, i) = f1.y;
    result(2, i) = f1.z;
    result(3, i) = n1.x;
    result(4, i) = n1.y;
    result(5, i) = n1.z;

    result(6, i)  = f2.x;
    result(7, i)  = f2.y;
    result(8, i)  = f2.z;
    result(9, i)  = n2.x;
    result(10, i) = n2.y;
    result(11, i) = n2.z;
  }
}

} // namespace blast
