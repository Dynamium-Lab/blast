#pragma once

namespace blast {


real obj_function(const Array& x, const Matrix &robot_cartesian_positions, const Box &box);
real obj_function(const Array& x, const Matrix &robot_cartesian_positions, const World* world);
real obj_function_gradient(const Array& x, const Matrix &robot_cartesian_positions, const World* world, Array* gradient);


} // namespace blast