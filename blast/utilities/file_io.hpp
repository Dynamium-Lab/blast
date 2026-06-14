#pragma once
#include <blast>

namespace blast {
inline host_fn void print_to_csv(const Matrix& m, const std::string& filename) {
  std::ofstream file;
  file.open(filename);
  for (u32 i = 0; i < m.rows; i++) {
    for (u32 j = 0; j < m.cols - 1; j++)
      file << m(i, j) << ",";
    file << m(i, m.cols - 1) << std::endl;
  }
  file.close();
}

inline host_fn void print_to_csv(const Trajectory& traj, const std::string& filename) {
  std::ofstream file;
  file.open(filename);

  // print header
  {
    file << "t";
    for (u32 i = 0; i < traj.pos.rows; i++)
      file << ",p[" << i << "]";
    for (u32 i = 0; i < traj.vel.rows; i++)
      file << ",v[" << i << "]";
    for (u32 i = 0; i < traj.acc.rows; i++)
      file << ",a[" << i << "]";
    file << std::endl;
  }

  // print data
  for (u32 i = 0; i < traj.t.size; i++) {
    file << traj.t[i];
    for (u32 j = 0; j < traj.pos.rows; j++)
      file << "," << traj.pos(j, i);
    for (u32 j = 0; j < traj.vel.rows; j++)
      file << "," << traj.vel(j, i);
    for (u32 j = 0; j < traj.acc.rows; j++)
      file << "," << traj.acc(j, i);
    file << std::endl;
  }
  file.close();
}

inline host_fn void print_to_csv(const std::vector<Trajectory>& traj, const std::string& filename) {
  std::ofstream file;
  file.open(filename);

  real start_time = 0.0;

  file << traj[0].t[0];
  for (u32 j = 0; j < traj[0].pos.rows; j++)
    file << "," << traj[0].pos(j, 0);
  for (u32 j = 0; j < traj[0].vel.rows; j++)
    file << "," << traj[0].vel(j, 0);
  for (u32 j = 0; j < traj[0].acc.rows; j++)
    file << "," << traj[0].acc(j, 0);
  file << std::endl;

  for (int task_id = 0; task_id < traj.size(); task_id++) {
    // int points_more = (int) std::ceil(res[task_id].x.back() * 1000.0) + 1;


    for (u32 i = 1; i < traj[task_id].t.size; i++) {
      file << traj[task_id].t[i] + start_time;
      for (u32 j = 0; j < traj[task_id].pos.rows; j++)
        file << "," << traj[task_id].pos(j, i);
      for (u32 j = 0; j < traj[task_id].vel.rows; j++)
        file << "," << traj[task_id].vel(j, i);
      for (u32 j = 0; j < traj[task_id].acc.rows; j++)
        file << "," << traj[task_id].acc(j, i);
      file << std::endl;
    }

    start_time += traj[task_id].t[traj[task_id].t.size - 1];
  }
  file.close();
}

inline host_fn blast::Matrix read_csv_matrix(const std::string& filename, const char* csv_sep = ",") {
  std::ifstream file(filename);
  std::cout << "Reading from file: " << filename << std::endl;
  if (!file.is_open()) {
    Assert(false);
    std::cerr << "ERROR: File is unavailable" << std::endl;
    return blast::Matrix(0, 0);
  }

  std::string line;
  u32         num_rows = 0, num_cols = 0;

  std::getline(file, line); // discard first line
  // Determine number of rows and columns
  if (std::getline(file, line)) {
    num_cols = std::count(line.begin(), line.end(), *csv_sep) + 1;
    num_rows++;
  }
  while (std::getline(file, line)) {
    num_rows++;
  }

  // Rewind file to read again
  file.clear();
  file.seekg(0);

  // Create matrix with detected dimensions
  blast::Matrix pos(num_rows, num_cols);

  std::getline(file, line); // discard first line
  // Read the data into the matrix
  for (u32 i = 0; i < num_rows; i++) {
    std::getline(file, line);
    std::istringstream iss(line);
    for (u32 j = 0; j < num_cols; j++) {
      std::string value;
      if (std::getline(iss, value, *csv_sep)) {
        pos(i, j) = std::stod(value);
      }
    }
  }
  return pos;
}

inline host_fn blast::Matrix read_csv_matrix_no_header(const std::string& filename, const char* csv_sep = ",") {
  std::ifstream file(filename);
  std::cout << "Reading from file: " << filename << std::endl;
  if (!file.is_open()) {
    Assert(false);
    std::cerr << "ERROR: File is unavailable" << std::endl;
    return blast::Matrix(0, 0);
  }

  std::string line;
  u32         num_rows = 0, num_cols = 0;

  // Determine number of rows and columns
  if (std::getline(file, line)) {
    num_cols = std::count(line.begin(), line.end(), *csv_sep) + 1;
    num_rows++;
  }
  while (std::getline(file, line)) {
    num_rows++;
  }

  // Rewind file to read again
  file.clear();
  file.seekg(0);

  // Create matrix with detected dimensions
  blast::Matrix pos(num_rows, num_cols);

  // Read the data into the matrix
  for (u32 i = 0; i < num_rows; i++) {
    std::getline(file, line);
    std::istringstream iss(line);
    for (u32 j = 0; j < num_cols; j++) {
      std::string value;
      if (std::getline(iss, value, *csv_sep)) {
        pos(i, j) = std::stod(value);
      }
    }
  }

  return pos;
}

inline host_fn Trajectory read_csv_trajectory_no_header(const std::string& filename, const char* csv_sep = ",") {
  std::ifstream file(filename);
  std::cout << "Reading from file: " << filename << std::endl;
  if (!file.is_open()) {
    Assert(false);
    std::cerr << "ERROR: File is unavailable" << std::endl;
    return {0, 0};
  }

  std::string line;
  u32         num_rows = 0, num_cols = 0;

  // Determine number of rows and columns
  if (std::getline(file, line)) {
    num_cols = std::count(line.begin(), line.end(), *csv_sep) + 1;
    num_rows++;
  }
  while (std::getline(file, line)) {
    num_rows++;
  }

  // Rewind file to read again
  file.clear();
  file.seekg(0);

  u32        num_joints = (num_cols - 1) / 3;
  Trajectory trajectory(num_rows, num_joints);

  // Read the data into the matrix
  for (u32 i = 0; i < num_rows; i++) {
    std::getline(file, line);
    std::istringstream iss(line);
    u32                col = 0;

    // Gets time
    {
      std::string value;
      if (std::getline(iss, value, *csv_sep))
        trajectory.t[i] = std::stod(value);
    }

    // Gets position
    for (col = 0; col < num_joints; col++) {
      std::string value;
      if (std::getline(iss, value, *csv_sep)) {
        trajectory.pos(col, i) = std::stod(value);
      }
    }

    // Gets velocity
    for (col = 0; col < num_joints; col++) {
      std::string value;
      if (std::getline(iss, value, *csv_sep)) {
        trajectory.vel(col, i) = std::stod(value);
      }
    }

    // Gets acceleration
    for (col = 0; col < num_joints; col++) {
      std::string value;
      if (std::getline(iss, value, *csv_sep)) {
        trajectory.acc(col, i) = std::stod(value);
      }
    }
  }

  return trajectory;
}

} // namespace blast
