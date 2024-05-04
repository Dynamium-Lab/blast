#include "blast.h"
#include "blast_math.h"

namespace blast {

struct Wolf1  {
    blast::Array x;
};

// Solves the collision distance problem with only one box at a time, returns best fitness score
real collision_gwo(const Matrix &robot_cartesian_positions, Box box, const int n_wolves, const int n_iterations) {
    // Initialization of GWO parameters
    const auto n_dimensions = 2;
    Array x_alpha(n_dimensions);
    Array x_beta(n_dimensions);
    Array x_delta(n_dimensions);

    real alpha_fit = INF_REAL;
    real beta_fit = INF_REAL;
    real delta_fit = INF_REAL;
    std::vector<Wolf1> wolf(n_wolves);
    // Initialization of wolves' positions
    for (int i = 4; i < n_wolves; i++) {
        wolf[i].x.resize(2);
        real fraction = (n_dimensions == 1) ? i / (n_dimensions) : i / (n_dimensions - 1);
        for (int j = 0; j < n_dimensions; j++) {
            wolf[i].x[j] = fraction;
            wolf[i].x[j] = 0.5*get_random() + 0.5;
        }
    }
    wolf[0].x.resize(2);
    wolf[0].x[0] = 0;
    wolf[0].x[1] = 0;
    
    wolf[1].x.resize(2);
    wolf[1].x[0] = 0;
    wolf[1].x[1] = 1;

    wolf[2].x.resize(2);
    wolf[2].x[0] = 1;
    wolf[2].x[1] = 0;

    wolf[3].x.resize(2);
    wolf[3].x[0] = 1;
    wolf[3].x[1] = 1;

    // GWO Computation
    for (int k = 0; k < n_iterations; k++) {
        real a = 2.0 - k * (2.0 / n_iterations);
        for (int i = 0; i < n_wolves; i++) {
            // Evaluating fitness and deciding which wolves are Alpha, Beta and Delta
            real fitness = obj_function(wolf[i].x, robot_cartesian_positions, box);
            // Alpha
            x_alpha = (fitness < alpha_fit) ? wolf[i].x : x_alpha;
            alpha_fit = (fitness < alpha_fit) ? fitness : alpha_fit; 
            
            // Beta
            x_beta = (fitness > alpha_fit && fitness < beta_fit) ? wolf[i].x : x_beta;
            beta_fit = (fitness > alpha_fit && fitness < beta_fit) ? fitness : beta_fit; 

            // Delta
            x_delta = (fitness > alpha_fit && fitness > beta_fit && fitness < delta_fit) ? wolf[i].x : x_delta;
            delta_fit = (fitness > alpha_fit && fitness > beta_fit && fitness < delta_fit) ? fitness : delta_fit; 
            
            // Search for the Prey (Diverging the wolves from each other)
            real a_alpha = a * ((2*(0.5 * get_random() + 0.5))-1);
            real a_beta = a * ((2*(0.5 * get_random() + 0.5))-1);
            real a_delta = a * ((2*(0.5 * get_random() + 0.5))-1);
                
            // Emphasising exploration 
            // Effects of obstacles to approaching prey in nature
            // Works like a weight to make it harder for the wolf to reach the prey
            real c_alpha = 2* (0.5 * get_random() + 0.5);
            real c_beta = 2* (0.5 * get_random() + 0.5);
            real c_delta = 2* (0.5 * get_random() + 0.5);
                
            for (int j = 0; j < n_dimensions; j++) {
                // Search Agent updates its position
                real d_alpha = std::abs(c_alpha * x_alpha[j] - wolf[i].x[j]);;
                real d_beta = std::abs(c_beta * x_beta[j] - wolf[i].x[j]);;
                real d_delta = std::abs(c_delta * x_delta[j] - wolf[i].x[j]);

                real X1 = (x_alpha[j] - a_alpha * d_alpha);
                real X2 = (x_beta[j] - a_beta * d_beta);
                real X3 = (x_delta[j] - a_delta * d_delta);

                wolf[i].x[j]= (X1 + X2 + X3) / 3;
            }
            wolf[i].x = clamp(wolf[i].x, 0, 1);
        }
    }
    real best_f = obj_function(x_alpha, robot_cartesian_positions, box);
    return best_f;
}

// Calls collision_gwo to solve collision distance problem one box at a time, one member at a time. Returns best fitness score
real test_collision_gwo_box(const Matrix &robot_cartesian_positions, World* world, const int n_wolves, const int n_iterations) {
    int n_caps = robot_cartesian_positions.rows/3 - 1;
    int n_points = robot_cartesian_positions.cols;
    real min_dist = INF_REAL;
    real temp_dist;
    Matrix temp(6, n_points);
    for (int j = 0; j < n_caps; j++) {
        for (int i = 0; i < n_points; i++) {
            temp(0, i) = robot_cartesian_positions(j*3, i); 
            temp(1, i) = robot_cartesian_positions(j*3+1, i); 
            temp(2, i) = robot_cartesian_positions(j*3+2, i); 
            temp(3, i) = robot_cartesian_positions(j*3+3, i); 
            temp(4, i) = robot_cartesian_positions(j*3+4, i); 
            temp(5, i) = robot_cartesian_positions(j*3+5, i); 
        } 
        for (int i = 0; i < world->boxes.size(); i++) {
            temp_dist = collision_gwo(temp, world->boxes[i], n_wolves, n_iterations);
            min_dist = temp_dist < min_dist ? temp_dist : min_dist;
        }
    }
    return min_dist;
}

// Solves the collision distance problem with the full world, returns best fitness score
real collision_gwo(const Matrix &robot_cartesian_positions, const World* world, const int n_wolves, const int n_iterations) {
    // Initialization of GWO parameters
    const auto n_dimensions = 2;
    Array x_alpha(n_dimensions);
    Array x_beta(n_dimensions);
    Array x_delta(n_dimensions);

    real alpha_fit = INF_REAL;
    real beta_fit = INF_REAL;
    real delta_fit = INF_REAL;
    std::vector<Wolf1> wolf(n_wolves);
    // Initialization of wolves' positions
    for (int i = 4; i < n_wolves; i++) {
        wolf[i].x.resize(2);
        real fraction = (n_dimensions == 1) ? i / (n_dimensions) : i / (n_dimensions - 1);
        for (int j = 0; j < n_dimensions; j++) {
            wolf[i].x[j] = fraction;
            wolf[i].x[j] = 0.5*get_random() + 0.5;
        }
    }
    wolf[0].x.resize(2);
    wolf[0].x[0] = 0;
    wolf[0].x[1] = 0;

    wolf[1].x.resize(2);
    wolf[1].x[0] = 0;
    wolf[1].x[1] = 1;

    wolf[2].x.resize(2);
    wolf[2].x[0] = 1;
    wolf[2].x[1] = 0;

    wolf[3].x.resize(2);
    wolf[3].x[0] = 1;
    wolf[3].x[1] = 1;

    // GWO Computation
    for (int k = 0; k < n_iterations; k++) {
        real a = 2.0 - k * (2.0 / n_iterations);
        for (int i = 0; i < n_wolves; i++) {
            // Evaluating fitness and deciding which wolves are Alpha, Beta and Delta
            real fitness = obj_function(wolf[i].x, robot_cartesian_positions, world);
            // Alpha
            x_alpha = (fitness < alpha_fit) ? wolf[i].x : x_alpha;
            alpha_fit = (fitness < alpha_fit) ? fitness : alpha_fit; 
            
            // Beta
            x_beta = (fitness > alpha_fit && fitness < beta_fit) ? wolf[i].x : x_beta;
            beta_fit = (fitness > alpha_fit && fitness < beta_fit) ? fitness : beta_fit; 

            // Delta
            x_delta = (fitness > alpha_fit && fitness > beta_fit && fitness < delta_fit) ? wolf[i].x : x_delta;
            delta_fit = (fitness > alpha_fit && fitness > beta_fit && fitness < delta_fit) ? fitness : delta_fit; 
            
            // Hunting
            // Approaching the prey 
            
            // Search for the Prey (Diverging the wolves from each other)
            real a_alpha = a * ((2*(0.5 * get_random() + 0.5))-1);
            real a_beta = a * ((2*(0.5 * get_random() + 0.5))-1);
            real a_delta = a * ((2*(0.5 * get_random() + 0.5))-1);
                
            // Emphasising exploration 
            // Effects of obstacles to approaching prey in nature
            // Works like a weight to make it harder for the wolf to reach the prey
            real c_alpha = 2* (0.5 * get_random() + 0.5);
            real c_beta = 2* (0.5 * get_random() + 0.5);
            real c_delta = 2* (0.5 * get_random() + 0.5);
                
            for (int j = 0; j < n_dimensions; j++) {
                // Search Agent updates its position
                real d_alpha = std::abs(c_alpha * x_alpha[j] - wolf[i].x[j]);;
                real d_beta = std::abs(c_beta * x_beta[j] - wolf[i].x[j]);;
                real d_delta = std::abs(c_delta * x_delta[j] - wolf[i].x[j]);

                real X1 = (x_alpha[j] - a_alpha * d_alpha);
                real X2 = (x_beta[j] - a_beta * d_beta);
                real X3 = (x_delta[j] - a_delta * d_delta);

                wolf[i].x[j]= (X1 + X2 + X3) / 3;
            }
            wolf[i].x = clamp(wolf[i].x, 0, 1);
        }
    }
    real best_f = obj_function(x_alpha, robot_cartesian_positions, world);
    return best_f;
}

// Calls collision_gwo to solve collision distance problem with the full world, one member at a time. Returns best fitness score
real test_collision_gwo_world_1caps(const Matrix &robot_cartesian_positions, const World* world, const int n_wolves, const int n_iterations) {
    int n_caps = robot_cartesian_positions.rows / 3 - 1;
    real distmin = INF_REAL;
    real current_dist;
    Matrix temp(6, robot_cartesian_positions.cols);
    for (int j = 0; j < n_caps; j++) {
        for (int i = 0; i < robot_cartesian_positions.cols; i++) {
            temp(0, i) = robot_cartesian_positions(j * 3, i);
            temp(1, i) = robot_cartesian_positions(j * 3 + 1, i);
            temp(2, i) = robot_cartesian_positions(j * 3 + 2, i);
            temp(3, i) = robot_cartesian_positions(j * 3 + 3, i);
            temp(4, i) = robot_cartesian_positions(j * 3 + 4, i);
            temp(5, i) = robot_cartesian_positions(j * 3 + 5, i);
        }
        current_dist = collision_gwo(temp, world, n_wolves, n_iterations);
        distmin = (distmin < 0) ? (current_dist > distmin ? current_dist : distmin) :
                                (current_dist < distmin ? current_dist : distmin);
    }
    return distmin;
}

// Calls collision_gwo to solve collision distance problem with the full world, all members at once. Returns best fitness score
real test_collision_gwo_world_full_robot(const Matrix &robot_cartesian_positions, World* world, int n_wolves, int n_iterations) {
    real distmin = collision_gwo(robot_cartesian_positions, world, n_wolves, n_iterations);
    return distmin;
}

} // namespace blast