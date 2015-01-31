
// Definition of class Society

#ifndef SOCIETY_H_INCLUDED
#define SOCIETY_H_INCLUDED

#include <vector>   // prototype for container vector
#include <string>   // prototype for string type
#include "Particle.h"

    /*  Topology of the society:
            The society is arranged in a shape of grid torus.
            The neighborhood contains five particles: ones below, above, to the right, to the left and itself.
        In the case, a 7*7 grid torus is used   */
class Society
{
public:

    // Constructor
    Society();
    Society(int dimension, const long double search_space[][2], const long double initialization_space[][2]);

    /*  Entrance of program */
    //  Test with seven tranditional benchmarks
    void one_one_start();            //  Test one bound handling technique with one function
    void all_all_start();            //  Test all bound handling technique with all function
    void all_one_start();            //  Test all bound handling technique with one function
    void one_all_start();            //  Test all bound handling technique with one function
    //  Test with twenty five CEC2005 benchmarks (Need comment or uncomment several codes)
    void cec2005_one_one_start();    //  Test one bound handling technique with one function
    void cec2005_all_one_start();    //  Test all bound handling technique with all function
    //  Test with flat landscape
    void flat_landscape_start();     //  Test all bound handling technique with flat landscape

    //-------------------------------------------------------------------------
    //  Update
    void update_lbest();    // Update local best position of all particles in the society
    void update_individual_best();   // Update the individual best position of all particles in the society
    void update_velocity_and_position();    // Update the velocity and position of all particles in the society
    void update_velocity_and_position_no_bound();   // For no bounds benchmarks, update velocity and position

    void update_lbest_flat_landscape();
    void update_individual_best_flat_landscape();

    //  Evaluate
    long double evaluate_member(const Particle* member) const;   // Evaluate the position of designated member
    long double evaluate_position(long double* position) const;   // Evaluate the designated position

    //  Display
    void display_society() const;
    int gbest_id() const;
    void display_gbest() const;
    void display_gbest_evaluation() const;
    void record_basic_info() const;
    void record_step() const;
    void record() const;

    //  Set
    void set_technique(int chosen_bound_handling_technique);    // Set bound handling technique.
    void set_technique();   // Set bound handling technique.
    void set_velocity_updating_method(int chosen_velocity_updating_method); // Set velocity update function.
    void set_velocity_updating_method();    // Set velocity update function.
    void set_benchmark_function(int chosen_benchmark_function); // Set benchmarks.
    void set_benchmark_function();  // Set benchmarks.
    void set_cec2005_benchmark_function();
    void set_dimension(int dimension);  // Set dimensions.
    void set_dimension();   // Set dimensions.
    void set_times(int times);  // Set number of generations. (for calculate Inertia Weight of old velocity)
    void set_total_times(int total_times);  // Set total number of generations.
    void set_total_times(); // Set total number of generations.
    void set_total_experimental_times();
    void set_society();     // Set social members with data members.
    void modify_spaces();  // Modify the search space and initializing space of social members.
    void random_shuffle_society();  // Random shuffle all particles' position and velocity

    //  Retrieve
    int get_technique() const;
    int get_velocity_updating_method() const;
    int get_benchmark_function() const;
    int get_dimension() const;
    int get_times() const;
    int get_total_times() const;
    string get_benchmark_function_name() const;
    string get_bound_handling_name() const;

private:

    //  Basic rules for velociety and position updating
    long double velocity_iterator(long double old_velocity, long double old_position, long double lbest, long double individual_best);
    long double position_iterator(long double old_velocity, long double old_position);
    long double velocity_clamp_operation(long double velocity, int dimension_id) const;   // Clamp used to limit maximal velocity

    /*  Benchmarks  */
    //  Seven traditional benchmark functions
    long double sphere_function(long double* position) const;
    long double rosenbrock_function(long double* position) const;
    long double rastrigin_function(long double* position) const;
    long double griewank_function(long double* position) const;
    long double ackley_function(long double* position) const;
    long double michalewicz_function(long double* position) const;
    long double schwefel_function(long double* position) const;
    long double cec2005_functions(long double* position) const;

    /*  Boundary Handling Techniques (Four types) */
    //  A. Change the personal/global(local) best update
    bool infinity(long double position, int dimension_id);

    //  B. Reposition infeasible particles
    pair<long double,long double> nearest_zero(long double position, long double velocity, int dimension_id);
    pair<long double,long double> nearest_deterministic_back(long double position, long double velocity, int dimension_id);
    pair<long double,long double> nearest_random_back(long double position, long double velocity, int dimension_id);
    long double nearest_unmodified(long double position, int dimension_id);   //  Return position

    pair<long double,long double> reflect_zero(long double position, long double velocity, int dimension_id);
    pair<long double,long double> reflect_deterministic_back(long double position, long double velocity, int dimension_id);
    pair<long double,long double> reflect_random_back(long double position, long double velocity, int dimension_id);
    long double reflect_unmodified(long double position, int dimension_id);   //  Return position

    pair<long double,long double> random_zero(long double position, long double velocity, int dimension_id);
    pair<long double,long double> random_deterministic_back(long double position, long double velocity, int dimension_id);
    pair<long double,long double> random_random_back(long double position, long double velocity, int dimension_id);
    long double random_unmodified(long double position, int dimension_id);    //  Return position

    //  C. Preventing infeasibility
    long double hyperbolic(long double old_position, long double velocity, int dimension_id);  //  Return velocity
    long double random_forth(long double old_position, long double velocity, int dimension_id);    //  Return velocity

    //  D. Other approaches
    pair<long double,long double> bounded_mirror(long double position, long double velocity, int dimension_id) ;
    long double period(long double position, int dimension_id);

    //  Create random numbers
    long double unifRand() const;      // Return uniform random value between 0-1 inclusively.
    long double unifRand(long double a, long double b) const;    // Return uniform random value between a-b inclusively.
    //  Regulate input from the customer
    int input_converter(string your_input) const;

    //-------------------------------------------------------------------------------
    //  Data members
    vector<Particle*> society;
    int chosen_bound_handling_technique;
    int chosen_velocity_updating_method;
    int chosen_benchmark_function;
    int dimension;  // dimension of the particles.
    int times;  // Used in one velocity updating method, to calculate Inertia Weight W
    int total_times;    // Used in one velocity updating method, to calculate Inertia Weights W
    int experimental_times;
    int total_experimental_times;

};

#endif // SOCIETY_H_INCLUDED
