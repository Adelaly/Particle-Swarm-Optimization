
// Declaration of class Particle

#ifndef PARTICLE_H_INCLUDED
#define PARTICLE_H_INCLUDED

#include <vector>   // prototype for container vector
#include <utility>  // prototype for pair

using namespace std;


class Particle
{
public:

    //  Constructor
    Particle(int dimension, const long double search_space_input[][2], const long double initialization_space_input[][2]);

    //  Modify data
    void set_position_at(int at, long double value);
    void set_velocity_at(int at, long double value);
    void set_individual_best(const long double* better_individual_best);
    void set_local_best(const long double* better_local_best);
    void set_search_space(long double search_upper_bound, long double search_lower_bound);
    void set_initialization_space(long double initialization_upper_bound, long double initialization_lower_bound);
    void set_out_of_bound(bool out_of_bound);    // Only used in Bound Handling Technique-----Infinity
    void random_shuffle_particle();

    //  Retrieve data
    int get_dimension() const;
    long double get_upper_bound_at(int at) const;
    long double get_lower_bound_at(int at) const;
    long double get_initialization_upper_bound(int at) const;
    long double get_initialization_lower_bound(int at) const;
    long double get_position_at(int at) const;
    long double get_velocity_at(int at) const;
    long double get_individual_best_at(int at) const;
    long double get_local_best_at(int at) const;
    long double* get_position() const;
    long double* get_velocity() const;
    long double* get_individual_best() const;
    long double* get_local_best() const;
    bool get_out_of_bound() const;   // Only used in Bound Handling Technique-----Infinity and Infinity_clamp

    //  Display data
    void display() const;
    void display_individual_best() const;

private:

    //  Create random numbers
    long double unifRand();  // Return uniform random value between 0-1 inclusively.
    long double unifRand(long double a, long double b);    // Return uniform random value between a-b inclusively.

    // Data members
    int dimension;
    vector<pair<long double, long double> > search_space; //  Search space
    vector<pair<long double, long double> > initialization_space; //  Initialization space
    long double* position;    //  Position
    long double* velocity;    //  Velocity
    long double* individual_best; //  Individual best(or Personal best)
    long double* lbest;   //  Local best
    bool out_of_bound;  // Only used in Bound Handling Technique------Infinitiy

};

#endif // PARTICLE_H_INCLUDED
