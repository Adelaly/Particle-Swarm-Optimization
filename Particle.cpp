
// Member-function definitions for class Particle

#include "Particle.h"

#include <iostream>
#include <vector>   // prototype for container vector
#include <utility>  // prototype for pair
#include <cstdlib>  // prototype for rand

using namespace std;

    // Constructor
    Particle::Particle(int dimension, const long double search_space_input[][2], const long double initialization_space_input[][2])
    {
        this->dimension = dimension;
        position = new long double[dimension];
        velocity = new long double[dimension];
        individual_best = new long double[dimension];
        lbest = new long double[dimension];

        for(int i=0; i<dimension; i++)
        {
            //  Search Space
            long double lower_bound = search_space_input[i][0];
            long double upper_bound = search_space_input[i][1];
            pair<long double, long double> space1 = make_pair(lower_bound, upper_bound);
            search_space.push_back(space1);

            //  Initialization Space
            long double initialization_lower_bound = initialization_space_input[i][0];
            long double initialization_upper_bound = initialization_space_input[i][1];
            pair<long double, long double> space2 = make_pair(initialization_lower_bound, initialization_upper_bound);
            initialization_space.push_back(space2);

            //  Half-diff methodology (to initialize velocity)
            long double random_position = unifRand(initialization_lower_bound, initialization_upper_bound);
            long double random_other_position = unifRand(lower_bound, upper_bound);

            //  Position
            position[i] = random_position;
            //  Velocity
            velocity[i] = 0.5*(random_position-random_other_position);
            //  Individual best
            individual_best[i] = random_position;
            //  Local best
            lbest[i] = random_position;

            //  Other
            out_of_bound = false;   // Only used in Bound Handling Technique-----Infinity
        }
    }

    // Set functions
    void Particle::set_position_at(int at, long double value)
    {
        position[at] = value;
    }

    void Particle::set_velocity_at(int at, long double value)
    {
        velocity[at] = value;
    }

    void Particle::set_individual_best(const long double* better_individual_best)
    {
        for(int i=0; i<dimension; i++)
            individual_best[i] = better_individual_best[i];
    }

    void Particle::set_local_best(const long double* better_local_best)
    {
        for(int i=0; i<dimension; i++)
            lbest[i] = better_local_best[i];
    }

    void Particle::set_search_space(long double search_upper_bound, long double search_lower_bound)
    {
        for(int i=0; i<dimension; i++)
        {
            (search_space[i]).first = search_lower_bound;
            (search_space[i]).second = search_upper_bound;
        }
    }

    void Particle::set_initialization_space(long double initialization_upper_bound, long double initialization_lower_bound)
    {
        for(int i=0; i<dimension; i++)
        {
            (initialization_space[i]).first = initialization_lower_bound;
            (initialization_space[i]).second = initialization_upper_bound;
        }
    }

    void Particle::set_out_of_bound(bool out_of_bound)  // Only used in Bound Handling Technique-----Infinity
    {
        this->out_of_bound = out_of_bound;
    }

    void Particle::random_shuffle_particle()        //  random shuffle positions and velocities of the particle
    {
        for(int i=0; i<dimension; i++)
        {
            //  Half-diff methodology (to initialize velocity)
            long double random_position = unifRand((initialization_space[i]).first, (initialization_space[i]).second);
            long double random_other_position = unifRand((search_space[i]).first, (search_space[i]).second);

            //  Position
            position[i] = random_position;
            //  Velocity
            velocity[i] = 0.5*(random_position-random_other_position);
            //  Individual best
            individual_best[i] = random_position;
            //  Local best
            lbest[i] = random_position;

            //  Other
            out_of_bound = false;   // Only used in Bound Handling Technique-----Infinity
        }
    }

    // Retrieve functions
    int Particle::get_dimension() const
    {
        return dimension;
    }

    long double Particle::get_upper_bound_at(int at) const
    {
        return (search_space[at]).second;
    }

    long double Particle::get_lower_bound_at(int at) const
    {
        return (search_space[at]).first;
    }

    long double Particle::get_initialization_upper_bound(int at) const
    {
        return (initialization_space[at]).second;
    }

    long double Particle::get_initialization_lower_bound(int at) const
    {
        return (initialization_space[at]).first;
    }

    long double Particle::get_position_at(int at) const
    {
        return position[at];
    }

    long double Particle::get_velocity_at(int at) const
    {
        return velocity[at];
    }

    long double Particle::get_individual_best_at(int at) const
    {
        return individual_best[at];
    }

    long double Particle::get_local_best_at(int at) const
    {
        return lbest[at];
    }

    long double* Particle::get_position() const
    {
        return position;
    }

    long double* Particle::get_velocity() const
    {
        return velocity;
    }

    long double* Particle::get_local_best() const
    {
        return lbest;
    }

    long double* Particle::get_individual_best() const
    {
        return individual_best;
    }

    bool Particle::get_out_of_bound() const  // Only used in Bound Handling Technique-----Infinity
    {
        return out_of_bound;
    }

    // Display data
    void Particle::display() const
    {
        cout << "position vector: " << endl;
        for(int i=0; i<dimension; i++)
            cout << position[i] << " ";

        cout << "\nvelocity vector: " << endl;
        for(int i=0; i<dimension; i++)
            cout << velocity[i] << " ";

        cout << "\nindividual best: " << endl;
        for(int i=0; i<dimension; i++)
            cout << individual_best[i] << " ";

        cout << "\nlocal best: " << endl;
        for(int i=0; i<dimension; i++)
            cout << lbest[i] << " ";
        cout << "\n\n\n";
    }

    void Particle::display_individual_best() const
    {
        for(int i=0; i<dimension; i++)
        {
            cout << individual_best[i] << " ";
            if((i+1)%5==0)
                cout << endl;
        }
    }

    // Create random numbers
    long double Particle::unifRand()
    {
        return rand() / (long double)(RAND_MAX);
    }

    long double Particle::unifRand(long double a, long double b)
    {
        return (b-a)*unifRand() + a;
    }
