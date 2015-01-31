
// Member-function definitions for class Society

# include "Society.h"
# include "cec2005_benchmarks.h"    // CEC2005 fuctions
# include "rand.h"                  // Random number generator for CEC2005 functions
# include "Particle.h"

# include <iostream>
# include <vector>   // prototype for vector type
# include <utility>  // prototype for pair type
# include <cstdlib>  // prototype for rand
# include <cmath>    // prototupe for fmod, pow
# include <sstream>  // for converting int to string
# include <string>   // prototype for string type, it included in sstream
# include <limits>   // prototupe for numeric_limit
# include <fstream>  // prototype for ofstream
# include <iomanip>  // prototype for setw


using namespace std;

// Const variables
// grid_leagth is the const variable for modifying topology of the society
const int grid_length = 7;

// Const variables served as constraint factors to velocity.
//  X acts as friction (to movement of the particles); C1 and C2 control
//   the relative attraction to the global best and personal best.
const long double X = 0.729843788;
const long double C1 = 2.05;
const long double C2 = 2.05;
const long double velocity_clamp_value = 10000;

// Const Variables served as constraint factors to velocity.
//  W acts as Inertia Weight which is not a global variable for data's security
//  W shifts from 0.9 to 04 linearly according to the iteration number of the search.
//  A_C1 and A_C2 are acceleration coefficients corresponding to local best and individual best position
const long double Inertia_Weight_Upper = 0.9;
const long double Inertia_Weight_Lower = 0.4;
const long double A_C1 = 2.0;
const long double A_C2 = 2.0;

const int the_record_step = 100;

    /* The topology of the society is arranged in a shape of grid torus.
       The neighborhood contains five particles:
        ones below, above, to the right, to the left and itself. */
    // A 7*7 grid torus
    Society::Society()
    {
        chosen_benchmark_function = 0;
        times = 0;
        total_experimental_times = 1;
    }

    Society::Society(int dimension, const long double search_space[][2], const long double initialization_space[][2])
    {
        for(int i=0; i<grid_length*grid_length; i++)
        {
            Particle* member = new Particle(dimension,search_space,initialization_space);
            society.push_back(member);
        }

        chosen_benchmark_function = 0;
        times = 0;
        total_experimental_times = 1;
    }

    //  Test one technique and Test one benchmark function(Only seven traditional benchmarks)
    void Society::one_one_start()
    {
        //  Basic initialization
        set_velocity_updating_method();
        set_technique();
        set_benchmark_function();
        set_dimension();
        set_total_times();
        set_total_experimental_times();

        set_society();  //  Creat swams, initialize spaces and initialize velocities and positions

        for(int h=1; h<=total_experimental_times; h++)
        {
            experimental_times = h;
            update_lbest();
            record_basic_info();

            //Working
            for(int i=0; i<total_times; i++)
            {
                if(times%the_record_step == 0)
                    record_step();
                cout << "times " << i << endl;
                update_velocity_and_position();
                update_individual_best();
                update_lbest();
            }
            display_society();
            display_gbest();
            record_step();

            random_shuffle_society();   //  Initialize velocities and positions
        }
    }

    //  Test all techniques and Test all benchmark functions(Only seven traditional benchmarks)
    void Society::all_all_start()
    {
        //  Basic initialization
        set_velocity_updating_method();
        set_dimension();
        set_total_times();
        set_total_experimental_times();

        set_society();  //  Creat swams, but not initialize spaces and therefore velocities and positions

        for(int i=1; i<=18; i++)
        {
            set_technique(i);

            for(int j=1; j<=7; j++)
            {
                set_benchmark_function(j);
                modify_spaces();   //  Initialize spaces

                for(int h=1; h<=total_experimental_times; h++)
                {
                    experimental_times = h;
                    random_shuffle_society();   // Initialize velocities and positions

                    update_lbest();
                    record_basic_info();
                    //  Working
                    for(int k=0; k<total_times; k++)
                    {
                        if(times%the_record_step == 0)
                            record_step();
                        update_velocity_and_position();
                        update_individual_best();
                        update_lbest();
                    }
                    display_gbest_evaluation();
                    record_step();
                }
            }
        }
    }

    //  Test all techniques and  Test one benchmark function(Only seven traditional benchmarks)
    void Society::all_one_start()
    {
        //  Basic initialization
        set_velocity_updating_method();
        set_benchmark_function();
        set_dimension();
        set_total_times();
        set_total_experimental_times();

        set_society();  //  Creat swarm, initialize spaces and initialize velocities and positions

        for(int i=1; i<=18; i++)
        {
            set_technique(i);

            for(int h=1; h<=total_experimental_times; h++)
            {
                experimental_times = h;

                update_lbest();
                record_basic_info();
                //  Working
                for(int k=0; k<total_times; k++)
                {
                    if(times%the_record_step == 0)
                        record_step();
                    update_velocity_and_position();
                    update_individual_best();
                    update_lbest();
                }
                display_gbest_evaluation();
                //record();
                record_step();
                random_shuffle_society();   //  Initialize velocities and positions
            }
        }
    }

    //  Test one bound handling technique and Test all benchmark functions(Only seven traditional benchmarks)
    void Society::one_all_start()
    {
        //  Basic initialization
        set_velocity_updating_method();
        set_technique();
        set_dimension();
        set_total_times();
        set_total_experimental_times();

        set_society();  //  Creat swarm, but not initialize spaces and therefore velocities and positions

        for(int j=1; j<=7; j++)
        {
            set_benchmark_function(j);
            modify_spaces();   //  Initialize spaces

            for(int h=1; h<=total_experimental_times; h++)
            {
                experimental_times = h;
                random_shuffle_society();   //  Initialize velocities and positions

                update_lbest();
                record_basic_info();
                //  Working
                for(int k=0; k<total_times; k++)
                {
                    if(times%the_record_step == 0)
                        record_step();
                    update_velocity_and_position();
                    update_individual_best();
                    update_lbest();
                }
                display_gbest_evaluation();
                //record();
                record_step();
            }
        }
    }

    // Test one technique and Test one chosen CEC2005 benchmark function
    void Society::cec2005_one_one_start()
    {
        // Basic Initialization
        set_velocity_updating_method();
        set_technique();
        set_cec2005_benchmark_function();
        cout << "\nNotice: dimension can only be 2, 10, 30 or 50.\n";
        set_dimension();
        set_total_times();
        set_total_experimental_times();

        set_society();  //  Creat swarm, initialize spaces and initialize velocities and positions

        // Initialize function variables
        // Need uncomment or comment at the place
        /* Assign nfunc in the begining */
        // 1 for f1-f14, 10 for f15-f25
        nfunc = 1;
        //nfunc = 10;
        nreal = dimension;
        if (dimension!=2 && dimension!=10 && dimension!=30 && dimension!=50)
        {
            cerr << "\n Wrong value of 'nreal(dimension)' entered, only 2, 10, 30, 50 variables are supported" << endl;
            exit(0);
        }

        seed = rand()/(long double)RAND_MAX;
        /* Call these routines to initialize random number generator */
        /* require for computing noise in some test problems */
        randomize();
        initrandomnormaldeviate();

        /* nreal and nfunc need to be initialized before calling these routines */
        /* Routine to allocate memory to global variables */
        allocate_memory();

        /* Routine the initalize global variables */
        initialize();

        /* For test problems 15 to 25, we need to calculate a normalizing quantity */
        /* The line below should be uncommented only for functions 15 to 25 */
        //calc_benchmark_norm();    /* Comment this line for functions 1 to 14 */

        for(int h=1; h<=total_experimental_times; h++)
        {
            experimental_times = h;

            update_lbest();
            record_basic_info();
            //Working
            for(int i=0; i<total_times; i++)
            {
                //cout << "times " << i << endl;
                if(times%the_record_step == 0)
                    record_step();
                if(chosen_benchmark_function==107||chosen_benchmark_function==125)
                    update_velocity_and_position_no_bound();
                else
                    update_velocity_and_position();
                update_individual_best();
                update_lbest();
            }
            display_society();
            display_gbest();
            record_step();

            random_shuffle_society();   //  Initialize velocities and position
        }

        /* Routine to free the memory allocated at run time */
        free_memory();
    }

     // Test all technique and Test one chosen CEC2005 benchmark function
    void Society::cec2005_all_one_start()
    {
        // Basic Initialization
        set_velocity_updating_method();
        set_cec2005_benchmark_function();
        cout << "\nNotice: dimension can only be 2, 10, 30 or 50.\n";
        set_dimension();
        set_total_times();
        set_total_experimental_times();

        set_society();  //  Creat swarm, initialize spaces and initialize velocities and positions

        // Initialize function variables
        // Need uncomment or comment at the place
        /* Assign nfunc in the begining */
        // 1 for f1-f14, 10 for f15-f25
        nfunc = 1;
        //nfunc = 10;
        nreal = dimension;
        if (dimension!=2 && dimension!=10 && dimension!=30 && dimension!=50)
        {
            cerr << "\n Wrong value of 'nreal(dimension)' entered, only 2, 10, 30, 50 variables are supported" << endl;
            exit(0);
        }

        seed = rand()/(long double)RAND_MAX;
        /* Call these routines to initialize random number generator */
        /* require for computing noise in some test problems */
        randomize();
        initrandomnormaldeviate();

        /* nreal and nfunc need to be initialized before calling these routines */
        /* Routine to allocate memory to global variables */
        allocate_memory();

        /* Routine the initalize global variables */
        initialize();

        /* For test problems 15 to 25, we need to calculate a normalizing quantity */
        /* The line below should be uncommented only for functions 15 to 25 */
        //calc_benchmark_norm();    /* Comment this line for functions 1 to 14 */

        for(int i=1; i<=18; i++)
        {
            set_technique(i);

            for(int h=1; h<=total_experimental_times; h++)
            {
                experimental_times = h;

                update_lbest();
                record_basic_info();
                //Working
                for(int i=0; i<total_times; i++)
                {
                    if(times%the_record_step == 0)
                        record_step();
                    if(chosen_benchmark_function==107||chosen_benchmark_function==125)
                        update_velocity_and_position_no_bound();
                    else
                        update_velocity_and_position();
                    update_individual_best();
                    update_lbest();
                }
                display_society();
                display_gbest();
                record_step();

                random_shuffle_society();   //  Initialize velocities and position
            }
        }

        /* Routine to free the memory allocated at run time */
        free_memory();
    }

    void Society::flat_landscape_start()
    {
        //  Basic initialization
        set_velocity_updating_method();
        set_total_times();
        set_total_experimental_times();

        //  Creat swarm, initialize spaces and initialize velocities and positions
        set_dimension(30);
        long double search_space[30][2];
        long double initialization_space[30][2];
        for(int i=0; i<30; i++)
        {
            search_space[i][0] = -100;
            search_space[i][1] = 100;
            initialization_space[i][0] = -100;
            initialization_space[i][1] = 100;
        }
        for(int i=0; i<grid_length*grid_length; i++)
        {
            Particle* member = new Particle(30,search_space,initialization_space);
            society.push_back(member);
        }
        for(int i=1; i<=18; i++)
        {
            set_technique(i);

            for(int h=1; h<=total_experimental_times; h++)
            {
                update_lbest_flat_landscape();

                stringstream ss;
                ss << h;
                string flat_test_times = ss.str();
                string file_name = "output_data/flat_landscape/" + get_bound_handling_name()
                                    + "_" + flat_test_times + ".txt";
                ofstream output(file_name.c_str(), ios::out);
                if(!output)
                {
                    cerr << "\nFile could not be opened for writing." << endl;
                    exit(0);
                }

                int count[30][20];
                for(int m=0; m<30; m++)
                {
                    for(int n=0; n<20; n++)
                        count[m][n] = 0;
                }

                //  Working
                for(int k=0; k<total_times; k++)
                {
                    update_velocity_and_position();
                 //cout << "\nok here 8.1\n";
                    update_individual_best_flat_landscape();
                // cout << "\nok here 8.2\n";
                    update_lbest_flat_landscape();
                    //  Recording
               // cout << "\nok here 8.3\n";
                    for(int l=0; l<grid_length*grid_length; l++)
                    {
                        for(int m=0; m<30; m++)
                        {
                            //cout << "\nok here 8.4\n";
                            long double pos = (society[l])->get_position_at(m);
                           // cout << "\nok here 8.5\n";
                            if(pos<100 && pos>=-100)
                                count[m][(int)((pos+100)/10)]++;
                            else if(pos==100)
                                count[m][19]++;
                        }
                    }
                }
               // cout << "\nok here 9\n";
                //  Recording
                for(int m=0; m<30; m++)
                {
                    output << "\ndimension: " << m+1 << endl;
                    for(int n=0; n<20; n++)
                        output << (long double)(count[m][n])/(total_times*49) << endl;
                }
                output.close();

                random_shuffle_society();   //  Initialize velocities and positions
               // cout << "\nok here 10\n";
                cout << "\n\nFinished one experiment for one technique!\n" << endl;
            }
        }
    }

    //  Update
    //  Update local best position of all particles in the society
    void Society::update_lbest()
    {
        const int population = society.size();
        for(int i=0; i<population; i++)
        {
            int neighborhood[5];
            // To the left
            if(i%grid_length==0)
            {
                neighborhood[0] = i+grid_length-1;
            }
            else
            {
                neighborhood[0] = i-1;
            }
            // To the right
            if((i+1)%grid_length==0)
            {
                neighborhood[1] = i-grid_length+1;
            }
            else
            {
                neighborhood[1] = i+1;
            }
            // Above
            if(i-grid_length<0)
            {
                neighborhood[2] = i-grid_length+grid_length*grid_length;
            }
            else
            {
                neighborhood[2] = i-grid_length;
            }
            // Below
            if(i+grid_length>grid_length*grid_length-1)
            {
                neighborhood[3] = i+grid_length-grid_length*grid_length;
            }
            else
            {
                neighborhood[3] = i+grid_length;
            }
            // Itself
            neighborhood[4] = i;

            //  Find the individual in the local society having the best position in history.
            int id = neighborhood[0];

            long double assessment =  evaluate_position((*(society[id])).get_individual_best());
            for(int j=1; j<5; j++)
            {
                long double temp = evaluate_position((*(society[neighborhood[j]])).get_individual_best());
                if(temp <= assessment)
                {
                    id = neighborhood[j];
                    assessment = temp;
                }
            }
            // Modify local best position;
            (*(society[i])).set_local_best((*(society[id])).get_individual_best());
        }
    }

    void Society::update_lbest_flat_landscape()
    {
        const int population = society.size();
        for(int i=0; i<population; i++)
        {
            int neighborhood[5];
            // To the left
            if(i%grid_length==0)
            {
                neighborhood[0] = i+grid_length-1;
            }
            else
            {
                neighborhood[0] = i-1;
            }
            // To the right
            if((i+1)%grid_length==0)
            {
                neighborhood[1] = i-grid_length+1;
            }
            else
            {
                neighborhood[1] = i+1;
            }
            // Above
            if(i-grid_length<0)
            {
                neighborhood[2] = i-grid_length+grid_length*grid_length;
            }
            else
            {
                neighborhood[2] = i-grid_length;
            }
            // Below
            if(i+grid_length>grid_length*grid_length-1)
            {
                neighborhood[3] = i+grid_length-grid_length*grid_length;
            }
            else
            {
                neighborhood[3] = i+grid_length;
            }
            // Itself
            neighborhood[4] = i;

            //  Find the individual in the local society having the best position in history.
            int chosen_lbest_id = rand() % 5;

            // Modify local best position;
            (*(society[i])).set_local_best((*(society[neighborhood[chosen_lbest_id]])).get_individual_best());
        }
    }

    // Update the individual best position of all particles in the society
    void Society::update_individual_best()
    {
        const int population = society.size();
        for(int i=0; i<population; i++)
        {
            if((*(society[i])).get_out_of_bound() == true)
                continue;
            if(evaluate_member(society[i]) <= evaluate_position((*(society[i])).get_individual_best()))
                (*(society[i])).set_individual_best((*(society[i])).get_position());
        }
    }

    void Society::update_individual_best_flat_landscape()
    {
        const int population = society.size();
        for(int i=0; i<population; i++)
        {
            if((*(society[i])).get_out_of_bound() == true)
                continue;
            if(rand()%2 == 0)
                (*(society[i])).set_individual_best((*(society[i])).get_position());
        }
    }

    // Update velocity and position of all particles
    void Society::update_velocity_and_position()
    {
        const int population = society.size();
        for(int i=0; i<population; i++)
        {
            const long double* lbest = (*(society[i])).get_local_best();
            const long double* individual_best = (*(society[i])).get_individual_best();
            const long double* old_position = (*(society[i])).get_position();
            const long double* old_velocity = (*(society[i])).get_velocity();
            bool is_out_of_bound = false;   //  For techniques Infinity and Infinity_Clamp
            for(int j=0; j<dimension; j++)
            {
                // Calculate velocity and position
                long double velocity;
                long double position;

                // Boundary Handling Techniques
                pair<long double,long double> position_velocity;
                switch(get_technique())
                {
                    case 1:     //  Bounded mirror
                        velocity = velocity_iterator(old_velocity[j], old_position[j], lbest[j], individual_best[j]);
                        velocity = velocity_clamp_operation(velocity,j);
                        position = position_iterator(old_position[j], velocity);
                        position_velocity = bounded_mirror(position,velocity, j);
                        break;
                    case 2:     //  Period
                        velocity = velocity_iterator(old_velocity[j], old_position[j], lbest[j], individual_best[j]);
                        velocity = velocity_clamp_operation(velocity,j);
                        position = period(position_iterator(old_position[j], velocity), j);
                        position_velocity = make_pair(position, velocity);
                        break;
                    case 3:     //  Reflect zero
                        velocity = velocity_iterator(old_velocity[j], old_position[j], lbest[j], individual_best[j]);
                        velocity = velocity_clamp_operation(velocity,j);
                        position = position_iterator(old_position[j], velocity);
                        position_velocity = reflect_zero(position, velocity, j);
                        break;
                    case 4:     //  Hyperbolic
                        velocity = velocity_iterator(old_velocity[j], old_position[j], lbest[j], individual_best[j]);
                        velocity = velocity_clamp_operation(velocity,j);
                        velocity = hyperbolic(old_position[j], velocity, j);    //  Normalize the velocity
                        velocity = velocity_clamp_operation(velocity,j);
                        position_velocity.second = velocity;
                        position_velocity.first = position_iterator(old_position[j], velocity);
                        break;
                    case 5:     //  Reflect unmodified
                        velocity = velocity_iterator(old_velocity[j], old_position[j], lbest[j], individual_best[j]);
                        velocity = velocity_clamp_operation(velocity,j);
                        position = position_iterator(old_position[j], velocity);
                        position_velocity.second = velocity;
                        position_velocity.first = reflect_unmodified(position, j);
                        break;
                    case 6:     //  Random forth
                        velocity = velocity_iterator(old_velocity[j], old_position[j], lbest[j], individual_best[j]);
                        velocity = velocity_clamp_operation(velocity,j);
                        velocity = random_forth(old_position[j], velocity, j);
                        velocity = velocity_clamp_operation(velocity,j);
                        position_velocity.second = velocity;
                        position_velocity.first = position_iterator(old_position[j], velocity);
                        break;
                    case 7:     //  Reflect deterministic back
                        velocity = velocity_iterator(old_velocity[j], old_position[j], lbest[j], individual_best[j]);
                        velocity = velocity_clamp_operation(velocity,j);
                        position = position_iterator(old_position[j], velocity);
                        position_velocity = reflect_deterministic_back(position, velocity, j);
                        break;
                    case 8:     //  Reflect random back
                        velocity = velocity_iterator(old_velocity[j], old_position[j], lbest[j], individual_best[j]);
                        velocity = velocity_clamp_operation(velocity,j);
                        position = position_iterator(old_position[j], velocity);
                        position_velocity = reflect_random_back(position, velocity, j);
                        break;
                    case 9:     //  Nearest zero
                        velocity = velocity_iterator(old_velocity[j], old_position[j], lbest[j], individual_best[j]);
                        velocity = velocity_clamp_operation(velocity,j);
                        position = position_iterator(old_position[j], velocity);
                        position_velocity = nearest_zero(position, velocity, j);
                        break;
                    case 10:    //  Nearest deterministic back
                        velocity = velocity_iterator(old_velocity[j], old_position[j], lbest[j], individual_best[j]);
                        velocity = velocity_clamp_operation(velocity,j);
                        position = position_iterator(old_position[j], velocity);
                        position_velocity = nearest_deterministic_back(position, velocity, j);
                        break;
                    case 11:    //  Nearest random back
                        velocity = velocity_iterator(old_velocity[j], old_position[j], lbest[j], individual_best[j]);
                        velocity = velocity_clamp_operation(velocity,j);
                        position = position_iterator(old_position[j], velocity);
                        position_velocity = nearest_random_back(position, velocity, j);
                        break;
                    case 12:    //  Nearest unmodifed
                        velocity = velocity_iterator(old_velocity[j], old_position[j], lbest[j], individual_best[j]);
                        velocity = velocity_clamp_operation(velocity,j);
                        position = position_iterator(old_position[j], velocity);
                        position_velocity.second = velocity;
                        position_velocity.first = nearest_unmodified(position, j);
                        break;
                    case 13:    //  Infinity
                        velocity = velocity_iterator(old_velocity[j], old_position[j], lbest[j], individual_best[j]);
                        position = position_iterator(old_position[j], velocity);
                        if(infinity(position, j) == true)
                            is_out_of_bound = true;
                        position_velocity = make_pair(position, velocity);
                        break;
                    case 14:    //  Infinity clamp
                        velocity = velocity_iterator(old_velocity[j], old_position[j], lbest[j], individual_best[j]);
                        velocity = velocity_clamp_operation(velocity, j);
                        position = position_iterator(old_position[j], velocity);
                        if(infinity(position, j) == true)
                            is_out_of_bound = true;
                        position_velocity = make_pair(position, velocity);
                        break;
                    case 15:    //  Random zero
                        velocity = velocity_iterator(old_velocity[j], old_position[j], lbest[j], individual_best[j]);
                        velocity = velocity_clamp_operation(velocity,j);
                        position = position_iterator(old_position[j], velocity);
                        position_velocity = random_zero(position, velocity, j);
                        break;
                    case 16:    //  Random deterministic back
                        velocity = velocity_iterator(old_velocity[j], old_position[j], lbest[j], individual_best[j]);
                        velocity = velocity_clamp_operation(velocity,j);
                        position = position_iterator(old_position[j], velocity);
                        position_velocity = random_deterministic_back(position, velocity, j);
                        break;
                    case 17:    //  Random random back
                        velocity = velocity_iterator(old_velocity[j], old_position[j], lbest[j], individual_best[j]);
                        velocity = velocity_clamp_operation(velocity,j);
                        position = position_iterator(old_position[j], velocity);
                        position_velocity = random_random_back(position, velocity, j);
                        break;
                    case 18:    //  Random unmodifed
                        velocity = velocity_iterator(old_velocity[j], old_position[j], lbest[j], individual_best[j]);
                        velocity = velocity_clamp_operation(velocity,j);
                        position = position_iterator(old_position[j], velocity);
                        position_velocity.second = velocity;
                        position_velocity.first = random_unmodified(position, j);
                        break;

                    default:
                        cout << "Error Input" << endl;
                        exit(0);
                        break;
                }
                // Update position and velocity at j-dimension of i-th particle
                (*(society[i])).set_position_at(j, position_velocity.first);
                (*(society[i])).set_velocity_at(j, position_velocity.second);
                (*(society[i])).set_out_of_bound(is_out_of_bound);  //  For technique Infinity and Infinity_clamp
            }
        }
        times++; // One generation pass
    }

    // For no bounds function, update velocity and position of all particles
    void Society::update_velocity_and_position_no_bound()
    {
        const int population = society.size();
        for(int i=0; i<population; i++)
        {
            const long double* lbest = (*(society[i])).get_local_best();
            const long double* individual_best = (*(society[i])).get_individual_best();
            const long double* old_position = (*(society[i])).get_position();
            const long double* old_velocity = (*(society[i])).get_velocity();
            for(int j=0; j<dimension; j++)
            {
                // Calculate velocity and position
                long double velocity =  velocity_iterator(old_velocity[j], old_position[j], lbest[j], individual_best[j]);;
                velocity = velocity_clamp_operation(velocity,j);

                // Boundary Handling Techniques
                (*(society[i])).set_position_at(j, position_iterator(old_position[j], velocity));
                (*(society[i])).set_velocity_at(j, velocity);
            }
        }
        times++;    // One generation pass
    }

    //  Evaluate
    //  Evaluate the position of designated member
    long double Society::evaluate_member(const Particle* member) const
    {
        if(chosen_benchmark_function == 1)  //  Sphere function
            return sphere_function(member->get_position());
        else if(chosen_benchmark_function == 2) //  Rosenbrock function
            return rosenbrock_function(member->get_position());
        else if(chosen_benchmark_function == 3) //  Rastringin function
            return rastrigin_function(member->get_position());
        else if(chosen_benchmark_function == 4) //  Sphere function
            return griewank_function(member->get_position());
        else if(chosen_benchmark_function == 5) //  Rosenbrock function
            return ackley_function(member->get_position());
        else if(chosen_benchmark_function == 6)
            return michalewicz_function(member->get_position());
        else if(chosen_benchmark_function == 7)
            return schwefel_function(member->get_position());
        else if(chosen_benchmark_function >= 101 && chosen_benchmark_function <=125)
            return cec2005_functions(member->get_position());
        else
        {
            cerr << "\n Error Input " << endl;
            exit(0);
        }
    }

    // Evaluate the designated position
    long double Society::evaluate_position(long double* position) const
    {
        if(chosen_benchmark_function == 1)  //  Sphere function
            return sphere_function(position);
        else if(chosen_benchmark_function == 2) //  Rosenbrock function
            return rosenbrock_function(position);
        else if(chosen_benchmark_function == 3) //  Rastringin function
            return rastrigin_function(position);
        else if(chosen_benchmark_function == 4) //  Sphere function
            return griewank_function(position);
        else if(chosen_benchmark_function == 5) //  Rosenbrock function
            return ackley_function(position);
        else if(chosen_benchmark_function == 6)
            return michalewicz_function(position);
        else if(chosen_benchmark_function == 7)
            return schwefel_function(position);
        else if(chosen_benchmark_function >= 101 && chosen_benchmark_function <=125)
            return cec2005_functions(position);
        else
        {
            cerr << "\n Error Input " << endl;
            exit(0);
        }
    }

    //  Display
    void Society::display_society() const
    {
        int population = society.size();
        for(int i=0; i<population; i++)
        {
            (*(society[i])).display();
            cout << endl << endl << endl;
        }
    }

    int Society::gbest_id() const
    {
        const int size = society.size();
        int best_particle_id = 0;
        long double best_score = evaluate_position((*(society[0])).get_individual_best());
        for(int i=1; i<size; i++)
        {
            long double score = evaluate_position((*(society[i])).get_individual_best());
            if(score <= best_score)
            {
                best_score = score;
                best_particle_id = i;
            }
        }
        return best_particle_id;
    }

    void Society::display_gbest() const
    {
        cout << " display gbest: " << endl;
        (*(society[gbest_id()])).display();
        cout << " evaluation is: " << evaluate_position((*(society[gbest_id()])).get_individual_best());
    }

    void Society::display_gbest_evaluation() const
    {
        cout << "Technique_id: " << get_technique() << "\tBenchmark id: " << get_benchmark_function() << endl;
        const int gb_id = gbest_id();
        cout << "Best evaluation is " << evaluate_position((*(society[gb_id])).get_individual_best()) << endl;
        (society[gb_id])->display_individual_best();
         cout << "\n\n";
    }

    void Society::record_basic_info() const
    {
        string technique_name = get_bound_handling_name();
        string function_name = get_benchmark_function_name();

        stringstream ss1,ss2,ss3;
        ss1 << dimension;
        string dim = ss1.str();
        ss2 << total_times;
        string generation = ss2.str();
        ss3 << experimental_times;
        string exp_times =  ss3.str();

        string file_name = "output_data/" + function_name + "/" + technique_name + "_" + function_name + "_" + dim + "_" + generation + "_" + exp_times + ".txt";
        ofstream output(file_name.c_str(), ios::out);
        if(!output)
        {
            cerr << "\nFile could not be opened for writing." << endl;
            exit(0);
        }
        output << "Technique name: " << technique_name << endl
            <<"Function name: " << function_name << endl
            << "Dimension: " << dimension << endl
            << "Iteration: " << total_times << endl
            << "Experimental times: " << experimental_times << endl
            << "Velocity update strategy: ";
        switch(chosen_velocity_updating_method)
        {
            case 1:
                output << "V=X*[V+c1*e1(lbest-pos)+c2*e2(pbest-pos])" << endl;
                break;
            case 2:
                output << "V=w*V+ac1*e1(lbest-pos)+ac2*e2(pbest-pos)" << endl;
                break;
            default:
                cerr << "\nWrong velocity input.\n";
                exit(0);
        }
        output << "V_Clamp: " << velocity_clamp_value << endl << endl
            << "Iteration\t" << "Evaluaton" << endl;
        output.close();
        cout << "\n\nBasic record done!\n";
    }

    void Society::record_step() const
    {
        string technique_name = get_bound_handling_name();
        string function_name = get_benchmark_function_name();

        stringstream ss1,ss2,ss3;
        ss1 << dimension;
        string dim = ss1.str();
        ss2 << total_times;
        string generation = ss2.str();
        ss3 << experimental_times;
        string exp_times =  ss3.str();

        string file_name = "output_data/" + function_name + "/" + technique_name + "_" + function_name + "_" + dim + "_" + generation + "_" + exp_times + ".txt";
        ofstream output(file_name.c_str(), ios::app);
        if(!output)
        {
            cerr << "\nFile could not be opened for writing." << endl;
            exit(0);
        }
        output << times << "\t" << scientific
            << evaluate_position((*(society[gbest_id()])).get_individual_best()) << endl;
        output.close();
        cout << "\nRecord done\n";
    }

    void Society::record() const
    {
        string technique_name = get_bound_handling_name();
        string function_name = get_benchmark_function_name();

        stringstream ss1,ss2,ss3;
        ss1 << dimension;
        string dim = ss1.str();
        ss2 << total_times;
        string generation = ss2.str();
        ss3 << experimental_times;
        string exp_times =  ss3.str();

        string file_name = "output_data/" + function_name + "/" + technique_name + "_" + function_name + "_" + dim + "_" + generation + "_" + exp_times + ".txt";
        ofstream output(file_name.c_str(), ios::out);
        if(!output)
        {
            cerr << "\nFile could not be opened for writing." << endl;
            exit(0);
        }
        output << "Technique name: " << technique_name << endl
            <<"Function name: " << function_name << endl
            << "Dimension: " << dimension << endl
            << "Iteration: " << total_times << endl
            << "Experimental times: " << experimental_times << endl
            << "Velocity update strategy: ";
        switch(chosen_velocity_updating_method)
        {
            case 1:
                output << "V=X*[V+c1*e1(lbest-pos)+c2*e2(pbest-pos])" << endl;
                break;
            case 2:
                output << "V=w*V+ac1*e1(lbest-pos)+ac2*e2(pbest-pos)" << endl;
                break;
            default:
                cerr << "\nWrong velocity input.\n";
                exit(0);
        }
        output << "V_Clamp: " << velocity_clamp_value << endl
             << "Best Eval: " << scientific << evaluate_position((*(society[gbest_id()])).get_individual_best()) << endl;
        output.close();
        cout << "\n\nRecord done!\n";
    }

    //  Set
    void Society::set_technique(int chosen_bound_handling_technique)
    {
        this->chosen_bound_handling_technique = chosen_bound_handling_technique;
    }

    void Society::set_technique()
    {
        cout << "*Please select bound handling technique, and input its index.\n\n"
                << "(1) bounded_mirror\n(2) period\n(3) reflect_zero\n(4) hyperbolic\n(5) reflect_unmodified"
                << "\n(6) random_forth\n(7) reflect_deterministic_back\n(8) reflect_random_back"
                << "\n(9) nearest_zero\n(10) nearest_deterministic_back\n(11) nearest_random_back"
                << "\n(12) nearest_unmodified\n(13) infinity\n(14) infinity_clamp\n(15) random_zero"
                << "\n(16) random_deterministic_back\n(17) random_random_back\n(18) random_unmodified"
                << "\nExit? Please input others."
                << "\n\n"
                << "Your input is ";
        string your_input;
        cin >> your_input;
        int selection = input_converter(your_input);
        if(selection >= 1 && selection <= 18)
            chosen_bound_handling_technique = selection;
        else
            exit(0);
        cout << endl << endl;
    }

    void Society::set_velocity_updating_method(int chosen_velocity_updating_method)
    {
        this->chosen_velocity_updating_method = chosen_velocity_updating_method;
    }

    void Society::set_velocity_updating_method()
    {
        cout << "*Please select velocity updating method, and input its index.\n\n"
                << "(1) Velocity = X * [Old Velocity + C1 * e1 * (Local Best - Old Position) + "
                << "C2 * e2 * (Individual Best - Old Position)]\n"
                << "X = 0.729843788\nC1 = 2.05\nC2 = 2.05\ne1,e2 uniform random value of 0 to 1\n\n"
                << "(2) Velocity = W * Old Velocity + A_C1 * e1 *(Local Best - Old Position) + "
                << "A_C2 * e2 * (Individual Best - Old Position)]\n"
                << "W shifts linearly from " << Inertia_Weight_Upper << " to " << Inertia_Weight_Lower << "."
                << "\nA_C1 = 2.0\nA_C2 = 2.0\ne1,e2 uniform random value of 0 to 1"
                << "\n\n Exit? Please input others."
                << "\n\n"
                << "Your input is ";
        string your_input;
        cin >> your_input;
        int selection = input_converter(your_input);
        if(selection >= 1 && selection <= 2)
            chosen_velocity_updating_method = selection;
        else
            exit(0);
        cout << endl << endl;
    }

    void Society::set_benchmark_function(int chosen_benchmark_function)
    {
        this->chosen_benchmark_function = chosen_benchmark_function;
    }

    void Society::set_benchmark_function()
    {
        cout << "*Please Select Benchmark Function, and input its index.\n\n"
                <<"(1) Sphere function\n"<< "(2) Rosenbrock function\n"
                <<"(3) Rastrigin function\n" << "(4) Griewank function\n"
                <<"(5) Ackley function\n" << "(6) Michalewicz function\n"
                <<"(7) Schwefel function"
                <<"\n Exit? Please input others."
                <<"\n\n"
                << "Your input is ";
        string your_input;
        cin >> your_input;
        int selection = input_converter(your_input);
        if(selection >= 1 && selection <= 7)
            chosen_benchmark_function = selection;
        else
            exit(0);
        cout << endl << endl;
    }

    void Society::set_cec2005_benchmark_function()
    {
        cout << "*Please Select CEC2005 Benchmark Function Suite F1-F25, and input 1 to 25.\n\n"
                <<" Exit? Please input others."
                <<"\n\n"
                << "Your input is ";
        string your_input;
        cin >> your_input;
        int selection = input_converter(your_input);
        if(selection >= 1 && selection <= 25)
            chosen_benchmark_function = selection+100;
        else
            exit(0);
        cout << endl << endl;
    }

    void Society::set_dimension(int dimension)
    {
        this->dimension = dimension;
    }

    void Society::set_dimension()
    {
        cout << "*Please dicide how many dimension should be dealt with. Your input is ";
        bool finish_set;
        do{
            string your_input;
            cin >> your_input;
            int selection = input_converter(your_input);
            if(selection >= 1 && selection <= 200)
            {
                finish_set = true;
                dimension = selection;
            }
            else
            {
                finish_set = false;
                cout << "Input must be legal and less than 200. If you want dimension larger than 200, "
                        <<"please change the code in function set_society and modify_society.\n"
                        <<"Your input is ";
            }
        }while(!finish_set);
        cout << endl << endl;
    }

    void Society::set_times(int times)
    {
        this->times = times;
    }

    void Society::set_total_times(int total_times)
    {
        this->total_times = total_times;
    }

    void Society::set_total_times()
    {
        cout << "*Please decide how many iteraton should go on. Your input is ";
        bool finish_set;
        do{
            string your_input;
            cin >> your_input;
            int selection = input_converter(your_input);
            if(selection >= 0)
            {
                finish_set = true;
                total_times = selection;
            }
            else
            {
                finish_set = false;
                cout << "Input must be legal and not less than 0. Your input is ";
            }
        }while(!finish_set);
        cout << endl << endl;
    }

    void Society::set_total_experimental_times()
    {
        cout << "*Please decide how many experimental times should go on. Your input is ";
        bool finish_set;
        do{
            string your_input;
            cin >> your_input;
            int selection = input_converter(your_input);
            if(selection >= 0)
            {
                finish_set = true;
                total_experimental_times = selection;
            }
            else
            {
                finish_set = false;
                cout << "Input must be legal and not less than 0. Your input is ";
            }
        }while(!finish_set);
        cout << endl << endl;
    }

    //  Default method for creating function with respect to benchmark function
    void Society::set_society()
    {
        long double search_upper_bound;
        long double search_lower_bound;
        long double initialize_upper_bound;
        long double initialize_lower_bound;
        switch(chosen_benchmark_function)
        {
            case 1:     //  Sphere function
                search_lower_bound = -100;
                search_upper_bound = 100;
                initialize_lower_bound = 50;
                initialize_upper_bound = 100;
                break;
            case 2:     //  Rosenbrock function
                search_lower_bound = -30;
                search_upper_bound = 30;
                initialize_lower_bound = 15;
                initialize_upper_bound = 30;
                break;
            case 3:
                search_lower_bound = -5.12;
                search_upper_bound = 5.12;
                initialize_lower_bound = 2.56;
                initialize_upper_bound = 5.12;
                break;
            case 4:     //  Sphere function
                search_lower_bound = -600;
                search_upper_bound = 600;
                initialize_lower_bound = 300;
                initialize_upper_bound = 600;
                break;
            case 5:     //  Rosenbrock function
                search_lower_bound = -32;
                search_upper_bound = 32;
                initialize_lower_bound = 16;
                initialize_upper_bound = 32;
                break;
            case 6:
                search_lower_bound = 0;
                search_upper_bound = 3.14;
                initialize_lower_bound = 2.355;
                initialize_upper_bound = 3.14;
                break;
            case 7:
                search_lower_bound = -500;
                search_upper_bound = 500;
                initialize_lower_bound = -250;
                initialize_upper_bound = 250;
                break;
            case 101:
            case 102:
            case 103:
            case 104:
            case 105:
            case 106:
            case 114:
                search_lower_bound = -100;
                search_upper_bound = 100;
                initialize_lower_bound = -100;
                initialize_upper_bound = 100;
                break;
            case 107:
                search_lower_bound = numeric_limits<long double>::min();
                search_upper_bound = numeric_limits<long double>::max();
                initialize_lower_bound = 0;
                initialize_upper_bound = 600;
                break;
            case 108:
                search_lower_bound = -32;
                search_upper_bound = 32;
                initialize_lower_bound = -32;
                initialize_upper_bound = 32;
                break;
            case 109:
            case 110:
            case 113:
            case 115:
            case 116:
            case 117:
            case 118:
            case 119:
            case 120:
            case 121:
            case 122:
            case 123:
            case 124:
                search_lower_bound = -5;
                search_upper_bound = 5;
                initialize_lower_bound = -5;
                initialize_upper_bound = 5;
                break;
            case 111:
                search_lower_bound = -0.5;
                search_upper_bound = 0.5;
                initialize_lower_bound = -0.5;
                initialize_upper_bound = 0.5;
                break;
            case 112:
                search_lower_bound = -PI;
                search_upper_bound = PI;
                initialize_lower_bound = -PI;
                initialize_upper_bound = PI;
                break;
            case 125:
                search_lower_bound = numeric_limits<long double>::min();
                search_upper_bound = numeric_limits<long double>::max();
                initialize_lower_bound = 2;
                initialize_upper_bound = 5;
                break;
            default:
                search_lower_bound = 0;
                search_upper_bound = 0;
                initialize_lower_bound = 0;
                initialize_upper_bound = 0;
                break;
        }

        long double search_space[200][2];
        long double initialization_space[200][2];
        if(dimension > 200)
        {
            cout << "Dimension is too large. Please change the code in function set_society and modify_society\n"
                    << "Or make the input of dimension less than 200." << endl;
            exit(0);
        }
        for(int i=0; i<dimension; i++)
        {
            search_space[i][0] = search_lower_bound;
            search_space[i][1] = search_upper_bound;
            initialization_space[i][0] = initialize_lower_bound;
            initialization_space[i][1] = initialize_upper_bound;
        }
        for(int i=0; i<grid_length*grid_length; i++)
        {
            Particle* member = new Particle(dimension,search_space,initialization_space);
            society.push_back(member);
        }
    }

    void Society::modify_spaces()
    {
        long double search_upper_bound;
        long double search_lower_bound;
        long double initialize_upper_bound;
        long double initialize_lower_bound;
        switch(chosen_benchmark_function)
        {
            case 1:     //  Sphere function
                search_lower_bound = -100;
                search_upper_bound = 100;
                initialize_lower_bound = 50;
                initialize_upper_bound = 100;
                break;
            case 2:     //  Rosenbrock function
                search_lower_bound = -30;
                search_upper_bound = 30;
                initialize_lower_bound = 15;
                initialize_upper_bound = 30;
                break;
            case 3:
                search_lower_bound = -5.12;
                search_upper_bound = 5.12;
                initialize_lower_bound = 2.56;
                initialize_upper_bound = 5.12;
                break;
            case 4:     //  Sphere function
                search_lower_bound = -600;
                search_upper_bound = 600;
                initialize_lower_bound = 300;
                initialize_upper_bound = 600;
                break;
            case 5:     //  Rosenbrock function
                search_lower_bound = -32;
                search_upper_bound = 32;
                initialize_lower_bound = 16;
                initialize_upper_bound = 32;
                break;
            case 6:
                search_lower_bound = 0;
                search_upper_bound = 3.14;
                initialize_lower_bound = 2.355;
                initialize_upper_bound = 3.14;
                break;
            case 7:
                search_lower_bound = -500;
                search_upper_bound = 500;
                initialize_lower_bound = -250;
                initialize_upper_bound = 250;
                break;
            default:
                cout << "Error Input" << endl;
                exit(0);
                break;
        }

        for(int i=0; i<grid_length*grid_length; i++)
        {
            (society[i])->set_search_space(search_upper_bound,search_lower_bound);
            (society[i])->set_initialization_space(initialize_upper_bound,initialize_lower_bound);
        }
    }

    void Society::random_shuffle_society()
    {
        times = 0;
        const int size = society.size();
        for(int i=0; i<size; i++)
            (*(society[i])).random_shuffle_particle();
    }

    int Society::get_technique() const
    {
        return chosen_bound_handling_technique;
    }

    int Society::get_velocity_updating_method() const
    {
        return chosen_velocity_updating_method;
    }

    int Society::get_dimension() const
    {
        return dimension;
    }

    int Society::get_times() const
    {
        return times;
    }

    int Society::get_total_times() const
    {
        return total_times;
    }

    int Society::get_benchmark_function() const
    {
        return chosen_benchmark_function;
    }

    string Society::get_benchmark_function_name() const
    {
        if(chosen_benchmark_function==1)
            return "Sphere";
        else if(chosen_benchmark_function==2)
            return "Rosenbrock";
        else if(chosen_benchmark_function==3)
            return "Rastringin";
        else if(chosen_benchmark_function==4)
            return "Griewank";
        else if(chosen_benchmark_function==5)
            return "Ackley";
        else if(chosen_benchmark_function==6)
            return "Michalewicz";
        else if(chosen_benchmark_function==7)
            return "Schwefel";
        else if(chosen_benchmark_function>=101&&chosen_benchmark_function<=125)
        {
            int intValue = chosen_benchmark_function - 100;
            stringstream ss;
            ss << intValue;
            string str = ss.str();
            return "F" + str;
        }
    }

    string Society::get_bound_handling_name() const
    {
        if(chosen_bound_handling_technique==1)
            return "Bounded_mirror";
        else if(chosen_bound_handling_technique==2)     //  Period
            return "Period";
        else if(chosen_bound_handling_technique==3)
            return "Reflect_zero";
        else if(chosen_bound_handling_technique==4)
            return "Hyperboloic";
        else if(chosen_bound_handling_technique==5)
            return "Reflect_unmodifed";
        else if(chosen_bound_handling_technique==6)
            return "Random_forth";
        else if(chosen_bound_handling_technique==7)
            return "Reflect_deterministic_back";
        else if(chosen_bound_handling_technique==8)
            return "Reflect_random_back";
        else if(chosen_bound_handling_technique==9)
            return "Nearest_zero";
        else if(chosen_bound_handling_technique==10)
            return "Nearest_deterministic_back";
        else if(chosen_bound_handling_technique==11)    //  Reflect random back
            return "Nearest_random_back";
        else if(chosen_bound_handling_technique==12)
            return "Nearest_unmodifed";
        else if(chosen_bound_handling_technique==13)
            return "Infinity";
        else if(chosen_bound_handling_technique==14)
            return "Infinity_clamp";
        else if(chosen_bound_handling_technique==15)
            return "Random_zero";
        else if(chosen_bound_handling_technique==16)
            return "Random_deterministic_back";
        else if(chosen_bound_handling_technique==17)
            return "Random_random_back";
        else if(chosen_bound_handling_technique==18)
            return "Random_unmodified";
        else
        {
            cerr << "\n Error of chosen_benchmark_function " << endl;
            exit(0);
        }
    }
    //-------------------------------------------------
    //  Private functions

    //  Basic rules for velocity and position updating
    long double Society::velocity_iterator(long double old_velocity, long double old_position, long double lbest, long double individual_best)
    {
        long double e1 = unifRand();
        long double e2 = unifRand();

        // Calculate velocity
        long double velocity;
        switch(get_velocity_updating_method())
        {
            case 1:
                velocity = X * (old_velocity + C1*e1*(lbest-old_position) +
                            C2*e2*(individual_best-old_position));
                break;
            case 2:
            {
                long double W = Inertia_Weight_Upper - (Inertia_Weight_Upper - Inertia_Weight_Lower) *
                            (get_times() / get_total_times());
                velocity = W * old_velocity + A_C1 * e1 * (lbest - old_position) +
                            A_C2 * e2 * (individual_best - old_position);
                break;
            }

            default:
                cout << "Error Input" << endl;
                exit(0);
                break;
        }
        return velocity;
    }

    long double Society::position_iterator(long double new_velocity, long double old_position)
    {
        return old_position + new_velocity;
    }

    long double Society::velocity_clamp_operation(long double velocity, int dimension_id) const
    {
        const long double upper_bound = (*(society[0])).get_upper_bound_at(dimension_id);
        const long double lower_bound = (*(society[0])).get_lower_bound_at(dimension_id);
        if(velocity > velocity_clamp_value * (upper_bound-lower_bound))
                return velocity_clamp_value * (upper_bound-lower_bound);
        else if( velocity < (-1)*velocity_clamp_value * (upper_bound-lower_bound))
                return (-1)*velocity_clamp_value * (upper_bound-lower_bound);
        return velocity;
    }

    //  Benchmark functions

    long double Society::sphere_function(long double* position) const
    {
        const int count = dimension;
        long double function_value = 0;
        for(int i=0; i<count; i++)
            function_value+=position[i]*position[i];
        return function_value;
    }

    long double Society::rosenbrock_function(long double* position) const
    {
        const int count = dimension;
        long double function_value = 0;
        for(int i=0; i<count-1; i++)
            function_value+=100 * (pow(position[i+1]-pow(position[i],2), 2)) + pow(1-position[i],2);
        return function_value;
    }

    long double Society::rastrigin_function(long double* position) const
    {
        const int count = dimension;
        long double function_value = 10 * count;
        for(int i=0; i<count; i++)
            function_value+=pow(position[i],2) - 10*cos(2*PI*position[i]);
        return function_value;
    }

    long double Society::griewank_function(long double* position) const
    {
        const int count = dimension;
        long double sum_part = 0;
        long double product_part = 1;
        for(int i=0; i<count; i++)
        {
            sum_part+=(position[i]*position[i])/4000;
            product_part*=cos(position[i]/sqrt(i+1));
        }
        long double function_value = sum_part - product_part + 1;
        return function_value;
    }

    long double Society::ackley_function(long double* position) const
    {
        const int count = dimension;
        long double first_part = 0;
        long double second_part = 0;
        for(int i=0; i<count; i++)
        {
            first_part+=position[i]*position[i];
            second_part+=cos(2*PI*position[i]);
        }
        long double function_value = -20*exp(-0.2*sqrt(first_part/count))
                                    - exp(second_part/count) + 20 + exp(1);
        return function_value;
    }

    long double Society::michalewicz_function(long double* position) const
    {
        const long double m = 10;
        const int count = dimension;
        long double function_value = 0;
        for(int i=0; i<count; i++)
            function_value-=sin(position[i])*pow(sin((i+1)*pow(position[i],2)/PI), 2*m);
        return function_value;
    }

    long double Society::schwefel_function(long double* position) const
    {
        const int count = dimension;
        long double function_value = 0;
        for(int i=0; i<count; i++)
            function_value-=position[i] * sin(sqrt(abs(position[i])));
        return function_value;
    }

    long double Society::cec2005_functions(long double* position) const
    {
        return calc_benchmark_func(position);
    }

    //  Bound handling techniques

    //  A. Change the personal/global(local) best update
    bool Society::infinity(long double position, int dimension_id)
    {
        long double upper_bound = (*(society[0])).get_upper_bound_at(dimension_id);
        long double lower_bound = (*(society[0])).get_lower_bound_at(dimension_id);
        if(position > upper_bound || position < lower_bound)
            return true;
        else
            return false;
    }

    //  B. Reposition infeasible particles
    pair<long double,long double> Society::nearest_zero(long double position, long double velocity, int dimension_id)
    {
        long double upper_bound = (*(society[0])).get_upper_bound_at(dimension_id);
        long double lower_bound = (*(society[0])).get_lower_bound_at(dimension_id);

        if(position > upper_bound)
        {
            position = upper_bound;
            velocity = 0;
        }
        else if(position < lower_bound)
        {
            position = lower_bound;
            velocity = 0;
        }
        return make_pair(position, velocity);
    }

    pair<long double,long double> Society::nearest_deterministic_back(long double position, long double velocity, int dimension_id)
    {
        long double upper_bound = (*(society[0])).get_upper_bound_at(dimension_id);
        long double lower_bound = (*(society[0])).get_lower_bound_at(dimension_id);

        if(position > upper_bound)
        {
            position = upper_bound;
            const long double deterministic_back_constant = 0.5;
            velocity = (-1) * deterministic_back_constant * velocity;
        }
        else if(position < lower_bound)
        {
            position = lower_bound;
            const long double deterministic_back_constant = 0.5;
            velocity = (-1) * deterministic_back_constant * velocity;
        }
        return make_pair(position, velocity);
    }

    pair<long double,long double> Society::nearest_random_back(long double position, long double velocity, int dimension_id)
    {
        long double upper_bound = (*(society[0])).get_upper_bound_at(dimension_id);
        long double lower_bound = (*(society[0])).get_lower_bound_at(dimension_id);

        if(position > upper_bound)
        {
            position = upper_bound;
            velocity = (-1) * unifRand() * velocity;
        }
        else if(position < lower_bound)
        {
            position = lower_bound;
            velocity = (-1) * unifRand() * velocity;
        }
        return make_pair(position, velocity);
    }

    long double Society::nearest_unmodified(long double position, int dimension_id)
    {
        long double upper_bound = (*(society[0])).get_upper_bound_at(dimension_id);
        long double lower_bound = (*(society[0])).get_lower_bound_at(dimension_id);

        if(position > upper_bound)
        {
            position = upper_bound;
        }
        else if(position < lower_bound)
        {
            position = lower_bound;
        }
        return position;
    }

    pair<long double,long double> Society::reflect_zero(long double position, long double velocity, int dimension_id)
    {
        long double upper_bound = (*(society[0])).get_upper_bound_at(dimension_id);
        long double lower_bound = (*(society[0])).get_lower_bound_at(dimension_id);

        while(position > upper_bound || position < lower_bound)
        {
            if(position > upper_bound)
                position = upper_bound - (position - upper_bound);
            else
                position = lower_bound + (lower_bound - position);
            velocity = 0;
        }
        return make_pair(position, velocity);
    }

    pair<long double,long double> Society::reflect_deterministic_back(long double position, long double velocity, int dimension_id)
    {
        long double upper_bound = (*(society[0])).get_upper_bound_at(dimension_id);
        long double lower_bound = (*(society[0])).get_lower_bound_at(dimension_id);

        while(position > upper_bound || position < lower_bound)
        {
            if(position > upper_bound)
                position = upper_bound - (position - upper_bound);
            else
                position = lower_bound + (lower_bound - position);
            const long double deterministic_back_constant = 0.5;
            velocity = (-1) * deterministic_back_constant * velocity;
        }
        return make_pair(position, velocity);
    }

    pair<long double,long double> Society::reflect_random_back(long double position, long double velocity, int dimension_id)
    {
        long double upper_bound = (*(society[0])).get_upper_bound_at(dimension_id);
        long double lower_bound = (*(society[0])).get_lower_bound_at(dimension_id);

        while(position > upper_bound || position < lower_bound)
        {
            if(position > upper_bound)
                position = upper_bound - (position - upper_bound);
            else
                position = lower_bound + (lower_bound - position);
            velocity = (-1) * unifRand() * velocity;
        }
        return make_pair(position, velocity);
    }

    long double Society::reflect_unmodified(long double position, int dimension_id)
    {
        long double upper_bound = (*(society[0])).get_upper_bound_at(dimension_id);
        long double lower_bound = (*(society[0])).get_lower_bound_at(dimension_id);

        while(position > upper_bound || position < lower_bound)
        {
            if(position > upper_bound)
                position = upper_bound - (position - upper_bound);
            else
                position = lower_bound + (lower_bound - position);
        }
        return position;
    }

    pair<long double,long double> Society::random_zero(long double position, long double velocity, int dimension_id)
    {
        long double upper_bound = (*(society[0])).get_upper_bound_at(dimension_id);
        long double lower_bound = (*(society[0])).get_lower_bound_at(dimension_id);

        if(position > upper_bound || position < lower_bound)
        {
            position = unifRand(lower_bound, upper_bound);
            velocity = 0;
        }
        return make_pair(position, velocity);
    }

    pair<long double,long double> Society::random_deterministic_back(long double position, long double velocity, int dimension_id)
    {
        long double upper_bound = (*(society[0])).get_upper_bound_at(dimension_id);
        long double lower_bound = (*(society[0])).get_lower_bound_at(dimension_id);

        if(position > upper_bound || position < lower_bound)
        {
            position = unifRand(lower_bound, upper_bound);
            const long double deterministic_back_constant = 0.5;
            velocity = (-1) * deterministic_back_constant * velocity;
        }
        return make_pair(position, velocity);
    }

    pair<long double,long double> Society::random_random_back(long double position, long double velocity, int dimension_id)
    {
        long double upper_bound = (*(society[0])).get_upper_bound_at(dimension_id);
        long double lower_bound = (*(society[0])).get_lower_bound_at(dimension_id);

        if(position > upper_bound || position < lower_bound)
        {
            position = unifRand(lower_bound, upper_bound);
            velocity = (-1) * unifRand() * velocity;
        }
        return make_pair(position, velocity);
    }

    long double Society::random_unmodified(long double position, int dimension_id)
    {
        long double upper_bound = (*(society[0])).get_upper_bound_at(dimension_id);
        long double lower_bound = (*(society[0])).get_lower_bound_at(dimension_id);

        if(position > upper_bound || position < lower_bound)
        {
            position = unifRand(lower_bound, upper_bound);
        }
        return position;
    }

    //  C. Preventing infeasibility
    long double Society::hyperbolic(long double old_position, long double velocity, int dimension_id)
    {
        long double upper_bound = (*(society[0])).get_upper_bound_at(dimension_id);
        long double lower_bound = (*(society[0])).get_lower_bound_at(dimension_id);

        if(velocity > 0)
        {
            velocity = velocity / (1 + velocity/(upper_bound - old_position));
        }
        else
        {
            velocity = velocity / (1 - velocity/(old_position - lower_bound));
        }
        return velocity;
    }

    long double Society::random_forth(long double old_position, long double velocity, int dimension_id)
    {
        long double upper_bound = (*(society[0])).get_upper_bound_at(dimension_id);
        long double lower_bound = (*(society[0])).get_lower_bound_at(dimension_id);

        if(old_position + velocity > upper_bound)
        {
            velocity = unifRand(0, upper_bound - old_position);
        }
        else if(old_position + velocity < lower_bound)
        {
            velocity = (-1) * unifRand(0, old_position - lower_bound);
        }
        return velocity;
    }

    //  D. Other approaches
    pair<long double,long double> Society::bounded_mirror(long double position, long double velocity, int dimension_id)
    {
        long double upper_bound = (*(society[0])).get_upper_bound_at(dimension_id);
        long double lower_bound = (*(society[0])).get_lower_bound_at(dimension_id);

        if(position > upper_bound || position < lower_bound)
        {
            // Period of image of bounded mirror is 2*(upper_bound-lower_bound)
            // The most trivial period is from lower_bound to 2*upper_bound-lower_bound
            long double numerator = position-lower_bound;
            long double denominator = 2*(upper_bound-lower_bound);
            long double remainder = fmod(numerator, denominator);
            long double projection;

            if(remainder>=0)
                projection = lower_bound + remainder;
            else
                projection = (2*upper_bound-lower_bound) + remainder;

            if(projection>=lower_bound&&projection<upper_bound)   // Projection in the first half period
            {
                position = projection;
            }
            else    // Projection in the second half period
            {
                position = 2*upper_bound - projection;
                velocity = (-1) * velocity;
            }
        }
        return make_pair(position,velocity);
    }

    long double Society::period(long double position, int dimension_id)
    {
        long double upper_bound = (*(society[0])).get_upper_bound_at(dimension_id);
        long double lower_bound = (*(society[0])).get_lower_bound_at(dimension_id);

        if(position > upper_bound || position < lower_bound)
        {
            long double numerator = position-lower_bound;
            long double denominator = upper_bound-lower_bound;
            long double remainder = fmod(numerator, denominator);
            long double projection;

            if(remainder>=0)
                projection = lower_bound + remainder;
            else
                projection = upper_bound + remainder;

            position = projection;
        }
        return position;
    }

    //  Create random numbers
    long double Society::unifRand() const
    {
        return rand() / (long double)(RAND_MAX);
    }

    long double Society::unifRand(long double a, long double b) const
    {
        return (b-a)*unifRand() + a;
    }

    int Society::input_converter(string your_input) const
    {
        bool illegal = false;
        const int size = your_input.size();
        for(int i=0; i<size; i++)
        {
            if(your_input[i]-'0' < 0 || your_input[i]-'0' > 9)
            {
                illegal = true;
                break;
            }
        }

        if(illegal)
            return -1;
        else
        {
            int total_number = 0;
            for(int i=0; i<size; i++)
                total_number += (your_input[i]-'0')*pow(10,size-i-1);
            return total_number;
        }
    }
