

//  Experiemental program
//  For Bound Handling Techniques in Particles Swarm Optimization

//  Director: Dr. Jun Zhang
//  Author: Qili Huang

//------------------------------------------------------------------
// Main function

#include "Particle.h"
#include  "Society.h"
#include "cec2005_benchmarks.h"
#include "rand.h"
#include "data_analysis.h"

#include <iostream>
#include <ctime>    // prototype for function time
#include <cstdlib>  // prototype for function srand
#include <string>   // prototype for string type


using namespace std;


int main()
{

    srand(time(0));

    Society society;

    cout << "\n---------------------------------\n";
    cout << "Welcome! Here is testing bound handling techniques in PSO.\n\n";
    string first_button;
    cout << "Please choose benchmark suites!\n"
        <<"(1) Seven Trandtional Benchmarks?\n --> Input 1.\n\n"
        <<"(2) CEC2005 Benchmarks, F1-F25, but only for 2,10,30,50 Dim.?\n --> Input 2.\n\n"
        <<"(3) Flat landscape?\n --> Input 3.\n\n"
        <<" Exit? --> Input others.\n\n"
        <<"Notice: Need uncomment or comment a few sentences for CEC2005 benchmarks\n\n"
        <<"Your input is ";
    cin >> first_button;
    cout << "\n\n";
    if(first_button == "1")
    {
        string second_button;
        cout << "(1)Test one bound handling technique with one function? Please input 1.\n"
            << "(2)Test one bound handling technique with all function? Please input 2.\n"
            << "(3)Test all bound handling technique with one function? Please input 3.\n"
            << "(4)Test all bound handling technique with all function? Please input 4.\n"
            << " Exit? Please input others.\n\n";
        cout << "Your input is ";
        cin >> second_button;
        cout << "\n\n";
        if(second_button == "1")
            society.one_one_start();
        else if(second_button == "2")
            society.one_all_start();
        else if(second_button == "3")
            society.all_one_start();
        else if(second_button == "4")
            society.all_all_start();
        else
            exit(0);
    }
    else if(first_button == "2")
    {
        string second_button;
        cout << "(1)Test one bound handling technique with one CEC2005 function? Please input 1.\n"
            << "(2)Test all bound handling technique with one CEC2005 function? Please input 2.\n"
            << " Exit? Please input others.\n\n";
        cout << "Your input is ";
        cin >> second_button;
        cout << "\n\n";
        if(second_button == "1")
            society.cec2005_one_one_start();
        else if(second_button == "2")
            society.cec2005_all_one_start();
        else
            exit(0);
        cout << "\n\nTest Over! \n";
    }
    else if(first_button == "3")
    {
        cout << "Here is flat landscape test for all bound handling techniques.\n";
        society.flat_landscape_start();
        cout << "\n\nTest Over! \n";
    }
    else
        exit(0);

    //statistic_analysis();
    //statistic_final_generation();
    //statistic_lg();
    //files_rearrangement();
    //print_technique_name();
    //print_generation_id();
    //random_choose_flat_landscape_distribution();
    //mean_stder_flat_landscape_distribution();
    //printInterval();
    return 0;
}

