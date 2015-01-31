
//  Source file for analysis codes

#include "data_analysis.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <vector>

using namespace std;


class Statistic
{
public:
    void SetValues(const long double *p, int count)
    {
        if(!(value.empty()))
            value.clear();
        amount = count;
        for(int i = 0; i < amount; i++)
            value.push_back(p[i]);
    }

    long double CalculateMean()
    {
        long double sum = 0;
        for(int i = 0; i < amount; i++)
            sum += value[i];
        return (sum / amount);
    }

    long double CalculateMedian()
    {
        vector<long double> temp (value);
        sort(temp.begin(), temp.end());
        if(amount%2 == 1)
            return temp[amount/2];
        else
            return (temp[amount/2 - 1] + temp[amount/2])/2;
    }

    long double CalculateMinElement()
    {
        return *(min_element(value.begin(), value.end()));
    }

    long double CalculateMaxElement()
    {
        return *(max_element(value.begin(), value.end()));
    }

    long double CalculateVariane()
    {
        mean = CalculateMean();

        long double temp = 0;
        for(int i = 0; i < amount; i++)
        {
             temp += (value[i] - mean) * (value[i] - mean) ;
        }
        return temp / amount;
    }

    long double CalculateSampleVariane()
    {
        mean = CalculateMean();
        long double temp = 0;
        for(int i = 0; i < amount; i++)
        {
             temp += (value[i] - mean) * (value[i] - mean) ;
        }
        return temp / (amount - 1);
    }

    long double GetStandardDeviation()
    {
        return sqrt(CalculateVariane());
    }

    long double GetSampleStandardDeviation()
    {
        return sqrt(CalculateSampleVariane());
    }

private:
    int amount;
    vector<long double> value;
    long double mean;

};


string get_benchmark_function_name(int number)
{
    if(number==1)
        return "Sphere";
    else if(number==2)
        return "Rosenbrock";
    else if(number==3)
        return "Rastringin";
    else if(number==4)
        return "Griewank";
    else if(number==5)
        return "Ackley";
    else if(number==6)
        return "Michalewicz";
    else if(number==7)
        return "Schwefel";
    else if(number>=101&&number<=125)
    {
        int intValue = number - 100;
        stringstream ss;
        ss << intValue;
        string str = ss.str();
        return "F" + str;
    }
    else
    {
        cerr << "\nError input of number(function type)\n";
        exit(0);
    }
}

 string get_bound_handling_name(int number)
{
        if(number==1)
            return "Bounded_mirror";
        else if(number==2)     //  Period
            return "Period";
        else if(number==3)
            return "Reflect_zero";
        else if(number==4)
            return "Hyperboloic";
        else if(number==5)
            return "Reflect_unmodifed";
        else if(number==6)
            return "Random_forth";
        else if(number==7)
            return "Reflect_deterministic_back";
        else if(number==8)
            return "Reflect_random_back";
        else if(number==9)
            return "Nearest_zero";
        else if(number==10)
            return "Nearest_deterministic_back";
        else if(number==11)    //  Reflect random back
            return "Nearest_random_back";
        else if(number==12)
            return "Nearest_unmodifed";
        else if(number==13)
            return "Infinity";
        else if(number==14)
            return "Infinity_clamp";
        else if(number==15)
            return "Random_zero";
        else if(number==16)
            return "Random_deterministic_back";
        else if(number==17)
            return "Random_random_back";
        else if(number==18)
            return "Random_unmodified";
        else
        {
            cerr << "\nError input of number(technique)\n";
            exit(1);
        }
}

void statistic_analysis(int id[], int count)
{
    for(int i=0; i<count; i++)
    {
        string function_name = get_benchmark_function_name(id[i]);
        for(int j=1; j<=18; j++)
        {
            string technique_name = get_bound_handling_name(j);

            // Read Experimental Data
            long double experiment_data[50][61];
            for(int k=0; k<50; k++)
            {
                stringstream ss;
                ss << k+1;
                string exp_times = ss.str();
                string input_file_name = "output_data/" + function_name + "/" + technique_name
                                        + "_" + function_name + "_30_6000_" + exp_times + ".txt";

                ifstream input(input_file_name.c_str(), ios::in);
                if(!input)
                {
                    cerr << "\nFail to open file for reading.\n";
                    exit(0);
                }
                //  Sift text header
                const int size = 256;
                char temp[size];
                for(int m=0; m<9; m++)
                    input.getline(temp, size);
                //  Read
                for(int l=0; l<61; l++)
                {
                    int first_item;
                    input >> first_item >> experiment_data[k][l];
                }
                input.close();
            }

            //  Statistic Operation and Write
            string output_file_name1 = "output_data/" + function_name + "/data_analysis/best, worst/"
                                        + technique_name + "_" + function_name + "_best_worst.txt";
            string output_file_name2 = "output_data/" + function_name + "/data_analysis/mean, standard error/"
                                        + technique_name + "_" + function_name + "_mean_standard_error.txt";
            string output_file_name3 = "output_data/" + function_name + "/data_analysis/median/"
                                        + technique_name + "_" + function_name + "_median.txt";

            ofstream output1(output_file_name1.c_str(), ios::out);
            ofstream output2(output_file_name2.c_str(), ios::out);
            ofstream output3(output_file_name3.c_str(), ios::out);
            if(!output1)
            {
                cerr << "\nFail to open file1 for writing.\n";
                exit(0);
            }
            if(!output2)
            {
                cerr << "\nFail to open file2 for writing.\n";
                exit(0);
            }
            if(!output3)
            {
                cerr << "\nFail to open file3 for writing.\n";
                exit(0);
            }

            //  Write text header
            output1 << "Technique name: " << technique_name << endl
                    <<"Function name: " << function_name << endl
                    << "Dimension: " << 30 << endl
                    << "Iteration: " << 6000 << endl
                    << "Experimental times: " << 50 << endl
                    << "Velocity update strategy: " << "V=X*[V+c1*e1(lbest-pos)+c2*e2(pbest-pos])" << endl
                    << "V_Clamp: " << 10000 << endl << endl
                    << "Iteration\t" << "Best\t" << "Worst" << endl;
            output2 << "Technique name: " << technique_name << endl
                    <<"Function name: " << function_name << endl
                    << "Dimension: " << 30 << endl
                    << "Iteration: " << 6000 << endl
                    << "Experimental times: " << 50 << endl
                    << "Velocity update strategy: " << "V=X*[V+c1*e1(lbest-pos)+c2*e2(pbest-pos])" << endl
                    << "V_Clamp: " << 10000 << endl << endl
                    << "Iteration\t" << "Mean\t" << "Standard Error" << endl;
            output3 << "Technique name: " << technique_name << endl
                    <<"Function name: " << function_name << endl
                    << "Dimension: " << 30 << endl
                    << "Iteration: " << 6000 << endl
                    << "Experimental times: " << 50 << endl
                    << "Velocity update strategy: " << "V=X*[V+c1*e1(lbest-pos)+c2*e2(pbest-pos])" << endl
                    << "V_Clamp: " << 10000 << endl << endl
                    << "Iteration\t" << "Median\t" << endl;

            Statistic pso_statistic;
            for(int k=0; k<61; k++)
            {
                long double sample[50];
                for(int l=0; l<50; l++)
                {
                    sample[l] = experiment_data[l][k];
                }
                //  At any a generation.
                pso_statistic.SetValues(sample, 50);
                output1 << k*100 << "\t" << scientific << pso_statistic.CalculateMinElement()
                        << "\t" << pso_statistic.CalculateMaxElement() << endl;
                output2 << k*100 << "\t" << scientific << pso_statistic.CalculateMean()
                        << "\t" << pso_statistic.GetStandardDeviation() << endl;
                output3 << k*100 << "\t" << scientific << pso_statistic.CalculateMedian() << endl;
            }
            output1.close();
            output2.close();
            output3.close();
        }
    }
   cout << "\nStatistic Step Done!\n\n";
}

void statistic_analysis()
{
    cout << "*Welcome use statistic analysis function!\n"
            <<"\nThis is for:\n\t(1) mean, standard error; \n\t(2) median; \n\t(3) best and worst.\n"
            <<"\n(This part usually will take less than a minute.)\n";

    cout << "\n\nAnalyze all functions? Input y/Y, otherwise input others. ";
    string all;

    cin >> all;
    if(all == "Y" || all == "y")
    {
        int files[] = {1,2,3,4,5,6,7,101,102,103,104,105,106,108,109,110,111,112,113,114};
        statistic_analysis(files, 20);
    }
    else
    {
        cout << "\n\nPlease input number of objects(functions) you need analyze. ";
        int number;
        cin >> number;
        cout << "\n\nPlease input id of objects. \n(1-7 for tradiational benchmarks, and 101-114 for F1-F14 and F7 is excluded.)\n\nYour input is:\n";
        int *files = new int[number];
        for(int i=0; i<number; i++)
            cin >> files[i];
        statistic_analysis(files, number);
    }
}

void statistic_final_generation(int id[], int count)
{
    string output_file_name4 = "output_data/data_analysis_for_all_function/best, worst/best_worst_final.txt";
    string output_file_name5 = "output_data/data_analysis_for_all_function/mean, standard error/mean_standard_error_final.txt";
    string output_file_name6 = "output_data/data_analysis_for_all_function/median/median_final.txt";
    string output_file_name7 = "output_data/data_analysis_for_all_function/best, worst/best_worst_final_onlydata.txt";
    string output_file_name8 = "output_data/data_analysis_for_all_function/mean, standard error/mean_standard_error_final_onlydata.txt";
    string output_file_name9 = "output_data/data_analysis_for_all_function/median/median_final_onlydata.txt";
    string output_file_name10 = "output_data/data_analysis_for_all_function/best, worst/best_worst_final_sort.txt";
    string output_file_name11 = "output_data/data_analysis_for_all_function/mean, standard error/mean_standard_error_final_sort.txt";
    string output_file_name12 = "output_data/data_analysis_for_all_function/median/median_final_sort.txt";

    ofstream output4(output_file_name4.c_str(), ios::out);
    ofstream output5(output_file_name5.c_str(), ios::out);
    ofstream output6(output_file_name6.c_str(), ios::out);
    ofstream output7(output_file_name7.c_str(), ios::out);
    ofstream output8(output_file_name8.c_str(), ios::out);
    ofstream output9(output_file_name9.c_str(), ios::out);
    ofstream output10(output_file_name10.c_str(), ios::out);
    ofstream output11(output_file_name11.c_str(), ios::out);
    ofstream output12(output_file_name12.c_str(), ios::out);
    if(!output4)
    {
        cerr << "\nFail to open file4: " << output_file_name4 << " for writing.\n";
        exit(0);
    }
    if(!output5)
    {
        cerr << "\nFail to open file5: " << output_file_name5 << " for writing.\n";
        exit(0);
    }
    if(!output6)
    {
        cerr << "\nFail to open file6: " << output_file_name6 << " for writing.\n";
        exit(0);
    }
    if(!output7)
    {
        cerr << "\nFail to open file7: " << output_file_name7 << " for writing.\n";
        exit(0);
    }
    if(!output8)
    {
        cerr << "\nFail to open file8: " << output_file_name8 << " for writing.\n";
        exit(0);
    }
    if(!output9)
    {
        cerr << "\nFail to open file9: " << output_file_name9 << " for writing.\n";
        exit(0);
    }
    if(!output10)
    {
        cerr << "\nFail to open file10: " << output_file_name10 << " for writing.\n";
        exit(0);
    }
    if(!output11)
    {
        cerr << "\nFail to open file11: " << output_file_name11 << " for writing.\n";
        exit(0);
    }
    if(!output12)
    {
        cerr << "\nFail to open file12: " << output_file_name12 << " for writing.\n";
        exit(0);
    }

    //  Write text header
    output4 <<"Readme: Here are all the best and worst data for any a Function with any a Technique.\n\n"
            << "Dimension: " << 30 << endl
            << "Iteration: " << 6000 << endl
            << "Experimental times: " << 50 << endl
            << "Velocity update strategy: " << "V=X*[V+c1*e1(lbest-pos)+c2*e2(pbest-pos])" << endl
            << "V_Clamp: " << 10000 << endl
            << "Best\t" << "Worst";
    output5 <<"Readme: Here are all the mean and standard error data for any a Function with any a Technique.\n\n"
            << "Dimension: " << 30 << endl
            << "Iteration: " << 6000 << endl
            << "Experimental times: " << 50 << endl
            << "Velocity update strategy: " << "V=X*[V+c1*e1(lbest-pos)+c2*e2(pbest-pos])" << endl
            << "V_Clamp: " << 10000 << endl
            << "Mean\t" << "Standard Error";
    output6 <<"Readme: Here are all the median data for any a Function with any a Technique.\n\n"
            << "Dimension: " << 30 << endl
            << "Iteration: " << 6000 << endl
            << "Experimental times: " << 50 << endl
            << "Velocity update strategy: " << "V=X*[V+c1*e1(lbest-pos)+c2*e2(pbest-pos])" << endl
            << "V_Clamp: " << 10000 << endl
            << "Median\t";

    for(int i=0; i<count; i++)
    {
        string function_name = get_benchmark_function_name(id[i]);

        output4 << "\n\n\n" << i+1 << "." << function_name << endl;
        output5 << "\n\n\n" << i+1 << "." << function_name << endl;
        output6 << "\n\n\n" << i+1 << "." << function_name << endl;
        output7 << "\n\n\n" << i+1 << "." << function_name << endl;
        output8 << "\n\n\n" << i+1 << "." << function_name << endl;
        output9 << "\n\n\n" << i+1 << "." << function_name << endl;
        output10 << "\n\n\n" << i+1 << "." << function_name << endl;
        output11 << "\n\n\n" << i+1 << "." << function_name << endl;
        output12 << "\n\n\n" << i+1 << "." << function_name << endl;

        long double final_generation_data1[18][2];
        long double final_generation_data2[18][2];
        long double final_generation_data3[18];
        for(int j=1; j<=18; j++)
        {
            string technique_name = get_bound_handling_name(j);

            string input_file_name1 = "output_data/" + function_name + "/data_analysis/best, worst/"
                                    + technique_name + "_" + function_name + "_best_worst.txt";
            string input_file_name2 = "output_data/" + function_name + "/data_analysis/mean, standard error/"
                                    + technique_name + "_" + function_name + "_mean_standard_error.txt";
            string input_file_name3 = "output_data/" + function_name + "/data_analysis/median/"
                                    + technique_name + "_" + function_name + "_median.txt";

            ifstream input1(input_file_name1.c_str(), ios::in);
            ifstream input2(input_file_name2.c_str(), ios::in);
            ifstream input3(input_file_name3.c_str(), ios::in);
            if(!input1)
            {
                cerr << "\nFail to open file: " << input_file_name1 << " for reading.\n";
                exit(0);
            }
            if(!input2)
            {
                cerr << "\nFail to open file: " << input_file_name2 << " for reading.\n";
                exit(0);
            }
            if(!input3)
            {
                cerr << "\nFail to open file: " << input_file_name3 << " for reading.\n";
                exit(0);
            }
            //  Sift text header
            const int size = 256;
            char temp[size];
            for(int m=0; m<9; m++)
            {
                input1.getline(temp, size);
                input2.getline(temp, size);
                input3.getline(temp, size);
            }
            //  Read. First 60 lines are overwritten.
            for(int l=0; l<61; l++)
            {
                int first_item;
                input1 >> first_item >> final_generation_data1[j-1][0] >> final_generation_data1[j-1][1];
                input2 >> first_item >> final_generation_data2[j-1][0] >> final_generation_data2[j-1][1];
                input3 >> first_item >> final_generation_data3[j-1];
            }
            input1.close();
            input2.close();
            input3.close();
        }

        string output_file_name1 = "output_data/" + function_name + "/data_analysis/best, worst/only final generation/"
                                    + function_name + "_best_worst_final.txt";
        string output_file_name2 = "output_data/" + function_name + "/data_analysis/mean, standard error/only final generation/"
                                    + function_name + "_mean_standard_error_final.txt";
        string output_file_name3 = "output_data/" + function_name + "/data_analysis/median/only final generation/"
                                    + function_name + "_median_final.txt";
        ofstream output1(output_file_name1.c_str(), ios::out);
        ofstream output2(output_file_name2.c_str(), ios::out);
        ofstream output3(output_file_name3.c_str(), ios::out);
        if(!output1)
        {
            cerr << "\nFail to open file1: " << output_file_name1 << " for writing.\n";
            exit(0);
        }
        if(!output2)
        {
            cerr << "\nFail to open file2: " << output_file_name2 << " for writing.\n";
            exit(0);
        }
        if(!output3)
        {
            cerr << "\nFail to open file3: " << output_file_name3 << " for writing.\n";
            exit(0);
        }

        //  Write text header
        output1 <<"Function name: " << function_name << endl
                << "Dimension: " << 30 << endl
                << "Iteration: " << 6000 << endl
                << "Experimental times: " << 50 << endl
                << "Velocity update strategy: " << "V=X*[V+c1*e1(lbest-pos)+c2*e2(pbest-pos])" << endl
                << "V_Clamp: " << 10000 << endl << endl
                << "Best\t" << "Worst" << endl;
        output2 <<"Function name: " << function_name << endl
                << "Dimension: " << 30 << endl
                << "Iteration: " << 6000 << endl
                << "Experimental times: " << 50 << endl
                << "Velocity update strategy: " << "V=X*[V+c1*e1(lbest-pos)+c2*e2(pbest-pos])" << endl
                << "V_Clamp: " << 10000 << endl << endl
                << "Mean\t" << "Standard Error" << endl;
        output3 <<"Function name: " << function_name << endl
                << "Dimension: " << 30 << endl
                << "Iteration: " << 6000 << endl
                << "Experimental times: " << 50 << endl
                << "Velocity update strategy: " << "V=X*[V+c1*e1(lbest-pos)+c2*e2(pbest-pos])" << endl
                << "V_Clamp: " << 10000 << endl << endl
                << "Median\t" << endl;

        for(int j=1; j<=18; j++)
        {
            output1 << "(" << j << ")" <<get_bound_handling_name(j) << "\n" << scientific
                    << final_generation_data1[j-1][0] << "\t" << final_generation_data1[j-1][1] << endl;
            output2 << "(" << j << ")" <<get_bound_handling_name(j) << "\n" << scientific
                    << final_generation_data2[j-1][0] << "\t" << final_generation_data2[j-1][1] << endl;
            output3 << "(" << j << ")" <<get_bound_handling_name(j) << "\n" << scientific
                    << final_generation_data3[j-1] << endl;

            output4 << "(" << j << ")" <<get_bound_handling_name(j) << "\n" << scientific
                    << final_generation_data1[j-1][0] << "\t" << final_generation_data1[j-1][1] << endl;
            output5 << "(" << j << ")" <<get_bound_handling_name(j) << "\n" << scientific
                    << final_generation_data2[j-1][0] << "\t" << final_generation_data2[j-1][1] << endl;
            output6 << "(" << j << ")" <<get_bound_handling_name(j) << "\n" << scientific
                    << final_generation_data3[j-1] << endl;

            output7 << scientific << final_generation_data1[j-1][0] << "\t" << final_generation_data1[j-1][1] << endl;
            output8 << scientific << final_generation_data2[j-1][0] << "\t" << final_generation_data2[j-1][1] << endl;
            output9 << scientific << final_generation_data3[j-1] << endl;

        }
        output1.close();
        output2.close();
        output3.close();

        long double temp1[18];
        long double temp2[18];

        for(int j=0; j<18; j++)
        {
            temp1[j] = final_generation_data1[j][0];
            temp2[j] = final_generation_data1[j][1];
        }
        output10 << "\tmin\tmax\nbest\t" << scientific
                << *(min_element(temp1, temp1+18)) << "\t" << *(max_element(temp1, temp1+18)) << endl;
        output10 << "\tmin\tmax\nworst\t" << scientific
                << *(min_element(temp2, temp2+18)) << "\t" << *(max_element(temp2, temp2+18)) << endl;
        for(int j=0; j<18; j++)
        {
            temp1[j] = final_generation_data2[j][0];
            temp2[j] = final_generation_data2[j][1];
        }
        output11 << "\tmin\tmax\nmean\t" << scientific
                << *(min_element(temp1, temp1+18)) << "\t" << *(max_element(temp1, temp1+18)) << endl;
        output11 << "\tmin\tmax\nstandard error\t" << scientific
                << *(min_element(temp2, temp2+18)) << "\t" << *(max_element(temp2, temp2+18)) << endl;
        for(int j=0; j<18; j++)
            temp1[j] = final_generation_data3[j];
        output12 << "\tmin\tmax\nmedian\t" << scientific
                << *(min_element(temp1, temp1+18)) << "\t" << *(max_element(temp1, temp1+18)) << endl;
    }
    output4.close();
    output5.close();
    output6.close();
    output7.close();
    output8.close();
    output9.close();
    output10.close();
    output11.close();
    output12.close();
    cout << "\nStatistic Step Done!\n\n";
}

void statistic_final_generation()
{
    cout << "\n*Welcome use statistic final function!\n"
            <<"\nThis is for final generation's values:\n\t(1) mean, standard error; \n\t(2) median; \n\t(3) best and worst.\n"
            <<"\n(This part usually will take less than a minute.)\n";

    cout << "\n\nAnalyze all functions? Input y/Y, otherwise input others. ";
    string all;

    cin >> all;
    if(all == "Y" || all == "y")
    {
        int files[] = {1,2,3,4,5,6,7,101,102,103,104,105,106,108,109,110,111,112,113,114};
        statistic_final_generation(files, 20);
    }
    else
    {
        cout << "\n\nPlease input number of objects(functions) you need analyze. ";
        int number;
        cin >> number;
        cout << "\n\nPlease input id of objects. \n(1-7 for tradiational benchmarks, and 101-114 for F1-F14 and F7 is excluded.)\n\nYour input is:\n";
        int *files = new int[number];
        for(int i=0; i<number; i++)
            cin >> files[i];
        statistic_final_generation(files, number);
    }
}

void statistic_lg(int id[], int count)
{
    for(int i=0; i<count; i++)
    {
        string function_name = get_benchmark_function_name(id[i]);

        string output_all_file_name = "output_data/" + function_name + "/data_analysis/median/log10/all_lg.txt";
        string output_all_file_name2 = "output_data/" + function_name + "/data_analysis/median/log10/all_onlydata_lg.txt";

        ofstream output_all(output_all_file_name.c_str(), ios::out);
        ofstream output_all2(output_all_file_name2.c_str(), ios::out);
        if(!output_all)
        {
            cerr << "\nFail to open file: " << output_all_file_name << " for writing.";
            exit(0);
        }
        if(!output_all2)
        {
            cerr << "\nFail to open file: " << output_all_file_name2 << " for writing.";
        }
        output_all << "\n\nFunction name: " << function_name << "\n\n\n\n";

        long double temp_all[18][61];
        //  These functions' median of 50 times of experiment will be no less than 0.
        if(id[i]==1 || id[i]==2 || id[i]==3 || id[i]==4 || id[i] ==5 || id[i]==106 || id[i]==111)
        {
            for(int j=1; j<=18; j++)
            {
                //  Open output text file
                string technique_name = get_bound_handling_name(j);

                output_all << "(" << j << ")." << technique_name << endl << endl;

                string output_file_name = "output_data/" + function_name + "/data_analysis/median/log10/"
                                            + technique_name + "_" + function_name + "_median_lg.txt";
                ofstream output(output_file_name.c_str(), ios::out);
                if(!output)
                {
                    cerr << "\nFail to open file for writing.\n";
                    exit(0);
                }
                //  Output text header
                output << "Function: " << function_name << "'s median theoretically is no less than 0.\n"
                        << "(1).Calculate lg(smallest positive median), if the median is equal to zero;\n"
                        << "(2).Calculate lg(smallest positive number), if it is in Ackley function and also"
                        << " as a result of precision of E in the computer, the median can be slightly less than zero.\n"
                        << "\n\n\nlog10(Median)\t\n\n";

                string input_file_name = "output_data/" + function_name + "/data_analysis/median/"
                                        + technique_name + "_" + function_name + "_median.txt";
                ifstream input(input_file_name.c_str(), ios::in);
                if(!input)
                {
                    cerr << "\nFail to open file: " << input_file_name << " for reading.\n";
                    exit(0);
                }

                //  Sieve useless lines out
                const int size = 256;
                char temp[size];
                for(int m=0; m<9; m++)
                    input.getline(temp, size);

                long double temp1[61];
                long double smallest_positive_number;
                for(int l=0; l<61; l++)
                {
                    int first_item;
                    input >> first_item >> temp1[l];
                    if(temp1[l] <= 0)
                    {
                        if(temp1[l-1] > 0)
                            smallest_positive_number = temp1[l-1];
                        output << scientific << log10(smallest_positive_number) << endl;
                        output_all << scientific << log10(smallest_positive_number) << endl;
                        temp_all[j-1][l] = log10(smallest_positive_number);
                        continue;
                    }
                    output << scientific << log10(temp1[l]) << endl;
                    output_all << scientific << log10(temp1[l]) << endl;
                    temp_all[j-1][l] = log10(temp1[l]);
                }
                output_all << "\n\n\n\n";
                input.close();
                output.close();
            }
        }

        //  These functions' median of 50 times of experiment can be either the positive number or negative number
        else
        {
            for(int j=1; j<=18; j++)
            {
                //  Open output text file
                string technique_name = get_bound_handling_name(j);

                output_all << "(" << j << ")." << technique_name << endl << endl;

                string output_file_name = "output_data/" + function_name + "/data_analysis/median/log10/"
                                            + technique_name + "_" + function_name + "_median_lg.txt";
                ofstream output(output_file_name.c_str(), ios::out);
                if(!output)
                {
                    cerr << "\nFail to open file: "<< output_file_name << " for writing.\n";
                    exit(0);
                }
                //  Output text header
                output << "Function: " << function_name << "'s median can be either the positive number or negative number.\n"
                        << "1. Calculate log10(median) directly, if median is bigger than 1;\n"
                        << "2. Set lg(median) as 0, if median is less than 1 and bigger than -1;\n"
                        << "3. Calculate the -log10(abs(median)), if median is less than -1.\n"
                        << "\n\nlog10(Median)\n\n";

                string input_file_name = "output_data/" + function_name + "/data_analysis/median/"
                                        + technique_name + "_" + function_name + "_median.txt";
                ifstream input(input_file_name.c_str(), ios::in);
                if(!input)
                {
                    cerr << "\nFail to open file for reading.\n";
                    exit(0);
                }

                //  Sieve useless lines out
                const int size = 256;
                char temp[size];
                for(int m=0; m<9; m++)
                    input.getline(temp, size);

                long double temp1;
                for(int l=0; l<61; l++)
                {
                    int first_item;
                    input >> first_item >> temp1;
                    if(temp1>1)
                    {
                        output << scientific << log10(temp1) << endl;
                        output_all << scientific << log10(temp1) << endl;
                        temp_all[j-1][l] = log10(temp1);
                    }
                    else if(temp1<=1 && temp1>=-1)
                    {
                        output << scientific << 0 << endl;
                        output_all << scientific << 0 << endl;
                        temp_all[j-1][l] = 0;
                    }
                    else
                    {
                        output << scientific << -log10(abs(temp1)) << endl;
                        output_all << scientific << -log10(abs(temp1)) << endl;
                        temp_all[j-1][l] = -log10(abs(temp1));
                    }
                }
                output_all << "\n\n\n\n";
                input.close();
                output.close();
            }
        }
        for(int l=0; l<61; l++)
        {
            output_all2 << l*100 << "\t";
            for(int j=0; j<17; j++)
                output_all2 << scientific << temp_all[j][l] << "\t";
            output_all2 << scientific << temp_all[17][l] << endl;
        }
        output_all.close();
        output_all2.close();
    }
    cout << "\nStatistic Step Done!\n";
}


void statistic_lg()
{

    cout << "\n*Welcome use statistic lg function!\n"
            <<"\nThis is for all generation's values:\n\t(1) mean, standard error;\n"
            <<"\n(This part usually will take less than a minute.)\n";

    cout << "\n\nAnalyze all functions? Input y/Y, otherwise input others. ";
    string all;

    cin >> all;
    if(all == "Y" || all == "y")
    {
        int files[] = {1,2,3,4,5,6,7,101,102,103,104,105,106,108,109,110,111,112,113,114};
        statistic_lg(files, 20);
    }
    else
    {
        cout << "\n\nPlease input number of objects(functions) you need analyze. ";
        int number;
        cin >> number;
        cout << "\n\nPlease input id of objects. \n(1-7 for tradiational benchmarks, and 101-114 for F1-F14 and F7 is excluded.)\n\nYour input is:\n";
        int *files = new int[number];
        for(int i=0; i<number; i++)
            cin >> files[i];
        statistic_lg(files, number);
    }
}

void files_rearrangement()
{
    string file_name;
    string technique_name;
    string function_name;
    string exp_times;

    //file_name = "output_data/Sphere_other15exp_";
    file_name = "output_data/F12/2/";

    function_name = get_benchmark_function_name(112);
    for(int i=7; i<=18; i++)
    {
        technique_name = get_bound_handling_name(i);
        for(int k=1; k<=25; k++)
        {
            stringstream ss1;
            ss1 << k;
            exp_times = ss1.str();

            string path_in = file_name + "/" + technique_name + "_" + function_name + "_30_6000_" + exp_times + ".txt";
            ifstream input(path_in.c_str(), ios::in);
            if(!input)
            {
                cerr << "\nFail to open the file for reading.\n";
                exit(0);
            }

            const int size = 256;
            char temp[size];
            for(int l=0; l<9; l++)
                input.getline(temp, size);
            int generation[100];
            long double evals[100];
            for(int m=0; m<61; m++)
                input >> generation[m] >> evals[m];

            stringstream ss2;
            ss2 << k+25;
            exp_times = ss2.str();

            string path_out = "output_data/F12/" + technique_name + "_" + function_name + "_30_6000_" + exp_times + ".txt";
            ofstream output(path_out.c_str(), ios::out);
            if(!output)
            {
                cerr << "\nFail to open file for writing.\n";\
                exit(0);
            }
            output << "Technique name: " << technique_name << endl
                    <<"Function name: " << function_name << endl
                    << "Dimension: " << 30 << endl
                    << "Iteration: " << 6000 << endl
                    << "Experimental times: " << exp_times << endl
                    << "Velocity update strategy: " << "V=X*[V+c1*e1(lbest-pos)+c2*e2(pbest-pos])" << endl
                    << "V_Clamp: " << 10000 << endl << endl
                    << "Iteration\t" << "Evaluaton" << endl;
            for(int m=0; m<61; m++)
                output << generation[m] << "\t" << scientific << evals[m] << endl;
        }
        cout << "\ndone!\n\n";
    }
}

void print_technique_name()
{
    ofstream output("output_data/data_analysis_for_all_function/all_technique_name.txt", ios::out);
    if(!output)
    {
        cerr << "\nFail to open file for writing.\n";
        exit(0);
    }
    for(int i=1; i<=18; i++)
        output << get_bound_handling_name(i) << endl;
    output.close();
    cout << "\nPrint done.\n";
}

void print_generation_id()
{
    ofstream output("output_data/data_analysis_for_all_function/generation_id.txt", ios::out);
    if(!output)
    {
        cerr << "\nFail to open file for writing.\n";
        exit(0);
    }
    for(int i=0; i<61; i++)
        output << i*100 << endl;
    output.close();
    cout << "\nPrint done.\n";
}

void random_choose_flat_landscape_distribution()
{
    cout << "\nWelcome use this function to random choose flat landscape distribution.\n";

    string file_name_all = "output_data/flat_landscape/random_chosen_distribution/all_onlydata.txt";
    ofstream all_output(file_name_all.c_str(), ios::out);

    long double all_temp[18][20];
    for(int i=1; i<=18; i++)
    {
        string technique_name = get_bound_handling_name(i);

        int at = rand() % 10 + 1;
        stringstream ss;
        ss << at;
        string file_id = ss.str();
        string file_name = "output_data/flat_landscape/" + technique_name + "_" + file_id + ".txt";
        ifstream input(file_name.c_str(), ios::in);
        if(!input)
        {
            cerr << "\nFail to open file: " << file_name << " for reading.\n";
            exit(0);
        }

        at = rand() % 30 + 1;
        ss.str(string());
        ss << at;
        string dimension_id = ss.str();
        string output_file_name = "output_data/flat_landscape/random_chosen_distribution/"
                                    + technique_name + "_" + file_id + "_dimension_" + dimension_id + ".txt";
        ofstream output(output_file_name.c_str(), ios::out);
        if(!output)
        {
            cerr << "\nFail to open the file: " << output_file_name << " for writing.\n";
            exit(0);
        }

        char temp[256];
        for(int j=0; j<at-1; j++)
            for(int k=0; k<22; k++)
                input.getline(temp, 256);
        input.getline(temp, 256);
        input .getline(temp, 256);

        long double p_temp;
        for(int j=0; j<20; j++)
        {
            input >> p_temp;
            all_temp[i-1][j] = p_temp;
            output << p_temp << endl;
        }

        input.close();
        output.close();
    }

    for(int j=0; j<20; j++)
    {
        for(int i=0; i<18; i++)
            all_output << scientific << all_temp[i][j] << "\t";
        all_output << endl;
    }
    all_output.close();

    cout << "\nStatistic Step Done!\n";
}

void printInterval()
{
    string file_name = "output_data/flat_landscape/random_chosen_distribution/interval.txt";
    ofstream output(file_name.c_str(), ios::out);
    for(int i=0; i<20; i++)
        output << -100+i*10 << endl;
    output << endl;
    for(int i=0; i<20; i++)
        output << i+1 << endl;
    output.close();
    cout << "\nStatistic Step Done\n";
}

void mean_stder_flat_landscape_distribution()
{
    cout << "\nWelcome use this function to analysis the flat landscape distribution of particles.\n";

    long double all_temp[18][20][2];
    for(int i=0; i<18; i++)
    {
        string technique_name = get_bound_handling_name(i+1);

        int at = rand() % 30 + 1;
        stringstream ss;
        ss << at;
        string dimension_id = ss.str();

        long double sample_temp[10][20];
        for(int t=0; t<10; t++)
        {
            ss.str(string());
            ss << t+1;
            string file_id = ss.str();
            string file_name = "output_data/flat_landscape/" + technique_name + "_" + file_id + ".txt";
            ifstream input(file_name.c_str(), ios::in);
            if(!input)
            {
                cerr << "\nFail to open file: " << file_name << " for reading.\n";
                exit(0);
            }

            char temp[256];
            for(int j=0; j<at-1; j++)
                for(int k=0; k<22; k++)
                    input.getline(temp, 256);
            input.getline(temp, 256);
            input .getline(temp, 256);

            long double p_temp;
            for(int j=0; j<20; j++)
            {
                input >> p_temp;
                sample_temp[t][j] = p_temp;
            }

            input.close();
        }

        string output_file_name = "output_data/flat_landscape/mean_standard_error_distribution/"
                                    + technique_name + "_dimension_" + dimension_id + ".txt";
        ofstream output(output_file_name.c_str(), ios::out);
        if(!output)
        {
            cerr << "\nFail to open the file: " << output_file_name << " for writing.\n";
            exit(0);
        }

        long double temp[10];
        Statistic st;
        for(int j=0; j<20; j++)
        {
            for(int t=0; t<10; t++)
                temp[t] = sample_temp[t][j];
            st.SetValues(temp, 10);
            all_temp[i][j][0] = st.CalculateMean();
            all_temp[i][j][1] = st.GetStandardDeviation();
            output << scientific << all_temp[i][j][0] << "\t" << all_temp[i][j][1] << endl;
        }
        output.close();
    }

    string file_name_all = "output_data/flat_landscape/mean_standard_error_distribution/all_onlydata.txt";
    string file_name_mean = "output_data/flat_landscape/mean_standard_error_distribution/all_mean_onlydata.txt";
    ofstream all_output(file_name_all.c_str(), ios::out);
    ofstream all_mean_output(file_name_mean.c_str(), ios::out);
    if(!all_output)
    {
        cerr << "\nFail to open file: " << file_name_all << " for writing.\n";
        exit(0);
    }
    if(!all_mean_output)
    {
        cerr << "\nFail to open file: " << file_name_mean << " for writing.\n";
        exit(0);
    }

    for(int j=0; j<20; j++)
    {
        for(int i=0; i<18; i++)
        {
            all_output << scientific << all_temp[i][j][0] << "\t" << all_temp[i][j][1] << "\t";
            all_mean_output << scientific << all_temp[i][j][0] << "\t";
        }
        all_output << endl;
        all_mean_output << endl;
    }
    all_output.close();
    all_mean_output.close();

    cout << "\nStatistic Step Done!\n";
}
