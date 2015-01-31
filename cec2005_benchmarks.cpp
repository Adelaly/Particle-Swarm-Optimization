
// Definitions

# include "cec2005_benchmarks.h"
# include "rand.h"

# include <iostream>
# include <cstdlib>
# include <cmath>
# include <fstream> // fstream
# include <ctime> // prototype for time
# include <cstdlib> // prototype for srand, rand()


using namespace std;

/* Global variables that you are required to initialize */
int nreal;                    /* number of real variables */
int nfunc;                    /* number of basic functions */
long double bound;            /* required for plotting the function profiles for nreal=2 */
int density;                /* density of grid points for plotting for nreal=2 */

/* Global variables being used in evaluation of various functions */
/* These are initalized in file def2 */
long double C;
long double global_bias;
long double *trans_x;
long double *basic_f;
long double *temp_x1;
long double *temp_x2;
long double *temp_x3;
long double *temp_x4;
long double *weight;
long double *sigma;
long double *lambda;
long double *bias;
long double *norm_x;
long double *norm_f;
long double **o;
long double **g;
long double ***l;



//==================================================================
//==================================================================

/* 1 Auxillary function definitions */

/* Function to return the maximum of two variables */
long double maximum (long double a, long double b)
{
    if (a>b)
    {
        return(a);
    }
    return (b);
}

/* Function to return the minimum of two variables */
long double minimum (long double a, long double b)
{
    if (a<b)
    {
        return (a);
    }
    return (b);
}

/* Function to return the modulus of a vector */
long double modulus (long double *x, int n)
{
    int i;
    long double res;
    res = 0.0;
    for (i=0; i<n; i++)
    {
        res += x[i]*x[i];
    }
    return (sqrt(res));
}

/* Function to return the dot product of two vecors */
long double dot (long double *a, long double *b, int n)
{
    int i;
    long double res;
    res = 0.0;
    for (i=0; i<n; i++)
    {
        res += a[i]*b[i];
    }
    return (res);
}

/* Function to return the mean of n variables */
long double mean (long double *x, int n)
{
    int i;
    long double res;
    res = 0.0;
    for (i=0; i<n; i++)
    {
        res += x[i];
    }
    return (res/(long double)n);
}


//=====================================================================
//=====================================================================

/* 2 Basic funcion definitions */

/* Code to evaluate ackley's function */
long double calc_ackley (long double *x)
{
    int i;
    long double sum1, sum2, res;
    sum1 = 0.0;
    sum2 = 0.0;
    for (i=0; i<nreal; i++)
    {
        sum1 += x[i]*x[i];
        sum2 += cos(2.0*PI*x[i]);
    }
    sum1 = -0.2*sqrt(sum1/nreal);
    sum2 /= nreal;
    res = 20.0 + E - 20.0*exp(sum1) - exp(sum2);
    return (res);
}

/* Code to evaluate rastrigin's function */
long double calc_rastrigin (long double *x)
{
    int i;
    long double res;
    res = 0.0;
    for (i=0; i<nreal; i++)
    {
        res += (x[i]*x[i] - 10.0*cos(2.0*PI*x[i]) + 10.0);
    }
    return (res);
}

/* Code to evaluate weierstrass's function */
long double calc_weierstrass (long double *x)
{
    int i, j;
    long double res;
    long double sum;
    long double a, b;
    int k_max;
    a = 0.5;
    b = 3.0;
    k_max = 20;
    res = 0.0;
    for (i=0; i<nreal; i++)
    {
        sum = 0.0;
        for (j=0; j<=k_max; j++)
        {
            sum += pow(a,j)*cos(2.0*PI*pow(b,j)*(x[i]+0.5));
        }
        res += sum;
    }
    return (res);
}

/* Code to evaluate griewank's function */
long double calc_griewank (long double *x)
{
    int i;
    long double s, p;
    long double res;
    s = 0.0;
    p = 1.0;
    for (i=0; i<nreal; i++)
    {
        s += x[i]*x[i];
        p *= cos(x[i]/sqrt(1.0+i));
    }
    res = 1.0 + s/4000.0 - p;
    return (res);
}

/* code to evaluate sphere function */
long double calc_sphere (long double *x)
{
    int i;
    long double res;
    res = 0.0;
    for (i=0; i<nreal; i++)
    {
        res += x[i]*x[i];
    }
    return (res);
}

/* Code to evaluate schwefel's function */
long double calc_schwefel (long double *x)
{
    int i, j;
    long double sum1, sum2;
    sum1 = 0.0;
    for (i=0; i<nreal; i++)
    {
        sum2 = 0.0;
        for (j=0; j<=i; j++)
        {
            sum2 += x[j];
        }
        sum1 += sum2*sum2;
    }
    return (sum1);
}

/* Code to evaluate rosenbrock's function */
long double calc_rosenbrock (long double *x)
{
    int i;
    long double res;
    res = 0.0;
    for (i=0; i<nreal-1; i++)
    {
        res += 100.0*pow((x[i]*x[i]-x[i+1]),2.0) + 1.0*pow((x[i]-1.0),2.0);
    }
    return (res);
}

/* Code to evaluate schaffer's function and rounding-off variables */
long double nc_schaffer (long double x, long double y)
{
    int i;
    int a;
    long double b;
    long double res;
    long double temp1, temp2;
    long double t1[2], t2[2];
    t1[0] = x;
    t1[1] = y;
    for (i=0; i<2; i++)
    {
        if (fabs(t1[i]) >= 0.5)
        {
            res = 2.0*t1[i];
            a = res;
            b = fabs(res-a);
            if (b<0.5)
            {
                t2[i] = a/2.0;
            }
            else
            {
                if (res<=0.0)
                {
                    t2[i] = (a-1.0)/2.0;
                }
                else
                {
                    t2[i] = (a+1.0)/2.0;
                }
            }
        }
        else
        {
            t2[i] = t1[i];
        }
    }
    temp1 = pow((sin(sqrt(pow(t2[0],2.0)+pow(t2[1],2.0)))),2.0);
    temp2 = 1.0 + 0.001*(pow(t2[0],2.0)+pow(t2[1],2.0));
    res = 0.5 + (temp1-0.5)/(pow(temp2,2.0));
    return (res);
}

/* Code to evaluate rastrigin's function and rounding-off variables */
long double nc_rastrigin (long double *x)
{
    int i;
    int a;
    long double b;
    long double res;
    for (i=0; i<nreal; i++)
    {
        if (fabs(x[i]) >= 0.5)
        {
            res = 2.0*x[i];
            a = res;
            b = fabs(res-a);
            if (b<0.5)
            {
                temp_x4[i] = a/2.0;
            }
            else
            {
                if (res<=0.0)
                {
                    temp_x4[i] = (a-1.0)/2.0;
                }
                else
                {
                    temp_x4[i] = (a+1.0)/2.0;
                }
            }
        }
        else
        {
            temp_x4[i] = x[i];
        }
    }
    res = 0.0;
    for (i=0; i<nreal; i++)
    {
        res += (temp_x4[i]*temp_x4[i] - 10.0*cos(2.0*PI*temp_x4[i]) + 10.0);
    }
    return (res);
}

//===========================================================================
//===========================================================================


/* Code to allocate memory to global variables being used in evaluation of functions */
void allocate_memory ()
{
    int i, j, k;
    norm_x = new long double[nreal];
    norm_f = new long double[nfunc];
    trans_x = new long double[nreal];
    basic_f = new long double[nfunc];
    temp_x1 = new long double[nreal];
    temp_x2 = new long double[nreal];
    temp_x3 = new long double[nreal];
    temp_x4 = new long double[nreal];
    weight = new long double[nfunc];
    sigma =  new long double[nfunc];
    lambda = new long double[nfunc];
    bias = new long double[nfunc];
    o = new long double*[nfunc];
    l = new long double**[nfunc];
    g = new long double*[nreal];
    for (i=0; i<nfunc; i++)
    {
        o[i] = new long double[nreal];
        l[i] = new long double*[nreal];
    }
    for (i=0; i<nreal; i++)
    {
        g[i] = new long double[nreal];
    }
    for (i=0; i<nfunc; i++)
    {
        for (j=0; j<nreal; j++)
        {
            l[i][j] = new long double[nreal];
        }
    }
    /* Do some trivial (common) initialization here itself */
	C = 2000.0;
    for (i=0; i<nreal; i++)
    {
        norm_x[i] = 5.0;
        trans_x[i] = 0.0;
        temp_x1[i] = 0.0;
        temp_x2[i] = 0.0;
        temp_x3[i] = 0.0;
        temp_x4[i] = 0.0;
        for (j=0; j<nreal; j++)
        {
            if (i==j)
            {
                g[i][j]=1.0;
            }
            else
            {
                g[i][j]=0.0;
            }
        }
    }
    for (i=0; i<nfunc; i++)
    {
        basic_f[i] = 0.0;
        norm_f[i] = 0.0;
        weight[i] = 1.0/(long double)nfunc;
        sigma[i] = 1.0;
        lambda[i] = 1.0;
        bias[i] = 100.0*(long double)i;
        for (j=0; j<nreal; j++)
        {
            o[i][j] = 0.0;
            for (k=0; k<nreal; k++)
            {
                if (j==k)
                {
                    l[i][j][k] = 1.0;
                }
                else
                {
                    l[i][j][k] = 0.0;
                }
            }
        }
    }
    return;
}

/* Code to transform a variable vector based on function index 'count' */
void transform (long double *x, int count)
{
    int i, j;
    for (i=0; i<nreal; i++)
    {
        temp_x1[i] = x[i] - o[count][i];
    }
    for (i=0; i<nreal; i++)
    {
        temp_x2[i] = temp_x1[i]/lambda[count];
    }
    for (j=0; j<nreal; j++)
    {
        temp_x3[j] = 0.0;
        for (i=0; i<nreal; i++)
        {
            temp_x3[j] += g[i][j]*temp_x2[i];
        }
    }
    for (j=0; j<nreal; j++)
    {
        trans_x[j] = 0.0;
        for (i=0; i<nreal; i++)
        {
            trans_x[j] += l[count][i][j]*temp_x3[i];
        }
    }
    return;
}

/* Code to transform a vector (with elements 5.0) based on function index 'count' */
void transform_norm (int count)
{
    int i, j;
    for (i=0; i<nreal; i++)
    {
        temp_x2[i] = 5.0/lambda[count];
    }
    for (j=0; j<nreal; j++)
    {
        temp_x3[j] = 0.0;
        for (i=0; i<nreal; i++)
        {
            temp_x3[j] += g[i][j]*temp_x2[i];
        }
    }
    for (j=0; j<nreal; j++)
    {
        trans_x[j] = 0.0;
        for (i=0; i<nreal; i++)
        {
            trans_x[j] += l[count][i][j]*temp_x3[i];
        }
    }
    return;
}

/* Code to compute the weights for a variable vector */
void calc_weight (long double *x)
{
    int i, j;
    long double sum;
    long double max;
    max = -INF;
    for (i=0; i<nfunc; i++)
    {
        sum = 0.0;
        for (j=0; j<nreal; j++)
        {
            sum += (x[j]-o[i][j])*(x[j]-o[i][j]);
        }
        weight[i] = exp(-(sum)/(2.0*nreal*sigma[i]*sigma[i]));
        max = maximum(max,weight[i]);
    }
    sum = 0.0;
    for (i=0; i<nfunc; i++)
    {
        if (weight[i]!=max)
        {
            weight[i] *= (1.0 - pow(max,10.0));
        }
        sum += weight[i];
    }
    if (sum==0.0)
    {
        for (i=0; i<nfunc; i++)
        {
            weight[i] = 1.0/(long double)nfunc;
        }
    }
    else
    {
        for (i=0; i<nfunc; i++)
        {
            weight[i] /= sum;
        }
    }
    return;
}

/* Code to free the allocated memory */
void free_memory()
{
    int i, j;
    delete [] norm_x;
    delete [] norm_f;
    delete [] trans_x;
    delete [] basic_f;
    delete [] temp_x1;
    delete [] temp_x2;
    delete [] temp_x3;
    delete [] temp_x4;
    delete [] weight;
    delete [] sigma;
    delete [] lambda;
    delete [] bias;
    for (i=0; i<nfunc; i++)
    {
        for (j=0; j<nreal; j++)
        {
            delete [] (l[i][j]);
        }
    }
    for (i=0; i<nfunc; i++)
    {
        delete [] (o[i]);
        delete [] (l[i]);
    }
    for (i=0; i<nreal; i++)
    {
        delete [] (g[i]);
    }
    delete [] o;
    delete [] l;
    delete [] g;
    return;
}

//======================================================================
//======================================================================
/* Test problem specific variable definitions */

/*  Initialize function definitions */
/*  Benchmark function definitions */
/*  Various benchmark function definitions */
# ifdef f1
void initialize()
{
    int i, j;
    fstream input;
    input.open("input_data/sphere_func_data.txt", ios::in);
    if(!input)
    {
        cerr << "File could not be opened for reading" << endl;
        exit(1);
    }
    for (i=0; i<nfunc; i++)
    {
        for (j=0; j<nreal; j++)
        {
            input >> o[i][j];
            o[i][j] -= 1.0;
            cout << scientific << "\nO[" << i+1 << "][" << j+1 << "] = " << o[i][j];
        }
    }
    input.close();

    bias[0] = -450.0;
    return;
}

long double calc_benchmark_func(long double *x)
{
    long double res;
    transform (x, 0);
    basic_f[0] = calc_sphere (trans_x);
    res = basic_f[0] + bias[0];
    return (res);
}
#endif

# ifdef f2
void initialize()
{
    int i, j;
    fstream input;
    input.open("input_data/schwefel_102_data.txt", ios::in);
    if(!input)
    {
        cerr << "File could not be opened for reading" << endl;
        exit(1);
    }
    for (i=0; i<nfunc; i++)
    {
        for (j=0; j<nreal; j++)
        {
            input >> o[i][j];
            o[i][j] -= 1.0;
            cout << scientific << "\nO[" << i+1 << "][" << j+1 << "] = " << o[i][j];
        }
    }
    input.close();

    bias[0] = -450.0;
    return;
}

long double calc_benchmark_func(long double *x)
{
    long double res;
    transform (x, 0);
    basic_f[0] = calc_schwefel (trans_x);
    res = basic_f[0] + bias[0];
    return (res);
}
#endif

# ifdef f3
void initialize()
{
    int i, j;
    fstream input;
    if (nreal==2)    input.open("input_data/elliptic_M_D2.txt", ios::in);
    if (nreal==10)   input.open("input_data/elliptic_M_D2.txt", ios::in);
    if (nreal==30)   input.open("input_data/elliptic_M_D2.txt", ios::in);
    if (nreal==50)   input.open("input_data/elliptic_M_D2.txt", ios::in);
    if(!input)
    {
        cerr << "File could not be opened for reading" << endl;
        exit(1);
    }
    for (i=0; i<nreal; i++)
    {
        for (j=0; j<nreal; j++)
        {
            input >> g[i][j];
            cout << scientific << "\nM[" << i+1 << "][" << j+1 << "] = " << g[i][j];
        }
    }
    input.close();

    input.open("input_data/high_cond_elliptic_rot_data.txt", ios::in);
    if(!input)
    {
        cerr << "File could not be opened for reading" << endl;
        exit(1);
    }
    for (i=0; i<nfunc; i++)
    {
        for (j=0; j<nreal; j++)
        {
            input >> o[i][j];
            cout << scientific << "\nO[" << i+1 << "][" << j+1 << "] = " << o[i][j];
        }
    }
    input.close();

    bias[0] = -450.0;
    return;
}

long double calc_benchmark_func(long double *x)
{
    int i;
    long double res;
    transform (x, 0);
    basic_f[0] = 0.0;
    for (i=0; i<nreal; i++)
    {
        basic_f[0] += trans_x[i]*trans_x[i]*pow(1.0e6,i/(nreal-1.0));
    }
    res = basic_f[0] + bias[0];
    return (res);
}
#endif

# ifdef f4
void initialize()
{
    int i, j;
    fstream input;
    input.open("input_data/schwefel_102_data.txt", ios::in);
    if(!input)
    {
        cerr << "File could not be opened for reading" << endl;
        exit(1);
    }
    for (i=0; i<nfunc; i++)
    {
        for (j=0; j<nreal; j++)
        {
            input >> o[i][j];
            o[i][j] -= 1.0;
            cout << scientific << "\nO[" << i+1 << "][" << j+1 << "] = " << o[i][j];
        }
    }
    input.close();

    bias[0] = -450.0;
    return;
}

long double calc_benchmark_func(long double *x)
{
    long double res;
    transform (x, 0);
    basic_f[0] = calc_schwefel(trans_x)*(1.0 + 0.4*fabs(randomnormaldeviate()));
    res = basic_f[0] + bias[0];
    return (res);
}
#endif

# ifdef f5
long double **A;
long double *B;

void initialize()
{
    int i, j;
    int index;
    fstream input;
    A = (long double **)malloc(nreal*sizeof(long double));
    for (i=0; i<nreal; i++)
    {
        A[i] = (long double *)malloc(nreal*sizeof(long double));
    }
    B = (long double *)malloc(nreal*sizeof(long double));
    input.open("input_data/schwefel_206_data.txt", ios::in);
    for (i=0; i<nfunc; i++)
    {
        for (j=0; j<nreal; j++)
        {
            input >> o[i][j];
            cout << scientific << "\nO[" << i+1 << "][" << j+1 << "] = " << o[i][j];
        }
        if(nreal<100)
        {
            long double temp;
            for(j=nreal; j<100; j++)
            {
                input >> temp;
            }
        }
    }
    for (i=0; i<nreal; i++)
    {
        for (j=0; j<nreal; j++)
        {
            input >> A[i][j];
            cout << scientific << "\nA[" << i+1 << "][" << j+1 << "] = " << A[i][j];
        }
        if(nreal<100)
        {
            long double temp;
            for(j=nreal; j<100; j++)
            {
                input >> temp;
            }
        }
    }
    input.close();

    if (nreal%4==0)
    {
        index = nreal/4;
    }
    else
    {
        index = nreal/4 + 1;
    }
    for (i=0; i<index; i++)
    {
        o[0][i] = -100.0;
    }
    index = (3*nreal)/4 - 1;
    for (i=index; i<nreal; i++)
    {
        o[0][i] = 100.0;
    }
    for (i=0; i<nreal; i++)
    {
        B[i] = 0.0;
        for (j=0; j<nreal; j++)
        {
            B[i] += A[i][j]*o[0][j];
        }
    }
    bias[0] = -310.0;
    return;
}

long double calc_benchmark_func(long double *x)
{
    int i, j;
    long double res;
    basic_f[0] = -INF;
    for (i=0; i<nreal; i++)
    {
        res=0.0;
        for (j=0; j<nreal; j++)
        {
            res += A[i][j]*x[j];
        }
        res = fabs(res-B[i]);
        if (basic_f[0] < res)
        {
            basic_f[0] = res;
        }
    }
    res = basic_f[0] + bias[0];
    return (res);
}
#endif

# ifdef f6
void initialize()
{
    int i, j;
    fstream input;
    input.open("input_data/rosenbrock_func_data.txt", ios::in);
    if(!input)
    {
        cerr << "File could not be opened for reading" << endl;
        exit(1);
    }
    for (i=0; i<nfunc; i++)
    {
        for (j=0; j<nreal; j++)
        {
            input >> o[i][j];
            o[i][j] -= 1.0;
            cout << scientific << "\nO[" << i+1 << "][" << j+1 << "] = " << o[i][j];
        }
    }
    input.close();

    bias[0] = 390.0;
    return;
}

long double calc_benchmark_func(long double *x)
{
    long double res;
    transform (x, 0);
    basic_f[0] = calc_rosenbrock(trans_x);
    res = basic_f[0] + bias[0];
    return (res);
}
#endif

# ifdef f7
void initialize()
{
    int i, j;
    fstream input;
    if (nreal==2)    input.open("input_data/griewank_M_D2.txt", ios::in);
    if (nreal==10)   input.open("input_data/griewank_M_D10.txt", ios::in);
    if (nreal==30)   input.open("input_data/griewank_M_D30.txt", ios::in);
    if (nreal==50)   input.open("input_data/griewank_M_D50.txt", ios::in);
    if(!input)
    {
        cerr << "File could not be opened for reading" << endl;
        exit(1);
    }
    for (i=0; i<nreal; i++)
    {
        for (j=0; j<nreal; j++)
        {
            input >> g[i][j];
            cout << scientific << "\nM[" << i+1 << "][" << j+1 << "] = " << g[i][j];
        }
    }
    input.close();

    input.open("input_data/griewank_func_data.txt", ios::in);
    if(!input)
    {
        cerr << "File could not be opened for reading" << endl;
        exit(1);
    }
    for (i=0; i<nfunc; i++)
    {
        for (j=0; j<nreal; j++)
        {
            input >> o[i][j];
            cout << scientific << "\nO[" << i+1 << "][" << j+1 << "] = " << o[i][j];
        }
    }
    input.close();

    bias[0] = -180.0;
    return;
}

long double calc_benchmark_func(long double *x)
{
    long double res;
    transform (x, 0);
    basic_f[0] = calc_griewank(trans_x);
    res = basic_f[0] + bias[0];
    return (res);
}
#endif

# ifdef f8
void initialize()
{
    int i, j;
    int index;
    fstream input;
    if (nreal==2)    input.open("input_data/ackley_M_D2.txt", ios::in);
    if (nreal==10)   input.open("input_data/ackley_M_D10.txt", ios::in);
    if (nreal==30)   input.open("input_data/ackley_M_D30.txt", ios::in);
    if (nreal==50)   input.open("input_data/ackley_M_D50.txt", ios::in);
    if(!input)
    {
        cerr << "File could not be opened for reading" << endl;
        exit(1);
    }
    for (i=0; i<nreal; i++)
    {
        for (j=0; j<nreal; j++)
        {
            input >> g[i][j];
            cout << scientific << "\nM[" << i+1 << "][" << j+1 << "] = " << g[i][j];
        }
    }
    input.close();

    input.open("input_data/ackley_func_data.txt", ios::in);
    if(!input)
    {
        cerr << "File could not be opened for reading" << endl;
        exit(1);
    }
    for (i=0; i<nfunc; i++)
    {
        for (j=0; j<nreal; j++)
        {
            input >> o[i][j];
            cout << scientific << "\nO[" << i+1 << "][" << j+1 << "] = " << o[i][j];
        }
    }
    input.close();

    index = nreal/2;
    for (i=1; i<=index; i++)
    {
        o[0][2*i-2] = -32.0;
    }
    bias[0] = -140.0;
    return;
}

long double calc_benchmark_func(long double *x)
{
    long double res;
    transform (x, 0);
    basic_f[0] = calc_ackley(trans_x);
    res = basic_f[0] + bias[0];
    return (res);
}
#endif

# ifdef f9
void initialize()
{
    int i, j;
    fstream input;
    input.open("input_data/rastrigin_func_data.txt", ios::in);
    if(!input)
    {
        cerr << "File could not be opened for reading" << endl;
        exit(1);
    }
    for (i=0; i<nfunc; i++)
    {
        for (j=0; j<nreal; j++)
        {
            input >> o[i][j];
            cout << scientific << "\nO[" << i+1 << "][" << j+1 << "] = " << o[i][j];
        }
    }
    input.close();

    bias[0] = -330.0;
    return;
}

long double calc_benchmark_func(long double *x)
{
    long double res;
    transform (x, 0);
    basic_f[0] = calc_rastrigin(trans_x);
    res = basic_f[0] + bias[0];
    return (res);
}
#endif

# ifdef f10
void initialize()
{
    int i, j;
    fstream input;
    if (nreal==2)    input.open("input_data/rastrigin_M_D2.txt", ios::in);
    if (nreal==10)   input.open("input_data/rastrigin_M_D10.txt", ios::in);
    if (nreal==30)   input.open("input_data/rastrigin_M_D30.txt", ios::in);
    if (nreal==50)   input.open("input_data/rastrigin_M_D50.txt", ios::in);
    if(!input)
    {
        cerr << "File could not be opened for reading" << endl;
        exit(1);
    }
    for (i=0; i<nreal; i++)
    {
        for (j=0; j<nreal; j++)
        {
            input >> g[i][j];
            cout << scientific << "\nM[" << i+1 << "][" << j+1 << "] = " << g[i][j];
        }
    }
    input.close();

    input.open("input_data/rastrigin_func_data.txt", ios::in);
    if(!input)
    {
        cerr << "File could not be opened for reading" << endl;
        exit(1);
    }
    for (i=0; i<nfunc; i++)
    {
        for (j=0; j<nreal; j++)
        {
            input >> o[i][j];
            cout << scientific << "\nO[" << i+1 << "][" << j+1 << "] = " << o[i][j];
        }
    }
    input.close();

    bias[0] = -330.0;
    return;
}

long double calc_benchmark_func(long double *x)
{
    long double res;
    transform (x, 0);
    basic_f[0] = calc_rastrigin(trans_x);
    res = basic_f[0] + bias[0];
    return (res);
}
#endif

# ifdef f11
void initialize()
{
    int i, j;
    fstream input;
    if (nreal==2)    input.open("input_data/weierstrass_M_D2.txt", ios::in);
    if (nreal==10)   input.open("input_data/weierstrass_M_D10.txt", ios::in);
    if (nreal==30)   input.open("input_data/weierstrass_M_D30.txt", ios::in);
    if (nreal==50)   input.open("input_data/weierstrass_M_D50.txt", ios::in);
    if(!input)
    {
        cerr << "File could not be opened for reading" << endl;
        exit(1);
    }
    for (i=0; i<nreal; i++)
    {
        for (j=0; j<nreal; j++)
        {
            input >> g[i][j];
            cout << scientific << "\nM[" << i+1 << "][" << j+1 << "] = " << g[i][j];
        }
    }
    input.close();

    input.open("input_data/weierstrass_data.txt", ios::in);
    if(!input)
    {
        cerr << "File could not be opened for reading" << endl;
        exit(1);
    }
    for (i=0; i<nfunc; i++)
    {
        for (j=0; j<nreal; j++)
        {
            input >> o[i][j];
            cout << scientific << "\nO[" << i+1 << "][" << j+1 << "] = " << o[i][j];
        }
    }
    input.close();

    bias[0] = 90.0;
    return;
}

long double calc_benchmark_func(long double *x)
{
    int i;
    long double res;
    for (i=0; i<nreal; i++)
    {
        norm_x[i] = 0.0;
    }
    transform (x, 0);
    basic_f[0] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    res = basic_f[0] + bias[0];
    return (res);
}
#endif

# ifdef f12
long double **A;
long double **B;
long double *alpha;

void initialize()
{
    int i, j;
    fstream input;
    A = (long double **)malloc(nreal*sizeof(long double));
    B = (long double **)malloc(nreal*sizeof(long double));
    alpha = (long double *)malloc(nreal*sizeof(long double));
    for (i=0; i<nreal; i++)
    {
        A[i] = (long double *)malloc(nreal*sizeof(long double));
        B[i] = (long double *)malloc(nreal*sizeof(long double));
    }
    input.open("input_data/schwefel_213_data.txt", ios::in);
    if(!input)
    {
        cerr << "File could not be opened for reading" << endl;
        exit(1);
    }
    /* Reading A */
    for (i=0; i<nreal; i++)
    {
        for (j=0; j<nreal; j++)
        {
            input >> A[i][j];
            cout << scientific << "\nA[" << i+1 << "][" << j+1 << "] = " << A[i][j];
        }
    }
    if (i!=100)
    {
        long double temp;
        for (i=nreal; i<100; i++)
        {
            input >> temp;
        }
    }
    cout << endl << endl;
    /* Reading B */
    for (i=0; i<nreal; i++)
    {
        for (j=0; j<nreal; j++)
        {
            input >> B[i][j];
            cout << scientific << "\nB[" << i+1 << "][" << j+1 << "] = " << B[i][j];
        }
    }
    if (i!=100)
    {
        long double temp;
        for (i=nreal; i<100; i++)
        {
            input >> temp;
        }
    }
    cout << endl << endl;
    /* Reading alpha */
    for (i=0; i<nreal; i++)
    {
        input >> alpha[i];
        cout << scientific << "\nalpha[" << i+1 << "] = " << alpha[i];
    }
    cout << endl << endl;
    input.close();

    bias[0] = -460.0;
    return;
}

long double calc_benchmark_func(long double *x)
{
    long double res;
    long double sum1, sum2;
    int i, j;
    basic_f[0] = 0.0;
    for (i=0; i<nreal; i++)
    {
        sum1 = 0.0;
        sum2 = 0.0;
        for (j=0; j<nreal; j++)
        {
            sum1 += A[i][j]*sin(alpha[j]) + B[i][j]*cos(alpha[j]);
            sum2 += A[i][j]*sin(x[j]) + B[i][j]*cos(x[j]);
        }
        basic_f[0] += pow((sum1-sum2),2.0);
    }
    res = basic_f[0] + bias[0];
    return (res);
}
#endif

# ifdef f13
void initialize()
{
    int i, j;
    fstream input;
    input.open("input_data/EF8F2_func_data.txt", ios::in);
    if(!input)
    {
        cerr << "File could not be opened for reading" << endl;
        exit(1);
    }
    for (i=0; i<nfunc; i++)
    {
        for (j=0; j<nreal; j++)
        {
            input >> o[i][j];
            o[i][j] -= 1.0;
            cout << scientific << "\nO[" << i+1 << "][" << j+1 << "] = " << o[i][j];
        }
    }
    cout << endl << endl;
    input.close();

    bias[0] = -130.0;
    return;
}

long double calc_benchmark_func(long double *x)
{
    int i;
    long double temp;
    long double res;
    transform(x,0);
    res = 0.0;
    for (i=0; i<nreal-1; i++)
    {
        temp = 100.0*pow((trans_x[i]*trans_x[i]-trans_x[i+1]),2.0) + 1.0*pow((trans_x[i]-1.0),2.0);
        res += (temp*temp)/4000.0 - cos(temp) + 1.0;
    }
    temp = 100.0*pow((trans_x[nreal-1]*trans_x[nreal-1]-trans_x[0]),2.0) + 1.0*pow((trans_x[nreal-1]-1.0),2.0);
    res += (temp*temp)/4000.0 - cos(temp) + 1.0 + bias[0];
    return (res);
}
#endif

# ifdef f14
void initialize()
{
    int i, j;
    fstream input;
    if (nreal==2)    input.open("input_data/E_ScafferF6_M_D2.txt", ios::in);
    if (nreal==10)   input.open("input_data/E_ScafferF6_M_D10.txt", ios::in);
    if (nreal==30)   input.open("input_data/E_ScafferF6_M_D30.txt", ios::in);
    if (nreal==50)   input.open("input_data/E_ScafferF6_M_D50.txt", ios::in);
    if(!input)
    {
        cerr << "File could not be opened for reading" << endl;
        exit(1);
    }
    for (i=0; i<nfunc; i++)
    {
        for (j=0; j<nreal; j++)
        {
            input >> g[i][j];
            cout << scientific << "\nM[" << i+1 << "][" << j+1 << "] = " << g[i][j];
        }
    }
    cout << endl << endl;
    input.close();

    input.open("input_data/E_ScafferF6_func_data.txt", ios::in);
    if(!input)
    {
        cerr << "File could not be opened for reading" << endl;
        exit(1);
    }
    for (i=0; i<nfunc; i++)
    {
        for (j=0; j<nreal; j++)
        {
            input >> o[i][j];
            cout << scientific << "\nO[" << i+1 << "][" << j+1 << "] = " << o[i][j];
        }
    }
    cout << endl << endl;
    input.close();

    bias[0] = -300.0;
    return;
}

long double calc_benchmark_func(long double *x)
{
    int i;
    long double temp1, temp2;
    long double res;
    transform(x,0);
    res = 0.0;
    for (i=0; i<nreal-1; i++)
    {
        temp1 = pow((sin(sqrt(pow(trans_x[i],2.0)+pow(trans_x[i+1],2.0)))),2.0);
        temp2 = 1.0 + 0.001*(pow(trans_x[i],2.0)+pow(trans_x[i+1],2.0));
        res += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
    }
    temp1 = pow((sin(sqrt(pow(trans_x[nreal-1],2.0)+pow(trans_x[0],2.0)))),2.0);
    temp2 = 1.0 + 0.001*(pow(trans_x[nreal-1],2.0)+pow(trans_x[0],2.0));
    res += 0.5 + (temp1-0.5)/(pow(temp2,2.0)) + bias[0];
    return (res);
}
#endif

# ifdef f15
void initialize()
{
    int i, j;
    fstream input;
    input.open("input_data/hybrid_func1_data.txt", ios::in);
    if(!input)
    {
        cerr << "File could not be opened for reading" << endl;
        exit(1);
    }
    for (i=0; i<nfunc; i++)
    {
        for (j=0; j<nreal; j++)
        {
            input >> o[i][j];
            cout << scientific << "\nO[" << i+1 << "][" << j+1 << "] = " << o[i][j];
        }
    }
    cout << endl << endl;
    input.close();

    lambda[0] = 1.0;
    lambda[1] = 1.0;
    lambda[2] = 10.0;
    lambda[3] = 10.0;
    lambda[4] = 1.0/12.0;
    lambda[5] = 1.0/12.0;
    lambda[6] = 5.0/32.0;
    lambda[7] = 5.0/32.0;
    lambda[8] = 1.0/20.0;
    lambda[9] = 1.0/20.0;
    global_bias = 120.0;
    return;
}

void calc_benchmark_norm()
{
    int i;
    transform_norm (0);    norm_f[0] = calc_rastrigin(trans_x);
    transform_norm (1);    norm_f[1] = calc_rastrigin(trans_x);
    for (i=0; i<nreal; i++)
    {
        norm_x[i] = 0.0;
    }
    transform_norm (2);    norm_f[2] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform_norm (3);    norm_f[3] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform_norm (4);    norm_f[4] = calc_griewank(trans_x);
    transform_norm (5);    norm_f[5] = calc_griewank(trans_x);
    transform_norm (6);    norm_f[6] = calc_ackley(trans_x);
    transform_norm (7);    norm_f[7] = calc_ackley(trans_x);
    transform_norm (8);    norm_f[8] = calc_sphere(trans_x);
    transform_norm (9);    norm_f[9] = calc_sphere(trans_x);
    return;
}

long double calc_benchmark_func(long double *x)
{
    int i;
    long double res;
    transform (x, 0);    basic_f[0] = calc_rastrigin(trans_x);
    transform (x, 1);    basic_f[1] = calc_rastrigin(trans_x);
    for (i=0; i<nreal; i++)
    {
        norm_x[i] = 0.0;
    }
    transform (x, 2);    basic_f[2] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform (x, 3);    basic_f[3] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform (x, 4);    basic_f[4] = calc_griewank(trans_x);
    transform (x, 5);    basic_f[5] = calc_griewank(trans_x);
    transform (x, 6);    basic_f[6] = calc_ackley(trans_x);
    transform (x, 7);    basic_f[7] = calc_ackley(trans_x);
    transform (x, 8);    basic_f[8] = calc_sphere(trans_x);
    transform (x, 9);    basic_f[9] = calc_sphere(trans_x);
    for (i=0; i<nfunc; i++)
    {
        basic_f[i] *= C/norm_f[i];
    }
    calc_weight(x);
    res = global_bias;
    for (i=0; i<nfunc; i++)
    {
        res += weight[i]*(basic_f[i]+bias[i]);
    }
    return (res);
}
#endif

# ifdef f16
void initialize()
{
    int i, j, k;
    fstream input;
    input.open("input_data/hybrid_func1_data.txt", ios::in);
    if(!input)
    {
        cerr << "File could not be opened for reading" << endl;
        exit(1);
    }
    for (i=0; i<nfunc; i++)
    {
        for (j=0; j<nreal; j++)
        {
            input >> o[i][j];
            cout << scientific << "\nO[" << i+1 << "][" << j+1 << "] = " << o[i][j];
        }
    }
    cout << endl << endl;
    input.close();

    if (nreal==2)    input.open("input_data/hybrid_func1_M_D2.txt", ios::in);
    if (nreal==10)   input.open("input_data/hybrid_func1_M_D10.txt", ios::in);
    if (nreal==30)   input.open("input_data/hybrid_func1_M_D30.txt", ios::in);
    if (nreal==50)   input.open("input_data/hybrid_func1_M_D50.txt", ios::in);
    if(!input)
    {
        cerr << "File could not be opened for reading" << endl;
        exit(1);
    }
    for (i=0; i<nfunc; i++)
    {
        for (j=0; j<nreal; j++)
        {
            for (k=0; k<nreal; k++)
            {
                input >> l[i][j][k];
                cout << scientific << " M[" << i+1 << "][" << j+1 << "][" << k+1 << "] = " << l[i][j][k] << endl;
            }
        }
    }
    input.close();

    lambda[0] = 1.0;
    lambda[1] = 1.0;
    lambda[2] = 10.0;
    lambda[3] = 10.0;
    lambda[4] = 1.0/12.0;
    lambda[5] = 1.0/12.0;
    lambda[6] = 5.0/32.0;
    lambda[7] = 5.0/32.0;
    lambda[8] = 1.0/20.0;
    lambda[9] = 1.0/20.0;
    global_bias = 120.0;
    return;
}

void calc_benchmark_norm()
{
    int i;
    transform_norm (0);    norm_f[0] = calc_rastrigin(trans_x);
    transform_norm (1);    norm_f[1] = calc_rastrigin(trans_x);
    for (i=0; i<nreal; i++)
    {
        norm_x[i] = 0.0;
    }
    transform_norm (2);    norm_f[2] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform_norm (3);    norm_f[3] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform_norm (4);    norm_f[4] = calc_griewank(trans_x);
    transform_norm (5);    norm_f[5] = calc_griewank(trans_x);
    transform_norm (6);    norm_f[6] = calc_ackley(trans_x);
    transform_norm (7);    norm_f[7] = calc_ackley(trans_x);
    transform_norm (8);    norm_f[8] = calc_sphere(trans_x);
    transform_norm (9);    norm_f[9] = calc_sphere(trans_x);
    return;
}

long double calc_benchmark_func(long double *x)
{
    int i;
    long double res;
    transform (x, 0);    basic_f[0] = calc_rastrigin(trans_x);
    transform (x, 1);    basic_f[1] = calc_rastrigin(trans_x);
    for (i=0; i<nreal; i++)
    {
        norm_x[i] = 0.0;
    }
    transform (x, 2);    basic_f[2] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform (x, 3);    basic_f[3] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform (x, 4);    basic_f[4] = calc_griewank(trans_x);
    transform (x, 5);    basic_f[5] = calc_griewank(trans_x);
    transform (x, 6);    basic_f[6] = calc_ackley(trans_x);
    transform (x, 7);    basic_f[7] = calc_ackley(trans_x);
    transform (x, 8);    basic_f[8] = calc_sphere(trans_x);
    transform (x, 9);    basic_f[9] = calc_sphere(trans_x);
    for (i=0; i<nfunc; i++)
    {
        basic_f[i] *= C/norm_f[i];
    }
    calc_weight(x);
    res = global_bias;
    for (i=0; i<nfunc; i++)
    {
        res += weight[i]*(basic_f[i]+bias[i]);
    }
    return (res);
}
#endif

# ifdef f17
void initialize()
{
    int i, j, k;
    fstream input;
    input.open("input_data/hybrid_func1_data.txt", ios::in);
    if(!input)
    {
        cerr << "File could not be opened for reading" << endl;
        exit(1);
    }
    for (i=0; i<nfunc; i++)
    {
        for (j=0; j<nreal; j++)
        {
            input >> o[i][j];
            cout << scientific << "\nO[" << i+1 << "][" << j+1 << "] = " << o[i][j];
        }
    }
    cout << endl << endl;
    input.close();

    if (nreal==2)    input.open("input_data/hybrid_func1_M_D2.txt", ios::in);
    if (nreal==10)   input.open("input_data/hybrid_func1_M_D10.txt", ios::in);
    if (nreal==30)   input.open("input_data/hybrid_func1_M_D30.txt", ios::in);
    if (nreal==50)   input.open("input_data/hybrid_func1_M_D50.txt", ios::in);
    if(!input)
    {
        cerr << "File could not be opened for reading" << endl;
        exit(1);
    }
    for (i=0; i<nfunc; i++)
    {
        for (j=0; j<nreal; j++)
        {
            for (k=0; k<nreal; k++)
            {
                input >> l[i][j][k];
                cout << scientific << " M[" << i+1 << "][" << j+1 << "][" << k+1 << "] = " << l[i][j][k] << endl;
            }
        }
    }
    input.close();

    lambda[0] = 1.0;
    lambda[1] = 1.0;
    lambda[2] = 10.0;
    lambda[3] = 10.0;
    lambda[4] = 1.0/12.0;
    lambda[5] = 1.0/12.0;
    lambda[6] = 5.0/32.0;
    lambda[7] = 5.0/32.0;
    lambda[8] = 1.0/20.0;
    lambda[9] = 1.0/20.0;
    global_bias = 120.0;
    return;
}

void calc_benchmark_norm()
{
    int i;
    transform_norm (0);    norm_f[0] = calc_rastrigin(trans_x);
    transform_norm (1);    norm_f[1] = calc_rastrigin(trans_x);
    for (i=0; i<nreal; i++)
    {
        norm_x[i] = 0.0;
    }
    transform_norm (2);    norm_f[2] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform_norm (3);    norm_f[3] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform_norm (4);    norm_f[4] = calc_griewank(trans_x);
    transform_norm (5);    norm_f[5] = calc_griewank(trans_x);
    transform_norm (6);    norm_f[6] = calc_ackley(trans_x);
    transform_norm (7);    norm_f[7] = calc_ackley(trans_x);
    transform_norm (8);    norm_f[8] = calc_sphere(trans_x);
    transform_norm (9);    norm_f[9] = calc_sphere(trans_x);
    return;
}

long double calc_benchmark_func(long double *x)
{
    int i;
    long double res;
    transform (x, 0);    basic_f[0] = calc_rastrigin(trans_x);
    transform (x, 1);    basic_f[1] = calc_rastrigin(trans_x);
    for (i=0; i<nreal; i++)
    {
        norm_x[i] = 0.0;
    }
    transform (x, 2);    basic_f[2] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform (x, 3);    basic_f[3] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform (x, 4);    basic_f[4] = calc_griewank(trans_x);
    transform (x, 5);    basic_f[5] = calc_griewank(trans_x);
    transform (x, 6);    basic_f[6] = calc_ackley(trans_x);
    transform (x, 7);    basic_f[7] = calc_ackley(trans_x);
    transform (x, 8);    basic_f[8] = calc_sphere(trans_x);
    transform (x, 9);    basic_f[9] = calc_sphere(trans_x);
    for (i=0; i<nfunc; i++)
    {
        basic_f[i] *= C/norm_f[i];
    }
    calc_weight(x);
    res = 0.0;
    for (i=0; i<nfunc; i++)
    {
        res += weight[i]*(basic_f[i]+bias[i]);
    }
    res = res*(1.0 + 0.2*fabs(randomnormaldeviate())) + global_bias;
    return (res);
}
#endif

# ifdef f18
void initialize()
{
    int i, j, k;
    fstream input;
    input.open("input_data/hybrid_func2_data.txt", ios::in);
    if(!input)
    {
        cerr << "File could not be opened for reading" << endl;
        exit(1);
    }
    for (i=0; i<nfunc; i++)
    {
        for (j=0; j<nreal; j++)
        {
            input >> o[i][j];
            cout << scientific << "\nO[" << i+1 << "][" << j+1 << "] = " << o[i][j];
        }
    }
    cout << endl << endl;
    input.close();

    if (nreal==2)    input.open("input_data/hybrid_func2_M_D2.txt", ios::in);
    if (nreal==10)   input.open("input_data/hybrid_func2_M_D10.txt", ios::in);
    if (nreal==30)   input.open("input_data/hybrid_func2_M_D30.txt", ios::in);
    if (nreal==50)   input.open("input_data/hybrid_func2_M_D50.txt", ios::in);
    if(!input)
    {
        cerr << "File could not be opened for reading" << endl;
        exit(1);
    }
    for (i=0; i<nfunc; i++)
    {
        for (j=0; j<nreal; j++)
        {
            for (k=0; k<nreal; k++)
            {
                input >> l[i][j][k];
                cout << scientific << " M[" << i+1 << "][" << j+1 << "][" << k+1 << "] = " << l[i][j][k] << endl;
            }
        }
    }
    input.close();

    for (i=0; i<nreal; i++)
    {
        o[nfunc-1][i] = 0.0;
    }
    sigma[0] = 1.0;
    sigma[1] = 2.0;
    sigma[2] = 1.5;
    sigma[3] = 1.5;
    sigma[4] = 1.0;
    sigma[5] = 1.0;
    sigma[6] = 1.5;
    sigma[7] = 1.5;
    sigma[8] = 2.0;
    sigma[9] = 2.0;
    lambda[0] = 5.0/16.0;
    lambda[1] = 5.0/32.0;
    lambda[2] = 2.0;
    lambda[3] = 1.0;
    lambda[4] = 1.0/10.0;
    lambda[5] = 1.0/20.0;
    lambda[6] = 20.0;
    lambda[7] = 10.0;
    lambda[8] = 1.0/6.0;
    lambda[9] = 1.0/12.0;
    global_bias = 10.0;
    return;
}

void calc_benchmark_norm()
{
    int i;
    transform_norm (0);    norm_f[0] = calc_ackley(trans_x);
    transform_norm (1);    norm_f[1] = calc_ackley(trans_x);
    transform_norm (2);    norm_f[2] = calc_rastrigin(trans_x);
    transform_norm (3);    norm_f[3] = calc_rastrigin(trans_x);
    transform_norm (4);    norm_f[4] = calc_sphere(trans_x);
    transform_norm (5);    norm_f[5] = calc_sphere(trans_x);
    for (i=0; i<nreal; i++)
    {
        norm_x[i] = 0.0;
    }
    transform_norm (6);    norm_f[6] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform_norm (7);    norm_f[7] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform_norm (8);    norm_f[8] = calc_griewank(trans_x);
    transform_norm (9);    norm_f[9] = calc_griewank(trans_x);
    return;
}

long double calc_benchmark_func(long double *x)
{
    int i;
    long double res;
    transform (x, 0);    basic_f[0] = calc_ackley(trans_x);
    transform (x, 1);    basic_f[1] = calc_ackley(trans_x);
    transform (x, 2);    basic_f[2] = calc_rastrigin(trans_x);
    transform (x, 3);    basic_f[3] = calc_rastrigin(trans_x);
    transform (x, 4);    basic_f[4] = calc_sphere(trans_x);
    transform (x, 5);    basic_f[5] = calc_sphere(trans_x);
    for (i=0; i<nreal; i++)
    {
        norm_x[i] = 0.0;
    }
    transform (x, 6);    basic_f[6] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform (x, 7);    basic_f[7] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform (x, 8);    basic_f[8] = calc_griewank(trans_x);
    transform (x, 9);    basic_f[9] = calc_griewank(trans_x);
    for (i=0; i<nfunc; i++)
    {
        basic_f[i] *= C/norm_f[i];
    }
    calc_weight(x);
    res = global_bias;
    for (i=0; i<nfunc; i++)
    {
        res += weight[i]*(basic_f[i]+bias[i]);
    }
    return (res);
}
#endif

# ifdef f19
void initialize()
{
    int i, j, k;
    fstream input;
    input.open("input_data/hybrid_func2_data.txt", ios::in);
    if(!input)
    {
        cerr << "File could not be opened for reading" << endl;
        exit(1);
    }
    for (i=0; i<nfunc; i++)
    {
        for (j=0; j<nreal; j++)
        {
            input >> o[i][j];
            cout << scientific << "\nO[" << i+1 << "][" << j+1 << "] = " << o[i][j];
        }
    }
    cout << endl << endl;
    input.close();

    if (nreal==2)    input.open("input_data/hybrid_func2_M_D2.txt", ios::in);
    if (nreal==10)   input.open("input_data/hybrid_func2_M_D10.txt", ios::in);
    if (nreal==30)   input.open("input_data/hybrid_func2_M_D30.txt", ios::in);
    if (nreal==50)   input.open("input_data/hybrid_func2_M_D50.txt", ios::in);
    if(!input)
    {
        cerr << "File could not be opened for reading" << endl;
        exit(1);
    }
    for (i=0; i<nfunc; i++)
    {
        for (j=0; j<nreal; j++)
        {
            for (k=0; k<nreal; k++)
            {
                input >> l[i][j][k];
                cout << scientific << " M[" << i+1 << "][" << j+1 << "][" << k+1 << "] = " << l[i][j][k] << endl;
            }
        }
    }
    input.close();

    for (i=0; i<nreal; i++)
    {
        o[nfunc-1][i] = 0.0;
    }
    sigma[0] = 0.1;
    sigma[1] = 2.0;
    sigma[2] = 1.5;
    sigma[3] = 1.5;
    sigma[4] = 1.0;
    sigma[5] = 1.0;
    sigma[6] = 1.5;
    sigma[7] = 1.5;
    sigma[8] = 2.0;
    sigma[9] = 2.0;
    lambda[0] = 0.5/32.0;
    lambda[1] = 5.0/32.0;
    lambda[2] = 2.0;
    lambda[3] = 1.0;
    lambda[4] = 1.0/10.0;
    lambda[5] = 1.0/20.0;
    lambda[6] = 20.0;
    lambda[7] = 10.0;
    lambda[8] = 1.0/6.0;
    lambda[9] = 1.0/12.0;
    global_bias = 10.0;
    return;
}

void calc_benchmark_norm()
{
    int i;
    transform_norm (0);    norm_f[0] = calc_ackley(trans_x);
    transform_norm (1);    norm_f[1] = calc_ackley(trans_x);
    transform_norm (2);    norm_f[2] = calc_rastrigin(trans_x);
    transform_norm (3);    norm_f[3] = calc_rastrigin(trans_x);
    transform_norm (4);    norm_f[4] = calc_sphere(trans_x);
    transform_norm (5);    norm_f[5] = calc_sphere(trans_x);
    for (i=0; i<nreal; i++)
    {
        norm_x[i] = 0.0;
    }
    transform_norm (6);    norm_f[6] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform_norm (7);    norm_f[7] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform_norm (8);    norm_f[8] = calc_griewank(trans_x);
    transform_norm (9);    norm_f[9] = calc_griewank(trans_x);
    return;
}

long double calc_benchmark_func(long double *x)
{
    int i;
    long double res;
    transform (x, 0);    basic_f[0] = calc_ackley(trans_x);
    transform (x, 1);    basic_f[1] = calc_ackley(trans_x);
    transform (x, 2);    basic_f[2] = calc_rastrigin(trans_x);
    transform (x, 3);    basic_f[3] = calc_rastrigin(trans_x);
    transform (x, 4);    basic_f[4] = calc_sphere(trans_x);
    transform (x, 5);    basic_f[5] = calc_sphere(trans_x);
    for (i=0; i<nreal; i++)
    {
        norm_x[i] = 0.0;
    }
    transform (x, 6);    basic_f[6] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform (x, 7);    basic_f[7] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform (x, 8);    basic_f[8] = calc_griewank(trans_x);
    transform (x, 9);    basic_f[9] = calc_griewank(trans_x);
    for (i=0; i<nfunc; i++)
    {
        basic_f[i] *= C/norm_f[i];
    }
    calc_weight(x);
    res = global_bias;
    for (i=0; i<nfunc; i++)
    {
        res += weight[i]*(basic_f[i]+bias[i]);
    }
    return (res);
}
#endif

# ifdef f20
void initialize()
{
    int i, j, k;
    int index;
    fstream input;
    input.open("input_data/hybrid_func2_data.txt", ios::in);
    if(!input)
    {
        cerr << "File could not be opened for reading" << endl;
        exit(1);
    }
    for (i=0; i<nfunc; i++)
    {
        for (j=0; j<nreal; j++)
        {
            input >> o[i][j];
            cout << scientific << "\nO[" << i+1 << "][" << j+1 << "] = " << o[i][j];
        }
    }
    cout << endl << endl;
    input.close();

    index = nreal/2;
    for (i=1; i<=index; i++)
    {
        o[0][2*i-1] = 5.0;
    }

    if (nreal==2)    input.open("input_data/hybrid_func2_M_D2.txt", ios::in);
    if (nreal==10)   input.open("input_data/hybrid_func2_M_D10.txt", ios::in);
    if (nreal==30)   input.open("input_data/hybrid_func2_M_D30.txt", ios::in);
    if (nreal==50)   input.open("input_data/hybrid_func2_M_D50.txt", ios::in);
    if(!input)
    {
        cerr << "File could not be opened for reading" << endl;
        exit(1);
    }
    for (i=0; i<nfunc; i++)
    {
        for (j=0; j<nreal; j++)
        {
            for (k=0; k<nreal; k++)
            {
                input >> l[i][j][k];
                cout << scientific << " M[" << i+1 << "][" << j+1 << "][" << k+1 << "] = " << l[i][j][k] << endl;
            }
        }
    }
    input.close();

    for (i=0; i<nreal; i++)
    {
        o[nfunc-1][i] = 0.0;
    }
    sigma[0] = 1.0;
    sigma[1] = 2.0;
    sigma[2] = 1.5;
    sigma[3] = 1.5;
    sigma[4] = 1.0;
    sigma[5] = 1.0;
    sigma[6] = 1.5;
    sigma[7] = 1.5;
    sigma[8] = 2.0;
    sigma[9] = 2.0;
    lambda[0] = 5.0/16.0;
    lambda[1] = 5.0/32.0;
    lambda[2] = 2.0;
    lambda[3] = 1.0;
    lambda[4] = 1.0/10.0;
    lambda[5] = 1.0/20.0;
    lambda[6] = 20.0;
    lambda[7] = 10.0;
    lambda[8] = 1.0/6.0;
    lambda[9] = 1.0/12.0;
    global_bias = 10.0;
    return;
}

void calc_benchmark_norm()
{
    int i;
    transform_norm (0);    norm_f[0] = calc_ackley(trans_x);
    transform_norm (1);    norm_f[1] = calc_ackley(trans_x);
    transform_norm (2);    norm_f[2] = calc_rastrigin(trans_x);
    transform_norm (3);    norm_f[3] = calc_rastrigin(trans_x);
    transform_norm (4);    norm_f[4] = calc_sphere(trans_x);
    transform_norm (5);    norm_f[5] = calc_sphere(trans_x);
    for (i=0; i<nreal; i++)
    {
        norm_x[i] = 0.0;
    }
    transform_norm (6);    norm_f[6] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform_norm (7);    norm_f[7] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform_norm (8);    norm_f[8] = calc_griewank(trans_x);
    transform_norm (9);    norm_f[9] = calc_griewank(trans_x);
    return;
}

long double calc_benchmark_func(long double *x)
{
    int i;
    long double res;
    transform (x, 0);    basic_f[0] = calc_ackley(trans_x);
    transform (x, 1);    basic_f[1] = calc_ackley(trans_x);
    transform (x, 2);    basic_f[2] = calc_rastrigin(trans_x);
    transform (x, 3);    basic_f[3] = calc_rastrigin(trans_x);
    transform (x, 4);    basic_f[4] = calc_sphere(trans_x);
    transform (x, 5);    basic_f[5] = calc_sphere(trans_x);
    for (i=0; i<nreal; i++)
    {
        norm_x[i] = 0.0;
    }
    transform (x, 6);    basic_f[6] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform (x, 7);    basic_f[7] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform (x, 8);    basic_f[8] = calc_griewank(trans_x);
    transform (x, 9);    basic_f[9] = calc_griewank(trans_x);
    for (i=0; i<nfunc; i++)
    {
        basic_f[i] *= C/norm_f[i];
    }
    calc_weight(x);
    res = global_bias;
    for (i=0; i<nfunc; i++)
    {
        res += weight[i]*(basic_f[i]+bias[i]);
    }
    return (res);
}
#endif

# ifdef f21
void initialize()
{
    int i, j, k;
    fstream input;
    input.open("input_data/hybrid_func3_data.txt", ios::in);
    if(!input)
    {
        cerr << "File could not be opened for reading" << endl;
        exit(1);
    }
    for (i=0; i<nfunc; i++)
    {
        for (j=0; j<nreal; j++)
        {
            input >> o[i][j];
            cout << scientific << "\nO[" << i+1 << "][" << j+1 << "] = " << o[i][j];
        }
    }
    cout << endl << endl;
    input.close();

    if (nreal==2)    input.open("input_data/hybrid_func3_M_D2.txt", ios::in);
    if (nreal==10)   input.open("input_data/hybrid_func3_M_D10.txt", ios::in);
    if (nreal==30)   input.open("input_data/hybrid_func3_M_D30.txt", ios::in);
    if (nreal==50)   input.open("input_data/hybrid_func3_M_D50.txt", ios::in);
    if(!input)
    {
        cerr << "File could not be opened for reading" << endl;
        exit(1);
    }
    for (i=0; i<nfunc; i++)
    {
        for (j=0; j<nreal; j++)
        {
            for (k=0; k<nreal; k++)
            {
                input >> l[i][j][k];
                cout << scientific << " M[" << i+1 << "][" << j+1 << "][" << k+1 << "] = " << l[i][j][k] << endl;
            }
        }
    }
    input.close();

    sigma[0] = 1.0;
    sigma[1] = 1.0;
    sigma[2] = 1.0;
    sigma[3] = 1.0;
    sigma[4] = 1.0;
    sigma[5] = 2.0;
    sigma[6] = 2.0;
    sigma[7] = 2.0;
    sigma[8] = 2.0;
    sigma[9] = 2.0;
    lambda[0] = 1.0/4.0;
    lambda[1] = 1.0/20.0;
    lambda[2] = 5.0;
    lambda[3] = 1.0;
    lambda[4] = 5.0;
    lambda[5] = 1.0;
    lambda[6] = 50.0;
    lambda[7] = 10.0;
    lambda[8] = 1.0/8.0;
    lambda[9] = 1.0/40.0;
    global_bias = 360.0;
    return;
}

void calc_benchmark_norm()
{
    int i;
    long double temp1, temp2, temp;
    transform_norm (0);
    norm_f[0] = 0.0;
    for (i=0; i<nreal-1; i++)
    {
        temp1 = pow((sin(sqrt(pow(trans_x[i],2.0)+pow(trans_x[i+1],2.0)))),2.0);
        temp2 = 1.0 + 0.001*(pow(trans_x[i],2.0)+pow(trans_x[i+1],2.0));
        norm_f[0] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
    }
    temp1 = pow((sin(sqrt(pow(trans_x[nreal-1],2.0)+pow(trans_x[0],2.0)))),2.0);
    temp2 = 1.0 + 0.001*(pow(trans_x[nreal-1],2.0)+pow(trans_x[0],2.0));
    norm_f[0] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
    transform_norm (1);
    norm_f[1] = 0.0;
    for (i=0; i<nreal-1; i++)
    {
        temp1 = pow((sin(sqrt(pow(trans_x[i],2.0)+pow(trans_x[i+1],2.0)))),2.0);
        temp2 = 1.0 + 0.001*(pow(trans_x[i],2.0)+pow(trans_x[i+1],2.0));
        norm_f[1] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
    }
    temp1 = pow((sin(sqrt(pow(trans_x[nreal-1],2.0)+pow(trans_x[0],2.0)))),2.0);
    temp2 = 1.0 + 0.001*(pow(trans_x[nreal-1],2.0)+pow(trans_x[0],2.0));
    norm_f[1] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
    transform_norm (2);    norm_f[2] = calc_rastrigin(trans_x);
    transform_norm (3);    norm_f[3] = calc_rastrigin(trans_x);
    transform_norm(4);
    norm_f[4] = 0.0;
    for (i=0; i<nreal-1; i++)
    {
        temp = 100.0*pow((trans_x[i]*trans_x[i]-trans_x[i+1]),2.0) + 1.0*pow((trans_x[i]-1.0),2.0);
        norm_f[4] += (temp*temp)/4000.0 - cos(temp) + 1.0;
    }
    temp = 100.0*pow((trans_x[nreal-1]*trans_x[nreal-1]-trans_x[0]),2.0) + 1.0*pow((trans_x[nreal-1]-1.0),2.0);
    norm_f[4] += (temp*temp)/4000.0 - cos(temp) + 1.0;
    transform_norm(5);
    norm_f[5] = 0.0;
    for (i=0; i<nreal-1; i++)
    {
        temp = 100.0*pow((trans_x[i]*trans_x[i]-trans_x[i+1]),2.0) + 1.0*pow((trans_x[i]-1.0),2.0);
        norm_f[5] += (temp*temp)/4000.0 - cos(temp) + 1.0;
    }
    temp = 100.0*pow((trans_x[nreal-1]*trans_x[nreal-1]-trans_x[0]),2.0) + 1.0*pow((trans_x[nreal-1]-1.0),2.0);
    norm_f[5] += (temp*temp)/4000.0 - cos(temp) + 1.0;
    for (i=0; i<nreal; i++)
    {
        norm_x[i] = 0.0;
    }
    transform_norm (6);    norm_f[6] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform_norm (7);    norm_f[7] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform_norm (8);    norm_f[8] = calc_griewank(trans_x);
    transform_norm (9);    norm_f[9] = calc_griewank(trans_x);
    return;
}

long double calc_benchmark_func(long double *x)
{
    int i;
    long double temp1, temp2, temp;
    long double res;
    transform (x, 0);
    basic_f[0] = 0.0;
    for (i=0; i<nreal-1; i++)
    {
        temp1 = pow((sin(sqrt(pow(trans_x[i],2.0)+pow(trans_x[i+1],2.0)))),2.0);
        temp2 = 1.0 + 0.001*(pow(trans_x[i],2.0)+pow(trans_x[i+1],2.0));
        basic_f[0] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
    }
    temp1 = pow((sin(sqrt(pow(trans_x[nreal-1],2.0)+pow(trans_x[0],2.0)))),2.0);
    temp2 = 1.0 + 0.001*(pow(trans_x[nreal-1],2.0)+pow(trans_x[0],2.0));
    basic_f[0] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
    transform (x, 1);
    basic_f[1] = 0.0;
    for (i=0; i<nreal-1; i++)
    {
        temp1 = pow((sin(sqrt(pow(trans_x[i],2.0)+pow(trans_x[i+1],2.0)))),2.0);
        temp2 = 1.0 + 0.001*(pow(trans_x[i],2.0)+pow(trans_x[i+1],2.0));
        basic_f[1] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
    }
    temp1 = pow((sin(sqrt(pow(trans_x[nreal-1],2.0)+pow(trans_x[0],2.0)))),2.0);
    temp2 = 1.0 + 0.001*(pow(trans_x[nreal-1],2.0)+pow(trans_x[0],2.0));
    basic_f[1] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
    transform (x, 2);    basic_f[2] = calc_rastrigin(trans_x);
    transform (x, 3);    basic_f[3] = calc_rastrigin(trans_x);
    transform (x, 4);
    basic_f[4] = 0.0;
    for (i=0; i<nreal-1; i++)
    {
        temp = 100.0*pow((trans_x[i]*trans_x[i]-trans_x[i+1]),2.0) + 1.0*pow((trans_x[i]-1.0),2.0);
        basic_f[4] += (temp*temp)/4000.0 - cos(temp) + 1.0;
    }
    temp = 100.0*pow((trans_x[nreal-1]*trans_x[nreal-1]-trans_x[0]),2.0) + 1.0*pow((trans_x[nreal-1]-1.0),2.0);
    basic_f[4] += (temp*temp)/4000.0 - cos(temp) + 1.0;
    transform(x, 5);
    basic_f[5] = 0.0;
    for (i=0; i<nreal-1; i++)
    {
        temp = 100.0*pow((trans_x[i]*trans_x[i]-trans_x[i+1]),2.0) + 1.0*pow((trans_x[i]-1.0),2.0);
        basic_f[5] += (temp*temp)/4000.0 - cos(temp) + 1.0;
    }
    temp = 100.0*pow((trans_x[nreal-1]*trans_x[nreal-1]-trans_x[0]),2.0) + 1.0*pow((trans_x[nreal-1]-1.0),2.0);
    basic_f[5] += (temp*temp)/4000.0 - cos(temp) + 1.0;
    for (i=0; i<nreal; i++)
    {
        norm_x[i] = 0.0;
    }
    transform (x, 6);    basic_f[6] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform (x, 7);    basic_f[7] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform (x, 8);    basic_f[8] = calc_griewank(trans_x);
    transform (x, 9);    basic_f[9] = calc_griewank(trans_x);
    for (i=0; i<nfunc; i++)
    {
        basic_f[i] *= C/norm_f[i];
    }
    calc_weight(x);
    res = global_bias;
    for (i=0; i<nfunc; i++)
    {
        res += weight[i]*(basic_f[i]+bias[i]);
    }
    return (res);
}
#endif

# ifdef f22
void initialize()
{
    int i, j, k;
    fstream input;
    input.open("input_data/hybrid_func3_data.txt", ios::in);
    if(!input)
    {
        cerr << "File could not be opened" << endl;
        exit(1);
    }
    for (i=0; i<nfunc; i++)
    {
        for (j=0; j<nreal; j++)
        {
            input >> o[i][j];
            cout << scientific << "\nO[" << i+1 << "][" << j+1 << "] = " << o[i][j];
        }
    }
    cout << endl << endl;
    input.close();

    if (nreal==2)    input.open("input_data/hybrid_func3_HM_D2.txt", ios::in);
    if (nreal==10)   input.open("input_data/hybrid_func3_HM_D10.txt", ios::in);
    if (nreal==30)   input.open("input_data/hybrid_func3_HM_D30.txt", ios::in);
    if (nreal==50)   input.open("input_data/hybrid_func3_HM_D50.txt", ios::in);
    if(!input)
    {
        cerr << "File could not be opened for reading" << endl;
        exit(1);
    }
    for (i=0; i<nfunc; i++)
    {
        for (j=0; j<nreal; j++)
        {
            for (k=0; k<nreal; k++)
            {
                input >> l[i][j][k];
                cout << scientific << " M[" << i+1 << "][" << j+1 << "][" << k+1 << "] = " << l[i][j][k] << endl;
            }
        }
    }
    input.close();

    sigma[0] = 1.0;
    sigma[1] = 1.0;
    sigma[2] = 1.0;
    sigma[3] = 1.0;
    sigma[4] = 1.0;
    sigma[5] = 2.0;
    sigma[6] = 2.0;
    sigma[7] = 2.0;
    sigma[8] = 2.0;
    sigma[9] = 2.0;
    lambda[0] = 1.0/4.0;
    lambda[1] = 1.0/20.0;
    lambda[2] = 5.0;
    lambda[3] = 1.0;
    lambda[4] = 5.0;
    lambda[5] = 1.0;
    lambda[6] = 50.0;
    lambda[7] = 10.0;
    lambda[8] = 1.0/8.0;
    lambda[9] = 1.0/40.0;
    global_bias = 360.0;
    return;
}

void calc_benchmark_norm()
{
    int i;
    long double temp1, temp2, temp;
    transform_norm (0);
    norm_f[0] = 0.0;
    for (i=0; i<nreal-1; i++)
    {
        temp1 = pow((sin(sqrt(pow(trans_x[i],2.0)+pow(trans_x[i+1],2.0)))),2.0);
        temp2 = 1.0 + 0.001*(pow(trans_x[i],2.0)+pow(trans_x[i+1],2.0));
        norm_f[0] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
    }
    temp1 = pow((sin(sqrt(pow(trans_x[nreal-1],2.0)+pow(trans_x[0],2.0)))),2.0);
    temp2 = 1.0 + 0.001*(pow(trans_x[nreal-1],2.0)+pow(trans_x[0],2.0));
    norm_f[0] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
    transform_norm (1);
    norm_f[1] = 0.0;
    for (i=0; i<nreal-1; i++)
    {
        temp1 = pow((sin(sqrt(pow(trans_x[i],2.0)+pow(trans_x[i+1],2.0)))),2.0);
        temp2 = 1.0 + 0.001*(pow(trans_x[i],2.0)+pow(trans_x[i+1],2.0));
        norm_f[1] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
    }
    temp1 = pow((sin(sqrt(pow(trans_x[nreal-1],2.0)+pow(trans_x[0],2.0)))),2.0);
    temp2 = 1.0 + 0.001*(pow(trans_x[nreal-1],2.0)+pow(trans_x[0],2.0));
    norm_f[1] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
    transform_norm (2);    norm_f[2] = calc_rastrigin(trans_x);
    transform_norm (3);    norm_f[3] = calc_rastrigin(trans_x);
    transform_norm(4);
    norm_f[4] = 0.0;
    for (i=0; i<nreal-1; i++)
    {
        temp = 100.0*pow((trans_x[i]*trans_x[i]-trans_x[i+1]),2.0) + 1.0*pow((trans_x[i]-1.0),2.0);
        norm_f[4] += (temp*temp)/4000.0 - cos(temp) + 1.0;
    }
    temp = 100.0*pow((trans_x[nreal-1]*trans_x[nreal-1]-trans_x[0]),2.0) + 1.0*pow((trans_x[nreal-1]-1.0),2.0);
    norm_f[4] += (temp*temp)/4000.0 - cos(temp) + 1.0;
    transform_norm(5);
    norm_f[5] = 0.0;
    for (i=0; i<nreal-1; i++)
    {
        temp = 100.0*pow((trans_x[i]*trans_x[i]-trans_x[i+1]),2.0) + 1.0*pow((trans_x[i]-1.0),2.0);
        norm_f[5] += (temp*temp)/4000.0 - cos(temp) + 1.0;
    }
    temp = 100.0*pow((trans_x[nreal-1]*trans_x[nreal-1]-trans_x[0]),2.0) + 1.0*pow((trans_x[nreal-1]-1.0),2.0);
    norm_f[5] += (temp*temp)/4000.0 - cos(temp) + 1.0;
    for (i=0; i<nreal; i++)
    {
        norm_x[i] = 0.0;
    }
    transform_norm (6);    norm_f[6] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform_norm (7);    norm_f[7] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform_norm (8);    norm_f[8] = calc_griewank(trans_x);
    transform_norm (9);    norm_f[9] = calc_griewank(trans_x);
    return;
}

long double calc_benchmark_func(long double *x)
{
    int i;
    long double temp1, temp2, temp;
    long double res;
    transform (x, 0);
    basic_f[0] = 0.0;
    for (i=0; i<nreal-1; i++)
    {
        temp1 = pow((sin(sqrt(pow(trans_x[i],2.0)+pow(trans_x[i+1],2.0)))),2.0);
        temp2 = 1.0 + 0.001*(pow(trans_x[i],2.0)+pow(trans_x[i+1],2.0));
        basic_f[0] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
    }
    temp1 = pow((sin(sqrt(pow(trans_x[nreal-1],2.0)+pow(trans_x[0],2.0)))),2.0);
    temp2 = 1.0 + 0.001*(pow(trans_x[nreal-1],2.0)+pow(trans_x[0],2.0));
    basic_f[0] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
    transform (x, 1);
    basic_f[1] = 0.0;
    for (i=0; i<nreal-1; i++)
    {
        temp1 = pow((sin(sqrt(pow(trans_x[i],2.0)+pow(trans_x[i+1],2.0)))),2.0);
        temp2 = 1.0 + 0.001*(pow(trans_x[i],2.0)+pow(trans_x[i+1],2.0));
        basic_f[1] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
    }
    temp1 = pow((sin(sqrt(pow(trans_x[nreal-1],2.0)+pow(trans_x[0],2.0)))),2.0);
    temp2 = 1.0 + 0.001*(pow(trans_x[nreal-1],2.0)+pow(trans_x[0],2.0));
    basic_f[1] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
    transform (x, 2);    basic_f[2] = calc_rastrigin(trans_x);
    transform (x, 3);    basic_f[3] = calc_rastrigin(trans_x);
    transform (x, 4);
    basic_f[4] = 0.0;
    for (i=0; i<nreal-1; i++)
    {
        temp = 100.0*pow((trans_x[i]*trans_x[i]-trans_x[i+1]),2.0) + 1.0*pow((trans_x[i]-1.0),2.0);
        basic_f[4] += (temp*temp)/4000.0 - cos(temp) + 1.0;
    }
    temp = 100.0*pow((trans_x[nreal-1]*trans_x[nreal-1]-trans_x[0]),2.0) + 1.0*pow((trans_x[nreal-1]-1.0),2.0);
    basic_f[4] += (temp*temp)/4000.0 - cos(temp) + 1.0;
    transform(x, 5);
    basic_f[5] = 0.0;
    for (i=0; i<nreal-1; i++)
    {
        temp = 100.0*pow((trans_x[i]*trans_x[i]-trans_x[i+1]),2.0) + 1.0*pow((trans_x[i]-1.0),2.0);
        basic_f[5] += (temp*temp)/4000.0 - cos(temp) + 1.0;
    }
    temp = 100.0*pow((trans_x[nreal-1]*trans_x[nreal-1]-trans_x[0]),2.0) + 1.0*pow((trans_x[nreal-1]-1.0),2.0);
    basic_f[5] += (temp*temp)/4000.0 - cos(temp) + 1.0;
    for (i=0; i<nreal; i++)
    {
        norm_x[i] = 0.0;
    }
    transform (x, 6);    basic_f[6] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform (x, 7);    basic_f[7] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform (x, 8);    basic_f[8] = calc_griewank(trans_x);
    transform (x, 9);    basic_f[9] = calc_griewank(trans_x);
    for (i=0; i<nfunc; i++)
    {
        basic_f[i] *= C/norm_f[i];
    }
    calc_weight(x);
    res = global_bias;
    for (i=0; i<nfunc; i++)
    {
        res += weight[i]*(basic_f[i]+bias[i]);
    }
    return (res);
}
#endif

# ifdef f23
void initialize()
{
    int i, j, k;
    fstream input;
    input.open("input_data/hybrid_func3_data.txt", ios::in);
    if(!input)
    {
        cerr << "File could not be opened for reading" << endl;
        exit(1);
    }
    for (i=0; i<nfunc; i++)
    {
        for (j=0; j<nreal; j++)
        {
            input >> o[i][j];
            cout << scientific << "\nO[" << i+1 << "][" << j+1 << "] = " << o[i][j];
        }
    }
    cout << endl << endl;
    input.close();

    if (nreal==2)    input.open("input_data/hybrid_func3_M_D2.txt", ios::in);
    if (nreal==10)   input.open("input_data/hybrid_func3_M_D10.txt", ios::in);
    if (nreal==30)   input.open("input_data/hybrid_func3_M_D30.txt", ios::in);
    if (nreal==50)   input.open("input_data/hybrid_func3_M_D50.txt", ios::in);
    if(!input)
    {
        cerr << "File could not be opened for reading" << endl;
        exit(1);
    }
    for (i=0; i<nfunc; i++)
    {
        for (j=0; j<nreal; j++)
        {
            for (k=0; k<nreal; k++)
            {
                input >> l[i][j][k];
                cout << scientific << " M[" << i+1 << "][" << j+1 << "][" << k+1 << "] = " << l[i][j][k] << endl;
            }
        }
    }
    input.close();

    sigma[0] = 1.0;
    sigma[1] = 1.0;
    sigma[2] = 1.0;
    sigma[3] = 1.0;
    sigma[4] = 1.0;
    sigma[5] = 2.0;
    sigma[6] = 2.0;
    sigma[7] = 2.0;
    sigma[8] = 2.0;
    sigma[9] = 2.0;
    lambda[0] = 1.0/4.0;
    lambda[1] = 1.0/20.0;
    lambda[2] = 5.0;
    lambda[3] = 1.0;
    lambda[4] = 5.0;
    lambda[5] = 1.0;
    lambda[6] = 50.0;
    lambda[7] = 10.0;
    lambda[8] = 1.0/8.0;
    lambda[9] = 1.0/40.0;
    global_bias = 360.0;
    return;
}

void calc_benchmark_norm()
{
    int i;
    long double temp1, temp2, temp;
    transform_norm (0);
    norm_f[0] = 0.0;
    for (i=0; i<nreal-1; i++)
    {
        temp1 = pow((sin(sqrt(pow(trans_x[i],2.0)+pow(trans_x[i+1],2.0)))),2.0);
        temp2 = 1.0 + 0.001*(pow(trans_x[i],2.0)+pow(trans_x[i+1],2.0));
        norm_f[0] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
    }
    temp1 = pow((sin(sqrt(pow(trans_x[nreal-1],2.0)+pow(trans_x[0],2.0)))),2.0);
    temp2 = 1.0 + 0.001*(pow(trans_x[nreal-1],2.0)+pow(trans_x[0],2.0));
    norm_f[0] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
    transform_norm (1);
    norm_f[1] = 0.0;
    for (i=0; i<nreal-1; i++)
    {
        temp1 = pow((sin(sqrt(pow(trans_x[i],2.0)+pow(trans_x[i+1],2.0)))),2.0);
        temp2 = 1.0 + 0.001*(pow(trans_x[i],2.0)+pow(trans_x[i+1],2.0));
        norm_f[1] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
    }
    temp1 = pow((sin(sqrt(pow(trans_x[nreal-1],2.0)+pow(trans_x[0],2.0)))),2.0);
    temp2 = 1.0 + 0.001*(pow(trans_x[nreal-1],2.0)+pow(trans_x[0],2.0));
    norm_f[1] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
    transform_norm (2);    norm_f[2] = calc_rastrigin(trans_x);
    transform_norm (3);    norm_f[3] = calc_rastrigin(trans_x);
    transform_norm(4);
    norm_f[4] = 0.0;
    for (i=0; i<nreal-1; i++)
    {
        temp = 100.0*pow((trans_x[i]*trans_x[i]-trans_x[i+1]),2.0) + 1.0*pow((trans_x[i]-1.0),2.0);
        norm_f[4] += (temp*temp)/4000.0 - cos(temp) + 1.0;
    }
    temp = 100.0*pow((trans_x[nreal-1]*trans_x[nreal-1]-trans_x[0]),2.0) + 1.0*pow((trans_x[nreal-1]-1.0),2.0);
    norm_f[4] += (temp*temp)/4000.0 - cos(temp) + 1.0;
    transform_norm(5);
    norm_f[5] = 0.0;
    for (i=0; i<nreal-1; i++)
    {
        temp = 100.0*pow((trans_x[i]*trans_x[i]-trans_x[i+1]),2.0) + 1.0*pow((trans_x[i]-1.0),2.0);
        norm_f[5] += (temp*temp)/4000.0 - cos(temp) + 1.0;
    }
    temp = 100.0*pow((trans_x[nreal-1]*trans_x[nreal-1]-trans_x[0]),2.0) + 1.0*pow((trans_x[nreal-1]-1.0),2.0);
    norm_f[5] += (temp*temp)/4000.0 - cos(temp) + 1.0;
    for (i=0; i<nreal; i++)
    {
        norm_x[i] = 0.0;
    }
    transform_norm (6);    norm_f[6] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform_norm (7);    norm_f[7] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform_norm (8);    norm_f[8] = calc_griewank(trans_x);
    transform_norm (9);    norm_f[9] = calc_griewank(trans_x);
    return;
}

long double calc_benchmark_func(long double *x)
{
    int i;
    long double temp1, temp2, temp;
    long double res;
    int a;
    long double b;
    for (i=0; i<nreal; i++)
    {
        if (fabs(x[i]-o[0][i]) >= 0.5)
        {
            res = 2.0*x[i];
            a = res;
            b = fabs(res-a);
            if (b<0.5)
            {
                temp_x4[i] = a/2.0;
            }
            else
            {
                if (res<=0.0)
                {
                    temp_x4[i] = (a-1.0)/2.0;
                }
                else
                {
                    temp_x4[i] = (a+1.0)/2.0;
                }
            }
        }
        else
        {
            temp_x4[i] = x[i];
        }
    }
    transform (temp_x4, 0);
    basic_f[0] = 0.0;
    for (i=0; i<nreal-1; i++)
    {
        temp1 = pow((sin(sqrt(pow(trans_x[i],2.0)+pow(trans_x[i+1],2.0)))),2.0);
        temp2 = 1.0 + 0.001*(pow(trans_x[i],2.0)+pow(trans_x[i+1],2.0));
        basic_f[0] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
    }
    temp1 = pow((sin(sqrt(pow(trans_x[nreal-1],2.0)+pow(trans_x[0],2.0)))),2.0);
    temp2 = 1.0 + 0.001*(pow(trans_x[nreal-1],2.0)+pow(trans_x[0],2.0));
    basic_f[0] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
    transform (temp_x4, 1);
    basic_f[1] = 0.0;
    for (i=0; i<nreal-1; i++)
    {
        temp1 = pow((sin(sqrt(pow(trans_x[i],2.0)+pow(trans_x[i+1],2.0)))),2.0);
        temp2 = 1.0 + 0.001*(pow(trans_x[i],2.0)+pow(trans_x[i+1],2.0));
        basic_f[1] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
    }
    temp1 = pow((sin(sqrt(pow(trans_x[nreal-1],2.0)+pow(trans_x[0],2.0)))),2.0);
    temp2 = 1.0 + 0.001*(pow(trans_x[nreal-1],2.0)+pow(trans_x[0],2.0));
    basic_f[1] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
    transform (temp_x4, 2);    basic_f[2] = calc_rastrigin(trans_x);
    transform (temp_x4, 3);    basic_f[3] = calc_rastrigin(trans_x);
    transform (temp_x4, 4);
    basic_f[4] = 0.0;
    for (i=0; i<nreal-1; i++)
    {
        temp = 100.0*pow((trans_x[i]*trans_x[i]-trans_x[i+1]),2.0) + 1.0*pow((trans_x[i]-1.0),2.0);
        basic_f[4] += (temp*temp)/4000.0 - cos(temp) + 1.0;
    }
    temp = 100.0*pow((trans_x[nreal-1]*trans_x[nreal-1]-trans_x[0]),2.0) + 1.0*pow((trans_x[nreal-1]-1.0),2.0);
    basic_f[4] += (temp*temp)/4000.0 - cos(temp) + 1.0;
    transform(temp_x4, 5);
    basic_f[5] = 0.0;
    for (i=0; i<nreal-1; i++)
    {
        temp = 100.0*pow((trans_x[i]*trans_x[i]-trans_x[i+1]),2.0) + 1.0*pow((trans_x[i]-1.0),2.0);
        basic_f[5] += (temp*temp)/4000.0 - cos(temp) + 1.0;
    }
    temp = 100.0*pow((trans_x[nreal-1]*trans_x[nreal-1]-trans_x[0]),2.0) + 1.0*pow((trans_x[nreal-1]-1.0),2.0);
    basic_f[5] += (temp*temp)/4000.0 - cos(temp) + 1.0;
    for (i=0; i<nreal; i++)
    {
        norm_x[i] = 0.0;
    }
    transform (temp_x4, 6);    basic_f[6] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform (temp_x4, 7);    basic_f[7] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform (temp_x4, 8);    basic_f[8] = calc_griewank(trans_x);
    transform (temp_x4, 9);    basic_f[9] = calc_griewank(trans_x);
    for (i=0; i<nfunc; i++)
    {
        basic_f[i] *= C/norm_f[i];
    }
    calc_weight(temp_x4);
    res = global_bias;
    for (i=0; i<nfunc; i++)
    {
        res += weight[i]*(basic_f[i]+bias[i]);
    }
    return (res);
}
#endif

# if defined f24 || defined f25
void initialize()
{
    int i, j, k;
    fstream input;
    input.open("input_data/hybrid_func4_data.txt", ios::in);
    if(!input)
    {
        cerr << "File could not be opened for reading" << endl;
        exit(1);
    }
    for (i=0; i<nfunc; i++)
    {
        for (j=0; j<nreal; j++)
        {
            input >> o[i][j];
            cout << scientific << "\nO[" << i+1 << "][" << j+1 << "] = " << o[i][j];
        }
    }
    cout << endl << endl;
    input.close();

    if (nreal==2)    input.open("input_data/hybrid_func4_M_D2.txt", ios::in);
    if (nreal==10)   input.open("input_data/hybrid_func4_M_D10.txt", ios::in);
    if (nreal==30)   input.open("input_data/hybrid_func4_M_D30.txt", ios::in);
    if (nreal==50)   input.open("input_data/hybrid_func4_M_D50.txt", ios::in);
    if(!input)
    {
        cerr << "File could not be opened for reading" << endl;
        exit(1);
    }
    for (i=0; i<nfunc; i++)
    {
        for (j=0; j<nreal; j++)
        {
            for (k=0; k<nreal; k++)
            {
                input >> l[i][j][k];
                cout << scientific << " M[" << i+1 << "][" << j+1 << "][" << k+1 << "] = " << l[i][j][k] << endl;
            }
        }
    }
    input.close();

    for (i=0; i<nfunc; i++)
    {
        sigma[i] = 2.0;
    }
    lambda[0] = 10.0;
    lambda[1] = 1.0/4.0;
    lambda[2] = 1.0;
    lambda[3] = 5.0/32.0;
    lambda[4] = 1.0;
    lambda[5] = 1.0/20.0;
    lambda[6] = 1.0/10.0;
    lambda[7] = 1.0;
    lambda[8] = 1.0/20.0;
    lambda[9] = 1.0/20.0;
    global_bias = 260.0;
    return;
}

void calc_benchmark_norm()
{
    int i;
    long double temp1, temp2, temp;
    for (i=0; i<nreal; i++)
    {
        norm_x[i] = 0.0;
    }
    transform_norm (0);    norm_f[0] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform_norm (1);
    norm_f[1] = 0.0;
    for (i=0; i<nreal-1; i++)
    {
        temp1 = pow((sin(sqrt(pow(trans_x[i],2.0)+pow(trans_x[i+1],2.0)))),2.0);
        temp2 = 1.0 + 0.001*(pow(trans_x[i],2.0)+pow(trans_x[i+1],2.0));
        norm_f[1] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
    }
    temp1 = pow((sin(sqrt(pow(trans_x[nreal-1],2.0)+pow(trans_x[0],2.0)))),2.0);
    temp2 = 1.0 + 0.001*(pow(trans_x[nreal-1],2.0)+pow(trans_x[0],2.0));
    norm_f[1] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
    transform_norm (2);
    norm_f[2] = 0.0;
    for (i=0; i<nreal-1; i++)
    {
        temp = 100.0*pow((trans_x[i]*trans_x[i]-trans_x[i+1]),2.0) + 1.0*pow((trans_x[i]-1.0),2.0);
        norm_f[2] += (temp*temp)/4000.0 - cos(temp) + 1.0;
    }
    temp = 100.0*pow((trans_x[nreal-1]*trans_x[nreal-1]-trans_x[0]),2.0) + 1.0*pow((trans_x[nreal-1]-1.0),2.0);
    norm_f[2] += (temp*temp)/4000.0 - cos(temp) + 1.0;
    transform_norm (3);    norm_f[3] = calc_ackley(trans_x);
    transform_norm (4);    norm_f[4] = calc_rastrigin(trans_x);
    transform_norm (5);    norm_f[5] = calc_griewank(trans_x);
    transform_norm (6);
    norm_f[6] = 0.0;
    for (i=0; i<nreal-1; i++)
    {
        norm_f[6] += nc_schaffer(trans_x[i], trans_x[i+1]);
    }
    norm_f[6] += nc_schaffer(trans_x[nreal-1], trans_x[0]);
    transform_norm(7);    norm_f[7] = nc_rastrigin(trans_x);
    transform_norm (8);
    norm_f[8] = 0.0;
    for (i=0; i<nreal; i++)
    {
        norm_f[8] += trans_x[i]*trans_x[i]*pow(1.0e6,i/(nreal-1.0));
    }
    transform_norm (9);    norm_f[9] = calc_sphere(trans_x)*(1.0 + 0.1*fabs(randomnormaldeviate()));
    return;
}

long double calc_benchmark_func(long double *x)
{
    int i;
    long double temp1, temp2, temp;
    long double res;
    for (i=0; i<nreal; i++)
    {
        norm_x[i] = 0.0;
    }
    /* First function */
    transform (x, 0);    basic_f[0] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);

    /* Second function */
    transform (x, 1);
    basic_f[1] = 0.0;
    for (i=0; i<nreal-1; i++)
    {
        temp1 = pow((sin(sqrt(pow(trans_x[i],2.0)+pow(trans_x[i+1],2.0)))),2.0);
        temp2 = 1.0 + 0.001*(pow(trans_x[i],2.0)+pow(trans_x[i+1],2.0));
        basic_f[1] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
    }
    temp1 = pow((sin(sqrt(pow(trans_x[nreal-1],2.0)+pow(trans_x[0],2.0)))),2.0);
    temp2 = 1.0 + 0.001*(pow(trans_x[nreal-1],2.0)+pow(trans_x[0],2.0));
    basic_f[1] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));

    /* Third Function */
    transform (x, 2);
    basic_f[2] = 0.0;
    for (i=0; i<nreal-1; i++)
    {
        temp = 100.0*pow((trans_x[i]*trans_x[i]-trans_x[i+1]),2.0) + 1.0*pow((trans_x[i]-1.0),2.0);
        basic_f[2] += (temp*temp)/4000.0 - cos(temp) + 1.0;
    }
    temp = 100.0*pow((trans_x[nreal-1]*trans_x[nreal-1]-trans_x[0]),2.0) + 1.0*pow((trans_x[nreal-1]-1.0),2.0);
    basic_f[2] += (temp*temp)/4000.0 - cos(temp) + 1.0;

    transform (x, 3);    basic_f[3] = calc_ackley(trans_x);
    transform (x, 4);    basic_f[4] = calc_rastrigin(trans_x);
    transform (x, 5);    basic_f[5] = calc_griewank(trans_x);

    /* Seventh Function */
    transform (x, 6);
    basic_f[6] = 0.0;
    for (i=0; i<nreal-1; i++)
    {
        basic_f[6] += nc_schaffer(trans_x[i], trans_x[i+1]);
    }
    basic_f[6] += nc_schaffer(trans_x[nreal-1], trans_x[0]);

    transform (x, 7);    basic_f[7] = nc_rastrigin(trans_x);

    transform (x, 8);
    basic_f[8] = 0.0;
    for (i=0; i<nreal; i++)
    {
        basic_f[8] += trans_x[i]*trans_x[i]*pow(1.0e6,i/(nreal-1.0));
    }
    transform (x, 9);    basic_f[9] = (calc_sphere(trans_x))*(1.0 + 0.1*fabs(randomnormaldeviate()));
    for (i=0; i<nfunc; i++)
    {
        basic_f[i] *= C/norm_f[i];
    }
    calc_weight(x);
    res = global_bias;
    for (i=0; i<nfunc; i++)
    {
        res += weight[i]*(basic_f[i]+bias[i]);
    }
    return (res);
}
#endif

