
//  Source file for CEC2005 benchmark functions

/* Global variable and function declarations */

#ifndef CEC2005_BENCHMARKS_H_INCLUDED
#define CEC2005_BENCHMARKS_H_INCLUDED


/* Global Constants */
# define INF 1.0e99
# define EPS 1.0e-14
# define E  2.7182818284590452353602874713526625
# define PI 3.1415926535897932384626433832795029

/* Uncomment one of the lines corresponding to a particular function */
/* Function identifier */
 # define f1
// # define f2
// # define f3
// # define f4
// # define f5
// # define f6
// # define f7
// # define f8
// # define f9
// # define f10
// # define f11
// # define f12
// # define f13
// # define f14
// # define f15
// # define f16
// # define f17
// # define f18
 //# define f19
// # define f20
// # define f21
// # define f22
// # define f23
// # define f24
// # define f25

/* Global variables that you are required to initialize */
extern int nreal;                    /* number of real variables */
extern int nfunc;                    /* number of basic functions */

/* Auxillary function declarations */
long double maximum (long double, long double);
long double minimum (long double, long double);
long double modulus (long double*, int);
long double dot (long double*, long double*, int);
long double mean (long double*, int);

/* Basic funcion declarations */
long double calc_ackley (long double*);
long double calc_rastrigin (long double*);
long double calc_weierstrass (long double*);
long double calc_griewank (long double*);
long double calc_sphere (long double*);
long double calc_schwefel (long double*);
long double calc_rosenbrock (long double *x);
long double nc_schaffer (long double, long double);
long double nc_rastrigin (long double*);

/* Utility function declarations */
void allocate_memory();
void initialize();
void transform (long double*, int);
void transform_norm (int);
void calc_weight (long double*);
void free_memory();

/* Benchmark function declaration */
long double calc_benchmark_func (long double*);
void calc_benchmark_norm();

#endif // CEC2005_BENCHMARKS_H_INCLUDED
