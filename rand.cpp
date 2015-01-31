
// Definition of random number generation routines

# include "rand.h"

# include <iostream>
# include <cstdlib>
# include <cmath>

using namespace std;

/* Variable declarations for the random number generator */
long double seed;
long double oldrand[55];
int jrand;
int rndcalcflag;
long double rndx1, rndx2;

/* Get seed number for random and start it up */
void randomize()
{
      int j1;
      for(j1=0; j1<=54; j1++)
      {
            oldrand[j1] = 0.0;
      }
      jrand=0;
      warmup_random (seed);
      return;
}

/* Get randomize off and running */
void warmup_random (long double seed)
{
      int j1, ii;
      long double new_random, prev_random;
      oldrand[54] = seed;
      new_random = 0.000000001;
      prev_random = seed;
      for(j1=1; j1<=54; j1++)
      {
            ii = (21*j1)%54;
            oldrand[ii] = new_random;
            new_random = prev_random-new_random;
            if(new_random<0.0)
            {
                  new_random += 1.0;
            }
            prev_random = oldrand[ii];
      }
      advance_random ();
      advance_random ();
      advance_random ();
      jrand = 0;
      return;
}

/* Create next batch of 55 random numbers */
void advance_random ()
{
      int j1;
      long double new_random;
      for(j1=0; j1<24; j1++)
      {
            new_random = oldrand[j1]-oldrand[j1+31];
            if(new_random<0.0)
            {
                  new_random = new_random+1.0;
            }
            oldrand[j1] = new_random;
      }
      for(j1=24; j1<55; j1++)
      {
            new_random = oldrand[j1]-oldrand[j1-24];
            if(new_random<0.0)
            {
                  new_random = new_random+1.0;
            }
            oldrand[j1] = new_random;
      }
}

/* Fetch a single random number between 0.0 and 1.0 */
long double randomperc()
{
      jrand++;
      if(jrand>=55)
      {
            jrand = 1;
            advance_random();
      }
      return((long double)oldrand[jrand]);
}

/* Fetch a single random integer between low and high including the bounds */
int rnd (int low, int high)
{
    int res;
    if (low >= high)
    {
        res = low;
    }
    else
    {
        res = low + (randomperc()*(high-low+1));
        if (res > high)
        {
            res = high;
        }
    }
    return (res);
}

/* Fetch a single random real number between low and high including the bounds */
long double rndreal (long double low, long double high)
{
    return (low + (high-low)*randomperc());
}

/* Initialize the randome generator for normal distribution */
void initrandomnormaldeviate()
{
    rndcalcflag = 1;
    return;
}

/* Return the noise value */
long double noise (long double mu, long double sigma)
{
    return((randomnormaldeviate()*sigma) + mu);
}

/* Compute the noise */
long double randomnormaldeviate()
{
    long double t;
    if(rndcalcflag)
    {
        rndx1 = sqrt(- 2.0*log(randomperc()));
        t = 6.2831853072*randomperc();
        rndx2 = sin(t);
        rndcalcflag = 0;
        return(rndx1*cos(t));
    }
    else
    {
        rndcalcflag = 1;
        return(rndx1*rndx2);
    }
}
