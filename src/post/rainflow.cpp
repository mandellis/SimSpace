#include "rainflow.h"
#include <cstdlib>
#include <iostream>
#include <fstream>

#include <vector>
#include <math.h>
#include <stdio.h>
#include <cmath>

#include <string.h>
#include <stdlib.h>
#include <conio.h>

#include <time.h>
#include <sys/types.h>
#include <sys/timeb.h>

/********************************************************************
 * Copyright:   (c) 2007-2008 Ladisk <www.fs.uni-lj.si/ladisk>
 * Author:      Primoz Cermelj <primoz.cermelj@gmail.com>
 *********************************************************************/
//! ----------------------
//! function: constructor
//! details:
//! ----------------------
rainflow::rainflow(QObject *parent): QObject(parent)
{
   ;
}

//! ------------------------
//! function: constructor I
//! details:
//! ------------------------
//rainflow::rainflow(GeometryTag loc, QObject *parent): QObject(parent)
//{
//    ;
//}

//! ------------------------
//! function: constructor II
//! details:
//! ------------------------
//rainflow::rainflow(GeometryTag loc, fatigueModel fm, QObject *parent): QObject(parent),
//   myLoc(loc),myFatigueModel(fm)
//{
 //   ;
//}
/*-------------------------------------------------------------
 * rf3
 *
 * Performs rainflow analysis without time analysis and returns
 * the actual number of rows in the output rf matrix - the
 * actual size of the data in array_out is (Nr x 3).
 *
 * Based on Adam Nieslony's rainflow.c for Matlab.
 *-------------------------------------------------------------*/
int rainflow::rf3(double *array_ext,	// (in) an array of turning points (see sig2ext on how to get these)
        int nr,				// (in) length of the array_ext (number of rows of the vector)
        double *array_out)	// (out) output matrix of size nr x 3; the columns are:
                            // 		cycles amplitude, cycles mean value, number of
                            // 		cycles (0.5 or 1.0). This array must be allocated
                            // 		apriori.
{
  double a[512], ampl, mean;
  int index, j, cNr, tot_num;

  tot_num = nr;

  // Init array_out to zero
  for (index=0; index<nr*3; index++)
  {
    array_out[index] = 0.0;
  }

  j = -1;
  cNr = 1;
  for (index=0; index<tot_num; index++)
  {
    a[++j]=*array_ext++;
    while ( (j >= 2) && (fabs(a[j-1]-a[j-2]) <= fabs(a[j]-a[j-1])) ) {
      ampl = fabs( (a[j-1]-a[j-2])/2 );
      switch(j){
        case 0: { break; }
        case 1: { break; }
        case 2: {
          mean=(a[0]+a[1])/2;
          a[0]=a[1];
          a[1]=a[2];
          j=1;
          if (ampl > 0) {
            *array_out ++= ampl;
            *array_out ++= mean;
            *array_out ++= 0.50;
          }
          break;
        }
        default: {
          mean = (a[j-1]+a[j-2])/2;
          a[j-2] = a[j];
          j = j-2;
          if (ampl > 0) {
            *array_out ++= ampl;
            *array_out ++= mean;
            *array_out ++= 1.00;
            cNr++;
          }
          break;
        }
      }
    }
  }
  for (index=0; index<j; index++)
  {
    ampl = fabs(a[index]-a[index+1])/2;
    mean = (a[index]+a[index+1])/2;
    if (ampl > 0)
    {
      *array_out ++= ampl;
      *array_out ++= mean;
      *array_out ++= 0.50;
    }
  }

  return (tot_num - 1 - cNr);
}


/*-------------------------------------------------------------
 * rf5
 *
 * Performs rainflow analysis with time analysis.
 *
 * Inputs:
 *      array_ext   an array of turning points (see sig2ext on how to get these)
 *      nr          length of the array_ext and  number of rows in array_out
 *      array_t     an array of time values
 *
 * Outputs:
 *      cnr         the actual number of rows in the rf matrix
 *      array_out   matrix of length nr x 5 where the result will be returned;
 *                  this array must be pre-allocated apriori.
 *
 * Based on Adam Nieslony's rainflow.c for Matlab.
 *-------------------------------------------------------------*/
int rainflow::rf5(double *array_ext, int nr, double *array_t, double *array_out){
  double a[512], t[512], ampl, mean, period, atime;
  int index, j, cNr, tot_num;

  tot_num = nr;

  // Init array_out to zero
  for (index=0; index<nr*5; index++) {
    array_out[index] = 0.0;
  }

  j = -1;
  cNr = 1;
  for (index=0; index<tot_num; index++) {
    a[++j]=*array_ext++;
    t[j]=*array_t++;
    while ( (j >= 2) && (fabs(a[j-1]-a[j-2]) <= fabs(a[j]-a[j-1])) ) {
      ampl=fabs( (a[j-1]-a[j-2])/2 );
      switch(j) {
        case 0: { break; }
        case 1: { break; }
        case 2: {
          mean=(a[0]+a[1])/2;
          period=(t[1]-t[0])*2;
          atime=t[0];
          a[0]=a[1];
          a[1]=a[2];
          t[0]=t[1];
          t[1]=t[2];
          j=1;
          if (ampl > 0) {
            *array_out++=ampl;
            *array_out++=mean;
            *array_out++=0.50;
            *array_out++=atime;
            *array_out++=period;
          }
          break;
        }
        default: {
          mean=(a[j-1]+a[j-2])/2;
          period=(t[j-1]-t[j-2])*2;
          atime=t[j-2];
          a[j-2]=a[j];
          t[j-2]=t[j];
          j=j-2;
          if (ampl > 0) {
            *array_out++=ampl;
            *array_out++=mean;
            *array_out++=1.00;
            *array_out++=atime;
            *array_out++=period;
            cNr++;
          }
          break;
        }
      }
    }
  }
  for (index=0; index<j; index++) {
    ampl=fabs(a[index]-a[index+1])/2;
    mean=(a[index]+a[index+1])/2;
    period=(t[index+1]-t[index])*2;
    atime=t[index];
    if (ampl > 0){
      *array_out++=ampl;
      *array_out++=mean;
      *array_out++=0.50;
      *array_out++=atime;
      *array_out++=period;
    }
  }

  return (tot_num - 1 - cNr);
}


/*-------------------------------------------------------------
 * sig2ext
 *
 * Searches local extrema from time course (signal) sig.
 *
 * Inputs:
 *      sig         an array of n time points
 *      time_sig    an array of n delta time points (pass NULL if this array
 *                  is to be assumed in the form of 0, 1, 2,....)
 *      n           number of points (sig, time_sig)
 *      clsn        number of classes (pass -1 if not to be used, i.e., no
 *                  divisions into classes)
 *      ext         (output) extrema found on sig
 *      exttime     (output) time values corresponding to ext;
 *                  if time was NULL, a dt=1 is assumed.
 * Outputs:
 *      np          number of extrema (number of points on the output)
 *
 * Based on Adam Nieslony's sig2ext.m (Matlab function).
 *-------------------------------------------------------------*/
int rainflow::sig2ext(double *sig, double *time_sig, int n, int clsn,
            double *ext, double *exttime)
{
    int i, have_time = 0;
    double smax, smin;
    double *w1;
    int *w;
    int np;

    if (time_sig != NULL)
    {
        have_time = 1;
    }
    if (clsn != -1)
    {
        clsn--;

        smin = this->arr_min(sig, n, 0);
        smax = this->arr_max(sig, n, NULL);
        for (i=0; i<n; i++)
        {
            sig[i] = round((sig[i]-smin)*clsn/(smax-smin))*(smax-smin)/clsn + smin;
        }
    }

    if (have_time != 1)
    {
        time_sig = NNEW(double, n);
        for (i=0; i<n; i++)
        {
            time_sig[i] = i;
        }
    }

    w = NNEW(int, n);
    w1 = this->diff(sig, n);
    for (i=1; i<n; i++)
    {
        if (w1[i-1]*w1[i] <= 0) w[i] = 1;
    }
    w[0] = 1;
    w[n-1] = 1;
    np = this->repl(ext, w, n, sig);
    np = this->repl(exttime, w, n, time_sig);

    free(w1);
    RENEW(w, int, np);
    w1 = this->diff(ext, np);
    for (i=1; i<np; i++)
    {
        if (!(w1[i-1]==0 && w1[i]==0)) w[i] = 1;
    }
    w[0] = 1;
    w[np-1] = 1;
    np = this->repl(ext, w, np, ext);
    np = this->repl(exttime, w, np, exttime);

    for (i=1; i<np; i++)
    {
        if (ext[i-1] != ext[i]) w[i] = 1;
    }
    w[0] = 1;
    np = this->repl(ext, w, np, ext);
    RENEW(w1, double, np-1);
    for (i=0; i<np-1; i++)
    {
        w1[i] = 0.5*( exttime[i+1] - exttime[i]);
        exttime[i] = exttime[i] + w1[i]*(!w[i+1]);
    }
    np = this->repl(exttime, w, np, exttime);

    if (np > 2)
    {
        free(w1);
        RENEW(w, int, np);
        w1 = this->diff(ext, np);
        for (i=1; i<np; i++){
            if (w1[i-1]*w1[i] < 0) w[i] = 1;
        }
        w[0] = 1;
        w[np-1] = 1;
        np = this->repl(ext, w, np, ext);
        np = this->repl(exttime, w, np, exttime);
    }

    free(w1);
    free(w);
    if (have_time != 1) free(time_sig);

    return np;
}

/*-------------------------------------------------------------
 * arr_min
 *
 * Returns minimum of a double vector. If the pos is given as a
 * non-NULL pointer, the position of the min value is returned
 * as well.
 *-------------------------------------------------------------*/
double rainflow::arr_min(double *sig, int n, int *pos)
{
    int i, ind=0;

    double min_val = sig[0];
    for (i=1; i<n; i++)
    {
        if (sig[i] < min_val)
        {
            min_val = sig[i];
            ind = i;
        }
    }
    if (pos != NULL) *pos = ind;
    return min_val;
}


/*-------------------------------------------------------------
 * arr_max
 *
 * Returns maximum of a double vector. If the pos is given as a
 * non-NULL pointer, the position of the max value is returned
 * as well.
 *-------------------------------------------------------------*/
double rainflow::arr_max(double *sig, int n, int *pos){
    int i, ind=0;

    double max_val = sig[0];
    for (i=1; i<n; i++)
    {
        if (sig[i] > max_val)
        {
            max_val = sig[i];
            ind = i;
        }
    }
    if (pos != NULL) *pos = ind;
    return max_val;
}

/*-------------------------------------------------------------
 * diff
 *
 * Returns a new, dynamically allocated vector with the length
 * n-1 and values defined as:
 *      v[2]-v[1]
 *      v[3]-v[2]
 *      ....
 *
 * Make sure the vector is cleared with free when no longer
 * needed.
 *-------------------------------------------------------------*/
double *rainflow::diff(double *vec, int n)
{
    // The length of vec_out is n-1!
    int i;
    double *vec_out;

    vec_out = NNEW(double, n-1);
    for (i=0; i<(n-1); i++)
    {
        vec_out[i] = vec[i+1] - vec[i];
    }
    return vec_out;
}


/*-------------------------------------------------------------
 * repl
 *
 * Replaces x with x_repl where filt is 1. x and filt are of the length n
 * while x_repl can be n or greater. x is replaced inplace. It returns
 * the number of terms replaced.
 *-------------------------------------------------------------*/
int rainflow::repl(double *x, int *filt, int n, double *x_repl)
{
    int i, j;
    j = 0;
    for (i=0; i<n; i++){
        if (filt[i] > 0){
            x[j] = x_repl[i];
            j++;
        }
    }
    return j;
}

//! ------------------
//! function: perform
//! details:
//! ------------------
bool rainflow::perform(std::map<int, std::vector<double>> strainDistTimeHistory, std::map<int,double> &damageDist)
{
    cout<<"rainflow::perform()->____function called____"<<endl;

    //!------------------------------------------
    //! deltaeps/2= espF*(2N)^c+sigmaF/E*(2N)^b
    //!              [0]     [1]   [2] [3]   [4]
    //!------------------------------------------
    double epsF = myFatigueModel.coeffs[0];
    double c = myFatigueModel.coeffs[1];
    double sigmaF = myFatigueModel.coeffs[2];
    double E = myFatigueModel.coeffs[3];
    double b = myFatigueModel.coeffs[4];

    for(std::map<int, std::vector<double>>::iterator it = strainDistTimeHistory.begin(); it!= strainDistTimeHistory.end(); it++)
    {
        int nodeID = it->first;
        std::vector<double> timeHistory = it->second;
        size_t timeSize = timeHistory.size();
        double *tH= NNEW(double, timeSize);
        for (int i=0; i<(timeSize); i++)
        {
            tH[i] = timeHistory[i];
        }

        double damage=0.0;
        double *ext=NNEW(double,timeSize*3);
        int np = this->rf3(tH,timeSize,ext);

        for(int i=0; i<np; i++)
        {
            double deltaEps = ext[3*i];
            //if(i==0) cout<<"rainflow::damage_index()->____ext value____"<<ext[3*i]<<"___"<<ext[3*i+1]<<"___"<<ext[3*i+2]<<endl;
            if(deltaEps>tol2)
            {
                double NN = solve_exact(deltaEps,epsF,c,sigmaF,E,b)+1e-12;
                //double NN = solve(deltaEps,epsF,c,sigmaF,E,b);
                damage+= 1/NN;
            }
        }
        //double damage = this->damage_index(tH,timeSize);
        damageDist.insert(std::make_pair(nodeID,damage));
    }
    return true;
}

//! -------------------------------------------------
//! function: solveEquation
//! details:  returns the number of cycles @ failure
//! -------------------------------------------------
double rainflow::solve(double eps, double epsF, double c, double sigmaF, double E, double b)
{
    //cout<<"rainflow::solve()->____function called____"<<endl;
    double N,Nn,initialRoot;

    switch(myFatigueModel.type)
    {
    case fatigueModel_BCM:
    {
        /*
        //!------------------------------------------------------------------
        //! Bisection Method: choose initial guess for Newton Raphson Method
        //! details:
        //!------------------------------------------------------------------
        int loopCounter=0;
        int initialRoot=0;
        int N1=1;    //initial value for f(N1)
        int N2=1e9;  //initial value for f(N2)
        int N3=0;
        int maxIter=5; // max iteration of Bisection Method
        int rootRange=0.5; // criteria for acceptance of initial guess for bisection method

        double f1 = eps-epsF*pow(2*N1,c)-(sigmaF/E)*pow(2*N1,b);
        double f2 = eps-epsF*pow(2*N2,c)-(sigmaF/E)*pow(2*N2,b);

        cout<<"rainflow::solve()->____bisection method, eps = "<<eps<<"____"<<endl;
        //cout<<"rainflow::solve()->____bisection method, f1 = "<<f1<<"____"<<endl;
        //cout<<"rainflow::solve()->____bisection method, f2 = "<<f2<<"____"<<endl;

        if(f1*f2<0)
        {
            while(1)
            {
                loopCounter++;
                cout<<"rainflow::solve()->____bisection method, enter while loppCounter = "<<loopCounter<<"____"<<endl;

                if(loopCounter<maxIter)
                {
                    N3=(N1+N2)/2;
                    double f3=eps-epsF*pow(2*N3,c)-(sigmaF/E)*pow(2*N3,b);
                    cout<<"rainflow::solve()->____bisection method, f3 = "<<f3<<"____"<<endl;

                    if(f3<rootRange && f3>-rootRange)
                    {
                        initialRoot=N3;
                        break;
                    }
                    if(f1*f3 < 0)
                    {
                        N2=N3;
                    }
                    if(f1*f3 > 0)
                    {
                        N1=N3;
                    }
                }
                else
                {
                    initialRoot=N3;
                    break;
                }
            }
        }
        else
        {
           initialRoot=N1;
        }
        cout<<"rainflow::solve()->____bisection method, initial calculated Root for NR = "<<initialRoot<<"____"<<endl;
        */

        //!---------------------------------------------------
        //! Newton Raphson Method: calculate root of equation
        //! details:
        //!---------------------------------------------------

        //! estimate initial root
        double level1=3.53e-3;
        double level2=3.567e-3;
        double level3=3.61e-3;
        double level4=3.651e-3;
        double level5=3.72e-3;
        double level6=3.8e-3;
        double level7=3.9e-3;
        double level8=4.0e-3;
        double level9=4.16e-3;
        double level10=4.37e-3;
        double level11=4.7e-3;
        double level12=9.64e-3;

        //add more level
        if(eps>tol2 && eps<level1) initialRoot=60000;
        else if(eps>level1 && eps<level2) initialRoot=55000;
        else if(eps>level2 && eps<level3) initialRoot=50000;
        else if(eps>level3 && eps<level4) initialRoot=45000;
        else if(eps>level4 && eps<level5) initialRoot=40000;
        else if(eps>level5 && eps<level6) initialRoot=33500;
        else if(eps>level6 && eps<level7) initialRoot=28000;
        else if(eps>level7 && eps<level8) initialRoot=22500;
        else if(eps>level8 && eps<level9) initialRoot=18500;
        else if(eps>level9 && eps<level10) initialRoot=14000;
        else if(eps>level10 && eps<level11) initialRoot=7500;
        else if(eps>level11 && eps<level12) initialRoot=500;
        else if(eps>level12) initialRoot=50;
        N=initialRoot;
        //cout<<"rainflow::solve()->____N = "<<N<<"____"<<endl;

        int loopCounter=0;
        while(loopCounter<20)
        {
            double f=eps-epsF*pow(2*N,c)-(sigmaF/E)*pow(2*N,b);
            double df=-2*c*epsF*pow(2*N,c-1)-2*b*(sigmaF/E)*pow(2*N,b-1);
            //cout<<"rainflow::solve()->____f = "<<f<<"____"<<endl;
            //cout<<"rainflow::solve()->____df = "<<df<<"____"<<endl;

            if(df!=0.0 /*|| fabs(df)>1e12*/)
            {
                Nn=N-f/df;
                //cout<<"rainflow::solve()->____NN = "<<Nn<<"____"<<endl;

                double err = abs((Nn-N));
                //cout<<"rainflow::solve()->____err = "<<err<<"____"<<endl;
                if(err<maxErr)
                {
                    break;
                }
                else
                {
                    N=Nn;
                    loopCounter++;
                    //cout<<"rainflow::solve()->____loop counter= "<<loopCounter<<"____"<<endl;
                }
            }
            else
            {
                //cout<<"rainflow::solve()->____break____"<<endl;
                break;
            }
        }
    }
        break;

    default:
        break;
    }
    return N;
}

//! ----------------------
//! function: solve_exact
//! details:
//! ----------------------
double rainflow::solve_exact(double eps, double epsF, double c, double sigmaF, double E, double b)
{
    //cout<<"rainflow::solve_exact()->____function called____"<<endl;
    double N;
    switch(myFatigueModel.type)
    {
    case fatigueModel_BCM:
    {
        /*
        double den = c-b; if(den == 0) den = 1e-12;
        double x = (1/(c-b))*log((eps/epsF)*(sigmaF/E));
        */
        double den = epsF+sigmaF/E*exp(b/c);
        //N = 0.5*eps/den;
        N = 0.5*pow(eps/den,1/c);
    }
        break;
    }
    return N;
}
