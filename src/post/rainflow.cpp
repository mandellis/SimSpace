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

#include <QVector>
#include <QList>

#include <time.h>
#include <sys/types.h>
#include <sys/timeb.h>

#define MAX 51000000

using namespace std;

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
rainflow::rainflow(GeometryTag loc, QObject *parent): QObject(parent),myLoc(loc)
{
    ;
}

//! ------------------------
//! function: constructor II
//! details:
//! ------------------------
rainflow::rainflow(GeometryTag loc, fatigueModel fm, QObject *parent): QObject(parent),
    myLoc(loc),myFatigueModel(fm)
{
    ;
}

//! -------------------------------------------------
//! function: solveEquation
//! details:  returns the number of cycles @ failure
//! -------------------------------------------------
double rainflow::solve(double eps, double epsF, double c, double sigmaF, double E, double b)
{
    cout<<"rainflow::solve()->____function called____"<<endl;
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
        cout<<"rainflow::solve()->____function called____"<<endl;
        double N;
        switch(myFatigueModel.type)
        {
        case fatigueModel_BCM:
        {
            double den = c-b; if(den == 0) den = 1e-12;
            double x = (1/(c-b))*log((eps/epsF)*(sigmaF/E));
            N = 0.5*exp(x);
        }
            break;
        }
        return N;
}

//! -----------------------
//! function: damage_index
//! details:
//! -----------------------
double rainflow::damage_index(const std::vector<double> &y)
{
    cout<<"rainflow::damage_index()->____function called____"<<endl;

    //!------------------------------------------
    //! deltaeps/2= espF*(2N)^c+sigmaF/E*(2N)^b
    //!              [0]     [1]   [2] [3]   [4]
    //!------------------------------------------
    double epsF = myFatigueModel.coeffs.at(0);
    double c = myFatigueModel.coeffs.at(1);
    double sigmaF = myFatigueModel.coeffs.at(2);
    double E = myFatigueModel.coeffs.at(3);
    double b = myFatigueModel.coeffs.at(4);

    double D=0.;
    double deltaEps;
    std::vector<double> B = this->rainflow_engine(y);
    size_t NbData = B.size();
    for(size_t i=0; i<NbData; i++)
    {
        deltaEps = B[i];
        if(deltaEps>tol2)
        {
            //double NN_ = solve_exact(deltaEps,epsF,c,sigmaF,E,b);
            double NN = solve(deltaEps,epsF,c,sigmaF,E,b);
            D+= 1/NN;
        }
    }
    return D;
}

//! --------------------------
//! function: rainflow_engine
//! details:
//! --------------------------
std::vector<double> rainflow::rainflow_engine(std::vector<double> y)
{
    cout<<"rainflow::rainflow_engine()->____function called____"<<endl;
    double sum;
    double ymax;
    double mina,maxa;
    double X,Y;

    int kv;
    int hold;
    int i,j,k,n;
    int num;
    int nkv;
    int last_a;

    std::vector<double> B;
    std::vector<double> a;

    double slope1;
    double slope2;

    ymax=0.;

    nkv=0;
    k=0;
    //	a[k]=y[k];
    a.push_back(y.at(k));

    k=1;    // ??

    int NP = int(y.size());
    for(i=1;i<(NP-1);i++)
    {
        slope1=(y.at(i)-y.at(i-1));
        slope2=(y.at(i+1)-y.at(i));

        if((slope1*slope2)<=0 && fabs(slope1)>0)
        {
            a.push_back(y.at(i));
            k++;
        }
    }
    a.push_back(y.at(NP-1));
    k++;

    last_a=k-1;
    hold=last_a;

    mina = 1.0e10;
    maxa = -mina;
    for(i=0;i<=last_a;i++)
    {
        if(a.at(i)<mina)
        {
            mina=a.at(i);
        }
        if(a.at(i)>maxa)
        {
            maxa=a.at(i);
        }
    }

    num=int(maxa-mina)+1;

    n=0;
    i=0;
    j=1;

    sum=0;

    kv=0;

    //int LLL=last_a;

    QVector<double> row;
    row.resize(4);

    while(1)
    {
        Y=(fabs(a.at(i)-a.at(i+1)));
        X=(fabs(a.at(j)-a.at(j+1)));
        if(X>=Y && Y>0 && Y<1.0e+10)
        {
            if(Y>ymax)
            {
                ymax=Y;
            }
            if(i==0)
            {
                n=0;
                sum+=0.5;

                row.insert(3,a.at(i+1));
                row.insert(2,a.at(i));
                row.insert(1,0.5);
                row.insert(0,Y);

                B.push_back(Y);
                kv++;
                a.erase(a.begin());
                last_a--;

                i=0;
                j=1;
            }
            else
            {
                sum+=1;
                row.insert(3,a.at(i+1));
                row.insert(2,a.at(i));
                row.insert(1,1);
                row.insert(0,Y);

                B.push_back(Y);

                kv++;
                n=0;

                a.erase(a.begin()+(i+1));
                a.erase(a.begin()+i);

                last_a-=2;

                i=0;
                j=1;
            }

            nkv++;

            /*
            if(nkv==3000)
            {
                double ratio = fabs((last_a)/double(LLL));
                nkv=0;
            }
            */
        }
        else
        {
            i++;
            j++;
        }

        if((j+1)>last_a)
        {
            break;
        }
    }


    for(i=0;i<(last_a);i++)
    {
        Y=(fabs(a.at(i)-a.at(i+1)));
        if(Y>0. && Y<1.0e+20)
        {
            sum+=0.5;

            row.insert(3,a.at(i+1));
            row.insert(2,a.at(i));
            row.insert(1,0.5);
            row.insert(0,Y*.5);
            //B.push_back(row);
            B.push_back(Y);

            kv++;

            if(Y>ymax) ymax=Y;
        }
    }
    return B;
}

//! ------------------
//! function: perform
//! details:
//! ------------------
bool rainflow::perform(std::map<int, std::vector<double>> strainDistTimeHistory, std::map<int,double> &damageDist)
{
    cout<<"rainflow::perform()->____function called____"<<endl;

    for(std::map<int, std::vector<double>>::iterator it = strainDistTimeHistory.begin(); it!= strainDistTimeHistory.end(); it++)
    {
        int nodeID = it->first;
        std::vector<double> timeHistory = it->second;
        double damage = this->damage_index(timeHistory);
        damageDist.insert(std::make_pair(nodeID,damage));
    }
    return true;
}
