//! custom includes
#include "posttools.h"
#include "cubicequation.h"
#include "postengine.h"

//! Qt
#include <QList>

//! C++
#include <iostream>
#include <fstream>
using namespace std;

postTools::postTools()
{
    ;
}

//! -----------------------------------------
//! 1) global iteration number
//! 2) step
//! 3) increment
//! 4) attempt
//! 5) actual step time
//! 6) actual total time
//! 7) average force
//! 8) largest residual force
//! 9) largest increment of displacement
//! 10) largest correction to displacement
//! 11) substep converged
//! -----------------------------------------
/*
//! ----------------------------------------------------------------------------------------
//! function: getStepSubStepByTime
//! details:
//! ----------------------------------------------------------------------------------------
bool postTools::getStepSubStepByTime(const QString &SolverOutputFile_path, double requiredTime, int &foundStep, int &foundSubStep)
{
    QList<double> timeList;
    QList<int> stepList;

    QList<int> subStepList;
    //QList<int> incrementList;

    int git;        //! global iteration number
    int stepNumber; //! step number
    int inc;        //! increment
    int att;        //! attempt
    double ast;     //! actual step time
    double time;    //! actual total time
    double avf;     //! average force
    double lrs;     //! largest residual force
    double lid;     //! largest increment of displacement
    double lcd;     //! largest correction to displacement
    bool con;       //! substep converged

    std::string val;
    ifstream is(SolverOutputFile_path.toStdString());
    if(!is.fail())
    {
        int k=0;
        //int subStep = 0;
        std::getline(is,val);

        for(;!is.eof();)
        {
            int N = sscanf(val.c_str(),"%d%d%d%d%lf%lf%lf%lf%lf%lf%d",&git,&stepNumber,&inc,&att,&ast,&time,&avf,&lrs,&lid,&lcd,&con);
            if(N!=11)
            {
                if(k==0)
                {
                    cout<<"LINE ERROR AT FIRST LINE"<<endl;
                    //! line error at first line
                    foundStep = -1;
                    foundSubStep = -1;
                    return false;
                }
                else
                {
                    cout<<"LINE ERROR"<<endl;
                    //! return the previously read data
                    return true;
                }
            }
            else    //! read exactly 10 data
            {
                if(con==true)
                {
                    //cout<<"convergence found"<<endl;
                    timeList<<time;

                    //! should be equivalent
                    stepList<<stepNumber;
                    subStepList<<inc;
                }
            }
            k++;
            std::getline(is,val);
        }

        //! diagnostic
        for(int i=0;i<timeList.length();i++)
        {
            cout<<"(Time = "<<timeList.at(i)<<", Step: "<<stepList.at(i)<<", subStep = "<<subStepList.at(i)<<")"<<endl;
        }

        //! end diagnostic
        //! process the solution info
        bool found = false;
        for(int i=1; i<=timeList.length()-1; i++)
        {
            double lower = timeList.at(i-1);
            double upper = timeList.at(i);

            if(requiredTime>=lower && requiredTime<upper)
            {
                cout<<"time point within an inner interval"<<endl;
                //! introduce average ... to do
                foundStep = stepList.at(i-1);
                foundSubStep = subStepList.at(i-1);
                found = true;
                break;
            }
            else if(requiredTime>0.0 && requiredTime<lower)
            {
                //! if the 0<=t<first time => force return the first substep
                cout<<"required time within first interval"<<endl;
                foundStep = 1;
                foundSubStep = 1;
                found = true;
                break;
            }
            else if(requiredTime==0.0)
            {
                //! if the retuiredTime==0 force return the last
                cout<<"returning the last"<<endl;
                foundStep = stepList.last();
                foundSubStep = subStepList.last();
                found = true;
                break;
            }
            else if(time == upper)
            {
                foundStep = stepList.at(i);
                foundSubStep = subStepList.at(i);
                found = true;
                break;
            }

        }
        if(!found)
        {
            //! the required time is outside the interval of available times
            //! bool SolverOuputToSolutionDataSet::getStepSubStepByTime(const QString &SolverOutputFile_path, double requiredTime, int &foundStep, int &foundSubStep)
            cout<<"error: required time greater than maximum"<<endl;
            foundStep = -1;
            foundSubStep = -1;
            return false;
        }
        else
        {
            //! return the last
            is.close();
            return true;
        }
    }
    else
    {
        //! file error
        foundStep = -1;
        foundSubStep = -1;
        if(is.is_open()) is.close();
        return false;
    }
}

//! ------------------------------
//! function: getStepSubStepBySet
//! details:
//! ------------------------------
bool  postTools::getStepSubStepBySet(const QString &SolverOutputFile_path,
                                     int requiredSet,
                                     int &foundStep,
                                     int &foundSubStep,
                                     double &foundAnalysisTime)
{
    QList<double> timeList;
    QList<int> stepList;
    QList<int> subStepList;

    int git;    //! global iteration number
    int stepNumber; //! step number
    int inc; //! increment
    int att;    //! attempt
    double ast; //! actual step time
    double time; //! actual total time
    double avf; //! average force
    double lrs; //! largest residual force
    double lid; //! largest increment of displacement
    double lcd; //! largest correction to displacement
    bool con;   //! substep converged

    std::string val;
    ifstream is(SolverOutputFile_path.toStdString());
    if(!is.fail())
    {
        int k=0;
        int subStep = 0;
        std::getline(is,val);
        for(;!is.eof();)
        {
            int N = sscanf(val.c_str(),"%d%d%d%d%lf%lf%lf%lf%lf%lf%d",&git,&stepNumber,&inc,&att,&ast,&time,&avf,&lrs,&lid,&lcd,&con);
            if(N!=11 && k==0)
            {
                //! error in line format
                foundAnalysisTime = -1.0;
                foundStep = -1;
                foundSubStep = -1;
                return false;
            }
            else
            {
                k++;
                if(con==true)
                {
                    timeList<<time;
                    stepList<<stepNumber;
                    subStep++;
                    subStepList<<subStep;
                }
            }
            std::getline(is,val);
        }

        //! process the solution data
        if(requiredSet<=stepList.length())
        {
            //! "-1" since the list of analysis times are stored in a QList
            foundAnalysisTime = timeList.at(requiredSet-1);
            foundStep = stepNumber;
            foundSubStep = subStepList.at(requiredSet-1);
            return true;

        }
        else
        {
            //! the required set does not exist
            foundAnalysisTime = -1.0;
            foundStep = -1;
            foundSubStep = -1;
            return false;
        }
    }
    else
    {
        //! error in reading file
        foundAnalysisTime = -1.0;
        foundStep = -1;
        foundSubStep = -1;
        return false;
    }
}
*/
//! ----------------------------------------------------------------
//! function: principalStressFromStressTensor
//! details:  return the principal stresses from minimum to maximum
//! ----------------------------------------------------------------
#include <algorithm>
QList<double> postTools::principalComponents(const QList<double> &sik)
{
    double s11 = sik[0];
    double s22 = sik[1];
    double s33 = sik[2];
    double s12 = sik[3];
    double s23 = sik[4];
    double s31 = sik[5];

    double I1 = s11+s22+s33;
    double I2 = s11*s22+s22*s33+s33*s11-pow(s12,2)-pow(s23,2)-pow(s31,2);
    double I3 = s11*s22*s33-s11*pow(s23,2)-s22*pow(s31,2)-s33*pow(s12,2)+2*s12*pow(s23,2);

    double x[3];

    //! ------------------------------------------------------------------------------------------------------------
    //! http://math.ivanovo.ac.ru/dalgebra/Khashin/poly/index.html
    //!
    //! In the case of three real roots function returns the number 3, the roots themselves back in x[0],x[1],x[2].
    //!
    //! Remark 1. Roots are not necessarily ordered!
    //! If two roots are match, the function returns 2 and in the array x still there are a three numbers.
    //! If the function returns 1, then x[0] is a real root and x[1]±i*x[2] is a pair of coplex-conjugated roots.
    //! Remark 2. Due to rounding errors pair of complex conjugate roots with a very small imaginary part can
    //! sometimes be a real root of multiplicity 2. For example, for equation x3 - 5x2 + 8x - 4 = 0 with roots 1,2,2
    //! we obtain roots     1.0, 2.0±i*9.6e-17. If the absolute value of the imaginary part of the root not greater
    //! than 1e-14, the SolveP3 itself replaces a pair on one valid double root, but the user must still be aware
    //! of the possibility of such a situation.
    //! ------------------------------------------------------------------------------------------------------------

    int NbRoots = SolveP3(x,-I1,I2,-I3);

    if(NbRoots==1)
    {
        x[1]=x[0];
        x[2]=x[0];
    }

    //! sort roots
    std::vector<double> vec;
    vec.push_back(x[0]);
    vec.push_back(x[1]);
    vec.push_back(x[2]);
    //! within triad, from the smallest to the largest
    std::sort(vec.begin(), vec.end());

    QList<double> values;

    for(int i=0; i<3; i++)
    {
        //cout<<" x["<<i<<"] = "<<x[i];
        //! from the largest to the smaller within the list
        values.push_back(vec.at(2-i));
    }
    cout<<endl;

    return values;
}

//! -----------------------------------
//! function: getStepSubStepByTimeDTM
//! details:
//! -----------------------------------
bool postTools::getStepSubStepByTimeDTM(QMap<double,QVector<int>> discreteTimeMap,double analysisTime,
                                         int &foundStep,
                                         int &foundSubStep)
{
    QMap<double,QVector<int>>::const_iterator mapIt;
    QVector<int> t;
    double curAnalysisTime;
    if(analysisTime = 0.0)
    {
        foundStep = discreteTimeMap.last().at(1);
        foundSubStep = discreteTimeMap.last().at(2);
        return true;
    }
    for(mapIt = discreteTimeMap.cbegin(); mapIt!= discreteTimeMap.cend(); ++mapIt)
    {
        curAnalysisTime = mapIt.key();
        t = mapIt.value();
        if (curAnalysisTime == analysisTime)
        {
            foundStep = t.at(1);
            foundSubStep = t.at(2);
            return true;
        }
    }
    return false;
}

//! ---------------------------------------------
//! function: getStepSubStepByTimeDTM
//! details:
//! ---------------------------------------------
bool postTools::getStepSubStepBySetDTM(QMap<double, QVector<int>> discreteTimeMap,
                                        int setNumber,
                                        double &analysisTime,
                                        int &foundStep,
                                        int &foundSubStep)
{
    QMap<double,QVector<int>>::const_iterator mapIt;
    int curSetNumber;
    double curTime;
    int curStep,curSubstep;

    QVector<int> t;
    for(mapIt = discreteTimeMap.cbegin(); mapIt!= discreteTimeMap.cend(); ++mapIt)
    {
        curTime = mapIt.key();
        t = mapIt.value();
        curSetNumber = t.at(0);
        curStep = t.at(1);
        curSubstep = t.at(2);
        if(setNumber == curSetNumber)
        {
            analysisTime = curTime;
            foundStep = curStep;
            foundSubStep = curSubstep;
            return true;
        }
    }
    return false;
}
