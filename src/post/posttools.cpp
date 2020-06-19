//! ----------------
//! custom includes
//! ----------------
#include "posttools.h"
#include "cubicequation.h"
#include "postengine.h"

//! ---
//! Qt
//! ---
#include <QList>

//! ----
//! C++
//! ----
#include <algorithm>
#include <iostream>
#include <fstream>
using namespace std;

postTools::postTools()
{
    ;
}

//! ------------------------------------------------------------------
//! function: principalComponents
//! details:  return the principal components from minimum to maximum
//! ------------------------------------------------------------------
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
    //cout<<endl;

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
    if(analysisTime == 0.0)
    {
        foundStep = discreteTimeMap.last().at(1);
        foundSubStep = discreteTimeMap.last().at(2);
        return true;
    }
    for(mapIt = discreteTimeMap.cbegin(); mapIt!= discreteTimeMap.cend(); ++mapIt)
    {
        curAnalysisTime = mapIt.key();
        //cout<<"postTools::getStepSubStepbyTimeDTM->____curAnalysisTime"<<curAnalysisTime<<", analysisTime "<<analysisTime<<endl;
        t = mapIt.value();
        if (curAnalysisTime == analysisTime)
        {
            foundStep = t.at(1);
            foundSubStep = t.at(2);
            cout<<"postTools::getStepSubStepbyTimeDTM->____step"<<foundStep<<", substep "<<foundSubStep<<endl;

            return true;
        }
    }
    return false;
}

//! ---------------------------------------------
//! function: getStepSubStepBySetDTM
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
