#ifndef SOLUTIONINFO_H
#define SOLUTIONINFO_H

//! ---
//! Qt
//! ---
#include <QVariant>
#include <QMetaType>
#include <iostream>
#include <fstream>

struct solutionInfo
{
    //! ---------------------------------------------------------
    //! 1) global iteration number
    //! 2) step number
    //! 3) increment
    //! 4) attempt
    //! 5) actual step time
    //! 6) actual total time
    //! 7) average (force/temp)
    //! 8) largest residual (force/flux)
    //! 9) largest increment of DOF (displacement/temperature)
    //! 10) largest correction to DOF (displacement/temperature)
    //! 11) substep converged
    //! ---------------------------------------------------------

public:

    //! ------------
    //! constructor
    //! ------------
    solutionInfo(){;}

    //! -----------------
    //! copy constructor
    //! -----------------
    solutionInfo(const solutionInfo &other)
    {
        globalIterationNb=other.globalIterationNb;
        stepNumber=other.stepNumber;
        increment=other.increment;
        attempt=other.attempt;
        steptime=other.steptime;
        time=other.time;
        average=other.average;
        largestResidual=other.largestResidual;
        largestDOFIncrement=other.largestDOFIncrement;
        largestDOFCorrection=other.largestDOFCorrection;
        subStepConverged=other.subStepConverged;
    }

    //! ------------
    //! operator ==
    //! ------------
    bool operator==(const solutionInfo& rhs)
    {
        return (globalIterationNb == rhs.globalIterationNb &&
        stepNumber == rhs.stepNumber &&
        increment == rhs.increment &&
        attempt == rhs.attempt &&
        steptime == rhs.steptime &&
        time == rhs.time &&
        average == rhs.average &&
        largestResidual == rhs.largestResidual &&
        largestDOFIncrement == rhs.largestDOFIncrement &&
        largestDOFCorrection == rhs.largestDOFCorrection &&
        subStepConverged == rhs.subStepConverged);
    }

    //! -----------
    //! operator =
    //! -----------
    solutionInfo operator=(const solutionInfo &other)
    {
        globalIterationNb=other.globalIterationNb;
        stepNumber=other.stepNumber;
        increment=other.increment;
        attempt=other.attempt;
        steptime=other.steptime;
        time=other.time;
        average=other.average;
        largestResidual=other.largestResidual;
        largestDOFIncrement=other.largestDOFIncrement;
        largestDOFCorrection=other.largestDOFCorrection;
        subStepConverged=other.subStepConverged;
    }

    int globalIterationNb;
    int stepNumber;
    int increment;
    int attempt;
    double steptime;
    double time;
    double average;
    double largestResidual;
    double largestDOFIncrement;
    double largestDOFCorrection;
    bool subStepConverged;

    //! -----------------
    //! function: toList
    //! details:  helper
    //! -----------------
    QList<QVariant> toList() const
    {
        QList<QVariant> RV;
        QVariant data;

        data.setValue(globalIterationNb); RV<<data;         //! [0]
        data.setValue(stepNumber); RV<<data;                //! [1]
        data.setValue(increment); RV<<data;                 //! [2]
        data.setValue(attempt); RV<<data;                   //! [3]
        data.setValue(steptime); RV<<data;                  //! [4]
        data.setValue(time); RV<<data;                      //! [5]
        data.setValue(average); RV<<data;                   //! [6]
        data.setValue(largestResidual); RV<<data;           //! [7]
        data.setValue(largestDOFIncrement); RV<<data;       //! [8]
        data.setValue(largestDOFCorrection); RV<<data;      //! [9]
        data.setValue(subStepConverged); RV<<data;          //! [10]

        return RV;
    }

    //! ------
    //! write
    //! ------
    void write(std::ofstream &out) const
    {
        out<<globalIterationNb<<std::endl;
        out<<stepNumber<<std::endl;
        out<<increment<<std::endl;
        out<<attempt<<std::endl;
        out<<steptime<<std::endl;
        out<<time<<std::endl;
        out<<average<<std::endl;
        out<<largestResidual<<std::endl;
        out<<largestDOFIncrement<<std::endl;
        out<<largestDOFCorrection<<std::endl;
        if(subStepConverged==false) out<<0<<std::endl;
        else out<<1<<std::endl;
    }

    //! -----
    //! read
    //! -----
    void read(std::ifstream &in)
    {
        in>>globalIterationNb;
        in>>stepNumber;
        in>>increment;
        in>>attempt;
        in>>steptime;
        in>>time;
        in>>average;
        in>>largestResidual;
        in>>largestDOFIncrement;
        in>>largestDOFCorrection;
        int val;
        in>>val;
        if(val==0) subStepConverged = false;
        else subStepConverged = true;
    }
};

Q_DECLARE_METATYPE(solutionInfo)

#endif // SOLUTIONINFO_H
