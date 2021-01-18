#ifndef POSTTOOLS_H
#define POSTTOOLS_H

//! Custom
#include "postengine.h"

//! Qt
#include <QString>
//#include "mathtools.h"

//! ---------------------------------------------
//! a collection of static functions and methods
//! ---------------------------------------------

class postTools
{

public:

    postTools();

    //! ------------------------------------------------------------------------------------------------------
    //! enter "requiredTime" (analysis time): get the corresponding pair (step/substep)
    //! current policy:
    //! 1) if requiredTime = 0.0 return the last available simulation time in the results
    //! 2) if requiredTime > time of the last available result return that
    //! 3) if requiredTime stands within an interval return floor: no average is performed at the moment
    //! 4) if requiredTime is smaller than zero (should never occur) return (step, substep) = (-1,-1) (error)
    //! ------------------------------------------------------------------------------------------------------
    /*static bool getStepSubStepByTime(const QString &SolverOutputFile_path, double requiredTime, int &foundStep, int &foundSubStep);

    //! -----------------------------------------------------------------------------------------------------
    //! entering the "requiredSet" it returns analysis time, step and substep
    //! -----------------------------------------------------------------------------------------------------
    static bool getStepSubStepBySet(const QString &SolverOutputFile_path,
                                    int requiredSet,
                                    int &foundStep,
                                    int &foundSubStep,
                                    double &foundAnalysisTime);*/

    //! ---------------------------------------------------------
    //! similar to the previous, but it uses a discrete time map
    //! ---------------------------------------------------------
    static bool getStepSubStepByTimeDTM(std::map<double, std::vector<int>> discreteTimeMap, double analysisTime,
                                         int &foundStep,
                                         int &foundSubStep);

    //! ---------------------------------------------------------
    //! similar to the previous, but it uses a discrete time map
    //! ---------------------------------------------------------
    static bool getStepSubStepBySetDTM(std::map<double, std::vector<int>> discreteTimeMap, int setNumber, double &analysisTime, int &foundStep, int &foundSubStep);

    //! ---------------------------------------------------------
    //! similar to the previous, but it uses a discrete time map
    //! ---------------------------------------------------------
    static bool postTools::getSetBySubStepByStepDTM(std::map<double, std::vector<int>> discreteTimeMap,
                                           int &setNumber,
                                           double &analysisTime,
                                           int step,
                                           int subStep);

    static void principalComponents(double *sik, double *values);
};

#endif // POSTTOOLS_H
