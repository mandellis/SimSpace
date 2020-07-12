#ifndef CCXCONSOLETOFILE_H
#define CCXCONSOLETOFILE_H

//! ---
//! Qt
//! ---
#include <QString>

//! ----------------
//! custom includes
//! ----------------
#include "solutioninfo.h"
#include "runterminationdata.h"
#include "qprogressindicator.h"

//! ----
//! C++
//! ----
#include <iostream>
#include <fstream>
using namespace std;

class CCXconsoleToFile
{
public:

    //! constructor
    CCXconsoleToFile();

    //! perform
    static runTerminationData perform(QString myTargetFileName,
                                      QString mySourceFileName,
                                      int analysisType,
                                      QList<solutionInfo> &listSolInfo,
                                      bool &simulationError,
                                      QProgressIndicator *aProgressIndicator=Q_NULLPTR);
};

#endif // CCXCONSOLETOFILE_H
