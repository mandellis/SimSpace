//! ----------------
//! custom includes
//! ----------------
#include "inputfilegenerator.h"
#include <qprogressindicator.h>
#include <qprogressevent.h>

//! ---
//! Qt
//! ---
#include <QString>

//! ----
//! CCX
//! ----
#include "qextendedstandarditem.h"
#include <simulationdatabase.h>
#include "writesolverfileclass.h"

//! ----
//! C++
//! ----
#include <iostream>
using namespace std;

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
inputFileGenerator::inputFileGenerator(QObject *parent):QThread(parent)
{
    ;
}

//! -------------------------------
//! function: setProgressIndicator
//! details:
//! -------------------------------
void inputFileGenerator::setProgressIndicator(QProgressIndicator *aProgressIndicator)
{
    myProgressIndicator = aProgressIndicator;
}

//! ------------------------
//! function: setParameters
//! details:
//! ------------------------
void inputFileGenerator::setParameters(std::vector<void*> parameters)
{
    myParameters = parameters;
    myInputFileName = QString(*(QString*)(myParameters[2]));
}

//! --------------
//! function: run
//! details:
//! --------------
void inputFileGenerator::run()
{
    cout<<"inputFileGenerator::run()->____run called. Thread: "<<QThread::currentThreadId()<<"____"<<endl;
    bool isDone = false;
    int solverTarget = 0;
    switch(solverTarget)
    {
    case 0:
    {
        //! ---------------------
        //! write CCX input file
        //! ---------------------
        //isDone = this->writeCCX();
        cout<<"inputFileGenerator::writeCCX()->____function called. Thread: "<<QThread::currentThreadId()<<"____"<<endl;
        if(myParameters.size() != 3) emit inputFileWritten(false);

        //! -----------------------
        //! unpack parameters:
        //! - simulation data base
        //! - simulation root
        //! - file name
        //! -----------------------
        simulationDataBase *sDB = static_cast<simulationDataBase*>(myParameters[0]);
        QExtendedStandardItem *simulationRootItem = static_cast<QExtendedStandardItem*>(myParameters[1]);
        writeSolverFileClass CCXInputGenerator(sDB,simulationRootItem);

        CCXInputGenerator.setProgressIndicator(myProgressIndicator);
        CCXInputGenerator.setName(myInputFileName);
        bool isDone = CCXInputGenerator.perform();
        emit inputFileWritten(isDone);
    }
        break;

    case 1:
    {
        //! ------------------------
        //! put here another solver
        //! ------------------------
        isDone = false;
    }
        break;

    default:
    {
        //! --------------
        //! other solvers
        //! --------------
        isDone = false;
    }
        break;
    }
    emit inputFileWritten(isDone);
}

//! -------------------
//! function: writeCCX
//! details:
//! -------------------
bool inputFileGenerator::writeCCX()
{
    cout<<"inputFileGenerator::writeCCX()->____function called. Thread: "<<QThread::currentThreadId()<<"____"<<endl;
    if(myParameters.size() != 3) return false;

    //! -----------------------
    //! unpack parameters:
    //! - simulation data base
    //! - simulation root
    //! - file name
    //! -----------------------
    simulationDataBase *sDB = static_cast<simulationDataBase*>(myParameters[0]);
    QExtendedStandardItem *simulationRootItem = static_cast<QExtendedStandardItem*>(myParameters[1]);
    writeSolverFileClass CCXInputGenerator(sDB,simulationRootItem);

    CCXInputGenerator.setProgressIndicator(myProgressIndicator);
    CCXInputGenerator.setName(myInputFileName);
    bool isDone = CCXInputGenerator.perform();
    return isDone;
}
