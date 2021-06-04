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
//! OF
//! ----
#include "ofwrite.h"

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
    cout<<"inputFileGenerator::inputFileGeneratorn()->____function called. Thread: "<<QThread::currentThreadId()<<"____"<<endl;
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

    QExtendedStandardItem *simulationRootItem = static_cast<QExtendedStandardItem*>(myParameters[1]);
    SimulationNodeClass *curNode = simulationRootItem->data(Qt::UserRole).value<SimulationNodeClass*>();
    int solverTarget = 0;
    if(curNode->getType()==SimulationNodeClass::nodeType_CFDAnalysis) solverTarget = 1;

    switch(solverTarget)
    {
    case 0:
    {
        //! ---------------------
        //! write CCX input file
        //! ---------------------
        isDone = this->writeCCX();
        if(myParameters.size() != 3) emit inputFileWritten(false);

        /*
         * run() should be executed in a thread different from the main one.
         * (that's what we would like to do). From documentation it appears
         * that a SLOT, defined in the class and called in within run() cannot
         * be executed in the new thread, but will stand in the old. It appears
         * also that, for obtaining this result, the SLOT should be defined in
         * an independent object (a worker). Actually even if writeCCX() is defined
         * outside run, its thread Id is different from the thread Id of the main
         * thread. Output example:
         *
         * inputFileGenerator::inputFileGeneratorn()->____function called. Thread: 0000000000001D0C____
         * inputFileGenerator::writeCCX()->____function called. Thread: 0000000000004EE8____
         * /
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
        */
    }
        break;

    case 1:
    {
        //! ------------------------
        //! put here another solver
        //! ------------------------
        isDone = this->writeOF();
        if(myParameters.size() != 3) emit inputFileWritten(false);
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

bool inputFileGenerator::writeOF()
{
    cout<<"inputFileGenerator::writeOF()->____function called. Thread: "<<QThread::currentThreadId()<<"____"<<endl;
    if(myParameters.size() != 3) return false;

    //! -----------------------
    //! unpack parameters:
    //! - simulation data base
    //! - simulation root
    //! - file name
    //! -----------------------
    simulationDataBase *sDB = static_cast<simulationDataBase*>(myParameters[0]);
    QExtendedStandardItem *simulationRootItem = static_cast<QExtendedStandardItem*>(myParameters[1]);
    ofwrite OFInputGenerator(sDB,simulationRootItem);

    OFInputGenerator.setProgressIndicator(myProgressIndicator);
    OFInputGenerator.setName(myInputFileName);
    bool isDone = OFInputGenerator.perform();
    return isDone;
}
