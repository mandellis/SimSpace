//! ----------------
//! custom includes
//! ----------------
#include "ccxconsoletofile.h"
#include "src/utils/tools.h"
#include "solutioninfo.h"
#include "src/utils/ccout.h"
#include "qsimulationstatusevent.h"
#include "qprogressevent.h"
#include "src/utils/global.h"

//! ---
//! Qt
//! ---
#include <QApplication>

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
CCXconsoleToFile::CCXconsoleToFile()
{
    ;
}

//! -----------------------------------------------------------------
//! function: perform
//! details:  return values
//!         - final read time (double)
//!         - solution info (QList<solutionInfo)
//!         - simulation error (bool)
//!           If a simulation error has occurred the final read time
//!           is the last converged substep
//! -----------------------------------------------------------------
runTerminationData CCXconsoleToFile::perform(QString myTargetFileName,
                                             QString mySourceFileName,
                                             int analysisType,
                                             QList<solutionInfo> &listSolInfo,
                                             bool &simulationError,
                                             QProgressIndicator *aProgressIndicator)
{
    cout<<"CCXconsoleToFile::perform()->____function called____"<<endl;

    //! ---------------
    //! open the files
    //! ---------------
    ifstream inputFile(mySourceFileName.toStdString());
    if(inputFile.is_open()==false) return runTerminationData();

    ofstream outFile;

    bool writeOutputFile = (myTargetFileName.isEmpty() || myTargetFileName.isNull())? false:true;
    if(writeOutputFile) outFile.open(myTargetFileName.toStdString());

    //! element of the list
    solutionInfo aSolutionInfo;

    //! simulation error
    double lastConvergedTime = 0.0;
    int lastConvergedStepNumber = 0;
    int lastConvergedSubStepNumber = 0;

    double lastConvergedTime_Old = 0.0;
    int lastConvergedStepNumber_Old = 0;
    int lastConvergedSubStepNumber_Old = 0;

    std::string val;
    std::getline(inputFile,val);

    QString sv = QString::fromStdString(val);

    int globalIterationNb = 0;
    for(;;)
    {
        //! check termination
        if(inputFile.eof()==true)
        {
            //cerr<<"end of file"<<endl;
            break;
        }

        lastConvergedTime_Old = lastConvergedTime;
        lastConvergedStepNumber_Old = lastConvergedStepNumber;
        lastConvergedSubStepNumber_Old = lastConvergedSubStepNumber;

        int stepNumber;
        if(sv.startsWith(" STEP"))
        {
            //! line " STEP            <int>" found
            sscanf(sv.toStdString().c_str()," STEP            %d",&stepNumber);
        }

        bool isConverged = false;
        bool endOfFile = false;

        if(sv.contains("increment") && sv.contains("attempt"))
        {
            int attempt;
            int increment;
            sscanf(sv.toStdString().c_str()," increment %d attempt %d ",&increment,&attempt);

            //! -------------------------
            //! increment size  <double>
            //! -------------------------
            std::getline(inputFile,val);

            double incrementSize;
            sscanf(sv.toStdString().c_str()," increment size= %lf",&incrementSize);

            //! -----------------------------------
            //! sum of previous increment <double>
            //! -----------------------------------
            std::getline(inputFile,val);

            double incrementsSum;
            sscanf(sv.toStdString().c_str()," sum of previous increments=%lf",&incrementsSum);

            //! --------------------------
            //! actual step time <double>
            //! --------------------------
            std::getline(inputFile,val);

            //!cerr<<val<<endl;
            double actualStepTime;
            sscanf(val.c_str()," actual step time=%lf",&actualStepTime);

            //! ---------------------------
            //! actual total time <double>
            //! ---------------------------
            std::getline(inputFile,val);

            //!cerr<<val<<endl;
            double actualTotalTime;
            sscanf(val.c_str()," actual total time=%lf",&actualTotalTime);

            std::getline(inputFile,val);
            sv = QString::fromStdString(val);

            for(;;)
            {
                if(sv.contains("iteration"))
                {
                    //!cerr<<val<<endl;                                                                //! iteration <n>
                    std::getline(inputFile,val);
                    sv = QString::fromStdString(val);

                    for(;;)
                    {
                        if(sv.contains("average"))
                        {
                            //! --------------------------------------------
                            //! 1) global iteration number
                            //! 2) step number
                            //! 3) increment
                            //! 4) attempt
                            //! 5) actual step time
                            //! 6) actual total time
                            //! 7) average (force/flux)
                            //! 8) largest residual force/flux
                            //! 9) largest increment of displacement/temp
                            //! 10) largest correction to displacement/temp
                            //! 11) substep converged
                            //! --------------------------------------------
                            if(writeOutputFile)
                            {
                                outFile<<++globalIterationNb<<"\t"; //! [1]
                                outFile<<stepNumber<<"\t";          //! [2]
                                outFile<<increment<<"\t";           //! [3]
                                outFile<<attempt<<"\t";             //! [4]
                                outFile<<actualStepTime<<"\t";      //! [5]
                                outFile<<actualTotalTime<<"\t";     //! [6]
                            }

                            aSolutionInfo.globalIterationNb = globalIterationNb;
                            aSolutionInfo.time = actualTotalTime;

                            //!cerr<<val<<endl;                                             //! average force/flux <double>
                            double average;
                            if(analysisType==0) sscanf(val.c_str()," average force= %lf",&average);
                            if(analysisType==1) sscanf(val.c_str()," average flux= %lf",&average);

                            if(writeOutputFile) outFile<<average<<"\t";    //! [7] write the largest residual force <double>
                            aSolutionInfo.average = average;

                            std::getline(inputFile,val);
                            //!cerr<<val<<endl;                                             //! time avg. force/flux <double>

                            std::getline(inputFile,val);
                            //!cerr<<val<<endl;                                            //! largest residual force/flux <double>

                            int node, dof;
                            double largestRes;
                            if(analysisType==0) sscanf(val.c_str()," largest residual force= %lf in node %d and dof %d",&largestRes,&node,&dof);
                            if(analysisType==1) sscanf(val.c_str()," largest residual flux= %lf in node %d and dof %d",&largestRes,&node,&dof);

                            if(writeOutputFile) outFile<<largestRes<<"\t";       //! [8] write the largest residual force/flux <double>
                            aSolutionInfo.largestResidual = largestRes;

                            std::getline(inputFile,val);
                            //!cerr<<val<<endl;                                             //! largest increment of DOF (displacement/temperature)

                            double largestDOFlIncrement;
                            if(analysisType==0) sscanf(val.c_str()," largest increment of disp=%lf",&largestDOFlIncrement);
                            if(analysisType==1) sscanf(val.c_str()," largest increment of temp=%lf",&largestDOFlIncrement);

                            if(writeOutputFile) outFile<<largestDOFlIncrement<<"\t";       //! [9] largest increment of DOF (displacement/temperature)
                            aSolutionInfo.largestDOFIncrement = largestDOFlIncrement;

                            std::getline(inputFile,val);
                            //!cerr<<val<<endl;                                            //! largest correction to disp

                            double largestDOFCorrection;
                            if(analysisType==0) sscanf(val.c_str()," largest correction to disp= %lf in node %d and dof %d",&largestDOFCorrection,&node,&dof);
                            if(analysisType==1) sscanf(val.c_str()," largest correction to temp= %lf in node %d and dof %d",&largestDOFCorrection,&node,&dof);

                            if(writeOutputFile) outFile<<largestDOFCorrection<<"\t";       //! [10]    largest correction to DOF (displacement/temperature)
                            aSolutionInfo.largestDOFCorrection = largestDOFCorrection;

                            //! --------------------------------------------
                            //! search for "convergence" or "no convergence"
                            //! --------------------------------------------
                            std::getline(inputFile,val);
                            sv = QString::fromStdString(val);

                            for(;;)
                            {
                                if(sv.contains("convergence")
                                        || sv.contains("no convergence")
                                        || sv.contains("divergence allowed"))
                                {
                                    //!cerr<<val<<endl;
                                    if(sv.contains("no convergence") || sv.contains("divergence allowed"))
                                    {

                                        isConverged=false;
                                        if(writeOutputFile) outFile<<0<<"\n";   //! [11]    substep converged
                                        aSolutionInfo.subStepConverged = false;
                                        listSolInfo<<aSolutionInfo;
                                        if(sv.contains("divergence allowed"))
                                        {
                                            //cout<<"FOUND \"divergence allowed\""<<endl;
                                            //! post an event "Divergence allowed"
                                            if(aProgressIndicator!=Q_NULLPTR)
                                            {
                                                QProgressEvent *e = new QProgressEvent(QProgressEvent_Update,0,100,0,"Divergence allowed",
                                                                                       QProgressEvent_None,-1,-1,-1,"Running CCX");
                                                QApplication::postEvent(aProgressIndicator,e);
                                                QApplication::processEvents();
                                            }
                                            //if(targetWidgetSimulationManager!=NULL)
                                            //{
                                            //    QSimulationStatusEvent *event = new QSimulationStatusEvent(stepNumber,increment,actualTotalTime,0,"Divergence allowed");
                                            //    QApplication::postEvent(targetWidgetSimulationManager,static_cast<QEvent*>(event));
                                            //    QApplication::processEvents();
                                            //}
                                        }
                                        else if(sv.contains("no convergence"))
                                        {
                                            if(aProgressIndicator!=Q_NULLPTR)
                                            {
                                                QProgressEvent *e = new QProgressEvent(QProgressEvent_Update,0,100,0,"No convergence",
                                                                                       QProgressEvent_None,-1,-1,-1,"Running CCX");
                                                QApplication::postEvent(aProgressIndicator,e);
                                                QApplication::processEvents();
                                            }
                                            /*
                                            //! post an event "No convergence"
                                            if(targetWidgetSimulationManager!=NULL)
                                            {
                                                QSimulationStatusEvent *event = new QSimulationStatusEvent(stepNumber,increment,actualTotalTime,0,"No convergence");
                                                QApplication::postEvent(targetWidgetSimulationManager,static_cast<QEvent*>(event));
                                                QApplication::processEvents();
                                            }
                                            */
                                        }
                                    }
                                    else
                                    {
                                        isConverged=true;
                                        if(writeOutputFile) outFile<<1<<"\n";   //! [11]    substep converged
                                        aSolutionInfo.subStepConverged = true;
                                        listSolInfo<<aSolutionInfo;

                                        lastConvergedTime = actualTotalTime;
                                        lastConvergedStepNumber = stepNumber;
                                        lastConvergedSubStepNumber = increment;

                                        //! post an event "Sub step converged"
                                        if(aProgressIndicator!=Q_NULLPTR)
                                        {
                                            QProgressEvent *e = new QProgressEvent(QProgressEvent_Update,0,100,0,"Converged",
                                                                                   QProgressEvent_None,-1,-1,-1,"Running CCX");
                                            QApplication::postEvent(aProgressIndicator,e);
                                            QApplication::processEvents();
                                        }
                                    }
                                    break;
                                }
                                std::getline(inputFile,val);
                                endOfFile = inputFile.eof();
                                if(endOfFile==true) break;
                                sv = QString::fromStdString(val);
                            }
                        }
                        if(isConverged==true) break;
                        endOfFile = inputFile.eof();
                        if(endOfFile==true) break;

                        std::getline(inputFile,val);
                        sv = QString::fromStdString(val);
                        if(sv.contains("*ERROR"))
                        {
                            cout<<"FOUND ERROR"<<endl;
                            simulationError = true;
                        }
                    }
                }
                if(isConverged==true) break;
                endOfFile = inputFile.eof();
                if(endOfFile==true)
                {
                    break;
                }

                if(simulationError==true) break;

                std::getline(inputFile,val);
                sv = QString::fromStdString(val);
            }
        }
        std::getline(inputFile,val);
        endOfFile = inputFile.eof();
        if(endOfFile==true) break;
        sv = QString::fromStdString(val);

        if(simulationError==true) break;
    }
    inputFile.close();
    if(writeOutputFile) outFile.close();

    runTerminationData rtd;
    double lastAvailableTime, lastAvailableStepNumber, lastAvailableSubStepNumber;

    if(simulationError==true)
    {
        cout<<"SOLUTION DID NOT CONVERGE"<<endl;
        if(listSolInfo.isEmpty())
        {
            //! the solver has failed at the very beginning for some reason
            //! error at the very beginning
            lastAvailableStepNumber = 0;
            lastAvailableSubStepNumber = 0;
            lastAvailableTime = 0;

            //! return value
            rtd.lastAvailableTime = lastAvailableTime;
            rtd.lastAvailableStep = lastAvailableStepNumber;
            rtd.lastAvailableSubStep = lastAvailableSubStepNumber;

            if(aProgressIndicator!=Q_NULLPTR)
            {
                QProgressEvent *e = new QProgressEvent(QProgressEvent_Reset,-1,-1,-1,"Solver failed to start",
                                                       QProgressEvent_None,-1,-1,-1,"");
                QApplication::postEvent(aProgressIndicator,e);
                QApplication::processEvents();
            }
        }
        else
        {
            int N = listSolInfo.length();
            if(N==1)
            {
                //! the solver has failed at the first substep
                //! error at the very first substep number
                lastAvailableStepNumber = 0;
                lastAvailableSubStepNumber = 0;
                lastAvailableTime = 0;

                //! return value
                rtd.lastAvailableTime = lastAvailableTime;
                rtd.lastAvailableStep = lastAvailableStepNumber;
                rtd.lastAvailableSubStep = lastAvailableSubStepNumber;

                if(aProgressIndicator!=Q_NULLPTR)
                {
                    QProgressEvent *e = new QProgressEvent(QProgressEvent_Reset,-1,-1,-1,"Error at first time step",
                                                           QProgressEvent_None,-1,-1,-1,"");
                    QApplication::postEvent(aProgressIndicator,e);
                    QApplication::processEvents();
                }
            }
            else if(N>1)
            {
                //! at least one time substep has been calculated
                //! error at a subsequent time step
                lastAvailableTime = lastConvergedTime_Old;
                lastAvailableStepNumber = lastConvergedStepNumber_Old;
                lastAvailableSubStepNumber = lastConvergedSubStepNumber_Old;

                if(aProgressIndicator!=Q_NULLPTR)
                {
                    QProgressEvent *e = new QProgressEvent(QProgressEvent_Reset,-1,-1,-1,"Solver did not converge",
                                                           QProgressEvent_None,-1,-1,-1,"");
                    QApplication::postEvent(aProgressIndicator,e);
                    QApplication::processEvents();
                }

                cout<<"SOLVER FAILED AT TIME: "<<lastAvailableTime<<endl;
                cout<<"LAST AVAILABLE STEP NR. "<<lastAvailableStepNumber<<endl;
                cout<<"LAST AVAILABLE SUBSTEP NR. "<<lastAvailableSubStepNumber<<endl;

                //! return value
                rtd.lastAvailableTime = lastAvailableTime;
                rtd.lastAvailableStep = lastAvailableStepNumber;
                rtd.lastAvailableSubStep = lastAvailableSubStepNumber;
            }
        }
    }
    else
    {
        //! ----------------
        //! return value
        //! ----------------
        rtd.lastAvailableTime = lastConvergedTime;
        rtd.lastAvailableStep = lastConvergedStepNumber;
        rtd.lastAvailableSubStep = lastConvergedSubStepNumber;
    }
    cout<<"------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"- SimulationManager::configureAndStartPostEngine()->____DATA AVAILABILITY AT THE END OF THE RUN____"<<endl;
    cout<<"- Last available time: "<<rtd.lastAvailableTime<<endl;
    cout<<"- Last available step: "<<rtd.lastAvailableStep<<endl;
    cout<<"- Last available substep: "<<rtd.lastAvailableSubStep<<endl;
    cout<<"------------------------------------------------------------------------------------------------------"<<endl;

    //! ----------------------------------------------
    //! update the progress bar of the solver manager
    //! ----------------------------------------------
    if(!listSolInfo.isEmpty())
    {
        if(listSolInfo.last().subStepConverged==true)
        {
            if(aProgressIndicator!=Q_NULLPTR)
            {
                int theProgress = int(100*listSolInfo.last().time);
                QProgressEvent *e = new QProgressEvent(QProgressEvent_Update,-1,-1,theProgress,"Substep converged",
                                                       QProgressEvent_None,-1,-1,-1,"Running CCX");
                QApplication::postEvent(aProgressIndicator,e);
                QApplication::processEvents();
            }
        }
    }
    return rtd;
}
