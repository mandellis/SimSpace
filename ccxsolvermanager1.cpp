//! ----------------
//! custom includes
//! ----------------
#include "ccxsolvermanager1.h"
#include "qprogressindicator.h"
#include <qccxsolvermessageevent.h>
#include <qconsoleevent.h>
#include <ccxsolvermessage.h>

//! ---
//! Qt
//! ---
#include <QProcess>
#include <QApplication>

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
CCXSolverManager1::CCXSolverManager1(QObject *parent):QThread(parent)
{
    myNbProcessors = 2;

    //! ---------------------------------------------------------------
    //! retrieve the solution monitor panel and the simulation manager
    //! ---------------------------------------------------------------
    QList<QWidget*> widgets = QApplication::allWidgets();
    for(int i=0; i<widgets.length(); i++)
    {
        if(widgets.at(i)->objectName()=="worksheetViewer")
        {
            mySolutionMonitorPanel = widgets.at(i);
            break;
        }
    }
    for(int i=0; i<widgets.length(); i++)
    {
        if(widgets.at(i)->objectName()=="simmanager")
        {
            mySimulationManager = widgets.at(i);
            break;
        }
    }
}

//! -----------------------
//! function: setInputFile
//! details:
//! -----------------------
void CCXSolverManager1::setInputFile(const QString &inputFile)
{
    myInputFileName = inputFile;
}

//! --------------------------
//! function: setNbProcessors
//! details:
//! --------------------------
void CCXSolverManager1::setNbProcessors(int NbProcessors)
{
    if(NbProcessors<1) NbProcessors = 1;
    myNbProcessors = NbProcessors;
}

//! -------------------------------
//! function: setProgressIndicator
//! details:
//! -------------------------------
void CCXSolverManager1::setProgressIndicator(QProgressIndicator *aProgressIndicator)
{
    myProgressIndicator = aProgressIndicator;
}

//! --------------
//! function: run
//! details:
//! --------------
void CCXSolverManager1::run()
{
    //! -----------------
    //! create a process
    //! -----------------
    myProcess = new QProcess(this);
    connect(myProcess,SIGNAL(readyReadStandardOutput()),this,SLOT(retrieveSolverData()));

    //! ---------------------------------------
    //! calculix solver: program name and args
    //! ---------------------------------------
    QString program;
    bool CCX_uses_SPOOLES = false;
    if(CCX_uses_SPOOLES)
    {
        program =  QString::fromLatin1(std::getenv("CCX_MT_SPOOLES"));
        if(program.isEmpty()) { ; }
        program.append("/ccx_2.16_spooles.exe");
    }
    else
    {
        program = QString::fromLatin1(std::getenv("CCX_MT_PARDISO"));
        if(program.isEmpty()) { ; }
        program.append("/ccx_215Pardiso.exe");
    }

    //! --------------------------
    //! "default" input file name
    //! --------------------------
    QString inputFileName("input");

    //! -----------------------------------------
    //! solution data directory: "initial" value
    //! -----------------------------------------
    QString solutionDataDir = myInputFileName;
    QString tmp = solutionDataDir.split("/").last();

    //! -------------------------------------
    //! solution data directory: final value
    //! -------------------------------------
    solutionDataDir.chop(tmp.length()+1);

    QList<QString> CCXarguments;
    CCXarguments<<inputFileName;

    //! ----------------
    //! number of cores
    //! ----------------
    QProcessEnvironment processEnv;
    processEnv.insert("omp_num_threads",QString("%1").arg(myNbProcessors));
    myProcess->setWorkingDirectory(solutionDataDir);
    myProcess->setProcessEnvironment(processEnv);
    myProcess->start(program,CCXarguments);
    myProcess->waitForFinished(-1);

    emit CCXRunFinished();
}

//! --------------------
//! function: copy_file
//! details:
//! --------------------
void CCXSolverManager1::copy_file(const char* srce_file, const char* dest_file)
{
    std::ifstream srce(srce_file, std::ios::binary);
    std::ofstream dest(dest_file, std::ios::binary);
    dest<<srce.rdbuf();
}

//! -----------------------------------------------------------
//! function: retrieveSolverData
//! details:  1) send the CCX solver messages to the worksheet
//!           2) append the CCX output to rawSolverData.txt
//! -----------------------------------------------------------
void CCXSolverManager1::retrieveSolverData()
{
    //cout<<"SolutionWorker::retrieveSolverData()->____function called on input file: "<<myInputFileName.toStdString()<<"____"<<endl;

    QString textBlock = myProcess->readAllStandardOutput();
    myCCXLog.append(textBlock);

    //! -------------------------------------------------
    //! send the solver output to the simulation manager
    //! -------------------------------------------------
    if(mySimulationManager!=Q_NULLPTR)
    {
        CCXSolverMessage msg(myCCXLog);
        QCCXSolverMessageEvent *e = new QCCXSolverMessageEvent(msg);
        QApplication::postEvent(mySimulationManager,e);
        QApplication::processEvents();
    }
    //! -------------------------------------------------
    //! send the solver messages to the worksheet viewer
    //! -------------------------------------------------
    if(mySolutionMonitorPanel!=Q_NULLPTR)
    {
        QConsoleEvent *e = new QConsoleEvent(textBlock);
        QApplication::postEvent(mySolutionMonitorPanel,e);
    }

    //! -----------------------------------------------
    //! files RawSolverData.txt RawSolverData_copy.txt
    //! -----------------------------------------------
    QString rawSolverDataFileName = myInputFileName;
    QString rawSolverDataFileName_copy = myInputFileName;

    int nbchar = rawSolverDataFileName.split("/").last().length();

    rawSolverDataFileName.chop(nbchar);
    rawSolverDataFileName_copy.chop(nbchar);
    rawSolverDataFileName.append("RawSolverOutput.txt");
    rawSolverDataFileName_copy.append("RawSolverOutput_copy.txt");

    QString lockFileName = myInputFileName;
    lockFileName.chop(nbchar);
    lockFileName.append("lock");

    //! -------------------------------------------------------------------------------
    //! The solver output is written into the RawSolverOutput.txt, which contains
    //! the convergence info. That file should be read by the convergence data
    //! viewer (convergenceDataViewer.cpp).
    //! The "write" operation, by this class, is triggered by the solver,
    //! while the "read" operation, by the convergence data viewer, triggered by
    //! the timeout() signal of a QTimer belonging to the "simulationManager.cpp".
    //! In such a way a multiple access to the same file, in write and read mode,
    //! could verify, and that should be avoided. In order to avoid this occurrence
    //! the following technique is used.
    //!
    //! [1] The solutionWorker emit readyReadStandardOuput() => the solution
    //!     worker append the block of text to the RawSolverOutput, using
    //!     the "retrieveSolverData" slot
    //! [2] If a "lock" file does not exists, a lock file is created,
    //!     and a copy of the RawSolverData.txt is made: "RawSolverOuput_copy.txt"
    //! [3] "RawSolverOutput_copy.txt" is closed, and the "lock" file is removed
    //!
    //! At the same time, in an asynchronous ways, the convergence data viewer
    //! tries to read the "RawSolverOutput_copy.txt", by checking the presence
    //! of the "lock" file. If the "lock" file exixts, it does not process
    //! the solver output (since it is in use by the solutionWorker). Otherwise,
    //! if the lock file does exists, the convergence data viewer creates a new
    //! one, then opens "RawSolverOutput_copy.txt" and processes the data, then closes
    //! both the files, finally removes the "lock".
    //! While this is occurring the "RawSolverData.txt" can be indipendently
    //! updated by the panel is updating (convergence info are not lost):
    //! only the copy process can not be performed. This has no consequences
    //! because panel update will be performed, in the worst case, at the end
    //! of the solution process.
    //! -------------------------------------------------------------------------------

    //! ----------------------------------
    //! append the text block to the file
    //! ----------------------------------
    cout<<"SolutionWorker::retrieveSolverData()->____writing a block of text into the raw solver output file____"<<endl;
    FILE *rawSolverData = fopen(rawSolverDataFileName.toStdString().c_str(),"a");
    fprintf(rawSolverData, textBlock.toStdString().c_str());
    fclose(rawSolverData);

    //! --------------------------------------------------------------------------
    //! the solver output can be written into the RawSolverOutput.txt file
    //! if that file is not opened by for processing by the simulationManager.
    //! The file status depends on the presence of the lock file.
    //! If the lock file does not exists create it, the write the CCX solver log,
    //! finally close and remove the lock
    //! --------------------------------------------------------------------------
    if(!QFile::exists(lockFileName))
    {
        cout<<"SolutionWorker::retrieveSolverData()->____lock file created____"<<endl;
        FILE *lockFile = fopen(lockFileName.toStdString().c_str(),"w");
        copy_file(rawSolverDataFileName.toStdString().c_str(), rawSolverDataFileName_copy.toStdString().c_str());
        fclose(lockFile);
        bool lockRemoved = QFile::remove(lockFileName);
        if(lockRemoved) cout<<"SolutionWorker::retrieveSolverData()->____lock file removed____"<<endl;
    }
    else
    {
        cout<<"SolutionWorker::retrieveSolverData()->____\"RawSolverOutput.txt\" already opened by simulationManager____"<<endl;
    }
}
