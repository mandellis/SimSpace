//! ----------------
//! custom includes
//! ----------------
#include "ccxsolvermanager.h"
#include "qprogressindicator.h"
#include <qprogressevent.h>

//! ---
//! Qt
//! ---
#include <QProcess>
#include <QApplication>

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
CCXSolverManager::CCXSolverManager(QObject *parent):QThread(parent)
{
    myNbProcessors = 2;
}

//! -----------------------
//! function: setInputFile
//! details:
//! -----------------------
void CCXSolverManager::setInputFile(const QString &inputFile)
{
    myInputFileName = inputFile;
}

//! --------------------------
//! function: setNbProcessors
//! details:
//! --------------------------
void CCXSolverManager::setNbProcessors(int NbProcessors)
{
    if(NbProcessors<1) NbProcessors = 2;
    myNbProcessors = NbProcessors;
}

//! -------------------------------
//! function: setProgressIndicator
//! details:
//! -------------------------------
void CCXSolverManager::setProgressIndicator(QProgressIndicator *aProgressIndicator)
{
    myProgressIndicator = aProgressIndicator;
}

//! --------------
//! function: run
//! details:
//! --------------
void CCXSolverManager::run()
{
    //! -----------------
    //! create a process
    //! -----------------
    QProcess process;

    //! --------------------------
    //! this send messages around
    //! --------------------------
    QString file = myInputFileName;
    Worker aWorker(&process,file);
    connect(&process,SIGNAL(readyReadStandardOutput()),&aWorker,SLOT(retrieveSolverData()));
    connect(&process,SIGNAL(finished(int)),&aWorker,SLOT(removeRawSolverData_copy()));
    connect(myProgressIndicator,SIGNAL(requestStopCCX()),&process,SLOT(kill()));

    //! ---------------------------------------
    //! calculix solver: program name and args
    //! ---------------------------------------
    QString program;
    bool CCX_uses_SPOOLES = false;
    if(CCX_uses_SPOOLES)
    {
        program = QString::fromLatin1(std::getenv("CCX_MT_SPOOLES"));
        if(program.isEmpty())
        {
            emit CCXRunFinished();
            return;
        }
        program.append("/ccx_2.16_spooles.exe");
    }
    else
    {
        program = QString::fromLatin1(std::getenv("CCX_MT_PARDISO"));
        if(program.isEmpty())
        {
            emit CCXRunFinished();
            return;
        }
        //program.append("/ccx_215Pardiso.exe");
        program.append("/ccx.exe");
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

    process.setWorkingDirectory(solutionDataDir);
    process.setProcessEnvironment(processEnv);
    process.start(program,CCXarguments);
    process.waitForFinished(-1);

    //! -------------
    //! run finished
    //! -------------
    emit CCXRunFinished();
}
