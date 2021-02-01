#ifndef CCXSOLVERMANAGER_H
#define CCXSOLVERMANAGER_H

#include <QWidget>
#include <QThread>
#include <QString>

#include "ccxsolvermessage.h"
#include "qccxsolvermessageevent.h"
#include "qconsoleevent.h"
#include <QApplication>
#include <QProcess>

class QProgressIndicator;

class CCXSolverManager: public QThread
{
    Q_OBJECT

private:

    void run() Q_DECL_OVERRIDE;

public:

    CCXSolverManager(QObject *parent = 0);
    void setNbProcessors(int NbProcessors=2);
    void setInputFile(const QString &inputFile);
    void setProgressIndicator(QProgressIndicator *aProgressIndicator);

private:

    QWidget *mySimulationManager;
    QWidget *mySolutionMonitorPanel;
    QString myInputFileName;
    QProgressIndicator *myProgressIndicator;
    QString myCCXLog;
    int myNbProcessors;

private slots:

public slots:

signals:

    void CCXRunFinished();
};


//! ---------------------
//! retrieve solver info
//! ---------------------
class Worker: public QObject
{
    Q_OBJECT

public:

    //! ----------------------
    //! function: constructor
    //! ----------------------
    Worker(QProcess *aProcess, QString inputFile, QString &alog = QString(""), QObject *parent = 0):
        process(aProcess),
        myInputFileName(inputFile),
        myLog(alog),
        QObject(parent)
    {
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

private:

    QWidget *mySimulationManager;
    QWidget *mySolutionMonitorPanel;
    QString myInputFileName;
    QString myLog;
    QProcess* process;

private:

    //! --------------------
    //! function: copy_file
    //! --------------------
    void copy_file(const char* srce_file, const char* dest_file)
    {
        std::ifstream srce(srce_file, std::ios::binary);
        std::ofstream dest(dest_file, std::ios::binary);
        dest<<srce.rdbuf();
    }

private slots:

    //! -----------------------------------
    //! function: removeRawSolverData_copy
    //! details:
    //! -----------------------------------
    void removeRawSolverData_copy()
    {
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

        QFile file(rawSolverDataFileName_copy);
        if(file.exists())
        {
            bool removed = file.remove();
            if(removed==true) cout<<"Worker::removeRawSolverData_copy()->____file removed____"<<endl;
            cout<<"Worker::removeRawSolverData_copy()->____file NOT removed____"<<endl;
        }
    }

    //! -----------------------------
    //! function: retrieveSolverData
    //! details:
    //! -----------------------------
    void retrieveSolverData()
    {
        cout<<"Worker::retrieveSolverData()->____get called from: "<<QThread::currentThreadId()<<"____"<<endl;
        QString textBlock = process->readAllStandardOutput();
        myLog.append(textBlock);

        //! -------------------------------------------------
        //! send the solver output to the simulation manager
        //! -------------------------------------------------
        if(mySimulationManager!=Q_NULLPTR)
        {
            CCXSolverMessage msg(myLog);
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
            this->copy_file(rawSolverDataFileName.toStdString().c_str(), rawSolverDataFileName_copy.toStdString().c_str());
            fclose(lockFile);
            bool lockRemoved = QFile::remove(lockFileName);
            if(lockRemoved) cout<<"SolutionWorker::retrieveSolverData()->____lock file removed____"<<endl;
        }
        else
        {
            cout<<"SolutionWorker::retrieveSolverData()->____\"RawSolverOutput.txt\" already opened by simulationManager____"<<endl;
        }
    }
};

#endif // CCXSOLVERMANAGER_H
