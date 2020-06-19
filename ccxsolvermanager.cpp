//! ----------------
//! custom includes
//! ----------------
#include "ccxsolvermanager.h"
#include "solutionworker.h"
#include "qprogressevent.h"

//! ---
//! Qt
//! ---
#include <QVBoxLayout>
#include <QProgressBar>
#include <QPushButton>
#include <QMessageBox>

//! ------------------------
//! function: constructor I
//! details:
//! ------------------------
CCXSolverManager::CCXSolverManager(const QString &inputFileName, QWidget *parent):QDialog(parent),
    myInputFile(inputFileName)
{
    cout<<"CCXSolverManager::CCXSolverManager()->____constructor called____"<<endl;
    cout<<"CCXSolverManager::CCXSolverManager()->____launching input file: "<<inputFileName.toStdString()<<"____"<<endl;

    //! ---------------------------
    //! create the graphic content
    //! ---------------------------
    this->createContent();

    //! ------------------------------------------------------
    //! create the solution worker: do not specify the parent
    //! ------------------------------------------------------
    SolutionWorker *theSolutionWorker = new SolutionWorker(myInputFile);
    theSolutionWorker->moveToThread(&workerThread);

    //! --------------------------------------------------------
    //! delete the solution worker when the thread has finished
    //! not sure about it is needed or not ...
    //! --------------------------------------------------------
    connect(&workerThread, &QThread::finished, theSolutionWorker, &QObject::deleteLater);

    //! ------------------------------------------------------------
    //! start the worker when CCXSolverManager::operate() is called
    //! ------------------------------------------------------------
    connect(this, &CCXSolverManager::operate, theSolutionWorker, &SolutionWorker::perform);
    connect(theSolutionWorker, SIGNAL(resultReady(int)), this, SLOT(emitSolutionProcessFinished(int)));

    //! -----------------------------------------------------------------------
    //! close the progress bar when the solution workers emit resultReady(int)
    //! -----------------------------------------------------------------------
    connect(theSolutionWorker, SIGNAL(resultReady(int)), this, SLOT(close()));

    //! -----------------------------------
    //! start the event loop of the thread
    //! -----------------------------------
    workerThread.start();
}

//! ---------------------
//! function: destructor
//! details:
//! ---------------------
CCXSolverManager::~CCXSolverManager()
{
    workerThread.quit();
    workerThread.wait();
}

//! -----------------------
//! function: setInputFile
//! details:
//! -----------------------
void CCXSolverManager::setInputFile(const QString &inputFile)
{
    myInputFile = inputFile;
}

//! ----------------------
//! function: eventFilter
//! details:
//! ----------------------
bool CCXSolverManager::eventFilter(QObject *object, QEvent *event)
{
    if(event->type()==QProgressEvent::type())
    {
        QProgressEvent *e = static_cast<QProgressEvent*>(event);
        switch(e->action())
        {
        case QProgressEvent_Init:
        {
            int min = e->min();
            int max = e->max();
            prgbar->setRange(min,max);
            prgbar->setValue(min);
        }
            break;
        case QProgressEvent_Update:
        {
            prgbar->setValue(e->val());
        }
            break;
        case QProgressEvent_Reset:
        {
            prgbar->setValue(e->min());
        }
            break;
        }
        //! the event stops here
        return true;
    }
    return QObject::eventFilter(object, event);
}

//! ------------------------
//! function: createContent
//! details:
//! ------------------------
void CCXSolverManager::createContent()
{
    this->setWindowTitle("Solver manager");
    this->setObjectName("solverManager");
    this->installEventFilter(this);
    this->setModal(false);

    //! -------------
    //! progress bar
    //! -------------
    QVBoxLayout *vl = new QVBoxLayout(this);
    prgbar = new QProgressBar(this);
    vl->addWidget(prgbar);

    QString progressBarStyleSheet = "QProgressBar"
                                    "{"
                                    "border: 2px solid grey;"
                                    "border-radius: 5px;"
                                    "text-align: center;"
                                    "}";

    prgbar->setStyleSheet(progressBarStyleSheet);
    prgbar->setFixedWidth(300);

    //! ------------
    //! stop button
    //! ------------
    myStopButton = new QPushButton(this);
    myStopButton->setText("Stop");
    myStopButton->setFixedWidth(100);

    //! ------------------------------
    //! horizontal layout for buttons
    //! ------------------------------
    QHBoxLayout *hl = new QHBoxLayout();
    hl->addWidget(myStopButton);

    vl->addLayout(hl);
    this->setLayout(vl);

    connect(myStopButton,SIGNAL(pressed()),this,SLOT(terminateRun()));

    //! -------------------------------
    //! disable the "X" (close) button
    //! -------------------------------
    setWindowFlags(((this->windowFlags() | Qt::CustomizeWindowHint) & ~Qt::WindowCloseButtonHint));

    //! ----------------
    //! show the dialog
    //! ----------------
    this->show();
}

//! -----------------------
//! function: terminateRun
//! details:
//! -----------------------
void CCXSolverManager::terminateRun()
{
    cout<<"CCXSolverManager::terminateRun()->____function called____"<<endl;

    //! ------------------------------------------
    //! a brute force terminate: kill the process
    //! to be modified ... to do ...
    //! ------------------------------------------
    workerThread.terminate();

    QString solutionFilePath = myInputFile;
    QString relativeSolutionFileName = solutionFilePath.split("/").last();
    solutionFilePath.chop(relativeSolutionFileName.length()+1);

    //! --------------------------------------------------------
    //! probably the solutionReady signal should not be emitted
    //! when the user has stopped the solver ...
    //! --------------------------------------------------------
    emit solutionReady(solutionFilePath);

    emit stopTimer();
    this->close();
}

//! --------------------------------------
//! function: emitSolutionProcessFinished
//! details:
//! --------------------------------------
void CCXSolverManager::emitSolutionProcessFinished(int code)
{
    cout<<"CCXSolverManager::emitSolutionProcessFinished()->____emitting signal solutionReady()____"<<endl;
    if(code!=0)
    {
        if(code==9999)
        {
            QMessageBox::critical(this, QString("Solver manager"),QString("The CCX solver is not installed on this computer"),QMessageBox::Ok);
            emit stopTimer();
            //this->close();
        }
        else
        {
            QMessageBox::critical(this, QString("Solver manager"),QString("Error: CCX cannot start (code: %1) Check your model").arg(code),QMessageBox::Ok);
            emit stopTimer();
            //this->close();
        }
    }
    else
    {
        //! ----------------------------------------------------
        //! emit the signal solution ready only if the solution
        //! (also a partial solution) has been obtained
        //! ----------------------------------------------------
        QString solutionFilePath = myInputFile;
        QString relativeSolutionFileName = solutionFilePath.split("/").last();
        solutionFilePath.chop(relativeSolutionFileName.length()+1);
        //cout<<"CCXSolverManager::emitSolutionProcessFinished()->____where is solver directory: "<<solutionFilePath.toStdString()<<"____"<<endl;
        emit solutionReady(solutionFilePath);
        emit stopTimer();
        //this->close();
    }
}

//! ------------------------
//! function: closeManager
//! details:  to be deleted
//! ------------------------
void CCXSolverManager::closeManager(int code)
{
    cout<<"@------------------------------------"<<endl;
    cout<<"@- closing manager with code: "<<code<<endl;
    cout<<"@------------------------------------"<<endl;
    this->close();
}
