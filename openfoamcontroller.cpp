//! ----------------
//! custom includes
//! ----------------
#include "openfoamcontroller.h"
#include "openfoamreader.h"
#include "qprogressevent.h"
#include "qprogressindicator.h"
#include "qextendedstandarditem.h"

//! -------
//! system
//! -------
#include <Windows.h>

//! ---
//! Qt
//! ---
#include <QTimer>

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
openFoamController::openFoamController(/*const QString &sourceDirPath,
                                       const QString &targetDirPath,*/
                                       int fileMode,
                                       QProgressIndicator *aProgressIndicator,
                                       QObject *parent): QObject(parent)
{
    myProgressIndicator = aProgressIndicator;

    //! ---------------------
    //! the target directory
    //! ---------------------
    //myTargetDirectory = targetDirPath;

    //! --------------
    //! Worker thread
    //! --------------
    OpenFoamReader *OFReader = new OpenFoamReader(/*sourceDirPath,targetDirPath,*/fileMode);
    OFReader->setProgressIndicator(myProgressIndicator);
    OFReader->moveToThread(&workerThread);
/*
#ifdef COSTAMP_VERSION
    OFReader->setTimeFolders(myTimeFolders);
#endif*/
    //! -----------------------------------
    //! connection for starting the thread
    //! -----------------------------------
    disconnect(this,SIGNAL(operate(SimulationNodeClass*)),OFReader,SLOT(perform(SimulationNodeClass*)));
    connect(this,SIGNAL(operate(SimulationNodeClass*)),OFReader,SLOT(perform(SimulationNodeClass*)));
    connect(this,SIGNAL(operate(SimulationNodeClass*)),this,SLOT(lockNode(SimulationNodeClass*)));

    //! ------------------------------------------------------------------------------------
    //! [*] from Qt documentation:
    //! "From Qt 4.8 onwards, it is possible to deallocate objects that live in a thread
    //! that has just ended, by connecting the finished() signal to QObject::deleteLater()"
    //! But: this is not able to stop the thread
    //! ------------------------------------------------------------------------------------
    disconnect(&workerThread, &QThread::finished, OFReader, &QObject::deleteLater);
    connect(&workerThread, &QThread::finished, OFReader, &QObject::deleteLater);
    //! -----------------------------------
    //! connection for stopping the thread
    //! this is able to stop the thread
    //! -----------------------------------
    disconnect(OFReader,SIGNAL(conversionFinished()),this,SLOT(stopThread()));
    connect(OFReader,SIGNAL(conversionFinished()),this,SLOT(stopThread()));
    //! ------------------------------------------------
    //! programmatically stop the thread - experimental
    //! ------------------------------------------------
    if(myProgressIndicator!=NULL) connect(myProgressIndicator,SIGNAL(stopPressed()),this,SLOT(stopThread()));
    //! ---------------------
    //! start the event loop
    //! ---------------------
    workerThread.start();
}

//! ---------------------
//! function: stopThread
//! details:
//! ---------------------
void openFoamController::stopThread()
{
    if(workerThread.isRunning())
    {
        cout<<"openFoamController::stopThread()->____stopping the thread____"<<endl;
        workerThread.quit();
        workerThread.terminate();
        Sleep(1000);
        myProgressIndicator->hide();

        //! ----------------
        //! unlock the node
        //! ----------------
        if(myNode!=NULL)
        {
            for(int i=0; i<myNode->getPropertyItems().size(); i++) myNode->getPropertyItems().at(i)->setEditable(true);
        }
    }

    //! ------------------------------------------------------------
    //! check uncomplete: remove the uncomplete translations
    //! search for a .lck file, remove it and the corresponding file
    //! -------------------------------------------------------------

    QDir resultsDataDir(myTargetDirectory);
    resultsDataDir.setFilter(QDir::NoDotAndDotDot | QDir::Files);
    QList<QString> fileList = resultsDataDir.entryList();
    QList<QString> filesToRemove;
    for(int i=0; i<fileList.length(); i++)
    {
        QString fileName = fileList.at(i);
        QString fileName1 = fileName.split(".").last();
        if(fileName1.split(".").last()=="lck")
        {
            cout<<"____"<<fileName.toStdString()<<"____"<<endl;
            cout<<"____"<<fileName1.toStdString()<<"____"<<endl;
            filesToRemove<<fileName;
            fileName.chop(4);
            filesToRemove<<fileName;
            exit(1);
        }
    }
    int found = false;
    for(int i=0; i<filesToRemove.length(); i++)
    {
        found = true;
        cout<<"____removing files: "<<filesToRemove.at(i).toStdString()<<"____"<<endl;
    }
    if(found) exit(1);

    //! --------------------------------------
    //! notify conversion finished or stopped
    //! --------------------------------------
    emit taskFinishedOrStopped();
}

//! ----------------------------
//! function: lockNode
//! details:  lock the controls
//! ----------------------------
void openFoamController::lockNode(SimulationNodeClass *node)
{
    myNode = node;
    for(int i=0; i<node->getPropertyItems().size(); i++) node->getPropertyItems().at(i)->setEditable(false);
}
