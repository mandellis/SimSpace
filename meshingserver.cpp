//! ----------------
//! custom includes
//! ----------------
#include "meshingserver.h"
#include "meshworker.h"
#include <qprogressindicator.h>
#include <qprogressevent.h>

//! ---
//! Qt
//! ---
#include <QTimer>
#include <QApplication>

//! ------
//! nglib
//! ------
namespace nglib
{
#include <nglib.h>
}
using namespace nglib;

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
MeshingServer::MeshingServer(QProgressIndicator *aProgressIndicator, QObject *parent): QObject(parent),myIsMeshingRunning(false)
{
    //! -----------------------
    //! the progress indicator
    //! -----------------------
    myProgressIndicator = aProgressIndicator;
    myProgressIndicator->setSecondaryBarVisible(true);
}

//! -----------------------------------
//! function: constructor
//! details:  not called at the moment
//! -----------------------------------
MeshingServer::MeshingServer(meshDataBase *mDB,
                             const QList<int> &listOfBodies,
                             bool isVolume,
                             QProgressIndicator *aProgressIndicator,
                             QObject *parent): QObject(parent),myIsMeshingRunning(false)
{
    //! -----------------------
    //! the progress indicator
    //! -----------------------
    myProgressIndicator = aProgressIndicator;
    myProgressIndicator->setSecondaryBarVisible(true);
    //if(listOfBodies.length()==1) myProgressIndicator->setPrimaryBarVisible(false);
    //        else myProgressIndicator->setPrimaryBarVisible(true);

    //! -----------------------
    //! create the mesh worker
    //! -----------------------
    meshWorker *theMeshWorker = new meshWorker(mDB,listOfBodies,isVolume,myProgressIndicator);
    theMeshWorker->setProgressIndicator(myProgressIndicator);
    theMeshWorker->moveToThread(&workerThread);

    //! ------------------------------------------------------------------------------------
    //! [*] from Qt documentation:
    //! "From Qt 4.8 onwards, it is possible to deallocate objects that live in a thread
    //! that has just ended, by connecting the finished() signal to QObject::deleteLater()"
    //! ------------------------------------------------------------------------------------
    disconnect(&workerThread, &QThread::finished, theMeshWorker, &QObject::deleteLater);
    disconnect(this, &MeshingServer::operate, theMeshWorker, &meshWorker::perform);
    disconnect(theMeshWorker, &meshWorker::resultReady, this, &MeshingServer::handleResult);

    connect(this, &MeshingServer::operate, theMeshWorker, &meshWorker::perform);

    //! ---------------------------------------------
    //! [*] but: this is not able to stop the thread
    //! ---------------------------------------------
    connect(&workerThread, &QThread::finished, theMeshWorker, &QObject::deleteLater);

    //! ----------------------------------------
    //! instead this is able to stop the thread
    //! ----------------------------------------
    connect(theMeshWorker, &meshWorker::resultReady, this, &MeshingServer::handleResult);

    //! ---------------------
    //! start the event loop
    //! ---------------------
    workerThread.start();

    //! -------------
    //! experimental
    //! -------------
    disconnect(myProgressIndicator,SIGNAL(abortNetgenSTLPressed()),this,SLOT(abortNetgenSTL()));
    connect(myProgressIndicator,SIGNAL(abortNetgenSTLPressed()),this,SLOT(abortNetgenSTL()));
}

//! ---------------------
//! function: isRunning
//! details:  get status
//! ---------------------
bool MeshingServer::isMeshingRunning()
{
    return myIsMeshingRunning;
}

//! ---------------
//! function: init
//! details:
//! ---------------
void MeshingServer::init(meshDataBase *mDB,
                         const QList<int> &listOfBodies,
                         bool isVolume,
                         QProgressIndicator *aProgressIndicator)
{
    //! ----------------
    //! status variable
    //! ----------------
    myIsMeshingRunning = false;

    //! -----------------------
    //! the progress indicator
    //! -----------------------
    myProgressIndicator = aProgressIndicator;
    if(myProgressIndicator!=Q_NULLPTR)
    {
        myProgressIndicator->setSecondaryBarVisible(true);
        cout<<"MeshingServer::init()->____"<<myProgressIndicator->getText().toStdString()<<"____"<<endl;
    }

    //! -----------------------------------
    //! create the mesh worker - no parent
    //! -----------------------------------
    meshWorker *theMeshWorker = new meshWorker(mDB,listOfBodies,isVolume,myProgressIndicator);
    theMeshWorker->moveToThread(&workerThread);

    //! ----------------------------------------------------------
    //! the "abortPressed()" signal is emitted by the main thread
    //! ----------------------------------------------------------
    //disconnect(myProgressIndicator,SIGNAL(abortPressed()),this,SLOT(abortMeshing()));
    //connect(myProgressIndicator,SIGNAL(abortPressed()),this,SLOT(abortMeshing()));

    disconnect(myProgressIndicator,SIGNAL(abortNetgenSTLPressed()),this,SLOT(abortSTLNetgen()));
    connect(myProgressIndicator,SIGNAL(abortNetgenSTLPressed()),this,SLOT(abortSTLNetgen()));

    //! ------------------------------------------------------------------------------------
    //! [*] from Qt documentation:
    //! "From Qt 4.8 onwards, it is possible to deallocate objects that live in a thread
    //! that has just ended, by connecting the finished() signal to QObject::deleteLater()"
    //! ------------------------------------------------------------------------------------

    //! ---------------------------------------------
    //! [*] but: this is not able to stop the thread
    //! ---------------------------------------------
    disconnect(&workerThread, &QThread::finished, theMeshWorker, &QObject::deleteLater);
    connect(&workerThread, &QThread::finished, theMeshWorker, &QObject::deleteLater);

    //! --------------------------------------------
    //! ... instead this is able to stop the thread
    //! --------------------------------------------
    disconnect(theMeshWorker, &meshWorker::resultReady, this, &MeshingServer::handleResult);
    connect(theMeshWorker, &meshWorker::resultReady, this, &MeshingServer::handleResult);

    //! --------------------------------------
    //! this for starting the meshing process
    //! --------------------------------------
    disconnect(this, &MeshingServer::operate, theMeshWorker, &meshWorker::perform);
    connect(this, &MeshingServer::operate, theMeshWorker, &meshWorker::perform);

    //! ---------------------
    //! start the event loop
    //! ---------------------
    if(!workerThread.isRunning()) workerThread.start();

    //! ---------------------------
    //! update the status variable
    //! ---------------------------
    myIsMeshingRunning = true;

    //! -----------------------------------------------------------------
    //! Retrieve the netgen progress information
    //! This timer is used for retrieving at a constant rate the Netgen
    //! status, by calling the nglib function "Ng_getMeshingStatus"
    //! -----------------------------------------------------------------
    myNetgenInquireTimer = new QTimer(this);
    myNetgenInquireTimer->setInterval(125);
    connect(myNetgenInquireTimer,SIGNAL(timeout()),this,SLOT(showNetgenProgress()));
    connect(theMeshWorker,SIGNAL(requestStoppingNetgenEnquireTimer()),this,SLOT(stopNetgenEnquireTime()));
    connect(theMeshWorker,SIGNAL(requestStartingNetgenEnquireTimer()),this,SLOT(startNetgenEnquireTime()));
}

//! -----------------------
//! function: handleResult
//! details:
//! -----------------------
void MeshingServer::handleResult()
{
    cout<<"MeshingServer::handleResult()->____stopping thread____"<<endl;
    this->stopThread();

    myIsMeshingRunning = false;
    bool isSuccessfull = true;
    emit meshingFinished(isSuccessfull);
}

//! -------------------------
//! function: abortSTLNetgen
//! details:
//! -------------------------
void MeshingServer::abortSTLNetgen()
{
    cout<<"MeshingServer::abortSTLNetgen()->____function called____"<<endl;

    myNetgenInquireTimer->stop();

    //! ------------
    //! stop Netgen
    //! ------------
    Ng_Terminate();

    this->stopThread();
    myIsMeshingRunning = false;
    bool isSuccessfull = false;
    emit meshingFinished(isSuccessfull);
}

//! --------------------------------------------------------------------
//! function: stopThread
//! details:  see note [*] => actually when connecting to "deleteLater"
//!           the thread do not terminate (the timer continues to run)
//! --------------------------------------------------------------------
void MeshingServer::stopThread()
{
    cout<<"MeshingServer::stopThread()->____function called____"<<endl;

    //! -------------------------------------------------------------------------
    //! use this interface test different ways for stopping the thread correctly
    //! -------------------------------------------------------------------------
    if(workerThread.isRunning()) workerThread.exit();

    //! -----------
    //! diagnostic
    //! -----------
    myNetgenInquireTimer->stop();
}

//! -----------------------
//! function: showNetgenProgress
//! details:
//! -----------------------
void MeshingServer::showNetgenProgress()
{
    static std::string oldmsg;

    if(myProgressIndicator!=Q_NULLPTR)
    {
        if(myProgressIndicator->isVisible())
        {
            double progress = Ng_getMeshingStatus().percent;
            std::string msg = Ng_getMeshingStatus().task;
            if(msg!="") oldmsg = msg;

            /*
            if(msg =="Surface meshing" || msg == "Volume meshing" || msg =="Volume optimization" ||
                    msg == "Find points" || msg == "Combine Improve" || msg == "Split Improve" ||
                    msg == "Swap Improve" || msg =="Swap Improve Surface" || msg =="Edge meshing" ||
                    msg == "Optimizing surface" || msg =="Setting local mesh size (elements per edge)" ||
                    msg == "Setting local mesh size (edge curvature)" || msg =="Setting local mesh size (face curvature)" ||
                    msg == "Setting local mesh size (close edges)" || msg == "Smooth Mesh" ||
                    msg == "Smooth Mesh Jacobian")
            */

            //! ----------------------
            //! post a progress event
            //! ----------------------
            QProgressEvent *progressEvent = new QProgressEvent(QProgressEvent_None,0,9999,0,
                                                               QString::fromStdString(oldmsg),
                                                               QProgressEvent_Update,0,100,int(progress),"Netgen meshing");
            QApplication::postEvent(myProgressIndicator,progressEvent);
            QApplication::processEvents();
        }
    }
}

//! --------------------------------
//! function: stopNetgenEquireTimer
//! details:
//! --------------------------------
void MeshingServer::stopNetgenEnquireTime()
{
    myNetgenInquireTimer->stop();
}

//! ---------------------------------
//! function: startNetgenEquireTimer
//! details:
//! ---------------------------------
void MeshingServer::startNetgenEnquireTime()
{
    myNetgenInquireTimer->start();
}
