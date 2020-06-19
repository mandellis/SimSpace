//! ----------------
//! custom includes
//! ----------------
#include "meshworker.h"
#include <qprogressevent.h>
#include <qprogressindicator.h>

//! ---
//! Qt
//! ---
#include <QTimer>

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
meshWorker::meshWorker(meshDataBase *mDB, QList<int> bodies, bool isVolume, QProgressIndicator *aProgressIndicator, QObject *parent):
    QObject(parent),
    myListOfBodies(bodies),
    myIsVolume(isVolume)
{
    //! -----------------------
    //! the progress indicator
    //! -----------------------
    myProgressIndicator = aProgressIndicator;
    cout<<"meshWorker::meshWorker()->____"<<myProgressIndicator->getText().toStdString()<<"____"<<endl;

    //! -----------------
    //! the mesher class
    //! -----------------
    theMesher = new MesherClass(mDB,bodies,isVolume,this);
    theMesher->setProgressIndicator(myProgressIndicator);
    connect(theMesher,SIGNAL(meshingFinished()),this,SLOT(emitResultReady()));

    //! ----------------
    //! netgen specific
    //! ----------------
    connect(theMesher,SIGNAL(requestStartingNetgenEnquireTimer()),this,SLOT(emitRequestStartingNetgenEnquireTimer()));
    connect(theMesher,SIGNAL(requestStoppingNetgenEnquireTimer()),this,SLOT(emitRequestStoppingNetgenEnquireTimer()));

    //! ----------------------------
    //! timer - diagnostic purposes
    //! ----------------------------
    QTimer *timer = new QTimer(this);
    timer->setInterval(250);
    connect(timer,SIGNAL(timeout()),this,SLOT(update()));
    timer->start();
}

//! -----------------------------------------------------------------
//! function: update
//! details:  a message for checking if the thread is running or not
//! -----------------------------------------------------------------
void meshWorker::update()
{
    static int ct = 0;
    ct +=250;
    cout<<"meshWorker::update()->____Meshing thread is running. Elapsed time: "<<ct<<"____"<<endl;
}

//! -------------------------------
//! function: setProgressIndicator
//! details:
//! -------------------------------
void meshWorker::setProgressIndicator(QProgressIndicator *aProgressIndicator)
{
    myProgressIndicator = aProgressIndicator;
}

//! ------------------
//! function: perform
//! details:
//! ------------------
void meshWorker::perform()
{
    cout<<"meshWorker::perform()->____perform called____"<<endl;
    theMesher->generateMesh();
}

//! --------------------------
//! function: emitResultReady
//! details:
//! --------------------------
void meshWorker::emitResultReady()
{
    emit resultReady();
}

//! ------------------------------------------------
//! function: emitRequestStartingNetgenEnquireTimer
//! details:
//! ------------------------------------------------
void meshWorker::emitRequestStartingNetgenEnquireTimer()
{
    emit requestStartingNetgenEnquireTimer();
}

//! ------------------------------------------------
//! function: emitRequestStoppingNetgenEnquireTimer
//! details:
//! ------------------------------------------------
void meshWorker::emitRequestStoppingNetgenEnquireTimer()
{
    emit requestStoppingNetgenEnquireTimer();
}
