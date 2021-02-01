//! ----------------
//! custom includes
//! ----------------
#include "interpolatorcontroller.h"
#include "qprogressevent.h"
#include "qprogressindicator.h"
#include "qextendedstandarditem.h"

//! -------
//! system
//! -------
#include <Windows.h>

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
interpolatorController::interpolatorController(const std::vector<Mapper3DClass::nodeSource> &vecSourceNodes,
                                               bool remap, int NbRemapSteps,
                                               int Nx, int Ny, int Nz,
                                               QProgressIndicator *aProgressIndicator,
                                               QObject* parent):
    QObject(parent)
{
    //! --------------
    //! Worker thread
    //! --------------
    Mapper3DClass *mapper = new Mapper3DClass();
    mapper->setProgressIndicator(aProgressIndicator);
    mapper->setSource(vecSourceNodes);
    mapper->setRemap(remap);
    mapper->setRemappingSteps(NbRemapSteps);
    mapper->setNbBuckets(Nx, Ny, Nz);
    mapper->moveToThread(&workerThread);

    //! -----------------------------------
    //! connection for starting the thread
    //! -----------------------------------
    disconnect(this,SIGNAL(operate(SimulationNodeClass*,int,occHandle(MeshVS_DataSource))),
               mapper,SLOT(interpolate(int,occHandle(MeshVS_DataSource))));
    connect(this,SIGNAL(operate(SimulationNodeClass*,int, occHandle(MeshVS_DataSource))),
            mapper,SLOT(interpolate(int,occHandle(MeshVS_DataSource))));

    //! ------------------------------
    //! lock the controls of the item
    //! ------------------------------
    connect(this,SIGNAL(operate(SimulationNodeClass*,int,occHandle(MeshVS_DataSource))),
            this,SLOT(lockNode(SimulationNodeClass*)));

    //! ------------------------------------------------------------------------------------
    //! [*] from Qt documentation:
    //! "From Qt 4.8 onwards, it is possible to deallocate objects that live in a thread
    //! that has just ended, by connecting the finished() signal to QObject::deleteLater()"
    //! But: this is not able to stop the thread
    //! ------------------------------------------------------------------------------------
    disconnect(&workerThread,&QThread::finished, mapper, &QObject::deleteLater);
    connect(&workerThread,&QThread::finished, mapper, &QObject::deleteLater);

    //! -----------------------------------
    //! connection for stopping the thread
    //! this is able to stop the thread
    //! -----------------------------------
    disconnect(mapper,SIGNAL(interpolationFinished()),this,SLOT(stopThread()));
    connect(mapper,SIGNAL(interpolationFinished()),this,SLOT(stopThread()));

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
void interpolatorController::stopThread()
{
    if(workerThread.isRunning())
    {
        cout<<"openFoamController::stopThread()->____trying to stop the thread____"<<endl;
        workerThread.terminate();
        Sleep(500);
        myProgressIndicator->hide();

        //! ----------------
        //! unlock the node
        //! ----------------
        if(myNode!=NULL)
        {
            for(int i=0; i<myNode->getPropertyItems().size(); i++) myNode->getPropertyItems().at(i)->setEditable(true);
        }
    }
}

//! ----------------------------
//! function: lockNode
//! details:  lock the controls
//! ----------------------------
void interpolatorController::lockNode(SimulationNodeClass *node)
{
    myNode = node;
    for(int i=0; i<node->getPropertyItems().size(); i++) node->getPropertyItems().at(i)->setEditable(false);
}
