#ifndef INTERPOLATORCONTROLLER_H
#define INTERPOLATORCONTROLLER_H

//! ---
//! Qt
//! ---
#include <QObject>
#include <QThread>

//! ----------------
//! custom includes
//! ----------------
#include "simulationnodeclass.h"
#include "mapper3dclass.h"

//! ----
//! OCC
//! ----
#include "occhandle.h"
#include <MeshVS_DataSource.hxx>

class QProgressIndicator;

class interpolatorController: public QObject
{
    Q_OBJECT

public:

    QThread workerThread;
    QProgressIndicator *myProgressIndicator;
    SimulationNodeClass *myNode;

    //! constructor
    interpolatorController(const std::vector<Mapper3DClass::nodeSource> &vecSourceNodes,
                           bool remap, int NbRemapSteps,
                           int Nx, int Ny, int Nz,
                           QProgressIndicator *aProgressIndicator=NULL,
                           QObject* parent=0);


    //! set progress indicator
    void setProgressIndicator(QProgressIndicator *aProgressIndicator) { myProgressIndicator = aProgressIndicator; }

private slots:

    //! for stopping the thread
    void stopThread();

    //! lock node
    void lockNode(SimulationNodeClass *node);

signals:

    //! for starting the thread
    void operate(SimulationNodeClass *node, int theAlgo, const occHandle(MeshVS_DataSource) &theTargetMesh);

    //! task finished or stopped - experimental
    void taskFinishedOrStopped();
};

#endif // INTERPOLATORCONTROLLER_H
