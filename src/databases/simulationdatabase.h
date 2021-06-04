#ifndef SIMULATIONDATABASE_H
#define SIMULATIONDATABASE_H

//! ----------------
//! custom includes
//! ----------------
#include <meshdatabase.h>
#include "simulationnodeclass.h"
#include "qextendedstandarditem.h"
#include "src/main/mydefines.h"
#include <indexedmapofmeshdatasources.h>

//! ----
//! OCC
//! ----
#include <TopoDS_Shape.hxx>
#include <TopTools_ListOfShape.hxx>
#include <TCollection_AsciiString.hxx>

//! ---
//! Qt
//! ---
#include <QMap>
#include <QObject>

//! ----
//! C++
//! ----
#include <vector>
#include <iostream>

class simulationDataBase: public meshDataBase
{
    Q_OBJECT

protected:

    //! the remote point root
    QExtendedStandardItem *RemotePointsRootItem;

    //! the connections root item
    QExtendedStandardItem *ConnectionsRootItem;

    //! the named selection root item
    QExtendedStandardItem *NamedSelectionRootItem;

    //! create the connections root node
    void createConnectionsRootNode();

    //! create the named selection root node
    void createNamedSelectionRootNode();

    //! create the standard item model
    virtual void createStandardModel();


public:

    //! create the thermal analysis root node
    void createThermalAnalysisRootNode();

    //! create the structural analysis root node
    void createStructuralAnalysisRootNode();

    //! create a combined analysis root
    void createCombinedAnalysisRootNode();

    //! create CFD analysis root
    void createCFDAnalysisRootNode();

    //! create particles in fields root
    void createParticlesInFieldsRootNode();

public slots:

    //! create the remote points root node
    void createRemotePointRoot();

public:

    //! get the model
    inline virtual QStandardItemModel* getModel() const override { return myModel; }

public:

    //! constructor - default
    //simulationDataBase(QObject *parent=0);

    //! constructor
    simulationDataBase(const TopoDS_Shape &shape = TopoDS_Shape(), const QString& theFilePath="", QObject *parent=0);

    //! constructor - from file
    simulationDataBase(const QList<SimulationNodeClass*> listOfNodes, const QString &archiveFileName, QObject *parent=0);

    //! destructor
    ~simulationDataBase();

public:

    QMap<int,IndexedMapOfMeshDataSources> MapItemBoudaryConditionDS;
};

#endif // SIMULATIONDATABASE_H
