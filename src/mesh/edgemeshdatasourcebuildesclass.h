#ifndef EDGEMESHDATASOURCEBUILDESCLASS_H
#define EDGEMESHDATASOURCEBUILDESCLASS_H

//! redefinition of opencascade Handle() macro
#include "occhandle.h"

//! Qt
#include <QObject>

//! custom includes
#include <meshdatabase.h>
#include <mesh.h>
#include <ng_meshvs_datasource1d.h>

//! OCC
#include <NCollection_Array1.hxx>

class edgeMeshDataSourceBuildesClass: public QObject
{

    Q_OBJECT

public:

    //! constructor
    edgeMeshDataSourceBuildesClass(meshDataBase *theMeshDataBase, QObject *parent=0);

    //! constructor

    //! perform on a body
    bool performOne(const TopoDS_Shape &theShape, NCollection_Array1<opencascade::handle<Ng_MeshVS_DataSource1D>> &edgeMeshArray);

    //! perform on all the bodies
    void perform();

    //! experimental
    QMap<int, QList<mesh::meshPoint>> getMap() {return edgeNr_to_listOfPoints;}

private:

    //! the mesh database
    meshDataBase *myDB;

    //! experimental
    QMap<int, QList<mesh::meshPoint>> edgeNr_to_listOfPoints;
};

#endif // EDGEMESHDATASOURCEBUILDESCLASS_H
