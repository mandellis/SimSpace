#ifndef RAYINTERSECTMESH_H
#define RAYINTERSECTMESH_H

//! ---
//! Qt
//! ---
#include <QObject>

//! ----
//! OCC
//! ----
#include <ng_meshvs_datasourceface.h>

//! ----------------
//! custom includes
//! ----------------
#include "occhandle.h"

//! ----
//! C++
//! ----
#include <memory>

class rayIntersectMesh: public QObject
{
    Q_OBJECT

public:

    rayIntersectMesh(QObject *parent=0);
    rayIntersectMesh(const occHandle(MeshVS_DataSource) &aMeshDS, QObject *parent =0);
    rayIntersectMesh(const occHandle(MeshVS_DataSource) &aMeshDS, double *origin, double *direction, QObject *parent =0);
    void setRay(double *origin, double *direction);
    void setMesh(const occHandle(MeshVS_DataSource) &aMeshDS);

private:

    double *myDirection;
    double *myOrigin;
    occHandle(MeshVS_DataSource) myMeshDS;

    std::map<meshElement2D,int> mapOfTriangles;

private:

    bool planeMeshIntersection(occHandle(Ng_MeshVS_DataSourceFace) &aSliceMesh);
};

#endif // RAYINTERSECTMESH_H
