#ifndef MESHSLICER_H
#define MESHSLICER_H

//! ----------------
//! custom includes
//! ----------------
#include "occhandle.h"
#include <MeshVS_DataSource.hxx>
#include <ng_meshvs_datasource3d.h>

class Ng_MeshVS_DataSourceFace;

class meshSlicer
{
public:

    //! constructor
    meshSlicer(const occHandle(MeshVS_DataSource) &aMeshDS = occHandle(MeshVS_DataSource)());

    //! set mesh data source
    void setMeshDataSource(const occHandle(MeshVS_DataSource) &aMeshDS) { myMeshDS = aMeshDS; }

    //! perform
    bool perform(double a, double b, double c, double d, occHandle(TColStd_HPackedMapOfInteger) &hiddenElementIDs);

    //! plane mesh intersection
    bool planeMeshIntersection(occHandle(Ng_MeshVS_DataSourceFace) &aSliceMesh, double a, double b,double c, double d);

private:

    occHandle(MeshVS_DataSource) myMeshDS;
};

#endif // MESHSLICER_H
