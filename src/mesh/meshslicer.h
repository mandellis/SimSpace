#ifndef MESHSLICER_H
#define MESHSLICER_H

//! ----------------
//! custom includes
//! ----------------
#include "occhandle.h"
#include <MeshVS_DataSource.hxx>
#include <ng_meshvs_datasource3d.h>

class meshSlicer
{
public:

    //! constructor
    meshSlicer(const occHandle(MeshVS_DataSource) &aMeshDS = occHandle(MeshVS_DataSource)());

    //! set mesh data source
    void setMeshDataSource(const occHandle(MeshVS_DataSource) &aMeshDS) { myMeshDS = aMeshDS; }

    //! perform
    bool perform(double a, double b, double c, double d, occHandle(TColStd_HPackedMapOfInteger) &hiddenElementIDs);

private:

    occHandle(MeshVS_DataSource) myMeshDS;
};

#endif // MESHSLICER_H
