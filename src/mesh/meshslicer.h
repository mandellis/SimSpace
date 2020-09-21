#ifndef MESHSLICER_H
#define MESHSLICER_H

//! ----------------
//! custom includes
//! ----------------
#include "occhandle.h"
#include <ng_meshvs_datasource3d.h>


class meshSlicer
{
public:

    //! constructor
    meshSlicer(const occHandle(Ng_MeshVS_DataSource3D) &aVolumeMeshDS = occHandle(Ng_MeshVS_DataSource3D)()):myVolumeMeshDS(aVolumeMeshDS)
    {;}

    //! set mesh data source
    void setMeshDataSource(const occHandle(Ng_MeshVS_DataSource3D) &aVolumeMeshDS) { myVolumeMeshDS = aVolumeMeshDS; }

    //! perform
    bool perform(double a, double b, double c, double d, occHandle(Ng_MeshVS_DataSource3D) &slicedVolumeMeshDS);

private:

    occHandle(Ng_MeshVS_DataSource3D) myVolumeMeshDS;

};

#endif // MESHSLICER_H
