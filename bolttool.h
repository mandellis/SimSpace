#ifndef BOLTTOOL_H
#define BOLTTOOL_H

//! ----------------
//! custom includes
//! ----------------
#include "occhandle.h"
#include <ng_meshvs_datasource3d.h>
#include <ng_meshvs_datasourceface.h>

//! ----
//! C++
//! ----
#include <vector>

class boltTool
{
public:

    //! ------------
    //! constructor
    //! ------------
    boltTool(const occHandle(Ng_MeshVS_DataSource3D) &aVolumeMeshDS);

    //! -----------------------------
    //! slice volume mesh with plane
    //! -----------------------------
    //bool sliceMeshWithPlane(double a, double b, double c, double d,
    //                        occHandle(MeshVS_DataSource) &slicedMeshDS,
    //                        std::vector<std::pair<int,std::string>> &vecCCXFaceDefs);
    bool sliceMeshWithPlane(double a, double b, double c, double d,
                            occHandle(MeshVS_DataSource) &slicedMeshDS,
                            std::vector<std::pair<int,int>> &vecCCXFaceDefs);

private:

    occHandle(Ng_MeshVS_DataSource3D) myVolumeMesh;
};

#endif // BOLTTOOL_H
