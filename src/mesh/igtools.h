#ifndef IGLTOOLS_H
#define IGLTOOLS_H

//! C++
#include <vector>

//! redefinition of opencascade Handle() macro
#include "occhandle.h"

//! Eigen
#include <Eigen/Dense>

//! custom includes
#include <ng_meshvs_datasourceface.h>

class iglTools
{
public:

    iglTools(){;}

    //! -----------------
    //! OCCMeshToIglMesh
    //! -----------------
    static void OCCMeshToIglMesh(const occHandle(Ng_MeshVS_DataSourceFace) &aFaceMesh,
                                 Eigen::MatrixXd &V,
                                 Eigen::MatrixXi &F);

    //! -----------------
    //! intersect meshes
    //! -----------------
    static bool iglIntersectMeshes(const occHandle(Ng_MeshVS_DataSourceFace) &aFaceMesh1,
                                   const occHandle(Ng_MeshVS_DataSourceFace) &aFaceMesh2,
                                   std::vector<int> &badNodes1,
                                   std::vector<int> &badNodes2);

    //! -------------------------
    //! check self intersections
    //! -------------------------
    static bool iglCheckSelfIntersection(const occHandle(Ng_MeshVS_DataSourceFace) &aFaceMesh,
                                         std::vector<int> &badNodes);

    //! ----------------------
    //! smooth scalar on mesh
    //! ----------------------
    static void smoothScalarOnMesh(const occHandle(Ng_MeshVS_DataSourceFace) &aFaceMesh,
                                   QMap<int,double> &scalarField,
                                   int mode = 1);
};

#endif // IGLTOOLS_H
