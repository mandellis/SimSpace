#ifndef MESHSLICER_H
#define MESHSLICER_H

//#define TEST_SLICER

//! ----------------
//! custom includes
//! ----------------
#include "occhandle.h"
#include <ng_meshvs_datasource3d.h>

//! ------
//! Eigen
//! ------
#include <Eigen/Dense>

//! ----
//! igl
//! ----
//#include <libigl/include/igl/copyleft/cgal/intersect_with_half_space.h>

class meshSlicer
{
public:

    //! constructor
    meshSlicer(const occHandle(Ng_MeshVS_DataSource3D) &aVolumeMeshDS = occHandle(Ng_MeshVS_DataSource3D)());

    //! set mesh data source
    void setMeshDataSource(const occHandle(Ng_MeshVS_DataSource3D) &aVolumeMeshDS);

    //! perform
    bool perform(double a, double b, double c, double d, occHandle(Ng_MeshVS_DataSource3D) &slicedVolumeMeshDS);

private:

    occHandle(Ng_MeshVS_DataSource3D) myVolumeMeshDS;

#ifdef TEST_SLICER
    Eigen::MatrixXd V;
    Eigen::MatrixXi T;
#endif
};

#endif // MESHSLICER_H
