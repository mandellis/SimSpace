#ifndef POINTTOMESHDISTANCE_H
#define POINTTOMESHDISTANCE_H

#include "occhandle.h"
#include <ng_meshvs_datasourceface.h>
#include <libigl/include/igl/embree/EmbreeIntersector.h>

//! ------
//! Eigen
//! ------
#include <Eigen/Dense>

//! ----
//! C++
//! ----
#include <iostream>

using namespace std;

class pointToMeshDistance
{
public:

    pointToMeshDistance(const occHandle(Ng_MeshVS_DataSourceFace) &aMeshDS = occHandle(Ng_MeshVS_DataSourceFace)());

    void init(const occHandle(Ng_MeshVS_DataSourceFace) &aMeshDS);

    bool distance(double *P, double *dir, float *d);

private:

    igl::embree::EmbreeIntersector myEmbree;
    Eigen::MatrixXd myV;
    Eigen::MatrixXi myF;
};

#endif // POINTTOMESHDISTANCE_H
