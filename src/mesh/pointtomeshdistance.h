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

    bool distance(double *P, double *dir, double *d);

private:

    Eigen::MatrixXd line_mesh_intersection(const Eigen::MatrixXd & V_source, const Eigen::MatrixXd  & N_source);

private:

    //Eigen::MatrixXd V;
    //Eigen::MatrixXi F;
    igl::embree::EmbreeIntersector myEmbree;
};

#endif // POINTTOMESHDISTANCE_H
