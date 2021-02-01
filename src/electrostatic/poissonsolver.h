#ifndef POISSONSOLVER_H
#define POISSONSOLVER_H

//! ------
//! Eigen
//! ------
#include <ext/eigen/Eigen/Dense>
#include <ext/eigen/Eigen/Core>
#include <ext/eigen/Eigen/Sparse>

//! ---
//! Qt
//! ---
#include <QObject>

//! ----
//! C++
//! ----
#include <memory>
#include <vector>
#include <map>
#include <iostream>
using namespace std;

//! ------------
//! opencascade
//! ------------
#include "occhandle.h"

class Ng_MeshVS_DataSource3D;
class Ng_MeshVS_DataSourceFace;

class PoissonSolver: public QObject
{
    Q_OBJECT

public:

    //! constructor: initialize with an occ volume mesh
    PoissonSolver(const occHandle(Ng_MeshVS_DataSource3D) &occVolumeMesh, QObject* parent=0);

    //! define patches
    void definePatches(const std::map<int,occHandle(Ng_MeshVS_DataSourceFace)> &patches);

    //! init potential on boundary
    void initPotentialOnBoundary(double val=0.0);

    //! set boundary condition/conditions
    void setBoundaryConditions(const std::map<int,double> &facesPotential, bool reinit);
    void setBoundaryCondition(int faceNr, double val, bool reinit);

    //! solve
    void solve(Eigen::VectorXd &Z, const Eigen::VectorXd &density);

private:

    //! mesh of the domain
    occHandle(Ng_MeshVS_DataSource3D) myDomainMesh;

    Eigen::VectorXd phi;    // electrostatic potential
    Eigen::MatrixXd V;      // vertexes coords
    Eigen::MatrixXi T;      // tetrahedra definitions
    Eigen::VectorXi all;    // all indexes
    Eigen::VectorXi in;     // interior indexes
    Eigen::VectorXi b;      // surface boundary indexes
    Eigen::VectorXd bc;     // boundary condition values

    std::map<int,Eigen::VectorXi> myPatches;  // patches
};

#endif // POISSONSOLVER_H
