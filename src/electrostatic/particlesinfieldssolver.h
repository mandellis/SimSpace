#ifndef PARTICLESINFIELDSSOLVER_H
#define PARTICLESINFIELDSSOLVER_H

//! ----------------
//! custom includes
//! ----------------
#include "particlesemitter.h"

//! ----
//! OCC
//! ----
#include <TopoDS_Face.hxx>

//! ----
//! C++
//! ----
#include <memory>
#include <utility>
#include <vector>
#include <iostream>
using namespace std;

//! ---
//! Qt
//! ---
#include <QObject>
#include <QStandardItem>
#include <QModelIndex>
#include <QStandardItemModel>
#include <QTreeView>


//! ------
//! Eigen
//! ------
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Sparse>

//! -------
//! libigl
//! -------
#include <libigl/include/igl/AABB.h>

//! ----
//! C++
//! ----
#include <string>
#include <iostream>
using namespace std;

class simulationDataBase;
class Ng_MeshVS_DataSource3D;
class PoissonSolver;

class particlesInFieldsSolver: public QObject
{
    Q_OBJECT

public:

    //! constructor
    particlesInFieldsSolver(simulationDataBase *sDB, QStandardItem *simulationRoot, QObject *parent=0);

    //! init
    bool init(QStandardItem *simulationRoot);

    //! create a particle emitter
    void createEmitter(const TopoDS_Face &aFace, const occHandle(Ng_MeshVS_DataSourceFace) &aFaceMesh, double intensity);

    //! run
    void run();

private:

    simulationDataBase *mySDB;
    QStandardItem *mySimulationRoot;
    Eigen::MatrixXi T;                          // mesh elements
    Eigen::MatrixXd V;                          // mesh points
    Eigen::SparseMatrix<double> M;              // mass matrix
    Eigen::SparseMatrix<double> Minv;           // inverse of
    Eigen::SparseMatrix<double> G;              // gradient operator
    igl::AABB<Eigen::MatrixXd,3> mySearchTree;  // search tree for point location

    Eigen::VectorXd myVolumes;                  // tetrahedra volumes
    Eigen::VectorXd myDensity;                  // density at node
    Eigen::VectorXd myElectrostaticPotential;   // electrostatic potential
    Eigen::MatrixXd myElectricField;            // electric field

    //std::shared_ptr<PoissonSolver> myPoissonSolver;
    PoissonSolver *myPoissonSolver;
    std::vector<shared_ptr<particlesEmitter>> myEmitters;

private:

    int myNbParticles;
    std::map<int,shared_ptr<particle>> myParticles;
    double myFinalTime;
    double myTimeStep;

private:

    bool getCellOfAPoint(double x, double y, double z, int& tetNumber);
    void computeChargeDensityAtNodes();
    void computeElectricField();
    void getElectricFieldAtPoint(double x, double y, double z, double &Ex, double &Ez, double &Ey);
    void moveParticles();
    std::string getTimeStamp();
};

#endif // PARTICLESINFIELDSSOLVER_H