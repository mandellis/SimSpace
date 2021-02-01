//! ----------------
//! custom includes
//! ----------------
#include "poissonsolver.h"
#include "ng_meshvs_datasource3d.h"
#include "ng_meshvs_datasource2d.h"
#include "ng_meshvs_datasourceface.h"

//! ------------
//! opencascade
//! ------------
#include <TColStd_MapIteratorOfMapOfInteger.hxx>
#include <TColStd_Array1OfInteger.hxx>
#include <TColStd_IndexedMapOfInteger.hxx>

//! ----
//! igl
//! ----
#include <ext/libigl/include/igl/cotmatrix.h>
#include <ext/libigl/include/igl/massmatrix.h>
#include <ext/libigl/include/igl/slice.h>
#include <ext/libigl/include/igl/min_quad_with_fixed.h>

//! ---------------------------------------------------------
//! function: constructor
//! details:  initialize with an occ volume mesh data source
//! ---------------------------------------------------------
PoissonSolver::PoissonSolver(const occHandle(Ng_MeshVS_DataSource3D) &occVolumeMesh, QObject *parent):
    myDomainMesh(occVolumeMesh),
    QObject(parent)
{
    //! ------------------------------
    //! all interior vertices indexes
    //! ------------------------------
    int NbPoints = occVolumeMesh->GetAllNodes().Extent();
    V.resize(NbPoints,3);
    for(int row=1; row<=NbPoints; row++)
    {
        const std::vector<double> &P = occVolumeMesh->getNodeCoordinates(row);
        for(int col=1; col<=3; col++)
        {
            V(row-1,col-1) = P[col];
            all(row-1) = row-1;
        }
    }

    //! ---------------------------
    //! definition of the elements
    //! ---------------------------
    int NbElements = occVolumeMesh->GetAllElements().Extent();
    T.resize(NbElements,4);
    int NbNodes = -1, buf[4];
    TColStd_Array1OfInteger nodeIDs(*buf,1,4);
    TColStd_IndexedMapOfInteger nodesMap = occVolumeMesh->myNodesMap;
    for(int i=1; i<=NbElements; i++)
    {
        int globalElementID = occVolumeMesh->myElementsMap.FindIndex(i);
        occVolumeMesh->GetNodesByElement(globalElementID,nodeIDs,NbNodes);

        int row = i-1;
        T(row,0) = nodesMap.FindIndex(nodeIDs(1))-1;
        T(row,1) = nodesMap.FindIndex(nodeIDs(2))-1;
        T(row,2) = nodesMap.FindIndex(nodeIDs(3))-1;
        T(row,3) = nodesMap.FindIndex(nodeIDs(4))-1;
    }

    //! ----------------
    //! surface indexes
    //! ----------------
    if(occVolumeMesh->myFaceToElements.isEmpty()) occVolumeMesh->buildFaceToElementConnectivity();
    occHandle(Ng_MeshVS_DataSource2D) surfaceMesh = new Ng_MeshVS_DataSource2D(occVolumeMesh);
    int NbPointsSurface = surfaceMesh->GetAllNodes().Extent();
    b.resize(NbPointsSurface);
    bc.resize(NbPointsSurface);

    int n=0;
    for(TColStd_MapIteratorOfPackedMapOfInteger it(surfaceMesh->GetAllNodes()); it.More(); it.Next())
    {
        int globalNodeID = it.Key();
        b(n) = occVolumeMesh->myNodesMap.FindIndex(globalNodeID)-1;
        n++;
    }

    //! -----------------
    //! interior indexes
    //! -----------------
    n = 0;
    in.resize(NbPoints-NbPointsSurface);
    for(TColStd_MapIteratorOfPackedMapOfInteger it(occVolumeMesh->GetAllNodes()); it.More(); it.Next())
    {
        int globalNodeID = it.Key();
        if(surfaceMesh->GetAllNodes().Contains(globalNodeID)) continue;
        in(n) = occVolumeMesh->myNodesMap.FindIndex(globalNodeID)-1;
    }
}

//! ---------------------------------
//! function: definePatches
//! details:  call after constructor
//! ---------------------------------
void PoissonSolver::definePatches(const std::map<int,occHandle(Ng_MeshVS_DataSourceFace)> &patches)
{
    int NbInserted = 0;
    std::map<int,int> alreadyInsertedIDs;
    for(std::map<int, occHandle(Ng_MeshVS_DataSourceFace)>::const_iterator it = patches.cbegin(); it !=patches.cend(); it++)
    {
        std::pair<int,occHandle(Ng_MeshVS_DataSourceFace)> apair = *it;
        int faceNr = apair.first;
        const occHandle(Ng_MeshVS_DataSourceFace) &faceMeshDS = apair.second;
        std::pair<int,Eigen::VectorXi> aPatch;
        aPatch.first = faceNr;
        Eigen::VectorXi ids;

        std::vector<int> eigenIDsForPatch;
        for(TColStd_MapIteratorOfPackedMapOfInteger it(faceMeshDS->GetAllNodes()); it.More(); it.Next())
        {
            int globalNodeID = it.Key();
            int candidateID = myDomainMesh->myNodesMap.FindIndex(globalNodeID)-1;
            if(alreadyInsertedIDs.count(candidateID)!=0) continue;
            eigenIDsForPatch.push_back(candidateID);
            alreadyInsertedIDs.insert(std::make_pair(NbInserted,candidateID));
            NbInserted++;
        }
        int NbIDsForPatch = int(eigenIDsForPatch.size());
        ids.resize(NbIDsForPatch);
        for(int i=0; i<NbIDsForPatch; i++) ids(i) = eigenIDsForPatch[i];

        aPatch.second = ids;
        myPatches.insert(aPatch);
    }
}

//! ----------------
//! function: solve
//! details:
//! ----------------
void PoissonSolver::solve(Eigen::VectorXd &Z, const Eigen::VectorXd &density)
{
    const double e0 = 8.854e-12;

    //! ---------------------------------------------------
    //! Construct and slice up Laplace - Beltrami operator
    //! ---------------------------------------------------
    Eigen::SparseMatrix<double> L;
    igl::cotmatrix(V,T,L);

    //! -------------------------------------------
    //! Inner points - could be put out from here?
    //! -------------------------------------------
    Eigen::SparseMatrix<double> L_in_in,L_in_b;
    igl::slice(L,in,in,L_in_in);
    igl::slice(L,in,b,L_in_b);

    //! ----------
    //! Solve PDE
    //! ----------
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> psolver(-L_in_in);
    //Eigen::VectorXd Z_in = psolver.solve(L_in_b*bc+density/e0);
    Z = psolver.solve(L_in_b*bc+density/e0);
}

//! --------------------------------
//! function: setBoundaryConditions
//! details:
//! --------------------------------
void PoissonSolver::setBoundaryConditions(const std::map<int,double> &facesPotential, bool reinit)
{
    if(reinit) this->initPotentialOnBoundary(0.0);
    for(std::map<int,double>::const_iterator it = facesPotential.cbegin(); it!=facesPotential.cend(); it++)
    {
        const std::pair<int,double> &apair = *it;
        int faceNr = apair.first;
        double val = apair.second;
        const Eigen::VectorXi &b_ = myPatches.at(faceNr);
        for(int i=0; i<b_.rows(); i++) bc(b_(i))=val;
    }
}

//! -------------------------------
//! function: setBoundaryCondition
//! details:
//! -------------------------------
void PoissonSolver::setBoundaryCondition(int faceNr, double val, bool reinit)
{
    if(reinit) this->initPotentialOnBoundary(0.0);
    const Eigen::VectorXi &b_ = myPatches.at(faceNr);
    for(int i=0; i<b_.rows(); i++) bc(b_(i))=val;
}

//! ----------------------------------
//! function: initPotentialOnBoundary
//! details:
//! ----------------------------------
void PoissonSolver::initPotentialOnBoundary(double val)
{
    for(int i=0; i<bc.rows(); i++) bc(i) = val;
}
