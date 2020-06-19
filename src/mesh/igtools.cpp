//! ----------------
//! custom includes
//! ----------------
#include "igtools.h"

#include <libigl/include/igl/cotmatrix.h>
#include <libigl/include/igl/massmatrix.h>
#include <libigl/include/igl/hessian_energy.h>
#include <libigl/include/igl/copyleft/cgal/remesh_self_intersections.h>
#include <libigl/include/igl/copyleft/cgal/SelfIntersectMesh.h>
#include <libigl/include/igl/copyleft/cgal/intersect_other.h>
#include "C:/CGAL-4.14/build_MSVC14.0/include/CGAL/compiler_config.h"

//! ----
//! OCC
//! ----
#include <TColStd_PackedMapOfInteger.hxx>
#include <TColStd_MapIteratorOfPackedMapOfInteger.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <MeshVS_EntityType.hxx>

//! ---------------------------------------------------------------------------
//! function: OCCMeshToIglMesh
//! details:  converts the opencascade face mesh data source into Eigen format
//! ---------------------------------------------------------------------------
void iglTools::OCCMeshToIglMesh(const occHandle(Ng_MeshVS_DataSourceFace) &aFaceMesh,
                                Eigen::MatrixXd &V,
                                Eigen::MatrixXi &F)
{
    //! --------------
    //! adding points
    //! --------------
    int NN = aFaceMesh->GetAllNodes().Extent();
    V.resize(NN,3);

    TColStd_MapIteratorOfPackedMapOfInteger it;
    int NbNodes;
    MeshVS_EntityType type;
    double buf[3];
    TColStd_Array1OfReal coords(*buf,1,3);
    for(it.Initialize(aFaceMesh->GetAllNodes());it.More();it.Next())
    {
        int globalID = it.Key();
        int localNodeID = aFaceMesh->myNodesMap.FindIndex(globalID);
        aFaceMesh->GetGeom(globalID,false,coords,NbNodes,type);
        int pos = localNodeID-1;
        V(pos,0) = coords(1);
        V(pos,1) = coords(2);
        V(pos,2) = coords(3);
    }

    //! ----------------
    //! adding elements
    //! ----------------
    int NE = aFaceMesh->GetAllElements().Extent();
    F.resize(NE,3);

    int bufn[3];
    TColStd_Array1OfInteger nodeIDs(*bufn,1,3);     //! only triangles
    for(it.Initialize(aFaceMesh->GetAllElements());it.More();it.Next())
    {
        int globalElementID = it.Key();
        int localElementID = aFaceMesh->myElementsMap.FindIndex(globalElementID);
        aFaceMesh->GetNodesByElement(globalElementID,nodeIDs,NbNodes);

        int pos = localElementID-1;
        for(int i=1; i<=NbNodes; i++)
        {
            int curGlobalNodeID = nodeIDs.Value(i);
            int curLocalNodeID = aFaceMesh->myNodesMap.FindIndex(curGlobalNodeID);
            F(pos,i-1) = curLocalNodeID-1;
        }
    }
}

//! -----------------------------------------------------------------------
//! function: intersectMeshes
//! details:  interface to INTERSECT_OTHER - see the following doc
//!
//! INTERSECT_OTHER Given a triangle mesh (VA,FA) and another mesh (VB,FB)
//! find all pairs of intersecting faces. Note that self-intersections are
//! ignored.
//!
//! Inputs:
//!
//!   VA  #V by 3 list of vertex positions
//!   FA  #F by 3 list of triangle indices into VA
//!   VB  #V by 3 list of vertex positions
//!   FB  #F by 3 list of triangle indices into VB
//!   params   whether to detect only and then whether to only find first
//!   intersection
//!
//! Outputs:
//!
//!   IF  #intersecting face pairs by 2 list of intersecting face pairs,
//!   indexing FA and FB
//!   VVAB  #VVAB by 3 list of vertex positions
//!   FFAB  #FFAB by 3 list of triangle indices into VVA
//!   JAB  #FFAB list of indices into [FA;FB] denoting birth triangle
//!   IMAB  #VVAB list of indices stitching duplicates (resulting from
//!     mesh intersections) together
//! -----------------------------------------------------------------------
bool iglTools::iglIntersectMeshes(const occHandle(Ng_MeshVS_DataSourceFace) &aFaceMesh1,
                                  const occHandle(Ng_MeshVS_DataSourceFace) &aFaceMesh2,
                                  std::vector<int> &badNodes1,
                                  std::vector<int> &badNodes2)
{
    //! ----------------------------------------------------
    //! first mesh - vertices coords, triangles definitions
    //! ----------------------------------------------------
    Eigen::MatrixXd V1;
    Eigen::MatrixXi F1;
    iglTools::OCCMeshToIglMesh(aFaceMesh1,V1,F1);

    //! -----------------------------------------------------
    //! second mesh - vertices coords, triangles definitions
    //! -----------------------------------------------------
    Eigen::MatrixXd V2;
    Eigen::MatrixXi F2;
    iglTools::OCCMeshToIglMesh(aFaceMesh2,V2,F2);

    Eigen::MatrixXi IF; //! intersecting faces
    bool first_only = false;

    //! --------------
    //! static method
    //! --------------
    igl::copyleft::cgal::intersect_other(V1,F1,V2,F2,first_only,IF);
    if(IF.rows()==0) return false;

    //! ------------------------------------------------------------------------
    //! else process data: retrieve the bad nodes of the first and second mesh
    //! ------------------------------------------------------------------------
    std::vector<int> bbadNodes1, bbadNodes2;
    for(int i=0; i<IF.rows(); i++)
    {
        bbadNodes1.push_back(F1(IF(i,0),0));
        bbadNodes1.push_back(F1(IF(i,0),1));
        bbadNodes1.push_back(F1(IF(i,0),2));
        bbadNodes2.push_back(F2(IF(i,1),0));
        bbadNodes2.push_back(F2(IF(i,1),1));
        bbadNodes2.push_back(F2(IF(i,1),2));
    }

    //! -------------------------
    //! remove duplicated values
    //! -------------------------
    std::set<int> B1(bbadNodes1.begin(), bbadNodes1.end());
    for(int n=0; n<B1.size(); n++)
    {
        int x = *std::next(B1.begin(), n);
        badNodes1.push_back(x);
    }

    std::set<int> B2(bbadNodes2.begin(), bbadNodes2.end());
    for(int n=0; n<B2.size(); n++)
    {
        int x = *std::next(B2.begin(), n);
        badNodes2.push_back(x);
    }
    return true;
}

//! ------------------------------------------------------------------------------
//! function: iglCheckSelfIntersection
//! details:  interface to SelfIntersectMesh
//!           return the number of intersecting pairs
//!
//! uses igl::copyleft::cgal::remesh_self_intersections(V,F,params,VV,FF,IF,H,IM)
//!
//! Given a triangle mesh (V,F) compute a new mesh (VV,FF) which is the same
//! as (V,F) except that any self-intersecting triangles in (V,F) have been
//! subdivided (new vertices and face created) so that the self-intersection
//! contour lies exactly on edges in (VV,FF). New vertices will appear in
//! original faces or on original edges. New vertices on edges are "merged"
//! only across original faces sharing that edge. This means that if the input
//! triangle mesh is a closed manifold the output will be too.
//!
//! Inputs:
//!   V  #V by 3 list of vertex positions
//!   F  #F by 3 list of triangle indices into V
//!   params  struct of optional parameters
//! Outputs:
//!   VV  #VV by 3 list of vertex positions
//!   FF  #FF by 3 list of triangle indices into VV
//!   IF  #intersecting face pairs by 2 list of intersecting face pairs,
//!     indexing F
//!   J  #FF list of indices into F denoting birth triangle
//!   IM  #VV list of indices into VV of unique vertices.
//! ------------------------------------------------------------------------------
//#include <libigl/include/igl/writeOBJ.h>
//#include <libigl/include/igl/writeSTL.h>
bool iglTools::iglCheckSelfIntersection(const occHandle(Ng_MeshVS_DataSourceFace) &aFaceMesh, std::vector<int> &badNodes)
{
    cout<<"iglTools::iglCheckSelfIntersection()->____function called____"<<endl;

    //! ---------------------------------------------
    //! the mesh to be checked for self intersection
    //! ---------------------------------------------
    Eigen::MatrixXd V(1,3);
    Eigen::MatrixXi F(1,3);
    iglTools::OCCMeshToIglMesh(aFaceMesh,V,F);

    Eigen::MatrixXd SV;
    Eigen::MatrixXi SF,IF;
    Eigen::VectorXi J,IM;

    igl::copyleft::cgal::RemeshSelfIntersectionsParam params;
    params.detect_only = true;

    igl::copyleft::cgal::remesh_self_intersections(V,F,params,SV,SF,IF,J,IM);

    //! ---------------------------
    //! print the "birth" triangle
    //! ---------------------------
    //cout<<"____Dimension of vector containing birth triangles____"<<endl;
    //cout<<"____("<<J.rows()<<", "<<J.cols()<<")____"<<endl;

    //! --------------------------------
    //! collect the vector of bad nodes
    //! --------------------------------
    //!cout<<"iglTools::iglCheckSelfIntersection()->____Number of bad triangles pairs: "<<IF.rows()<<"____"<<endl;

    std::vector<int> bbadNodes;
    for(int i=0; i<IF.rows(); i++)
    {
        //!cout<<"iglTools::iglCheckSelfIntersection()->_____bad triangle pair("<<IF(i,0)<<", "<<IF(i,1)<<")_____"<<endl;
        bbadNodes.push_back(F(IF(i,0),0));
        bbadNodes.push_back(F(IF(i,0),1));
        bbadNodes.push_back(F(IF(i,0),2));
        bbadNodes.push_back(F(IF(i,1),0));
        bbadNodes.push_back(F(IF(i,1),1));
        bbadNodes.push_back(F(IF(i,1),2));
    }

    //! -------------------------
    //! remove duplicated values
    //! -------------------------
    std::set<int> B(bbadNodes.begin(), bbadNodes.end());
    for(int n=0; n<B.size(); n++)
    {
        int x = *std::next(B.begin(), n);
        badNodes.push_back(x);
    }

    if(IF.rows()>=1) return true;   //! at least one intersection
    return false;
}

//! -------------------------------------------------------------------------------------------
//! function: smoothScalarOnMesh
//! details:  https://github.com/libigl/libigl/blob/master/tutorial/712_DataSmoothing/main.cpp
//!           https://libigl.github.io/tutorial/#data-smoothing
//! -------------------------------------------------------------------------------------------
void iglTools::smoothScalarOnMesh(const occHandle(Ng_MeshVS_DataSourceFace) &aFaceMesh,
                                  QMap<int,double> &scalarField,
                                  int mode)
{
    typedef Eigen::SparseMatrix<double> SparseMat;

    //! --------------------------
    //! convert into Eigen format
    //! --------------------------
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    //iglTools::OCCMeshToIglMesh(aFaceMesh,V,F);
    //! --------------
    //! adding points
    //! --------------
    int NN = aFaceMesh->GetAllNodes().Extent();
    V.resize(NN,3);

    TColStd_MapIteratorOfPackedMapOfInteger it;
    int NbNodes;
    MeshVS_EntityType type;
    double buf[3];
    TColStd_Array1OfReal coords(*buf,1,3);
    for(it.Initialize(aFaceMesh->GetAllNodes());it.More();it.Next())
    {
        int globalID = it.Key();
        int localNodeID = aFaceMesh->myNodesMap.FindIndex(globalID);
        aFaceMesh->GetGeom(globalID,false,coords,NbNodes,type);
        int pos = localNodeID-1;
        V(pos,0) = coords(1);
        V(pos,1) = coords(2);
        V(pos,2) = coords(3);
    }

    //! ----------------
    //! adding elements
    //! ----------------
    int NE = aFaceMesh->GetAllElements().Extent();
    F.resize(NE,3);

    int bufn[3];
    TColStd_Array1OfInteger nodeIDs(*bufn,1,3);     //! only triangles
    for(it.Initialize(aFaceMesh->GetAllElements());it.More();it.Next())
    {
        int globalElementID = it.Key();
        int localElementID = aFaceMesh->myElementsMap.FindIndex(globalElementID);
        aFaceMesh->GetNodesByElement(globalElementID,nodeIDs,NbNodes);

        int pos = localElementID-1;
        for(int i=1; i<=NbNodes; i++)
        {
            int curGlobalNodeID = nodeIDs.Value(i);
            int curLocalNodeID = aFaceMesh->myNodesMap.FindIndex(curGlobalNodeID);
            F(pos,i-1) = curLocalNodeID-1;
        }
    }

    //! -------------------------
    //! convert the scalar field
    //! -------------------------
    Eigen::VectorXd znoisy(scalarField.size());
    for(QMap<int,double>::iterator it = scalarField.begin(); it!=scalarField.end(); it++)
    {
        int globalNodeID = it.key();
        double value = it.value();
        znoisy(globalNodeID-1) = value;
    }

    //! -------------------------
    //! convert the scalar field
    //! -------------------------
    Eigen::VectorXd *zsmoothed;
    switch(mode)
    {
    case 0:
    {
        //! ------------------------------------------
        //! Constructing the squared Laplacian energy
        //! ------------------------------------------
        SparseMat L, M;
        igl::cotmatrix(V, F, L);
        igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);
        Eigen::SimplicialLDLT<SparseMat> solver(M);
        SparseMat MinvL = solver.solve(L);
        SparseMat QL = L.transpose()*MinvL;

        //! ---------------------------------
        //! Solve to find Laplacian-smoothed
        //! ---------------------------------
        const double al = 8e-4;
        Eigen::SimplicialLDLT<SparseMat> lapSolver(al*QL + (1.-al)*M);
        Eigen::VectorXd zl = lapSolver.solve(al*M*znoisy);
        zsmoothed = &zl;
    }
        break;

    case 1:
    {
        //! ----------------------------------------
        //! Constructing the squared Hessian energy
        //! ----------------------------------------
        SparseMat M;
        igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);
        Eigen::SimplicialLDLT<SparseMat> solver(M);

        SparseMat QH;
        igl::hessian_energy(V, F, QH);
        const double ah = 0.05;
        /*
        Eigen::SimplicialLDLT<SparseMat> hessSolver(ah*QH + (1.-ah)*M);
        Eigen::VectorXd zh;
        //! -----------------------------------------
        //! Solve to find Hessian-smoothed solutions
        //! -----------------------------------------
        for(int n=1; n<=10; n++)
        {
            cout<<"____smoothing step: "<<n<<"____"<<endl;
            zh = hessSolver.solve(ah*M*znoisy);
            znoisy = zh;
        }
        zsmoothed = &zh;
        */
        //! -----------------------------------------
        //! Solve to find Hessian-smoothed solutions
        //! -----------------------------------------
        Eigen::SimplicialLDLT<SparseMat> hessSolver(ah*QH + (1.-ah)*M);
        Eigen::VectorXd zh = hessSolver.solve(ah*M*znoisy);
        zsmoothed = &zh;
    }
        break;
    }

    for(int index = 0; index<zsmoothed->rows(); index++)
    {
        int globalNodeID = index+1;
        double value = (*zsmoothed)(index);
        scalarField.insert(globalNodeID,value);
    }
}
