//! ----------------
//! custom includes
//! ----------------
#include "pointtomeshdistance.h"

//! ------
//! Eigen
//! ------
#include <Eigen/Dense>

//! -------
//! libigl
//! -------
#include <libigl/include/igl/embree/line_mesh_intersection.h>
#include <libigl/include/igl/Hit.h>
#include <libigl/include/igl/per_vertex_normals.h>

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
pointToMeshDistance::pointToMeshDistance(const occHandle(Ng_MeshVS_DataSourceFace) &aMeshDS)
{
    this->init(aMeshDS);
}

//! ---------------
//! function: init
//! details:
//! ---------------
void pointToMeshDistance::init(const occHandle(Ng_MeshVS_DataSourceFace) &aMeshDS)
{
    if(aMeshDS.IsNull()) return;
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    //! --------------------------------------------
    //! adding points to Eigen mesh data structure
    //! -------------------------------------------
    int NN = aMeshDS->GetAllNodes().Extent();
    V.resize(NN,3);
    TColStd_MapIteratorOfPackedMapOfInteger it(aMeshDS->GetAllNodes());
    int NbNodes;
    MeshVS_EntityType aType;
    double buf[3];
    TColStd_Array1OfReal coords(*buf,1,3);
    for(int localNodeID = 1; localNodeID<=NN; localNodeID++, it.Next())
    {
        int globalNodeID = it.Key();
        aMeshDS->GetGeom(globalNodeID,false,coords,NbNodes,aType);
        int pos = localNodeID-1;
        V(pos,0) = coords(1);
        V(pos,1) = coords(2);
        V(pos,2) = coords(3);
    }

    //! ---------------------------------------------
    //! adding elements to Eigen mesh data structure
    //! ---------------------------------------------
    int NE = aMeshDS->GetAllElements().Extent();
    F.resize(NE,3);

    int bufn[3];
    TColStd_Array1OfInteger nodeIDs(*bufn,1,3);
    int NbElements = aMeshDS->GetAllElements().Extent();
    TColStd_MapIteratorOfPackedMapOfInteger itt(aMeshDS->GetAllElements().Extent());
    for(int localElementID = 1; localElementID<=NbElements; localElementID++, itt.Next())
    {
        int globalElementID = itt.Key();
        aMeshDS->GetNodesByElement(globalElementID,nodeIDs,NbNodes);
        int pos = localElementID-1;
        for(int i=1; i<=NbNodes; i++)
        {
            int curGlobalNodeID = nodeIDs.Value(i);
            int curLocalNodeID = aMeshDS->myNodesMap.FindIndex(curGlobalNodeID);
            F(pos,i-1) = curLocalNodeID-1;
        }
    }

    //! ------------
    //! init embree
    //! ------------
    myEmbree.init(V.cast<float>(),F.cast<int>());
}

//! -------------------------
//! function: signedDistance
//! details:
//! -------------------------
bool pointToMeshDistance::distance(double *P, double *dir, double *d)
{
    Eigen::MatrixXd V_source; V_source.resize(1,3);
    Eigen::MatrixXd N_source; N_source.resize(1,3);
    V_source(0,0) = P[0]; V_source(0,1) = P[1]; V_source(0,2) = P[2];
    N_source(0,0) = dir[0]; N_source(0,1) = dir[1]; N_source(0,2) = dir[2];

    Eigen::MatrixXd intersectionPoint = this->line_mesh_intersection(V_source,N_source);
    double xi = intersectionPoint(0);
    double yi = intersectionPoint(1);
    double zi = intersectionPoint(2);
    if(xi == 0 && yi == 0 && zi == -1) return false;

    *d = sqrt(pow(xi-P[0],2)+pow(yi-P[1],2)+pow(zi-P[2],2));
    return true;
}

//! ---------------------------------------------------------
//! function: line_mesh_intersection
//! details:  it works on multiple (source, direction) pairs
//! ---------------------------------------------------------
Eigen::MatrixXd pointToMeshDistance::line_mesh_intersection(const Eigen::MatrixXd & V_source, const Eigen::MatrixXd  & N_source)
{
  const double tol = 0.00001;

  Eigen::MatrixXd ray_pos = V_source;
  Eigen::MatrixXd ray_dir = N_source;

  // Allocate matrix for the result
  Eigen::MatrixXd R;
  R.resize(V_source.rows(), 3);

  // Shoot rays from the source to the target
  for (unsigned i=0; i<ray_pos.rows(); ++i)
  {
    igl::Hit A,B;

    // Shoot ray A
    Eigen::RowVector3d A_pos = ray_pos.row(i) + tol * ray_dir.row(i);
    Eigen::RowVector3d A_dir = -ray_dir.row(i);

    bool A_hit = myEmbree.intersectBeam(A_pos.cast<float>(), A_dir.cast<float>(),A);

    Eigen::RowVector3d B_pos = ray_pos.row(i) - tol * ray_dir.row(i);
    Eigen::RowVector3d B_dir = ray_dir.row(i);

    bool B_hit = myEmbree.intersectBeam(B_pos.cast<float>(), B_dir.cast<float>(),B);

    int choice = -1;

    if (A_hit && ! B_hit) choice = 0;
    else if (!A_hit && B_hit) choice = 1;
    else if (A_hit && B_hit) choice = A.t > B.t;

    Eigen::RowVector3d temp;

    if (choice == -1) temp << -1, 0, 0;
    else if (choice == 0) temp << A.id, A.u, A.v;
    else if (choice == 1) temp << B.id, B.u, B.v;

    R.row(i) = temp;
  }
  return R;
}
