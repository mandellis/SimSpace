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
    //cout<<"pointToMeshDistance::pointToMeshDistance()->____contructor called____"<<endl;
    this->init(aMeshDS);
}

//! ---------------
//! function: init
//! details:
//! ---------------
void pointToMeshDistance::init(const occHandle(Ng_MeshVS_DataSourceFace) &aMeshDS)
{
    if(aMeshDS.IsNull()) return;

    //cout<<"pointToMeshDistance::init()->____function called____"<<endl;
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    int NN = aMeshDS->GetAllNodes().Extent();
    V.resize(NN,3);

    TColStd_MapIteratorOfPackedMapOfInteger it;
    int NbNodes;
    MeshVS_EntityType type;
    double buf[3];
    TColStd_Array1OfReal coords(*buf,1,3);
    for(it.Initialize(aMeshDS->GetAllNodes());it.More();it.Next())
    {
        int globalID = it.Key();
        int localNodeID = aMeshDS->myNodesMap.FindIndex(globalID);
        aMeshDS->GetGeom(globalID,false,coords,NbNodes,type);
        int pos = localNodeID-1;
        V(pos,0) = coords(1);
        V(pos,1) = coords(2);
        V(pos,2) = coords(3);
    }

    //! ----------------
    //! adding elements
    //! ----------------
    int NE = aMeshDS->GetAllElements().Extent();
    F.resize(NE,3);

    int bufn[3];
    TColStd_Array1OfInteger nodeIDs(*bufn,1,3);     //! only triangles
    for(it.Initialize(aMeshDS->GetAllElements());it.More();it.Next())
    {
        int globalElementID = it.Key();
        int localElementID = aMeshDS->myElementsMap.FindIndex(globalElementID);
        aMeshDS->GetNodesByElement(globalElementID,nodeIDs,NbNodes);

        int pos = localElementID-1;
        for(int i=1; i<=NbNodes; i++)
        {
            int curGlobalNodeID = nodeIDs.Value(i);
            int curLocalNodeID = aMeshDS->myNodesMap.FindIndex(curGlobalNodeID);
            F(pos,i-1) = curLocalNodeID-1;
        }
    }

    //! -------------------------------
    //! store the mesh and init embree
    //! -------------------------------
    myV = V;
    myF = F;
    myEmbree.init(V.cast<float>(),F.cast<int>());
}

//! -------------------------
//! function: signedDistance
//! details:
//! -------------------------
bool pointToMeshDistance::distance(double *P, double *dir, float *d)
{
    Eigen::RowVector3f origin (P[0],P[1],P[2]);
    Eigen::RowVector3f direction(dir[0],dir[1],dir[2]);
    igl::Hit hit;
    double tnear = 1e-3;
    bool isDone = myEmbree.intersectRay(origin.cast<float>(),direction.cast<float>(),hit,tnear);
    if(isDone==false) *d = 0.0;
    else *d = hit.t;
    return isDone;
}
