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

    //! --------------------------------------------
    //! adding points to Eigen mesh data structure
    //! -------------------------------------------
    int NN = aMeshDS->GetAllNodes().Extent();
    //cout<<"____"<<NN<<"____"<<endl;
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
    //cout<<"____"<<NE<<"____"<<endl;
    F.resize(NE,3);

    int bufn[3];
    TColStd_Array1OfInteger nodeIDs(*bufn,1,3);
    int NbElements = aMeshDS->GetAllElements().Extent();
    TColStd_MapIteratorOfPackedMapOfInteger itt(aMeshDS->GetAllElements());
    int mask[3] {1,2,3};
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
            //F(pos,mask[i-1]-1) = curLocalNodeID-1;
        }
    }

    //! ------------
    //! init embree
    //! ------------
    //cout<<"pointToMeshDistance::init()->____initializing embree____"<<endl;
    myEmbree.init(V.cast<float>(),F.cast<int>());
    //cout<<"pointToMeshDistance::init()->____embree initialized____"<<endl;
}

//! -------------------------
//! function: signedDistance
//! details:
//! -------------------------
bool pointToMeshDistance::distance(double *P, double *dir, double *d)
{
    //cout<<"pointToMeshDistance::distance()->____function called____"<<endl;
    Eigen::MatrixXd V_source;
    V_source.resize(1,3);
    V_source(0,0) = P[0]; V_source(0,1) = P[1]; V_source(0,2) = P[2];

    Eigen::MatrixXd N_source;
    N_source.resize(1,3);
    N_source(0,0) = dir[0]; N_source(0,1) = dir[1]; N_source(0,2) = dir[2];

    Eigen::MatrixXd ray_pos = V_source;
    Eigen::MatrixXd ray_dir = N_source;
    std::vector<igl::Hit> hits;
    int num_rays;
    bool isEmpty = myEmbree.intersectRay(ray_pos.cast<float>(),ray_dir.cast<float>(),hits,num_rays);
    //cout<<"____intersect? "<<(isEmpty==true? "N":"Y")<<"____"<<endl;
    //cout<<"____num rays: "<<num_rays<<"____"<<endl;
    //cout<<"____hits nb:  "<<hits.size()<<"____"<<endl;

    if(isEmpty==false && hits.size()>=2)
    {
        std::sort(hits.begin(),hits.end(),[](igl::Hit a, igl::Hit b) {return a.t > b.t;});
        //const igl::Hit &theHit = hits[1];
        //cout<<"____"<<theHit.t<<"____"<<endl;
        *d = hits[1].t;
        return true;
    }
    return false;
}
