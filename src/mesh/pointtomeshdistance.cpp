//! ----------------
//! custom includes
//! ----------------
#include "pointtomeshdistance.h"

//! ------
//! Eigen
//! ------
#include <ext/eigen/Eigen/Dense>

//! -------
//! libigl
//! -------
#include <ext/libigl/include/igl/embree/line_mesh_intersection.h>
#include <ext/libigl/include/igl/Hit.h>

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

    cout<<"pointToMeshDistance::init()->____function called____"<<endl;
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

//! -------------------
//! function: distance
//! details:
//! -------------------
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

//! -------------------------
//! function: gapMeasurement
//! details:
//! -------------------------
bool pointToMeshDistance::distance(double *P, double *dir, double firstConeAmplitude, int NbCones, float *d)
{
    const double PI = 3.1415926534;

    //! ---------------------------------------------
    //! vector of gaps measured along different rays
    //! ---------------------------------------------
    std::vector<float> gaps;

    Eigen::RowVector3f origin (P[0],P[1],P[2]);
    Eigen::RowVector3f direction(dir[0],dir[1],dir[2]);
    igl::Hit hit;
    double tnear = 1e-3;
    bool isDone = myEmbree.intersectRay(origin.cast<float>(),direction.cast<float>(),hit,tnear);
    if(isDone==true) gaps.push_back(hit.t);

    //! ---------------------------------------------------------------------
    //! create a start direction which will be rotated along the node normal
    //! This vector will be rotated by steps around the node normal, creating
    //! a cone of rays
    //! ---------------------------------------------------------------------
    double angleDeviation = firstConeAmplitude;
    for(int pass = 1; pass<=NbCones; pass++)
    {
        double v_0[3];
        if(dir[2]!=0)
        {
            v_0[0] = dir[0];
            v_0[1] = dir[1];
            v_0[2] = (cos(pass*angleDeviation)-dir[0]*dir[0]-dir[1]*dir[1])/dir[2];
        }
        else if(dir[1]!=0)
        {
            v_0[0] = (cos(pass*angleDeviation)-dir[0]*dir[0]-dir[2]*dir[2])/dir[1];
            v_0[1] = dir[1];
            v_0[2] = dir[2];
        }
        else if(dir[0]!=0)
        {
            v_0[0] = (cos(pass*angleDeviation)-dir[1]*dir[1]-dir[2]*dir[2])/dir[0];
            v_0[1] = dir[1];
            v_0[2] = dir[2];
        }

        double l = sqrt(pow(v_0[0],2)+pow(v_0[1],2)+pow(v_0[2],2));
        v_0[0] /= l;
        v_0[1] /= l;
        v_0[2] /= l;

        int N = 16;                     // cast 16 rays from P
        double dangle = 2*PI/N;
        for(int i = 0; i<N; i++)
        {
            double dir_rot[3];
            this->rotateVector(v_0,dir,dangle,dir_rot);
            Eigen::RowVector3f dir_rotated(dir_rot[0],dir_rot[1],dir_rot[2]);
            isDone = myEmbree.intersectRay(origin.cast<float>(),dir_rotated.cast<float>(),hit,tnear);
            if(isDone==true) gaps.push_back(hit.t);
            v_0[0] = dir_rot[0];
            v_0[1] = dir_rot[1];
            v_0[2] = dir_rot[2];
        }
    }
    if(gaps.empty()==true)
    {
        *d = 0;
        return false;
    }
    double gap = *std::min(gaps.begin(),gaps.end());
    *d = gap;
    return true;
}

//! ---------------------------------------------------------------
//! function: rotateVector
//! details: rotate vector v about direction k by an angle "angle"
//! ---------------------------------------------------------------
void pointToMeshDistance::rotateVector(double *v, double *k, double angle, double *r)
{
    //! ---------------------
    //! i       j       k
    //! k[0]    k[1]    k[2]
    //! v[0]    v[1]    v[2]
    //! ---------------------
    double sx = k[1]*v[2]-k[2]*v[1];
    double sy = k[2]*v[0]-k[0]*v[2];
    double sz = k[0]*v[1]-k[1]*v[0];

    double dot = v[0]*k[0]+v[1]*k[1]+v[2]*k[2];
    if(dot>1) dot = 1;
    if(dot<-1) dot = -1;
    r[0] = v[0]*cos(angle)+ sx*sin(angle) + k[0]*dot*(1-cos(angle));
    r[1] = v[1]*cos(angle)+ sy*sin(angle) + k[1]*dot*(1-cos(angle));
    r[2] = v[2]*cos(angle)+ sz*sin(angle) + k[2]*dot*(1-cos(angle));

    //! ----------
    //! normalize
    //! ----------
    double l = sqrt(pow(r[0],2)+pow(r[1],2)+pow(r[2],2));
    r[0] /=l;
    r[1] /=l;
    r[2] /=l;
}
