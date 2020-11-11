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
    //int mask[3] {1,2,3};
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

    myV = V;
    myF = F;

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
    //FILE *fp = fopen("D:/box.txt","a");
    //cout<<"pointToMeshDistance::distance()->____function called____"<<endl;
    Eigen::RowVector3f origin (P[0],P[1],P[2]);
    Eigen::RowVector3f direction(dir[0],dir[1],dir[2]);

    std::vector<igl::Hit> hits;
    int numRays;
    bool isEmpty = myEmbree.intersectRay((origin+0.001*direction).cast<float>(),direction.cast<float>(),hits,numRays);
    if(isEmpty==true)
    {
        *d = 0.0;
        return false;
    }

    bool isHitValid = false;
    float distance = 0.0;
    for(int i=0; i<hits.size(); i++)
    {
        const igl::Hit &curHit = hits[i];
        int id = curHit.id;    // triangle id
        float u = curHit.u;    // displacement from the first vertex
        float v = curHit.v;    // displacement from the second vertex
        //cout<<"____triangle ID: "<<id<<"____"<<endl;
        //cout<<"____("<<u<<", "<<v<<")____"<<endl;
        if(u>1.0 || v>1.0) exit(9999);
        int A = myF(id,0);
        int B = myF(id,1);
        int C = myF(id,2);
        double xA = myV(A,0); double yA = myV(A,1); double zA = myV(A,2);
        double xB = myV(B,0); double yB = myV(B,1); double zB = myV(B,2);
        double xC = myV(C,0); double yC = myV(C,1); double zC = myV(C,2);

        //! -----------------------------------
        //! barycentric coordinates
        //! P=u∗A+v∗B+(1-u-v)∗C
        //! translation is added (+xA,+yA,+zA)  ??????
        //! -----------------------------------
        double xP = u*xA+v*xB+(1-u-v)*xC;//+xA;
        double yP = u*yA+v*yB+(1-u-v)*yC;//+yA;
        double zP = u*zA+v*zB+(1-u-v)*zC;//+zA;
        //fprintf(fp,"%lf\t%lf\t%lf\n",xP,yP,zP);

        //cout<<"____("<<xP<<", "<<yP<<", "<<zP<<")____"<<endl;

        float dot = (xP-P[0])*dir[0]+(yP-P[1])*dir[1]+(zP-P[2])*dir[2];
        dot /= sqrt(pow(xP-P[0],2)+pow(yP-P[1],2)+pow(zP-P[2],2))*1.0;  // 1.0 is the modulus of n
        if(dot<-1) dot = -1;
        if(dot>1) dot = 1;
        float angle = std::acos(dot);
        //cout<<"____angle: "<<angle*180.0/3.1415926534<<"____"<<endl;
        if(angle<=5.0*3.1415926534/180.0)
        {
            distance = sqrt(pow(xP-P[0],2)+pow(yP-P[1],2)+pow(zP-P[2],2));
            isHitValid = true;
            break;
        }
        else
        {
            distance = 0.0;
        }
    }
    //fclose(fp);

    if(isHitValid==true)
    {
        *d = distance;
        cout<<"____Distance: "<<distance<<"____"<<endl;
        return true;
    }
    return false;
}
