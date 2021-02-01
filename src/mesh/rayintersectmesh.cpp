//! ----------------
//! custom includes
//! ----------------
#include "rayintersectmesh.h"
#include "polygon.h"
#include <MeshVS_HArray1OfSequenceOfInteger.hxx>
#include <MeshVS_EntityType.hxx>
#include <meshelementbycoords.h>
#include "src/utils/hash_c.h"

//! ----
//! OCC
//! ----
#include <TColStd_MapIteratorOfPackedMapOfInteger.hxx>
#include <TColStd_SequenceOfInteger.hxx>

//! ----
//! C++
//! ----
#include <vector>

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
rayIntersectMesh::rayIntersectMesh(QObject* parent):QObject(parent)
{
    ;
}

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
rayIntersectMesh::rayIntersectMesh(const occHandle(MeshVS_DataSource) &aMeshDS, QObject* parent):QObject(parent)
{
    myMeshDS = aMeshDS;
}

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
rayIntersectMesh::rayIntersectMesh(const occHandle(MeshVS_DataSource) &aMeshDS, double *origin, double *direction, QObject *parent):QObject(parent)
{
    myMeshDS = aMeshDS;
    myOrigin = origin;
    myDirection = direction;
}

//! -----------------
//! function: setRay
//! details:
//! ------------------
void rayIntersectMesh::setRay(double *origin, double *direction)
{
    for(int i=0; i<3; i++)
    {
        myDirection[i] = direction[i];
        myOrigin[i] = origin[i];
    }
}

//! ------------------
//! function: setMesh
//! details:
//! ------------------
void rayIntersectMesh::setMesh(const occHandle(MeshVS_DataSource) &aMeshDS)
{
    myMeshDS = aMeshDS;
}

//! ----------------------------
//! function: retrieveTriangles
//! details:
//! ----------------------------
bool rayIntersectMesh::planeMeshIntersection(occHandle(Ng_MeshVS_DataSourceFace) &aSliceMesh)
{
    struct segment
    {
        int v1,v2;
        segment(int av1 =0, int av2=0){ v1=av1; v2=av2; }
        segment(const segment &aSegment){ v1 = aSegment.v1; v2 = aSegment.v2; }
        bool operator == (const segment &aSegment)
        {
            if(v1==aSegment.v1 && v2==aSegment.v2) return true;
            if(v1==aSegment.v2 && v2==aSegment.v1) return true;
            return false;
        }
        segment operator = (const segment &aSegment)
        {
            v1 = aSegment.v1; v2 = aSegment.v2;
            return *this;
        }
        void sort() { if(v2<v1) { int x=v2; v2=v1; v1=x; } }
        bool operator < (const segment &aSegment) const
        {
            //sort the node IDs of the segments
            int vv1, vv2; if(v2<v1) { vv1=v2; vv2=v1; }
            int vvv1, vvv2; if(aSegment.v2<aSegment.v1) { vvv1=v2; vvv2=v1; }

            //! -------------
            //! hash combine
            //! -------------
            size_t seed1=0;
            size_t seed2=0;
            hash_c<int>(seed1,vv1); hash_c<int>(seed1,vv2);
            hash_c<int>(seed2,vvv1); hash_c<int>(seed2,vvv2);
            if(seed1<seed2) return true; return false;
        }
    };

    std::vector<meshElementByCoords> vecSliceElements;
    std::map<segment,double*> alreadyVisitedSegments;
    TColStd_MapIteratorOfPackedMapOfInteger it(myMeshDS->GetAllElements());
    for(int localElementID = 1; localElementID<=myMeshDS->GetAllElements().Extent(); localElementID++, it.Next())
    {
        meshElementByCoords aSliceElement;

        int globalElementID = it.Key();
        int NbNodes;
        MeshVS_EntityType aType;
        double buf[60];
        TColStd_Array1OfReal coords(*buf,1,60);
        myMeshDS->GetGeom(globalElementID,true,coords,NbNodes,aType);
        std::vector<polygon::Point> aCloud;
        for(int n=0; n<NbNodes; n++)
        {
            int s = 3*n;
            polygon::Point aP(coords(s+1),coords(s+2),coords(s+3));
            aCloud.push_back(aP);
        }
        //! ------------------
        //! test intersection
        //! ------------------
        double a = myDirection[0];
        double b = myDirection[1];
        double c = myDirection[2];
        double d = -a*myOrigin[0]-b*myOrigin[1]-c*myOrigin[2];
        bool intersect = polygon::testPolygonPlaneIntersection(aCloud,a,b,c,d);
        if(intersect == false) continue;

        //! ---------------
        //! scan the faces
        //! ---------------
        occHandle(MeshVS_HArray1OfSequenceOfInteger) topology;
        myMeshDS->Get3DGeom(globalElementID,NbNodes,topology);

        int NbFaces = topology->Length();
        for(int f = 1; f<=NbFaces; f++)
        {
            //! ---------------------------------
            //! scan the segments of the element
            //! optimization ... to do
            //! ---------------------------------
            int NbNodesOfElement, aBuf[20];
            TColStd_Array1OfInteger nodeIDs(*aBuf,1,20);
            myMeshDS->GetNodesByElement(globalElementID,nodeIDs,NbNodesOfElement);

            TColStd_SequenceOfInteger aSeq = topology->Value(f);
            int N = aSeq.Length();  // number of segments
            for(int i = 0; i<N; i++)
            {
                int index1 = aSeq.Value(i+1)+1;
                int globalNodeID1 = nodeIDs(index1);
                int index2 = aSeq.Value((i+1)%N+1)+1;
                int globalNodeID2 = nodeIDs(index2);

                std::map<segment,double*>::iterator it = alreadyVisitedSegments.find(segment(globalNodeID1,globalNodeID2));
                if(it != alreadyVisitedSegments.end()) continue;

                //! -------------------------------------------------------------
                //! try to find an intersection: test plane segment intersection
                //! -------------------------------------------------------------
                myMeshDS->GetGeom(globalNodeID1,false,coords,NbNodes,aType);
                double xP0 = coords(1);
                double yP0 = coords(2);
                double zP0 = coords(3);

                myMeshDS->GetGeom(globalNodeID2,false,coords,NbNodes,aType);
                double xP1 = coords(1);
                double yP1 = coords(2);
                double zP1 = coords(3);

                //! -------------------------------------------
                //! parametric equation of the (P1-P0) segment
                //! -------------------------------------------
                double L = sqrt(pow(xP1-xP0,2)+pow(yP1-yP0,2)+pow(zP1-zP0,2));
                if(L<1e-6) continue;
                double kx = (xP1-xP0)/L;
                double ky = (yP1-yP0)/L;
                double kz = (zP1-zP0)/L;
                double den = a*kx+b*ky*c*kz;
                if(fabs(den)<1e-6) continue; // no intersection
                double t = -(a*xP0+b*yP0*c*zP0+d)/(a*kx+b*ky*c*kz);
                //if(t==1 && t==0) continue;  // invalid intersections

                //! ----------------------------
                //! valid point of intersection
                //! ----------------------------
                double xI = xP0+kx*t;
                double yI = yP0+ky*t;
                double zI = zP0+kz*t;
                mesh::meshPoint aP(xI,yI,zI);
                aSliceElement.pointList.push_back(aP);
                double intersection[3]{xI, yI, zI};
                alreadyVisitedSegments.insert(std::make_pair(segment(globalNodeID1,globalNodeID2),intersection));
            }
            if(aSliceElement.pointList.size()<3) continue;
            switch(aSliceElement.pointList.size())
            {
            case 3: aSliceElement.type = TRIG; break;
            case 4: aSliceElement.type = QUAD; break;
            case 5: aSliceElement.type = PENTA; break;
            case 6: aSliceElement.type = TRIG6; break;
            case 7: aSliceElement.type = EPTA; break;
            case 8: aSliceElement.type = QUAD8; break;
            }
            vecSliceElements.push_back(aSliceElement);
        }
    }
    if(vecSliceElements.size()==0) return false;
    aSliceMesh = new Ng_MeshVS_DataSourceFace(vecSliceElements,true,true);
    if(aSliceMesh.IsNull()) return false;
    return true;
}
