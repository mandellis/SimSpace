//! ----------------
//! custom includes
//! ----------------
#include <meshslicer.h>
#include <polygon.h>
#include <meshelementbycoords.h>
#include <mesh.h>

//! ----
//! OCC
//! ----
#include <TColStd_MapIteratorOfPackedMapOfInteger.hxx>
#include <TColStd_Array1OfInteger.hxx>
#include <TColStd_HPackedMapOfInteger.hxx>
#include <MeshVS_EntityType.hxx>
#include <ng_meshvs_datasourceface.h>

//! ----
//! C++
//! ----
#include <vector>
#include <iostream>
#include <ppl.h>
using namespace std;

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
meshSlicer::meshSlicer(const occHandle(MeshVS_DataSource) &aMeshDS):myMeshDS(aMeshDS)
{
    ;
}

//! ------------------
//! function: perform
//! details:
//! ------------------
bool meshSlicer::perform(double a, double b, double c, double d, occHandle(TColStd_HPackedMapOfInteger) &hiddenElementIDs)
{
    //cout<<"meshSlicer::perform()->____function called____"<<endl;
    if(myMeshDS.IsNull()) return false;

    TColStd_PackedMapOfInteger amap;
    int NbNodes, bufi[20];
    TColStd_Array1OfInteger nodeIDs(*bufi,1,20);
    MeshVS_EntityType aType;
    double buf[3];
    TColStd_Array1OfReal coords(*buf,1,3);

    TColStd_MapIteratorOfPackedMapOfInteger it(myMeshDS->GetAllElements());
    int NbElements = myMeshDS->GetAllElements().Extent();
    for(int localElementID = 1; localElementID<=NbElements; localElementID++, it.Next())
    {
        int globalElementID = it.Key();
        myMeshDS->GetNodesByElement(globalElementID,nodeIDs,NbNodes);
        int k=0;
        for(int i=1; i<=NbNodes; i++)
        {
            int globalNodeID = nodeIDs(i);
            myMeshDS->GetGeom(globalNodeID,false,coords,NbNodes,aType);
            double distance = (a*coords(1)+b*coords(2)+c*coords(3)+d)/sqrt(a*a+b*b+c*c);
            if(distance>=0) break;
            k++;
        }
        if(k==NbNodes) amap.Add(globalElementID);
    }
    if(hiddenElementIDs.IsNull()) hiddenElementIDs = new TColStd_HPackedMapOfInteger;
    hiddenElementIDs->ChangeMap() = amap;
    if(amap.IsEmpty()) return false;
    return true;
}

//! --------------------------------
//! function: planeMeshIntersection
//! details:
//! --------------------------------
bool meshSlicer::planeMeshIntersection(occHandle(Ng_MeshVS_DataSourceFace) &aSliceMesh, double a, double b, double c, double d)
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
        cout<<"meshSlicer::planeMeshIntersection()->____scanning element: "<<localElementID<<"____"<<endl;
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
            cout<<"______scanning face nr: "<<f<<"____"<<endl;

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
                double t = -(a*xP0+b*yP0*c*zP0+d)/den;
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
