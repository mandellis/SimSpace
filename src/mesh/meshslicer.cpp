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
/*
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

    bool toBeAdded = false;
    TColStd_MapIteratorOfPackedMapOfInteger it(myMeshDS->GetAllElements());
    int NbElements = myMeshDS->GetAllElements().Extent();
    for(int localElementID = 1; localElementID<=NbElements; localElementID++, it.Next())
    {
        int globalElementID = it.Key();
        myMeshDS->GetNodesByElement(globalElementID,nodeIDs,NbNodes);
        toBeAdded = true;
        int k=0;
        for(int i=1; i<=NbNodes; i++)
        {
            int globalNodeID = nodeIDs(i);
            myMeshDS->GetGeom(globalNodeID,false,coords,NbNodes,aType);
            double distance = polygon::pointPlaneDistance(polygon::Point(coords(1),coords(2),coords(3)),a,b,c,d);
            //if(distance>=0) { toBeAdded = false; break; }
            if(distance>=0) { break; }
            k++;
        }
        //if(toBeAdded == true) amap.Add(globalElementID);
        if(k==NbNodes) amap.Add(globalElementID);
    }
    if(hiddenElementIDs.IsNull()) hiddenElementIDs = new TColStd_HPackedMapOfInteger;
    hiddenElementIDs->ChangeMap() = amap;
    if(amap.IsEmpty()) return false;
    return true;
}
*/

//! ---------------------------
//! function: perform
//! details:  try optimization
//! ---------------------------
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
