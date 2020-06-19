#ifndef RENUMBERINGTOOL_H
#define RENUMBERINGTOOL_H

#include <mesh.h>
#include <ng_meshvs_datasource2d.h>
#include <ng_meshvs_datasource3d.h>
#include <ng_meshvs_datasourceface.h>
#include <hash_c.h>
#include "occhandle.h"

#include <map>
#include <vector>

#include <TColStd_MapIteratorOfPackedMapOfInteger.hxx>
#include <TColStd_PackedMapOfInteger.hxx>
#include <TColStd_IndexedMapOfInteger.hxx>
#include <MeshVS_EntityType.hxx>
#include <TColStd_Array1OfReal.hxx>

template <class S, class T>
void renumberMeshNodes(S* sourceDS, T* targetDS)
{
    //! ----------------------------------------
    //! prepare the information for renumbering
    //! ----------------------------------------
    std::map<size_t,int> mapNodeCoords_globalNodeID;
    std::vector<int> globalNodeIDs;
    std::vector<mesh::tolerantPoint> vecTolerantPoints;

    int NbNodes;
    double buf[3];
    TColStd_Array1OfReal coords(*buf,1,3);
    MeshVS_EntityType aType;
    for(TColStd_MapIteratorOfPackedMapOfInteger it(sourceDS->GetAllNodes()); it.More(); it.Next())
    {
        int globalNodeID = it.Key();
        bool isDone = sourceDS->GetGeom(globalNodeID,false,coords,NbNodes,aType);
        if(isDone==false) continue;

        //! ---------------------------------------------------------
        //! vector of the points (x,y,z) and vector of their indices
        //! ---------------------------------------------------------
        vecTolerantPoints.push_back(mesh::tolerantPoint(coords(1),coords(2),coords(3)));
        globalNodeIDs.push_back(globalNodeID);

        size_t seed = 0;
        for(int j=1; j<=3; j++) hash_c<double>(seed, coords(j));
        std::pair<size_t,int> aPair;
        aPair.first = seed;
        aPair.second = globalNodeID;
        mapNodeCoords_globalNodeID.insert(aPair);
    }

    if(strcmp(typeid(T).name(),"Ng_MeshVS_DataSourceFace")==0)
        occHandle(Ng_MeshVS_DataSourceFace)::DownCast(targetDS)->renumberNodes(mapNodeCoords_globalNodeID,vecTolerantPoints,globalNodeIDs);
    if(strcmp(typeid(T).name(),"Ng_MeshVS_DataSource2D")==0)
        occHandle(Ng_MeshVS_DataSourceFace)::DownCast(targetDS)->renumberNodes(mapNodeCoords_globalNodeID,vecTolerantPoints,globalNodeIDs);
}

#endif // RENUMBERINGTOOL_H
