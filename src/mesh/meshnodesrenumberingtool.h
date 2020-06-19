#ifndef MESHNODESRENUMBERINGTOOL_H
#define MESHNODESRENUMBERINGTOOL_H

#include <mesh.h>
#include <ng_meshvs_datasource2d.h>
#include <ng_meshvs_datasource3d.h>
#include <ng_meshvs_datasourceface.h>
#include <hash_c.h>
#include "occhandle.h"

#include <map>
#include <vector>

#include <MeshVS_DataSource.hxx>
#include <TColStd_MapIteratorOfPackedMapOfInteger.hxx>
#include <TColStd_PackedMapOfInteger.hxx>
#include <TColStd_IndexedMapOfInteger.hxx>
#include <MeshVS_EntityType.hxx>
#include <TColStd_Array1OfReal.hxx>


class meshNodesRenumberingTool
{
public:

    meshNodesRenumberingTool(const occHandle(MeshVS_DataSource) &sourceDS=occHandle(MeshVS_DataSource)());
    void setSource(const occHandle(MeshVS_DataSource) &sourceDS);

    template <class T>
    void perform(T* targetDS)
    {
        if(strcmp(typeid(T).name(),"class Ng_MeshVS_DataSourceFace")==0)
        {
            occHandle(Ng_MeshVS_DataSourceFace)::DownCast(targetDS)->renumberNodes(mapNodeCoords_globalNodeID,vecTolerantPoints,globalNodeIDs);
        }
        if(strcmp(typeid(T).name(),"class Ng_MeshVS_DataSource2D")==0)
        {
            occHandle(Ng_MeshVS_DataSource2D)::DownCast(targetDS)->renumberNodes(mapNodeCoords_globalNodeID,vecTolerantPoints,globalNodeIDs);
        }
    }

private:

    std::map<size_t,int> mapNodeCoords_globalNodeID;
    std::vector<int> globalNodeIDs;
    std::vector<mesh::tolerantPoint> vecTolerantPoints;
};

#endif // MESHNODESRENUMBERINGTOOL_H
