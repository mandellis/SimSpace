//! ----------------
//! custom includes
//! ----------------
#include "occmeshtoccxmesh.h"

//! ----
//! OCC
//! ----
#include <MeshVS_EntityType.hxx>
#include <MeshVS_DataSource.hxx>
#include <TColStd_MapIteratorOfPackedMapOfInteger.hxx>

//! --------------------------
//! function: performCCXtoOCC
//! details:
//! --------------------------
std::map<int,int> OCCMeshToCCXmesh::performCCXtoOCC(GeometryTag loc, meshDataBase *mDB)
{
    int bodyIndex = loc.parentShapeNr;
    int subTopologyIndex = loc.subTopNr;
    TopAbs_ShapeEnum type = loc.subShapeType;

    std::map<int,int> map;

    if(type==TopAbs_SOLID)
    {
        int offset = 0;
        for(int k=1; k<bodyIndex; k++)
        {
            //if(mDB->MapOfIsActive.value(k)==Property::SuppressionStatus_Active)
            if(!mDB->ArrayOfMeshDS.value(k).IsNull())
                offset = offset+mDB->ArrayOfMeshDS.value(k)->GetAllNodes().Extent();
        }

        occHandle(MeshVS_DataSource) theCurMeshDS = mDB->ArrayOfMeshDS.value(bodyIndex);
        //occHandle(Ng_MeshVS_DataSource3D) theCurMeshDS = occHandle(Ng_MeshVS_DataSource3D)::DownCast(mDB->ArrayOfMeshDS.value(bodyIndex));
        if(theCurMeshDS.IsNull()) return map;
        for(TColStd_MapIteratorOfPackedMapOfInteger anIter(theCurMeshDS->GetAllNodes()); anIter.More(); anIter.Next())
        {
            int nodeID = anIter.Key()+offset;
            map.insert(std::make_pair(nodeID,anIter.Key()));
        }
    }
    if(type==TopAbs_FACE)
    {
        int offset = 0;
        for(int k=1; k<bodyIndex; k++)
        {
            if(!mDB->ArrayOfMeshDS.value(k).IsNull())
            offset = offset+mDB->ArrayOfMeshDS.value(k)->GetAllNodes().Extent();
        }

        occHandle(MeshVS_DataSource) theCurMeshDS = mDB->ArrayOfMeshDSOnFaces.getValue(bodyIndex,subTopologyIndex);
        if(theCurMeshDS.IsNull()) return map;
        for(TColStd_MapIteratorOfPackedMapOfInteger anIter(theCurMeshDS->GetAllNodes()); anIter.More(); anIter.Next())
        {
            int nodeID = anIter.Key()+offset;
            map.insert(std::make_pair(nodeID,anIter.Key()));
        }
    }
    if(type==TopAbs_EDGE)
    {
        //! to do
    }
    if(type==TopAbs_VERTEX)
    {
        //! to do
    }
    return map;
}

//! --------------------------
//! function: performOCCtoCCX
//! details:
//! --------------------------
std::map<int,int> OCCMeshToCCXmesh::performOCCtoCCX(GeometryTag loc, meshDataBase *mDB)
{
    int bodyIndex = loc.parentShapeNr;
    int subTopologyIndex = loc.subTopNr;
    TopAbs_ShapeEnum type = loc.subShapeType;

    std::map<int,int> map;

    if(type==TopAbs_SOLID)
    {
        int offset = 0;
        for(int k=1; k<bodyIndex; k++)
        {
            //if(mDB->MapOfIsActive.value(k)==Property::SuppressionStatus_Active)
            if(!mDB->ArrayOfMeshDS.value(k).IsNull())
                offset = offset+mDB->ArrayOfMeshDS.value(k)->GetAllNodes().Extent();
        }

        occHandle(MeshVS_DataSource) theCurMeshDS = mDB->ArrayOfMeshDS.value(bodyIndex);
        //occHandle(Ng_MeshVS_DataSource3D) theCurMeshDS = occHandle(Ng_MeshVS_DataSource3D)::DownCast(mDB->ArrayOfMeshDS.value(bodyIndex));
        if(theCurMeshDS.IsNull()) return map;
        for(TColStd_MapIteratorOfPackedMapOfInteger anIter(theCurMeshDS->GetAllNodes()); anIter.More(); anIter.Next())
        {
            int nodeID = anIter.Key()+offset;
            map.insert(std::make_pair(anIter.Key(),nodeID));
        }
    }
    if(type==TopAbs_FACE)
    {
        int offset = 0;
        for(int k=1; k<bodyIndex; k++)
        {
            offset = offset+mDB->ArrayOfMeshDS.value(k)->GetAllNodes().Extent();
        }

        occHandle(MeshVS_DataSource) theCurMeshDS = mDB->ArrayOfMeshDSOnFaces.getValue(bodyIndex,subTopologyIndex);
        if(theCurMeshDS.IsNull()) return map;
        for(TColStd_MapIteratorOfPackedMapOfInteger anIter(theCurMeshDS->GetAllNodes()); anIter.More(); anIter.Next())
        {
            int nodeID = anIter.Key()+offset;
            map.insert(std::make_pair(anIter.Key(),nodeID));
        }
    }
    if(type==TopAbs_EDGE)
    {
        //! to do
    }
    if(type==TopAbs_VERTEX)
    {
        //! to do
    }
    return map;
}
