//! custom includes
#include "occmeshtoccxmesh.h"

//! OCC
#include <MeshVS_EntityType.hxx>
#include <MeshVS_DataSource.hxx>
#include <TColStd_MapIteratorOfPackedMapOfInteger.hxx>

QMap<int,int> OCCMeshToCCXmesh::perform(GeometryTag loc, meshDataBase *mDB)
{
    int bodyIndex = loc.parentShapeNr;
    int subTopologyIndex = loc.subTopNr;
    TopAbs_ShapeEnum type = loc.subShapeType;

    QMap<int,int> map;

    if(type==TopAbs_SOLID)
    {
        //cout<<"creating map for solids"<<endl;
        int offset = 0;
        for(int k=1; k<bodyIndex; k++)
        {
            //if(mDB->MapOfIsActive.value(k)==Property::SuppressionStatus_Active)
            if(!mDB->ArrayOfMeshDS.value(k).IsNull())
                offset = offset+mDB->ArrayOfMeshDS.value(k)->GetAllNodes().Extent();
        }
        //cout<<"____node offset calculated"<<endl;
        //cout<<"____bodyIndex: "<<bodyIndex<<endl;
        //cout<<mDB->bodyMap.size()<<endl;

        //occHandle(MeshVS_DataSource) theCurMeshDS = mDB->ArrayOfMeshDS.value(bodyIndex);
        occHandle(Ng_MeshVS_DataSource3D) theCurMeshDS = occHandle(Ng_MeshVS_DataSource3D)::DownCast(mDB->ArrayOfMeshDS.value(bodyIndex));
        if(theCurMeshDS.IsNull()) cout<<"____INVALID MESH____"<<endl;
        cout<<"____"<<theCurMeshDS->GetAllNodes().Extent()<<"____"<<endl;
        cout<<"____"<<theCurMeshDS->GetAllElements().Extent()<<"____"<<endl;
        TColStd_MapIteratorOfPackedMapOfInteger anIter;
        for(anIter.Initialize(theCurMeshDS->GetAllNodes()); anIter.More(); anIter.Next())
        {
            int nodeID = anIter.Key()+offset;
            map.insert(nodeID,anIter.Key());
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

        TColStd_MapIteratorOfPackedMapOfInteger anIter;
        for(anIter.Initialize(theCurMeshDS->GetAllNodes()); anIter.More(); anIter.Next())
        {
            int nodeID = anIter.Key()+offset;
            map.insert(nodeID,anIter.Key());
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
