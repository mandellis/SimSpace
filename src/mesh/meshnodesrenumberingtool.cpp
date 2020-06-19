#include "meshnodesrenumberingtool.h"

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
meshNodesRenumberingTool::meshNodesRenumberingTool(const occHandle(MeshVS_DataSource) &sourceDS)
{
    if(sourceDS.IsNull()==false) this->setSource(sourceDS);
}

//! --------------------
//! function: setSource
//! details:
//! --------------------
void meshNodesRenumberingTool::setSource(const occHandle(MeshVS_DataSource) &sourceDS)
{
    //! ----------------------------------------
    //! prepare the information for renumbering
    //! ----------------------------------------
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
}
/*
//! ------------------
//! function: perform
//! details:
//! ------------------
template <class T>
void meshNodesRenumberingTool::perform(T* targetDS)
{
    cout<<typeid(T).name()<<endl; exit(10);
    if(strcmp(typeid(T).name(),"Ng_MeshVS_DataSourceFace")==0)
    {
        occHandle(Ng_MeshVS_DataSourceFace)::DownCast(targetDS)->renumberNodes(mapNodeCoords_globalNodeID,vecTolerantPoints,globalNodeIDs);
    }
    if(strcmp(typeid(T).name(),"Ng_MeshVS_DataSource2D")==0)
    {
        occHandle(Ng_MeshVS_DataSource2D)::DownCast(targetDS)->renumberNodes(mapNodeCoords_globalNodeID,vecTolerantPoints,globalNodeIDs);
    }
}
*/
