//! ----------------
//! custom includes
//! ----------------
#include "bolttool.h"
#include <meshelement2d.h>
#include <meshelementbycoords.h>
#include <polygon.h>
#include "src/utils/hash_c.h"

//! ----
//! C++
//! ----
#include <algorithm>
#include <vector>
#include <map>
#include <set>
#include <utility>

//! ----
//! OCC
//! ----
#include <MeshVS_HArray1OfSequenceOfInteger.hxx>
#include <TColStd_MapIteratorOfPackedMapOfInteger.hxx>
#include <TColStd_Array1OfInteger.hxx>

//! ------------
//! constructor
//! ------------
boltTool::boltTool(const occHandle(Ng_MeshVS_DataSource3D) &aVolumeMeshDS):myVolumeMesh(aVolumeMeshDS)
{
    ;
}

//! -----------------------------
//! function: sliceMeshWithPlane
//! details:
//! -----------------------------
bool boltTool::sliceMeshWithPlane(double a, double b, double c, double d,
                                  occHandle(MeshVS_DataSource) &slicedMeshDS,
                                  std::vector<std::pair<int,int>> &vecCCXFaceDefs)
{
    cout<<"boltTool::sliceMeshWithPlane()->____function called____"<<endl;

    std::vector<meshElementByCoords> volumeElementsList;

    //! -----------------------------------------------
    //! search for the intersected volume elements
    //! prepare also the map for renumbering the nodes
    //! -----------------------------------------------
    std::map<size_t,int> mapNodeCoordsToNodeGlobalID;
    std::vector<mesh::tolerantPoint> vecTolerantPoints;
    std::vector<int> vecGlobalNodeIDs;

    for(TColStd_MapIteratorOfPackedMapOfInteger it = myVolumeMesh->GetAllElements(); it.More(); it.Next())
    {
        meshElementByCoords aVolumeMeshElement;
        std::vector<polygon::Point> aPointCloud;

        int globalElementID = it.Key();
        int NbNodes, nbuf[20];
        TColStd_Array1OfInteger nodeIDs(*nbuf,1,20);
        myVolumeMesh->GetNodesByElement(globalElementID,nodeIDs,NbNodes);

        aVolumeMeshElement.ID = globalElementID;
        for(int i=1; i<=NbNodes; i++)
        {
            int globalNodeID = nodeIDs(i);
            int localNodeID = myVolumeMesh->myNodesMap.FindIndex(globalNodeID);
            std::vector<double> pc = myVolumeMesh->getNodeCoordinates(localNodeID);
            double x = pc[0];
            double y = pc[1];
            double z = pc[2];
            aPointCloud.push_back(polygon::Point(x,y,z));
            aVolumeMeshElement.pointList<<mesh::meshPoint(x,y,z,globalNodeID);

            //! ----------------------------
            //! information for renumbering
            //! ----------------------------
            size_t seed = 0;
            hash_c<double>(seed,x);
            hash_c<double>(seed,y);
            hash_c<double>(seed,z);
            //std::pair<size_t,int> apair;
            //apair.first = seed;
            //apair.second = globalNodeID;
            //mapNodeCoordsToNodeGlobalID.insert(apair);
            mapNodeCoordsToNodeGlobalID.insert(std::make_pair(seed,globalNodeID));
            const double tolerance = 1e-6;
            vecTolerantPoints.push_back(mesh::tolerantPoint(x,y,z,tolerance));
            vecGlobalNodeIDs.push_back(globalNodeID);
        }
        switch(NbNodes)
        {
        case 4: aVolumeMeshElement.type = TET; break;
        case 5: aVolumeMeshElement.type = PYRAM; break;
        case 6: aVolumeMeshElement.type = PRISM; break;
        case 8: aVolumeMeshElement.type = HEXA; break;
        case 10: aVolumeMeshElement.type = TET10; break;
        case 13: aVolumeMeshElement.type = PYRAM13; break;
        case 15: aVolumeMeshElement.type = PRISM15; break;
        case 20: aVolumeMeshElement.type = HEXA20; break;
        }
        bool intersect = polygon::testPolygonPlaneIntersection(aPointCloud,a,b,c,d);
        if(!intersect) continue;
        volumeElementsList.push_back(aVolumeMeshElement);
    }
    if(volumeElementsList.size()==0) return false;
    cout<<"tag00"<<volumeElementsList.size()<<endl;

    occHandle(Ng_MeshVS_DataSource3D) volumeSlicedMesh = new Ng_MeshVS_DataSource3D(volumeElementsList,false,false);

    //! ----------------------------------------------
    //! test: visualization of the sliced volume mesh
    //! ----------------------------------------------
    //slicedMeshDS = volumeSlicedMesh;

    //! ---------------------------
    //! build the CCX connectivity
    //! ---------------------------
    std::map<meshElement2D,std::vector<std::pair<int,int>>> CCXFaceConnectivity;
    volumeSlicedMesh->buildCCXFaceToElementConnectivity(CCXFaceConnectivity);
    cout<<"tag00"<<endl;

    //! -----------------------------------------------
    //! corresponding surface mesh
    //! here the CCX connectivity should also be built
    //! -----------------------------------------------
    volumeSlicedMesh->buildFaceToElementConnectivity();
    cout<<"tag01"<<endl;

    occHandle(Ng_MeshVS_DataSource2D) surfaceSlicedMesh = new Ng_MeshVS_DataSource2D(volumeSlicedMesh);

    //! ------------------------------------------------------------------
    //! test: visualization of the surface mesh of the sliced volume mesh
    //! ------------------------------------------------------------------
    //slicedMeshDS = surfaceSlicedMesh;
    cout<<"tag03"<<endl;

    //! -----------------------------------------------------------------
    //! build the "two layers" mesh: the overall surface mesh is needed,
    //! in order to discard the "lateral" surface elements
    //! -----------------------------------------------------------------
    if(myVolumeMesh->myFaceToElements.isEmpty()) myVolumeMesh->buildFaceToElementConnectivity();
    occHandle(Ng_MeshVS_DataSource2D) overallSurfaceMesh = new Ng_MeshVS_DataSource2D(myVolumeMesh);

    //! ----------------------------------------
    //! use a map because faster when accessing
    //! ----------------------------------------
    std::map<meshElement2D,int> serviceMap;
    cout<<"tag00"<<endl;
    int h=0;
    for(TColStd_MapIteratorOfPackedMapOfInteger it(overallSurfaceMesh->GetAllElements()); it.More(); it.Next())
    {
        int NbNodes, nbuf[8];
        TColStd_Array1OfInteger nodeIDs(*nbuf,1,8);
        overallSurfaceMesh->GetNodesByElement(it.Key(),nodeIDs,NbNodes);
        meshElement2D aMeshElement2D;
        for(int i=1; i<=NbNodes; i++) aMeshElement2D.nodeIDs<<nodeIDs(i);

        //! ----------------
        //! unessential ...
        //! ----------------
        aMeshElement2D.ID = it.Key();
        switch(NbNodes)
        {
        case 3: aMeshElement2D.type = TRIG; break;
        case 4: aMeshElement2D.type = QUAD; break;
        case 6: aMeshElement2D.type = TRIG6; break;
        case 8: aMeshElement2D.type = QUAD8; break;
        }
        //std::pair<meshElement2D,int> apair;
        //apair.first = aMeshElement2D;
        //apair.second = ++h;
        //serviceMap.insert(apair);
        serviceMap.insert(std::make_pair(aMeshElement2D,++h));
    }
    cout<<"tag01"<<endl;

    //! ------------------------------------------
    //! "Two layer mesh" S2
    //! S surface mesh of the whole volume
    //! S1 surface mesh of the sliced volume mesh
    //! Si = S intersection S1
    //! S2 = S1 - Si
    //! ------------------------------------------
    cout<<"tag03"<<endl;

    std::vector<meshElementByCoords> twoLayersMeshElements;
    for(TColStd_MapIteratorOfPackedMapOfInteger it(surfaceSlicedMesh->GetAllElements()); it.More(); it.Next())
    {
        meshElementByCoords aMeshElement;

        //! ------------------------------------------------
        //! check if the current surface element is present
        //! in the "serviceMap" (if Yes, discard it)
        //! ------------------------------------------------
        int NbNodes, nbuf[8];
        TColStd_Array1OfInteger nodeIDs(*nbuf,1,8);
        surfaceSlicedMesh->GetNodesByElement(it.Key(),nodeIDs,NbNodes);
        meshElement2D aMeshElement2D;
        for(int i=1; i<=NbNodes; i++)
        {
            int globalNodeID = nodeIDs(i);
            aMeshElement2D.nodeIDs<<globalNodeID;

            int localNodeID = myVolumeMesh->myNodesMap.FindIndex(globalNodeID);
            const std::vector<double> &pc = myVolumeMesh->getNodeCoordinates(localNodeID);
            aMeshElement.pointList<<mesh::meshPoint(pc[0],pc[1],pc[2],globalNodeID);
        }

        //! ------------------------------
        //! unessential for meshElement2D
        //! ------------------------------
        aMeshElement2D.ID = it.Key();
        aMeshElement.ID = it.Key();
        switch(NbNodes)
        {
        case 3:
        {
            aMeshElement2D.type = TRIG;
            aMeshElement.type = TRIG;
        }
            break;
        case 4:
        {
            aMeshElement2D.type = QUAD;
            aMeshElement.type = QUAD;
        }
            break;
        case 6:
        {
            aMeshElement2D.type = TRIG6;
            aMeshElement.type = TRIG6;
        }
            break;
        case 8:
        {
            aMeshElement2D.type = QUAD8;
            aMeshElement.type = QUAD8;
        }
            break;
        }

        //! ------------------------------------------------
        //! if the surface element lays on the surface mesh
        //! discard it
        //! ------------------------------------------------
        if(serviceMap.count(aMeshElement2D)!=0) continue;

        twoLayersMeshElements.push_back(aMeshElement);
    }
    cout<<"tag04"<<endl;

    occHandle(Ng_MeshVS_DataSourceFace) twoLayerMeshFaceDS = new Ng_MeshVS_DataSourceFace(twoLayersMeshElements,false,false);

    //! --------------------------------------------------------------
    //! segment to element connectivity of the "two layers" face mesh
    //! --------------------------------------------------------------
    twoLayerMeshFaceDS->computeFreeMeshSegments();

    const QMap<mesh::meshSegment,QList<int>> &segmentToElement = twoLayerMeshFaceDS->mySegmentToElement;

    //! -------------------------------------------
    //! test: visualization of the two layers mesh
    //! -------------------------------------------
    //slicedMeshDS = twoLayerMeshFaceDS;

    //! --------------------------------------
    //! build a "single element" surface mesh
    //! --------------------------------------
    int seed = 1;
    int globalElementID = twoLayerMeshFaceDS->myElementsMap.FindKey(seed);

    int NbNodes, nbuf[8];
    TColStd_Array1OfInteger nodeIDs(*nbuf,1,8);
    twoLayerMeshFaceDS->GetNodesByElement(globalElementID,nodeIDs,NbNodes);
    meshElementByCoords aMeshElement;
    aMeshElement.ID = globalElementID;
    for(int i=1; i<=NbNodes; i++)
    {
        int globalNodeID = nodeIDs(i);
        int localNodeID = twoLayerMeshFaceDS->myNodesMap.FindIndex(globalNodeID);
        const std::vector<double> &pc = twoLayerMeshFaceDS->getNodeCoordinates(localNodeID);
        aMeshElement.pointList<<mesh::meshPoint(pc[0],pc[1],pc[2],globalNodeID);
    }
    cout<<"tag05"<<endl;

    switch(NbNodes)
    {
    case 3: aMeshElement.type = TRIG; break;
    case 4: aMeshElement.type = QUAD; break;    
    case 6: aMeshElement.type = TRIG6; break;
    case 8: aMeshElement.type = QUAD8; break;
    }

    //! ---------------------------------------------------
    //! test: visualization of the "one element" face mesh
    //! ---------------------------------------------------
    //slicedMeshDS = oneElementMeshFaceDS;

    //! --------------------------------------------------------------
    //! this is the vector of mesh elements that will be filled,
    //! used at the end for building the "one layer" mesh data source
    //! --------------------------------------------------------------
    std::vector<meshElementByCoords> vecOneLayerMeshElements;
    vecOneLayerMeshElements.push_back(aMeshElement);

    occHandle(Ng_MeshVS_DataSourceFace) oneElementMeshFaceDS = new Ng_MeshVS_DataSourceFace(vecOneLayerMeshElements,false,false);
    oneElementMeshFaceDS->computeFreeMeshSegments();
    //cout<<"____initial number of boundary segments (should be 3 or 4): "<<oneElementMeshFaceDS->myBoundarySegments.length()<<"____"<<endl;

    //! -------------------------------------------
    //! extend to the adjacent surface elements up
    //! to the boundary of the single mesh layer
    //! -------------------------------------------
    std::vector<int> alreadyAdded;
    alreadyAdded.push_back(globalElementID);

    int added = 1;
    while(added != 0)
    {
        added = 0;
        int NbBoundarySegments = oneElementMeshFaceDS->myBoundarySegments.length();
        //cout<<"____number of boundary segments of the \"partial\" mesh: "<<NbBoundarySegments<<"____"<<endl;

        for(int i=0 ; i<NbBoundarySegments; i++)
        {
            const mesh::meshSegment &aSegment = oneElementMeshFaceDS->myBoundarySegments.at(i);

            //! ---------------------------------------
            //! list of elements attached to a segment
            //! ---------------------------------------
            QList<int> attachedElements = segmentToElement.value(aSegment);
            //if(aSegment.nodeIDs.size()==2)
            //    cout<<"____exploring segment: ("<<aSegment.nodeIDs[0]<<", "<<aSegment.nodeIDs[1]<<")____"<<endl;
            //else if(aSegment.nodeIDs.size()==3)
            //    cout<<"____exploring segment: ("<<aSegment.nodeIDs[0]<<", "<<aSegment.nodeIDs[1]<<", "<<aSegment.nodeIDs[2]<<")____"<<endl;
            //cout<<"____number of attached elements: "<<attachedElements.length()<<"____"<<endl;

            //! --------------------------------------------------------
            //! there is not an element attached to the current segment
            //! --------------------------------------------------------
            if(attachedElements.isEmpty()) continue;

            //! --------------------------------------------------------------
            //! there is at least one element attacted to the current segment
            //! --------------------------------------------------------------
            int NbAttachedElements = attachedElements.length();
            for(int j=0; j<NbAttachedElements; j++)
            {
                int curGlobalElementID = attachedElements.at(j);

                //! --------------------------------------------
                //! define a candidate element for beeing added
                //! --------------------------------------------
                int NbNodes, buf[8];
                TColStd_Array1OfInteger nodeIDs(*buf,1,8);
                twoLayerMeshFaceDS->GetNodesByElement(curGlobalElementID,nodeIDs,NbNodes);

                meshElementByCoords me;
                me.ID = curGlobalElementID;
                for(int j=1; j<=NbNodes; j++)
                {
                    int globalNodeID = nodeIDs(j);
                    int localNodeID = twoLayerMeshFaceDS->myNodesMap.FindIndex(globalNodeID);
                    const std::vector<double> &pc = twoLayerMeshFaceDS->getNodeCoordinates(localNodeID);
                    me.pointList<<mesh::meshPoint(pc[0],pc[1],pc[2],globalNodeID);
                }
                switch(NbNodes)
                {
                case 3: me.type = TRIG; break;
                case 4: me.type = QUAD; break;
                case 6: me.type = TRIG6; break;
                case 8: me.type = QUAD8; break;
                }

                //! ---------------------------------------------
                //! check if this element has already been added
                //! ---------------------------------------------
                if(std::find(alreadyAdded.begin(),alreadyAdded.end(),curGlobalElementID)==alreadyAdded.end())
                {
                    vecOneLayerMeshElements.push_back(me);
                    alreadyAdded.push_back(curGlobalElementID);
                    added++;
                }
            }
        }

        //! ---------------
        //! build the mesh
        //! ---------------
        //cout<<"____number of piled up elements: "<<vecOneLayerMeshElements.size()<<"____"<<endl;
        oneElementMeshFaceDS = new Ng_MeshVS_DataSourceFace(vecOneLayerMeshElements,false,false);
        oneElementMeshFaceDS->computeFreeMeshSegments();
    }

    //! ------------------------------------------------------------
    //! for visualization - final check: show the final "mini mesh"
    //! ------------------------------------------------------------
    slicedMeshDS = oneElementMeshFaceDS;
    cout<<"tag06"<<endl;

    //!-----------------------------------------------------------------------
    //! note: the vector of elements attached to the mesh element 2D always
    //! contains one element, since along the code we worked onto the surface
    //! mesh of a volume mesh slice (so ".at(0)" is used)
    //! ----------------------------------------------------------------------
    for(TColStd_MapIteratorOfPackedMapOfInteger it(oneElementMeshFaceDS->GetAllElements()); it.More(); it.Next())
    {
        int globalElementID = it.Key();
        meshElement2D aMeshElement2D;
        int NbNodes, nbuf[8];
        TColStd_Array1OfInteger nodeIDs(*nbuf,1,8);
        oneElementMeshFaceDS->GetNodesByElement(globalElementID,nodeIDs,NbNodes);
        for(int i=1; i<=NbNodes; i++) aMeshElement2D.nodeIDs<<nodeIDs(i);

        //! --------------------------
        //! unessential for this task
        //! --------------------------
        aMeshElement2D.ID = globalElementID;
        switch(NbNodes)
        {
        case 3: aMeshElement2D.type = TRIG; break;
        case 4: aMeshElement2D.type = QUAD; break;
        case 6: aMeshElement2D.type = TRIG6; break;
        case 8: aMeshElement2D.type = QUAD8; break;
        }

        if(CCXFaceConnectivity.at(aMeshElement2D).size()!=1)
        {
            cerr<<"____more than one volume element (i.e. 2) attached to the current face____"<<endl;
            return false;
        }
        std::map<meshElement2D,std::vector<std::pair<int,int>>>::iterator it1 = CCXFaceConnectivity.find(aMeshElement2D);
        if(it1==CCXFaceConnectivity.end())
        {
            cerr<<"____mesh element 2D not found____"<<endl;
            return false;
        }
        std::pair<int,int> *volumeMeshGlobalID_CCXFace = &CCXFaceConnectivity.at(aMeshElement2D).at(0);
        volumeMeshGlobalID_CCXFace->first = volumeSlicedMesh->myElementsMap.FindKey(volumeMeshGlobalID_CCXFace->first);
        vecCCXFaceDefs.push_back(*volumeMeshGlobalID_CCXFace);
        //cout<<(*volumeMeshGlobalID_CCXFace).first<<", "<<(*volumeMeshGlobalID_CCXFace).second<<endl;
    }
    cout<<"tag07"<<endl;

    return true;
}
