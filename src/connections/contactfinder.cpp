//! ----------------
//! custom includes
//! ----------------
#include "contactfinder.h"
#include "global.h"
#include "geometrydatabase.h"
#include "mesh.h"
#include "ng_meshvs_datasource3d.h"
#include "ng_meshvs_datasourceface.h"
#include "polygon.h"
#include "qprogressindicator.h"
#include "qprogressevent.h"
#include "extendedrwstl.h"
#include "stlapiwriter.h"
#include "meshtools.h"
#include "triangletotriangleintersection.h"

//! ----
//! OCC
//! ----
#include <MeshVS_DataSource.hxx>
#include <TColStd_MapIteratorOfPackedMapOfInteger.hxx>
#include <TColStd_Array1OfInteger.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <BRepTools.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <TopAbs_ShapeEnum.hxx>

//! ---
//! Qt
//! ---
#include <QApplication>
#include <QThread>

//! --------
//! C++ STL
//! --------
#include <set>

#include <algorithm>
#include <chrono>

//! ------------------------
//! function: contactFinder
//! details:
//! ------------------------
contactFinder::contactFinder(const std::vector<std::pair<GeometryTag,GeometryTag>> &bodyPairs,
                             geometryDataBase *gDB,
                             QProgressIndicator* aProgressIndicator):
    myGDB(gDB),
    myProgressIndicator(aProgressIndicator)
{
    bool rebuildTessellation = false;
    std::map<int,occHandle(MeshVS_DataSource)> mapOfSurfaceMeshDataSources;    // di bellezza
    prepareSTLTessellation_useDisk(bodyPairs,
                                   rebuildTessellation,
                                   mapOfSurfaceMeshDataSources,
                                   m_subMeshes);

    //! ---------------------------------------
    //! rebuild the surface mesh for each body
    //! ---------------------------------------
    //std::map<int,occHandle(Ng_MeshVS_DataSourceFace)> mapSurfaceMeshDS;
    for(std::map<int,std::map<int,occHandle(Ng_MeshVS_DataSourceFace)>>::const_iterator it = m_subMeshes.cbegin(); it!=m_subMeshes.cend(); it++)
    {
        int bodyIndex = it->first;
        const std::map<int,occHandle(Ng_MeshVS_DataSourceFace)>& subMeshes = it->second;
        occHandle(Ng_MeshVS_DataSourceFace) aSurfaceMeshDS;
        std::vector<occHandle(Ng_MeshVS_DataSourceFace)> vecMeshDS;
        for(std::map<int,occHandle(Ng_MeshVS_DataSourceFace)>::const_iterator it1 = subMeshes.cbegin(); it1!=subMeshes.cend(); it1++)
        {
            //int faceNr = it->first;
            occHandle(Ng_MeshVS_DataSourceFace) aMeshDS = it1->second;
            if(aMeshDS.IsNull()) continue;
            vecMeshDS.push_back(aMeshDS);
        }
        occHandle(Ng_MeshVS_DataSourceFace) overall = new Ng_MeshVS_DataSourceFace(vecMeshDS);
        if(overall.IsNull()) continue;
        mapSurfaceMeshDS.insert(std::make_pair(bodyIndex,overall));
        //cout<<"____"<<bodyIndex<<" NE: "<<overall->GetAllElements().Extent()<<" "<<"NN: "<<overall->GetAllNodes().Extent()<<"____"<<endl;
    }

    //! -------------------
    //! init the raycsters
    //! -------------------
    initRayCaster(mapSurfaceMeshDS);

    //! -----------------------------------------------------------
    //! key => body index;
    //! value => std::map<int,int> => key: globalElementID, faceNr
    //! -----------------------------------------------------------
    //std::map<int,std::map<int,int>> MM;
    for(std::map<int,std::map<int,occHandle(Ng_MeshVS_DataSourceFace)>>::const_iterator it= m_subMeshes.cbegin(); it!=m_subMeshes.cend();it++)
    {
        int bodyIndex = it->first;
        const std::map<int,occHandle(Ng_MeshVS_DataSourceFace)> sm = it->second;
        std::map<int,int> M;
        for(std::map<int,occHandle(Ng_MeshVS_DataSourceFace)>::const_iterator it1 = sm.cbegin(); it1!=sm.cend(); it1++)
        {
            int faceNr = it1->first;
            const occHandle(Ng_MeshVS_DataSourceFace) &aMeshDS = it1->second;
            if(aMeshDS.IsNull()) continue;
            for(TColStd_MapIteratorOfPackedMapOfInteger it2(aMeshDS->GetAllElements()); it2.More(); it2.Next())
                M.insert(std::make_pair(it2.Key(),faceNr));
        }
        MM.insert(std::make_pair(bodyIndex,M));
    }
    cout<<"____tag01____"<<endl;

}

//! -------------------------------
//! function: setProgressIndicator
//! details:
//! -------------------------------
void contactFinder::setProgressIndicator(QProgressIndicator *aProgressIndicator)
{
    myProgressIndicator = aProgressIndicator;
}

//! ------------------------
//! function: initRayCaster
//! details:
//! ------------------------
#include "igtools.h"
void contactFinder::initRayCaster(const std::map<int,occHandle(Ng_MeshVS_DataSourceFace)> &mapSurfaceMeshes)
{
    cout<<"contactFinder::initRayCaster()->____function called____"<<endl;

    for(std::map<int,occHandle(Ng_MeshVS_DataSourceFace)>::const_iterator itm = mapSurfaceMeshes.cbegin(); itm!=mapSurfaceMeshes.cend(); itm++)
    {
        int bodyIndex = itm->first;
        const occHandle(Ng_MeshVS_DataSourceFace) &aMeshDS = itm->second;
        if(aMeshDS.IsNull()) continue;

        cout<<"____working on mesh of body: "<<bodyIndex<<"____"<<endl;

        //! -------------------
        //! nodes and elements
        //! -------------------
        Eigen::MatrixXd V;
        Eigen::MatrixXi F;
        //iglTools::OCCMeshToIglMesh(aMeshDS,V,F);

        //! --------------------------
        //! convert into eigen format
        //! --------------------------
        int NN = aMeshDS->GetAllNodes().Extent();
        V.resize(NN,3);

        TColStd_MapIteratorOfPackedMapOfInteger it;
        int NbNodes;
        MeshVS_EntityType type;
        double buf[3];
        TColStd_Array1OfReal coords(*buf,1,3);
        for(it.Initialize(aMeshDS->GetAllNodes());it.More();it.Next())
        {
            int globalID = it.Key();
            int localNodeID = aMeshDS->myNodesMap.FindIndex(globalID);
            aMeshDS->GetGeom(globalID,false,coords,NbNodes,type);
            int pos = localNodeID-1;
            V(pos,0) = coords(1);
            V(pos,1) = coords(2);
            V(pos,2) = coords(3);
            //cout<<"____"<<pos<<"("<<V(pos,0)<<", "<<V(pos,1)<<", "<<V(pos,2)<<")____"<<endl;
        }

        //! ----------------
        //! adding elements
        //! ----------------
        int NE = aMeshDS->GetAllElements().Extent();
        F.resize(NE,3);
        int bufn[3];
        TColStd_Array1OfInteger nodeIDs(*bufn,1,3);     //! only triangles
        for(it.Initialize(aMeshDS->GetAllElements());it.More();it.Next())
        {
            int globalElementID = it.Key();
            int localElementID = aMeshDS->myElementsMap.FindIndex(globalElementID);
            aMeshDS->GetNodesByElement(globalElementID,nodeIDs,NbNodes);
            int pos = localElementID-1;
            for(int i=1; i<=NbNodes; i++)
            {
                int curGlobalNodeID = nodeIDs.Value(i);
                int curLocalNodeID = aMeshDS->myNodesMap.FindIndex(curGlobalNodeID);
                F(pos,i-1) = curLocalNodeID-1;
                //if(i==NbNodes)cout<<F(pos,i-1)<<endl;
                //else cout<<F(pos,i-1)<<" ";
            }
        }


        //! ------------
        //! init embree
        //! ------------
        std::shared_ptr<igl::embree::EmbreeIntersector> anEmbree(new igl::embree::EmbreeIntersector);
        anEmbree->init(V.cast<float>(),F.cast<int>());

        mapOfRayCasters.insert(std::make_pair(bodyIndex,anEmbree));
    }
}

/*
//! ----------------------------------------------------------------
//! function: AngularCriterion
//! details:  input a vector in 3D space with components (xA,yA,zA)
//!           input a vector in 3D space with components (xB,yB,zB)
//! ----------------------------------------------------------------
bool contactFinder::AngularCriterion(double xA, double yA, double zA, double xB, double yB, double zB)
{
    const double pi = 3.141592654;
    double dot = xA*xB+yA*yB+zA*zB;
    if(dot>1.0) dot = 1.0;
    if(dot<-1.0) dot = -1.0;
    double angle = 180.0*acos(dot)/pi;
    if(angle<=180 && angle>=180-30) return true;
    else return false;
}
*/
//! ------------------
//! function: perform
//! details:
//! ------------------
bool contactFinder::perform(std::vector<std::pair<std::vector<GeometryTag>,std::vector<GeometryTag>>> &allContactsPairs,
                            double tolerance,
                            int grouping)
{
    if(myVectorOfBodyPairs.empty()) return false;       //! bodies not set

    //! ------------------
    //! define the result
    //! ------------------
    std::vector<std::pair<std::vector<GeometryTag>,std::vector<GeometryTag>>> allContactPairs;

    //! -----------------
    //! has been stopped
    //! -----------------
    bool stopped = false;

    //! ----------------
    //! disable "Abort"
    //! ----------------
    myProgressIndicator->disableStop();

    //! ------------------------
    //! set the status: running
    //! ------------------------
    Global::status().code = 1;

    //! ---------------------
    //! create an init event
    //! ---------------------
    if(myProgressIndicator!=nullptr)
    {
                int NbPairs = int(myVectorOfBodyPairs.size());
        QString task = "Creating automatic connections";
        QProgressEvent *aProgressEvent = new QProgressEvent(QProgressEvent_Init,0,NbPairs-1,0,"Searching for contact pairs",
                                                            QProgressEvent_Init,0,100,0,task);
        QApplication::postEvent(myProgressIndicator,aProgressEvent);
        QApplication::processEvents();
        QThread::msleep(500);
    }

    //! --------------------------------------
    //! iterate over the vector of body pairs
    //! --------------------------------------
    int done = 0;
    double TotalTimeForSearch = 0;
    for(std::vector<std::pair<GeometryTag,GeometryTag>>::const_iterator it = myVectorOfBodyPairs.cbegin(); it!=myVectorOfBodyPairs.cend(); it++)
    {
        //! --------------
        //! enable "Stop"
        //! --------------
       if(myProgressIndicator!=nullptr) myProgressIndicator->enableStop();

        //! -----------------
        //! check the status
        //! -----------------
        if(Global::status().code ==0)
        {
            Global::status().code = 1;
            if(myProgressIndicator!=nullptr)
            {
                QProgressEvent* aProgressEvent = new QProgressEvent(QProgressEvent_Reset,0,100,0,"Process interrupted by the user",
                                                                    QProgressEvent_Reset,0,100,0,"");
                QApplication::postEvent(myProgressIndicator,aProgressEvent);
                QApplication::processEvents();
            }
            stopped = true;
            return stopped;
        }


        //! ---------------------------------------------
        //! define the result: a vector of mesh pairs
        //! indexed as the input vector of geometry tags
        //! ---------------------------------------------
        //std::vector<std::pair<std::vector<GeometryTag>,std::vector<GeometryTag>>> allContactPairs;
        for(std::vector<std::pair<GeometryTag,GeometryTag>>::iterator it =  myVectorOfBodyPairs.begin(); it!= myVectorOfBodyPairs.end(); it++)
        {
            const std::pair<GeometryTag,GeometryTag> aBodyPair = *it;
            int sourceIndex = aBodyPair.first.parentShapeNr;
            int targetIndex = aBodyPair.second.parentShapeNr;
            cout<<"____WORKING ON BODY PAIR: ("<<sourceIndex<<", "<<targetIndex<<")____"<<endl;

            std::shared_ptr<igl::embree::EmbreeIntersector> rayCaster = mapOfRayCasters.at(targetIndex);                    //! target mesh
            const std::map<int,occHandle(Ng_MeshVS_DataSourceFace)>& faceDataSources = m_subMeshes.at(sourceIndex);         //! source mesh

            int start = getTimePoint();

            //! iterazione sulle facce del corpo source
            for(std::map<int,occHandle(Ng_MeshVS_DataSourceFace)>::const_iterator it1=faceDataSources.cbegin(); it1!=faceDataSources.cend(); it1++)
            {
                std::vector<GeometryTag> tagsSource,tagsTarget;

                int faceNr = it1->first;
                const occHandle(Ng_MeshVS_DataSourceFace) &aMeshDS = it1->second;
                if(aMeshDS.IsNull()) continue;
                bool foundIntersectedFace = false;
                int faceNr_Intercepted = -1;

                //! analisi di una singola faccia del corpo source:
                //! termina appena trovata un elemento della faccia target
                //! (e quindi la faccia target) intersecato
                for(TColStd_MapIteratorOfPackedMapOfInteger it2(aMeshDS->GetAllElements()); it2.More(); it2.Next())
                {
                    int NbNodes, globalElementID = it2.Key();
                    MeshVS_EntityType aType;
                    double buf[60];
                    TColStd_Array1OfReal coords(*buf,1,60);
                    if(aMeshDS->GetGeom(globalElementID,true,coords,NbNodes,aType)) continue;
                    double nx,ny,nz;
                    if(!aMeshDS->GetNormal(globalElementID,20,nx,ny,nz)) continue;
                    std::vector<polygon::Point> aSurfacePolygon;
                    for(int n=0; n<NbNodes; n++)
                    {
                        int s = 3*n;
                        aSurfacePolygon.push_back(polygon::Point(coords(s),coords(s+1),coords(s+2)));
                    }
                    polygon::Point C;
                    double area = 0;
                    polygon::getPolygonCenter(aSurfacePolygon,C,area);

                    //! -----------
                    //! cast a ray
                    //! -----------
                    Eigen::RowVector3f origin(C.x,C.y,C.z);
                    Eigen::RowVector3f direction0(nx,ny,nz);    // this normal is directed towards the exterior
                    igl::Hit hit;
                    bool isDone = false;
                    try { isDone = rayCaster->intersectRay(origin.cast<float>(),direction0.cast<float>(),hit,0,1000); }
                    catch(...) { continue; }
                    if(fabs(hit.t)>tolerance) continue;
                    int interceptedElement = hit.id+1;
                    int globalElementID_intercepted = mapSurfaceMeshDS.at(targetIndex)->myElementsMap.FindKey(interceptedElement);
                    const std::map<int,int> &M = MM.at(targetIndex);
                    std::map<int,int>::const_iterator it3 = M.find(globalElementID_intercepted);
                    if(it3==M.cend()) continue;
                    faceNr_Intercepted = it3->second;
                    foundIntersectedFace = true;
                    break;
                }

                if(foundIntersectedFace)
                {
                    GeometryTag aSourceTag;
                    aSourceTag.isParent = false;
                    aSourceTag.parentShapeNr = sourceIndex;
                    aSourceTag.subTopNr = faceNr;
                    aSourceTag.subShapeType = TopAbs_FACE;
                    tagsSource.push_back(aSourceTag);

                    GeometryTag aTargetTag;
                    aTargetTag.isParent = false;
                    aTargetTag.parentShapeNr = targetIndex;
                    aTargetTag.subTopNr = faceNr_Intercepted;
                    aTargetTag.subShapeType = TopAbs_FACE;
                    tagsTarget.push_back(aTargetTag);

                    std::pair<std::vector<GeometryTag>,std::vector<GeometryTag>> apair_;
                    apair_.first = tagsSource;
                    apair_.second = tagsTarget;
                    allContactPairs.push_back(apair_);
                }
            }

        }

        if(stopped)
        {
            cout<<"contactFinder::perform()->____process interrupted by the used during mesh proximity analysis____"<<endl;
            Global::status().code = 1;

            //! reset the progress indicator
            QProgressEvent* aProgressEvent = new QProgressEvent(QProgressEvent_Reset,0,100,0,"Process interrupted by the user",
                                                QProgressEvent_Reset,0,100,0,"");
            QApplication::postEvent(myProgressIndicator,aProgressEvent);
            QApplication::processEvents();
            return stopped;
        }

        //! post an update event
        if(myProgressIndicator!=nullptr)
        {
            QProgressEvent* aProgressEvent = new QProgressEvent(QProgressEvent_Update,0,9999,done,"Contact pairs found",
                                                                QProgressEvent_None,0,9999,0,"Creating automatic connections");
            QApplication::postEvent(myProgressIndicator,aProgressEvent);
            QApplication::processEvents();
            done++;
        }
    }

    cout<<"@------------------------------------@"<<endl;
    cout<<"@- total time for search [s] "<<TotalTimeForSearch<<endl;
    cout<<"@------------------------------------@"<<endl;

    if(myProgressIndicator!=nullptr)
    {
        QProgressEvent *aProgressEvent = new QProgressEvent(QProgressEvent_Reset,0,100,0,"Contact pairs found",
                                                            QProgressEvent_Reset,0,100,0,"");
        QApplication::postEvent(myProgressIndicator,aProgressEvent);
        QApplication::processEvents();
        QThread::msleep(500);
    }

    std::vector<std::pair<std::vector<GeometryTag>,std::vector<GeometryTag>>> reordered;
    this->groupBy(allContactsPairs,reordered,grouping);
    allContactsPairs.clear();       //! unessential: left here for showing the intention
    allContactsPairs = reordered;

    //! return true: the process has not been stopped
    return true;
}

//! -----------------------
//! function: getTimePoint
//! details:
//! -----------------------
int contactFinder::getTimePoint()
{
    std::chrono::time_point<std::chrono::system_clock> now = std::chrono::system_clock::now();
    auto duration = now.time_since_epoch();
    auto millis = std::chrono::duration_cast<std::chrono::milliseconds>(duration).count();
    return millis;
}

//! --------------------------------------------------------------------------------------------
//! function: prepareSTLTessellation_useDisk
//! details:  given a vector of body pairs (actually of GeometryTag) returns:
//!           1. a map whose index is a body index, and whose value is the STL mesh data source
//!           2. a map whose index is a body index, and whose value is a map, whose index is
//!              an STL triangle hash code and whose value the corresponding face number
//! --------------------------------------------------------------------------------------------
#include <experimental/filesystem>
namespace fs = experimental::filesystem;
int contactFinder::prepareSTLTessellation_useDisk(const std::vector<std::pair<GeometryTag,GeometryTag>> &vectorOfBodyPairs,
                                                  bool rebuildTessellation,
                                                  std::map<int,occHandle(MeshVS_DataSource)> &mapOfSurfaceMeshDataSources,
                                                  std::map<int,std::map<int,occHandle(Ng_MeshVS_DataSourceFace)>> &mapOfMapOfFaceMeshDataSources)
{
    cout<<"contactFinder::prepareSTLTessellation_useDisk()->____function called____"<<endl;

    std::set<int> parentShapeNumbers;
    for(std::vector<std::pair<GeometryTag,GeometryTag>>::const_iterator it = vectorOfBodyPairs.cbegin();
        it != vectorOfBodyPairs.cend(); it++)
    {
        parentShapeNumbers.insert(it->first.parentShapeNr);
        parentShapeNumbers.insert(it->second.parentShapeNr);
    }

    //! -------------------------
    //! rebuild the tessellation
    //! -------------------------
    rebuildTessellation = true;
    if(rebuildTessellation)
    {
        for(std::set<int>::iterator it = parentShapeNumbers.begin(); it!= parentShapeNumbers.end(); it++)
        {
            int bodyIndex = *it;
            //const body &curBody = myGDB->MapOfBodyTopologyMap
            //        myGDB->mapOfBodies.at(bodyIndex);
            const TopoDS_Shape &curShape = myGDB->bodyMap.value(bodyIndex);
            if(curShape.IsNull()) continue;
            double ad = 0.3; double ld = 0.1;
            BRepTools::Clean(curShape);
            BRepMesh_IncrementalMesh aMesher(curShape,ld,true,ad,true,false);
        }
    }
    cout<<"contactFinder::prepareSTLTessellation_useDisk()->____tag00____"<<endl;

    std::string tempDir = getPathOfExecutable()+"\\supportFilesDir";
    if(fs::exists(fs::path(tempDir)))
    {
        fs::remove_all(fs::path(tempDir));
        fs::create_directories(fs::path(tempDir));
    }
    else fs::create_directories(fs::path(tempDir));
    cout<<"contactFinder::prepareSTLTessellation_useDisk()->____tag01____"<<endl;

    //! -------------------------------------------------
    //! convert the tessellations into mesh data sources
    //! -------------------------------------------------
    int count = 0;
    for(std::set<int>::iterator it = parentShapeNumbers.begin(); it!= parentShapeNumbers.end(); it++)
    {
        int bodyIndex = *it;
        //const body &curBody = myGDB->mapOfBodies.at(bodyIndex);
        //const TopoDS_Shape &curShape = curBody.shape();
        const TopoDS_Shape &curShape = myGDB->bodyMap.value(bodyIndex);

        if(curShape.IsNull()) continue;

        std::string fileName = tempDir+"\\STLMesh_"+std::to_string(++count)+".stl";
        QString filenameQt =QString::fromStdString(fileName);
        //! ----------------------------
        //! write the extended stl file
        //! ----------------------------
        STLAPIWriter::WriteExtended(curShape,fileName.c_str());

        //! ---------------------------
        //! read the extended stl file
        //! ---------------------------
        //int NbFaces = myGDB->mapOfBodies.at(bodyIndex).faceMap().Extent();
        int NbFaces = myGDB->MapOfBodyTopologyMap.value(bodyIndex).faceMap.Extent();
        occHandle(Ng_MeshVS_DataSource2D) anSTLMesh;

        std::vector<occHandle(Ng_MeshVS_DataSourceFace)> faceSubMeshes;
        for(int i=0; i<=NbFaces; i++) faceSubMeshes.push_back(occHandle(Ng_MeshVS_DataSourceFace)());
        MeshTools::toSTLMesh(curShape,filenameQt,anSTLMesh,faceSubMeshes);

        //! ----------------------------------------------
        //! fill the map of the surface mesh data sources
        //! ----------------------------------------------
        mapOfSurfaceMeshDataSources.insert(std::make_pair(bodyIndex,anSTLMesh));

        //! -----------------------------------------------
        //! fill the map of the face mesh stl data sources
        //! -----------------------------------------------
        std::map<int,occHandle(Ng_MeshVS_DataSourceFace)> mapOfFaceMeshDataSources;
        for(int faceNr = 0; faceNr<=NbFaces; faceNr++)
        {
            occHandle(Ng_MeshVS_DataSourceFace) aFaceMeshDS = faceSubMeshes[faceNr];
            if(aFaceMeshDS.IsNull()) continue;
            mapOfFaceMeshDataSources.insert(std::make_pair(faceNr,aFaceMeshDS));
        }
        mapOfMapOfFaceMeshDataSources.insert(std::make_pair(bodyIndex,mapOfFaceMeshDataSources));
    }
    cout<<"contactFinder::prepareSTLTessellation_useDisk()->____tag02____"<<endl;

    //! remov the support files directory
    //fs::remove_all(fs::path(tempDir));
    cout<<"contactFinder::prepareSTLTessellation_useDisk()->____tag03____"<<endl;

    return int(mapOfSurfaceMeshDataSources.size());
    cout<<"contactFinder::prepareSTLTessellation_useDisk()->____exiting____"<<endl;

}


//! -------------------
//! function: groupBy
//! details:
//! -------------------
void contactFinder::groupBy(const std::vector<std::pair<std::vector<GeometryTag>,std::vector<GeometryTag>>> &allContactsPairs,
                             std::vector<std::pair<std::vector<GeometryTag>,std::vector<GeometryTag>>> &reordered,
                             int grouping)
{
    switch(grouping)
    {
    //! ---------------
    //! by bodies pair
    //! ---------------
    case 0:
    {
        cout<<"____FILTERING BY BODY PAIR____"<<endl;
        for(std::vector<std::pair<GeometryTag,GeometryTag>>::iterator it = myVectorOfBodyPairs.begin(); it!=myVectorOfBodyPairs.end(); ++it)
        {
            std::pair<std::vector<GeometryTag>,std::vector<GeometryTag>> reorderedPair;
            std::vector<GeometryTag> tagsMasterReordered;
            std::vector<GeometryTag> tagsSlaveReordered;

            std::pair<GeometryTag,GeometryTag> &curPair = *it;
            int masterBodyIndex = curPair.first.parentShapeNr;
            int slaveBodyIndex = curPair.second.parentShapeNr;

            int Nb = int(allContactsPairs.size());
            for(int i=0; i<Nb; i++)
            {
                std::pair<std::vector<GeometryTag>,std::vector<GeometryTag>> aContactPair = allContactsPairs[i];
                std::vector<GeometryTag> tagsMaster = aContactPair.first;
                std::vector<GeometryTag> tagsSlave = aContactPair.second;
                GeometryTag aMtag = tagsMaster.at(0);
                GeometryTag aStag = tagsSlave.at(0);
                if(aMtag.parentShapeNr == masterBodyIndex && aStag.parentShapeNr == slaveBodyIndex)
                {
                    tagsMasterReordered.push_back(aMtag);
                    tagsSlaveReordered.push_back(aStag);
                }
            }
            if(tagsMasterReordered.size()>0 && tagsSlaveReordered.size()>0)
            {
                reorderedPair.first = tagsMasterReordered;
                reorderedPair.second = tagsSlaveReordered;
                reordered.push_back(reorderedPair);
            }
        }
    }
        break;

    //! ----------------
    //! by master faces
    //! ----------------
    case 1:
    {
        std::map<GeometryTag,std::vector<GeometryTag>> supportMap;
        int Nb = int(allContactsPairs.size());
        for(int i=0; i<Nb; i++)
        {
            std::pair<std::vector<GeometryTag>,std::vector<GeometryTag>> aContact = allContactsPairs[i];
            GeometryTag masterTag = aContact.first[0];
            std::map<GeometryTag,std::vector<GeometryTag>>::iterator it = supportMap.find(masterTag);
            if(it==supportMap.end())
            {
                std::pair<GeometryTag,std::vector<GeometryTag>> element;
                element.first = masterTag;
                element.second = aContact.second;
                supportMap.insert(element);
            }
            else
            {
                //(*it).second.append(aContact.second);
                (*it).second.insert((*it).second.end(),aContact.second.begin(),aContact.second.end());
            }
        }

        for(std::map<GeometryTag,std::vector<GeometryTag>>::iterator it = supportMap.begin(); it!=supportMap.end(); it++)
        {
            std::pair<GeometryTag,std::vector<GeometryTag>> aPair = *it;
            std::vector<GeometryTag> masterTags;
            masterTags.push_back(aPair.first);
            std::vector<GeometryTag> slaveTags;
            slaveTags = aPair.second;

            std::pair<std::vector<GeometryTag>,std::vector<GeometryTag>> aPair_;
            aPair_.first = masterTags;
            aPair_.second = slaveTags;
            reordered.push_back(aPair_);
        }

        //! --------
        //! testing
        //! --------
        cout<<"@____AGGREGATION____@"<<endl;
        for(std::map<GeometryTag,std::vector<GeometryTag>>::iterator it = supportMap.begin(); it!=supportMap.end(); ++it)
        {
            GeometryTag mt = (*it).first;
            std::vector<GeometryTag> sts = (*it).second;
            cout<<"____master: ("<<mt.parentShapeNr<<", "<<mt.subTopNr<<")____"<<endl;
            for(int i=0; i<sts.size(); i++)
            {
                cout<<"_______slave: ("<<sts[i].parentShapeNr<<", "<<sts[i].subTopNr<<")_______"<<endl;
            }
        }
        cout<<"@____END AGGREGATION____@"<<endl;
    }
        break;

    //! ---------------
    //! by slave faces
    //! ---------------
    case 2:
    {
        //! -----------------------------------------
        //! key => slave face; value => master faces
        //! -----------------------------------------
        std::map<GeometryTag,std::vector<GeometryTag>> supportMap;
        int Nb = int(allContactsPairs.size());
        for(int i=0; i<Nb; i++)
        {
            std::pair<std::vector<GeometryTag>,std::vector<GeometryTag>> aContact = allContactsPairs[i];
            GeometryTag slaveTag = aContact.second[0];
            std::map<GeometryTag,std::vector<GeometryTag>>::iterator it = supportMap.find(slaveTag);
            if(it==supportMap.end())
            {
                std::pair<GeometryTag,std::vector<GeometryTag>> element;
                element.first = slaveTag;
                element.second = aContact.first;
                supportMap.insert(element);
            }
            else
            {
                //(*it).second.append(aContact.first);
                it->second.insert(it->second.end(),aContact.first.begin(),aContact.first.end());
            }
        }

        for(std::map<GeometryTag,std::vector<GeometryTag>>::iterator it = supportMap.begin(); it!=supportMap.end(); it++)
        {
            std::pair<GeometryTag,std::vector<GeometryTag>> aPair = *it;
            std::vector<GeometryTag> slaveTags;
            slaveTags.push_back(aPair.first);
            std::vector<GeometryTag> masterTags;
            masterTags = aPair.second;

            std::pair<std::vector<GeometryTag>,std::vector<GeometryTag>> aPair_;
            aPair_.first = masterTags;
            aPair_.second = slaveTags;
            reordered.push_back(aPair_);
        }
    }
        break;

    //! ---------------
    //! by master body
    //! ---------------
    case 3:
    {
        cout<<"____FILTERING BY MASTER BODY____"<<endl;
        //! ----------------------------------------
        //! key => index of master body
        //! value => index of the array of contacts
        //! ----------------------------------------
        std::map<int,std::vector<int>> supportMap;
        int Nb = int(allContactsPairs.size());

        for(int i=0; i<Nb; i++)
        {
            std::pair<std::vector<GeometryTag>,std::vector<GeometryTag>> aContact = allContactsPairs[i];
            //! ----------------------------------------------------------------
            //! by definition the "masterTags" vector contains only one element
            //! ----------------------------------------------------------------
            std::vector<GeometryTag> masterTags = aContact.first;

            int masterBodyNumber = masterTags[0].parentShapeNr;
            std::map<int,std::vector<int>>::iterator it = supportMap.find(masterBodyNumber);
            if(it == supportMap.end())
            {
                std::pair<int,std::vector<int>> element;
                element.first = masterBodyNumber;
                std::vector<int> indexList;
                indexList.push_back(i);
                element.second = indexList;
                supportMap.insert(element);
            }
            else { (*it).second.push_back(i); }
        }

        for(std::map<int,std::vector<int>>::iterator it = supportMap.begin(); it!=supportMap.end(); it++)
        {

            std::vector<int> listOfIndexOfContactPairs = (*it).second;
            std::vector<GeometryTag> masterTags, slaveTags;
            for(int k=0; k<listOfIndexOfContactPairs.size(); k++)
            {
                int indexOfContactPair = listOfIndexOfContactPairs[k];
                std::pair<std::vector<GeometryTag>,std::vector<GeometryTag>> apair = allContactsPairs[indexOfContactPair];
                masterTags.insert(masterTags.end(),apair.first.begin(),apair.first.end());
                slaveTags.insert(slaveTags.end(),apair.second.begin(),apair.second.end());
            }
            if(masterTags.size()==0) { cerr<<"____EMPTY MASTER TAGS____"<<endl; exit(1); }
            if(slaveTags.size()==0) { cerr<<"____EMPTY SLAVE TAGS____"<<endl; exit(2); }

            //! --------------------------
            //! define a new contact pair
            //! --------------------------
            std::pair<std::vector<GeometryTag>,std::vector<GeometryTag>> aContact;
            aContact.first = masterTags;
            aContact.second = slaveTags;
            reordered.push_back(aContact);
        }
    }
        break;

    //! ----------------
    //! By slave bodies
    //! ----------------
    case 4:
    {
        cout<<"____FILTERING BY SLAVE BODY____"<<endl;
        //! ----------------------------------------
        //! key => index of slave body
        //! value => index of the array of contacts
        //! ----------------------------------------
        std::map<int,std::vector<int>> supportMap;
        int Nb = int(allContactsPairs.size());
        for(int i=0; i<Nb; i++)
        {
            std::pair<std::vector<GeometryTag>,std::vector<GeometryTag>> aContact = allContactsPairs[i];

            //! ---------------------------------------------------------------
            //! by definition the "slaveTags" vector contains only one element
            //! ---------------------------------------------------------------
            std::vector<GeometryTag> slaveTags = aContact.second;

            int slaveBodyNumber = slaveTags[0].parentShapeNr;
            std::map<int,std::vector<int>>::iterator it = supportMap.find(slaveBodyNumber);
            if(it == supportMap.end())
            {
                std::pair<int,std::vector<int>> element;
                element.first = slaveBodyNumber;
                std::vector<int> indexList;
                indexList.push_back(i);
                element.second = indexList;
                supportMap.insert(element);
            }
            else { (*it).second.push_back(i); }
        }

        for(std::map<int,std::vector<int>>::iterator it = supportMap.begin(); it!=supportMap.end(); it++)
        {
            std::vector<int> listOfIndexOfContactPairs = (*it).second;

            std::vector<GeometryTag> masterTags, slaveTags;
            for(int k=0; k<listOfIndexOfContactPairs.size(); k++)
            {
                int indexOfContactPair = listOfIndexOfContactPairs[k];
                std::pair<std::vector<GeometryTag>,std::vector<GeometryTag>> apair = allContactsPairs[indexOfContactPair];
                masterTags.insert(masterTags.end(),apair.first.begin(),apair.first.end());
                slaveTags.insert(slaveTags.end(),apair.second.begin(),apair.second.end());
            }

            if(masterTags.size()==0) { cerr<<"____EMPTY MASTER TAGS____"<<endl; exit(1); }
            if(slaveTags.size()==0) { cerr<<"____EMPTY SLAVE TAGS____"<<endl; exit(2); }

            //! --------------------------
            //! define a new contact pair
            //! --------------------------
            std::pair<std::vector<GeometryTag>,std::vector<GeometryTag>> aContact;
            aContact.first = masterTags;
            aContact.second = slaveTags;
            reordered.push_back(aContact);
        }
    }
        break;

    //! -----------------------
    //! no grouping ("case 5")
    //! -----------------------
    default:
    {
        reordered = allContactsPairs;
    }
        break;
    }
}



//! ------------------------------
//! function: getPathOfExecutable
//! details:
//! ------------------------------
#include <Windows.h>
#include <codecvt>
std::string contactFinder::getPathOfExecutable()
{
    using convert_typeX = std::codecvt_utf8<wchar_t>;
    std::wstring_convert<convert_typeX, wchar_t> converterX;
    wchar_t buff[MAX_PATH+10]={0};
    GetModuleFileNameW(nullptr, buff, MAX_PATH);
    std::wstring fulllocation(buff);
    int last_slash = (int)fulllocation.find_last_of(L"\\");
    fulllocation = fulllocation.substr(0, last_slash);
    return converterX.to_bytes(fulllocation);
}
