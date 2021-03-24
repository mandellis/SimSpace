//! ----------------
//! custom includes
//! ----------------
#include "surfacemeshtofacemeshes.h"
#include <StlMesh_Mesh.hxx>
#include <ng_meshvs_datasourceface.h>
#include <ng_meshvs_datasource2d.h>
#include <extendedrwstl.h>
#include <stlapiwriter.h>
#include <stldoctor.h>

//! ---
//! Qt
//! ---
#include <QList>
#include <QMap>

//! ----
//! OCC
//! ----
#include <TopExp.hxx>
#include <TopAbs_ShapeEnum.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TColStd_MapIteratorOfPackedMapOfInteger.hxx>
#include <MeshVS_EntityType.hxx>
#include <BRepTools.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopExp.hxx>
#include <algorithm>

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
surfaceMeshToFaceMeshes::surfaceMeshToFaceMeshes(const occHandle(Ng_MeshVS_DataSource2D) &aSurfaceMeshDS, const TopoDS_Shape &aShape)
{
    myMeshDS = aSurfaceMeshDS;
    myShape = aShape;
}

//! --------------------------------------------------------------
//! function: initRayCaster
//! details:  init the embree ray caster using a mesh data source
//! --------------------------------------------------------------
void surfaceMeshToFaceMeshes::initRayCaster(const occHandle(Ng_MeshVS_DataSource2D) &aSurfaceMeshDS)
{
    cout<<"surfaceMeshToFaceMeshes::initRayCaster()->____function called____"<<endl;

    //! -------------
    //! sanity check
    //! -------------
    if(aSurfaceMeshDS.IsNull()) return;
    if(aSurfaceMeshDS->GetAllElements().Extent()<1) return;
    if(aSurfaceMeshDS->GetAllNodes().Extent()<3) return;

    //! -------------------
    //! nodes and elements
    //! -------------------
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    //! --------------------------
    //! convert into eigen format
    //! --------------------------
    int NN = aSurfaceMeshDS->GetAllNodes().Extent();
    V.resize(NN,3);

    TColStd_MapIteratorOfPackedMapOfInteger it;
    int NbNodes;
    MeshVS_EntityType type;
    double buf[3];
    TColStd_Array1OfReal coords(*buf,1,3);
    for(it.Initialize(aSurfaceMeshDS->GetAllNodes());it.More();it.Next())
    {
        int globalID = it.Key();
        int localNodeID = aSurfaceMeshDS->myNodesMap.FindIndex(globalID);
        aSurfaceMeshDS->GetGeom(globalID,false,coords,NbNodes,type);
        int pos = localNodeID-1;
        V(pos,0) = coords(1);
        V(pos,1) = coords(2);
        V(pos,2) = coords(3);
    }

    //! ----------------
    //! adding elements
    //! ----------------
    int NE = aSurfaceMeshDS->GetAllElements().Extent();
    F.resize(NE,3);

    int bufn[3];
    TColStd_Array1OfInteger nodeIDs(*bufn,1,3);     //! only triangles
    for(it.Initialize(aSurfaceMeshDS->GetAllElements());it.More();it.Next())
    {
        int globalElementID = it.Key();
        int localElementID = aSurfaceMeshDS->myElementsMap.FindIndex(globalElementID);
        aSurfaceMeshDS->GetNodesByElement(globalElementID,nodeIDs,NbNodes);

        int pos = localElementID-1;
        for(int i=1; i<=NbNodes; i++)
        {
            int curGlobalNodeID = nodeIDs.Value(i);
            int curLocalNodeID = aSurfaceMeshDS->myNodesMap.FindIndex(curGlobalNodeID);
            F(pos,i-1) = curLocalNodeID-1;
        }
    }

    //! -------------------------------
    //! store the mesh and init embree
    //! -------------------------------
    myV = V;
    myF = F;
    myEmbree.init(V.cast<float>(),F.cast<int>());

    cout<<"surfaceMeshToFaceMeshes::initRayCaster()->____raycaster initialized____"<<endl;
}

//! -------------------
//! function: setShape
//! details:
//! -------------------
void surfaceMeshToFaceMeshes::setShape(const TopoDS_Shape &aShape)
{
    myShape = aShape;
}

//! ------------------
//! function: setMesh
//! details:
//! ------------------
void surfaceMeshToFaceMeshes::setMesh(const opencascade::handle<Ng_MeshVS_DataSource2D> &aMeshDS)
{
    myMeshDS = aMeshDS;
}

//! ------------------------------
//! function; getPathOfExecutable
//! details:  helper
//! ------------------------------
std::string surfaceMeshToFaceMeshes::getPathOfExecutable()
{
    LPWSTR buffer;
    GetModuleFileName(NULL, buffer, MAX_PATH);
    std::string aString;
    aString.reserve(wcslen(buffer));
    for (;*buffer; buffer++) aString += (char)*buffer;
    std::string::size_type pos = aString.find_last_of("\\/");
    if (pos == std::string::npos) return "";
    return aString.substr(0, pos);
}

//! ------------------------------------------------------------------------------
//! function: buildShapeTessellation
//! details:  retrieve from the shape the map (faceNr, stl face mesh data source)
//!           and init the ray caster with the STL surface mesh
//! ------------------------------------------------------------------------------
bool surfaceMeshToFaceMeshes::buildShapeTessellation(const TopoDS_Shape &aShape,
                                                     occHandle(Ng_MeshVS_DataSource2D) &surfaceMeshDS_STL,
                                                     std::map<int,occHandle(Ng_MeshVS_DataSourceFace)> &mapOfFaceMeshDS_STL)
{
    cout<<"surfaceMeshToFaceMeshes::buildShapeTessellation()->____function called____"<<endl;

    //! -----------------------------------------
    //! write the shape tessellation on disk
    //! (before generate a quality tessellation)
    //! -----------------------------------------
    //std::string fileName = this->getPathOfExecutable()+"\\tessellationExtended.stl";
    std::string fileName("D:/tessellationExtended.stl");
    cout<<"surfaceMeshToFaceMeshes::buildShapeTessellation()->____tessellation file path: "<<fileName<<"____"<<endl;
    STLAPIWriter aWriter;

    BRepTools::Clean(aShape);
    bool isInParallel = true;
    bool isRelative = true;
    BRepMesh_IncrementalMesh tessellator(aShape,0.01,isRelative,0.1,isInParallel);
    aWriter.WriteExtended(aShape,fileName.c_str());

    QString outputFileName;
    bool healTessellation = false;
    if(healTessellation)
    {
        //! ------------------------------------------------------
        //! heal the tessellation with one of the available tools
        //! ------------------------------------------------------
        STLdoctor aSTLDoctor;
        outputFileName = QString("D:/tessellationExtendedAndCorrected.stl");
        cout<<"surfaceMeshToFaceMeshes::buildShapeTessellation()->____output file path: "<<fileName<<"____"<<endl;
        //aSTLDoctor.performMeshFix(QString::fromStdString(fileName),outputFileName);
        aSTLDoctor.perform(QString::fromStdString(fileName),outputFileName);
    }

    //! ---------------------------
    //! read the extended stl file
    //! ---------------------------
    TopTools_IndexedMapOfShape faceMap;
    TopExp::MapShapes(aShape,TopAbs_FACE,faceMap);
    std::vector<occHandle(StlMesh_Mesh)> vecFaceStlMesh;
    int NbFaces = faceMap.Extent();
    for(int faceNr=0; faceNr<=NbFaces; faceNr++)
    {
        occHandle(StlMesh_Mesh) aStlMesh = new StlMesh_Mesh();
        vecFaceStlMesh.push_back(aStlMesh);
    }
    QList<QMap<int,int>> maps_localToGlobal_nodeIDs;
    QList<QMap<int,int>> maps_localToGlobal_elementIDs;

    occHandle(StlMesh_Mesh) overallSTL;
    if(healTessellation)
    {
        overallSTL = ExtendedRWStl::ReadExtendedSTLAscii(outputFileName,
                                                         vecFaceStlMesh,
                                                         maps_localToGlobal_nodeIDs,
                                                         maps_localToGlobal_elementIDs);
    }
    else
    {
         overallSTL = ExtendedRWStl::ReadExtendedSTLAscii(QString::fromStdString(fileName),
                                                                                 vecFaceStlMesh,
                                                                                 maps_localToGlobal_nodeIDs,
                                                                                 maps_localToGlobal_elementIDs);
    }


    if(overallSTL.IsNull()) return false;
    if(overallSTL->NbTriangles()<1 || overallSTL->NbVertices()<3) return false;

    //! -----------------------------------
    //! build the surface mesh data source
    //! -----------------------------------
    surfaceMeshDS_STL = new Ng_MeshVS_DataSource2D(overallSTL);
    if(surfaceMeshDS_STL.IsNull()) return false;

    //! -----------------------------------------------------
    //! generate face mesh data sources using the stl meshes
    //! and scan each face mesh data source: build a map of
    //! stl triangles (each triangle is a triad of integers)
    //! versus the face number
    //! -----------------------------------------------------
    for(int faceNr = 1; faceNr<=NbFaces; faceNr++)
    {
        occHandle(StlMesh_Mesh) curStlMesh = vecFaceStlMesh.at(faceNr);
        if(curStlMesh.IsNull())
        {
            mapOfFaceMeshDS_STL.insert(std::make_pair(faceNr,occHandle(Ng_MeshVS_DataSourceFace)()));
            continue;
        }
        occHandle(Ng_MeshVS_DataSourceFace) faceMeshDS = new Ng_MeshVS_DataSourceFace(curStlMesh,maps_localToGlobal_nodeIDs,maps_localToGlobal_elementIDs,faceNr);
        mapOfFaceMeshDS_STL.insert(std::make_pair(faceNr,faceMeshDS));
        //cout<<"face nr: "<<faceNr<<" nodes: "<<faceMeshDS->GetAllNodes().Extent()<<" elements: "<<faceMeshDS->GetAllNodes().Extent()<<"____"<<endl;
    }
    cout<<"surfaceMeshToFaceMeshes::buildShapeTessellation()->____exiting function____"<<endl;
    return true;
}

//! ------------------
//! function: perform
//! details:
//! ------------------
bool surfaceMeshToFaceMeshes::perform(std::map<int,occHandle(Ng_MeshVS_DataSourceFace)> &theFaceMeshDataSources)
{
    cout<<"surfaceMeshToFaceMeshes::perform()->____function called____"<<endl;

    //! -------------------------
    //! point to face number map
    //! -------------------------
    std::map<mesh::meshPoint,std::vector<int>> pointsToFaceNumber;

    //! -----------------------------------
    //! map (faceNr, STL face data source)
    //! -----------------------------------
    occHandle(Ng_MeshVS_DataSource2D) STLMeshDS;
    std::map<int,occHandle(Ng_MeshVS_DataSourceFace)> mapOfFaceMeshDS_STL;
    bool done = this->buildShapeTessellation(myShape,STLMeshDS,mapOfFaceMeshDS_STL);
    if(done == false)
    {
        cout<<"surfaceMeshToFaceMeshes::perform()->____tessellation not done____"<<endl;
        return false;
    }

    //! --------------------
    //! init the ray caster
    //! --------------------
    this->initRayCaster(STLMeshDS);

    //! -----------------------------------------
    //! build a map (meshElement2D, face number)
    //! -----------------------------------------
    std::map<meshElement2D,int> mapTriangleFaceNr;
    for(std::map<int,occHandle(Ng_MeshVS_DataSourceFace)>::iterator it = mapOfFaceMeshDS_STL.begin(); it!=mapOfFaceMeshDS_STL.end(); it++)
    {
        int faceNr = it->first;
        const occHandle(Ng_MeshVS_DataSourceFace) &aSTLFaceMeshDS = it->second;
        for(TColStd_MapIteratorOfPackedMapOfInteger it_(aSTLFaceMeshDS->GetAllElements()); it_.More(); it_.Next())
        {
            int globalElementID = it_.Key();
            int NbNodes, buf[8];
            TColStd_Array1OfInteger nodeIDs(*buf,1,8);
            aSTLFaceMeshDS->GetNodesByElement(globalElementID,nodeIDs,NbNodes);
            meshElement2D aMeshElement;
            aMeshElement.ID = globalElementID;
            for(int n=1; n<=NbNodes; n++) aMeshElement.nodeIDs<<nodeIDs(n);
            mapTriangleFaceNr.insert(std::make_pair(aMeshElement,faceNr));
        }
    }

    //! -----------------------------------------------------------------------------
    //! iterate over the elements of the FEM surface mesh
    //! 1. find the center of the triangular element
    //! 2. shift outside the surface mesh by a fixed amount
    //! 3. cast a ray and found the intersected STL triangle
    //! 4. since STL triangles are tagged and "mapTriangleFaceNr" has been generated
    //!    the face number is known
    //! 5. fill the map "elementsOfTheFaces" (face number, vector of face elements)
    //! -----------------------------------------------------------------------------
    int NbElements_FEM_found = 0;
    int NbElements_FEM_NotFound = 0;
    int NbPoints_FEM_mesh = myMeshDS->GetAllElements().Extent();

    //! ----------------------------------------------------------
    //! the final result
    //! - elements of the faces
    //! - elements that cannot be associated using this algorithm
    //! ----------------------------------------------------------
    std::map<int,std::vector<meshElementByCoords>> elementsOfTheFaces;
    std::vector<meshElementByCoords> unassociatedElements;

    for(TColStd_MapIteratorOfPackedMapOfInteger it(myMeshDS->GetAllElements()); it.More(); it.Next())
    {
        int globalElementID = it.Key();
        int NbNodes;
        double bufd[24];
        TColStd_Array1OfReal coords(*bufd,1,24);
        MeshVS_EntityType aType;
        myMeshDS->GetGeom(globalElementID,true,coords,NbNodes,aType);

        //! ----------------------
        //! center of the element
        //! ----------------------
        double Cx, Cy,Cz;
        Cx = Cy = Cz = 0;
        for(int n=0; n<NbNodes; n++)
        {
            int s = 3*n;
            Cx += coords(s+1);
            Cy += coords(s+2);
            Cz += coords(s+3);
        }
        Cx /= 3.0; Cy /= 3.0; Cz /= 3.0;

        //! ---------------------------------------
        //! normal at the current FEM mesh element
        //! ---------------------------------------
        double nx,ny,nz;
        myMeshDS->GetNormal(globalElementID,10,nx,ny,nz);

        //! --------------------------------------------------
        //! shift the point outside the center of the element
        //! cast a ray from outside towards the STL mesh
        //! . compute a proper shift
        //! --------------------------------------------------
        Eigen::RowVector3f origin(Cx,Cy,Cz);
        Eigen::RowVector3f direction0(nx,ny,nz);    // this normal is directed towards the exterior
        Eigen::RowVector3f direction1(-nx,-ny,-nz); // this normal is directed towards the interior
        igl::Hit hit0, hit1;
        bool isDone0 = myEmbree.intersectRay(origin.cast<float>(),direction0.cast<float>(),hit0,0.0,100.0);
        bool isDone1 = myEmbree.intersectRay(origin.cast<float>(),direction1.cast<float>(),hit1,0.0,100.0);

        int intersectedElement = -1;    // dummy initialization
        if(isDone0 == true && isDone1 == false)
        {
            intersectedElement = hit0.id;
        }
        else if(isDone0 == false && isDone1 == true)
        {
            intersectedElement = hit1.id;
        }
        else if(isDone0 == true && isDone1 == true)
        {
            if(fabs(hit0.t)<=fabs(hit1.t)) intersectedElement = hit0.id;
            else intersectedElement = hit1.id;
        }
        else
        {
            //! -----------------------------------------------------------
            //! this element cannot be associated through ray intersection
            //! because the ray inersection with the an stl triangle does
            //! not exist
            //! -----------------------------------------------------------
            NbElements_FEM_NotFound++;
            meshElementByCoords aMeshElement;
            aMeshElement.ID = globalElementID;
            int bufn[8], NbNodes__;
            TColStd_Array1OfInteger nodeIDs(*bufn,1,8);
            myMeshDS->GetNodesByElement(globalElementID,nodeIDs,NbNodes__);
            for(int n=0; n<NbNodes__;n++)
            {
                int s = 3*n;
                mesh::meshPoint aMeshPoint(coords(s+1),coords(s+2),coords(s+3),nodeIDs(n+1));
                aMeshElement.pointList<<aMeshPoint;
            }
            unassociatedElements.push_back(aMeshElement);
            continue;       // must continue!
        }

        int localElementID_STL = intersectedElement+1;                                      // local ID of the intersected element within the OCC framework
        int globalElementID_STL = STLMeshDS->myElementsMap.FindKey(localElementID_STL);     // global ID of the intersected element within the OCC framework

        //cout<<"surfaceMeshToFaceMeshes::perform()->____intersected element: ("<<localElementID_STL<<", "<<globalElementID_STL<<")____"<<endl;

        int NbNodes_, buf[8];
        TColStd_Array1OfInteger nodeIDs(*buf,1,8);
        STLMeshDS->GetNodesByElement(globalElementID_STL,nodeIDs,NbNodes_);
        meshElement2D aMeshElement;
        aMeshElement.ID = globalElementID_STL;
        for(int n=1; n<=NbNodes_; n++) aMeshElement.nodeIDs<<nodeIDs(n);

        std::map<meshElement2D,int>::iterator ite = mapTriangleFaceNr.find(aMeshElement);
        if(ite==mapTriangleFaceNr.end())
        {
            //! -----------------------------------------------------
            //! the mesh triangle cannot be associated with any face
            //! for some strange reason - should never occur
            //! -----------------------------------------------------
            meshElementByCoords aMeshElement;
            aMeshElement.ID = globalElementID;
            int NbNodes_, bufn[8];
            TColStd_Array1OfInteger nodeIDs_(*bufn,1,8);
            myMeshDS->GetNodesByElement(globalElementID,nodeIDs_,NbNodes_);

            for(int n=0; n<NbNodes_;n++)
            {
                int s = 3*n;
                mesh::meshPoint aMeshPoint(coords(s+1),coords(s+2),coords(s+3),nodeIDs_(n+1));
                aMeshElement.pointList<<aMeshPoint;
            }
            unassociatedElements.push_back(aMeshElement);

            NbElements_FEM_NotFound++;
            cout<<"surfaceMeshToFaceMeshes::perform()->____raycast failure: the current element will be reassigned using brute force____"<<endl;
            continue;
        }

        //! --------------------------------------------------------
        //! the mesh element can be associated with the face number
        //! --------------------------------------------------------
        NbElements_FEM_found++;

        int faceNr = ite->second;
        meshElementByCoords aMeshElementByCoords;
        aMeshElementByCoords.ID = globalElementID;

        myMeshDS->GetNodesByElement(globalElementID,nodeIDs,NbNodes);
        for(int n=0; n<NbNodes; n++)
        {
            int globalNodeID = nodeIDs(n+1);
            int s = 3*n;
            mesh::meshPoint aMeshPoint(coords(s+1),coords(s+2),coords(s+3),globalNodeID);
            aMeshElementByCoords.pointList<<aMeshPoint;
        }
        std::map<int,std::vector<meshElementByCoords>>::iterator itt = elementsOfTheFaces.find(faceNr);
        if(itt==elementsOfTheFaces.end())
        {
            std::vector<meshElementByCoords> v { aMeshElementByCoords };
            elementsOfTheFaces.insert(std::make_pair(faceNr,v));
        }
        else itt->second.push_back(aMeshElementByCoords);
    }

    cout<<"surfaceMeshToFaceMeshes::perform()->____stage 1 - heuristic____"<<endl;
    cout<<"surfaceMeshToFaceMeshes::perform()->____number of surface mesh elements "<<NbPoints_FEM_mesh<<"____"<<endl;
    cout<<"surfaceMeshToFaceMeshes::perform()->____number of assigned elements: "<<NbElements_FEM_found<<"____"<<endl;
    cout<<"surfaceMeshToFaceMeshes::perform()->____number of unassigned elements: "<<NbElements_FEM_NotFound<<"____"<<endl;

    //! ---------------------------------
    //! handling "lost" surface elements
    //! 1. search by segments
    //! ---------------------------------
    if(unassociatedElements.empty()==false)
    {
        cout<<"surfaceMeshToFaceMeshes::perform()->____stage 2 - brute force____"<<endl;

        //! ---------------------------------
        //! mesh segments to face number map
        //! ---------------------------------
        cout<<"surfaceMeshToFaceMeshes::perform()->____generating face to segments map____"<<endl;
        std::map<mesh::meshSegment,std::vector<int>> segmentFaceNrMap;
        for(std::map<int,occHandle(Ng_MeshVS_DataSourceFace)>::iterator it = theFaceMeshDataSources.begin(); it!=theFaceMeshDataSources.end(); it++)
        {
            int faceNr = it->first;
            occHandle(Ng_MeshVS_DataSourceFace) aFaceMeshDS = it->second;

            if(aFaceMeshDS->mySegmentToElement.isEmpty()) aFaceMeshDS->computeFreeMeshSegments();
            const QMap<mesh::meshSegment,QList<int>> &segmentsToElement = aFaceMeshDS->mySegmentToElement;
            for(QMap<mesh::meshSegment,QList<int>>::const_iterator it_ = segmentsToElement.cbegin(); it_!=segmentsToElement.cend(); it_++)
            {
                mesh::meshSegment aSegment = it_.key();
                std::map<mesh::meshSegment,std::vector<int>>::iterator it__ = segmentFaceNrMap.find(aSegment);
                if(it__==segmentFaceNrMap.end())
                {
                    std::vector<int> v { faceNr };
                    segmentFaceNrMap.insert(std::make_pair(aSegment,v));
                }
                else it__->second.push_back(faceNr);
            }
        }
        cout<<"surfaceMeshToFaceMeshes::perform()->____face to segments map generated____"<<endl;

        //! ---------------------------------------
        //! start reassigning using segment search
        //! ---------------------------------------
        cout<<"surfaceMeshToFaceMeshes::perform()->____start reassociating: "<<unassociatedElements.size()<<" lost elements____"<<endl;
        int NbReassigned = 0;
        for(std::vector<meshElementByCoords>::iterator it = unassociatedElements.begin(); it!=unassociatedElements.end(); /*it++*/)
        {
            const meshElementByCoords &curElement = *it;
            cout<<"surfaceMeshToFaceMeshes::perform()->____trying to assign element ID: "<<curElement.ID<<"____"<<endl;
            int NbNodes = curElement.pointList.length();

            //! ---------------------------------------------
            //! store the segments of the unassigned element
            //! ---------------------------------------------
            std::vector<mesh::meshSegment> vecSegmentsOfElement;    // segments of the unassigned element
            for(int i=1; i<=NbNodes; i++)
            {
                mesh::meshSegment aMeshSegment;
                int firstIndex = (i-1)%NbNodes;
                int secondIndex = i%NbNodes;
                aMeshSegment.nodeIDs<<curElement.pointList[firstIndex].ID;
                aMeshSegment.nodeIDs<<curElement.pointList[secondIndex].ID;
                vecSegmentsOfElement.push_back(aMeshSegment);
            }

            //! --------------------------------------------------------------------
            //! check the position of the current segment into the segmentFaceNrMap
            //! => find one of the possibile faces of the element
            //! --------------------------------------------------------------------
            bool currentElementReassigned = false;
            for(std::vector<mesh::meshSegment>::iterator it_ = vecSegmentsOfElement.begin(); it_ !=vecSegmentsOfElement.end(); it_++)
            {
                const mesh::meshSegment &curSegment = *it_;
                std::map<mesh::meshSegment,std::vector<int>>::iterator its = segmentFaceNrMap.find(curSegment);

                cout<<"surfaceMeshToFaceMeshes::perform()->____check position of segment ("<<curSegment.nodeIDs[0]<<", "<<
                      curSegment.nodeIDs[1]<<")";

                if(its!=segmentFaceNrMap.end())
                {
                    cout<<": found____"<<endl;
                    const std::vector<int> &possibleFaceNrs = its->second;
                    int faceNr = possibleFaceNrs[0];
                    std::map<int,std::vector<meshElementByCoords>>::iterator ite = elementsOfTheFaces.find(faceNr);
                    if(ite == elementsOfTheFaces.end())
                    {
                        std::vector<meshElementByCoords> v {curElement};
                        elementsOfTheFaces.insert(std::make_pair(faceNr,v));
                    }
                    else ite->second.push_back(curElement);

                    //! ------------------------------------------------------------------
                    //! the element has been reassigned to the face "faceNr"
                    //! before exiting add the other two segments to the segmentFaceNrMap
                    //! ------------------------------------------------------------------
                    for(size_t n=0; n<vecSegmentsOfElement.size(); n++)
                    {
                        const mesh::meshSegment &aSegment = vecSegmentsOfElement[n];
                        if(aSegment == curSegment)
                        {
                            cout<<"____found the already added segment: continue____"<<endl;
                        }
                        std::map<mesh::meshSegment,std::vector<int>>::iterator itss = segmentFaceNrMap.find(aSegment);
                        if(itss==segmentFaceNrMap.end())
                        {
                            cout<<"____adding____"<<endl;
                            std::vector<int> v { faceNr };
                            segmentFaceNrMap.insert(std::make_pair(v,faceNr));
                        }
                        else
                        {
                            cout<<"____piling____"<<endl;
                            itss->second.push_back(faceNr);
                        }
                    }

                    NbReassigned++;
                    currentElementReassigned = true;
                    break;
                }
                else
                {
                    cout<<": not found____"<<endl;
                }
            }
            if(currentElementReassigned) it = unassociatedElements.erase(it);
            else it++;
        }

        //! -------------------------------------------------------------------
        //! if needed try to associate using another strategy
        //! use std::map<mesh::meshPoint,std::vector<int>> pointsToFaceNumber;
        //! -------------------------------------------------------------------
        for(std::vector<meshElementByCoords>::iterator it = unassociatedElements.begin(); it!=unassociatedElements.end(); /*it++*/)
        {
            if(pointsToFaceNumber.empty())
            {
                cout<<"surfaceMeshToFaceMeshes::perform()->____generating mesh point to face number map____"<<endl;
                for(std::map<int,std::vector<meshElementByCoords>>::iterator it1 = elementsOfTheFaces.begin(); it1 != elementsOfTheFaces.end(); it1++)
                {
                    int faceNumber = it1->first;
                    const std::vector<meshElementByCoords> &vecMeshElements = it1->second;
                    for(std::vector<meshElementByCoords>::const_iterator it2 =vecMeshElements.cbegin(); it2!=vecMeshElements.cend(); it2++)
                    {
                        const meshElementByCoords &aMeshElement = *it2;
                        for(int j = 1; j<aMeshElement.pointList.length(); j++)
                        {
                            const mesh::meshPoint &aMeshPoint = aMeshElement.pointList[j];
                            std::map<mesh::meshPoint,std::vector<int>>::iterator it3 = pointsToFaceNumber.find(aMeshPoint);
                            if(it3==pointsToFaceNumber.end())
                            {
                                std::vector<int> v { faceNumber };
                                pointsToFaceNumber.insert(std::make_pair(aMeshPoint, v));
                            }
                            else
                            {
                                if(std::find(it3->second.begin(), it3->second.end(), faceNumber)==it3->second.end())
                                    it3->second.push_back(faceNumber);
                            }
                        }
                    }
                }
                cout<<"surfaceMeshToFaceMeshes::perform()->____mesh point to face number map generated____"<<endl;
            }

            bool reassigned = false;
            const meshElementByCoords &curPendingElement = *it;
            for(int j=0; j<curPendingElement.pointList.length(); j++)
            {
                const mesh::meshPoint &aMeshPoint = curPendingElement.pointList[j];
                std::map<mesh::meshPoint,std::vector<int>>::iterator itf = pointsToFaceNumber.find(aMeshPoint);
                if(itf != pointsToFaceNumber.end())
                {
                    std::vector<int> vecfaceNr = itf->second;
                    int faceNr = vecfaceNr[0];  // take the first...
                    std::map<int,std::vector<meshElementByCoords>>::iterator ite = elementsOfTheFaces.find(faceNr);
                    if(ite==elementsOfTheFaces.end())
                    {
                        std::vector<meshElementByCoords> v {curPendingElement};
                        elementsOfTheFaces.insert(std::make_pair(faceNr,v));
                    }
                    else ite->second.push_back(curPendingElement);

                    //! -------------------------------------
                    //! other ... to complete
                    //! 1. update the pointsToFaceNumber map
                    //! -------------------------------------
                    for(int k=0; k<curPendingElement.pointList.length(); k++)
                    {
                        const mesh::meshPoint &aPoint = curPendingElement.pointList[k];
                        if(aPoint==aMeshPoint) continue;
                        std::map<mesh::meshPoint,std::vector<int>>::iterator itff = pointsToFaceNumber.find(aPoint);
                        if(itff==pointsToFaceNumber.end())
                        {
                            std::vector<int> v { faceNr };
                            pointsToFaceNumber.insert(std::make_pair(aPoint,v));
                        }
                        else itff->second.push_back(faceNr);
                    }

                    it = unassociatedElements.erase(it);
                    reassigned = true;
                    NbReassigned++;
                    break;
                }
                else it++;
            }
        }

        cout<<"surfaceMeshToFaceMeshes::perform()->____# "<<unassociatedElements.size()<<" have not been reassigned____"<<endl;
    }
    if(elementsOfTheFaces.empty()) exit(3000);
    //! ---------------------------------
    //! build the face mesh data sources
    //! ---------------------------------
    for(std::map<int,std::vector<meshElementByCoords>>::iterator it = elementsOfTheFaces.begin(); it!=elementsOfTheFaces.end(); it++)
    {
        int faceNr = it->first;
        const std::vector<meshElementByCoords> &elements = it->second;
        occHandle(Ng_MeshVS_DataSourceFace) aFaceMeshDS = new Ng_MeshVS_DataSourceFace(elements,false,false);
        theFaceMeshDataSources.insert(std::make_pair(faceNr,aFaceMeshDS));
    }
    if(unassociatedElements.size()==0) return false;
    return true;
}
