#include "surfacemeshtofacemeshes.h"
#include <extendedrwstl.h>
#include <stlapiwriter.h>
#include <QList>
#include <QMap>
#include <StlMesh_Mesh.hxx>
#include <TopExp.hxx>
#include <TopAbs_ShapeEnum.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TColStd_MapIteratorOfPackedMapOfInteger.hxx>
#include <MeshVS_EntityType.hxx>
#include <ng_meshvs_datasourceface.h>
#include <ng_meshvs_datasource2d.h>

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

    //! ----------------------------------------------------------------
    //! convert the shape into a tessellation with tags written on disk
    //! ----------------------------------------------------------------
    std::string fileName("D:/support.stl");
    STLAPIWriter aWriter;
    aWriter.WriteExtended(aShape,fileName.c_str());

    //! ------------------------------------------------
    //! init a vector of (NbFaces+1) empty StlMesh_Mesh
    //! ------------------------------------------------
    TopTools_IndexedMapOfShape faceMap;
    TopExp::MapShapes(aShape,TopAbs_FACE,faceMap);
    std::vector<occHandle(StlMesh_Mesh)> vecFaceStlMesh;
    int NbFaces = faceMap.Extent();
    for(int faceNr=0; faceNr<=NbFaces; faceNr++)
    {
        occHandle(StlMesh_Mesh) aStlMesh = new StlMesh_Mesh();
        vecFaceStlMesh.push_back(aStlMesh);
    }

    //! -----------------------------------------------------------------
    //! read the tessellation from disk - build the overall surface mesh
    //! and the face mesh data sources
    //! -----------------------------------------------------------------
    QString surfaceMeshFilePath("D:/support.stl");
    QList<QMap<int,int>> maps_localToGlobal_nodeIDs;
    QList<QMap<int,int>> maps_localToGlobal_elementIDs;
    occHandle(StlMesh_Mesh) overallSTL = ExtendedRWStl::ReadExtendedSTLAscii(surfaceMeshFilePath,
                                                                             vecFaceStlMesh,
                                                                             maps_localToGlobal_nodeIDs,
                                                                             maps_localToGlobal_elementIDs);

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

    //! -----------------------------------
    //! map (faceNr, STL face data source)
    //! -----------------------------------
    occHandle(Ng_MeshVS_DataSource2D) STLMeshDS;
    std::map<int,occHandle(Ng_MeshVS_DataSourceFace)> mapOfFaceMeshDS_STL;
    bool done = this->buildShapeTessellation(myShape,STLMeshDS,mapOfFaceMeshDS_STL);
    if(done == false) return false;

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

    //! -----------------------------------------------------
    //! iterate over the elements of the FEM surface mesh
    //! 1. find the center of the triangular element
    //! 2. shift outside the surface mesh by a fixed amount
    //! 3. cast a ray and found the intersected STL triangle
    //! 4.
    //! -----------------------------------------------------
    int NbPoints_FEM_meshFound = 0;
    int NbPoints_FEM_meshNotFound = 0;
    int NbPoints_FEM_mesh = myMeshDS->GetAllElements().Extent();

    std::map<int,std::vector<meshElementByCoords>> elementsOfTheFaces;
    for(TColStd_MapIteratorOfPackedMapOfInteger it(myMeshDS->GetAllElements()); it.More(); it.Next())
    {
        int globalElementID = it.Key();
        int NbNodes;
        double buf[24];
        TColStd_Array1OfReal coords(*buf,1,24);
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
        Eigen::RowVector3f origin0(Cx,Cy,Cz);
        Eigen::RowVector3f direction0(nx,ny,nz);
        igl::Hit hit0;
        bool isDone = myEmbree.intersectRay(origin0.cast<float>(),direction0.cast<float>(),hit0,0.0);
        double eps = 10.0;    // at this stage this is an arbitrary value
        if(isDone)
        {
            double distanceFromIntersection = hit0.v;
            eps = distanceFromIntersection/2.0;
        }

        double x = Cx+nx*eps;
        double y = Cy+ny*eps;
        double z = Cz+nz*eps;
        Eigen::RowVector3f origin (x,y,z);
        Eigen::RowVector3f direction(-nx,-ny,-nz);
        igl::Hit hit;
        double tnear = 0.0;
        isDone = myEmbree.intersectRay(origin.cast<float>(),direction.cast<float>(),hit,tnear);
        if(isDone == true)
        {
            int intersectedElement = hit.id;                                                    // ID of the intersected element within the eigen framework
            int localElementID_STL = intersectedElement+1;                                      // local ID of the intersected element within the OCC framework
            int globalElementID_STL = STLMeshDS->myElementsMap.FindKey(localElementID_STL);     // global ID of the intersected element within the OCC framework

            //cout<<"surfaceMeshToFaceMeshes::perform()->____intersected element: "<<intersectedElement<<"____"<<endl;
            //cout<<"surfaceMeshToFaceMeshes::perform()->____intersected element: "<<localElementID_STL<<"____"<<endl;
            //cout<<"surfaceMeshToFaceMeshes::perform()->____intersected element: "<<globalElementID_STL<<"____"<<endl;

            int NbNodes, buf[8];
            TColStd_Array1OfInteger nodeIDs(*buf,1,8);
            STLMeshDS->GetNodesByElement(globalElementID_STL,nodeIDs,NbNodes);
            meshElement2D aMeshElement;
            aMeshElement.ID = globalElementID_STL;
            for(int n=1; n<=NbNodes; n++) aMeshElement.nodeIDs<<nodeIDs(n);

            std::map<meshElement2D,int>::iterator ite = mapTriangleFaceNr.find(aMeshElement);
            if(ite==mapTriangleFaceNr.end())
            {
                //! -----------------------------------------------------
                //! the mesh triangle cannot be associated with any face
                //! -----------------------------------------------------
                NbPoints_FEM_meshNotFound++;
                continue;
            }
            //! --------------------------------------------------------
            //! the mesh element can be associated with the face number
            //! --------------------------------------------------------
            NbPoints_FEM_meshFound++;

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
    }
    cout<<"surfaceMeshToFaceMeshes::perform()->____total number of the surface mesh: "<<NbPoints_FEM_mesh<<"____"<<endl;
    cout<<"surfaceMeshToFaceMeshes::perform()->____found: "<<NbPoints_FEM_meshFound<<"____"<<endl;
    cout<<"surfaceMeshToFaceMeshes::perform()->____found: "<<NbPoints_FEM_meshNotFound<<"____"<<endl;

    //! ---------------------------------
    //! build the face mesh data sources
    //! ---------------------------------
    for(std::map<int,std::vector<meshElementByCoords>>::iterator it = elementsOfTheFaces.begin(); it!=elementsOfTheFaces.end(); it++)
    {
        int faceNr = it->first;
        cout<<"surfaceMeshToFaceMeshes::perform()->____building the face mesh data source of face nr: "<<faceNr<<"____"<<endl;
        const std::vector<meshElementByCoords> &elements = it->second;
        occHandle(Ng_MeshVS_DataSourceFace) aFaceMeshDS = new Ng_MeshVS_DataSourceFace(elements,false,false);
        theFaceMeshDataSources.insert(std::make_pair(faceNr,aFaceMeshDS));
    }
    return true;
}
