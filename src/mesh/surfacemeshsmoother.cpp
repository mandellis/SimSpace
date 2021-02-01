//! ----------------
//! custom includes
//! ----------------
#include "surfacemeshsmoother.h"
#include "qprogressevent.h"
#include "src/utils/global.h"
#include "userMessage.h"
#include <mesh.h>
#include <src/utils/meshtools.h>
#include <tetqualityclass.h>

//! ----
//! OCC
//! ----
#include <TColStd_MapIteratorOfPackedMapOfInteger.hxx>
#include <MeshVS_EntityType.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <TopoDS.hxx>
#include <StlMesh_Mesh.hxx>
#include <TopExp.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <MeshVS_DataSource.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <BRepTools.hxx>

//! ----
//! C++
//! ----
#include <vector>
#include <chrono>
#include <random>

//! ---
//! Qt
//! ---
#include <QApplication>
#include <QMap>
#include <QList>
#include <QThread>

//! ----------------
//! custom includes
//! ----------------
#include <mesh.h>
#include <ng_meshvs_datasourceface.h>
#include <ng_meshvs_datasource2d.h>
#include <ng_meshvs_datasource3d.h>
#include <stlapiwriter.h>
#include <extendedrwstl.h>

//! ------
//! Eigen
//! ------
#include <ext/eigen/Eigen/Dense>

//! -------
//! libigl
//! -------
#include <ext/libigl/include/igl/embree/EmbreeIntersector.h>

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
surfaceMeshSmoother::surfaceMeshSmoother(QProgressIndicator *aProgressIndicator, QObject *parent):QObject(parent)
{
    this->setProgressIndicator(aProgressIndicator);
}

//! -------------------------------
//! function: setProgressIndicator
//! details:
//! -------------------------------
void surfaceMeshSmoother::setProgressIndicator(QProgressIndicator *aProgressIndicator)
{
    myProgressIndicator = aProgressIndicator;
    if(myProgressIndicator!=Q_NULLPTR)
    {
        disconnect(this,SIGNAL(requestDisableStop()),myProgressIndicator,SLOT(disableStop()));
        disconnect(this,SIGNAL(requestEnableStop()),myProgressIndicator,SLOT(enableStop()));
        connect(this,SIGNAL(requestDisableStop()),myProgressIndicator,SLOT(disableStop()));
        connect(this,SIGNAL(requestEnableStop()),myProgressIndicator,SLOT(enableStop()));
    }
}

//! ------------------
//! function: perform
//! details:
//! ------------------
bool surfaceMeshSmoother::perform(meshDataBase *mDB, int bodyIndex, const double geometricTolerance)
{
    cout<<"surfaceMeshSmoother::perform()->____function called____"<<endl;

    //! ---------------------------------------------------------
    //! disable stop during this preliminary steps
    //! - note: cannot enable/disable the progress indicator
    //! directly, since it lives in the main thread: use signals
    //! ---------------------------------------------------------
    emit requestDisableStop();

    //! --------
    //! message
    //! --------
    userMessage um;
    um.isDone = false;
    um.message ="";

    //! --------------------------------
    //! retrieve the shape tessellation
    //! --------------------------------
    TopoDS_Shape aShape = mDB->bodyMap.value(bodyIndex);

    //! -----------------------------------
    //! remesh the shape with high quality
    //! so it "overlaps" the cad shape
    //! -----------------------------------
    cout<<"********************************"<<endl;
    cout<<" generating new CAD tessellator "<<endl;

    bool isInParallel = true;
    bool isRelative = true;

    BRepTools::Clean(aShape);
    BRepMesh_IncrementalMesh tessellator(aShape,0.01,isRelative,0.1,isInParallel);

    cout<<" CAD tessellator generated "<<endl;
    cout<<"********************************"<<endl;

    occHandle(Ng_MeshVS_DataSource2D) surfaceMeshDS_STL;
    std::map<int,occHandle(Ng_MeshVS_DataSourceFace)> mapOfFaceMeshDS_STL;
    bool isDone = this->buildShapeTessellation(aShape,surfaceMeshDS_STL,mapOfFaceMeshDS_STL);

    //! ------------------
    //! init a ray caster
    //! ------------------
    isDone = this->initRayCaster(surfaceMeshDS_STL);

    //! --------------------------------------------------------
    //! retrieve the surface and volume mesh from the data base
    //! compute the "best" node normals for the surface mesh
    //! --------------------------------------------------------
    occHandle(MeshVS_DataSource) aMesh = mDB->ArrayOfMeshDS2D.value(bodyIndex);
    if(aMesh.IsNull()) return false;
    occHandle(Ng_MeshVS_DataSource2D) surfaceMesh = occHandle(Ng_MeshVS_DataSource2D)::DownCast(aMesh);

    //! ----------------------------------------------
    //! a very small shift of the surface mesh points
    //! ----------------------------------------------
    //this->perturbMesh(surfaceMesh);

    surfaceMesh->computeNormalAtNodes();

    //! ------------------------
    //! set the status: running
    //! ------------------------
    Global::status().code = 1;

    //! -----------------------------------------------------------
    //! init a progress bar with the number of surface mesh points
    //! -----------------------------------------------------------
    QString task("Projecting mesh points on geometry");
    QProgressEvent *e;
    if(myProgressIndicator!=Q_NULLPTR)
    {
        int NbPoints = surfaceMesh->GetAllNodes().Extent();
        e = new QProgressEvent(QProgressEvent_None,-1,-1,-1,"Smoothing mesh",QProgressEvent_Init,1,NbPoints,0,task);
        QApplication::postEvent(myProgressIndicator,e);
        QApplication::processEvents();
        QThread::msleep(1000);
    }

    //! ------------
    //! enable stop
    //! ------------
    emit requestEnableStop();

    std::map<int,std::vector<double>> updatedPoints;
    int NbProgress = 0;
    int NbUpdated = 0;
    int NbMeshPoints = surfaceMesh->GetAllNodes().Extent();
    for(TColStd_MapIteratorOfPackedMapOfInteger it(surfaceMesh->GetAllNodes()); it.More(); it.Next())
    {
        //! ------------------
        //! check termination
        //! ------------------
        if(Global::status().code == 0)
        {
            cout<<"*********************************"<<endl;
            cout<<" process interrupted by the user "<<endl;
            cout<<"*********************************"<<endl;
            updatedPoints.clear();
            break;
        }

        int globalNodeID = it.Key();
        int NbNodes;
        double buf[3];
        TColStd_Array1OfReal coords(*buf,1,3);
        MeshVS_EntityType aType;
        surfaceMesh->GetGeom(globalNodeID,false,coords,NbNodes,aType);

        QList<double> nodeNormal;
        surfaceMesh->myNodeNormals.value(globalNodeID,nodeNormal);
        double nx = nodeNormal[0];
        double ny = nodeNormal[1];
        double nz = nodeNormal[2];

        //! -----------------------------------------------------------------------------
        //! from the current mesh point cast two rays, one towards the exterior (ray 0),
        //! one towards the interior (ray 1)
        //! -----------------------------------------------------------------------------
        Eigen::RowVector3f P(coords(1),coords(2),coords(3));

        Eigen::RowVector3f direction0(nx,ny,nz);
        igl::Hit hit0;
        bool isDone0 = myEmbree.intersectRay(P.cast<float>(),direction0.cast<float>(),hit0,0.0);

        Eigen::RowVector3f direction1(-nx,-ny,-nz);
        igl::Hit hit1;
        bool isDone1 = myEmbree.intersectRay(P.cast<float>(),direction1.cast<float>(),hit1,0.0);

        bool intersectionOK = false;
        std::vector<double> P_intersection;
        if(isDone0 == true && isDone1 == false)
        {
            intersectionOK = true;
            float s = hit0.t;
            double x = P[0]+nx*s;
            double y = P[1]+ny*s;
            double z = P[2]+nz*s;
            P_intersection.push_back(x);
            P_intersection.push_back(y);
            P_intersection.push_back(z);
            NbUpdated++;
            //cout<<"____is done 0____"<<endl;
        }
        else if(isDone0 == false && isDone1 == true)
        {
            intersectionOK = true;
            float s = hit1.t;
            double x = P[0]-nx*s;
            double y = P[1]-ny*s;
            double z = P[2]-nz*s;
            P_intersection.push_back(x);
            P_intersection.push_back(y);
            P_intersection.push_back(z);
            NbUpdated++;
            //cout<<"____is done 1____"<<endl;
        }
        else if(isDone0 == true && isDone1 == true)
        {
            intersectionOK = true;
            float s0 = hit0.t;
            float s1 = hit1.t;
            if(fabs(s0)<=fabs(s1))
            {
                double x = P[0]+nx*s0;
                double y = P[1]+ny*s0;
                double z = P[2]+nz*s0;
                P_intersection.push_back(x);
                P_intersection.push_back(y);
                P_intersection.push_back(z);
                NbUpdated++;
                //cout<<"____is done 0 & 1 => choosing 0____"<<endl;
            }
            else
            {
                double x = P[0]-nx*s1;
                double y = P[1]-ny*s1;
                double z = P[2]-nz*s1;
                P_intersection.push_back(x);
                P_intersection.push_back(y);
                P_intersection.push_back(z);
                NbUpdated++;
                //cout<<"____is done 0 & 1 => choosing 1____"<<endl;
            }
        }
        if(intersectionOK==false)
        {
            //vcout<<"____no intersection found____"<<endl;
            //continue; // no intersection found
        }
        else
        {
            updatedPoints.insert(std::make_pair(globalNodeID,P_intersection));
        }

        //! --------------------------
        //! send progress information
        //! --------------------------
        NbProgress++;
        if(myProgressIndicator!=Q_NULLPTR)
        {
            e = new QProgressEvent(QProgressEvent_None,-1,-1,-1,"Smoothing mesh",QProgressEvent_Update,1,NbMeshPoints,NbProgress,task);
            QApplication::postEvent(myProgressIndicator,e);
            QApplication::processEvents();
        }
    }

    cout<<"***************************************************"<<endl;
    cout<<" "<<updatedPoints.size()<<" projected onto the initial tessellation"<<endl;
    cout<<" "<<NbUpdated<<" projected onto the initial tessellation"<<endl;
    cout<<"***************************************************"<<endl;

    //! --------------------------
    //! no point has been shifted
    //! --------------------------
    if(updatedPoints.size()==0) return false;

    emit requestDisableStop();

    //! ---------------------------------------------------------------------
    //! iterate over the map of the displaced points change the surface mesh
    //! before displacing points store the old positions
    //! ---------------------------------------------------------------------
    std::map<int,std::vector<double>> oldPositions;
    for(std::map<int,std::vector<double>>::iterator it = updatedPoints.begin(); it != updatedPoints.end(); it++)
    {
        double Pold[3];
        this->getPointCoordinate(surfaceMesh,it->first,Pold);
        oldPositions.insert(std::make_pair(it->first,std::vector<double> { Pold[0], Pold[1], Pold[2]} ));
    }

    for(std::map<int,std::vector<double>>::iterator it = updatedPoints.begin(); it != updatedPoints.end(); /*it++*/)
    {
        int globalNodeID = it->first;
        const std::vector<double> &P_updated = it->second;
        int localNodeID = surfaceMesh->myNodesMap.FindIndex(globalNodeID);

        double Pold[3];
        this->getPointCoordinate(surfaceMesh,globalNodeID,Pold);
        double oldToNewDistance = sqrt(pow(Pold[0]-P_updated[0],2)+pow(Pold[1]-P_updated[1],2)+pow(Pold[2]-P_updated[2],2));
        if(oldToNewDistance>=geometricTolerance)    // the projected point is too far from its original position
        {
            it = updatedPoints.erase(it);
            continue;
        }
        surfaceMesh->changeNodeCoords(localNodeID,P_updated);
        it++;
    }

    //! -------------------------------------------------------------
    //! check triangles: iterate over the surface elementes attached
    //! to the current node, and check their validity
    //! -------------------------------------------------------------
    cout<<"**********************************************************"<<endl;
    cout<<" start checking surface element validity after the shifts "<<endl;

    std::map<int,std::vector<int>> nodeToSurfaceElementConnectivityMap;
    MeshTools::buildPointToElementConnectivity(surfaceMesh,nodeToSurfaceElementConnectivityMap);
    for(std::map<int,std::vector<double>>::iterator itup = updatedPoints.begin(); itup != updatedPoints.end(); /*itup++*/)
    {
        int globalNodeID_toCheck = itup->first;
        bool isSurfaceElementOK = true;
        const std::vector<int> &attachedSurfaceElements = nodeToSurfaceElementConnectivityMap.at(globalNodeID_toCheck);
        for(std::vector<int>::const_iterator itt = attachedSurfaceElements.cbegin(); itt!=attachedSurfaceElements.cend(); itt++)
        {
            int globalElementID_attached = *itt;
            int NbNodes;
            MeshVS_EntityType aType;
            double buf[24];
            TColStd_Array1OfReal coords(*buf,1,24);
            surfaceMesh->GetGeom(globalElementID_attached,true,coords,NbNodes,aType);
            double x0 = coords(1); double y0 = coords(2); double z0 = coords(3);
            double x1 = coords(4); double y1 = coords(5); double z1 = coords(6);
            double x2 = coords(6); double y2 = coords(7); double z2 = coords(8);
            double l10x = x1-x0; double l10y = y1-y0; double l10z = z1-z0;
            double l20x = x2-x0; double l20y = y2-y0; double l20z = z2-z0;
            double det = (l10y*l20z-l10z*l20y)+(l10z*l20x-l10x*l20z)+(l10x*l20y-l10y*l20x);
            if(det<=0)
            {
                cerr<<" an invalid surface element has been generated "<<endl;
                isSurfaceElementOK = false;
                itup = updatedPoints.erase(itup);
                break;
            }
            else
            {
                itup++;
            }
        }
        if(isSurfaceElementOK==false)
        {
            std::map<int,std::vector<double>>::iterator itold = oldPositions.find(globalNodeID_toCheck);
            std::vector<double> oldPositionOfPoint = itold->second;
            int localNodeID_toCheck = surfaceMesh->myNodesMap.FindIndex(globalNodeID_toCheck);
            surfaceMesh->changeNodeCoords(localNodeID_toCheck,oldPositionOfPoint);
        }
    }

    cout<<" surface element check done "<<endl;
    cout<<"**********************************************************"<<endl;

    //! ------------------------------------------------------------------
    //! if the volume mesh has been generated it is possible to check the
    //! consequences of the nodal displacements on the volume elements
    //! ------------------------------------------------------------------
    std::map<int,std::vector<int>> nodeToVolumeElementMap;
    occHandle(MeshVS_DataSource) m = mDB->ArrayOfMeshDS.value(bodyIndex);
    if(m.IsNull()==false)
    {
        occHandle(Ng_MeshVS_DataSource3D) volumeMesh = occHandle(Ng_MeshVS_DataSource3D)::DownCast(m);
        MeshTools::buildPointToElementConnectivity(volumeMesh,nodeToVolumeElementMap);

        for(std::map<int,std::vector<double>>::iterator it_ = updatedPoints.begin(); it_!=updatedPoints.end(); it_++)
        {
            int globalNodeID_ = it_->first;
            const std::vector<double> &P_updated = it_->second;

            //! -------------------------------------------------------
            //! store the old point coordinates before moving the node
            //! -------------------------------------------------------
            double P[3];
            this->getPointCoordinate(volumeMesh,globalNodeID_,P);
            std::vector<double> P_old { P[0], P[1], P[2] };

            //! ------------------
            //! displace the node
            //! ------------------
            int localNodeID_ = volumeMesh->myNodesMap.FindIndex(globalNodeID_);
            std::vector<double> Pnew { P_updated[0], P_updated[1], P_updated[2] };
            volumeMesh->changeNodeCoords(localNodeID_,Pnew);

            //! ---------------------------------------------------
            //! check if the displacement does not invert elements
            //! ---------------------------------------------------
            bool isShiftOK = true;
            const std::vector<int> &attachedElements = nodeToVolumeElementMap.find(globalNodeID_)->second;
            for(std::vector<int>::const_iterator ite = attachedElements.cbegin(); ite != attachedElements.cend(); ite++)
            {
                int globalElementID_attached = *ite;
                int NbNodesOfElement;
                MeshVS_EntityType aType;
                double bufd[24];
                TColStd_Array1OfReal coordsElement(*bufd,1,24);
                volumeMesh->GetGeom(globalElementID_attached,true,coordsElement,NbNodesOfElement,aType);
                TetQualityClass aTet(coordsElement(1),coordsElement(2),coordsElement(3),
                                     coordsElement(4),coordsElement(5),coordsElement(6),
                                     coordsElement(7),coordsElement(8),coordsElement(9),
                                     coordsElement(10),coordsElement(11),coordsElement(12));
                if(aTet.Volume()<=0)
                {
                    cerr<<"____shifting the node: "<<globalNodeID_<<" has made the volume element: "<<globalElementID_attached<<" not valid____"<<endl;
                    isShiftOK = false;
                    break;
                }
            }
            if(isShiftOK==false)
            {
                std::map<int,std::vector<double>>::iterator pos = updatedPoints.find(globalNodeID_);
                it_ = updatedPoints.erase(pos);
                volumeMesh->changeNodeCoords(localNodeID_,P_old);
                int localNodeID_surface = surfaceMesh->myNodesMap.FindIndex(globalNodeID_);
                surfaceMesh->changeNodeCoords(localNodeID_surface,P_old);
            }
            else it_++;
        }
        mDB->ArrayOfMeshDS.insert(bodyIndex,volumeMesh);
    }
    mDB->ArrayOfMeshDS2D.insert(bodyIndex, surfaceMesh);

    //! ----------------------------------
    //! update the face mesh data sources
    //! ----------------------------------
    int NbFaces = mDB->MapOfBodyTopologyMap.value(bodyIndex).faceMap.Extent();
    for(int faceNr = 1; faceNr<=NbFaces; faceNr++)
    {
        occHandle(MeshVS_DataSource) m = mDB->ArrayOfMeshDSOnFaces.getValue(bodyIndex,faceNr);
        if(m.IsNull()) continue;
        occHandle(Ng_MeshVS_DataSourceFace) aFaceMeshDS = occHandle(Ng_MeshVS_DataSourceFace)::DownCast(m);
        for(std::map<int,std::vector<double>>::iterator it_ = updatedPoints.begin(); it_!=updatedPoints.end(); it_++)
        {
            int globalNodeID_toShift = it_->first;
            int localNodeID_toShift = aFaceMeshDS->myNodesMap.FindIndex(globalNodeID_toShift);
            const std::vector<double> &Pnew = it_->second;
            aFaceMeshDS->changeNodeCoords(localNodeID_toShift,Pnew);
        }
    }
    return true;
}

//! ---------------------------------
//! function: buildShapeTessellation
//! details:
//! ---------------------------------
bool surfaceMeshSmoother::buildShapeTessellation(const TopoDS_Shape &aShape,
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

//! ------------------------
//! function: initRayCaster
//! details:
//! ------------------------
bool surfaceMeshSmoother::initRayCaster(const occHandle(Ng_MeshVS_DataSource2D) &aSurfaceMeshDS)
{
    cout<<"surfaceMeshToFaceMeshes::initRayCaster()->____function called____"<<endl;

    //! -------------
    //! sanity check
    //! -------------
    if(aSurfaceMeshDS.IsNull()) return false;
    if(aSurfaceMeshDS->GetAllElements().Extent()<1) return false;
    if(aSurfaceMeshDS->GetAllNodes().Extent()<3) return false;

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

    //! ------------
    //! init embree
    //! ------------
    myEmbree.init(V.cast<float>(),F.cast<int>());

    cout<<"surfaceMeshToFaceMeshes::initRayCaster()->____raycaster initialized____"<<endl;
    return true;
}

//! ----------------------
//! function: perturbMesh
//! details:
//! ----------------------
void surfaceMeshSmoother::perturbMesh(occHandle(Ng_MeshVS_DataSource2D) &surfaceMesh)
{
    cout<<"surfaceMeshSmoother::perturbMesh()->____function called____"<<endl;
    std::random_device rd;
    std::mt19937::result_type seed = rd() ^ (
                (std::mt19937::result_type)
                std::chrono::duration_cast<std::chrono::seconds>(
                    std::chrono::system_clock::now().time_since_epoch()
                    ).count() +
                (std::mt19937::result_type)
                std::chrono::duration_cast<std::chrono::microseconds>(
                    std::chrono::high_resolution_clock::now().time_since_epoch()
                    ).count() );

    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> distrib(-1e-5, 1e-5);

    for(TColStd_MapIteratorOfPackedMapOfInteger it(surfaceMesh->GetAllNodes()); it.More(); it.Next())
    {
        int globalNodeID = it.Key();
        int localNodeID = surfaceMesh->myNodesMap.FindIndex(globalNodeID);
        double P[3];
        this->getPointCoordinate(surfaceMesh,globalNodeID,P);
        double x = P[0]+distrib(gen);
        double y = P[1]+distrib(gen);
        double z = P[2]+distrib(gen);
        std::vector<double> perturbedPoint { x, y, z };
        surfaceMesh->changeNodeCoords(localNodeID,perturbedPoint);
    }
    cout<<"surfaceMeshSmoother::perturbMesh()->____exiting function: mesh point perturbed____"<<endl;
}

//! ------------------------------
//! function: getPointCoordinates
//! details:  helper
//! ------------------------------
bool surfaceMeshSmoother::getPointCoordinate(const occHandle(MeshVS_DataSource) &aMeshDS, int globalNodeID, double *P)
{
    MeshVS_EntityType aType;
    double buf[3];
    TColStd_Array1OfReal coords(*buf,1,3);
    int NbNodes;
    bool isDone = aMeshDS->GetGeom(globalNodeID,false,coords,NbNodes,aType);
    for(int n=0; n<3; n++) P[0] = coords(n+1);
    return isDone;
}
