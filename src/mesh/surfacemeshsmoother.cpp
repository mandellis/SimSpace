//! ----------------
//! custom includes
//! ----------------
#include "surfacemeshsmoother.h"
#include "qprogressevent.h"
#include "global.h"
#include "userMessage.h"
#include <mesh.h>
#include <nettgentool2.h>

//! ----
//! OCC
//! ----
#include <TColStd_MapIteratorOfPackedMapOfInteger.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <MeshVS_EntityType.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <gp_Pnt.hxx>
#include <TopExp_Explorer.hxx>
#include <TopAbs_ShapeEnum.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <GeomAdaptor_Surface.hxx>
#include <Geom_Surface.hxx>
#include <TopAbs_State.hxx>
#include <BRepClass_FaceClassifier.hxx>
#include <BRep_Tool.hxx>
#include <Extrema_ExtAlgo.hxx>
#include <TopExp.hxx>
#include <TopTools_IndexedMapOfShape.hxx>

//! ----
//! C++
//! ----
#include <vector>
#include <chrono>

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
#include "OCCface.h"

//! ------
//! Eigen
//! ------
#include <Eigen/Dense>

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
bool surfaceMeshSmoother::perform(meshDataBase *mDB, int bodyIndex, bool rebuildVolumeMesh)
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

    occHandle(MeshVS_DataSource) aMesh = mDB->ArrayOfMeshDS2D.value(bodyIndex);
    if(aMesh.IsNull()) return false;
    occHandle(Ng_MeshVS_DataSource2D) surfaceMesh = occHandle(Ng_MeshVS_DataSource2D)::DownCast(aMesh);

    //! --------------------------------------------------------
    //! face map of the shape: user for constructing an OCCFace
    //! --------------------------------------------------------
    const TopoDS_Shape &curBody = mDB->bodyMap.value(bodyIndex);
    TopTools_IndexedMapOfShape faceMap;
    TopExp::MapShapes(curBody,TopAbs_FACE,faceMap);
    std::vector<TopoDS_Face> vecFaces;
    int NbGeometryFaces = faceMap.Extent();
    for(int faceNr = 1; faceNr<=NbGeometryFaces; faceNr++)
    {
        TopoDS_Face curFace = TopoDS::Face(faceMap.FindKey(faceNr));
        vecFaces.push_back(curFace);
        if(curFace.IsNull())
        {
            cout<<"surfaceMeshSmoother::perform()->____the face nr: "<<faceNr<<" is null: jumping over it____"<<endl;
            continue;
        }
    }

    //! -------------------------
    //! this supports projection
    //! -------------------------
    OCCFace SurfaceBody(vecFaces);

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
        e = new QProgressEvent(QProgressEvent_None,-1,-1,-1,"Smoothing mesh",QProgressEvent_Init,0,NbPoints,0,task);
        QApplication::postEvent(myProgressIndicator,e);
        QApplication::processEvents();
        QThread::msleep(1000);
    }

    //! ------------
    //! enable stop
    //! ------------
    emit requestEnableStop();

    int count = 0;
    int projected = 0;
    std::map<int,mesh::meshPoint> mapOfProjectedPoints;
    for(TColStd_MapIteratorOfPackedMapOfInteger it(surfaceMesh->GetAllNodes()); it.More(); it.Next())
    {
        count++;
        if(Global::status().code==0)
        {
            if(myProgressIndicator!=Q_NULLPTR)
            {
                e = new QProgressEvent(QProgressEvent_Reset,0,100,0,"Process interrupted by the user",
                                       QProgressEvent_Reset,0,100,0,"",Q_NULLPTR);
                QApplication::postEvent(myProgressIndicator,e);
                QApplication::processEvents();
            }
            Global::status().code = 1;  //?
            return false;
        }

        int globalNodeID = it.Key();
        int localNodeID = surfaceMesh->myNodesMap.FindIndex(globalNodeID);
        std::vector<double> aP = surfaceMesh->getNodeCoordinates(localNodeID);
        double theProjectedPoint[3];
        bool isDone = SurfaceBody.pointProjection(aP.data(), theProjectedPoint);
        if(isDone == true)
        {
            projected++;

            //! --------------------------------------------
            //! store the projected mesh point into the map
            //! --------------------------------------------
            std::pair<int,mesh::meshPoint> mapElement;
            mapElement.first = globalNodeID;
            mapElement.second = mesh::meshPoint(theProjectedPoint[0],theProjectedPoint[1],theProjectedPoint[2],globalNodeID);
            mapOfProjectedPoints.insert(mapElement);

            std::vector<double> newCoords { theProjectedPoint[0], theProjectedPoint[1], theProjectedPoint[2]};

            //! ---------------------------------
            //! change the position of the point
            //! ---------------------------------
            surfaceMesh->changeNodeCoords(localNodeID,newCoords);
        }
        //! --------------------
        //! update the progress
        //! --------------------
        if(myProgressIndicator!=Q_NULLPTR)
        {
            if(count%100 == 0)
            {
                QString msg = QString("Smoothing the surface mesh");
                e = new QProgressEvent(QProgressEvent_None,-1,-1,-1,msg,QProgressEvent_Update,-1,-1,count,task);
                QApplication::postEvent(myProgressIndicator,e);
                QApplication::processEvents();
            }
        }
    }

    emit requestDisableStop();

    //! -----------------------------------------------
    //! recompute element and nodal normals and insert
    //! the new surface mesh into the database
    //! -----------------------------------------------
    surfaceMesh->computeNormalAtElements();
    surfaceMesh->computeNormalAtNodes();
    mDB->ArrayOfMeshDS2D.insert(bodyIndex,surfaceMesh);

    //! ---------------------------------------
    //! update also the face mesh data sources
    //! ---------------------------------------
    if(myProgressIndicator!=Q_NULLPTR)
    {
        e = new QProgressEvent(QProgressEvent_None,-1,-1,-1,"Smoothing mesh",QProgressEvent_Init,1,NbGeometryFaces,1,task);
        QApplication::postEvent(myProgressIndicator,e);
        QApplication::processEvents();
    }

    for(int faceNr = 1; faceNr<= NbGeometryFaces; faceNr++)
    {
        occHandle(MeshVS_DataSource) aDS = mDB->ArrayOfMeshDSOnFaces.getValue(bodyIndex,faceNr);
        if(aDS.IsNull()) continue;
        occHandle(Ng_MeshVS_DataSourceFace) curFaceMeshDS = occHandle(Ng_MeshVS_DataSourceFace)::DownCast(aDS);
        for(TColStd_MapIteratorOfPackedMapOfInteger it(curFaceMeshDS->GetAllNodes()); it.More(); it.Next())
        {
            int globalNodeID = it.Key();
            std::map<int,mesh::meshPoint>::iterator it1 = mapOfProjectedPoints.find(globalNodeID);
            if(it1==mapOfProjectedPoints.end()) continue;
            mesh::meshPoint aP = (*it1).second;
            int localNodeID = curFaceMeshDS->myNodesMap.FindIndex(globalNodeID);
            std::vector<double> P{aP.x,aP.y,aP.z};
            curFaceMeshDS->changeNodeCoords(localNodeID,P);
        }

        //! ------------------------------------------
        //! recompute normal at elements and at nodes
        //! insert the face mesh DS into the database
        //! ------------------------------------------
        curFaceMeshDS->computeNormalAtElements();
        curFaceMeshDS->computeNormalAtNodes();
        mDB->ArrayOfMeshDSOnFaces.setValue(bodyIndex,faceNr,curFaceMeshDS);

        if(myProgressIndicator!=Q_NULLPTR)
        {
            QString msg = QString("Smoothing surface mesh");
            e = new QProgressEvent(QProgressEvent_None,-1,-1,-1,msg,QProgressEvent_Update,-1,-1,faceNr,task);
            QApplication::postEvent(myProgressIndicator,e);
            QApplication::processEvents();
        }
    }

    cout<<"@-----------------------------------"<<endl;
    cout<<"@- projected points: "<<100.0*double(projected)/double(surfaceMesh->GetAllNodes().Extent())<<" %"<<endl;
    cout<<"@-----------------------------------"<<endl;

    //! ----------------------------------------------
    //! check if there is a volume mesh to be adapted
    //! (can use also the flag "toBeUpdated")
    //! ----------------------------------------------
    aMesh = mDB->ArrayOfMeshDS.value(bodyIndex);
    if(!aMesh.IsNull())
    {
        if(rebuildVolumeMesh==true)
        {
            occHandle(Ng_MeshVS_DataSource3D) curVolumeMesh = occHandle(Ng_MeshVS_DataSource3D)::DownCast(aMesh);
            for(TColStd_MapIteratorOfPackedMapOfInteger it(surfaceMesh->GetAllNodes()); it.More(); it.Next())
            {
                int globalNodeID = it.Key();
                std::map<int,mesh::meshPoint>::iterator it1 = mapOfProjectedPoints.find(globalNodeID);
                if(it1==mapOfProjectedPoints.end()) continue;
                mesh::meshPoint aP = (*it1).second;
                std::vector<double> P {aP.x,aP.y,aP.z};
                int localNodeID = curVolumeMesh->myNodesMap.FindKey(globalNodeID);
                curVolumeMesh->changeNodeCoords(localNodeID,P);
            }
            mDB->ArrayOfMeshDS.insert(bodyIndex,curVolumeMesh);
        }
        else
        {
            NettgenTool2 netgenMesher;
            netgenMesher.setMeshDataBase(mDB);
            netgenMesher.setProgressIndicator(myProgressIndicator);
            occHandle(Ng_MeshVS_DataSource3D) volMesh;
            um = netgenMesher.meshVolumeFromSurfaceMesh(bodyIndex,surfaceMesh,volMesh);
            if(um.isDone==true) mDB->ArrayOfMeshDS.insert(bodyIndex,volMesh);
            else
            {
                //! -----------------------------------------------------
                //! reset all the mesh data sources for the current body
                //! -----------------------------------------------------
                mDB->ArrayOfMeshDS.insert(bodyIndex,occHandle(Ng_MeshVS_DataSource3D)());
                mDB->ArrayOfMeshDS2D.insert(bodyIndex,occHandle(Ng_MeshVS_DataSource2D)());
                for(int faceNr = 1; faceNr<= NbGeometryFaces; faceNr++)
                    mDB->ArrayOfMeshDSOnFaces.setValue(bodyIndex,faceNr,occHandle(Ng_MeshVS_DataSourceFace)());
                mDB->ArrayOfMeshIsToBeUdpdated.insert(bodyIndex,true);
            }
        }
    }
    else
    {
        cout<<"@ --------------------------"<<endl;
        cout<<"@ - there is no volume mesh "<<endl;
        cout<<"@ --------------------------"<<endl;
        return false;
    }
    return true;
}
