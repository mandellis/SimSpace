//! ----------------
//! custom includes
//! ----------------
#include "tethex.h"
#include "tetsplitter.h"
#include <meshelementbycoords.h>
#include "global.h"
#include "qprogressevent.h"

//! ---
//! Qt
//! ---
#include <QApplication>
#include <QThread>

//! ----
//! OCC
//! ----
#include <TColStd_PackedMapOfInteger.hxx>
#include <TColStd_MapIteratorOfPackedMapOfInteger.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <MeshVS_EntityType.hxx>

//! ----
//! C++
//! ----
#include <vector>
#include <map>
#include <iostream>
using namespace std;

//! ------------------
//! function: perform
//! details:
//! ------------------
userMessage TetHex::perform(int bodyIndex)
{
    cout<<"TetHex::perform()->____function called____"<<endl;

    //! -------------------------------------------------
    //! retrieve the initial tet mesh from the data base
    //! -------------------------------------------------
    occHandle(Ng_MeshVS_DataSource3D) aTetMesh = occHandle(Ng_MeshVS_DataSource3D)::DownCast(myMDB->ArrayOfMeshDS.value(bodyIndex));

    //! -------------------
    //! define the results
    //! -------------------
    occHandle(Ng_MeshVS_DataSource3D) HexaMeshDS;
    occHandle(Ng_MeshVS_DataSource2D) QuadMeshDS;

    //! ------------------------------------------
    //! generate the hexahedral mesh
    //! insert the data source into the data base
    //! ------------------------------------------
    cout<<"____tag00____"<<endl;
    userMessage um = makeHexa(aTetMesh,HexaMeshDS);

    if(!um.isDone)
    {
        um.isDone = false;
        um.message = QString("TetHex: error in generating the volume mesh data source");
        return um;
    }
    //myMDB->ArrayOfMeshDS.insert(bodyIndex,HexaMeshDS);

    //! ------------------------------------------------------------
    //! generate the surface mesh data source
    //! insert the quad surface mesh data source into the data base
    //! ------------------------------------------------------------
    if(HexaMeshDS->myFaceToElements.isEmpty()) HexaMeshDS->buildFaceToElementConnectivity();
    QuadMeshDS = new Ng_MeshVS_DataSource2D(HexaMeshDS);

    //! -------------
    //! sanity check
    //! -------------
    if(QuadMeshDS.IsNull())
    {
        um.isDone = false;
        um.message = QString("TetHex: error in generating the surface mesh data source");
        Global::status().myMessages->appendMessage(um);
        return um;
    }

    //! --------------------------------------------
    //! generate the map for the renumbering
    //! hash(x,y,z) [size_t] <=> globalNodeID [int]
    //! --------------------------------------------
    std::map<size_t,int> mapNodeCoords_globalNodeID;
    std::vector<int> globalNodeIDs;
    std::vector<mesh::tolerantPoint> vecTolerantPoints;

    for(TColStd_MapIteratorOfPackedMapOfInteger it(QuadMeshDS->GetAllNodes()); it.More(); it.Next())
    {
        int globalNodeID = it.Key();
        int localNodeID = QuadMeshDS->myNodesMap.FindIndex(globalNodeID);
        const std::vector<double> &aPoint = QuadMeshDS->getNodeCoordinates(localNodeID);

        //! ---------------------------------------------------------
        //! vector of the points (x,y,z) and vector of their indices
        //! ---------------------------------------------------------
        vecTolerantPoints.push_back(mesh::tolerantPoint(aPoint[0],aPoint[1],aPoint[2]));
        globalNodeIDs.push_back(globalNodeID);

        size_t seed = 0;
        for(int j=0; j<3; j++) hash_c<double>(seed, aPoint[j]);
        std::pair<size_t,int> aPair;
        aPair.first = seed;
        aPair.second = globalNodeID;
        mapNodeCoords_globalNodeID.insert(aPair);
    }

    //! ----------------------
    //! post a progress event
    //! ----------------------
    int NbFaces = myMDB->MapOfBodyTopologyMap.value(bodyIndex).faceMap.Extent();
    QProgressEvent *aProgressEvent = new QProgressEvent(QProgressEvent_None,0,999,0,"Building face mesh data sources",
                                                       QProgressEvent_Init,1,NbFaces,0,"TetHex conversion");
    if(myProgressIndicator!=Q_NULLPTR)
    {
        QApplication::postEvent(myProgressIndicator,aProgressEvent);
        QApplication::processEvents();
    }

    //! ------------------------------------
    //! generate the face mesh data sources
    //! ------------------------------------
    for(int faceNr = 1; faceNr<=NbFaces; faceNr++)
    {
        //! ----------------------
        //! post a progress event
        //! ----------------------
        aProgressEvent = new QProgressEvent(QProgressEvent_None,0,999,0,"Building face mesh data sources",
                                            QProgressEvent_Init,0,NbFaces-1,0,"TetHex conversion");
        if(myProgressIndicator!=Q_NULLPTR)
        {
            QApplication::postEvent(myProgressIndicator,aProgressEvent);
            QApplication::processEvents();
        }

        //! ------------------
        //! check termination
        //! ------------------
        if(Global::status().code==0)
        {
            um.isDone = false;
            um.message = "TextHex: process interrupted by the user while generating the face mesh data sources";
            Global::status().myMessages->appendMessage(um);
            return um;
        }

        //! --------------------------------------------------------------
        //! retrieve the current face mesh data source from the data base
        //! This is composed by surface elements that are not quads
        //! --------------------------------------------------------------
        occHandle(Ng_MeshVS_DataSourceFace) curFaceMeshDS = occHandle(Ng_MeshVS_DataSourceFace)::DownCast(myMDB->ArrayOfMeshDSOnFaces.getValue(bodyIndex,faceNr));
        if(curFaceMeshDS.IsNull()) continue;

        //! ---------------------------------------------------------------
        //! iterate over the triangles and split each one into three quads
        //! ---------------------------------------------------------------
        std::vector<meshElementByCoords> quadDecomposition;
        for(TColStd_MapIteratorOfPackedMapOfInteger it = curFaceMeshDS->GetAllElements(); it.More(); it.Next())
        {
            int globalElementID = it.Key();

            double buf[12];
            TColStd_Array1OfReal coords(*buf,1,12);
            MeshVS_EntityType type;
            int NbNodes;
            if(curFaceMeshDS->GetGeom(globalElementID,true,coords,NbNodes,type))
            {
                meshElementByCoords aMeshElement;
                for(int n=0; n<NbNodes; n++)
                {
                    int s = 3*n;
                    double x = coords(s+1);
                    double y = coords(s+2);
                    double z = coords(s+3);
                    mesh::meshPoint aP(x,y,z);
                    aMeshElement.pointList<<aP;
                    aMeshElement.type = QUAD;
                }
                std::vector<meshElementByCoords> threeQuads = tetSplitter::splitTri(aMeshElement);
                for(int k=0; k<3; k++) quadDecomposition.push_back(threeQuads[k]);
            }
        }

        //! -----------------------------------------
        //! call the face mesh data base constructor
        //! -----------------------------------------
        occHandle(Ng_MeshVS_DataSourceFace) quadFaceMeshDS = new Ng_MeshVS_DataSourceFace(quadDecomposition);
        if(!quadFaceMeshDS.IsNull())
        {
            //! -------------
            //! disable stop
            //! -------------
            if(myProgressIndicator!=Q_NULLPTR)
            {
                //myProgressIndicator->disableStop();
                this->setStopButtonEnabled(false);
            }

            //! -------------------------
            //! renumber face mesh nodes
            //! -------------------------
            //mapNodeCoords_globalNodeID.clear();
            quadFaceMeshDS->renumberNodes(mapNodeCoords_globalNodeID,vecTolerantPoints,globalNodeIDs);

            //! -------------------------------------------------
            //! put the face mesh data sources into the database
            //! -------------------------------------------------
            myMDB->ArrayOfMeshDSOnFaces.setValue(bodyIndex,faceNr,quadFaceMeshDS);

            //! ------------
            //! enable stop
            //! ------------
            if(myProgressIndicator!=Q_NULLPTR)
            {
                this->setStopButtonEnabled(true);
                //myProgressIndicator->enableStop();
            }
        }    
    }

    //! ----------------------
    //! post a progress event
    //! ----------------------
    if(myProgressIndicator!=Q_NULLPTR)
    {
        aProgressEvent = new QProgressEvent(QProgressEvent_None,0,100,0,"TetHex conversion done",QProgressEvent_Reset,0,0,0,myTask);
        QApplication::postEvent(myProgressIndicator,aProgressEvent);
        QApplication::processEvents();
    }

    //! --------------------------------------
    //! insert the results into the data base
    //! --------------------------------------
    myMDB->ArrayOfMeshDS.insert(bodyIndex,HexaMeshDS);
    myMDB->ArrayOfMeshDS2D.insert(bodyIndex,QuadMeshDS);

    um.isDone = true;
    um.message = QString("TetHex: convertion done on body ")+myMDB->MapOfBodyNames.value(bodyIndex);

    cout<<"TetHex::perform()->____exiting function____"<<endl;

    Global::status().myMessages->appendMessage(um);
    return um;
}

//! -------------------
//! function: makeHexa
//! details:
//! -------------------
userMessage TetHex::makeHexa(const occHandle(Ng_MeshVS_DataSource3D) &aTetMesh, occHandle(Ng_MeshVS_DataSource3D) &HexaMeshDS)
{
    //! --------------------
    //! init a user message
    //! --------------------
    userMessage um(false,"");

    //! -------------
    //! sanity check
    //! -------------
    if(aTetMesh.IsNull())
    {
        um.message = "TetHex: the input mesh is null";
        return um;
    }
    if(aTetMesh->GetAllElements().Extent()<1 || aTetMesh->GetAllNodes().Extent()<4)
    {
        um.message = "TetHex: the input mesh is not valid";
        return um;
    }

    //! -----------------------
    //! the number of elements
    //! -----------------------
    int NbSteps = aTetMesh->GetAllElements().Extent();
    cout<<"userMessage::perform()->____number of elements: "<<NbSteps<<"____"<<endl;

    //! ----------------------
    //! send a progress event
    //! ----------------------
    if(myProgressIndicator!=Q_NULLPTR)
    {
        QProgressEvent *pe = new QProgressEvent(QProgressEvent_None,0,100,0,"Starting TetHex",QProgressEvent_Init,0,NbSteps-1,0,myTask);
        QApplication::postEvent(myProgressIndicator,pe);
        QApplication::processEvents();
    }

    //! -----------------------------------------
    //! all the hexa element after tet splitting
    //! -----------------------------------------
    std::vector<meshElementByCoords> allHexaElements;

    //! -------------------------------
    //! iterate over the mesh elements
    //! -------------------------------
    int step = 0;   // the progress
    int NbNodes;
    double buf[12];
    TColStd_Array1OfReal coords(*buf,1,12);
    MeshVS_EntityType aType;
    for(TColStd_MapIteratorOfPackedMapOfInteger it(aTetMesh->GetAllElements()); it.More(); it.Next())
    {
        //! ------------------
        //! check termination
        //! ------------------
        if(Global::status().code ==0)
        {
            um.isDone = false;
            um.message = "TetHex: conversion stopped by the user";
            return um;
        }

        int globalElementID = it.Key();
        aTetMesh->GetGeom(globalElementID,true,coords,NbNodes,aType);

        QList<mesh::meshPoint> pointList;
        for(int i=0; i<NbNodes; i++)
        {
            int s = 3*i;
            double x = coords(s+1);
            double y = coords(s+2);
            double z = coords(s+3);

            mesh::meshPoint aP(x,y,z);
            pointList<<aP;
        }

        meshElementByCoords curMeshElement;
        curMeshElement.pointList<<pointList;
        std::vector<meshElementByCoords> hexaElements = tetSplitter::splitTet(curMeshElement);
        for(int i=0; i<4; i++) allHexaElements.push_back(hexaElements[i]);

        //! --------------------
        //! send progress event
        //! --------------------
        if(myProgressIndicator!=Q_NULLPTR)
        {
            if(NbSteps%25==0)
            {
                QString msg = QString("Converting tetrahedra to hexahedra");
                QProgressEvent *pe = new QProgressEvent(QProgressEvent_None,0,100,0,msg,QProgressEvent_Update,0,NbSteps-1,step,myTask);
                QApplication::postEvent(myProgressIndicator,pe);
                QApplication::processEvents();
                step++;
            }
        }
    }

    if(allHexaElements.size()<4)
    {
        um.isDone = false;
        um.message = QString("Invalid mesh");
        return um;
    }

    //! --------------------------
    //! call the mesh constructor
    //! --------------------------
    HexaMeshDS = new Ng_MeshVS_DataSource3D(allHexaElements);

    /*
    //! ----------------------
    //! post a progress event
    //! ----------------------
    if(myProgressIndicator!=Q_NULLPTR)
    {
        QProgressEvent *pe = new QProgressEvent(QProgressEvent_None,0,100,0,"TetHex conversion done",QProgressEvent_Reset,0,0,0,task);
        QApplication::postEvent(myProgressIndicator,pe);
        QApplication::processEvents();
    }
    */

    //! -------------------------------------
    //! post a message into the tabular view
    //! -------------------------------------
    um.isDone = true;
    um.message = QString("TetHex conversion done");
    Global::status().myMessages->appendMessage(um);
    return um;
}

//! ---------------------------------------------------------------
//! function: setStopButtonEnabled
//! details:  enable/disable the progress indicator stop button
//!           here an event is posted to the application main loop
//! ---------------------------------------------------------------
void TetHex::setStopButtonEnabled(bool isEnabled)
{
    QProgressEvent *pe;
    if(isEnabled) pe = new QProgressEvent(QProgressEvent_EnableStop,0,0,0,"",QProgressEvent_None,0,0,0,myTask);
    else pe = new QProgressEvent(QProgressEvent_DisableStop,0,0,0,"",QProgressEvent_None,0,0,0,myTask);
    QApplication::postEvent(myProgressIndicator,pe);
}
