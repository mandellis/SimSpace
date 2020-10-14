//#define OLD_VERSION
#define MESHING_STATUS_MESSAGE_PERSISTANCE 500
#define USE_MESHFIX

//! ----------------
//! custom includes
//! ----------------
#include <mesherclass.h>
#include <ng_mesher2.h>
#include <occmesher.h>
#include <meshtools.h>
#include <tetgenmesher.h>
#include <edgemeshdatasourcebuildesclass.h>
#include <ng_meshvs_datasource3d.h>
#include "testtools.h"
#include "qprogressevent.h"
#include "qprogressindicator.h"
#include "tools.h"
#include "ccout.h"
#include "extendedrwstl.h"
#include "stldoctor.h"
#include "stlapiwriter.h"
#include <prismaticlayer.h>
#include <prismaticlayerparameters.h>
#include "meshtools.h"
#include "geomtoolsclass.h"
#include <tetwildmesher.h>
#include <surfacemeshsmoother.h>
//#include <renumberingtool.h>
#include <meshnodesrenumberingtool.h>

//! -------------
//! experimental
//! -------------
#include <meshselectnlayers.h>

//! ----
//! OCC
//! ----
#include <TopoDS_Solid.hxx>
#include <TopExp.hxx>
#include <BRepTools.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <StlMesh_Mesh.hxx>
#include <StlMesh.hxx>
#include <RWStl.hxx>
#include <NCollection_Vector.hxx>
#include <TColStd_MapIteratorOfMapOfInteger.hxx>

//! ---
//! Qt
//! ---
#include <QApplication>

//! ---------------------
//! function: destructor
//! details:
//! ---------------------
MesherClass::~MesherClass()
{
    cout<<"MesherClass::~MesherClass()->____DESTRUCTOR CALLED____"<<endl;
    //! -------------------------
    //! remove the support files
    //! -------------------------
    this->removeSupportFiles(true);
}

//! ------------------------------
//! function: constructor
//! details:  default constructor
//! ------------------------------
MesherClass::MesherClass(QObject *parent):QObject(parent)
{
    //! -------------------------
    //! remove the support files
    //! -------------------------
    this->removeSupportFiles(true);    
}

//! -----------------------------------------------------------
//! function: init
//! details:  used with the constructor without initialization
//! -----------------------------------------------------------
void MesherClass::init(meshDataBase *mDB, const QList<int> &listOfBodies, bool isVolume)
{
    myBodyList = listOfBodies;
    myMeshDB = mDB;
    myIsVolume = isVolume;

    //! --------------------------------------------------
    //! the Netgen mesher is incorporated as class member
    //! the Tetgen mesher is incorporated as class member
    //! the Express mesh is incorporated as class member ... to do ...
    //! --------------------------------------------------
    myNetgenMesher = new NettgenTool2(this);
    myTetMesher = new TetgenMesher(myMeshDB,this);

    //! -------------------------
    //! remove the support files
    //! -------------------------
    this->removeSupportFiles(true);
}

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
MesherClass::MesherClass(meshDataBase *mDB, const QList<int> &listOfBodies, bool isVolume, QObject *parent):
    QObject(parent)
{
    myBodyList = listOfBodies;
    myMeshDB = mDB;
    myIsVolume = isVolume;

    //! ------------------------------------------------------
    //! the Netgen mesher is incorporated as class member
    //! the Netgen stl mesher is incorporated as class member
    //! the Tetgen mesher is incorporated as class member
    //! ------------------------------------------------------
    myNetgenMesher = new NettgenTool2(this);
    myTetMesher = new TetgenMesher(myMeshDB,this);

    //! -------------------------------------------------------------
    //! new design: the Express mesh is incorporated as class member
    //! -------------------------------------------------------------
    //! ...

    //! -------------------------
    //! remove the support files
    //! -------------------------
    this->removeSupportFiles(true);
}

//! ----------------------------------------------------------------
//! function: setProgressIndicator
//! details:  set the progress indicator for all the mesher engines
//! ----------------------------------------------------------------
void MesherClass::setProgressIndicator(QProgressIndicator *aProgressIndicator)
{
    if(aProgressIndicator!= Q_NULLPTR)
    {
        myProgressIndicator = aProgressIndicator;
        disconnect(this,SIGNAL(requestDisableStop()),myProgressIndicator,SLOT(disableStop()));
        disconnect(this,SIGNAL(requestDisableStop()),myProgressIndicator,SLOT(enableStop()));
        connect(this,SIGNAL(requestDisableStop()),myProgressIndicator,SLOT(disableStop()));
        connect(this,SIGNAL(requestEnableStop()),myProgressIndicator,SLOT(enableStop()));

        myTetMesher->setProgressIndicator(myProgressIndicator);
        myNetgenMesher->setProgressIndicator(myProgressIndicator);
    }

    //! ----------------
    //! netgen specific
    //! ----------------
    disconnect(myNetgenMesher,SIGNAL(requestStartingNetgenEnquireTimer()),this,SLOT(handleRequestStartingNetgenEnquireTimer()));
    disconnect(myNetgenMesher,SIGNAL(requestStoppingNetgenEnquireTimer()),this,SLOT(handleRequestStoppingNetgenEnquireTimer()));

    connect(myNetgenMesher,SIGNAL(requestStartingNetgenEnquireTimer()),this,SLOT(handleRequestStartingNetgenEnquireTimer()));
    connect(myNetgenMesher,SIGNAL(requestStoppingNetgenEnquireTimer()),this,SLOT(handleRequestStoppingNetgenEnquireTimer()));
}

//! --------------------------
//! function: postUpdateEvent
//! details:
//! --------------------------
void MesherClass::postUpdateEvent(QProgressEvent *aProgressEvent)
{
    if(myProgressIndicator!=Q_NULLPTR)
    {
        cout<<"MesherClass::postUpdateEvent()->____posting an update event____"<<endl;
        QApplication::postEvent(myProgressIndicator, aProgressEvent);
        QApplication::processEvents();
    }
}

//! ----------------------------------
//! function: generateMesh
//! details:  iterate over the bodies
//! ----------------------------------
#define EXPERIMENTAL_BLOCK
void MesherClass::generateMesh()
{
    cout<<"MesherClass::generateMesh()->____function called____"<<endl;

    //! ---------------------------------------------------------
    //! number of meshes
    //! progress bar: one meshing step, 2D or 3D, for each body
    //! ---------------------------------------------------------
    int NbSteps = myBodyList.length();

    //if(NbSteps==1) myProgressIndicator->setPrimaryBarVisible(false);

    //! --------------------------------------------
    //! send an "Init" progress to the progress bar
    //! --------------------------------------------
    int done=0;
    QProgressEvent *pe = new QProgressEvent(QProgressEvent_Init,0,NbSteps-1,0,"",QProgressEvent_Init,0,100,0,0);
    postUpdateEvent(pe);

    for(int i=0; i<myBodyList.length(); i++)
    {
        //! -----------------------------------------------
        //! at every iteration check if "Stop" was pressed
        //! -----------------------------------------------
        QApplication::processEvents();

        int code = Global::status().code;
        if(code==1)
        {
            int bodyIndex = myBodyList.at(i);

            //! ---------------------------------
            //! a part of the meshing parameters
            //! ---------------------------------
            meshParam mp = myMeshDB->getMeshingParameters(bodyIndex);

            //! -------------------------------------------
            //! define the surface and volume data sources
            //! -------------------------------------------
            occHandle(MeshVS_DataSource) mainMesh2D;
            occHandle(Ng_MeshVS_DataSource3D) mainMesh3D;

            //! ----------------------------------------------------------------
            //! this paramter determines which surface meshing strategy is used
            //! ----------------------------------------------------------------
            bool useBRep = myMeshDB->mapOfUseBRep.value(bodyIndex);
            bool isMeshingInMemory = myMeshDB->MapOfIsMeshingRunningInMemory.value(bodyIndex);
            bool isSecondOrder = myMeshDB->ArrayOfMeshOrder.value(bodyIndex);

            //! --------------------------------------------------------------
            //! generate the surface mesh. In case of prismatic mesh generate
            //! the prismatic 3D elements
            //! --------------------------------------------------------------
            if(myIsVolume==false)
            {
                cout<<"@____start generating the surface mesh____"<<endl;

                //! --------------------
                //! update the progress
                //! --------------------
                QProgressEvent *pe = new QProgressEvent();
                pe->setVal(done);
                pe->setMessage(QString("Start meshing body %1").arg(myMeshDB->MapOfBodyNames.value(bodyIndex)));
                this->postUpdateEvent(pe);

                Sleep(MESHING_STATUS_MESSAGE_PERSISTANCE);

                //! ------------------------
                //! the surface mesh engine
                //! ------------------------
                Property::meshEngine2D MeshEngine2D = myMeshDB->ArrayOfMesh2DEngine.value(bodyIndex);

                //! ---------------------------------------------------------------
                //! set the "to be updated" flag: value = "true" since the meshing
                //! process is not considered as finished without a volume mesh
                //! ---------------------------------------------------------------
                myMeshDB->ArrayOfMeshIsToBeUdpdated.insert(bodyIndex,true);

                //! ---------------------------------------------------------------------
                //! useBRep ON  => mesh CAD a geometry (use the boundary representation)
                //! useBRep OFF => mesh from a tessellation of the boundary
                //! ---------------------------------------------------------------------
                if(useBRep)
                {
                    switch(MeshEngine2D)
                    {
                    case Property::meshEngine2D_Netgen:
                    {
                        //! --------------------------------------------------------------
                        //! the progress indicator knows which is the running mesh engine
                        //! --------------------------------------------------------------
                        QProgressEvent *pe = new QProgressEvent(QProgressEvent_Init,0,NbSteps-1,0,
                                                                "Meshing surface using Netgen",
                                                                QProgressEvent_Init,0,100,0,"Netgen meshing");
                        this->postUpdateEvent(pe);
                        userMessage mr = Netgen_generateSurfaceMesh(bodyIndex,mainMesh2D);
                        Global::status().myMessages->appendMessage(mr);
                    }
                        break;

                    case Property::meshEngine2D_OCC_ExpressMesh:
                    {
                        QProgressEvent *pe = new QProgressEvent(QProgressEvent_Init,0,NbSteps-1,0,"Meshing surface using Express Mesh",QProgressEvent_Init,0,100,0);
                        postUpdateEvent(pe);
                        userMessage mr = this->ExMesh_generateSurfaceMesh(mp,bodyIndex,mainMesh2D);
                        Global::status().myMessages->appendMessage(mr);
                    }
                        break;
                    }
                }
                else
                {
                    switch(MeshEngine2D)
                    {
                    case Property::meshEngine2D_Netgen_STL:
                    {
                        //! --------------------------------------------------------------
                        //! the progress indicator knows which is the running mesh engine
                        //! --------------------------------------------------------------
                        QProgressEvent *pe = new QProgressEvent(QProgressEvent_Init,0,NbSteps-1,0,
                                                                "Meshing surface using Netgen STL",
                                                                QProgressEvent_Init,0,100,0,"Netgen meshing");
                        postUpdateEvent(pe);

                        userMessage mr = Netgen_STL_generateSurfaceMesh(bodyIndex,mainMesh2D);
                        Global::status().myMessages->appendMessage(mr);

                        //! ----------------------------------------------------------
                        //! experimental: project the surface nodes onto the geometry
                        //! ----------------------------------------------------------
                        bool meshCorrectionOn = myMeshDB->mapOfGeometryCorrection.value(bodyIndex);
                        if(meshCorrectionOn)
                        {
                            surfaceMeshSmoother *aSurfaceMeshSmoother = new surfaceMeshSmoother(myProgressIndicator,this);
                            aSurfaceMeshSmoother->perform(myMeshDB,bodyIndex);
                        }
                    }
                        break;
                    }
                }

                //! ---------------------------------------------------------
                //! after the surface mesh has been generated, generate
                //! the boundary mesh, prismatic or tetrahedral, if required
                //! ---------------------------------------------------------
                bool isInflationDefinedOnBody = myMeshDB->HasPrismaticFaces.value(bodyIndex);
                if(isInflationDefinedOnBody)
                {
                    occHandle(Ng_MeshVS_DataSourceFace) theLastInflatedMesh;
                    QList<occHandle(Ng_MeshVS_DataSourceFace)> listOfInflatedMeshes;
                    userMessage mr = this->PrismaticLayers_generatePrismaticMesh(bodyIndex,theLastInflatedMesh,listOfInflatedMeshes);

                    // for diagnostic
                    myMeshDB->ArrayOfMeshDS2D.insert(bodyIndex,theLastInflatedMesh);

                    Global::status().myMessages->appendMessage(mr);
                }

                //! ---------------------------------------------------------------
                //! set the "to be updated" flag: value = "true" since the meshing
                //! process is not considered as finished without a volume mesh
                //! ---------------------------------------------------------------
                myMeshDB->ArrayOfMeshIsToBeUdpdated.insert(bodyIndex, true);

                //! --------------------
                //! update the progress
                //! --------------------
                pe = new QProgressEvent();
                pe->setVal(done);
                pe->setMessage(QString("Surface mesh of body %1 done").arg(myMeshDB->MapOfBodyNames.value(bodyIndex)));
                postUpdateEvent(pe);
                Sleep(MESHING_STATUS_MESSAGE_PERSISTANCE);
            }
            //! -------------------------
            //! generate the volume mesh
            //! -------------------------
            else
            {
                //! --------------------
                //! update the progress
                //! --------------------
                QProgressEvent *pe = new QProgressEvent();
                pe->setVal(done);
                pe->setMessage(QString("Start meshing body %1").arg(myMeshDB->MapOfBodyNames.value(bodyIndex)));
                postUpdateEvent(pe);
                Sleep(MESHING_STATUS_MESSAGE_PERSISTANCE);

                //! ------------------------------------
                //! the surface and volume mesh engines
                //! ------------------------------------
                Property::meshEngine2D MeshEngine2D = myMeshDB->ArrayOfMesh2DEngine.value(bodyIndex);
                Property::meshEngine3D MeshEngine3D = myMeshDB->ArrayOfMesh3DEngine.value(bodyIndex);

                //! -------------------------------------------------
                //! set the initial update flag for the current body
                //! -------------------------------------------------
                myMeshDB->ArrayOfMeshIsToBeUdpdated.insert(bodyIndex,true);

                switch(MeshEngine2D)
                {
                //! -------------------------
                //! Netgen as surface mesher
                //! -------------------------
                case Property::meshEngine2D_Netgen:
                {
                    switch(MeshEngine3D)
                    {
                    case Property::meshEngine3D_Netgen:
                    {
                        userMessage mr = this->Netgen_generateVolumeMesh(bodyIndex,mainMesh2D,mainMesh3D);
                        if(mr.isDone) myMeshDB->ArrayOfMeshIsToBeUdpdated.insert(bodyIndex,false);
                        Global::status().myMessages->appendMessage(mr);
                    }
                        break;

                    case Property::meshEngine3D_Tetgen:
                    {
                        userMessage mr = this->Netgen_generateSurfaceMesh(bodyIndex,mainMesh2D);
                        if(mr.isDone)
                        {
                            //! -----------------------------------------------
                            //! the surface mesh has been generated by Netgen:
                            //! It's a quality mesh: do not add SPs on surface
                            //! -----------------------------------------------
                            int preserveSurfaceMesh = 1;
                            userMessage mr = this->Tetgen_generateVolumeMesh(bodyIndex,preserveSurfaceMesh,isMeshingInMemory,done);
                            if(mr.isDone) myMeshDB->ArrayOfMeshIsToBeUdpdated.insert(bodyIndex,false);
                        }
                        Global::status().myMessages->appendMessage(mr);
                    }
                        break;
                    }
                }
                    break;

                case Property::meshEngine2D_OCC_ExpressMesh:
                {
                    switch (MeshEngine3D)
                    {
                    case Property::meshEngine3D_Netgen:
                    {
                        userMessage mr = this->ExMesh_generateVolumeMesh(mp,bodyIndex,mainMesh2D,mainMesh3D,done);
                        if(mr.isDone) myMeshDB->ArrayOfMeshIsToBeUdpdated.insert(bodyIndex,false);
                        Global::status().myMessages->appendMessage(mr);
                    }
                        break;

                    case Property::meshEngine3D_Tetgen:
                    {
                        userMessage mr = this->ExMesh_generateSurfaceMesh(mp,bodyIndex,mainMesh2D);
                        if(mr.isDone)
                        {
                            //! ----------------------------------------------------
                            //! the surface mesh has been generated by Express mesh
                            //! It's supposed to be a quality mesh BUT the surface
                            //! mesh is actually not suitable for FEM/FVM: add SPs
                            //! ----------------------------------------------------
                            int preserveSurfaceMesh = 0;
                            userMessage mr = this->Tetgen_generateVolumeMesh(bodyIndex,preserveSurfaceMesh,isMeshingInMemory,done);
                            if(mr.isDone) myMeshDB->ArrayOfMeshIsToBeUdpdated.insert(bodyIndex,false);
                        }
                        Global::status().myMessages->appendMessage(mr);
                    }
                        break;
                    }
                }
                    break;

                //! --------------------------
                //! Netgen STL surface mesher
                //! --------------------------
                case Property::meshEngine2D_Netgen_STL:
                {
                    switch(MeshEngine3D)
                    {
                    case Property::meshEngine3D_Netgen_STL:
                    {
                        userMessage mr = this->Netgen_STL_generateVolumeMesh(bodyIndex,mainMesh2D,mainMesh3D);
                        if(mr.isDone) myMeshDB->ArrayOfMeshIsToBeUdpdated.insert(bodyIndex,false);
                        Global::status().myMessages->appendMessage(mr);

                        //! ----------------------------------------------------------
                        //! experimental: project the surface nodes onto the geometry
                        //! ----------------------------------------------------------
                        bool meshCorrectionOn = myMeshDB->mapOfGeometryCorrection.value(bodyIndex);
                        if(meshCorrectionOn)
                        {
                            surfaceMeshSmoother *aSurfaceMeshSmoother = new surfaceMeshSmoother(myProgressIndicator,this);
                            aSurfaceMeshSmoother->perform(myMeshDB,bodyIndex);
                        }
                    }
                    break;

                    case Property::meshEngine3D_Tetgen:
                    {
                        //! -------------------------------------------
                        //! generate the surface mesh using Netgen STL
                        //! -------------------------------------------
                        userMessage mr = this->Netgen_STL_generateSurfaceMesh(bodyIndex,mainMesh2D);
                        if(mr.isDone)
                        {
                            //! -----------------------------------------------
                            //! the surface mesh has been generated by Netgen:
                            //! preserve this quality surface mesh
                            //! -----------------------------------------------
                            int preserveSurfaceMesh = 0;
                            userMessage mr = this->Tetgen_generateVolumeMesh(bodyIndex,preserveSurfaceMesh,isMeshingInMemory,done);
                            if(mr.isDone)
                            {
                                myMeshDB->ArrayOfMeshIsToBeUdpdated.insert(bodyIndex,false);
                            }
                            Global::status().myMessages->appendMessage(mr);
                        }
                    }
                        break;
                    }
                    break;
                }

                //! ----------------------------------------------------------
                //! NULL surface mesher: this option is automatically handled
                //! within the interface when the volume mesher is Tetgen BR
                //! or Tetwild, since they only need the tessellation
                //! ----------------------------------------------------------
                case Property::meshEngine2D_NULL:
                {
                    if(MeshEngine3D == Property::meshEngine3D_Tetgen_BR)
                    {
                        //! ---------------------------------
                        //! do preliminary work on .stl mesh
                        //! ---------------------------------
                        QString inputExtendedSTL_S = this->processSTL(bodyIndex);

                        //! ----------------------------------------------------------------------------------
                        //! read the extended .stl file and return the surface mesh and face mesh datasources
                        //! ----------------------------------------------------------------------------------
                        occHandle(Ng_MeshVS_DataSource2D) anSTL_MeshVS_DataSource;
                        int NbGeometryFaces = myMeshDB->MapOfBodyTopologyMap.value(bodyIndex).faceMap.Extent();

                        std::vector<occHandle(Ng_MeshVS_DataSourceFace)> array_STL_meshDS;
                        for(int i=0; i<=NbGeometryFaces; i++) array_STL_meshDS.push_back(occHandle(Ng_MeshVS_DataSourceFace)());
                        const TopoDS_Shape &shape = myMeshDB->bodyMap.value(bodyIndex);
                        bool isDone = MeshTools::toSTLMesh1(shape, inputExtendedSTL_S, anSTL_MeshVS_DataSource, array_STL_meshDS);

                        if(isDone)
                        {
                            //! ---------------------------------------
                            //! put the surface mesh into the database
                            //! ---------------------------------------
                            mainMesh2D = anSTL_MeshVS_DataSource;
                            myMeshDB->ArrayOfMeshDS2D.insert(bodyIndex,mainMesh2D);

                            //! ------------------------------------------------
                            //! put the face mesh datasources into the database
                            //! ------------------------------------------------
                            for(int faceNr = 0; faceNr<=NbGeometryFaces; faceNr++)
                            {
                                if(array_STL_meshDS.at(faceNr).IsNull() == true) continue;
                                myMeshDB->ArrayOfMeshDSOnFaces.setValue(bodyIndex,faceNr,array_STL_meshDS.at(faceNr));
                            }

                            //! ---------------------------------------------
                            //! the main calls for generating the volume mesh
                            //! the input tassellation will not be respected
                            //! Tetgen will add points to the surface mesh
                            //! ---------------------------------------------
                            userMessage mr;
                            int preserveSurfaceMesh = 0;
                            bool runInMemory = myMeshDB->MapOfIsMeshingRunningInMemory.value(bodyIndex);
                            mr = Tetgen_generateVolumeMesh(bodyIndex, preserveSurfaceMesh, runInMemory, done);
                            Global::status().myMessages->appendMessage(mr);
                        }
                        else
                        {
                            userMessage mr(false,QString("Cannot generate the volume mesh with Tetgen: error in building the PLC"));
                            myMeshDB->ArrayOfMeshIsToBeUdpdated.insert(bodyIndex, true);
                            Global::status().myMessages->appendMessage(mr);
                        }

                        //! ----------------------------------------------------------
                        //! experimental: project the surface nodes onto the geometry
                        //! ----------------------------------------------------------
                        //surfaceMeshSmoother::perform(myMeshDB,bodyIndex);

                        //! ---------------------------------------------
                        //! remove the Tetgen support files, if existing
                        //! ---------------------------------------------
                        this->removeSupportFiles();
                    }
                    if(MeshEngine3D == Property::meshEngine3D_TetWild)
                    {
                        bool runInMemory = myMeshDB->MapOfIsMeshingRunningInMemory.value(bodyIndex);

                        //! -------------------------------------------------------
                        //! retrieve all the TetWild parameters from the data base
                        //! only the ones which are of interests will be used
                        //! -------------------------------------------------------
                        int envelopeSizingType = myMeshDB->mapOfEnvelopeSizingType.value(bodyIndex);
                        double relativeEnvelopeSize = myMeshDB->mapOfRelativeEnvelopeSize.value(bodyIndex);
                        double absoluteEnvelopeSize = myMeshDB->mapOfAbsoluteEnvelopeSize.value(bodyIndex);

                        int idealLenghtSizing = myMeshDB->mapOfIdealLengthSizingType.value(bodyIndex);
                        double relativeLength = myMeshDB->mapOfIdealLengthRelativeSize.value(bodyIndex);
                        double absoluteLength = myMeshDB->mapOfIdealLengthAbsoluteSize.value(bodyIndex);

                        //! --------------------------------------
                        //! clean the tessellation from the shape
                        //! --------------------------------------
                        if(runInMemory==false)
                        {
                            //! ---------------------------------------------------------------------------------------
                            //! TetWild run on disk
                            //! The TetWild command line accects "relative" envelope sizg and "relative" initial ideal
                            //! legnth, by consequence, if into the interface "absolute" values switch are active, the
                            //! absolute envelope sizing and absolute ideal length should be converted into relative
                            //! values, according to their internal definition within the mesher
                            //! ---------------------------------------------------------------------------------------
                            BRepTools::Clean(myMeshDB->bodyMap.value(bodyIndex));
                            QString inputSoup = this->processTetWildSTL(bodyIndex);
                            tetWildMesher aTetWildMesher;

                            //! ---------------------
                            //! test some parameters
                            //! ---------------------
                            double envelopeSize;    //! in TetWild command line this parameter is relative
                            double idealLength;     //! in TetWild command line this parameter is relative

                            if(envelopeSizingType==0) envelopeSize = relativeEnvelopeSize;
                            else
                            {
                                //! ---------------------------------------------------
                                //! find the diagonal of the bounding box of the shape
                                //! ---------------------------------------------------
                                const TopoDS_Shape &aShape = myMeshDB->bodyMap.value(bodyIndex);
                                double L1,L2,L3,D;
                                GeomToolsClass::getBoundingBox(aShape,L1,L2,L3);
                                D = sqrt(L1*L1+L2*L2+L3*L3);
                                envelopeSize = absoluteEnvelopeSize/D;
                            }
                            if(idealLenghtSizing==0) idealLength = relativeLength;
                            else
                            {
                                //! ---------------------------------------------------
                                //! find the diagonal of the bounding box of the shape
                                //! ---------------------------------------------------
                                const TopoDS_Shape &aShape = myMeshDB->bodyMap.value(bodyIndex);
                                double L1,L2,L3,D;
                                GeomToolsClass::getBoundingBox(aShape,L1,L2,L3);
                                D = sqrt(L1*L1+L2*L2+L3*L3);
                                idealLength = absoluteLength/D;
                            }
                            aTetWildMesher.setParameters(idealLength,envelopeSize);

                            bool isDone = aTetWildMesher.perform_onDisk(inputSoup);
                            if(isDone)
                            {
                                bool isMeshReady = aTetWildMesher.retrieveMeshDataSourceFromDisk(mainMesh3D);

                                //! insert a valid or null 3D data source into the mesh data base
                                myMeshDB->ArrayOfMeshDS.insert(bodyIndex,mainMesh3D);

                                //! ------------------------------------------------------------------
                                //! insert a valid or null 2D mesh data source into the data base
                                //! using the following constructor: it builds the surface mesh using
                                //! the connectivity info (they must be generated here [*])
                                //! ------------------------------------------------------------------
                                mainMesh3D->buildFaceToElementConnectivity();
                                mainMesh2D = new Ng_MeshVS_DataSource2D(mainMesh3D);

                                cout<<"@____surface mesh: "<<mainMesh2D->GetAllElements().Extent()<<" elements; "<<mainMesh2D->GetAllNodes().Extent()<<" nodes"<<endl;

                                myMeshDB->ArrayOfMeshDS2D.insert(bodyIndex,mainMesh2D);

                                if(isMeshReady) myMeshDB->ArrayOfMeshIsToBeUdpdated.insert(bodyIndex, false);
                                else myMeshDB->ArrayOfMeshIsToBeUdpdated.insert(bodyIndex, true);
                                userMessage mr(true,QString("Tetwild has successfully generated the volume mesh"));
                                Global::status().myMessages->appendMessage(mr);
                            }
                            else
                            {
                                myMeshDB->ArrayOfMeshIsToBeUdpdated.insert(bodyIndex, true);
                                userMessage mr(false,QString("Cannot generate the volume mesh with Tetwild"));
                                Global::status().myMessages->appendMessage(mr);
                            }
                        }
                        else
                        {
                            //! --------------------------------------------
                            //! TetWild run in memory - not implemented yet
                            //! --------------------------------------------
                            myMeshDB->ArrayOfMeshIsToBeUdpdated.insert(bodyIndex, false);
                        }
                    }
                }
                    break;
                }
            }

            //! ---------------------------------------------------
            //! this will convert a tet mesh into a full hexa mesh
            //! ---------------------------------------------------
            if(myMeshDB->ArrayOfMeshType.value(bodyIndex) == 1)
            {
                TetHex aTetHex(myMeshDB);
                aTetHex.setProgressIndicator(myProgressIndicator);
                userMessage mr = aTetHex.perform(bodyIndex);

                if(mr.isDone) myMeshDB->ArrayOfMeshIsToBeUdpdated.insert(bodyIndex, false);
            }

            /*
            //! -------------------------------------------------------------
            //! for testing element to element connectivity - can be removed
            //! -------------------------------------------------------------
            char fname[128];
            sprintf(fname,"D:\\elToEl_%d.txt",bodyIndex);
            FILE *file = fopen(fname,"w");
            cout<<endl<<"@____printing element to element connectivity____"<<endl;
            const occHandle(Ng_MeshVS_DataSource3D) &m = occHandle(Ng_MeshVS_DataSource3D)::DownCast(myMeshDB->ArrayOfMeshDS.value(bodyIndex));
            std::map<int,std::vector<int>> cmap;
            m->buildElementToElementConnectivity(cmap);
            for(std::map<int,std::vector<int>>::iterator it = cmap.begin(); it!=cmap.end(); ++it)
            {
                std::pair<int,std::vector<int>> aPair = *it;
                int elementID = aPair.first;
                std::vector<int> attachedElements = aPair.second;
                fprintf(file,"%d\t=>\t",elementID);
                //cout<<"____"<<elementID<<"____"<<endl;
                int Nb = int(attachedElements.size());
                for(int i=0; i<Nb-1; i++) fprintf(file,"%d\t",attachedElements[i]);
                fprintf(file,"%d\n",attachedElements[Nb-1]);
            }
            cout<<"@____end of printing element to element connectivity____"<<endl<<endl;
            fclose(file);
            //! ------------
            //! end testing
            //! ------------

            //! -----------------------------------------------------------------------
            //! moreover: build a layer of elements using the face numbers is vecFaces
            //! -----------------------------------------------------------------------
            std::vector<int> vecFaces{1,6};
            QList<occHandle(Ng_MeshVS_DataSourceFace)> listOfFaces;
            for(std::vector<int>::iterator it = vecFaces.begin(); it!=vecFaces.end(); it++)
            {
                int faceNr = *it;
                occHandle(Ng_MeshVS_DataSourceFace) curFace = occHandle(Ng_MeshVS_DataSourceFace)::DownCast(myMeshDB->ArrayOfMeshDSOnFaces.getValue(bodyIndex,faceNr));
                listOfFaces<<curFace;
            }
            occHandle(Ng_MeshVS_DataSourceFace) mergedFaceMeshDS = new Ng_MeshVS_DataSourceFace(listOfFaces);
            occHandle(Ng_MeshVS_DataSource3D) aVolumeMesh = occHandle(Ng_MeshVS_DataSource3D)::DownCast(myMeshDB->ArrayOfMeshDS.value(1));
            meshSelectNLayers selector(aVolumeMesh);
            selector.setSteps(2);
            selector.setStartingSurfaceElements(mergedFaceMeshDS);
            occHandle(Ng_MeshVS_DataSource3D) selectedElements;
            bool isDone = selector.getElements(selectedElements);
            if(isDone) myMeshDB->ArrayOfMeshDS.insert(1,selectedElements);
            cout<<"____number of selected elements: "<<selectedElements->GetAllElements().Extent()<<"____"<<endl;
            //!-------------
            //! end testing
            //! ------------
            */
            //! ------------------------------------------
            //! end of the surface/volume meshing process
            //! ------------------------------------------

            if(isSecondOrder)
            {
                //! ----------------------------------
                //! a tool for renumbering mesh nodes
                //! ----------------------------------
                meshNodesRenumberingTool aTool;
                if(myIsVolume)
                {
                    //! ----------------------------------------
                    //! convert the volume mesh into n-th order
                    //! ----------------------------------------
                    occHandle(Ng_MeshVS_DataSource3D) firstOrderVolumeMesh = occHandle(Ng_MeshVS_DataSource3D)::DownCast(myMeshDB->ArrayOfMeshDS.value(bodyIndex));
                    occHandle(Ng_MeshVS_DataSource3D) secondOrderVolumeMeshDS = new Ng_MeshVS_DataSource3D(firstOrderVolumeMesh,2);
                    myMeshDB->ArrayOfMeshDS.insert(bodyIndex,secondOrderVolumeMeshDS);

                    //! ------------------------------------------------
                    //! build the n-th order surface mesh - topological
                    //! no need of renumbering
                    //! -----------------------------------------------
                    //secondOrderVolumeMeshDS->buildFaceToElementConnectivity();
                    //occHandle(Ng_MeshVS_DataSource2D) secondOrderSurfaceMeshDS = new Ng_MeshVS_DataSource2D(secondOrderVolumeMeshDS);
                    //myMeshDB->ArrayOfMeshDS2D.insert(bodyIndex,secondOrderSurfaceMeshDS);

                    //! -----------------------------------------
                    //! convert the surface mesh into n-th order
                    //! -----------------------------------------
                    occHandle(Ng_MeshVS_DataSource2D) firstOrderSurfaceMesh = occHandle(Ng_MeshVS_DataSource2D)::DownCast(myMeshDB->ArrayOfMeshDS2D.value(bodyIndex));
                    occHandle(Ng_MeshVS_DataSource2D) secondOrderSurfaceMeshDS = new Ng_MeshVS_DataSource2D(firstOrderSurfaceMesh,2);
                    aTool.setSource(secondOrderVolumeMeshDS);
                    aTool.perform<Ng_MeshVS_DataSource2D>(secondOrderSurfaceMeshDS.get());

                    /*
                    //! ---------------
                    //! renumber nodes
                    //! ---------------
                    std::map<size_t,int> mapNodeCoords_globalNodeID;
                    std::vector<int> globalNodeIDs;
                    std::vector<mesh::tolerantPoint> vecTolerantPoints;
                    for(TColStd_MapIteratorOfPackedMapOfInteger it(secondOrderVolumeMeshDS->GetAllNodes()); it.More(); it.Next())
                    {
                        int globalNodeID = it.Key();
                        int localNodeID = secondOrderVolumeMeshDS->myNodesMap.FindIndex(globalNodeID);
                        const std::vector<double> &aPoint = secondOrderVolumeMeshDS->getNodeCoordinates(localNodeID);

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

                    secondOrderSurfaceMeshDS->renumberNodes(mapNodeCoords_globalNodeID,vecTolerantPoints,globalNodeIDs);
                    */
                    myMeshDB->ArrayOfMeshDS2D.insert(bodyIndex,secondOrderSurfaceMeshDS);

                    /*
                    //! ------------------------------------------
                    //! prepare information for renumbering faces
                    //! ------------------------------------------
                    std::map<size_t,int> mapNodeCoords_globalNodeID;
                    std::vector<int> globalNodeIDs;
                    std::vector<mesh::tolerantPoint> vecTolerantPoints;
                    mapNodeCoords_globalNodeID.clear();
                    globalNodeIDs.clear();
                    vecTolerantPoints.clear();
                    for(TColStd_MapIteratorOfPackedMapOfInteger it(secondOrderSurfaceMeshDS->GetAllNodes()); it.More(); it.Next())
                    {
                        int globalNodeID = it.Key();
                        int localNodeID = secondOrderSurfaceMeshDS->myNodesMap.FindIndex(globalNodeID);
                        const std::vector<double> &aPoint = secondOrderSurfaceMeshDS->getNodeCoordinates(localNodeID);

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
                    */

                    //! ----------------------------------
                    //! build the face mesh data sources
                    //! and renumber nodes
                    //! ----------------------------------
                    int NbFaces = myMeshDB->MapOfBodyTopologyMap.value(bodyIndex).faceMap.Extent();
                    for(int faceNr = 1; faceNr<=NbFaces; faceNr++)
                    {
                        if(myMeshDB->ArrayOfMeshDSOnFaces.getValue(bodyIndex,faceNr).IsNull()) continue;

                        occHandle(Ng_MeshVS_DataSourceFace) aFaceFirstOrderFaceMesh =
                                occHandle(Ng_MeshVS_DataSourceFace)::DownCast(myMeshDB->ArrayOfMeshDSOnFaces.getValue(bodyIndex,faceNr));
                        occHandle(Ng_MeshVS_DataSourceFace) aSecondOrderFaceMesh =
                                new Ng_MeshVS_DataSourceFace(aFaceFirstOrderFaceMesh,2);

                        aTool.perform<Ng_MeshVS_DataSourceFace>(aSecondOrderFaceMesh.get());
                        //aSecondOrderFaceMesh->renumberNodes(mapNodeCoords_globalNodeID,vecTolerantPoints,globalNodeIDs);
                        myMeshDB->ArrayOfMeshDSOnFaces.setValue(bodyIndex,faceNr,aSecondOrderFaceMesh);
                    }

                    //! ----------------
                    //! curved elements
                    //! ----------------
                    bool curvedElements = myMeshDB->mapOfIsElementStraight.value(bodyIndex);
                    if(curvedElements)
                    {
                        ;
                    }
                }
                else
                {
                    //! -----------------------------------------------------------------------
                    //! volume mesh is not required: convert the surface mesh into n-th order
                    //! -----------------------------------------------------------------------
                    occHandle(Ng_MeshVS_DataSource2D) firstOrderSurfaceMesh = occHandle(Ng_MeshVS_DataSource2D)::DownCast(myMeshDB->ArrayOfMeshDS2D.value(bodyIndex));
                    occHandle(Ng_MeshVS_DataSource2D) secondOrderSurfaceMeshDS = new Ng_MeshVS_DataSource2D(firstOrderSurfaceMesh,2);
                    myMeshDB->ArrayOfMeshDS2D.insert(bodyIndex,secondOrderSurfaceMeshDS);

                    /*
                    //! ------------------------------------------
                    //! prepare information for renumbering faces
                    //! ------------------------------------------
                    std::map<size_t,int> mapNodeCoords_globalNodeID;
                    std::vector<int> globalNodeIDs;
                    std::vector<mesh::tolerantPoint> vecTolerantPoints;
                    for(TColStd_MapIteratorOfPackedMapOfInteger it(secondOrderSurfaceMeshDS->GetAllNodes()); it.More(); it.Next())
                    {
                        int globalNodeID = it.Key();
                        int localNodeID = secondOrderSurfaceMeshDS->myNodesMap.FindIndex(globalNodeID);
                        const std::vector<double> &aPoint = secondOrderSurfaceMeshDS->getNodeCoordinates(localNodeID);

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
                    */
                    aTool.setSource(secondOrderSurfaceMeshDS);

                    //! ----------------------------------
                    //! build the face mesh data sources
                    //! and renumber nodes
                    //! ----------------------------------
                    int NbFaces = myMeshDB->MapOfBodyTopologyMap.value(bodyIndex).faceMap.Extent();
                    for(int faceNr = 1; faceNr<=NbFaces; faceNr++)
                    {
                        occHandle(Ng_MeshVS_DataSourceFace) aFaceFirstOrderFaceMesh =
                                occHandle(Ng_MeshVS_DataSourceFace)::DownCast(myMeshDB->ArrayOfMeshDSOnFaces.getValue(bodyIndex,faceNr));
                        occHandle(Ng_MeshVS_DataSourceFace) aSecondOrderFaceMesh =
                                new Ng_MeshVS_DataSourceFace(aFaceFirstOrderFaceMesh,2);

                        //aSecondOrderFaceMesh->renumberNodes(mapNodeCoords_globalNodeID,vecTolerantPoints,globalNodeIDs);
                        aTool.perform<Ng_MeshVS_DataSourceFace>(aSecondOrderFaceMesh.get());
                        myMeshDB->ArrayOfMeshDSOnFaces.setValue(bodyIndex,faceNr,aSecondOrderFaceMesh);
                    }
                }
            }



#ifdef CREATE1D_DATASOURCES

            //! for the moment the generation of the edge mesh datasources
            //! has been limited to the surface mesh step, because of an unknown origin error
            //! occurring when generating the volume mesh
            QMap<int, QList<meshPoint>> edgeNr_to_listOfPoints;
            if(myIsVolume==false)
            {
                //! build 1D meshes. A 1D mesh if generated from the surface mesh
                if(!mainMesh2D.IsNull())
                {
                    edgeMeshDataSourceBuildesClass edgeMeshBuilder(myMeshDB);
                    edgeMeshBuilder.perform();
                }
                else
                {
                    cout<<"MesherClass::MesherClass()->____surface mesh is NULL: cannot generate 1D meshes____"<<endl;
                }
            }
            else
            {
                cout<<"MesherClass::MesherClass()->____edge mesh datasources are not generated in case of volume mesh____"<<endl;
                ccout("MesherClass::MesherClass()->____edge mesh datasources are not generated in case of volume mesh____");
            }
#endif

        }
        else
        {
            QProgressEvent *pe = new QProgressEvent();
            pe->setVal(done);
            pe->setMessage(QString("Interrupt request received"));
            postUpdateEvent(pe);
            Sleep(MESHING_STATUS_MESSAGE_PERSISTANCE*2);

            cout<<"\\-------------------------------------------\\"<<endl;
            cout<<"\\- meshing process interrupted by the user -\\"<<endl;
            cout<<"\\-------------------------------------------\\"<<endl;
            break;
        }
        done++;
    }

    //! ------------------------
    //! send event meshing done
    //! ------------------------
    pe = new QProgressEvent();
    pe->setVal(done);
    pe->setMessage(QString("Meshing done"));
    postUpdateEvent(pe);
    Sleep(MESHING_STATUS_MESSAGE_PERSISTANCE*2);

    //! ------------------------------------------
    //! update progress - remove the progress bar
    //! ------------------------------------------
    pe = new QProgressEvent(QProgressEventAction::QProgressEvent_Reset,0,100,0,"",QProgressEventAction::QProgressEvent_Reset,0,100,0,"");
    postUpdateEvent(pe);

    //! ---------------------
    //! emit meshingFinished
    //! ---------------------
    emit meshingFinished();

    //! --------------------------
    //! re-init the global status
    //! --------------------------
    Global::status().code = 1;
    Global::status().task="";
}

//! -------------------------------------
//! function: Netgen_generateSurfaceMesh
//! details:
//! -------------------------------------
userMessage MesherClass::Netgen_generateSurfaceMesh(int bodyIndex, occHandle(MeshVS_DataSource) &mainMesh2D)
{
    myNetgenMesher->setMeshDataBase(myMeshDB);
    occHandle(Ng_MeshVS_DataSource3D) dummyDS;
    bool isVolume = false;
    userMessage mr = myNetgenMesher->performBrep(bodyIndex,isVolume,mainMesh2D,dummyDS);
    if(mr.isDone)
    {
        myMeshDB->ArrayOfMeshDS2D.insert(bodyIndex,mainMesh2D);
        myMeshDB->ArrayOfMeshDS.insert(bodyIndex,dummyDS);
    }
    return mr;
}

//! -----------------------------------------
//! function: Netgen_STL_generateSurfaceMesh
//! details:
//! -----------------------------------------
userMessage MesherClass::Netgen_STL_generateSurfaceMesh(int bodyIndex,
                                                        occHandle(MeshVS_DataSource) &mainMesh2D)
{
    //! ------------------------------------------------
    //! do work on the .stl mesh used for visualization
    //! ------------------------------------------------
    QString inputExtendedSTL_S = this->processSTL(bodyIndex);

    //! ----------------------------------------------------------------------------
    //! read the extended (with connectivity info), healed and simplified .stl file
    //! and return the surface mesh and face mesh datasources
    //! ----------------------------------------------------------------------------
    occHandle(Ng_MeshVS_DataSource2D) anSTL_MeshVS_DataSource;

    //! ---------------------------------
    //! vector of face mesh data sources
    //! ---------------------------------
    std::vector<occHandle(Ng_MeshVS_DataSourceFace)> array_STL_meshDS;
    int NbGeometryFaces = myMeshDB->MapOfBodyTopologyMap.value(bodyIndex).faceMap.Extent();
    for(int i=0; i<=NbGeometryFaces; i++) array_STL_meshDS.push_back(occHandle(Ng_MeshVS_DataSourceFace)());
    const TopoDS_Shape &shape = myMeshDB->bodyMap.value(bodyIndex);
    bool isDone = MeshTools::toSTLMesh1(shape, inputExtendedSTL_S, anSTL_MeshVS_DataSource, array_STL_meshDS);

    if(!isDone)
    {
        cout<<"MesherClass::Netgen_STL_generateSurfaceMesh()->____cannot start the meshing process: the .stl mesh input cannot be generated___"<<endl;
        userMessage mr(false,"cannot start the meshing process: the .stl mesh input cannot be generated");
        return mr;
    }

    //! ------------------------------------------------------------------
    //! put the surface and face mesh datasources read from the support
    //! file into the database: they will be read by the tool
    //! ------------------------------------------------------------------
    mainMesh2D = anSTL_MeshVS_DataSource;
    myMeshDB->ArrayOfMeshDS2D.insert(bodyIndex,mainMesh2D);
    for(int faceNr = 0; faceNr<=NbGeometryFaces; faceNr++)
    {
        if(!array_STL_meshDS.at(faceNr).IsNull())
        {
            //cerr<<"____face nr: "<<faceNr<<" nodes: "<<array_STL_meshDS.at(faceNr)->GetAllNodes().Extent()<<"____"<<endl;
            myMeshDB->ArrayOfMeshDSOnFaces.setValue(bodyIndex,faceNr,array_STL_meshDS.at(faceNr));
        }
        else
        {
            cerr<<"____face nr: "<<faceNr<<"  NULL____"<<endl;
        }
    }

    //! ----------------------------------------------
    //! the main call for generating the surface mesh
    //! ----------------------------------------------
    myNetgenMesher->setMeshDataBase(myMeshDB);
    occHandle(Ng_MeshVS_DataSource3D) dummyDS;
    bool isVolume = false;
    userMessage mr = myNetgenMesher->perform(bodyIndex,isVolume,mainMesh2D,dummyDS);
    if(mr.isDone)
    {
        myMeshDB->ArrayOfMeshDS2D.insert(bodyIndex,mainMesh2D);
        myMeshDB->ArrayOfMeshDS.insert(bodyIndex,dummyDS);
    }

    //! --------------------------------------------------
    //! clear the tessellation if the meshing process has
    //! not been successfull
    //! --------------------------------------------------
    if(mr.isDone==false)
    {
        for(int faceNr = 0; faceNr<=NbGeometryFaces; faceNr++)
            myMeshDB->ArrayOfMeshDSOnFaces.setValue(bodyIndex,faceNr,occHandle(MeshVS_DataSource)());
    }

    //! ----------------------
    //! remove supports files
    //! ----------------------
    this->removeSupportFiles();
    return mr;
}

//! -------------------------------------
//! function: ExMesh_generateSurfaceMesh
//! details:
//! -------------------------------------
userMessage MesherClass::ExMesh_generateSurfaceMesh(const meshParam &mp,
                                                      int bodyIndex,
                                                      occHandle(MeshVS_DataSource) &mainMesh2D)
{
    TopoDS_Shape shape = myMeshDB->bodyMap.value(bodyIndex);
    occMesher theOCCMesher(shape,mp);
    theOCCMesher.setProgressIndicator(myProgressIndicator);

    //! -----------------------------------------------------------------
    //! "done" passed to ::performSurface(...) activate the progress bar
    //! -----------------------------------------------------------------
    bool deleteNetgenPointers = true;
    userMessage mr = theOCCMesher.performSurface(mainMesh2D,mp.isSecondOrder,mp.isElementStraight,deleteNetgenPointers);

    if(mr.isDone)
    {
        myMeshDB->ArrayOfMeshDS2D.insert(bodyIndex,mainMesh2D);
        QMap<int,occHandle(Ng_MeshVS_DataSourceFace)> faceDataSources = theOCCMesher.getFaceDataSources();
        for(QMap<int,occHandle(Ng_MeshVS_DataSourceFace)>::iterator it=faceDataSources.begin(); it!=faceDataSources.end(); ++it)
        {
            int faceNr = it.key();
            const occHandle(Ng_MeshVS_DataSourceFace) &aFaceMeshDS = it.value();
            myMeshDB->ArrayOfMeshDSOnFaces.setValue(bodyIndex,faceNr,aFaceMeshDS);
        }
    }
    return mr;
}

//! ----------------------------------
//! function: STL_generateSurfaceMesh
//! details:  unused?
//! ----------------------------------
bool MesherClass::STL_generateSurfaceMesh(const meshParam &mp, int bodyIndex, occHandle(MeshVS_DataSource) &mainMesh2D)
{
    cout<<"MesherClass::generateMesh()->____start generating the surface mesh using STL____"<<endl;

    Q_UNUSED(mp)

    QString inputExtendedSTL_S =  processSTL(bodyIndex);

    //! -------------------------------------------------------------
    //! the MeshTool::toSTLMesh returns the surface mesh data source
    //! and the array of face mesh datasources
    //! -------------------------------------------------------------
    occHandle(Ng_MeshVS_DataSource2D) anSTL_meshDS;
    int NbGeometryFaces = myMeshDB->MapOfBodyTopologyMap.value(bodyIndex).faceMap.Extent();
    const TopoDS_Shape &shape = myMeshDB->bodyMap.value(bodyIndex);

    std::vector<occHandle(Ng_MeshVS_DataSourceFace)> array_STL_meshDS;
    for(int i=0; i<=NbGeometryFaces; i++) array_STL_meshDS.push_back(occHandle(Ng_MeshVS_DataSourceFace)());
    bool isDone = MeshTools::toSTLMesh1(shape, inputExtendedSTL_S, anSTL_meshDS, array_STL_meshDS);

    if(!isDone) return false;

    mainMesh2D = anSTL_meshDS;
    myMeshDB->ArrayOfMeshDS2D.insert(bodyIndex,mainMesh2D);

    //! ----------------------------------
    //! extract the face mesh datasources
    //! ----------------------------------
    int Nfaces = myMeshDB->MapOfBodyTopologyMap.value(bodyIndex).faceMap.Extent();
    for(int faceNr = 0; faceNr<=Nfaces; faceNr++)
    {
        if(array_STL_meshDS.at(faceNr).IsNull())
        {
            myMeshDB->ArrayOfMeshDSOnFaces.setValue(bodyIndex,faceNr,occHandle(Ng_MeshVS_DataSourceFace)());
            continue;
        }
        if(array_STL_meshDS.at(faceNr)->GetAllElements().Extent()<1)
        {
            myMeshDB->ArrayOfMeshDSOnFaces.setValue(bodyIndex,faceNr,occHandle(Ng_MeshVS_DataSourceFace)());
            continue;
        }
        myMeshDB->ArrayOfMeshDSOnFaces.setValue(bodyIndex,faceNr,array_STL_meshDS.at(faceNr));
    }
    return true;
}

/*
//! --------------------------------------------------
//! function: PrismaticLayers_generateTetBoundaryMesh
//! details:
//! --------------------------------------------------
userMessage MesherClass::PrismaticLayers_generateTetBoundaryMesh(int bodyIndex)
{
    cout<<"MesherClass::PrismaticLayers_generateTetBoundaryMesh()->____function called____"<<endl;

    //! ---------------
    //! meshing result
    //! ---------------
    userMessage mr;

    //! ------------------------
    //! disable the Stop button
    //! ------------------------
    if(myProgressIndicator!=Q_NULLPTR) emit requestDisableStop();

    //! ---------------------------------------
    //! initialize the prismatic layer builder
    //! ---------------------------------------
    prismaticLayer prismaticLayerBuilder(myMeshDB);
    if(myProgressIndicator!=Q_NULLPTR) prismaticLayerBuilder.setProgressIndicator(myProgressIndicator);
    prismaticLayerBuilder.setBody(bodyIndex);

    //! ---------------------------
    //! set the progress indicator
    //! ---------------------------
    if(myProgressIndicator!=Q_NULLPTR) prismaticLayerBuilder.setProgressIndicator(myProgressIndicator);

    //! ------------------------------------------------------
    //! retrieve the list of the prismatic faces for the body
    //! ------------------------------------------------------
    const QList<int> &prismaticFacesOnBody = myMeshDB->prismaticFaces.value(bodyIndex);
    prismaticLayerBuilder.setPrismaticFaces(prismaticFacesOnBody);

    //! ---------------------------------------------------------------
    //! merge the prismatic face meshes into a unique mesh data source
    //! ---------------------------------------------------------------
    bool isDone = prismaticLayerBuilder.mergeFaceMeshes();
    if(!isDone)
    {
        cout<<"MesherClass::PrismaticLayers_generatePrismaticMesh()->____error in merging the prismatic face mesh data sources____"<<endl;
        mr.isDone = false;
        mr.message = QString("Boundary mesh generation cannot be started: error in merging the face mesh data sources");
        return mr;
    }

    //! ------------------------------------------------------------
    //! retrieve and set the prismatic mesh parameters for the body
    //! set the parameters for the prismatic layers generation
    //! ------------------------------------------------------------
    const prismaticLayerParameters &parameters = myMeshDB->prismaticMeshParameters.value(bodyIndex);
    prismaticLayerBuilder.setParameters(parameters);

    //! -----------------------------
    //! tet boundary layer generator
    //! -----------------------------
    occHandle(Ng_MeshVS_DataSource3D) meshAtWalls;
    occHandle(Ng_MeshVS_DataSourceFace) lastInflatedMeshDS;
    prismaticLayerBuilder.generateTetLayers(meshAtWalls,lastInflatedMeshDS);
    myMeshDB->ArrayOfMeshDS.insert(bodyIndex,meshAtWalls);

    //! -----------------------
    //! enable the Stop button
    //! -----------------------
    if(myProgressIndicator!=Q_NULLPTR) emit requestEnableStop();

    mr.isDone = true;
    mr.message = QString("The prismatic mesh has been generated");
    return mr;
}
*/

//! ------------------------------------------------
//! function: PrismaticLayers_generatePrismaticMesh
//! details:
//! ------------------------------------------------
userMessage MesherClass::PrismaticLayers_generatePrismaticMesh(int bodyIndex,
                                                               occHandle(Ng_MeshVS_DataSourceFace) &theLastInflatedMesh,
                                                               QList<occHandle(Ng_MeshVS_DataSourceFace)> &listOfInflatedMeshes)
{    
    cout<<"MesherClass::PrismaticLayers_generatePrismaticMesh()->____function called____"<<endl;

    //! ---------------
    //! meshing result
    //! ---------------
    userMessage mr;

    //! ------------------------
    //! disable the Stop button
    //! ------------------------
    if(myProgressIndicator!=Q_NULLPTR) emit requestDisableStop();

    //! ---------------------------------------
    //! initialize the prismatic layer builder
    //! ---------------------------------------
    prismaticLayer prismaticLayerBuilder(myMeshDB);
    if(myProgressIndicator!=Q_NULLPTR) prismaticLayerBuilder.setProgressIndicator(myProgressIndicator);
    prismaticLayerBuilder.setBody(bodyIndex);

    //! --------------------------------------------------------------------------
    //!  retrieve the list of the prismatic faces for the body (list of face IDs)
    //! --------------------------------------------------------------------------
    const QList<int> &prismaticFacesOnBody = myMeshDB->prismaticFaces.value(bodyIndex);
    std::vector<int> prismaticFaces = prismaticFacesOnBody.toVector().toStdVector();
    prismaticLayerBuilder.setPrismaticFaces(prismaticFaces);

    //! ---------------------------------------------------------------
    //! merge the prismatic face meshes into a unique mesh data source
    //! ---------------------------------------------------------------
    bool isDone = prismaticLayerBuilder.mergeFaceMeshes();
    if(!isDone)
    {
        cout<<"MesherClass::PrismaticLayers_generatePrismaticMesh()->____error in merging the prismatic face mesh data sources____"<<endl;
        mr.isDone = false;
        mr.message = QString("Inflation cannot start: error in merging the face mesh data sources");
        return mr;
    }

    //! ------------------------------------------------------------
    //! retrieve and set the prismatic mesh parameters for the body
    //! set the parameters for the prismatic layers generation
    //! ------------------------------------------------------------
    const prismaticLayerParameters &parameters = myMeshDB->prismaticMeshParameters.value(bodyIndex);
    prismaticLayerBuilder.setParameters(parameters);

    occHandle(Ng_MeshVS_DataSource3D) prismaticMeshDS;
    isDone = prismaticLayerBuilder.generateMeshAtWalls(prismaticMeshDS,theLastInflatedMesh,listOfInflatedMeshes);

    //! -----------------------------------------------------
    //! the "prismaticMeshDS" is inserted in "ArrayOfMeshDS"
    //! as a volume mesh
    //! -----------------------------------------------------
    if(isDone) myMeshDB->ArrayOfMeshDS.insert(bodyIndex,prismaticMeshDS);
    else myMeshDB->ArrayOfMeshDS.insert(bodyIndex,occHandle(MeshVS_DataSource)());

    //! -----------------------
    //! enable the Stop button
    //! -----------------------
    if(myProgressIndicator!=Q_NULLPTR) emit requestEnableStop();

    mr.isDone = true;
    mr.message = QString("The prismatic mesh has been generated");
    return mr;
}

//! ------------------------------------
//! function: Netgen_generateVolumeMesh
//! details:
//! ------------------------------------
userMessage MesherClass::Netgen_generateVolumeMesh(int bodyIndex,
                                                   occHandle(MeshVS_DataSource) &mainMesh2D,
                                                   occHandle(MeshVS_DataSource) &mainMesh3D)
{
    userMessage mr;
    bool isInflationDefinedOnBody = myMeshDB->HasPrismaticFaces.value(bodyIndex);
    if(!isInflationDefinedOnBody)
    {
        //! ------------------------
        //! Netgen as volume mesher
        //! ------------------------
        myNetgenMesher->setMeshDataBase(myMeshDB);
        bool isVolume = true;
        userMessage mr = myNetgenMesher->performBrep(bodyIndex,isVolume,mainMesh2D,mainMesh3D);
        if(!mr.isDone) return mr;

        //! ------------------------------------------
        //! record the volume and the surface mesh
        //! here the face mesh data sources have been
        //! already put into the mesh data base
        //! ------------------------------------------
        myMeshDB->ArrayOfMeshDS2D.insert(bodyIndex,mainMesh2D);
        myMeshDB->ArrayOfMeshDS.insert(bodyIndex,mainMesh3D);
        return mr;
    }
    else
    {
        //! -------------------------
        //! check algorithm pre/post
        //! -------------------------
        const prismaticLayerParameters &parameters = myMeshDB->prismaticMeshParameters.value(bodyIndex);
        int algo = parameters.generationAlgorithm;

        if(algo ==0)
        {
            //! --------------------------
            //! pre-inflation algorithm
            //! generate the surface mesh
            //! --------------------------
            mr = this->Netgen_generateSurfaceMesh(bodyIndex,mainMesh2D);
            if(!mr.isDone)
            {
                mr.message = QString("Error in generating the surface mesh");
                return mr;
            }
            //! -----------------------------------------------
            //! generate the prismatic elements: the function
            //! PrismaticLayers_generatePrismaticMesh puts the
            //! boundary mesh data source into the database
            //! -----------------------------------------------
            QList<occHandle(Ng_MeshVS_DataSourceFace)> listOfInflatedMeshes;
            occHandle(Ng_MeshVS_DataSourceFace) theLastInflatedMesh;
            mr = this->PrismaticLayers_generatePrismaticMesh(bodyIndex,theLastInflatedMesh,listOfInflatedMeshes);

            if(mr.isDone==false)
            {
                mr.isDone = false;
                mr.message = QString("Error in generating the prismatic meshes");
                return mr;
            }

            //! ----------------------------------------------------
            //! retrieve the prismatic mesh from the mesh data base
            //! ----------------------------------------------------
            occHandle(Ng_MeshVS_DataSource3D) prismatic3DMesh = occHandle(Ng_MeshVS_DataSource3D)::DownCast(myMeshDB->ArrayOfMeshDS.value(bodyIndex));
            occHandle(Ng_MeshVS_DataSource3D) postInflationVolumeMesh;

            //! --------------------------------------------------------------
            //! generate the volume mesh from the post inflation surface mesh
            //! --------------------------------------------------------------
            mr = myNetgenMesher->meshVolumeFromSurfaceMesh(bodyIndex,theLastInflatedMesh,postInflationVolumeMesh);
            if(mr.isDone == false)
            {
                mr.isDone = false;
                mr.message = QString("Error in generating the post inflation volume mesh");
                return mr;
            }

            //! ---------------------------------------------------------
            //! record the post inflation mesh - for diagnostic purposes
            //! ---------------------------------------------------------
            //myMeshDB->ArrayOfMeshDS.insert(bodyIndex,postInflationVolumeMesh);

            //! -------------------------------------------------------------
            //! sum the prismatic mesh and the volume mesh
            //! it uses the 3D mesh constructor from a list of meshElementByCoords
            //! which perform the automatic node renumbering
            //! -------------------------------------------------------------
            std::vector<meshElementByCoords> elementList1, elementList2;
            MeshTools::toListOf3DElements(prismatic3DMesh,elementList1);
            MeshTools::toListOf3DElements(postInflationVolumeMesh,elementList2);
            elementList1.reserve(elementList1.size()+elementList2.size());
            elementList1.insert(elementList1.end(),elementList2.begin(),elementList2.end());
            mainMesh3D = new Ng_MeshVS_DataSource3D(elementList1);

            //! -----------------------------------------
            //! use the surface mesh built topologically
            //! -----------------------------------------
            if(!mainMesh3D.IsNull())
            {
                const occHandle(Ng_MeshVS_DataSource3D) &curMeshDS3D = occHandle(Ng_MeshVS_DataSource3D)::DownCast(mainMesh3D);
                if(curMeshDS3D->myFaceToElements.isEmpty()) curMeshDS3D->buildFaceToElementConnectivity();
                mainMesh2D = new Ng_MeshVS_DataSource2D(curMeshDS3D);
            }

            //! -----------------------------------------------------
            //! insert the main surface and volume mesh data sources
            //! -----------------------------------------------------
            myMeshDB->ArrayOfMeshDS2D.insert(bodyIndex, mainMesh2D);
            myMeshDB->ArrayOfMeshDS.insert(bodyIndex, mainMesh3D);

            //! ----------------------------------------------------------------
            //! add the face mesh datasources
            //! The nodes of the face mesh data source nodes must be renumbered
            //! using the new globalNodeIDs. To this aim:
            //! - generate a map (pointHash, globalNodeID) using the surface or
            //!   volume mesh
            //! - retrieve the "old" face mesh data source from the mesh data
            //!   base and use Ng_MeshVS_DataSourceFace::renumberNodes()
            //! ----------------------------------------------------------------
            std::map<size_t,int> mapNodeCoords_globalNodeID;
            occHandle(Ng_MeshVS_DataSource2D) curSurfaceMeshDS = occHandle(Ng_MeshVS_DataSource2D)::DownCast(myMeshDB->ArrayOfMeshDS2D.value(bodyIndex));
            const std::vector<mesh::tolerantPoint> &vecTolerantPoints = curSurfaceMeshDS->buildTolerantPoints(1e-6);
            std::vector<int> vecGlobalNodeIDs;

            for(TColStd_MapIteratorOfPackedMapOfInteger it(curSurfaceMeshDS->GetAllNodes()); it.More(); it.Next())
            {
                int globalNodeID = it.Key();
                vecGlobalNodeIDs.push_back(it.Key());
                int localNodeID = curSurfaceMeshDS->myNodesMap.FindIndex(globalNodeID);
                const std::vector<double> &aPoint = curSurfaceMeshDS->getNodeCoordinates(localNodeID);

                size_t seed = 0;
                for(int j=0; j<3; j++) hash_c<double>(seed, aPoint[j]);
                std::pair<size_t,int> aPair;
                aPair.first = seed;
                aPair.second = globalNodeID;
                mapNodeCoords_globalNodeID.insert(aPair);
            }

            //! ----------------------------------------
            //! iterate over the face mesh data sources
            //! ----------------------------------------
            int NbGeometryFaces = myMeshDB->MapOfBodyTopologyMap.value(bodyIndex).faceMap.Extent();
            for(int faceNr = 1; faceNr<=NbGeometryFaces; faceNr++)
            {
                //! -------------------
                //! renumber the nodes
                //! -------------------
                if(myMeshDB->ArrayOfMeshDSOnFaces.getValue(bodyIndex,faceNr).IsNull())
                {
                    cout<<"- the face mesh data source nr: "<<faceNr<<" is null. Jumping over it"<<endl;
                    continue;
                }
                occHandle(Ng_MeshVS_DataSourceFace) curFaceMeshDS = occHandle(Ng_MeshVS_DataSourceFace)::DownCast(myMeshDB->ArrayOfMeshDSOnFaces.getValue(bodyIndex,faceNr));
                curFaceMeshDS->renumberNodes(mapNodeCoords_globalNodeID,vecTolerantPoints,vecGlobalNodeIDs);

                //! -------------------------------------------------
                //! put the face mesh data sources into the database
                //! -------------------------------------------------
                myMeshDB->ArrayOfMeshDSOnFaces.setValue(bodyIndex,faceNr,curFaceMeshDS);
            }

            mr.isDone = true;
            mr.message = QString("Netgen volume mesh generated");
            return mr;
        }
        else
        {
            //! -------------------------
            //! post-inflation algorithm
            //! Netgen as volume mesher
            //! -------------------------
            myNetgenMesher->setMeshDataBase(myMeshDB);
            bool isVolume = true;
            userMessage mr = myNetgenMesher->performBrep(bodyIndex,isVolume,mainMesh2D,mainMesh3D);
            if(mr.isDone)
            {
                myMeshDB->ArrayOfMeshDS2D.insert(bodyIndex,mainMesh2D);
                myMeshDB->ArrayOfMeshDS.insert(bodyIndex,mainMesh3D);
            }

            //! ----------------------------------------------------
            //! build the prismatic mesh elements. The function
            //! PrismaticLayers_generatePrismaticMesh_post puts the
            //! boundary mesh data source into the database
            //! ----------------------------------------------------
            occHandle(Ng_MeshVS_DataSource3D) prismaticMeshDS;
            mr = this->PrismaticLayers_generatePrismaticMesh_post(bodyIndex,prismaticMeshDS);
            if(mr.isDone==false)
            {
                mr.isDone = false;
                mr.message = QString("Error in generating the prismatic mesh");
                return mr;
            }

            //! -------------------------------------------------------------------
            //! sum the prismatic mesh and the volume mesh
            //! it uses the 3D mesh constructor from a list of meshElementByCoords
            //! which perform an automatic node renumbering
            //! -------------------------------------------------------------------
            std::vector<meshElementByCoords> elementList1, elementList2;
            MeshTools::toListOf3DElements(prismaticMeshDS,elementList1);
            MeshTools::toListOf3DElements(occHandle(Ng_MeshVS_DataSource3D)::DownCast(mainMesh3D),elementList2);
            elementList1.reserve(elementList1.size()+elementList2.size());
            elementList1.insert(elementList1.end(),elementList2.begin(),elementList2.end());

            occHandle(Ng_MeshVS_DataSource3D) summedMeshDS = new Ng_MeshVS_DataSource3D(elementList1);
            summedMeshDS->buildFaceToElementConnectivity();

            mainMesh3D = summedMeshDS;
            myMeshDB->ArrayOfMeshDS.insert(bodyIndex,mainMesh3D);

            //! ---------------------------------
            //! surface mesh built topologically
            //! ---------------------------------
            mainMesh2D = new Ng_MeshVS_DataSource2D(summedMeshDS);
            myMeshDB->ArrayOfMeshDS2D.insert(bodyIndex,mainMesh2D);

            mr.isDone = true;
            mr.message = QString("The prismatic mesh has been succesfully generated");
            return mr;
        }
    }
}

//! ------------------------------------
//! function: Tetgen_generateVolumeMesh
//! details:
//! ------------------------------------
userMessage MesherClass::Tetgen_generateVolumeMesh(int bodyIndex, int preserveSurfaceMesh, bool runInMemory, int done)
{
    cout<<"MesherClass::Tetgen_generateVolumeMesh()->____function called____"<<endl;

    userMessage mr;
    bool hasPrismaticFaces = myMeshDB->HasPrismaticFaces.value(bodyIndex);

    if(!hasPrismaticFaces)  // standard meshing without boundary mesh
    {
        //! ---------------------------------------------------
        //! start meshing volume: initialize the tetgen mesher
        //! ---------------------------------------------------
        if(runInMemory==false)
        {
            cout<<"MesherClass::Tetgen_generateVolumeMesh()->____Tetgen will run on disk____"<<endl;

            //! -----------------------------------------------------------------
            //! arrayOfMeshDS is used by the Tetgen mesher for building the PLC
            //! the face mesh data sources are retrieved from the mesh database,
            //! so they should have been generated in advance
            //! ------------------------------------------------------------------
            int NbGeometryFaces = myMeshDB->MapOfBodyTopologyMap.value(bodyIndex).faceMap.Extent();
            NCollection_Array1<occHandle(Ng_MeshVS_DataSourceFace)> inputArrayOfFaceMeshDS(0,NbGeometryFaces);
            for(int faceNr=0; faceNr<=NbGeometryFaces; faceNr++)
            {
                occHandle(Ng_MeshVS_DataSourceFace) aFaceMeshDS = occHandle(Ng_MeshVS_DataSourceFace)::DownCast(myMeshDB->ArrayOfMeshDSOnFaces.getValue(bodyIndex,faceNr));
                inputArrayOfFaceMeshDS.SetValue(faceNr,aFaceMeshDS);
            }

            //! -------------
            //! set switches
            //! -------------
            myTetMesher->setSwitches(preserveSurfaceMesh, bodyIndex);

            //! -----------------------------------------------------------------
            //! tetgenArrayOfFaceMeshDS, tetgenVolumeMeshDS, tetgenSurfaceMeshDS
            //! are the returning data
            //! -----------------------------------------------------------------
            NCollection_Array1<occHandle(Ng_MeshVS_DataSourceFace)> tetgenArrayOfFaceMeshDS(0,NbGeometryFaces);
            occHandle(Ng_MeshVS_DataSource3D) tetgenVolumeMeshDS;
            occHandle(Ng_MeshVS_DataSource2D) tetgenSurfaceMeshDS;

            //! -----------------------------------------------------------------------
            //! this will generate <bodyName>1.ele, <bodyName>1.face, <bodyName>1.edge
            //! the file name is internally found by tetMesher from the mesh database
            //! -----------------------------------------------------------------------
            int exitCode = myTetMesher->performOnDisk1(inputArrayOfFaceMeshDS,bodyIndex,tetgenVolumeMeshDS,tetgenSurfaceMeshDS,tetgenArrayOfFaceMeshDS,done);

            if(exitCode==0)
            {
                //! ------------
                //! "0" success
                //! ------------
                myMeshDB->ArrayOfMeshDS.insert(bodyIndex,tetgenVolumeMeshDS);
                myMeshDB->ArrayOfMeshDS2D.insert(bodyIndex,tetgenSurfaceMeshDS);

                //! ------------------------------------------------
                //! put the face mesh datasources into the database
                //! ------------------------------------------------
                occHandle(Ng_MeshVS_DataSourceFace) aFaceMeshDS;
                for(int faceNr = 0; faceNr<=NbGeometryFaces; faceNr++)
                {
                    aFaceMeshDS = tetgenArrayOfFaceMeshDS.Value(faceNr);
                    myMeshDB->ArrayOfMeshDSOnFaces.setValue(bodyIndex,faceNr,aFaceMeshDS);
                }

                mr.isDone = true;
                mr.message = QString("Tetgen has successfully generated the volume mesh");
                return mr;
            }
            else if(exitCode==-1)
            {
                mr.isDone = false;
                mr.message = QString("Tetgen failure. Error in building PLC for body index: %1").arg(bodyIndex);
                cout<<mr.message.toStdString()<<endl;
                return mr;
            }
            else if(exitCode==-2)
            {
                mr.isDone = false;
                mr.message = QString("Tetgen stopped while generating the face mesh data sources").arg(exitCode);
                cout<<mr.message.toStdString()<<endl;
                return mr;
            }
            else
            {
                mr.isDone = false;
                mr.message = QString("Tetgen failure. Tetgen exit code: %1____").arg(exitCode);
                cout<<mr.message.toStdString()<<endl;
                return mr;
            }
        }
        else    // use tetgen library
        {
            cout<<"MesherClass::Tetgen_generateVolumeMesh()->____Tetgen will run in memory____"<<endl;

            //! --------------------------------------
            //! create the PLC using the surface mesh
            //! --------------------------------------
            QList<int> invalidFaceTags;
            bool isPLCValid = myTetMesher->buildPLC(bodyIndex,invalidFaceTags);
            if(isPLCValid==false)
            {
                mr.isDone = false;
                mr.message = QString("Tetgen failure. Error in building PLC for body index: %1").arg(bodyIndex);
                cout<<mr.message.toStdString()<<endl;
                return mr;
            }
            //! ---------------------------
            //! set the meshing parameters
            //! ---------------------------
            myTetMesher->setSwitches(preserveSurfaceMesh, bodyIndex);

            tetgenio meshOut;
            bool isVolDone = myTetMesher->perform(&meshOut, bodyIndex, false);
            if(isVolDone==false)
            {
                mr.isDone = false;
                mr.message = QString("Tetgen meshing failure on body index: %1").arg(bodyIndex);
                cout<<mr.message.toStdString()<<endl;
                return mr;
            }

            cout<<"MesherClass::generateMesh()->____Tetgen mesher: building the tetgen mesh datasources____"<<endl;

            //! -------------------------
            //! volume mesh data sources
            //! -------------------------
            occHandle(Ng_MeshVS_DataSource3D) aMeshDS = new Ng_MeshVS_DataSource3D(meshOut);
            myMeshDB->ArrayOfMeshDS.insert(bodyIndex,aMeshDS);
            cout<<"MesherClass::generateMesh()->____volume mesh OK____"<<endl;

            //! -------------------------
            //! surface mesh data source
            //! -------------------------
            aMeshDS->buildFaceToElementConnectivity();
            occHandle(Ng_MeshVS_DataSource2D) aMeshDS2D = new Ng_MeshVS_DataSource2D(aMeshDS);

            //occHandle(Ng_MeshVS_DataSource2D) aMeshDS2D = new Ng_MeshVS_DataSource2D(meshOut);
            //myMeshDB->ArrayOfMeshDS2D.insert(bodyIndex,aMeshDS2D);
            //cout<<"MesherClass::generateMesh()->____surface mesh OK____"<<endl;

            //! -------------------------------
            //! generate the face data sources
            //! -------------------------------
            int NbGeometryFaces = myMeshDB->MapOfBodyTopologyMap.value(bodyIndex).faceMap.Extent();

            //! --------------------------------
            //! init the secondary progress bar
            //! --------------------------------
            if(myProgressIndicator!=Q_NULLPTR)
            {
                QProgressEvent *pe = new QProgressEvent();
                pe->setVal(done);
                pe->setMessage("Generating the face mesh data sources");
                pe->setAction1(QProgressEventAction::QProgressEvent_Init);
                pe->setRange1(0, NbGeometryFaces);
                QApplication::postEvent(myProgressIndicator,pe);
            }

            cout<<"MesherClass::generateMesh()->____start generating face mesh data sources____"<<endl;
            for(int faceNr=1; faceNr<=NbGeometryFaces;faceNr++)
            {
                occHandle(Ng_MeshVS_DataSourceFace) aFaceMeshDS = new Ng_MeshVS_DataSourceFace(&meshOut,faceNr);
                myMeshDB->ArrayOfMeshDSOnFaces.setValue(bodyIndex,faceNr,aFaceMeshDS);

                //! -------------------------------------------
                //! post update events on face mesh generation
                //! -------------------------------------------
                if(myProgressIndicator!=Q_NULLPTR && faceNr%5==0)
                {
                    QProgressEvent *pe = new QProgressEvent();
                    pe->setVal(done);
                    pe->setMessage(QString("Face mesh %1").arg(faceNr));
                    pe->setVal1(faceNr);
                    QApplication::postEvent(myProgressIndicator,pe);
                }
            }

            mr.isDone = true;
            mr.message = QString("MesherClass::Tetgen_generateVolumeMesh()->____Tetgen has successfully generated the volume mesh____");
            return mr;
        }
    }
    else    // meshing with boundary mesh
    {
        //! **********************************************************
        //!
        //! the block of code for the generation of the boundary mesh
        //! is in PrismaticLayers_generatePrismaticMesh
        //!
        //! **********************************************************
        occHandle(Ng_MeshVS_DataSourceFace) theLastInflatedMesh;
        QList<occHandle(Ng_MeshVS_DataSourceFace)> listOfInflatedMeshes;
        mr = this->PrismaticLayers_generatePrismaticMesh(bodyIndex,theLastInflatedMesh,listOfInflatedMeshes);
        if(mr.isDone==false)
        {
            mr.isDone = false;
            mr.message = "Error in generating the prismatic 3D mesh";
            return mr;
        }

        //! ----------------------------------------------------------------
        //! the volume mesh generated starting from the last displaced mesh
        //! ----------------------------------------------------------------
        occHandle(Ng_MeshVS_DataSource3D) tetgenVolumeMeshDS;

        if(runInMemory==false)  // tetgen run on disk
        {
            cout<<"____Tetgen will run on disk____"<<endl;

            //! ---------------------------------------------------
            //! start meshing volume: initialize the tetgen mesher
            //! ---------------------------------------------------
            NCollection_Array1<occHandle(Ng_MeshVS_DataSourceFace)> inputArrayOfFaceMeshDS(0,0);
            inputArrayOfFaceMeshDS.SetValue(0,theLastInflatedMesh);

            //! -------------
            //! set switches
            //! -------------
            preserveSurfaceMesh = 1;
            myTetMesher->setSwitches(preserveSurfaceMesh, bodyIndex);

            //! -----------------------------------------------------------------------
            //! this will generate <bodyName>1.ele, <bodyName>1.face, <bodyName>1.edge
            //! the file name is internally found by tetMesher from the mesh database
            //! -----------------------------------------------------------------------
            int exitCode = myTetMesher->performOnDisk1(inputArrayOfFaceMeshDS, bodyIndex, tetgenVolumeMeshDS);
            if(exitCode!=0)
            {
                mr.isDone = false;
                mr.message = QString("Tetgen volume mesh failure %1").arg(exitCode);
                cout<<mr.message.toStdString()<<endl;
                return mr;
            }

            //! ------------
            //! "0" success
            //! ------------
            // record the post inflation mesh: for testing purposes
            //myMeshDB->ArrayOfMeshDS.insert(bodyIndex,tetgenVolumeMeshDS);
        }
        else  // tetgen run in memory
        {
            //! --------------------------------------
            //! create the PLC using the surface mesh
            //! --------------------------------------
            //QList<int> invalidFaceTags;
            //bool isPLCValid = myTetMesher->buildPLC(bodyIndex,invalidFaceTags);
            bool isPLCValid = myTetMesher->buildPLC(theLastInflatedMesh);
            if(isPLCValid==false)
            {
                mr.isDone = false;
                mr.message = QString("Tetgen failure. Error in building PLC for body index: %1").arg(bodyIndex);
                cout<<mr.message.toStdString()<<endl;
                return mr;
            }

            //! ---------------------------
            //! set the meshing parameters
            //! ---------------------------
            int preserveSurfaceMesh = true;
            myTetMesher->setSwitches(preserveSurfaceMesh, bodyIndex);

            tetgenio meshOut;
            bool isVolDone = myTetMesher->perform(&meshOut, bodyIndex, false);
            if(isVolDone==false)
            {
                mr.isDone = false;
                mr.message = QString("Tetgen meshing failure on body index: %1").arg(bodyIndex);
                cout<<mr.message.toStdString()<<endl;
                return mr;
            }
            tetgenVolumeMeshDS = new Ng_MeshVS_DataSource3D(meshOut);
        }

        //! --------------------------
        //! volume mesh not generated
        //! --------------------------
        if(tetgenVolumeMeshDS.IsNull())
        {
            mr.isDone = false;
            mr.message = QString("Cannot merge meshes for body %1. The interior mesh is null").arg(bodyIndex);
            return mr;
        }

        //! -------------------------------------------------
        //! merge the volume meshes
        //! retrieve the prismatic 3D mesh from the database
        //! -------------------------------------------------
        occHandle(Ng_MeshVS_DataSource3D) prismatic3DMesh =
                occHandle(Ng_MeshVS_DataSource3D)::DownCast(myMeshDB->ArrayOfMeshDS.value(bodyIndex));

        std::vector<meshElementByCoords> elementList1, elementList2;
        MeshTools::toListOf3DElements(prismatic3DMesh,elementList1);
        MeshTools::toListOf3DElements(tetgenVolumeMeshDS,elementList2);
        elementList1.reserve(elementList1.size()+elementList2.size());
        elementList1.insert(elementList1.end(),elementList2.begin(),elementList2.end());
        occHandle(Ng_MeshVS_DataSource3D) mainMesh3D = new Ng_MeshVS_DataSource3D(elementList1,true,true);
        if(mainMesh3D.IsNull())
        {
            mr.isDone = false;
            mr.message = QString("Error in merging the prismatic and interionr mesh for body %1").arg(bodyIndex);
            return mr;
        }
        //myMeshDB->ArrayOfMeshDS.insert(bodyIndex,mainMesh3D);

        //! for testing purposes - show the interior mesh
        myMeshDB->ArrayOfMeshDS.insert(bodyIndex,tetgenVolumeMeshDS);

        mr.isDone = true;
        mr.message = QString("Prismatic 3D mesh and interior mesh for body %1 successfully merged").arg(bodyIndex);
        return mr;
    }
}

//! ----------------------------------------
//! function: Netgen_STL_generateVolumeMesh
//! details:
//! ----------------------------------------
userMessage MesherClass::Netgen_STL_generateVolumeMesh(int bodyIndex,
                                                       occHandle(MeshVS_DataSource) &mainMesh2D,
                                                       occHandle(MeshVS_DataSource) &mainMesh3D)
{
    cout<<"MesherClass::Netgen_STL_generateVolumeMesh()->____function called____"<<endl;

    //! -------------------
    //! the meshing result
    //! -------------------
    userMessage mr;

    //! ---------------------------------
    //! do preliminary work on .stl mesh
    //! ---------------------------------
    QString inputExtendedSTL_S = this->processSTL(bodyIndex);

    //! ----------------------------------------------------------------------------------
    //! read the extended .stl file and return the surface mesh and face mesh datasources
    //! ----------------------------------------------------------------------------------
    occHandle(Ng_MeshVS_DataSource2D) anSTL_MeshVS_DataSource;
    int NbGeometryFaces = myMeshDB->MapOfBodyTopologyMap.value(bodyIndex).faceMap.Extent();

    std::vector<occHandle(Ng_MeshVS_DataSourceFace)> array_STL_meshDS;
    for(int i=0; i<=NbGeometryFaces; i++) array_STL_meshDS.push_back(occHandle(Ng_MeshVS_DataSourceFace)());
    const TopoDS_Shape &shape = myMeshDB->bodyMap.value(bodyIndex);
    bool isDone = MeshTools::toSTLMesh1(shape, inputExtendedSTL_S, anSTL_MeshVS_DataSource, array_STL_meshDS);

    if(!isDone)
    {
        mr.isDone = false;
        mr.message = "Error in writing the .stl file";
        return mr;
    }

    //! -------------------------------------------------------
    //! put the .stl face mesh data sources into the database
    //! (since the prismatic layer builder reads the mesh data
    //! sources from the mesh data base)
    //! -------------------------------------------------------
    mainMesh2D = anSTL_MeshVS_DataSource;
    myMeshDB->ArrayOfMeshDS2D.insert(bodyIndex,mainMesh2D);
    for(int faceNr = 0; faceNr<=NbGeometryFaces; faceNr++)
    {
        if(!array_STL_meshDS.at(faceNr).IsNull())
            myMeshDB->ArrayOfMeshDSOnFaces.setValue(bodyIndex,faceNr,array_STL_meshDS.at(faceNr));
    }

    //! ------------------------------------------
    //! check if inflation is defined on the body
    //! ------------------------------------------
    bool isInflationDefinedOnBody = myMeshDB->HasPrismaticFaces.value(bodyIndex);
    if(!isInflationDefinedOnBody)
    {
        //! --------------------------------------------
        //! the main call for generating the volume mesh
        //! --------------------------------------------
        myNetgenMesher->setMeshDataBase(myMeshDB);
        bool isVolume = true;
        mr = myNetgenMesher->perform(bodyIndex,isVolume,mainMesh2D,mainMesh3D);

        //! ----------------------------------------------------------
        //! put the surface and volume mesh into the mesh data source
        //! (NULL or valid)
        //! ----------------------------------------------------------
        myMeshDB->ArrayOfMeshDS2D.insert(bodyIndex,mainMesh2D);
        myMeshDB->ArrayOfMeshDS.insert(bodyIndex,mainMesh3D);

        //! -----------------------------------------------
        //! clear the .stl mesh data sources from the view
        //! after a meshing process failure
        //! -----------------------------------------------
        if(!mr.isDone)
        {
            for(int faceNr = 0; faceNr<=NbGeometryFaces; faceNr++)
                myMeshDB->ArrayOfMeshDSOnFaces.setValue(bodyIndex,faceNr,occHandle(MeshVS_DataSource)());
            mr.isDone = false;
            mr.message = "Netgen error in generating the volume mesh from STL input";
            return mr;
        }

        //! ---------------------------------------------------------
        //! put the face mesh data sources into the data base
        //! OK - done: the data sources are set directly by the tool
        //! ---------------------------------------------------------

        //! ----------------------
        //! remove supports files
        //! ----------------------
        this->removeSupportFiles();
        return mr;
    }
    else
    {
        //! -------------------------
        //! check algorithm pre/post
        //! -------------------------
        const prismaticLayerParameters &parameters = myMeshDB->prismaticMeshParameters.value(bodyIndex);
        int algo = parameters.generationAlgorithm;
        if(algo==0)
        {
            //! --------------------------------
            //! generate the prismatic elements
            //! --------------------------------
            occHandle(Ng_MeshVS_DataSourceFace) theLastInflatedMesh;
            mr = this->Netgen_STL_generateSurfaceMesh(bodyIndex,mainMesh2D);
            if(!mr.isDone)
            {
                cout<<"MesherClass::Netgen_generateVolumeMesh()->____error in generating the surface mesh____"<<endl;
                mr.message = QString("Error in generating the surface mesh before inflation");
                return mr;
            }
            QList<occHandle(Ng_MeshVS_DataSourceFace)> listOfInflatedMeshes;
            mr = this->PrismaticLayers_generatePrismaticMesh(bodyIndex,theLastInflatedMesh,listOfInflatedMeshes);
            if(mr.isDone==false) return mr;

            occHandle(Ng_MeshVS_DataSource3D) prismatic3DMesh = occHandle(Ng_MeshVS_DataSource3D)::DownCast(myMeshDB->ArrayOfMeshDS.value(bodyIndex));
            occHandle(Ng_MeshVS_DataSource3D) postInflationVolumeMesh;

            mr = myNetgenMesher->meshVolumeFromSurfaceMesh(bodyIndex,theLastInflatedMesh,postInflationVolumeMesh);
            if(mr.isDone == false)
            {
                mr.isDone = false;
                mr.message = QString("Error in generating the post inflation volume mesh");
                return mr;
            }

            //! ---------------------------------------------------------
            //! record the post inflation mesh - for diagnostic purposes
            //! ---------------------------------------------------------
            //myMeshDB->ArrayOfMeshDS.insert(bodyIndex,postInflationVolumeMesh);

            //! -------------------------------------------------------------
            //! sum the prismatic mesh and the volume mesh
            //! it uses the 3D mesh constructor from a list of meshElementByCoords
            //! which perform an automatic node renumbering
            //! -------------------------------------------------------------
            std::vector<meshElementByCoords> elementList1, elementList2;
            MeshTools::toListOf3DElements(prismatic3DMesh,elementList1);
            MeshTools::toListOf3DElements(postInflationVolumeMesh,elementList2);
            elementList1.reserve(elementList1.size()+elementList2.size());
            elementList1.insert(elementList1.end(),elementList2.begin(),elementList2.end());
            mainMesh3D = new Ng_MeshVS_DataSource3D(elementList1);

            //! -----------------------------------------
            //! use the surface mesh built topologically
            //! -----------------------------------------
            if(!mainMesh3D.IsNull())
            {
                const occHandle(Ng_MeshVS_DataSource3D) &curMeshDS3D = occHandle(Ng_MeshVS_DataSource3D)::DownCast(mainMesh3D);
                if(curMeshDS3D->myFaceToElements.isEmpty()) curMeshDS3D->buildFaceToElementConnectivity();
                mainMesh2D = new Ng_MeshVS_DataSource2D(curMeshDS3D);
            }

            mr.isDone = true;
            mr.message = QString("Netgen volume mesh from stl geometry generated");

            myMeshDB->ArrayOfMeshDS2D.insert(bodyIndex, mainMesh2D);
            myMeshDB->ArrayOfMeshDS.insert(bodyIndex, mainMesh3D);

            cout<<"\\-------------------------------------------\\"<<endl;
            cout<<"\\- preparing the map of new global node IDs \\"<<endl;
            cout<<"\\-------------------------------------------\\"<<endl;

            //! ----------------------------------------------------------------
            //! put the face mesh data sources into the data base
            //! OK - done: the data sources are set directly by the tool
            //! The nodes of the face mesh data source nodes must be renumbered
            //! using the new globalNodeIDs. To this aim:
            //! - generate a map (pointHash, globalNodeID) using the surface or
            //!   volume mesh
            //! - retrieve the old mesh data source from the mesh data base and
            //!   use Ng_MeshVS_DataSourceFace::renumberNodes()
            //! ----------------------------------------------------------------
            std::map<size_t,int> mapNodeCoords_globalNodeID;
            occHandle(Ng_MeshVS_DataSource2D) curSurfaceMeshDS = occHandle(Ng_MeshVS_DataSource2D)::DownCast(myMeshDB->ArrayOfMeshDS2D.value(bodyIndex));
            std::vector<mesh::tolerantPoint> vecTolerantPoints = curSurfaceMeshDS->buildTolerantPoints(1e-6);
            std::vector<int> vecGlobalNodeIDs;

            for(TColStd_MapIteratorOfPackedMapOfInteger it(curSurfaceMeshDS->GetAllNodes()); it.More(); it.Next())
            {
                int globalNodeID = it.Key();
                vecGlobalNodeIDs.push_back(globalNodeID);
                int localNodeID = curSurfaceMeshDS->myNodesMap.FindIndex(globalNodeID);
                const std::vector<double> &aPoint = curSurfaceMeshDS->getNodeCoordinates(localNodeID);

                size_t seed = 0;
                for(int j=0; j<3; j++) hash_c<double>(seed, aPoint[j]);
                //std::pair<size_t,int> aPair;
                //aPair.first = seed;
                //aPair.second = globalNodeID;
                mapNodeCoords_globalNodeID.insert(std::make_pair(seed,globalNodeID));
            }

            cout<<"\\-----------------------------------------------------\\"<<endl;
            cout<<"\\- start nodal renumbering the face mesh data sources \\"<<endl;
            cout<<"\\-----------------------------------------------------\\"<<endl;

            //! ----------------------------------------
            //! iterate over the face mesh data sources
            //! ----------------------------------------
            for(int faceNr = 1; faceNr<=NbGeometryFaces; faceNr++)
            {
                cout<<"- renumbering face data source nr: "<<faceNr<<endl;

                //! -------------------
                //! renumber the nodes
                //! -------------------
                if(myMeshDB->ArrayOfMeshDSOnFaces.getValue(bodyIndex,faceNr).IsNull())
                {
                    cout<<"- the face mesh data source nr: "<<faceNr<<" is null. Jumping over it"<<endl;
                    continue;
                }
                occHandle(Ng_MeshVS_DataSourceFace) curFaceMeshDS = occHandle(Ng_MeshVS_DataSourceFace)::DownCast(myMeshDB->ArrayOfMeshDSOnFaces.getValue(bodyIndex,faceNr));
                curFaceMeshDS->renumberNodes(mapNodeCoords_globalNodeID,vecTolerantPoints,vecGlobalNodeIDs);

                //! -------------------------------------------------
                //! put the face mesh data sources into the database
                //! -------------------------------------------------
                myMeshDB->ArrayOfMeshDSOnFaces.setValue(bodyIndex,faceNr,curFaceMeshDS);
            }
            cout<<"\\-----------------------------------------------------\\"<<endl;
            cout<<"\\-         renumbering successfully done              \\"<<endl;
            cout<<"\\-----------------------------------------------------\\"<<endl;

            this->removeSupportFiles();
            return mr;
        }
        else
        {
            //! --------------------------------------------
            //! the main call for generating the volume mesh
            //! --------------------------------------------
            myNetgenMesher->setMeshDataBase(myMeshDB);
            bool isVolume = true;
            mr = myNetgenMesher->perform(bodyIndex,isVolume,mainMesh2D,mainMesh3D);

            //! ----------------------------------------------------------
            //! put the surface and volume mesh into the mesh data source
            //! (NULL or valid)
            //! ----------------------------------------------------------
            myMeshDB->ArrayOfMeshDS2D.insert(bodyIndex,mainMesh2D);
            myMeshDB->ArrayOfMeshDS.insert(bodyIndex,mainMesh3D);

            //! -----------------------------------------------
            //! clear the .stl mesh data sources from the view
            //! after a meshing process failure
            //! -----------------------------------------------
            if(!mr.isDone)
            {
                for(int faceNr = 0; faceNr<=NbGeometryFaces; faceNr++)
                    myMeshDB->ArrayOfMeshDSOnFaces.setValue(bodyIndex,faceNr,occHandle(MeshVS_DataSource)());
                mr.isDone = false;
                mr.message = "Netgen error in generating the volume mesh from STL input";
                return mr;
            }

            //! ---------------------------------------------------------
            //! put the face mesh data sources into the data base
            //! OK - done: the data sources are set directly by the tool
            //! ---------------------------------------------------------

            //! ----------------------
            //! remove supports files
            //! ----------------------
            this->removeSupportFiles();

            mr.isDone = true;
            mr.message = QString("Not supported");
            return mr;
        }
    }
}

//! -----------------------------
//! function: removeSupportFiles
//! details:
//! -----------------------------
void MesherClass::removeSupportFiles(bool remove)
{
    remove = true;
    if(remove)
    {
        //! ----------------------
        //! remove supports files
        //! ----------------------
        QString workingDir = tools::getWorkingDir();
        QString inputStandardSTL = workingDir+"/inputStandardSTLMesh.stl";
        QString inputExtendedSTL = workingDir+"/inputExtendedSTLMesh.stl";
        QString inputExtendedSTL_H = workingDir+"/inputExtendedSTLMesh_H.stl";
        QString inputExtendedSTL_S = workingDir+"/inputExtendedSTLMesh_S.stl";
        QString inputExtendedSTL_SS = workingDir+"/inputExtendedSTLMesh_SH.stl";

        QFile file;
        file.setFileName(inputStandardSTL);
        if(file.exists()) file.remove();
        file.setFileName(inputExtendedSTL);
        if(file.exists()) file.remove();
        file.setFileName(inputExtendedSTL_H);
        if(file.exists()) file.remove();
        file.setFileName(inputExtendedSTL_S);
        if(file.exists()) file.remove();
        file.setFileName(inputExtendedSTL_SS);
        if(file.exists()) file.remove();
    }
}

//! -------------------------------------------------
//! function: processTetWildSTL
//! details:  generate the triangle soup for TetWild
//! -------------------------------------------------
QString MesherClass::processTetWildSTL(int bodyIndex)
{
    //! --------------------------------------
    //! clean the tessellation from the shape
    //! --------------------------------------
    BRepTools::Clean(myMeshDB->bodyMap.value(bodyIndex));

    //! -------------------------------------------
    //! retrieve the parameters of the discretizer
    //! -------------------------------------------
    Property::meshEngine2D surfaceTessellator = myMeshDB->mapOfDiscretizer.value(bodyIndex);
    double angularDeflection = myMeshDB->mapOfAngularDeflection.value(bodyIndex);
    double linearDeflection = myMeshDB->mapOfLinearDeflection.value(bodyIndex);

    //! -----------------------------
    //! generate the "triangle soup"
    //! -----------------------------
    switch(surfaceTessellator)
    {
    case Property::meshEngine2D_OCC_STL:
    {
        bool isInParallel = true;
        bool isRelative = true;
        BRepMesh_IncrementalMesh stlOccMesher(myMeshDB->bodyMap.value(bodyIndex),linearDeflection,isRelative,angularDeflection,isInParallel);
        Q_UNUSED(stlOccMesher)
    }
        break;

    case Property::meshEngine2D_OCC_ExpressMesh:
    {
        occMesher anExpressMeshMesher(myMeshDB->bodyMap.value(bodyIndex));
        double minSize = myMeshDB->ArrayOfMinBodyElementSize.value(bodyIndex);
        double maxSize = myMeshDB->ArrayOfMaxBodyElementSize.value(bodyIndex);

        meshParam mp;
        mp.minElementSize = minSize;
        mp.maxElementSize = maxSize;

        bool isDone = anExpressMeshMesher.perform(angularDeflection, linearDeflection, minSize, maxSize);
        if(!isDone)
        {
            bool isInParallel = false;
            bool isRelative = true;
            BRepMesh_IncrementalMesh stlOccMesher(myMeshDB->bodyMap.value(bodyIndex),linearDeflection,isRelative,angularDeflection,isInParallel);
            Q_UNUSED(stlOccMesher)
        }
    }
        break;
    }

    QString inputStandardSTL = tools::getWorkingDir()+"/inputStandardSTLMesh.stl";
    const TopoDS_Shape &shape = myMeshDB->bodyMap.value(bodyIndex);
    STLAPIWriter aWriter;
    aWriter.Write(shape,inputStandardSTL.toStdString().c_str());
    return inputStandardSTL;
}

//! ---------------------
//! function: processSTL
//! details:
//! ---------------------
QString MesherClass::processSTL(int bodyIndex)
{
    cout<<"MesherClass::processSTL()->____function called____"<<endl;

    //! --------------------------
    //! disable the "Stop" button
    //! --------------------------
    if(myProgressIndicator!=Q_NULLPTR) emit requestDisableStop();

    //! ---------------------------------------------
    //! remove the previous support files if present
    //! ---------------------------------------------
    this->removeSupportFiles();

    //! -------------------------------------------
    //! retrieve the parameters of the discretizer
    //! -------------------------------------------
    Property::meshEngine2D surfaceDiscretizer = myMeshDB->mapOfDiscretizer.value(bodyIndex);

    double angularDeflection = myMeshDB->mapOfAngularDeflection.value(bodyIndex);
    double linearDeflection = myMeshDB->mapOfLinearDeflection.value(bodyIndex);
    bool isHealingOn = myMeshDB->MapOfIsBodyHealingOn.value(bodyIndex);
    bool isSimplificationOn = myMeshDB->MapOfIsMeshSimplificationOn.value(bodyIndex);
    int defeaturingBy = myMeshDB->MapOfMeshSimplificationBy.value(bodyIndex);
    double defeaturingParameterValue = myMeshDB->MapOfBodyDefeaturingParameterValue.value(bodyIndex);
    bool featurePreserve = myMeshDB->mapOfFeaturePreserving.value(bodyIndex);

    //! --------------------------------------
    //! clean the tessellation from the shape
    //! --------------------------------------
    BRepTools::Clean(myMeshDB->bodyMap.value(bodyIndex));

    //! --------------------------
    //! generate the tessellation
    //! --------------------------
    switch(surfaceDiscretizer)
    {
    case Property::meshEngine2D_OCC_STL:
    {
        cout<<"\\-------------------------------------------------------------------------------------------------------\\"<<endl;
        cout<<"\\ MesherClass::processSTL()->____using BRepIncrementalMesh mesh for generating the visualization mesh____"<<endl;
        cout<<"\\ linear deflection: "<<linearDeflection<<endl;
        cout<<"\\ angular deflection: "<<angularDeflection<<endl;
        cout<<"\\-------------------------------------------------------------------------------------------------------\\"<<endl;

        bool isInParallel = true;
        bool isRelative = true;
        BRepMesh_IncrementalMesh stlOccMesher(myMeshDB->bodyMap.value(bodyIndex),linearDeflection, isRelative, angularDeflection, isInParallel);

        Q_UNUSED(stlOccMesher)
    }
        break;

    case Property::meshEngine2D_OCC_ExpressMesh:
    {
        cout<<"\\--------------------------------------------------------------------------------------------\\"<<endl;
        cout<<"\\ MesherClass::processSTL()->____using Express mesh for generating the visualization mesh____"<<endl;
        cout<<"\\ linear deflection: "<<linearDeflection<<endl;
        cout<<"\\ angular deflection: "<<angularDeflection<<endl;
        cout<<"\\--------------------------------------------------------------------------------------------\\"<<endl;

        occMesher anExpressMeshMesher(myMeshDB->bodyMap.value(bodyIndex));

        //! --------------------------
        //! read also the body sizing
        //! --------------------------
        double minSize = myMeshDB->ArrayOfMinBodyElementSize.value(bodyIndex);
        double maxSize = myMeshDB->ArrayOfMaxBodyElementSize.value(bodyIndex);

        meshParam mp;
        mp.minElementSize = minSize;
        mp.maxElementSize = maxSize;

        cout<<"\\------------------------------------------"<<endl;
        cout<<"\\ linear deviation: "<<linearDeflection<<endl;
        cout<<"\\ angular deflection: "<<angularDeflection<<endl;
        cout<<"\\ max size: "<<maxSize<<endl;
        cout<<"\\ min size: "<<minSize<<endl;
        cout<<"\\------------------------------------------"<<endl;

        bool isDone = anExpressMeshMesher.perform(angularDeflection, linearDeflection, minSize, maxSize);
        if(!isDone)
        {
            cout<<"MesherClass::processSTL()->____Express mesh failure in generating the mesh for visualization. Use a standard .stl____"<<endl;

            bool isInParallel = true;
            bool isRelative = true;
            BRepMesh_IncrementalMesh stlOccMesher(myMeshDB->bodyMap.value(bodyIndex),linearDeflection, isRelative, angularDeflection,isInParallel);
            Q_UNUSED(stlOccMesher)
        }
    }
        break;
    }    

    //! -----------------------------------------------------
    //! write a standard .stl file - for diagnostic purposes
    //! write the extended .stl file
    //! -----------------------------------------------------
    QString inputStandardSTL = tools::getWorkingDir()+"/inputStandardSTLMesh.stl";
    QString inputExtendedSTL = tools::getWorkingDir()+"/inputExtendedSTLMesh.stl";

    const TopoDS_Shape &shape = myMeshDB->bodyMap.value(bodyIndex);
    STLAPIWriter aWriter;
    aWriter.Write(shape,inputStandardSTL.toStdString().c_str());
    aWriter.WriteExtended(shape,inputExtendedSTL.toStdString().c_str());

    if(isHealingOn==false && isSimplificationOn==false) return inputExtendedSTL;

    //! -------------------------------------------------
    //! prepare the data for the feature preserve option
    //! -------------------------------------------------
    std::vector<mesh::tolerantPoint> fixedPoints;
    if(featurePreserve==true)
    {
        occHandle(Ng_MeshVS_DataSource2D) surfaceMeshDS;
        std::vector<occHandle(Ng_MeshVS_DataSourceFace)> vecFaceMeshDS;
        int NbFaces = myMeshDB->MapOfBodyTopologyMap.value(bodyIndex).faceMap.Extent();
        for(int i=0; i<=NbFaces; i++) vecFaceMeshDS.push_back(occHandle(Ng_MeshVS_DataSourceFace)());
        //cout<<"____"<<NbFaces<<"____"<<endl;

        MeshTools::toSTLMesh1(shape,inputExtendedSTL,surfaceMeshDS,vecFaceMeshDS);

        QList<occHandle(Ng_MeshVS_DataSourceFace)> listOfFaces;
        for(int i=0; i<myMeshDB->featuredTags.size(); i++)
        {
            //cout<<"____feature tag: ("<<myMeshDB->featuredTags[i].parentShapeNr<<", "<<myMeshDB->featuredTags[i].subTopNr<<")____"<<endl;
            const GeometryTag &aTag = myMeshDB->featuredTags[i];
            if(aTag.parentShapeNr != bodyIndex) continue;
            occHandle(Ng_MeshVS_DataSourceFace) aFaceDS = vecFaceMeshDS.at(aTag.subTopNr);
            if(aFaceDS.IsNull()) continue;
            listOfFaces<<aFaceDS;
        }

        occHandle(Ng_MeshVS_DataSourceFace) BCMeshDS;
        if(listOfFaces.isEmpty()==false)
        {
            BCMeshDS = new Ng_MeshVS_DataSourceFace(listOfFaces);
            BCMeshDS->computeFreeMeshSegments();
            BCMeshDS->myBoundarySegments;
            cout<<"____number mesh points on BC boundaries: "<<BCMeshDS->myBoundaryPoints.length()<<"____"<<endl;

            for(QList<int>::iterator it = BCMeshDS->myBoundaryPoints.begin(); it!=BCMeshDS->myBoundaryPoints.end(); it++)
            {
                int globalNodeID = *it;
                int NbNodes;
                double c[3];
                TColStd_Array1OfReal coords(*c,1,3);
                MeshVS_EntityType aType;
                BCMeshDS->GetGeom(globalNodeID,false,coords,NbNodes,aType);
                //cout<<"____boundary point ("<<coords(1)<<", "<<coords(2)<<", "<<coords(3)<<")____"<<endl;
                double tolerance = 1e-3;
                fixedPoints.push_back(mesh::tolerantPoint(coords(1),coords(2),coords(3),tolerance));
            }
        }
    }

    //! ---------------------------------------------------
    //! a name for the corrected (ADMesh.exe) output file
    //! a name for the simplified (Simplify.h) output file
    //! ---------------------------------------------------
    QString inputExtendedSTL_H = tools::getWorkingDir()+"/inputExtendedSTLMesh_H.stl";
    QString inputExtendedSTL_S = tools::getWorkingDir()+"/inputExtendedSTLMesh_S.stl";
    QString inputExtendedSTL_SH = tools::getWorkingDir()+"/inputExtendedSTLMesh_SH.stl";

    //! ----------------------------------------------
    //! experimental - suppress healing and/or ADMesh
    //! ----------------------------------------------
    if(isHealingOn)
    {
        //! -------------------------------
        //! ADMesh => .stl mesh correction
        //! -------------------------------
        //STLdoctor aSTLDoctor;
        //aSTLDoctor.perform(inputExtendedSTL,inputExtendedSTL_H);

        //! -------------------------
        //! correction using MeshFix
        //! -------------------------
        STLdoctor aSTLDoctor;
        //aSTLDoctor.setProgressIndicator(myProgressIndicator);
        aSTLDoctor.performMeshFix(inputExtendedSTL, inputExtendedSTL_H);
    }
    else
    {
        inputExtendedSTL_H = inputExtendedSTL;
    }

    if(isSimplificationOn)
    {
        //! --------------------
        //! update the progress
        //! --------------------
        if(myProgressIndicator!=Q_NULLPTR)
        {
            QProgressEvent *progressEvent = new QProgressEvent();
            progressEvent->setMessage("Start mesh based defeaturing");
            this->postUpdateEvent(progressEvent);
            Sleep(250);
        }

        int method;
        double surfaceMeshSizeReduction;
        double pairDistance;        

        switch(defeaturingBy)
        {
        case 0: //! using triangle count
        {
            method = 0;
            STLdoctor aSTLDoctor;
            surfaceMeshSizeReduction = defeaturingParameterValue;

            //! feature preserving code here
            aSTLDoctor.simplifyMesh(inputExtendedSTL_H,inputExtendedSTL_S,method,defeaturingParameterValue,-1,fixedPoints);
        }
            break;

        case 1: //! using pair distance
        {
            method = 1;
            STLdoctor aSTLDoctor;
            pairDistance = defeaturingParameterValue;

            //! feature preserving code here
            aSTLDoctor.simplifyMesh(inputExtendedSTL_H,inputExtendedSTL_S,method,-1,defeaturingParameterValue,fixedPoints);
        }
            break;
        }

        //! ----------------------------------------------------------------
        //! last mesh correction using MeshFix at the end of simplification
        //! (in principle the tessellation should be watertight)
        //! ----------------------------------------------------------------
        STLdoctor aSTLDoctor;
        aSTLDoctor.performMeshFix(inputExtendedSTL_S, inputExtendedSTL_SH);

        //! --------------------
        //! update the progress
        //! --------------------
        if(myProgressIndicator!=Q_NULLPTR)
        {
            QProgressEvent *progressEvent = new QProgressEvent();
            progressEvent->setMessage("Mesh based defeaturing done. Start meshing");
            this->postUpdateEvent(progressEvent);
            Sleep(250);
        }
    }
    else
    {
        inputExtendedSTL_SH = inputExtendedSTL_H;
    }

    //! -------------------------
    //! enable the "Stop" button
    //! -------------------------
    if(myProgressIndicator!=Q_NULLPTR) emit requestEnableStop();

    return inputExtendedSTL_SH;
}

//! -----------------------------------------------------
//! function: ExMesh_generateVolumeMesh
//! details:  it works only with Netgen as volume mesher
//! -----------------------------------------------------
userMessage MesherClass::ExMesh_generateVolumeMesh(const meshParam &mp,
                                                     int bodyIndex,
                                                     occHandle(MeshVS_DataSource) &mainMesh2D,
                                                     occHandle(MeshVS_DataSource) &mainMesh3D,
                                                     int done)
{
    cout<<"____meshing surface with EXPRESS MESH____"<<endl;

    //! ----------------------------------
    //! begin generating the surface mesh
    //! ----------------------------------
    const TopoDS_Shape &shape = myMeshDB->bodyMap.value(bodyIndex);
    occMesher theOCCMesher(shape,mp);
    theOCCMesher.setProgressIndicator(myProgressIndicator);
    theOCCMesher.setVolumeMesher(Property::meshEngine3D_Netgen);

    userMessage mr = theOCCMesher.performVolume(mainMesh3D,mainMesh2D,mp.isSecondOrder,true,done);

    if(mr.isDone)
    {
        myMeshDB->ArrayOfMeshDS.insert(bodyIndex,mainMesh3D);
        myMeshDB->ArrayOfMeshDS2D.insert(bodyIndex,mainMesh2D);
        QMap<int,occHandle(Ng_MeshVS_DataSourceFace)> faceMeshDS = theOCCMesher.getFaceDataSources();
        for(QMap<int,occHandle(Ng_MeshVS_DataSourceFace)>::iterator it = faceMeshDS.begin(); it != faceMeshDS.end(); ++it)
        {
            int faceNr = it.key();
            const occHandle(Ng_MeshVS_DataSourceFace) &curFaceMeshDS = it.value();
            myMeshDB->ArrayOfMeshDSOnFaces.setValue(bodyIndex, faceNr, curFaceMeshDS);
        }
    }
    return mr;
}

//! -----------------------
//! function: abortMeshing
//! details:
//! -----------------------
void MesherClass::abort()
{
    cout<<"MesherClass::abort()->____abort received____"<<endl;

    //! ---------------------------------------------------------------------
    //! interrupt all the meshing processes followinig the current in action
    //! ---------------------------------------------------------------------
    cout<<"MesherClass::abort()->____resetting all the mesh data sources____"<<endl;

    //! --------------------------------
    //! reset ALL the mesh data sources
    //! --------------------------------
    for(QMap<int,TopoDS_Shape>::iterator it = myMeshDB->bodyMap.begin(); it!=myMeshDB->bodyMap.end(); ++it)
    {
        int bodyIndex = it.key();
        myMeshDB->ArrayOfMeshDS.insert(bodyIndex,occHandle(MeshVS_DataSource)());
        myMeshDB->ArrayOfMeshDS2D.insert(bodyIndex,occHandle(MeshVS_DataSource)());
        int NbFaces = myMeshDB->MapOfBodyTopologyMap.value(bodyIndex).faceMap.Extent();
        int NbEdges = myMeshDB->MapOfBodyTopologyMap.value(bodyIndex).edgeMap.Extent();
        for(int faceNr=0; faceNr<=NbFaces; faceNr++)
        {
            myMeshDB->ArrayOfMeshDSOnFaces.setValue(bodyIndex,faceNr,occHandle(MeshVS_DataSource)());
        }
        for(int edgeNr=1; edgeNr<=NbEdges; edgeNr++)
        {
            myMeshDB->ArrayOfMeshDSOnEdges.setValue(bodyIndex,edgeNr,occHandle(MeshVS_DataSource)());
        }
    }
    //! --------------------
    //! set the update flag
    //! --------------------
    for(QMap<int,TopoDS_Shape>::iterator it = myMeshDB->bodyMap.begin(); it!=myMeshDB->bodyMap.end(); ++it)
    {
        int bodyIndex = it.key();
        myMeshDB->ArrayOfMeshIsToBeUdpdated.insert(bodyIndex,true);
    }
    cerr<<"MesherClass::abort()->____EMITTING \"ABORT RECEIVED\"____"<<endl;

    //! ------------------------------------------
    //! update progress - remove the progress bar
    //! ------------------------------------------
    QProgressEvent *pe = new QProgressEvent();
    pe->setMessage(QString("Abort request received"));
    this->postUpdateEvent(pe);

    pe = new QProgressEvent(QProgressEventAction::QProgressEvent_Reset,0,100,0,"",QProgressEventAction::QProgressEvent_Reset,0,100,0);
    this->postUpdateEvent(pe);

    //emit meshingFinished();

    //! -------------------------
    //! remove the support files
    //! -------------------------
    this->removeSupportFiles(true);
}

//! --------------------------------------------------
//! function: handleRequestStartingNetgenEnquireTimer
//! details:
//! --------------------------------------------------
void MesherClass::handleRequestStartingNetgenEnquireTimer()
{
    emit requestStartingNetgenEnquireTimer();
}

//! --------------------------------------------------
//! function: handleRequestStoppingNetgenEnquireTimer
//! details:
//! --------------------------------------------------
void MesherClass::handleRequestStoppingNetgenEnquireTimer()
{
    emit requestStoppingNetgenEnquireTimer();
}

//! -------------------------------------------------------
//! function: PrismaticLayers_generatePrismaticMesh_post
//! details:  the volume mesh must be generated in advance
//! -------------------------------------------------------
userMessage MesherClass::PrismaticLayers_generatePrismaticMesh_post(int bodyIndex,
                                                                    occHandle(Ng_MeshVS_DataSource3D) &prismaticMeshDS)
{
    cout<<"MesherClass::PrismaticLayers_generatePrismaticMesh_post()->____function called____"<<endl;

    //! ---------------
    //! meshing result
    //! ---------------
    userMessage mr;

    //! ------------------------
    //! disable the Stop button
    //! ------------------------
    if(myProgressIndicator!=Q_NULLPTR) emit requestDisableStop();

    //! ----------------------------------------------------------------------
    //! initialize the prismatic layer builder: mesh data base and body index
    //! ----------------------------------------------------------------------
    prismaticLayer prismaticLayerBuilder(myMeshDB);
    prismaticLayerBuilder.setBody(bodyIndex);

    //! ---------------------------
    //! initialize the prismatic layer builder: set the progress indicator
    //! ---------------------------
    if(myProgressIndicator!=Q_NULLPTR) prismaticLayerBuilder.setProgressIndicator(myProgressIndicator);

    //! ------------------------------------------------------
    //! initialize the prismatic layer builder
    //! retrieve the list of the prismatic faces for the body
    //! ------------------------------------------------------
    const QList<int> &prismaticFacesOnBody = myMeshDB->prismaticFaces.value(bodyIndex);
    std::vector<int> prismaticFaces = prismaticFacesOnBody.toVector().toStdVector();
    prismaticLayerBuilder.setPrismaticFaces(prismaticFaces);

    //! -------------------------------------------------------
    //! initialize the prismatic layer builder
    //! set the parameters for the prismatic layers generation
    //! -------------------------------------------------------
    const prismaticLayerParameters &parameters = myMeshDB->prismaticMeshParameters.value(bodyIndex);
    prismaticLayerBuilder.setParameters(parameters);

    //! -----------------------------------------------------------------------
    //! displace the surface mesh and "compress" the pre-inflation volume mesh
    //! which should have be present into the mesh database
    //! -----------------------------------------------------------------------
    QList<occHandle(Ng_MeshVS_DataSourceFace)> theInflatedMeshes;
    occHandle(Ng_MeshVS_DataSource3D) preInflationVolumeMeshDS = occHandle(Ng_MeshVS_DataSource3D)::DownCast(myMeshDB->ArrayOfMeshDS.value(bodyIndex));
    bool isDone = prismaticLayerBuilder.inflateMeshAndCompress(theInflatedMeshes,preInflationVolumeMeshDS);
    if(isDone==false)
    {
        mr.isDone = false;
        mr.message = QString("Inflation stopped");
        return mr;
    }

    //! -----------
    //! diagnostic
    //! -----------
    cout<<"MesherClass::PrismaticLayers_generatePrismaticMesh()->____number of inflated meshes: "<<theInflatedMeshes.size()<<"____"<<endl;
    cout<<"@----------------------------------------------------------@"<<endl;
    for(int i=0; i<theInflatedMeshes.length(); i++)
    {
        int NN = theInflatedMeshes.at(i)->GetAllNodes().Extent();
        int NE = theInflatedMeshes.at(i)->GetAllElements().Extent();
        cout<<"@____surface layer nr "<<i<<". Nb Nodes: "<<NN<<"; Nb Elements: "<<NE<<"____"<<endl;
    }
    cout<<"@----------------------------------------------------------@"<<endl;

    //! ---------------------------------------------
    //! generate the prismatic 3D mesh from the list
    //! of the displaced meshes
    //! ---------------------------------------------
    cout<<"MesherClass::PrismaticLayers_generatePrismaticMesh()->____start generating prismatic 3D elements____"<<endl;
    isDone = prismaticLayerBuilder.buildPrismaticElements(theInflatedMeshes,prismaticMeshDS);

    /* disabilito per un momento
    //! -----------------------------------------------------
    //! the "prismaticMeshDS" is inserted in "ArrayOfMeshDS"
    //! as a volume mesh
    //! -----------------------------------------------------
    if(isDone)
    {
        myMeshDB->ArrayOfMeshDS.insert(bodyIndex,prismaticMeshDS);
    }
    else
    {
        myMeshDB->ArrayOfMeshDS.insert(bodyIndex,occHandle(MeshVS_DataSource)());
    }
    */

    //! -----------------------
    //! enable the Stop button
    //! -----------------------
    if(myProgressIndicator!=Q_NULLPTR) emit requestEnableStop();

    //! ---------------------------------------------------------------------------------
    //! check the result: return "true" if at least the first displaced mesh is not null
    //! ---------------------------------------------------------------------------------
    if(theInflatedMeshes.at(1).IsNull())
    {
        mr.isDone = false;
        mr.message = QString("Error in generaring the prismatic mesh");
        return mr;
    }
    mr.isDone = true;
    mr.message = QString("The prismatic mesh has been generated");
    return mr;
}
