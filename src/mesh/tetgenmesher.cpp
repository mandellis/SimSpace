//#define TETGEN_PATH "D:\\tetgen1.5.0\\build\\Release\\tetgen.exe"
#define TETGEN_PATH "C:\\Tetgen\\tetgen.exe"

//! ----------------
//! custom includes
//! ----------------
#include <tetgenmesher.h>
#include <ng_meshvs_datasource2d.h>
#include "meshtools.h"
#include <ng_meshvs_datasource3d.h>
#include <ng_meshvs_datasourceface.h>
#include "se_exception.h"
#include "mydefines.h"
#include "ccout.h"
#include "tools.h"
#include "qprogressevent.h"

//! ----
//! OCC
//! ----
#include <MeshVS_DataSource.hxx>
#include <TColStd_PackedMapOfInteger.hxx>
#include <TColStd_MapIteratorOfPackedMapOfInteger.hxx>
#include <TColStd_IndexedMapOfInteger.hxx>
#include <TopAbs_ShapeEnum.hxx>
#include <TopoDS_Face.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopTools_MapIteratorOfMapOfShape.hxx>

//! ---
//! Qt
//! ---
#include <QMessageBox>
#include <QDir>
#include <QApplication>
#include <QProcess>
#include <QWidget>

//! ----
//! C++
//! ----
#include <iostream>
using namespace std;

//! -------
//! global
//! -------
#include <global.h>

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
TetgenMesher::TetgenMesher(meshDataBase *mDB, QObject *parent):QObject(parent)
{
    myMeshDB = mDB;
    tetgenProcess = new QProcess(this);
    this->createTetgenDirs();
}

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
TetgenMesher::TetgenMesher(QObject *parent):QObject(parent)
{
    tetgenProcess = new QProcess(this);
    this->createTetgenDirs();
}

//! ----------------------
//! function: setDataBase
//! details:
//! ----------------------
void TetgenMesher::init(meshDataBase *mDB)
{
    myMeshDB = mDB;
}

//! -------------------------------
//! function: setProgressIndicator
//! details:
//! -------------------------------
void TetgenMesher::setProgressIndicator (QProgressIndicator *aProgressIndicator)
{
    if(aProgressIndicator!=Q_NULLPTR)
    {
        myProgressIndicator = aProgressIndicator;
        disconnect(this,SIGNAL(requestEnableStop()),myProgressIndicator,SLOT(enableStop()));
        disconnect(this,SIGNAL(requestDisableStop()),myProgressIndicator,SLOT(disableStop()));
        connect(this,SIGNAL(requestEnableStop()),myProgressIndicator,SLOT(enableStop()));
        connect(this,SIGNAL(requestDisableStop()),myProgressIndicator,SLOT(disableStop()));
    }
}

#ifndef PLCNEW
//! -------------------------------------------------------
//! function: buildPLC
//! details:  build the PLC reprentation of a TopoDS_Shape
//! -------------------------------------------------------
bool TetgenMesher::buildPLC(int bodyIndex, QList<int> &invalidFaceTags, bool saveTetgenFiles)
{
    cout<<"------------------------------"<<endl;
    cout<<"- TETGEN IS BUILDING THE PLC -"<<endl;
    cout<<"------------------------------"<<endl;

    myPLC.initialize();
    myPLC.firstnumber = 1;

    //! ----------------------------
    //! the overall 2D surface mesh
    //! ----------------------------
    occHandle(Ng_MeshVS_DataSource2D) meshDS2D = occHandle(Ng_MeshVS_DataSource2D)::DownCast(myMeshDB->ArrayOfMeshDS2D.value(bodyIndex));
    if(!meshDS2D.IsNull())
    {
        //! ------------------------------
        //! number of surface mesh points
        //! ------------------------------
        myPLC.numberofpoints = meshDS2D->GetAllNodes().Extent();

        cout<<"TetgenMesher::init()->____number of PLC nodes: "<<myPLC.numberofpoints<<"____"<<endl;

        //! ------------
        //! points list
        //! ------------
        myPLC.pointlist = new REAL [myPLC.numberofpoints*3];          //! list of REAL
        myPLC.pointmarkerlist = new int[myPLC.numberofpoints*3];      //! list of int

        //! --------------------------------------
        //! nodes coordinates from the datasource
        //! --------------------------------------
        TColStd_PackedMapOfInteger nodeMap = meshDS2D->GetAllNodes();
        TColStd_MapIteratorOfPackedMapOfInteger nodeIter;
        Standard_Real aCoordsBuf[3];
        TColStd_Array1OfReal Coords(*aCoordsBuf,1,3);

        int S=0;
        Standard_Integer nbNodes;
        MeshVS_EntityType aType;
        for(nodeIter.Initialize(nodeMap);nodeIter.More();nodeIter.Next())
        {
            int nodeID = nodeIter.Key();
            if(!meshDS2D->GetGeom(nodeID,Standard_False,Coords,nbNodes,aType))
            {
                cout<<"TetgenMesher()->TetgenMesher::init()->____nodeID "<<nodeID<<" not found____"<<endl;
                continue;
            }

            for(int j=1;j<=3;j++)
            {
                myPLC.pointlist[S]=Coords(j);
                //! ----------------------------------
                //! tag the point as a boundary point
                //! ----------------------------------
                myPLC.pointmarkerlist[S]=1;
                S++;
            }
        }

        //! -----------------------------------------
        //! each surface element is taken as a facet
        //! -----------------------------------------
        int NSE = meshDS2D->GetAllElements().Extent();
        myPLC.numberoffacets = NSE;

        cout<<"TetgenMesher::init()->____number of PLC facets: "<<NSE<<"____"<<endl;

        myPLC.facetmarkerlist = new int[NSE];
        myPLC.facetlist = new tetgenio::facet[NSE];

        int NFaces = myMeshDB->MapOfBodyTopologyMap.value(bodyIndex).faceMap.Extent();

        cout<<"TetgenMesher::init()->____number of geometry faces: "<<NFaces<<"____"<<endl;

        //! -----------------------------------------------
        //! build the PLC using the face mesh data sources
        //! -----------------------------------------------
        for(int faceNr=0, S=0; faceNr<=NFaces; faceNr++)
        {
            //!cout<<"TetgenMesher::init()->____Working on face Nr: "<<faceNr<<"____"<<endl;

            const occHandle(Ng_MeshVS_DataSourceFace) &curFaceDS = occHandle(Ng_MeshVS_DataSourceFace)::DownCast(myMeshDB->ArrayOfMeshDSOnFaces.getValue(bodyIndex,faceNr));

            //! -------------------------------
            //! one of the face meshes is null
            //! -------------------------------
            if((curFaceDS.IsNull() || curFaceDS->GetAllElements().Extent()==0) && faceNr==0)
            {
                //!cout<<"TetgenMesher::init()->____the fake ace Nr: "<<faceNr<<" has no elements: jumping over this face"<<endl;

                invalidFaceTags<<faceNr;
                continue;
            }
            if((curFaceDS.IsNull() || curFaceDS->GetAllElements().Extent()==0) && faceNr!=0)
            {
                //!cout<<"TetgenMesher::init()->____the face Nr: "<<faceNr<<" has no elements: jumping over this face"<<endl;

                invalidFaceTags<<faceNr;
                continue;
            }
            else
            {
                const occHandle(Ng_MeshVS_DataSourceFace) &faceMeshDS =
                        occHandle(Ng_MeshVS_DataSourceFace)::DownCast(myMeshDB->ArrayOfMeshDSOnFaces.getValue(bodyIndex,faceNr));

                int NbNodes;
                Standard_Integer aNodeBuf[3];
                TColStd_Array1OfInteger nodeIDs(*aNodeBuf,1,3);

                //! ------------------------------------------------------
                //! scan all the elements (triangles) of the surface mesh
                //! el is the local (face) element number
                //! ------------------------------------------------------
                for(int el = 1; el<=faceMeshDS->GetAllElements().Extent(); el++)
                {
                    //cout<<"_____Adding facet Nr: "<<el<<" ->";
                    tetgenio::facet *facet = &myPLC.facetlist[el-1+S];

                    //! the "facet" is a triangle
                    facet->numberofpolygons = 1;
                    facet->numberofholes = 0;
                    facet->holelist = NULL;

                    myPLC.facetmarkerlist[el-1+S] = faceNr;

                    facet->polygonlist = new tetgenio::polygon[1];
                    tetgenio::polygon *polygon = &facet->polygonlist[0];

                    int globalElementID =faceMeshDS->myElementsMap.FindKey(el);

                    //! ------------------------
                    //! use the overall 2D mesh
                    //! ------------------------
                    bool nodesOfTheElementFound = meshDS2D->GetNodesByElement(globalElementID,nodeIDs,NbNodes);

                    if(nodesOfTheElementFound!=true)
                    {
                        //! 2D mesh not consistent
                        return false;
                    }

                    polygon->numberofvertices = NbNodes;
                    polygon->vertexlist = new int[NbNodes];

                    for(int j=0; j<NbNodes; j++)
                    {
                        polygon->vertexlist[j] = nodeIDs(j+1);
                    }
                }
                S = S + faceMeshDS->GetAllElements().Extent();
            }
        }

        //! -----------
        //! diagnostic
        //! -----------
        cout<<"TetgenMesher::init()->____the following faces have not a mesh: ";
        for(int i=0; i<invalidFaceTags.length()-1; i++) cout<<invalidFaceTags.at(i)<<", ";
        cout<<invalidFaceTags.last()<<endl;

        cout<<"TetgenMesher::init()->____The PLC has been built____"<<endl;

        if(saveTetgenFiles)
        {
            cout<<"TetgenMesher::init()->____Tetgen is writing the PLC on disk____"<<endl;
            char fp[512];
            QDir curDir(myPLCDir);
            QString fileName = myMeshDB->MapOfBodyNames.value(bodyIndex);
            QString filePath = curDir.absolutePath()+"/"+fileName;
            sprintf(fp,filePath.toStdString().c_str());

            cout<<"TetgenMesher::init->___saving PLC nodes : "<<filePath.toStdString()<<"____"<<endl;
            myPLC.save_nodes(fp);

            cout<<"TetgenMesher::init->___saving PLC facets: "<<filePath.toStdString()<<"____"<<endl;
            myPLC.save_poly(fp);

            cout<<"TetgenMesher::init()->____nodes and facets written____"<<endl;
        }
        return true;
    }
    else
    {
        cout<<"TetgenMesher::init()->____error in generating the PLC: the input surface mesh is null____"<<endl;
        return false;
    }
}
#endif

//! -------------------------------------------------------------------------
//! function: SEfunc
//! details:  delocalize perform in order to catch also unhandled exceptions
//! -------------------------------------------------------------------------
void TetgenMesher::SEFunc(tetgenio *meshOut)
{
    __try
    {
        //! ------------------------------
        //! setup the switches
        //! p - takes as an input .poly
        //! Y - preserve the surface mesh
        //! V - verbose
        //! q - quality meshing
        //! a - maximum volume constraint
        //! ------------------------------
        cout<<"TetgenMesher::SEFunc()->____try: tetrahedralize called____"<<endl;
        tetrahedralize(mySwitches, &myPLC, meshOut);
        cout<<"TetgenMesher::SEFunc()->__try: MESHING DONE____"<<endl;
    }
    __finally
    {
        cout<<"TetgenMesher::SEFunc()->____in finally____"<<endl;
    }
}

void trans_func_T(unsigned int u, EXCEPTION_POINTERS* pExp)
{
    Q_UNUSED(pExp)
    Q_UNUSED(u)
    printf( "In trans_func.\n" );
    throw SE_Exception();
}

//! -------------------------
//! function: perform
//! details:  build the mesh
//! -------------------------
bool TetgenMesher::perform(tetgenio *meshOut, int bodyIndex, bool saveMesh)
{
    cout<<"TetgenMesher::perform()->____perform called____"<<endl;
    try
    {
        _set_se_translator(trans_func_T);
        cout<<"TetgenMesher::perform()->____Tetgen running in memory switches: "<<mySwitches<<"____"<<endl;
        this->SEFunc(meshOut);
    }
    catch(SE_Exception e)
    {
        Q_UNUSED(e)
        cout<<"TetgenMesher::perform()->____Caught a __try exception with SE_Exception"<<endl;
        return false;
    }
    catch(int tetgenError)
    {
        QString errorMsg;
        switch(tetgenError)
        {
        case 1: errorMsg = "Tetgen runs out of memory"; break;
        case 2: errorMsg = "Please report this bug to Hang.Si@wias-berlin.de.\n"
                           "Include the message above, your input data set, and the exact\n"
                           "command line you used to run this program"; break;
        case 3: errorMsg = "A self-intersection was detected\n"
                           "Hint: use -d option to detect all self-intersections"; break;
        case 4: errorMsg = "A very small input feature size was detected\n"
                           "Hint: use -T option to set a smaller tolerance\n"
                           " If you want to ignore this possible error, set a smaller tolerance\n"
                           "by the -T switch, default is 1e−8"; break;
        case 5: errorMsg = "Two very close input facets were detected. Program stopped.\n"
                           "Hint: use -Y option to avoid adding Steiner points in boundary"; break;
        case 10: errorMsg = "An input error was detected. Program stopped"; break;
        default: break;
        }
        return false;
    }

    //! --------------
    //! save the mesh
    //! --------------
    if(saveMesh)
    {
        QDir curDir(myPLCMesh);
        QString outFileName = myMeshDB->MapOfBodyNames.value(bodyIndex);
        QString outFilePath = curDir.absolutePath()+"/"+outFileName+"_out";
        char fpout[512];
        sprintf(fpout,outFilePath.toStdString().c_str());
        ccout(QString("TetgenMesher::perform()->____saving Tetgen final mesh: ")+outFilePath+"_____");
        meshOut->save_nodes(fpout);
        meshOut->save_elements(fpout);
        meshOut->save_faces(fpout);
    }
    return true;
}

//! ---------------------------------------------------------
//! function: setSwitches
//! details:  set meshing parameters
//!           When defining the switch string do not use "-"
//! ---------------------------------------------------------
void TetgenMesher::setSwitches(int preserveSurfaceMesh, int bodyIndex)
{
    double tolerance = TOLERANCE;
    double L = myMeshDB->ArrayOfMaxBodyElementSize.value(bodyIndex);
    double maxVolSize = (4.0/3.0)*3.14159*pow(L,3);

    switch(preserveSurfaceMesh)
    {
    case 0:
    {
        //! ------------------------------------------------
        //! surface elements can be changed
        //! O: optimization level 10, 7 edge/face flips,
        //!    vertex smoothing, vertex insertion/deletion:
        //!    Attach "O10/7"
        //! ------------------------------------------------
        sprintf(mySwitches,"pq1.4/0.0a%.2lfT%.2eV",maxVolSize,tolerance);
    }
        break;

    case 1:
    {
        //! ------------------------------------------------
        //! 1: optimization level 10, 7 edge/face flips,
        //!    vertex smoothing, vertex insertion/deletion:
        //!    Attach "O10/7"
        //! ------------------------------------------------
        sprintf(mySwitches,"pYq1.4/0.0a%.2lfT%.2eV",maxVolSize,tolerance);
    }
        break;
    }
}

//! ---------------------------
//! function: createTetgenDirs
//! details:
//! ---------------------------
void TetgenMesher::createTetgenDirs()
{
    //! ---------------------------------------------------
    //! create the directories for storing the tetgen data
    //! ---------------------------------------------------
    QDir curDir;
    curDir.cd(tools::getWorkingDir());

    ccout("TetgenMesher::creatingTetgenDirs->____creating the Tetgen directories____");

    if(curDir.cd("Tetgen PLC")==false)
    {
        curDir.mkdir("Tetgen PLC");
        curDir.cd("Tetgen PLC");
        myPLCDir = curDir.absolutePath();
        curDir.cdUp();
    }
    else
    {
        curDir.cd("Tetgen PLC");
        myPLCDir = curDir.absolutePath();
        tools::clearDir(curDir.absolutePath());
        curDir.cdUp();
    }

    if(curDir.cd("Tetgen mesh")==false)
    {
        curDir.mkdir("Tetgen mesh");
        curDir.cd("Tetgen mesh");
        myPLCMesh =curDir.absolutePath();
        curDir.cdUp();
    }
    else
    {
        curDir.cd("Tetgen mesh");
        myPLCMesh =curDir.absolutePath();
        tools::clearDir(curDir.absolutePath());
        curDir.cdUp();
    }
}

//!-----------------------------------------------------------------
//! function: createTetgenSupportFilesDir
//! details:  a directory for support files when meshing using disk
//!-----------------------------------------------------------------
void TetgenMesher::createTetgenSupportFilesDir(bool clearPreviousContent)
{
    QDir curDir;
    curDir.cd(tools::getWorkingDir());

    ccout("TetgenMesher::createTetgenSupportFilesDir->____creating the Tetgen directories____");
    if(curDir.cd("TetgenSupportFiles")==false)
    {
        curDir.mkdir("TetgenSupportFiles");
        curDir.cd("TetgenSupportFiles");
        myTetgenSupportFilesDir =curDir.absolutePath();
        curDir.cdUp();
    }
    else
    {
        curDir.cd("TetgenSupportFiles");
        myTetgenSupportFilesDir =curDir.absolutePath();
        if(clearPreviousContent==true) tools::clearDir(curDir.absolutePath());
        curDir.cdUp();
    }
}

//! ------------------------------------------
//! function: retrieveMeshDataSourcesFromDisk
//! details:
//! ------------------------------------------
bool TetgenMesher::retrieveMeshDataSourcesFromDisk(int bodyIndex,
                                                   occHandle(Ng_MeshVS_DataSource3D) &volumeMeshDS,
                                                   occHandle(Ng_MeshVS_DataSource2D) &surfaceMeshDS,
                                                   NCollection_Array1<occHandle(Ng_MeshVS_DataSourceFace)> &arrayOfFaceMeshDS,
                                                   int done)
{
    //! ----------------------------------------------------------------------------------
    //! use the body name and append the body index (two bodies could have the same name)
    //! ----------------------------------------------------------------------------------
    QString bodyName = myMeshDB->MapOfBodyNames.value(bodyIndex)+QString("_%1").arg(bodyIndex);

    QString eleFileName = myTetgenSupportFilesDir+"/"+bodyName+".1.ele";
    QString nodeFileName = myTetgenSupportFilesDir+"/"+bodyName+".1.node";
    QString faceFileName = myTetgenSupportFilesDir+"/"+bodyName+".1.face";
    QString edgeFileName = myTetgenSupportFilesDir+"/"+bodyName+".1.edge";

    cout<<"TetgenMesher::retrieveMeshDataSourcesFromDisk()->____retrieving the mesh datasources for body: "<<bodyName.toStdString()<<"____"<<endl;

    //! --------------------------------
    //! init the secondary progress bar
    //! --------------------------------
    int N_events = myMeshDB->MapOfBodyTopologyMap.value(bodyIndex).faceMap.Extent()+2;

    if(myProgressIndicator!=Q_NULLPTR)
    {
        QProgressEvent *pe = new QProgressEvent();
        pe->setVal(done);
        pe->setMessage("Generating the face mesh data sources");
        pe->setAction1(QProgressEventAction::QProgressEvent_Init);
        pe->setRange1(0, N_events);
        QApplication::postEvent(myProgressIndicator,pe);
        QApplication::processEvents();
    }

    //! -------------------------------------------------------------
    //! test if the .ele file exixts - the control is redundant
    //! since already included in Ng_MeshVS_DataSource3D constructor
    //! -------------------------------------------------------------
    FILE *eleFile = fopen(eleFileName.toStdString().c_str(),"r");
    if(eleFile==NULL)
    {
        cout<<"TetgenMesher::retrieveMeshDataSourcesFromDisk()->____cannot retrieve the volume mesh datasource for body: "<<bodyName.toStdString()<<"____"<<endl;
        return false;
    }
    fclose(eleFile);
    volumeMeshDS = new Ng_MeshVS_DataSource3D(eleFileName,nodeFileName);

    //! -----------------------------------------------
    //! generate the face to element connectivity data
    //! -----------------------------------------------
    cout<<"TetgenMesher::retrieveMeshDataSourcesFromDisk()->____generating the face to element connectivity data____"<<endl;
    this->setStopButtonEnabled(false);
    volumeMeshDS->buildFaceToElementConnectivity();

    //! --------------------------------------------
    //! post an event: the surface mesh data source
    //! --------------------------------------------
    if(myProgressIndicator!=Q_NULLPTR)
    {
        QProgressEvent *pe = new QProgressEvent();
        pe->setVal(done);
        pe->setMessage(QString("Processing volume mesh"));
        pe->setVal1(1);
        QApplication::postEvent(myProgressIndicator,pe);
        QApplication::processEvents();
    }

    if(volumeMeshDS.IsNull())
    {
        cout<<"TetgenMesher::retrieveMeshDataSourcesFromDisk()->____error in generating the volume mesh____"<<endl;
        return false;
    }
    cout<<"TetgenMesher::retrieveMeshDataSourcesFromDisk()->____volume mesh OK____"<<endl;

    //! ------------------------------
    //! test if the .face file exists
    //! ------------------------------
    cout<<"TetgenMesher::retrieveMeshDataSourcesFromDisk()->____retrieving the overall surface mesh____"<<endl;
    FILE *faceFile = fopen(faceFileName.toStdString().c_str(),"r");
    if(faceFile==NULL)
    {
        cout<<"TetgenMesher::retrieveMeshDataSourcesFromDisk()->____cannot retrieve the surface mesh datasource for body: "<<bodyName.toStdString()<<"____"<<endl;
        return false;
    }
    fclose(faceFile);

    surfaceMeshDS = new Ng_MeshVS_DataSource2D(faceFileName,nodeFileName);
    if(surfaceMeshDS.IsNull())
    {
        cout<<"TetgenMesher::retrieveMeshDataSourcesFromDisk()->____error in generating the surface mesh____"<<endl;
        return false;
    }

    //! --------------------------------------------
    //! post an event: the surface mesh data source
    //! --------------------------------------------
    if(myProgressIndicator!=Q_NULLPTR)
    {
        QProgressEvent *pe = new QProgressEvent();
        pe->setVal(done);
        pe->setMessage(QString("Processing volume mesh"));
        pe->setVal1(2);
        QApplication::postEvent(myProgressIndicator,pe);
        QApplication::processEvents();
    }

    cout<<"TetgenMesher::retrieveMeshDataSourcesFromDisk()->____surface mesh OK____"<<endl;

    //! -----------------------------------------------------
    //! generate the face mesh data sources
    //! (the existance of the files has been checked before)
    //! Important note: even if the fake "0" face is used
    //! for building the PLC, tetgen will not generate the
    //! mesh data source for that face, since considered a
    //! "correction", and assumed to be of "small" area.
    //! For that reason the cycle starts from "1", and the
    //! "0" face explicitly is nullified
    //! -----------------------------------------------------

    //QList<QList<double>> listOfMeshPoints = TetgenMesher::readNodeFile(nodeFileName);
    //QList<QList<int>> listOfFaces = TetgenMesher::readFaceFile(faceFileName);

    int NbGeometryFaces = myMeshDB->MapOfBodyTopologyMap.value(bodyIndex).faceMap.Extent();

    //! -------------------
    //! read the node file
    //! -------------------
    QMap<int,mesh::meshPoint> indexedMapOfMeshPoints;
    this->readNodeFile(nodeFileName,indexedMapOfMeshPoints);

    this->setStopButtonEnabled(true);

    occHandle(Ng_MeshVS_DataSourceFace) aFaceMeshDS;
    for(int faceNr = 1; faceNr<=NbGeometryFaces; faceNr++)
    {
        if(Global::status().code == 0) return false;    // process interrupted by the user

        cout<<"TetgenMesher::retrieveMeshDataSourcesFromDisk()->____retrieving the face mesh for face nr. "<<faceNr<<"____"<<endl;

        //aFaceMeshDS = new Ng_MeshVS_DataSourceFace(faceFileName,nodeFileName,faceNr);

        aFaceMeshDS = new Ng_MeshVS_DataSourceFace(faceFileName,indexedMapOfMeshPoints,faceNr);
        //aFaceMeshDS = new Ng_MeshVS_DataSourceFace(listOfMeshPoints,listOfFaces,faceNr);  //cesere

        arrayOfFaceMeshDS.SetValue(faceNr,aFaceMeshDS);

        //! -------------------------------------------
        //! post update events on face mesh generation
        //! -------------------------------------------
        if(faceNr%50==0 && myProgressIndicator!=Q_NULLPTR)
        {
            QProgressEvent *pe = new QProgressEvent(QProgressEvent_None,0,9999,0,"Generate face mesh data sources",
                                                    QProgressEvent_Update,0,9999,faceNr,"Tetgen building face mesh datasources from disk");
            QApplication::postEvent(myProgressIndicator,pe);
            QApplication::processEvents();
        }
    }

    //! ---------------------
    //! nullify the "0" face
    //! ---------------------
    cout<<"TetgenMesher::retrieveMeshDataSourcesFromDisk()->____nullifying the \"0\" face____"<<endl;
    arrayOfFaceMeshDS.SetValue(0,occHandle(Ng_MeshVS_DataSourceFace)());

    return true;
}

//! -------------------------------------------
//! function: retrieveVolumeMeshDataSourceOnly
//! details:
//! -------------------------------------------
bool TetgenMesher::retrieveVolumeMeshDataSourceOnly(int bodyIndex, occHandle(Ng_MeshVS_DataSource3D) &volumeMeshDS)
{
    //! ----------------------------------------------------------------------------------
    //! use the body name and append the body index (two bodies could have the same name)
    //! ----------------------------------------------------------------------------------
    QString bodyName = myMeshDB->MapOfBodyNames.value(bodyIndex)+QString("_%1").arg(bodyIndex);

    QString eleFileName = myTetgenSupportFilesDir+"/"+bodyName+".1.ele";
    QString nodeFileName = myTetgenSupportFilesDir+"/"+bodyName+".1.node";
    QString faceFileName = myTetgenSupportFilesDir+"/"+bodyName+".1.face";
    QString edgeFileName = myTetgenSupportFilesDir+"/"+bodyName+".1.edge";

    cout<<"TetgenMesher::retrieveMeshDataSourcesFromDisk()->____retrieving the volume mesh datasources for body: "<<bodyName.toStdString()<<"____"<<endl;

    //! -----------------------------------------------------------------
    //! test if the .ele file exixts - the control is actually redundant
    //! since already included in Ng_MeshVS_DataSource3D constructor
    //! -----------------------------------------------------------------
    FILE *eleFile = fopen(eleFileName.toStdString().c_str(),"r");
    if(eleFile==NULL)
    {
        cout<<"TetgenMesher::retrieveMeshDataSourcesFromDisk()->____cannot retrieve the volume mesh datasource for body: "<<bodyName.toStdString()<<"____"<<endl;
        return false;
    }
    fclose(eleFile);
    volumeMeshDS = new Ng_MeshVS_DataSource3D(eleFileName,nodeFileName);
    return true;
}

#ifdef PLCNEW
//! -------------------------------------------------------------
//! function: buildPLC
//! details:  build the PLC reprentation of a TopoDS_Shape
//!           this version does not use the overall surface mesh
//! -------------------------------------------------------------
bool TetgenMesher::buildPLC(int bodyIndex, QList<int> &invalidFaceTags, bool saveTetgenFiles)
{
    ccout("TetgenMesher::init()->____function called____");

    if(myMeshDB!=NULL)
    {
        cout<<"------------------------------"<<endl;
        cout<<"- TETGEN IS BUILDING THE PLC -"<<endl;
        cout<<"------------------------------"<<endl;

        ccout(QString("------------------------------"));
        ccout(QString("- TETGEN IS BUILDING THE PLC -"));
        ccout(QString("------------------------------"));

        QList<QVector<double>> pointList;

        myPLC.initialize();
        myPLC.firstnumber = 1;

        int NbFacets = 0;
        int NbTopologyFaces = myMeshDB->MapOfBodyTopologyMap.value(bodyIndex).faceMap.Extent();

        cout<<"TetgenMesher::buildPLC()->____number of topology faces: "<<NbTopologyFaces<<"____"<<endl;
        ccout(QString("TetgenMesher::buildPLC()->____number of topology faces: %1____").arg(NbTopologyFaces));

        int S = 0;
        for(int faceNr = 0; faceNr<=NbTopologyFaces; faceNr++)
        {
            const occHandle(Ng_MeshVS_DataSourceFace) &curFaceMeshDS = occHandle(Ng_MeshVS_DataSourceFace)::DownCast(myMeshDB->ArrayOfMeshDSOnFaces.getValue(bodyIndex,faceNr));
            if(!curFaceMeshDS.IsNull())
            {
                if(curFaceMeshDS->GetAllElements().Extent()>0)
                {
                    cout<<"TetgenMesher::buildPLC()->____working (1) on face nr: "<<faceNr<<"____"<<endl;
                    ccout(QString("TetgenMesher::buildPLC()->____working (1) on face nr: %1____").arg(faceNr));

                    TColStd_PackedMapOfInteger eMap = curFaceMeshDS->GetAllElements();
                    TColStd_MapIteratorOfPackedMapOfInteger eIt;
                    MeshVS_EntityType type;
                    int NbNodes;
                    double buf[24];
                    TColStd_Array1OfReal coords(*buf,1,24);
                    for(eIt.Initialize(eMap); eIt.More(); eIt.Next(), S++)
                    {
                        int ID = eIt.Key();
                        curFaceMeshDS->GetGeom(ID,true,coords,NbNodes,type);
                        for(int k=0; k<NbNodes; k++)
                        {
                            QVector<double> P;
                            P.push_back(coords(3*k+1));
                            P.push_back(coords(3*k+2));
                            P.push_back(coords(3*k+3));
                            if(!pointList.contains(P))
                            {
                                pointList<<P;
                            }
                        }
                    }
                    NbFacets = NbFacets + curFaceMeshDS->GetAllElements().Extent();
                }
            }
            else
            {
                cout<<"TetgenMesher::buildPLC()->____working (1) on face nr: "<<faceNr<<". The face has a null mesh: jumping over it____"<<endl;
                ccout(QString("TetgenMesher::buildPLC()->____working (1) on face nr: %1. The face has a null mesh: jumping over it____").arg(faceNr));

                invalidFaceTags<<faceNr;
            }
        }

        myPLC.numberoffacets = NbFacets;
        myPLC.facetlist = new tetgenio::facet[NbFacets];
        myPLC.facetmarkerlist = new int[NbFacets];

        ccout(QString("____Number of facets %1").arg(NbFacets));

        S = 0;
        for(int faceNr=0; faceNr<=NbTopologyFaces; faceNr++)
        {
            const occHandle(Ng_MeshVS_DataSourceFace) &curFaceMeshDS = occHandle(Ng_MeshVS_DataSourceFace)::DownCast(myMeshDB->ArrayOfMeshDSOnFaces.getValue(bodyIndex,faceNr));
            if(!curFaceMeshDS.IsNull())
            {
                if(curFaceMeshDS->GetAllElements().Extent()>0)
                {
                    cout<<"TetgenMesher::buildPLC()->____working (2) on face nr: "<<faceNr<<"____"<<endl;
                    ccout(QString("____working on face nr: %1____").arg(faceNr));

                    TColStd_PackedMapOfInteger eMap = curFaceMeshDS->GetAllElements();
                    TColStd_MapIteratorOfPackedMapOfInteger eIt;

                    int el = 1;
                    for(eIt.Initialize(eMap);eIt.More();eIt.Next(), el++)
                    {
                        int ID = eIt.Key();
                        int NbNodes;
                        double buf[24];
                        MeshVS_EntityType type;
                        TColStd_Array1OfReal coords(*buf,1,24);
                        curFaceMeshDS->GetGeom(ID,true,coords,NbNodes,type);

                        myPLC.facetmarkerlist[el-1+S] = faceNr;

                        tetgenio::facet *facet = &myPLC.facetlist[el-1+S];
                        facet->numberofpolygons = 1;
                        facet->numberofholes = 0;
                        facet->holelist = NULL;
                        facet->polygonlist = new tetgenio::polygon[1];
                        tetgenio::polygon *polygon = &facet->polygonlist[0];

                        polygon->numberofvertices = NbNodes;
                        polygon->vertexlist = new int[NbNodes];

                        int index;
                        for(int j=0; j<NbNodes; j++)
                        {
                            QVector<double> P;
                            P.push_back(coords(3*j+1));
                            P.push_back(coords(3*j+2));
                            P.push_back(coords(3*j+3));

                            polygon->vertexlist[j] = pointList.indexOf(P)+1;
                            index = polygon->vertexlist[j];
                        }
                        //cout<<"____"<<el-1+S<<"____"<<NbNodes<<"____"<<"____"<<index<<"____"<<endl;
                    }
                }
                S = S+curFaceMeshDS->GetAllElements().Extent();
            }
        }

        //! --------------------
        //! surface mesh points
        //! --------------------
        myPLC.numberofpoints = pointList.length();
        myPLC.pointlist = new REAL [myPLC.numberofpoints*3];
        myPLC.pointmarkerlist = new int[myPLC.numberofpoints*3];

        ccout(QString("____Number of points: %1____").arg(pointList.length()));

        for(int S=0; S<pointList.length(); S++)
        {
            const QVector<double> &P = pointList.at(S);
            for(int j=0;j<3;j++)
            {
                myPLC.pointlist[3*S+j]=P.at(j);
            }
            //! ----------------------------------
            //! tag the point as a boundary point
            //! ----------------------------------
            myPLC.pointmarkerlist[S]=1;
        }

        cout<<"------------------------------"<<endl;
        cout<<"- THE PLC HAS BEEN BUILT     -"<<endl;
        cout<<"------------------------------"<<endl;

        ccout(QString("------------------------------"));
        ccout(QString("- THE PLC HAS BEEN BUILT     -"));
        ccout(QString("------------------------------"));

        if(saveTetgenFiles)
        {
            cout<<"------------------------------"<<endl;
            cout<<"- TETGEN IS WRITING THE PlC  -"<<endl;
            cout<<"------------------------------"<<endl;
            ccout("------------------------------");
            ccout("- TETGEN IS WRITING THE PlC  -");
            ccout("------------------------------");

            char fp[512];
            QDir curDir(myPLCDir);
            QString fileName = myMeshDB->mapOfBodyNames.value(bodyIndex);
            QString filePath = curDir.absolutePath()+"/"+fileName;
            sprintf(fp,filePath.toStdString().c_str());

            ccout(QString("TetgenMesher::init->___saving PLC nodes : ").append(filePath).append("____"));
            myPLC.save_nodes(fp);

            ccout(QString("TetgenMesher::init->___saving PLC facets: ").append(filePath).append("____"));
            myPLC.save_poly(fp);

            cout<<"------------------------------"<<endl;
            cout<<"- NODES AND FACES WRITTEN    -"<<endl;
            cout<<"------------------------------"<<endl;
            ccout("------------------------------");
            ccout("- NODES AND FACES WRITTEN    -");
            ccout("------------------------------");
        }
        return true;
    }
    else
    {
        cout<<"------------------------------------------------------"<<endl;
        cout<<"- CANNOT BUILD THE PLC: THE MESH DATABASE IS NULL    -"<<endl;
        cout<<"------------------------------------------------------"<<endl;
        ccout("------------------------------------------------------");
        ccout("- CANNOT BUILD THE PLC: THE MESH DATABASE IS NULL    -");
        ccout("------------------------------------------------------");
        return false;
    }
}
#endif

//! --------------------------------
//! function: redirectStandardOuput
//! details:
//! --------------------------------
void TetgenMesher::redirectTetgenOutput()
{
    std::string message = tetgenProcess->readAllStandardOutput().toStdString();
    QProgressEvent *progressEvent;
    if(QString::fromStdString(message).contains(QString("Initializing memorypools.")))
    {
        progressEvent = new QProgressEvent(QProgressEvent_None,0,9999,0,"",QProgressEvent_Update,0,9999,1,"Tetgen meshing running on disk",this);
    }
    if(QString::fromStdString(message).contains(QString("Initializing robust predicates.")))
    {
        progressEvent = new QProgressEvent(QProgressEvent_None,0,9999,0,"",QProgressEvent_Update,0,9999,2,"Tetgen meshing running on disk",this);
    }
    if(QString::fromStdString(message).contains(QString("Delaunizing vertices...")))
    {
        progressEvent = new QProgressEvent(QProgressEvent_None,0,9999,0,"",QProgressEvent_Update,0,9999,3,"Tetgen meshing running on disk",this);
    }
    if(QString::fromStdString(message).contains(QString("Creating surface mesh ...")))
    {
        progressEvent = new QProgressEvent(QProgressEvent_None,0,9999,0,"",QProgressEvent_Update,0,9999,4,"Tetgen meshing running on disk",this);
    }
    if(QString::fromStdString(message).contains(QString("Recovering boundaries...")))
    {
        progressEvent = new QProgressEvent(QProgressEvent_None,0,9999,0,"",QProgressEvent_Update,0,9999,5,"Tetgen meshing running on disk",this);
    }
    if(QString::fromStdString(message).contains(QString("Removing exterior tetrahedra ...")))
    {
        progressEvent = new QProgressEvent(QProgressEvent_None,0,9999,0,"",QProgressEvent_Update,0,9999,6,"Tetgen meshing running on disk",this);
    }
    if(QString::fromStdString(message).contains(QString("Recovering Delaunayness...")))
    {
        progressEvent = new QProgressEvent(QProgressEvent_None,0,9999,0,"",QProgressEvent_Update,0,9999,7,"Tetgen meshing running on disk",this);
    }
    if(QString::fromStdString(message).contains(QString("Refining mesh...")))
    {
        progressEvent = new QProgressEvent(QProgressEvent_None,0,9999,0,"",QProgressEvent_Update,0,9999,8,"Tetgen meshing running on disk",this);
    }
    if(QString::fromStdString(message).contains(QString("Optimizing mesh...")))
    {
        progressEvent = new QProgressEvent(QProgressEvent_None,0,9999,0,"",QProgressEvent_Update,0,9999,9,"Tetgen meshing running on disk",this);
    }
    if(QString::fromStdString(message).contains(QString("Writing nodes.")))
    {
        progressEvent = new QProgressEvent(QProgressEvent_None,0,9999,0,"",QProgressEvent_Update,0,9999,10,"Tetgen meshing running on disk",this);
    }
    if(QString::fromStdString(message).contains(QString("Writing elements.")))
    {
        progressEvent = new QProgressEvent(QProgressEvent_None,0,9999,0,"",QProgressEvent_Update,0,9999,11,"Tetgen meshing running on disk",this);
    }
    if(QString::fromStdString(message).contains(QString("Writing faces.")))
    {
        progressEvent = new QProgressEvent(QProgressEvent_None,0,9999,0,"",QProgressEvent_Update,0,9999,12,"Tetgen meshing running on disk",this);
    }
    if(QString::fromStdString(message).contains(QString("Writing edges.")))
    {
        progressEvent = new QProgressEvent(QProgressEvent_None,0,9999,0,"",QProgressEvent_Update,0,9999,13,"Tetgen meshing running on disk",this);
    }

    if(progressEvent!=Q_NULLPTR)
    {
        QApplication::postEvent(myProgressIndicator,progressEvent);
        QApplication::processEvents();
    }
}

//! ----------------------------------------------------------------------------------------
//! function: performOnDisk1
//! details:  the array of face datasources is used as an input for building the Tetgen PLC
//!           the last three arguments are the returning values
//! ----------------------------------------------------------------------------------------
int TetgenMesher::performOnDisk1(const NCollection_Array1<occHandle(Ng_MeshVS_DataSourceFace>) &arrayOfFaceDS,
                                 int bodyIndex,
                                 occHandle(Ng_MeshVS_DataSource3D) &tetgenMesh3D,
                                 occHandle(Ng_MeshVS_DataSource2D) &tetgenMesh2D,
                                 NCollection_Array1<occHandle(Ng_MeshVS_DataSourceFace)> &tetgenArrayOfFaceDS,
                                 int done)
{
    cout<<"TetgenMesher::performOnDisk1()->____function called____"<<endl;

    //! -----------------------
    //! retrieve the body name
    //! -----------------------
    cout<<myMeshDB->MapOfBodyNames.value(bodyIndex).toStdString()<<endl;
    QString bodyName = myMeshDB->MapOfBodyNames.value(bodyIndex)+QString("_%1").arg(bodyIndex);

    //! ----------------------------------------------
    //! create the Tetgen directory for support files
    //! ----------------------------------------------
    bool clearPreviousContent = true;
    this->createTetgenSupportFilesDir(clearPreviousContent);

    //! ----------------------------
    //! send an Init progress event
    //! ----------------------------
    QProgressEvent *progressEvent = new QProgressEvent(QProgressEvent_None,0,9999,0,"",
                                                       QProgressEvent_Init,0,13,0,"Tetgen meshing running on disk",this);
    QApplication::postEvent(myProgressIndicator,progressEvent);
    QApplication::processEvents();

    //! ----------------------------------------------------------------
    //! the PLC is written using the two separate files <bodyName>.node
    //! and <bodyName>.poly00
    //! ----------------------------------------------------------------
    QString nodeFilePath = myTetgenSupportFilesDir+"/"+bodyName+".node";
    QString polyFilePath = myTetgenSupportFilesDir+"/"+bodyName+".poly";

    cout<<"TetgenMesher::performOnDisk1()->____reading node file: "<<nodeFilePath.toStdString()<<"____"<<endl;
    cout<<"TetgenMesher::performOnDisk1()->____reading poly file: "<<polyFilePath.toStdString()<<"____"<<endl;

    //! -----------------------------------------------------
    //! build the PLC: generate the .node and the .poly file
    //! uses the external tool MeshTools
    //! if a the progess indicator is passed as argument the
    //! process can be stopped by the user
    //! -----------------------------------------------------
    bool isPLCDone = MeshTools::buildPLC(arrayOfFaceDS,nodeFilePath,polyFilePath,myProgressIndicator);
    if(!isPLCDone) return -1;

    //! ---------------------------
    //! configure and start Tetgen
    //! ---------------------------
    QString program = QString(TETGEN_PATH);
    tetgenProcess->setProgram(program);

    //! ------------------
    //! set the arguments
    //! ------------------
    QList<QString> arguments;
    //arguments<<QString::fromLatin1(mySwitches)<<polyFilePath;
    arguments<<QString("-")+QString::fromLatin1(mySwitches)<<polyFilePath;

    cout<<"TetgenMesher::performOnDisk1()->____arguments: "<<mySwitches<<" "<<polyFilePath.toStdString()<<endl;

    tetgenProcess->setArguments(arguments);

    //! --------------------------
    //! set the working directory
    //! --------------------------
    cout<<"TetgenMesher::performOnDisk()->____setting the output directory: "<<myTetgenSupportFilesDir.toStdString()<<"____"<<endl;
    tetgenProcess->setWorkingDirectory(myTetgenSupportFilesDir);

    cout<<"TetgenMesher::performOnDisk1()->____Tetgen meshing started____"<<endl;

    //! --------------------------------------------------------------
    //! this connection allows the creation of Tetgen progress events
    //! --------------------------------------------------------------
    disconnect(tetgenProcess,SIGNAL(readyReadStandardOutput()),this,SLOT(redirectTetgenOutput()));
    connect(tetgenProcess,SIGNAL(readyReadStandardOutput()),this,SLOT(redirectTetgenOutput()));

    tetgenProcess->start();
    int exitCode = tetgenProcess->exitCode();
    tetgenProcess->waitForFinished(-1);

    //! ------------------------------------
    //! generate the face mesh data sources
    //! ------------------------------------
    bool isDone = this->retrieveMeshDataSourcesFromDisk(bodyIndex,tetgenMesh3D,tetgenMesh2D,tetgenArrayOfFaceDS,done);
    if(isDone == false) exitCode = -2;
    return exitCode;
}

//! ----------------------------------------------------------------------------
//! function: performOnDisk1
//! details:  this is an overloaded function. It returns only the volume mesh.
//!           It is used when generating prismatic meshes at the boundary, and
//!           at only the interior the volume mesh is of interest
//! ----------------------------------------------------------------------------
int TetgenMesher::performOnDisk1(const NCollection_Array1<occHandle(Ng_MeshVS_DataSourceFace)> &arrayOfFaceDS,
                                 int bodyIndex,
                                 occHandle(Ng_MeshVS_DataSource3D) &tetgenMesh3D)
{
    cout<<"TetgenMesher::performOnDisk1()->____(overload) function called____"<<endl;

    //! -----------------------
    //! retrieve the body name
    //! -----------------------
    QString bodyName = myMeshDB->MapOfBodyNames.value(bodyIndex)+QString("_%1").arg(bodyIndex);

    //! ----------------------------------------------
    //! create the Tetgen directory for support files
    //! ----------------------------------------------
    bool clearPreviousContent = true;
    this->createTetgenSupportFilesDir(clearPreviousContent);

    //! ----------------------------------------------------------------
    //! the PLC is written using the two separate files <bodyName>.node
    //! and <bodyName>.poly
    //! ----------------------------------------------------------------
    QString nodeFilePath = myTetgenSupportFilesDir+"/"+bodyName+".node";
    QString polyFilePath = myTetgenSupportFilesDir+"/"+bodyName+".poly";

    cout<<"TetgenMesher::performOnDisk1()->____reading node file: "<<nodeFilePath.toStdString()<<"____"<<endl;
    cout<<"TetgenMesher::performOnDisk1()->____reading poly file: "<<polyFilePath.toStdString()<<"____"<<endl;

    //! -----------------------------------------------------
    //! build the PLC: generate the .node and the .poly file
    //! uses the external tool MeshTools
    //! -----------------------------------------------------
    bool isPLCDone = MeshTools::buildPLC(arrayOfFaceDS,nodeFilePath,polyFilePath);
    if(isPLCDone==false)
    {
        cout<<"TetgenMesher::performOnDisk1()->____error in generating the PLC____"<<endl;
        return -1;
    }

    //! ---------------------------
    //! configure and start Tetgen
    //! ---------------------------
    QString program = QString(TETGEN_PATH);
    cout<<"TetgenMesher::performOnDisk1()->____program: "<<program.toStdString()<<"____"<<endl;
    tetgenProcess->setProgram(program);

    //! ------------------
    //! set the arguments
    //! ------------------
    QList<QString> arguments;
    //arguments<<QString::fromLatin1(mySwitches)<<polyFilePath;
    arguments<<QString("-")+QString::fromLatin1(mySwitches)<<polyFilePath;

    cout<<"TetgenMesher::performOnDisk1()->____arguments: "<<mySwitches<<" "<<polyFilePath.toStdString()<<endl;

    tetgenProcess->setArguments(arguments);

    //! --------------------------
    //! set the working directory
    //! --------------------------
    tetgenProcess->setWorkingDirectory(myTetgenSupportFilesDir);

    cout<<"TetgenMesher::performOnDisk()->____setting the output directory: "<<myTetgenSupportFilesDir.toStdString()<<"____"<<endl;
    cout<<"TetgenMesher::performOnDisk1()->____Tetgen meshing started____"<<endl;

    disconnect(tetgenProcess,SIGNAL(readyReadStandardOutput()),this,SLOT(redirectTetgenOutput()));
    connect(tetgenProcess,SIGNAL(readyReadStandardOutput()),this,SLOT(redirectTetgenOutput()));
    tetgenProcess->start();
    int exitCode = tetgenProcess->exitCode();
    tetgenProcess->waitForFinished(-1);

    if(exitCode!=0)
    {
        cout<<"____process terminated by the user____"<<endl;
        return exitCode;
    }

    //! ------------------------------------
    //! generate the volume mesh datasource
    //! ------------------------------------
    this->retrieveVolumeMeshDataSourceOnly(bodyIndex,tetgenMesh3D);
    return exitCode;
}

//! -----------------------
//! function: readNodeFile
//! details:
//! -----------------------
QList<QList<double>> TetgenMesher::readNodeFile(const QString &nodeFileName)
{
    QList<QList<double>> listOfAllNodes;
    if(nodeFileName.isEmpty()) return QList<QList<double>>();
    FILE *nodeFile = fopen(nodeFileName.toStdString().c_str(),"r");
    if(nodeFile==NULL) return QList<QList<double>>();

    cout<<"TetgenMesher::readNodeFile()->____Start reading node file: "<<nodeFileName.toStdString()<<"____"<<endl;

    //! ----------------------------------------------------------------------------------------
    //! First line: <# of points> <dimension (3)> <# of attributes> <boundary markers (0 or 1)>
    //! ----------------------------------------------------------------------------------------
    int NbPoints;
    int dimension;
    int NbAttributes;
    int boundaryMarkerTag;
    int pointNumber;
    double x,y,z;
    int tag;
    fscanf(nodeFile,"%d%d%d%d",&NbPoints,&dimension,&NbAttributes,&boundaryMarkerTag);

    cout<<"TetgenMesher::readNodeFile()->____Start reading mesh points___"<<endl;

    for(int line=1; line<=NbPoints; line++)
    {
        fscanf(nodeFile,"%d%lf%lf%lf%d",&pointNumber,&x,&y,&z,&tag);
        //if(tag!=0)
        //if(tag==faceNr)
        {
            QList<double> aPoint;
            aPoint<<x<<y<<z;
            listOfAllNodes<<aPoint;
        }
    }

    cout<<"TetgenMesher::readNodeFile()->____closing the .node file____"<<endl;

    fclose(nodeFile);
    return listOfAllNodes;
}

//! -----------------------
//! function: readFaceFile
//! details:
//! -----------------------
QList<QList<int>> TetgenMesher::readFaceFile(const QString &faceFileName)
{
    QList<QList<int>> listOfFaces;
    if(faceFileName.isNull()) return QList<QList<int>>();

    FILE *faceFile = fopen(faceFileName.toStdString().c_str(),"r");
    if(faceFile==NULL)
    {
        cerr<<"TetgenMesher::readFaceFile()->____the .face file cannot be opened____"<<endl;
        return QList<QList<int>>();
    }

    cout<<"TetgenMesher::readFaceFile->()___Start reading face file: "<<faceFileName.toStdString()<<"____"<<endl;

    //! --------------------
    //! read the first line
    //! --------------------
    int NbElements;
    int boundaryTag;
    fscanf(faceFile,"%d%d",&NbElements,&boundaryTag);

    for(int line=1; line<=NbElements; line++)
    {
        int number;
        int n1,n2,n3,n4,n5,n6;
        int faceTag;
        QList<int> faceDefinition;  //! list of nodes followed by face tag number

        if(5==fscanf(faceFile,"%d%d%d%d%d",&number,&n1,&n2,&n3,&faceTag))
        {
            faceDefinition<<number<<n1<<n2<<n3<<faceTag;
        }
        else if(8==fscanf(faceFile,"%d%d%d%d%d%d%d%d",&number,&n1,&n2,&n3,&n4,&n5,&n6,&faceTag))
        {
            faceDefinition<<number<<n1<<n2<<n3<<n3<<n4<<n4<<n6<<faceTag;
        }
        listOfFaces<<faceDefinition;
    }

    cout<<"TetgenMesher::readFaceFile->()____closing the .face file____"<<endl;

    //! ---------------
    //! close the file
    //! ---------------
    fclose(faceFile);
    return listOfFaces;
}


//! --------------------------------------------------------------------------------------
//! function: readNodeFile
//! details:  return value: the nodes ot the whole 3D mesh are stored into the .node file
//! --------------------------------------------------------------------------------------
bool TetgenMesher::readNodeFile(const QString &nodeFileName, QMap<int,mesh::meshPoint> &indexedMapOfAllMeshPoints)
{
    cout<<"TetgenMesher::readNodeFile()->____Start reading node file: "<<nodeFileName.toStdString()<<"____"<<endl;

    //! -----------------------------------------------------------
    //! the .node file contains the nodes of the whole volume mesh
    //! -----------------------------------------------------------
    FILE *nodeFile = fopen(nodeFileName.toStdString().c_str(),"r");
    if(nodeFile==NULL)
    {
        cout<<"TetgenMesher::readNodeFile()->____the .node file cannot be opened____"<<endl;
        return false;
    }

    //! ----------------------------------------------------------------------------------------
    //! First line: <# of points> <dimension (3)> <# of attributes> <boundary markers (0 or 1)>
    //! ----------------------------------------------------------------------------------------
    int NbPoints;
    int dimension;
    int NbAttributes;
    int boundaryMarkerTag;
    int pointNumber;
    double x,y,z;
    int tag;
    fscanf(nodeFile,"%d%d%d%d",&NbPoints,&dimension,&NbAttributes,&boundaryMarkerTag);

    //! -----------------
    //! console messages
    //! -----------------
    cout<<"TetgenMesher::readNodeFile()->____Start reading mesh points___"<<endl;

    //! ------------------------------------
    //! "line" is actually the globalNodeID
    //! ------------------------------------
    for(int line=1; line<=NbPoints; line++)
    {
        fscanf(nodeFile,"%d%lf%lf%lf%d",&pointNumber,&x,&y,&z,&tag);
        mesh::meshPoint aMeshPoint(x,y,z,line);
        indexedMapOfAllMeshPoints.insert(line,aMeshPoint);
    }

    cout<<"TetgenMesher::readNodeFile()->____mesh points read____"<<endl;
    fclose(nodeFile);
    return true;
}

//! -------------------------------
//! function: setStopButtonEnabled
//! details:
//! -------------------------------
void TetgenMesher::setStopButtonEnabled(bool isEnabled)
{
    QProgressEvent *pe;
    if(isEnabled) pe = new QProgressEvent(QProgressEvent_EnableStop,0,0,0,"",QProgressEvent_None,0,0,0,"Tetgen meshing running on disk");
    else pe = new QProgressEvent(QProgressEvent_DisableStop,0,0,0,"",QProgressEvent_None,0,0,0,"Tetgen meshing running on disk");
    QApplication::postEvent(myProgressIndicator,pe);
}
