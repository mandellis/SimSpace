//! ----------------
//! custom includes
//! ----------------
#include "ng_mesher2.h"
#include "ng_meshvs_datasource2d.h"
#include "ng_meshvs_datasource3d.h"
#include "ng_meshvs_datasourceface.h"
#include "src/main/mydefines.h"
#include "src/utils/tools.h"
#include "src/utils/ccout.h"
#include "qprogressevent.h"

//! ----
//! OCC
//! ----
#include <Bnd_Box.hxx>
#include <BRepBndLib.hxx>
#include <BRepGProp.hxx>
#include <BRep_Tool.hxx>
#include <GProp_GProps.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Edge.hxx>
#include <TopExp_Explorer.hxx>
#include <TopExp.hxx>
#include <TopAbs_ShapeEnum.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TCollection_AsciiString.hxx>
#include <TColStd_ListOfInteger.hxx>
#include <TColStd_HArray2OfInteger.hxx>
#include <TColStd_HArray2OfReal.hxx>

//! ---
//! Qt
//! ---
#include <QMessageBox>
#include <QString>
#include <QDir>

//! ----
//! C++
//! ----
#include <iostream>

//! nglib
#include <ngexception.hpp>

const double theDeflection = 0.01;

//! -----------------------
//! function: trans_func_N
//! details:
//! -----------------------
void trans_func_N(unsigned int u, EXCEPTION_POINTERS* pExp)
{
    Q_UNUSED(pExp)
    Q_UNUSED(u)
    printf( "In trans_func.\n" );
    throw SE_Exception();
}

//! --------------------------------
//! function: initDefaultParameters
//! details:
//! --------------------------------
void NetgenMesher::initDefaultParameters()
{
    //! Netgen - nglib meshing parameters
    mp.uselocalh = 1;                   // Enable/Disable usage of local mesh size modifiers
    mp.maxh = 100;                      // Maximum allowed mesh size
    mp.minh= 1;                         // Minimum allowed ->global<- mesh size allowed
    mp.fineness = 0.5;                  // Mesh density: 0...1 (0 => coarse; 1 => fine)
    mp.grading = 0.3;                   // Mesh grading: 0...1 (0 => uniform mesh; 1 => aggressive local grading)
    mp.elementsperedge=4;               // Number of elements to generate per edge
    mp.elementspercurve=2;              // Elements to generate per curvature radius
    mp.closeedgeenable=0;               // Enable/Disable mesh refinement at close edges
    mp.closeedgefact=0;                 // Factor to use for refinement at close edges (larger => finer)
    mp.minedgelenenable=0;              // Enable/Disable user defined minimum edge length for edge subdivision
    mp.minedgelen=1e-2;                 // Minimum edge length to use while subdividing the edges (default = 1e-4)
    mp.second_order=0;                  // Generate second-order surface and volume elements
    mp.quad_dominated=0;                // Creates a Quad-dominated mesh
    mp.optsurfmeshenable=1;             // Enable/Disable automatic surface mesh optimization
    mp.optvolmeshenable=1;              // Enable/Disable automatic volume mesh optimization
    mp.optsteps_2d=2;                   // Number of optimize steps to use for 2-D mesh optimization
    mp.optsteps_3d=2;                   // Number of optimize steps to use for 3-D mesh optimization
    mp.check_overlapping_boundary = 1;
    mp.check_overlap = 1;
    mp.invert_tets = 0;
    mp.invert_trigs = 0;
}

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
NetgenMesher::NetgenMesher(QObject *parent):QObject(parent)
{
    //! ------------------------
    //! init default parameters
    //! ------------------------
    this->initDefaultParameters();

    //! -----------------------
    //! the progress indicator
    //! -----------------------
    myProgressIndicator = Q_NULLPTR;
}

//! ------------------------
//! function: constructor I
//! details:
//! ------------------------
NetgenMesher::NetgenMesher(meshParam meshParams, QObject *parent):QObject(parent)
{
    //! ---------------------
    //! init the netgen mesh
    //! ---------------------
    Ng_Init();
    myNgMesh = Ng_NewMesh();

    //! ------------------------
    //! init default parameters
    //! ------------------------
    this->initDefaultParameters();

    //! -----------------------
    //! the progress indicator
    //! -----------------------
    myProgressIndicator = Q_NULLPTR;

    //! -----------------------
    //! set meshing parameters
    //! -----------------------
    mp.grading = meshParams.grading;
    mp.minh= meshParams.minElementSize;
    mp.maxh= meshParams.maxElementSize;
}

//! -------------------------
//! function: constructor II
//! details:
//! -------------------------
NetgenMesher::NetgenMesher(TopoDS_Shape theShape, meshParam meshParams, QObject *parent):
    QObject(parent),
    myTopoDS_Shape(theShape)
{
    cout<<"NetgenMesher::NetgenMesher()->____constructor II called____"<<endl;

    //! ---------------------
    //! init the netgen mesh
    //! ---------------------
    Ng_Init();
    myNgMesh = Ng_NewMesh();

    //! ------------------------
    //! init default parameters
    //! ------------------------
    this->initDefaultParameters();

    //! -----------------------
    //! the progress indicator
    //! -----------------------
    myProgressIndicator = Q_NULLPTR;

    //! -----------------------------
    //! TopoDS_Shape to Ng_OCC_Shape
    //! -----------------------------
    Ng_OCC_Shape occShape = (Ng_OCC_Shape)(&myTopoDS_Shape);

    //! --------------------------------
    //! Ng_OCC_Shape && Ng_OCC_Geometry
    //! --------------------------------
    myNg_OCC_geometry = Ng_UseOCCGeometry(occShape);

    //! -------------------
    //! set edge mesh size
    //! -------------------
    QMap<int, double> edgeSizeMap = meshParams.edgeSize;
    for (QMap<int, double>::iterator anIterE = edgeSizeMap.begin(); anIterE!= edgeSizeMap.end(); ++anIterE)
    {
        int aKey = anIterE.key();
        double aValue = anIterE.value();
        Ng_OCC_SetEdgeSize(myNg_OCC_geometry, aKey, aValue);
    }

    //! -------------------------
    //! set local face mesh size
    //! -------------------------
    QMap<int, double> faceSizeMap = meshParams.faceSize;
    for (QMap<int, double>::iterator anIterF = faceSizeMap.begin(); anIterF!= faceSizeMap.end(); ++anIterF)
    {
        int aKey = anIterF.key();
        double aValue = anIterF.value();
        Ng_OCC_SetFaceSize(myNg_OCC_geometry, aKey, aValue);
    }

    //! ----------------------
    //! set local vertex size
    //! ----------------------
    QMap<int, double> vertexSizeMap = meshParams.vertexSize;
    QMap<int, double> pinballSizeMap = meshParams.pinballSize;

    if(!vertexSizeMap.isEmpty())
    {
        QMap<int, double>::iterator anIterV;
        QMap<int, double>::iterator anIterPin;

        for (anIterV = vertexSizeMap.begin(), anIterPin = pinballSizeMap.begin();
             anIterV!=vertexSizeMap.end() && anIterPin!=pinballSizeMap.end(); ++anIterV, ++anIterPin)
        {
            //! vertex number
            int aKey = anIterV.key();
            //! element size
            double aValue = anIterV.value();
            //!  - region of influence
            double pinball = anIterPin.value();
            TopTools_IndexedMapOfShape vmap;
            TopExp::MapShapes(theShape,TopAbs_VERTEX,vmap);
            gp_Pnt V = BRep_Tool::Pnt(TopoDS::Vertex(vmap.FindKey(aKey)));
            double C[3];
            C[0] = V.X();
            C[1] = V.Y();
            C[2] = V.Z();
            //!cout<<"Element size: "<<aValue<<" Point(x, y, z) = "<<"("<<C[0]<<", "<<C[1]<<", "<<C[2]<<")"<<endl;
            Ng_RestrictMeshSizeSphere(myNgMesh, C, pinball, aValue);
        }
    }

    //! -----------------------
    //! set meshing parameters
    //! -----------------------
    mp.grading = meshParams.grading;
    mp.minh= meshParams.minElementSize;
    mp.maxh= meshParams.maxElementSize;
    Ng_OCC_SetLocalMeshSize(myNg_OCC_geometry, myNgMesh, &mp, theDeflection);
    mp.quad_dominated = meshParams.isSurfaceQuad;
}

//! ---------------------
//! function: destructor
//! details:
//! ---------------------
NetgenMesher::~NetgenMesher()
{
    cout<<"NetgenMesher::~NetgenMesher()->____DESTRUCTOR CALLED____"<<endl;
}

//! -------------------------------
//! function: setProgressIndicator
//! details:
//! -------------------------------
void NetgenMesher::setProgressIndicator (QProgressIndicator *theProgressIndicator)
{
    myProgressIndicator = theProgressIndicator;
    //connect(myProgressIndicator,SIGNAL(abortPressed()),this,SLOT(abort()));
}

//! -----------------------------------------------------------------------
//! function: performVolume
//! details:  requires constructor II - with shape and meshing paramenters
//! -----------------------------------------------------------------------
userMessage NetgenMesher::performVolume(occHandle(MeshVS_DataSource)& mainMeshVS_DataSource3D,
                                          occHandle(MeshVS_DataSource)& mainMeshVS_DataSource2D,
                                          bool isSecondOrder,
                                          bool isElementStraight)
{
    //! ----------------------------------------------------------------------
    //! generate first the surface mesh, but do not delete the Netgen pointer
    //! ----------------------------------------------------------------------
    bool deleteNgPointers = false;
    userMessage mr = this->performSurface(mainMeshVS_DataSource2D,isSecondOrder,isElementStraight,deleteNgPointers);
    if(mr.isDone)
    {
        try
        {
            //! generate the volume mesh
            _set_se_translator(trans_func_N);
            this->SEFunc_GenerateVolume(myNgMesh, &mp);
        }
        catch(netgen::NgException)
        {
            //! for external use
            throw NG_ERROR;

            this->deleteNgPointers();
            userMessage mr(false,"Netgen: the volume mesh cannot be generated");
            return mr;
        }
        catch(SE_Exception)
        {
            //! for external use
            throw SE_Exception();

            this->deleteNgPointers();
            userMessage mr(false,"Netgen has generated an system error while generating the volume mesh");
            return mr;
        }

        //! ----------------------
        //! check if second order
        //! ----------------------
        if(isSecondOrder==true)
        {
            if(isElementStraight) Ng_Generate_SecondOrder(myNgMesh);
            else Ng_OCC_Generate_SecondOrder(myNg_OCC_geometry, myNgMesh);
        }

        //! ----------------------------
        //! resulting mesh data sources
        //! ----------------------------
        mainMeshVS_DataSource2D = new Ng_MeshVS_DataSource2D(myNgMesh);
        mainMeshVS_DataSource3D = new Ng_MeshVS_DataSource3D(myNgMesh);

        if(mainMeshVS_DataSource2D.IsNull())
        {
            this->deleteNgPointers();
            userMessage mr(true, "Cannot convert Netgen pointer into surface mesh data source");
            return mr;
        }

        if(mainMeshVS_DataSource3D.IsNull())
        {
            this->deleteNgPointers();
            userMessage mr(true, "Cannot convert Netgen pointer into volume mesh data source");
            return mr;
        }

        //! ------------------------------------
        //! generate the face mesh data sources
        //! ------------------------------------
#ifdef GENERATE_FACE_MESH_DATASOURCES
        this->generateFaceMeshDataSources();
#endif
        //! ---------------------------
        //! delete the Netgen pointers
        //! ---------------------------
        this->deleteNgPointers();

        userMessage mr(true, "Volume mesh generated");
        return mr;
    }
    return mr;
}

//! -------------------------
//! function: performSurface
//! details:
//! -------------------------
userMessage NetgenMesher::performSurface(occHandle(MeshVS_DataSource)& mainMeshVS_DataSource2D,
                                           bool isSecondOrder,
                                           bool isElementStraight,
                                           bool deleteNgPointers)
{
    //! -----------------------
    //! generate the edge mesh
    //! -----------------------
    Ng_Result ng_res;
    ng_res=Ng_OCC_GenerateEdgeMesh(myNg_OCC_geometry,myNgMesh,&mp);

    try
    {
        _set_se_translator(trans_func_N);
        this->SEFunc_GenerateSurface(&mp);
    }
    catch(SE_Exception e)
    {
        //! --------------------------------------------------------------------
        //! unhandled system exception generated by the surface meshing process
        //! --------------------------------------------------------------------
        this->deleteNgPointers();
        userMessage mr(false,"Netgen system error while generating the surface mesh");
        return mr;
    }
    catch(netgen::NgException)
    {
        //! --------------
        //! handled error
        //! --------------
        this->deleteNgPointers();
        userMessage mr(false,"Netgen error while generating the surface mesh");
        return mr;
    }

    //! ----------------------
    //! generate second order
    //! ----------------------
    if(isSecondOrder==true)
    {
        if(isElementStraight) Ng_Generate_SecondOrder(myNgMesh);
        else Ng_OCC_Generate_SecondOrder(myNg_OCC_geometry, myNgMesh);
    }

    mainMeshVS_DataSource2D = new Ng_MeshVS_DataSource2D(myNgMesh);
    if(mainMeshVS_DataSource2D.IsNull())
    {
        this->deleteNgPointers();
        userMessage mr(false,"Cannot convert Netgen pointer into surface mesh data source");
        return mr;
    }

    //! ------------------------------------
    //! generate the face mesh data sources
    //! ------------------------------------
#ifdef GENERATE_FACE_MESH_DATASOURCES
    this->generateFaceMeshDataSources();
#endif

    //! ----------------------------------------
    //! delete the Netgen pointers if requested
    //! ----------------------------------------
    if(deleteNgPointers) this->deleteNgPointers();

    userMessage mr(true,"Surface mesh sucessfully generated");
    return mr;
}

//! ----------------------------------
//! function: SEFunct_generateSurface
//! details:
//! ----------------------------------
void NetgenMesher::SEFunc_GenerateSurface(Ng_Meshing_Parameters *mp)
{
    __try
    {
        cout<<"NetgenMesher::SEFunc_GenerateSurface()->____try____"<<endl;
        Ng_OCC_GenerateSurfaceMesh(myNg_OCC_geometry,myNgMesh,mp);
    }
    __finally
    {
        cout<<"NetgenMesher::SEFunc_GenerateSurface()->____finally____"<<endl;
    }
}

//! --------------------------------------------------------------------
//! function: performVolume
//! details:  generate a volume mesh from a netgen surface mesh
//!           generated externally. The input Netgen mesh must contains
//!           face descriptors
//! --------------------------------------------------------------------
userMessage NetgenMesher::performVolume(meshParam meshingParameters,
                                          Ng_Mesh *mesh,
                                          occHandle(MeshVS_DataSource) &mainMeshVS_DataSource3D)
{
    cout<<"NetgenMesher::performVolume()->____function called____"<<endl;

    //! ------------------------------------------------------------------------
    //! change the default setting: correct surface elements generated by EMESH
    //! ------------------------------------------------------------------------
    mp.invert_trigs = 1;
    mp.check_overlap = 1;

    //! ---------------------------------------------
    //! change the default Netgen meshing parameters
    //! ---------------------------------------------
    mp.grading = meshingParameters.grading;
    mp.uselocalh = 1;
    mp.minh= meshingParameters.minElementSize;
    mp.maxh= meshingParameters.maxElementSize;
    mp.optsteps_3d = 3;

    try
    {
        _set_se_translator(trans_func_N);
        this->SEFunc_GenerateVolume(mesh, &mp);
    }
    catch(SE_Exception e)
    {
        cout<<"____Netgen has generated a system error while generating volume mesh____"<<endl;
        userMessage mr(false,"Netgen has generated a system error while generating volume mesh");
        return mr;
    }
    catch(netgen::NgException)
    {
        cout<<"____Netgen error while generating volume mesh____"<<endl;
        mainMeshVS_DataSource3D = new Ng_MeshVS_DataSource3D();
        userMessage mr(false,"Netgen error while generating volume mesh");
        return mr;
    }

    //! -------------------------------------------
    //! try generating the volume mesh data source
    //! -------------------------------------------
    mainMeshVS_DataSource3D = new Ng_MeshVS_DataSource3D(mesh);
    if(mainMeshVS_DataSource3D.IsNull())
    {
        cout<<"____Cannot convert the netgen pointer into a volume mesh data source____"<<endl;
        mainMeshVS_DataSource3D = new Ng_MeshVS_DataSource3D();
        userMessage mr(false,"Cannot convert the netgen pointer into a volume mesh data source");
        return mr;
    }

    //! ---------------------------
    //! delete the Netgen pointers
    //! ---------------------------
    //this->deleteNgPointers();

    userMessage mr(true,"Netgen volume mesh successfully generated");
    return mr;
}

//! ---------------------------------
//! function: SEFunct_generateVolume
//! details:
//! ---------------------------------
void NetgenMesher::SEFunc_GenerateVolume(Ng_Mesh *NgMesh, Ng_Meshing_Parameters *mp)
{
    __try
    {
        Ng_GenerateVolumeMesh(NgMesh, mp);
        cout<<"____in try: volume meshing done____"<<endl;
    }
    __finally
    {
        cout<<"____in finally____"<<endl;
    }
}

//! ---------------------------
//! function: deleteNgPointers
//! details:
//! ---------------------------
void NetgenMesher::deleteNgPointers()
{
    if(myNgMesh!=NULL)
    {
        cout<<"NetgenMesher::deleteNgPointers()->____delete Netgen mesh pointer____"<<endl;
        Ng_DeleteMesh(myNgMesh);
        myNgMesh = NULL;
    }
    if(myNg_OCC_geometry!=NULL)
    {
        cout<<"NetgenMesher::deleteNgPointers()->____delete Netgen OCC geometry____"<<endl;
        Ng_OCC_DeleteGeometry(myNg_OCC_geometry);
        myNg_OCC_geometry = NULL;
    }
}

//! ----------------------------------------------
//! function: generateFaceMeshDataSources
//! details:  generate the face mesh data sources
//! ----------------------------------------------
void NetgenMesher::generateFaceMeshDataSources(int done)
{
    TopTools_IndexedMapOfShape faceMap;
    TopExp::MapShapes(myTopoDS_Shape,TopAbs_FACE,faceMap);
    int NbFaces = faceMap.Extent();

    //! --------------------------------
    //! init the secondary progress bar
    //! --------------------------------
    QProgressEvent *pe = new QProgressEvent();
    pe->setVal(done);
    pe->setMessage("Generating the face mesh data sources");
    pe->setAction1(QProgressEventAction::QProgressEvent_Init);
    pe->setRange1(0, NbFaces);
    if(myProgressIndicator) QApplication::postEvent(myProgressIndicator,pe);

    for(int faceNr=1; faceNr<=NbFaces; faceNr++)
    {
        const occHandle(Ng_MeshVS_DataSourceFace) &curFaceMeshDS = new Ng_MeshVS_DataSourceFace(myNgMesh, faceNr);
        faceMeshDS.insert(faceNr, curFaceMeshDS);

        //! -------------------------------------------
        //! post update events on face mesh generation
        //! -------------------------------------------
        if(faceNr%1==0)
        {
            QProgressEvent *pe = new QProgressEvent();
            pe->setVal(done);
            pe->setMessage(QString("Face mesh %1").arg(faceNr));
            pe->setVal1(faceNr);
            if(myProgressIndicator) QApplication::postEvent(myProgressIndicator,pe);
        }
    }
}

//! --------------------
//! function: configure
//! details:
//! --------------------
void NetgenMesher::configure(const TopoDS_Shape &shape, meshParam theMeshingParams)
{
    //! ---------------------
    //! init the netgen mesh
    //! ---------------------
    Ng_Init();
    myNgMesh = Ng_NewMesh();

    //! ---------------
    //! init the shape
    //! ---------------
    myTopoDS_Shape = shape;

    //! -----------------------------
    //! TopoDS_Shape to Ng_OCC_Shape
    //! -----------------------------
    Ng_OCC_Shape occShape = (Ng_OCC_Shape)(&myTopoDS_Shape);

    //! --------------------------------
    //! Ng_OCC_Shape && Ng_OCC_Geometry
    //! --------------------------------
    myNg_OCC_geometry = Ng_UseOCCGeometry(occShape);

    //! -------------------
    //! set edge mesh size
    //! -------------------
    QMap<int, double> edgeSizeMap = theMeshingParams.edgeSize;
    for (QMap<int, double>::iterator anIterE = edgeSizeMap.begin(); anIterE!= edgeSizeMap.end(); ++anIterE)
    {
        int aKey = anIterE.key();
        double aValue = anIterE.value();
        Ng_OCC_SetEdgeSize(myNg_OCC_geometry, aKey, aValue);
    }

    //! -------------------------
    //! set local face mesh size
    //! -------------------------
    QMap<int, double> faceSizeMap = theMeshingParams.faceSize;
    for (QMap<int, double>::iterator anIterF = faceSizeMap.begin(); anIterF!= faceSizeMap.end(); ++anIterF)
    {
        int aKey = anIterF.key();
        double aValue = anIterF.value();
        Ng_OCC_SetFaceSize(myNg_OCC_geometry, aKey, aValue);
    }

    //! ----------------------
    //! set local vertex size
    //! ----------------------
    QMap<int, double> vertexSizeMap = theMeshingParams.vertexSize;
    QMap<int, double> pinballSizeMap = theMeshingParams.pinballSize;

    if(!vertexSizeMap.isEmpty())
    {
        QMap<int, double>::iterator anIterV;
        QMap<int, double>::iterator anIterPin;

        for (anIterV = vertexSizeMap.begin(), anIterPin = pinballSizeMap.begin();
             anIterV!=vertexSizeMap.end() && anIterPin!=pinballSizeMap.end(); ++anIterV, ++anIterPin)
        {
            //! vertex number
            int aKey = anIterV.key();

            //! element size
            double aValue = anIterV.value();

            //!  - region of influence
            double pinball = anIterPin.value();
            TopTools_IndexedMapOfShape vmap;
            TopExp::MapShapes(myTopoDS_Shape,TopAbs_VERTEX,vmap);
            gp_Pnt V = BRep_Tool::Pnt(TopoDS::Vertex(vmap.FindKey(aKey)));
            double C[3];
            C[0] = V.X();
            C[1] = V.Y();
            C[2] = V.Z();
            //!cout<<"Element size: "<<aValue<<" Point(x, y, z) = "<<"("<<C[0]<<", "<<C[1]<<", "<<C[2]<<")"<<endl;
            Ng_RestrictMeshSizeSphere(myNgMesh, C, pinball, aValue);
        }
    }

    //! -----------------------
    //! set meshing parameters
    //! -----------------------
    mp.grading = theMeshingParams.grading;
    mp.minh= theMeshingParams.minElementSize;
    mp.maxh= theMeshingParams.maxElementSize;

    const double theDeflection = 0.01;
    Ng_OCC_SetLocalMeshSize(myNg_OCC_geometry, myNgMesh, &mp, theDeflection);
    mp.quad_dominated = theMeshingParams.isSurfaceQuad;

    /*
    //! ------------------------
    //! test new nglib function
    //! ------------------------
    Ng_Meshing_Parameters tmpmp;
    tmpmp.grading = 0.2;
    tmpmp.minh = 0.1;
    tmpmp.maxh = 5;
    tmpmp.uselocalh = 1;
    std::map<int,double> ems;
    std::map<int,double> fms;
    Ng_OCC_MeshSizingOnEdges meshSizingOnEdges = (Ng_OCC_MeshSizingOnEdges)(&ems);
    Ng_OCC_MeshSizingOnFaces meshSizingOnFaces = (Ng_OCC_MeshSizingOnFaces)(&fms);
    void *occShape1 = (void*)(&myTopoDS_Shape);

    std::vector<Ng_meshPoint> meshPoints;
    void *meshPointsP = (void*)(&meshPoints);
    std::vector<Ng_surfaceMeshElement> meshElements;
    void *meshElementsP = (void*)(&meshElements);
    std::map<int,std::vector<int>> faceMeshDefinition;
    void *fmd = (void*)(&faceMeshDefinition);

    bool isDone = false;
    isDone = Ng_OCC_buildSurfaceMesh(occShape1,
                                     &tmpmp,
                                     &meshSizingOnEdges,
                                     &meshSizingOnFaces,
                                     meshPointsP,
                                     meshElementsP,
                                     fmd);

    //! write into a file
    if(isDone)
    {
        FILE *f = fopen("D:/netgenMeshTest.txt","w");
        fprintf(f,"Nodes: %d\n",meshPoints.size());
        for(int n=0; n<meshPoints.size(); n++)
        {
            const Ng_meshPoint &P = meshPoints.at(n);
            fprintf(f,"%d\t%lf\t%lf\t%lf\n",P.nodeID,P.x,P.y,P.z);
        }
        fprintf(f,"Elements: %d\n",meshElements.size());
        for(int n=0; n<meshElements.size(); n++)
        {
            const Ng_surfaceMeshElement &element = meshElements.at(n);
            for(int k=0; k<element.nodeIDs.size()-1; k++) fprintf(f,"%d\t",element.nodeIDs.at(k));
            fprintf(f,"%d\n",element.nodeIDs.at(element.nodeIDs.size()-1));
        }
        fclose(f);
    }
    //! -------------------------------
    //! end testing new nglib function
    //! -------------------------------
    */
}

//! ----------------
//! function: abort
//! details:
//! ----------------
void NetgenMesher::abort()
{
    ;
}
