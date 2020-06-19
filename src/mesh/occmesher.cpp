//! custom includes
#include "occMesher.h"
#include "meshtools.h"
#include <ng_mesher2.h>
#include <ng_meshvs_datasource3d.h>
#include "ccout.h"
#include "qprogressevent.h"

//! OCC
#include <TopoDS_Shape.hxx>
#include <QMShape_Parameters.hxx>
#include <QMData_MeshParameters.hxx>
#include <MeshVS_DataSource.hxx>
#include <OMFDS_Mesh.hxx>
#include <OMFTools.hxx>
#include <OMFVS_DataSource.hxx>
#include <TopExp.hxx>
#include <TopAbs_ShapeEnum.hxx>

//! -------------------------------------
//! deviation and deflection angle [deg]
//! -------------------------------------
#define DEFLECTION 15
#define DEVIATION_ANGLE 15

//! -----------------------
//! function: trans_func_O
//! details:
//! -----------------------
void trans_func_O(unsigned int u, EXCEPTION_POINTERS* pExp)
{
    Q_UNUSED(pExp)
    Q_UNUSED(u)
    printf( "In trans_func.\n" );
    throw SE_Exception();
}

//! -------------------------
//! function: constructor II
//! details:
//! -------------------------
occMesher::occMesher(const TopoDS_Shape &theShape, meshParam meshParams, QObject *parent):QObject(parent),
    myTopoDS_Shape(theShape),
    myMeshParameters(meshParams),
    myProgressIndicator(Q_NULLPTR)
{
    const double PI = 3.1415926535897932384626433832795;

    //! ------------------------------------------------------------
    //! init the mesh parameters
    //! use Express mesh for getting autocomputed parameters, then
    //! copy values into the private member
    //! ------------------------------------------------------------
    QMShape_Parameters shapeBasedMeshingParameters;
    shapeBasedMeshingParameters.Load(myTopoDS_Shape);
    shapeBasedMeshingParameters.Compute();
    myEMeshMeshParameters = shapeBasedMeshingParameters.GetParameters();

    //! ----------------------------------
    //! override some of the autocomputed
    //! ----------------------------------
    myEMeshMeshParameters.SetMinElemSize(meshParams.minElementSize);
    myEMeshMeshParameters.SetMaxElemSize(meshParams.maxElementSize);
    myEMeshMeshParameters.SetDeflection(DEFLECTION*PI/180);
    myEMeshMeshParameters.SetDeviationAngle(DEVIATION_ANGLE*PI/180);
    myEMeshMeshParameters.SetCheckQuadInteriorForDeflection(Standard_True);
    // experimental feature according to the manual. It is not connected
    //myEMeshMeshParameters.SetCheckDelaunay(Standard_True);

    //! ---------------
    //! init the shape
    //! ---------------
    myTessellator.Init(myTopoDS_Shape);

    //! -------------------------
    //! set algorithm parameters
    //! -------------------------
    myTessellator.ToCheckSelfInter(true);
    myTessellator.ToCheckMutualInter(true);
    myTessellator.ToCheckClosed(true);
    myTessellator.SetQuadTreeMinDepth(2);
    myTessellator.SetQuadTreeMaxDepth(128);
    myTessellator.SetParallel(true);

    //! ----------------------------
    //! load the meshing parameters
    //! ----------------------------
    myTessellator.LoadParameters(myEMeshMeshParameters);
}

//! ------------------------
//! function: constructor I
//! details:
//! ------------------------
occMesher::occMesher(const TopoDS_Shape &aShape, QObject *parent):
    myTopoDS_Shape(aShape), QObject(parent), myProgressIndicator(Q_NULLPTR)
{
    //! ---------------------------
    //! init the meshing paramters
    //! ---------------------------
    QMShape_Parameters shapeBasedMeshingParameters;
    shapeBasedMeshingParameters.Load(myTopoDS_Shape);
    shapeBasedMeshingParameters.Compute();

    //! -----
    //! init
    //! -----
    myTessellator.Init(myTopoDS_Shape);

    //! -------------------------------------
    //! set the meshing algorithm parameters
    //! -------------------------------------
    myTessellator.ToCheckSelfInter(Standard_True);
    myTessellator.ToCheckMutualInter(Standard_True);
    myTessellator.ToCheckClosed(Standard_True);
    myTessellator.SetQuadTreeMinDepth(2);
    myTessellator.SetQuadTreeMaxDepth(128);
    myTessellator.SetParallel(Standard_True);

    //! mesh method for faces
    //bool isSurfaceQuad = meshParam.isSurfaceQuad;
    //bool isVolumeHexa = false;          //! unused: MeshExpress is a surface mesher
    //this->setMeshMethod(isSurfaceQuad,isVolumeHexa);

    //! ---------------------
    //! surface element type
    //! ---------------------
    myEMeshMeshParameters.SetMeshElemType(QMData_MeshParameters::Tri);
}

//! ---------------------
//! function: destructor
//! details:
//! ---------------------
occMesher::~occMesher()
{
    cout<<"occMesher::~occMesher()->____DESTRUCTOR CALLED____"<<endl;
}

//! -----------------------------------
//! function: setBodyMeshingParameters
//! details:
//! -----------------------------------
void occMesher::setBodyMeshingParameters(meshParam mp)
{
    const double PI = 3.1415926535897932384626433832795;

    //! ----------------------------------------------------------
    //! Netgen mesing parameters (for generating the volume mesh)
    //! ----------------------------------------------------------
    myMeshParameters.minElementSize = mp.minElementSize;
    myMeshParameters.maxElementSize = mp.maxElementSize;
    //myMeshParameters.grading = mp.grading;

    //! ------------------------------
    //! copy some of the autocomputed
    //! ------------------------------
    QMShape_Parameters shapeBasedMeshingParameters;
    shapeBasedMeshingParameters.Load(myTopoDS_Shape);
    shapeBasedMeshingParameters.Compute();
    myEMeshMeshParameters = shapeBasedMeshingParameters.GetParameters();

    //! ----------------------------------
    //! override some of the autocomputed
    //! ----------------------------------
    cout<<"occMesher::setBodyMeshingParameters()->____overriding autocomputed parameters____"<<endl;
    myEMeshMeshParameters.SetMinElemSize(mp.minElementSize);
    myEMeshMeshParameters.SetMaxElemSize(mp.maxElementSize);
    myEMeshMeshParameters.SetDeflection(DEFLECTION*PI/180.0);
    myEMeshMeshParameters.SetDeviationAngle(DEVIATION_ANGLE*PI/180.0);
    myEMeshMeshParameters.SetCheckDelaunay(Standard_True);
    myEMeshMeshParameters.SetCheckQuadInteriorForDeflection(Standard_True);       //!29/03/2018

    //! ----------------------------
    //! load the meshing parameters
    //! ----------------------------
    myTessellator.LoadParameters(myEMeshMeshParameters);
}

//! ------------------------
//! function: setMeshMethod
//! details:
//! ------------------------
void occMesher::setMeshMethod(int isSurfaceQuad, int isVolumeHexa)
{
    Q_UNUSED(isVolumeHexa);
    if(isSurfaceQuad)
    {
        myEMeshMeshParameters.SetMeshElemType(QMData_MeshParameters::Quad);
        cout<<"occMesher::setMeshMethod()->____set quad surface on body____"<<endl;
    }
    else
    {
        myEMeshMeshParameters.SetMeshElemType(QMData_MeshParameters::Tri);
    }
}

//! ---------------------------------
//! function: perform
//! details:  check if used... to do
//! ---------------------------------
Standard_Boolean occMesher::perform(NCollection_Handle<NCollection_List<Standard_Integer>> &badFacesIndices)
{
    //! --------------------------
    //! reload meshing parameters
    //! --------------------------
    myTessellator.LoadParameters(myEMeshMeshParameters);

    //! diagnostic - can be removed
    cout<<"occMesher::perform()->____Minimum size: "<<myEMeshMeshParameters.MinElemSize()<<"____"<<endl;
    cout<<"occMesher::perform()->____Maximum size: "<<myEMeshMeshParameters.MaxElemSize()<<"____"<<endl;
    cout<<"occMesher::perform()->____Deflection: "<<myEMeshMeshParameters.Deflection()<<"____"<<endl;
    cout<<"occMesher::perform()->____Deviation angle: "<<myEMeshMeshParameters.DeviationAngle()<<"____"<<endl;
    //! end diagnostic

    myTessellator.Perform();
    badFacesIndices = myTessellator.BadFacesIndices();

    cout<<"occMesher::perform()->____done____"<<endl;
    return myTessellator.IsDone();
}

//! ---------------------------------------------------
//! function: perform
//! details:  just for generating a visualization mesh
//! ---------------------------------------------------
Standard_Boolean occMesher::perform(double angularDeflection, double linearDeflection, double minSize, double maxSize)
{
    //! ---------------------------
    //! init the meshing paramters
    //! ---------------------------
    QMShape_Parameters shapeBasedMeshingParameters;
    shapeBasedMeshingParameters.Load(myTopoDS_Shape);
    shapeBasedMeshingParameters.Compute();

    cout<<"____Autocomputed Express Mesh parameters____"<<endl;
    cout<<"____Angular deflection: "<<shapeBasedMeshingParameters.GetParameters().DeviationAngle()<<"____"<<endl;
    cout<<"____Linear deflection: "<<shapeBasedMeshingParameters.GetParameters().Deflection()<<"____"<<endl;
    cout<<"____Minimum size: "<<shapeBasedMeshingParameters.GetParameters().MinElemSize()<<"____"<<endl;
    cout<<"____Maximum size: "<<shapeBasedMeshingParameters.GetParameters().MaxElemSize()<<"____"<<endl;

    //! ------------------
    //! copy autocomputed
    //! ------------------
    myEMeshMeshParameters = shapeBasedMeshingParameters;

    //! ------------------------------
    //! compare autocomputed and user
    //! ------------------------------
    if(angularDeflection<shapeBasedMeshingParameters.GetParameters().DeviationAngle())
    {
        cout<<"____Overriding autocomputed angular deflection____"<<endl;
        myEMeshMeshParameters.SetDeviationAngle(angularDeflection);
    }
    if(linearDeflection<shapeBasedMeshingParameters.GetParameters().Deflection())
    {
        cout<<"____Overriding autocomputed linear deflection____"<<endl;
        myEMeshMeshParameters.SetDeflection(linearDeflection);
    }
    if(minSize<shapeBasedMeshingParameters.GetParameters().MinElemSize())
    {
        cout<<"____Overriding autocomputed minimum element size____"<<endl;
        myEMeshMeshParameters.SetMinElemSize(minSize);
    }
    if(maxSize<shapeBasedMeshingParameters.GetParameters().MaxElemSize())
    {
        cout<<"____Overriding autocomputed maximum element size____"<<endl;
        myEMeshMeshParameters.SetMinElemSize(maxSize);
    }
    //! --------------------
    //! load the parameters
    //! --------------------
    myTessellator.LoadParameters(myEMeshMeshParameters);

    //! -------------------------------------------------
    //! perform: the final mesh is stored into the shape
    //! -------------------------------------------------
    myTessellator.Perform();

    return myTessellator.IsDone();
}

//! ---------------------------------
//! function: SEFunc_GenerateSurface
//! details:
//! ---------------------------------
void occMesher::SEFunc_GenerateSurface()
{
    __try
    {
        cout<<"occMesher::SEFunc_GenerateSurface()->____in try____"<<endl;
        myTessellator.Perform();
    }
    __finally
    {
        cout<<"occMesher::SEFunc_GenerateSurface()->____in finally____"<<endl;
    }
}

//! -------------------------
//! function: performSurface
//! details:
//! -------------------------
userMessage occMesher::performSurface(occHandle(MeshVS_DataSource)& theMeshVS_DataSource,
                                        NCollection_Handle<NCollection_List<Standard_Integer>> &badFacesIndices)
{
    cout<<"____tag01____"<<endl;
    //! --------------------
    //! load the parameters
    //! --------------------
    myTessellator.LoadParameters(myEMeshMeshParameters);
    cout<<"____tag02____"<<endl;

    try
    {
        cout<<"____tag03____"<<endl;

        _set_se_translator(trans_func_O);
        cout<<"____tag04____"<<endl;

        this->SEFunc_GenerateSurface();
    }
    catch(int)
    {
        userMessage mr(false,"Express mesh cannot generate the surface mesh");
        return mr;
    }
    catch(SE_Exception e)
    {
        userMessage mr(false,"EMESH has generated a system error");
        return mr;
    }

    //! --------------------------
    //! retrieve the data sources
    //! --------------------------
    badFacesIndices = myTessellator.BadFacesIndices();

    bool isComputeNormals = true;
    bool isGroup = false;
    try
    {
        bool isDone = OMFTools::TriangulationToMesh(myTopoDS_Shape,myOMFDS_Mesh,isComputeNormals,isGroup);
        if(isDone)
        {
            theMeshVS_DataSource = new OMFVS_DataSource(myOMFDS_Mesh, Standard_True);
            QString message("EMESH the surface mesh has been generated");
            userMessage mr(true,message);
            return mr;
        }
        else
        {
            QString message("EMESH not all the face mesh data sources have been generated");
            userMessage mr(true,message);
            cout<<message.toStdString()<<endl;
            return mr;
        }
    }
    catch(...)
    {
        QString message("EMESH error in generating the face mesh datasources");
        userMessage mr(false,message);
        return mr;
    }
}

//! -------------------------
//! function: performSurface
//! details:
//! -------------------------
userMessage occMesher::performSurface(occHandle(MeshVS_DataSource)&  mainMesh2D,
                                        bool isSecondOrder,
                                        bool isElementStraight,
                                        bool deleteNgPointers,
                                        int done)
{
    Q_UNUSED(isElementStraight)
    Q_UNUSED(deleteNgPointers)

    //! ---------------------------------
    //! load the Express Mesh parameters
    //! ---------------------------------
    myTessellator.LoadParameters(myEMeshMeshParameters);

    cout<<"____tag08: min size: "<<myEMeshMeshParameters.MinElemSize()<<"____"<<endl;
    cout<<"____tag08: max size: "<<myEMeshMeshParameters.MaxElemSize()<<"____"<<endl;
    try
    {
        _set_se_translator(trans_func_O);
        this->SEFunc_GenerateSurface();
    }
    catch(int)
    {
        userMessage mr(false,"Express mesh cannot generate the surface mesh");
        return mr;
    }
    catch(SE_Exception e)
    {
        userMessage mr(false,"Express mesh has generated a system error");
        return mr;
    }

    Standard_Boolean isComputeNormals = Standard_True;
    Standard_Boolean isGroup = Standard_False;
    Standard_Boolean isDone;

    try
    {
        cout<<"____tag09____"<<endl;

        isDone = OMFTools::TriangulationToMesh(myTopoDS_Shape,myOMFDS_Mesh,isComputeNormals,isGroup);

        cout<<"____tag10____"<<endl;
    }
    catch(...)
    {
        userMessage mr(false,"Express mesh error when retrieving the mesh data source from the shape");
        return mr;
    }

    if(isDone)
    {
        occHandle(MeshVS_DataSource) theMeshVS_DataSource = new OMFVS_DataSource(myOMFDS_Mesh, Standard_True);
        NCollection_Handle<NCollection_List<Standard_Integer>> badFacesIndices = myTessellator.BadFacesIndices();

        //! ---------------------------------------------------------------
        //! convert into Netgen format
        //! the resulting Netgen mesh HAS face descriptors, and Netgen can
        //! use it for generating the volume mesh
        //! ---------------------------------------------------------------
        Ng_Mesh *NgMesh = Ng_NewMesh();
        NgMesh = MeshTools::OCCDSToNegtenSurfaceMesh(myTopoDS_Shape, badFacesIndices);

        if(NgMesh==NULL)
        {
            userMessage mr(false,"Cannot convert the Express mesh into a Netgen pointer");
            return mr;
        }

        //! ---------------------------------
        //! set the main mesh 2D data source
        //! ---------------------------------
        mainMesh2D = new Ng_MeshVS_DataSource2D(NgMesh);

        if(mainMesh2D.IsNull())
        {
            Ng_DeleteMesh(NgMesh);
            Ng_Exit();
            userMessage mr(false,"Cannot generate the surface mesh data source from Netgen pointer");
            return mr;
        }

        //! --------------------------------------------------------------------
        //! check if second order: the mesh midside nodes (2nd order) are added
        //! using Netgen function. Superparametric elements are generater
        //! --------------------------------------------------------------------
        if(isSecondOrder)
        {
            ccout("occMesher()->____converting surface mesh into 2nd order____");
            Ng_Generate_SecondOrder(NgMesh);
            mainMesh2D = new Ng_MeshVS_DataSource2D(NgMesh);
        }

        //! -----------------------------------
        //! generate the face mesh datasources
        //! -----------------------------------
#ifdef GENERATE_FACE_MESH_DATASOURCES
        this->generateFaceMeshDataSources(NgMesh, done);
#endif

        //! -----------------------------------------------------
        //! delete the Netgen pointer generated inside MeshTools
        //! -----------------------------------------------------
        //if(deleteNgPointers)
        //{
        //    cout<<"occMesher::performSurface()->____deleting the Netgen pointer____"<<endl;
        //    if (myNgMesh!=NULL) Ng_DeleteMesh(myNgMesh);
        //}

        userMessage mr(true,"Surface mesh sucessfully generated");
        return mr;
    }
    else
    {
        userMessage mr(false,"EMESH error in generating the face mesh datasources");
        return mr;
    }
}

//! -------------------------------------------------------------------------
//! function: generateFaceMeshDataSources
//! details:  generate the face mesh data sources from a Netgen mesh pointer
//! -------------------------------------------------------------------------
void occMesher::generateFaceMeshDataSources(Ng_Mesh *aNgMesh, int done)
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

    for(int faceNr=1; faceNr<=NbFaces;faceNr++)
    {
        const occHandle(Ng_MeshVS_DataSourceFace) &curFaceMeshDS = new Ng_MeshVS_DataSourceFace(aNgMesh, faceNr);
        faceDataSources.insert(faceNr,curFaceMeshDS);

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

//! --------------------------
//! function: deleteNgPointer
//! details:
//! --------------------------
//void occMesher::deleteNgPointers()
//{
//    if(myNgMesh!=NULL) Ng_DeleteMesh(myNgMesh);
//}

//! ------------------------
//! function: performVolume
//! details:
//! ------------------------
#include <tetgenmesher.h>
userMessage occMesher::performVolume(occHandle(MeshVS_DataSource)&  mainMesh3D,
                                       occHandle(MeshVS_DataSource)&  mainMesh2D,
                                       bool isSecondOrder,
                                       bool isElementStraight,
                                       int done)
{
    isSecondOrder = false;

    //! ----------------------------------------------------
    //! second order elements are generated by adding nodes
    //! ----------------------------------------------------
    isElementStraight = true;

    //! ----------------------------------------------------
    //! generate an Express mesh (first order) surface mesh
    //! ----------------------------------------------------
    userMessage mr;

    //! ---------------------------------
    //! load the Express Mesh parameters
    //! ---------------------------------
    myTessellator.LoadParameters(myEMeshMeshParameters);

    try
    {
        _set_se_translator(trans_func_O);
        this->SEFunc_GenerateSurface();
    }
    catch(int)
    {
        userMessage mr(false,"occMesher::performVolume: Express mesh cannot generate the surface mesh");
        return mr;
    }
    catch(SE_Exception e)
    {
        userMessage mr(false,"occMesher::performVolume: Express mesh has generated a system error");
        return mr;
    }

    Standard_Boolean isComputeNormals = Standard_True;
    Standard_Boolean isGroup = Standard_False;
    Standard_Boolean isDone;

    try
    {
        isDone = OMFTools::TriangulationToMesh(myTopoDS_Shape,myOMFDS_Mesh,isComputeNormals,isGroup);
    }
    catch(...)
    {
        userMessage mr(false,"occMesher::performVolume: Express mesh error when retrieving the mesh data source from the shape");
        return mr;
    }

    if(!isDone)
    {
        userMessage mr(false,"occMesher::performVolume: Express mesh error in generating the face mesh datasources");
        return mr;
    }

    occHandle(MeshVS_DataSource) theMeshVS_DataSource = new OMFVS_DataSource(myOMFDS_Mesh, Standard_True);
    NCollection_Handle<NCollection_List<Standard_Integer>> badFacesIndices = myTessellator.BadFacesIndices();

    //! ----------------------------------------------------------------
    //! convert into Netgen format
    //! the resulting Netgen mesh is first order, HAS face descriptors,
    //! and Netgen can use it for generating the volume mesh
    //! ----------------------------------------------------------------
    Ng_Mesh *NgMesh = MeshTools::OCCDSToNegtenSurfaceMesh(myTopoDS_Shape, badFacesIndices);
    if(NgMesh==NULL)
    {
        userMessage mr(false,"occMesher::performVolume: cannot convert the Express mesh into a Netgen pointer");
        return mr;
    }

    //! ---------------------------------
    //! set the main mesh 2D data source
    //! ---------------------------------
    mainMesh2D = new Ng_MeshVS_DataSource2D(NgMesh);

    if(mainMesh2D.IsNull())
    {
        Ng_DeleteMesh(NgMesh);
        userMessage mr(false,"occMesher::performVolume: cannot generate the surface mesh data source from Netgen pointer");
        return mr;
    }

    //! -----------------------------------
    //! generate the face mesh datasources
    //! -----------------------------------
#ifdef GENERATE_FACE_MESH_DATASOURCES
    this->generateFaceMeshDataSources(NgMesh, done);
#endif

    //! -------------------------------
    //! end of surface mesh generation
    //! -------------------------------
    switch(myVolumeMeshEngine)
    {
    case Property::meshEngine3D_Netgen:
    {
        //! --------------------------------------------------------
        //! the following constructor inits Netgen without geometry
        //! the constructor internally the following parameters
        //! mp.invert_trigs = 1
        //! mp.check_overlap = 1
        //! mp.optsteps_3d = 1
        //! --------------------------------------------------------
        NetgenMesher aNetgenMesher(myMeshParameters);
        mr = aNetgenMesher.performVolume(myMeshParameters,NgMesh,mainMesh3D);

        if(!mr.isDone) return mr;

        /*
        //! -------------------------------------------------------------
        //! get the face mesh data sources generated previously
        //! Note: the face mesh data sources are internally generated by
        //!       the method "performSurface"
        //! -------------------------------------------------------------
        QMap<int,occHandle(Ng_MeshVS_DataSourceFace)> faceMesDS = this->getFaceDataSources();
        for(QMap<int,occHandle(Ng_MeshVS_DataSourceFace)>::iterator it = faceMesDS.begin(); it!= faceMesDS.end(); ++it)
        {
            int faceNr = it.key();
            const occHandle(Ng_MeshVS_DataSourceFace) &curFaceDS = it.value();
            faceMesDS.insert(faceNr,curFaceDS);
        }
        */

        //! ----------------------
        //! check if second order
        //! ----------------------
        if(myMeshParameters.isSecondOrder) Ng_Generate_SecondOrder(NgMesh);

        //! -----------------------
        //! build the data sources
        //! -----------------------
        mainMesh3D = new Ng_MeshVS_DataSource3D(NgMesh);
        if(mainMesh3D.IsNull())
        {
            Ng_DeleteMesh(NgMesh);
            userMessage mr(false,"Cannot generate the surface mesh from Netgen pointer");
            return mr;
        }

        mainMesh2D = new Ng_MeshVS_DataSource2D(NgMesh);
        if(mainMesh2D.IsNull())
        {
            Ng_DeleteMesh(NgMesh);
            userMessage mr(false,"Cannot generate the surface volume from Netgen pointer");
            return mr;
        }

#ifdef GENERATE_FACE_MESH_DATASOURCES
        this->generateFaceMeshDataSources(NgMesh, done);
        /*
        //! ------------------------------------
        //! generate the face mesh data sources
        //! ------------------------------------
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
            occHandle(Ng_MeshVS_DataSourceFace) curFaceMeshDS = new Ng_MeshVS_DataSourceFace(NgMesh, faceNr);
            faceDataSources.insert(faceNr,curFaceMeshDS);

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
        */
#endif

    }
        break;
    }
    return mr;
}
