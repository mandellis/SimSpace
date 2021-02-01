#include <netgentools.h>
#include <ng_meshvs_datasourceface.h>
#include <ng_meshvs_datasource2d.h>
#include <ng_meshvs_datasource3d.h>

#include "src/utils/ccout.h"
#include "se_exception.h"
#include <ngexception.hpp>
#include "qprogressevent.h"

#include <TColStd_MapIteratorOfPackedMapOfInteger.hxx>

//! -------------------------
//! function: trans_func_STL
//! details:
//! -------------------------
void trans_func_STL(unsigned int u, EXCEPTION_POINTERS* pExp)
{
    Q_UNUSED(pExp)
    Q_UNUSED(u)
    printf("In trans_func_STL\n");
    throw SE_Exception();
}

//! ---------------------
//! function: destructor
//! details:
//! ---------------------
NetgenTools::~NetgenTools()
{
    cout<<"NetgenTools::~NetgenTools()->___destructor called____"<<endl;
}

//! -------------------------------
//! function: setProgressIndicator
//! details:
//! -------------------------------
void NetgenTools::setProgressIndicator (QProgressIndicator *theProgressIndicator)
{
    myProgressIndicator = theProgressIndicator;
    disconnect(myProgressIndicator,SIGNAL(stopPressed()),this,SLOT(abort()));
    connect(myProgressIndicator,SIGNAL(stopPressed()),this,SLOT(abort()));
}

//! ---------------------------
//! function: buildSTLBoundary
//! details:
//! ---------------------------
bool NetgenTools::buildSTLBoundary(int bodyIndex, Ng_STL_Geometry *anSTLgeo, Ng_Mesh *aNgMesh)
{
    cout<<"NetgenTools::buildSTLBoundary()->____function called____"<<endl;

    //! ----------------------------------------------------
    //! tags of the faces which do not have a triangulation
    //! ----------------------------------------------------
    myInvalidFaceTags.clear();

    //! -------------------------
    //! number of geometry faces
    //! -------------------------
    int Nb = myMDB->MapOfBodyTopologyMap.value(bodyIndex).faceMap.Extent();

    //! --------------------------------------
    //! add stl triangles to the netgen input
    //! --------------------------------------
    cout<<"NetgenTools::buildSTLBoundary()->____start adding triangles____"<<endl;

    int offset = 0;
    std::vector<int> faceTags;
    for(int faceNr = 0; faceNr<=Nb; faceNr++)
    {
        const occHandle(MeshVS_DataSource) &meshDS = myMDB->ArrayOfMeshDSOnFaces.getValue(bodyIndex,faceNr);
        if(meshDS.IsNull())
        {
            cout<<"NetgenTools::buildSTLBoundary()->____face nr: "<<faceNr<<" has not a valid mesh. Jumping over it____"<<endl;
            myInvalidFaceTags<<faceNr;
            offset++;
            continue;
        }

        //! -------------------------------------------------------
        //! only valid faces (i.e. having an .stl mesh) are tagged
        //! -------------------------------------------------------
        faceTags.push_back(faceNr);

        const occHandle(Ng_MeshVS_DataSourceFace) &aFaceMeshDS = occHandle(Ng_MeshVS_DataSourceFace)::DownCast(meshDS);
        TColStd_PackedMapOfInteger eMap = aFaceMeshDS->GetAllElements();
        for(TColStd_MapIteratorOfPackedMapOfInteger eMapIt(eMap);eMapIt.More();eMapIt.Next())
        {
            int globalElementID = eMapIt.Key();
            int NbNodes;
            MeshVS_EntityType type;
            Standard_Real buf[9];
            TColStd_Array1OfReal coord(*buf,1,9);
            aFaceMeshDS->GetGeom(globalElementID,true,coord,NbNodes,type);
            double p1[3], p2[3], p3[3];
            p1[0] = coord(1); p1[1] = coord(2); p1[2] = coord(3);
            p2[0] = coord(4); p2[1] = coord(5); p2[2] = coord(6);
            p3[0] = coord(7); p3[1] = coord(8); p3[2] = coord(9);

            Ng_STL_AddTriangle(anSTLgeo,p1,p2,p3);
        }
    }

    //! ----------------------------------------------
    //! all the mesh segments to be used for defining
    //! the netgen .stl edges
    //! ----------------------------------------------
    //QList<Ng_MeshVS_DataSourceFace::meshSegment> allSegments;
    QList<mesh::meshSegment> allSegments;

    //! ----------------------------------
    //! retrieve the contour of the faces
    //! ----------------------------------
    for(int faceNr=0; faceNr<=Nb; faceNr++)
    {
        const occHandle(MeshVS_DataSource) &curMeshDS = myMDB->ArrayOfMeshDSOnFaces.getValue(bodyIndex,faceNr);
        if(curMeshDS.IsNull()) continue;

        const occHandle(Ng_MeshVS_DataSourceFace) &curFaceMeshDS = occHandle(Ng_MeshVS_DataSourceFace)::DownCast(curMeshDS);
        if(curFaceMeshDS->myBoundarySegments.isEmpty()) curFaceMeshDS->computeFreeMeshSegments();

        //QList<Ng_MeshVS_DataSourceFace::meshSegment> meshSegments = curFaceMeshDS->myBoundarySegments;
        QList<mesh::meshSegment> meshSegments = curFaceMeshDS->myBoundarySegments;
        int NbMeshSegments = meshSegments.length();

        //cout<<"NetgenTools::buildSTLBoundary()->____face nr: "<<faceNr<<" Number of segments: "<<NbMeshSegments<<"____"<<endl;

        for(int i=0; i<NbMeshSegments; i++)
        {
            //! -----------------------------------------------
            //! take a segment mesh segment of the 1D boundary
            //! -----------------------------------------------
            //Ng_MeshVS_DataSourceFace::meshSegment aMeshSegment =  meshSegments.at(i);
            mesh::meshSegment aMeshSegment =  meshSegments.at(i);
            if(!allSegments.contains(aMeshSegment)) allSegments<<aMeshSegment;
        }
    }

    Ng_Result res;
    try
    {
        res = Ng_STL_InitSTLGeometry(anSTLgeo);
    }
    catch(...)
    {
        //if(anSTLgeo!=NULL)
        {
            cout<<"NetgenTools::buildSTLBoundary->____delete STL geometry pointer____"<<endl;
            //Ng_STL_DeleteGeometry(anSTLgeo);
            cout<<"NetgenTools::buildSTLBoundary->____STL geometry pointer not deleted____"<<endl;
        }
        //if(aNgMesh!=NULL)
        {
            cout<<"NetgenTools::buildSTLBoundary->____delete Netgen mesh pointer____"<<endl;
            Ng_DeleteMesh(aNgMesh);
            cout<<"NetgenTools::buildSTLBoundary->____Netgen mesh pointer deleted____"<<endl;
        }
        return false;
    }

    if(res==NG_SURFACE_INPUT_ERROR)
    {
        if(anSTLgeo!=NULL)
        {
            cout<<"NetgenTools::buildSTLBoundary->____delete STL geometry pointer____"<<endl;
            Ng_STL_DeleteGeometry(anSTLgeo);
            cout<<"NetgenTools::buildSTLBoundary->____STL geometry pointer deleted____"<<endl;
        }
        if(aNgMesh!=NULL)
        {
            cout<<"NetgenTools::performSurface()->____delete Netgen mesh pointer____"<<endl;
            Ng_DeleteMesh(aNgMesh);
            cout<<"NetgenTools::performSurface()->____Netgen mesh pointer deleted____"<<endl;
        }
        return false;
    }
    /*
    //! -----------------------------------------------------
    //! init Netgen STL geometry. Ng_STL_InitSTLGeometry():
    //! 1.  calls STLGeometry constructor
    //! 2.  read all the triangles. At this step the points
    //!     within the netgen stl data structure are defined
    //! -----------------------------------------------------    
    Ng_Result res;
    try
    {
        _set_se_translator(trans_func_STL);
        res = SEFunc_InitSTLGeometry(anSTLgeo);
    }
    catch(SE_Exception e)
    {
        Q_UNUSED(e)

        cout<<"NetgenTools::buildSTLBoundary()->____system error when generating the STL geometry____"<<endl;

        if(aNgMesh!=NULL) Ng_DeleteMesh(aNgMesh);
        if(anSTLgeo!=NULL) Ng_STL_DeleteGeometry(anSTLgeo);
        return false;
    }
    catch(netgen::NgException)
    {
        cout<<"NetgenTools::buildSTLBoundary()->____catching a NgException____"<<endl;

        //if(aNgMesh!= NULL) Ng_DeleteMesh(aNgMesh);
        //if(anSTLgeo!= NULL) Ng_STL_DeleteGeometry(anSTLgeo);
        return false;
    }
    if(res==NG_SURFACE_INPUT_ERROR)
    {
        //! cesere
        return false;
    }
    */

    //! ---------------------
    //! add all the segments
    //! ---------------------
    std::vector<std::pair<std::vector<double>,std::vector<double>>> vectorOfAllSegments;
    occHandle(Ng_MeshVS_DataSource2D) surfaceMeshDS = occHandle(Ng_MeshVS_DataSource2D)::DownCast(myMDB->ArrayOfMeshDS2D.value(bodyIndex));
    for(int j=0; j<allSegments.length(); j++)
    {
        //Ng_MeshVS_DataSourceFace::meshSegment aMeshSeg = allSegments.at(j);
        mesh::meshSegment aMeshSeg = allSegments.at(j);

        //! -------------------
        //! nodes of a segment
        //! -------------------
        MeshVS_EntityType type;
        int NN;

        //! -----------------------------
        //! first point: get coordinates
        //! -----------------------------
        int globalNodeID = aMeshSeg.nodeIDs.at(0);
        double buf[3];
        TColStd_Array1OfReal coords1(*buf,1,3);
        surfaceMeshDS->GetGeom(globalNodeID,false,coords1,NN,type);

        //! -----------------------------------
        //! define the first point of the edge
        //! -----------------------------------
        std::vector<double> point1;
        point1.push_back(coords1(1));
        point1.push_back(coords1(2));
        point1.push_back(coords1(3));

        //! ------------------------------
        //! second point: get coordinates
        //! ------------------------------
        globalNodeID = aMeshSeg.nodeIDs.at(1);
        double buf1[3];
        TColStd_Array1OfReal coords2(*buf1,1,3);
        surfaceMeshDS->GetGeom(globalNodeID,false,coords2,NN,type);

        //! ------------------------------------
        //! define the second point of the edge
        //! ------------------------------------
        std::vector<double> point2;
        point2.push_back(coords2(1));
        point2.push_back(coords2(2));
        point2.push_back(coords2(3));

        double P1[3],P2[3];
        P1[0] = coords1(1); P1[1] = coords1(2); P1[2] = coords1(3);
        P2[0] = coords2(1); P2[1] = coords2(2); P2[2] = coords2(3);

        std::pair<std::vector<double>,std::vector<double>> seg;
        std::vector<double> fP{ P1[0], P2[0], P1[2] };
        std::vector<double> sP{ P1[1], P2[1], P1[0] };

        seg.first = fP;
        seg.second = sP;
        vectorOfAllSegments.push_back(seg);

        Ng_STL_AddEdge(anSTLgeo,P1,P2);
    }

    //res = Ng_STL_MakeEdges(anSTLgeo,aNgMesh,&mp);
    Ng_STL_AddEdgeProgrammatically(anSTLgeo,aNgMesh,vectorOfAllSegments);

    //void *tess = (void*)(&tessellation);    // for testing new nglib function
    //Ng_STL_buildSurfaceMesh(tess);          // for testing new nglib function

    //! -------
    //! sizing
    //! -------
    mp.uselocalh = 1;
    mp.maxh = myMDB->ArrayOfMaxBodyElementSize.value(bodyIndex);
    mp.minh = myMDB->ArrayOfMinBodyElementSize.value(bodyIndex);
    mp.grading = myMDB->ArrayOfGradingValue.value(bodyIndex);
    mp.optsurfmeshenable = 2;
    mp.optsteps_2d = myMDB->ArrayOfSmoothingSteps.value(bodyIndex);
    mp.optvolmeshenable = 2;
    mp.optsteps_3d = myMDB->ArrayOfSmoothingSteps.value(bodyIndex);

    cout<<"NetgenTools::buildSTLBoundary->____netgen stl geometry generated____"<<endl;
    return true;
}

//! ---------------------------------
//! function: SEFunc_InitSTLGeometry
//! details:
//! ---------------------------------
Ng_Result NetgenTools::SEFunc_InitSTLGeometry(Ng_STL_Geometry *geo)
{
    Ng_Result res;
    __try
    {
        res = Ng_STL_InitSTLGeometry(geo);
        if(res==NG_OK)
        {
            cout<<"NetgenTools::SEFunc_InitSTLGeometry()->___in try: NG_OK____"<<endl;
        }
        else
        {
            cout<<"NetgenTools::SEFunc_InitSTLGeometry()->___in try NG_SURFACE_INPUT_ERROR____"<<endl;
        }
    }
    __finally
    {
        //! -----------
        //! diagnostic
        //! -----------
        if(res==NG_OK) cout<<"NetgenTools::SEFunc_InitSTLGeometry()->___in finally: NG_OK____"<<endl;
        else cout<<"NetgenTools::SEFunc_InitSTLGeometry()->___in finally: NG_SURFACE_INPUT_ERROR____"<<endl;
    }
    return res;
}

//! -------------------------
//! function: performSurface
//! details:
//! -------------------------
#include <memory>
userMessage NetgenTools::performSurface(int bodyIndex,
                                          occHandle(MeshVS_DataSource)& mainMeshVS_DataSource2D,
                                          bool isSecondOrder,
                                          bool isElementStraight,
                                          int done,
                                          bool deletePointers)
{
    cout<<"NetgenTools::performSurface()->____function called____"<<endl;

    //! ------------
    //! Netgen mesh
    //! ------------
    Ng_Mesh *aNgMesh = Ng_NewMesh();

    //! --------------------
    //! Netgen STL geometry
    //! --------------------
    Ng_STL_Geometry *anSTLgeo = Ng_STL_NewGeometry();

    //! ------------------------------------------------------------
    //! build the stl boundary: this fill also myInvalidFaceTags
    //! i.e. the list of geometry faces which have not an .stl mesh
    //! ------------------------------------------------------------
    bool isDone = this->buildSTLBoundary(bodyIndex,anSTLgeo,aNgMesh);
    if(!isDone)
    {
        //! ------------------------------------------------------
        //! the .stl boundary cannot be generated for some reason
        //! ------------------------------------------------------
        cout<<"NetgenTools::performSurface()->____error in generating the netgel stl boundary____"<<endl;

        //Ng_STL_DeleteGeometry(anSTLgeo);

        userMessage mr(false,"Cannot start the Netgen STL: the STL boundary cannot be built");
        return mr;
    }

    Ng_Result ng_res;
    try
    {
        _set_se_translator(trans_func_STL);
        this->SEFunc_GenerateSurface(anSTLgeo,aNgMesh,&mp,ng_res);
    }
    catch(SE_Exception e)
    {
        //! ----------------------------------------------------------------------
        //! unhandled system exception generated by the surface meshing algorithm
        //! ----------------------------------------------------------------------
        cout<<"NetgenTools::performSurface()->____unhandled exception: delete pointers____"<<endl;

        if(anSTLgeo!=NULL)
        {
            cout<<"NetgenTools::performSurface()->____delete STL geometry pointer____"<<endl;
            //Ng_STL_DeleteGeometry(anSTLgeo);
            cout<<"NetgenTools::performSurface()->____STL geometry pointer not deleted____"<<endl;
        }
        if(aNgMesh!=NULL)
        {
            cout<<"NetgenTools::performSurface()->____delete Netgen mesh pointer____"<<endl;
            Ng_DeleteMesh(aNgMesh);
            cout<<"NetgenTools::performSurface()->____Netgen mesh pointer deleted____"<<endl;
        }

        userMessage mr(false,"Netgen: unhandled exception while generating the surface mesh");
        return mr;
    }
    catch(netgen::NgException)
    {
        //! --------------
        //! handled error
        //! --------------
        cout<<"NetgenTools::performSurface()->____Netgen exception: delete pointers____"<<endl;

        if(anSTLgeo!=NULL)
        {
            cout<<"NetgenTools::performSurface()->____delete STL geometry pointer____"<<endl;
            //Ng_STL_DeleteGeometry(anSTLgeo);
            cout<<"NetgenTools::performSurface()->____STL geometry pointer not deleted____"<<endl;
        }
        if(aNgMesh!=NULL)
        {
            cout<<"NetgenTools::performSurface()->____delete Netgen mesh pointer____"<<endl;
            Ng_DeleteMesh(aNgMesh);
            cout<<"NetgenTools::performSurface()->____Netgen mesh pointer deleted____"<<endl;
        }
        userMessage mr(false,"Netgen error while generating the surface mesh");
        return mr;
    }

    //! ------------------------------------------------------------------------
    //! generate second order - always straight for the moment
    //! (additional reconstruction capabilities should be implemented for this)
    //! ------------------------------------------------------------------------
    if(isSecondOrder==true)
    {
        if(isElementStraight) Ng_STL_Generate_SecondOrder(anSTLgeo,aNgMesh);
        else Ng_STL_Generate_SecondOrder(anSTLgeo,aNgMesh);
    }

    mainMeshVS_DataSource2D = new Ng_MeshVS_DataSource2D(aNgMesh);
    if(mainMeshVS_DataSource2D.IsNull())
    {
        userMessage mr(false,"Error in generating the surface mesh data source");
        return mr;
    }
    myMDB->ArrayOfMeshDS2D.insert(bodyIndex,mainMeshVS_DataSource2D);

    cout<<"NetgenTools::performSurface()->____deleting the pointers____"<<endl;
    if(deletePointers==true)
    {
        cout<<"NetgenTools::performSurface()->___deleting the netgen mesh____"<<endl;
        Ng_DeleteMesh(aNgMesh);
        cout<<"NetgenTools::performSurface()->____netgen mesh deleted____"<<endl;

        cout<<"NetgenTools::performSurface()->____deleting the netgen stl geometry____"<<endl;
        //Ng_STL_DeleteGeometry(anSTLgeo);
        cout<<"NetgenTools::performSurface()->____netgen stl geometry deleted____"<<endl;
    }
    userMessage mr(true,"Surface mesh sucessfully generated");
    return mr;
}

//! ----------------------------------
//! function: SEFunct_generateSurface
//! details:
//! ----------------------------------
void NetgenTools::SEFunc_GenerateSurface(Ng_STL_Geometry *geo, Ng_Mesh *mesh, Ng_Meshing_Parameters *mp, Ng_Result &res)
{
    __try
    {
        cout<<"NetgenTools::SEFunc_GenerateSurface()->____try____"<<endl;
        res = Ng_STL_GenerateSurfaceMesh(geo,mesh,mp);
    }
    __finally
    {
        cout<<"NetgenTools::SEFunc_GenerateSurface()->____finally____"<<endl;
    }
}

//! ---------------------------------
//! function: SEFunct_generateVolume
//! details:
//! ---------------------------------
void NetgenTools::SEFunc_GenerateVolume(Ng_Mesh *mesh, Ng_Meshing_Parameters *mp, Ng_Result &res)
{
    __try
    {
        cout<<"NetgenTools::SEFunc_GenerateVolume()->____try____"<<endl;
        res = Ng_GenerateVolumeMesh(mesh,mp);
    }
    __finally
    {
        cout<<"NetgenTools::SEFunc_GenerateVolume()->____finally____"<<endl;
    }
}

//! ------------------------
//! function: performVolume
//! details:
//! ------------------------
userMessage NetgenTools::performVolume(int bodyIndex,
                                         occHandle(MeshVS_DataSource)& mainMeshVS_DataSource3D,
                                         occHandle(MeshVS_DataSource)& mainMeshVS_DataSource2D,
                                         bool isSecondOrder,
                                         bool isElementStraight,
                                         int done)
{
    Q_UNUSED(isElementStraight)
    userMessage mr;

    //! --------------------------------
    //! first generate the surface mesh
    //! --------------------------------
    Ng_Mesh *aNgMesh = Ng_NewMesh();
    Ng_STL_Geometry *anSTLgeo = Ng_STL_NewGeometry();

    //! -----------------------
    //! build the stl boundary
    //! -----------------------
    bool isDone = this->buildSTLBoundary(bodyIndex,anSTLgeo,aNgMesh);
    if(!isDone)
    {
        userMessage mr(false,"Cannot start the Netgen STL process");
        return mr;
    }

    Ng_Result ng_res;
    try
    {
        _set_se_translator(trans_func_STL);
        this->SEFunc_GenerateSurface(anSTLgeo,aNgMesh,&mp,ng_res);
    }
    catch(SE_Exception e)
    {
        //! ----------------------------------------------------------------------
        //! unhandled system exception generated by the surface meshing algorithm
        //! ----------------------------------------------------------------------
        cout<<"NetgenTools::performSurface()->____unhandled exception: delete pointers____"<<endl;

        if(aNgMesh!=NULL)
        {
            cout<<"NetgenTools::performSurface()->____deleting the Netgen mesh pointer____"<<endl;
            Ng_DeleteMesh(aNgMesh);
            cout<<"NetgenTools::performSurface()->____Netgen mesh pointer deleted____"<<endl;
        }
        if(anSTLgeo!=NULL)
        {
            cout<<"NetgenTools::performSurface()->____deleting the Netgen geomeetry mesh pointer____"<<endl;
            //Ng_STL_DeleteGeometry(anSTLgeo);
            cout<<"NetgenTools::performSurface()->____Netgen geometry mesh pointer NOT deleted____"<<endl;
        }

        userMessage mr(false,"Netgen system error while generating the surface mesh");
        return mr;
    }
    catch(netgen::NgException)
    {
        //! --------------
        //! handled error
        //! --------------
        cout<<"NetgenTools::performSurface()->____handled exception: delete pointers____"<<endl;

        if(aNgMesh!=NULL)
        {
            cout<<"NetgenTools::performSurface()->____deleting the Netgen mesh pointer____"<<endl;
            Ng_DeleteMesh(aNgMesh);
            cout<<"NetgenTools::performSurface()->____Netgen mesh pointer deleted____"<<endl;
        }
        if(anSTLgeo!=NULL)
        {
            cout<<"NetgenTools::performSurface()->____deleting the Netgen geomeetry mesh pointer____"<<endl;
            //Ng_STL_DeleteGeometry(anSTLgeo);
            cout<<"NetgenTools::performSurface()->____Netgen geometry mesh pointer NOT deleted____"<<endl;
        }

        userMessage mr(false,"Netgen error while generating the surface mesh");
        return mr;
    }

    if(ng_res == NG_ERROR) exit(1);

    mainMeshVS_DataSource2D = new Ng_MeshVS_DataSource2D(aNgMesh);
    if(mainMeshVS_DataSource2D.IsNull())
    {
        if(aNgMesh!=NULL)
        {
            cout<<"NetgenTools::performSurface()->____deleting the Netgen mesh pointer____"<<endl;
            Ng_DeleteMesh(aNgMesh);
            cout<<"NetgenTools::performSurface()->____Netgen mesh pointer deleted____"<<endl;
        }
        if(anSTLgeo!=NULL)
        {
            cout<<"NetgenTools::performSurface()->____deleting the Netgen geomeetry mesh pointer____"<<endl;
            Ng_STL_DeleteGeometry(anSTLgeo);
            cout<<"NetgenTools::performSurface()->____Netgen geometry mesh pointer deleted____"<<endl;
        }
        userMessage mr(false,"NgTool cannot convert the Netgen pointer into a surface mesh datasource");
        return mr;
    }

    try
    {
        Ng_Result ng_res;
        _set_se_translator(trans_func_STL);
        this->SEFunc_GenerateVolume(aNgMesh,&mp,ng_res);
    }
    catch(SE_Exception e)
    {
        //! ----------------------------------------------------------------------
        //! unhandled system exception generated by the surface meshing algorithm
        //! ----------------------------------------------------------------------        
        cout<<"NetgenTools::performSurface()->____unhandled exception: delete pointers____"<<endl;

        if(aNgMesh!=NULL)
        {
            cout<<"NetgenTools::performSurface()->____deleting the Netgen mesh pointer____"<<endl;
            Ng_DeleteMesh(aNgMesh);
            cout<<"NetgenTools::performSurface()->____Netgen mesh pointer deleted____"<<endl;
        }
        if(anSTLgeo!=NULL)
        {
            cout<<"NetgenTools::performSurface()->____deleting the Netgen geomeetry mesh pointer____"<<endl;
            //Ng_STL_DeleteGeometry(anSTLgeo);
            cout<<"NetgenTools::performSurface()->____Netgen geometry mesh pointer NOT deleted____"<<endl;
        }

        userMessage mr(false,"Netgen system exception while generating the volume mesh");
        return mr;
    }
    catch(netgen::NgException)
    {
        //! --------------
        //! handled error
        //! --------------
        cout<<"NetgenTools::performSurface()->____unhandled exception: delete pointers____"<<endl;

        if(aNgMesh!=NULL)
        {
            cout<<"NetgenTools::performSurface()->____deleting the Netgen mesh pointer____"<<endl;
            Ng_DeleteMesh(aNgMesh);
            cout<<"NetgenTools::performSurface()->____Netgen mesh pointer deleted____"<<endl;
        }
        if(anSTLgeo!=NULL)
        {
            cout<<"NetgenTools::performSurface()->____deleting the Netgen geomeetry mesh pointer____"<<endl;
            //Ng_STL_DeleteGeometry(anSTLgeo);
            cout<<"NetgenTools::performSurface()->____Netgen geometry mesh pointer NOT deleted____"<<endl;
        }

        userMessage mr(false,"Netgen error while generating the volume mesh");
        return mr;
    }

    //! -------------------------------------------------------------------------
    //! generate second order - always straight for the moment
    //! (additional "reconstruction capabilities should be implemented for this)
    //! -------------------------------------------------------------------------
    if(isSecondOrder==true)
    {
        if(isElementStraight) Ng_STL_Generate_SecondOrder(anSTLgeo,aNgMesh);
        else Ng_STL_Generate_SecondOrder(anSTLgeo,aNgMesh);
    }

    //! ----------------------------
    //! resulting mesh data sources
    //! 2D mesh data source
    //! ----------------------------
    mainMeshVS_DataSource2D = new Ng_MeshVS_DataSource2D(aNgMesh);
    if(mainMeshVS_DataSource2D.IsNull())
    {
        cout<<"NetgenTools::performVolume()->____cannot generate the surface mesh data source. Delete Metgen pointers____"<<endl;

        if(aNgMesh!=NULL) Ng_DeleteMesh(aNgMesh);
        if(anSTLgeo!=NULL) Ng_STL_DeleteGeometry(anSTLgeo);

        userMessage mr(false,"NetgenTools cannot convert the Netgen pointer into a surface mesh data source");
        return mr;
    }
    myMDB->ArrayOfMeshDS2D.insert(bodyIndex,mainMeshVS_DataSource2D);

    //! ----------------------------
    //! resulting mesh data sources
    //! 3D mesh data source
    //! ----------------------------
    mainMeshVS_DataSource3D = new Ng_MeshVS_DataSource3D(aNgMesh);
    if(mainMeshVS_DataSource3D.IsNull())
    {
        cout<<"NetgenTools::performVolume()->____cannot generate the volume mesh data source. Delete Metgen pointers____"<<endl;

        if(aNgMesh!=NULL)
        {
            cout<<"NetgenTools::performSurface()->____deleting the Netgen mesh pointer____"<<endl;
            Ng_DeleteMesh(aNgMesh);
            cout<<"NetgenTools::performSurface()->____Netgen mesh pointer deleted____"<<endl;
        }
        if(anSTLgeo!=NULL)
        {
            cout<<"NetgenTools::performSurface()->____deleting the Netgen geomeetry mesh pointer____"<<endl;
            //Ng_STL_DeleteGeometry(anSTLgeo);
            cout<<"NetgenTools::performSurface()->____Netgen geometry mesh pointer NOT deleted____"<<endl;
        }

        userMessage mr(false,"NetgenTools cannot convert the Netgen pointer into a surface mesh data source");
        return mr;
    }
    myMDB->ArrayOfMeshDS.insert(bodyIndex,mainMeshVS_DataSource3D);

    //! ---------------------------
    //! delete the Netgen pointers
    //! ---------------------------
    if(aNgMesh!=NULL)
    {
        cout<<"NetgenTools::performSurface()->____deleting the Netgen mesh pointer____"<<endl;
        Ng_DeleteMesh(aNgMesh);
        cout<<"NetgenTools::performSurface()->____Netgen mesh pointer deleted____"<<endl;
    }
    if(anSTLgeo!=NULL)
    {
        cout<<"NetgenTools::performSurface()->____deleting the Netgen geomeetry mesh pointer____"<<endl;
        //Ng_STL_DeleteGeometry(anSTLgeo);
        cout<<"NetgenTools::performSurface()->____Netgen geometry mesh pointer NOT deleted____"<<endl;
    }

    mr.isDone=true;
    mr.message=QString("Surface mesh sucessfully generated");
    return mr;
}

//! --------------------------------
//! function: initDefaultParameters
//! details:
//! --------------------------------
void NetgenTools::initDefaultParameters()
{
    //! -------------------------
    //! nglib meshing parameters
    //! -------------------------
    mp.uselocalh = 1;                   // Enable/Disable usage of local mesh size modifiers
    mp.maxh = 100;                      // Maximum allowed mesh size
    mp.minh= 1;                         // Minimum allowed ->global<- mesh size allowed
    mp.fineness = 0.5;                  // Mesh density: 0...1 (0 => coarse; 1 => fine)
    mp.grading = 0.3;                   // Mesh grading: 0...1 (0 => uniform mesh; 1 => aggressive local grading)
    mp.elementsperedge=3;               // Number of elements to generate per edge
    mp.elementspercurve=3;              // Elements to generate per curvature radius
    mp.closeedgeenable=0;               // Enable/Disable mesh refinement at close edges
    mp.closeedgefact=0;                 // Factor to use for refinement at close edges (larger => finer)
    mp.minedgelenenable=0;              // Enable/Disable user defined minimum edge length for edge subdivision
    mp.minedgelen=1e-2;                 // Minimum edge length to use while subdividing the edges (default = 1e-4)
    mp.second_order=0;                  // Generate second-order surface and volume elements
    mp.quad_dominated=0;                // Creates a Quad-dominated mesh
    mp.optsurfmeshenable=1;             // Enable/Disable automatic surface mesh optimization
    mp.optvolmeshenable=1;              // Enable/Disable automatic volume mesh optimization
    mp.optsteps_2d=2;                   // Number of optimize steps to use for 2-D mesh optimization
    mp.optsteps_3d=1;                   // Number of optimize steps to use for 3-D mesh optimization
    mp.check_overlapping_boundary = 1;
    mp.check_overlap = 1;
    mp.invert_tets = 0;
    mp.invert_trigs = 0;
}

//! ----------------
//! function: abort
//! details:
//! ----------------
void NetgenTools::abort()
{
    if(myIsRunning==true)
    {
        cout<<"NetgenTools::abort()->____abort received: stopping Netgen STL mesh engine____"<<endl;
        //Ng_Terminate();
        myIsRunning = false;
    }
}
