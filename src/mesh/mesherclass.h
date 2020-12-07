#ifndef MESHERCLASS_H
#define MESHERCLASS_H

//#define CREATE1D_DATASOURCES

//! ----------------
//! custom includes
//! ----------------
#include <meshdatabase.h>
#include <ng_meshvs_datasourceface.h>
#include "qprogressevent.h"
#include "qprogressindicator.h"
#include <userMessage.h>
#include <tetgenmesher.h>
#include <nettgentool2.h>
#include <tethex.h>

//! ---
//! Qt
//! ---
#include <QObject>
#include <QList>

//! ----
//! C++
//! ----
#include <vector>

//! -------
//! global
//! -------
#include "global.h"

class MesherClass: public QObject
{

    Q_OBJECT

public:

    //! -------------
    //! constructors
    //! -------------
    MesherClass(meshDataBase *mDB, const QList<int> &listOfBodies, bool isVolume, QObject *parent=0);
    MesherClass(QObject *parent=0);

    ~MesherClass();

    void init(meshDataBase *mDB, const QList<int> &listOfBodies, bool isVolume);

private:

    meshDataBase *myMeshDB;
    QList<int> myBodyList;
    bool myIsVolume;
    QProgressIndicator *myProgressIndicator;
    int mainProgress;

private:

    //! --------------------------------------------------------------------------------------
    //! The Tetgen volume mesher (the Tetgen function have been incapsulated into the class).
    //! The Tetgen mesher can work in memory or on disk
    //! --------------------------------------------------------------------------------------
    TetgenMesher *myTetMesher;

    //! ---------------------------------------------------------------------------------------
    //! The Netgen surface abd volume mesher. It uses a completely modified version of the
    //! netgen original nglib.dll. In summary two functions have been added to above library,
    //! taking as an input a BRep file the other an .stl tessellation and the mesh controls.
    //! Both return the meshing results, i.e. the surface, volume, and the face submeshes.
    //! Internally they use shared (smart) pointers, in order to handle memory deallocation
    //! in case of failure or external interruption. The Netgen mesher allows user interaction
    //! ---------------------------------------------------------------------------------------
    NettgenTool2 *myNetgenMesher;

    //! ---------------------------------------------------------------------------------
    //! The old Netgen mesher, which uses the (only slighyly) modified nglib funtions,
    //! communicating the one with the other using mesh and geometry netgen pointers.
    //! At the moment whithin the mesherClass only a NetgenMesher is used for generating
    //! a volume mesh starting from a post inflation mesh. Also this function should be
    //! eliminated in the future, since it does not use shared pointers
    //! ---------------------------------------------------------------------------------
    //NetgenMesher *theNetgenMesher;

public:

    //! set the database
    void setDataBase(meshDataBase *mDB){ myMeshDB = mDB;}

    //! set the bodies to mesh
    void setBodyList(const QList<int> bodyList){ myBodyList = bodyList;}

    //! generate mesh
    void generateMesh();

    //! post a progress event
    void postUpdateEvent(QProgressEvent *aProgressEvent);

    //! set progress indicator
    void setProgressIndicator(QProgressIndicator *aProgressIndicator);

private:

    //! ----------------------------------
    //! surface mesh generation functions
    //! ----------------------------------
    userMessage Netgen_generateSurfaceMesh(int bodyIndex, occHandle(MeshVS_DataSource) &mainMesh2D);
    userMessage Netgen_STL_generateSurfaceMesh(int bodyIndex, occHandle(MeshVS_DataSource) &mainMesh2D);
    userMessage ExMesh_generateSurfaceMesh(const meshParam &mp, int bodyIndex, occHandle(MeshVS_DataSource) &mainMesh2D);

    bool STL_generateSurfaceMesh(const meshParam &mp, int bodyIndex, occHandle(MeshVS_DataSource) &mainMesh2D);

    //! ------------------------------------
    //! new functions/methods - volume mesh
    //! ------------------------------------
    userMessage Netgen_generateVolumeMesh(int bodyIndex, occHandle(MeshVS_DataSource) &mainMesh2D, occHandle(MeshVS_DataSource) &mainMesh3D);
    userMessage Netgen_STL_generateVolumeMesh(int bodyIndex, occHandle(MeshVS_DataSource) &mainMesh2D, occHandle(MeshVS_DataSource) &mainMesh3D);
    userMessage Tetgen_generateVolumeMesh(int bodyIndex, int preserveSurfaceMesh, bool runInMemory, int done);
    userMessage ExMesh_generateVolumeMesh(const meshParam &mp, int bodyIndex, occHandle(MeshVS_DataSource) &mainMesh2D, occHandle(MeshVS_DataSource) &mainMesh3D, int done);

    //! ----------------------------
    //! prismatic layers - pre algo
    //! ----------------------------
    userMessage PrismaticLayers_generatePrismaticMesh(int bodyIndex,
                                                      occHandle(Ng_MeshVS_DataSourceFace) &theLastInflatedMesh,
                                                      QList<occHandle(Ng_MeshVS_DataSourceFace)> &listOfInflatedMeshes);

    //! -----------------------------
    //! prismatic layers - post algo
    //! -----------------------------
    userMessage PrismaticLayers_generatePrismaticMesh_post(int bodyIndex,
                                                           occHandle(Ng_MeshVS_DataSource3D) &prismaticMeshDS);

    //! --------------------------
    //! tetrahedral boundary mesh
    //! --------------------------
    //userMessage PrismaticLayers_generateTetBoundaryMesh(int bodyIndex);

    //! -------------------------
    //! remove the support files
    //! -------------------------
    void removeSupportFiles(bool remove=true);

    //! processMesh
    bool processMesh(int bodyIndex);

    //! ----------------
    //! processSTL file
    //! ----------------
    QString processSTL(int bodyIndex);

    //! --------------------------------------------
    //! generate an .stl file as input for TetWild
    //! --------------------------------------------
    QString processTetWildSTL(int bodyIndex);

public slots:

    //! --------------------------
    //! abort the meshing process
    //! --------------------------
    void abort();

private slots:

    //! --------------------------
    //! connected to NettgenTool2
    //! --------------------------
    void handleRequestStartingNetgenEnquireTimer();
    void handleRequestStoppingNetgenEnquireTimer();

signals:

    void meshingFinished();

    //! ---------------------------
    //! connected to meshingServer
    //! ---------------------------
    void requestStoppingNetgenEnquireTimer();
    void requestStartingNetgenEnquireTimer();

    void requestEnableStop();
    void requestDisableStop();
};

#endif // MESHERCLASS_H
