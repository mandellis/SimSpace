#ifndef TETGENMESHER_H
#define TETGENMESHER_H

#define TOLERANCE 1.0E-10

//#define PLCNEW

//! ---
//! Qt
//! ---
#include <QObject>
#include <QString>
#include <QProcess>
#include <QMap>

//! ----------------
//! custom includes
//! ----------------
#include "mydefines.h"
#include <mesh.h>
#include <meshdatabase.h>
#include <ng_meshvs_datasource3d.h>
#include <ng_meshvs_datasource2d.h>
#include <ng_meshvs_datasourceface.h>
#include "qprogressindicator.h"

//! ----
//! OCC
//! ----
#include <TopoDS_Shape.hxx>
#include <TopoDS_Solid.hxx>
#include <TopoDS.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>

//! -------
//! Tetgen
//! -------
#include <tetgen.h>

//! ----
//! C++
//! ----
#include <stdio.h>
#include <windows.h>
#include <eh.h>

class TetgenMesher: public QObject
{
    Q_OBJECT

public:

    //! -----------------------------------
    //! constructor without initialization
    //! -----------------------------------
    TetgenMesher(QObject*parent=0);

    //! -------------------------------------------------------
    //! constructor with initialization - requires setDataBase
    //! -------------------------------------------------------
    TetgenMesher(meshDataBase *mDB, QObject*parent=0);

    //! ----------------------
    //! set the mesh database
    //! ----------------------
    void init(meshDataBase *mDB);

    //! --------------------------------
    //! build the PLC running in memory
    //! --------------------------------
#ifdef PLCNEW
    bool buildPLC(int bodyIndex, QList<int> &invalidFaceTags, bool saveTetgenFiles = false);
#endif

#ifndef PLCNEW
    bool buildPLC(int bodyIndex, QList<int> &invalidFaceTags, bool saveTetgenFiles = false);
#endif

    bool buildPLC(const occHandle(MeshVS_DataSource) &aSurfaceMesh);

    //! --------------------
    //! set Tetgen switches
    //! --------------------
    void setSwitches(int preserveSurfaceMesh, int bodyIndex);

    //! --------
    //! perform
    //! --------
    bool perform(tetgenio *meshOut, int bodyIndex, bool saveMesh=true);

    //! -----------------------
    //! set progress indicator
    //! -----------------------
    void setProgressIndicator (QProgressIndicator *aProgressIndicator);

    //! --------------------
    //! perform on disk new
    //! --------------------
    int performOnDisk1(const NCollection_Array1<occHandle(Ng_MeshVS_DataSourceFace>) &inputArrayOfFaceDS,
                       int bodyIndex,
                       occHandle(Ng_MeshVS_DataSource3D) &tetgenMesh3D,
                       occHandle(Ng_MeshVS_DataSource2D) &tetgenMesh2D,
                       NCollection_Array1<occHandle(Ng_MeshVS_DataSourceFace)> &tetgenArrayOfFaceDS,
                       int done =0);

    //! --------------------
    //! perform on disk new
    //! --------------------
    int performOnDisk1(const NCollection_Array1<occHandle(Ng_MeshVS_DataSourceFace)> &inputArrayOfFaceDS,
                       int bodyIndex,
                       occHandle(Ng_MeshVS_DataSource3D) &tetgenMesh3D);


    //! ----------------------------------------------------------
    //! create a directory for support files (meshing using disk)
    //! ----------------------------------------------------------
    void createTetgenSupportFilesDir(bool clearPreviousContent=true);

    //! --------------------------------------------------
    //! delete the content of the support files directory
    //! --------------------------------------------------
    void deleteTetgenSupportFilesDir() {;}

    //! -----------------------------------------
    //! retrieve the mesh data sources from disk
    //! -----------------------------------------
    bool retrieveMeshDataSourcesFromDisk(int bodyIndex,
                                         occHandle(Ng_MeshVS_DataSource3D) &volumeMeshDS,
                                         occHandle(Ng_MeshVS_DataSource2D) &surfaceMeshDS,
                                         NCollection_Array1<occHandle(Ng_MeshVS_DataSourceFace)> &arrayOfFaceMeshDS,
                                         int done = 0);

    //! ----------------------------------------------------
    //! retrieve from disk only the volume mesh data source
    //! ----------------------------------------------------
    bool retrieveVolumeMeshDataSourceOnly(int bodyIndex, occHandle(Ng_MeshVS_DataSource3D) &volumeMeshDS);


private:

    //! --------------------------------
    //! create tetgen support directory
    //! --------------------------------
    void createTetgenDirs();

    //! --------------------------
    //! system exception handling
    //! --------------------------
    void SEFunc(tetgenio *meshOut);

private:

    meshDataBase *myMeshDB;
    tetgenio myPLC;
    QProcess *tetgenProcess;

    char mySwitches[128];

    QString myPLCDir;
    QString myPLCMesh;
    QString myTetgenSupportFilesDir;

    //! the progress indicator
    QProgressIndicator* myProgressIndicator;

    //! get the tetgen support files dir location
    QString getTetgenSupportFilesDir() const { return myTetgenSupportFilesDir; }

    //! enable/disable stop button
    void setStopButtonEnabled(bool isEnabled);

private slots:

    void redirectTetgenOutput();

public:

    //! read the .face file
    QList<QList<int>> readFaceFile(const QString &faceFileName);

    //! read the .node file
    QList<QList<double>> readNodeFile(const QString &nodeFileName);

    //! read the .node file
    bool readNodeFile(const QString &nodeFileName, QMap<int,mesh::meshPoint> &indexedMapOfAllMeshPoints);

    //! get the shell process
    QProcess* getProcess() { return tetgenProcess; }

signals:

    //! ---------------------------------------------------------
    //! signals for driving the button of the progress indicator
    //! ---------------------------------------------------------
    void requestEnableStop();
    void requestDisableStop();
};

void trans_func_T(unsigned int, EXCEPTION_POINTERS*);

#endif // TETGENMESHER_H
