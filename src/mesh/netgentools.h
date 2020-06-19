#ifndef NETGENTOOLS_H
#define NETGENTOOLS_H

//! ----------------
//! custom includes
//! ----------------
#include <meshdatabase.h>
#include <userMessage.h>
#include <ng_meshvs_datasourceface.h>
#include "qprogressindicator.h"

//! ----
//! C++
//! ----
#include <stdio.h>
#include <windows.h>
#include <eh.h>

//! ---
//! Qt
//! ---
#include <QObject>
#include <QMap>

namespace nglib
{
    #include <nglib.h>
}
using namespace nglib;

class NetgenTools: public QObject
{
    Q_OBJECT

private:

    meshDataBase *myMDB;
    Ng_Meshing_Parameters mp;
    bool myIsRunning;
    QProgressIndicator *myProgressIndicator;

public:

    //! --------------
    //! constructor I
    //! --------------
    NetgenTools(QObject *parent=0):QObject(parent)
    {
        this->initDefaultParameters();
        myProgressIndicator = Q_NULLPTR;
        myIsRunning = false;
    }

    //! ---------------
    //! constructor II
    //! ---------------
    NetgenTools(meshDataBase *mDB, QObject *parent=0):QObject(parent),
        myMDB(mDB)
    {
        this->initDefaultParameters();
        myProgressIndicator = Q_NULLPTR;
        myIsRunning = false;
    }

    //! destructor
    ~NetgenTools();

    //! -----------------------
    //! set progress indicator
    //! -----------------------
    void setProgressIndicator (QProgressIndicator *theProgressIndicator);

    //! set mesh database
    void setMeshDataBase(meshDataBase *mDB) { myMDB = mDB; }

    //! ----------------
    //! perform surface
    //! ----------------
    userMessage performSurface(int bodyIndex,
                                 occHandle(MeshVS_DataSource)& mainMeshVS_DataSource2D,
                                 bool isSecondOrder=false,
                                 bool isElementStraight = false,
                                 int done=0,
                                 bool deletePointers=true);

    //! ---------------
    //! perform volume
    //! ---------------
    userMessage performVolume(int bodyIndex,
                                occHandle(MeshVS_DataSource)& mainMeshVS_DataSource3D,
                                occHandle(MeshVS_DataSource)& mainMeshVS_DataSource2D,
                                bool isSecondOrder=false,
                                bool isElementStraight=true,
                                int done=0);

    //! ---------------
    //! set is running
    //! ---------------
    void setIsRunning(bool isRunning) { myIsRunning = isRunning; }

    //! ---------------
    //! get is running
    //! ---------------
    bool isRunning() { return myIsRunning; }

    //! get face mesh data sources
    //QMap<int,occHandle(Ng_MeshVS_DataSourceFace)> getFaceMeshDataSources() { return faceMeshDS; }

private:

    //! ----------------------------------
    //! build the STL boundary with edges
    //! ----------------------------------
    bool buildSTLBoundary(int bodyIndex, Ng_STL_Geometry *anSTLgeo, Ng_Mesh *aNgMesh);

    //! -------------
    //! SE functions
    //! -------------
    void SEFunc_GenerateSurface(Ng_STL_Geometry *geo, Ng_Mesh *mesh, Ng_Meshing_Parameters *mp, Ng_Result &res);
    void SEFunc_GenerateVolume(Ng_Mesh *mesh, Ng_Meshing_Parameters *mp, Ng_Result &res);
    Ng_Result SEFunc_InitSTLGeometry(Ng_STL_Geometry *geo);

    //! init the default Meshing parameters
    void initDefaultParameters();

    //! face mesh data sources
    QMap<int,occHandle(Ng_MeshVS_DataSourceFace)> faceMeshDS;

    //! generate face mesh data sources
    void generateFaceMeshDataSources(int done=0);

    //! tags of invalid faces (which have not an .stl mesh)
    QList<int> myInvalidFaceTags;

public slots:

    void abort();
};

#endif // NETGENTOOLS_H
