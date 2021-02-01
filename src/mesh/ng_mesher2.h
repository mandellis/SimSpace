#ifndef NG_MESHER2_H
#define NG_MESHER2_H

#define EXPERIMENTAL_SURFACE

//! ----------------
//! custom includes
//! ----------------
#include <meshdatabase.h>
#include "src/main/mydefines.h"
#include "se_exception.h"
#include <userMessage.h>
#include "qprogressindicator.h"

//! ----
//! OCC
//! ----
#include <BRep_Tool.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>
#include <TopExp.hxx>
#include <MeshVS_Mesh.hxx>
#include <MeshVS_DataSource.hxx>
#include <ng_meshvs_datasourceface.h>
#include <NCollection_Array1.hxx>

//! ---
//! Qt
//! ---
#include <QObject>
#include <QMap>

//! ----
//! C++
//! ----
#include <stdio.h>
#include <windows.h>
#include <eh.h>

//! Netgen
namespace nglib
{
    #include <nglib.h>
}
using namespace nglib;

class NetgenMesher: public QObject
{

    Q_OBJECT

private:

    //! NgMesh
    Ng_Mesh *myNgMesh;

    //! TopoDS_Shape
    TopoDS_Shape myTopoDS_Shape;

    //! Ng_OCC_Geometry
    Ng_OCC_Geometry *myNg_OCC_geometry;

    //! Netgen meshing parameters
    Ng_Meshing_Parameters mp;

    //! the progress indicator
    QProgressIndicator* myProgressIndicator;

public:

    //! constructor
    NetgenMesher(QObject *parent = 0);

    //! configure - for constructor without parameters
    void configure(const TopoDS_Shape &shape, meshParam theMeshingParams);

    //! contructor I
    NetgenMesher(meshParam meshParams, QObject *parent=0);

    //! constructor II
    NetgenMesher(TopoDS_Shape theShape, meshParam mp, QObject *parent=0);

    //! destructor
    ~NetgenMesher();

    //! -----------------------
    //! set progress indicator
    //! -----------------------
    void setProgressIndicator (QProgressIndicator *theProgressIndicator);

    //! perform surface
    userMessage performSurface(occHandle(MeshVS_DataSource)& mainMeshVS_DataSource2D,
                                 bool isSecondOrder=false,
                                 bool isElementStraight = false,
                                 bool deleteNgPointers=true);

    //! perform volume
    userMessage performVolume(occHandle(MeshVS_DataSource)& mainMeshVS_DataSource3D,
                                occHandle(MeshVS_DataSource)& mainMeshVS_DataSource2D,
                                bool isSecondOrder=false,
                                bool isElementStraight=false);

    //! perform volume
    userMessage performVolume(meshParam meshingParameters,
                                Ng_Mesh *mesh,
                                occHandle(MeshVS_DataSource) &mainMeshVS_DataSource3D);

    //! generate face mesh data sources
    void generateFaceMeshDataSources(int done=0);

    //! delete netgen pointers
    void deleteNgPointers();

private:

    //! init the default parameters
    void initDefaultParameters();

    //! system error functions
    void SEFunc_GenerateSurface(Ng_Meshing_Parameters *mp);
    void SEFunc_GenerateVolume(Ng_Mesh *NgMesh, Ng_Meshing_Parameters *mp);

    //! face mesh data sources
    QMap<int,occHandle(Ng_MeshVS_DataSourceFace)> faceMeshDS;

public:

    //! get face mesh data sources
    QMap<int,occHandle(Ng_MeshVS_DataSourceFace)> getFaceDataSources() { return faceMeshDS; }

public slots:

    //! abort
    void abort();
};

void trans_func_N(unsigned int, EXCEPTION_POINTERS*);

#endif // NG_MESHER2_H
