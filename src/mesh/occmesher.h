#ifndef OCCMESHER_H
#define OCCMESHER_H

//! OCC
#include <TopoDS_Shape.hxx>
#include <QMShape_Tessellator.hxx>
#include <QMShape_Parameters.hxx>
#include <QMData_MeshParameters.hxx>
#include <OMFVS_DataSource.hxx>
#include <OMFQM_IMesh.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <ng_meshvs_datasourceface.h>

//! Qt
#include <QObject>
#include <QMap>
#include <QList>

//! custom includes
#include "mydefines.h"
#include "se_exception.h"
#include <userMessage.h>
#include "qprogressindicator.h"
#include "property.h"

//! C++
#include <stdio.h>
#include <windows.h>
#include <eh.h>

//! Netgen
namespace nglib
{
#include <nglib.h>
}
using namespace nglib;


class OMFDS_Mesh;
class MeshVS_DataSource;

class occMesher: public QObject
{

    Q_OBJECT

public:

    //! constructor I
    occMesher(const TopoDS_Shape &aShape, QObject* parent=0);

    //! constructor II
    occMesher(const TopoDS_Shape &theShape, meshParam meshParams, QObject *parent=0);

    //! destructor
    ~occMesher();

    //! set global meshing parameters
    void setBodyMeshingParameters(meshParam mp);

    //! set the mesh method
    void setMeshMethod(int isSurfaceQuad, int isVolumeHexa);

    //! perform - check if used ... to do
    Standard_Boolean perform(NCollection_Handle<NCollection_List<Standard_Integer>> &badFacesIndices);

    //! generates a mesh for visualization
    Standard_Boolean perform(double angularDeflection, double linearDeflection, double minSize, double maxSize);

    //! perform surface
    userMessage performSurface(occHandle(MeshVS_DataSource)& theMeshVS_DataSource,
                                    NCollection_Handle<NCollection_List<Standard_Integer>> &badFacesIndices);

    //! perform surface
    userMessage performSurface(occHandle(MeshVS_DataSource)& mainMeshVS_DataSource2D,
                                 bool isSecondOrder=false,
                                 bool isElementStraight = false,
                                 bool deleteNgPointers=true,
                                 int done = 0);

    //! set the volume mesher
    void setVolumeMesher(Property::meshEngine3D theVolumeMeshEngine) { myVolumeMeshEngine = theVolumeMeshEngine; }

    //! perform volume
    userMessage performVolume(occHandle(MeshVS_DataSource)& mainMeshVS_DataSource3D,
                                occHandle(MeshVS_DataSource)& mainMeshVS_DataSource2D,
                                bool isSecondOrder=false,
                                bool isElementStraight=false,
                                int done=0);

    //! delete netgen pointers
    //void deleteNgPointers();

    //! the mesh DS
    occHandle(OMFDS_Mesh) myOMFDS_Mesh;

private:

    //! geometry and topology content
    TopoDS_Shape myTopoDS_Shape;

    //! the OCC mesher class
    QMShape_Tessellator myTessellator;

    //! the Express mesh meshing parameters
    QMData_MeshParameters myEMeshMeshParameters;

    //! meshing parameters (interface: both for Express Mesh and Netgen)
    meshParam myMeshParameters;

    //! mesh engine 3D
    Property::meshEngine3D myVolumeMeshEngine;

    //! Netgen mesh pointer
    //Ng_Mesh *myNgMesh;

    //! the progress indicator
    QProgressIndicator* myProgressIndicator;

    //! the mesh with faces sub-meshes
    occHandle(OMFDS_Mesh) myMeshDS_withSubMeshes;

    //! SE functions
    void SEFunc_GenerateSurface();
    void SEFunc_GenerateVolume();

    //! bad faces indices
    NCollection_Handle<NCollection_List<Standard_Integer>> myBadFacesIndices;

    //! face mesh data sources
    QMap<int,occHandle(Ng_MeshVS_DataSourceFace)> faceDataSources;

    //! generate the face mesh data sources
    void generateFaceMeshDataSources(Ng_Mesh *aNgMesh, int done=0);

public:

    NCollection_Handle<NCollection_List<Standard_Integer>> getBadFaces() { return myBadFacesIndices; }

    //! getFaceDataSources
    QMap<int,occHandle(Ng_MeshVS_DataSourceFace)> getFaceDataSources() { return faceDataSources; }

    //! set the progress indicator
    void setProgressIndicator(QProgressIndicator *progressIndicator) { myProgressIndicator = progressIndicator; }

    //! get the Netgen mesh pointer
    //Ng_Mesh* getNetgenMeshPointer() { return myNgMesh; }
};

void trans_func_O(unsigned int, EXCEPTION_POINTERS*);

#endif // occMesher_H
