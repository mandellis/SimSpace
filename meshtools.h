#ifndef MESHTOOLS_H
#define MESHTOOLS_H

//! ----
//! C++
//! ----
#include <map>

//! ----
//! OCC
//! ----
#include <MeshVS_Mesh.hxx>
#include <MeshVS_DataSource.hxx>
#include <Poly_Triangulation.hxx>
#include <TopoDS_Shape.hxx>
#include <NCollection_Handle.hxx>
#include <NCollection_List.hxx>
#include <NCollection_Vector.hxx>
#include <NCollection_Array1.hxx>
#include <gp_Pnt.hxx>

#include <TopoDS_Face.hxx>
#include <TopoDS_Edge.hxx>

//! ---
//! Qt
//! ---
#include <QList>
#include <qhash.h>

//! ----------------
//! custom includes
//! ----------------
#include <ng_meshvs_datasourceface.h>
#include <ng_meshvs_datasource2d.h>
#include <meshdatabase.h>
#include <mesh.h>

namespace nglib
{
    #include <nglib.h>
}
using namespace nglib;

class QProgressIndicator;

class MeshTools
{
public:

    MeshTools();
    static bool toPolyTriangulation(const occHandle(Ng_MeshVS_DataSourceFace) &faceMeshDS, occHandle(Poly_Triangulation) &aTriangulation);
    //static bool toPolyTriangulation(const occHandle(MeshVS_DataSource) &faceMeshDS, occHandle(Poly_Triangulation) &aTriangulation);

    //! ---------------------------------------------------------------------
    //! extract the STL visualization mesh from the shape (it uses the disk)
    //! ---------------------------------------------------------------------
    static bool toSTLMesh1(const TopoDS_Shape &shape,
                          const QString &surfaceMeshFilePath,
                          occHandle(Ng_MeshVS_DataSource2D) &surfaceMeshDS,
                          std::vector<occHandle(Ng_MeshVS_DataSourceFace)> &arrayOfFaceSTL_MeshVS_DataSource,
                          QProgressIndicator *aProgressIndicator = Q_NULLPTR,
                          int done=0);

    //! ------------------------------------
    //! store mesh to shape for shaded view
    //! ------------------------------------
    static bool storeMeshToShape(const TopoDS_Shape &shape, occHandle(Poly_Triangulation) &theTriagulation);

    //! -----------------------------------------------------------------------
    //! convert the Express mesh stored within the shape into a Netgen pointer
    //! -----------------------------------------------------------------------
    static Ng_Mesh* OCCDSToNegtenSurfaceMesh(const TopoDS_Shape &theShape,
                                             const NCollection_Handle<NCollection_List<Standard_Integer>> &badFacesIndices);

    //! ------------------------------------------------------------------------------------
    //! generate a tetgen PLC from the array of mesh datasources, and write it on the disk
    //! into two separate files .node, .poly
    //! ------------------------------------------------------------------------------------
    static bool buildPLC(const NCollection_Array1<occHandle(Ng_MeshVS_DataSourceFace)> &arrayOfFaceDS,
                         const QString &nodeFilePath,
                         const QString &polyFilePath,
                         QProgressIndicator *progressIndicator = Q_NULLPTR);

    /*
    static bool buildColoredMesh(const occHandle(MeshVS_DataSource) &theMeshVS_DataSource,
                                 const QMap<int,double> &res,
                                 occHandle(MeshVS_Mesh) &aColoredMesh,
                                 double min = -100,
                                 double max = 100,
                                 int numberOfLevels = 10,
                                 bool showEdges = false,
                                 bool autoscale = true);

                                 */
    static bool buildColoredMesh(const occHandle(MeshVS_DataSource) &theMeshVS_DataSource,
                                 const std::map<int,double> &res,
                                 occHandle(MeshVS_Mesh) &aColoredMesh,
                                 double min = -100,
                                 double max = 100,
                                 int numberOfLevels = 10,
                                 bool showEdges = false,
                                 bool autoscale = true);

    static bool buildIsoStrip(const occHandle(MeshVS_DataSource) &theMeshVS_DataSource,
                              const std::map<int,double> &res,
                              const std::map<int,gp_Vec> &displacementMap,
                              double scale,
                              double min,
                              double max,
                              int NbLevels,
                              occHandle(MeshVS_Mesh) &aColoredMesh,
                              bool showEdges = false);

    static bool buildIsoSurfaces(const occHandle(MeshVS_DataSource) &theMeshVS_DataSource,
                                 const std::map<int,double> &res,
                                 double min,
                                 double max,
                                 int NbLevels,
                                 occHandle(MeshVS_Mesh) &aColoredMesh,
                                 bool showEdges = false);

    static bool buildIsoStrip(const occHandle(MeshVS_DataSource) &theMeshVS_DataSource,
                              const std::map<int,double> &res,
                              double min,
                              double max,
                              int NbLevels,
                              occHandle(MeshVS_Mesh) &aColoredMesh,
                              bool showEdges = false);

    static bool arrayOfFaceDataSourcesToExtendedStlFile(const NCollection_Array1<occHandle(Ng_MeshVS_DataSourceFace)> &arrayOfFaceMeshDS,
                                                        const QString &extendedStlFileName);

    static void sortPointsOnEdge(const QList<mesh::meshPoint> &pointsOfTheEdge,
                                 const TopoDS_Edge &anEdge,
                                 QList<mesh::meshPoint> &sortedPointsOfTheEdge);

    // commented since unused
    //static bool faceMeshDataSourceToListOfNodesElements(const occHandle(MeshVS_DataSource) &aMeshDS,
    //                                                    QList<mesh::meshPoint> &meshPointList,
    //                                                    QList<meshElement2D> &meshElementList);

    //! ---------------------------------------------
    //! Netgen pointer from surface mesh data source
    //! ---------------------------------------------
    static Ng_Mesh* surfaceMeshDSToNetgenMesh(const occHandle(MeshVS_DataSource) &aSurfaceMeshDS);

    //! -------------------------------
    //! filter volume elements by type
    //! -------------------------------
    static void filterVolumeElementsByType(const occHandle(MeshVS_DataSource) &inputMesh,
                                           ElemType aType,
                                           occHandle(MeshVS_DataSource) &outputMesh);

    static void toListOf3DElements(const occHandle(Ng_MeshVS_DataSource3D) &inputMesh, std::vector<meshElementByCoords> &elements);
    static std::map<GeometryTag,std::vector<occHandle(MeshVS_Mesh)>> groupMeshes(const std::map<GeometryTag,occHandle(MeshVS_Mesh)> &mapOfMeshes);
    static occHandle(MeshVS_DataSource) mergeMesh(const occHandle(MeshVS_DataSource) &mesh1, const occHandle(MeshVS_DataSource) &mesh2);
    static void computeAngleDefectMap(const occHandle(Ng_MeshVS_DataSourceFace) &aMeshDS, std::map<int,double> &mapOfAngleDefect);

};

class point
{
public:

    point();
    point(double x, double y, double z) {mX =x; mY = y; mZ = z;}

    double X() const {return mX;}
    double Y() const {return mY;}
    double Z() const {return mZ;}

private:

    double mX, mY, mZ;
};

//! -----------------------------------------------------------------
//! from Qt documentation
//! A QHash's key type has additional requirements other than being
//! an assignable data type: it must provide operator==(), and there
//! must also be a qHash() function in the type's namespace that
//! returns a hash value for an argument of the key's type.
//! -----------------------------------------------------------------
inline bool operator==(const point &p1, const point &p2)
{
    if(fabs(p1.X()-p2.X())<=1e-8 && fabs(p1.Y()-p2.Y())<=1e-8 && fabs(p1.Z()-p2.Z())<=1e-8)
        return true;
    else return false;
}

inline unsigned int qHash(const point &key, uint seed)
{
    char keystringc[128];
    sprintf(keystringc,"%.12e_%.12e_%.12e",key.X(),key.Y(),key.Z());
    QString keystring;
    keystring.fromLatin1(keystringc);
    return qHash(keystring, seed);
}

#endif // MESHTOOLS_H
