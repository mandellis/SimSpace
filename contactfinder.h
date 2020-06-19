#ifndef CONTACTFINDER_H
#define CONTACTFINDER_H

//! ---------------------------------------------------
//! note: the version compiled using polygon functions
//! is sligthly faster
//! ---------------------------------------------------
//!#define USE_POLYGON

//! ----
//! OCC
//! ----
#include <TopoDS_Shape.hxx>
#include <TopoDS_Solid.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS.hxx>

//! ----------------
//! custom includes
//! ----------------
#include <geometrytag.h>
#include <meshdatabase.h>
#include <polygon.h>
#include <mesh.h>
#include <geometrytag.h>

//! ----
//! C++
//! ----
#include <vector>
#include <map>

//! ---------
//! typedefs
//! ---------
typedef std::vector<polygon::Point> triangle;


class Ng_MeshVS_DataSource2D;
class Ng_MeshVS_DataSourceFace;
class MeshVS_DataSource;
class QProgressIndicator;

class contactFinder
{
public:

    //! ------------
    //! constructor
    //! ------------
    contactFinder();

    //! --------------
    //! set data base
    //! --------------
    void setDataBase(meshDataBase *mDB);

    //! -------------
    //! setBodyPairs
    //! -------------
    void setBodyPairs(const std::vector<std::pair<GeometryTag,GeometryTag>> &vectorOfBodyPairs);

    //! ----------------------------------------------------------------
    //! proximity calculation: source - reference; target - inner cycle
    //! ----------------------------------------------------------------
    bool sourceNodesToTargetElement1(const occHandle(Ng_MeshVS_DataSource2D) &aSourceMesh,
                                     const occHandle(Ng_MeshVS_DataSource2D) &aTargetMesh,
                                     const std::map<int,occHandle(MeshVS_DataSource)> &mapOfFaceMeshDSForBodyTarget,
                                     const std::map<int,occHandle(MeshVS_DataSource)> &mapOfFaceMeshDSForBodySource,
                                     double tolerance,
                                     const std::map<size_t,int> &trianglesMapForBodyTarget,
                                     const std::map<size_t,int> &trianglesMapForBodySource,
                                     std::vector<std::pair<int,int>> &vecSourceTargetFaceNr);

    //! ----------------------------------
    //! prepareSTLTessellation - use disk
    //! ----------------------------------
    int prepareSTLTessellation_useDisk(const std::vector<std::pair<GeometryTag,GeometryTag>> &vectorOfBodyPairs,
                                       bool rebuildTessellation,
                                       std::map<int,occHandle(MeshVS_DataSource)> &mapOfSurfaceMeshDataSources,
                                       std::map<int,std::map<int,occHandle(MeshVS_DataSource)>> &mapOfMapOfFaceMeshDataSources);


    //! ----------------------
    //! set angular criterion
    //! ----------------------
    void setAngularCriterion(double angle);

    //! ----------------------------
    //! check the angular criterion
    //! ----------------------------
    bool AngularCriterion(double xA, double yA, double zA, double xB, double yB, double zB);

    //! -------------------
    //! area of a triangle
    //! -------------------
    double triangleArea(const polygon::Point &V1, const polygon::Point &V2, const polygon::Point &V3)
    {
        //! ----------
        //!	x1	y1	1
        //!	x2	y2	1
        //!	x3	y3	1
        //! ----------
        double x1 = V1.x; double y1 = V1.y;
        double x2 = V2.x; double y2 = V2.y;
        double x3 = V3.x; double y3 = V3.y;
        return (x1*(y2-y3)-y1*(x2-x3)+1.0*(x2*y3-y2*x3));
    }

    //! --------
    //! perform
    //! --------
    bool perform(const std::vector<std::pair<GeometryTag,GeometryTag>> &vectorOfBodyPairs,
                 std::vector<std::pair<QVector<GeometryTag>,QVector<GeometryTag>>> &allContactsPairs,
                 double tolerance,
                 int grouping);

    //! ---------
    //! group by
    //! ---------
    void groupBy(const std::vector<std::pair<QVector<GeometryTag>, QVector<GeometryTag>>> &allContactsPairs,
                  std::vector<std::pair<QVector<GeometryTag>, QVector<GeometryTag> > > &reordered,
                  int grouping);

    //! ------------------------------------
    //! is point projection inside triangle
    //! ------------------------------------
    bool isPointProjectionInsideTriangle(const std::vector<polygon::Point> &aTriangle,
                                         const polygon::Point &P,
                                         polygon::Point &PP,
                                         double &distance);


    //! --------------
    //! triangle hash
    //! --------------
    size_t triangleHash(const triangle &aTriangle);

    //! -----------------------
    //! set progress indicator
    //! -----------------------
    void setProgressIndicator(QProgressIndicator *aProgressIndicator);

private:

    meshDataBase *myMDB;
    std::vector<std::pair<GeometryTag,GeometryTag>> myVectorOfBodyPairs;
    occHandle(Ng_MeshVS_DataSource2D) myMasterSurfaceMesh;
    occHandle(Ng_MeshVS_DataSource2D) mySlaveSurfaceMesh;

    GeometryTag findFaceTag(const TopoDS_Face &aFace);

    double myAngularCriterion;
    QProgressIndicator *myProgressIndicator;

    //! experimental
    std::map<int,occHandle(MeshVS_DataSource)> myMapOfSTLMeshDataSources;


private:

    //! ---------------------
    //! performance analysis
    //! ---------------------
    int getTimePoint();

    //! experimental
    void TriangleToTriangleTest(const triangle &trig1, const triangle &trig2, bool &intersecting, double &distance);
};

#endif // CONTACTFINDER_H
