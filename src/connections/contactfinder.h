#ifndef CONTACTFINDER_H
#define CONTACTFINDER_H

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
#include "occhandle.h"
#include "geometrytag.h"
#include <meshdatabase.h>
#include "polygon.h"
#include "mesh.h"
#include "geometrytag.h"

//! --------
//! C++ STL
//! --------
#include <vector>
#include <map>
#include <string>

#include <memory>

//! -------
//! libigl
//! -------
#include <ext/libigl/include/igl/embree/EmbreeIntersector.h>

//! ---------
//! typedefs
//! ---------
typedef std::vector<polygon::Point> triangle;

class Ng_MeshVS_DataSource2D;
class Ng_MeshVS_DataSourceFace;
class MeshVS_DataSource;
class QProgressIndicator;
class geometryDataBase;

class contactFinder
{

private:

    geometryDataBase *myGDB;
    std::map<int,std::map<int,occHandle(Ng_MeshVS_DataSourceFace)>> m_subMeshes;
    std::map<int,std::map<int,int>> MM;
    std::map<int,occHandle(Ng_MeshVS_DataSourceFace)> mapSurfaceMeshDS;
    std::vector<std::pair<GeometryTag,GeometryTag>> myVectorOfBodyPairs;
    std::map<int,std::shared_ptr<igl::embree::EmbreeIntersector>> mapOfRayCasters;
    QProgressIndicator *myProgressIndicator;

public:

    contactFinder(const std::vector<std::pair<GeometryTag,GeometryTag>> &bodyPairs,
                  geometryDataBase *gDB,
                  QProgressIndicator *aProgressIndicator=nullptr);

    void setProgressIndicator(QProgressIndicator *aProgressIndicator);

    int prepareSTLTessellation_useDisk(const std::vector<std::pair<GeometryTag,GeometryTag>> &vectorOfBodyPairs,
                                       bool rebuildTessellation,
                                       std::map<int,occHandle(MeshVS_DataSource)> &mapOfSurfaceMeshDataSources,
                                       std::map<int,std::map<int,occHandle(Ng_MeshVS_DataSourceFace)>> &mapOfMapOfFaceMeshDataSources);

    //bool AngularCriterion(double xA, double yA, double zA, double xB, double yB, double zB);

    bool perform(std::vector<std::pair<std::vector<GeometryTag>,std::vector<GeometryTag>>> &allContactsPairs,
                 double tolerance,
                 int grouping);

    void groupBy(const std::vector<std::pair<std::vector<GeometryTag>, std::vector<GeometryTag>>> &allContactsPairs,
                  std::vector<std::pair<std::vector<GeometryTag>, std::vector<GeometryTag> > > &reordered,
                  int grouping);

private:

    void initRayCaster(const std::map<int, occHandle(Ng_MeshVS_DataSourceFace)> &mapSurfaceMeshes);
    int getTimePoint();
    std::string getPathOfExecutable();
};

#endif // CONTACTFINDER_H
