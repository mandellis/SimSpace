#ifndef SURFACEMESHTOFACEMESHES_H
#define SURFACEMESHTOFACEMESHES_H

//! ----------------
//! custom includes
//! ----------------
#include <meshelement2d.h>
#include <occhandle.h>

//! ------------
//! opencascade
//! ------------
#include <TopoDS_Shape.hxx>

//! ----
//! C++
//! ----
#include <map>
#include <iostream>
using namespace std;

//! -------
//! libigl
//! -------
#include <libigl/include/igl/embree/EmbreeIntersector.h>

//! --------
//! Windows
//! --------
#include <Windows.h>

class Ng_MeshVS_DataSourceFace;
class Ng_MeshVS_DataSource2D;

class surfaceMeshToFaceMeshes
{
public:

    surfaceMeshToFaceMeshes(const occHandle(Ng_MeshVS_DataSource2D) &aMeshDS = occHandle(Ng_MeshVS_DataSource2D)(), const TopoDS_Shape &aShape = TopoDS_Shape());
    void setShape(const TopoDS_Shape &aShape);
    void setMesh(const occHandle(Ng_MeshVS_DataSource2D) &aMeshDS);
    bool perform(std::map<int,occHandle(Ng_MeshVS_DataSourceFace)> &theFaceMeshDataSources);

private:

    bool buildShapeTessellation(const TopoDS_Shape &aShape,
                                occHandle(Ng_MeshVS_DataSource2D) &surfaceMeshDS_STL,
                                std::map<int, occHandle(Ng_MeshVS_DataSourceFace)> &mapOfFaceMeshDS_STL);      // fills the myMeshElementToFaceNumberMap

    void initRayCaster(const occHandle(Ng_MeshVS_DataSource2D) &aMeshDS);

    std::string getPathOfExecutable();

private:

    occHandle(Ng_MeshVS_DataSource2D) myMeshDS;
    Eigen::MatrixXd myV;
    Eigen::MatrixXi myF;

    igl::embree::EmbreeIntersector myEmbree;

    TopoDS_Shape myShape;

};

#endif // SURFACEMESHTOFACEMESHES_H
