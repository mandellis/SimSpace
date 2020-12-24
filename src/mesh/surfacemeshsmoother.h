#ifndef SURFACEMESHSMOOTHER_H
#define SURFACEMESHSMOOTHER_H


//! ----------------
//! custom includes
//! ----------------
#include "occhandle.h"
#include "meshdatabase.h"

//! ----
//! OCC
//! ----
#include <TopoDS_Shape.hxx>
#include "qprogressindicator.h"

//! ---
//! Qt
//! ---
#include <QObject>

//! ----
//! C++
//! ----
#include <map>

//! -------
//! libigl
//! -------
#include <libigl/include/igl/embree/EmbreeIntersector.h>

class MeshVS_DataSource;
class Ng_MeshVS_DataSource2D;

class surfaceMeshSmoother: public QObject
{
    Q_OBJECT

public:

    //! constructor
    surfaceMeshSmoother(QProgressIndicator *aProgressIndicator = Q_NULLPTR, QObject *parent=0);

    //! perform
    bool perform(meshDataBase *mDB, int bodyIndex, const double geometricTolerance = 5.0);

    //! set progress indicator
    void setProgressIndicator(QProgressIndicator *aProgressIndicator);

private:

    //! progress indicator
    QProgressIndicator *myProgressIndicator;

    //! build shape tessellation
    bool buildShapeTessellation(const TopoDS_Shape &aShape,
                                occHandle(Ng_MeshVS_DataSource2D) &surfaceMeshDS_STL,
                                std::map<int,occHandle(Ng_MeshVS_DataSourceFace)> &mapOfFaceMeshDS_STL);

    //! init a raycaster
    bool initRayCaster(const occHandle(Ng_MeshVS_DataSource2D) &aSurfaceMeshDS);

    //! perturn mesh points
    void perturbMesh(occHandle(Ng_MeshVS_DataSource2D) &surfaceMesh);

    //! helper
    bool getPointCoordinate(const occHandle(MeshVS_DataSource) &aMeshDS, int globalNodeID, double *P);

private:

    //! the ray caster
    igl::embree::EmbreeIntersector myEmbree;

signals:

    void requestEnableStop();
    void requestDisableStop();
};

#endif // SURFACEMESHSMOOTHER_H
