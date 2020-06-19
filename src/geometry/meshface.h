#ifndef MESHFACE_H
#define MESHFACE_H

//! ----------------
//! custom includes
//! ----------------
#include "geometryface.h"
#include <ng_meshvs_datasourceface.h>

class meshFace: public geometryFace
{
public:

    meshFace(const opencascade::handle<Ng_MeshVS_DataSourceFace> &theFaceMeshDS = occHandle(Ng_MeshVS_DataSourceFace)());
    virtual double area() override;
    virtual bool pointProjection(double *aPoint, double *theProjectedPoint, double eps=0.0) override;
    void setGeometry(const opencascade::handle<Ng_MeshVS_DataSourceFace> &theFaceMeshDS);

    void randomPointOn(double *point);

private:

    opencascade::handle<Ng_MeshVS_DataSourceFace> myFaceMeshDS;

};

#endif // MESHFACE_H
