#ifndef QMORPH_H
#define QMORPH_H

#include "occhandle.h"

//! ----
//! OCC
//! ----
#include <TopoDS_Face.hxx>

//! ----------------
//! custom includes
//! ----------------
#include <ng_meshvs_datasourceface.h>
#include <mesh.h>

//! ----
//! C++
//! ----
#include <map>

class QMorph
{
public:

    QMorph(const TopoDS_Face &aFace, const occHandle(Ng_MeshVS_DataSourceFace) &trigMesh);

private:

    TopoDS_Face myFace;
    occHandle(Ng_MeshVS_DataSourceFace) myTrigMesh;

    QList<mesh::meshSegment> myFront;

    enum meshEdgeStatus
    {
        state00,
        state01,
        state10,
        state11
    };

    void setMeshEdgeState();

};

#endif // QMORPH_H
