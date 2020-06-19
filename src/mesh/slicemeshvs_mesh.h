#ifndef SILCEMESHVS_MESH_H
#define SILCEMESHVS_MESH_H

#include <Standard_DefineHandle.hxx>
#include <MeshVS_Mesh.hxx>
#include <AIS_KindOfInteractive.hxx>

class SliceMeshVS_Mesh: public MeshVS_Mesh
{
public:

    Standard_EXPORT SliceMeshVS_Mesh();

    //! Returns Object as the type of Interactive Object
    Standard_EXPORT virtual AIS_KindOfInteractive Type() const Standard_OVERRIDE;

    //! Returns index 15 for this object
    Standard_EXPORT virtual Standard_Integer Signature() const Standard_OVERRIDE;

    DEFINE_STANDARD_RTTIEXT(SliceMeshVS_Mesh,MeshVS_Mesh)
};

DEFINE_STANDARD_HANDLE(SliceMeshVS_Mesh, MeshVS_Mesh)

#endif // SILCEMESHVS_MESH_H
