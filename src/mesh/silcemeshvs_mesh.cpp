#include <slicemeshvs_mesh.h>
#include "ais_customsignatures.h"
#include <AIS_KindOfInteractive.hxx>

IMPLEMENT_STANDARD_HANDLE(SliceMeshVS_Mesh,MeshVS_Mesh)
IMPLEMENT_STANDARD_RTTIEXT(SliceMeshVS_Mesh,MeshVS_Mesh)

SliceMeshVS_Mesh::SliceMeshVS_Mesh():MeshVS_Mesh()
{
    ;
}

AIS_KindOfInteractive SliceMeshVS_Mesh::Type() const
{
    return AIS_KOI_None;
}

Standard_Integer SliceMeshVS_Mesh::Signature() const
{
    return CUSTOM_SIGNATURE_SLICED_MESH;
}
