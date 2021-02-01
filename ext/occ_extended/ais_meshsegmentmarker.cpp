#include "ais_meshsegmentmarker.h"
#include "ais_customsignatures.h"

IMPLEMENT_STANDARD_HANDLE(AIS_MeshSegmentMarker,AIS_Shape)
IMPLEMENT_STANDARD_RTTIEXT(AIS_MeshSegmentMarker,AIS_Shape)

AIS_MeshSegmentMarker::AIS_MeshSegmentMarker(const TopoDS_Shape &shape): AIS_Shape(shape)
{
    ;
}

AIS_KindOfInteractive AIS_MeshSegmentMarker::Type() const
{
    return AIS_KOI_Shape;
}

Standard_Integer AIS_MeshSegmentMarker::Signature() const
{
    return CUSTOM_SIGNATURE_MESH_SEGMENT_MARKER;
}
