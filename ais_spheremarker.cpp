#include "AIS_SphereMarker.h"
#include "ais_customsignatures.h"

IMPLEMENT_STANDARD_HANDLE(AIS_SphereMarker,AIS_Shape)
IMPLEMENT_STANDARD_RTTIEXT(AIS_SphereMarker,AIS_Shape)

AIS_SphereMarker::AIS_SphereMarker(const TopoDS_Shape &shape): AIS_Shape(shape)
{
    ;
}

AIS_KindOfInteractive AIS_SphereMarker::Type() const
{
    return AIS_KOI_Shape;
}

Standard_Integer AIS_SphereMarker::Signature() const
{
    return CUSTOM_SIGNATURE_SPHERE_MARKER;
}
