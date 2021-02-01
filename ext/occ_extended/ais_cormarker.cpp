#include "ais_cormarker.h"
#include "ais_customsignatures.h"

IMPLEMENT_STANDARD_HANDLE(AIS_CORMarker,AIS_Shape)
IMPLEMENT_STANDARD_RTTIEXT(AIS_CORMarker,AIS_Shape)

AIS_CORMarker::AIS_CORMarker(const TopoDS_Shape &shape): AIS_Shape(shape)
{
    ;
}

AIS_KindOfInteractive AIS_CORMarker::Type() const
{
    return AIS_KOI_Shape;
}

Standard_Integer AIS_CORMarker::Signature() const
{
    return CUSTOM_SIGNATURE_COR_MARKER;
}
