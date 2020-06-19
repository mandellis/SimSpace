#include "ais_curvedarrowmarker.h"
#include "ais_customsignatures.h"

IMPLEMENT_STANDARD_HANDLE(AIS_CurvedArrowMarker,AIS_Shape)
IMPLEMENT_STANDARD_RTTIEXT(AIS_CurvedArrowMarker,AIS_Shape)

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
AIS_CurvedArrowMarker::AIS_CurvedArrowMarker(const TopoDS_Shape &shape): AIS_Shape(shape)
{
    ;
}

//! ---------------
//! function: Type
//! details:
//! ---------------
AIS_KindOfInteractive AIS_CurvedArrowMarker::Type() const
{
    return AIS_KOI_Shape;
}

//! --------------------
//! function: Signature
//! details:
//! --------------------
Standard_Integer AIS_CurvedArrowMarker::Signature() const
{
    return CUSTOM_SIGNATURE_CURVED_ARROW;
}
