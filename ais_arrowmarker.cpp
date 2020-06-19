//! custom includes
#include "ais_arrowmarker.h"
#include "ais_customsignatures.h"

//! OCC
#include <TopoDS_Shape.hxx>

IMPLEMENT_STANDARD_HANDLE(AIS_ArrowMarker,AIS_Shape)
IMPLEMENT_STANDARD_RTTIEXT(AIS_ArrowMarker,AIS_Shape)

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
AIS_ArrowMarker::AIS_ArrowMarker(const TopoDS_Shape &shape): AIS_Shape(shape)
{
    ;
}

//! ---------------
//! function: Type
//! details:
//! ---------------
AIS_KindOfInteractive AIS_ArrowMarker::Type() const
{
    return AIS_KOI_Shape;
}

//! --------------------
//! function: Signature
//! details:
//! --------------------
Standard_Integer AIS_ArrowMarker::Signature() const
{
    return CUSTOM_SIGNATURE_ARROW_MARKER;
}
