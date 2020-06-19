//! custom includes
#include "ais_doublearrowmarker.h"
#include "ais_customsignatures.h"

//! OCC
#include <TopoDS_Shape.hxx>

IMPLEMENT_STANDARD_HANDLE(AIS_DoubleArrowMarker,AIS_Shape)
IMPLEMENT_STANDARD_RTTIEXT(AIS_DoubleArrowMarker,AIS_Shape)

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
AIS_DoubleArrowMarker::AIS_DoubleArrowMarker(const TopoDS_Shape &shape): AIS_Shape(shape)
{
    ;
}

//! ---------------
//! function: Type
//! details:
//! ---------------
AIS_KindOfInteractive AIS_DoubleArrowMarker::Type() const
{
    return AIS_KOI_Shape;
}

//! --------------------
//! function: Signature
//! details:
//! --------------------
Standard_Integer AIS_DoubleArrowMarker::Signature() const
{
    return CUSTOM_SIGNATURE_DOUBLE_ARROW_MARKER;
}
