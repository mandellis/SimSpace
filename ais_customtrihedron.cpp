//! custom includes
#include "ais_customtrihedron.h"
#include "ais_customsignatures.h"

//! OCC
#include <TopoDS_Shape.hxx>

IMPLEMENT_STANDARD_HANDLE(AIS_CustomTrihedron,AIS_ColoredShape)
IMPLEMENT_STANDARD_RTTIEXT(AIS_CustomTrihedron,AIS_ColoredShape)

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
AIS_CustomTrihedron::AIS_CustomTrihedron(const TopoDS_Shape &shape): AIS_ColoredShape(shape)
{
    ;
}

//! ---------------
//! function: Type
//! details:
//! ---------------
AIS_KindOfInteractive AIS_CustomTrihedron::Type() const
{
    return AIS_KOI_Shape;
}

//! --------------------
//! function: Signature
//! details:
//! --------------------
Standard_Integer AIS_CustomTrihedron::Signature() const
{
    return CUSTOM_SIGNATURE_CUSTOM_TRIHEDRON_MARKER;
}
