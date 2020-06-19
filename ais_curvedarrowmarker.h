#ifndef AIS_CURVEDARROWMARKER_H
#define AIS_CURVEDARROWMARKER_H

//! OCC
#include <Standard_DefineHandle.hxx>
#include <AIS_Shape.hxx>
#include <AIS_KindOfInteractive.hxx>
#include <TopoDS_Shape.hxx>

class AIS_CurvedArrowMarker: public AIS_Shape
{
public:

    Standard_EXPORT AIS_CurvedArrowMarker(const TopoDS_Shape &shape);

    //! Returns Object as the type of Interactive Object
    Standard_EXPORT virtual AIS_KindOfInteractive Type() const Standard_OVERRIDE;

    //! Returns index 14 for "AIS_CurvedArrow"
    Standard_EXPORT virtual Standard_Integer Signature() const Standard_OVERRIDE;

    DEFINE_STANDARD_RTTIEXT(AIS_CurvedArrowMarker,AIS_Shape)
};

DEFINE_STANDARD_HANDLE(AIS_CurvedArrowMarker,AIS_Shape)

#endif // AIS_CURVEDARROWMARKER_H
