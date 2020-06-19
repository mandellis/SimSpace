#ifndef AIS_ARROWMARKER_H
#define AIS_ARROWMARKER_H

//! OCC
#include <Standard_DefineHandle.hxx>
#include <AIS_Shape.hxx>
#include <AIS_KindOfInteractive.hxx>
#include <TopoDS_Shape.hxx>

class AIS_ArrowMarker: public AIS_Shape
{
public:

    Standard_EXPORT AIS_ArrowMarker(const TopoDS_Shape &shape);

    //! Returns Object as the type of Interactive Object
    Standard_EXPORT virtual AIS_KindOfInteractive Type() const Standard_OVERRIDE;

    //! Returns index 11 for "Arrow marker"
    Standard_EXPORT virtual Standard_Integer Signature() const Standard_OVERRIDE;

    DEFINE_STANDARD_RTTIEXT(AIS_ArrowMarker,AIS_Shape)
};

DEFINE_STANDARD_HANDLE(AIS_ArrowMarker,AIS_Shape)

#endif // AIS_ARROWMARKER_H
