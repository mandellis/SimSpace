#ifndef AIS_CORMARKER_H
#define AIS_CORMARKER_H

#include <Standard_DefineHandle.hxx>
#include <AIS_Shape.hxx>
#include <AIS_KindOfInteractive.hxx>
#include <TopoDS_Shape.hxx>

class AIS_CORMarker: public AIS_Shape
{
public:

    Standard_EXPORT AIS_CORMarker(const TopoDS_Shape &shape);

    //! Returns Object as the type of Interactive Object
    Standard_EXPORT virtual AIS_KindOfInteractive Type() const Standard_OVERRIDE;

    //! Returns index 8 for "Center of rotation"
    Standard_EXPORT virtual Standard_Integer Signature() const Standard_OVERRIDE;

    DEFINE_STANDARD_RTTIEXT(AIS_CORMarker,AIS_Shape)
};

DEFINE_STANDARD_HANDLE(AIS_CORMarker, AIS_Shape)

#endif // AIS_CORMARKER_H
