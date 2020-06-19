#ifndef AIS_SPHEREMARKER_H
#define AIS_SPHEREMARKER_H

#include <Standard_DefineHandle.hxx>
#include <AIS_Shape.hxx>
#include <AIS_KindOfInteractive.hxx>
#include <TopoDS_Shape.hxx>

class AIS_SphereMarker: public AIS_Shape
{
public:

    Standard_EXPORT AIS_SphereMarker(const TopoDS_Shape &shape);

    //! Returns Object as the type of Interactive Object
    Standard_EXPORT virtual AIS_KindOfInteractive Type() const Standard_OVERRIDE;

    //! Returns index 10 for "AIS_SphereMarker"
    Standard_EXPORT virtual Standard_Integer Signature() const Standard_OVERRIDE;

    DEFINE_STANDARD_RTTIEXT(AIS_SphereMarker,AIS_Shape)
};

DEFINE_STANDARD_HANDLE(AIS_SphereMarker, AIS_Shape)

#endif // AIS_SPHEREMARKER_H
