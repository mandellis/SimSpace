#ifndef AIS_CUSTOMTRIHEDRON_H
#define AIS_CUSTOMTRIHEDRON_H

//! OCC
#include <Standard_DefineHandle.hxx>
#include <AIS_ColoredShape.hxx>
#include <AIS_KindOfInteractive.hxx>
#include <TopoDS_Shape.hxx>

class AIS_CustomTrihedron: public AIS_ColoredShape
{
public:

    Standard_EXPORT AIS_CustomTrihedron(const TopoDS_Shape &shape);

    //! Returns Object as the type of Interactive Object
    Standard_EXPORT virtual AIS_KindOfInteractive Type() const Standard_OVERRIDE;

    //! Returns index 13 for "AIS_CustomTrihedron"
    Standard_EXPORT virtual Standard_Integer Signature() const Standard_OVERRIDE;

    DEFINE_STANDARD_RTTIEXT(AIS_CustomTrihedron,AIS_ColoredShape)
};

DEFINE_STANDARD_HANDLE(AIS_CustomTrihedron,AIS_ColoredShape)

#endif // AIS_CUSTOMTRIHEDRON_H
