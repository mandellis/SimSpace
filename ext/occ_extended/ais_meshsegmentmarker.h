#ifndef AIS_MESHSEGMENTMARKER_H
#define AIS_MESHSEGMENTMARKER_H

#include <Standard_DefineHandle.hxx>
#include <AIS_Shape.hxx>
#include <AIS_KindOfInteractive.hxx>
#include <TopoDS_Shape.hxx>

class AIS_MeshSegmentMarker: public AIS_Shape
{
public:

    Standard_EXPORT AIS_MeshSegmentMarker(const TopoDS_Shape &shape);

    //! Returns Object as the type of Interactive Object
    Standard_EXPORT virtual AIS_KindOfInteractive Type() const Standard_OVERRIDE;

    //! Returns index 9 for "Mesh segment"
    Standard_EXPORT virtual Standard_Integer Signature() const Standard_OVERRIDE;

    DEFINE_STANDARD_RTTIEXT(AIS_MeshSegmentMarker,AIS_Shape)
};

DEFINE_STANDARD_HANDLE(AIS_MeshSegmentMarker, AIS_Shape)


#endif // AIS_MESHSEGMENTMARKER_H
