#ifndef AIS_EXTENDEDSHAPE_H
#define AIS_EXTENDEDSHAPE_H

//! redefinition of opencascade Handle() macro
#include "occhandle.h"

#include <Standard_DefineHandle.hxx>
#include <AIS_ColoredShape.hxx>
#include <Quantity_Color.hxx>
#include <Prs3d_Presentation.hxx>
#include <PrsMgr_PresentationManager3d.hxx>
#include <Prs3d_Projector.hxx>

class AIS_ExtendedShape: public AIS_ColoredShape
{

private:

    //! The geometry and topology content
    TopoDS_Shape myTopoDSShape;

    //! The name of the shape
    Standard_CString myName;

    //! index of the shape
    int myIndex;

    //! is visible
    bool isObjectVisible;

public:

    DEFINE_STANDARD_RTTIEXT(AIS_ExtendedShape,AIS_ColoredShape)

    //! constructor
    AIS_ExtendedShape(const TopoDS_Shape &);

    virtual void Compute(const Handle_PrsMgr_PresentationManager3d &,
                         const Handle_Prs3d_Presentation& aPresentation,
                         const Standard_Integer aMode);

    virtual void Compute(const Handle_Prs3d_Projector &aProjector, const Handle_Prs3d_Presentation &aPresentation);

    //! set the visibility flag
    void setShapeVisibility(Standard_Boolean isVisible) { isObjectVisible = isVisible; }

    //! get the visibility flag
    Standard_Boolean isVisible() { return isObjectVisible; } const

    //! set index
    void setIndex(int index) { myIndex = index; }

    //! get index
    int index() {return myIndex; }const

    //! get name
    Standard_CString name() { return myName; } const

    //! set name
    void setName(Standard_CString name) {myName = name; }
};

DEFINE_STANDARD_HANDLE(AIS_ExtendedShape, AIS_ColoredShape)

#endif // AIS_ExtendedShape_H
