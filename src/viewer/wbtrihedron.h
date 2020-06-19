#ifndef WBTRIHEDRON_H
#define WBTRIHEDRON_H

#include <Standard.hxx>
#include <Standard_DefineHandle.hxx>
#include <AIS_InteractiveObject.hxx>
#include "occhandle.h"
#include <Prs3d_Presentation.hxx>
#include <PrsMgr_PresentationManager3d.hxx>
#include <SelectMgr_Selection.hxx>

#include <Poly_Triangulation.hxx>
#include <TopLoc_Location.hxx>

class WBTrihedron: public AIS_InteractiveObject
{

public:

    /*
    enum selectionMode
    {
        selectionMode_wholeSelect,
        selectionMode_singleAxisSelect
    };
    */
public:

    //! constructor
    Standard_EXPORT WBTrihedron(double size = 10);

    //! set selection mode
    //void setSelectionMode(selectionMode aSelectionMode);

    //DEFINE_STANDARD_RTTIEXT(WBTrihedron,AIS_InteractiveObject)

protected:


private:

    double mySize;
    occHandle(Poly_Triangulation) polyTriangulationAxisZ;
    occHandle(Poly_Triangulation) polyTriangulationTipAxisZ;
    occHandle(Poly_Triangulation) polyTriangulationTipHollowDiskZ;
    occHandle(Poly_Triangulation) polyTriangulationAxisY;
    occHandle(Poly_Triangulation) polyTriangulationTipAxisY;
    occHandle(Poly_Triangulation) polyTriangulationTipHollowDiskY;
    occHandle(Poly_Triangulation) polyTriangulationAxisX;
    occHandle(Poly_Triangulation) polyTriangulationTipAxisX;
    occHandle(Poly_Triangulation) polyTriangulationTipHollowDiskX;

    //TopLoc_Location zLoc;

private:

    void Compute (const occHandle(PrsMgr_PresentationManager3d)& /*thePresentationManager*/,
                  const occHandle(Prs3d_Presentation)& thePresentation,
                  const Standard_Integer /*theMode*/) Standard_OVERRIDE;


    void ComputeSelection (const occHandle(SelectMgr_Selection)& theSel,
                           const Standard_Integer aMode) Standard_OVERRIDE;

};

//DEFINE_STANDARD_HANDLE(WBTrihedron,AIS_InteractiveObject)

#endif // WBTRIHEDRON_H
