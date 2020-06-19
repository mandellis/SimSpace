#ifndef RESULTPRESENTATION_H
#define RESULTPRESENTATION_H

#include <string>
#include <iostream>
#include <fstream>
using namespace std;

struct resultPresentation
{
    enum combinedView
    {
        combinedView_resultOnly,
        combinedView_undeformedWireFrame,
        combinedView_undeformedModel
    };

    //! ----------
    //! view mode
    //! ----------
    combinedView theCombinedView;

    //! ------------------------------
    //! result mesh is shown deformed
    //! ------------------------------
    bool isDeformedView;
    double theScale;

    //! ---------------------
    //! mesh edge visibility
    //! ---------------------
    bool isMeshVisible;

    //! ---------------------------------------------------
    //! constructor
    //! true scale/no undeformed wireframe/no mesh visible
    //! ---------------------------------------------------
    resultPresentation(combinedView aCombinedView=combinedView_resultOnly,
                       bool anIsDeformedView=false,
                       double aScale=1.0,
                       bool aMeshVisible=false):
        theCombinedView(aCombinedView),
        isDeformedView(anIsDeformedView),
        theScale(aScale),
        isMeshVisible(aMeshVisible)
    {;}

    //! -----------------
    //! copy constructor
    //! -----------------
    resultPresentation(const resultPresentation &aRP)
    {
        theCombinedView = aRP.theCombinedView;
        isDeformedView = aRP.isDeformedView;
        theScale = aRP.theScale;
        isMeshVisible = aRP.isMeshVisible;
    }

    //! -----------
    //! operator =
    //! -----------
    resultPresentation operator = (const resultPresentation &aRP)
    {
        theCombinedView = aRP.theCombinedView;
        isDeformedView = aRP.isDeformedView;
        theScale = aRP.theScale;
        isMeshVisible = aRP.isMeshVisible;
        return *this;
    }
};

#endif // RESULTPRESENTATION_H
