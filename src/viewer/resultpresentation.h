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
        combinedView_meshVisible,
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

    //! ---------------------------------------------------
    //! constructor
    //! true scale/no undeformed wireframe/no mesh visible
    //! ---------------------------------------------------
    resultPresentation(combinedView aCombinedView=combinedView_resultOnly,
                       bool anIsDeformedView=false,
                       double aScale=1.0):
        theCombinedView(aCombinedView),
        isDeformedView(anIsDeformedView),
        theScale(aScale)
    {;}

    //! -----------------
    //! copy constructor
    //! -----------------
    resultPresentation(const resultPresentation &aRP)
    {
        theCombinedView = aRP.theCombinedView;
        isDeformedView = aRP.isDeformedView;
        theScale = aRP.theScale;
    }

    //! -----------
    //! operator =
    //! -----------
    resultPresentation operator = (const resultPresentation &aRP)
    {
        theCombinedView = aRP.theCombinedView;
        isDeformedView = aRP.isDeformedView;
        theScale = aRP.theScale;
        return *this;
    }
};

#endif // RESULTPRESENTATION_H
