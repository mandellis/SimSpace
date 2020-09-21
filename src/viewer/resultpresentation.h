#ifndef RESULTPRESENTATION_H
#define RESULTPRESENTATION_H

#include <string>
#include <iostream>
#include <fstream>
using namespace std;

struct resultPresentation
{
    //! ----------
    //! view mode
    //! ----------
    enum combinedView
    {
        combinedView_resultOnly,
        combinedView_meshVisible,
        combinedView_undeformedWireFrame,
        combinedView_undeformedModel
    };

    combinedView theCombinedView;
    bool useExteriorMeshForVolumeResults;
    double theScale;

    //! ------------
    //! constructor
    //! ------------
    resultPresentation(combinedView aCombinedView = combinedView_resultOnly,
                       double aScale = 1.0,
                       bool aUseExteriorMeshForVolumeResults = true):
        theCombinedView(aCombinedView),
        useExteriorMeshForVolumeResults(aUseExteriorMeshForVolumeResults),
        theScale(aScale)
    {;}

    //! -----------------
    //! copy constructor
    //! -----------------
    resultPresentation(const resultPresentation &aRP)
    {
        theCombinedView = aRP.theCombinedView;
        useExteriorMeshForVolumeResults = aRP.useExteriorMeshForVolumeResults;
        theScale = aRP.theScale;
    }

    //! -----------
    //! operator =
    //! -----------
    resultPresentation operator = (const resultPresentation &aRP)
    {
        theCombinedView = aRP.theCombinedView;
        useExteriorMeshForVolumeResults = aRP.useExteriorMeshForVolumeResults;
        theScale = aRP.theScale;
        return *this;
    }

    //! ------------
    //! operator ==
    //! ------------
    bool operator == (const resultPresentation &aRP)
    {
        if(theCombinedView != aRP.theCombinedView) return false;
        if(useExteriorMeshForVolumeResults != aRP.useExteriorMeshForVolumeResults) return false;
        if(theScale != aRP.theScale) return false;
        return true;
    }

    //! ------------
    //! operator !=
    //! ------------
    bool operator != (const resultPresentation &aRP)
    {
        if(*this == aRP) return false;
        return true;
    }
};

#endif // RESULTPRESENTATION_H
