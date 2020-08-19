#ifndef RESULTPRESENTATION_H
#define RESULTPRESENTATION_H

struct resultPresentation
{
    //! ----------------------------
    //! type of result presentation
    //! ----------------------------
    enum typeOfPresentation
    {
        typeOfPresentation_isostrips,
        typeOfPresentation_isosurfaces,
        typeOfPresentation_isolines,
        typeOfPresentation_nodalresults
    };

    //! --------------
    //! combined view
    //! --------------
    enum combinedView
    {
        combinedView_resultOnly,
        combinedView_meshVisible,
        combinedView_undeformedWireFrame,
        combinedView_undeformedModel
    };

    combinedView theCombinedView;
    typeOfPresentation theTypeOfPresentation;
    bool useExteriorMeshForVolumeResults;
    double theScale;

    //! ------------
    //! constructor
    //! ------------
    resultPresentation(combinedView aCombinedView = combinedView_resultOnly,
                       typeOfPresentation aTypeOfPresentation = typeOfPresentation_isostrips,
                       double aScale = 1.0,
                       bool aUseExteriorMeshForVolumeResults = true):
        theCombinedView(aCombinedView),
        useExteriorMeshForVolumeResults(aUseExteriorMeshForVolumeResults),
        theScale(aScale),
        theTypeOfPresentation(aTypeOfPresentation)
    {;}

    //! -----------------
    //! copy constructor
    //! -----------------
    resultPresentation(const resultPresentation &aRP)
    {
        theTypeOfPresentation = aRP.theTypeOfPresentation;
        theCombinedView = aRP.theCombinedView;
        useExteriorMeshForVolumeResults = aRP.useExteriorMeshForVolumeResults;
        theScale = aRP.theScale;
    }

    //! -----------
    //! operator =
    //! -----------
    resultPresentation operator = (const resultPresentation &aRP)
    {
        theTypeOfPresentation = aRP.theTypeOfPresentation;
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
        if(theTypeOfPresentation != aRP.theTypeOfPresentation) return false;
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
