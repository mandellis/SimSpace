#ifndef PRISMATICLAYERSPARAMETER_H
#define PRISMATICLAYERSPARAMETER_H

//! ----------------
//! custom includes
//! ----------------
#include "property.h"

//! ---
//! Qt
//! ---
#include <QMetaType>

enum prismaticLayer_sizing
{
    prismaticLayer_sizing_FirstLayerThickness,
    prismaticLayer_sizing_TotalThickness
};
Q_DECLARE_METATYPE(prismaticLayer_sizing)

struct prismaticLayerParameters
{
    prismaticLayer_sizing typeOfSizing;
    int NbLayers;
    double firstLayerThickness;
    double totalThickness;
    double expRatio;

    //! ---------------
    //! new parameters
    //! ---------------
    double curvatureSensitivityForShrink;
    int NbGuidingVectorSmoothingSteps;
    int NbLayerThicknessSmoothingSteps;
    double curvatureSensitivityForGuidingVectorsSmoothing;
    double curvatureSensitivityForThicknessSmoothing;
    //! ----------------------
    //! end of new parameters
    //! ----------------------

    bool lockBoundary;

    bool checkSelfIntersections;
    bool checkMutualIntersections;
    int generationAlgorithm;
    int boundaryMeshType;               // "0" Hybrid "1" Tetrahedral

    //! ----------------------
    //! operator = assignment
    //! ----------------------
    prismaticLayerParameters operator = (const prismaticLayerParameters &rhs)
    {
        typeOfSizing = rhs.typeOfSizing;
        NbLayers = rhs.NbLayers;
        firstLayerThickness = rhs.firstLayerThickness;
        totalThickness = rhs.totalThickness;
        expRatio = rhs.expRatio;

        //! ---------------
        //! new parameters
        //! ---------------
        curvatureSensitivityForShrink = rhs.curvatureSensitivityForShrink;
        NbGuidingVectorSmoothingSteps = rhs.NbGuidingVectorSmoothingSteps;
        NbLayerThicknessSmoothingSteps = rhs.NbLayerThicknessSmoothingSteps;
        curvatureSensitivityForGuidingVectorsSmoothing = rhs.curvatureSensitivityForGuidingVectorsSmoothing;
        curvatureSensitivityForThicknessSmoothing = rhs.curvatureSensitivityForThicknessSmoothing;
        //! ----------------------
        //! end of new parameters
        //! ----------------------

        checkSelfIntersections = rhs.checkSelfIntersections;
        checkMutualIntersections = rhs.checkMutualIntersections;
        lockBoundary = rhs.lockBoundary;
        generationAlgorithm = rhs.generationAlgorithm;
        boundaryMeshType = rhs.boundaryMeshType;

        return *this;
    }
};
Q_DECLARE_METATYPE(prismaticLayerParameters)

#endif // PRISMATICLAYERSPARAMETER_H
