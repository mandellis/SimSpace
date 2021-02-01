#ifndef GRAPHICSTOOLS_H
#define GRAPHICSTOOLS_H

//! ----------------
//! custom includes
//! ----------------
#include "simulationnodeclass.h"
#include "ext/occ_extended/ais_colorscaleextended.h"

//! ----
//! OCC
//! ----
#include <Quantity_NameOfColor.hxx>

class graphicsTools
{
public:

    graphicsTools();

    static Quantity_NameOfColor getModelFeatureColor(SimulationNodeClass::nodeType theType);
    static void createColorBox(double min, double max, int Nintervals, occHandle(AIS_ColorScaleExtended) &aColorBox);
};

#endif // GRAPHICSTOOLS_H
