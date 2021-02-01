#ifndef MARKERBUILDER_H
#define MARKERBUILDER_H

#include "simulationnodeclass.h"
#include <geometrydatabase.h>

class markerBuilder
{
public:

    static bool addMarker(SimulationNodeClass *node, geometryDataBase *gDB);
};

#endif // MARKERBUILDER_H
