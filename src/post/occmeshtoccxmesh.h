#ifndef OCCMESHTOCCXMESH_H
#define OCCMESHTOCCXMESH_H

#include <meshdatabase.h>
#include "mydefines.h"

class OCCMeshToCCXmesh
{

public:

    static std::map<int,int> perform(GeometryTag loc, meshDataBase *mDB);
};

#endif // OCCMESHTOCCXMESH_H
