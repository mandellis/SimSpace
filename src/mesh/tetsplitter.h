#ifndef TETSPLITTER_H
#define TETSPLITTER_H

#include <mesh.h>
#include <meshelementbycoords.h>
#include <meshelement2d.h>

class tetSplitter
{
public:

    tetSplitter(){;}

    static std::vector<meshElementByCoords> splitTet(const meshElementByCoords &aTet);
    static std::vector<meshElementByCoords> splitTri(const meshElementByCoords &aTri);

private:

    meshElementByCoords myTet;

};

#endif // TETSPLITTER_H
