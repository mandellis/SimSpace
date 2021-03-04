#ifndef POLYHEDRON_H
#define POLYHEDRON_H

#include <polygon.h>
#include <vector>
using namespace std;

namespace polyhedron
{
// begin namespace

typedef std::vector<polygon::Point> face;

class cell
{
public:

    //! constructor
    cell(const std::vector<face> &faces = std::vector<face>());

    //! number of faces
    int getNbFaces() const { return int(myFaces.size()); }

    //! return a face
    face getFace(int faceNr) const { return myFaces[faceNr]; }

    //! return the centroid
    bool getCentroid(double &xc, double &yc, double &zc) const;

    //! return the volume
    double getVolume() const;

private:

    std::vector<face> myFaces;
    std::vector<face> myTriangulation;
};

// end namespace
}

#endif // POLYHEDRON_H
