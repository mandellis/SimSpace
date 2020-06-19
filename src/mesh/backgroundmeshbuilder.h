#ifndef BACKGROUNDMESHBUILDER_H
#define BACKGROUNDMESHBUILDER_H

#include <TopoDS_Shape.hxx>
#include <vectortool.h>

class backgroundMeshBuilder
{
public:

    //! constructor
    backgroundMeshBuilder(const std::vector<TopoDS_Shape> &listOfShape);

    //! set the shape
    void setShape(const std::vector<TopoDS_Shape> &listOfShapes) { myListOfShape = listOfShapes; }

private:

    std::vector<TopoDS_Shape> myListOfShape;
};

#endif // BACKGROUNDMESHBUILDER_H
