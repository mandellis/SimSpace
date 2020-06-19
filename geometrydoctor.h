#ifndef GEOMETRYDOCTOR_H
#define GEOMETRYDOCTOR_H

#include <TopoDS_Shape.hxx>
#include <TopTools_IndexedMapOfShape.hxx>

class geometryDoctor
{
public:
    geometryDoctor(TopoDS_Shape shape);
    void BuildFMap();
    void Analyze();
    TopoDS_Shape Repair(double tolerance,Standard_Boolean fixsmalledges,
                        Standard_Boolean fixspotstripfaces,
                        Standard_Boolean sewfaces,
                        Standard_Boolean makesolids);

    TopTools_IndexedMapOfShape fmap, emap, vmap, somap, shmap, wmap;

private:
    TopoDS_Shape myShape;
};

#endif // GEOMETRYDOCTOR_H
