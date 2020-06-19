#include "backgroundmeshbuilder.h"
#include "geomtoolsclass.h"
#include <Bnd_Box.hxx>
#include <BRepBndLib.hxx>

backgroundMeshBuilder::backgroundMeshBuilder(const std::vector<TopoDS_Shape> &listOfShape):
    myListOfShape(listOfShape)
{
    /*
    //! --------------------------------
    //! retrieve the shape bounding box
    //! --------------------------------
    Bnd_Box boundingBox;
    BRepBndLib::Add(shape, boundingBox);
    double Xmin,Ymin,Zmin,Xmax,Ymax,Zmax;
    boundingBox.Get(Xmin,Ymin,Zmin,Xmax,Ymax,Zmax);
    double delta_X = fabs(Xmax-Xmin);
    double delta_Y = fabs(Ymax-Ymin);
    double delta_Z = fabs(Zmax-Zmin);

    //! --------------
    //! enlarge a bit
    //! --------------
    Xmin = Xmin-delta_X*0.02;
    Xmax = Xmax+delta_X*0.02;
    Ymin = Ymin-delta_Y*0.02;
    Ymax = Ymax+delta_Y*0.02;
    Zmin = Zmin-delta_Z*0.02;
    Zmax = Zmax-delta_Z*0.02;
    */
}
