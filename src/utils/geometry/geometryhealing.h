#ifndef GEOMETRYHEALING_H
#define GEOMETRYHEALING_H

//! redefinition of opencascade Handle() macro
#include "occhandle.h"

#include <QObject>

//! OCC
#include <TopoDS_Shape.hxx>

//! minimum edge
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>

class geometryHealing: public QObject
{

    Q_OBJECT

public:

    geometryHealing(QObject *parent=0);
    geometryHealing(const TopoDS_Shape &aShape, QObject *parent=0);
    inline void setShape(const TopoDS_Shape &aShape) { myTopoDS_Shape = aShape; }

    void perform(double Prec);

    void removeDegenerateEdges(TopoDS_Shape &shape);
    void fixFaces(TopoDS_Shape &shape);

    inline TopoDS_Shape getShape() { return myTopoDS_Shape; }

private:

    //! geometry
    TopoDS_Shape myTopoDS_Shape;

private:

    //! minimum edge lenght
    double getMinimumEdge(const TopoDS_Shape &aShape);
};


#endif // GEOMETRYHEALING_H
