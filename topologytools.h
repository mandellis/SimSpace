#ifndef TOPOLOGYTOOLS_H
#define TOPOLOGYTOOLS_H

//! custom includes
#include "mydefines.h"
#include <geometrydatabase.h>

//! Qt
#include <QList>
#include <QMap>

//! OCC
#include <TopoDS_Shape.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Compound.hxx>

class TopologyTools
{
public:

    TopologyTools();

    static QVector<GeometryTag> generateLocationPairs(geometryDataBase *gDB, const ListOfShape& scope);

    //! ------------------------------------------
    //! retrieve the main shapes within a compoud
    //! ------------------------------------------
    static TopoDS_Shape exploreCompound(const TopoDS_Compound &aComp,
                                        QList<TopoDS_Shape> &csolids,
                                        QList<TopoDS_Shape> &solids,
                                        QList<TopoDS_Shape> &shells,
                                        QList<TopoDS_Shape> &faces,
                                        QList<TopoDS_Shape> &wires,
                                        QList<TopoDS_Shape> &edges,
                                        bool load3DBodies=true,
                                        bool loadSurfaceBodies=true,
                                        bool loadLineBodies=false);

    //! --------------------------------
    //! make shell from a list of faces
    //! --------------------------------
    static TopoDS_Shell makeShellFromFaces(const QList<TopoDS_Face> &aFaceList);
};

#endif // TOPOLOGYTOOLS_H
