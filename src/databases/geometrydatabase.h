#ifndef GEOMETRYDATABASE_H
#define GEOMETRYDATABASE_H

//! ----------------
//! custom includes
//! ----------------
#include "simulationnodeclass.h"
#include "qextendedstandarditem.h"
#include "src/registeredMetatypes/topods_shape_reg.h"

//! ---
//! Qt
//! ---
#include <QStandardItemModel>
#include <QStandardItem>
#include <QObject>
#include <QList>

//! ----
//! OCC
//! ----
#include <TopoDS_Shape.hxx>
#include <NCollection_Array1.hxx>
#include <TColStd_HArray1OfInteger.hxx>
#include <TColStd_HArray2OfInteger.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TColStd_Array1OfAsciiString.hxx>
#include <TColStd_Array1OfBoolean.hxx>
#include <TopTools_SequenceOfShape.hxx>
#include <TopTools_ListOfShape.hxx>
#include <NCollection_Map.hxx>

#include <Bnd_Box.hxx>
#include <BRepBndLib.hxx>

//! ----------------
//! custom includes
//! ----------------
#include "qprogressevent.h"
#include "qprogressindicator.h"
#include "src/utils/global.h"

class TopoDS_Solid;
class TopoDS_Shell;
class TopoDS_Wire;

//! -------------------------------
//! struct variable: topology maps
//! 1) indexed map of shells
//! 2) indexed map of faces
//! 3) indexed map of edges
//! 4) indexed map of vertices
//! for each body (3D, 2D, 1D)
//! -------------------------------
struct TopologyMap
{
    TopTools_IndexedMapOfShape csolidMap;
    TopTools_IndexedMapOfShape solidMap;
    TopTools_IndexedMapOfShape shellMap;
    TopTools_IndexedMapOfShape faceMap;
    TopTools_IndexedMapOfShape wireMap;
    TopTools_IndexedMapOfShape edgeMap;
    TopTools_IndexedMapOfShape vertexMap;
};
Q_DECLARE_METATYPE (TopologyMap)

//! ----------------
//! body properties
//! ----------------
struct geometryProperties_
{
    Standard_Real area;
    Standard_Real volume;
    Standard_Real xB,yB,zB;     //! center of mass coordinates
    Standard_Real Ixx,Iyy,Izz;  //! principal moments of inertia
};


class geometryDataBase: public QObject
{
    Q_OBJECT

protected:

    //! the file location
    QString fullSourceFileName;

    //! ----------------------------------------------------------------
    //! Geometry content: this is the shape provided by the CAD reader
    //! "<CAD reader>.oneShape()": myTopoDS_Shape is a TopAbs_COMPOUND
    //! (highest entry in the topology)
    //! ----------------------------------------------------------------
    TopoDS_Shape myTopoDS_Shape;

    //! the model containing the simulation nodes
    QStandardItemModel *myModel;

    //! the invisible root item
    QExtendedStandardItem *myInvisibleRootItem;

    //! the model root
    QExtendedStandardItem *myRootItem;

    //! the "Import" item
    QExtendedStandardItem *myImportNode;

    //! the "Geometry" root item
    QExtendedStandardItem *GeometryRootItem;

    //! the "Coordinate systems" root item
    QExtendedStandardItem *CoordinateSystemsRootitem;

protected:

    virtual void createStandardModel();

    //! create geometry nodes
    void createGeometryNodes();

    //! create "Coordinate systems" root item
    void createCoordinateSystemsRoot();

    //! get the topology map
    TopologyMap getTopologyMaps(const TopoDS_Shape &aShape);

public:

    //! return the model
    inline virtual QStandardItemModel *getModel() const { return myModel; }

    //! transfer names
    void transferNames();

    //! get full file name
    inline QString getSourceFilePath() const { return fullSourceFileName; }

    //! progress indicator
    QProgressIndicator *myProgressIndicator;

public:

    //! ------------------------------------
    //! constructor I - with initialization
    //! ------------------------------------
    geometryDataBase(const TopoDS_Shape& theShape = TopoDS_Shape(), const QString &theFilePath ="", QObject *parent=0);

    //! ---------------------------------
    //! constructor II - build from file
    //! ---------------------------------
    geometryDataBase(const QList<SimulationNodeClass*> listOfNodes, const QString &archiveFileName, QObject *parent=0);

    //! destructor
    virtual ~geometryDataBase();

    //! -------------------------------------------------------------
    //! 1) 3D body = free solid (not belonging to a composite solid)
    //! 2) 2D body = free shell, free TopoDS_Face
    //! 3) 1D body = free TopoDS_Wire, free TopoDS_Edge
    //!
    //! experimental - trying to replace the OOC data structures
    //! using something more general
    //! ---------------------------------------------------------
    QMap<int,TopoDS_Shape> bodyMap;

    //! ---------------------------------------------------------
    //! Each one of the previous bodies is numbered in the map
    //! nb={1, 2, ..., Nb}. The nb-th body has a topology map
    //! associated with it: if the body is 3D the map will start
    //! with shells, if is 2D the map will start with faces, ...
    //! ---------------------------------------------------------
    QMap<int,TopologyMap>  MapOfBodyTopologyMap;

    //! ----------------
    //! The parts names
    //! ----------------
    QMap<int,QString> MapOfBodyNames;

    //! ------------------------------
    //! activation/suppression status
    //! ------------------------------
    QMap<int,bool> MapOfIsActive;

    //! ----------------------------------------------------------------------------------------------------------
    //! int => body index within the body map
    //! int => face index of the body within the face topology map
    //! QList<std::pair<int,int>> the list of the twin faces: first index = parent body; second index = twin face
    //! ----------------------------------------------------------------------------------------------------------
    QMap<int,QMap<int,QList<std::pair<int,int>>>> bodyTwinFaces;

    //! ----------------------------------------------------------------------------------------------------------
    //! int => body index within the body map
    //! int => edge index of the body within the edge topology map
    //! QList<std::pair<int,int>> the list of the twin edges: first index = parent body; second index = twin edge
    //! ----------------------------------------------------------------------------------------------------------
    QMap<int,QMap<int,QList<std::pair<int,int>>>> bodyTwinEdges;

    //! ---------------------------------------
    //! Return the number of 3D, 2D, 1D bodies
    //! ---------------------------------------
    inline int N3D() { return myN3D; }
    inline int N2D() { return myN2D; }
    inline int N1D() { return myN1D; }

    //! -----------------
    //! obtain the shape
    //! -----------------
    inline void shape(TopoDS_Shape &theShape) { theShape = myTopoDS_Shape; }
    inline TopoDS_Shape shape() const { return myTopoDS_Shape; }

    //! --------------
    //! set the shape
    //! --------------
    virtual void update(const TopoDS_Shape &aShape);

    //! ---------------
    //! reset database
    //! ---------------
    virtual void resetDataBase();

    //! -----------------------
    //! set progress indicator
    //! -----------------------
    void setProgressIndicator(QProgressIndicator *aProgressIndicator) { myProgressIndicator = aProgressIndicator; }

private:

    //! ---------------------
    //! create geometry root
    //! ---------------------
    void createGeometryRoot();

    //! -----------------------------
    //! create the "Model" root item
    //! -----------------------------
    void createModelRootItem();

    //! -------------------------
    //! create the "Import" item
    //! -------------------------
    void createImportItem();

    //! -----------
    //! build maps
    //! -----------
    void buildMaps(const TopoDS_Shape &shape);

    //! ----------------------------
    //! Number of 3D, 2D, 1D bodies
    //! ----------------------------
    int myN3D, myN2D, myN1D;

    //! -------------------------------
    //! Print a summary of what loaded
    //! -------------------------------
    void printSummary();

    //! -----------------
    //! get bounding box
    //! -----------------
    void getBoundingBox(const TopoDS_Shape &shape, double &L1, double &L2, double &L3);
    //double boundingBox(const TopoDS_Shape &shape);

    //! ----------------------
    //! tolerance for healing
    //! ----------------------
    double toleranceForHealing(const TopoDS_Shape &shape, int typeOfTolerance);

public:

    //! get the number of the subshape contained in a shape
    bool getSubShapeNr(const TopoDS_Shape &aSubShape, int &mainShapeIndex, int &subShapeIndex, TopAbs_ShapeEnum &subShapeType) const;
    bool getSubShapeNr1(const TopoDS_Shape &aSubShape, int &mainShapeIndex, int &subShapeIndex, TopAbs_ShapeEnum &subShapeType) const;

    //! build the interactive objects
    void buildIOs();

    //! --------------------------
    //! function: boundingBox
    //! details:  return the diag
    //! --------------------------
    double boundingBox(const TopoDS_Shape &shape)
    {
        Bnd_Box boundingBox;
        BRepBndLib::Add(shape, boundingBox);
        Standard_Real Xmin,Ymin,Zmin,Xmax,Ymax,Zmax;
        boundingBox.Get(Xmin,Ymin,Zmin,Xmax,Ymax,Zmax);
        double L1 = fabs(Xmax-Xmin);
        double L2 = fabs(Ymax-Ymin);
        double L3 = fabs(Zmax-Zmin);
        double diag = sqrt(pow(L1,2)+pow(L2,2)+pow(L3,2));
        return diag;
    }

public:

    virtual void removeBodyFromDataBase(int bodyIndex);

signals:

};

#endif // GEOMETRYDATABASE_H
