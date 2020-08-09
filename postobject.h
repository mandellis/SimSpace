#ifndef POSTOBJECT_H
#define POSTOBJECT_H

//! ----------------
//! custom includes
//! ----------------
#include "ais_colorscaleextended.h"
#include "mapofmeshdatasources.h"
#include <meshdatabase.h>
#include <simulationdata.h>

//! ----
//! OCC
//! ----
#include <MeshVS_Mesh.hxx>
#include <MeshVS_DataSource.hxx>
#include <MeshVS_DeformedDataSource.hxx>
#include <NCollection_TListIterator.hxx>
#include <MeshVS_NodalColorPrsBuilder.hxx>
#include <MeshVS_DataSource.hxx>
#include <TColStd_PackedMapOfInteger.hxx>
#include <TColStd_MapIteratorOfPackedMapOfInteger.hxx>
#include <gp_Vec.hxx>

//! --------------------------------------
//! replacement for opencascade::handle<>
//! --------------------------------------
#include "occhandle.h"

//! ---
//! Qt
//! ---
#include <QMetaType>
#include <QString>
#include <QMap>
#include <QVector>
#include <Quantity_Color.hxx>

class postObject
{

public:

    //! --------------------
    //! default constructor
    //! --------------------
    postObject();

    //! ------------
    //! constructor
    //! ------------
    postObject(const QMap<GeometryTag, QList<QMap<int, double>>> &resMap);

    //! ------------
    //! constructor
    //! ------------
    postObject(const QMap<GeometryTag, QList<QMap<int, double>>> &resMap, const QVector<GeometryTag> &aVecLoc);

    //! ------------
    //! constructor
    //! ------------
    postObject(const QMap<GeometryTag,QList<QMap<int,double>>> &resMap, const QVector<GeometryTag> &aVecLoc, const QString& aName);

    //! ------------
    //! constructor
    //! ------------
    postObject(const QMap<GeometryTag,QList<QMap<int,double>>> &resMap,
               const QVector<GeometryTag> &aVecLoc,
               const QMap<GeometryTag,QMap<int,gp_Vec>> mapOfMapOfNodalDiplacements,
               const QString& aName, bool aShowSolidMeshAsSurface);

    //! ----------------------
    //! constructor from file
    //! ----------------------
    postObject(ifstream &file);

    //! -----------
    //! destructor
    //! -----------
    ~postObject()
    {
        //cout<<"postObject:~postObject()->____DESTRUCTOR CALLED____"<<endl;
    }

    //! -----------------
    //! copy constructor
    //! -----------------
    postObject(const postObject &other)
    {
        theMeshes = other.theMeshes;                    // check if the operator = id defined for the class
        theMeshDataSources = other.theMeshDataSources;  // check if the operator = id defined for the class
        AISColorScale = other.AISColorScale;            // check if the operator = id defined for the class
        name = other.name;
        myVecLoc = other.myVecLoc;
        theData = other.theData;
        myMapOfNodalDisplacements = other.myMapOfNodalDisplacements;
        mySolutionDataComponent = other.mySolutionDataComponent;
        myIsAutoscale = other.myIsAutoscale;
        myMin = other.myMin;
        myMax = other.myMax;
        myNbLevels = other.myNbLevels;
        myShowSolidMeshAsSurface = other.myShowSolidMeshAsSurface;
    }

    //! -------------------------------------------------
    //! set the map of the vectorial nodal displacements
    //! -------------------------------------------------
    void setMapOfNodalDisplacements(const QMap<GeometryTag,QMap<int,gp_Vec>> &mapDisplMap)
    {
        myMapOfNodalDisplacements = mapDisplMap;
    }

    //! -------------------------------------------------
    //! get the map of the vectorial nodal displacements
    //! -------------------------------------------------
    QMap<GeometryTag,QMap<int,gp_Vec>> getMapOfNodalDisplacements()
    {
        return myMapOfNodalDisplacements;
    }

private:

    //! -------------------------------------------------
    //! in case of results on a volume mesh tells
    //! if plot only the surface mesh or the volume mesh
    //! -------------------------------------------------
    bool myShowSolidMeshAsSurface;

    //! ---------
    //! the data
    //! ---------
    QMap<GeometryTag,QList<QMap<int,double>>> theData;

    //! -------------------------------------
    //! MeshVS_Mesh objects (colored meshes)
    //! -------------------------------------
    QMap<GeometryTag,occHandle(MeshVS_Mesh)> theMeshes;
    QMap<GeometryTag,occHandle(MeshVS_DeformedDataSource)> theMeshDataSources;
    QMap<GeometryTag,occHandle(MeshVS_DeformedDataSource)> theMeshDataSourcesForView;

    //! ----------------
    //! the color scale
    //! ----------------
    occHandle(AIS_ColorScaleExtended) AISColorScale;

    //! ----------------------------
    //! the name of the post object
    //! ----------------------------
    QString name;

    //! ----------
    //! locations
    //! ----------
    QVector<GeometryTag> myVecLoc;

    //! --------------------------------------
    //! map of map of the nodal displacements
    //! --------------------------------------
    QMap<GeometryTag,QMap<int,gp_Vec>> myMapOfNodalDisplacements;

    //! ----------------------------
    //! the solution data component
    //! ----------------------------
    int mySolutionDataComponent;
    bool myIsAutoscale;
    double myMin,myMax;
    int myNbLevels;

    //! ------------------
    //! deformation scale
    //! ------------------
    double myScale;

private:

    //! reset meshes
    void resetMeshes();

    //! get min max
    std::pair<double,double> getMinMax(int component);

public:

    //! set name
    void setName(const QString &aName) {name = aName; }

    //! get name
    QString getName() const { return name;}

    //! get data
    QMap<GeometryTag,QList<QMap<int,double>>> getData() { return theData; }

    //! get locations
    QVector<GeometryTag> getLocations() const { return myVecLoc; }

    //! NbMeshes
    int NbMeshes() const { return theData.size(); }

    //! write
    void write(std::ofstream &file);

    //! read: read only the colors
    void read(std::ifstream &file);

    //! write mesh
    void writeMesh(ofstream &file, const occHandle(MeshVS_DataSource) &theMeshDS);

    //! rebuild "theMeshes" with data
    //void buildMeshIO(const mapOfMeshDataSources &aMapOfMeshDataSources,
    //                 double min=-1e20, double max=1e20, int Nlevels=9, bool autoscale=true, int component=0);

    void buildMeshIO(double min=-1e20, double max=1e20, int Nlevels=9, bool autoscale=true, int component=0, double deformationScale = 1.0);

    //! init
    void init(meshDataBase *mDB, int component=0);

    //! update mesh
    //void updateView(bool showMeshEdges);

    //! clone
    postObject clone(const postObject &other);

    //! get colored meshes
    QMap<GeometryTag,occHandle(MeshVS_Mesh>) getColoredMeshes() const { return theMeshes; }

    //! get the mesh data sources
    QMap<GeometryTag,occHandle(MeshVS_DeformedDataSource)> getMeshDataSources() const { return theMeshDataSources; }

    //! is empty
    bool isEmpty() const{ return theData.isEmpty(); }

    //! update scaled view (read myScale)
    void updateScaledView();

public:

    //! get the color box
    occHandle(AIS_ColorScaleExtended) getColorBox() const { return AISColorScale; }

    //! set scale
    void setScale(double scale) { myScale = scale; }

    //! get scale
    double getScale() { return myScale; }

    //! get solution data component
    int getSolutionDataComponent() {return mySolutionDataComponent; }

    //! get min/max
    double getMin() { return myMin; }
    double getMax() { return myMax; }

    //! get auto number of levels
    int getNbLevels() { return myNbLevels; }

    //! experimental
    void updateMapping(int mapping);

private:

    void writeIntoStream(ofstream &os, const opencascade::handle<MeshVS_DataSource> &aMeshDS);
    bool readMeshFromStream(ifstream &stream, occHandle(MeshVS_DataSource) &aMeshDS);

    //! -------
    //! helper
    //! -------
    int hueFromValue(int theValue,int theMin,int theMax)
    {
        int aMinLimit (0), aMaxLimit (230);
        int aHue = aMaxLimit;
        if (theMin!=theMax) aHue = (int)(aMaxLimit -(aMaxLimit-aMinLimit)*(theValue-theMin)/(theMax - theMin));
        aHue = std::min (std::max (aMinLimit, aHue), aMaxLimit);
        return aHue;
    }
};

Q_DECLARE_METATYPE(postObject)


typedef std::shared_ptr<postObject> sharedPostObject;
Q_DECLARE_METATYPE(sharedPostObject)

#endif // POSTOBJECT_H
