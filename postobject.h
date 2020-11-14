
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
#include <TColStd_HPackedMapOfInteger.hxx>

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
//#include <QVector>
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
    postObject(const std::map<GeometryTag, std::vector<std::map<int, double> > > &resMap);

    //! ------------
    //! constructor
    //! ------------
    postObject(const std::map<GeometryTag, std::vector<std::map<int, double>>> &resMap, const std::vector<GeometryTag> &aVecLoc);

    //! ------------
    //! constructor
    //! ------------
    postObject(const std::map<GeometryTag,std::vector<std::map<int,double>>> &resMap, const std::vector<GeometryTag> &aVecLoc, const QString& aName);

    //! ------------
    //! constructor
    //! ------------
    postObject(const std::map<GeometryTag,std::vector<std::map<int,double>>> &resMap,
               const std::vector<GeometryTag> &aVecLoc,
               const std::map<GeometryTag,std::map<int,gp_Vec>> mapOfMapOfNodalDiplacements,
               const QString& aName, bool useSurfaceMeshForVolumeResults);

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
        myUseSurfaceMeshForVolumeResults = other.myUseSurfaceMeshForVolumeResults;
    }

    //! -------------------------------------------------
    //! set the map of the vectorial nodal displacements
    //! -------------------------------------------------
    void setMapOfNodalDisplacements(const std::map<GeometryTag,std::map<int,gp_Vec>> &mapDisplMap);

    //! -------------------------------------------------
    //! get the map of the vectorial nodal displacements
    //! -------------------------------------------------
    std::map<GeometryTag,std::map<int,gp_Vec>> getMapOfNodalDisplacements();

private:

    //! -------------------------------------------------
    //! in case of results on a volume mesh tells
    //! if plot only the surface mesh or the volume mesh
    //! -------------------------------------------------
    bool myUseSurfaceMeshForVolumeResults;

    //! ---------
    //! the data
    //! ---------
    std::map<GeometryTag,std::vector<std::map<int,double>>> theData;

    //! -----------------------------------
    //! mesh objects and mesh data sources
    //! -----------------------------------
    std::map<GeometryTag,occHandle(MeshVS_Mesh)> theMeshes;
    std::map<GeometryTag,occHandle(MeshVS_DeformedDataSource)> theMeshDataSources;
    std::map<GeometryTag,occHandle(MeshVS_DeformedDataSource)> theMeshDataSourcesForView;

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
    std::vector<GeometryTag> myVecLoc;

    //! --------------------------------------
    //! map of map of the nodal displacements
    //! --------------------------------------
    std::map<GeometryTag,std::map<int,gp_Vec>> myMapOfNodalDisplacements;

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

    //! ----------------
    //! hidden elements
    //! ----------------
    std::map<GeometryTag,occHandle(TColStd_HPackedMapOfInteger)> myMapOfHiddenElements;

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
    std::map<GeometryTag,std::vector<std::map<int,double>>> getData() { return theData; }

    //! get locations
    std::vector<GeometryTag> getLocations() const { return myVecLoc; }

    //! NbMeshes
    int NbMeshes() const { return (int)theData.size(); }

    //! write
    void write(std::ofstream &file);

    //! read: read only the colors
    void read(std::ifstream &file);

    //! write mesh
    void writeMesh(ofstream &file, const occHandle(MeshVS_DataSource) &theMeshDS);

    //! build mesh IO
    bool buildMeshIO(double min=-1e20, double max=1e20, int Nlevels=10, bool autoscale=true, int component=0, double deformationScale = 1.0);

    //! init - OBSOLETE
    void init(meshDataBase *mDB);
    void setMeshDataSources(const std::map<GeometryTag,occHandle(MeshVS_DataSource)> &aMapOfMeshDS);

    //! get colored meshes
    std::map<GeometryTag,occHandle(MeshVS_Mesh>) getColoredMeshes() const { return theMeshes; }

    //! get the mesh data sources
    std::map<GeometryTag,occHandle(MeshVS_DeformedDataSource)> getMeshDataSources() const { return theMeshDataSources; }

    //! get the mesh data sources for view
    std::map<GeometryTag,occHandle(MeshVS_DeformedDataSource)> getMeshDataSourcesForView() const { return theMeshDataSourcesForView; }

    //! is empty
    bool isEmpty() const{ return (theData.size()==0? true:false); }

    //! update scaled view (internally read myScale)
    void updateScaledView();

    //! compute hidden elements - unused - but do not remove
    void computeHiddenElements(const std::map<int,std::vector<double>> &mapOfClipPlanes);


public:

    //! get the color box
    occHandle(AIS_ColorScaleExtended) getColorBox() const { return AISColorScale; }

    //! get scale
    double getScale() { return myScale; }

    //! get solution data component
    int getSolutionDataComponent() {return mySolutionDataComponent; }

    //! get min/max
    double getMin() { return myMin; }
    double getMax() { return myMax; }

    //! get auto number of levels
    int getNbLevels() { return myNbLevels; }

    //! is autoscale
    bool IsAutoscale() { return myIsAutoscale; }

    //! experimental
    void updateMapping(int mapping);

    //! set scale
    void setScale(double scale) { myScale = scale; }

    //! set mode
    void setMode (bool useSurfaceMeshForVolumeResults) { myUseSurfaceMeshForVolumeResults = useSurfaceMeshForVolumeResults; }

private:

    void writeIntoStream(ofstream &os, const occHandle(MeshVS_DataSource) &aMeshDS);
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
