#ifndef MESHDATABASE_H
#define MESHDATABASE_H

#include <QMap>
//#include <map>
#include "src/utils/hash_c.h"
//#include <utility>

struct doubleKey
{
    int key1;
    int key2;

    doubleKey(int akey1, int akey2) { key1= akey1; key2= akey2;}
    doubleKey(const doubleKey &other) { key1 = other.key1; key2 = other.key2; }
    inline doubleKey* operator = (const doubleKey& other) { key1 = other.key1; key2 = other.key2; return this; }

    inline bool operator == (const doubleKey &other) const { if(key1==other.key1 && key2==other.key2) return true; return false; }

    inline bool operator < (const doubleKey& other) const
    {
        size_t seed1 = 0;
        hash_c<int>(seed1,key1);
        hash_c<int>(seed1,key2);

        size_t seed2 = 0;
        hash_c<int>(seed2,other.key1);
        hash_c<int>(seed2,other.key2);

        if(seed1<seed2) return true;
        return false;
    }
};

/* an equivalent definition
template <class Q>
struct Q2DMap
{
    QMap<doubleKey,Q> innerMap;
    void setValue(int akey1, int akey2, Q value)
    {
        doubleKey dk(akey1, akey2);
        innerMap.insert(dk,value);
    }
    inline Q getValue(int akey1, int akey2) { return innerMap.value(doubleKey(akey1,akey2)); }
    inline void clear() { innerMap.clear(); }
};
*/

/*
//! -------------------------
//! use C++ standard library
//! -------------------------
template <class T>
class Q2DMap: public std::map<doubleKey,T>
{
public:

    void setValue(int aKey1, int aKey2, T value)
    {
        doubleKey dk(aKey1,aKey2);
        std::map<doubleKey,T>::iterator it = this->find(dk);
        if(it==this->end())
        {
            std::pair<doubleKey,T> aPair;
            aPair.first = dk;
            aPair.second = T;
            this->insert(aPair);
        }
        else it->second = value;
    }
    T getValue(int aKey1, int aKey2)
    {
        doubleKey dk(aKey1,aKey2);
        std::map<doubleKey,T>::iterator it = this->find(dk);
        return it->second;
    }
    void remove(const doubleKey &aKey)
    {
        this->erase(aKey);
    }
};
*/

template <class T>
class Q2DMap:public QMap<doubleKey,T>
{
public:

    void setValue(int aKey1, int aKey2, T value)
    {
        doubleKey dk(aKey1,aKey2);
        this->insert(dk, value);
    }
    T getValue(int aKey1, int aKey2)
    {
        doubleKey dk(aKey1,aKey2);
        return this->value(dk);
    }
};

//! ----------------
//! custom includes
//! ----------------
#include <geometrydatabase.h>
#include "src/utils/myenumvariables.h"
#include "qextendedstandarditem.h"
#include "property.h"
#include "src/main/mydefines.h"
#include "prismaticlayerparameters.h"
#include <ng_meshvs_datasourceface.h>

//! ----
//! OCC
//! ----
#include <TColStd_Array1OfReal.hxx>
#include <TColStd_Array2OfReal.hxx>
#include <NCollection_Array1.hxx>
#include <NCollection_Array2.hxx>
#include <TColStd_Array2OfReal.hxx>
#include <TColStd_Array1OfBoolean.hxx>
#include <TColStd_Array2OfBoolean.hxx>
#include <TopTools_SequenceOfShape.hxx>
#include <TopTools_ListOfShape.hxx>
#include <MeshVS_DataSource.hxx>

//! ----
//! C++
//! ----
#include <vector>

//! ---
//! Qt
//! ---
#include <QObject>
#include <QVector>

//! ----------------
//! custom includes
//! ----------------
#include "geometrytag.h"

class TopoDS_Shape;
class MeshVS_DataSource;

class meshDataBase: public geometryDataBase
{

    Q_OBJECT

protected:

    //! the mesh root node
    QExtendedStandardItem *MeshRootItem;

    //! create the mesh root node
    void createMeshRootNode();

    //! create the standard model
    virtual void createStandardModel();

private:

    //! initialize the array of settings
    void initializeMeshArrays();

public:

    //! return the model
    virtual inline QStandardItemModel *getModel() const override { return myModel; }

public:

    //! constructor I - default, not used at the momoment
    //meshDataBase(QObject *parent=0);

    //! constructor II - constructor with initialization
    explicit meshDataBase(const TopoDS_Shape &shape = TopoDS_Shape(), const QString &theFilePath ="", QObject *parent=0);

    //! constructor III - constructor from a list of nodes
    explicit meshDataBase(const QList<SimulationNodeClass*> listOfNodes, const QString &archiveFileName, QObject *parent=0);

    //! destructor
    virtual ~meshDataBase();

    //! update
    virtual void update(const TopoDS_Shape &aShape);

    //! reset database
    virtual void resetDataBase();

    //! -------------------
    //! submesh generation
    //! -------------------
    bool subMeshGeneration = true;

    //! ------------------------------
    //! face/edge/vertex element size
    //! ------------------------------
    Q2DMap<double> MapOfElementSizeOnFace;
    Q2DMap<double> MapOfElementSizeOnEdge;
    Q2DMap<int> MapOfNumberOfDivisionOnEdge;
    Q2DMap<bool> MapOfSizingTypeOnEdge;
    Q2DMap<double> MapOfElementSizeOnVertex;
    Q2DMap<double> MapOfVertexPinball;

    //! ----------------
    //! prismatic faces
    //! ----------------
    QMap<int,bool> HasPrismaticFaces;
    QMap<int,bool> defaultHasPrismaticFaces;

    QMap<int,QList<int>> prismaticFaces;        //! key: body index; value: QList<int> list of faces having inflation defined on
    QMap<int,prismaticLayerParameters> prismaticMeshParameters;     //! key: body index; value:  QList<prismaticLayerParameters> prismatic layers parameters on that body

    //! ------------------------------------
    //! face/edge/vertex modification flags
    //! ------------------------------------
    Q2DMap<bool> MapOfIsFaceModified;
    Q2DMap<bool> MapOfIsEdgeModified;
    Q2DMap<bool> MapOfIsVertexModified;

    //! ---------------------------------
    //! healing/defeaturing and defaults
    //! ---------------------------------
    QMap<int,Property::meshEngine2D> mapOfDiscretizer;
    QMap<int,double> mapOfAngularDeflection;
    QMap<int,double> mapOfLinearDeflection;
    QMap<int,double> mapOfMinFaceSize;              //! property of the surface tassellator
    QMap<int,double> mapOfMaxFaceSize;              //! property of the surface tassellator
    QMap<int,bool> MapOfIsBodyDefeaturingOn;
    QMap<int,bool> MapOfIsBodyHealingOn;
    QMap<int,bool> MapOfIsMeshSimplificationOn;
    QMap<int,int> MapOfMeshSimplificationBy;
    QMap<int,double> MapOfBodyDefeaturingParameterValue;

    QMap<int,bool> defaultIsBodyDefeaturingOn;
    QMap<int,bool> defaultIsBodyHealingOn;
    QMap<int,bool> defaultIsMeshSimplificationOn;
    QMap<int,double> defaultMeshDefeaturingParameterValue;
    QMap<int,int> defaultMeshDefeaturingBy;
    QMap<int,double> defaultMinFaceSize;            //! default property of the surface tassellator
    QMap<int,double> defaultMaxFaceSize;            //! default property of the surface tassellator

    //! ----------------------------------
    //! meshing process in memory/on disk
    //! ----------------------------------
    QMap<int,bool> MapOfIsMeshingRunningInMemory;
    QMap<int,bool> defaultIsMeshingRunningInMemory;

    //! ------------------------
    //! feature preserving flag
    //! ------------------------
    QMap<int,bool> mapOfFeaturePreserving;
    QMap<int,bool> defaultMapOfFeaturePreserving;
    std::vector<GeometryTag> featuredTags;

    //! ----------------------------------------------------
    //! geometry correction after patch independent meshing
    //! ----------------------------------------------------
    QMap<int,bool> mapOfGeometryCorrection;
    QMap<int,bool> defaultMapOfGeometryCorrection;

    //! ------------------------------
    //! new TetWild mesher parameters
    //! ------------------------------
    QMap<int,int> mapOfEnvelopeSizingType;
    QMap<int,double> mapOfRelativeEnvelopeSize;
    QMap<int,double> mapOfAbsoluteEnvelopeSize;

    QMap<int,int> mapOfIdealLengthSizingType;
    QMap<int,double> mapOfIdealLengthRelativeSize;
    QMap<int,double> mapOfIdealLengthAbsoluteSize;

    QMap<int,int> defaultMapOfEnvelopeSizingType;
    QMap<int,double> defaultMapOfRelativeEnvelopeSize;
    QMap<int,double> defaultMapOfAbsoluteEnvelopeSize;

    QMap<int,int> defaultMapOfIdealLengthSizingType;
    QMap<int,double> defaultMapOfIdealLengthRelativeSize;
    QMap<int,double> defaultMapOfIdealLengthAbsoluteSize;

    //! ----------------------
    //! The mesh data sources
    //! ----------------------
    QMap<int,occHandle(MeshVS_DataSource)> ArrayOfMeshDS;
    QMap<int,occHandle(MeshVS_DataSource)> ArrayOfMeshDS2D;
    Q2DMap<occHandle(MeshVS_DataSource)> ArrayOfMeshDSOnFaces;
    Q2DMap<occHandle(MeshVS_DataSource)> ArrayOfMeshDSOnEdges;
    
    //! experimental
    Q2DMap<occHandle(MeshVS_DataSource)> ArrayOfPrismaticMeshes;
    Q2DMap<occHandle(MeshVS_DataSource)> ArrayOfNonPrismaticMeshes;
    
    //! ------------------------------
    //!  mesh parameters and defaults
    //! ------------------------------
    QMap<int,bool> mapOfUseBRep;
    QMap<int,Property::meshEngine2D> ArrayOfMesh2DEngine;
    QMap<int,Property::meshEngine3D> ArrayOfMesh3DEngine;
    QMap<int,Property::meshOrder> ArrayOfMeshOrder;
    QMap<int,double> ArrayOfGradingValue;
    QMap<int,double> ArrayOfMinBodyElementSize;
    QMap<int,double> ArrayOfMaxBodyElementSize;
    QMap<int,bool> mapOfIsElementStraight;
    QMap<int,int> ArrayOfSmoothingSteps;

    QMap<int,bool> defaultMapOfUseBRep;
    QMap<int,Property::meshEngine2D> default2DMeshEngine;
    QMap<int,Property::meshEngine3D> default3DMeshEngine;
    QMap<int,Property::meshOrder> defaultMeshOrder;
    QMap<int,double> defaultVolumeGrading;
    QMap<int,double> defaultVolumeMinElementSize;
    QMap<int,double> defaultVolumeMaxElementSize;
    QMap<int,bool> defaultIsElementStraight;
    QMap<int,int> defaultSmoothingSteps;

    QMap<int,Property::meshEngine2D> defaultMapOfDiscretizer;
    QMap<int,double> defaultMapOfAngularDeflection;
    QMap<int,double> defaultMapOfLinearDeflection;

    //! -----------------
    //! mesh update flag
    //! -----------------
    QMap<int,bool> ArrayOfMeshIsToBeUdpdated;

    //! ------------------------------------------------------------
    //! the type of mesh for each body - surface mesh & volume mesh
    //! ------------------------------------------------------------
    //QMap<int,Property::meshType_Surface> ArrayOfSurfaceMeshType;
    //QMap<int,Property::meshType_Volume> ArrayOfVolumeMeshType;
    //QMap<int,Property::meshType_Surface> defaultSurfaceMeshType;
    //QMap<int,Property::meshType_Volume> defaultVolumeMeshType;

    QMap<int,int> ArrayOfMeshType;
    QMap<int,int> ArrayOfDefaultMeshType;

    //! --------------------------------------
    //! set/calculate the meshing paramenters
    //! --------------------------------------
    void calculateDefaultMeshingParameters();

    //! ----------------------------------
    //! apply the default mesh parameters
    //! ----------------------------------
    void applyDefaultMeshingParameters();

    //! -----------------------
    //! get meshing parameters
    //! -----------------------
    meshParam getMeshingParameters(int bodyIndex);

    void NmaxEntity(int &NMaxFace, int &NMaxEdge, int &NMaxVertex);
    void shapeSize(const TopoDS_Shape &aShape, double &diag, double &minlength);

public:

    //! ------------------------------------------------------
    //! remote body from database - acts on mesh data sources
    //! ------------------------------------------------------
    virtual void removeBodyFromDataBase(int bodyIndex)
    {
        geometryDataBase::removeBodyFromDataBase(bodyIndex);

        //! -------------------------------
        //! remove volume and surface mesh
        //! -------------------------------
        ArrayOfMeshDS.remove(bodyIndex);
        ArrayOfMeshDS.insert(bodyIndex,occHandle(MeshVS_DataSource)());
        ArrayOfMeshDS2D.remove(bodyIndex);
        ArrayOfMeshDS2D.insert(bodyIndex,occHandle(MeshVS_DataSource)());

        //! -------------------
        //! flag to be updated
        //! -------------------
        ArrayOfMeshIsToBeUdpdated.insert(bodyIndex,true);

        //! -------------------
        //! remove face meshes
        //! -------------------
        for(QMap<int,TopologyMap>::iterator it = MapOfBodyTopologyMap.begin(); it != MapOfBodyTopologyMap.end(); it++)
        {
            int curBodyIndex = it.key();
            if(curBodyIndex == bodyIndex)
            {
                TopTools_IndexedMapOfShape faceMap = it.value().faceMap;
                for(int faceNr=1; faceNr<=faceMap.Extent(); faceNr++)
                {
                    doubleKey aDoubleKey(bodyIndex,faceNr);
                    ArrayOfMeshDSOnFaces.remove(aDoubleKey);
                    ArrayOfMeshDSOnFaces.insert(aDoubleKey,occHandle(MeshVS_DataSource)());
                }
            }
        }
        //! -------------------
        //! remove edge meshes
        //! -------------------
        for(QMap<int,TopologyMap>::iterator it = MapOfBodyTopologyMap.begin(); it != MapOfBodyTopologyMap.end(); it++)
        {
            int curBodyIndex = it.key();
            if(curBodyIndex == bodyIndex)
            {
                TopTools_IndexedMapOfShape edgeMap = it.value().edgeMap;
                for(int edgeNr=1; edgeNr<=edgeMap.Extent(); edgeNr++)
                {
                    doubleKey aDoubleKey(bodyIndex,edgeNr);
                    ArrayOfMeshDSOnEdges.remove(aDoubleKey);
                    ArrayOfMeshDSOnEdges.insert(aDoubleKey,occHandle(MeshVS_DataSource)());
                }
            }
        }
    }

public slots:

};

#endif // MESHDATABASE_H
