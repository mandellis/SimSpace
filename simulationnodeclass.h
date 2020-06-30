#ifndef SIMULATIONNODECLASS_H
#define SIMULATIONNODECLASS_H

//! ---
//! Qt
//! ---
#include <QObject>
#include <QString>
#include <QVector>
#include <QMap>
#include <QMetaType>

//! ----
//! OCC
//! ----
#include <TopTools_Array1OfListOfShape.hxx>

//! ----------------
//! custom includes
//! ----------------
#include "listofshape.h"
#include "myenumvariables.h"
#include "property.h"

class QStandardItemModel;
class QStandardItem;
class load;
class CustomTableModel;
class QExtendedStandardItem;

class SimulationNodeClass: public QObject
{
    Q_OBJECT

public:

    enum nodeType
    {
        nodeType_NULL,                         //! NULL node
        nodeType_root,                         //! root of the tree (not the invisible item root)

        nodeType_import,                       //! geometry importing options

        nodeType_meshControl,                  //! the family of mesh controls nodes
        nodeType_meshGrading,                  //! grading Netgen
        nodeType_meshMinElementSize,           //! min body element size
        nodeType_meshMaxElementSize,           //! max body element size
        nodeType_meshBodyMeshControl,          //! grading, min, max
        nodeType_meshBodyMeshMethod,           //! 2D/3D engine, 2D/3D elements, mesh order
        nodeType_meshMethod,                   //! mesh method
        nodeType_meshFaceSize,                 //! face size
        nodeType_meshEdgeSize,                 //! edge size
        nodeType_meshVertexSize,               //! size around a point
        nodeType_meshPrismaticLayer,           //! prismatic layer
        nodeType_meshGenericSizing,            //! generic sizing control
        nodeType_meshMeshType,                 //! type of mesh (full tri/tet, full quad/hexa, quad/hexa dominant)
        nodeType_meshMeshMetric,               //! mesh metric

        nodeType_connection,                    //! the family or the root node
        nodeType_connectionPair,                //! a connection pair
        nodeType_connectionGroup,               //! a connection group

        nodeType_geometry,                      //! the family or the root node
        nodeType_geometryBody,                  //! a body
        nodeType_geometryPart,                  //! a part made of bodies
        nodeType_repairTool,                    //! a repair tool

        //! -----------------
        //! named selections
        //! -----------------
        nodeType_namedSelection,                //! the family of named selection
        nodeType_namedSelectionGeometry,        //! named selection made of geometry
        nodeType_namedSelectionNodal,           //! nodal selection made of nodes
        nodeType_namedSelectionElement,         //! named selection made of element

        //! --------------------------
        //! thermal environment stuff
        //! --------------------------
        nodeType_thermalAnalysis,
        nodeType_thermalAnalysisSettings,
        nodeType_thermalAnalysisTemperature,
        nodeType_thermalAnalysisConvection,
        nodeType_thermalAnalysisThermalFlux,
        nodeType_thermalAnalysisThermalFlow,
        nodeType_thermalAnalysisThermalPower,
        nodeType_thermalAnalysisRadiation,
        nodeType_thermalAnalysisAdiabaticWall,

        //! --------------------------------
        //! thermal analysis solution nodes
        //! --------------------------------
        nodeType_thermalAnalysisSolution,
        nodeType_thermalAnalysisSolutionInformation,
        nodeType_solutionThermalTemperature,
        nodeType_solutionThermalFlux,

        //! -----------------------------
        //! structural environment stuff
        //! -----------------------------
        nodeType_structuralAnalysis,
        nodeType_structuralAnalysisSettings,
        nodeType_structuralAnalysisBoundaryContidion_FixedSupport,
        nodeType_structuralAnalysisBoundaryCondition_CylindricalSupport,
        nodeType_structuralAnalysisBoundaryCondition_FrictionlessSupport,
        nodeType_structuralAnalysisBoundaryCondition_CompressionOnlySupport,
        nodeType_structuralAnalysisBoundaryCondition_Displacement,
        nodeType_structuralAnalysisBoundaryCondition_Moment,
        nodeType_structuralAnalysisBoundaryCondition_Pressure,
        nodeType_structuralAnalysisBoundaryCondition_Acceleration,
        nodeType_structuralAnalysisBoundaryCondition_RotationalVelocity,
        nodeType_structuralAnalysisBoundaryCondition_RemoteForce,
        nodeType_structuralAnalysisBoundaryCondition_Force,
        nodeType_structuralAnalysisBoundaryCondition_RemoteDisplacement,
        nodeType_structuralAnalysisBoundaryCondition_RemoteRotation,
        nodeType_structuralAnalysisBoltPretension,
        nodeType_structuralAnalysisThermalCondition,
        nodeType_structuralAnalysisBoundaryCondition_ImportedTemperatureDistribution,
        nodeType_pointMass,
        nodeType_mapper,                            //! mapper3D
        nodeType_importedBodyScalar,                //! imported scalar
        nodeType_modelChange,                       //! model change

#ifdef COSTAMP_VERSION
        nodeType_processParameters,
        nodeType_processParametersPressure,
        nodeType_processParametersClosureForce,
#endif
        nodeType_coordinateSystems,                 //! the family or the root node
        nodeType_coordinateSystem,                  //! a coordinate system, 64
        nodeType_coordinateSystem_global,           //! the global coordinate system, 65

        nodeType_StructuralAnalysisSolution,                        //! the family of the post-processing items
        nodeType_StructuralAnalysisSolutionInformation,             //! the solution information
        nodeType_solutionStructuralNodalDisplacement,                         //! nodal displacement
        nodeType_solutionStructuralStress,                    //! stress
        nodeType_solutionStructuralTotalStrain,               //! total strain
        nodeType_solutionStructuralMechanicalStrain,          //! mechanical strain
        nodeType_solutionStructuralThermalStrain,             //! thermal strain (as difference)
        nodeType_solutionStructuralTemperature,               //! temperature
        nodeType_solutionStructuralEquivalentPlasticStrain,   //! equivalent plastic strain
        nodeType_solutionStructuralNodalForces,               //! nodal forces
        nodeType_solutionStructuralFatigueTool,               //! fatigue tool
        nodeType_solutionStructuralContact,                   //! contact

#ifdef COSTAMP_VERSION
        nodeType_timeStepBuilder,                   //! time step builder
#endif
        nodeType_remotePointRoot,                   //! the family of remote points
        nodeType_remotePoint,                       //! a remote point

        nodeType_OpenFoamScalarData,                 //! open foam scalar data translator
        nodeType_postObject,                         //! a generic post object

        nodeType_command,                             //! a snippet

        //! ------------------
        //! combined analysis
        //! ------------------
        nodeType_combinedAnalysis,
        nodeType_combinedAnalysisSettings,
        nodeType_combinedAnalysisSolution,
        nodeType_combinedAnalysisSolutionInformation,

        //! --------------------
        //! particles in fields
        //! --------------------
        nodeType_particlesInFieldsAnalysis,
        nodeType_particlesInFieldsAnalysisSettings,
        nodeType_particlesInFieldsSolution,
        nodeType_particlesInFieldsSolutionInformation,
        nodeType_electrostaticPotential,
        nodeType_magneticField,
        nodeType_particlesInFieldsParticlePack
    };

    Q_ENUM(nodeType)

public:

    //! constuctor I
    SimulationNodeClass(QObject *parent=0);

    //! constructor II
    SimulationNodeClass(const QString &aName, const nodeType &aType,
                        const QVector<Property> &props, QObject *parent=0);

    //! copy constructor
    SimulationNodeClass(const SimulationNodeClass &other);

    //! destructor
    ~SimulationNodeClass();

    //! type - return the type in a QString format for serialization
    const QString type() const;

    //! setType - set the type - for serialization
    void setType1(QString type);

    //! type to number
    static int typeToInt(const QString &type);

    //! has tabular data
    bool hasTabularData();

    //! get tabular data
    CustomTableModel* getTabularDataModel();

    //! set valid
    void setValid(bool isValid = true)
    {
        Property prop_isValid("Is valid",isValid,Property::PropertyGroup_Validity);
        if(this->getPropertyItem("Is valid")==NULL) this->addProperty(prop_isValid);
        else this->replaceProperty("Is valid",prop_isValid);
    }

    //! add time tag
    void addTimeTag();

private:

    //! the node name
    QString myNodeName;

    //! the type
    nodeType myNodeType;

    //! the values
    QVector<Property> myProps;

    //! the model
    QStandardItemModel *myNodeModel;

    //! the root of the model
    QExtendedStandardItem *myNodeRootItem;

    //! the tabular data
    CustomTableModel *theTabularDataModel;

    //! the separators
    QExtendedStandardItem *itemDefinition;
    QExtendedStandardItem *itemBoundingBox;
    QExtendedStandardItem *itemProperties;
    QExtendedStandardItem *itemStatistics;
    QExtendedStandardItem *itemMaterial;
    QExtendedStandardItem *itemGraphicProperties;
    QExtendedStandardItem *itemScope;
    QExtendedStandardItem *itemStepControls;
    QExtendedStandardItem *itemSolverControl;
    QExtendedStandardItem *itemLighting;
    QExtendedStandardItem *itemAutoDetection;
    QExtendedStandardItem *itemTransparency;
    QExtendedStandardItem *itemDisplay;
    QExtendedStandardItem *itemAdvanced;
    QExtendedStandardItem *itemOrigin;
    QExtendedStandardItem *itemDirectionalData;
    QExtendedStandardItem *itemTransformations;
    QExtendedStandardItem *itemInformation;
    QExtendedStandardItem *itemValidity;
    QExtendedStandardItem *itemEnvironmentTemperature;
    QExtendedStandardItem *itemMeshViewOptions;
    QExtendedStandardItem *itemSizing;
    QExtendedStandardItem *itemGraphicObject;
    QExtendedStandardItem *itemOutputSettings;
    QExtendedStandardItem *itemBackground;
    QExtendedStandardItem *itemSolutionInformation;
    QExtendedStandardItem *itemConvergenceCriteria;
    QExtendedStandardItem *itemTimeIncrementation;
    QExtendedStandardItem *itemCutBack;
    QExtendedStandardItem *itemLineSearch;
    QExtendedStandardItem *itemHidden;
    QExtendedStandardItem *itemColorBox;
    QExtendedStandardItem *itemDefeaturing;
    QExtendedStandardItem *itemMethod;
    //QExtendedStandardItem *itemShrink;
    QExtendedStandardItem *itemSmoothing;

    QExtendedStandardItem *itemTetWild;
    QExtendedStandardItem *itemMeshDataSources;
    QExtendedStandardItem *itemIdentifier;

    QExtendedStandardItem *itemConvection;
    QExtendedStandardItem *itemCouplingTime;
    QExtendedStandardItem *itemEmitter;
    QExtendedStandardItem *itemPosition;

    Property::loadDefinition old_componentX;
    Property::loadDefinition old_componentY;
    Property::loadDefinition old_componentZ;
    Property::loadDefinition old_magnitude;
    QVector<double> old_direction;

private:

    //! create separators
    void createSeparators();

private slots:

    //! create the node model
    void createNodeModel();

public:

    Property::loadDefinition getOldXLoadDefinition() {return old_componentX;}
    Property::loadDefinition getOldYLoadDefinition() {return old_componentY;}
    Property::loadDefinition getOldZLoadDefinition() {return old_componentZ;}
    Property::loadDefinition getOldMagnitude() {return old_magnitude;}
    QVector<double> getOldDirection() {return old_direction;}

    void updateOldLoadDefinition(int componentNumber, Property::loadDefinition loadDef)
    {
        switch(componentNumber)
        {
        case 1: old_componentX = loadDef; break;
        case 2: old_componentY = loadDef; break;
        case 3: old_componentZ = loadDef; break;
        }
    }
    void updateOldMagnitude(Property::loadDefinition oldMagnitude) { old_magnitude = oldMagnitude; }
    void updateOldDirection(const QVector<double> &oldDirection) {old_direction = oldDirection;}

    //! get family
    nodeType getFamily();

    //! set name
    inline void setName(QString &name) { myNodeName = name; }

    //! get name
    inline QString getName() const { return myNodeName; }

    //! set type
    inline void setType(nodeType aType) { myNodeType = aType; }

    //! get type
    inline nodeType getType() const { return myNodeType ;}

    //! get the node model
    inline QStandardItemModel* getModel() const { return myNodeModel; }

    //! get the vector of properties
    inline QVector<Property> getProperties() const { return myProps; }

    //! get items
    QVector<QExtendedStandardItem*> getPropertyItems() const;

    //! get one item - retrieve an item by the name of the property
    QExtendedStandardItem* getPropertyItem(const QString &propertyName, QObject *caller=0) const;

    //! add the property
    void addProperty(Property theProp, int position = -1);

    //! remove a property
    bool removeProperty(const QString &keyPropertyName);

    //! replace a property
    bool replaceProperty(const QString &keyPropertyName, Property theProp);

    //! add tabular data
    void createTabularData(const QVector<load> &vecLoad, bool addFirstRow = true);

    //! set empty tabular data
    //void initTabularData();

    //! ------------------------
    //! get value of a property
    //! ------------------------
    template<class T>
    inline T getPropertyValue(const QString &propertyName)
    {
        return this->getPropertyItem(propertyName)->data(Qt::UserRole).value<Property>().getData().value<T>();
    }

    //! is analysis root
    bool isAnalysisRoot();

    //! is analysis settings
    bool isAnalysisSettings();

    //! is simulaion set up node
    bool isSimulationSetUpNode();

    //! is Solution node
    bool isSolution();

    //! is analysis result
    bool isAnalysisResult();

    //! is solution information
    bool isSolutionInformation();
};

Q_DECLARE_METATYPE(SimulationNodeClass)

#endif // SIMULATIONNODECLASS_H
