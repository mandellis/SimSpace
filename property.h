#ifndef PROPERTY_H
#define PROPERTY_H

//! ---
//! Qt
//! ---
#include <QString>
#include <QVariant>
#include <QMetaType>
#include <QObject>          //! this is needed for the Q_GADGET macro
#include <QMetaEnum>        //! enum variables handling
#include <QDataStream>
#include <QStandardItem>

//! ----------------
//! custom includes
//! ----------------
#include "listofshape.h"

class Property
{
    Q_GADGET

public:

    enum PropertyGroup
    {
        PropertyGroup_NULL,
        PropertyGroup_Display,
        PropertyGroup_Lighting,
        PropertyGroup_AutoDetection,
        PropertyGroup_Transparency,
        PropertyGroup_Definition,
        PropertyGroup_BoundingBox,
        PropertyGroup_Properties,
        PropertyGroup_Statistics,
        PropertyGroup_Material,
        PropertyGroup_GraphicProperties,
        PropertyGroup_Scope,
        PropertyGroup_StepControls,
        PropertyGroup_SolverControls,
        PropertyGroup_Advanced,
        PropertyGroup_Origin,
        PropertyGroup_DirectionalData,
        PropertyGroup_Transformations,
        PropertyGroup_Information,
        PropertyGroup_Validity,
        PropertyGroup_EnvironmentTemperature,
        PropertyGroup_MeshViewOptions,
        PropertyGroup_Sizing,
        PropertyGroup_GraphicObjects,
        PropertyGroup_OutputSettings,
        PropertyGroup_Background,
        PropertyGroup_SolutionInfo,
        PropertyGroup_ConvergenceCriteria,
        PropertyGroup_TimeIncrementation,
        PropertyGroup_CutBack,
        PropertyGroup_LineSearch,
        PropertyGroup_Hidden,
        PropertyGroup_ColorBox,
        PropertyGroup_Defeaturing,
        PropertyGroup_Method,
        PropertyGroup_Smoothing,
        PropertyGroup_TetWild,
        PropertyGroup_MeshDataSources,
        PropertyGroup_Identifier,
        PropertyGroup_CouplingTime,
        PropertyGroup_Emitter,
        PropertyGroup_Position
    };
    Q_ENUM(PropertyGroup)

    enum meshEngine2D
    {
        meshEngine2D_ProgramControlled,
        meshEngine2D_Netgen,
        meshEngine2D_OCC_ExpressMesh,
        meshEngine2D_OCC_STL,
        meshEngine2D_Netgen_STL,
        meshEngine2D_NULL   // no surface mesher when using Netgen with a PLC generated using a surface tessellator
    } me2d;
    Q_ENUM(meshEngine2D)

    enum meshEngine3D
    {
        meshEngine3D_Netgen,
        meshEngine3D_Tetgen,
        meshEngine3D_Netgen_STL,
        meshEngine3D_Tetgen_BR,
        meshEngine3D_TetWild
    } me3d;
    Q_ENUM(meshEngine3D)

    enum meshOrder
    {
        meshOrder_First,
        meshOrder_Second
    } mor;
    Q_ENUM(meshOrder)

    enum meshType_Surface
    {
        meshType_Surface_AllTrig,
        meshType_Surface_QuadDominant
    } mts;
    Q_ENUM(meshType_Surface)

    enum meshType_Volume
    {
        meshType_Volume_AllTet,
        meshType_Volume_HexaDominant
    } mtv;
    Q_ENUM(meshType_Volume)

    enum ScopingMethod
    {
        ScopingMethod_GeometrySelection,
        ScopingMethod_NamedSelection,
        ScopingMethod_RemotePoint,
        ScopingMethod_Automatic,
        ScopingMethod_MeshSelection
    } sm;
    Q_ENUM(ScopingMethod)

    enum SuppressionStatus
    {
        SuppressionStatus_Active,
        SuppressionStatus_Suppressed
    } ss;
    Q_ENUM(SuppressionStatus)

    enum autoTimeStepping
    {
        autoTimeStepping_ProgramControlled,
        autoTimeStepping_ON,
        autoTimeStepping_OFF
    } ats;
    Q_ENUM(autoTimeStepping)

    enum integrationScheme
    {
        integrationScheme_full,
        integrationScheme_reduced
    } is;
    Q_ENUM(integrationScheme)

    enum boltStatusDefinedBy
    {
        boltStatusDefinedBy_load,
        boltStatusDefinedBy_adjustment,
        boltStatusDefinedBy_open,
        boltStatusDefinedBy_lock
    } bs;
    Q_ENUM(boltStatusDefinedBy)

    enum solverType
    {
        solverType_programControlled,
        solverType_direct,
        solverType_iterative
    } st;
    Q_ENUM(solverType)

    enum elementControl
    {
        elementControl_programControlled,
        elementControl_manual
    } ec;
    Q_ENUM(elementControl)

    enum DOFfreedom
    {
        DOFfreedom_fixed,
        DOFfreedom_free
    } df;
    Q_ENUM(DOFfreedom)

    enum typeOfValue
    {
        typeOfValue_constant,
        typeOfValue_tabular,
        typeOfValue_function,
        typeOfValue_free,
        typeOfValue_undefined
    } tv;
    Q_ENUM(typeOfValue)

    enum defineBy
    {
        defineBy_components,
        defineBy_vector,
        defineBy_normal,
        defineBy_geometrySelection,         //! for coordinate systems
        defineBy_globalCoordinates,         //! for coordinate systems
    } db;
    Q_ENUM(defineBy)

    enum solverEngine
    {
        solverEngine_CCX,
        solverEngine_Ansys,
        solverEngine_OpenFoam
    } se;
    Q_ENUM(solverEngine)

    enum contactType
    {
        contactType_bonded,
        contactType_frictional,
        contactType_frictionless,
        contactType_tied
    } ct;
    Q_ENUM(contactType)

    enum contactBehavior
    {
        contactBehavior_asymmetric,
        contactBehavior_symmetric,
    } cb;
    Q_ENUM(contactBehavior)

    enum contactFormulation
    {
        contactFormulation_penalty,
        contactFormulation_MPC,
    } cf;
    Q_ENUM(contactFormulation)

    enum overpressureFunction
    {
        overpressureFunction_linear,
        overpressureFunction_exponential,
        overpressureFunction_tied,
    } of;
    Q_ENUM(overpressureFunction)

    enum modelChangeActivationStatus
    {
        modelChangeActivationStatus_Remove = -1,
        modelChangeActivationStatus_Inactive = 0,
        modelChangeActivationStatus_Add = 1
    } mcas;
    Q_ENUM(modelChangeActivationStatus)

    enum loadType
    {
        //! force and remote force
        loadType_forceX,
        loadType_forceY,
        loadType_forceZ,
        loadType_forceMagnitude,

        //! rotational velocity
        loadType_rotationalVelocityX,
        loadType_rotationalVelocityY,
        loadType_rotationalVelocityZ,
        loadType_rotationalVelocityMagnitude,

        //! moment
        loadType_momentX,
        loadType_momentY,
        loadType_momentZ,
        loadType_momentMagnitude,

        //! acceleration
        loadType_accelerationX,
        loadType_accelerationY,
        loadType_accelerationZ,
        loadType_accelerationMagnitude,

        //! displacement and remote displacement
        loadType_displacementX,
        loadType_displacementY,
        loadType_displacementZ,
        loadType_displacementMagnitude,
        loadType_displacementFree,          // maybe can be removed

        //! remote rotation
        loadType_remoteRotationX,
        loadType_remoteRotationY,
        loadType_remoteRotationZ,
        loadType_remoteRotationMagnitude,
        loadType_remoteRotationFree,        // maybe can be removed

        loadType_pressureMagnitude,

        //! --------------------
        //! thermal environment
        //! --------------------
        loadType_temperatureMagnitude,
        loadType_thermalFlowMagnitude,
        loadType_thermalFluxMagnitude,
        loadType_thermalPowerMagnitude,
        loadType_thermalConvectionFilmCoefficientMagnitude,
        loadType_thermalConvectionReferenceTemperatureMagnitude,

        loadType_thermalConditionTemperature,

        //! -----------------------
        //! in "Analysis settings"
        //! -----------------------
        loadType_stepNumber,
        loadType_stepEndTime,
        loadType_solverType,
        loadType_fieldParameters,
        loadType_fluxConvergence,
        loadType_solutionConvergence,
        loadType_timeIncrementation,
        loadType_timeIncrementationParameters,
        loadType_cutBack,
        loadType_cutBackParameters,
        loadType_lineSearch,
        loadType_lineSearchParameters,
        loadType_outputSettings,
        loadType_storeResultsAt,

        loadType_NLgeom,

        loadType_autoTimeStepping,
        loadType_time,
        loadType_none,

        loadType_boltStatusDefinedBy,   //! data in the first column in table for bolt
        loadType_boltForce,             //! data in the second column in table for bolt
        loadType_boltAdjustment,        //! data in the third column in table for bolt

        loadType_modelChange,

        //! --------------
        //! analysis type
        //! --------------
        loadType_analysisType,
        loadType_timeIntegration
    } lt;
    Q_ENUM(loadType)

    // -------------
    // experimental
    // -------------
    //enum boltDefineBy
    //{
    //    boltDefineBy_load,                     //! for bolt pretension
    //    boltDefinedBy_adjustment,              //! for bolt pretension
    //    boltDefineBy_open,                     //! for bolt pretension
    //    boltDefineBy_lock                      //! for bolt pretension
    //} bdf;
    //Q_ENUM(boltDefineBy)

    enum loadDefinition
    {
        loadDefinition_constant,        //! - 0
        loadDefinition_tabularData,     //! - 1
        loadDefinition_function,        //! - 2
        loadDefinition_free             //! - 3
    } ld;
    Q_ENUM(loadDefinition)

    enum typeOfTransformation
    {
        typeOfTransformation_offsetX,
        typeOfTransformation_offsetY,
        typeOfTransformation_offsetZ,
        typeOfTransformation_rotationX,
        typeOfTransformation_rotationY,
        typeOfTransformation_rotationZ,
        typeOfTransformation_none
    } tt;
    Q_ENUM(typeOfTransformation)

    enum solutionInformation
    {
        solutionInformation_solveRequired,
        solutionInformation_done,
        solutionInformation_running,
        solutionInformation_interrupted,
        solutionInformation_failed
    } sinf;
    Q_ENUM(solutionInformation)

    enum meshMethod
    {
        meshMethod_automatic,       //0
        meshMethod_NetgenTetgen,    //1
        meshMethod_NetgenNetgen,    //2
        meshMethod_EMeshTetgen,     //3
        meshMethod_EMeshNetgen      //4
    } mm;
    Q_ENUM(meshMethod)

    enum analysisType
    {
        analysisType_structural,
        analysisType_thermal,
        analysisType_modal,
        analysisType_frequencyResponse,
        analysisType_uncoupledTemperatureDisplacement,
        analysisType_coupledTemperatureDisplacement
    } ant;
    Q_ENUM(analysisType)

    enum timeIntegration
    {
        timeIntegration_steadyState,
        timeIntegration_transient
    } ti;
    Q_ENUM(timeIntegration)

public:

    //! constructor
    Property();

    //! constructor I
    Property(QString name, QVariant data, PropertyGroup group);

    //! copy constructor
    Property(const Property &other);

    //! operator ==
    bool operator==(const Property &rhs) const
    {
        if(propertyName != rhs.propertyName) return false;
        if(propertyGroup != rhs.propertyGroup) return false;
        return true;
    }

private:

    void init(const QVariant &data);

private:

    QString propertyName;
    QVariant propertyValue;
    PropertyGroup propertyGroup;

public:

    //! set/get name
    inline void setName(const QString &name) { propertyName = name; }
    inline QString getName() const { return propertyName; }

    //! set/get data
    inline void setData(const QVariant &data);
    inline QVariant getData() const { return propertyValue; }

    //! set/get group
    inline void setGroup(const PropertyGroup &aGroup) { propertyGroup=aGroup; }
    inline PropertyGroup getGroup() const { return propertyGroup; }

    const QString getPropertyGroup() const;
    void setPropertyGroup(const QString &thePropertyGroup);

    const Property::PropertyGroup getPropertyGroupEnum() const;

    const QString getPropertyValue(const QString &propName) const;
    void setPropertyValue(const QString &enumName, const QString &thePropertyValue);

public:

    //static QMultiMap<QString, QString> propertyMap();
    static QMap<QString, QString> propertyMap();

    static void readProperty(std::ifstream &in, Property &prop);
    static void readVoid(ifstream &in, QStandardItem *item);

    static void writeVoid(ofstream& outFile, void *p);
    static void writeProperty(ofstream& out, const Property &prop);
    static void writeProperty1(ofstream& out, const Property &prop);
};


Q_DECLARE_METATYPE(Property)

#endif // PROPERTY_H
