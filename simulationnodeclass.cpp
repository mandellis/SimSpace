//! ----------------
//! custom includes
//! ----------------
#include "SimulationNodeClass.h"
#include "qextendedstandarditem.h"
#include "load.h"
#include "property.h"
#include "customtablemodel.h"
#include "detailviewer.h"
#include "simulationmanager.h"

//! ---
//! Qt
//! ---
#include <QStandardItemModel>
#include <QThread>

//! ----
//! C++
//! ----
#include <chrono>

//!--------------------------------------------------------------//
//! function: constructor I /'ìtrrrrèpì985                                      //
//! details:  tahnks to Gilda for comments above                 //
//!--------------------------------------------------------------//
SimulationNodeClass::SimulationNodeClass(QObject *parent):QObject(parent)
{
    //! --------------
    //! the node name
    //! --------------
    myNodeName = "";

    //! ---------
    //! the type
    //! ---------
    myNodeType = nodeType_NULL;

    //! ------------------------------
    //! an empty vector of properties
    //! ------------------------------
    myProps.empty();

    //! ---------------
    //! the node model
    //! ---------------
    this->createNodeModel();

    //! -----------------
    //! add the time tag
    //! -----------------
    //this->addTimeTag();
}

//!--------------------------
//! function: constructor II
//! details:
//! -------------------------
SimulationNodeClass::SimulationNodeClass(const QString &aName, const nodeType &aType, const QVector<Property> &props, QObject *parent):
    QObject(parent),
    myProps(props),
    myNodeType(aType),
    myNodeName(aName)
{
    //cout<<"SimulationNodeClass::SimulationNodeClass()->____constructor II called____"<<endl;

    //! ---------------
    //! the node model
    //! ---------------
    this->createNodeModel();

    //! -----------------
    //! add the time tag
    //! -----------------
    //this->addTimeTag();

    //! ---------------------------------
    //! initialize some load definitions
    //! ---------------------------------
    if(aType==nodeType_structuralAnalysisBoundaryCondition_Displacement ||
            aType ==nodeType_structuralAnalysisBoundaryCondition_RemoteDisplacement ||
            aType ==nodeType_structuralAnalysisBoundaryCondition_RemoteRotation)
    {
        old_componentX = Property::loadDefinition_free;
        old_componentY = Property::loadDefinition_free;
        old_componentZ = Property::loadDefinition_free;
        old_magnitude = Property::loadDefinition_constant;
    }
    else
    {
        old_componentX = Property::loadDefinition_constant;
        old_componentY = Property::loadDefinition_constant;
        old_componentZ = Property::loadDefinition_constant;
        old_magnitude = Property::loadDefinition_constant;
    }

    //! -------------------------
    //! initialize a "Direction"
    //! -------------------------
    old_direction.push_back(0); old_direction.push_back(0); old_direction.push_back(0);
    old_direction.push_back(1); old_direction.push_back(0); old_direction.push_back(0);
}

//!---------------------------
//! function: createNodeModel
//! details:
//! --------------------------
void SimulationNodeClass::createNodeModel()
{
    //cout<<"SimulationNodeClass::createNodeModel()->____function called____"<<endl;
    myNodeModel = new QStandardItemModel(this);
    myNodeModel->setColumnCount(2);
    myNodeModel->setHorizontalHeaderLabels(QStringList()<<"Property"<<"Value");
    myNodeRootItem = static_cast<QExtendedStandardItem*>(myNodeModel->invisibleRootItem());

    QVariant data;

    //! ----------------------
    //! create the separators
    //! ----------------------
    this->createSeparators();

    //! ------------------------------------------------------
    //! create and attach the items containing the properties
    //! ------------------------------------------------------
    for(int i=0; i<myProps.length(); i++) this->addProperty(myProps.at(i));

    //! ------------------------------------------------------------------
    //! additional items needing a special treatment
    //! the values of all the properties created here are used as default
    //! ------------------------------------------------------------------
    QExtendedStandardItem *item = this->getPropertyItem("Define by");
    if(item!=Q_NULLPTR)
    {
        //! create node-specific properties
        Property::defineBy theDefineBy = item->data(Qt::UserRole).value<Property>().getData().value<Property::defineBy>();
        if(theDefineBy==Property::defineBy_vector)
        {
            //! --------------------------------------------------------------------------------
            //! check first if the node does not have "Direction" and "Magnitude" defined
            //! this must be done here, since when reloading from disk a duplicated "Direction"
            //! and a duplicated "Magnitude" would be created
            //! --------------------------------------------------------------------------------
            if(this->getPropertyItem("Direction")==Q_NULLPTR && myNodeType!=SimulationNodeClass::nodeType_structuralAnalysisBoltPretension)
            {
                //! ----------------------------------------------------------
                //! The direction is identified by the three cosines
                //! Also an origin is set, in order to place a marker
                //! (an arrow) correctly positioned on the model. Here vec[6]
                //! Xorigin = vec[0], Yorigin = vec[1], Zorigin = vec[2]
                //! cos_alfa = vec[3], cos_beta = vec[4], cos_gamma = vec[5]
                //! ----------------------------------------------------------

                //!cout<<"SimulationNodeClass::createNodeModel()->____adding 'Direction' for a 'Define by vector'____"<<endl;
                QVector<double> vec;
                //! the definition of the origin point
                vec.append(0.0); vec.append(0.0); vec.append(0.0);
                //! the definition of the direction => global X axis direction by default
                vec.append(1.0); vec.append(0.0); vec.append(0.0);

                //QVariant data;
                data.setValue(vec);
                Property property_Direction("Direction",data,Property::PropertyGroup_Definition);
                this->addProperty(property_Direction);
            }

            if(this->getPropertyItem("Magnitude")==Q_NULLPTR && myNodeType!=SimulationNodeClass::nodeType_structuralAnalysisBoltPretension)
            {
                //!cout<<"SimulationNodeClass::createNodeModel()->____adding 'Magnitude' for a 'Define by vector'____"<<endl;
                Property::loadDefinition theLoadDefinition;
                if(this->getType()==nodeType_structuralAnalysisBoundaryCondition_Displacement)
                {
                    theLoadDefinition = Property::loadDefinition_free;
                }
                else
                {
                    //theLoadDefinition = loadDefinition_tabularData;
                    theLoadDefinition = Property::loadDefinition_constant;
                }
                //QVariant data;
                data.setValue(theLoadDefinition);
                Property property_magnitude("Magnitude",data,Property::PropertyGroup_Definition);
                this->addProperty(property_magnitude);
            }
        }
        else if (theDefineBy==Property::defineBy_components)
        {
            cout<<"____a load defined by components has been created____"<<endl;
        }
        else if(theDefineBy==Property::defineBy_normal)
        {
            cout<<"____a pressure has been created: only \"Normal to\" is shown____"<<endl;
        }
    }
}

//!---------------------------
//! function: add a property
//! details:  helper function
//! --------------------------
void SimulationNodeClass::addProperty(Property theProp, int position)
{
    QVariant data;

    //! this is the item "label" (first column)
    QExtendedStandardItem *itemLabel = new QExtendedStandardItem();
    QVariant keyName;
    keyName.setValue(theProp.getName());

    itemLabel->setData(keyName,Qt::DisplayRole);
    data.setValue(keyName);
    itemLabel->setData(data,Qt::UserRole);

    //! this is the item containing the property (second column)
    data.setValue(theProp);
    QExtendedStandardItem *anItem = new QExtendedStandardItem();
    anItem->setData(data, Qt::UserRole);

    //!anItem->setSelectable(true);

    //! a row with two columns is appended
    QList<QStandardItem*> list;
    list<<itemLabel<<anItem;

    Property::PropertyGroup thePropertyGroup = theProp.getGroup();

    switch(thePropertyGroup)
    {
    case Property::PropertyGroup_BoundingBox:
        if(position == -1)itemBoundingBox->appendRow(list);
        else itemBoundingBox->insertRow(position,list);
        break;
    case Property::PropertyGroup_Definition:
        if(position == -1)itemDefinition->appendRow(list);
        else itemDefinition->insertRow(position,list);
        break;
    case Property::PropertyGroup_GraphicProperties:
        if(position == -1)itemGraphicProperties->appendRow(list);
        else itemGraphicProperties->insertRow(position,list);
        break;
    case Property::PropertyGroup_Material:
        if(position == -1)itemMaterial->appendRow(list);
        else itemMaterial->insertRow(position,list);
        break;
    case Property::PropertyGroup_Properties:
        if(position == -1)itemProperties->appendRow(list);
        else itemProperties->insertRow(position,list);
        break;
    case Property::PropertyGroup_Scope:
        if(position == -1)itemScope->appendRow(list);
        else itemScope->insertRow(position,list);
        break;
    case Property::PropertyGroup_Statistics:
        if(position == -1)itemStatistics->appendRow(list);
        else itemStatistics->insertRow(position,list);
        break;
    case Property::PropertyGroup_StepControls:
        if(position == -1)itemStepControls->appendRow(list);
        else itemStepControls->insertRow(position,list);
        break;
    case Property::PropertyGroup_SolverControls:
        if(position == -1)itemSolverControl->appendRow(list);
        else itemSolverControl->insertRow(position,list);
        break;
    case Property::PropertyGroup_Lighting:
        if(position ==-1)itemLighting->appendRow(list);
        else itemLighting->insertRow(position,list);
        break;
    case Property::PropertyGroup_AutoDetection:
        if(position ==-1)itemAutoDetection->appendRow(list);
        else itemAutoDetection->insertRow(position,list);
        break;
    case Property::PropertyGroup_Transparency:
        if(position == -1)itemTransparency->appendRow(list);
        else itemTransparency->insertRow(position,list);
        break;
    case Property::PropertyGroup_Display:
        if(position == -1)itemDisplay->appendRow(list);
        else itemDisplay->insertRow(position,list);
        break;
    case Property::PropertyGroup_Advanced:
        if(position == -1)itemAdvanced->appendRow(list);
        else itemAdvanced->insertRow(position,list);
        break;
    case Property::PropertyGroup_Origin:
        if(position == -1)itemOrigin->appendRow(list);
        else itemOrigin->insertRow(position,list);
        break;
    case Property::PropertyGroup_DirectionalData:
        if(position == -1)itemDirectionalData->appendRow(list);
        else itemDirectionalData->insertRow(position,list);
        break;
    case Property::PropertyGroup_Transformations:
        if(position == -1)itemTransformations->appendRow(list);
        else itemTransformations->insertRow(position,list);
        break;
    case Property::PropertyGroup_Information:
        if(position == -1)itemInformation->appendRow(list);
        else itemInformation->insertRow(position,list);
        break;
    case Property::PropertyGroup_Validity:
        if(position == -1)itemValidity->appendRow(list);
        else itemValidity->insertRow(position,list);
        break;
    case Property::PropertyGroup_EnvironmentTemperature:
        if(position == -1)itemEnvironmentTemperature->appendRow(list);
        else itemEnvironmentTemperature->insertRow(position,list);
        break;
    case Property::PropertyGroup_MeshViewOptions:
        if(position == -1)itemMeshViewOptions->appendRow(list);
        else itemMeshViewOptions->insertRow(position,list);
        break;
    case Property::PropertyGroup_Sizing:
        if(position == -1)itemSizing->appendRow(list);
        else itemSizing->insertRow(position,list);
        break;
    case Property::PropertyGroup_GraphicObjects:
        if(position == -1)itemGraphicObject->appendRow(list);
        else itemGraphicObject->insertRow(position,list);
        break;
    case Property::PropertyGroup_OutputSettings:
        if(position == -1)itemOutputSettings->appendRow(list);
        else itemOutputSettings->insertRow(position,list);
        break;
    case Property::PropertyGroup_Background:
        if(position == -1) itemBackground->appendRow(list);
        else itemBackground->insertRow(position,list);
        break;
    case Property::PropertyGroup_SolutionInfo:
        if(position == -1) itemSolutionInformation->appendRow(list);
        else itemSolutionInformation->insertRow(position,list);
        break;
    case Property::PropertyGroup_ConvergenceCriteria:
        if(position == -1) itemConvergenceCriteria->appendRow(list);
        else itemConvergenceCriteria->insertRow(position,list);
        break;
    case Property::PropertyGroup_TimeIncrementation:
        if(position == -1) itemTimeIncrementation->appendRow(list);
        else itemTimeIncrementation->insertRow(position,list);
        break;
    case Property::PropertyGroup_CutBack:
        if(position == -1) itemCutBack->appendRow(list);
        else itemCutBack->insertRow(position,list);
        break;
    case Property::PropertyGroup_LineSearch:
        if(position == -1) itemLineSearch->appendRow(list);
        else itemLineSearch->insertRow(position,list);
        break;
    case Property::PropertyGroup_Hidden:
        if(position == -1) itemHidden->appendRow(list);
        else itemHidden->insertRow(position,list);
        break;
    case Property::PropertyGroup_ColorBox:
        if(position == -1) itemColorBox->appendRow(list);
        else itemColorBox->insertRow(position,list);
        break;
    case Property::PropertyGroup_Defeaturing:
        if(position == -1) itemDefeaturing->appendRow(list);
        else itemDefeaturing->insertRow(position,list);
        break;
    case Property::PropertyGroup_Method:
        if(position == -1) itemMethod->appendRow(list);
        else itemMethod->insertRow(position,list);
        break;

    case Property::PropertyGroup_Smoothing:
        if(position == -1) itemSmoothing->appendRow(list);
        else itemSmoothing->insertRow(position,list);
        break;

    case Property::PropertyGroup_TetWild:
        if(position == -1) itemTetWild->appendRow(list);
        else itemTetWild->insertRow(position,list);
        break;
    case Property::PropertyGroup_MeshDataSources:
        if(position == -1) itemMeshDataSources->appendRow(list);
        else itemMeshDataSources->insertRow(position,list);
        break;
    case Property::PropertyGroup_Identifier:
        if(position == -1) itemIdentifier->appendRow(list);
        else itemIdentifier->insertRow(position,list);
        break;
    case Property::PropertyGroup_CouplingTime:
        if(position == -1) itemCouplingTime->appendRow(list);
        else itemCouplingTime->insertRow(position,list);
        break;
    case Property::PropertyGroup_Emitter:
        if(position == -1) itemEmitter->appendRow(list);
        else itemEmitter->insertRow(position,list);
        break;
    case Property::PropertyGroup_Position:
        if(position == -1) itemPosition->appendRow(list);
        else itemPosition->insertRow(position,list);
        break;
    }
}

//!---------------------------------------------------------
//! function: createTabularData
//! details:  build the tabular data from a vector of loads
//! --------------------------------------------------------
void SimulationNodeClass::createTabularData(const QVector<load> &vecLoad, bool addFirstRow)
{
    cout<<"SimulationNodeClass::createTabularData()->____function called____"<<endl;
    theTabularDataModel = new CustomTableModel(vecLoad, addFirstRow, this);
    cout<<"SimulationNodeClass::createTabularData()->____exiting____"<<endl;
}

//!-----------------------------------------
//! function: getTabularDataModel
//! details:  return the tabular data model
//! ----------------------------------------
CustomTableModel *SimulationNodeClass::getTabularDataModel()
{
    return theTabularDataModel;
}

//!----------------------------
//! function: copy constructor
//! details:
//!----------------------------
SimulationNodeClass::SimulationNodeClass(const SimulationNodeClass &other)
{
    myNodeName = other.myNodeName;
    myNodeType = other.myNodeType;

    QVector<Property> props = other.getProperties();
    for(QVector<Property>::iterator it = props.begin(); it!=props.end(); ++it)
    {
        const Property &prop = *it;
        myProps.push_back(prop);
    }
    this->createNodeModel();
}

//!----------------------
//! function: destructor
//! details:
//!----------------------
SimulationNodeClass::~SimulationNodeClass()
{
    ;
}

//!---------------------
//! function: getFamily
//! details:
//! --------------------
SimulationNodeClass::nodeType SimulationNodeClass::getFamily()
{
    nodeType RV;
    switch(this->myNodeType)
    {
    case nodeType_root:
        RV = nodeType_root;
        break;

    case nodeType_coordinateSystem:
    case nodeType_coordinateSystems:
    case nodeType_coordinateSystem_global:
        RV = nodeType_coordinateSystems;
        break;

    case nodeType_thermalAnalysis:
    case nodeType_thermalAnalysisSettings:
    case nodeType_thermalAnalysisTemperature:
    case nodeType_thermalAnalysisConvection:
    case nodeType_thermalAnalysisThermalFlux:
    case nodeType_thermalAnalysisThermalFlow:
    case nodeType_thermalAnalysisThermalPower:
    case nodeType_thermalAnalysisRadiation:
    case nodeType_thermalAnalysisAdiabaticWall:
        RV = nodeType_thermalAnalysis;
        break;


    case nodeType_thermalAnalysisSolution:
    case nodeType_thermalAnalysisSolutionInformation:
    case nodeType_solutionThermalFlux:
    case nodeType_solutionThermalTemperature:
        RV = nodeType_thermalAnalysisSolution;
        break;

    case nodeType_meshBodyMeshControl:
    case nodeType_meshBodyMeshMethod:
    case nodeType_meshEdgeSize:
    case nodeType_meshFaceSize:
    case nodeType_meshVertexSize:
    case nodeType_meshGrading:
    case nodeType_meshMinElementSize:
    case nodeType_meshMaxElementSize:
    case nodeType_meshGenericSizing:
    case nodeType_meshControl:
    case nodeType_meshPrismaticLayer:
    case nodeType_meshMeshType:
    case nodeType_meshMeshMetric:
        RV = nodeType_meshControl;
        break;

    case nodeType_structuralAnalysis:
    case nodeType_structuralAnalysisSettings:
    case nodeType_structuralAnalysisBoundaryContidion_FixedSupport:
    case nodeType_structuralAnalysisBoundaryCondition_CylindricalSupport:
    case nodeType_structuralAnalysisBoundaryCondition_FrictionlessSupport:
    case nodeType_structuralAnalysisBoundaryCondition_CompressionOnlySupport:
    case nodeType_structuralAnalysisBoundaryCondition_Acceleration:
    case nodeType_structuralAnalysisBoundaryCondition_RotationalVelocity:
    case nodeType_structuralAnalysisBoundaryCondition_RemoteForce:
    case nodeType_structuralAnalysisBoundaryCondition_Force:
    case nodeType_structuralAnalysisBoundaryCondition_Moment:
    case nodeType_structuralAnalysisBoundaryCondition_Pressure:
    case nodeType_structuralAnalysisBoundaryCondition_Displacement:
    case nodeType_structuralAnalysisBoundaryCondition_RemoteDisplacement:
    case nodeType_structuralAnalysisBoundaryCondition_RemoteRotation:
    case nodeType_structuralAnalysisThermalCondition:
    case nodeType_structuralAnalysisBoltPretension:
    case nodeType_structuralAnalysisBoundaryCondition_ImportedTemperatureDistribution:
    case nodeType_mapper:
    case nodeType_importedBodyScalar:
    case nodeType_OpenFoamScalarData:
    case nodeType_modelChange:
#ifdef COSTAMP_VERSION
    case nodeType_timeStepBuilder:
#endif
        RV = nodeType_structuralAnalysis;
        break;

    case nodeType_geometryBody:
    case nodeType_geometryPart:
    case nodeType_geometry:
    case nodeType_pointMass:
    case nodeType_repairTool:
        RV = nodeType_geometry;
        break;

    case nodeType_connection:
    case nodeType_connectionPair:
    case nodeType_connectionGroup:
        RV = nodeType_connection;
        break;

    case nodeType_namedSelection:
    case nodeType_namedSelectionGeometry:
    case nodeType_namedSelectionElement:
    case nodeType_namedSelectionNodal:
        RV = nodeType_namedSelection;
        break;

    case nodeType_StructuralAnalysisSolution:
    case nodeType_StructuralAnalysisSolutionInformation:
    case nodeType_solutionStructuralNodalDisplacement:
    case nodeType_solutionStructuralStress:
    case nodeType_solutionStructuralTotalStrain:
    case nodeType_solutionStructuralThermalStrain:
    case nodeType_solutionStructuralMechanicalStrain:
    case nodeType_solutionStructuralFatigueTool:
    case nodeType_solutionStructuralTemperature:
    case nodeType_solutionStructuralEquivalentPlasticStrain:
    case nodeType_solutionStructuralNodalForces:
    case nodeType_solutionStructuralContact:
    case nodeType_solutionStructuralGamma:
    case nodeType_solutionStructuralReactionForce:
        RV = nodeType_StructuralAnalysisSolution;
        break;

    case nodeType_combinedAnalysis:
    case nodeType_combinedAnalysisSettings:
        RV = nodeType_combinedAnalysis;
        break;

   case nodeType_particlesInFieldsAnalysis:
   case nodeType_particlesInFieldsAnalysisSettings:
   case nodeType_electrostaticPotential:
        RV = nodeType_particlesInFieldsAnalysis;
        break;

    case nodeType_particlesInFieldsSolution:
    case nodeType_particlesInFieldsSolutionInformation:
        RV = nodeType_particlesInFieldsSolution;
        break;

    case nodeType_combinedAnalysisSolution:
    case nodeType_combinedAnalysisSolutionInformation:
        RV = nodeType_combinedAnalysisSolution;
        break;

    case nodeType_remotePointRoot:
    case nodeType_remotePoint:
        RV = nodeType_remotePointRoot;
        break;

    case nodeType_meshMethod:
        RV = nodeType_meshControl;
        break;

    case nodeType_import:
        RV = nodeType_import;
        break;

    case nodeType_postObject:
        RV = nodeType_postObject;
        break;

#ifdef COSTAMP_VERSION
    case nodeType_processParametersClosureForce:
    case nodeType_processParametersPressure:
        RV = nodeType_processParameters;
        break;
#endif

    default:
        break;
    }
    return RV;
}

//!----------------------------
//! function: getPropertyItems
//! details:  helper function
//!----------------------------
QVector<QExtendedStandardItem*> SimulationNodeClass::getPropertyItems() const
{
    QVector<QExtendedStandardItem*> propertyItems;
    for(int k=0;k<myNodeRootItem->rowCount();k++)
    {
        //! the child item is the "separator" item
        QExtendedStandardItem *childItem = static_cast<QExtendedStandardItem*>(myNodeRootItem->child(k,0));
        for(int j=0;j<childItem->rowCount();j++)
        {
            QExtendedStandardItem *propertyItem = static_cast<QExtendedStandardItem*>(childItem->child(j,1));
            propertyItems.append(propertyItem);
        }
    }
    return propertyItems;
}

//!-------------------------------------------------------------
//! function: getPropertyItem
//! details:  helper function - return the tree item containing
//!           the property with key name "propertyName"
//!-------------------------------------------------------------
QExtendedStandardItem* SimulationNodeClass::getPropertyItem(const QString &propertyName, QObject *caller) const
{
    if(caller!=Q_NULLPTR)
        cout<<"SimulationNodeClass::getPropertyItem()->____called in "<<caller->objectName().toStdString()<<"____"<<endl;
    for(int k=0;k<myNodeRootItem->rowCount();k++)
    {
        QExtendedStandardItem *definitionItem = static_cast<QExtendedStandardItem*>(myNodeRootItem->child(k,0));
        for(int j=0;j<definitionItem->rowCount();j++)
        {
            QExtendedStandardItem *propertyItem = static_cast<QExtendedStandardItem*>(definitionItem->child(j,1));
            if(propertyItem->data(Qt::UserRole).value<Property>().getName() == propertyName)
            {
                //!cout<<"SimulationNodeClass::getPropertyItem()->____Property ->"<<propertyName.toStdString()<<"<- found____"<<endl;
                return propertyItem;
            }
        }
    }
    //cout<<"SimulationNodeClass::getPropertyItem()->____Property key name \""<<propertyName.toStdString()<<"\" NOT found____"<<endl;
    return Q_NULLPTR;
}

//!-----------------------------------------------------------
//! function: createSeparators
//! details:  create the graphic separators for the tree view
//!-----------------------------------------------------------
void SimulationNodeClass::createSeparators()
{
    //! -----------------------------------------------------------------
    //! define all the possible separators - contained in Qt::UserRole+1
    //! -----------------------------------------------------------------
    itemDefinition = new QExtendedStandardItem();
    itemDefinition->setData("Definition", Qt::DisplayRole);
    itemBoundingBox = new QExtendedStandardItem();
    itemBoundingBox->setData("Bounding box", Qt::DisplayRole);
    itemProperties = new QExtendedStandardItem();
    itemProperties->setData("Properties", Qt::DisplayRole);
    itemStatistics = new QExtendedStandardItem();
    itemStatistics->setData("Statistics", Qt::DisplayRole);
    itemMaterial = new QExtendedStandardItem();
    itemMaterial->setData("Material", Qt::DisplayRole);
    itemGraphicProperties = new QExtendedStandardItem();
    itemGraphicProperties->setData("Graphic properties", Qt::DisplayRole);
    itemScope = new QExtendedStandardItem();
    itemScope->setData("Scope", Qt::DisplayRole);
    itemStepControls = new QExtendedStandardItem();
    itemStepControls->setData("Step controls", Qt::DisplayRole);
    itemSolverControl = new QExtendedStandardItem();
    itemSolverControl->setData("Solver controls", Qt::DisplayRole);
    itemLighting = new QExtendedStandardItem();
    itemLighting->setData("Lighting", Qt::DisplayRole);
    itemAutoDetection = new QExtendedStandardItem();
    itemAutoDetection->setData("Auto detection", Qt::DisplayRole);
    itemTransparency = new QExtendedStandardItem();
    itemTransparency->setData("Transparency", Qt::DisplayRole);
    itemDisplay = new QExtendedStandardItem();
    itemDisplay->setData("Display", Qt::DisplayRole);
    itemAdvanced = new QExtendedStandardItem();
    itemAdvanced->setData("Advanced", Qt::DisplayRole);
    itemOrigin = new QExtendedStandardItem();
    itemOrigin->setData("Origin", Qt::DisplayRole);
    itemDirectionalData = new QExtendedStandardItem();
    itemDirectionalData->setData("Directional data",Qt::DisplayRole);
    itemTransformations = new QExtendedStandardItem();
    itemTransformations->setData("Transformations", Qt::DisplayRole);
    itemInformation = new QExtendedStandardItem();
    itemInformation->setData("Information", Qt::DisplayRole);
    itemValidity = new QExtendedStandardItem();
    itemValidity->setData("Status", Qt::DisplayRole);
    itemEnvironmentTemperature = new QExtendedStandardItem();
    itemEnvironmentTemperature->setData("Environment temperature", Qt::DisplayRole);
    itemMeshViewOptions = new QExtendedStandardItem();
    itemMeshViewOptions->setData("Mesh view options", Qt::DisplayRole);
    itemSizing = new QExtendedStandardItem();
    itemSizing->setData("Sizing",Qt::DisplayRole);
    itemGraphicObject = new QExtendedStandardItem();
    itemGraphicObject->setData("Graphic object",Qt::DisplayRole);
    itemOutputSettings = new QExtendedStandardItem();
    itemOutputSettings->setData("Output settings",Qt::DisplayRole);
    itemBackground = new QExtendedStandardItem();
    itemBackground->setData("Background",Qt::DisplayRole);
    itemSolutionInformation= new QExtendedStandardItem();
    itemSolutionInformation->setData("Solution information",Qt::DisplayRole);
    itemConvergenceCriteria = new QExtendedStandardItem();
    itemConvergenceCriteria->setData("Convergence criteria",Qt::DisplayRole);
    itemTimeIncrementation = new QExtendedStandardItem();
    itemTimeIncrementation->setData("Time incrementation",Qt::DisplayRole);
    itemCutBack = new QExtendedStandardItem();
    itemCutBack->setData("Cutback",Qt::DisplayRole);
    itemLineSearch = new QExtendedStandardItem();
    itemLineSearch->setData("Line search",Qt::DisplayRole);
    itemHidden = new QExtendedStandardItem();
    itemHidden->setData("Hidden",Qt::DisplayRole);
    itemColorBox = new QExtendedStandardItem();
    itemColorBox->setData("Color box",Qt::DisplayRole);
    itemDefeaturing = new QExtendedStandardItem();
    itemDefeaturing->setData("Defeaturing",Qt::DisplayRole);
    itemMethod = new QExtendedStandardItem();
    itemMethod->setData("Method",Qt::DisplayRole);
    itemSmoothing = new QExtendedStandardItem();
    itemSmoothing->setData("Smoothing",Qt::DisplayRole);
    itemTetWild = new QExtendedStandardItem();
    itemTetWild->setData("Experimental mesher parameters",Qt::DisplayRole);
    itemMeshDataSources = new QExtendedStandardItem();
    itemMeshDataSources->setData("Mesh data sources",Qt::DisplayRole);
    itemIdentifier = new QExtendedStandardItem();
    itemIdentifier->setData("Identifier",Qt::DisplayRole);
    itemConvection = new QExtendedStandardItem();
    itemConvection->setData("Convection",Qt::DisplayRole);
    itemCouplingTime = new QExtendedStandardItem();
    itemCouplingTime->setData("Couplig time",Qt::DisplayRole);
    itemEmitter = new QExtendedStandardItem();
    itemEmitter->setData("Emitter",Qt::DisplayRole);
    itemPosition = new QExtendedStandardItem();
    itemPosition->setData("Position",Qt::DisplayRole);

    itemDefinition->setData("separator",Qt::UserRole+1);
    itemBoundingBox->setData("separator",Qt::UserRole+1);
    itemProperties->setData("separator",Qt::UserRole+1);
    itemStatistics->setData("separator",Qt::UserRole+1);
    itemMaterial->setData("separator",Qt::UserRole+1);
    itemGraphicProperties->setData("separator",Qt::UserRole+1);
    itemScope->setData("separator",Qt::UserRole+1);
    itemStepControls->setData("separator",Qt::UserRole+1);
    itemSolverControl->setData("separator",Qt::UserRole+1);
    itemLighting->setData("separator",Qt::UserRole+1);
    itemAutoDetection->setData("separator",Qt::UserRole+1);
    itemTransparency->setData("separator",Qt::UserRole+1);
    itemDisplay->setData("separator",Qt::UserRole+1);
    itemAdvanced->setData("separator",Qt::UserRole+1);
    itemOrigin->setData("separator",Qt::UserRole+1);
    itemDirectionalData->setData("separator",Qt::UserRole+1);
    itemTransformations->setData("separator",Qt::UserRole+1);
    itemInformation->setData("separator",Qt::UserRole+1);
    itemValidity->setData("separator",Qt::UserRole+1);
    itemEnvironmentTemperature->setData("separator",Qt::UserRole+1);
    itemMeshViewOptions->setData("separator",Qt::UserRole+1);
    itemSizing->setData("separator",Qt::UserRole+1);
    itemGraphicObject->setData("separator",Qt::UserRole+1);
    itemOutputSettings->setData("separator",Qt::UserRole+1);
    itemBackground->setData("separator", Qt::UserRole+1);
    itemSolutionInformation->setData("separator",Qt::UserRole+1);
    itemConvergenceCriteria->setData("separator",Qt::UserRole+1);
    itemTimeIncrementation->setData("separator",Qt::UserRole+1);
    itemCutBack->setData("separator",Qt::UserRole+1);
    itemLineSearch->setData("separator",Qt::UserRole+1);
    itemHidden->setData("separator",Qt::UserRole+1);
    itemColorBox->setData("separator",Qt::UserRole+1);
    itemDefeaturing->setData("separator",Qt::UserRole+1);
    itemMethod->setData("separator",Qt::UserRole+1);
    itemSmoothing->setData("separator",Qt::UserRole+1);
    itemTetWild->setData("separator",Qt::UserRole+1);
    itemMeshDataSources->setData("separator",Qt::UserRole+1);
    itemIdentifier->setData("separator",Qt::UserRole+1);
    itemConvection->setData("separator",Qt::UserRole+1);
    itemCouplingTime->setData("separator",Qt::UserRole+1);
    itemEmitter->setData("separator",Qt::UserRole+1);
    itemPosition->setData("separator",Qt::UserRole+1);

    //! ------------------------
    //! generate the separators
    //! ------------------------
    nodeType theNodeType = this->getType();

    switch(theNodeType)
    {
    case nodeType_geometry:
    case nodeType_geometryBody:
        if(theNodeType==nodeType_geometryBody || theNodeType ==nodeType_geometryPart)
        {
            myNodeRootItem->appendRow(itemGraphicProperties);
        }
        myNodeRootItem->appendRow(itemScope);

        myNodeRootItem->appendRow(itemDefinition);
        if(theNodeType==nodeType_geometryBody || theNodeType ==nodeType_geometryPart)
        {
            myNodeRootItem->appendRow(itemMaterial);
            myNodeRootItem->appendRow(itemStatistics);
        }
        myNodeRootItem->appendRow(itemBoundingBox);
        myNodeRootItem->appendRow(itemProperties);
        myNodeRootItem->appendRow(itemAdvanced);
        myNodeRootItem->appendRow(itemMeshDataSources);
        break;

    case nodeType_coordinateSystem:
    case nodeType_coordinateSystem_global:
        myNodeRootItem->appendRow(itemOrigin);
        myNodeRootItem->appendRow(itemDirectionalData);
        myNodeRootItem->appendRow(itemDefinition);
        if(theNodeType!=nodeType_coordinateSystem_global)
        {
            myNodeRootItem->appendRow(itemTransformations);
        }
        break;

    case nodeType_meshControl:
        myNodeRootItem->appendRow(itemStatistics);
        myNodeRootItem->appendRow(itemMeshViewOptions);
        myNodeRootItem->appendRow(itemSizing);
        myNodeRootItem->appendRow(itemAdvanced);
        break;

    case nodeType_structuralAnalysisBoundaryCondition_Acceleration:
    case nodeType_structuralAnalysisBoundaryCondition_RotationalVelocity:
    case nodeType_structuralAnalysisBoundaryContidion_FixedSupport:
    case nodeType_structuralAnalysisBoundaryCondition_FrictionlessSupport:
    case nodeType_structuralAnalysisBoundaryCondition_CylindricalSupport:
    case nodeType_structuralAnalysisBoundaryCondition_CompressionOnlySupport:
    case nodeType_structuralAnalysisBoundaryCondition_RemoteForce:
    case nodeType_structuralAnalysisBoundaryCondition_Force:
    case nodeType_structuralAnalysisBoundaryCondition_Displacement:
    case nodeType_structuralAnalysisBoundaryCondition_RemoteDisplacement:
    case nodeType_structuralAnalysisBoundaryCondition_RemoteRotation:
    case nodeType_structuralAnalysisBoundaryCondition_Moment:
    case nodeType_structuralAnalysisBoltPretension:
    case nodeType_structuralAnalysisBoundaryCondition_Pressure:
    case nodeType_structuralAnalysisThermalCondition:
        myNodeRootItem->appendRow(itemScope);
        myNodeRootItem->appendRow(itemDefinition);
        myNodeRootItem->appendRow(itemAdvanced);
        myNodeRootItem->appendRow(itemHidden);
        myNodeRootItem->appendRow(itemGraphicObject);
        myNodeRootItem->appendRow(itemMeshDataSources);
        break;

    case nodeType_pointMass:
        myNodeRootItem->appendRow(itemScope);
        myNodeRootItem->appendRow(itemDefinition);
        myNodeRootItem->appendRow(itemPosition);
        myNodeRootItem->appendRow(itemGraphicObject);
        myNodeRootItem->appendRow(itemMeshDataSources);
        break;

    case nodeType_structuralAnalysisBoundaryCondition_ImportedTemperatureDistribution:
        myNodeRootItem->appendRow(itemScope);
        myNodeRootItem->appendRow(itemDefinition);
        myNodeRootItem->appendRow(itemCouplingTime);
        break;

    case nodeType_thermalAnalysisTemperature:
    case nodeType_thermalAnalysisThermalFlux:
    case nodeType_thermalAnalysisThermalFlow:
    case nodeType_thermalAnalysisThermalPower:
    case nodeType_thermalAnalysisRadiation:
    case nodeType_thermalAnalysisAdiabaticWall:
        myNodeRootItem->appendRow(itemScope);
        myNodeRootItem->appendRow(itemDefinition);
        myNodeRootItem->appendRow(itemAdvanced);
        myNodeRootItem->appendRow(itemHidden);
        myNodeRootItem->appendRow(itemGraphicObject);
        myNodeRootItem->appendRow(itemMeshDataSources);
        break;

    case nodeType_thermalAnalysisConvection:
        myNodeRootItem->appendRow(itemScope);
        myNodeRootItem->appendRow(itemDefinition);
        myNodeRootItem->appendRow(itemConvection);
        myNodeRootItem->appendRow(itemAdvanced);
        myNodeRootItem->appendRow(itemHidden);
        myNodeRootItem->appendRow(itemGraphicObject);
        myNodeRootItem->appendRow(itemMeshDataSources);
        break;

    case nodeType_electrostaticPotential:
        myNodeRootItem->appendRow(itemScope);
        myNodeRootItem->appendRow(itemDefinition);
        myNodeRootItem->appendRow(itemEmitter);
        myNodeRootItem->appendRow(itemAdvanced);
        //myNodeRootItem->appendRow(itemHidden);
        //myNodeRootItem->appendRow(itemGraphicObject);
        myNodeRootItem->appendRow(itemMeshDataSources);
        break;

        //! -------------
        //! Mesh section
        //! -------------
    case nodeType_meshBodyMeshMethod:
    case nodeType_meshBodyMeshControl:
    case nodeType_meshFaceSize:
    case nodeType_meshEdgeSize:
    case nodeType_meshVertexSize:
        myNodeRootItem->appendRow(itemScope);
        myNodeRootItem->appendRow(itemDefinition);
        myNodeRootItem->appendRow(itemAdvanced);
        myNodeRootItem->appendRow(itemGraphicObject);
        break;

    case nodeType_meshMeshType:
        myNodeRootItem->appendRow(itemScope);
        myNodeRootItem->appendRow(itemDefinition);
        break;

    case nodeType_meshMeshMetric:
        myNodeRootItem->appendRow(itemScope);
        myNodeRootItem->appendRow(itemDefinition);
        break;

    case nodeType_meshPrismaticLayer:
        myNodeRootItem->appendRow(itemScope);
        myNodeRootItem->appendRow(itemDefinition);
        myNodeRootItem->appendRow(itemSmoothing);
        myNodeRootItem->appendRow(itemAdvanced);
        break;

        //! ------------------
        //! Analysis settings
        //! ------------------
    case nodeType_thermalAnalysisSettings:
    case nodeType_structuralAnalysisSettings:
    case nodeType_combinedAnalysisSettings:
        myNodeRootItem->appendRow(itemStepControls);
        myNodeRootItem->appendRow(itemSolverControl);
        myNodeRootItem->appendRow(itemConvergenceCriteria);
        myNodeRootItem->appendRow(itemTimeIncrementation);
        myNodeRootItem->appendRow(itemCutBack);
        myNodeRootItem->appendRow(itemLineSearch);
        myNodeRootItem->appendRow(itemOutputSettings);
        break;
    case nodeType_particlesInFieldsAnalysisSettings:
        myNodeRootItem->appendRow(itemStepControls);
        myNodeRootItem->appendRow(itemSolverControl);
        break;

    case nodeType_root:
        myNodeRootItem->appendRow(itemBackground);
        myNodeRootItem->appendRow(itemLighting);
        break;

        //! ------------
        //! Connections
        //! ------------
    case nodeType_connection:
        myNodeRootItem->appendRow(itemAutoDetection);
        myNodeRootItem->appendRow(itemTransparency);
        break;

    case nodeType_connectionGroup:
        myNodeRootItem->appendRow(itemAutoDetection);
        myNodeRootItem->appendRow(itemScope);
        myNodeRootItem->appendRow(itemDefinition);
        break;

    case nodeType_connectionPair:
        myNodeRootItem->appendRow(itemScope);
        myNodeRootItem->appendRow(itemDefinition);
        myNodeRootItem->appendRow(itemAdvanced);
        myNodeRootItem->appendRow(itemMeshDataSources);
        break;

        //! ----------------
        //! Named selection
        //! ----------------
    case nodeType_namedSelection:
        myNodeRootItem->appendRow(itemDisplay);
        break;

    case nodeType_namedSelectionGeometry:
        myNodeRootItem->appendRow(itemScope);
        myNodeRootItem->appendRow(itemDefinition);
        myNodeRootItem->appendRow(itemStatistics);
        break;

    case nodeType_namedSelectionElement:
        myNodeRootItem->appendRow(itemScope);
        myNodeRootItem->appendRow(itemDefinition);
        myNodeRootItem->appendRow(itemMeshDataSources);
        myNodeRootItem->appendRow(itemStatistics);
        break;

        //! ---------
        //! Solution
        //! ---------
    case nodeType_thermalAnalysisSolution:
    case nodeType_StructuralAnalysisSolution:
    case nodeType_combinedAnalysisSolution:
    case nodeType_particlesInFieldsSolution:
        myNodeRootItem->appendRow(itemInformation);
        break;


    case nodeType_structuralAnalysis:
    case nodeType_thermalAnalysis:
    case nodeType_combinedAnalysis:
        myNodeRootItem->appendRow(itemEnvironmentTemperature);
        break;

    case nodeType_remotePoint:
        myNodeRootItem->appendRow(itemScope);
        myNodeRootItem->appendRow(itemDefinition);
        myNodeRootItem->appendRow(itemAdvanced);
        myNodeRootItem->appendRow(itemGraphicObject);
        break;

    case nodeType_mapper:
        myNodeRootItem->appendRow(itemDefinition);
        myNodeRootItem->appendRow(itemColorBox);
        break;

    case nodeType_OpenFoamScalarData:
        myNodeRootItem->appendRow(itemDefinition);
        myNodeRootItem->appendRow(itemOutputSettings);
        break;

    case nodeType_importedBodyScalar:
        myNodeRootItem->appendRow(itemScope);
        myNodeRootItem->appendRow(itemDefinition);
        myNodeRootItem->appendRow(itemCouplingTime);
        myNodeRootItem->appendRow(itemGraphicObject);
        myNodeRootItem->appendRow(itemAdvanced);
        myNodeRootItem->appendRow(itemColorBox);
        break;

    case nodeType_postObject:
        myNodeRootItem->appendRow(itemScope);
        myNodeRootItem->appendRow(itemDefinition);
        myNodeRootItem->appendRow(itemGraphicObject);
        myNodeRootItem->appendRow(itemColorBox);
        break;

        //! -------------------------------
        //! structural environment results
        //! -------------------------------
    case nodeType_solutionStructuralNodalDisplacement:
    case nodeType_solutionStructuralStress:
    case nodeType_solutionStructuralTotalStrain:
    case nodeType_solutionStructuralMechanicalStrain:
    case nodeType_solutionStructuralThermalStrain:
    case nodeType_solutionStructuralFatigueTool:
    case nodeType_solutionStructuralTemperature:
    case nodeType_solutionStructuralEquivalentPlasticStrain:
    case nodeType_solutionStructuralNodalForces:
    case nodeType_solutionStructuralContact:
    case nodeType_solutionStructuralGamma:
    case nodeType_solutionStructuralReactionForce:
        myNodeRootItem->appendRow(itemScope);
        myNodeRootItem->appendRow(itemDefinition);
        myNodeRootItem->appendRow(itemGraphicObject);
        myNodeRootItem->appendRow(itemColorBox);
        break;

        //! ----------------------------
        //! thermal environment results
        //! ----------------------------
    case SimulationNodeClass::nodeType_solutionThermalTemperature:
    case SimulationNodeClass::nodeType_solutionThermalFlux:
        myNodeRootItem->appendRow(itemScope);
        myNodeRootItem->appendRow(itemDefinition);
        myNodeRootItem->appendRow(itemGraphicObject);
        myNodeRootItem->appendRow(itemColorBox);
        break;

        //! ---------------------
        //! Solution information
        //! ---------------------
    case nodeType_StructuralAnalysisSolutionInformation:
    case nodeType_thermalAnalysisSolutionInformation:
    case nodeType_combinedAnalysisSolutionInformation:
    case nodeType_particlesInFieldsSolutionInformation:
        myNodeRootItem->appendRow(itemSolutionInformation);
        myNodeRootItem->appendRow(itemHidden);
        break;

    case nodeType_repairTool:
        myNodeRootItem->appendRow(itemScope);
        myNodeRootItem->appendRow(itemDefinition);
        break;

    case nodeType_meshMethod:
        myNodeRootItem->appendRow(itemScope);
        myNodeRootItem->appendRow(itemMethod);
        myNodeRootItem->appendRow(itemDefinition);
        myNodeRootItem->appendRow(itemDefeaturing);
        myNodeRootItem->appendRow(itemTetWild);
        break;

    case nodeType_import:
        myNodeRootItem->appendRow(itemDefinition);
        break;

        //! -------------
        //! model change
        //! -------------
    case nodeType_modelChange:
        myNodeRootItem->appendRow(itemDefinition);
        myNodeRootItem->appendRow(itemScope);
        break;


#ifdef COSTAMP_VERSION
    case nodeType_timeStepBuilder:
        myNodeRootItem->appendRow(itemDefinition);
        break;

    case nodeType_processParametersClosureForce:
    case nodeType_processParametersPressure:
        myNodeRootItem->appendRow(itemScope);
        myNodeRootItem->appendRow(itemDefinition);
#endif
    }

    //! -----------------------
    //! for all the node types
    //! -----------------------
    myNodeRootItem->appendRow(itemValidity);
    myNodeRootItem->appendRow(itemIdentifier);
}

//! --------------------------
//! function: removeProperty
//! details:  helper function
//! --------------------------
bool SimulationNodeClass::removeProperty(const QString &keyPropertyName)
{
    QExtendedStandardItem *item = this->getPropertyItem(keyPropertyName);
    if(item!=NULL)
    {
        this->getModel()->removeRow(item->index().row(),item->index().parent());
        return true;
    }
    return false;
}

//! --------------------------
//! function: replaceProperty
//! details:
//! --------------------------
bool SimulationNodeClass::replaceProperty(const QString &keyPropertyName, Property theProp)
{
    QExtendedStandardItem *item= this->getPropertyItem(keyPropertyName);
    if(item!=NULL)
    {
        QVariant data;
        data.setValue(theProp);
        item->setData(data,Qt::UserRole);
        return true;
    }
    else
    {
        cerr<<"SimulationNodeClass::replaceProperty()->____error: the property: "<<keyPropertyName.toStdString()<<" cannot be replaced____"<<endl;
        return false;
    }
    return false;
}

//! -------------------------------------
//! function: type
//! details:  return the type as QString
//! -------------------------------------
const QString SimulationNodeClass::type() const
{
    const QMetaObject &mo = SimulationNodeClass::staticMetaObject;
    int index = mo.indexOfEnumerator("nodeType");
    QMetaEnum metaEnum = mo.enumerator(index);
    return metaEnum.valueToKey(myNodeType);
}

//! --------------------
//! function: typeToInt
//! details:
//! --------------------
int SimulationNodeClass::typeToInt(const QString& type)
{
    const QMetaObject &mo = SimulationNodeClass::staticMetaObject;
    int index = mo.indexOfEnumerator("nodeType");
    QMetaEnum metaEnum = mo.enumerator(index);
    return metaEnum.keysToValue(type.toStdString().c_str());
}

//! -------------------
//! function: setType1
//! details:
//! -------------------
void SimulationNodeClass::setType1(QString type)
{
    const QMetaObject &mo = SimulationNodeClass::staticMetaObject;
    int index = mo.indexOfEnumerator("nodeType");
    QMetaEnum metaEnum = mo.enumerator(index);
    int value = metaEnum.keyToValue(qPrintable(type));
    myNodeType = static_cast<nodeType>(value);
}

//! -------------------------
//! function: hasTabularData
//! details:
//! -------------------------
bool SimulationNodeClass::hasTabularData()
{
    switch(this->getType())
    {
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Displacement:
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteDisplacement:
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteRotation:
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RotationalVelocity:
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Acceleration:
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Force:
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteForce:
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Moment:
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Pressure:
    case SimulationNodeClass::nodeType_structuralAnalysisThermalCondition:
    case SimulationNodeClass::nodeType_structuralAnalysisBoltPretension:

    case SimulationNodeClass::nodeType_thermalAnalysisConvection:
    case SimulationNodeClass::nodeType_thermalAnalysisRadiation:
    case SimulationNodeClass::nodeType_thermalAnalysisThermalFlux:
    case SimulationNodeClass::nodeType_thermalAnalysisThermalFlow:
    case SimulationNodeClass::nodeType_thermalAnalysisTemperature:
    case SimulationNodeClass::nodeType_thermalAnalysisThermalPower:
        return true;
        break;
    default:
        break;
    }
    return false;
}

//! -------------------------------------------------
//! function: isAnalysisRoot
//! details:  true if the <this> is an analysis root
//! -------------------------------------------------
bool SimulationNodeClass::isAnalysisRoot()
{
    if(myNodeType == SimulationNodeClass::nodeType_structuralAnalysis ||
            myNodeType == SimulationNodeClass::nodeType_thermalAnalysis ||
            myNodeType == SimulationNodeClass::nodeType_combinedAnalysis ||
            myNodeType == SimulationNodeClass::nodeType_particlesInFieldsAnalysis)
    return true;
    return false;
}

//! -----------------------------
//! function: generateTimeString
//! details:  monothonic timer
//! -----------------------------
std::string SimulationNodeClass::generateTimeString()
{
    using namespace std::chrono;

    // get current time
    auto now = system_clock::now();

    // get number of milliseconds for the current second
    // (remainder after division into seconds)
    auto ms = duration_cast<milliseconds>(now.time_since_epoch()) % 1000;

    // convert to std::time_t in order to convert to std::tm (broken time)
    auto timer = system_clock::to_time_t(now);

    // convert to broken time
    std::tm bt = *std::localtime(&timer);

    std::ostringstream oss;
    oss << std::put_time(&bt, "%Y%m%d%H%M%S"); // HH:MM:SS
    oss << std::setfill('0') << std::setw(3) << ms.count();
    return oss.str();
}

//! ---------------------
//! function: addTimeTag
//! details:
//! ---------------------
void SimulationNodeClass::addTimeTag()
{
    //! ---------------------
    //! Time tag - wait 5 ms
    //! ---------------------
    QThread::msleep(5);
    QVariant data;
    data.setValue(QString::fromStdString(SimulationNodeClass::generateTimeString()));
    Property prop_timeTag("Time tag",data,Property::PropertyGroup_Identifier);
    this->addProperty(prop_timeTag);
}

//! --------------------------------
//! function: isSimulationSetUpNode
//! details:
//! --------------------------------
bool SimulationNodeClass::isSimulationSetUpNode()
{
    if(myNodeType == SimulationNodeClass::nodeType_electrostaticPotential ||
            myNodeType == SimulationNodeClass::nodeType_magneticField ||
            myNodeType == SimulationNodeClass::nodeType_particlesInFieldsParticlePack ||
            myNodeType == SimulationNodeClass::nodeType_thermalAnalysisAdiabaticWall ||
            myNodeType == SimulationNodeClass::nodeType_thermalAnalysisConvection ||
            myNodeType == SimulationNodeClass::nodeType_thermalAnalysisRadiation ||
            myNodeType == SimulationNodeClass::nodeType_thermalAnalysisTemperature ||
            myNodeType == SimulationNodeClass::nodeType_thermalAnalysisThermalFlow ||
            myNodeType == SimulationNodeClass::nodeType_thermalAnalysisThermalFlux ||
            myNodeType == SimulationNodeClass::nodeType_thermalAnalysisThermalPower ||
            myNodeType == SimulationNodeClass::nodeType_structuralAnalysisBoltPretension ||
            myNodeType == SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Acceleration ||
            myNodeType == SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_CylindricalSupport ||
            myNodeType == SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_CompressionOnlySupport ||
            myNodeType == SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Displacement ||
            myNodeType == SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Force ||
            myNodeType == SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_FrictionlessSupport ||
            myNodeType == SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Moment ||
            myNodeType == SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Pressure ||
            myNodeType == SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteDisplacement ||
            myNodeType == SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteForce ||
            myNodeType == SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteRotation ||
            myNodeType == SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RotationalVelocity ||
            myNodeType == SimulationNodeClass::nodeType_structuralAnalysisBoundaryContidion_FixedSupport ||
            myNodeType == SimulationNodeClass::nodeType_structuralAnalysisThermalCondition ||
            myNodeType == SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_ImportedTemperatureDistribution ||
            myNodeType == SimulationNodeClass::nodeType_mapper ||
            myNodeType == SimulationNodeClass::nodeType_modelChange ||
            myNodeType == SimulationNodeClass::nodeType_mapper ||
            myNodeType == SimulationNodeClass::nodeType_modelChange)
        return true;

#ifdef COSTAMP_VERSION
    if(myNodeType == SimulationNodeClass::nodeType_timeStepBuilder) return true;
#endif

    return false;
}

//! -------------------------------------
//! function: isChildSimulationSetUpNode
//! details:
//! -------------------------------------
bool SimulationNodeClass::isChildSimulationSetUpNode()
{
    if(myNodeType == SimulationNodeClass::nodeType_importedBodyScalar ||
            myNodeType == SimulationNodeClass::nodeType_OpenFoamScalarData
            )
        return true;
    return false;
}


//! -----------------------------
//! function: isAnalysisSettings
//! details:  helper
//! -----------------------------
bool SimulationNodeClass::isAnalysisSettings()
{
    if(myNodeType == SimulationNodeClass::nodeType_structuralAnalysisSettings ||
            myNodeType == SimulationNodeClass::nodeType_thermalAnalysisSettings ||
            myNodeType == SimulationNodeClass::nodeType_combinedAnalysisSettings ||
            myNodeType == SimulationNodeClass::nodeType_particlesInFieldsAnalysisSettings)
        return true;
    return false;
}

//! ---------------------
//! function: isSolution
//! details:  helper
//! ---------------------
bool SimulationNodeClass::isSolution()
{
    if(myNodeType == SimulationNodeClass::nodeType_StructuralAnalysisSolution ||
            myNodeType == SimulationNodeClass::nodeType_thermalAnalysisSolution ||
            myNodeType == SimulationNodeClass::nodeType_combinedAnalysisSolution ||
            myNodeType == SimulationNodeClass::nodeType_particlesInFieldsSolution)
        return true;
    return false;
}

//! --------------------------------
//! function: isSolutionInformation
//! details:  helper
//! --------------------------------
bool SimulationNodeClass::isSolutionInformation()
{
    if(myNodeType == SimulationNodeClass::nodeType_StructuralAnalysisSolutionInformation ||
            myNodeType == SimulationNodeClass::nodeType_thermalAnalysisSolutionInformation ||
            myNodeType == SimulationNodeClass::nodeType_combinedAnalysisSolutionInformation ||
            myNodeType == SimulationNodeClass::nodeType_particlesInFieldsSolutionInformation)
        return true;
    return false;
}

//! ---------------------------
//! function: isAnalysisResult
//! details:  helper
//! ---------------------------
bool SimulationNodeClass::isAnalysisResult()
{
    if(myNodeType == SimulationNodeClass::nodeType_solutionStructuralEquivalentPlasticStrain ||
            myNodeType == SimulationNodeClass::nodeType_solutionStructuralFatigueTool ||
            myNodeType == SimulationNodeClass::nodeType_solutionStructuralMechanicalStrain ||
            myNodeType == SimulationNodeClass::nodeType_solutionStructuralNodalDisplacement ||
            myNodeType == SimulationNodeClass::nodeType_solutionStructuralNodalForces ||
            myNodeType == SimulationNodeClass::nodeType_solutionStructuralStress ||
            myNodeType == SimulationNodeClass::nodeType_solutionStructuralFatigueTool ||
            myNodeType == SimulationNodeClass::nodeType_solutionStructuralTemperature ||
            myNodeType == SimulationNodeClass::nodeType_solutionStructuralThermalStrain ||
            myNodeType == SimulationNodeClass::nodeType_solutionStructuralTotalStrain ||
            myNodeType == SimulationNodeClass::nodeType_solutionThermalTemperature ||
            myNodeType == SimulationNodeClass::nodeType_solutionThermalFlux ||
            myNodeType == SimulationNodeClass::nodeType_solutionStructuralContact ||
            myNodeType == SimulationNodeClass::nodeType_postObject)
        return true;
    return false;
}
Q_DECLARE_OPAQUE_POINTER(SimulationNodeClass*)
