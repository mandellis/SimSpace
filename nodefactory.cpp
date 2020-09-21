//! ----------------
//! custom includes
//! ----------------
#include "nodefactory.h"
#include "listofshape.h"
#include "property.h"
#include "mydefines.h"
#include "topologytools.h"
#include "mydefines.h"
#include "tools.h"
#include "parser.h"
#include "prebuiltcontactoptions.h"
#include "load.h"
#include <occPreGLwidget.h>

#include "markers.h"
#include "handle_ais_doublearrowmarker_reg.h"
#include "handle_ais_trihedron_reg.h"
#include <indexedmapofmeshdatasources.h>        //! registered typedef for QMap<int,occHandle(MeshVS_DataSource)>
#include <connectionpairgenerationoptions.h>

//! ----
//! OCC
//! ----
#include <TopoDS.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Solid.hxx>
#include <gp_Dir.hxx>

//! ---
//! Qt
//! ---
#include <QDateTime>
#include <QMessageBox>

//! ----
//! C++
//! ----
#include <memory>

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
nodeFactory::nodeFactory()
{
    ;
}

//! -------------------------------------
//! function: nodeFromScratch
//! details:  create a node from scratch
//! -------------------------------------
SimulationNodeClass* nodeFactory::nodeFromScratch(SimulationNodeClass::nodeType type,
                                                  meshDataBase *mDB,
                                                  const occHandle(AIS_InteractiveContext)& CTX,
                                                  QVariant addOptions)
{
    cout<<"nodeFactory::nodeFromScratch()->____function called____"<<endl;

    //! ---------------------
    //! vector of properties
    //! ---------------------
    QVector<Property> vecProp;

    //! -----------------
    //! name of the node
    //! -----------------
    QString name;

    //! ----------
    //! node data
    //! ----------
    QVariant data;

    //! -------------------
    //! suppression status
    //! -------------------
    data.setValue(Property::SuppressionStatus_Active);
    Property prop_suppressed("Suppressed",data,Property::PropertyGroup_Definition);

    //! -----------------
    //! "Scoping method"
    //! -----------------
    data.setValue(Property::ScopingMethod_GeometrySelection);
    Property prop_scopingMethod("Scoping method",data,Property::PropertyGroup_Scope);

    ListOfShape scope;
    if(!CTX.IsNull()) for(CTX->InitSelected(); CTX->MoreSelected(); CTX->NextSelected()) scope.Append(CTX->SelectedShape());

    //! -------------------
    //! "Scope" and "Tags"
    //! -------------------
    std::vector<GeometryTag> vecLoc = TopologyTools::generateLocationPairs(mDB,scope);
    data.setValue(vecLoc);
    Property prop_scope("Geometry",data,Property::PropertyGroup_Scope);
    Property prop_tags("Tags",data,Property::PropertyGroup_Scope);

    switch(type)
    {
    case SimulationNodeClass::nodeType_pointMass:
    {
        name = "Point mass"; //bubi

        //! ------------
        //! under scope
        //! ------------
        vecProp.push_back(prop_scopingMethod);
        vecProp.push_back(prop_scope);
        vecProp.push_back(prop_tags);

        //! ---------------
        //! under position
        //! ---------------
        data.setValue(0.0);
        vecProp.push_back(Property("X coordinate",data,Property::PropertyGroup_Position));
        vecProp.push_back(Property("Y coordinate",data,Property::PropertyGroup_Position));
        vecProp.push_back(Property("Z coordinate",data,Property::PropertyGroup_Position));

        //! -----------------
        //! under definition
        //! -----------------
        vecProp.push_back(prop_suppressed);

        vecProp.push_back(Property("Mass",data,Property::PropertyGroup_Definition));
        vecProp.push_back(Property("Jx",data,Property::PropertyGroup_Definition));
        vecProp.push_back(Property("Jy",data,Property::PropertyGroup_Definition));
        vecProp.push_back(Property("Jz",data,Property::PropertyGroup_Definition));
    }
        break;
    //! -------------------------
    //! temperature distribution
    //! -------------------------
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_ImportedTemperatureDistribution:
    {
        name = "Imported body temperature";

        //! ------------
        //! under scope
        //! ------------
        vecProp.push_back(prop_scopingMethod);
        vecProp.push_back(prop_scope);
        vecProp.push_back(prop_tags);

        //! -----------------
        //! under definition
        //! -----------------
        vecProp.push_back(prop_suppressed);

        //! ------------------------------------------
        //! use a smart pointer for this "dummy" node
        //! ------------------------------------------
        SimulationNodeClass *aNode = new SimulationNodeClass("No selection",SimulationNodeClass::nodeType_NULL, QVector<Property>());
        //aNode->addTimeTag();

        QStandardItem *item = new QStandardItem();
        item->setData(aNode->getName(),Qt::DisplayRole);
        data.setValue(aNode);
        item->setData(data,Qt::UserRole);

        void *p = (void*)(item);
        data.setValue(p);
        vecProp.push_back(Property("Analysis",data,Property::PropertyGroup_Definition));

        //! ------------------------------------------
        //! use a smart pointer for this "dummy" node
        //! ------------------------------------------
        SimulationNodeClass *aNode_ = new SimulationNodeClass("No temperature result selected",SimulationNodeClass::nodeType_NULL, QVector<Property>());
        //aNode_->addTimeTag();

        QStandardItem *itemTemperatureDistribution = new QStandardItem();
        data.setValue(aNode_->getName());
        itemTemperatureDistribution->setData(data,Qt::DisplayRole);
        data.setValue(aNode_);
        itemTemperatureDistribution->setData(data,Qt::UserRole);

        p = (void*)(itemTemperatureDistribution);
        data.setValue(p);
        vecProp.push_back(Property("Imported body temperature",data,Property::PropertyGroup_Definition));

        //! --------------
        //! Coupling time
        //! --------------
        data.setValue(1);
        Property prop_couplingTimeStep("Structural time step",data,Property::PropertyGroup_CouplingTime);
        vecProp.push_back(prop_couplingTimeStep);
    }
        break;

    //! -------------
    //! fatigue tool
    //! -------------
    case SimulationNodeClass::nodeType_solutionStructuralFatigueTool:
    {
        name = "Fatigue tool";

        //! --------------------------------------------
        //! fatigue method
        //! 0 -> Basquin Coffin Manson - default option
        //! 1 -> [...]
        //! --------------------------------------------
        Property prop_fatigueMethod("Fatigue algo",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_fatigueMethod);

        /*
        //! ---------------------------
        //! fatigue curve coefficients
        //! ---------------------------
        double val0,val1,val2,val3;
        val0 = 1407.6;  //! rupture stress
        val1 = 0.6;     //! rupture strain
        val2 = 0.09;    //! Ramberg Osgood
        val3 = 0.07;    //! b

        data.setValue(val0);
        Property prop_coefficient_0("Rupture stress",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_coefficient_0);

        data.setValue(val1);
        Property prop_coefficient_1("Rupture strain",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_coefficient_1);

        data.setValue(val2);
        Property prop_coefficient_2("n'",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_coefficient_2);

        data.setValue(val3);
        Property prop_coefficient_3("b",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_coefficient_3);
        */

        data.setValue(0);
        Property prop_stressStrainSource("Stress/strain source",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_stressStrainSource);

        data.setValue(1);
        Property prop_component("Component",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_component);

        //! -----------------------
        //! triaxiality correction
        //! -----------------------
        data.setValue(0);
        Property prop_triaxialityCorrection("Triaxiality correction",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_triaxialityCorrection);

        //! -----------------
        //! number of cycles
        //! -----------------
        data.setValue(1000);
        Property prop_numberOfCycles("Number of cycles",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_numberOfCycles);

        //! ---------------------------------------------------
        //! type of scale: "0" => "Autoscale"; "1" => "Custom"
        //! ---------------------------------------------------
        int typeOfScale = 0;
        data.setValue(typeOfScale);
        Property prop_scaleType("Scale type",data,Property::PropertyGroup_ColorBox);
        vecProp.push_back(prop_scaleType);

        //! ------------------------------------------------
        //! Number of intervals: "0" => "Program controlled
        //! ------------------------------------------------
        int NbIntervals = INITIAL_NUMBER_OF_COLORBOX_LEVELS;
        data.setValue(NbIntervals);
        Property prop_NbIntervals("# intervals",data,Property::PropertyGroup_ColorBox);
        vecProp.push_back(prop_NbIntervals);

        //! -----------------------------------
        //! Mapping - 0 linear 1 - logarithmic
        //! -----------------------------------
        int mapping = 0;
        data.setValue(0);
        Property prop_mapping("Mapping",data,Property::PropertyGroup_ColorBox);
        vecProp.push_back(prop_mapping);

        //! ------------
        //! under scope
        //! ------------
        vecProp.push_back(prop_scopingMethod);
        vecProp.push_back(prop_scope);
        vecProp.push_back(prop_tags);
    }
        break;

        //! -------
        //! strain
        //! -------
    case SimulationNodeClass::nodeType_solutionStructuralTotalStrain:
    case SimulationNodeClass::nodeType_solutionStructuralMechanicalStrain:
    case SimulationNodeClass::nodeType_solutionStructuralThermalStrain:
    {
        //! read options
        int typeOfStrain = addOptions.toInt();

        //! prepare the name
        switch(typeOfStrain)
        {
        case 0: name = "Total equivalent strain"; break;
        case 1: name = "Total strain intensity"; break;
        case 2: name = "Total maximum principal strain"; break;
        case 3: name = "Total middle principal strain"; break;
        case 4: name = "Total minimum principal strain"; break;

        case 5: name = "Mechanical equivalent strain"; break;
        case 6: name = "Mechanical strain intensity"; break;
        case 7: name = "Mechanical maximum principal strain"; break;
        case 8: name = "Mechanical middle principal strain"; break;
        case 9: name = "Mechanical minimum principal strain"; break;

        case 10: name = "Thermal equivalent strain"; break;
        case 11: name = "Thermal strain intensity"; break;
        case 12: name = "Thermal maximum principal strain"; break;
        case 13: name = "Thermal middle principal strain"; break;
        case 14: name = "Thermal minimum principal strain"; break;
        }

        //! ------------------------------------------------------------------
        //! redefine tags and scope:
        //! initially the scope is "All bodies" for the structural diagnostic
        //! ------------------------------------------------------------------
        int NbBodies = mDB->bodyMap.size();
        ListOfShape scope;
        for(int i=1; i<= NbBodies; i++)
        {
            TopoDS_Solid aSolid = TopoDS::Solid(mDB->bodyMap.value(i));
            scope.Append(aSolid);
        }
        std::vector<GeometryTag> vecLocAllBodies = TopologyTools::generateLocationPairs(mDB, scope);

        data.setValue(vecLocAllBodies);
        Property prop_scopeAllBodies("Geometry",data,Property::PropertyGroup_Scope);
        Property prop_tagsAllBodies("Tags",data,Property::PropertyGroup_Scope);

        //! under scope
        vecProp.push_back(prop_scopingMethod);
        vecProp.push_back(prop_scopeAllBodies);
        vecProp.push_back(prop_tagsAllBodies);

        /*
        //! under scope
        vecProp.push_back(prop_scopingMethod);
        vecProp.push_back(prop_scope);
        vecProp.push_back(prop_tags);
        */

        //! under definition
        Property prop_type("Type ",typeOfStrain,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_type);

        //! time info
        data.setValue(0.0);
        Property prop_By("By",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_By);
        Property prop_displayTime("Display time",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_displayTime);

        //! mode
        data.setValue(0);
        Property prop_modeNumber("Mode number",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_modeNumber);

        //! type of scale: "0" => "Autoscale"; "1" => "Custom"
        int typeOfScale = 0;
        data.setValue(typeOfScale);
        Property prop_scaleType("Scale type",data,Property::PropertyGroup_ColorBox);
        vecProp.push_back(prop_scaleType);

        //! under definition
        //Property prop_type("Type ",typeOfStrain,Property::PropertyGroup_Definition);
        //vecProp.push_back(prop_type);

        //! Number of intervals: "0" => "Program controlled
        int NbIntervals = INITIAL_NUMBER_OF_COLORBOX_LEVELS;
        data.setValue(NbIntervals);
        Property prop_NbIntervals("# intervals",data,Property::PropertyGroup_ColorBox);
        vecProp.push_back(prop_NbIntervals);
    }
        break;

        //! -------
        //! stress
        //! -------
    case SimulationNodeClass::nodeType_solutionStructuralStress:
    {
        //! read options
        int typeOfStress = addOptions.toInt();

        //! prepare the name
        switch(typeOfStress)
        {
        case 0: name = "Equivalent stress"; break;
        case 1: name = "Stress intensity"; break;
        case 2: name = "Maximum principal stress"; break;
        case 3: name = "Middle principal stress"; break;
        case 4: name = "Minimum principal stress"; break;
        case 5: name = "Normal stress X"; break;
        case 6: name = "Normal stress Y"; break;
        case 7: name = "Normal stress Z"; break;
        case 8: name = "Shear stress XY"; break;
        case 9: name = "Shear stress YZ"; break;
        case 10: name = "Shear stress XZ"; break;
        }

        //! -------------------------------------------------------------------
        //! redefine tags and scope:
        //! initially the scope is "All bodies" for the structural diagnostic
        //! -------------------------------------------------------------------
        int NbBodies = mDB->bodyMap.size();
        ListOfShape scope;
        for(int i=1; i<= NbBodies; i++)
        {
            TopoDS_Solid aSolid = TopoDS::Solid(mDB->bodyMap.value(i));
            scope.Append(aSolid);
        }
        std::vector<GeometryTag> vecLocAllBodies = TopologyTools::generateLocationPairs(mDB, scope);

        data.setValue(vecLocAllBodies);
        Property prop_scopeAllBodies("Geometry",data,Property::PropertyGroup_Scope);
        Property prop_tagsAllBodies("Tags",data,Property::PropertyGroup_Scope);

        //! under scope
        vecProp.push_back(prop_scopingMethod);
        vecProp.push_back(prop_scopeAllBodies);
        vecProp.push_back(prop_tagsAllBodies);

        /*
        //! under scope
        vecProp.push_back(prop_scopingMethod);
        vecProp.push_back(prop_scope);
        vecProp.push_back(prop_tags);
        */

        //! under definition
        Property prop_type("Type ",typeOfStress,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_type);

        //! time info
        data.setValue(0.0);
        Property prop_By("By",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_By);
        Property prop_displayTime("Display time",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_displayTime);

        //! mode
        data.setValue(0);
        Property prop_modeNumber("Mode number",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_modeNumber);
        //Property prop_type("Type ",typeOfStress,Property::PropertyGroup_Definition);
        //vecProp.push_back(prop_type);

        //! type of scale: "0" => "Autoscale"; "1" => "Custom"
        int typeOfScale = 0;
        data.setValue(typeOfScale);
        Property prop_scaleType("Scale type",data,Property::PropertyGroup_ColorBox);
        vecProp.push_back(prop_scaleType);

        //! Number of intervals: "0" => "Program controlled
        int NbIntervals = INITIAL_NUMBER_OF_COLORBOX_LEVELS;
        data.setValue(NbIntervals);
        Property prop_NbIntervals("# intervals",data,Property::PropertyGroup_ColorBox);
        vecProp.push_back(prop_NbIntervals);
    }
        break;

        //! --------------------
        //! nodal displacements
        //! --------------------
    case SimulationNodeClass::nodeType_solutionStructuralNodalDisplacement:
    {
        //! read options
        int typeOfDisplacement = addOptions.toInt();

        //! prepare the name
        switch(typeOfDisplacement)
        {
        case 0: name = "Total displacement"; break;
        case 1: name = "Directional displacement X"; break;
        case 2: name = "Directional displacement Y"; break;
        case 3: name = "Directional displacement Z"; break;
        }

        for(std::vector<GeometryTag>::iterator it = vecLoc.begin(); it!=vecLoc.end(); ++it)
        {
            const GeometryTag curLoc = *it;
            int parentShapeNr = curLoc.parentShapeNr;
            QString suffix = mDB->MapOfBodyNames.value(parentShapeNr)+" ";
            name.append(suffix);
        }

        //! ------------------------------------------------------------------
        //! redefine tags and scope:
        //! initially the scope is "All bodies" for the structural diagnostic
        //! ------------------------------------------------------------------
        int NbBodies = mDB->bodyMap.size();
        ListOfShape scope;
        for(int i=1; i<= NbBodies; i++)
        {
            TopoDS_Solid aSolid = TopoDS::Solid(mDB->bodyMap.value(i));
            scope.Append(aSolid);
        }
        std::vector<GeometryTag> vecLocAllBodies = TopologyTools::generateLocationPairs(mDB, scope);

        data.setValue(vecLocAllBodies);
        Property prop_scopeAllBodies("Geometry",data,Property::PropertyGroup_Scope);
        Property prop_tagsAllBodies("Tags",data,Property::PropertyGroup_Scope);

        //! under scope
        vecProp.push_back(prop_scopingMethod);
        vecProp.push_back(prop_scopeAllBodies);
        vecProp.push_back(prop_tagsAllBodies);

        //! under definition
        //! the component (total, X, Y, Z) is provided by the calling function through the addOptions
        Property prop_type("Type ",typeOfDisplacement,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_type);

        //! ---------------------------------------------------------------------------
        //! "By" selector: it provides options "Time" and "Set"
        //! Here "By"="Time", so the "Display time" property is created
        //! If "By" is changed into "Set", the DetailViewer replace the "Display time"
        //! with "Set Number"
        //! 0 => "Time"
        //! 1 => "Set number"
        //! ---------------------------------------------------------------------------
        data.setValue(0.0);
        Property prop_By("By",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_By);
        Property prop_displayTime("Display time",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_displayTime);

        //! mode
        data.setValue(0);
        Property prop_modeNumber("Mode number",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_modeNumber);

        //! type of scale: "0" => "Autoscale"; "1" => "Custom"
        int typeOfScale = 0;
        data.setValue(typeOfScale);
        Property prop_scaleType("Scale type",data,Property::PropertyGroup_ColorBox);
        vecProp.push_back(prop_scaleType);

        //! Number of intervals: "0" => "Program controlled
        int NbIntervals = INITIAL_NUMBER_OF_COLORBOX_LEVELS;
        data.setValue(NbIntervals);
        Property prop_NbIntervals("# intervals",data,Property::PropertyGroup_ColorBox);
        vecProp.push_back(prop_NbIntervals);
     }
        break;

        //! -----------------
        //! connection group
        //! -----------------
    case SimulationNodeClass::nodeType_connectionGroup:
    {
        name = "Connection group";

        //! ------------
        //! under scope
        //! ------------
        vecProp.push_back(prop_scopingMethod);
        vecProp.push_back(prop_scope);
        vecProp.push_back(prop_tags);

        //! --------------------------------------------------------------
        //! distance tolerance for the automatic contact search algorithm
        //! --------------------------------------------------------------
        double tolerance = 0.1;         //! a default value
        Property prop_tolerance("Tolerance",tolerance,Property::PropertyGroup_AutoDetection);
        vecProp.push_back(prop_tolerance);

        //! ------------------
        //! angular criterion
        //! ------------------
        double critAngle = 30.0;
        Property prop_angularCriterion("Angular criterion",critAngle,Property::PropertyGroup_AutoDetection);
        vecProp.push_back(prop_angularCriterion);

        //! -----------------------
        //! grouping options:
        //! "0" => by bodies
        //! "1" => by master faces
        //! "2" => by slave faces
        //! "3" = >by master body
        //! "4" => none
        //! -----------------------
        data.setValue(0);
        Property prop_grouping("Grouping",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_grouping);
    }
        break;

        //! ---------------
        //! generic sizing
        //! ---------------
    case SimulationNodeClass::nodeType_meshGenericSizing:
    {
        name = "Sizing";

        //! under definition
        vecProp.push_back(prop_suppressed);

        //! under scope
        vecProp.push_back(prop_scopingMethod);
        vecProp.push_back(prop_scope);
        vecProp.push_back(prop_tags);
    }
        break;

        //! ----------------------
        //! open foam scalar data
        //! ----------------------
    case SimulationNodeClass::nodeType_OpenFoamScalarData:
    {
        name = "OpenFoam scalar field importer";

        QString sourceDir = tools::getWorkingDir();
        data.setValue(sourceDir);
        Property prop_sourceFile("Source directory",data,Property::PropertyGroup_Definition);

        //! get the current working dir: initialize the saving directory path
        QString targetDir = tools::getWorkingDir().append("/openfoam_files");
        data.setValue(targetDir);
        Property prop_savingDirectoryPath("Target directory",data,Property::PropertyGroup_Definition);

        //! under definition
        vecProp.push_back(prop_suppressed);
        vecProp.push_back(prop_sourceFile);
        vecProp.push_back(prop_savingDirectoryPath);

        //! under definition
        //data.setValue(false);
        //Property prop_translate("Translate",data,Property::PropertyGroup_Definition);
        //vecProp.push_back(prop_translate);
        //! props timeList
        QVector<double> timeList;
        data.setValue(timeList);
        Property prop_timeList("Time list",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_timeList);

        //! under ouptup settings
        data.setValue(1);
        Property prop_outputFiles("Split data",data,Property::PropertyGroup_OutputSettings);
        vecProp.push_back(prop_outputFiles);
    }
        break;

        //! ------------
        //! post object
        //! ------------
    case SimulationNodeClass::nodeType_postObject:
    {
        name = "Post object";

        //! under scope
        //vecProp.push_back(prop_scopingMethod);
        //vecProp.push_back(prop_scope);
        vecProp.push_back(prop_tags);

        //! under definition
        vecProp.push_back(prop_suppressed);

        //! ---------------------------------------------------
        //! type of scale: "0" => "Autoscale"; "1" => "Custom"
        //! ---------------------------------------------------
        int typeOfScale = 0;
        data.setValue(typeOfScale);
        Property prop_scaleType("Scale type",data,Property::PropertyGroup_ColorBox);
        vecProp.push_back(prop_scaleType);

        //! ------------------------------------------------
        //! Number of intervals: "0" => "Program controlled
        //! ------------------------------------------------
        int NbIntervals = INITIAL_NUMBER_OF_COLORBOX_LEVELS;
        data.setValue(NbIntervals);
        Property prop_NbIntervals("# intervals",data,Property::PropertyGroup_ColorBox);
        vecProp.push_back(prop_NbIntervals);
    }
        break;

        //! --------------------------------
        //! imported body scalar (OpenFoam)
        //! --------------------------------
    case SimulationNodeClass::nodeType_importedBodyScalar:
    {
        cout<<"nodeFactory::nodeFromScratch->____creating imported body scalar____"<<endl;
        name = "Imported body scalar";

        //! -------------------------------
        //! property "Step selection mode"
        //! 0 => all steps
        //! 1 => first step
        //! 2 => last step
        //! 3 => by number
        //! 4 => time history importer
        //! -------------------------------
        data.setValue(int(0));
        Property prop_stepSelectionMode("Step selection mode",data,Property::PropertyGroup_Definition);

        QString filePath ="";
        data.setValue(filePath);
        Property prop_sourceFile("Source file",data,Property::PropertyGroup_Definition);

        //! -------------------------
        //! "0" => "in pinball" algo
        //! "1" => shape functions
        //! -------------------------
        int algo = 0;
        data.setValue(algo);
        Property prop_interpolationAlgo("Algorithm",data,Property::PropertyGroup_Advanced);

        double pinballValue = 10.0;
        data.setValue(pinballValue);
        Property prop_pinball("Pinball",data,Property::PropertyGroup_Advanced);

        //! ----------
        //! remap try
        //! ----------
        bool remap = false;
        data.setValue(remap);
        Property prop_remapFlag("Remap",data,Property::PropertyGroup_Advanced);

        //! ----------------
        //! a default value
        //! ----------------
        int NbBuckets = 10;
        data.setValue(NbBuckets);
        Property prop_Xbuckets("X buckets",data,Property::PropertyGroup_Advanced);
        Property prop_Ybuckets("Y buckets",data,Property::PropertyGroup_Advanced);
        Property prop_Zbuckets("Z buckets",data,Property::PropertyGroup_Advanced);

        //! -----------------
        //! under definition
        //! -----------------
        vecProp.push_back(prop_suppressed);
        vecProp.push_back(prop_stepSelectionMode);
        vecProp.push_back(prop_sourceFile);

        //! ------------
        //! under scope
        //! ------------
        vecProp.push_back(prop_scopingMethod);
        vecProp.push_back(prop_scope);
        vecProp.push_back(prop_tags);

        //! -----------------
        //! under definition
        //! -----------------
        //data.setValue(false);
        //Property prop_generate("Generate",data,Property::PropertyGroup_Definition);
        //vecProp.push_back(prop_generate);

        //! ---------------
        //! under advanced
        //! ---------------
        vecProp.push_back(prop_interpolationAlgo);
        vecProp.push_back(prop_pinball);
        vecProp.push_back(prop_Xbuckets);
        vecProp.push_back(prop_Ybuckets);
        vecProp.push_back(prop_Zbuckets);
        vecProp.push_back(prop_remapFlag);

        //! -------------------
        //! depends on "Remap"
        //! -------------------
        if(remap==true)
        {
            //! a default value
            int NbRemappingSteps = 2;
            data.setValue(NbRemappingSteps);
            Property prop_remappingSteps("Remapping steps",data,Property::PropertyGroup_Advanced);
            vecProp.push_back(prop_remappingSteps);
        }
    }
        break;

        //! ----------------------
        //! mapper - interpolator
        //! ----------------------
    case SimulationNodeClass::nodeType_mapper:
    {
        //cout<<"nodeFactory::nodeFromScratch->____creating item mapper3D____"<<endl;
        name = "Mapper";

        //! add suppression
        vecProp.push_back(prop_suppressed);
    }
        break;

        //! ------------------------------------------
        //! fixed, frictionless, cylindrical supports
        //! ------------------------------------------
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryContidion_FixedSupport:
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_FrictionlessSupport:
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_CylindricalSupport:
    {
        //cout<<"nodeFactory::nodeFromScratch->____creating frictionless/fixed/cylindrical____"<<endl;
        switch(type)
        {
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryContidion_FixedSupport:
            name = "Fixed support";
            break;
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_FrictionlessSupport:
            name = "Frictionless support";
            break;
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_CylindricalSupport:
        {
            name = "Cylindrical support";
            Property prop_radial("Radial",QVariant(Property::DOFfreedom_fixed),Property::PropertyGroup_Definition);
            Property prop_axial("Axial",QVariant(Property::DOFfreedom_fixed),Property::PropertyGroup_Definition);
            Property prop_tangential("Tangential",QVariant(Property::DOFfreedom_fixed),Property::PropertyGroup_Definition);

            vecProp.push_back(prop_radial);
            vecProp.push_back(prop_axial);
            vecProp.push_back(prop_tangential);
        }
            break;
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_CompressionOnlySupport:
        {
            name = "Compression only support";
        }
            break;
        }
        //! -------------------
        //! under "Definition"
        //! -------------------
        vecProp.push_back(prop_suppressed);

        //! --------------
        //! under "Scope"
        //! --------------
        vecProp.push_back(prop_scopingMethod);
        vecProp.push_back(prop_scope);
        vecProp.push_back(prop_tags);
    }
        break;

    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_CompressionOnlySupport:
    {
        name = "Compression only support";

        //! ----------------------------------------------
        //! coefficient K: if "0" => "Program controlled"
        //! ----------------------------------------------
        data.setValue(0.0);
        Property prop_K("K",data,Property::PropertyGroup_Advanced);
        vecProp.push_back(prop_K);

        //! --------------------------------------------
        //! sigma infinity: if "0" => "Program controlled"
        //! --------------------------------------------
        data.setValue(0.0);
        Property prop_sigmaInf("Sigma infinity",data,Property::PropertyGroup_Advanced);
        vecProp.push_back(prop_sigmaInf);

        //! -------------------
        //! under "Definition"
        //! -------------------
        vecProp.push_back(prop_suppressed);

        //! --------------
        //! under "Scope"
        //! --------------
        vecProp.push_back(prop_scopingMethod);
        vecProp.push_back(prop_scope);
        vecProp.push_back(prop_tags);

    }
        break;

        //! ----------------
        //! bolt pretension
        //! ----------------
    case SimulationNodeClass::nodeType_structuralAnalysisBoltPretension:
    {
        cout<<"nodeFactory::nodeFromScratch->____creating item \"Bolt pretension\"____"<<endl;
        name = "Bolt pretension";

        //! under "Scope"
        vecProp.push_back(prop_scopingMethod);
        vecProp.push_back(prop_scope);
        vecProp.push_back(prop_tags);

        //! under "Definition"
        vecProp.push_back(prop_suppressed);

        //! -----------------------------------------
        //! coordinate system defining the bolt load
        //! -----------------------------------------
        void *CSGlobal = addOptions.value<void*>();
        data.setValue(CSGlobal);

        Property prop_coordinateSystem("Coordinate system",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_coordinateSystem);

        data.setValue(Property::boltStatusDefinedBy_load);
        Property prop_defineBy("Bolt status",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_defineBy);

        data.setValue(0.0);
        Property prop_load("Load",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_load);

        data.setValue(0.0);
        Property prop_adjustment("Adjustment",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_adjustment);

        //! -------
        //! hidden
        //! -------
        QVector<double> refPoint{0,0,0};
        data.setValue(refPoint);
        Property prop_refPoint("Reference point",data,Property::PropertyGroup_Hidden);
        vecProp.push_back(prop_refPoint);
    }
        break;

        //! -------------
        //! remote point
        //! -------------
    case SimulationNodeClass::nodeType_remotePoint:
    {
        cout<<"nodeFactory::nodeFromScratch->____creating item \"Remote point\"____"<<endl;
        name = "Remote point";

        Property prop_Xcoordinate;
        prop_Xcoordinate.setName("X coordinate");
        prop_Xcoordinate.setGroup(Property::PropertyGroup_Scope);
        Property prop_Ycoordinate;
        prop_Ycoordinate.setName("Y coordinate");
        prop_Ycoordinate.setGroup(Property::PropertyGroup_Scope);
        Property prop_Zcoordinate;
        prop_Zcoordinate.setName("Z coordinate");
        prop_Zcoordinate.setGroup(Property::PropertyGroup_Scope);

        Property prop_XAbscoord;
        prop_XAbscoord.setName("X abs coordinate");
        prop_XAbscoord.setGroup(Property::PropertyGroup_Scope);
        Property prop_YAbscoord;
        prop_YAbscoord.setName("Y abs coordinate");
        prop_YAbscoord.setGroup(Property::PropertyGroup_Scope);
        Property prop_ZAbscoord;
        prop_ZAbscoord.setName("Z abs coordinate");
        prop_ZAbscoord.setGroup(Property::PropertyGroup_Scope);

        if(addOptions==QVariant())
        {
            cout<<"____remote point created without options (default remote point settings)____"<<endl;

            //! coordinates in the current system of reference (initially the global)
            data.setValue(0.0);

            //! coordinates in the system of reference
            prop_Xcoordinate.setData(data);
            prop_Ycoordinate.setData(data);
            prop_Zcoordinate.setData(data);

            //! absolute coordinates
            prop_XAbscoord.setData(data);
            prop_YAbscoord.setData(data);
            prop_ZAbscoord.setData(data);
        }
        else
        {
            double xcoord = addOptions.value<QVector<double>>().at(0);
            double ycoord = addOptions.value<QVector<double>>().at(1);
            double zcoord = addOptions.value<QVector<double>>().at(2);

            cerr<<"____remote point created with options (assigned location)____"<<endl;
            cerr<<"____("<<xcoord<<", "<<ycoord<<", "<<zcoord<<")____"<<endl;

            //! coordinates
            data.setValue(xcoord);
            prop_Xcoordinate.setData(data);
            data.setValue(ycoord);
            prop_Ycoordinate.setData(data);
            data.setValue(zcoord);
            prop_Zcoordinate.setData(data);

            //! ------------------------------------------------
            //! absolute coordinates (global coordinate system)
            //! ------------------------------------------------
            data.setValue(xcoord);
            prop_XAbscoord.setData(data);
            data.setValue(ycoord);
            prop_YAbscoord.setData(data);
            data.setValue(zcoord);
            prop_ZAbscoord.setData(data);
        }

        std::vector<GeometryTag> vecLoc;
        data.setValue(vecLoc);
        Property prop_location("Location",data,Property::PropertyGroup_Scope);

        //! under "Scope"
        vecProp.push_back(prop_scopingMethod);
        vecProp.push_back(prop_scope);
        vecProp.push_back(prop_tags);

        vecProp.push_back(prop_Xcoordinate);
        vecProp.push_back(prop_Ycoordinate);
        vecProp.push_back(prop_Zcoordinate);

        vecProp.push_back(prop_XAbscoord);
        vecProp.push_back(prop_YAbscoord);
        vecProp.push_back(prop_ZAbscoord);

        vecProp.push_back(prop_location);

        //! under "Definition"
        vecProp.push_back(prop_suppressed);

        //! ------------------------------------------------------
        //! type of coupling: "0" => kinematic "1" => distributed
        //! ------------------------------------------------------
        data.setValue(1);
        Property prop_toc("Coupling",data,Property::PropertyGroup_Advanced);

        //! --------------------------------------------------------
        //! DOFs selection: "0" => program controlled "1" => manual
        //! --------------------------------------------------------
        data.setValue(1);
        Property prop_DOFselection("DOFs selection",data,Property::PropertyGroup_Advanced);

        //! ---------------------------------------------
        //! DOFs selectors "0" => inactive "1" => active
        //! ---------------------------------------------
        data.setValue(1);
        Property prop_DOFX("X component ",data,Property::PropertyGroup_Advanced);
        Property prop_DOFY("Y component ",data,Property::PropertyGroup_Advanced);
        Property prop_DOFZ("Z component ",data,Property::PropertyGroup_Advanced);
        Property prop_ROTX("X rotation",data,Property::PropertyGroup_Advanced);
        Property prop_ROTY("Y rotation",data,Property::PropertyGroup_Advanced);
        Property prop_ROTZ("Z rotation",data,Property::PropertyGroup_Advanced);

        //! under advanced
        vecProp.push_back(prop_toc);
        vecProp.push_back(prop_DOFselection);
        vecProp.push_back(prop_DOFX);
        vecProp.push_back(prop_DOFY);
        vecProp.push_back(prop_DOFZ);
        vecProp.push_back(prop_ROTX);
        vecProp.push_back(prop_ROTY);
        vecProp.push_back(prop_ROTZ);

        //! ----------------------------------------
        //! add a sphere marker as a graphic object
        //! ----------------------------------------
        AIS_SphereMarker_handle_reg sphereMarker = markers::buildSphereMarker(gp_Pnt(0,0,0),10);
        data.setValue(sphereMarker);
        Property prop_sphereMarker("Graphic object",data,Property::PropertyGroup_GraphicObjects);
        vecProp.push_back(prop_sphereMarker);
    }
        break;

        //! -------------
        //! contact pair
        //! -------------
    case SimulationNodeClass::nodeType_connectionPair:
    {
        cout<<"nodeFactory::nodeFromScratch()->____creating connection node____"<<endl;

        QVariant data;
        std::vector<GeometryTag> vecLocMaster;
        std::vector<GeometryTag> vecLocSlave;

        //! -------------------
        //! suppression status
        //! -------------------
        vecProp.push_back(prop_suppressed);

        const connectionPairGenerationOption &options = addOptions.value<connectionPairGenerationOption>();
        bool isManual = options.manual;
        if(isManual==false)
        {
            cout<<"____creating an automatic contact____"<<endl;
            //! -----
            //! name
            //! -----
            name = options.connectionPairName;
            cout<<"____name: "<<name.toStdString()<<"____"<<endl;

            //! ---------------
            //! scoping method
            //! ---------------
            data.setValue(Property::ScopingMethod_GeometrySelection);
            Property prop_scopingMethod("Scoping method",data,Property::PropertyGroup_Scope);
            vecProp.push_back(prop_scopingMethod);

            //! -------
            //! master
            //! -------
            vecLocMaster = std::vector<GeometryTag>();
            data.setValue(vecLocMaster);
            Property prop_scopeMaster("Master",data,Property::PropertyGroup_Scope);
            Property prop_MasterTags("Tags master",data,Property::PropertyGroup_Scope);
            vecProp.push_back(prop_scopeMaster);
            vecProp.push_back(prop_MasterTags);

            //! ------
            //! slave
            //! ------
            vecLocSlave = std::vector<GeometryTag>();
            data.setValue(vecLocSlave);
            Property prop_scopeSlave("Slave",data,Property::PropertyGroup_Scope);
            Property prop_SlaveTags("Tags slave",data,Property::PropertyGroup_Scope);
            vecProp.push_back(prop_scopeSlave);
            vecProp.push_back(prop_SlaveTags);

        }
        else
        {
            //! -----
            //! name
            //! -----
            name = "Contact pair";

            //! ---------------
            //! scoping method
            //! ---------------
            data.setValue(Property::ScopingMethod_GeometrySelection);
            Property prop_scopingMethod("Scoping method",data,Property::PropertyGroup_Scope);
            vecProp.push_back(prop_scopingMethod);

            //! -------
            //! master
            //! -------
            vecLocMaster = vecLoc;
            data.setValue(vecLocMaster);
            Property prop_scopeMaster("Master",data,Property::PropertyGroup_Scope);
            Property prop_MasterTags("Tags master",data,Property::PropertyGroup_Scope);
            vecProp.push_back(prop_scopeMaster);
            vecProp.push_back(prop_MasterTags);

            //! ------
            //! slave
            //! ------
            vecLocSlave = std::vector<GeometryTag>();
            data.setValue(vecLocSlave);
            Property prop_scopeSlave("Slave",data,Property::PropertyGroup_Scope);
            Property prop_SlaveTags("Tags slave",data,Property::PropertyGroup_Scope);
            vecProp.push_back(prop_scopeSlave);
            vecProp.push_back(prop_SlaveTags);
        }

        //! -----------------------------------------------------------------------
        //! the contact formualation: contact pair initially treated with lagrangian
        //! -----------------------------------------------------------------------
        Property::contactFormulation theContactFormulation = Property::contactFormulation_lagrange;
        data.setValue(theContactFormulation);
        Property prop_connectionFormulation("Formulation",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_connectionFormulation);

        //! ---------------------------------------------------------------
        //! the contact type: contact pair initially created as frictional
        //! ---------------------------------------------------------------
        Property::contactType theContactType = Property::contactType_frictional;
        data.setValue(theContactType);
        Property prop_connectionType("Type",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_connectionType);

        //! ---------------------------------------------
        //! the contact behavior - symmetric by default
        //! ---------------------------------------------
        Property::contactBehavior theContactBehavior = Property::contactBehavior_symmetric;
        data.setValue(theContactBehavior);
        Property prop_contactBehavior("Behavior",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_contactBehavior);

        //! -------------------------------------
        //! friction coefficient: 0.1 by default
        //! -------------------------------------
        data.setValue(0.1);
        Property prop_frictionCoefficient("Friction coefficient",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_frictionCoefficient);

        //! --------------------------------------------
        //! small sliding: "0" => inactive "1" => active
        //! --------------------------------------------
        data.setValue(int(1));
        Property prop_smallSliding("Small sliding",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_smallSliding);

        //! --------------------------------------------
        //! adjust to touch: "0" => inactive "1" => active
        //! --------------------------------------------
        data.setValue(int(1));
        Property prop_adjust("Adjust to touch",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_adjust);

        //! ----------------------
        //! Overpressure function
        //! ----------------------
        Property::overpressureFunction theOPfunction = Property::overpressureFunction_linear;
        data.setValue(theOPfunction);
        Property prop_OPfunction("Overpressure",data,Property::PropertyGroup_Advanced);
        vecProp.push_back(prop_OPfunction);

        //! ----------------------------------------------
        //! coefficient K: if "0" => "Program controlled"
        //! ----------------------------------------------
        data.setValue(0.0);
        Property prop_K("K",data,Property::PropertyGroup_Advanced);
        vecProp.push_back(prop_K);

        //! --------------------------------------------
        //! sigma infinity: if "0" => "Program controlled"
        //! --------------------------------------------
        data.setValue(0.0);
        Property prop_sigmaInf("Sigma infinity",data,Property::PropertyGroup_Advanced);
        vecProp.push_back(prop_sigmaInf);

        //! ----------------------------------------------------
        //! thermal conductance: if "0" => "Program controlled"
        //! ----------------------------------------------------
        data.setValue(0.0);
        Property prop_thermalCond("Thermal conductance",data,Property::PropertyGroup_Advanced);
        vecProp.push_back(prop_thermalCond);

        //! --------------------------------------------
        //! constant C0: if "0" => "Program controlled"
        //! --------------------------------------------
        data.setValue(0.0);
        Property prop_C0("C0",data,Property::PropertyGroup_Advanced);
        vecProp.push_back(prop_C0);
        //! --------------------------------------------
        //! constant P0: if "0" => "Program controlled"
        //! --------------------------------------------
        data.setValue(0.0);
        Property prop_P0("P0",data,Property::PropertyGroup_Advanced);
        vecProp.push_back(prop_P0);
        //! ------------------------------------------------
        //! constant lamdba: if "0" => "Program controlled"
        //! ------------------------------------------------
        data.setValue(0.0);
        Property prop_lambda("Lambda",data,Property::PropertyGroup_Advanced);
        vecProp.push_back(prop_lambda);

        //! -----------------------------------------------------------------------------------------
        //! parent time tag: this automatically establishes a parenthood with the "Connection group"
        //! -----------------------------------------------------------------------------------------
        data.setValue(addOptions.value<connectionPairGenerationOption>().timeTag);
        Property prop_parentTimeTag("Parent time tag",data,Property::PropertyGroup_Identifier);

        vecProp.push_back(prop_parentTimeTag);
    }
        break;

        //! ------------------------------------------------------
        //! coordinate system - this can be edited by the user
        //! the global coordinate system has been already created
        //! within the geometryDataBase class
        //! ------------------------------------------------------
    case SimulationNodeClass::nodeType_coordinateSystem:
    {
        name = "Coordinate system";

        //! -------------------
        //! "Define by" switch
        //! -------------------
        QVariant data;
        data.setValue(Property::defineBy_geometrySelection);

        //! warning: the property key name "Define by " includes space
        Property prop_defineBy("Define by ",data,Property::PropertyGroup_Origin);
        vecProp.push_back(prop_defineBy);

        //! -----------------
        //! directional data
        //! -----------------
        QVector<double> dirX{1,0,0}, dirY{0,1,0}, dirZ{0,0,1};
        QVariant datax,datay,dataz;

        datax.setValue(dirX);
        datay.setValue(dirY);
        dataz.setValue(dirZ);

        Property prop_X_axisData("X axis data",datax,Property::PropertyGroup_DirectionalData);
        Property prop_Y_axisData("Y axis data",datay,Property::PropertyGroup_DirectionalData);
        Property prop_Z_axisData("Z axis data",dataz,Property::PropertyGroup_DirectionalData);
        vecProp.push_back(prop_X_axisData);
        vecProp.push_back(prop_Y_axisData);
        vecProp.push_back(prop_Z_axisData);

        //! ----------------
        //! the base origin
        //! ----------------
        QVector<double> baseOriginCoords {0,0,0};
        data.setValue(baseOriginCoords);
        Property prop_baseOrigin("Base origin",data,Property::PropertyGroup_Transformations);
        vecProp.push_back(prop_baseOrigin);

        //! --------------------------
        //! the base directional data
        //! --------------------------
        QVector<QVector<double>> baseDirectionalData;

        QVector<double> xAxis{1,0,0}, yAxis{0,1,0}, zAxis{0,0,1};
        baseDirectionalData.push_back(xAxis); baseDirectionalData.push_back(yAxis); baseDirectionalData.push_back(zAxis);
        data.setValue(baseDirectionalData);
        Property prop_baseDirectionalData("Base directional data",data,Property::PropertyGroup_Transformations);
        vecProp.push_back(prop_baseDirectionalData);

        //! ----------------------------------------------
        //! here "Geometry" and "Tags" are under "Origin"
        //! ----------------------------------------------
        data.setValue(vecLoc);
        Property prop_geometry("Geometry",data,Property::PropertyGroup_Origin);
        Property prop_Tags("Tags",data,Property::PropertyGroup_Origin);
        vecProp.push_back(prop_geometry);
        vecProp.push_back(prop_Tags);

        //! -----------------------------------
        //! "Origin X", "Origin Y", "Origin Z"
        //! -----------------------------------
        Property prop_Xorigin("Origin X",0,Property::PropertyGroup_Origin);
        Property prop_Yorigin("Origin Y",0,Property::PropertyGroup_Origin);
        Property prop_Zorigin("Origin Z",0,Property::PropertyGroup_Origin);
        vecProp.push_back(prop_Xorigin);
        vecProp.push_back(prop_Yorigin);
        vecProp.push_back(prop_Zorigin);

        //! --------------------------------------
        //! add the trihedron as a graphic object
        //! --------------------------------------
        AIS_Trihedron_handle_reg theTrihedron = markers::buildTrihedron(baseOriginCoords,baseDirectionalData,20);
        data.setValue(theTrihedron);
        Property prop_trihedron("Graphic object",data,Property::PropertyGroup_GraphicObjects);
        vecProp.push_back(prop_trihedron);
    }
        break;

        //! ----------------------------------
        //! "Displacement" boundary condition
        //! ----------------------------------
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Displacement:
    {
        name = "Displacement";

        Property::defineBy theDefineby = Property::defineBy_components;
        data.setValue(theDefineby);
        Property prop_defineBy("Define by",data,Property::PropertyGroup_Definition);

        vecProp.push_back(prop_scopingMethod);
        vecProp.push_back(prop_scope);
        vecProp.push_back(prop_tags);
        vecProp.push_back(prop_suppressed);
        vecProp.push_back(prop_defineBy);

        //! ----------------------------------------------------------------------------------
        //! the "Displacement" is created from scratch setting the three components to "Free"
        //! ----------------------------------------------------------------------------------
        Property::loadDefinition load_componentX = Property::loadDefinition_free;
        data.setValue(load_componentX);
        Property prop_componentX("X component",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_componentX);

        Property::loadDefinition load_componentY = Property::loadDefinition_free;
        data.setValue(load_componentY);
        Property prop_componentY("Y component",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_componentY);

        Property::loadDefinition load_componentZ = Property::loadDefinition_free;
        data.setValue(load_componentZ);
        Property prop_componentZ("Z component",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_componentZ);        
    }
        break;

        //! ------------------------------------------------------------
        //! "Remote rotation" the rotation defined using a remote point
        //! ------------------------------------------------------------
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteRotation:
    {
        name = "Remote rotation";

        Property::defineBy theDefineby = Property::defineBy_components;
        data.setValue(theDefineby);
        Property prop_defineBy("Define by",data,Property::PropertyGroup_Definition);

        vecProp.push_back(prop_scopingMethod);
        vecProp.push_back(prop_scope);
        vecProp.push_back(prop_tags);
        vecProp.push_back(prop_suppressed);
        vecProp.push_back(prop_defineBy);

        //! -----------------------------------------
        //! "Remote displacement" boundary condition
        //! -----------------------------------------
        Property::loadDefinition load_rotationX = Property::loadDefinition_free;
        data.setValue(load_rotationX);
        Property prop_rotationX("X component",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_rotationX);

        Property::loadDefinition load_rotationY = Property::loadDefinition_free;
        data.setValue(load_rotationY);
        Property prop_rotationY("Y component",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_rotationY);

        Property::loadDefinition load_rotationZ = Property::loadDefinition_free;
        data.setValue(load_rotationZ);
        Property prop_rotationZ("Z component",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_rotationZ);

        //! ---------------------------------------------------------
        //! Property "Coupling": "0" => kinematic "1" => distributed
        //! ---------------------------------------------------------
        data.setValue(1);
        Property prop_Coupling("Coupling",data,Property::PropertyGroup_Advanced);

        //! -------------------------------------------------------
        //! DOFs selection "0" => program controlled "1" => manual
        //! -------------------------------------------------------
        data.setValue(1);
        Property prop_DOFselection("DOFs selection",data,Property::PropertyGroup_Advanced);

        //! ---------------------------------------------
        //! DOFs selectors "0" => inactive "1" => active
        //! ---------------------------------------------
        data.setValue(1);
        Property prop_DOFX("X component ",data,Property::PropertyGroup_Advanced);
        Property prop_DOFY("Y component ",data,Property::PropertyGroup_Advanced);
        Property prop_DOFZ("Z component ",data,Property::PropertyGroup_Advanced);
        Property prop_ROTX("X rotation",data,Property::PropertyGroup_Advanced);
        Property prop_ROTY("Y rotation",data,Property::PropertyGroup_Advanced);
        Property prop_ROTZ("Z rotation",data,Property::PropertyGroup_Advanced);

        //! -----------------
        //! under "Advanced"
        //! -----------------
        vecProp.push_back(prop_Coupling);
        vecProp.push_back(prop_DOFselection);
        vecProp.push_back(prop_DOFX);
        vecProp.push_back(prop_DOFY);
        vecProp.push_back(prop_DOFZ);
        vecProp.push_back(prop_ROTX);
        vecProp.push_back(prop_ROTY);
        vecProp.push_back(prop_ROTZ);

        //! ------------------------------------------
        //! "Hidden properties"
        //! coordinates of the reference node (0,0,0)
        //! ------------------------------------------
        QVector<double> referencePointCoords;
        referencePointCoords.clear();
        data.setValue(referencePointCoords);
        Property prop_referencePoint("Reference point",data,Property::PropertyGroup_Hidden);
        vecProp.push_back(prop_referencePoint);
    }
        break;

        //! -----------------------------------------
        //! "Remote displacement" boundary condition
        //! -----------------------------------------
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteDisplacement:
    {
        name = "Remote displacement";

        Property::defineBy theDefineby = Property::defineBy_components;
        data.setValue(theDefineby);
        Property prop_defineBy("Define by",data,Property::PropertyGroup_Definition);

        vecProp.push_back(prop_scopingMethod);
        vecProp.push_back(prop_scope);
        vecProp.push_back(prop_tags);
        vecProp.push_back(prop_suppressed);
        vecProp.push_back(prop_defineBy);

        //! -----------------------------------------------------------------------------------------
        //! the "Remote displacement" is created from scratch setting the three components to "Free"
        //! -----------------------------------------------------------------------------------------
        Property::loadDefinition load_componentX = Property::loadDefinition_free;
        data.setValue(load_componentX);
        Property prop_componentX("X component",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_componentX);

        Property::loadDefinition load_componentY = Property::loadDefinition_free;
        data.setValue(load_componentY);
        Property prop_componentY("Y component",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_componentY);

        Property::loadDefinition load_componentZ = Property::loadDefinition_free;
        data.setValue(load_componentZ);
        Property prop_componentZ("Z component",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_componentZ);

        //! ---------------------------------------------------------
        //! Property "Coupling": "0" => kinematic "1" => distributed
        //! ---------------------------------------------------------
        data.setValue(1);
        Property prop_Coupling("Coupling",data,Property::PropertyGroup_Advanced);

        //! -------------------------------------------------------
        //! DOFs selection "0" => program controlled "1" => manual
        //! -------------------------------------------------------
        data.setValue(1);
        Property prop_DOFselection("DOFs selection",data,Property::PropertyGroup_Advanced);

        //! ---------------------------------------------
        //! DOFs selectors "0" => inactive "1" => active
        //! ---------------------------------------------
        data.setValue(1);
        Property prop_DOFX("X component ",data,Property::PropertyGroup_Advanced);
        Property prop_DOFY("Y component ",data,Property::PropertyGroup_Advanced);
        Property prop_DOFZ("Z component ",data,Property::PropertyGroup_Advanced);
        Property prop_ROTX("X rotation",data,Property::PropertyGroup_Advanced);
        Property prop_ROTY("Y rotation",data,Property::PropertyGroup_Advanced);
        Property prop_ROTZ("Z rotation",data,Property::PropertyGroup_Advanced);

        //! -----------------
        //! under "Advanced"
        //! -----------------
        vecProp.push_back(prop_Coupling);
        vecProp.push_back(prop_DOFselection);
        vecProp.push_back(prop_DOFX);
        vecProp.push_back(prop_DOFY);
        vecProp.push_back(prop_DOFZ);
        vecProp.push_back(prop_ROTX);
        vecProp.push_back(prop_ROTY);
        vecProp.push_back(prop_ROTZ);

        //! ------------------------------------------
        //! "Hidden properties"
        //! coordinates of the reference node (0,0,0)
        //! ------------------------------------------
        QVector<double> referencePointCoords;
        referencePointCoords.clear();
        //referencePointCoords.push_back(0);
        //referencePointCoords.push_back(0);
        //referencePointCoords.push_back(0);
        data.setValue(referencePointCoords);
        Property prop_referencePoint("Reference point",data,Property::PropertyGroup_Hidden);
        vecProp.push_back(prop_referencePoint);
    }
        break;

        //! ------------------------------------------------
        //! Force/Remote force/Moment/Pressure/Acceleration
        //! ------------------------------------------------
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteForce:
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Force:
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Moment:
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Acceleration:
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RotationalVelocity:
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Pressure:
    case SimulationNodeClass::nodeType_structuralAnalysisThermalCondition:
    {
        //! -----------------
        //! generate a name
        //! -----------------
        switch(type)
        {
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Force: name = "Force"; break;
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteForce: name = "Remote force"; break;
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Moment: name = "Moment"; break;
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Acceleration: name = "Acceleration"; break;
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RotationalVelocity: name = "Rotational velocity"; break;
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Pressure: name = "Pressure"; break;
        case SimulationNodeClass::nodeType_structuralAnalysisThermalCondition: name = "Thermal condition"; break;
        }

        //! -------------------------------------------------------
        //! Add the property "Define by" for Force/Moment/Pressure
        //! "Pressure" has only the option "Normal to"
        //! "Thermal condition" has not the property "Define by"
        //! -------------------------------------------------------
        if(type == SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Force ||
                type == SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteForce ||
                type == SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Moment ||
                type == SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Acceleration ||
                type == SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RotationalVelocity)
        {
            Property::defineBy theDefineby = Property::defineBy_vector;
            data.setValue(theDefineby);
            Property prop_defineBy("Define by",data,Property::PropertyGroup_Definition);
            vecProp.push_back(prop_defineBy);
        }
        else if(type == SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Pressure)
        {
            Property::defineBy theDefineby = Property::defineBy_normal;
            data.setValue(theDefineby);
            Property prop_defineBy("Define by",data,Property::PropertyGroup_Definition);
            vecProp.push_back(prop_defineBy);

            data.setValue(Property::loadDefinition_constant);
            Property prop_loadMagnitude("Magnitude",data,Property::PropertyGroup_Definition);
            data.setValue(prop_loadMagnitude);
            vecProp.push_back(prop_loadMagnitude);
        }

        //! -----------------------------------
        //! advanced properties for "Coupling"
        //! -----------------------------------
        if(type == SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Force ||
                type == SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteForce ||
                type == SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Moment ||
                type == SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Moment)
        {
            //! ---------------------------------------------------------
            //! Property "Coupling": "0" => kinematic "1" => distributed
            //! ---------------------------------------------------------
            data.setValue(1);
            Property prop_Coupling("Coupling",data,Property::PropertyGroup_Advanced);

            //! -------------------------------------------------------
            //! DOFs selection "0" => program controlled "1" => manual
            //! -------------------------------------------------------
            data.setValue(1);
            Property prop_DOFselection("DOFs selection",data,Property::PropertyGroup_Advanced);

            //! ---------------------------------------------
            //! DOFs selectors "0" => inactive "1" => active
            //! ---------------------------------------------
            data.setValue(1);
            Property prop_DOFX("X component ",data,Property::PropertyGroup_Advanced);
            Property prop_DOFY("Y component ",data,Property::PropertyGroup_Advanced);
            Property prop_DOFZ("Z component ",data,Property::PropertyGroup_Advanced);
            Property prop_ROTX("X rotation",data,Property::PropertyGroup_Advanced);
            Property prop_ROTY("Y rotation",data,Property::PropertyGroup_Advanced);
            Property prop_ROTZ("Z rotation",data,Property::PropertyGroup_Advanced);

            //! under "Advanced"
            vecProp.push_back(prop_Coupling);
            vecProp.push_back(prop_DOFselection);
            vecProp.push_back(prop_DOFX);
            vecProp.push_back(prop_DOFY);
            vecProp.push_back(prop_DOFZ);
            vecProp.push_back(prop_ROTX);
            vecProp.push_back(prop_ROTY);
            vecProp.push_back(prop_ROTZ);

            //! ------------------------------------------
            //! "Hidden properties"
            //! coordinates of the reference node (0,0,0)
            //! ------------------------------------------
            QVector<double> referencePointCoords;
            referencePointCoords.clear();
            //referencePointCoords.push_back(0);
            //referencePointCoords.push_back(0);
            //referencePointCoords.push_back(0);
            data.setValue(referencePointCoords);
            Property prop_referencePoint("Reference point",data,Property::PropertyGroup_Hidden);
            vecProp.push_back(prop_referencePoint);
        }

        //! -------------------------------
        //! under "Scope" and "Definition"
        //! -------------------------------
        vecProp.push_back(prop_suppressed);
        vecProp.push_back(prop_scopingMethod);
        vecProp.push_back(prop_scope);
        vecProp.push_back(prop_tags);
    }
        break;

    //! ---------------------------------------
    //! thermal boundary condition convection
    //! ---------------------------------------
    case SimulationNodeClass::nodeType_thermalAnalysisConvection:
    {
        name = "Convection";

        //! -------------------------------
        //! under "Scope" and "Definition"
        //! -------------------------------
        vecProp.push_back(prop_suppressed);
        vecProp.push_back(prop_scopingMethod);
        vecProp.push_back(prop_scope);
        vecProp.push_back(prop_tags);

        //! -----------------------
        //! convection coefficient
        //! -----------------------
        data.setValue(Property::loadDefinition_constant);
        Property prop_loadConvectionCoefficient("Film coefficient",data,Property::PropertyGroup_Definition);
        data.setValue(prop_loadConvectionCoefficient);
        vecProp.push_back(prop_loadConvectionCoefficient);

        //! ----------------------
        //! reference temperature
        //! ----------------------
        data.setValue(Property::loadDefinition_constant);
        Property prop_loadReferenceTemperature("Reference temperature",data,Property::PropertyGroup_Definition);
        data.setValue(prop_loadReferenceTemperature);
        vecProp.push_back(prop_loadReferenceTemperature);
    }
        break;

    //! ---------------
    //! adiabatic wall
    //! ---------------
    case SimulationNodeClass::nodeType_thermalAnalysisAdiabaticWall:
    {
        name = "Adiabatic wall";

        //! -------------------------------
        //! under "Scope" and "Definition"
        //! -------------------------------
        vecProp.push_back(prop_suppressed);
        vecProp.push_back(prop_scopingMethod);
        vecProp.push_back(prop_scope);
        vecProp.push_back(prop_tags);
    }
        break;

    //! -------------------------------------------------
    //! thermal flux and thermal flow boundary condition
    //! -------------------------------------------------
    case SimulationNodeClass::nodeType_thermalAnalysisTemperature:
    case SimulationNodeClass::nodeType_thermalAnalysisThermalFlux:
    case SimulationNodeClass::nodeType_thermalAnalysisThermalFlow:
    case SimulationNodeClass::nodeType_thermalAnalysisThermalPower:
    {
        switch (type)
        {
        case SimulationNodeClass::nodeType_thermalAnalysisTemperature: name ="Temperature"; break;
        case SimulationNodeClass::nodeType_thermalAnalysisThermalFlux: name ="Thermal flux"; break;
        case SimulationNodeClass::nodeType_thermalAnalysisThermalFlow: name ="Thermal flow"; break;
        case SimulationNodeClass::nodeType_thermalAnalysisThermalPower: name = "Heat generation"; break;
        }

        //! -------------------------------
        //! under "Scope" and "Definition"
        //! -------------------------------
        vecProp.push_back(prop_suppressed);
        vecProp.push_back(prop_scopingMethod);
        vecProp.push_back(prop_scope);
        vecProp.push_back(prop_tags);

        //! -----------------------
        //! convection coefficient
        //! -----------------------
        data.setValue(Property::loadDefinition_constant);
        Property prop_loadMagnitude("Magnitude",data,Property::PropertyGroup_Definition);
        data.setValue(prop_loadMagnitude);
        vecProp.push_back(prop_loadMagnitude);
    }
        break;

        //! ------------
        //! Repair tool
        //! ------------
    case SimulationNodeClass::nodeType_repairTool:
    {
        name = "Repair tool";
        Property prop_fixSmallEdges("Fix small edges",false,Property::PropertyGroup_Definition);
        Property prop_sewFaces("Sew small faces",false,Property::PropertyGroup_Definition);

        //! under "Definition"
        vecProp.push_back(prop_scopingMethod);
        vecProp.push_back(prop_scope);
        vecProp.push_back(prop_tags);
        vecProp.push_back(prop_suppressed);
    }
        break;

        //! ------------
        //! Mesh method
        //! ------------
    case SimulationNodeClass::nodeType_meshMethod:
    {
        name = "Method";
        vecProp.push_back(prop_suppressed);
        vecProp.push_back(prop_scopingMethod);
        vecProp.push_back(prop_scope);
        vecProp.push_back(prop_tags);

        //! -------------------------------------------
        //! default node with BRep on - no defeaturing
        //! -------------------------------------------
        data.setValue(Property::meshEngine2D_Netgen);
        Property prop_surfaceMesher("Surface mesher",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_surfaceMesher);

        data.setValue(Property::meshEngine3D_Netgen);
        Property prop_volumeMesher("Volume mesher",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_volumeMesher);

        data.setValue(Property::meshOrder_First);
        Property prop_meshOrder("Mesh order",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_meshOrder);

        bool useBRep = true;
        data.setValue(useBRep);
        Property prop_useBRep("Patch conforming",data,Property::PropertyGroup_Method);
        vecProp.push_back(prop_useBRep);

        /*
        bool defeaturing = false;
        data.setValue(defeaturing);
        Property prop_defeaturing("Defeaturing",data,Property::PropertyGroup_Defeaturing);
        vecProp.push_back(prop_defeaturing);
        */

        //! ----------------------------------------------------
        //! default node with BRep off - defeatuging by "Level"
        //! ----------------------------------------------------
        /*
        bool useBRep = false;
        data.setValue(useBRep);
        Property prop_useBRep("Patch conforming",data,Property::PropertyGroup_Method);
        vecProp.pop_back(prop_useBRep);

        data.setValue(Property::meshEngine2D_OCC_STL);
        Property prop_surfaceDiscretizer("Tessellator",data,Property::PropertyGroup_Advanced);
        vecProp.push_back(prop_surfaceDiscretizer);

        double angularDeflection = 0.01;
        data.setValue(angularDeflection);
        Property prop_angularDeflection("Angular deflection",data,Property::PropertyGroup_Advanced);
        vecProp.push_back(prop_angularDeflection);

        double linearDeflection = 0.1;
        Property prop_linearDeflection("Linear deflection",data,Property::PropertyGroup_Advanced);
        vecProp.push_back(prop_linearDeflection);

        //! add healing ...

        bool defeaturing = false;
        data.setValue(defeaturing);
        Property prop_defeaturing("Defeaturing",data,Property::propertyGroup_Defeaturing);
        vecProp.push_back(prop_defeaturing);

        int by = 0;
        data.setValue(by);
        Property prop_by("By",data,Property::propertyGroup_Defeaturing);
        vecProp.push_back(prop_by);

        double level = 0.95;
        data.setValue(level);
        Property prop_level("Level",data,Property::PropertyGroup_Defeaturing);
        vecProp.pop_back(prop_level);
        */
    }
        break;

    //! ----------
    //! Mesh type
    //! ----------
    case SimulationNodeClass::nodeType_meshMeshType:
    {
        name = "Mesh type";
        vecProp.push_back(prop_suppressed);
        vecProp.push_back(prop_scopingMethod);
        vecProp.push_back(prop_scope);
        vecProp.push_back(prop_tags);

        //! -----------------------
        //! 0 => fully tetrahedral
        //! 1 => fully hexahedral
        //! 2 => hexa-dominant
        //! -----------------------
        data.setValue(0);
        Property prop_meshType("Mesh type",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_meshType);
    }
        break;

    //! ------------
    //! Mesh metric
    //! ------------
    case SimulationNodeClass::nodeType_meshMeshMetric:
    {
        name = "Mesh metric";
        vecProp.push_back(prop_scopingMethod);
        vecProp.push_back(prop_scope);
        vecProp.push_back(prop_tags);

        //! --------------------------------------------
        //! metric type - currently implemented for tet
        //! 0 => "Modified Mean Ratio (MMR)"
        //! 1 => "Modified Condition Number (MCN)"
        //! 2 => "Modified volume-length (MIVL)"
        //! --------------------------------------------
        data.setValue(0);
        Property prop_meshMetric("Metric type",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_meshMetric);

        histogramData hisData;
        data.setValue(hisData);
        Property prop_meshMetricData("Metric data",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_meshMetricData);
    }
        break;

    //! ------------
    //! Body sizing
    //! ------------
    case SimulationNodeClass::nodeType_meshBodyMeshControl:
    {
        name = "Body sizing";
        vecProp.push_back(prop_suppressed);
        vecProp.push_back(prop_scopingMethod);
        vecProp.push_back(prop_scope);
        vecProp.push_back(prop_tags);

        data.setValue(0.1);       //! some default values
        Property prop_minElSize("Min element size",data,Property::PropertyGroup_Definition);
        data.setValue(100);     //! some default values
        Property prop_maxElSize("Max element size",data,Property::PropertyGroup_Definition);
        data.setValue(0.3);     //! some default values
        Property prop_grading("Grading",data,Property::PropertyGroup_Definition);

        vecProp.push_back(prop_minElSize);
        vecProp.push_back(prop_maxElSize);
        vecProp.push_back(prop_grading);

        if(addOptions.isValid())
        {
            if(addOptions.canConvert<ListOfShape>()) data.setValue(addOptions);
            Property prop_graphObject("Colored shape",data,Property::PropertyGroup_GraphicObjects);
            vecProp.push_back(prop_graphObject);
        }
    }
        break;

        //! ------------
        //! Face sizing
        //! ------------
    case SimulationNodeClass::nodeType_meshFaceSize:
    {
        name = "Face sizing";
        vecProp.push_back(prop_suppressed);
        vecProp.push_back(prop_scopingMethod);
        vecProp.push_back(prop_scope);
        vecProp.push_back(prop_tags);

        data.setValue(10);       //! a default value
        Property prop_faceMeshSizing("Face sizing",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_faceMeshSizing);
    }
        break;

        //! ------------
        //! Edge sizing
        //! ------------
    case SimulationNodeClass::nodeType_meshEdgeSize:
    {
        name = "Edge sizing";
        vecProp.push_back(prop_suppressed);
        vecProp.push_back(prop_scopingMethod);
        vecProp.push_back(prop_scope);
        vecProp.push_back(prop_tags);

        //! --------------------------------------------------------------
        //! edge sizing type: 0 => sizing length 1 => number of divisions
        //! --------------------------------------------------------------
        data.setValue(0);
        Property prop_edgeSizingType("Sizing type",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_edgeSizingType);

        //! -----------------------------------------------------------------------
        //! a default sizing value: put both the mutually exclusive mesh controls.
        //! Only one will appear on the GUI, according to the "Sizing type"
        //! -----------------------------------------------------------------------
        data.setValue(double(10));
        Property prop_edgeMeshSizing("Element size",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_edgeMeshSizing);

        data.setValue(int(5));
        Property prop_NbDivisions("Number of divisions",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_NbDivisions);
    }
        break;

        //! --------------
        //! Vertex sizing
        //! --------------
    case SimulationNodeClass::nodeType_meshVertexSize:
    {
        name = "Vertex sizing";
        vecProp.push_back(prop_suppressed);
        vecProp.push_back(prop_scopingMethod);
        vecProp.push_back(prop_scope);
        vecProp.push_back(prop_tags);

        //! --------
        //! pinball
        //! --------
        data.setValue(10);
        Property prop_edgeSizingType("Pinball",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_edgeSizingType);

        //! -----------------------
        //! a default sizing value
        //! -----------------------
        data.setValue(10);
        Property prop_edgeMeshSizing("Element size",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_edgeMeshSizing);
    }
        break;

        //! ----------------
        //! Prismatic layer
        //! ----------------
    case SimulationNodeClass::nodeType_meshPrismaticLayer:
    {
        name = "Prismatic layer";
        vecProp.push_back(prop_suppressed);
        vecProp.push_back(prop_scopingMethod);
        vecProp.push_back(prop_scope);
        vecProp.push_back(prop_tags);

        //! ------------------------------------------------------------------
        //! generation of prismatic layer through
        //! (first layer thickness, total number of layers, Expansion ratio) => 0
        //! (total thickness, total number of layers, Expansion ratio) => 1
        //! ------------------------------------------------------------------

        //! --------------------------
        //! "Boundary scoping method"
        //! --------------------------
        data.setValue(Property::ScopingMethod_GeometrySelection);
        Property prop_boundaryScopingMethod("Boundary scoping method",data,Property::PropertyGroup_Definition);

        //! ------------------
        //! scope: "Boundary"
        //! ------------------
        ListOfShape boundaryScope;
        std::vector<GeometryTag> vecLoc = TopologyTools::generateLocationPairs(mDB,boundaryScope);
        data.setValue(vecLoc);

        Property prop_boundaryScope("Boundary",data,Property::PropertyGroup_Definition);

        //! --------------
        //! scope: "Tags"
        //! --------------
        Property prop_boundaryTags("Boundary tags",data,Property::PropertyGroup_Definition);

        int options = 0;
        data.setValue(options);
        Property prop_options("Options",data,Property::PropertyGroup_Definition);

        int numberOfLayers = 3;
        data.setValue(numberOfLayers);
        Property prop_NbLayers("Number of layers",data,Property::PropertyGroup_Definition);

        double firstLayerHeight = 1.0;
        data.setValue(firstLayerHeight);
        Property prop_firstLayerHeight("First layer height",data,Property::PropertyGroup_Definition);

        double growthRate = 1.2;
        data.setValue(growthRate);
        Property prop_growthRate("Expansion ratio",data,Property::PropertyGroup_Definition);

        //! ----------------------------------------
        //! algorithm for prismatic mesh generation
        //! 0 => "Pre" 1 => "Post"
        //! ----------------------------------------
        int algo = 0;
        data.setValue(algo);
        Property prop_generationAlgorithm("Algorithm",data,Property::PropertyGroup_Definition);

        //! -------------------
        //! Boundary mesh type
        //! 0 => "Hybrid"
        //! 1 => "Tetrahedral"
        //! -------------------
        int boundaryMeshType = 0;
        data.setValue(boundaryMeshType);
        Property prop_boundaryMeshType("Boundary mesh type",data,Property::PropertyGroup_Definition);

        //! ---------------
        //! new parameters
        //! ---------------
        double curvatureSensitivityForShrink = 100;
        data.setValue(curvatureSensitivityForShrink);
        Property prop_curvatureSensitivityForShrink("Curvature sensitivity",data,Property::PropertyGroup_Advanced);

        int NbGuidingVectorSmoothingSteps = 50;
        data.setValue(NbGuidingVectorSmoothingSteps);
        Property prop_NbGuidingVectorSmoothingSteps("Guiding vectors smoothing steps",data,Property::PropertyGroup_Smoothing);

        int NbLayerThicknessSmoothingSteps = 50;
        data.setValue(NbLayerThicknessSmoothingSteps);
        Property prop_NbLayerThicknessSmoothingSteps("Thickness smoothing steps",data,Property::PropertyGroup_Smoothing);

        double curvatureSensitivityForGuidingVectorsSmoothing = 5;
        data.setValue(curvatureSensitivityForGuidingVectorsSmoothing);
        Property prop_curvatureSensitivityForGuidingVectorsSmoothing("Guiding vector smoothing - curvature sensitivity",data,Property::PropertyGroup_Smoothing);

        double curvatureSensitivityForThicknessSmoothing = 50;
        data.setValue(curvatureSensitivityForThicknessSmoothing);
        Property prop_curvatureSensitivityForThickessSmoothing("Thickness smoothing - curvature sensitivity",data,Property::PropertyGroup_Smoothing);
        //! ----------------------
        //! end of new parameters
        //! ----------------------

        bool lockBoundary = true;
        data.setValue(lockBoundary);
        Property prop_lockBoundary("Lock boundary",data,Property::PropertyGroup_Advanced);

        bool checkSelfIntersections = true;
        data.setValue(checkSelfIntersections);
        Property prop_checkSelfIntersection("Check self intersections",data,Property::PropertyGroup_Advanced);

        bool checkMutualIntersections = true;
        data.setValue(checkMutualIntersections);
        Property prop_checkMutualIntersections("Check mutual intersections",data,Property::PropertyGroup_Advanced);

        vecProp.push_back(prop_boundaryScopingMethod);
        vecProp.push_back(prop_boundaryScope);
        vecProp.push_back(prop_boundaryTags);
        vecProp.push_back(prop_options);
        vecProp.push_back(prop_NbLayers);
        vecProp.push_back(prop_firstLayerHeight);
        vecProp.push_back(prop_growthRate);
        vecProp.push_back(prop_generationAlgorithm);
        vecProp.push_back(prop_boundaryMeshType);

        //! ---------------
        //! new parameters
        //! ---------------
        vecProp.push_back(prop_curvatureSensitivityForShrink);
        vecProp.push_back(prop_NbGuidingVectorSmoothingSteps);
        vecProp.push_back(prop_NbLayerThicknessSmoothingSteps);
        vecProp.push_back(prop_curvatureSensitivityForGuidingVectorsSmoothing);
        vecProp.push_back(prop_curvatureSensitivityForThickessSmoothing);
        //! ----------------------
        //! end of new parameters
        //! ----------------------

        vecProp.push_back(prop_lockBoundary);
        vecProp.push_back(prop_checkSelfIntersection);
        vecProp.push_back(prop_checkMutualIntersections);
    }
        break;

        //! ----------------
        //! Named selection
        //! ----------------
    case SimulationNodeClass::nodeType_namedSelectionGeometry:
    {
        name = "Named selection";
        vecProp.push_back(prop_scopingMethod);
        vecProp.push_back(prop_scope);
        vecProp.push_back(prop_tags);
    }
        break;

        //! -----------------------
        //! Mesh element selection
        //! -----------------------
    case SimulationNodeClass::nodeType_namedSelectionElement:
    {
        name = "Mesh element selection";
        //vecProp.push_back(prop_scopingMethod);
        vecProp.push_back(prop_scope);
        vecProp.push_back(prop_tags);

        //! --------------------------------------
        //! "Selection method"
        //! "0" elements numbers list
        //! "1" GUI selection - graphical picking
        //! --------------------------------------
        data.setValue(0);
        Property prop_selectionMethod("Selection method",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_selectionMethod);

        //! -----------------
        //! the element list
        //! -----------------
        QString elementList("");
        data.setValue(elementList);
        Property prop_elementList("Element list",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_elementList);
    }
        break;

        //! ------------
        //! Temperature
        //! ------------
    case SimulationNodeClass::nodeType_solutionStructuralTemperature:
    {
        name = "Temperature";
        //! ------------------------------------------------------------------
        //! redefine tags and scope:
        //! initially the scope is "All bodies" for the structural diagnostic
        //! ------------------------------------------------------------------
        int NbBodies = mDB->bodyMap.size();
        ListOfShape scope;
        for(int i=1; i<= NbBodies; i++)
        {
            TopoDS_Solid aSolid = TopoDS::Solid(mDB->bodyMap.value(i));
            scope.Append(aSolid);
        }
        std::vector<GeometryTag> vecLocAllBodies = TopologyTools::generateLocationPairs(mDB, scope);

        data.setValue(vecLocAllBodies);
        Property prop_scopeAllBodies("Geometry",data,Property::PropertyGroup_Scope);
        Property prop_tagsAllBodies("Tags",data,Property::PropertyGroup_Scope);

        //! under scope
        vecProp.push_back(prop_scopingMethod);
        vecProp.push_back(prop_scope);
        vecProp.push_back(prop_tags);

        int dummyType = 0;
        data.setValue(dummyType);
        Property prop_type("Type ",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_type);

        //! ---------------------------------------------------------------------------
        //! under definition
        //! "By" selector: it provides options "Time" and "Set"
        //! Here "By"="Time", so the "Display time" property is created
        //! If "By" is changed into "Set", the DetailViewer replace the "Display time"
        //! with "Set Number"
        //! 0 => "Time"
        //! 1 => "Set number"
        //! ---------------------------------------------------------------------------
        data.setValue(0.0);
        Property prop_By("By",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_By);
        Property prop_displayTime("Display time",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_displayTime);

        //! mode
        data.setValue(0);
        Property prop_modeNumber("Mode number",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_modeNumber);

        //! ---------------------------------------------------
        //! type of scale: "0" => "Autoscale"; "1" => "Custom"
        //! ---------------------------------------------------
        int typeOfScale = 0;
        data.setValue(typeOfScale);
        Property prop_scaleType("Scale type",data,Property::PropertyGroup_ColorBox);
        vecProp.push_back(prop_scaleType);

        //! ------------------------------------------------
        //! Number of intervals: "0" => "Program controlled
        //! ------------------------------------------------
        int NbIntervals = INITIAL_NUMBER_OF_COLORBOX_LEVELS;
        data.setValue(NbIntervals);
        Property prop_NbIntervals("# intervals",data,Property::PropertyGroup_ColorBox);
        vecProp.push_back(prop_NbIntervals);
    }
        break;

        //! --------------------------
        //! Equivalent plastic strain
        //! --------------------------
    case SimulationNodeClass::nodeType_solutionStructuralEquivalentPlasticStrain:
    {
        //! read options
        int typeOfEPS = addOptions.toInt();

        //! prepare the name
        switch(typeOfEPS)
        {
        case 0: name = "Equivalent plastic strain"; break;
        }

        for(std::vector<GeometryTag>::iterator it = vecLoc.begin(); it!=vecLoc.end(); ++it)
        {
            const GeometryTag curLoc = *it;
            int parentShapeNr = curLoc.parentShapeNr;
            QString suffix = mDB->MapOfBodyNames.value(parentShapeNr)+" ";
            name.append(suffix);
        }

        //! ------------------------------------------------------------------
        //! redefine tags and scope:
        //! initially the scope is "All bodies" for the structural diagnostic
        //! ------------------------------------------------------------------
        int NbBodies = mDB->bodyMap.size();
        ListOfShape scope;
        for(int i=1; i<= NbBodies; i++)
        {
            TopoDS_Solid aSolid = TopoDS::Solid(mDB->bodyMap.value(i));
            scope.Append(aSolid);
        }
        std::vector<GeometryTag> vecLocAllBodies = TopologyTools::generateLocationPairs(mDB, scope);

        data.setValue(vecLocAllBodies);
        Property prop_scopeAllBodies("Geometry",data,Property::PropertyGroup_Scope);
        Property prop_tagsAllBodies("Tags",data,Property::PropertyGroup_Scope);

        //! under scope
        vecProp.push_back(prop_scopingMethod);
        vecProp.push_back(prop_scope);
        vecProp.push_back(prop_tags);

        //! under definition
        Property prop_type("Type ",typeOfEPS,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_type);

        //! ---------------------------------------------------------------------------
        //! under definition
        //! "By" selector: it provides options "Time" and "Set"
        //! Here "By"="Time", so the "Display time" property is created
        //! If "By" is changed into "Set", the DetailViewer replace the "Display time"
        //! with "Set Number"
        //! 0 => "Time"
        //! 1 => "Set number"
        //! ---------------------------------------------------------------------------
        data.setValue(0.0);
        Property prop_By("By",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_By);
        Property prop_displayTime("Display time",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_displayTime);

        //! mode
        data.setValue(0);
        Property prop_modeNumber("Mode number",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_modeNumber);

        //! ---------------------------------------------------
        //! type of scale: "0" => "Autoscale"; "1" => "Custom"
        //! ---------------------------------------------------
        int typeOfScale = 0;
        data.setValue(typeOfScale);
        Property prop_scaleType("Scale type",data,Property::PropertyGroup_ColorBox);
        vecProp.push_back(prop_scaleType);

        //! ------------------------------------------------
        //! Number of intervals: "0" => "Program controlled
        //! ------------------------------------------------
        int NbIntervals = INITIAL_NUMBER_OF_COLORBOX_LEVELS;
        data.setValue(NbIntervals);
        Property prop_NbIntervals("# intervals",data,Property::PropertyGroup_ColorBox);
        vecProp.push_back(prop_NbIntervals);
    }
        break;

        //! -------------
        //! nodal forces
        //! -------------
    case SimulationNodeClass::nodeType_solutionStructuralNodalForces:
    {
        //! read options
        int typeOfForces = addOptions.toInt();

        //! prepare the name
        switch(typeOfForces)
        {
        case 0: name = "Total force"; break;
        case 1: name = "Directional force X"; break;
        case 2: name = "Directional force Y"; break;
        case 3: name = "Directional force Z"; break;
        }

        for(std::vector<GeometryTag>::iterator it = vecLoc.begin(); it!=vecLoc.end(); ++it)
        {
            const GeometryTag curLoc = *it;
            int parentShapeNr = curLoc.parentShapeNr;
            QString suffix = mDB->MapOfBodyNames.value(parentShapeNr)+" ";
            name.append(suffix);
        }

        //! under scope
        vecProp.push_back(prop_scopingMethod);
        vecProp.push_back(prop_scope);
        vecProp.push_back(prop_tags);

        //! under definition
        //! the component (total, X, Y, Z) is provided by the calling function through the addOptions
        Property prop_type("Type ",typeOfForces,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_type);

        //! ---------------------------------------------------------------------------
        //! under definition
        //! "By" selector: it provides options "Time" and "Set"
        //! Here "By"="Time", so the "Display time" property is created
        //! If "By" is changed into "Set", the DetailViewer replace the "Display time"
        //! with "Set Number"
        //! 0 => "Time"
        //! 1 => "Set number"
        //! ---------------------------------------------------------------------------
        data.setValue(0.0);
        Property prop_By("By",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_By);
        Property prop_displayTime("Display time",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_displayTime);

        //! mode
        data.setValue(0);
        Property prop_modeNumber("Mode number",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_modeNumber);

        //! ---------------------------------------------------
        //! type of scale: "0" => "Autoscale"; "1" => "Custom"
        //! ---------------------------------------------------
        int typeOfScale = 0;
        data.setValue(typeOfScale);
        Property prop_scaleType("Scale type",data,Property::PropertyGroup_ColorBox);
        vecProp.push_back(prop_scaleType);

        //! ------------------------------------------------
        //! Number of intervals: "0" => "Program controlled
        //! ------------------------------------------------
        int NbIntervals = INITIAL_NUMBER_OF_COLORBOX_LEVELS;
        data.setValue(NbIntervals);
        Property prop_NbIntervals("# intervals",data,Property::PropertyGroup_ColorBox);
        vecProp.push_back(prop_NbIntervals);
    }
        break;
    case SimulationNodeClass::nodeType_solutionStructuralReactionForce:
    {
        //! read options
        int typeOfForces = addOptions.toInt();

        //! prepare the name
        switch(typeOfForces)
        {
        case 0: name = "Total reaction force"; break;
        case 1: name = "Directional reaction force X"; break;
        case 2: name = "Directional reaction force Y"; break;
        case 3: name = "Directional reaction force Z"; break;
        }

        for(std::vector<GeometryTag>::iterator it = vecLoc.begin(); it!=vecLoc.end(); ++it)
        {
            const GeometryTag curLoc = *it;
            int parentShapeNr = curLoc.parentShapeNr;
            QString suffix = mDB->MapOfBodyNames.value(parentShapeNr)+" ";
            name.append(suffix);
        }

        //! under scope
        vecProp.push_back(prop_scopingMethod);
        vecProp.push_back(prop_scope);
        vecProp.push_back(prop_tags);

        //! under definition
        //! the component (total, X, Y, Z) is provided by the calling function through the addOptions
        Property prop_type("Type ",typeOfForces,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_type);

        //! ---------------------------------------------------------------------------
        //! under definition
        //! "By" selector: it provides options "Time" and "Set"
        //! Here "By"="Time", so the "Display time" property is created
        //! If "By" is changed into "Set", the DetailViewer replace the "Display time"
        //! with "Set Number"
        //! 0 => "Time"
        //! 1 => "Set number"
        //! ---------------------------------------------------------------------------
        data.setValue(0.0);
        Property prop_By("By",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_By);
        Property prop_displayTime("Display time",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_displayTime);

        //! mode
        data.setValue(0);
        Property prop_modeNumber("Mode number",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_modeNumber);

        //! ---------------------------------------------------
        //! type of scale: "0" => "Autoscale"; "1" => "Custom"
        //! ---------------------------------------------------
        int typeOfScale = 0;
        data.setValue(typeOfScale);
        Property prop_scaleType("Scale type",data,Property::PropertyGroup_ColorBox);
        vecProp.push_back(prop_scaleType);

        //! ------------------------------------------------
        //! Number of intervals: "0" => "Program controlled
        //! ------------------------------------------------
        int NbIntervals = INITIAL_NUMBER_OF_COLORBOX_LEVELS;
        data.setValue(NbIntervals);
        Property prop_NbIntervals("# intervals",data,Property::PropertyGroup_ColorBox);
        vecProp.push_back(prop_NbIntervals);
    }
        break;

        //! ----------------------
        //! contact (diagnostics)
        //! ----------------------
    case SimulationNodeClass::nodeType_solutionStructuralContact:
    {
        //! read options
        int aType = addOptions.toInt();

        //! prepare the name
        switch(aType)
        {
        case 0: name = "Contact pressure"; break;
        case 1: name = "Fritcional stress"; break;
        case 2: name = "Penetration"; break;
        case 3: name = "Sliding"; break;
        }

        for(std::vector<GeometryTag>::iterator it = vecLoc.begin(); it!=vecLoc.end(); ++it)
        {
            const GeometryTag curLoc = *it;
            int parentShapeNr = curLoc.parentShapeNr;
            QString suffix = mDB->MapOfBodyNames.value(parentShapeNr)+" ";
            name.append(suffix);
        }

        //! under scope
        vecProp.push_back(prop_scopingMethod);
        vecProp.push_back(prop_scope);
        vecProp.push_back(prop_tags);

        //! under definition
        //! the component (total, X, Y, Z) is provided by the calling function through the addOptions
        Property prop_type("Type ",aType,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_type);

        //! ---------------------------------------------------------------------------
        //! under definition
        //! "By" selector: it provides options "Time" and "Set"
        //! Here "By"="Time", so the "Display time" property is created
        //! If "By" is changed into "Set", the DetailViewer replace the "Display time"
        //! with "Set Number"
        //! 0 => "Time"
        //! 1 => "Set number"
        //! ---------------------------------------------------------------------------
        data.setValue(0.0);
        Property prop_By("By",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_By);
        Property prop_displayTime("Display time",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_displayTime);

        //! mode
        data.setValue(0);
        Property prop_modeNumber("Mode number",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_modeNumber);

        //! ---------------------------------------------------
        //! type of scale: "0" => "Autoscale"; "1" => "Custom"
        //! ---------------------------------------------------
        int typeOfScale = 0;
        data.setValue(typeOfScale);
        Property prop_scaleType("Scale type",data,Property::PropertyGroup_ColorBox);
        vecProp.push_back(prop_scaleType);

        //! ------------------------------------------------
        //! Number of intervals: "0" => "Program controlled
        //! ------------------------------------------------
        int NbIntervals = INITIAL_NUMBER_OF_COLORBOX_LEVELS;
        data.setValue(NbIntervals);
        Property prop_NbIntervals("# intervals",data,Property::PropertyGroup_ColorBox);
        vecProp.push_back(prop_NbIntervals);
    }
        break;

        //! -----------------------------
        //! thermal solution temperature
        //! -----------------------------
    case SimulationNodeClass::nodeType_solutionThermalTemperature:
    {
        name = "Temperature";

        //! under scope
        //vecProp.push_back(prop_scopingMethod);
        //vecProp.push_back(prop_scope);
        //vecProp.push_back(prop_tags);

        //! ------------------------------------------------------------------
        //! redefine tags and scope:
        //! initially the scope is "All bodies" for the structural diagnostic
        //! ------------------------------------------------------------------
        int NbBodies = mDB->bodyMap.size();
        ListOfShape scope;
        for(int i=1; i<= NbBodies; i++)
        {
            TopoDS_Solid aSolid = TopoDS::Solid(mDB->bodyMap.value(i));
            scope.Append(aSolid);
        }
        std::vector<GeometryTag> vecLocAllBodies = TopologyTools::generateLocationPairs(mDB, scope);

        data.setValue(vecLocAllBodies);
        Property prop_scopeAllBodies("Geometry",data,Property::PropertyGroup_Scope);
        Property prop_tagsAllBodies("Tags",data,Property::PropertyGroup_Scope);

        //! under scope
        vecProp.push_back(prop_scopingMethod);
        vecProp.push_back(prop_scopeAllBodies);
        vecProp.push_back(prop_tagsAllBodies);

        //! --------
        //! "Type "
        //! --------
        int dummyType = 0;
        data.setValue(dummyType);
        Property prop_type("Type ",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_type);

        //! ---------------------------------------------------------------------------
        //! under definition
        //! "By" selector: it provides options "Time" and "Set"
        //! Here "By"="Time", so the "Display time" property is created
        //! If "By" is changed into "Set", the DetailViewer replace the "Display time"
        //! with "Set Number"
        //! 0 => "Time"
        //! 1 => "Set number"
        //! ---------------------------------------------------------------------------
        data.setValue(0.0);
        Property prop_By("By",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_By);
        Property prop_displayTime("Display time",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_displayTime);

        //! mode
        data.setValue(0);
        Property prop_modeNumber("Mode number",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_modeNumber);

        //! ---------------------------------------------------
        //! type of scale: "0" => "Autoscale"; "1" => "Custom"
        //! ---------------------------------------------------
        int typeOfScale = 0;
        data.setValue(typeOfScale);
        Property prop_scaleType("Scale type",data,Property::PropertyGroup_ColorBox);
        vecProp.push_back(prop_scaleType);

        //! ------------------------------------------------
        //! Number of intervals: "0" => "Program controlled
        //! ------------------------------------------------
        int NbIntervals = INITIAL_NUMBER_OF_COLORBOX_LEVELS;
        data.setValue(NbIntervals);
        Property prop_NbIntervals("# intervals",data,Property::PropertyGroup_ColorBox);
        vecProp.push_back(prop_NbIntervals);
    }
        break;

    case SimulationNodeClass::nodeType_solutionStructuralGamma:
    {
        name = "Gamma";
        //! ------------------------------------------------------------------
        //! redefine tags and scope:
        //! initially the scope is "All bodies" for the structural diagnostic
        //! ------------------------------------------------------------------
        int NbBodies = mDB->bodyMap.size();
        ListOfShape scope;
        for(int i=1; i<= NbBodies; i++)
        {
            TopoDS_Solid aSolid = TopoDS::Solid(mDB->bodyMap.value(i));
            scope.Append(aSolid);
        }
        std::vector<GeometryTag> vecLocAllBodies = TopologyTools::generateLocationPairs(mDB, scope);

        data.setValue(vecLocAllBodies);
        Property prop_scopeAllBodies("Geometry",data,Property::PropertyGroup_Scope);
        Property prop_tagsAllBodies("Tags",data,Property::PropertyGroup_Scope);

        //! under scope
        vecProp.push_back(prop_scopingMethod);
        vecProp.push_back(prop_scopeAllBodies);
        vecProp.push_back(prop_tagsAllBodies);

        //! --------
        //! "Type "
        //! --------
        int dummyType = 0;
        data.setValue(dummyType);
        Property prop_type("Type ",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_type);

        //! ---------------------------------------------------------------------------
        //! under definition
        //! "By" selector: it provides options "Time" and "Set"
        //! Here "By"="Time", so the "Display time" property is created
        //! If "By" is changed into "Set", the DetailViewer replace the "Display time"
        //! with "Set Number"
        //! 0 => "Time"
        //! 1 => "Set number"
        //! ---------------------------------------------------------------------------
        data.setValue(0.0);
        Property prop_By("By",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_By);
        Property prop_displayTime("Display time",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_displayTime);

        //! mode
        data.setValue(0);
        Property prop_modeNumber("Mode number",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_modeNumber);
/*
        //! mode
        data.setValue(0);
        Property prop_modeNumber("Mode number",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_modeNumber);

        //! mode
        data.setValue(0);
        Property prop_modeNumber("Mode number",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_modeNumber);

        //! mode
        data.setValue(0);
        Property prop_modeNumber("Mode number",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_modeNumber);*/

        //! ---------------------------------------------------
        //! type of scale: "0" => "Autoscale"; "1" => "Custom"
        //! ---------------------------------------------------
        int typeOfScale = 0;
        data.setValue(typeOfScale);
        Property prop_scaleType("Scale type",data,Property::PropertyGroup_ColorBox);
        vecProp.push_back(prop_scaleType);

        //! ------------------------------------------------
        //! Number of intervals: "0" => "Program controlled
        //! ------------------------------------------------
        int NbIntervals = INITIAL_NUMBER_OF_COLORBOX_LEVELS;
        data.setValue(NbIntervals);
        Property prop_NbIntervals("# intervals",data,Property::PropertyGroup_ColorBox);
        vecProp.push_back(prop_NbIntervals);
    }
        break;

        //! ------------------------------
        //! thermal solution thermal flux
        //! ------------------------------
    case SimulationNodeClass::nodeType_solutionThermalFlux:
    {
        name = "Thermal flux";
        //! ------------------------------------------------------------------
        //! redefine tags and scope:
        //! initially the scope is "All bodies" for the structural diagnostic
        //! ------------------------------------------------------------------
        int NbBodies = mDB->bodyMap.size();
        ListOfShape scope;
        for(int i=1; i<= NbBodies; i++)
        {
            TopoDS_Solid aSolid = TopoDS::Solid(mDB->bodyMap.value(i));
            scope.Append(aSolid);
        }
        std::vector<GeometryTag> vecLocAllBodies = TopologyTools::generateLocationPairs(mDB, scope);

        data.setValue(vecLocAllBodies);
        Property prop_scopeAllBodies("Geometry",data,Property::PropertyGroup_Scope);
        Property prop_tagsAllBodies("Tags",data,Property::PropertyGroup_Scope);

        //! under scope
        vecProp.push_back(prop_scopingMethod);
        vecProp.push_back(prop_scopeAllBodies);
        vecProp.push_back(prop_tagsAllBodies);

        //! --------
        //! "Type "
        //! --------
        int dummyType = 0;
        data.setValue(dummyType);
        Property prop_type("Type ",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_type);

        //! ---------------------------------------------------------------------------
        //! under definition
        //! "By" selector: it provides options "Time" and "Set"
        //! Here "By"="Time", so the "Display time" property is created
        //! If "By" is changed into "Set", the DetailViewer replace the "Display time"
        //! with "Set Number"
        //! 0 => "Time"
        //! 1 => "Set number"
        //! ---------------------------------------------------------------------------
        data.setValue(0.0);
        Property prop_By("By",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_By);
        Property prop_displayTime("Display time",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_displayTime);

        //! mode
        data.setValue(0);
        Property prop_modeNumber("Mode number",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_modeNumber);

        //! ---------------------------------------------------
        //! type of scale: "0" => "Autoscale"; "1" => "Custom"
        //! ---------------------------------------------------
        int typeOfScale = 0;
        data.setValue(typeOfScale);
        Property prop_scaleType("Scale type",data,Property::PropertyGroup_ColorBox);
        vecProp.push_back(prop_scaleType);

        //! ------------------------------------------------
        //! Number of intervals: "0" => "Program controlled
        //! ------------------------------------------------
        int NbIntervals = INITIAL_NUMBER_OF_COLORBOX_LEVELS;
        data.setValue(NbIntervals);
        Property prop_NbIntervals("# intervals",data,Property::PropertyGroup_ColorBox);
        vecProp.push_back(prop_NbIntervals);
    }
        break;

        //! -------------
        //! model change
        //! -------------
    case SimulationNodeClass::nodeType_modelChange:
    {
        name = "Model change";

        //! ------------------------------------
        //! "0" model change on body/ies
        //! "1" model change on contact pairs
        //! ------------------------------------
        data.setValue(0);   //! body by default
        Property prop_itemType("Item type",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_itemType);

        vecProp.push_back(prop_scopingMethod);
        vecProp.push_back(prop_scope);
        vecProp.push_back(prop_tags);
    }
        break;

        //! ---------------------
        //! "Electrostatic wall"
        //! ---------------------
    case SimulationNodeClass::nodeType_electrostaticPotential:
    {
        name = "Electrostatic wall";

        vecProp.push_back(prop_scopingMethod);
        vecProp.push_back(prop_scope);
        vecProp.push_back(prop_tags);

        double wallPotential = 0;
        data.setValue(wallPotential);
        Property prop_wallPotential("Potential",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_wallPotential);

        int isEmitter = 0;
        data.setValue(isEmitter);
        Property prop_isEmitter("Emitter",data,Property::PropertyGroup_Definition);
        vecProp.push_back(prop_isEmitter);

    }
        break;
    }

    //! ----------------
    //! create the node
    //! ----------------
    SimulationNodeClass *node = new SimulationNodeClass(name,type,vecProp);

    //! ---------------
    //! add a time tag
    //! ---------------
    node->addTimeTag();

    return node;
}

//! -------------------------------------------
//! function: nodeFromFile
//! details:  create a node from file - unused
//! -------------------------------------------
SimulationNodeClass *nodeFactory::nodeFromFile(const QString &fileName)
{
    cout<<"nodeFactory::nodeFromFile->____function called____"<<endl;
    ifstream in;
    in.open(fileName.toStdString());

    in.close();
    return 0;
}
