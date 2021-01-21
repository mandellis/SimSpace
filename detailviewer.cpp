//! ----------------
//! custom includes
//! ----------------
#include "detailviewer.h"
#include "simulationnodeclass.h"
#include "generaldelegate.h"
#include "property.h"
#include "mydefines.h"
#include "tablewidget.h"
#include "qextendedstandarditem.h"
#include "customtablemodel.h"
#include "simulationmanager.h"
#include "qbackgroundevent.h"
#include "tools.h"
#include "tabulardatacolumns.h"
#include "maintreetools.h"
#include <indexedmapofmeshdatasources.h>

//! ---
//! Qt
//! ---
#include <QStandardItemModel>
#include <QVariant>
#include <QHeaderView>
#include <QApplication>
#include <QDir>
#include <QScrollBar>

//! ----
//! OCC
//! ----
#include <TopoDS_Face.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Vertex.hxx>
#include <gp_Trsf.hxx>
#include <gp_Pnt.hxx>
#include "geomtoolsclass.h"
#include "qextendedstandarditem.h"
#include <TopoDS.hxx>
#include <AIS_Trihedron.hxx>            //! for testing purposes
#include <Prs3d_Drawer.hxx>
#include <Prs3d_DatumAspect.hxx>
#include <Prs3d_LineAspect.hxx>
#include <Quantity_NameOfColor.hxx>
#include <Prs3d_ArrowAspect.hxx>
#include <Prs3d_Arrow.hxx>
#include <Aspect_WidthOfLine.hxx>
#include <Prs3d_ShadingAspect.hxx>
#include <Aspect_TypeOfFacingModel.hxx>
#include <Graphic3d_MaterialAspect.hxx>
#include <Graphic3d_NameOfMaterial.hxx>
#include <DsgPrs_DatumPrs.hxx>
#include <Geom_Surface.hxx>
#include <GeomAdaptor_Surface.hxx>
#include <Geom_Line.hxx>
#include <GeomAdaptor_Curve.hxx>
#include <GeomAbs_CurveType.hxx>
#include <TopoDS_Builder.hxx>
#include <GProp_GProps.hxx>
#include <BRepGProp.hxx>

//! ----
//! C++
//! ----
#include <sstream>

//! ----------------------------------------
//! function: mousePressEvent
//! details:  start editing after one click
//! ----------------------------------------
void DetailViewer::mousePressEvent(QMouseEvent *event)
{
    if (event->button() == Qt::LeftButton)
    {
        QModelIndex index = indexAt(event->pos());
        if (index.column() == 1)
        {
            edit(index);
        }
    }
    QTreeView::mousePressEvent(event);
}

//! ----------------------------
//! function: createConnections
//! details:
//! ----------------------------
void DetailViewer::createConnections()
{
    connect(myGeneralDelegate,SIGNAL(suppressionChanged(Property::SuppressionStatus)),this,SLOT(handleSuppressionPropertyChange(Property::SuppressionStatus)));
    connect(myGeneralDelegate,SIGNAL(visibilityChanged_()),this,SLOT(handleVisibilityChange()));

    connect(myGeneralDelegate,SIGNAL(elementControlChanged()),this,SLOT(handleElementControlChange()));
    connect(myGeneralDelegate,SIGNAL(scopingMethodChanged()),this,SLOT(handleScopingMethodChange()));

    //! -----------------------------------
    //! "updateTags" is a crictical method
    //! -----------------------------------
    connect(myGeneralDelegate,SIGNAL(scopeChanged()),this,SLOT(updateTags()));

    connect(myGeneralDelegate,SIGNAL(numberOfStepChanged()),this,SLOT(handleNumberOfStepChanged()));
    connect(myGeneralDelegate,SIGNAL(StepEndTimeChanged()),this,SLOT(handleStepEndTimeChanged()));
    connect(myGeneralDelegate,SIGNAL(solverTypeChanged()),this,SLOT(handleSolverTypeChanged()));
    connect(myGeneralDelegate,SIGNAL(currentStepNumberChanged()),this,SLOT(handleCurrentStepNumberChanged()));

    //! -------------------
    //! coordinate systems
    //! -------------------
    connect(myGeneralDelegate,SIGNAL(originChanged()),this,SLOT(handleOriginChanged()));
    connect(myGeneralDelegate,SIGNAL(originAndDirectionChanged()),this,SLOT(handleOriginAndDirectionChanged()));
    connect(myGeneralDelegate,SIGNAL(originChangedByValues()),this,SLOT(handleOriginChangedByValue()));
    connect(myGeneralDelegate,SIGNAL(transformationsChanged()),this,SLOT(applyTransformations()));

    connect(myGeneralDelegate,SIGNAL(defineBy_Changed()),this,SLOT(handleDefineBy_Changed()));
    connect(myGeneralDelegate,SIGNAL(defineByChanged()),this,SLOT(handleDefineByChanged()));

    connect(myGeneralDelegate,SIGNAL(loadDefinitionFilmCoefficientChanged(QString)),this,SLOT(handleFilmCoefficientLoadDefinitionChanged(QString)));
    connect(myGeneralDelegate,SIGNAL(loadDefinitionReferenceTemperatureChanged(QString)),this,SLOT(handleReferenceTemperatureLoadDefinitionChanged(QString)));
    connect(myGeneralDelegate,SIGNAL(loadDefinitionMagnitudeChanged(QString)),this,SLOT(handleMagnitudeLoadDefinitionChanged(QString)));
    connect(myGeneralDelegate,SIGNAL(loadDefinitionXChanged(QString)),this,SLOT(handleXLoadDefinitionChanged(QString)));
    connect(myGeneralDelegate,SIGNAL(loadDefinitionYChanged(QString)),this,SLOT(handleYLoadDefinitionChanged(QString)));
    connect(myGeneralDelegate,SIGNAL(loadDefinitionZChanged(QString)),this,SLOT(handleZLoadDefinitionChanged(QString)));

    connect(myGeneralDelegate,SIGNAL(autoTimeSteppingChanged()),this,SLOT(handleAutoTimeSteppingChanged()));
    connect(myGeneralDelegate,SIGNAL(timeDivisionChanged()),this,SLOT(handleTimeDivisionChanged()));
    connect(myGeneralDelegate,SIGNAL(meshNodesVisibilityChanged(bool)),this,SLOT(handleChangeMeshNodesVisibility(bool)));
    connect(myGeneralDelegate,SIGNAL(meshSmoothingChanged()),this,SLOT(handleMeshSmoothingChange()));

    connect(myGeneralDelegate,SIGNAL(requestStartInterpolator()),this,SLOT(handleRequestStartInterpolator()));
    connect(myGeneralDelegate,SIGNAL(requestStartOpenFoamScalarDataTranslator()),this,SLOT(handleRequestStartOpenFoamScalarDataTranslator()));

    connect(myGeneralDelegate,SIGNAL(globalMeshControlChanged()),this,SLOT(handleGlobalMeshControlChange()));

    connect(myGeneralDelegate,SIGNAL(remapFlagChanged()),this,SLOT(handleRemapFlagChanged()));
    connect(myGeneralDelegate,SIGNAL(stepSelectionModeChanged()),this,SLOT(handleStepSelectionModeChange()));

    connect(myGeneralDelegate,SIGNAL(backgroundChanged()),this,SLOT(postBackgroundChangeEvent()));
    connect(myGeneralDelegate,SIGNAL(NbThreadsChanged()),this,SLOT(handleNbThreadsChanged()));

    //connect(myGeneralDelegate,SIGNAL(byChanged()),this,SLOT(handleByChanged()));

    connect(myGeneralDelegate,SIGNAL(solutionComponentChanged()),this,SLOT(handleSolutionComponentChanged()));
    connect(myGeneralDelegate,SIGNAL(typeOfSizingChanged()),this,SLOT(handleTypeOfSizingChanged()));

    connect(myGeneralDelegate,SIGNAL(solutionInformationUpdateIntervalChanged()),this,SLOT(handleSolutionInformationUpdateIntervalChanged()));
    connect(myGeneralDelegate,SIGNAL(BoltStatusDefinedByChanged()),this,SLOT(handleBoltStatusDefinedByChanged()));
    connect(myGeneralDelegate,SIGNAL(BoltLoadChanged()),this,SLOT(handleBoltLoadChanged()));
    connect(myGeneralDelegate,SIGNAL(BoltAdjustmentChanged()),this,SLOT(handleBoltAdjustmentChanged()));

    //! ------
    //! bolts
    //! ------
    connect(myGeneralDelegate,SIGNAL(boltCSChanged()),this,SLOT(handleBoltCSChanged()));

    connect(myGeneralDelegate,SIGNAL(accelerationChanged()),this,SLOT(handleAccelerationChanged()));
    connect(myGeneralDelegate,SIGNAL(momentChanged()),this,SLOT(handleMomentChanged()));

    connect(myGeneralDelegate,SIGNAL(fluxConvegernceChanged()),this,SLOT(handleFluxConvergenceChanged()));
    connect(myGeneralDelegate,SIGNAL(solutionConvegernceChanged()),this,SLOT(handleSolutionConvergenceChanged()));
    connect(myGeneralDelegate,SIGNAL(fieldParametersChanged()),this,SLOT(updateFieldParameters()));
    connect(myGeneralDelegate,SIGNAL(timeIncrementationChanged()),this,SLOT(updateTimeIncrementationParameters()));
    connect(myGeneralDelegate,SIGNAL(cutbackFactorsChanged()),this,SLOT(updateCutBackFactors()));
    connect(myGeneralDelegate,SIGNAL(lineSearchChanged()),this,SLOT(updateLineSearch()));
    connect(myGeneralDelegate,SIGNAL(lineSearchParametersChanged()),this,SLOT(updateLineSearch()));

    connect(myGeneralDelegate,SIGNAL(outputControlsChanged()),this,SLOT(handleOutputControlsChanged()));
    connect(myGeneralDelegate,SIGNAL(storeResultsAtChanged()),this,SLOT(handleStoreResultsAtChanged()));

    //! --------------
    //! remote points
    //! --------------
    connect(myGeneralDelegate,SIGNAL(remotePointChangedByLocation()),this,SLOT(handleRemotePointChangedByLocation()));
    connect(myGeneralDelegate,SIGNAL(remotePointSystemOfReferenceChanged()),this,SLOT(handleRemotePointSystemOfReferenceChanged()));

    //! -----------------
    //! prismatic layers
    //! -----------------
    connect(myGeneralDelegate,SIGNAL(prismaticLayerOptionsChanged()),this,SLOT(handlePrismaticLayerOptions()));
    connect(myGeneralDelegate,SIGNAL(boundaryScopingMethodChanged()),this,SLOT(handleBoundaryScopingMethodChanged()));
    connect(myGeneralDelegate,SIGNAL(boundaryScopeChanged()),this,SLOT(updateBoudaryTags()));

    //! -----------------------------------------------------------------------
    //! at the moment this work has been transferred to the simulation manager
    //! -----------------------------------------------------------------------
    //connect(myGeneralDelegate,SIGNAL(DOFselectorChanged()),this,SLOT(handleDOFselectorChange()));
    //connect(myGeneralDelegate,SIGNAL(couplingChanged()),this,SLOT(handleCouplingChanged()));

    connect(myGeneralDelegate,SIGNAL(meshMethodChanged()),this,SLOT(handleMeshMethodChanged()));
    connect(myGeneralDelegate,SIGNAL(meshSimplificationChanged()),this,SLOT(handleMeshSimplificationChanged()));
    connect(myGeneralDelegate,SIGNAL(meshDefeaturingChanged()),this,SLOT(handleMeshDefeaturingChanged()));

    //! ---------------------------
    //! new mesh method definition
    //! ---------------------------
    connect(myGeneralDelegate,SIGNAL(BRepFlagChanged()),this,SLOT(handleBRepFlagChanged()));
    connect(myGeneralDelegate,SIGNAL(defeaturingFlagChanged()),this,SLOT(handleDefeaturingFlagChanged()));
    connect(myGeneralDelegate,SIGNAL(simplificationFlagChanged()),this,SLOT(handleSimplificationFlagChanged()));
    connect(myGeneralDelegate,SIGNAL(meshSimplificationByChanged()),this,SLOT(handleMeshSimplificationByChanged()));
    connect(myGeneralDelegate,SIGNAL(volumeMesherChanged()),this,SLOT(handleVolumeMesherChanged()));
    connect(myGeneralDelegate,SIGNAL(TessellatorChanged()),this,SLOT(handleTessellatorChanged()));

    //! ------------------------
    //! geometry import options
    //! ------------------------
    connect(myGeneralDelegate,SIGNAL(geometryHealingChanged()),this,SLOT(handleGeometryHealingChanged()));

#ifdef COSTAMP_VERSION
    connect(myGeneralDelegate,SIGNAL(requestStartTimeStepBuilder()),this,SLOT(handleRequestStartTimeStepBuilder()));
#endif

    //! -------------------
    //! TetWild parameters
    //! -------------------
    connect(myGeneralDelegate,SIGNAL(idealLengthChanged()),this,SLOT(handleIdealLengthChanged()));
    connect(myGeneralDelegate,SIGNAL(envelopeSizingChanged()),this,SLOT(handleEnvelopeSizingChanged()));

    //! --------------
    //! interpolation
    //! --------------
    connect(myGeneralDelegate,SIGNAL(interpolationAlgorithmChanged()),this,SLOT(handleInterpolationAlgorithmChanged()));

    //! -------------
    //! model change
    //! -------------
    connect(myGeneralDelegate,SIGNAL(modelChangeScopingMethodChanged()),this,SLOT(handleModelChangeScopingMethodChanged()));
    connect(myGeneralDelegate,SIGNAL(modelChangeActivationStatusChanged()),this,SLOT(handleModelChangeActivationStatusChanged()));

    //! -------------
    //! transparency
    //! -------------
    connect(myGeneralDelegate,SIGNAL(transparencyChanged()),this,SLOT(handleTransparencyChanged()));

    //! ---------------------------------------------
    //! selection method of a mesh element selection
    //! ---------------------------------------------
    connect(myGeneralDelegate,SIGNAL(selectionMethodChanged()),this,SLOT(handleSelectionMethodChanged()));
    connect(myGeneralDelegate,SIGNAL(elementListChanged()),this,SLOT(handleMeshElementListChanged()));

    connect(myGeneralDelegate,SIGNAL(meshMetricChanged()),this,SLOT(handleMeshMetricChanged()));

    //! -------------
    //! fatigue tool
    //! -------------
    connect(myGeneralDelegate,SIGNAL(fatigueAlgoChanged()),this,SLOT(handleFatigueAlgoChanged()));

    //! --------------
    //! analysis type
    //! --------------
    connect(myGeneralDelegate,SIGNAL(analysisTypeChangedChanged()),this,SLOT(handleAnalysisTypeChanged()));

    //! -----------------
    //! Time integration
    //! -----------------
    connect(myGeneralDelegate,SIGNAL(timeIntegrationChanged()),this,SLOT(handleTimeIntegrationChanged()));

    //! --------
    //! Emitter
    //! --------
    connect(myGeneralDelegate,SIGNAL(EmitterStatusChanged()),this,SLOT(handleEmitterStatusChanged()));
}

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
DetailViewer::DetailViewer(QWidget *parent): QTreeView(parent)
{
    cout<<"DetailViewer::DetailViewer()->____default constructor called____"<<endl;

    //! drawing style
    this->setAlternatingRowColors(true);

    //! delegate
    myGeneralDelegate = new GeneralDelegate(this);
    this->setItemDelegate(myGeneralDelegate);

    //! items with child(ren) have an expander
    this->setRootIsDecorated(true);

    //! hide the headers
    //this->header()->hide();

    //! the interactive context [...] unitialized

    //! create the connections
    this->createConnections();

    //! allow multiple selection
    this->setSelectionMode(QAbstractItemView::MultiSelection);

    //! horizontal scrollbar
    this->horizontalScrollBar()->setEnabled(true);
    this->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOn);

    //! set object name
    this->setObjectName("detailViewer");
}

//! ------------------------
//! function: constructor I
//! details:
//! ------------------------
DetailViewer::DetailViewer(const occHandle(AIS_InteractiveContext) &aCTX, QWidget *parent):QTreeView(parent),myCTX(aCTX)
{
    cout<<"DetailViewer::DetailViewer()->____CONSTRUCTOR I CALLED____"<<endl;

    //! delegate    
    myGeneralDelegate = new GeneralDelegate(myCTX, this);
    this->setItemDelegate(myGeneralDelegate);

    //! items with child(ren) have an expander
    this->setRootIsDecorated(true);

    //! hide the headers
    //this->header()->hide();

    //! create connections
    this->createConnections();

    //! allow multiple selection
    this->setSelectionMode(QAbstractItemView::MultiSelection);

    //! set object name
    this->setObjectName("detailViewer");
}

//! --------------------------------
//! function: setContext
//! details:  also for the delegate
//! --------------------------------
void DetailViewer::setContext(const occHandle(AIS_InteractiveContext) &aCTX)
{
    cout<<"DetailViewer::setContext()->____function called____"<<endl;
    if(aCTX.IsNull()==true)
    {
        cout<<"DetailViewer::setContext()->____NULL CONTEXT____"<<endl;
        return;
    }
    myCTX = aCTX;
    myGeneralDelegate->setContext(myCTX);
    cout<<"DetailViewer::setContext()->____context OK____"<<endl;
}

//! --------------------------------
//! function: setMeshContext
//! details:  also for the delegate
//! --------------------------------
void DetailViewer::setMeshContext(const occHandle(AIS_InteractiveContext) &aMeshCTX)
{
    cout<<"DetailViewer::setMeshContext()->____function called____"<<endl;
    if(aMeshCTX.IsNull()==true)
    {
        cout<<"DetailViewer::setMeshContext()->____NULL MESH CONTEXT____"<<endl;
        return;
    }
    myMeshCTX = aMeshCTX;
    myGeneralDelegate->setMeshContext(aMeshCTX);
    cout<<"DetailViewer::setMeshContext()->____mesh context OK____"<<endl;
}

//! ----------------------
//! function: setTheModel
//! details:
//! ----------------------
void DetailViewer::setTheModel(const QModelIndex &anIndex)
{
    //cout<<"DetailViewer::setTheModel()->____function called____"<<endl;

    //! ---------------------
    //! set the window title
    //! ---------------------
    QString windowTitle = QString("Details of ").append(anIndex.data(Qt::DisplayRole).toString());
    this->parentWidget()->setWindowTitle(windowTitle);

    QVariant data = anIndex.data(Qt::UserRole);
    SimulationNodeClass *aNode = data.value<SimulationNodeClass*>();
    this->setModel(aNode->getModel());

    //! -----------------------------------------------------
    //! warning: please, do not use Qt TreeView::expandAll()
    //! since it chashes (I don't know the reason). Use the
    //! DetailViewer::expandChildren() instead
    //! -----------------------------------------------------
    this->expandChildren();

    //! ----------------------------------------
    //! the source model index and current node
    //! ----------------------------------------
    myCurModelIndex = anIndex;
    myCurNode = data.value<SimulationNodeClass*>();

    SimulationNodeClass::nodeType nodeType = myCurNode->getType();
    switch (nodeType)
    {
    case SimulationNodeClass::nodeType_remotePoint:
    {
        /*
        //! ------------------------------
        //! hide the absolute coordinates
        //! ------------------------------
        this->setPropertyVisible("X abs coordinate",true);
        this->setPropertyVisible("Y abs coordinate",true);
        this->setPropertyVisible("Z abs coordinate",true);
        */
    }
        break;

    case SimulationNodeClass::nodeType_namedSelectionGeometry:
    {
        QExtendedStandardItem *item = myCurNode->getPropertyItem("Scoping method");
        if(item!=Q_NULLPTR) this->setRowHidden(item->row(),item->parent()->index(),true);
    }
        break;

    case SimulationNodeClass::nodeType_meshEdgeSize:
    {
        int sizingType = myCurNode->getPropertyItem("Sizing type")->data(Qt::UserRole).value<Property>().getData().toInt();
        if(sizingType==0)
        {
            QExtendedStandardItem *itemNbDivisions = myCurNode->getPropertyItem("Number of divisions");
            QModelIndex index = itemNbDivisions->index();
            this->setRowHidden(index.row(),index.parent(),true);
        }
        else
        {
            QExtendedStandardItem *itemElementSize = myCurNode->getPropertyItem("Element size");
            QModelIndex index = itemElementSize->index();
            this->setRowHidden(index.row(),index.parent(),true);
        }
    }
        break;

    case SimulationNodeClass::nodeType_structuralAnalysisBoltPretension:
    {
        QExtendedStandardItem *itemLoad = myCurNode->getPropertyItem("Load");
        QModelIndex indexLoad = itemLoad->index();
        QExtendedStandardItem *itemAdjustment = myCurNode->getPropertyItem("Adjustment");
        QModelIndex indexAdjustment = itemAdjustment->index();

        Property::boltStatusDefinedBy boltStatusDefinedBy = myCurNode->getPropertyValue<Property::boltStatusDefinedBy>("Bolt status");

        switch(boltStatusDefinedBy)
        {
        case Property::boltStatusDefinedBy_load:
        {
            this->setRowHidden(indexLoad.row(),indexLoad.parent(),false);
            this->setRowHidden(indexAdjustment.row(),indexAdjustment.parent(),true);

            //! ------------------------------------------
            //! the "Adjustment" control must be disabled
            //! ------------------------------------------
            QVariant data;
            data.setValue(QString("N/A"));
            SimulationNodeClass *nodeAnalysisSettings = myCurModelIndex.parent().child(0,0).data(Qt::UserRole).value<SimulationNodeClass*>();
            int currentStepNumber = nodeAnalysisSettings->getPropertyItem("Current step number")->data(Qt::UserRole).value<Property>().getData().toInt();

            SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
            int SC = mainTreeTools::calculateStartColumn(sm->myTreeView);
            int row = currentStepNumber;
            int col = SC+2;
            CustomTableModel *tabularDataModel = nodeAnalysisSettings->getTabularDataModel();
            QModelIndex index = tabularDataModel->makeIndex(row,col);
            tabularDataModel->setData(index,data,Qt::EditRole);
        }
            break;

        case Property::boltStatusDefinedBy_adjustment:
        {
            this->setRowHidden(indexLoad.row(),indexLoad.parent(),true);
            this->setRowHidden(indexAdjustment.row(),indexAdjustment.parent(),false);

            //! ------------------------------------
            //! the "Load" control must be disabled
            //! ------------------------------------
            QVariant data;
            data.setValue(QString("N/A"));
            SimulationNodeClass *nodeAnalysisSettings = myCurModelIndex.parent().child(0,0).data(Qt::UserRole).value<SimulationNodeClass*>();
            int currentStepNumber = nodeAnalysisSettings->getPropertyValue<int>("Current step number");

            SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
            int SC = mainTreeTools::calculateStartColumn(sm->myTreeView);

            int row = currentStepNumber;
            int col = SC+1;
            CustomTableModel *tabularDataModel = nodeAnalysisSettings->getTabularDataModel();
            QModelIndex index = tabularDataModel->makeIndex(row,col);
            tabularDataModel->setData(index,data,Qt::EditRole);
        }
            break;

        case Property::boltStatusDefinedBy_open:
        case Property::boltStatusDefinedBy_lock:
        {
            //! -----------------------------------------------------------
            //! both the "Load" and "Adjustment" controls must be disabled
            //! -----------------------------------------------------------
            this->setRowHidden(indexLoad.row(),indexLoad.parent(),false);
            this->setRowHidden(indexAdjustment.row(),indexAdjustment.parent(),false);
            QVariant data;
            data.setValue(QString("N/A"));
            SimulationNodeClass *nodeAnalysisSettings = myCurModelIndex.parent().child(0,0).data(Qt::UserRole).value<SimulationNodeClass*>();
            int currentStepNumber = nodeAnalysisSettings->getPropertyValue<int>("Current step number");

            SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
            int SC = mainTreeTools::calculateStartColumn(sm->myTreeView);
            int row = currentStepNumber;

            CustomTableModel *tabularDataModel = nodeAnalysisSettings->getTabularDataModel();
            QModelIndex indexLoad = tabularDataModel->makeIndex(row,SC+1);
            tabularDataModel->setData(indexLoad,data,Qt::EditRole);

            QModelIndex indexAdjustment = tabularDataModel->makeIndex(row,SC+2);
            tabularDataModel->setData(indexAdjustment,data,Qt::EditRole);
        }
            break;
        }
    }
        break;

    case SimulationNodeClass::nodeType_structuralAnalysisSettings:
    case SimulationNodeClass::nodeType_thermalAnalysisSettings:
    {
        //! --------------------------------------------------------------------------------------------------
        //! handle the visibility of some Analysis settings according to "Program controlled", "Custom", ecc
        //! --------------------------------------------------------------------------------------------------
        bool isHidden;
        QStandardItem *theParentItem;

        int fluxConvergence = myCurNode->getPropertyValue<int>("Flux convergence");
        QExtendedStandardItem *itemFluxConvergence = this->getNode()->getPropertyItem("Flux convergence");
        theParentItem = itemFluxConvergence->parent();
        if(fluxConvergence == 2 || fluxConvergence ==0) isHidden = true; else isHidden = false;
        for(int row=1; row<=6; row++) this->setRowHidden(row,theParentItem->index(),isHidden);

        int solutionConvergence = myCurNode->getPropertyValue<int>("Solution convergence");
        QExtendedStandardItem *itemSolutionCOnvergence = this->getNode()->getPropertyItem("Solution convergence");
        theParentItem = itemSolutionCOnvergence->parent();
        if(solutionConvergence == 2) isHidden = true; else isHidden = false;
        for(int row = 8; row<=9; row++) this->setRowHidden(row,theParentItem->index(),isHidden);

        QExtendedStandardItem *itemTimeIncrementation = this->getNode()->getPropertyItem("Time incrementation");
        theParentItem = itemTimeIncrementation->parent();
        int timeIncrementation = myCurNode->getPropertyValue<int>("Time incrementation");
        if(timeIncrementation == 0) isHidden = true; else isHidden = false;
        for(int row = 1; row<=theParentItem->rowCount(); row++) this->setRowHidden(row,theParentItem->index(),isHidden);

        QExtendedStandardItem *itemCutbackFactors = this->getNode()->getPropertyItem("Cutback factors");
        theParentItem = itemCutbackFactors->parent();
        int cutBackFactors = myCurNode->getPropertyValue<int>("Cutback factors");
        if(cutBackFactors == 0) isHidden = true; else isHidden = false;
        for(int row = 1; row<=theParentItem->rowCount(); row++) this->setRowHidden(row,theParentItem->index(),isHidden);

        QExtendedStandardItem *itemLineSearch = this->getNode()->getPropertyItem("Line search");
        theParentItem = itemLineSearch->parent();
        int lineSearch = myCurNode->getPropertyValue<int>("Line search");
        if(lineSearch == 0) isHidden = true; else isHidden = false;
        for(int row = 1; row<=theParentItem->rowCount(); row++) this->setRowHidden(row,theParentItem->index(),isHidden);
    }
        break;

    case SimulationNodeClass::nodeType_connectionPair:
    case SimulationNodeClass::nodeType_connectionGroup:
    {
        //! -------------------------
        //! hide the "Time tag" item
        //! -------------------------
        //this->setPropertyVisible("Time tag");

        //! ------------------------------------
        //! hide "Tags master" and "Tags slave"
        //! ------------------------------------
        QStandardItem *itemTagsMaster = myCurNode->getPropertyItem("Tags master");
        if(itemTagsMaster!=Q_NULLPTR)
        {
            int row =itemTagsMaster->row();
            QModelIndex separatorIndex = itemTagsMaster->parent()->index();
            this->setRowHidden(row,separatorIndex,true);
        }
        QStandardItem *itemTagsSlave = myCurNode->getPropertyItem("Tags slave");
        if(itemTagsSlave!=Q_NULLPTR)
        {
            int row =itemTagsSlave->row();
            QModelIndex separatorIndex = itemTagsSlave->parent()->index();
            this->setRowHidden(row,separatorIndex,true);
        }
    }
        break;

    case SimulationNodeClass::nodeType_meshPrismaticLayer:
    {
        this->setPropertyVisible("Tags");
        this->setPropertyVisible("Boundary tags");
    }
        break;
    }

    SimulationNodeClass::nodeType nodeFamily = myCurNode->getFamily();
    switch(nodeFamily)
    {
    case SimulationNodeClass::nodeType_meshControl:
    {
        //! ----------------------------
        //! hide "Tags" for a mesh item
        //! ----------------------------
        //QExtendedStandardItem *item = myCurNode->getPropertyItem("Tags");
        //if(item!=Q_NULLPTR) this->setRowHidden(item->row(),item->parent()->index(),true);
    }
        break;

    case SimulationNodeClass::nodeType_namedSelection:
    {
        //! ---------------------------------------
        //! hide "Tags" for a named selection item
        //! ---------------------------------------
        QExtendedStandardItem *item = myCurNode->getPropertyItem("Tags");
        if(item!=Q_NULLPTR) this->setRowHidden(item->row(),item->parent()->index(),true);
    }
        break;

        /*
    case SimulationNodeClass::nodeType_coordinateSystems:
    {
        //! -----------------------------------------
        //! hide "Tags" for a coordinate system item
        //! -----------------------------------------
        QExtendedStandardItem *item = myCurNode->getPropertyItem("Tags");
        if(item!=Q_NULLPTR) this->setRowHidden(item->row(),item->parent()->index(),true);
    }
        break;
        */
    }
}


//! ------------------------------------------------------------------------------------
//! function: setTheModel
//! details:  experimental overload - use for handling multiple contact pairs selection
//! ------------------------------------------------------------------------------------
void DetailViewer::setTheModel(SimulationNodeClass* aNode)
{
    cout<<"DetailViewer::setTheModel()->____function (overload) called____"<<endl;
    if(aNode==Q_NULLPTR) return;

    //! ---------------------
    //! set the window title
    //! ---------------------
    QString windowTitle = QString("Details of ").append(aNode->getName());
    this->parentWidget()->setWindowTitle(windowTitle);
    this->setModel(aNode->getModel());

    //! -----------------------------------------------------
    //! warning: please, do not use Qt TreeView::expandAll()
    //! since it chashes (I don't know the reason). Use the
    //! DetailViewer::expandChildren() instead
    //! -----------------------------------------------------
    this->expandChildren();

    //! -----------------------
    //! the source model index
    //! -----------------------
    myCurModelIndex = QModelIndex();

    //! -----------------------
    //! store the current node
    //! -----------------------
    myCurNode = aNode;

    //! -------------------------
    //! hide the "Time tag" item
    //! -------------------------
    //this->setPropertyVisible("Time tag");

    //! ------------------------------------
    //! hide "Tags master" and "Tags slave"
    //! ------------------------------------
    QStandardItem *itemTagsMaster = myCurNode->getPropertyItem("Tags master");
    if(itemTagsMaster!=Q_NULLPTR)
    {
        int row =itemTagsMaster->row();
        QModelIndex separatorIndex = itemTagsMaster->parent()->index();
        this->setRowHidden(row,separatorIndex,true);
    }
    QStandardItem *itemTagsSlave = myCurNode->getPropertyItem("Tags slave");
    if(itemTagsSlave!=Q_NULLPTR)
    {
        int row =itemTagsSlave->row();
        QModelIndex separatorIndex = itemTagsSlave->parent()->index();
        this->setRowHidden(row,separatorIndex,true);
    }

    cout<<"DetailViewer::setTheModel()->____exiting function (overload)____"<<endl;
}

//! --------------------------
//! function: clearTree
//! details:  clear tree view
//! --------------------------
void DetailViewer::clearTree()
{
    if(this->model()!=Q_NULLPTR)
    static_cast<QStandardItemModel*>(this->model())->clear();
    this->parentWidget()->setWindowTitle("Detail viewer");
}

//! ------------------------------------------
//! function: handleSuppressionPropertyChange
//! details:  to do ... beta
//! ------------------------------------------
#include <occPreGLwidget.h>
void DetailViewer::handleSuppressionPropertyChange(Property::SuppressionStatus newSuppressionStatus)
{
    QString test = newSuppressionStatus==Property::SuppressionStatus_Active? "ACTIVATION":"SUPPRESSION";
    cout<<"DetailViewer::handleSuppressionPropertyChange()->____request "<<test.toStdString()<<"____"<<endl;

    occPreGLWidget *occViewer = qobject_cast<occPreGLWidget*>(this->sender());
    if(occViewer==Q_NULLPTR)
    {
        cerr<<"____REQUEST COMING FROM THE MAIN TREE____"<<endl;
    }
    else
    {
        cerr<<"____REQUEST COMING FROM VIEWER____"<<endl;
    }
    emit requestChangeNodeSuppressionStatus(newSuppressionStatus);
}

//! ---------------------------
//! function: handleVisibility
//! details:  to do ... beta
//! ---------------------------
void DetailViewer::handleVisibilityChange()
{
    QExtendedStandardItem *item_visibility = myCurNode->getPropertyItem("Visible");
    if(item_visibility!=Q_NULLPTR)
    {
        bool isVisibile = item_visibility->data(Qt::UserRole).value<Property>().getData().toBool();
        cout<<"DetailViewer::handleVisibilityChange()->____the new visibility is "<<isVisibile<<"____"<<endl;
        if(isVisibile==true) emit requestHandleVisibilityChange(true);
        else emit requestHandleVisibilityChange(false);
    }
}

//! -------------------------------------------------------------
//! function: handleElementControlChange
//! details:  "Element control" refers to the integration method
//! -------------------------------------------------------------
void DetailViewer::handleElementControlChange()
{
    emit requestChangeElementControl();
}

//! -----------------------------------
//! function: handleScopingMehodChange
//! details:
//! -----------------------------------
void DetailViewer::handleScopingMethodChange()
{
    cout<<"DetailViewer::handleScopingMethodChange()->____function called____"<<endl;

    //myCurNode->getModel()->blockSignals(true);
    this->connectToSimulationManager(false);

    QExtendedStandardItem* item = myCurNode->getPropertyItem("Scoping method");
    if(item!=Q_NULLPTR)
    {
        Property::ScopingMethod theScopingMethod = item->data(Qt::UserRole).value<Property>().getData().value<Property::ScopingMethod>();
        switch(theScopingMethod)
        {
        case Property::ScopingMethod_GeometrySelection:
        {
            myCurNode->removeProperty("Remote points");
            if(myCurNode->getType()==SimulationNodeClass::nodeType_connectionPair)
            {
                myCurNode->removeProperty("Master");
                myCurNode->removeProperty("Tags master");
                myCurNode->removeProperty("Slave");
                myCurNode->removeProperty("Tags slave");

                myCurNode->removeProperty("Master mesh data sources");
                myCurNode->removeProperty("Slave mesh data sources");

                QVariant data;
                std::vector<GeometryTag> vecLoc;
                data.setValue(vecLoc);

                Property prop_master("Master",data,Property::PropertyGroup_Scope);
                Property prop_slave("Slave",data,Property::PropertyGroup_Scope);

                Property prop_tagsMaster("Tags master",data,Property::PropertyGroup_Scope);
                Property prop_tagsSlave("Tags slave",data,Property::PropertyGroup_Scope);

                myCurNode->addProperty(prop_master);
                myCurNode->addProperty(prop_tagsMaster);
                myCurNode->addProperty(prop_slave);
                myCurNode->addProperty(prop_tagsSlave);

                //! ----------------------------------------------------------------------------------------
                //! The private slot "updateTags()" is called when the delegate emits "scopeChanged()"
                //! Here the "Scoping method" is changed, so "updateTags()" is not called: call it manually
                //! ----------------------------------------------------------------------------------------
                this->updateTags();
            }
            else
            {
                myCurNode->removeProperty("Named selection");
                std::vector<GeometryTag> curTags = myCurNode->getPropertyValue<std::vector<GeometryTag>>("Tags");
                if(curTags.size()!=0)
                {
                    QVariant data;
                    data.setValue(curTags);
                    myCurNode->addProperty(Property("Geometry",data,Property::PropertyGroup_Scope),1);
                }
                else
                {
                    myCurNode->removeProperty("Tags");

                    std::vector<GeometryTag> vecLoc;
                    QVariant data;
                    data.setValue(vecLoc);
                    Property property_scope("Geometry",data,Property::PropertyGroup_Scope);
                    Property property_tags("Tags",data,Property::PropertyGroup_Scope);
                    myCurNode->addProperty(property_scope,1);
                    myCurNode->addProperty(property_tags,2);

                    //! ----------------------------------------------------------------------------------------
                    //! The private slot "updateTags()" is called when the delegate emits "scopeChanged()"
                    //! Here the "Scoping method" is changed, so "updateTags()" is not called: call it manually
                    //! ----------------------------------------------------------------------------------------
                    this->updateTags();
                }
            }
        }
            break;

        case Property::ScopingMethod_NamedSelection:
        {
            myCurNode->removeProperty("Remote points");
            if(myCurNode->getType()==SimulationNodeClass::nodeType_connectionPair)
            {
                myCurNode->removeProperty("Master");
                myCurNode->removeProperty("Slave");
                myCurNode->removeProperty("Tags master");
                myCurNode->removeProperty("Tags slave");

                myCurNode->removeProperty("Master mesh data sources");
                myCurNode->removeProperty("Slave mesh data sources");

                SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
                QExtendedStandardItem *itemNSRoot = sm->getTreeItem(SimulationNodeClass::nodeType_namedSelection);
                QExtendedStandardItem *itemNSempty = static_cast<QExtendedStandardItem*>(itemNSRoot->child(0,0));
                void *p = (void*)itemNSempty;

                QVariant data;
                data.setValue(p);
                Property prop_master("Master",data,Property::PropertyGroup_Scope);
                Property prop_slave("Slave",data,Property::PropertyGroup_Scope);

                std::vector<GeometryTag> vecLocs;
                data.setValue(vecLocs);

                Property prop_tagsMaster("Tags master",data,Property::PropertyGroup_Scope);
                Property prop_tagsSlave("Tags slave",data,Property::PropertyGroup_Scope);

                myCurNode->addProperty(prop_master);
                myCurNode->addProperty(prop_tagsMaster);

                myCurNode->addProperty(prop_slave);
                myCurNode->addProperty(prop_tagsSlave);

                //! ----------------------------------------------------------------------------------------
                //! The private slot "updateTags()" is called when the delegate emits "scopeChanged()"
                //! Here the "Scoping method" is changed, so "updateTags()" is not called: call it manually
                //! ----------------------------------------------------------------------------------------
                this->updateTags();
            }
            else
            {
                myCurNode->removeProperty("Geometry");
                myCurNode->removeProperty("Tags");

                SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
                QExtendedStandardItem *itemNSRoot = sm->getTreeItem(SimulationNodeClass::nodeType_namedSelection);
                QExtendedStandardItem *itemNSempty = static_cast<QExtendedStandardItem*>(itemNSRoot->child(0,0));
                void *p = (void*)itemNSempty;
                QVariant data;
                data.setValue(p);

                Property prop_namedSelection("Named selection",data,Property::PropertyGroup_Scope);
                myCurNode->addProperty(prop_namedSelection,1);
                std::vector<GeometryTag> vecLoc;
                data.setValue(vecLoc);
                Property prop_tags("Tags",data,Property::PropertyGroup_Scope);
                myCurNode->addProperty(prop_tags,2);

                //! ----------------------------------------------------------------------------------------
                //! The private slot "updateTags()" is called when the delegate emits "scopeChanged()"
                //! Here the "Scoping method" is changed, so "updateTags()" is not called: call it manually
                //! ----------------------------------------------------------------------------------------
                this->updateTags();
            }
        }
            break;

        case Property::ScopingMethod_RemotePoint:
        {
            cout<<"DetailViewer::handleScopingMethodChange()->____updating \"Tags\" from \"Remote point\"____"<<endl;
            //! ---------------------------------------
            //! remove "Geometry" or "Named selection"
            //! ---------------------------------------
            myCurNode->removeProperty("Geometry");
            myCurNode->removeProperty("Named selection");
            //myCurNode->removeProperty("Tags");

            //! -------------------
            //! add "Remote point"
            //! -------------------
            QVariant data;

            //! scan the Remote point branch
            SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
            QStandardItem *itemRemotePointRoot = sm->getTreeItem(SimulationNodeClass::nodeType_remotePointRoot);

            //! this "if" statement is unessential, since the delegate does not return
            //! any editor if the remote point root has not been created
            if(itemRemotePointRoot!=Q_NULLPTR)
            {
                QStandardItem *itemEmptyRemotePoint = itemRemotePointRoot->child(0,0);
                void *p = (void*)itemEmptyRemotePoint;
                data.setValue(p);
                Property prop_remotePoint("Remote points",data,Property::PropertyGroup_Scope);
                myCurNode->addProperty(prop_remotePoint,1);
            }
            //! ----------------------------------------------------------------------------------------
            //! The private slot "updateTags()" is called when the delegate emits "scopeChanged()"
            //! Here the "Scoping method" is changed, so "updateTags()" is not called: call it manually
            //! ----------------------------------------------------------------------------------------
            this->updateTags();
        }
            break;

        case Property::ScopingMethod_Automatic:
        {
            //! --------------------------------------
            //! this holds only for a connection pair
            //! --------------------------------------
            QVariant data;
            IndexedMapOfMeshDataSources emptyMapOfMeshDS;   //cese
            data.setValue(emptyMapOfMeshDS);
            if(myCurNode->getPropertyItem("Master mesh data sources")==Q_NULLPTR)
            {
                Property prop_masterMeshDS("Master mesh data sources",data,Property::PropertyGroup_MeshDataSources);
                myCurNode->addProperty(prop_masterMeshDS);
            }
            if(myCurNode->getPropertyItem("Slave mesh data sources")==Q_NULLPTR)
            {
                Property prop_slaveMeshDS("Slave mesh data sources",data,Property::PropertyGroup_MeshDataSources);
                myCurNode->addProperty(prop_slaveMeshDS);
            }
        }
            break;
        }
    }
    //myCurNode->getModel()->blockSignals(false);
    this->connectToSimulationManager(true);
}

//! -----------------------------------------------------------
//! function: updateTags
//! details:  for "Tags"/"Slave tags"/"Master tags" copies the
//!           content (std::vector<GeometryTag>) into the property
//! -----------------------------------------------------------
#include "markers.h"
#include "markerbuilder.h"
#include <datasourcebuilder.h>
#include <indexedmapofmeshdatasources.h>
#include "maintreetools.h"
void DetailViewer::updateTags()
{
    cout<<"DetailViewer::updateTags()->____function called____"<<endl;

    //myCurNode->getModel()->blockSignals(true);

    //! -----------------------------------------------------------------------
    //! it is necessary to treat the "Named selection" and "Coordinate system"
    //! separately, since they don't have the property "Scoping method"
    //! -----------------------------------------------------------------------
    SimulationNodeClass::nodeType nodeType = myCurNode->getType();
    SimulationNodeClass::nodeType theFamily = myCurNode->getFamily();

    //! ---------------------------------------------------
    //! NOT a "Named selection", NOR a "Coordinate system"
    //! ---------------------------------------------------
    if(theFamily!=SimulationNodeClass::nodeType_namedSelection && nodeType!=SimulationNodeClass::nodeType_coordinateSystem)
    {
        Property::ScopingMethod scopingMethod = myCurNode->getPropertyValue<Property::ScopingMethod>("Scoping method");
        switch(scopingMethod)
        {
        case Property::ScopingMethod_GeometrySelection:
        {
            //! -------------------------------------------
            //! scope defined by direct geometry selection
            //! -------------------------------------------
            QExtendedStandardItem *itemGeometry = myCurNode->getPropertyItem("Geometry");
            QExtendedStandardItem *itemLocation = myCurNode->getPropertyItem("Location");
            QExtendedStandardItem *itemMaster = myCurNode->getPropertyItem("Master");
            QExtendedStandardItem *itemSlave = myCurNode->getPropertyItem("Slave");

            //! -----------------------------------------------------
            //! extract the scope from "Geometry" or from "Location"
            //! actually the case "Location" should never occur
            //! -----------------------------------------------------
            QVariant data;
            std::vector<GeometryTag> vecLoc;
            if(itemGeometry!=Q_NULLPTR || itemLocation!=Q_NULLPTR)
            {
                cout<<"DetailViewer::updateTags()->____Reading scope from \"Geometry\" or \"Location\"____"<<endl;
                std::vector<GeometryTag> vecLoc = itemGeometry->data(Qt::UserRole).value<Property>().getData().value<std::vector<GeometryTag>>();
                data.setValue(vecLoc);

                Property prop_tags("Tags",data,Property::PropertyGroup_Scope);
                myCurNode->replaceProperty("Tags",prop_tags);
            }
            else if(itemMaster!=Q_NULLPTR && itemSlave!=Q_NULLPTR) //! "&&" is unessential, since master is always present with slave
            {
                cout<<"DetailViewer::updateTags()->____Reading scope from \"Master\"____"<<endl;
                vecLoc = itemMaster->data(Qt::UserRole).value<Property>().getData().value<std::vector<GeometryTag>>();
                data.setValue(vecLoc);
                Property prop_masterTags("Tags master",data,Property::PropertyGroup_Scope);
                myCurNode->replaceProperty("Tags master",prop_masterTags);

                cout<<"DetailViewer::updateTags()->____Reading scope from \"Slave\"____"<<endl;
                vecLoc = itemSlave->data(Qt::UserRole).value<Property>().getData().value<std::vector<GeometryTag>>();
                data.setValue(vecLoc);

                Property prop_slaveTags("Tags slave",data,Property::PropertyGroup_Scope);
                myCurNode->replaceProperty("Tags slave",prop_slaveTags);
            }
        }
        break;

        case Property::ScopingMethod_NamedSelection:
        {
            //! -----------------------------------
            //! scope defined by "Named selection"
            //! -----------------------------------
            QExtendedStandardItem *itemGeometry = myCurNode->getPropertyItem("Named selection");
            QExtendedStandardItem *itemMaster = myCurNode->getPropertyItem("Master");
            QExtendedStandardItem *itemSlave = myCurNode->getPropertyItem("Slave");

            QExtendedStandardItem *itemBoundaryNamedSelection = myCurNode->getPropertyItem("Boundary named selection");
            if(itemBoundaryNamedSelection!=Q_NULLPTR)
            {
                void *p = itemBoundaryNamedSelection->data(Qt::UserRole).value<Property>().getData().value<void*>();

                QExtendedStandardItem *itemNS = static_cast<QExtendedStandardItem*>(p);

                SimulationNodeClass *nodeNS = itemNS->data(Qt::UserRole).value<SimulationNodeClass*>();
                std::vector<GeometryTag> vecLoc = nodeNS->getPropertyValue<std::vector<GeometryTag>>("Tags");

                QVariant data;
                data.setValue(vecLoc);
                Property prop_tags("Boundary tags",data,Property::PropertyGroup_Scope);
                myCurNode->replaceProperty("Boundary tags",prop_tags);
            }

            if(itemGeometry!=Q_NULLPTR)
            {
                void *p = itemGeometry->data(Qt::UserRole).value<Property>().getData().value<void*>();

                QExtendedStandardItem *itemNS = static_cast<QExtendedStandardItem*>(p);

                SimulationNodeClass *nodeNS = itemNS->data(Qt::UserRole).value<SimulationNodeClass*>();
                std::vector<GeometryTag> vecLoc = nodeNS->getPropertyValue<std::vector<GeometryTag>>("Tags");

                QVariant data;
                data.setValue(vecLoc);
                Property prop_tags("Tags",data,Property::PropertyGroup_Scope);
                myCurNode->replaceProperty("Tags",prop_tags);
            }
            if(itemMaster!=Q_NULLPTR && itemSlave!=Q_NULLPTR)
            {
                QExtendedStandardItem *itemNS;
                SimulationNodeClass *nodeNS;
                std::vector<GeometryTag> vecLoc;

                //! -----------------------------------------------------------------
                //! [1] rebuild the pairs (parent shape, child shape) for the Master
                //! -----------------------------------------------------------------
                void *p = itemMaster->data(Qt::UserRole).value<Property>().getData().value<void*>();
                itemNS = static_cast<QExtendedStandardItem*>(p);
                nodeNS = itemNS->data(Qt::UserRole).value<SimulationNodeClass*>();

                //! see: SimulationNodeClass* nodeFactory::nodeFromScratch Ref. [1]
                vecLoc = nodeNS->getPropertyValue<std::vector<GeometryTag>>("Tags");

                QVariant data;
                data.setValue(vecLoc);
                Property prop_masterTags("Tags master",data,Property::PropertyGroup_Scope);
                myCurNode->replaceProperty("Tags master",prop_masterTags);

                //! -----------------------------------------------------------------
                //! [2] rebuild the pairs (parent shape, child shape) for the Slave
                //! -----------------------------------------------------------------
                p = itemSlave->data(Qt::UserRole).value<Property>().getData().value<void*>();
                itemNS = static_cast<QExtendedStandardItem*>(p);
                nodeNS = itemNS->data(Qt::UserRole).value<SimulationNodeClass*>();

                //! see: SimulationNodeClass* nodeFactory::nodeFromScratch Ref. [1]
                vecLoc = nodeNS->getPropertyValue<std::vector<GeometryTag>>("Tags");
                data.setValue(vecLoc);
                Property prop_slaveTags("Tags slave",data,Property::PropertyGroup_Scope);
                myCurNode->replaceProperty("Tags slave",prop_slaveTags);
            }
        }
        break;

        case Property::ScopingMethod_RemotePoint:
        {
            cout<<"DetailViewer::updateTags()->____UPDATING TAGS BY REMOTE POINT____"<<endl;

            //! --------------------------------
            //! scope defined by "Remote point"
            //! --------------------------------
            void *p = myCurNode->getPropertyItem("Remote points")->data(Qt::UserRole).value<Property>().getData().value<void*>();
            QStandardItem *curItemRemotePointInItem = (QStandardItem*) p;
            std::vector<GeometryTag> vecLoc =curItemRemotePointInItem->data(Qt::UserRole).value<SimulationNodeClass*>()->getPropertyValue<std::vector<GeometryTag>>("Tags");
            QVariant data;
            data.setValue(vecLoc);
            Property prop("Tags",data,Property::PropertyGroup_Scope);
            myCurNode->replaceProperty("Tags",prop);
        }
            break;
        }
    }
    else
    {
        //! --------------------------------------------------------------
        //! for the "Named selection" the scope is always defined through
        //! a direct geometry selection
        //! retrieve the scope from "Geometry" or from "Location"
        //! --------------------------------------------------------------
        QVariant data;
        switch(nodeType)
        {
        case SimulationNodeClass::nodeType_namedSelectionGeometry:
        {
            cout<<"DetailViewer::updateTags()->____handling \"Geometry\"____"<<endl;
            std::vector<GeometryTag> vecLoc = myCurNode->getPropertyValue<std::vector<GeometryTag>>("Geometry");
            data.setValue(vecLoc);
            Property prop_tags("Tags",data,Property::PropertyGroup_Scope);
            myCurNode->replaceProperty("Tags", prop_tags);
        }
            break;

        case SimulationNodeClass::nodeType_coordinateSystem:
        {
            cout<<"DetailViewer::updateTags()->____handling \"Coordinate system\"____"<<endl;
            QExtendedStandardItem *item = myCurNode->getPropertyItem("Geometry");
            if(item==Q_NULLPTR)
            {
                item = myCurNode->getPropertyItem("Location");
            }
            std::vector<GeometryTag> vecLoc = item->data(Qt::UserRole).value<Property>().getData().value<std::vector<GeometryTag>>();
            data.setValue(vecLoc);
            Property prop_tags("Tags",data,Property::PropertyGroup_Origin);
            myCurNode->replaceProperty("Tags", prop_tags);
        }
            break;
        }
    }

    //! ---------------------------------------------------------------------
    //! after updating the "Tags", if the item is a "Force", "Remote force",
    //! "Remote displacement", "Moment" update also the "Reference node"
    //! ---------------------------------------------------------------------
    //SimulationNodeClass::nodeType nodeType = myCurNode->getType();
    if(nodeType == SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteForce ||
            nodeType == SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteDisplacement ||
            nodeType == SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteRotation ||
            nodeType == SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Force ||
            nodeType == SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Moment)
    {
        SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
        simulationDataBase *sdb = sm->getDataBase();
        std::vector<GeometryTag> vecLoc = myCurNode->getPropertyValue<std::vector<GeometryTag>>("Tags");
        QList<double> newReferencePoint = GeomToolsClass::calculateCentroid(sdb,vecLoc);
        //cout<<"DetailViewer::updateTags()->____updating reference point ("<<newReferencePoint.at(0)<<", "<<newReferencePoint.at(1)<<", "<<newReferencePoint.at(2)<<")____"<<endl;
        QVariant data;
        data.setValue(newReferencePoint.toVector());
        Property prop_referencePoint("Reference point",data,Property::PropertyGroup_Hidden);
        myCurNode->replaceProperty("Reference point",prop_referencePoint);
    }

    emit requestChangeColor();

    SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
    markerBuilder::addMarker(myCurNode,sm->getDataBase());

    //myCurNode->getModel()->blockSignals(false);

    cout<<"DetailViewer::updateTags()->____update tags exiting____"<<endl;
}

//! ------------------------------------
//! function: handleNumberOfStepChanged
//! details:
//! ------------------------------------
void DetailViewer::handleNumberOfStepChanged()
{
    emit requestResizeTabularData();
}

//! ------------------------------------
//! function: handleAnalysisTypeChanged
//! details:
//! ------------------------------------
void DetailViewer::handleAnalysisTypeChanged()
{
    Property::analysisType analysisType = myCurNode->getPropertyValue<Property::analysisType>("Analysis type");
    int tableRow = myCurNode->getPropertyValue<int>("Current step number");
    CustomTableModel *tabData = myCurNode->getTabularDataModel();
    QVariant data;
    data.setValue(analysisType);
    tabData->setDataRC(data,tableRow,TABULAR_DATA_ANALYSIS_TYPE_COLUMN,Qt::EditRole);
}

//! ---------------------------------------
//! function: handleTimeIntegrationChanged
//! details:
//! ---------------------------------------
void DetailViewer::handleTimeIntegrationChanged()
{
    Property::timeIntegration timeIntegration = myCurNode->getPropertyValue<Property::timeIntegration>("Static/Transient");
    int tableRow = myCurNode->getPropertyValue<int>("Current step number");
    QVariant data;
    data.setValue(timeIntegration);
    CustomTableModel *tabData = myCurNode->getTabularDataModel();
    tabData->setDataRC(data,tableRow,TABULAR_DATA_TIME_INTEGRATION_COLUMN,Qt::EditRole);
}

//! -----------------------------------
//! function: handleStepEndTimeChanged
//! details:
//! -----------------------------------
void DetailViewer::handleStepEndTimeChanged()
{
    //! the value of the "Step end time"
    double newStepEndTime = myCurNode->getPropertyValue<double>("Step end time");
    int tableRow = myCurNode->getPropertyValue<int>("Current step number");
    CustomTableModel *tabData = myCurNode->getTabularDataModel();
    QVariant data;
    data.setValue(newStepEndTime);
    tabData->setDataRC(data,tableRow,TABULAR_DATA_STEP_END_TIME_COLUMN,Qt::EditRole);
    //tabData->setDataRC(newStepEndTime,tableRow,TABULAR_DATA_STEP_END_TIME_COLUMN,Qt::EditRole);
}

//! ----------------------------------
//! function: handleSolverTypeChanged
//! details:
//! ----------------------------------
void DetailViewer::handleSolverTypeChanged()
{
    Property::solverType theSolverType = myCurNode->getPropertyValue<Property::solverType>("Solver type");
    int row = myCurNode->getPropertyValue<int>("Current step number");
    QVariant data;
    data.setValue(theSolverType);
    myCurNode->getTabularDataModel()->setDataRC(data,row,TABULAR_DATA_SOLVER_TYPE_COLUMN);
}

//! -----------------------------------------
//! function: handleCurrentStepNumberChanged
//! details:
//! -----------------------------------------
void DetailViewer::handleCurrentStepNumberChanged()
{
    cout<<"DetailViewer::handleCurrentStepNumberChanged()->____function called____"<<endl;
    QVariant data,data1;

    CustomTableModel *tabModel = myCurNode->getTabularDataModel();
    int currentStepNumber = myCurNode->getPropertyValue<int>("Current step number");

    //! ---------------------
    //! time division:
    //! "Step end time"
    //! "Auto time stepping"
    //! "Solver type"
    //! ---------------------
    double stepEndTime = tabModel->dataRC(currentStepNumber,1,Qt::EditRole).toDouble();
    QVector<int> N = tabModel->dataRC(currentStepNumber,3,Qt::EditRole).value<QVector<int>>();
    Property::solverType theSolverType = tabModel->dataRC(currentStepNumber,TABULAR_DATA_SOLVER_TYPE_COLUMN,Qt::EditRole).value<Property::solverType>();

    Property::timeIntegration timeIntegration = tabModel->dataRC(currentStepNumber,TABULAR_DATA_TIME_INTEGRATION_COLUMN,Qt::EditRole).value<Property::timeIntegration>();
    data.setValue(timeIntegration);
    myCurNode->replaceProperty("Static/Transient",Property("Static/Transient",data,Property::PropertyGroup_StepControls));

    //! -----------------
    //! Field parameters
    //! -----------------
    QVector<double> fieldParameters = tabModel->dataRC(currentStepNumber,4,Qt::EditRole).value<QVector<double>>();

    //! "Flux convergence" (5-th column of the table)
    int fluxConvergence = tabModel->dataRC(currentStepNumber,5,Qt::EditRole).toInt();

    //! "--Value" (or "--q_alpha_u")
    double Value = fieldParameters.at(3);

    //! "--q_alpha_0"
    double q_alpha_0 = fieldParameters.at(2);

    //! "--R_alpha_n"
    double R_alpha_n = fieldParameters.at(0);

    //! "--R_alpha_P"
    double R_alpha_P = fieldParameters.at(4);

    //! "--R_alpha_l"
    double R_alpha_l = fieldParameters.at(7);

    //! "--epsilon_alpha"
    double epsilon_alpha = fieldParameters.at(5);

    //! "Solution convergence" (6-th column of the table)
    int solutionConvergence = tabModel->dataRC(currentStepNumber,6,Qt::EditRole).toInt();

    //! "--C_alpha_n"
    double C_alpha_n = fieldParameters.at(1);

    //! "--C_alpha_epsilon"
    double C_alpha_epsilon = fieldParameters.at(6);

    //! -------------------------------------
    //! remove the "Step end time" item
    //! remove the "Solver type" item
    //! remove the "Auto time stepping" item
    //! -------------------------------------
    myCurNode->removeProperty("Step end time");
    myCurNode->removeProperty("Solver type");
    myCurNode->removeProperty("Auto time stepping");

    //! ---------------------------------
    //! rebuild the "Step end time" item
    //! ---------------------------------
    data1.setValue(stepEndTime);
    Property prop_StepEndTime("Step end time",data1,Property::PropertyGroup_StepControls);
    myCurNode->addProperty(prop_StepEndTime,2);

    //! --------------------------------------
    //! rebuild the "Auto time stepping" item
    //! --------------------------------------
    if(N[3]==0)
    {
        myCurNode->removeProperty("Initial substeps");
        myCurNode->removeProperty("Maximum substeps");
        myCurNode->removeProperty("Minimum substeps");
        myCurNode->removeProperty("Number of substeps");
        data1.setValue(Property::autoTimeStepping_ProgramControlled);
        Property prop_autoTimeStepping("Auto time stepping",data1,Property::PropertyGroup_StepControls);
        myCurNode->addProperty(prop_autoTimeStepping,3);
    }
    else if(N[3]==1)
    {
        myCurNode->removeProperty("Initial substeps");
        myCurNode->removeProperty("Maximum substeps");
        myCurNode->removeProperty("Minimum substeps");
        myCurNode->removeProperty("Number of substeps");
        data1.setValue(Property::autoTimeStepping_ON);
        Property prop_autoTimeStepping("Auto time stepping",data1,Property::PropertyGroup_StepControls);
        myCurNode->addProperty(prop_autoTimeStepping,3);

        data.setValue(N[0]);
        Property prop_NdivIni("Initial substeps",data,Property::PropertyGroup_StepControls);
        myCurNode->addProperty(prop_NdivIni);
        data.setValue(N[1]);
        Property prop_NdivMin("Minimum substeps",data,Property::PropertyGroup_StepControls);
        myCurNode->addProperty(prop_NdivMin);
        data.setValue(N[2]);
        Property prop_NdivMax("Maximum substeps",data,Property::PropertyGroup_StepControls);
        myCurNode->addProperty(prop_NdivMax);
    }
    else if(N[3]==2)
    {
        myCurNode->removeProperty("Initial substeps");
        myCurNode->removeProperty("Maximum substeps");
        myCurNode->removeProperty("Minimum substeps");
        myCurNode->removeProperty("Number of substeps");
        data1.setValue(Property::autoTimeStepping_OFF);
        Property prop_autoTimeStepping("Auto time stepping",data1,Property::PropertyGroup_StepControls);
        myCurNode->addProperty(prop_autoTimeStepping,3);

        data.setValue(N[0]);
        Property prop_Ndiv("Number of substeps",data,Property::PropertyGroup_StepControls);
        myCurNode->addProperty(prop_Ndiv);
    }

    //! -------------------------------
    //! rebuild the "Solver type" item
    //! -------------------------------
    data1.setValue(theSolverType);
    Property prop_solverType("Solver type",data1,Property::PropertyGroup_SolverControls);
    myCurNode->addProperty(prop_solverType);

    //! ---------------------------------
    //! remove the "Convergence criteria"
    //! ---------------------------------
    myCurNode->removeProperty("Flux convergence");
    myCurNode->removeProperty("--Value");
    myCurNode->removeProperty("--q_alpha_0");
    myCurNode->removeProperty("--R_alpha_n");
    myCurNode->removeProperty("--R_alpha_P");
    myCurNode->removeProperty("--R_alpha_l");
    myCurNode->removeProperty("--epsilon_alpha");
    myCurNode->removeProperty("Solution convergence");
    myCurNode->removeProperty("--C_alpha_n");
    myCurNode->removeProperty("--C_alpha_epsilon");

    //! ----------------------------------
    //! rebuild the "Convergence criteria"
    //! ----------------------------------
    data.setValue(fluxConvergence);
    Property prop_fluxConvergence("Flux convergence",data,Property::PropertyGroup_ConvergenceCriteria);
    myCurNode->addProperty(prop_fluxConvergence);

    data.setValue(Value);
    Property prop_Value("--Value",data,Property::PropertyGroup_ConvergenceCriteria);
    myCurNode->addProperty(prop_Value);

    data.setValue(q_alpha_0);
    Property prop_q_alpha_0("--q_alpha_0",data,Property::PropertyGroup_ConvergenceCriteria);
    myCurNode->addProperty(prop_q_alpha_0);

    data.setValue(R_alpha_n);
    Property prop_R_alpha_n("--R_alpha_n",data,Property::PropertyGroup_ConvergenceCriteria);
    myCurNode->addProperty(prop_R_alpha_n);

    data.setValue(R_alpha_P);
    Property prop_R_alpha_P("--R_alpha_P",data,Property::PropertyGroup_ConvergenceCriteria);
    myCurNode->addProperty(prop_R_alpha_P);

    data.setValue(R_alpha_l);
    Property prop_R_alpha_l("--R_alpha_l",data,Property::PropertyGroup_ConvergenceCriteria);
    myCurNode->addProperty(prop_R_alpha_l);

    data.setValue(epsilon_alpha);
    Property prop_epsilon_alpha("--epsilon_alpha",data,Property::PropertyGroup_ConvergenceCriteria);
    myCurNode->addProperty(prop_epsilon_alpha);

    data.setValue(solutionConvergence);
    Property prop_solutionConvergence("Solution convergence",data,Property::PropertyGroup_ConvergenceCriteria);
    myCurNode->addProperty(prop_solutionConvergence);

    data.setValue(C_alpha_n);
    Property prop_C_alpha_n("--C_alpha_n",data,Property::PropertyGroup_ConvergenceCriteria);
    myCurNode->addProperty(prop_C_alpha_n);

    data.setValue(C_alpha_epsilon);
    Property prop_C_alpha_epsilon("--C_alpha_epsilon",data,Property::PropertyGroup_ConvergenceCriteria);
    myCurNode->addProperty(prop_C_alpha_epsilon);    

    //! -------------------------------------------------------------------------------------------
    //! if "Flux convergence" is "Remove" or "Program controlled" hide all controls, else show all
    //! -------------------------------------------------------------------------------------------
    bool isHidden2;
    if(fluxConvergence==0 || fluxConvergence==1) isHidden2 = true; else isHidden2 = false;
    QStandardItem *theSeparator2 = myCurNode->getModel()->itemFromIndex(myCurNode->getPropertyItem("Flux convergence")->index().parent());
    for(int row=1; row<=6; row++) this->setRowHidden(row,theSeparator2->index(),isHidden2);

    //! -------------------------------------------------------------------------------------------
    //! if "Flux convergence" is "Remove" or "Program controlled" hide all controls, else show all
    //! -------------------------------------------------------------------------------------------
    bool isHidden3;
    if(solutionConvergence==0 || solutionConvergence==1) isHidden3 = true; else isHidden3 = false;
    QStandardItem *theSeparator3 = myCurNode->getModel()->itemFromIndex(myCurNode->getPropertyItem("Solution convergence")->index().parent());
    for(int row=8; row<=9; row++) this->setRowHidden(row,theSeparator3->index(),isHidden3);

    //! -------------------------------
    //! time incrementation parameters
    //! -------------------------------
    myCurNode->removeProperty("Time incrementation");
    myCurNode->removeProperty("I_0");
    myCurNode->removeProperty("I_R");
    myCurNode->removeProperty("I_P");
    myCurNode->removeProperty("I_C");
    myCurNode->removeProperty("I_L");
    myCurNode->removeProperty("I_G");
    myCurNode->removeProperty("I_S");
    myCurNode->removeProperty("I_A");
    myCurNode->removeProperty("I_J");
    myCurNode->removeProperty("I_T");

    int timeIncrementationType = tabModel->dataRC(currentStepNumber,8,Qt::EditRole).toInt();
    data.setValue(timeIncrementationType);
    Property prop_timeIncrementationType("Time incrementation",data,Property::PropertyGroup_TimeIncrementation);
    myCurNode->addProperty(prop_timeIncrementationType);

    QVector<int> timeIncrementationParameters = tabModel->dataRC(currentStepNumber,7,Qt::EditRole).value<QVector<int>>();

    int I_0 = timeIncrementationParameters.at(0);
    data.setValue(I_0);
    Property prop_I_0("I_0",data,Property::PropertyGroup_TimeIncrementation);
    myCurNode->addProperty(prop_I_0);

    int I_R = timeIncrementationParameters.at(1);
    data.setValue(I_R);
    Property prop_I_R("I_R",data,Property::PropertyGroup_TimeIncrementation);
    myCurNode->addProperty(prop_I_R);

    int I_P = timeIncrementationParameters.at(2);
    data.setValue(I_P);
    Property prop_I_P("I_P",data,Property::PropertyGroup_TimeIncrementation);
    myCurNode->addProperty(prop_I_P);

    int I_C = timeIncrementationParameters.at(3);
    data.setValue(I_C);
    Property prop_I_C("I_C",data,Property::PropertyGroup_TimeIncrementation);
    myCurNode->addProperty(prop_I_C);

    int I_L = timeIncrementationParameters.at(4);
    data.setValue(I_L);
    Property prop_I_L("I_L",data,Property::PropertyGroup_TimeIncrementation);
    myCurNode->addProperty(prop_I_L);

    int I_G = timeIncrementationParameters.at(5);
    data.setValue(I_G);
    Property prop_I_G("I_G",data,Property::PropertyGroup_TimeIncrementation);
    myCurNode->addProperty(prop_I_G);

    int I_S = timeIncrementationParameters.at(6);
    data.setValue(I_S);
    Property prop_I_S("I_S",data,Property::PropertyGroup_TimeIncrementation);
    myCurNode->addProperty(prop_I_S);

    int I_A = timeIncrementationParameters.at(7);
    data.setValue(I_A);
    Property prop_I_A("I_A",data,Property::PropertyGroup_TimeIncrementation);
    myCurNode->addProperty(prop_I_A);

    int I_J = timeIncrementationParameters.at(8);
    data.setValue(I_J);
    Property prop_I_J("I_J",data,Property::PropertyGroup_TimeIncrementation);
    myCurNode->addProperty(prop_I_J);

    int I_T = timeIncrementationParameters.at(9);
    data.setValue(I_T);
    Property prop_I_T("I_T",data,Property::PropertyGroup_TimeIncrementation);
    myCurNode->addProperty(prop_I_T);

    //! ----------------------------------------------------------------------
    //! if "Time incrementation" is "Custom" hide all controls, else show all
    //! ----------------------------------------------------------------------
    bool isHidden;
    if(timeIncrementationType==0) isHidden = true; else isHidden = false;
    QStandardItem *theSeparator = myCurNode->getModel()->itemFromIndex(myCurNode->getPropertyItem("Time incrementation")->index().parent());
    for(int row=1; row<theSeparator->rowCount(); row++) this->setRowHidden(row,theSeparator->index(),isHidden);

    //! ---------------------------------
    //! cutback factors and cutback type
    //! ---------------------------------
    int cutBackType = tabModel->dataRC(currentStepNumber,9,Qt::EditRole).toInt();

    data.setValue(cutBackType);
    Property prop_cutBackType("Cutback factors",data,Property::PropertyGroup_CutBack);
    myCurNode->removeProperty("Cutback factors");
    myCurNode->addProperty(prop_cutBackType);

    QVector<double> cutBackParameters = tabModel->dataRC(currentStepNumber,10,Qt::EditRole).value<QVector<double>>();

    double D_f = cutBackParameters.at(0);
    double D_C = cutBackParameters.at(1);
    double D_B = cutBackParameters.at(2);
    double D_A = cutBackParameters.at(3);
    double D_S = cutBackParameters.at(4);
    double D_H = cutBackParameters.at(5);
    double D_D = cutBackParameters.at(6);
    double W_G = cutBackParameters.at(7);

    myCurNode->removeProperty("D_f");
    myCurNode->removeProperty("D_C");
    myCurNode->removeProperty("D_B");
    myCurNode->removeProperty("D_A");
    myCurNode->removeProperty("D_S");
    myCurNode->removeProperty("D_H");
    myCurNode->removeProperty("D_D");
    myCurNode->removeProperty("W_G");

    data.setValue(D_f);
    Property property_D_f("D_f",data,Property::PropertyGroup_CutBack);
    myCurNode->addProperty(property_D_f);

    data.setValue(D_C);
    Property property_D_C("D_C",data,Property::PropertyGroup_CutBack);
    myCurNode->addProperty(property_D_C);

    data.setValue(D_B);
    Property property_D_B("D_B",data,Property::PropertyGroup_CutBack);
    myCurNode->addProperty(property_D_B);

    data.setValue(D_A);
    Property property_D_A("D_A",data,Property::PropertyGroup_CutBack);
    myCurNode->addProperty(property_D_A);

    data.setValue(D_S);
    Property property_D_S("D_S",data,Property::PropertyGroup_CutBack);
    myCurNode->addProperty(property_D_S);

    data.setValue(D_H);
    Property property_D_H("D_H",data,Property::PropertyGroup_CutBack);
    myCurNode->addProperty(property_D_H);

    data.setValue(D_D);
    Property property_D_D("D_D",data,Property::PropertyGroup_CutBack);
    myCurNode->addProperty(property_D_D);

    data.setValue(W_G);
    Property property_W_G("W_G",data,Property::PropertyGroup_CutBack);
    myCurNode->addProperty(property_W_G);

    //! ----------------------------------------------------------------------
    //! if "Cutback factors" is "Custom" hide all controls, else show all
    //! ----------------------------------------------------------------------
    bool isHidden1;
    if(cutBackType==0) isHidden1 = true; else isHidden1 = false;
    QStandardItem *theSeparator1 = myCurNode->getModel()->itemFromIndex(myCurNode->getPropertyItem("Cutback factors")->index().parent());
    for(int row=1; row<theSeparator1->rowCount(); row++) this->setRowHidden(row,theSeparator1->index(),isHidden);

    //! ---------------------------------------
    //! line search and line search parameters
    //! ---------------------------------------
    myCurNode->removeProperty("Line search");
    myCurNode->removeProperty("Min value");
    myCurNode->removeProperty("Max value");

    int lineSearch = tabModel->dataRC(currentStepNumber,11,Qt::EditRole).toInt();
    QVector<double> lineSearchParameters = tabModel->dataRC(currentStepNumber,12,Qt::EditRole).value<QVector<double>>();

    data.setValue(lineSearch);
    Property property_lineSearch("Line search",data,Property::PropertyGroup_LineSearch);
    myCurNode->addProperty(property_lineSearch);

    data.setValue(lineSearchParameters.at(0));
    Property property_minLineSearchValue("Min value",data,Property::PropertyGroup_LineSearch);
    myCurNode->addProperty(property_minLineSearchValue);

    data.setValue(lineSearchParameters.at(1));
    Property property_maxLineSearchValue("Max value",data,Property::PropertyGroup_LineSearch);
    myCurNode->addProperty(property_maxLineSearchValue);

    //! ------------------------------------------------------------------
    //! handle the visibility of the "Min value" and "Max value" controls
    //! ------------------------------------------------------------------
    bool isHidden4 = false;
    if(lineSearch==0) isHidden4 = true;
    QStandardItem *item = myCurNode->getPropertyItem("Min value");
    this->setRowHidden(item->index().row(),item->index().parent(),isHidden4);
    item = myCurNode->getPropertyItem("Max value");
    this->setRowHidden(item->index().row(),item->index().parent(),isHidden4);

    //! ----------------
    //! Output controls
    //! ----------------
    QVector<int> outputControls = tabModel->dataRC(currentStepNumber,TABULAR_DATA_OUPUT_CONTROLS_COLUMN,Qt::EditRole).value<QVector<int>>();
    data.setValue(outputControls.at(0));
    Property prop_stress("Stress",data,Property::PropertyGroup_OutputSettings);
    myCurNode->replaceProperty("Stress",prop_stress);

    data.setValue(outputControls.at(1));
    Property prop_strain("Strain",data,Property::PropertyGroup_OutputSettings);
    myCurNode->replaceProperty("Strain",prop_strain);

    data.setValue(outputControls.at(2));
    Property prop_reactionForces("Reaction forces",data,Property::PropertyGroup_OutputSettings);
    myCurNode->replaceProperty("Reaction forces",prop_reactionForces);

    data.setValue(outputControls.at(3));
    Property prop_contactData("Contact data",data,Property::PropertyGroup_OutputSettings);
    myCurNode->replaceProperty("Contact data",prop_contactData);
}

//! -----------------------------------------------------------------------
//! function: updateDetailViewerFromTabularData
//! details:  refresh the detail viewer following a change in tabular data
//! -----------------------------------------------------------------------
void DetailViewer::updateDetailViewerFromTabularData(QModelIndex topLeftIndex, QModelIndex bottomRightIndex, QVector<int> roles)
{
    //! ---------------
    //! metodo critico
    //! ---------------
    Q_UNUSED(roles)
    Q_UNUSED(bottomRightIndex)

    static int i;
    cout<<"DetailViewer::updateDetailViewerFromTabularData()->____function called "<<i++<<"____"<<endl;

    CustomTableModel *tabularDataModel;

    //! -----------------------------------
    //! retrieve the "Current step number"
    //! -----------------------------------
    SimulationNodeClass *nodeAnalysisSettings = Q_NULLPTR;
    if(myCurNode->isAnalysisRoot()) nodeAnalysisSettings = myCurModelIndex.child(0,0).data(Qt::UserRole).value<SimulationNodeClass*>();
    if(myCurNode->isAnalysisSettings()) nodeAnalysisSettings = myCurModelIndex.data(Qt::UserRole).value<SimulationNodeClass*>();
    if(myCurNode->isSimulationSetUpNode()) nodeAnalysisSettings = myCurModelIndex.parent().child(0,0).data(Qt::UserRole).value<SimulationNodeClass*>();
    if(myCurNode->isSolution()) nodeAnalysisSettings = myCurModelIndex.parent().child(0,0).data(Qt::UserRole).value<SimulationNodeClass*>();
    if(myCurNode->isSolutionInformation()) nodeAnalysisSettings = myCurModelIndex.parent().parent().child(0,0).data(Qt::UserRole).value<SimulationNodeClass*>();
    if(myCurNode->isAnalysisResult()) nodeAnalysisSettings = myCurModelIndex.parent().parent().child(0,0).data(Qt::UserRole).value<SimulationNodeClass*>();
    if(myCurNode->isChildSimulationSetUpNode()) nodeAnalysisSettings = myCurModelIndex.parent().parent().child(0,0).data(Qt::UserRole).value<SimulationNodeClass*>();
    if(myCurNode->isNephewSimulationSetUpNode()) nodeAnalysisSettings = myCurModelIndex.parent().parent().parent().child(0,0).data(Qt::UserRole).value<SimulationNodeClass*>();

    if(nodeAnalysisSettings==Q_NULLPTR)
    {
        cerr<<"@---------------------------------------------------------@"<<endl;
        cerr<<"@- Analysis settings not found- check interface behavior -@"<<endl;
        cerr<<"@---------------------------------------------------------@"<<endl;
        return;
    }

    int currentStepNumber = nodeAnalysisSettings->getPropertyValue<int>("Current step number");
    tabularDataModel = nodeAnalysisSettings->getTabularDataModel();

    //! -------------------------------------------------------------------------------------------------------------------------
    //! Once changed, the node model emits the Qt signal "itemChanged()", which is connected with the SLOT "handleItemChange()".
    //! That SLOT function modifies the tabular data, putting into the table the values set by the DetailViewer controls.
    //! Instead the function defined in this part of code modifies the DetailViewer content (i.e. the node model) using the
    //! tabular data. In order to avoid the ping-pong effect first disconnect "itemChanged()" from "handleItemChange()", then,
    //! after the changes have been applied, reconnect (this is done at the end [*])
    //! -------------------------------------------------------------------------------------------------------------------------
    SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
    //QStandardItemModel *nodeModel = myCurNode->getModel();
    //disconnect(nodeModel,SIGNAL(itemChanged(QStandardItem*)),sm,SLOT(handleItemChange(QStandardItem*)));

    myCurNode->getModel()->blockSignals(true);

    int row = topLeftIndex.row();
    if(row == currentStepNumber)
    {
        SimulationNodeClass::nodeType nodeType = myCurNode->getType();
        switch(nodeType)
        {
        case SimulationNodeClass::nodeType_structuralAnalysisSettings:
        case SimulationNodeClass::nodeType_thermalAnalysisSettings:
        case SimulationNodeClass::nodeType_combinedAnalysisSettings:
        {
            int column = topLeftIndex.column();
            QVariant data = tabularDataModel->dataRC(row,column);
            switch(column)
            {
            case TABULAR_DATA_STEP_END_TIME_COLUMN:
            {
                double stepEndTime = data.toDouble();
                data.setValue(stepEndTime);
                myCurNode->replaceProperty("Step end time",Property("Step end time",data,Property::PropertyGroup_StepControls));
            }
                break;

            case TABULAR_DATA_SOLVER_TYPE_COLUMN:
            {
                Property::solverType theSolverType = data.value<Property::solverType>();
                data.setValue(theSolverType);
                myCurNode->replaceProperty("Solver type",Property("Solver type",data,Property::PropertyGroup_SolverControls));
            }
                break;

            case TABULAR_DATA_ANALYSIS_TYPE_COLUMN:
            {
                int analysisType = data.toInt();
                data.setValue(analysisType);
                myCurNode->replaceProperty("Analysis type",Property("Analysis type",data,Property::PropertyGroup_StepControls));
            }
                break;

            case TABULAR_DATA_TIME_INTEGRATION_COLUMN:
            {
                Property::timeIntegration timeIntegration = data.value<Property::timeIntegration>();
                data.setValue(timeIntegration);
                myCurNode->replaceProperty("Static/Transient",Property("Static/Transient",data,Property::PropertyGroup_StepControls));
            }
                break;
            }
        }
            break;

        case SimulationNodeClass::nodeType_structuralAnalysisBoltPretension:
        {
            cout<<"DetailViewer::updateDetailViewerFromTabularData()->____updating \"Bolt pretension\"____"<<endl;

            int SC = mainTreeTools::calculateStartColumn(sm->myTreeView);

            int col_boltStatus = SC;
            int col_boltLoad = SC+1;
            int col_boltAdjustment = SC+2;

            //cout<<"DetailViewer::updateDetailViewerFromTabularData()->____using data in columns {"<<SC<<", "<<SC+1<<", "<<SC+2<<"}____"<<endl;
            QVariant boltStatus = tabularDataModel->dataRC(row,col_boltStatus,Qt::EditRole);
            //if(boltStatus.canConvert<Property::defineBy>()) exit(1);
            QVariant boltLoad = tabularDataModel->dataRC(row,col_boltLoad,Qt::EditRole);
            //if(boltLoad.canConvert<Property::defineBy>()) exit(2);
            QVariant boltAdjustment = tabularDataModel->dataRC(row,col_boltAdjustment,Qt::EditRole);
            //if(boltAdjustment.canConvert<Property::defineBy>()) exit(3);

            //! -----------------------------------
            //! define the property to be replaced
            //! -----------------------------------
            Property prop_boltStatus("Bolt status",boltStatus,Property::PropertyGroup_Definition);
            Property prop_boltLoad("Load",boltLoad,Property::PropertyGroup_Definition);
            Property prop_boltAdjustment("Adjustment",boltAdjustment,Property::PropertyGroup_Definition);

            myCurNode->replaceProperty("Bolt status",prop_boltStatus);
            myCurNode->replaceProperty("Load",prop_boltLoad);
            myCurNode->replaceProperty("Adjustment",prop_boltAdjustment);

            //! ----------------------------------------------------------------
            //! handle the visibility of the properties "Load" and "Adjustment"
            //! ----------------------------------------------------------------
            QExtendedStandardItem *itemLoad = myCurNode->getPropertyItem("Load");
            QModelIndex indexLoad = itemLoad->index();
            QExtendedStandardItem *itemAdjustment = myCurNode->getPropertyItem("Adjustment");
            QModelIndex indexAdjustment = itemAdjustment->index();

            switch(boltStatus.value<Property::boltStatusDefinedBy>())
            {
            case Property::boltStatusDefinedBy_load:
            {
                //! hide the control for the "Adjustment" property and show the control for the "Load" property
                cout<<"DetailViewer::updateDetailViewerFromTabularData()->____bolt status defined by \"Load\"____"<<endl;
                if(!this->isRowHidden(indexAdjustment.row(),indexAdjustment.parent())) this->setRowHidden(indexAdjustment.row(),indexAdjustment.parent(),true);
                if(this->isRowHidden(indexLoad.row(),indexLoad.parent())) this->setRowHidden(indexLoad.row(),indexLoad.parent(),false);
            }
                break;

            case Property::boltStatusDefinedBy_adjustment:
            {
                //! hide the control for the "Load" property and show the control for the "Adjustment" property
                cout<<"DetailViewer::updateDetailViewerFromTabularData()->____bolt status defined by \"Adjustment\"____"<<endl;
                if(this->isRowHidden(indexAdjustment.row(),indexAdjustment.parent())) this->setRowHidden(indexAdjustment.row(),indexAdjustment.parent(),false);
                if(!this->isRowHidden(indexLoad.row(),indexLoad.parent())) this->setRowHidden(indexLoad.row(),indexLoad.parent(),true);
            }
                break;

            case Property::boltStatusDefinedBy_open:
            case Property::boltStatusDefinedBy_lock:
            {
                //! hide the controls for the "Load" and "Adjustment" properties
                cout<<"DetailViewer::updateDetailViewerFromTabularData()->____bolt status defined by \"Open\" or \"Lock\"____"<<endl;
                if(!this->isRowHidden(indexAdjustment.row(),indexAdjustment.parent())) this->setRowHidden(indexAdjustment.row(),indexAdjustment.parent(),true);
                if(!this->isRowHidden(indexLoad.row(),indexLoad.parent())) this->setRowHidden(indexLoad.row(),indexLoad.parent(),true);
            }
                break;
            }
        }
            break;

        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Force:
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteForce:
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Moment:
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RotationalVelocity:
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Acceleration:
        {
            cout<<"____function called for a vectorial quantity____"<<endl;
            Property::defineBy theDefineBy = myCurNode->getPropertyValue<Property::defineBy>("Define by");
            switch(theDefineBy)
            {
            case Property::defineBy_normal:
            case Property::defineBy_vector:
            {
                cout<<"____property defined by direction: updating magnitude____"<<endl;

                //! diagnostic - can be removed
                SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
                int col = mainTreeTools::calculateStartColumn(sm->myTreeView);

                int row = currentStepNumber;
                double magnitude = tabularDataModel->dataRC(row,col,Qt::EditRole).toDouble();
                cout<<"____step nr: "<<row<<" tab col= "<<col<<" read val: "<<magnitude<<"____"<<endl;
                //! end diagnostic

                QVariant data;
                data.setValue(Property::loadDefinition_tabularData);
                Property prop_magnitude("Magnitude",data,Property::PropertyGroup_Definition);
                myCurNode->replaceProperty("Magnitude",prop_magnitude);
            }
                break;

            case Property::defineBy_components:
            {
                cout<<"____property defined by components: updating 3 values____"<<endl;

                QList<QString> propNames;
                propNames<<"X component"<<"Y component"<<"Z component";
                for(int i=0; i<3; i++)
                {
                    //! ----------------------------
                    //! diagnostic - can be removed
                    //! ----------------------------
                    SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
                    int col = mainTreeTools::calculateStartColumn(sm->myTreeView)+i;
                    int row = currentStepNumber;
                    double componentValue = tabularDataModel->dataRC(row,col,Qt::EditRole).toDouble();
                    cout<<"____step nr: "<<row<<" tab col= "<<col<<" "<<propNames.at(i).toStdString()<<" read val: "<<componentValue<<"____"<<endl;
                    //! ---------------
                    //! end diagnostic
                    //! ---------------

                    QVariant data;
                    data.setValue(Property::loadDefinition_tabularData);
                    QString propName = propNames.at(i);
                    Property prop_component(propName,data,Property::PropertyGroup_Definition);
                    myCurNode->replaceProperty(propName,prop_component);
                }
            }
                break;
            }

            //! ---------------------------------------------------------------------------
            //! the following for updating the marker when the data changes occur in table
            //! ---------------------------------------------------------------------------
            if(nodeType==SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Moment) this->handleMomentChanged();
            if(nodeType==SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Acceleration) this->handleAccelerationChanged();
            // add here the other cases...
        }
            break;

        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Displacement:
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteDisplacement:
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteRotation:
        {
            cout<<"____function called for \"Displacement\"____"<<endl;
            Property::defineBy theDefineBy = myCurNode->getPropertyItem("Define by")->data(Qt::UserRole).value<Property>().getData().value<Property::defineBy>();
            switch(theDefineBy)
            {
            case Property::defineBy_vector:
            case Property::defineBy_normal:
            {
                cout<<"____property defined by direction: updating magnitude____"<<endl;

                //! diagnostic - can be removed
                SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
                int col = mainTreeTools::calculateStartColumn(sm->myTreeView);

                int row = currentStepNumber;
                double magnitude = tabularDataModel->dataRC(row,col,Qt::EditRole).toDouble();
                cout<<"____step nr: "<<row<<" tab col= "<<col<<" read val: "<<magnitude<<"____"<<endl;
                //! end diagnostic

                QVariant data;
                data.setValue(Property::loadDefinition_tabularData);
                Property prop_magnitude("Magnitude",data,Property::PropertyGroup_Definition);
                myCurNode->replaceProperty("Magnitude",prop_magnitude);
            }
                break;

            case Property::defineBy_components:
            {
                cout<<"____property defined by components: updating values____"<<endl;

                //! -----------------------------------------------------------
                //! change the components: take into account the "Free" option
                //! -----------------------------------------------------------
                Property::loadDefinition theLoadDefinitionXcomponent = myCurNode->getPropertyValue<Property::loadDefinition>("X component");
                Property::loadDefinition theLoadDefinitionYcomponent = myCurNode->getPropertyValue<Property::loadDefinition>("Y component");
                Property::loadDefinition theLoadDefinitionZcomponent = myCurNode->getPropertyValue<Property::loadDefinition>("Z component");
                bool c0 = theLoadDefinitionXcomponent != Property::loadDefinition_free? true:false;
                bool c1 = theLoadDefinitionYcomponent != Property::loadDefinition_free? true:false;
                bool c2 = theLoadDefinitionZcomponent != Property::loadDefinition_free? true:false;

                //! --------------------------
                //! only one component active
                //! --------------------------
                if((c0==true && c1==false && c2 == false) || (c0==false && c1==true && c2 == false) || (c0==false && c1==false && c2 == true))
                {
                    QString propertyName;

                    if(c0 == true) propertyName = "X component";
                    if(c1 == true) propertyName = "Y component";
                    if(c2 == true) propertyName = "Z component";

                    //! diagnostic - can be removed
                    //int col = sm->calculateStartColumn();
                    SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
                    int col = mainTreeTools::calculateStartColumn(sm->myTreeView);
                    int row = currentStepNumber;
                    double componentValue = tabularDataModel->dataRC(row,col,Qt::EditRole).toDouble();
                    cout<<"____step nr: "<<row<<" tab col= "<<col<<" "<<propertyName.toStdString()<<" read val: "<<componentValue<<"____"<<endl;
                    //! end diagnostic

                    QVariant data;
                    data.setValue(Property::loadDefinition_tabularData);
                    Property prop_displacementComponent(propertyName,data,Property::PropertyGroup_Definition);
                    myCurNode->replaceProperty(propertyName,prop_displacementComponent);
                }
                //! --------------------------
                //! two components are active
                //! --------------------------
                if((c0==true && c1==true && c2 == false) || (c0==true && c1==false && c2 == true) || (c0==false && c1==true && c2 == true))
                {
                    if(c0==true && c1==true && c2 == false)
                    {
                        //! diagnostic - can be removed
                        int row = currentStepNumber;
                        //int col = sm->calculateStartColumn();
                        SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
                        int col = mainTreeTools::calculateStartColumn(sm->myTreeView);
                        double componentValue1 = tabularDataModel->dataRC(row,col,Qt::EditRole).toDouble();
                        double componentValue2 = tabularDataModel->dataRC(row,col+1,Qt::EditRole).toDouble();
                        cout<<"____step nr: "<<row<<" tab col= "<<col<<" "<<" X component read val: "<<componentValue1<<"____"<<endl;
                        cout<<"____step nr: "<<row<<" tab col= "<<col+1<<" "<<" Y component read val: "<<componentValue2<<"____"<<endl;
                        //! end diagnostic

                        QVariant data;
                        data.setValue(Property::loadDefinition_tabularData);
                        Property prop_displacementComponentX("X component",data,Property::PropertyGroup_Definition);
                        myCurNode->replaceProperty("X component",prop_displacementComponentX);
                        Property prop_displacementComponentY("Y component",data,Property::PropertyGroup_Definition);
                        myCurNode->replaceProperty("Y component",prop_displacementComponentY);
                    }
                    if(c0==true && c1==false && c2 == true)
                    {
                        //! diagnostic - can be removed
                        int row = currentStepNumber;
                        //int col = sm->calculateStartColumn();
                        SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));//bubi
                        int col = mainTreeTools::calculateStartColumn(sm->myTreeView);
                        double componentValue1 = tabularDataModel->dataRC(row,col,Qt::EditRole).toDouble();
                        double componentValue2 = tabularDataModel->dataRC(row,col+1,Qt::EditRole).toDouble();
                        cout<<"____step nr: "<<row<<" tab col= "<<col<<" "<<" X component read val: "<<componentValue1<<"____"<<endl;
                        cout<<"____step nr: "<<row<<" tab col= "<<col+1<<" "<<" Z component read val: "<<componentValue2<<"____"<<endl;
                        //! end diagnostic

                        QVariant data;
                        data.setValue(Property::loadDefinition_tabularData);
                        Property prop_displacementComponentX("X component",data,Property::PropertyGroup_Definition);
                        myCurNode->replaceProperty("X component",prop_displacementComponentX);
                        Property prop_displacementComponentZ("Z component",data,Property::PropertyGroup_Definition);
                        myCurNode->replaceProperty("Z component",prop_displacementComponentZ);
                    }
                    if(c0==false && c1==true && c2 == true)
                    {
                        //! diagnostic - can be removed
                        int row = currentStepNumber;
                        //int col = sm->calculateStartColumn();
                        SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));//bubi
                        int col = mainTreeTools::calculateStartColumn(sm->myTreeView);
                        double componentValue1 = tabularDataModel->dataRC(row,col,Qt::EditRole).toDouble();
                        double componentValue2 = tabularDataModel->dataRC(row,col+1,Qt::EditRole).toDouble();
                        cout<<"____step nr: "<<row<<" tab col= "<<col<<" "<<" Y component read val: "<<componentValue1<<"____"<<endl;
                        cout<<"____step nr: "<<row<<" tab col= "<<col+1<<" "<<" Z component read val: "<<componentValue2<<"____"<<endl;
                        //! end diagnostic

                        QVariant data;
                        data.setValue(Property::loadDefinition_tabularData);
                        Property prop_displacementComponentY("Y component",data,Property::PropertyGroup_Definition);
                        myCurNode->replaceProperty("Y component",prop_displacementComponentY);
                        Property prop_displacementComponentZ("Z component",data,Property::PropertyGroup_Definition);
                        myCurNode->replaceProperty("Z component",prop_displacementComponentZ);
                    }
                }
                //! ----------------------------
                //! three components are active
                //! ----------------------------
                if(c0==true && c1==true && c2==true)
                {
                    QList<QString> propNames;
                    propNames<<"X component"<<"Y component"<<"Z component";
                    for(int i=0; i<3; i++)
                    {
                        //! diagnostic - can be removed
                        //int col = sm->calculateStartColumn()+i;
                        SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));//bubi
                        int col = mainTreeTools::calculateStartColumn(sm->myTreeView)+1;
                        int row = currentStepNumber;
                        double componentValue = tabularDataModel->dataRC(row,col,Qt::EditRole).toDouble();
                        cout<<"____step nr: "<<row<<" tab col= "<<col<<" "<<propNames.at(i).toStdString()<<" read val: "<<componentValue<<"____"<<endl;
                        //! end diagnostic

                        QVariant data;
                        data.setValue(Property::loadDefinition_tabularData);
                        QString propName = propNames.at(i);
                        Property prop_component(propName,data,Property::PropertyGroup_Definition);
                        myCurNode->replaceProperty(propName,prop_component);
                    }
                }
            }
                break;
            }
        }
            break;
        case SimulationNodeClass::nodeType_structuralAnalysisThermalCondition:
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Pressure:
        case SimulationNodeClass::nodeType_thermalAnalysisConvection:
        case SimulationNodeClass::nodeType_thermalAnalysisTemperature:
        case SimulationNodeClass::nodeType_thermalAnalysisThermalFlow:
        case SimulationNodeClass::nodeType_thermalAnalysisThermalFlux:
        case SimulationNodeClass::nodeType_thermalAnalysisThermalPower:
        {
            //cout<<"DetailViewer::updateDetailViewerFromTabularData()->____function called for \"Pressure\"____"<<endl;
            QVariant data;
            data.setValue(Property::loadDefinition_tabularData);
            Property prop_magnitude("Magnitude",data,Property::PropertyGroup_Definition);
            myCurNode->replaceProperty("Magnitude",prop_magnitude);
        }
            break;

        case SimulationNodeClass::nodeType_modelChange:
        {
            //cout<<"DetailViewer::updateDetailViewerFromTabularData()->____function called for \"Model change\"____"<<endl;

            SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
            int col = mainTreeTools::calculateStartColumn(sm->myTreeView);
            int row = currentStepNumber;
            Property::modelChangeActivationStatus val = tabularDataModel->dataRC(row,col,Qt::EditRole).value<Property::modelChangeActivationStatus>();

            QVariant data;
            data.setValue(val);
            Property prop_activationStatus("Activation status",data,Property::PropertyGroup_Definition);
            myCurNode->replaceProperty("Activation status",prop_activationStatus);
        }
            break;
        }
    }

    //! -----------------
    //! [*] reconnection
    //! -----------------
    //if(myCurNode->getType()!=SimulationNodeClass::nodeType_structuralAnalysisSettings &&
    //        myCurNode->getType()!=SimulationNodeClass::nodeType_thermalAnalysisSettings)
    //    connect(nodeModel,SIGNAL(itemChanged(QStandardItem*)),sm,SLOT(handleItemChange(QStandardItem*)));

    myCurNode->getModel()->blockSignals(false);

    cout<<"DetailViewer::updateDetailViewerFromTabularData()->____exiting____"<<endl;
}

//! --------------------------------
//! function: selectionHasChanged()
//! details:
//! --------------------------------
void DetailViewer::selectionHasChanged()
{
    //cout<<"DetailViewer::selectionHasChanged()->____function called____"<<endl;
    emit requestHandleSelectionChanged();
}

//! ---------------------------------
//! function: handleDefineBy_Changed
//! details:
//! ---------------------------------
void DetailViewer::handleDefineBy_Changed()
{
    cout<<"DetailViewer::handleDefineBy_Changed()->____function called____"<<endl;
    QExtendedStandardItem *itemDefineBy_ = myCurNode->getPropertyItem("Define by ");
    if(itemDefineBy_!=Q_NULLPTR)     //! it's safe
    {
        QExtendedStandardItem *itemGeometry = myCurNode->getPropertyItem("Geometry");
        QExtendedStandardItem *itemLocation = myCurNode->getPropertyItem("Location");
        Property::defineBy theDefineBy = itemDefineBy_->data(Qt::UserRole).value<Property>().getData().value<Property::defineBy>();

        if(theDefineBy==Property::defineBy_geometrySelection)
        {
            if(itemGeometry==Q_NULLPTR)
            {
                //! add
                QVariant data;
                std::vector<GeometryTag> vecLoc;
                data.setValue(vecLoc);
                Property prop_geometry("Geometry",data,Property::PropertyGroup_Origin);
                myCurNode->addProperty(prop_geometry);
            }
            if(itemLocation!=Q_NULLPTR)
            {
                //! remove
                myCurNode->getModel()->removeRow(itemLocation->index().row(),itemLocation->parent()->index());
            }
        }
        else    //! theDefineBy = defineBy_globalCoordinates
        {
            if(itemGeometry!=Q_NULLPTR)
            {
                //! remove
                myCurNode->getModel()->removeRow(itemGeometry->index().row(),itemGeometry->parent()->index());
            }
            if(itemLocation==Q_NULLPTR)
            {
                //! add
                QVariant data;
                std::vector<GeometryTag> vecLoc;
                data.setValue(vecLoc);
                Property prop_location("Location",data,Property::PropertyGroup_Origin);
                myCurNode->addProperty(prop_location);
            }
        }
    }
}

//! ----------------------------------------------------------------
//! function: handleOriginAndDirectionChanged
//! details:  handle the case in which the origin and the direction
//!           of the CS have been changed by a geometry selection
//! ----------------------------------------------------------------
void DetailViewer::handleOriginAndDirectionChanged()
{
    cout<<"DetailViewer::handleOriginAndDirectionChanged()->____function called____"<<endl;

    //! -------------------------------------------
    //! this function needs the geometry data base
    //! -------------------------------------------
    SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
    geometryDataBase *gdb = sm->getDataBase();

    QVector<double> Origin;
    QVector<QVector<double>> directionalData;

    QExtendedStandardItem *itemTags = myCurNode->getPropertyItem("Tags");
    std::vector<GeometryTag> vecLoc = itemTags->data(Qt::UserRole).value<Property>().getData().value<std::vector<GeometryTag>>();

    if(vecLoc.size()!=0)
    {
        //! block the node signals
        this->connectToSimulationManager(true);

        //! ------------------------------------------
        //! use the last location in case of multiple
        //! ------------------------------------------
        int bodyIndex = vecLoc[vecLoc.size()-1].parentShapeNr;
        int subShapeIndex =  vecLoc[vecLoc.size()-1].subTopNr;

        TopoDS_Face theFace = TopoDS::Face(gdb->MapOfBodyTopologyMap.value(bodyIndex).faceMap.FindKey(subShapeIndex));
        GeomAbs_SurfaceType surfaceType;
        GeomToolsClass::getFaceType(theFace,surfaceType);

        gp_Pnt P;
        gp_Ax1 theDirection;        //! main axis (Z-axis) of the local coordinate system

        switch(surfaceType)
        {
        case GeomAbs_Plane:
        {
            GeomToolsClass::getPlanarFaceInfo(theFace,theDirection,P);
            cout<<"DetailViewer::handleOriginAndDirectionChanged()->____the new origin ("<<P.X()<<" ,"<<P.Y()<<", "<<P.Z()<<")____"<<endl;
        }
            break;
            //! other cases ... to do
        }
        //! -----------------------------------------------------------------
        //! create a coordinate system through origin and direction (Z-axis)
        //! the X and Y axes are automatically set
        //! -----------------------------------------------------------------
        gp_Ax2 theCoordinateSystem(P,theDirection.Direction());

        //! ----------------------------
        //! update the directional data
        //! ----------------------------
        gp_Dir ex = gp_Ax2().XDirection();
        gp_Dir ey = gp_Ax2().YDirection();
        gp_Dir ez = gp_Ax2().Direction();

        gp_Dir ex_loc = theCoordinateSystem.XDirection();
        gp_Dir ey_loc = theCoordinateSystem.YDirection();
        gp_Dir ez_loc = theCoordinateSystem.Direction();

        double A11,A12,A13,A21,A22,A23,A31,A32,A33;
        A11 = ex_loc.Dot(ex); A12 = ex_loc.Dot(ey); A13 = ex_loc.Dot(ez);
        A21 = ey_loc.Dot(ex); A22 = ey_loc.Dot(ey); A23 = ey_loc.Dot(ez);
        A31 = ez_loc.Dot(ex); A32 = ez_loc.Dot(ey); A33 = ez_loc.Dot(ez);

        QVector<double> X_AxisData, Y_AxisData, Z_AxisData;
        X_AxisData.append(A11); X_AxisData.append(A12); X_AxisData.append(A13);
        Y_AxisData.append(A21); Y_AxisData.append(A22); Y_AxisData.append(A23);
        Z_AxisData.append(A31); Z_AxisData.append(A32); Z_AxisData.append(A33);

        QVariant datax,datay,dataz,data;
        datax.setValue(X_AxisData);
        datay.setValue(Y_AxisData);
        dataz.setValue(Z_AxisData);
        Property prop_datax("X axis data",datax,Property::PropertyGroup_DirectionalData);
        Property prop_datay("Y axis data",datay,Property::PropertyGroup_DirectionalData);
        Property prop_dataz("Z axis data",dataz,Property::PropertyGroup_DirectionalData);

        directionalData.push_back(X_AxisData);
        directionalData.push_back(Y_AxisData);
        directionalData.push_back(Z_AxisData);

        data.setValue(prop_datax);
        myCurNode->replaceProperty("X axis data",prop_datax);
        myCurNode->replaceProperty("Y axis data",prop_datay);
        myCurNode->replaceProperty("Z axis data",prop_dataz);

        //! ------------------
        //! update the origin
        //! ------------------
        Property prop_Xorigin("Origin X",P.X(),Property::PropertyGroup_Origin);
        Property prop_Yorigin("Origin Y",P.Y(),Property::PropertyGroup_Origin);
        Property prop_Zorigin("Origin Z",P.Z(),Property::PropertyGroup_Origin);

        Origin.push_back(P.X());
        Origin.push_back(P.Y());
        Origin.push_back(P.Z());

        myCurNode->replaceProperty("Origin X",prop_Xorigin);
        myCurNode->replaceProperty("Origin Y",prop_Yorigin);
        myCurNode->replaceProperty("Origin Z",prop_Zorigin);
        //! end update the origin and the directional data

        //! copy the new origin coordinates and the new directional data
        //! into the base origin cell and base directional data cell
        this->updateBaseData();

        //! unblock the node signals
        this->connectToSimulationManager(true);

        //! request the OCC viewer to display the new trihedron
        emit requestDisplayTrihedron(Origin,directionalData);
    }
}

//! ------------------------------------------------------------------------
//! function: handleOriginChanged
//! details:  the origin of the CS has been changed by a geometry selection
//! ------------------------------------------------------------------------
void DetailViewer::handleOriginChanged()
{
    cout<<"DetailViewer::handleOriginChanged()->____function called____"<<endl;
    QVector<double> Origin;
    QExtendedStandardItem *itemScope = myCurNode->getPropertyItem("Location");
    if(itemScope==Q_NULLPTR) return;

    this->connectToSimulationManager(false);

    //! ------------------
    //! update the origin
    //! ------------------
    gp_Pnt P;
    std::vector<GeometryTag> vecLoc = itemScope->data(Qt::UserRole).value<Property>().getData().value<std::vector<GeometryTag>>();

    if(vecLoc.size()==0) return;

    SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
    geometryDataBase *gdb = sm->getDataBase();

    //! ------------------------------------------
    //! use the last location in case of multiple
    //! ------------------------------------------
    const GeometryTag &loc = vecLoc[vecLoc.size()-1];
    int bodyIndex = loc.parentShapeNr;
    int subShapeNr = loc.subTopNr;

    switch(loc.subShapeType)
    {
    case TopAbs_FACE:
    {
        TopoDS_Face theFace = TopoDS::Face(gdb->MapOfBodyTopologyMap.value(bodyIndex).faceMap.FindKey(subShapeNr));
        GeomAbs_SurfaceType surfaceType;
        GeomToolsClass::getFaceType(theFace,surfaceType);
        gp_Ax1 theDirectionTemp;

        switch(surfaceType)
        {
        case GeomAbs_Plane:
            GeomToolsClass::getPlanarFaceInfo(theFace,theDirectionTemp,P);
            break;

            //! other cases ... to do ...
        }
    }
        break;

    case TopAbs_EDGE:
    {
        TopoDS_Edge theEdge = TopoDS::Edge(gdb->MapOfBodyTopologyMap.value(bodyIndex).edgeMap.FindKey(subShapeNr));
        GeomAbs_CurveType curveType;
        GeomToolsClass::getCurveType(theEdge, curveType);
        switch(curveType)
        {
        case GeomAbs_Circle:
        {
            ;
        }
            break;

        }
    }
        break;

    case TopAbs_VERTEX:
    {
        TopoDS_Vertex V = TopoDS::Vertex(gdb->MapOfBodyTopologyMap.value(bodyIndex).edgeMap.FindKey(subShapeNr));
        P = BRep_Tool::Pnt(V);
    }
        break;
    }

    //! ----------------------------------------
    //! update "Origin X" "Origin Y" "Origin Z"
    //! ----------------------------------------
    Property prop_Xorigin("Origin X",P.X(),Property::PropertyGroup_Origin);
    Property prop_Yorigin("Origin Y",P.Y(),Property::PropertyGroup_Origin);
    Property prop_Zorigin("Origin Z",P.Z(),Property::PropertyGroup_Origin);

    Origin.push_back(P.X()); Origin.push_back(P.Y()); Origin.push_back(P.Z());

    QVariant data;
    data.setValue(prop_Xorigin);
    myCurNode->getPropertyItem(("Origin X"))->setData(data, Qt::UserRole);
    data.setValue(prop_Yorigin);
    myCurNode->getPropertyItem(("Origin Y"))->setData(data, Qt::UserRole);
    data.setValue(prop_Zorigin);
    myCurNode->getPropertyItem(("Origin Z"))->setData(data, Qt::UserRole);

    //! get the directional data
    QVector<QVector<double>> directionalData = this->getDirectionalData();

    //! copy the new origin coordinates and the new directional data
    //! into the base origin cell and base directional data cell
    this->updateBaseData();

    //! unlock signals
    this->connectToSimulationManager(true);

    //! request the OCC viewer to display the new trihedron
    emit requestDisplayTrihedron(Origin, directionalData);
}

//! ----------------------------------------------------
//! function: handleOriginChangedByValue()
//! details:  the origin of a CS has been changed by a
//!           direct specifications of the coordinates
//! ----------------------------------------------------
void DetailViewer::handleOriginChangedByValue()
{
    this->connectToSimulationManager(false);

    //! recover the coordinates of the origin
    QVector<double> Origin = this->getOrigin();

    //! recover the directional data (not modified here, already present in the node)
    QVector<QVector<double>> directionalData = this->getDirectionalData();

    //! copy the new origin coordinates and the new directional data
    //! into the base origin cell and base directional data cell
    this->updateBaseData();

    this->connectToSimulationManager(true);

    //! request the OCC viewer to display a new trihedron
    emit requestDisplayTrihedron(Origin, directionalData);
}

//! ------------------------------------------------------
//! function: handleRemotePointSystemOfReferenceChanged()
//! details:
//! ------------------------------------------------------
void DetailViewer::handleRemotePointSystemOfReferenceChanged()
{
    cout<<"handleRemotePointSystemOfReferenceChanged()->____function called____"<<endl;
    cout<<"handleRemotePointSystemOfReferenceChanged()->____calculating new coordinates____"<<endl;
    /*
    //! --------------------------------
    //! update the absolute coordinates
    //! --------------------------------
    void *p = myCurNode->getPropertyItem("Coordinate system")->data(Qt::UserRole).value<Property>().getData().value<void*>();
    QExtendedStandardItem *item = (QExtendedStandardItem*)p;
    SimulationNodeClass *CS = item->data(Qt::UserRole).value<SimulationNodeClass*>();
    if(CS->getType()!=SimulationNodeClass::nodeType_coordinateSystem_global)
    {
        QVector<double> CS_origin = CS->getPropertyValue<QVector<double>>("Base origin");
        QVector<double> X_axisData = CS->getPropertyValue<QVector<double>>("X axis data");
        QVector<double> Y_axisData = CS->getPropertyValue<QVector<double>>("Y axis data");
        QVector<double> Z_axisData = CS->getPropertyValue<QVector<double>>("Z axis data");

        cout<<"handleRemotePointSystemOfReferenceChanged()->____CS origin: ("<<CS_origin.at(0)<<", "<<CS_origin.at(1)<<", "<<CS_origin.at(2)<<")____"<<endl;
        cout<<"handleRemotePointSystemOfReferenceChanged()->____CS directional data____"<<endl;
        cout<<"____("<<X_axisData.at(0)<<", "<<X_axisData.at(1)<<", "<<X_axisData.at(2)<<")____"<<endl;
        cout<<"____("<<Y_axisData.at(0)<<", "<<Y_axisData.at(1)<<", "<<Y_axisData.at(2)<<")____"<<endl;
        cout<<"____("<<Z_axisData.at(0)<<", "<<Z_axisData.at(1)<<", "<<Z_axisData.at(2)<<")____"<<endl;

        double xO = CS_origin.at(0);
        double yO = CS_origin.at(1);
        double zO = CS_origin.at(2);

        double xP = myCurNode->getPropertyValue<double>("X coordinate");
        double yP = myCurNode->getPropertyValue<double>("Y coordinate");
        double zP = myCurNode->getPropertyValue<double>("Z coordinate");

        double xnew = xO + (xP-xO)*X_axisData.at(0) + (yP-yO)*Y_axisData.at(0) + (zP-zO)*Z_axisData.at(0);
        double ynew = yO + (xP-xO)*X_axisData.at(1) + (yP-yO)*Y_axisData.at(1) + (zP-zO)*Z_axisData.at(1);
        double znew = xO + (xP-xO)*X_axisData.at(2) + (yP-yO)*Y_axisData.at(2) + (zP-zO)*Z_axisData.at(2);

        QVariant data;
        data.setValue(xnew);
        Property prop_Xcoord("X abs coordinate",data,Property::PropertyGroup_Scope);

        data.setValue(ynew);
        Property prop_Ycoord("Y abs coordinate",data,Property::PropertyGroup_Scope);

        data.setValue(znew);
        Property prop_Zcoord("Z abs coordinate",data,Property::PropertyGroup_Scope);

        myCurNode->replaceProperty("X abs coordinate",prop_Xcoord);
        myCurNode->replaceProperty("Y abs coordinate",prop_Ycoord);
        myCurNode->replaceProperty("Z abs coordinate",prop_Zcoord);
    }
    else
    {
        double xP = myCurNode->getPropertyValue<double>("X coordinate");
        double yP = myCurNode->getPropertyValue<double>("Y coordinate");
        double zP = myCurNode->getPropertyValue<double>("Z coordinate");

        QVariant data;
        data.setValue(xP);
        Property prop_Xcoord("X abs coordinate",data,Property::PropertyGroup_Scope);

        data.setValue(yP);
        Property prop_Ycoord("Y abs coordinate",data,Property::PropertyGroup_Scope);

        data.setValue(zP);
        Property prop_Zcoord("Z abs coordinate",data,Property::PropertyGroup_Scope);

        myCurNode->replaceProperty("X abs coordinate",prop_Xcoord);
        myCurNode->replaceProperty("Y abs coordinate",prop_Ycoord);
        myCurNode->replaceProperty("Z abs coordinate",prop_Zcoord);
    }
    */
}

//! ---------------------------------------------
//! function: handleRemotePointChangedByLocation
//! details:
//! ---------------------------------------------
void DetailViewer::handleRemotePointChangedByLocation()
{
    cout<<"DetailViewer::handleRemotePointLocationChanged()->____function called____"<<endl;
    gp_Pnt P;
    gp_Ax1 theAxis;
    QExtendedStandardItem *itemLocation = myCurNode->getPropertyItem("Location");
    if(itemLocation!=Q_NULLPTR)
    {
        //! --------------------------------------------
        //! the updated point P: calculate the location
        //! using the last entry in "Location"
        //! --------------------------------------------
        const std::vector<GeometryTag> &vecLoc = itemLocation->data(Qt::UserRole).value<Property>().getData().value<std::vector<GeometryTag>>();
        if(vecLoc.size()!=0)
        {
            SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
            geometryDataBase *gdb = sm->getDataBase();

            if(vecLoc.size()==1)
            {
                //! --------------------
                //! one entity selected
                //! --------------------
                cout<<"DetailViewer::handleRemotePointLocationChanged()->____one single entity selected____"<<endl;

                const GeometryTag &loc = vecLoc[vecLoc.size()-1];
                int bodyIndex = loc.parentShapeNr;
                int subShapeNr = loc.subTopNr;

                cout<<"DetailViewer::handleRemotePointLocationChanged()->____"<<bodyIndex<<", "<<subShapeNr<<"____"<<endl;
                cout<<"DetailViewer::handleRemotePointLocationChanged()->____shape type: "<<(loc.subShapeType==TopAbs_FACE? "FACE____":"other than face____")<<endl;

                switch(loc.subShapeType)
                {
                case TopAbs_FACE:
                {
                    //! ----------------------------------------------
                    //! important node: do not use constant reference
                    //! const TopoDS_Face &theFace = ...
                    //! ----------------------------------------------
                    TopoDS_Face theFace = TopoDS::Face(gdb->MapOfBodyTopologyMap.value(bodyIndex).faceMap.FindKey(subShapeNr));

                    opencascade::handle<Geom_Surface> surface = BRep_Tool::Surface(theFace);

                    GeomAbs_SurfaceType surfaceType;
                    GeomToolsClass::getFaceType(theFace,surfaceType);

                    switch (surfaceType)
                    {
                    case GeomAbs_Plane: GeomToolsClass::getPlanarFaceInfo(theFace,theAxis,P); break;
                    case GeomAbs_Cylinder: GeomToolsClass::getCylindricalFaceInfo(theFace,theAxis,P); break;
                    case GeomAbs_Sphere: GeomToolsClass::getSphericalFaceInfo(theFace,theAxis,P); break;
                    case GeomAbs_Cone: GeomToolsClass::getConicalFaceInfo(theFace,theAxis,P); break;
                    default:
                        //! other cases ... to do ...
                        break;
                    }
                }
                    break;

                case TopAbs_EDGE:
                {
                    const TopoDS_Edge &theEdge = TopoDS::Edge(gdb->MapOfBodyTopologyMap.value(bodyIndex).edgeMap.FindKey(subShapeNr));
                    GeomAbs_CurveType curveType;
                    GeomToolsClass::getCurveType(theEdge,curveType);
                    switch(curveType)
                    {
                    case GeomAbs_Line: GeomToolsClass::getEdgeInfo(theEdge,theAxis,P); break;
                    case GeomAbs_Circle: GeomToolsClass::getCircleInfo(theEdge,theAxis,P); break;
                    case GeomAbs_Ellipse: GeomToolsClass::getEllipseInfo(theEdge,theAxis,P); break;
                    default:
                        //! other cases ... to do ...
                        break;
                    }
                }
                    break;

                case TopAbs_VERTEX:
                {
                    const TopoDS_Vertex &V = TopoDS::Vertex(gdb->MapOfBodyTopologyMap.value(bodyIndex).vertexMap.FindKey(subShapeNr));
                    P = BRep_Tool::Pnt(V);
                }
                    break;
                }
                cout<<"DetailViewer::handleRemotePointLocationChanged()->____Point("<<P.X()<<", "<<P.Y()<<", "<<P.Z()<<")____"<<endl;
            }
            else
            {
                //! ----------------------------------------------------------
                //! more than one entity selected
                //! if multiple entities are selected they have the same type
                //! ----------------------------------------------------------
                cout<<"DetailViewer::handleRemotePointLocationChanged()->____more than one entity selected____"<<endl;
                TopoDS_Compound aComp;
                TopoDS_Builder aBuilder;
                aBuilder.MakeCompound(aComp);

                bool selectionContainsVertices = false;
                bool selectionContainsFaces = false;
                bool selectionContainsShells = false;
                bool selectionContainsWires = false;
                bool selectionContainsEdges = false;

                for(std::vector<GeometryTag>::const_iterator anIt = vecLoc.cbegin(); anIt!=vecLoc.cend(); ++anIt)
                {
                    const GeometryTag &aLoc = *anIt;
                    TopAbs_ShapeEnum type = aLoc.subShapeType;
                    bool isParent = aLoc.isParent;
                    int parentShapeNr = aLoc.parentShapeNr;
                    int subShapeNr = aLoc.subTopNr;
                    if(isParent==true)
                    {
                        TopoDS_Shape aShape = sm->getDataBase()->bodyMap.value(parentShapeNr);
                        aBuilder.Add(aComp,aShape);
                    }
                    else
                    {
                        //! ----------------------------------------------------------------
                        //! if the selection is made of vertices, immediately exit from the
                        //! cycle, and compute the center of mass of the set of vertices
                        //! ----------------------------------------------------------------
                        if(type==TopAbs_VERTEX)
                        {
                            cout<<"DetailViewer::handleRemotePointLocationChanged()->____selection contain vertices____"<<endl;
                            selectionContainsVertices = true;
                            break;
                        }
                        switch(type)
                        {
                        case TopAbs_SHELL:
                        {
                            const TopoDS_Shape &aShape = sm->getDataBase()->MapOfBodyTopologyMap.value(parentShapeNr).shellMap.FindKey(subShapeNr);
                            aBuilder.Add(aComp,aShape);
                            selectionContainsShells = true;
                        }
                            break;

                        case TopAbs_FACE:
                        {
                            TopoDS_Shape aShape = sm->getDataBase()->MapOfBodyTopologyMap.value(parentShapeNr).faceMap.FindKey(subShapeNr);
                            aBuilder.Add(aComp,aShape);
                            selectionContainsFaces = true;
                        }
                            break;

                        case TopAbs_WIRE:
                        {
                            TopoDS_Shape aShape = sm->getDataBase()->MapOfBodyTopologyMap.value(parentShapeNr).wireMap.FindKey(subShapeNr);
                            aBuilder.Add(aComp,aShape);
                            selectionContainsWires = true;
                        }
                            break;

                        case TopAbs_EDGE:
                        {
                            TopoDS_Shape aShape = sm->getDataBase()->MapOfBodyTopologyMap.value(parentShapeNr).edgeMap.FindKey(subShapeNr);
                            aBuilder.Add(aComp,aShape);
                            selectionContainsEdges = true;
                        }
                            break;
                        }
                    }
                }
                if(selectionContainsEdges==true || selectionContainsWires==true)
                {
                    cout<<"DetailViewer::handleRemotePointLocationChanged()->____calculation on line bodies____"<<endl;
                    GProp_GProps props;
                    BRepGProp::LinearProperties(aComp,props);
                    P = props.CentreOfMass();
                    cout<<"DetailViewer::handleRemotePointLocationChanged()->____Point("<<P.X()<<", "<<P.Y()<<", "<<P.Z()<<")____"<<endl;
                }
                if(selectionContainsFaces==true || selectionContainsShells==true)
                {
                    cout<<"DetailViewer::handleRemotePointLocationChanged()->____calculation on surface bodies____"<<endl;
                    GProp_GProps props;
                    BRepGProp::SurfaceProperties(aComp,props);
                    P = props.CentreOfMass();
                    cout<<"DetailViewer::handleRemotePointLocationChanged()->____Point("<<P.X()<<", "<<P.Y()<<", "<<P.Z()<<")____"<<endl;
                }
                if(selectionContainsVertices==true)
                {
                    //! ------------------------------
                    //! iterate over the vertices and
                    //! compute the centroid
                    //! ------------------------------
                    double xB = 0, yB = 0, zB = 0;
                    int NbPoints = (int)vecLoc.size();
                    cout<<"DetailViewer::handleRemotePointLocationChanged()->____number of selected points: "<<NbPoints<<"____"<<endl;
                    for(std::vector<GeometryTag>::const_iterator anIt = vecLoc.cbegin(); anIt!=vecLoc.cend(); ++anIt)
                    {
                        const GeometryTag &aLoc = *anIt;
                        int parentShapeNr = aLoc.parentShapeNr;
                        int subShapeNr = aLoc.subTopNr;
                        TopoDS_Vertex aVertex = TopoDS::Vertex(sm->getDataBase()->MapOfBodyTopologyMap.value(parentShapeNr).vertexMap.FindKey(subShapeNr));
                        const gp_Pnt &aPoint = BRep_Tool::Pnt(aVertex);
                        xB = xB + aPoint.X();
                        yB = yB + aPoint.Y();
                        zB = zB + aPoint.Z();
                    }
                    xB = xB/NbPoints; yB = yB/NbPoints; zB = zB/NbPoints;
                    P.SetX(xB); P.SetY(yB); P.SetZ(zB);
                    cout<<"DetailViewer::handleRemotePointLocationChanged()->____Point("<<xB<<", "<<yB<<", "<<zB<<")____"<<endl;
                }
            }
        }

        //! --------------
        //! block signals
        //! --------------
        myCurNode->getModel()->blockSignals(true);

        //! ----------------------------------------------------
        //! update "X coordinate" "Y coordinate" "Z coordinate"
        //! ----------------------------------------------------
        QVariant data;
        data.setValue(P.X());
        Property prop_Xcoordinate("X coordinate",data,Property::PropertyGroup_Scope);
        data.setValue(P.Y());
        Property prop_Ycoordinate("Y coordinate",data,Property::PropertyGroup_Scope);
        data.setValue(P.Z());
        Property prop_Zcoordinate("Z coordinate",data,Property::PropertyGroup_Scope);

        myCurNode->replaceProperty("X coordinate",prop_Xcoordinate);
        myCurNode->replaceProperty("Y coordinate",prop_Ycoordinate);
        myCurNode->replaceProperty("Z coordinate",prop_Zcoordinate);

        //! --------------------------------------
        //! update the absolute coordinates using
        //! the current coordinate system
        //! --------------------------------------
        void *p = myCurNode->getPropertyValue<void*>("Coordinate system");
        QStandardItem *itemCS = static_cast<QStandardItem*>(p);
        SimulationNodeClass *nodeCS = itemCS->data(Qt::UserRole).value<SimulationNodeClass*>();
        SimulationNodeClass::nodeType type = nodeCS->getType();
        double xAbsCoord, yAbsCoord, zAbsCoord;
        switch(type)
        {
        case SimulationNodeClass::nodeType_coordinateSystem_global:
        {
            xAbsCoord = P.X();
            yAbsCoord = P.Y();
            zAbsCoord = P.Z();
        }
            break;

        case SimulationNodeClass::nodeType_coordinateSystem:
        {
            QVector<QVector<double>> directionalData = nodeCS->getPropertyValue<QVector<QVector<double>>>("Base directional data");
            xAbsCoord = directionalData.at(0).at(0)*P.X()+directionalData.at(1).at(0)*P.Y()+directionalData.at(2).at(0)*P.Z();
            yAbsCoord = directionalData.at(0).at(1)*P.X()+directionalData.at(1).at(1)*P.Y()+directionalData.at(2).at(1)*P.Z();
            xAbsCoord = directionalData.at(0).at(2)*P.X()+directionalData.at(1).at(2)*P.Y()+directionalData.at(2).at(2)*P.Z();
        }
            break;
        }

        data.setValue(xAbsCoord);
        Property prop_xAbsCoord("X abs coordinate",data,Property::PropertyGroup_Scope);
        myCurNode->replaceProperty("X abs coordinate",prop_xAbsCoord);
        data.setValue(yAbsCoord);
        Property prop_yAbsCoord("Y abs coordinate",data,Property::PropertyGroup_Scope);
        myCurNode->replaceProperty("Y abs coordinate",prop_yAbsCoord);
        data.setValue(zAbsCoord);
        Property prop_zAbsCoord("Z abs coordinate",data,Property::PropertyGroup_Scope);
        myCurNode->replaceProperty("Z abs coordinate",prop_zAbsCoord);

        //! ------------------
        //! update the marker
        //! ------------------
        SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
        cout<<"DetailViewer::handleRemotePointLocationChanged()->____updating the marker____"<<endl;
        emit requestHideAllMarkers(true);
        markerBuilder::addMarker(myCurNode,sm->getDataBase());

        //! ----------------
        //! unblock signals
        //! ----------------
        myCurNode->blockSignals(false);
    }
}

//! -------------------------------------------
//! function: DetailViewer::addTranformation()
//! details:
//! -------------------------------------------
bool DetailViewer::addPropertyToNode(const QString &name, Property::PropertyGroup thePropertyGroup)
{
    if(myCurNode!=Q_NULLPTR && myCurNode->getModel()->rowCount()!=0)
    {
        QVariant data;
        Property aProp(name,data,thePropertyGroup);
        myCurNode->addProperty(aProp);
        return true;
    }
    return false;
}

//! ---------------------------------------
//! function: deleteTransformationFromNode
//! details:
//! ---------------------------------------
void DetailViewer::deleteTransformationFromNode()
{
    if(myCurNode!=Q_NULLPTR) //! avoids crash when the a system of reference item does not exist
    {
        //!cout<<"DetailViewer::deleteTransformationFromNode()->____function called____"<<endl;
        QString name;
        //! remove by selecting the label or the value cell
        if(this->currentIndex().column()==0)name = this->currentIndex().data(Qt::DisplayRole).toString();
        else name = this->currentIndex().data(Qt::UserRole).value<Property>().getName();

        if(name == "Offset X" || name == "Offset Y" || name =="Offset Z" ||
           name == "Rotation X" || name =="Rotation Y" || name =="Rotation Z")
        {
            QModelIndex index = this->currentIndex();

            //! get the transformation to remove
            cout<<"removing: "<<this->currentIndex().data(Qt::UserRole).value<Property>().getName().toStdString()<<endl;

            //! remove the item from the tree
            myCurNode->getModel()->removeRow(index.row(),index.parent());

            //! insert code here for updating
            this->applyTransformations();
        }
    }
}

//! -------------------------------------------
//! function: getOrigin
//! details:  avoid writing four lines of code
//! -------------------------------------------
QVector<double> DetailViewer::getOrigin()
{
    cout<<"DetailViewer::getOrigin()->____function called____"<<endl;    
    double X_origin = myCurNode->getPropertyValue<double>("Origin X");
    double Y_origin = myCurNode->getPropertyValue<double>("Origin Y");
    double Z_origin = myCurNode->getPropertyValue<double>("Origin Z");
    QVector<double> Origin {X_origin, Y_origin, Z_origin};
    return Origin;
}

//! ---------------------------------------------
//! function: getDirectionalData
//! details:  avoid repeating four lines of code
//! ---------------------------------------------
QVector<QVector<double>> DetailViewer::getDirectionalData()
{
    QVector<double> dirDataX = myCurNode->getPropertyValue<QVector<double>>("X axis data");
    QVector<double> dirDataY = myCurNode->getPropertyValue<QVector<double>>("Y axis data");
    QVector<double> dirDataZ = myCurNode->getPropertyValue<QVector<double>>("Z axis data");
    QVector<QVector<double>> directionalData{dirDataX,dirDataY,dirDataZ};
    return directionalData;
}

void DetailViewer::updateCS()
{
    ;
}

//! -------------------------------------------------------------------
//! function: applyTransformation
//! details:  this apply a single transformation. The function is used
//!           by the "applyTransformations"
//! -------------------------------------------------------------------
void DetailViewer::applyTransformation(Property::typeOfTransformation theType, double delta)
{
    cout<<"DetailViewer::applyTransformation()->____function called____"<<endl;
    //! -----------------
    //! lock the signals
    //! -----------------
    myCurNode->getModel()->blockSignals(false);

    //! ----------------------
    //! the coordinate system
    //! ----------------------
    gp_Ax2 theCoordinateSystem;

    QVector<double> Origin = myCurNode->getPropertyValue<QVector<double>>("Base origin");
    double X = Origin.at(0); double Y = Origin.at(1); double Z = Origin.at(2);

    QVector<QVector<double>> baseDirData = myCurNode->getPropertyValue<QVector<QVector<double>>>("Base directional data");
    QVector<double> xAxisData = baseDirData.at(0);
    QVector<double> yAxisData = baseDirData.at(1);
    QVector<double> zAxisData = baseDirData.at(2);

    //! -----------------------------
    //! create the coordinate system
    //! -----------------------------
    theCoordinateSystem.SetLocation(gp_Pnt(X,Y,Z));
    theCoordinateSystem.SetDirection(gp_Dir(zAxisData.at(0),zAxisData.at(1),zAxisData.at(2)));
    theCoordinateSystem.SetXDirection(gp_Dir(xAxisData.at(0),xAxisData.at(1),xAxisData.at(2)));

    gp_Ax2 newCS;
    newCS = theCoordinateSystem;

    switch(theType)
    {
    case Property::typeOfTransformation_offsetX:
        newCS = theCoordinateSystem.Translated(gp_Vec(theCoordinateSystem.XDirection()).Multiplied(delta));
        break;
    case Property::typeOfTransformation_offsetY:
        newCS = theCoordinateSystem.Translated(gp_Vec(theCoordinateSystem.YDirection()).Multiplied(delta));
        break;
    case Property::typeOfTransformation_offsetZ:
        newCS = theCoordinateSystem.Translated(gp_Vec(theCoordinateSystem.Direction()).Multiplied(delta));
        break;
    case Property::typeOfTransformation_rotationX:
        newCS = theCoordinateSystem.Rotated(gp_Ax1(theCoordinateSystem.Location(),theCoordinateSystem.XDirection()),delta*M_PI/180);
        break;
    case Property::typeOfTransformation_rotationY:
        newCS = theCoordinateSystem.Rotated(gp_Ax1(theCoordinateSystem.Location(),theCoordinateSystem.YDirection()),delta*M_PI/180);
        break;
    case Property::typeOfTransformation_rotationZ:
        newCS = theCoordinateSystem.Rotated(gp_Ax1(theCoordinateSystem.Location(),theCoordinateSystem.Direction()),delta*M_PI/180);
        break;
    }

    //! -------------------------
    //! update the "Base origin"
    //! -------------------------
    QVariant data;
    QVector<double> updatedBaseOrigin;
    updatedBaseOrigin.push_back(newCS.Location().X());
    updatedBaseOrigin.push_back(newCS.Location().Y());
    updatedBaseOrigin.push_back(newCS.Location().Z());
    data.setValue(updatedBaseOrigin);
    Property prop_baseOrigin("Base origin",data,Property::PropertyGroup_Transformations);
    data.setValue(prop_baseOrigin);
    myCurNode->getPropertyItem("Base origin")->setData(data,Qt::UserRole);

    //! -----------------------------------
    //! update the "Base directional data"
    //! -----------------------------------
    xAxisData.replace(0,newCS.XDirection().X());
    xAxisData.replace(1,newCS.XDirection().Y());
    xAxisData.replace(2,newCS.XDirection().Z());
    yAxisData.replace(0,newCS.YDirection().X());
    yAxisData.replace(1,newCS.YDirection().Y());
    yAxisData.replace(2,newCS.YDirection().Z());
    zAxisData.replace(0,newCS.Direction().X());
    zAxisData.replace(1,newCS.Direction().Y());
    zAxisData.replace(2,newCS.Direction().Z());

    QVector<QVector<double>> updatedBaseDirectionalData;
    updatedBaseDirectionalData.push_back(xAxisData);
    updatedBaseDirectionalData.push_back(yAxisData);
    updatedBaseDirectionalData.push_back(zAxisData);
    data.setValue(updatedBaseDirectionalData);
    Property prop_baseDirectionalData("Base directional data",data,Property::PropertyGroup_Transformations);
    data.setValue(prop_baseDirectionalData);
    myCurNode->getPropertyItem("Base directional data")->setData(data,Qt::UserRole);

    //! -------------------
    //! unlock the signals
    //! -------------------
    myCurNode->getModel()->blockSignals(false);
}

//! -------------------------------
//! function: applyTransformations
//! details:
//! -------------------------------
void DetailViewer::applyTransformations()
{
    //! ------------------
    //! block the signals
    //! ------------------
    this->connectToSimulationManager(false);
    this->updateBaseData();
    int i;
    for(i=0; i<myCurNode->getModel()->rowCount(); i++)
    {
        if(myCurNode->getModel()->item(i,0)->data(Qt::DisplayRole).toString()=="Transformations") break;
    }

    int NumberOfTransformations = myCurNode->getModel()->item(i,0)->rowCount()-2;
    for(int k=0; k<NumberOfTransformations; k++)
    {
        QExtendedStandardItem *itemTrsf_label = static_cast<QExtendedStandardItem*>(myCurNode->getModel()->item(i,0)->child(k+2,0));
        QExtendedStandardItem *itemTrsf_value = static_cast<QExtendedStandardItem*>(myCurNode->getModel()->item(i,0)->child(k+2,1));
        QString name = itemTrsf_label->data(Qt::DisplayRole).toString();
        double value = itemTrsf_value->data(Qt::UserRole).value<Property>().getData().toDouble();
        if(name == "Offset X")
        {
            this->applyTransformation(Property::typeOfTransformation_offsetX,value);
        }
        if(name == "Offset Y")
        {
            this->applyTransformation(Property::typeOfTransformation_offsetY,value);
        }
        if(name == "Offset Z")
        {
            this->applyTransformation(Property::typeOfTransformation_offsetZ,value);
        }
        if(name == "Rotation X")
        {
            this->applyTransformation(Property::typeOfTransformation_rotationX,value);
        }
        if(name == "Rotation Y")
        {
            this->applyTransformation(Property::typeOfTransformation_rotationY,value);
        }
        if(name == "Rotation Z")
        {
            this->applyTransformation(Property::typeOfTransformation_rotationZ,value);
        }
    }

    //! -------------------
    //! unlock the signals
    //! -------------------
    this->connectToSimulationManager(true);

    //! ------------------------------------------------------------------------
    //! retrieve the current item (SimulationManager tree) and emit itemChanged
    //! ------------------------------------------------------------------------
    QModelIndex modelIndex = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"))->myTreeView->currentIndex();
    QStandardItem *item = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"))->getModel()->itemFromIndex(modelIndex);
    emit myCurNode->getModel()->itemChanged(item);

    //! ----------------------------------------
    //! request displayng the coordinate system
    //! ----------------------------------------
    QVector<double> O = myCurNode->getPropertyValue<QVector<double>>("Base origin");
    QVector<QVector<double>> DD = myCurNode->getPropertyValue<QVector<QVector<double>>>("Base directional data");
    this->requestDisplayTrihedron(O,DD);
}

//! -------------------------------------------------------------------------
//! function: updateBaseOriginCoords
//! details:  this function copies the "Origin X", "Origin Y", and "Origin Y"
//!           into the "base origin coordinates"
//! --------------------------------------------------------------------------
void DetailViewer::updateBaseOriginCoords()
{
    cout<<"DetailViewer::updateOriginCoords->____function called____"<<endl;

    //! -----------------
    //! lock the signals
    //! -----------------
    this->connectToSimulationManager(false);

    //! ------------------
    //! update the origin
    //! ------------------
    QVector<double> origin = this->getOrigin();
    double x = origin.at(0);
    double y = origin.at(1);
    double z = origin.at(2);

    QVector<double> baseOriginCoords {x,y,z};
    QVariant data;
    data.setValue(baseOriginCoords);
    Property prop_baseOrigin("Base origin",data,Property::PropertyGroup_Transformations);
    data.setValue(prop_baseOrigin);
    myCurNode->getPropertyItem("Base origin")->setData(data,Qt::UserRole);

    //! -------------------
    //! unlock the signals
    //! -------------------
    this->connectToSimulationManager(true);
}

//! -------------------------------------------------------------------------------
//! function: updateBaseDirectionalData
//! details:  this function copies the "X Axis data", "Y axis data", "Z axis data"
//!           into the "Base directional data"
//! -------------------------------------------------------------------------------
void DetailViewer::updateBaseDirectionalData()
{
    cout<<"DetailViewer::updateBaseDirectionalData()->____function called____"<<endl;

    //! -----------------
    //! lock the signals
    //! -----------------
    this->connectToSimulationManager(false);

    QVector<QVector<double>> dirData = this->getDirectionalData();
    QVariant data;
    data.setValue(dirData);
    Property prop_dirData("Base directional data",data,Property::PropertyGroup_Transformations);
    data.setValue(prop_dirData);
    myCurNode->getPropertyItem("Base directional data")->setData(data,Qt::UserRole);

    //! -------------------
    //! unlock the signals
    //! -------------------
    this->connectToSimulationManager(true);
}

//! ----------------------------------------------------------------------------------
//! function: updateBaseData()
//! details:  this function merges "updateBaseOrigin" and "updateBaseDirectionalData"
//! ----------------------------------------------------------------------------------
void DetailViewer::updateBaseData()
{
    cout<<"DetailViewer::updateBaseData()->____function called____"<<endl;
    this->updateBaseOriginCoords();
    this->updateBaseDirectionalData();
}

//! -----------------------------------------
//! function: getNode
//! details:  return the SimulationNodeClass
//! -----------------------------------------
SimulationNodeClass *DetailViewer::getNode()
{
    return myCurNode;
}

//! ------------------------------------------------------
//! function: handleDefineByChanged
//! details:  signal connected with the SimulationManager
//!           which handles the tabular data
//! ------------------------------------------------------
void DetailViewer::handleDefineByChanged()
{
    cout<<"DetailViewer::handleDefineByChanged()->____function called____"<<endl;

    //! -----------------------
    //! block the node signals
    //! -----------------------
    myCurNode->getModel()->blockSignals(true);
    //this->connectToSimulationManager(false);

    //! --------------------------------------------------------------
    //! additional switches are added/removed here, while the changes
    //! of the Tabular data must be done by the Simulation manager,
    //! which has direct access to the "Analysis settings" node
    //! --------------------------------------------------------------
    Property::loadDefinition old_componentX;
    Property::loadDefinition old_componentY;
    Property::loadDefinition old_componentZ;
    Property::loadDefinition old_magnitude;
    QVector<double> old_direction;

    Property::defineBy defineBy = myCurNode->getPropertyValue<Property::defineBy>("Define by");
    switch(defineBy)
    {
    case Property::defineBy_vector:
    {
        cout<<"DetailViewer::handleDefineByChanged()->____transition to vector: removing components and direction____"<<endl;

        //! ----------------------------------
        //! remove the components, if present
        //! ----------------------------------
        QExtendedStandardItem* itemXcomponent = myCurNode->getPropertyItem("X component");
        if(itemXcomponent!=Q_NULLPTR)
        {
            old_componentX = myCurNode->getOldXLoadDefinition();
            myCurNode->updateOldLoadDefinition(1,old_componentX);
            myCurNode->removeProperty("X component");
        }
        QExtendedStandardItem* itemYcomponent = myCurNode->getPropertyItem("Y component");
        if(itemYcomponent!=Q_NULLPTR)
        {
            old_componentY = myCurNode->getOldYLoadDefinition();
            myCurNode->updateOldLoadDefinition(2,old_componentY);
            myCurNode->removeProperty("Y component");
        }
        QExtendedStandardItem* itemZcomponent = myCurNode->getPropertyItem("Z component");
        if(itemZcomponent!=Q_NULLPTR)
        {
            old_componentZ = myCurNode->getOldZLoadDefinition();
            myCurNode->updateOldLoadDefinition(3,old_componentZ);
            myCurNode->removeProperty("Z component");
        }
        //! ----------------------------------
        //! remove the CS selector if present
        //! ----------------------------------
        QExtendedStandardItem *item_CS = myCurNode->getPropertyItem("Coordinate system");
        if(item_CS!=Q_NULLPTR)
        {
            myCurNode->removeProperty("Coordinate system");
        }
        //! ---------------------------------
        //! add the magnitude if not present
        //! ---------------------------------
        QExtendedStandardItem *item_magnitude = myCurNode->getPropertyItem("Magnitude");
        if(item_magnitude==Q_NULLPTR)
        {
            QVariant data;
            old_magnitude = myCurNode->getOldMagnitude();
            data.setValue(old_magnitude);
            Property property_magnitude("Magnitude",data,Property::PropertyGroup_Definition);
            myCurNode->addProperty(property_magnitude);
        }
        //! ---------------------------------
        //! add the direction if not present
        //! ---------------------------------
        QExtendedStandardItem *item_direction = myCurNode->getPropertyItem("Direction");
        if(item_direction==Q_NULLPTR)
        {
            QVariant data;
            old_direction = myCurNode->getOldDirection();
            data.setValue(old_direction);
            Property property_vectorDirection("Direction",data,Property::PropertyGroup_Definition);
            myCurNode->addProperty(property_vectorDirection);
        }
    }
        break;

    case Property::defineBy_components:
    {
        cout<<"DetailViewer::handleDefineByChanged()->____transition to components: removing magnitude and direction____"<<endl;

        //! --------------------------------
        //! remove the magnitude if present
        //! --------------------------------
        QExtendedStandardItem* itemMagnitude = myCurNode->getPropertyItem("Magnitude");
        if(itemMagnitude!=Q_NULLPTR)
        {
            //! store old value before removing
            old_magnitude = itemMagnitude->data(Qt::UserRole).value<Property>().getData().value<Property::loadDefinition>();
            myCurNode->updateOldMagnitude(old_magnitude);
            myCurNode->getModel()->removeRow(itemMagnitude->index().row(),itemMagnitude->parent()->index());
        }
        //! --------------------------------
        //! remove the direction if present
        //! --------------------------------
        QExtendedStandardItem *itemDirection = myCurNode->getPropertyItem("Direction");
        if(itemDirection!=Q_NULLPTR)
        {
            //! store the old values before removing
            old_direction = itemDirection->data(Qt::UserRole).value<Property>().getData().value<QVector<double>>();
            myCurNode->updateOldDirection(old_direction);
            myCurNode->getModel()->removeRow(itemDirection->index().row(),itemDirection->parent()->index());
        }
        //! -----------------------------------
        //! add the CS selector if not present
        //! -----------------------------------
        QExtendedStandardItem *item_CS = myCurNode->getPropertyItem("Coordinate system");
        if(item_CS==Q_NULLPTR)
        {
            SimulationManager *theSimulationManager = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));

            QExtendedStandardItem *itemCSroot = theSimulationManager->getTreeItem(SimulationNodeClass::nodeType_coordinateSystems);
            QExtendedStandardItem *itemGlobalCS = static_cast<QExtendedStandardItem*>(itemCSroot->child(0,0));
            void *itemGlobalCSvoid = (void*)itemGlobalCS;
            QVariant data;
            data.setValue(itemGlobalCSvoid);
            Property prop_coordinateSystem("Coordinate system",data,Property::PropertyGroup_Definition);
            myCurNode->addProperty(prop_coordinateSystem,1);
        }
        //! ----------------------------------
        //! add the components if not present
        //! ----------------------------------
        QExtendedStandardItem *item_componentX = myCurNode->getPropertyItem("X component");
        if(item_componentX==Q_NULLPTR)
        {
            QVariant data;
            old_componentX = myCurNode->getOldXLoadDefinition();
            data.setValue(old_componentX);
            Property property_Xcomponent("X component",data,Property::PropertyGroup_Definition);
            myCurNode->addProperty(property_Xcomponent);
        }
        QExtendedStandardItem *item_componentY = myCurNode->getPropertyItem("Y component");
        if(item_componentY==Q_NULLPTR)
        {
            QVariant data;
            old_componentY = myCurNode->getOldYLoadDefinition();
            data.setValue(old_componentY);
            Property property_Ycomponent("Y component",data,Property::PropertyGroup_Definition);
            myCurNode->addProperty(property_Ycomponent);
        }
        QExtendedStandardItem *item_componentZ = myCurNode->getPropertyItem("Z component");
        if(item_componentZ==Q_NULLPTR)
        {
            QVariant data;
            old_componentZ = myCurNode->getOldZLoadDefinition();
            data.setValue(old_componentZ);
            Property property_Zcomponent("Z component",data,Property::PropertyGroup_Definition);
            myCurNode->addProperty(property_Zcomponent);
        }
    }
        break;
    }

    //! ------------------------
    //! unlock the node signals
    //! ------------------------
    myCurNode->getModel()->blockSignals(false);
    //this->connectToSimulationManager(false);

    emit requestHandleTabularData();
}

//! --------------------------------------
//! function: handleLoadDefinitionChanged
//! details:
//! --------------------------------------
void DetailViewer::handleLoadDefinitionChanged()
{
    //cout<<"DetailViewer::handleLoadDefinitionChanged()->____function called____"<<endl;
    emit requestHandleLoadDefinitionChanged();
}

//! -------------------------------------------------
//! function: handleFilmCoefficientDefinitionChanged
//! details:  specific for convection
//! -------------------------------------------------
void DetailViewer::handleFilmCoefficientLoadDefinitionChanged(const QString& textData)
{
    emit requestHandleFilmCoefficientLoadDefinitionChanged(textData);
}

//! -------------------------------------------------
//! function: handleFilmCoefficientDefinitionChanged
//! details:  specific for convection
//! -------------------------------------------------
void DetailViewer::handleReferenceTemperatureLoadDefinitionChanged(const QString& textData)
{
    emit requestHandleReferenceTemperatureLoadDefinitionChanged(textData);
}

//! --------------------------------------------------------------------
//! function: handleMagnitudeLoadDefinitionChanged
//! details:  "Magnitude" is a property for magnitude of vectorial data
//! --------------------------------------------------------------------
void DetailViewer::handleMagnitudeLoadDefinitionChanged(const QString& textData)
{
    cout<<"DetailViewer::handleMagnitudeLoadDefinitionChanged()->____function called____"<<endl;
    emit requestHandleMagnitudeLoadDefinitionChanged(textData);
}

//! ---------------------------------------
//! function: handleXLoadDefinitionChanged
//! details:
//! ---------------------------------------
void DetailViewer::handleXLoadDefinitionChanged(const QString& textData)
{
    //cout<<"DetailViewer::handleXLoadDefinitionChanged()->____function called____"<<endl;
    //cout<<"DetailViewer::handleXLoadDefinitionChanged()->____editor text: "<<textData.toStdString()<<"____"<<endl;
    emit requestHandleXLoadDefinitionChanged(textData);
}

//! -----------------------------------------
//! function: handleYLoadDefinitionChanged()
//! details:
//! -----------------------------------------
void DetailViewer::handleYLoadDefinitionChanged(const QString& textData)
{
    //cout<<"DetailViewer::handleYLoadDefinitionChanged()->____function called____"<<endl;
    //cout<<"DetailViewer::handleYLoadDefinitionChanged()->____editor text: "<<textData.toStdString()<<"____"<<endl;
    emit requestHandleYLoadDefinitionChanged(textData);
}

//! -----------------------------------------
//! function: handleZLoadDefinitionChanged()
//! details:
//! -----------------------------------------
void DetailViewer::handleZLoadDefinitionChanged(const QString &textData)
{
    //cout<<"DetailViewer::handleZLoadDefinitionChanged()->____function called____"<<endl;
    //cout<<"DetailViewer::handleZLoadDefinitionChanged()->____editor text: "<<textData.toStdString()<<"____"<<endl;
    emit requestHandleZLoadDefinitionChanged(textData);
}

//! ------------------------------------------
//! function: handleAutoTimeSteppingChanged()
//! details:
//! ------------------------------------------
void DetailViewer::handleAutoTimeSteppingChanged()
{
    //cout<<"DetailViewer::handleAutoTimeSteppingChanged()->____function called____"<<endl;
    Property::autoTimeStepping theAutoTimeStepping =  myCurNode->getPropertyValue<Property::autoTimeStepping>("Auto time stepping");

    QExtendedStandardItem *item;
    int currentStepNumber = myCurNode->getPropertyValue<int>("Current step number");

    if(theAutoTimeStepping == Property::autoTimeStepping_ProgramControlled)
    {
        //! set the default values
        QVector<int> N;
        N.append(5); N.append(2); N.append(10); N.append(0);

        item = myCurNode->getPropertyItem("Number of substeps");
        if(item!=Q_NULLPTR)
        {
            myCurNode->getModel()->removeRow(item->index().row(),item->index().parent());
        }
        item = myCurNode->getPropertyItem("Initial substeps");
        if(item!=Q_NULLPTR)
        {
            myCurNode->getModel()->removeRow(item->index().row(),item->index().parent());
            item = myCurNode->getPropertyItem("Minimum substeps");
            myCurNode->getModel()->removeRow(item->index().row(),item->index().parent());
            item = myCurNode->getPropertyItem("Maximum substeps");
            myCurNode->getModel()->removeRow(item->index().row(),item->index().parent());
        }
        QVariant data;
        data.setValue(N);
        myCurNode->getTabularDataModel()->setDataRC(data,currentStepNumber,3);
    }
    else if(theAutoTimeStepping == Property::autoTimeStepping_ON)
    {

        //! remove the previous settings
        item = myCurNode->getPropertyItem("Number of substeps");
        if(item!=Q_NULLPTR)
        {
            myCurNode->getModel()->removeRow(item->index().row(),item->index().parent());
        }

        //! read from the tabData
        //int currentStepNumber = myCurNode->getPropertyItem("Current step number")->data(Qt::UserRole).value<Property>().getData().toInt();
        CustomTableModel *tabData = myCurNode->getTabularDataModel();
        QVector<int> N = tabData->dataRC(currentStepNumber,3,Qt::EditRole).value<QVector<int>>();

        //! -----
        QVariant data;
        data.setValue(N[0]);
        Property prop_NdivIni("Initial substeps",data,Property::PropertyGroup_StepControls);
        myCurNode->addProperty(prop_NdivIni);
        data.setValue(N[1]);
        Property prop_NdivMin("Minimum substeps",data,Property::PropertyGroup_StepControls);
        myCurNode->addProperty(prop_NdivMin);
        data.setValue(N[2]);
        Property prop_NdivMax("Maximum substeps",data,Property::PropertyGroup_StepControls);
        myCurNode->addProperty(prop_NdivMax);

        //! set the time step division policy
        N.remove(3);
        N.push_back(1); //! "1" is for "ON"
        data.setValue(N);
        tabData->setDataRC(data,currentStepNumber,3);
    }
    else if(theAutoTimeStepping == Property::autoTimeStepping_OFF)
    {
        item = myCurNode->getPropertyItem("Initial substeps");
        if(item!=Q_NULLPTR)
        {
            myCurNode->getModel()->removeRow(item->index().row(),item->index().parent());
            item = myCurNode->getPropertyItem("Minimum substeps");
            myCurNode->getModel()->removeRow(item->index().row(),item->index().parent());
            item = myCurNode->getPropertyItem("Maximum substeps");
            myCurNode->getModel()->removeRow(item->index().row(),item->index().parent());
        }

        //! read from the tabData
        //int currentStepNumber = myCurNode->getPropertyItem("Current step number")->data(Qt::UserRole).value<Property>().getData().toInt();
        CustomTableModel *tabData = myCurNode->getTabularDataModel();
        QVector<int> N = tabData->dataRC(currentStepNumber,3,Qt::EditRole).value<QVector<int>>();

        QVariant data;
        data.setValue(N[0]);
        Property prop_Ndiv("Number of substeps",data,Property::PropertyGroup_StepControls);
        myCurNode->addProperty(prop_Ndiv);

        //! set the time step division policy
        N.remove(3);
        N.push_back(2); //! "1" is for "ON"
        data.setValue(N);
        tabData->setDataRC(data,currentStepNumber,3);
    }
}

//! ----------------------------------
//! function: handleTimeStepDivisione
//! details:
//! ----------------------------------
void DetailViewer::handleTimeDivisionChanged()
{
    int currentStepNumber = myCurNode->getPropertyItem("Current step number")->data(Qt::UserRole).value<Property>().getData().toInt();
    Property::autoTimeStepping ats = myCurNode->getPropertyItem("Auto time stepping")->data(Qt::UserRole).value<Property>().getData().value<Property::autoTimeStepping>();

    int n0,n1,n2,n3;

    switch(ats)
    {
    case Property::autoTimeStepping_ProgramControlled:
    {
        n0 = 5;
        n1 = 2;
        n2 = 10;
        n3 = 0;
    }
        break;
    case Property::autoTimeStepping_ON:
    {
        n0 = myCurNode->getPropertyItem("Initial substeps")->data(Qt::UserRole).value<Property>().getData().toInt();
        n1 = myCurNode->getPropertyItem("Minimum substeps")->data(Qt::UserRole).value<Property>().getData().toInt();
        n2 = myCurNode->getPropertyItem("Maximum substeps")->data(Qt::UserRole).value<Property>().getData().toInt();
        n3 = 1;
    }
        break;
    case Property::autoTimeStepping_OFF:
        n0 = myCurNode->getPropertyItem("Number of substeps")->data(Qt::UserRole).value<Property>().getData().toInt();
        n1 = n2 = 0;
        n3 = 2;
        break;
    }
    QVector<int> N; N.append(n0); N.append(n1); N.append(n2); N.append(n3);
    QVariant data; data.setValue(N);
    myCurNode->getTabularDataModel()->setDataRC(data,currentStepNumber,3);
}

//! ------------------------------------------
//! function: handleChangeMeshNodesVisibility
//! details:
//! ------------------------------------------
void DetailViewer::handleChangeMeshNodesVisibility(bool meshNodesVisible)
{
    cout<<"DetailViewer::handleChangeMeshNodesVisibility->____function called: value: "<<meshNodesVisible<<endl;
    emit requestChangeMeshNodesVisibility(meshNodesVisible);
}

//! ------------------------------------
//! function: handleMeshSmoothingChange
//! details:
//! ------------------------------------
void DetailViewer::handleMeshSmoothingChange()
{
    emit requestInvalidateAllMeshes();
}

//! -----------------------------------------
//! function: handleRequestStartInterpolator
//! details:
//! -----------------------------------------
void DetailViewer::handleRequestStartInterpolator()
{
    cout<<"DetailViewer::handleRequestStartInterpolator()->____function called____"<<endl;
    emit startInterpolator();}

//! ----------------------------------------
//! function: handleGlobalMeshControlChange
//! details:
//! ----------------------------------------
void DetailViewer::handleGlobalMeshControlChange()
{
    cout<<"DetailViewer::handleglobalMeshControlChange->____function called____"<<endl;
    emit requestGlobalMeshControlChange();
}

//! -----------------------------------------------------------
//! function: handleRequestStartOpenFoamScalarDataTranslator()
//! -----------------------------------------------------------
void DetailViewer::handleRequestStartOpenFoamScalarDataTranslator()
{
    cout<<"DetailViewer::handleRequestStartOpenFoamScalarDataTranslator->____function called____"<<endl;
    emit startOpenFoamScalarDataTranslator();
}

//! ---------------------------------------------------
//! function: handleRemapFlagChanged()
//! details:  add/remove the control "Remapping steps"
//! ---------------------------------------------------
void DetailViewer::handleRemapFlagChanged()
{
    static int firstCall;
    firstCall++;
    static int remappingStepsOld;
    bool remapFlag = myCurNode->getPropertyItem("Remap")->data(Qt::UserRole).value<Property>().getData().toBool();
    if(remapFlag)
    {
        QVariant data;
        if(firstCall==1) data.setValue(1);
        else data.setValue(remappingStepsOld);
        Property prop_remappingSteps("Remapping steps",data,Property::PropertyGroup_Advanced);
        myCurNode->addProperty(prop_remappingSteps);
    }
    else
    {
        remappingStepsOld = myCurNode->getPropertyItem("Remapping steps")->data(Qt::UserRole).value<Property>().getData().toInt();
        myCurNode->removeProperty("Remapping steps");
    }
}

//! -----------------------------------------
//! function: handleStepSelectionModeChanged
//! details:
//! -----------------------------------------
void DetailViewer::handleStepSelectionModeChange()
{
    cout<<"DetailViewer::handleStepSelectionModeChange->____function called____"<<endl;
    int val = myCurNode->getPropertyItem("Step selection mode")->data(Qt::UserRole).value<Property>().getData().toInt();
    QVariant data;

    switch(val)
    {
    //! 0 => all steps 1 => first 2 => last 3 => by number
    case 0: case 1: case 2:
    {
        if(myCurNode->getPropertyItem("Source directory")!=Q_NULLPTR) myCurNode->removeProperty("Source directory");
        if(myCurNode->getPropertyItem("Step number")!=Q_NULLPTR) myCurNode->removeProperty("Step number");

        if(myCurNode->getPropertyItem("Source file")==Q_NULLPTR)
        {
            data.setValue(QString());
            Property prop_sourceFile("Source file",data,Property::PropertyGroup_Definition);
            myCurNode->addProperty(prop_sourceFile,2);
        }
    }
        break;

    //! step selection mode by number
    case 3:
    {
        if(myCurNode->getPropertyItem("Source directory")!=Q_NULLPTR) myCurNode->removeProperty("Source directory");
        if(myCurNode->getPropertyItem("Step number")==Q_NULLPTR)
        {
            data.setValue(1);
            Property prop_stepNumber("Step number",data,Property::PropertyGroup_Definition);
            myCurNode->addProperty(prop_stepNumber,2);
        }

        if(myCurNode->getPropertyItem("Source file")==Q_NULLPTR)
        {
            data.setValue(QString());
            Property prop_sourceFile("Source file",data,Property::PropertyGroup_Definition);
            myCurNode->addProperty(prop_sourceFile,2);
        }
    }
        break;

    //! automatic time history importer
    case 4:
    {
        if(myCurNode->getPropertyItem("Source file")!=Q_NULLPTR) myCurNode->removeProperty("Source file");
        if(myCurNode->getPropertyItem("Step number")==Q_NULLPTR)
        {
            data.setValue(1);
            Property prop_stepNumber("Step number",data,Property::PropertyGroup_Definition);
            myCurNode->addProperty(prop_stepNumber);
        }
        if(myCurNode->getPropertyItem("Source directory")==Q_NULLPTR)
        {
            QDir curDir;
            data.setValue(QString(curDir.currentPath()));
            Property prop_sourceDir("Source directory",data,Property::PropertyGroup_Definition);
            myCurNode->addProperty(prop_sourceDir,2);
        }
    }
        break;
    }
}

//! ------------------------------------
//! function: postBackgroundChangeEvent
//! details:
//! ------------------------------------
void DetailViewer::postBackgroundChangeEvent()
{
    cout<<"DetailViewer::postBackgroundChangeEvent()->____function called____"<<endl;
    QWidget *widget = tools::getWidgetByName("maingwindow");
    if(widget!=Q_NULLPTR)
    {
        int gradient = myCurNode->getPropertyItem("Gradient")->data(Qt::UserRole).value<Property>().getData().toInt();
        QVector<int> rgb1 = myCurNode->getPropertyItem("First color")->data(Qt::UserRole).value<Property>().getData().value<QVector<int>>();
        QVector<int> rgb2 = myCurNode->getPropertyItem("Second color")->data(Qt::UserRole).value<Property>().getData().value<QVector<int>>();

        QColor firstColor;
        firstColor.setRed(rgb1.at(0));
        firstColor.setGreen(rgb1.at(1));
        firstColor.setBlue(rgb1.at(2));

        QColor secondColor;
        secondColor.setRed(rgb2.at(0));
        secondColor.setGreen(rgb2.at(1));
        secondColor.setBlue(rgb2.at(2));

        cout<<"DetailViewer::postBackgroundChangeEvent()->____Gradient: "<<gradient<<" First color: "<<
              firstColor.name().toStdString()<<" Second color: "<<secondColor.name().toStdString()<<"____"<<endl;

        QBackgroundEvent *theBkgEvent = new QBackgroundEvent(gradient, firstColor, secondColor);

        cout<<"DetailViewer::postBackgroundChangeEvent()->____posting event: "<<theBkgEvent->type()<<"____"<<endl;

        QApplication::postEvent(widget,theBkgEvent);
    }
    else
    {
        cout<<"DetailViewer::postBackgroundChangeEvent()->____cannot post bkg event____"<<endl;
    }
}

//! ---------------------------------
//! function: handleNbThreadsChanged
//! details:
//! ---------------------------------
void DetailViewer::handleNbThreadsChanged()
{
    cout<<"DetailViewer::handleNbThreadsChanged()->____function called____"<<endl;
    int NbThreads = myCurNode->getPropertyItem("Number of threads")->data(Qt::UserRole).value<Property>().getData().toInt();
    cout<<"DetailViewer::handleNbThreadsChanged()->____"<<NbThreads<<"____"<<endl;
}

//! ---------------------------------------------------
//! function: handleByChanged
//! details:  add/remove controls to the detail viewer
//! ---------------------------------------------------
void DetailViewer::handleByChanged()
{
    cout<<"DetailViewer::handleByChanged()->____function called____"<<endl;
    myCurNode->getModel()->blockSignals(true);
    int by = myCurNode->getPropertyValue<int>("By");
    QVariant data;
    switch (by)
    {
    case 0:
    {
        myCurNode->removeProperty("Set number");
        data.setValue(double(0.0));
        Property prop_displayTime("Display time",data,Property::PropertyGroup_Definition);
        myCurNode->addProperty(prop_displayTime,2);
    }
        break;

    case 1:
    {
        myCurNode->removeProperty("Display time");
        data.setValue(int(1));
        Property prop_setNumber("Set number",data,Property::PropertyGroup_Definition);
        myCurNode->addProperty(prop_setNumber,2);
    }
        break;
    }
    myCurNode->getModel()->blockSignals(false);
}

//! ------------------------------------------------
//! function: handleByChanged
//! details:  a solution component has been changed
//!           the previous result must be deleted
//! ------------------------------------------------
void DetailViewer::handleSolutionComponentChanged()
{
    cout<<"DetailViewer::handleSolutionComponentChanged()->____function called____"<<endl;
    //! ...
}

//! -----------------------------------------
//! function: expandChildren
//! details:  expand children under a parent
//! -----------------------------------------
void DetailViewer::expandChildren(QModelIndex parent)
{
    for(int r=0; r<this->model()->rowCount(parent); ++r)
    {
        QModelIndex index = this->model()->index(r, 0, parent);
        // here the your applicable code
        if(this->model()->hasChildren(index) )
        {
            expandChildren(index);
            if(!this->isExpanded(index))
                this->expand(index);
        }
    }
}

//! ------------------------------------
//! function: handleTypeOfSizingChanged
//! details:
//! ------------------------------------
void DetailViewer::handleTypeOfSizingChanged()
{
    cout<<"DetailViewer::handleTypeOfSizingChanged()->____function called____"<<endl;
    QExtendedStandardItem *itemType = myCurNode->getPropertyItem("Sizing type");
    int typeOfSizing = itemType->data(Qt::UserRole).value<Property>().getData().toInt();

    QVariant data;
    if(myCurNode->getPropertyItem("Element size")==Q_NULLPTR)
    {
        cout<<"____Adding element size____"<<endl;
        data.setValue(double(10.0));
        Property prop_elementSize("Element size",data,Property::PropertyGroup_Definition);
        myCurNode->addProperty(prop_elementSize);
    }
    if(myCurNode->getPropertyItem("Number of divisions")==Q_NULLPTR)
    {
        cout<<"____Adding number of divisions____"<<endl;
        data.setValue(int(5));
        Property prop_NbDivisions("Number of divisions",data,Property::PropertyGroup_Definition);
        myCurNode->addProperty(prop_NbDivisions);
    }

    QExtendedStandardItem *itemElementSize = myCurNode->getPropertyItem("Element size");
    QExtendedStandardItem *itemNbDivisions = myCurNode->getPropertyItem("Number of divisions");
    QModelIndex indexElementSize = itemElementSize->index();
    QModelIndex indexNbDivisions = itemNbDivisions->index();

    switch(typeOfSizing)
    {
    case 0:
    {
        if(this->isRowHidden(indexElementSize.row(),indexElementSize.parent())) this->setRowHidden(indexElementSize.row(),indexElementSize.parent(),false);
        this->setRowHidden(indexNbDivisions.row(),indexNbDivisions.parent(),true);     
    }
        break;

    case 1:
    {
        if(this->isRowHidden(indexNbDivisions.row(),indexNbDivisions.parent())) this->setRowHidden(indexNbDivisions.row(),indexNbDivisions.parent(),false);
        this->setRowHidden(indexElementSize.row(),indexElementSize.parent(),true);
    }
        break;
    }
}

//! -------------------------------------------------------------
//! function: handleBoltStatusDefinedByChanged
//! details:  handle the visibility of the two controls "Load"
//!           and "Adjustment"
//! -------------------------------------------------------------
void DetailViewer::handleBoltStatusDefinedByChanged()
{
    cout<<"DetailViewer::handleBoltStatusDefinedByChanged()->____function called____"<<endl;
    QExtendedStandardItem *itemDefineBy = myCurNode->getPropertyItem("Bolt status");
    Property::boltStatusDefinedBy boltDefineBy = itemDefineBy->data(Qt::UserRole).value<Property>().getData().value<Property::boltStatusDefinedBy>();

    //QExtendedStandardItem *itemLoad = myCurNode->getPropertyItem("Load");
    //QExtendedStandardItem *itemAdjustment = myCurNode->getPropertyItem("Adjustment");

    SimulationNodeClass* nodeAnalysisSettings = myCurModelIndex.parent().child(0,0).data(Qt::UserRole).value<SimulationNodeClass*>();
    int currentStepNumber = nodeAnalysisSettings->getPropertyItem("Current step number")->data(Qt::UserRole).value<Property>().getData().toInt();
    int row = currentStepNumber;
    //int col = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"))->calculateStartColumn();
    SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
    int col = mainTreeTools::calculateStartColumn(sm->myTreeView);
    CustomTableModel *tabularDataModel = nodeAnalysisSettings->getTabularDataModel();
    QVariant data;

    switch(boltDefineBy)
    {
    case Property::boltStatusDefinedBy_load:
    {
        //! -----------------------------------------------
        //! show "Load" control, hide "Adjustment" control
        //! -----------------------------------------------
        this->setPropertyVisible("Load",true);
        this->setPropertyVisible("Adjustment",false);

        //! -------------------------------------------------------
        //! enable "Load" and disable "Adjustment" in tabular data
        //! -------------------------------------------------------
        double boltLoad = 0.0;
        data.setValue(boltLoad);
        QModelIndex indexLoad = tabularDataModel->makeIndex(row,col+1);
        tabularDataModel->setData(indexLoad,data,Qt::EditRole);

        data.setValue(QString("N/A"));
        QModelIndex indexAdjustment = tabularDataModel->makeIndex(row,col+2);
        tabularDataModel->setData(indexAdjustment,data,Qt::EditRole);
    }
        break;

    case Property::boltStatusDefinedBy_adjustment:
    {
        //! -----------------------------------------------
        //! show "Adjustment" control, hide "Load" control
        //! -----------------------------------------------
        this->setPropertyVisible("Adjustment",true);
        this->setPropertyVisible("Load",false);

        //! -------------------------------------------------------
        //! enable "Load" and disable "Adjustment" in tabular data
        //! -------------------------------------------------------
        double boltAdjustment =0.0;
        data.setValue(boltAdjustment);
        QModelIndex indexAdjustment = tabularDataModel->makeIndex(row,col+2);
        tabularDataModel->setData(indexAdjustment,data,Qt::EditRole);

        data.setValue(QString("N/A"));
        QModelIndex indexLoad = tabularDataModel->makeIndex(row,col+1);
        tabularDataModel->setData(indexLoad,data,Qt::EditRole);
    }
        break;

    case Property::boltStatusDefinedBy_open:
    case Property::boltStatusDefinedBy_lock:
    {
        //! --------------------------------------
        //! hide "Load" and "Adjustment" controls
        //! --------------------------------------
        this->setPropertyVisible("Adjustment",false);
        this->setPropertyVisible("Load",false);

        //! ------------------------------------------------
        //! disable "Load" and "Adjustment" in tabular data
        //! ------------------------------------------------
        data.setValue(QString("N/A"));
        QModelIndex indexLoad = tabularDataModel->makeIndex(row,col+1);
        tabularDataModel->setData(indexLoad,data,Qt::EditRole);

        QModelIndex indexAdjustment = tabularDataModel->makeIndex(row,col+2);
        tabularDataModel->setData(indexAdjustment,data,Qt::EditRole);
    }
        break;
    }
    data.setValue(boltDefineBy);
    QModelIndex indexBoltStatusDefinedBy = tabularDataModel->makeIndex(row,col);
    tabularDataModel->setData(indexBoltStatusDefinedBy,data,Qt::EditRole);
}

//! -------------------------------------------------------------
//! function: handleBoltLoadChanged
//! details:  transfer the value of "Load" into the tabular data
//! -------------------------------------------------------------
void DetailViewer::handleBoltLoadChanged()
{
    cout<<"DetailViewer::handleBoltLoadChanged()->____function called: transfer the \"Load\" value into the tabular data viewer____"<<endl;
    double load = myCurNode->getPropertyValue<double>("Load");
    //cout<<"____current value of the bolt load: "<<load<<"____"<<endl;

    //! --------------------------------------
    //! retrieve the "Analysis settings" item
    //! --------------------------------------
    SimulationNodeClass *nodeAnalysisSettings = myCurModelIndex.parent().child(0,0).data(Qt::UserRole).value<SimulationNodeClass*>();

    int currentTimeStep = nodeAnalysisSettings->getPropertyValue<int>("Current step number");
    CustomTableModel *tabularDataModel = nodeAnalysisSettings->getTabularDataModel();

    int row = currentTimeStep;
    SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
    int col = mainTreeTools::calculateStartColumn(sm->myTreeView)+1;

    QModelIndex indexLoad = tabularDataModel->makeIndex(row,col);
    QVariant data;
    data.setValue(load);
    tabularDataModel->setData(indexLoad,data,Qt::EditRole);
}

//! ----------------------------------------------------------------
//! function: handleBoltAdjustmentChanged
//! details:  transfer the "Adjustment" value into the tabular data
//! ----------------------------------------------------------------
void DetailViewer::handleBoltAdjustmentChanged()
{
    cout<<"DetailViewer::handleBoltAdjustmentChanged()->____function called: transfer the \"Adjustment\" value into the tabular data viewer____"<<endl;
    double adjustment = myCurNode->getPropertyValue<double>("Adjustment");

    //! --------------------------------------
    //! retrieve the "Analysis settings" item
    //! --------------------------------------
    SimulationNodeClass *nodeAnalysisSettings = myCurModelIndex.parent().child(0,0).data(Qt::UserRole).value<SimulationNodeClass*>();

    int currentTimeStep = nodeAnalysisSettings->getPropertyItem("Current step number")->data(Qt::UserRole).value<Property>().getData().toInt();
    CustomTableModel *tabularDataModel = nodeAnalysisSettings->getTabularDataModel();
    int row = currentTimeStep;
    //int col = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"))->calculateStartColumn()+2;
    SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
    int col = mainTreeTools::calculateStartColumn(sm->myTreeView)+2;

    QModelIndex indexAdjustment = tabularDataModel->makeIndex(row,col);
    QVariant data;
    data.setValue(adjustment);
    tabularDataModel->setData(indexAdjustment,data,Qt::EditRole);
}

//! ------------------------------------------------------------------
//! function: handleDOFselectorChange
//! details:  0 => Program controlled => remove all the DOFs controls
//!           1 => Manual => add all the controls
//! ------------------------------------------------------------------
void DetailViewer::handleDOFselectorChange()    // now disconnected
{
    cout<<"DetailViewer::handleDOFselectorChange()->____function called____"<<endl;
    int typeOfCoupling = myCurNode->getPropertyItem("Coupling")->data(Qt::UserRole).value<Property>().getData().toInt();
    int DOFsSelection = myCurNode->getPropertyItem("DOFs selection")->data(Qt::UserRole).value<Property>().getData().toInt();

    QStandardItem *itemXcomponent = myCurNode->getPropertyItem("X component ");
    QStandardItem *itemYcomponent = myCurNode->getPropertyItem("Y component ");
    QStandardItem *itemZcomponent = myCurNode->getPropertyItem("Z component ");
    QStandardItem *itemXrotation = myCurNode->getPropertyItem("X rotation");
    QStandardItem *itemYrotation = myCurNode->getPropertyItem("Y rotation");
    QStandardItem *itemZrotation = myCurNode->getPropertyItem("Z rotation");

    Property::PropertyGroup thePropertyGroup = Property::PropertyGroup_Advanced;

    QVariant data;
    switch(DOFsSelection)
    {
    case 0:
    {
        //! -----------------------------------------------
        //! "Kinematic"/"Distributed" "Program controlled"
        //! -----------------------------------------------
        myCurNode->removeProperty("X component ");
        myCurNode->removeProperty("Y component ");
        myCurNode->removeProperty("Z component ");
        myCurNode->removeProperty("X rotation");
        myCurNode->removeProperty("Y rotation");
        myCurNode->removeProperty("Z rotation");
    }
        break;

    case 1:
    {
        switch(typeOfCoupling)
        {
        //! ---------------------
        //! "Kinematic" "Manual"
        //! ---------------------
        case 0:
        {
            if(itemXcomponent==Q_NULLPTR)
            {
                data.setValue(1);
                Property prop_Xcomponent("X component ",data,thePropertyGroup);
                myCurNode->addProperty(prop_Xcomponent);
            }
            if(itemYcomponent==Q_NULLPTR)
            {
                data.setValue(1);
                Property prop_Ycomponent("Y component ",data,thePropertyGroup);
                myCurNode->addProperty(prop_Ycomponent);
            }
            if(itemZcomponent==Q_NULLPTR)
            {
                data.setValue(1);
                Property prop_Zcomponent("Z component ",data,thePropertyGroup);
                myCurNode->addProperty(prop_Zcomponent);
            }
            myCurNode->removeProperty("X rotation");
            myCurNode->removeProperty("Y rotation");
            myCurNode->removeProperty("Z rotation");
        }
            break;

        case 1:
        {
            //! -----------------------
            //! "Distributed" "Manual"
            //! -----------------------
            if(itemXcomponent==Q_NULLPTR)
            {
                data.setValue(1);
                Property prop_Xcomponent("X component ",data,thePropertyGroup);
                myCurNode->addProperty(prop_Xcomponent);
            }
            if(itemYcomponent==Q_NULLPTR)
            {
                data.setValue(1);
                Property prop_Ycomponent("Y component ",data,thePropertyGroup);
                myCurNode->addProperty(prop_Ycomponent);
            }
            if(itemZcomponent==Q_NULLPTR)
            {
                data.setValue(1);
                Property prop_Zcomponent("Z component ",data,thePropertyGroup);
                myCurNode->addProperty(prop_Zcomponent);
            }
            if(itemXrotation==Q_NULLPTR)
            {
                data.setValue(1);
                Property prop_Xrot("X rotation",data,thePropertyGroup);
                myCurNode->addProperty(prop_Xrot);
            }
            if(itemYrotation==Q_NULLPTR)
            {
                data.setValue(1);
                Property prop_Yrot("Y rotation",data,thePropertyGroup);
                myCurNode->addProperty(prop_Yrot);
            }
            if(itemZrotation==Q_NULLPTR)
            {
                data.setValue(1);
                Property prop_Zrot("Z rotation",data,thePropertyGroup);
                myCurNode->addProperty(prop_Zrot);
            }
        }
            break;
        }
    }
        break;
    }
}

//! -------------------------------------------------------------
//! function: handleCouplingChanged()
//! details:  0 => Kinematic 1 => Distributed
//! -------------------------------------------------------------
void DetailViewer::handleCouplingChanged()
{
    int typeOfCoupling = myCurNode->getPropertyValue<int>("Coupling");
    int DOFsSelection = myCurNode->getPropertyValue<int>("DOFs selection");

    QStandardItem *itemXcomponent = myCurNode->getPropertyItem("X component ");
    QStandardItem *itemYcomponent = myCurNode->getPropertyItem("Y component ");
    QStandardItem *itemZcomponent = myCurNode->getPropertyItem("Z component ");

    QStandardItem *itemXrotation = myCurNode->getPropertyItem("X rotation");
    QStandardItem *itemYrotation = myCurNode->getPropertyItem("Y rotation");
    QStandardItem *itemZrotation = myCurNode->getPropertyItem("Z rotation");

    QVariant data;

    switch(typeOfCoupling)
    {
    case 0:
    {
        //! ---------------------------------------------------
        //! "Kinematic" coupling does not have rotational DOFs
        //! ---------------------------------------------------
        myCurNode->removeProperty("X rotation");
        myCurNode->removeProperty("Y rotation");
        myCurNode->removeProperty("Z rotation");
        switch(DOFsSelection)
        {
        case 0:
        {
            //! -------------------------------------
            //! "Kinematic" and "Program controlled"
            //! -------------------------------------
            myCurNode->removeProperty("X component ");
            myCurNode->removeProperty("Y component ");
            myCurNode->removeProperty("Z component ");
        }
            break;

        case 1:
        {
            //! -------------------------
            //! "Kinematic" and "Manual"
            //! -------------------------
            if(itemXcomponent==Q_NULLPTR)
            {
                data.setValue(1);
                Property prop_Xcomponent("X component ",data,Property::PropertyGroup_Advanced);
                myCurNode->addProperty(prop_Xcomponent);
            }
            if(itemYcomponent==Q_NULLPTR)
            {
                data.setValue(1);
                Property prop_Ycomponent("Y component ",data,Property::PropertyGroup_Advanced);
                myCurNode->addProperty(prop_Ycomponent);
            }
            if(itemZcomponent==Q_NULLPTR)
            {
                data.setValue(1);
                Property prop_Zcomponent("Z component ",data,Property::PropertyGroup_Advanced);
                myCurNode->addProperty(prop_Zcomponent);
            }
            myCurNode->removeProperty("X rotation");
            myCurNode->removeProperty("Y rotation");
            myCurNode->removeProperty("Z rotation");
        }
            break;
        }
    }
        break;

    case 1:
    {
        switch(DOFsSelection)
        {
        case 0:
        {
            //! ---------------------------------------
            //! "Distributed" and "Program controlled"
            //! ---------------------------------------
            myCurNode->removeProperty("X component ");
            myCurNode->removeProperty("Y component ");
            myCurNode->removeProperty("Z component ");
            myCurNode->removeProperty("X rotation");
            myCurNode->removeProperty("Y rotation");
            myCurNode->removeProperty("Z rotation");
        }
            break;

        case 1:
        {
            //! ---------------------------
            //! "Distributed" and "Manual"
            //! ---------------------------
            if(itemXcomponent==Q_NULLPTR)
            {
                data.setValue(1);
                Property prop_Xcomponent("X component ",data,Property::PropertyGroup_Advanced);
                myCurNode->addProperty(prop_Xcomponent);
            }
            if(itemYcomponent==Q_NULLPTR)
            {
                data.setValue(1);
                Property prop_Ycomponent("Y component ",data,Property::PropertyGroup_Advanced);
                myCurNode->addProperty(prop_Ycomponent);
            }
            if(itemZcomponent==Q_NULLPTR)
            {
                data.setValue(1);
                Property prop_Zcomponent("Z component ",data,Property::PropertyGroup_Advanced);
                myCurNode->addProperty(prop_Zcomponent);
            }
            if(itemXrotation==Q_NULLPTR)
            {
                data.setValue(1);
                Property prop_Xrot("X rotation",data,Property::PropertyGroup_Advanced);
                myCurNode->addProperty(prop_Xrot);
            }
            if(itemYrotation==Q_NULLPTR)
            {
                data.setValue(1);
                Property prop_Yrot("Y rotation",data,Property::PropertyGroup_Advanced);
                myCurNode->addProperty(prop_Yrot);
            }
            if(itemZrotation==Q_NULLPTR)
            {
                data.setValue(1);
                Property prop_Zrot("Z rotation",data,Property::PropertyGroup_Advanced);
                myCurNode->addProperty(prop_Zrot);
            }
        }
            break;
        }
    }
        break;
    }
}

//! --------------------------------------------------------
//! function: handleFluxConvergenceChanged
//! details:  0 => Remove; 1 => On; 2 => Program controlled
//! --------------------------------------------------------
void DetailViewer::handleFluxConvergenceChanged()
{
    cout<<"DetailViewer::handleFluxConvergenceChanged()->____function called____"<<endl;

    QExtendedStandardItem *itemFluxConvergence = myCurNode->getPropertyItem("Flux convergence");
    int val = itemFluxConvergence->data(Qt::UserRole).value<Property>().getData().toInt();
    QVariant data;
    //! --------------------------------------------------------------
    //! val = 0 => Remove val = 1 => On val = 2 => Program controlled
    //! --------------------------------------------------------------
    QExtendedStandardItem *item_Value = myCurNode->getPropertyItem("--Value");
    QExtendedStandardItem *item_q_alpha_0 = myCurNode->getPropertyItem("--q_alpha_0");
    QExtendedStandardItem *item_R_alpha_n = myCurNode->getPropertyItem("--R_alpha_n");
    QExtendedStandardItem *item_R_alpha_P = myCurNode->getPropertyItem("--R_alpha_P");
    QExtendedStandardItem *item_R_alpha_l = myCurNode->getPropertyItem("--R_alpha_l");
    QExtendedStandardItem *item_epsilon_alpha = myCurNode->getPropertyItem("--epsilon_alpha");

    switch(val)
    {
    case 0:
    {
        //! -------------------
        //! remove
        //! -------------------
        if(item_Value!=Q_NULLPTR) myCurNode->removeProperty("--Value");                      //! row 1
        if(item_q_alpha_0!=Q_NULLPTR) myCurNode->removeProperty("--q_alpha_0");              //! row 2
        if(item_R_alpha_n!=Q_NULLPTR) myCurNode->removeProperty("--R_alpha_n");              //! row 3
        if(item_R_alpha_P!=Q_NULLPTR) myCurNode->removeProperty("--R_alpha_P");              //! row 4
        if(item_R_alpha_l!=Q_NULLPTR) myCurNode->removeProperty("--R_alpha_l");              //! row 5
        if(item_epsilon_alpha!=Q_NULLPTR) myCurNode->removeProperty("--epsilon_alpha");      //! row 6

        data.setValue(0.0);
        Property prop_Value("--Value",data,Property::PropertyGroup_ConvergenceCriteria);
        myCurNode->addProperty(prop_Value,1);
        myCurNode->getPropertyItem("--Value")->setEditable(false);
        this->setRowHidden(myCurNode->getPropertyItem("--Value")->index().row(),myCurNode->getPropertyItem("--Value")->index().parent(),true);

        data.setValue(0.0);
        Property prop_q_alpha_0("--q_alpha_0",data,Property::PropertyGroup_ConvergenceCriteria);
        myCurNode->addProperty(prop_q_alpha_0,2);
        myCurNode->getPropertyItem("--q_alpha_0")->setEditable(false);
        this->setRowHidden(myCurNode->getPropertyItem("--q_alpha_0")->index().row(),myCurNode->getPropertyItem("--q_alpha_0")->index().parent(),true);

        data.setValue(1e10);
        Property prop_R_alpha_n("--R_alpha_n",data,Property::PropertyGroup_ConvergenceCriteria);
        myCurNode->addProperty(prop_R_alpha_n,3);
        myCurNode->getPropertyItem("--R_alpha_n")->setEditable(false);
        this->setRowHidden(myCurNode->getPropertyItem("--R_alpha_n")->index().row(),myCurNode->getPropertyItem("--R_alpha_n")->index().parent(),true);

        data.setValue(1e10);
        Property prop_R_alpha_P("--R_alpha_P",data,Property::PropertyGroup_ConvergenceCriteria);
        myCurNode->addProperty(prop_R_alpha_P,4);
        myCurNode->getPropertyItem("--R_alpha_P")->setEditable(false);
        this->setRowHidden(myCurNode->getPropertyItem("--R_alpha_P")->index().row(),myCurNode->getPropertyItem("--R_alpha_P")->index().parent(),true);

        data.setValue(1e10);
        Property prop_c_alpha_i_max("--R_alpha_l",data,Property::PropertyGroup_ConvergenceCriteria);
        myCurNode->addProperty(prop_c_alpha_i_max,5);
        myCurNode->getPropertyItem("--R_alpha_l")->setEditable(false);
        this->setRowHidden(myCurNode->getPropertyItem("--R_alpha_l")->index().row(),myCurNode->getPropertyItem("--R_alpha_l")->index().parent(),true);

        data.setValue(1e-5);
        Property prop_epsilon_alpha("--epsilon_alpha",data,Property::PropertyGroup_ConvergenceCriteria);
        myCurNode->addProperty(prop_epsilon_alpha,6);
        myCurNode->getPropertyItem("--epsilon_alpha")->setEditable(false);
        this->setRowHidden(myCurNode->getPropertyItem("--epsilon_alpha")->index().row(),myCurNode->getPropertyItem("--epsilon_alpha")->index().parent(),true);
    }
        break;

    case 1:
    {   
        //! -------------------
        //! On
        //! -------------------
        if(item_Value!=Q_NULLPTR) myCurNode->removeProperty("--Value");
        if(item_q_alpha_0!=Q_NULLPTR) myCurNode->removeProperty("--q_alpha_0");
        if(item_R_alpha_n!=Q_NULLPTR) myCurNode->removeProperty("--R_alpha_n");
        if(item_R_alpha_n!=Q_NULLPTR) myCurNode->removeProperty("--R_alpha_P");
        if(item_R_alpha_l!=Q_NULLPTR) myCurNode->removeProperty("--R_alpha_l");
        if(item_epsilon_alpha!=Q_NULLPTR) myCurNode->removeProperty("--epsilon_alpha");

        data.setValue(0.0);
        Property prop_Value("--Value",data,Property::PropertyGroup_ConvergenceCriteria);
        myCurNode->addProperty(prop_Value,1);
        myCurNode->getPropertyItem("--Value")->setEditable(true);
        this->setRowHidden(myCurNode->getPropertyItem("--Value")->index().row(),myCurNode->getPropertyItem("--Value")->index().parent(),false);

        data.setValue(0.0);
        Property prop_q_alpha_0("--q_alpha_0",data,Property::PropertyGroup_ConvergenceCriteria);
        myCurNode->addProperty(prop_q_alpha_0,2);
        myCurNode->getPropertyItem("--q_alpha_0")->setEditable(true);
        this->setRowHidden(myCurNode->getPropertyItem("--q_alpha_0")->index().row(),myCurNode->getPropertyItem("--q_alpha_0")->index().parent(),false);

        data.setValue(0.001);
        Property prop_R_alpha_n("--R_alpha_n",data,Property::PropertyGroup_ConvergenceCriteria);
        myCurNode->addProperty(prop_R_alpha_n,3);
        myCurNode->getPropertyItem("--R_alpha_n")->setEditable(true);
        this->setRowHidden(myCurNode->getPropertyItem("--R_alpha_n")->index().row(),myCurNode->getPropertyItem("--R_alpha_n")->index().parent(),false);

        data.setValue(0.02);
        Property prop_R_alpha_P("--R_alpha_P",data,Property::PropertyGroup_ConvergenceCriteria);
        myCurNode->addProperty(prop_R_alpha_P,4);
        myCurNode->getPropertyItem("--R_alpha_P")->setEditable(true);
        this->setRowHidden(myCurNode->getPropertyItem("--R_alpha_P")->index().row(),myCurNode->getPropertyItem("--R_alpha_P")->index().parent(),false);

        data.setValue(1e-8);
        Property prop_c_alpha_i_max("--R_alpha_l",data,Property::PropertyGroup_ConvergenceCriteria);
        myCurNode->addProperty(prop_c_alpha_i_max,5);
        myCurNode->getPropertyItem("--R_alpha_l")->setEditable(true);
        this->setRowHidden(myCurNode->getPropertyItem("--R_alpha_l")->index().row(),myCurNode->getPropertyItem("--R_alpha_l")->index().parent(),false);

        data.setValue(1e-5);
        Property prop_epsilon_alpha("--epsilon_alpha",data,Property::PropertyGroup_ConvergenceCriteria);
        myCurNode->addProperty(prop_epsilon_alpha,6);
        myCurNode->getPropertyItem("--epsilon_alpha")->setEditable(true);
        this->setRowHidden(myCurNode->getPropertyItem("--epsilon_alpha")->index().row(),myCurNode->getPropertyItem("--epsilon_alpha")->index().parent(),false);
    }
        break;

    case 2:
    {
        //! -------------------
        //! Program controlled
        //! -------------------
        if(item_Value!=Q_NULLPTR) myCurNode->removeProperty("--Value");
        if(item_q_alpha_0!=Q_NULLPTR) myCurNode->removeProperty("--q_alpha_0");
        if(item_R_alpha_n!=Q_NULLPTR) myCurNode->removeProperty("--R_alpha_n");
        if(item_R_alpha_n!=Q_NULLPTR) myCurNode->removeProperty("--R_alpha_P");
        if(item_R_alpha_l!=Q_NULLPTR) myCurNode->removeProperty("--R_alpha_l");
        if(item_epsilon_alpha!=Q_NULLPTR) myCurNode->removeProperty("--epsilon_alpha");

        data.setValue(0.0);
        Property prop_Value("--Value",data,Property::PropertyGroup_ConvergenceCriteria);
        myCurNode->addProperty(prop_Value,1);
        myCurNode->getPropertyItem("--Value")->setEditable(false);
        this->setRowHidden(myCurNode->getPropertyItem("--Value")->index().row(),myCurNode->getPropertyItem("--Value")->index().parent(),true);

        data.setValue(0.0);
        Property prop_q_alpha_0("--q_alpha_0",data,Property::PropertyGroup_ConvergenceCriteria);
        myCurNode->addProperty(prop_q_alpha_0,2);
        myCurNode->getPropertyItem("--q_alpha_0")->setEditable(false);
        this->setRowHidden(myCurNode->getPropertyItem("--q_alpha_0")->index().row(),myCurNode->getPropertyItem("--q_alpha_0")->index().parent(),true);

        data.setValue(0.001);
        Property prop_R_alpha_n("--R_alpha_n",data,Property::PropertyGroup_ConvergenceCriteria);
        myCurNode->addProperty(prop_R_alpha_n,3);
        myCurNode->getPropertyItem("--R_alpha_n")->setEditable(false);
        this->setRowHidden(myCurNode->getPropertyItem("--R_alpha_n")->index().row(),myCurNode->getPropertyItem("--R_alpha_n")->index().parent(),true);

        data.setValue(0.02);
        Property prop_R_alpha_P("--R_alpha_P",data,Property::PropertyGroup_ConvergenceCriteria);
        myCurNode->addProperty(prop_R_alpha_P,4);
        myCurNode->getPropertyItem("--R_alpha_P")->setEditable(false);
        this->setRowHidden(myCurNode->getPropertyItem("--R_alpha_P")->index().row(),myCurNode->getPropertyItem("--R_alpha_P")->index().parent(),true);

        data.setValue(1e-8);
        Property prop_c_alpha_i_max("--R_alpha_l",data,Property::PropertyGroup_ConvergenceCriteria);
        myCurNode->addProperty(prop_c_alpha_i_max,5);
        myCurNode->getPropertyItem("--R_alpha_l")->setEditable(false);
        this->setRowHidden(myCurNode->getPropertyItem("--R_alpha_l")->index().row(),myCurNode->getPropertyItem("--R_alpha_l")->index().parent(),true);

        data.setValue(1e10);
        Property prop_epsilon_alpha("--epsilon_alpha",data,Property::PropertyGroup_ConvergenceCriteria);
        myCurNode->addProperty(prop_epsilon_alpha,6);
        myCurNode->getPropertyItem("--epsilon_alpha")->setEditable(false);
        this->setRowHidden(myCurNode->getPropertyItem("--epsilon_alpha")->index().row(),myCurNode->getPropertyItem("--epsilon_alpha")->index().parent(),true);
    }
        break;
    }

    //! -----------------------------------------------------
    //! update the field parameters at the current time step
    //! -----------------------------------------------------
    int  N = myCurNode->getPropertyItem("Current step number")->data(Qt::UserRole).value<Property>().getData().toInt();
    CustomTableModel *tabModel = myCurNode->getTabularDataModel();

    data = tabModel->dataRC(N,4,Qt::EditRole).value<QVariant>();
    if(data.canConvert<QVector<double>>()) cout<<"____can convert____"<<endl;
    QVector<double> fieldParameters = data.value<QVector<double>>();

    //! ------
    //! pos 3
    //! ------
    double Value = myCurNode->getPropertyItem("--Value")->data(Qt::UserRole).value<Property>().getData().toDouble();
    fieldParameters.replace(3,Value);

    //! ------
    //! pos 2
    //! ------
    double q_alpha_0 = myCurNode->getPropertyItem("--q_alpha_0")->data(Qt::UserRole).value<Property>().getData().toDouble();
    fieldParameters.replace(2,q_alpha_0);

    //! ------
    //! pos 0
    //! ------
    double R_alpha_n = myCurNode->getPropertyItem("--R_alpha_n")->data(Qt::UserRole).value<Property>().getData().toDouble();
    fieldParameters.replace(0,R_alpha_n);

    //! ------
    //! pos 4
    //! ------
    double R_alpha_P = myCurNode->getPropertyItem("--R_alpha_P")->data(Qt::UserRole).value<Property>().getData().toDouble();
    fieldParameters.replace(4,R_alpha_P);

    //! ------
    //! pos 7
    //! ------
    double R_alpha_l = myCurNode->getPropertyItem("--R_alpha_l")->data(Qt::UserRole).value<Property>().getData().toDouble();
    fieldParameters.replace(7,R_alpha_l);

    //! ------
    //! pos 5
    //! ------
    double epsilon_alpha = myCurNode->getPropertyItem("--epsilon_alpha")->data(Qt::UserRole).value<Property>().getData().toDouble();
    fieldParameters.replace(5,epsilon_alpha);

    data.setValue(fieldParameters);
    tabModel->setDataRC(data,N,4);

    //! -----------
    //! column # 5
    //! -----------
    int fluxConvergence = myCurNode->getPropertyItem("Flux convergence")->data(Qt::UserRole).value<Property>().getData().toInt();
    data.setValue(fluxConvergence);
    tabModel->setDataRC(data,N,5);

    //! -----------
    //! column # 6
    //! -----------
    int solutionConvergence = myCurNode->getPropertyItem("Solution convergence")->data(Qt::UserRole).value<Property>().getData().toInt();
    data.setValue(solutionConvergence);
    tabModel->setDataRC(data,N,6);
}

//! --------------------------------------------------------
//! function: handleSolutionConvergenceChanged
//! details:  0 => Remove; 1 => On; 2 => Program controlled
//! --------------------------------------------------------
void DetailViewer::handleSolutionConvergenceChanged()
{
    cout<<"DetailViewer::handleSolutionConvergenceChanged()->____function called____"<<endl;

    QExtendedStandardItem *itemSolutionConvergence = myCurNode->getPropertyItem("Solution convergence");
    int val = itemSolutionConvergence->data(Qt::UserRole).value<Property>().getData().toInt();
    QVariant data;
    //! --------------------------------------------------------------
    //! val = 0 => Remove val = 1 => On val = 2 => Program controlled
    //! --------------------------------------------------------------
    QExtendedStandardItem *item_C_alpha_n = myCurNode->getPropertyItem("--C_alpha_n");
    QExtendedStandardItem *item_C_alpha_epsilon = myCurNode->getPropertyItem("--C_alpha_epsilon");

    switch(val)
    {
    case 0:
    {
        //! -------------------
        //! Remove
        //! -------------------
        if(item_C_alpha_n!=Q_NULLPTR) myCurNode->removeProperty("--C_alpha_n");
        if(item_C_alpha_epsilon!=Q_NULLPTR) myCurNode->removeProperty("--C_alpha_epsilon");

        data.setValue(1e10);
        Property prop_C_alpha_n("--C_alpha_n",data,Property::PropertyGroup_ConvergenceCriteria);
        myCurNode->addProperty(prop_C_alpha_n,8);
        myCurNode->getPropertyItem("--C_alpha_n")->setEditable(false);
        this->setRowHidden(myCurNode->getPropertyItem("--C_alpha_n")->index().row(),myCurNode->getPropertyItem("--C_alpha_n")->index().parent(),true);

        data.setValue(1e10);
        Property prop_C_epsilon_alpha("--C_alpha_epsilon",data,Property::PropertyGroup_ConvergenceCriteria);
        myCurNode->addProperty(prop_C_epsilon_alpha,9);
        myCurNode->getPropertyItem("--C_alpha_epsilon")->setEditable(false);
        this->setRowHidden(myCurNode->getPropertyItem("--C_alpha_epsilon")->index().row(),myCurNode->getPropertyItem("--C_alpha_epsilon")->index().parent(),true);
    }
        break;

    case 1:
    {
        //! -------------------
        //! On
        //! -------------------
        if(item_C_alpha_n!=Q_NULLPTR) myCurNode->removeProperty("--C_alpha_n");
        if(item_C_alpha_epsilon!=Q_NULLPTR) myCurNode->removeProperty("--C_alpha_epsilon");

        data.setValue(0.02);
        Property prop_C_alpha_n("--C_alpha_n",data,Property::PropertyGroup_ConvergenceCriteria);
        myCurNode->addProperty(prop_C_alpha_n,8);
        myCurNode->getPropertyItem("--C_alpha_n")->setEditable(true);
        this->setRowHidden(myCurNode->getPropertyItem("--C_alpha_n")->index().row(),myCurNode->getPropertyItem("--C_alpha_n")->index().parent(),false);

        data.setValue(0.001);
        Property prop_C_epsilon_alpha("--C_alpha_epsilon",data,Property::PropertyGroup_ConvergenceCriteria);
        myCurNode->addProperty(prop_C_epsilon_alpha,9);
        myCurNode->getPropertyItem("--C_alpha_epsilon")->setEditable(true);
        this->setRowHidden(myCurNode->getPropertyItem("--C_alpha_epsilon")->index().row(),myCurNode->getPropertyItem("--C_alpha_epsilon")->index().parent(),false);
    }
        break;

    case 2:
    {
        //! -------------------
        //! Program controlled
        //! -------------------
        if(item_C_alpha_n!=Q_NULLPTR) myCurNode->removeProperty("--C_alpha_n");
        if(item_C_alpha_epsilon!=Q_NULLPTR) myCurNode->removeProperty("--C_alpha_epsilon");

        data.setValue(0.02);
        Property prop_C_alpha_n("--C_alpha_n",data,Property::PropertyGroup_ConvergenceCriteria);
        myCurNode->addProperty(prop_C_alpha_n,8);
        myCurNode->getPropertyItem("--C_alpha_n")->setEditable(false);
        this->setRowHidden(myCurNode->getPropertyItem("--C_alpha_n")->index().row(),myCurNode->getPropertyItem("--C_alpha_n")->index().parent(),true);

        data.setValue(0.001);
        Property prop_C_epsilon_alpha("--C_alpha_epsilon",data,Property::PropertyGroup_ConvergenceCriteria);
        myCurNode->addProperty(prop_C_epsilon_alpha,9);
        myCurNode->getPropertyItem("--C_alpha_epsilon")->setEditable(false);
        this->setRowHidden(myCurNode->getPropertyItem("--C_alpha_epsilon")->index().row(),myCurNode->getPropertyItem("--C_alpha_epsilon")->index().parent(),true);
    }
        break;
    }

    //! -----------------------------------------------------
    //! update the field parameters at the current time step
    //! -----------------------------------------------------
    int  N = myCurNode->getPropertyItem("Current step number")->data(Qt::UserRole).value<Property>().getData().toInt();
    CustomTableModel *tabModel = myCurNode->getTabularDataModel();

    cout<<"____current time step: "<<N<<"____"<<endl;
    data = tabModel->dataRC(N,4,Qt::EditRole).value<QVariant>();
    if(data.canConvert<QVector<double>>()) cout<<"____can convert____"<<endl;
    QVector<double> fieldParameters = data.value<QVector<double>>();

    //! ------
    //! pos 1
    //! ------
    double C_alpha_n = myCurNode->getPropertyItem("--C_alpha_n")->data(Qt::UserRole).value<Property>().getData().toDouble();
    fieldParameters.replace(1,C_alpha_n);

    //! ------
    //! pos 6
    //! ------
    double C_alpha_epsilon = myCurNode->getPropertyItem("--C_alpha_epsilon")->data(Qt::UserRole).value<Property>().getData().toDouble();
    fieldParameters.replace(6,C_alpha_epsilon);

    data.setValue(fieldParameters);
    tabModel->setDataRC(data,N,4);

    //! -----------
    //! column # 5
    //! -----------
    int fluxConvergence = myCurNode->getPropertyItem("Flux convergence")->data(Qt::UserRole).value<Property>().getData().toInt();
    data.setValue(fluxConvergence);
    tabModel->setDataRC(data,N,5);

    //! -----------
    //! column # 6
    //! -----------
    int solutionConvergence = myCurNode->getPropertyItem("Solution convergence")->data(Qt::UserRole).value<Property>().getData().toInt();
    data.setValue(solutionConvergence);
    tabModel->setDataRC(data,N,6);
}

//! --------------------------------
//! function: updateFieldParameters
//! details:
//! --------------------------------
void DetailViewer::updateFieldParameters()
{
    cout<<"DetailViewer::updateFieldParameters()->____function called____"<<endl;
    //! -----------------------------------------------------
    //! update the field parameters at the current time step
    //! -----------------------------------------------------
    int  N = myCurNode->getPropertyValue<int>("Current step number");
    CustomTableModel *tabModel = myCurNode->getTabularDataModel();

    cout<<"____current time step: "<<N<<"____"<<endl;
    QVariant data;
    data = tabModel->dataRC(N,4,Qt::EditRole).value<QVariant>();
    QVector<double> fieldParameters = data.value<QVector<double>>();

    //! ------
    //! pos 3
    //! ------
    double Value = myCurNode->getPropertyItem("--Value")->data(Qt::UserRole).value<Property>().getData().toDouble();
    fieldParameters.replace(3,Value);

    //! ------
    //! pos 2
    //! ------
    double q_alpha_0 = myCurNode->getPropertyItem("--q_alpha_0")->data(Qt::UserRole).value<Property>().getData().toDouble();
    fieldParameters.replace(2,q_alpha_0);

    //! ------
    //! pos 0
    //! ------
    double R_alpha_n = myCurNode->getPropertyItem("--R_alpha_n")->data(Qt::UserRole).value<Property>().getData().toDouble();
    fieldParameters.replace(0,R_alpha_n);

    //! ------
    //! pos 4
    //! ------
    double R_alpha_P = myCurNode->getPropertyItem("--R_alpha_P")->data(Qt::UserRole).value<Property>().getData().toDouble();
    fieldParameters.replace(4,R_alpha_P);

    //! ------
    //! pos 7
    //! ------
    double R_alpha_l = myCurNode->getPropertyItem("--R_alpha_l")->data(Qt::UserRole).value<Property>().getData().toDouble();
    fieldParameters.replace(7,R_alpha_l);

    //! ------
    //! pos 5
    //! ------
    double epsilon_alpha = myCurNode->getPropertyItem("--epsilon_alpha")->data(Qt::UserRole).value<Property>().getData().toDouble();
    fieldParameters.replace(5,epsilon_alpha);

    //! ------
    //! pos 1
    //! ------
    double C_alpha_n = myCurNode->getPropertyItem("--C_alpha_n")->data(Qt::UserRole).value<Property>().getData().toDouble();
    fieldParameters.replace(1,C_alpha_n);

    //! ------
    //! pos 6
    //! ------
    double C_alpha_epsilon = myCurNode->getPropertyItem("--C_alpha_epsilon")->data(Qt::UserRole).value<Property>().getData().toDouble();
    fieldParameters.replace(6,C_alpha_epsilon);

    data.setValue(fieldParameters);
    tabModel->setDataRC(data,N,4);

    //! -----------
    //! column # 5
    //! -----------
    int fluxConvergence = myCurNode->getPropertyItem("Flux convergence")->data(Qt::UserRole).value<Property>().getData().toInt();
    data.setValue(fluxConvergence);
    tabModel->setDataRC(data,N,5);

    //! -----------
    //! column # 6
    //! -----------
    int solutionConvergence = myCurNode->getPropertyItem("Solution convergence")->data(Qt::UserRole).value<Property>().getData().toInt();
    data.setValue(solutionConvergence);
    tabModel->setDataRC(data,N,6);
}

//! ---------------------------------------------
//! function: updateTimeIncrementationParameters
//! details:
//! ---------------------------------------------
void DetailViewer::updateTimeIncrementationParameters()
{
    cout<<"DetailViewer::updateTimeIncrementationParameters()->____function called____"<<endl;

    int  N = myCurNode->getPropertyValue<int>("Current step number");
    CustomTableModel *tabModel = myCurNode->getTabularDataModel();

    int timeIncrementationType = myCurNode->getPropertyValue<int>("Time incrementation");
    QVariant data;
    data = tabModel->dataRC(N,7,Qt::EditRole).value<QVariant>();
    QVector<int> timeIncrementationParameters = data.value<QVector<int>>();

    bool isHidden;
    int I_0,I_R,I_P,I_C,I_L,I_G,I_S,I_A,I_J,I_T;
    switch(timeIncrementationType)
    {
    case 0:
    {
        //! -----------------------------------------------------------------------
        //! put here the default CCX solver values or an "optimized" set of values
        //! for the current application
        //! -----------------------------------------------------------------------
        I_0 = 4;
        I_R = 8;
        I_P = 9;
        I_C = 16;
        I_L = 10;
        I_G = 4;
        I_S = -1;
        I_A = 5;
        I_J = -1;
        I_T = -1;

        myCurNode->removeProperty("I_0");
        myCurNode->removeProperty("I_R");
        myCurNode->removeProperty("I_P");
        myCurNode->removeProperty("I_C");
        myCurNode->removeProperty("I_L");
        myCurNode->removeProperty("I_G");
        myCurNode->removeProperty("I_S");
        myCurNode->removeProperty("I_A");
        myCurNode->removeProperty("I_J");
        myCurNode->removeProperty("I_T");

        data.setValue(I_0);
        Property property_I_0("I_0",data,Property::PropertyGroup_TimeIncrementation);
        myCurNode->addProperty(property_I_0);
        data.setValue(I_R);
        Property property_I_R("I_R",data,Property::PropertyGroup_TimeIncrementation);
        myCurNode->addProperty(property_I_R);
        data.setValue(I_P);
        Property property_I_P("I_P",data,Property::PropertyGroup_TimeIncrementation);
        myCurNode->addProperty(property_I_P);
        data.setValue(I_C);
        Property property_I_C("I_C",data,Property::PropertyGroup_TimeIncrementation);
        myCurNode->addProperty(property_I_C);
        data.setValue(I_L);
        Property property_I_L("I_L",data,Property::PropertyGroup_TimeIncrementation);
        myCurNode->addProperty(property_I_L);
        data.setValue(I_G);
        Property property_I_G("I_G",data,Property::PropertyGroup_TimeIncrementation);
        myCurNode->addProperty(property_I_G);
        data.setValue(I_S);
        Property property_I_S("I_S",data,Property::PropertyGroup_TimeIncrementation);
        myCurNode->addProperty(property_I_S);
        data.setValue(I_A);
        Property property_I_A("I_A",data,Property::PropertyGroup_TimeIncrementation);
        myCurNode->addProperty(property_I_A);
        data.setValue(I_J);
        Property property_I_J("I_J",data,Property::PropertyGroup_TimeIncrementation);
        myCurNode->addProperty(property_I_J);
        data.setValue(I_T);
        Property property_I_T("I_T",data,Property::PropertyGroup_TimeIncrementation);
        myCurNode->addProperty(property_I_T);

        isHidden = true;
    }
        break;

    case 1:
    {
        I_0 = myCurNode->getPropertyValue<int>("I_0");
        I_R = myCurNode->getPropertyValue<int>("I_R");
        I_P = myCurNode->getPropertyValue<int>("I_P");
        I_C = myCurNode->getPropertyValue<int>("I_C");
        I_L = myCurNode->getPropertyValue<int>("I_L");
        I_G = myCurNode->getPropertyValue<int>("I_G");
        I_S = myCurNode->getPropertyValue<int>("I_S");
        I_A = myCurNode->getPropertyValue<int>("I_A");
        I_J = myCurNode->getPropertyValue<int>("I_J");
        I_T = myCurNode->getPropertyValue<int>("I_T");

        isHidden = false;
    }
        break;
    }

    timeIncrementationParameters.replace(0,I_0);
    timeIncrementationParameters.replace(1,I_R);
    timeIncrementationParameters.replace(2,I_P);
    timeIncrementationParameters.replace(3,I_C);
    timeIncrementationParameters.replace(4,I_L);
    timeIncrementationParameters.replace(5,I_G);
    timeIncrementationParameters.replace(6,I_S);
    timeIncrementationParameters.replace(7,I_A);
    timeIncrementationParameters.replace(8,I_J);
    timeIncrementationParameters.replace(9,I_T);

    //! ------------------------
    //! modify the tabular data
    //! ------------------------
    data.setValue(timeIncrementationParameters);
    tabModel->setDataRC(data,N,7,Qt::EditRole);

    data.setValue(timeIncrementationType);
    tabModel->setDataRC(data,N,8,Qt::EditRole);

    this->setRowHidden(myCurNode->getPropertyItem("I_0")->index().row(),myCurNode->getPropertyItem("I_0")->index().parent(),isHidden);
    this->setRowHidden(myCurNode->getPropertyItem("I_R")->index().row(),myCurNode->getPropertyItem("I_R")->index().parent(),isHidden);
    this->setRowHidden(myCurNode->getPropertyItem("I_P")->index().row(),myCurNode->getPropertyItem("I_P")->index().parent(),isHidden);
    this->setRowHidden(myCurNode->getPropertyItem("I_C")->index().row(),myCurNode->getPropertyItem("I_C")->index().parent(),isHidden);
    this->setRowHidden(myCurNode->getPropertyItem("I_L")->index().row(),myCurNode->getPropertyItem("I_L")->index().parent(),isHidden);
    this->setRowHidden(myCurNode->getPropertyItem("I_G")->index().row(),myCurNode->getPropertyItem("I_G")->index().parent(),isHidden);
    this->setRowHidden(myCurNode->getPropertyItem("I_S")->index().row(),myCurNode->getPropertyItem("I_S")->index().parent(),isHidden);
    this->setRowHidden(myCurNode->getPropertyItem("I_A")->index().row(),myCurNode->getPropertyItem("I_A")->index().parent(),isHidden);
    this->setRowHidden(myCurNode->getPropertyItem("I_J")->index().row(),myCurNode->getPropertyItem("I_J")->index().parent(),isHidden);
    this->setRowHidden(myCurNode->getPropertyItem("I_T")->index().row(),myCurNode->getPropertyItem("I_T")->index().parent(),isHidden);
}

//! -------------------------------
//! function: updateCutBackFactors
//! details:
//! -------------------------------
void DetailViewer::updateCutBackFactors()
{
    cout<<"DetailViewer::updateCutBackFactors()->____function called____"<<endl;
    int  N = myCurNode->getPropertyValue<int>("Current step number");
    CustomTableModel *tabModel = myCurNode->getTabularDataModel();

    QVariant data;
    data = tabModel->dataRC(N,TABULAR_DATA_CUTBACK_PARAMETERS_COLUMN,Qt::EditRole).value<QVariant>();
    QVector<double> cutBackParameters = data.value<QVector<double>>();

    bool isHidden;
    double D_f,D_C,D_B,D_A,D_S,D_H,D_D,W_G;

    int cutBackType = myCurNode->getPropertyValue<int>("Cutback factors");
    switch(cutBackType)
    {
    case 0:
    {
        D_f = 0.25;
        D_C = 0.5;
        D_B = 0.75;
        D_A = 0.85;
        D_S = -1;   //! unused
        D_H = -1;   //! unused
        D_D = 1.5;
        W_G = -1;   //! unused

        myCurNode->removeProperty("D_f");
        myCurNode->removeProperty("D_C");
        myCurNode->removeProperty("D_B");
        myCurNode->removeProperty("D_A");
        myCurNode->removeProperty("D_S");
        myCurNode->removeProperty("D_H");
        myCurNode->removeProperty("D_D");
        myCurNode->removeProperty("W_G");

        data.setValue(D_f);
        Property property_D_f("D_f",data,Property::PropertyGroup_CutBack);
        myCurNode->addProperty(property_D_f);

        data.setValue(D_C);
        Property property_D_C("D_C",data,Property::PropertyGroup_CutBack);
        myCurNode->addProperty(property_D_C);

        data.setValue(D_B);
        Property property_D_B("D_B",data,Property::PropertyGroup_CutBack);
        myCurNode->addProperty(property_D_B);

        data.setValue(D_A);
        Property property_D_A("D_A",data,Property::PropertyGroup_CutBack);
        myCurNode->addProperty(property_D_A);

        data.setValue(D_S);
        Property property_D_S("D_S",data,Property::PropertyGroup_CutBack);
        myCurNode->addProperty(property_D_S);

        data.setValue(D_H);
        Property property_D_H("D_H",data,Property::PropertyGroup_CutBack);
        myCurNode->addProperty(property_D_H);

        data.setValue(D_D);
        Property property_D_D("D_D",data,Property::PropertyGroup_CutBack);
        myCurNode->addProperty(property_D_D);

        data.setValue(W_G);
        Property property_W_G("W_G",data,Property::PropertyGroup_CutBack);
        myCurNode->addProperty(property_W_G);

        isHidden = true;
    }
        break;

    case 1:
    {
        D_f = myCurNode->getPropertyValue<double>("D_f");
        D_C = myCurNode->getPropertyValue<double>("D_C");
        D_B = myCurNode->getPropertyValue<double>("D_B");
        D_A = myCurNode->getPropertyValue<double>("D_A");
        D_S = myCurNode->getPropertyValue<double>("D_S");
        D_H = myCurNode->getPropertyValue<double>("D_H");
        D_D = myCurNode->getPropertyValue<double>("D_D");
        W_G = myCurNode->getPropertyValue<double>("W_G");

        isHidden = false;
    }
        break;
    }

    this->setRowHidden(myCurNode->getPropertyItem("D_f")->index().row(),myCurNode->getPropertyItem("D_f")->index().parent(),isHidden);
    this->setRowHidden(myCurNode->getPropertyItem("D_C")->index().row(),myCurNode->getPropertyItem("D_C")->index().parent(),isHidden);
    this->setRowHidden(myCurNode->getPropertyItem("D_B")->index().row(),myCurNode->getPropertyItem("D_B")->index().parent(),isHidden);
    this->setRowHidden(myCurNode->getPropertyItem("D_A")->index().row(),myCurNode->getPropertyItem("D_A")->index().parent(),isHidden);
    this->setRowHidden(myCurNode->getPropertyItem("D_S")->index().row(),myCurNode->getPropertyItem("D_S")->index().parent(),isHidden);
    this->setRowHidden(myCurNode->getPropertyItem("D_H")->index().row(),myCurNode->getPropertyItem("D_H")->index().parent(),isHidden);
    this->setRowHidden(myCurNode->getPropertyItem("D_D")->index().row(),myCurNode->getPropertyItem("D_D")->index().parent(),isHidden);
    this->setRowHidden(myCurNode->getPropertyItem("W_G")->index().row(),myCurNode->getPropertyItem("W_G")->index().parent(),isHidden);

    cutBackParameters.replace(0,D_f);
    cutBackParameters.replace(1,D_C);
    cutBackParameters.replace(2,D_B);
    cutBackParameters.replace(3,D_A);
    cutBackParameters.replace(4,D_S);
    cutBackParameters.replace(5,D_H);
    cutBackParameters.replace(6,D_D);
    cutBackParameters.replace(7,W_G);

    //! ------------------------
    //! modify the tabular data
    //! ------------------------
    data.setValue(cutBackType);
    tabModel->setDataRC(data,N,TABULAR_DATA_CUTBACK_TYPE_COLUMN,Qt::EditRole);

    data.setValue(cutBackParameters);
    tabModel->setDataRC(data,N,TABULAR_DATA_CUTBACK_PARAMETERS_COLUMN,Qt::EditRole);
}

//! -----------------------------
//! function: update line search
//! details:
//! -----------------------------
void DetailViewer::updateLineSearch()
{
    cout<<"DetailViewer::updateLineSearch()->____function called____"<<endl;
    int  N = myCurNode->getPropertyItem("Current step number")->data(Qt::UserRole).value<Property>().getData().toInt();
    CustomTableModel *tabModel = myCurNode->getTabularDataModel();

    //! -----------------------
    //! line search parameters
    //! -----------------------
    QVariant data;
    data = tabModel->dataRC(N,TABULAR_DATA_LINE_SEARCH_PARAMETERS_COLUMN,Qt::EditRole).value<QVariant>();
    QVector<double> lineSearchParameters = data.value<QVector<double>>();

    int lineSearch = myCurNode->getPropertyValue<int>("Line search");
    switch(lineSearch)
    {
    //! ---------------------
    //! "Program controlled"
    //! ---------------------
    case 0:
    {
        cout<<"____updating using program controlled____"<<endl;

        //! -------------------------------------
        //! CCX defaults or application defaults
        //! -------------------------------------
        double minLineSearch = 0.25;
        double maxLineSearch = 1.01;

        //! update the Detail viewer controls
        data.setValue(minLineSearch);
        Property prop_minValue("Min value",data,Property::PropertyGroup_LineSearch);
        data.setValue(maxLineSearch);
        Property prop_maxValue("Max value",data,Property::PropertyGroup_LineSearch);

        myCurNode->removeProperty("Min value");
        myCurNode->removeProperty("Max value");

        myCurNode->addProperty(prop_minValue);
        myCurNode->addProperty(prop_maxValue);

        //! update the table
        QVector<double> newlineSearchParameters;
        newlineSearchParameters.push_back(0.25);
        newlineSearchParameters.push_back(1.01);
        data.setValue(newlineSearchParameters);
        tabModel->setDataRC(data,N,TABULAR_DATA_LINE_SEARCH_PARAMETERS_COLUMN,Qt::EditRole);
    }
        break;

    //! ---------
    //! "Custom"
    //! ---------
    case 1:
    {
        cout<<"____updating using custom____"<<endl;
        double minVal = myCurNode->getPropertyValue<double>("Min value");
        double maxVal = myCurNode->getPropertyValue<double>("Max value");
        lineSearchParameters.replace(0,minVal);
        lineSearchParameters.replace(1,maxVal);
        data.setValue(lineSearchParameters);
        tabModel->setDataRC(data,N,TABULAR_DATA_LINE_SEARCH_PARAMETERS_COLUMN,Qt::EditRole);
    }
        break;
    }

    //! ------------------
    //! handle visibility
    //! ------------------
    bool isHidden = true;
    if(lineSearch == 1) isHidden = false;
    QStandardItem *item = myCurNode->getPropertyItem("Min value");
    this->setRowHidden(item->index().row(),item->index().parent(),isHidden);
    item = myCurNode->getPropertyItem("Max value");
    this->setRowHidden(item->index().row(),item->index().parent(),isHidden);
}

//! --------------------------------------
//! function: handleOutputControlsChanged
//! details:
//! --------------------------------------
void DetailViewer::handleOutputControlsChanged()
{
    cout<<"DetailViewer::handleOutputControlsChanged()->____function called____"<<endl;

    //! row of the table - current step number
    int currentStep = myCurNode->getPropertyValue<int>("Current step number");

    //! ------------------------------------------------
    //! update the table using the detail viewer values
    //! ------------------------------------------------
    int stressFlag = myCurNode->getPropertyValue<int>("Stress");
    int strainFlag = myCurNode->getPropertyValue<int>("Strain");
    int reactionForcesFlag = myCurNode->getPropertyValue<int>("Reaction forces");
    int contactDataFlag = myCurNode->getPropertyValue<int>("Contact data");

    QVector<int> flags;
    flags.push_back(stressFlag);
    flags.push_back(strainFlag);
    flags.push_back(reactionForcesFlag);
    flags.push_back(contactDataFlag);

    //! -----------------
    //! update the table
    //! -----------------
    QVariant data;
    data.setValue(flags);
    myCurNode->getTabularDataModel()->setDataRC(data,currentStep,TABULAR_DATA_OUPUT_CONTROLS_COLUMN,Qt::EditRole);
}

//! --------------------------------------------
//! function: handleStoreResultsAtChanged()
//! details:  add/remove the "--Value" selector
//! --------------------------------------------
void DetailViewer::handleStoreResultsAtChanged()
{
    cout<<"DetailViewer::handleStoreResultsAtChanged()->____function called____"<<endl;

    //! -------------------------------------------------------------------------------------
    //! add/remove the control "--Recurrence rate" depending on the "Store results at" value
    //! -------------------------------------------------------------------------------------
    int storeResultsAt = myCurNode->getPropertyValue<int>("Store results at");
    cout<<"DetailViewer::handleStoreResultsAtChanged()->____"<<storeResultsAt<<"____"<<endl;

    switch(storeResultsAt)
    {
    case 0:
    case 1:
    {
        //! --------------------------------------------------
        //! "0" => "All time points" "1" => "Last time point"
        //! remove the "--Recurrence rate" selector
        //! --------------------------------------------------
        myCurNode->removeProperty("--Recurrence rate");
    }
        break;

    case 2:
    {
        //! -------------------------------------
        //! "2" => "Specified recurrence rate"
        //! add the "--Recurrence rate" selector
        //! -------------------------------------
        int recurrenceRate = 1;     //! set a default value for "--Recurrence rate"
        QVariant data;
        data.setValue(recurrenceRate);
        Property prop_recurrenceRate("--Recurrence rate",data,Property::PropertyGroup_OutputSettings);
        myCurNode->addProperty(prop_recurrenceRate);
    }
        break;
    }

    //! row of the table - current step number
    int currentStep = myCurNode->getPropertyValue<int>("Current step number");

    //! --------------------------------------------------
    //! update the table: the second value is FREQUENCY
    //! "All time points"           => in table (0,1)
    //! "Last time point"           => in table (1,9999)
    //! "Specified recurrence rate" => in table (2,value)
    //! --------------------------------------------------
    switch(storeResultsAt)
    {
    case 0: cout<<"____store results at \"All time points\"____"<<endl; break;
    case 1: cout<<"____store results at \"Last time point\"____"<<endl; break;
    case 2: cout<<"____store results at \"Specified recurrence rate\"____"<<endl; break;
    }

    QVector<int> flags1;
    flags1.push_back(storeResultsAt);
    int FREQUENCY;

    if(storeResultsAt == 0) FREQUENCY = 1;
    if(storeResultsAt == 1) FREQUENCY = 999999;
    if(storeResultsAt == 2) FREQUENCY = myCurNode->getPropertyValue<int>("--Recurrence rate");
    flags1.push_back(FREQUENCY);

    QVariant data; data.setValue(flags1);
    myCurNode->getTabularDataModel()->setDataRC(data,currentStep,TABULAR_DATA_STORE_RESULTS_AT_COLUMN,Qt::EditRole);
    cout<<"____writing pair: ("<<flags1.at(0)<<", "<<flags1.at(1)<<")____"<<endl;
}

//! ----------------------------------------------------------------
//! function: setPropertyVisible
//! details:  hide the row of the model containing a given property
//! ----------------------------------------------------------------
void DetailViewer::setPropertyVisible(const QString &propertyName, bool visible)
{
    QExtendedStandardItem *item = myCurNode->getPropertyItem(propertyName);
    if(item==Q_NULLPTR)
    {
        cerr<<"DetailViewer::setPropertyVisible()->____cannot make the property \""<<propertyName.toStdString()<<"\" "<<
              (visible == true? "visible": "hidden")<<" since not present____"<<endl;
    }
    this->setRowHidden(item->index().row(),item->parent()->index(),visible);
}

//! ---------------------------
//! function: expandSeparator
//! details:
//! ---------------------------
void DetailViewer::expandSeparator(const QString &separatorName)
{
    QStandardItemModel *nodeModel = myCurNode->getModel();
    for(int row=0; row<nodeModel->rowCount(); row++)
    {
       QStandardItem *childItem = nodeModel->item(row,0);
       if(childItem->data(Qt::UserRole+1).toString()==separatorName)
       {
           this->expand(childItem->index());
           break;
       }
    }
}

//! ------------------------------
//! function: collapseSeparator()
//! details:  collapse separator
//! ------------------------------
void DetailViewer::collapseSeparator(const QString &separatorName)
{
    QStandardItemModel *nodeModel = myCurNode->getModel();
    for(int row=0; row<nodeModel->rowCount(); row++)
    {
        QStandardItem *childItem = nodeModel->item(row,0);
        if(childItem->data(Qt::UserRole+1).toString()==separatorName)
        {
            this->collapse(childItem->index());
            break;
        }
    }
}

//! -------------------------------------------------------------
//! function: handleScaleChanged
//! details:  handle all the properties of the "Color box" group
//! -------------------------------------------------------------
void DetailViewer::handleColorBoxScaleChanged()
{
    //! --------------
    //! block signals
    //! --------------
    SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
    disconnect(myCurNode->getModel(),SIGNAL(itemChanged(QStandardItem*)),sm,SLOT(handleItemChange(QStandardItem*)));

    static int firstCall;
    firstCall++;
    cout<<"DetailViewer::handleColorBoxScaleChanged()->____function called: "<<firstCall<<"____"<<endl;
    static double min, max;
    if(firstCall==1)
    {
        min = -100;
        max = 100;
    }
    else
    {
        if(myCurNode->getPropertyItem("Min")!=Q_NULLPTR) min = myCurNode->getPropertyValue<double>("Min");
        if(myCurNode->getPropertyItem("Max")!=Q_NULLPTR) max = myCurNode->getPropertyValue<double>("Max");
    }
    int scaleType = myCurNode->getPropertyValue<int>("Scale type");
    switch(scaleType)
    {
    case 0:
    {
        //! -----------------
        //! 0 => "Autoscale"
        //! -----------------
        if(myCurNode->getPropertyItem("Min")!=Q_NULLPTR) min = myCurNode->getPropertyValue<double>("Min");
        if(myCurNode->getPropertyItem("Max")!=Q_NULLPTR) max = myCurNode->getPropertyValue<double>("Max");
        myCurNode->removeProperty("Min");
        myCurNode->removeProperty("Max");
    }
        break;

    case 1:
    {
        //! ---------------------------------------------------
        //! 1 => "Custom scale"
        //! automatically correct errors in "Min" "Max" values
        //! ---------------------------------------------------
        if(min>max) min = max-0.10*min;
        if(max<min) max = min+0.10*max;
        if(max==min) min = max-0.10*max;
        QVariant data;
        data.setValue(min);
        Property prop_min("Min",data,Property::PropertyGroup_ColorBox);
        data.setValue(max);
        Property prop_max("Max",data,Property::PropertyGroup_ColorBox);
        myCurNode->removeProperty("Min");
        myCurNode->removeProperty("Max");
        myCurNode->addProperty(prop_min,1);
        myCurNode->addProperty(prop_max,2);
    }
        break;
    }

    //! ----------------
    //! unblock signals
    //! ----------------
    connect(myCurNode->getModel(),SIGNAL(itemChanged(QStandardItem*)),sm,SLOT(handleItemChange(QStandardItem*)));
}

//! -----------------------------------------------------------
//! function: handleMeshDefeaturingChanged()
//! details:  enable/disable mesh healing/defeaturing controls
//! -----------------------------------------------------------
void DetailViewer::handleMeshDefeaturingChanged()
{
    cout<<"DetailViewer::handleMeshDefeaturingChanged()->____function called____"<<endl;
    static int wasMeshHealingOn;
    static int wasMeshSimplificationOn;
    static int meshSimplificationBy;
    static int level;
    static double pairDistance;

    if(myCurNode->getPropertyItem("Mesh healing")!=Q_NULLPTR) wasMeshHealingOn = myCurNode->getPropertyValue<int>("Mesh healing");
    if(myCurNode->getPropertyItem("Mesh simplification")!=Q_NULLPTR) wasMeshSimplificationOn = myCurNode->getPropertyValue<int>("Mesh simplification");
    if(myCurNode->getPropertyItem("By")!=Q_NULLPTR) meshSimplificationBy = myCurNode->getPropertyValue<int>("By");
    if(myCurNode->getPropertyItem("Level")) level = myCurNode->getPropertyValue<int>("Level");
    if(myCurNode->getPropertyItem("Pair distance")) pairDistance = myCurNode->getPropertyValue<double>("Pair distance");

    int meshDefeaturingStatus = myCurNode->getPropertyValue<int>("Mesh defeaturing");
    switch(meshDefeaturingStatus)
    {
    case 0:
    {
        //! remove all the properties (existing and not existing)
        QList<QString> propertiesToRemove;
        propertiesToRemove<<"Mesh healing"<<"Mesh simplification"<<"By"<<"Level"<<"Pair distance";
        for(int i=0; i<propertiesToRemove.length(); i++) myCurNode->removeProperty(propertiesToRemove.at(i));
    }
        break;

    case 1:
    {
        QVariant data;

        //! rebuild "Mesh healing"
        if(wasMeshHealingOn==1) data.setValue(1);
        else data.setValue(0);
        Property prop_meshHealing("Mesh healing",data,Property::PropertyGroup_Advanced);
        myCurNode->addProperty(prop_meshHealing,2);

        //! rebuild "Mesh simplification"
        if(wasMeshSimplificationOn==1) data.setValue(1);
        else data.setValue(0);
        Property prop_meshSimplification("Mesh simplification",data,Property::PropertyGroup_Advanced);
        myCurNode->addProperty(prop_meshSimplification,3);

        //! rebuild "By"
        if(wasMeshSimplificationOn==1)
        {
            data.setValue(meshSimplificationBy);
            Property prop_meshSimplificationBy("By",data,Property::PropertyGroup_Advanced);
            myCurNode->addProperty(prop_meshSimplificationBy,4);

            if(meshSimplificationBy==0)
            {
                //! rebuild mesh simplification by "Surface mesh size"
                data.setValue(level);
                Property prop_level("Level",data,Property::PropertyGroup_Advanced);
                myCurNode->addProperty(prop_level,5);
            }
            else
            {
                //! rebuild mesh simplification by "Pair distance"
                data.setValue(pairDistance);
                Property prop_pairDistance("Pair distance",data,Property::PropertyGroup_Advanced);
                myCurNode->addProperty(prop_pairDistance,5);
            }
        }
    }
        break;
    }
}

//! --------------------------------------------
//! function: handleMeshSimplificationChanged()
//! detail:
//! --------------------------------------------
void DetailViewer::handleMeshSimplificationChanged()
{
    cout<<"DetailViewer::handleMeshSimplificationChanged()->____function called____"<<endl;

    static int wasPairDistanceActive;
    static int wasLevelActive;
    static int level;
    static double pairDistance;

    QVariant data;
    int meshDefeaturingFlag = myCurNode->getPropertyValue<int>("Mesh simplification");
    switch(meshDefeaturingFlag)
    {
    case 0:
    {
        //! Off
        if(myCurNode->getPropertyItem("Pair distance")==Q_NULLPTR)
        {
            wasPairDistanceActive = 0;
        }
        else
        {
            wasPairDistanceActive = 1;
            pairDistance = myCurNode->getPropertyValue<double>("Pair distance");
        }
        if(myCurNode->getPropertyItem("Level")==Q_NULLPTR)
        {
            wasLevelActive = 0;
        }
        else
        {
            wasLevelActive = 1;
            level = myCurNode->getPropertyValue<int>("Level");
        }
        myCurNode->removeProperty("By");
        myCurNode->removeProperty("Pair distance");
        myCurNode->removeProperty("Level");
    }
        break;

    case 1:
    {
        //! On
        if(wasPairDistanceActive==1)
        {
            data.setValue(1);
            Property prop_meshDefeaturingBy("By",data,Property::PropertyGroup_Advanced);
            //myCurNode->addProperty(prop_meshDefeaturingBy,2);
            //myCurNode->addProperty(prop_meshDefeaturingBy,3);
            myCurNode->addProperty(prop_meshDefeaturingBy,4);

            data.setValue(pairDistance);
            Property prop_pairDistance("Pair distance",data,Property::PropertyGroup_Advanced);
            //myCurNode->addProperty(prop_pairDistance,3);
            //myCurNode->addProperty(prop_pairDistance,4);
            myCurNode->addProperty(prop_pairDistance,5);
        }
        if(wasLevelActive==1)
        {
            data.setValue(0);
            Property prop_meshDefeaturingBy("By",data,Property::PropertyGroup_Advanced);
            //myCurNode->addProperty(prop_meshDefeaturingBy,2);
            //myCurNode->addProperty(prop_meshDefeaturingBy,3);
            myCurNode->addProperty(prop_meshDefeaturingBy,4);

            data.setValue(level);
            Property prop_level("Level",data,Property::PropertyGroup_Advanced);
            //myCurNode->addProperty(prop_surfaceMeshSize,3);
            //myCurNode->addProperty(prop_surfaceMeshSize,4);
            myCurNode->addProperty(prop_level,5);
        }
    }
        break;
    }
}


//! ------------------------------------
//! function: handleMeshMethodChanged()
//! details:
//! ------------------------------------
void DetailViewer::handleMeshMethodChanged()
{
    cout<<"DetailViewer::handleMeshMethodChanged()->____function called____"<<endl;

    static Property::meshMethod method;
    static int meshDefeaturing;
    static int meshSimplification;
    static int meshHealing;
    static int by;
    static int level;
    static double pairDistance;
    static Property::meshEngine2D triangleSurfaceMesher;

    if(myCurNode->getPropertyItem("Mesh defeaturing")!=Q_NULLPTR) meshDefeaturing = myCurNode->getPropertyValue<int>("Mesh defeaturing");
    if(myCurNode->getPropertyItem("Mesh healing")!=Q_NULLPTR) meshHealing = myCurNode->getPropertyValue<int>("Mesh healing");
    if(myCurNode->getPropertyItem("Mesh simplification")!=Q_NULLPTR) meshSimplification = myCurNode->getPropertyValue<int>("Mesh simplification");
    if(myCurNode->getPropertyItem("Surface mesher")!=Q_NULLPTR) triangleSurfaceMesher = myCurNode->getPropertyValue<Property::meshEngine2D>("Surface mesher");
    if(myCurNode->getPropertyItem("By")!=Q_NULLPTR)
    {
        by = myCurNode->getPropertyValue<int>("By");
        if(by == 1) pairDistance = myCurNode->getPropertyValue<double>("Pair distance");
        else level = myCurNode->getPropertyValue<int>("Level");
    }

    Property::meshMethod mm = myCurNode->getPropertyValue<Property::meshMethod>("Method");
    switch(mm)
    {
    case Property::meshMethod_EMeshNetgen:
    case Property::meshMethod_EMeshTetgen:
    case Property::meshMethod_NetgenNetgen:
    case Property::meshMethod_NetgenTetgen:
    {
        method = mm;
        QList<QString> properties;
        properties<<"Mesh defeaturing"<<
                    "Surface mesher"<<
                    "Mesh healing"<<
                    "Mesh simplification"<<
                    "By"<<
                    "Pair distance"<<
                    "Level";
        for(int i=0; i<properties.length(); i++) myCurNode->removeProperty(properties.at(i));
    }
        break;

    case Property::meshMethod_automatic:
    {
        QVariant data;

        data.setValue(triangleSurfaceMesher);
        Property prop_triangleSurfaceMesher("Surface mesher",data,Property::PropertyGroup_Advanced);
        myCurNode->addProperty(prop_triangleSurfaceMesher,0);

        data.setValue(meshDefeaturing);
        Property prop_meshDefeaturing("Mesh defeaturing",data,Property::PropertyGroup_Advanced);
        myCurNode->addProperty(prop_meshDefeaturing,1);

        //! ---------------------------------------------------------------------------------
        //! add "By" "Pair distance" or "Surface mesh size" only if "Mesh simplification" is on
        //! ---------------------------------------------------------------------------------
        if(meshDefeaturing==1)
        {
            data.setValue(meshHealing);
            Property prop_meshHealing("Mesh healing",data,Property::PropertyGroup_Advanced);
            myCurNode->addProperty(prop_meshHealing,2);

            data.setValue(meshSimplification);
            Property prop_meshSimplification("Mesh simplification",data,Property::PropertyGroup_Advanced);
            myCurNode->addProperty(prop_meshSimplification,3);

            data.setValue(by);
            Property prop_by("By",data,Property::PropertyGroup_Advanced);
            myCurNode->addProperty(prop_by,4);

            if(by==1)
            {
                data.setValue(pairDistance);
                Property prop_pairDistance("Pair distance",data,Property::PropertyGroup_Advanced);
                myCurNode->addProperty(prop_pairDistance,5);
            }
            else
            {
                data.setValue(level);
                Property prop_surfaceMeshSize("Level",data,Property::PropertyGroup_Advanced);
                myCurNode->addProperty(prop_surfaceMeshSize,5);
            }
        }
    }
        break;
    }
}

//! ----------------------------------------
//! function: handlePrismaticLayerOptions()
//! details:
//! ----------------------------------------
void DetailViewer::handlePrismaticLayerOptions()
{
    cout<<"DetailViewer::handlePrismaticLayerOptions()->____function called____"<<endl;
    QVariant data;
    static double firstLayerHeight_old;
    static double totalThickness_old;

    int options = myCurNode->getPropertyValue<int>("Options");
    if(options==0) cout<<"____setting up \"First layer thickness\" options____"<<endl;
    else if(options==1) cout<<"____setting up \"Total thickness\" options____"<<endl;
    switch(options)
    {
    //! ------------------------------------------------------------
    //! generate prismatic layers by "First layer thickness" option
    //! ------------------------------------------------------------
    case 0:
    {
        if(myCurNode->getPropertyItem("Total thickness")!=Q_NULLPTR)
            totalThickness_old = myCurNode->getPropertyValue<double>("Total thickness");
        else totalThickness_old = 1.0;

        myCurNode->removeProperty("Total thickness");

        data.setValue(firstLayerHeight_old);
        Property prop_firstLayerHeight("First layer height",data,Property::PropertyGroup_Definition);
        myCurNode->addProperty(prop_firstLayerHeight,6);
    }
        break;
    //! ------------------------------------------------------
    //! generate prismatic layers by "Total thickness" option
    //! ------------------------------------------------------
    case 1:
    {
        if(myCurNode->getPropertyItem("First layer height")!=Q_NULLPTR)
            firstLayerHeight_old = myCurNode->getPropertyValue<double>("First layer height");
        else firstLayerHeight_old = 1.0;

        myCurNode->removeProperty("First layer height");

        data.setValue(totalThickness_old);
        Property prop_totalThickness("Total thickness",data,Property::PropertyGroup_Definition);
        myCurNode->addProperty(prop_totalThickness,6);
    }
        break;
    }
}

//! ---------------------------------------------
//! function: handleBoundaryScopingMethodChanged
//! details:
//! ---------------------------------------------
void DetailViewer::handleBoundaryScopingMethodChanged()
{
    cout<<"DetailViewer::handleBoundaryScopingMethodChanged()->____function called____"<<endl;
    QVariant data;
    Property::ScopingMethod scopingMethod = myCurNode->getPropertyValue<Property::ScopingMethod>("Boundary scoping method");
    switch(scopingMethod)
    {
    case Property::ScopingMethod_GeometrySelection:
    {
        std::vector<GeometryTag> tags = myCurNode->getPropertyValue<std::vector<GeometryTag>>("Boundary tags");
        //cout<<"____NUMBER OF FACES: "<<(int)tags.size()<<"____"<<endl;
        myCurNode->removeProperty("Boundary named selection");
        myCurNode->removeProperty("Boundary tags");
        data.setValue(tags);
        Property prop_namedSelection("Boundary",data,Property::PropertyGroup_Definition);
        myCurNode->addProperty(prop_namedSelection,2);
        Property prop_boundaryTags("Boundary tags",data,Property::PropertyGroup_Definition);
        myCurNode->addProperty(prop_boundaryTags,3);

        //! update tags?
    }
        break;

    case Property::ScopingMethod_NamedSelection:
    {
        myCurNode->removeProperty("Boundary");
        myCurNode->removeProperty("Boundary tags");

        SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
        QExtendedStandardItem *itemNSRoot = sm->getTreeItem(SimulationNodeClass::nodeType_namedSelection);
        QExtendedStandardItem *itemNSempty = static_cast<QExtendedStandardItem*>(itemNSRoot->child(0,0));
        void *p = (void*)itemNSempty;
        data.setValue(p);

        Property prop_namedSelection("Boundary named selection",data,Property::PropertyGroup_Definition);
        myCurNode->addProperty(prop_namedSelection,2);
        std::vector<GeometryTag> vecLoc;
        data.setValue(vecLoc);
        Property prop_boundaryTags("Boundary tags",data,Property::PropertyGroup_Definition);
        myCurNode->addProperty(prop_boundaryTags,3);

        //! update tags ... to do
    }
        break;
    }

}

//! -----------------------------
//! function: updateBoundaryTags
//! details:
//! -----------------------------
void DetailViewer::updateBoudaryTags()
{
    cout<<"DetailViewer::updateBoudaryTags()->____function called____"<<endl;

    QVariant data;
    Property::ScopingMethod scopingMethod =  myCurNode->getPropertyValue<Property::ScopingMethod>("Boundary scoping method");
    switch(scopingMethod)
    {
    case Property::ScopingMethod_GeometrySelection:
    {
        const std::vector<GeometryTag> &vecLoc = myCurNode->getPropertyValue<std::vector<GeometryTag>>("Boundary");
        data.setValue(vecLoc);
        Property prop_boundaryTags("Boundary tags",data,Property::PropertyGroup_Definition);
        myCurNode->replaceProperty("Boundary tags",prop_boundaryTags);
    }
        break;

    case Property::ScopingMethod_NamedSelection:
    {
        void *p = myCurNode->getPropertyValue<void*>("Boundary named selection");
        QStandardItem *itemNS = static_cast<QStandardItem*>(p);
        SimulationNodeClass *nodeNS = itemNS->data(Qt::UserRole).value<SimulationNodeClass*>();
        const std::vector<GeometryTag> &vecLoc = nodeNS->getPropertyValue<std::vector<GeometryTag>>("Tags");
        data.setValue(vecLoc);
        Property prop_boundaryTags("Boundary tags",data,Property::PropertyGroup_Definition);
        myCurNode->replaceProperty("Boundary tags",prop_boundaryTags);
    }
        break;
    }
}

//! ----------------------------------------------------
//! function: handleBoltCSChanged
//! details:  for the moment it changes the position of
//!           the double marker
//! ----------------------------------------------------
void DetailViewer::handleBoltCSChanged()
{
    //cout<<"DetailViewer::handleBoltCSChanged()->____function called____"<<endl;
    SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
    bool isDone = markerBuilder::addMarker(myCurNode,sm->getDataBase());
    if(isDone==true)
    {
        occHandle(AIS_Shape) theMarker = myCurNode->getPropertyValue<AIS_DoubleArrowMarker_handle_reg>("Graphic object");
        emit requestHideAllMarkers();
        theMarker->SetZLayer(Graphic3d_ZLayerId_TopOSD);
        myCTX->Display(theMarker,1,-1,true,false);
    }

    //! --------------
    //! block signals
    //! --------------
    myCurNode->getModel()->blockSignals(true);

    //! ---------------------------
    //! change the reference point
    //! ---------------------------
    void *cs = myCurNode->getPropertyValue<void*>("Coordinate system");
    QStandardItem *itemCS = static_cast<QStandardItem*>(cs);
    SimulationNodeClass *nodeCS = itemCS->data(Qt::UserRole).value<SimulationNodeClass*>();

    QVector<double> origin;
    if(nodeCS->getType()==SimulationNodeClass::nodeType_coordinateSystem_global)
    {
        for(int i=0; i<3; i++) origin.push_back(0);
    }
    else
    {
        QVector<double> baseOrigin = nodeCS->getPropertyValue<QVector<double>>("Base origin");
        for(int i=0; i<3; i++) origin.push_back(baseOrigin[i]);
    }
    //cout<<"____bolt: new reference point("<<origin[0]<<", "<<origin[1]<<", "<<origin[2]<<")____"<<endl;
    QVariant data;
    data.setValue(origin);
    Property prop_refPoint("Reference point",data,Property::PropertyGroup_Hidden);
    myCurNode->replaceProperty("Reference point",prop_refPoint);

    //! ----------------
    //! unblock signals
    //! ----------------
    myCurNode->getModel()->blockSignals(false);
}

//! ---------------------------------------------------------------
//! function: handleAccelerationChanged
//! details:  for the moment it changes the position of the marker
//! ---------------------------------------------------------------
void DetailViewer::handleAccelerationChanged()
{
    cout<<"DetailViewer::handleAccelerationChanged()->____function called____"<<endl;
    SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
    emit requestHideAllMarkers(true);
    bool isDone = markerBuilder::addMarker(myCurNode,sm->getDataBase());
    if(isDone==true)
    {
        occHandle(AIS_Shape) theMarker = myCurNode->getPropertyValue<AIS_ArrowMarker_handle_reg>("Graphic object");
        theMarker->SetZLayer(Graphic3d_ZLayerId_TopOSD);
        myCTX->Display(theMarker,1,-1,true,false);
    }
}

//! ---------------------------------------------------------------
//! function: handleMomentChanged
//! details:  for the moment it changes the position of the marker
//! ---------------------------------------------------------------
void DetailViewer::handleMomentChanged()
{
    cout<<"DetailViewer::handleMomentChanged()->____function called____"<<endl;
    SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
    emit requestHideAllMarkers(true);
    bool isDone = markerBuilder::addMarker(myCurNode,sm->getDataBase());
    if(isDone==true)
    {
        occHandle(AIS_Shape) theMarker = myCurNode->getPropertyValue<AIS_CurvedArrowMarker_handle_reg>("Graphic object");
        theMarker->SetZLayer(Graphic3d_ZLayerId_TopOSD);
        myCTX->Display(theMarker,1,-1,true,false);
    }
}

//! --------------------------------
//! function: handleBRepFlagChanged
//! details:
//! --------------------------------
void DetailViewer::handleBRepFlagChanged()
{
    cout<<"DetailViewer::handleBRepFlagChanged()->____function called____"<<endl;

    static Property::meshEngine2D discretizer;
    if(myCurNode->getPropertyItem("Tessellator")!=Q_NULLPTR) discretizer = myCurNode->getPropertyValue<Property::meshEngine2D>("Tessellator");
    else discretizer = Property::meshEngine2D_OCC_STL;

    static double angularDeflection;
    if(myCurNode->getPropertyItem("Angular deflection")!=Q_NULLPTR) angularDeflection = myCurNode->getPropertyValue<double>("Angular deflection");
    else angularDeflection = 0.5;

    static double linearDeflection;
    if(myCurNode->getPropertyItem("Linear deflection")!=Q_NULLPTR) linearDeflection = myCurNode->getPropertyValue<double>("Linear deflection");
    else linearDeflection = 0.1;

    static double minFaceSize;
    if(myCurNode->getPropertyItem("Min face size")!=Q_NULLPTR) minFaceSize = myCurNode->getPropertyValue<double>("Min face size");
    else minFaceSize = 0.1;

    static double maxFaceSize;
    if(myCurNode->getPropertyItem("Max face size")!=Q_NULLPTR) maxFaceSize = myCurNode->getPropertyValue<double>("Max face size");
    else maxFaceSize = 100;

    static bool defeaturing;
    if(myCurNode->getPropertyItem("Defeaturing")!=Q_NULLPTR) defeaturing = myCurNode->getPropertyValue<bool>("Defeaturing");
    else defeaturing = false;

    static bool healing;
    if(myCurNode->getPropertyItem("Healing")!=Q_NULLPTR) healing = myCurNode->getPropertyValue<bool>("Healing");
    else healing = false;

    static bool simplification;
    if(myCurNode->getPropertyItem("Simplification")!=Q_NULLPTR) simplification = myCurNode->getPropertyValue<bool>("Simplification");
    else simplification = true;

    static int by;
    if(myCurNode->getPropertyItem("By")!=Q_NULLPTR) by = myCurNode->getPropertyValue<int>("By");
    else by = 0;

    static double level;
    if(myCurNode->getPropertyItem("Level")!=Q_NULLPTR) level = myCurNode->getPropertyValue<double>("Level");
    else level = 0.90;

    static bool featurePreserve = false;
    if(myCurNode->getPropertyItem("Preserve boundary condition edges")!=Q_NULLPTR) featurePreserve = myCurNode->getPropertyValue<bool>("Preserve boundary conditions edges");

    static bool projectMeshPoints = false;
    if(myCurNode->getPropertyItem("Project points on geometry")!=Q_NULLPTR) projectMeshPoints = myCurNode->getPropertyValue<bool>("Project points on geometry");

    static double pairDistance;
    if(myCurNode->getPropertyItem("Pair distance")!=Q_NULLPTR) pairDistance = myCurNode->getPropertyValue<double>("Pair distance");
    else pairDistance = 0.1;

    static Property::meshEngine2D surfaceMesher_BRepOff;
    if(myCurNode->getPropertyItem("Surface mesher")!=Q_NULLPTR) surfaceMesher_BRepOff = myCurNode->getPropertyValue<Property::meshEngine2D>("Surface mesher");
    else surfaceMesher_BRepOff = Property::meshEngine2D_Netgen_STL;

    static Property::meshEngine2D surfaceMesher_BRepOn;
    if(myCurNode->getPropertyItem("Surface mesher")!=Q_NULLPTR) surfaceMesher_BRepOn = myCurNode->getPropertyValue<Property::meshEngine2D>("Surface mesher");
    else surfaceMesher_BRepOn = Property::meshEngine2D_Netgen;

    static Property::meshEngine3D volumeMesher_BRepOff;
    if(myCurNode->getPropertyItem("Volume mesher")!=Q_NULLPTR) volumeMesher_BRepOff = myCurNode->getPropertyValue<Property::meshEngine3D>("Volume mesher");
    else volumeMesher_BRepOff = Property::meshEngine3D_Netgen_STL;

    static Property::meshEngine3D volumeMesher_BRepOn;
    if(myCurNode->getPropertyItem("Volume mesher")!=Q_NULLPTR) volumeMesher_BRepOn = myCurNode->getPropertyValue<Property::meshEngine3D>("Volume mesher");
    else volumeMesher_BRepOn = Property::meshEngine3D_Netgen;

    static bool runInMemory;
    if(myCurNode->getPropertyItem("Run in memory")!=Q_NULLPTR) runInMemory = myCurNode->getPropertyValue<bool>("Run in memory");
    else runInMemory = false;

    //! -------------------
    //! TetWild parameters
    //! -------------------
    static int envelopeSizing;
    static double relativeSize;
    static double size;
    static int idealLengthSizing;
    static double relativeLength;
    static double length;

    QVariant data;

    //! ------------------------------------------------------
    //! this property is always defined for the "Method" item
    //! ------------------------------------------------------
    bool useBRep = myCurNode->getPropertyValue<bool>("Patch conforming");
    if(useBRep==true)
    {
        //! ------------------
        //! meshing uses BRep
        //! ------------------

        //! hide the "Defeaturing branch" to do ...

        //! -----------------------------
        //! store and remove
        //! remove "Tessellator"
        //! remove "Angular deflection"
        //! remove "Linear deflection"
        //! -----------------------------
        if(myCurNode->getPropertyItem("Discretizer")!=Q_NULLPTR) discretizer = myCurNode->getPropertyValue<Property::meshEngine2D>("Discretizer");
        if(myCurNode->getPropertyItem("Angular deflection")!=Q_NULLPTR) angularDeflection = myCurNode->getPropertyValue<double>("Angular deflection");
        if(myCurNode->getPropertyItem("Linear deflection")!=Q_NULLPTR) linearDeflection = myCurNode->getPropertyValue<double>("Linear deflection");
        if(myCurNode->getPropertyItem("Min face size")!=Q_NULLPTR) minFaceSize = myCurNode->getPropertyValue<double>("Min face size");
        if(myCurNode->getPropertyItem("Max face size")!=Q_NULLPTR) maxFaceSize = myCurNode->getPropertyValue<double>("Max face size");

        myCurNode->removeProperty("Tessellator");
        myCurNode->removeProperty("Angular deflection");
        myCurNode->removeProperty("Linear deflection");
        myCurNode->removeProperty("Min face size");
        myCurNode->removeProperty("Max face size");

        //! ------------------------
        //! store and remove
        //! remove "Defeaturing"
        //! remove "Healing"
        //! remove "Simplification"
        //! remove "By"
        //! remove "Level"
        //! remove "Pair distance"
        //! remote "Preserve boundary conditions edges"
        //! ------------------------
        if(myCurNode->getPropertyItem("Defeaturing")!=Q_NULLPTR) defeaturing = myCurNode->getPropertyValue<bool>("Defeaturing");
        if(myCurNode->getPropertyItem("Healing")!=Q_NULLPTR) healing = myCurNode->getPropertyValue<bool>("Healing");
        if(myCurNode->getPropertyItem("Simplification")!=Q_NULLPTR) simplification = myCurNode->getPropertyValue<bool>("Simplification");
        if(myCurNode->getPropertyItem("By")!=Q_NULLPTR) by = myCurNode->getPropertyValue<int>("By");
        if(myCurNode->getPropertyItem("Level")!=Q_NULLPTR) level = myCurNode->getPropertyValue<double>("Level");
        if(myCurNode->getPropertyItem("Pair distance")!=Q_NULLPTR) pairDistance = myCurNode->getPropertyValue<double>("Pair distance");
        if(myCurNode->getPropertyItem("Preserve boundary conditions edges")!=Q_NULLPTR) featurePreserve = myCurNode->getPropertyValue<bool>("Preserve boundary conditions edges");
        if(myCurNode->getPropertyItem("Project points on geometry")!=Q_NULLPTR) projectMeshPoints = myCurNode->getPropertyValue<bool>("Project points on geometry");

        myCurNode->removeProperty("Defeaturing");
        myCurNode->removeProperty("Healing");
        myCurNode->removeProperty("Simplification");
        myCurNode->removeProperty("By");
        myCurNode->removeProperty("Level");
        myCurNode->removeProperty("Pair distance");
        myCurNode->removeProperty("Preserve boundary conditions edges");
        myCurNode->removeProperty("Project points on geometry");
        myCurNode->removeProperty("Run in memory");

        //! ----------------------------------------
        //! store and remove the TetWild parameters
        //! ----------------------------------------
        if(myCurNode->getPropertyItem("Envelope sizing")!=Q_NULLPTR) envelopeSizing = myCurNode->getPropertyValue<int>("Envelope sizing");
        if(myCurNode->getPropertyItem("Relative size")!=Q_NULLPTR) relativeSize = myCurNode->getPropertyValue<double>("Relative size");
        if(myCurNode->getPropertyItem("Size")!=Q_NULLPTR) size = myCurNode->getPropertyValue<double>("Size");
        if(myCurNode->getPropertyItem("Ideal length sizing")!=Q_NULLPTR) idealLengthSizing = myCurNode->getPropertyValue<double>("Ideal length sizing");
        if(myCurNode->getPropertyItem("Relative length")!=Q_NULLPTR) relativeLength = myCurNode->getPropertyValue<double>("Relative length");
        if(myCurNode->getPropertyItem("Length")!=Q_NULLPTR) length = myCurNode->getPropertyValue<double>("Length");

        myCurNode->removeProperty("Envelope sizing");
        myCurNode->removeProperty("Relative size");
        myCurNode->removeProperty("Size");
        myCurNode->removeProperty("Ideal length sizing");
        myCurNode->removeProperty("Relative length");
        myCurNode->removeProperty("Length");

        //! ------------------------
        //! update the mesh engines
        //! ------------------------
        //data.setValue(surfaceMesher_BRepOn);
        data.setValue(Property::meshEngine2D_Netgen);
        Property prop_surfaceMesher("Surface mesher",data,Property::PropertyGroup_Definition);
        if(myCurNode->getPropertyItem("Surface mesher")!=Q_NULLPTR) myCurNode->replaceProperty("Surface mesher",prop_surfaceMesher);
        else myCurNode->addProperty(prop_surfaceMesher,1);

        //data.setValue(volumeMesher_BRepOff);
        data.setValue(Property::meshEngine3D_Netgen);
        Property prop_volumeMesher("Volume mesher",data,Property::PropertyGroup_Definition);
        myCurNode->replaceProperty("Volume mesher",prop_volumeMesher);
    }
    else
    {
        //! ----------------------------------------------------
        //! meshing process uses a surface tessellation as seed
        //! ----------------------------------------------------
        data.setValue(discretizer);
        Property prop_surfaceDiscretizer("Tessellator",data,Property::PropertyGroup_Method);
        myCurNode->addProperty(prop_surfaceDiscretizer,1);

        data.setValue(angularDeflection);
        Property prop_angularDeflection("Angular deflection",data,Property::PropertyGroup_Method);
        myCurNode->addProperty(prop_angularDeflection,2);

        data.setValue(linearDeflection);
        Property prop_linearDeflection("Linear deflection",data,Property::PropertyGroup_Method);
        myCurNode->addProperty(prop_linearDeflection,3);

        //! ---------------------------------------------------
        //! "Min face size" and "Max face size" apply only for
        //! the Express mesh tessellator
        //! ---------------------------------------------------
        if(myCurNode->getPropertyValue<Property::meshEngine2D>("Tessellator")==Property::meshEngine2D_OCC_ExpressMesh)
        {
            data.setValue(minFaceSize);
            Property prop_minFaceSize("Min face size",data,Property::PropertyGroup_Method);
            myCurNode->addProperty(prop_minFaceSize,4);

            data.setValue(maxFaceSize);
            Property prop_maxFaceSize("Max face size",data,Property::PropertyGroup_Method);
            myCurNode->addProperty(prop_maxFaceSize,5);
        }
        else
        {
            myCurNode->removeProperty("Min face size");
            myCurNode->removeProperty("Max face size");
        }

        data.setValue(defeaturing);
        Property prop_defeaturing("Defeaturing",data,Property::PropertyGroup_Defeaturing);
        myCurNode->addProperty(prop_defeaturing,0);

        if(defeaturing==true)
        {
            data.setValue(healing);
            Property prop_healing("Healing",data,Property::PropertyGroup_Defeaturing);
            myCurNode->addProperty(prop_healing,1);
            data.setValue(simplification);
            Property prop_simplification("Simplification",data,Property::PropertyGroup_Defeaturing);
            myCurNode->addProperty(prop_simplification,2);
            if(simplification==true)
            {
                data.setValue(by);
                Property prop_by("By",data,Property::PropertyGroup_Defeaturing);
                myCurNode->addProperty(prop_by,3);

                if(by==0)
                {
                    data.setValue(level);
                    Property prop_level("Level",data,Property::PropertyGroup_Defeaturing);
                    myCurNode->addProperty(prop_level,4);
                }
                else
                {
                    data.setValue(pairDistance);
                    Property prop_pairDistance("Pair distance",data,Property::PropertyGroup_Defeaturing);
                    myCurNode->addProperty(prop_pairDistance,4);
                }

                data.setValue(featurePreserve);
                myCurNode->addProperty(Property("Preserve boundary conditions edges",data,Property::PropertyGroup_Defeaturing),5);

                data.setValue(projectMeshPoints);
                myCurNode->addProperty(Property("Project points on geometry",data,Property::PropertyGroup_Defeaturing),6);
            }
        }

        //! ------------------------------------------------------
        //! since mesh based defeaturing is enabled the available
        //! surface mesher is "Netgen_STL" or "Tetgen BR"
        //! ------------------------------------------------------
        //! store
        surfaceMesher_BRepOn = myCurNode->getPropertyValue<Property::meshEngine2D>("Surface mesher");

        //! replace
        data.setValue(Property::meshEngine2D_Netgen_STL);
        //data.setValue(surfaceMesher_BRepOff);
        Property prop_surfaceMesher("Surface mesher",data,Property::PropertyGroup_Definition);
        myCurNode->removeProperty("Surface mesher");
        myCurNode->addProperty(prop_surfaceMesher,1);

        //! need to update the label
        //! to do ...

        //! store
        volumeMesher_BRepOn = myCurNode->getPropertyValue<Property::meshEngine3D>("Volume mesher");

        //! replace
        //data.setValue(volumeMesher_BRepOff);
        data.setValue(Property::meshEngine3D_Netgen_STL);
        Property prop_volumeMesher("Volume mesher",data,Property::PropertyGroup_Definition);
        myCurNode->removeProperty("Volume mesher");
        myCurNode->addProperty(prop_volumeMesher,2);

        //! need to update the label
        //! to do ...

        //! -----------------------
        //! handle "Run in memory"
        //! -----------------------
        switch(myCurNode->getPropertyValue<Property::meshEngine3D>("Volume mesher"))
        {
        case Property::meshEngine3D_Tetgen:
        case Property::meshEngine3D_Tetgen_BR:
        case Property::meshEngine3D_TetWild:
        {
            if(myCurNode->getPropertyItem("Run in memory")==Q_NULLPTR)
            {
                data.setValue(runInMemory);
                Property prop_runInMemory("Run in memory",data,Property::PropertyGroup_Definition);
                myCurNode->addProperty(prop_runInMemory);
            }
        }
            break;
        case Property::meshEngine3D_Netgen_STL:
        case Property::meshEngine3D_Netgen:
        {
            //! store and replace
            if(myCurNode->getPropertyItem("Run in memory")!=Q_NULLPTR) runInMemory = myCurNode->getPropertyValue<bool>("Run in memory");
            myCurNode->removeProperty("Run in memory");
        }
            break;
        }
    }    
}

//! -----------------------------------
//! function: handleTessellatorChanged
//! details:
//! -----------------------------------
void DetailViewer::handleTessellatorChanged()
{
    static double minFaceSize, maxFaceSize;

    if(myCurNode->getPropertyItem("Min face size")!=Q_NULLPTR) minFaceSize = myCurNode->getPropertyValue<double>("Min face size");
    else minFaceSize = 0.01;

    if(myCurNode->getPropertyItem("Max face size")!=Q_NULLPTR) maxFaceSize = myCurNode->getPropertyValue<double>("Max face size");
    else maxFaceSize = 1000;

    Property::meshEngine2D Tessellator = myCurNode->getPropertyValue<Property::meshEngine2D>("Tessellator");
    switch(Tessellator)
    {
    case Property::meshEngine2D_OCC_ExpressMesh:
    {
        QVariant data;
        data.setValue(minFaceSize);
        Property prop_minFaceSize("Min face size",data,Property::PropertyGroup_Method);
        myCurNode->addProperty(prop_minFaceSize,4);
        data.setValue(maxFaceSize);
        Property prop_maxFaceSize("Max face size",data,Property::PropertyGroup_Method);
        myCurNode->addProperty(prop_maxFaceSize,5);
    }
        break;
    case Property::meshEngine2D_OCC_STL:
    {
        minFaceSize = myCurNode->getPropertyValue<double>("Min face size");
        maxFaceSize = myCurNode->getPropertyValue<double>("Max face size");

        myCurNode->removeProperty("Min face size");
        myCurNode->removeProperty("Max face size");
    }
        break;
    }
 }

//! ----------------------------------------
//! function: handleDefeaturingFlagChanged
//! details:
//! ----------------------------------------
void DetailViewer::handleDefeaturingFlagChanged()
{
    cout<<"DetailViewer::handleDefeaturingFlagChanged()->____function called____"<<endl;

    QVariant data;

    static bool healing;
    if(myCurNode->getPropertyItem("Healing")==Q_NULLPTR) healing = true;
    else healing = myCurNode->getPropertyValue<bool>("Healing");

    static bool simplification;
    if(myCurNode->getPropertyItem("Simplification")==Q_NULLPTR) simplification = true;
    else simplification = myCurNode->getPropertyValue<bool>("Simplification");

    static int by;
    if(myCurNode->getPropertyItem("By")==Q_NULLPTR) by = 0;
    else by = myCurNode->getPropertyValue<int>("By");

    static double level;
    if(myCurNode->getPropertyItem("Level")==Q_NULLPTR) level = 0.95;
    else level = myCurNode->getPropertyValue<double>("Level");

    static double pairDistance;
    if(myCurNode->getPropertyItem("Pair distance")==Q_NULLPTR) pairDistance = 0.01;
    else pairDistance = myCurNode->getPropertyValue<double>("Pair distance");

    static bool featurePreserve;
    if(myCurNode->getPropertyItem("Preserve boundary conditions edges")==Q_NULLPTR) featurePreserve = false;
    else featurePreserve = myCurNode->getPropertyValue<bool>("Preserve boundary conditions edges");

    static bool projectPoints;
    if(myCurNode->getPropertyItem("Project points on geometry")==Q_NULLPTR) projectPoints = false;
    else projectPoints = myCurNode->getPropertyValue<bool>("Project points on geometry");

    bool defeaturing = myCurNode->getPropertyValue<bool>("Defeaturing");
    if(defeaturing)
    {
        data.setValue(healing);
        Property prop_healing("Healing",data,Property::PropertyGroup_Defeaturing);
        myCurNode->addProperty(prop_healing,1);

        data.setValue(simplification);
        Property prop_simplification("Simplification",data,Property::PropertyGroup_Defeaturing);
        myCurNode->addProperty(prop_simplification,2);

        if(simplification==true)
        {
            data.setValue(by);
            Property prop_by("By",data,Property::PropertyGroup_Defeaturing);
            myCurNode->addProperty(prop_by,3);

            if(by==0)
            {
                data.setValue(level);
                Property prop_level("Level",data,Property::PropertyGroup_Defeaturing);
                myCurNode->addProperty(prop_level,4);
            }
            else
            {
                data.setValue(pairDistance);
                Property prop_distance("Level",data,Property::PropertyGroup_Defeaturing);
                myCurNode->addProperty(prop_distance,4);
            }

            data.setValue(featurePreserve);
            myCurNode->addProperty(Property("Preserve boundary conditions edges",data,Property::PropertyGroup_Defeaturing),5);

            data.setValue(projectPoints);
            myCurNode->addProperty(Property("Project points on geometry",data,Property::PropertyGroup_Defeaturing),6);
        }
    }
    else
    {
        //! -----------------
        //! store and remove
        //! -----------------
        if(myCurNode->getPropertyItem("Healing")!=Q_NULLPTR) healing = myCurNode->getPropertyValue<bool>("Healing");
        myCurNode->removeProperty("Healing");
        if(myCurNode->getPropertyItem("Simplification")!=Q_NULLPTR) simplification = myCurNode->getPropertyValue<bool>("Simplification");
        myCurNode->removeProperty("Simplification");
        if(myCurNode->getPropertyItem("By")!=Q_NULLPTR) by = myCurNode->getPropertyValue<int>("By");
        myCurNode->removeProperty("By");
        if(myCurNode->getPropertyItem("Level")!=Q_NULLPTR) level = myCurNode->getPropertyValue<double>("Level");
        myCurNode->removeProperty("Level");
        if(myCurNode->getPropertyItem("Pair distance")!=Q_NULLPTR) pairDistance = myCurNode->getPropertyValue<double>("Pair distance");
        myCurNode->removeProperty("Pair distance");
        if(myCurNode->getPropertyItem("Preserve boundary conditions edges")!=Q_NULLPTR) featurePreserve = myCurNode->getPropertyValue<bool>("Preserve boundary conditions edges");
        myCurNode->removeProperty("Preserve boundary conditions edges");
        if(myCurNode->getPropertyItem("Project points on geometry")!=Q_NULLPTR) projectPoints = myCurNode->getPropertyValue<bool>("Project points on geometry");
        myCurNode->removeProperty("Project points on geometry");
    }
}

//! ------------------------------------------
//! function: handleSimplificationFlagChanged
//! details:
//! ------------------------------------------
void DetailViewer::handleSimplificationFlagChanged()
{
    QVariant data;

    static int by;
    if(myCurNode->getPropertyItem("By")==Q_NULLPTR) by = 0;
    else by = myCurNode->getPropertyValue<int>("By");

    static double level;
    if(myCurNode->getPropertyItem("Level")==Q_NULLPTR) level = 0.95;
    else level = myCurNode->getPropertyValue<double>("Level");

    static double pairDistance;
    if(myCurNode->getPropertyItem("Pair distance")==Q_NULLPTR) pairDistance = 0.01;
    else pairDistance = myCurNode->getPropertyValue<double>("Pair distance");

    static bool featurePreserve = false;
    if(myCurNode->getPropertyItem("Preserve boundary conditions edges")==Q_NULLPTR) featurePreserve = false;
    else featurePreserve = myCurNode->getPropertyValue<bool>("Preserve boundary conditions edges");

    static bool projectPoints = false;
    if(myCurNode->getPropertyItem("Project points on geometry")==Q_NULLPTR) projectPoints = false;
    else projectPoints = myCurNode->getPropertyValue<bool>("Project points on geometry");

    bool simplification = myCurNode->getPropertyValue<bool>("Simplification");
    if(simplification==true)
    {
        data.setValue(by);
        Property prop_by("By",data,Property::PropertyGroup_Defeaturing);
        myCurNode->addProperty(prop_by,3);
        if(by==0)
        {
            data.setValue(level);
            Property prop_level("Level",data,Property::PropertyGroup_Defeaturing);
            myCurNode->addProperty(prop_level,4);
        }
        else
        {
            data.setValue(pairDistance);
            Property prop_pairDistance("Pair distance",data,Property::PropertyGroup_Defeaturing);
            myCurNode->addProperty(prop_pairDistance,4);
        }
        data.setValue(featurePreserve);
        myCurNode->addProperty(Property("Preserve boundary conditions edges",data,Property::PropertyGroup_Defeaturing),5);

        data.setValue(projectPoints);
        myCurNode->addProperty(Property("Project points on geometry",data,Property::PropertyGroup_Defeaturing),6);
    }
    else
    {
        //! -----------------
        //! store and remove
        //! -----------------
        if(myCurNode->getPropertyItem("By")!=Q_NULLPTR) by = myCurNode->getPropertyValue<int>("By");
        if(myCurNode->getPropertyItem("Level")!=Q_NULLPTR) level = myCurNode->getPropertyValue<double>("Level");
        if(myCurNode->getPropertyItem("Pair distance")!=Q_NULLPTR) pairDistance = myCurNode->getPropertyValue<double>("Pair distance");
        if(myCurNode->getPropertyItem("Preserve boundary conditions edges")!=Q_NULLPTR) featurePreserve = myCurNode->getPropertyValue<bool>("Preserve boundary conditions edges");
        if(myCurNode->getPropertyItem("Project points on geometry")!=Q_NULLPTR) projectPoints = myCurNode->getPropertyValue<bool>("Project points on geometry");

        myCurNode->removeProperty("By");
        myCurNode->removeProperty("Level");
        myCurNode->removeProperty("Pair distance");
        myCurNode->removeProperty("Preserve boundary conditions edges");
        myCurNode->removeProperty("Project points on geometry");
    }
}

//! ----------------------------------------------
//! function: handleMeshSimplificationByChanged()
//! details:
//! ----------------------------------------------
void DetailViewer::handleMeshSimplificationByChanged()
{
    cout<<"DetailViewer::handleMeshSimplificationByChanged()->____function called____"<<endl;

    static double level;
    static double pairDistance;

    if(myCurNode->getPropertyItem("Level")!=Q_NULLPTR) level = myCurNode->getPropertyValue<double>("Level");
    else level = 0.95;

    if(myCurNode->getPropertyItem("Pair distance")!=Q_NULLPTR) pairDistance = myCurNode->getPropertyValue<double>("Pair distance");
    else pairDistance = 0.01;

    QVariant data;
    int by = myCurNode->getPropertyValue<int>("By");
    switch(by)
    {
    case 0: //! activate "Level"
    {
        pairDistance = myCurNode->getPropertyValue<double>("Pair distance");
        myCurNode->removeProperty("Pair distance");
        data.setValue(level);
        Property prop_level("Level",data,Property::PropertyGroup_Defeaturing);
        myCurNode->addProperty(prop_level,4);
    }
        break;

    case 1: //! activate "Pair distance"
    {
        level = myCurNode->getPropertyValue<double>("Level");
        myCurNode->removeProperty("Level");
        data.setValue(pairDistance);
        Property prop_pairDistance("Pair distance",data,Property::PropertyGroup_Defeaturing);
        myCurNode->addProperty(prop_pairDistance,4);
    }
        break;
    }
}

//! ------------------------------------
//! function: handleVolumeMesherChanged
//! details:
//! ------------------------------------
void DetailViewer::handleVolumeMesherChanged()
{
    //! ---------------------------------
    //! store the defeaturing properties
    //! ---------------------------------
    static bool isDefeaturingON = false;
    if(myCurNode->getPropertyItem("Defeaturing")!=Q_NULLPTR) isDefeaturingON = myCurNode->getPropertyValue<bool>("Defeaturing");

    static bool isHealingOn = false;
    static bool isSimplificationOn = false;
    static int by = 0;
    static double level = 0.90;
    static double pairDistance = 0.1;
    if(isDefeaturingON)
    {
        isHealingOn = myCurNode->getPropertyValue<bool>("Healing");
        isSimplificationOn = myCurNode->getPropertyValue<bool>("Simplification");
        if(isSimplificationOn)
        {
            by = myCurNode->getPropertyValue<int>("By");
            if(by==0)
            {
                level = myCurNode->getPropertyValue<double>("Level");
            }
            else
            {
                pairDistance = myCurNode->getPropertyValue<double>("Pair distance");
            }
        }
    }
    //! ------------------------------------------------------
    //! end store the defeaturing options - could be removed?
    //! ------------------------------------------------------

    QVariant data;

    static Property::meshEngine2D surfaceMesher;
    if(myCurNode->getPropertyItem("Surface mesher")!=Q_NULLPTR) surfaceMesher = myCurNode->getPropertyValue<Property::meshEngine2D>("Surface mesher");
    else surfaceMesher = Property::meshEngine2D_Netgen_STL;

    Property::meshEngine3D volumeMesher = myCurNode->getPropertyValue<Property::meshEngine3D>("Volume mesher");
    switch(volumeMesher)
    {
    case Property::meshEngine3D_Tetgen_BR:
    case Property::meshEngine3D_TetWild:
    {
        //! before removing the Surface mesh property check if it is defined
        if(myCurNode->getPropertyItem("Surface mesher")!=Q_NULLPTR)
        {
            //! store and remove
            //if(myCurNode->getPropertyValue<Property::meshEngine2D>("Surface mesher"))
                surfaceMesher = myCurNode->getPropertyValue<Property::meshEngine2D>("Surface mesher");
            myCurNode->removeProperty("Surface mesher");
        }
    }
        break;

    case Property::meshEngine3D_Netgen_STL:
    case Property::meshEngine3D_Tetgen:
    {
        if(myCurNode->getPropertyItem("Surface mesher")==Q_NULLPTR)
        {
        data.setValue(surfaceMesher);
        Property prop_surfaceMesher("Surface mesher",data,Property::PropertyGroup_Definition);
        myCurNode->addProperty(prop_surfaceMesher,1);
        }
        else
        {
            data.setValue(surfaceMesher);
            Property prop_surfaceMesher("Surface mesher",data,Property::PropertyGroup_Definition);
            myCurNode->replaceProperty("Surface mesher",prop_surfaceMesher);
        }
    }
        break;
    }

    //! -----------------------
    //! handle "Run in memory"
    //! -----------------------
    switch(volumeMesher)
    {
    case Property::meshEngine3D_Netgen_STL:
    case Property::meshEngine3D_Netgen:
    {
        myCurNode->removeProperty("Run in memory");
    }
        break;
    case Property::meshEngine3D_Tetgen_BR:
    case Property::meshEngine3D_Tetgen:
    case Property::meshEngine3D_TetWild:
    {
        if(myCurNode->getPropertyItem("Run in memory")==Q_NULLPTR)
        {
            data.setValue(false);
            Property prop_runInMemory("Run in memory",data,Property::PropertyGroup_Definition);
            myCurNode->addProperty(prop_runInMemory);
        }
    }
        break;
    }

    //! --------------------------------
    //! handle the "Defeaturing" branch
    //! --------------------------------
    switch(volumeMesher)
    {
    case Property::meshEngine3D_Tetgen:
    {
        bool useBRep = myCurNode->getPropertyValue<bool>("Patch conforming");
        if(useBRep==false)
        {
            //! -----------------------------
            //! add the "Defeaturing" branch
            //! -----------------------------
            if(myCurNode->getPropertyItem("Defeaturing")==Q_NULLPTR)
            {
                data.setValue(false);
                Property prop_defeaturing("Defeaturing",data,Property::PropertyGroup_Defeaturing);
                myCurNode->addProperty(prop_defeaturing,0);
            }
            //! ---------------------------------------
            //! remove the TewWild parameters settings
            //! ---------------------------------------
            myCurNode->removeProperty("Envelope sizing");
            myCurNode->removeProperty("Relative size");
            myCurNode->removeProperty("Size");
            myCurNode->removeProperty("Ideal length sizing");
            myCurNode->removeProperty("Relative length");
            myCurNode->removeProperty("Length");
        }
        else
        {
            //! ---------------------------------------
            //! remove the TewWild parameters settings
            //! (all, defined or not defined)
            //! ---------------------------------------
            myCurNode->removeProperty("Envelope sizing");
            myCurNode->removeProperty("Relative size");
            myCurNode->removeProperty("Size");
            myCurNode->removeProperty("Ideal length sizing");
            myCurNode->removeProperty("Relative length");
            myCurNode->removeProperty("Length");

            //! -------------------------------------
            //! remove the "Defeaturing" branch
            //! (all the properties, defined or not)
            //! -------------------------------------
            myCurNode->removeProperty("Defeaturing");
            myCurNode->removeProperty("Healing");
            myCurNode->removeProperty("Simplification");
            myCurNode->removeProperty("By");
            myCurNode->removeProperty("Level");
            myCurNode->removeProperty("Pair distance");
        }
    }
        break;

    case Property::meshEngine3D_Tetgen_BR:
    case Property::meshEngine3D_Netgen_STL:
    {
        //! -----------------------------
        //! add the "Defeaturing" branch
        //! -----------------------------
        if(myCurNode->getPropertyItem("Defeaturing")==Q_NULLPTR)
        {
            data.setValue(false);
            Property prop_defeaturing("Defeaturing",data,Property::PropertyGroup_Defeaturing);
            myCurNode->addProperty(prop_defeaturing,0);
        }
        //! ---------------------------------------
        //! remove the TewWild parameters settings
        //! ---------------------------------------
        myCurNode->removeProperty("Envelope sizing");
        myCurNode->removeProperty("Relative size");
        myCurNode->removeProperty("Size");
        myCurNode->removeProperty("Ideal length sizing");
        myCurNode->removeProperty("Relative length");
        myCurNode->removeProperty("Length");
    }
        break;

    case Property::meshEngine3D_TetWild:
    {
        //! -------------------------------------
        //! remove the "Defeaturing" branch
        //! (all the properties, defined or not)
        //! -------------------------------------
        /*
        myCurNode->removeProperty("Defeaturing");
        myCurNode->removeProperty("Healing");
        myCurNode->removeProperty("Simplification");
        myCurNode->removeProperty("By");
        myCurNode->removeProperty("Level");
        myCurNode->removeProperty("Pair distance");
        */
        //! hide the "Defeaturing" branch ... to do ...

        //! -------------------------------------------------
        //! TetWild parameters
        //!
        //! Envelope sizing: 0 => relative 1 => absolute
        //! Relative size (for option "0") => 0.0025
        //! Ideal length sizing: 0 => relative 1 => absolute
        //! Length (for option "0") => 0.05;
        //! -------------------------------------------------
        int envelopeSizing = 0;
        data.setValue(envelopeSizing);
        Property prop_envelopeSizingType("Envelope sizing", data, Property::PropertyGroup_TetWild);
        myCurNode->addProperty(prop_envelopeSizingType);
        double value = 0.001;
        data.setValue(value);
        Property prop_relativeSize("Relative size",data, Property::PropertyGroup_TetWild);
        myCurNode->addProperty(prop_relativeSize);
        int idealLengthSizing = 0;
        data.setValue(idealLengthSizing);
        Property prop_idealLength("Ideal length sizing",data, Property::PropertyGroup_TetWild);
        myCurNode->addProperty(prop_idealLength);
        double relativeLength = 0.05;
        data.setValue(relativeLength);
        Property prop_relativeLength("Relative length",data, Property::PropertyGroup_TetWild);
        myCurNode->addProperty(prop_relativeLength);
    }
        break;
    }
}



//! -------------------------------------------
//! function: handleGeometryHealingChanged
//! details:  tolerance = 0 => "Normal"
//!           tolerance = 1 => "Loose"
//!           tolerance = 2 => "User defined"
//! ------------------------------------------
void DetailViewer::handleGeometryHealingChanged()
{
    static int tolerance;
    if(myCurNode->getPropertyItem("Tolerance")!=Q_NULLPTR) tolerance = myCurNode->getPropertyValue<double>("Tolerance");
    else tolerance = 1;

    static double toleranceValue;
    if(myCurNode->getPropertyItem("Tolerance value")!=Q_NULLPTR) toleranceValue = myCurNode->getPropertyValue<double>("Tolerance value");
    else toleranceValue = 0.0001;

    QVariant data;
    bool isHealingActive = myCurNode->getPropertyValue<bool>("Geometry healing");
    if(isHealingActive)
    {
        data.setValue(tolerance);

        Property prop_tolerance("Tolerance",data,Property::PropertyGroup_Definition);

        if(myCurNode->getPropertyItem("Tolerance")==Q_NULLPTR) myCurNode->addProperty(prop_tolerance,5);
        else myCurNode->replaceProperty("Tolerance",prop_tolerance);

        //! --------------------------------------------------
        //! check what kind of "Tolerance" was defined before
        //! --------------------------------------------------
        switch(tolerance)
        {
        case 1: case 2: //! "Normal" or "Loose"
        {
            //! store "Tolerance value" and remove
            if(myCurNode->getPropertyItem("Tolerance value")!=Q_NULLPTR) toleranceValue = myCurNode->getPropertyValue<double>("Tolerance value");
            myCurNode->removeProperty("Tolerance value");
        }
            break;

        case 0: //! "User defined"
        {
            data.setValue(toleranceValue);
            Property prop_toleranceValue("Tolerance value",data,Property::PropertyGroup_Definition);
            myCurNode->addProperty(prop_toleranceValue,6);
        }
            break;
        }
    }
    else
    {
        //! store "Tolerance" and remove
        tolerance = myCurNode->getPropertyValue<double>("Tolerance");
        myCurNode->removeProperty("Tolerance");

        //! store "Tolerance value" if present, and remove
        if(myCurNode->getPropertyItem("Tolerance value")!=Q_NULLPTR) toleranceValue = myCurNode->getPropertyValue<double>("Tolerance value");
        myCurNode->removeProperty("Tolerance value");
    }
}

#ifdef COSTAMP_VERSION
void DetailViewer::handleRequestStartTimeStepBuilder()
{
    cout<<"_____REQUEST START TIME STEP BUILDER____"<<endl;
    emit requestStartTimeStepBuilder();
}
#endif

//! -----------------------------------
//! function: handleIdealLengthChanged
//! details:
//! -----------------------------------
void DetailViewer::handleIdealLengthChanged()
{
    static double absIdealLength, relIdealLength;

    myCurNode->getModel()->blockSignals(true);
    QVariant data;
    int idealLength = myCurNode->getPropertyValue<int>("Ideal length sizing");
    switch(idealLength)
    {
    case 0:
    {
        //! store and remove
        absIdealLength = myCurNode->getPropertyValue<double>("Absolute length");
        myCurNode->removeProperty("Absolute length");
        //! new property
        data.setValue(relIdealLength);
        Property prop_relIdealLength("Relative length",data,Property::PropertyGroup_TetWild);
        myCurNode->addProperty(prop_relIdealLength,3);
    }
        break;
    case 1:
    {
        //! store and remove
        relIdealLength = myCurNode->getPropertyValue<double>("Relative length");
        myCurNode->removeProperty("Relative length");
        //! new property
        data.setValue(absIdealLength);
        Property prop_absIdealLength("Absolute length",data,Property::PropertyGroup_TetWild);
        myCurNode->addProperty(prop_absIdealLength,3);
    }
        break;
    }
    myCurNode->getModel()->blockSignals(false);
}

//! --------------------------------------
//! function: handleEnvelopeSizingChanged
//! details:
//! --------------------------------------
void DetailViewer::handleEnvelopeSizingChanged()
{
    static double relSize,absSize;
    QVariant data;

    int envSizing = myCurNode->getPropertyValue<int>("Envelope sizing");
    switch(envSizing)
    {
    case 0:
    {
        //! store and remove
        absSize = myCurNode->getPropertyValue<double>("Absolute size");
        myCurNode->removeProperty("Absolute size");
        //! new property
        data.setValue(relSize);
        Property prop_relSize("Relative size",data,Property::PropertyGroup_TetWild);
        myCurNode->addProperty(prop_relSize,1);
    }
        break;
    case 1:
    {
        //! store and remove
        relSize = myCurNode->getPropertyValue<double>("Relative size");
        myCurNode->removeProperty("Relative size");
        //! new property
        data.setValue(absSize);
        Property prop_absSize("Absolute size",data,Property::PropertyGroup_TetWild);
        myCurNode->addProperty(prop_absSize,1);
    }
        break;
    }
}

//! ----------------------------------------------
//! function: handleInterpolationAlgorithmChanged
//! details:
//! ----------------------------------------------
void DetailViewer::handleInterpolationAlgorithmChanged()
{
    cout<<"DetailViewer::handleInterpolationAlgorithmChanged()->____function called____"<<endl;
    static double pinballVal;

    if(myCurNode->getPropertyItem("Pinball")!=Q_NULLPTR) pinballVal = myCurNode->getPropertyValue<double>("Pinball");
    else pinballVal = 5.0;
    int algo = myCurNode->getPropertyValue<int>("Algorithm");
    QVariant data;
    switch(algo)
    {
    //case 0:
    case 0: case 1:
    {
        if(myCurNode->getPropertyItem("Pinball")==Q_NULLPTR)
        {
            data.setValue(pinballVal);
            Property prop_pinball("Pinball",data,Property::PropertyGroup_Advanced);
            myCurNode->addProperty(prop_pinball,1);
        }
    }
        break;

    //case 1:
    case 2:
    {
        //! -----------------
        //! store and remove
        //! -----------------
        if(myCurNode->getPropertyItem("Pinball")!=Q_NULLPTR) pinballVal = myCurNode->getPropertyValue<double>("Pinball");
        myCurNode->removeProperty("Pinball");
    }
        break;
    }
}

//! ----------------------------------------------
//! function: removeModelBranch
//! details:  remove a branch from the node model
//! ----------------------------------------------
#include "nodefactory.h"
void DetailViewer::removeModelBranch(QStandardItem *theSeparator)
{
    cout<<"DetailViewer::removeModelBranch()->____function called____"<<endl;
    QStandardItemModel *theNodeModel = static_cast<QStandardItemModel*>(myCurNode->getModel());
    QModelIndex indexOfTheSeparator = theNodeModel->indexFromItem(theSeparator);
    if(!indexOfTheSeparator.isValid()) return;
    theNodeModel->removeRow(indexOfTheSeparator.row());
}

//! ------------------------------------------------
//! function: handleModelChangeScopingMethodChanged
//! details:
//! ------------------------------------------------
void DetailViewer::handleModelChangeScopingMethodChanged()
{
    cout<<"DetailViewer::handleModelChangeScopingMethodChanged()->____function called____"<<endl;

    QVariant data;
    int itemType = myCurNode->getPropertyValue<int>("Item type");
    switch(itemType)
    {
    case 0:
    {
        myCurNode->removeProperty("Contact pair");

        data.setValue(Property::ScopingMethod_GeometrySelection);
        Property prop_scopingMethod("Scoping method",data,Property::PropertyGroup_Scope);
        std::vector<GeometryTag> emptyTags;
        data.setValue(emptyTags);
        Property prop_scope("Geometry",data,Property::PropertyGroup_Scope);
        Property prop_tags("Tags",data,Property::PropertyGroup_Scope);

        myCurNode->addProperty(prop_scopingMethod);
        myCurNode->addProperty(prop_scope);
        myCurNode->addProperty(prop_tags);
    }
        break;

    case 1:
    {
        myCurNode->removeProperty("Scoping method");
        myCurNode->removeProperty("Geometry");
        myCurNode->removeProperty("Tags");

        SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
        QStandardItem *itemContactRoot = sm->getTreeItem(SimulationNodeClass::nodeType_connection);
        //cout<<"____number of children: "<<itemContactGroup->rowCount()<<"____"<<endl;
        QStandardItem *itemDummyContactPair = itemContactRoot->child(0,0);
        if(itemDummyContactPair!=Q_NULLPTR)
        {
            void *p = (void*)(itemDummyContactPair);
            data.setValue(p);
            Property prop_itemContactPair("Contact pair",data,Property::PropertyGroup_Scope);
            myCurNode->addProperty(prop_itemContactPair);
        }
    }
        break;
    }
}

//! ---------------------------------------------------
//! function: handleModelChangeActivationStatusChanged
//! details:
//! ---------------------------------------------------
void DetailViewer::handleModelChangeActivationStatusChanged()
{
    cout<<"DetailViewer::handleModelChangeActivationStatusChanged()->____function called____"<<endl;

    //! ----------------------------------------------------------------
    //! avoids ping-pong - reconnection at the end of this function [*]
    //! ----------------------------------------------------------------
    SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
    TableWidget *tableWidget = static_cast<TableWidget*>(tools::getWidgetByName("messagesAndLoadsWidget"));
    disconnect(tableWidget,SIGNAL(requestUpdateDetailViewer(QModelIndex,QModelIndex,QVector<int>)),
               this,SLOT(updateDetailViewerFromTabularData(QModelIndex,QModelIndex,QVector<int>)));

    //! --------------------------
    //! retrieve the tabular data
    //! --------------------------
    SimulationNodeClass *nodeAnalysisSettings = sm->getAnalysisSettingsNodeFromCurrentItem();

    CustomTableModel *tabData = nodeAnalysisSettings->getTabularDataModel();

    //! --------------------------------------
    //! retrieve the current step time number
    //! --------------------------------------
    int currentStepNumber = nodeAnalysisSettings->getPropertyValue<int>("Current step number");

    int row = currentStepNumber;
    int col = mainTreeTools::calculateStartColumn(sm->myTreeView);

    //! ----------------------
    //! generate the new data
    //! ----------------------
    QVariant data;
    Property::modelChangeActivationStatus as = myCurNode->getPropertyValue<Property::modelChangeActivationStatus>("Activation status");

    data.setValue(as);
    tabData->setDataRC(data,row,col,Qt::EditRole);

    //! -------------
    //! reconnection
    //! -------------
    connect(tableWidget,SIGNAL(requestUpdateDetailViewer(QModelIndex,QModelIndex,QVector<int>)),
            this,SLOT(updateDetailViewerFromTabularData(QModelIndex,QModelIndex,QVector<int>)));
}

//! ------------------------------------
//! function: handleTransparencyChanged
//! details:
//! ------------------------------------
void DetailViewer::handleTransparencyChanged()
{
    double aLevel = myCurNode->getPropertyValue<double>("Transparency");
    emit requestHandleTransparencyChanged(aLevel);
}

//! ---------------------------------------
//! function: handleSelectionMethodChanged
//! details:
//! ---------------------------------------
void DetailViewer::handleSelectionMethodChanged()//cesere
{
    this->connectToSimulationManager(false);
    //myCurNode->getModel()->blockSignals(true);

    QVariant data;
    int selectionMethod = myCurNode->getPropertyValue<int>("Selection method");
    switch(selectionMethod)
    {
    case 0:
    {
        data.setValue(QString(""));
        Property prop_elementList("Element list",data,Property::PropertyGroup_Definition);
        myCurNode->removeProperty("Mesh entities");
        myCurNode->addProperty(prop_elementList,1);
    }
        break;
    case 1:
    {
        data.setValue(std::vector<int>());
        myCurNode->removeProperty("Element list");
        Property prop_meshEntities("Mesh entities",data,Property::PropertyGroup_Definition);
        myCurNode->removeProperty("Element list");
        myCurNode->addProperty(prop_meshEntities,1);
    }
        break;
    }

    //myCurNode->getModel()->blockSignals(false);
    this->connectToSimulationManager(true);
}

//! ----------------------------------
//! function: addMeshElementSelection
//! details:
//! ----------------------------------
#include <meshelementbycoords.h>
void DetailViewer::handleMeshElementListChanged()
{
    cout<<"DetailViewer::addMeshElementSelection()->____function called____"<<endl;
    std::vector<int> elementList;
    QString str = myCurNode->getPropertyValue<QString>("Element list");
    std::stringstream ss(str.toStdString());
    while (ss.good())
    {
        std::string substr;
        getline(ss,substr,',');
        elementList.push_back(std::atoi(substr.c_str()));
        cout<<"____element: "<<atoi(substr.c_str())<<"____"<<endl;
    }
    std::vector<GeometryTag> tags = myCurNode->getPropertyValue<std::vector<GeometryTag>>("Geometry");
    int bodyIndex = tags[0].parentShapeNr;
    cout<<"____body index: "<<bodyIndex<<"____"<<endl;
    SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
    const occHandle(Ng_MeshVS_DataSource3D) &volumeMesh = occHandle(Ng_MeshVS_DataSource3D)::DownCast(sm->getDataBase()->ArrayOfMeshDS.value(bodyIndex));
    if(volumeMesh.IsNull())
    {
        cout<<"____NULL mesh____"<<endl;
        return;
    }
    std::vector<meshElementByCoords> listOfElements;
    for(int i = 0; i<elementList.size(); i++)
    {
        int globalElementID = elementList[i];
        int NbNodes, buf[20];
        TColStd_Array1OfInteger nodeIDs(*buf,1,20);
        volumeMesh->GetNodesByElement(globalElementID,nodeIDs,NbNodes);
        meshElementByCoords aMeshElement;
        aMeshElement.ID = globalElementID;
        for(int k = 1; k<=NbNodes; k++)
        {
            int globalNodeID = nodeIDs(k);
            int localNodeID = volumeMesh->myNodesMap.FindIndex(globalNodeID);
            const std::vector<double> &pc = volumeMesh->getNodeCoordinates(localNodeID);
            aMeshElement.pointList<<mesh::meshPoint(pc[0],pc[1],pc[2],globalNodeID);
        }
        switch(NbNodes)
        {
        case 4: aMeshElement.type = TET; break;
        case 5: aMeshElement.type = PYRAM; break;
        case 6: aMeshElement.type = PRISM; break;
        case 8: aMeshElement.type = HEXA; break;
        }
        listOfElements.push_back(aMeshElement);
    }
    //! --------------------------
    //! call the mesh constructor
    //! --------------------------
    occHandle(Ng_MeshVS_DataSource3D) miniMesh = new Ng_MeshVS_DataSource3D(listOfElements);
    cout<<"____elements: "<<miniMesh->GetAllElements().Extent()<<"____"<<endl;
    cout<<"____nodes: "<<miniMesh->GetAllNodes().Extent()<<"____"<<endl;
    //! -----------------------------------------
    //! store the mesh data source into the item
    //! -----------------------------------------
    QVariant data;
    data.setValue(miniMesh);
    Property prop_meshDS("Selected elements",data,Property::PropertyGroup_MeshDataSources);
    myCurNode->addProperty(prop_meshDS);
}

//! ----------------------------------
//! function: handleMeshMetricChanged
//! details:
//! ----------------------------------
void DetailViewer::handleMeshMetricChanged()
{
    cout<<"DetailViewer::handleMeshMetricChanged()->____function called____"<<endl;
    emit requestHandleMeshMetricChanged();
}

//! -----------------------------------
//! function: handleFatigueAlgoChanged
//! details:
//! -----------------------------------
void DetailViewer::handleFatigueAlgoChanged()
{
    QVariant data;
    int val = myCurNode->getPropertyValue<int>("Fatigue algo");
    switch(val)
    {
    //! ----------------------
    //! Basquin Coffin Manson
    //! ----------------------
    case 0:
    {
        //myCurNode->getModel()->blockSignals(true);
        if(myCurNode->getPropertyItem("Stress/strain source")==Q_NULLPTR)
        {
            data.setValue(1);   //! "Mechanical" by default
            myCurNode->addProperty(Property("Stress/strain source",data,Property::PropertyGroup_Definition),1);
        }
        //myCurNode->getModel()->blockSignals(false);
        if(myCurNode->getPropertyItem("Triaxiality correction")==Q_NULLPTR)
        {
            data.setValue(0);
            myCurNode->addProperty(Property("Triaxiality correction",data,Property::PropertyGroup_Definition),3);
        }
    }
        break;
    //! ------------------------------------
    //! Twice Yield equivalent stress range
    //! ------------------------------------
    case 1:
    {
        //myCurNode->getModel()->blockSignals(true);
        myCurNode->removeProperty("Stress/strain source");
        myCurNode->removeProperty("Triaxiality correction");
        //myCurNode->getModel()->blockSignals(false);

        //! -------------------------------------
        //! force the "Component" property value
        //! -------------------------------------
        data.setValue(0);
        myCurNode->replaceProperty("Component",Property("Component",data,Property::PropertyGroup_Definition));
    }
        break;
    }
}

//! -------------------------------------
//! function: handleEmitterStatusChanged
//! details:
//! -------------------------------------
void DetailViewer::handleEmitterStatusChanged()
{
    cout<<"DetailViewer::handleEmitterStatusChanged()->____function called____"<<endl;
    static unsigned int call;
    static double intensity;
    static double mass;
    static double charge;

    this->connectToSimulationManager(false);

    QVariant data;
    bool isEmitter = myCurNode->getPropertyValue<bool>("Emitter");
    if(isEmitter)
    {
        if(call<1)
        {
            // set default values
            mass = 9.11e-31;
            charge = 1.6e-19;
            intensity = 1e7;
        }
        data.setValue(intensity);
        myCurNode->addProperty(Property("Intensity",data,Property::PropertyGroup_Emitter));
        data.setValue(mass);
        myCurNode->addProperty(Property("Particle mass",data,Property::PropertyGroup_Emitter));
        data.setValue(charge);
        myCurNode->addProperty(Property("Electric charge",data,Property::PropertyGroup_Emitter));
    }
    else
    {
        intensity = myCurNode->getPropertyValue<double>("Intensity");
        mass = myCurNode->getPropertyValue<double>("Particle mass");
        charge = myCurNode->getPropertyValue<double>("Electric charge");
        myCurNode->removeProperty("Intensity");
        myCurNode->removeProperty("Particle mass");
        myCurNode->removeProperty("Electric charge");
    }
    this->connectToSimulationManager(true);
    call++;
}

//! -------------------------------------
//! function: connectToSimulationManager
//! details:  helper
//! -------------------------------------
void DetailViewer::connectToSimulationManager(bool toConnect)
{
    SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
    QStandardItemModel *nodeModel = myCurNode->getModel();
    if(toConnect)
    {
        disconnect(nodeModel,SIGNAL(itemChanged(QStandardItem*)),sm,SLOT(handleItemChange(QStandardItem*)));
        connect(nodeModel,SIGNAL(itemChanged(QStandardItem*)),sm,SLOT(handleItemChange(QStandardItem*)));
    }
    else
    {
        disconnect(nodeModel,SIGNAL(itemChanged(QStandardItem*)),sm,SLOT(handleItemChange(QStandardItem*)));
    }
}
