#ifndef DETAILVIEWER_H
#define DETAILVIEWER_H

//! ----------------
//! custom includes
//! ----------------
#include "occhandle.h"
#include "myenumvariables.h"
#include "property.h"

//! ---
//! Qt
//! ---
#include <QTreeView>
#include <QTableView>
#include <QModelIndex>
#include <QStandardItemModel>
#include <QWidget>
#include <QMouseEvent>

//! ----
//! OCC
//! ----
#include <AIS_InteractiveContext.hxx>


class GeneralDelegate;
class SimulationNodeClass;

class DetailViewer: public QTreeView
{
    Q_OBJECT

private:

    //! for mantaining a connection with the main tree
    QModelIndex myCurModelIndex;

    //! the current node
    SimulationNodeClass *myCurNode;

    //! start editing after one click
    void mousePressEvent(QMouseEvent *event);

    //! expand children
    void expandChildren(QModelIndex parent = QModelIndex());

    //! create connections
    void createConnections();

    //! remove a branch from node model
    void removeModelBranch(QStandardItem *theSeparator);

    //! connect/disconnect to simulation manager - helper4
    void connectToSimulationManager(bool toConnect);

public:

    //! constructor
    DetailViewer(QWidget *parent=0);

    //! constructor I
    DetailViewer(const occHandle(AIS_InteractiveContext) &aCTX, QWidget *parent=0);

    //! destructor
    virtual ~DetailViewer() { cout<<"DetailViewer::~DetailViewer()->____DESTRUCTOR CALLED____"<<endl; }

    //! delegate
    GeneralDelegate *myGeneralDelegate;

    //! access function
    SimulationNodeClass *getNode();

    //! get value of a property
    template<class T>
    inline T getPropertyValue(SimulationNodeClass *aNode, const QString &propertyName)
    {
        return aNode->getPropertyItem(propertyName)->data(Qt::UserRole).value<Property>().getData().value<T>();
    }

    //! hide item
    void setPropertyVisible(const QString &propertyName, bool visible=true);

    //! expand separator
    void expandSeparator(const QString &separatorName);

    //! collapse separator
    void collapseSeparator(const QString &separatorName);

    //! set current node - force it - experimental
    void setCurrentNode(SimulationNodeClass *aNode) { myCurNode = aNode; }

private:

    //! current multiple selection node - experimental
    SimulationNodeClass *myCurMultipleSelectionNode;

public:

    //! experimental - functions to be completed ...
    void setCurrentMultipleSelectionNode(SimulationNodeClass* aNode){ myCurMultipleSelectionNode = aNode; }
    SimulationNodeClass* getCurrentMultipleSelectionNode() { return myCurMultipleSelectionNode; }

public slots:

    //! set context
    void setContext(const occHandle(AIS_InteractiveContext) &aCTX);

    //! set mesh context
    void setMeshContext(const occHandle(AIS_InteractiveContext) &aMeshCTX);

    //! slot
    void setTheModel(const QModelIndex &anIndex);
    void setTheModel(SimulationNodeClass *aNode);

    //! clear the tree
    void clearTree();

    //! add property to node
    bool addPropertyToNode(const QString &name, Property::PropertyGroup thePropertyGroup);

    //! delete property from node
    void deleteTransformationFromNode();

    //! update "Field parameters"
    void updateFieldParameters();

private:

    //! interactive context
    occHandle(AIS_InteractiveContext) myCTX;

    //! interactive mesh context
    occHandle(AIS_InteractiveContext) myMeshCTX;

private slots:

    void handleSuppressionPropertyChange(Property::SuppressionStatus newSuppressionStatus);
    void handleVisibilityChange();
    void handleElementControlChange();
    void handleScopingMethodChange();
    void updateTags();

    void handleNumberOfStepChanged();
    void handleStepEndTimeChanged();
    void handleSolverTypeChanged();
    void handleCurrentStepNumberChanged();
    void handleAnalysisTypeChanged();
    void handleTimeIntegrationChanged();

    void selectionHasChanged();
    void handleOriginChanged();
    void handleOriginAndDirectionChanged();
    void handleOriginChangedByValue();

    void handleDefineBy_Changed();
    void handleDefineByChanged();

    QVector<double> getOrigin();
    QVector<QVector<double>> getDirectionalData();

    void updateBaseOriginCoords();
    void updateBaseDirectionalData();
    void updateBaseData();

    void applyTransformation(Property::typeOfTransformation theType, double delta);
    void applyTransformations();
    void updateCS();

    void handleLoadDefinitionChanged();
    void handleFilmCoefficientLoadDefinitionChanged(const QString& textData);
    void handleReferenceTemperatureLoadDefinitionChanged(const QString& textData);
    void handleMagnitudeLoadDefinitionChanged(const QString &textData);
    void handleXLoadDefinitionChanged(const QString &textData);
    void handleYLoadDefinitionChanged(const QString& textData);
    void handleZLoadDefinitionChanged(const QString& textData);

    void handleAutoTimeSteppingChanged();
    void handleTimeDivisionChanged();
    void updateDetailViewerFromTabularData(QModelIndex topLeftIndex, QModelIndex bottomRightIndex, QVector<int> roles);
    void handleChangeMeshNodesVisibility(bool meshNodesVisible);
    void handleMeshSmoothingChange();
    void handleRequestStartInterpolator();
    void handleRequestStartOpenFoamScalarDataTranslator();
    void handleGlobalMeshControlChange();

    void handleRemapFlagChanged();
    void handleStepSelectionModeChange();

    void postBackgroundChangeEvent();

    void handleNbThreadsChanged();
    //void handleByChanged();
    void handleSolutionComponentChanged();
    void handleTypeOfSizingChanged();
    void handleSolutionInformationUpdateIntervalChanged() { emit requestHandleSolutionInformationUpdateIntervalChanged(); }

    //! ------
    //! bolts
    //! ------
    void handleBoltStatusDefinedByChanged();
    void handleBoltLoadChanged();
    void handleBoltAdjustmentChanged();
    void handleBoltCSChanged();             // marker changed

    //! -------------
    //! acceleration
    //! -------------
    void handleAccelerationChanged();       // marker changed

    //! -------
    //! moment
    //! -------
    void handleMomentChanged();             // moment changed

    void handleFluxConvergenceChanged();
    void handleSolutionConvergenceChanged();
    void updateTimeIncrementationParameters();
    void updateCutBackFactors();
    void updateLineSearch();

    void handleOutputControlsChanged();
    void handleStoreResultsAtChanged();

    //! --------------
    //! remote points
    //! --------------
    void handleRemotePointChangedByLocation();
    void handleRemotePointSystemOfReferenceChanged();

    void handleMeshMethodChanged();
    void handleMeshDefeaturingChanged();
    void handleMeshSimplificationChanged();
    void handleMeshSimplificationByChanged();

    //! -----------------
    //! prismatic layers
    //! -----------------
    void handlePrismaticLayerOptions();
    void handleBoundaryScopingMethodChanged();
    void updateBoudaryTags();

    //! ------------
    //! mesh method
    //! ------------
    void handleBRepFlagChanged();
    void handleDefeaturingFlagChanged();
    void handleSimplificationFlagChanged();
    void handleVolumeMesherChanged();
    void handleTessellatorChanged();

    //! ---------------------------
    //! geometry importing options
    //! ---------------------------
    void handleGeometryHealingChanged();

    //! -------------
    //! model change
    //! -------------
    void handleModelChangeScopingMethodChanged();
    void handleModelChangeActivationStatusChanged();

    //! -----------------------
    //! mesh element selection
    //! -----------------------
    void handleSelectionMethodChanged();
    void handleMeshElementListChanged();

    //! --------------------
    //! mesh metric display
    //! --------------------
    void handleMeshMetricChanged();

    //! -------------
    //! fatigue algo
    //! -------------
    void handleFatigueAlgoChanged();

    //! --------
    //! emitter
    //! --------
    void handleEmitterStatusChanged();

#ifdef COSTAMP_VERSION
    void handleRequestStartTimeStepBuilder();
#endif

public slots:

    //! ---------------------------------------------------------------------------------
    //! at the moment the following task have been transferred to the simulation manager
    //! by consequence the methods are public
    //! ---------------------------------------------------------------------------------
    //! remote points, remote force/displacement
    void handleDOFselectorChange();
    void handleCouplingChanged();
    void handleColorBoxScaleChanged();

    //! -----------------------
    //! "By" time/step-substep
    //! -----------------------
    void handleByChanged();

    //! -------------------
    //! TetWild parameters
    //! -------------------
    void handleIdealLengthChanged();
    void handleEnvelopeSizingChanged();

    void handleInterpolationAlgorithmChanged();
    void handleTransparencyChanged();

signals:

    void requestChangeNodeSuppressionStatus(Property::SuppressionStatus newSuppressionStatus);
    void requestHandleVisibilityChange(bool isVisibile);
    void requestChangeElementControl();
    void requestChangeColor();

    //! connect to the "SimulationManager" class
    void requestResizeTabularData();

    //! connect these signals to the "SimulationManager" class
    void requestHandleStepEndTimeChanged();
    //void requestHandleSolverTypeChanged();
    void requestHandleCurrentStepNumberChanged();

    //! signal to be connected to the "DirectionSelector" class
    void requestHandleSelectionChanged();
    //! signal to be connected to the "myOccPreGLWidget" class
    void requestDisplayTrihedron(QVector<double>,QVector<QVector<double>>, int size=0);

    void requestHandleTabularData();
    void requestHandleLoadDefinitionChanged();
    void requestHandleFilmCoefficientLoadDefinitionChanged(const QString& textData);
    void requestHandleReferenceTemperatureLoadDefinitionChanged(const QString &textData);
    void requestHandleMagnitudeLoadDefinitionChanged(const QString& textData);
    void requestHandleXLoadDefinitionChanged(const QString& textData);
    void requestHandleYLoadDefinitionChanged(const QString& textData);
    void requestHandleZLoadDefinitionChanged(const QString& textData);
    void requestChangeMeshNodesVisibility(bool meshNodesVisible);
    void requestInvalidateAllMeshes();
    void startInterpolator();
    void startOpenFoamScalarDataTranslator();
    void requestGlobalMeshControlChange();
    void requestHandleSolutionInformationUpdateIntervalChanged();
    //void requestHandleBoltControls();     to be removed
    void requestHideAllMarkers(bool = false);
    void requestHandleTransparencyChanged(double aLevel);
    void requestHandleMeshMetricChanged();

#ifdef COSTAMP_VERSION
    void requestStartTimeStepBuilder();
#endif
};

#endif // DETAILVIEWER_H
