#ifndef GENERALDELEGATE_H
#define GENERALDELEGATE_H

//! ---
//! Qt
//! ---
#include <QItemDelegate>
#include <QModelIndex>
#include <QObject>
#include <QSize>
#include <QSpinBox>
#include <QWidget>
#include <QStyledItemDelegate>
#include <QPoint>

//! ----
//! OCC
//! ----
#include <AIS_InteractiveContext.hxx>
#include <TopTools_ListOfShape.hxx>

//! ----------------
//! custom includes
//! ----------------
#include "src/utils/myenumvariables.h"
#include "property.h"
#include "simulationnodeclass.h"
#include "src/main/simulationmanager.h"
#include "detailviewer.h"

class GeneralDelegate : public QStyledItemDelegate
{
    Q_OBJECT

public:

    //! constructor
    explicit GeneralDelegate(QWidget *parent=0);

    //! constuctor
    explicit GeneralDelegate(const occHandle(AIS_InteractiveContext) &aCTX, QWidget *parent=0);

    //! destructor
    ~GeneralDelegate()
    {
        //delete pointToClick;
    }

    //! Create Editor when we constructing the delegate
    QWidget* createEditor(QWidget *parent, const QStyleOptionViewItem &option, const QModelIndex &index) const override;

    //! Set the Editor
    void setEditorData(QWidget *editor, const QModelIndex &index) const override;

    //! When data are modified,  this model reflect the change
   void setModelData(QWidget *editor, QAbstractItemModel *model, const QModelIndex &index) const override;

    //! Give the widget the info on size and location
    void updateEditorGeometry(QWidget *editor, const QStyleOptionViewItem &option, const QModelIndex &index) const override;

    //! paint
    void paint(QPainter *painter, const QStyleOptionViewItem &opt, const QModelIndex &index) const override;

    //! size hint
    QSize sizeHint(const QStyleOptionViewItem & option, const QModelIndex & index) const override;

private:

    DetailViewer *myDetailViewer;    

    //! ---------------------------------------
    //! experimental - for multiple selections
    //! ---------------------------------------
    SimulationNodeClass* getCurrentMultipleSelectionNode()
    {
        //DetailViewer *dw = static_cast<DetailViewer*>(tools::getWidgetByName("detailViewer"));
        DetailViewer *dw = static_cast<DetailViewer*>(this->parent());
        SimulationNodeClass *aNode = dw->getCurrentMultipleSelectionNode();
        return aNode;
    }

    inline SimulationNodeClass* getCurrentNode() const
    {
        SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
        QModelIndex index = sm->myTreeView->currentIndex();
        SimulationNodeClass *node = index.data(Qt::UserRole).value<SimulationNodeClass*>();
        return node;
    }

    QStandardItem* getCurrenItem() const
    {
        SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
        QModelIndex index = sm->myTreeView->currentIndex();
        return sm->getModel()->itemFromIndex(index);
    }

private:

    //! the interactive context for geometry
    occHandle(AIS_InteractiveContext) myCTX;

    //! the interactive context for mesh
    occHandle(AIS_InteractiveContext) myMeshCTX;

public:

    //! set the interactive context
    void setContext(const occHandle(AIS_InteractiveContext) &aCTX);

    //! set the interactive mesh context
    void setMeshContext(const occHandle(AIS_InteractiveContext) &aMeshCTX);

private slots:

    void commitAndCloseShapeSelectorEditor();                           //! close the shape selector
    void commitAndCloseLineEdit();                                      //! close the QLineEdit
    void commitAndCloseComboBox();                                      //! close the QComboBox
    void commitAndCloseElementControlComboBox();                        //! close the QComboBox of the "Element control"
    void commitAndCloseScopingMethodComboBox();                         //! close the QComboBox of the "Scoping method"
    void commitAndCloseDefineByControlComboBox();                       //! close the QComboBox of the "Define by"
    void commitAndCloseNSSelector();                                    //! close the Named Selection selector
    void commitAndCloseCSSelector();                                    //! 02/01/2018
    void commitAndCloseSolverTypeComboBox();                            //! 06/12/2017

    void NumberOfStepSpinBoxClosed();                                   //! 25/11/2017
    void CurrentStepSpinBoxClosed();                                    //! 25/11/2017
    void StepEndTimeLineEditClosed();                                   //! 25/11/2017
    //void commitAndCloseLe();                                          //! 29/11/2017

    void commitAndCloseLe_filmCoefficient();
    void commitAndCloseLe_referenceTemperature();
    void commitAndCloseLe_magnitude();
    void commitAndCloseLe_Xcomponent();
    void commitAndCloseLe_Ycomponent();
    void commitAndCloseLe_Zcomponent();

    void commitAndCloseDirectionSelector();                             //! 11/12/2017
    void commitAndCloseComboBoxDefineBy_();                             //! 11/12/2017

    void emitOriginChanged();                                           //! 13/12/2017
    void emitOriginAndDirectionChanged();                               //! 15/12/2017
    void emitRemotePointChangedByLocation();                            //! 02/11/2018

    void commitAndCloseLineEditCSOrigin();                              //! 18/12/2017
    void commitAndCloseLineEdiCSTransformation();                       //! 18/12/2017

    void commitAndCloseAutoTimeStepping();                              //! 09/01/2018
    void commitAndCloseSubstepEditors();                                //! 09/01/2018
    void commitAndCloseVisibilityComboBox();                            //! 17/01/2018

    void commitAndCloseShowMeshNodes();                                 //! 05/03/2018
    void commitAndCloseSmoothingControl();                              //! 07/03/2018
    void commitAndCloseStepSelectionMode();                             //! 22/03/2018
    void commitAndCloseGenerate();                                      //! 22/03/2018
    //void commitAndCloseTranslate();                                     //! 09/05/2018
    void commitAndCloseSubmeshesControl();                              //! 25/03/2018
    void commitAndCloseRemapFlagControl();                              //! 21/04/2018
    //void commitAndCloseSurfaceMesherSelector();                         //! 27/04/2018
    void commitAndCloseMidsideNodesSelector();                          //! 27/04/2018
    void commitAndCloseInitialSizeSeedSelector();                       //! 01/05/2018
    void commitAndCloseRelevanceControl();                              //! 02/05/2018
    void commitAndCloseSplitFileModeSelector();                         //! 09/05/2018
    void commitAndCloseBackgroundControl();                             //! 06/06/2018
    void commitAndCloseCouplingSelector();                              //! 07/06/2018
    void commitAndCloseRPDOFselector();                                 //! 07/06/2018
    void commitAndCloseDOFswitchSelector();                             //! 07/06/2018
    void commitAndCloseNbThreads();                                     //! 08/06/2018
    void commitAndCloseSolutionInformation();                           //! 12/06/2018
    void commitAndCloseSolutionComponentSelector();                     //! 20/06/2018
    void commitAndCloseBySelector();                                    //! 22/06/2018
    void closeComboBox();                                               //! 26/06/2018
    void commitAndCloseContactToleranceEditor();                        //! 10/07/2018
    void commitAndCloseComboBoxTypeOfSizing();                          //! 21/08/2018
    void commitAndCloseStraightSidedElementsControl();                  //! 25/08/2018
    void commitAndCloseUpdateInterval();                                //! 11/09/2018
    void commitAndCloseBoltStatusDefinedBy();                           //! 16/09/2018
    void commitAndCloseBoltLoadControl();                               //! 19/09/2018
    void commitAndCloseBoltAdjustmentControl();                         //! 19/09/2018
    void commitAndCloseComboBoxFluxConvergence();                       //! 12/10/2018
    void commitAndCloseComboBoxSolutionConvergence();                   //! 12/10/2018
    void commitAndCloseFieldParametersLineEditEditor();                 //! 16/10/2018
    void commitAndCloseTimeIncrementationEditor();                      //! 17/10/2018
    void commitAndCloseTimeIncrementationControl();                     //! 17/10/2018
    void commitAndCloseCutBackParametersLineEdit();                     //! 18/10/2018
    void commitAndCloseCutBackEditor();                                 //! 18/10/2018
    void commitAndCloseLineSearchEditor();                              //! 19/10/2018
    void commitAndCloseLineParametersChanged();                         //! 19/10/2018
    void commitAndCloseSmallSlidingControl();                           //! 24/10/2018
    void commitAndCloseComboBoxForOutputSettings();                     //! 28/10/2018
    void commitAndCloseStoreResultsAt();                                //! 28/10/2018
    void commitAndCloseRemotePointsSelector();                          //! 08/11/2018
    void commitAndCloseScaleTypeSelector();                             //! 13/11/2018
    void commitAndCloseNumberOfInterval();                              //! 13/11/2018
    void commitAndCloseMinMaxControls();                                //! 13/11/2018
    void commitAndCloseMeshMethodSelector();                            //! 15/12/2018
    void commitAndCloseMeshSimplificationControl();                     //! 15/12/2018
    void commitAndCloseMeshHealingControl();                            //! 15/12/2018
    void commitAndCloseMeshDefeaturingParameterValueControl();          //! 15/12/2018
    void commitAndCloseMeshDefeaturingControl();                        //! 04/01/2019 - ON/OFF mesh defeaturing controls (healing, simplification)
    void commitAndClosePrismaticLayerOptions();                         //! 21/01/2019 - options for generatnig the prismatic layers
    void commitAndCloseBoundaryScopingMethod();                         //! 22/01/2019
    void commitAndCloseShapeSelectorEditor_boundaryPrismaticLayer();    //! 23/01/2019
    void commitAndCloseUseBRep();                                       //! 15/05/2019
    void commitAndCloseDefeaturingControl();                            //! 16/05/2019
    void commitAndCloseSimplificationControl();                         //! 17/05/2019
    void commitAndCloseVolumeMesher();                                  //! 19/05/2019
    void commitAndCloseGeometryHealing();                               //! 19/05/2019
    void commitAndCloseTessellator();                                   //! 09/06/2019
    void commitAndCloseAdjustControl();
    void displayLevelSliderValue(int pos);
    void commitAndCloseFileSelect();                                    //! 31/05/2019
    void commitAndCloseToleranceType();                                 //! 06/06/2019

#ifdef COSTAMP_VERSION
    void commitAndCloseTimeStepBuilderButton();                         //! 29/04/2019
    void commitAndCloseTimeStepBuilderFileSelector();                   //! 07/05/2019
#endif

    //void commitAndCloseShrinkFunction();                                //! 10/07/2019

    //! Tetwild
    void commitAndCloseEnvelopeSizing();                                //! 10/09/2019
    void commitAndCloseIdealLength();                                   //! 10/09/2019

    void commitAndCloseInterpolationAlgorithm();                        //! 30/09/2019
    void commitAndCloseItemType();                                      //! 08/01/2020
    void commitAndCloseContactPairSelector();                           //! 09/01/2020
    void commitAndCloseModelChangeActivationStatus();                   //! 13/01/2020
    void commitAndCloseLineEditTransparency();                          //! 28/02/2020
    void commitAndCloseSelectionMethod();                               //! 05/03/2020
    void commitAndCloseLineEditElementList();                           //! 05/03/2020
    void commitAndCloseMeshMetricCombobox();                            //! 11/04/2020
    void commitAndCloseFatigueAlgo();                                   //! 12/05/2020
    void commitAndCloseAnalysisType();                                  //! 19/05/2020
    void commitAndCloseTimeIntegrationChanged();                        //! 20/05/2020
    void commitAndCloseComboBoxEmitter();

signals:

    void suppressionChanged(Property::SuppressionStatus newSuppressioStatus) const;
    void visibilityChanged(bool newIsVisible) const;
    void visibilityChanged_();                                          //!17/01/2018
    void requestChangeNodeSuppressionStatus(Property::SuppressionStatus newSuppressionStatus);
    void elementControlChanged() const;
    void scopingMethodChanged() const;
    void scopeChanged() const;                              //! signal emitted at closing by the shape selector

    void numberOfStepChanged() const;                       //!
    void currentStepNumberChanged() const;
    void solverTypeChanged() const;
    void StepEndTimeChanged() const;

    void tabularDataChanged() const;
    void originChanged() const;
    void originChangedByValues();
    void originAndDirectionChanged();
    void transformationsChanged() const;
    void defineBy_Changed();                                //! the method for changing the origin of coordinate system coordinates has changed
    void defineByChanged();

    void loadDefinitionChanged();
    void loadDefinitionFilmCoefficientChanged(const QString &textData);
    void loadDefinitionReferenceTemperatureChanged(const QString &textData);
    void loadDefinitionMagnitudeChanged(const QString &textData);
    void loadDefinitionXChanged(const QString &textData);
    void loadDefinitionYChanged(const QString &textData);
    void loadDefinitionZChanged(const QString &textData);

    void autoTimeSteppingChanged();                         //! 09/01/2018
    void timeDivisionChanged();                             //! "

    void meshNodesVisibilityChanged(bool meshNodesVisible);                       //! 05/03/2018
    void meshSmoothingChanged();
    void stepSelectionModeChanged();
    void requestStartInterpolator();
    void requestStartOpenFoamScalarDataTranslator();
    void globalMeshControlChanged();
    void meshOrderChanged();

    //! openfoam scalar data translator: split mode => one single file => several files
    void splitModeChanged();
    void remapFlagChanged();

    void backgroundChanged();
    void DOFselectorChanged();
    void NbThreadsChanged();

    void solutionInformationChanged();
    void solutionComponentChanged();
    void byChanged();
    void automaticContactSearchToleranceChanged();
    void typeOfSizingChanged();
    void solutionInformationUpdateIntervalChanged();
    void BoltStatusDefinedByChanged();
    void BoltLoadChanged();
    void BoltAdjustmentChanged();
    void couplingChanged();
    void fluxConvegernceChanged();
    void solutionConvegernceChanged();
    void fieldParametersChanged();
    void timeIncrementationChanged();
    void cutbackFactorsChanged();
    void lineSearchChanged();
    void lineSearchParametersChanged();

    void outputControlsChanged();
    void storeResultsAtChanged();

    //! remote points
    void remotePointChangedByLocation();
    void remotePointSystemOfReferenceChanged();

    void scaleChanged();
    void colorBoxScaleChanged();
    void meshMethodChanged();
    void meshHealingChanged();
    void meshDefeaturingChanged();
    void meshSimplificationChanged();
    void meshSimplificationByChanged();
    void meshSimplificationParameterValueChanged();

    void prismaticLayerOptionsChanged();
    void boundaryScopingMethodChanged();
    void boundaryScopeChanged();

    void boltCSChanged();
    void accelerationChanged();
    void momentChanged();

    void BRepFlagChanged();
    void defeaturingFlagChanged();
    void simplificationFlagChanged();
    void volumeMesherChanged();
    void TessellatorChanged();

    //! ------------------
    //! importing options
    //! ------------------
    void geometryHealingChanged();

#ifdef COSTAMP_VERSION
    void requestStartTimeStepBuilder();
#endif

    //! -------------------
    //! TetWild parameters
    //! -------------------
    void idealLengthChanged();
    void envelopeSizingChanged();

    void interpolationAlgorithmChanged();
    void modelChangeScopingMethodChanged();
    void modelChangeActivationStatusChanged();

    void transparencyChanged();
    void selectionMethodChanged();
    void elementListChanged();

    void meshMetricChanged();
    void fatigueAlgoChanged();
    void analysisTypeChangedChanged();
    void timeIntegrationChanged();

    void EmitterStatusChanged();
};

#endif
