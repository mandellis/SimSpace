#ifdef COSTAMP_VERSION
#include "src/TimeStepBuilder/timestepbuilder.h"
#endif

//! ----------------
//! custom includes
//! ----------------
#include "global.h"
#include "simulationmanager.h"
#include "qextendedstandarditem.h"
#include "AIS_ExtendedShape.h"
#include "ais_colorscaleextended.h"
#include "mydefines.h"
#include "property.h"
#include "customtablemodel.h"
#include "load.h"
#include "geomtoolsclass.h"
#include "markers.h"
#include "deserializerclass.h"
#include "nodefactory.h"
#include "tools.h"
#include "openfoamreader.h"
#include "contextmenubuilder.h"
//#include "ccxsolvermanager.h"
#include "ccxconsoletofile.h"
#include "meshingserver.h"
#include "postengine.h"
#include "posttools.h"
#include "occmeshtoccxmesh.h"
#include "qprogressevent.h"
#include "frdreader.h"
#include "contactfinder.h"
#include "parser.h"
#include "runterminationdata.h"
#include "qconsoleevent.h"
#include "solutioninfo.h"
#include "modelloader.h"
#include "tablewidget.h"
#include "tabulardatacolumns.h"
#include <src/cliptool/clipTool.h>
#include "handle_ais_doublearrowmarker_reg.h"
#include "markerbuilder.h"
#include "maintreetools.h"
#include <connectionpairgenerationoptions.h>
#include <tabulardataviewerclass1.h>

#include <ng_mesher2.h>
#include <occmesher.h>
#include <edgemeshdatasourcebuildesclass.h>
#include <ng_meshvs_datasource3d.h>
#include <mesherclass.h>
#include <facedatasourcebuilder.h>
#include <elementtypes.h>
#include "meshvs_mesh_handle_reg.h"

#include "topologytools.h"
#include "meshtools.h"
#include "graphicstools.h"
#include "mapper3dclass.h"
#include "listofmesh.h"
#include "ccout.h"
#include "qsimulationstatusevent.h"
#include "tablewidget.h"
#include "generaldelegate.h"

#include "qccxsolvermessageevent.h"
#include "ccxsolvermanager1.h"
#include "ccxsolvermessage.h"
#include "ccxtools.h"

#include "inputfilegenerator.h"

//! ----
//! C++
//! ----
#include <set>
#include <algorithm>
#include <chrono>

//! -------
//! quazip
//! -------
#include <JlCompress.h>

//! ----
//! OCC
//! ----
#include <AIS_Shape.hxx>
#include <TopExp_Explorer.hxx>
#include <Geom_Surface.hxx>
#include <BRep_Tool.hxx>
#include <GeomAdaptor_Surface.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>
#include <TopExp.hxx>
#include <BRepTools.hxx>
#include <TCollection_AsciiString.hxx>
#include <MeshVS_NodalColorPrsBuilder.hxx>
#include <MeshVS_BuilderPriority.hxx>
#include <MeshVS_MeshPrsBuilder.hxx>
#include <Quantity_Color.hxx>
#include <MeshVS_Drawer.hxx>
#include <MeshVS_DrawerAttribute.hxx>
#include <MeshVS_DisplayModeFlags.hxx>
#include <Aspect_TypeOfColorScalePosition.hxx>
#include <TopTools_MapIteratorOfMapOfShape.hxx>
#include <BRep_Builder.hxx>

//! ---
//! Qt
//! ---
#include <QStandardItemModel>
#include <QStandardItem>
#include <QItemSelection>
#include <QHeaderView>
#include <QTreeView>
#include <QTableView>
#include <QVBoxLayout>
#include <QMenu>
#include <QMessageBox>
#include <QVector>
#include <QFileDialog>
#include <QProgressDialog>
#include <QLabel>
#include <QDirIterator>
#include <QFile>
#include <QDir>
#include <QApplication>
#include <QProcess>
#include <QTimer>
#include <QItemDelegate>

//! ----
//! C++
//! ----
#include <iostream>
using namespace std;


#include "wbtrihedron.h"
//! ----------------------
//! function: constructor
//! details:  default
//! ----------------------
SimulationManager::SimulationManager(QWidget *parent): QWidget(parent)
{
    cout<<"SimulationManager::SimulationManager()->____default constructor called____"<<endl;

    //! -------------------
    //! create the widgets
    //! -------------------
    this->createContent();

    //! -------------
    //! event filter
    //! -------------
    this->installEventFilter(this);

    //! ----------------------------------------------------
    //! serializer, deserializer, postEngine, meshingServer
    //! ----------------------------------------------------
    mySerializer = new serializerClass(this);
    myDeserializer = new deserializerClass(this);
    myPostEngine = new postEngine(this);
    QProgressIndicator *aProgressIndicator = static_cast<QProgressIndicator*>(tools::getWidgetByName("progressIndicator"));
    myMeshingServer = new MeshingServer(aProgressIndicator,this);
    myInputFileGenerator = new inputFileGenerator(this);
    myInputFileGenerator->setProgressIndicator(aProgressIndicator);

    //! -----------------
    //! status variables
    //! -----------------
    myIsMeshingRunning = false;
    myIsCalculationRunning = false;

    //! -----------------------------------------------------------------------------
    //! stop the timer when the solution is ready (present or not, converged or not)
    //! -----------------------------------------------------------------------------
    myTimer = new QTimer(this);
    connect(myTimer,SIGNAL(timeout()),this,SLOT(retrieveSolverInfo()));

    //! -------------------------------------
    //! create an empty simulation data base
    //! -------------------------------------
    //this->createSimulationDataBaseEmpty();
    //this->createSimulationNode(SimulationNodeClass::nodeType_timeStepBuilder);

    //! -----------------------------------
    //! init the current project directory
    //! -----------------------------------
    //myCurrentProjectDir = QString("");

    //! -----------------------------
    //! the current running analysis
    //! -----------------------------
    myCurrentRunningAnalysis = Q_NULLPTR;
}

//! ------------------------
//! function: constructor I
//! details:
//! ------------------------
SimulationManager::SimulationManager(const occHandle(AIS_InteractiveContext) &aCTX, QWidget *parent): QWidget(parent), myCTX(aCTX)
{
    cout<<"SimulationManager::SimulationManager()->____constructor I called____"<<endl;

    //! -------------------
    //! create the widgets
    //! -------------------
    this->createContent();

    //! -------------
    //! event filter
    //! -------------
    this->installEventFilter(this);

    mySerializer = new serializerClass(this);
    myDeserializer = new deserializerClass(this);
    myPostEngine = new postEngine(this);

    QProgressIndicator *aProgressIndicator = static_cast<QProgressIndicator*>(tools::getWidgetByName("progressIndicator"));
    myMeshingServer = new MeshingServer(aProgressIndicator,this);
    myInputFileGenerator = new inputFileGenerator(this);
    myInputFileGenerator->setProgressIndicator(aProgressIndicator);

    myTimer = new QTimer(this);

    //! -----------------
    //! status variables
    //! -----------------
    myIsMeshingRunning = false;
    myIsCalculationRunning = false;

    //! -----------------------------------
    //! init the current project directory
    //! -----------------------------------
    //myCurrentProjectDir = QString("");

    //! -----------------------------
    //! the current running analysis
    //! -----------------------------
    myCurrentRunningAnalysis =  Q_NULLPTR;
}

//! ----------------------------
//! function: setSelectionModel
//! details:
//! ----------------------------
void SimulationManager::setSelectionModel()
{
    //! -------------------------------------------
    //! retrieve the selection model from the tree
    //! -------------------------------------------
    mySelectionModel = myTreeView->selectionModel();
}

//! ----------------------
//! function: highlighter
//! details:
//! ----------------------
void SimulationManager::highlighter(QModelIndex modelIndex)
{
    //cout<<"SimulationManager::highlighter()->____function called____"<<endl;
    Q_UNUSED(modelIndex)

    //! -------------------------------
    //! set the active analysis branch
    //! -------------------------------
    this->setTheActiveAnalysisBranch();

    QList<QModelIndex> selectedIndexes = myTreeView->selectionModel()->selectedRows();

    //! ----------------------------------------------------------------
    //! try to handle multiple selections within a group - experimental - gildotta
    //! ----------------------------------------------------------------
    if(selectedIndexes.length()>1)
    {
        cout<<"@---------------------------------@"<<endl;
        cout<<"@- handling a multiple selection -@"<<endl;
        cout<<"@---------------------------------@"<<endl;
        bool areAllConnectionPairs = true;
        for(int i=0; i<selectedIndexes.length(); i++)
        {
            QModelIndex modelIndex = selectedIndexes.at(i);
            SimulationNodeClass *node = modelIndex.data(Qt::UserRole).value<SimulationNodeClass*>();
            if(node->getType()!=SimulationNodeClass::nodeType_connectionPair)
            {
                areAllConnectionPairs = false;
                break;
            }
        }
        if(areAllConnectionPairs==true)
        {
            cout<<"@-----------------------------------------------------@"<<endl;
            cout<<"@- handling a multiple selection of connection pairs -@"<<endl;
            cout<<"@-----------------------------------------------------@"<<endl;
            QMap<QString,QVariant> propNameToData;
            QMap<QString,Property::PropertyGroup> propNameToGroup;
            int NbIndexes = selectedIndexes.length();
            for(int i=0; i<NbIndexes; i++)
            {
                SimulationNodeClass *node = modelIndex.data(Qt::UserRole).value<SimulationNodeClass*>();
                QVector<Property> vecProp = node->getProperties();
                for(QVector<Property>::iterator it = vecProp.begin(); it!= vecProp.end(); it++)
                {
                    const Property &aProp = *it;
                    QString name = aProp.getName();
                    //cout<<"____name: "<<name.toStdString()<<"____"<<endl;
                    QVariant data = aProp.getData();
                    Property::PropertyGroup group = aProp.getPropertyGroupEnum();
                    propNameToData.insertMulti(name,data);
                    propNameToGroup.insertMulti(name,group);
                }
            }

            //! --------------------------------------------------------------
            //! check if a property is contained into the map NbIndexes times
            //! --------------------------------------------------------------
            QVector<Property> vecProps;
            QList<QString> names = propNameToData.keys();
            QList<QString> alreadyInsertedNames;
            for(int n=0; n<names.length(); n++)
            {
                QString aName = names.at(n);
                if(propNameToData.values(aName).length()==NbIndexes)
                {
                    if(alreadyInsertedNames.contains(aName)) continue;
                    Property::PropertyGroup aPropGroup = propNameToGroup.value(aName);
                    QVariant data = QVariant();
                    Property aProp(aName,data,aPropGroup);
                    vecProps.push_back(aProp);
                    alreadyInsertedNames<<aName;
                    cout<<"____adding property: "<<aProp.getName().toStdString()<<"____"<<endl;
                }
            }
            //! --------------------------------------------------------------------------
            //! now build a dummy node and force it as the current into the detail viewer
            //! --------------------------------------------------------------------------
            SimulationNodeClass *aDummyNode = new SimulationNodeClass("Multiple contact selection",SimulationNodeClass::nodeType_connectionPair,vecProps,this);
            connect(aDummyNode->getModel(),SIGNAL(itemChanged(QStandardItem*)),this,SLOT(handleItemChange(QStandardItem*)));
            DetailViewer *detailViewer = static_cast<DetailViewer*>(tools::getWidgetByName("detailViewer"));
            detailViewer->setCurrentMultipleSelectionNode(aDummyNode); // set the current multiple selection node
            detailViewer->setTheModel(aDummyNode);

            return;
        }
    }
    //! -----------------
    //! end experimental
    //! -----------------

    if(!selectedIndexes.isEmpty())
    {
        DetailViewer *detailViewer = static_cast<DetailViewer*>(tools::getWidgetByName("detailViewer"));
        detailViewer->setCurrentMultipleSelectionNode(Q_NULLPTR);

        //! ----------------------------------------------
        //! for the moment the detail viewer is activated
        //! only for the last element of the selection
        //! ----------------------------------------------
        const QModelIndex &index = selectedIndexes.last();

        if(index.isValid() && index.data(Qt::UserRole).isValid())
        {
            //! ------------------------------------------
            //! hide the double view port - experimental
            //! initially hide the double view port:
            //! requestShowDoubleViewPort(true) is called
            //! only if a connectionPair has been clicked
            //! ------------------------------------------
            emit requestShowDoubleViewPort(false);

            //! ---------------------------------------------
            //! transfer the node model to the detail viewer
            //! ---------------------------------------------
            emit updatedetailViewer(index);

            QStandardItem *curItem = myModel->itemFromIndex(index);
            SimulationNodeClass *theNode = curItem->data(Qt::UserRole).value<SimulationNodeClass*>();
            SimulationNodeClass::nodeType theNodeType = theNode->getType();
            cout<<"SimulationManager::highlighter()->____highlight "<<theNode->getName().toStdString()<<"____"<<endl;

            //! -------------------------------------------------
            //! build a context menu suitable for the node type
            //! this function has been removed from here and put
            //! into showContextMenu() for synch reasons. Line
            //! left here for documentation
            //! -------------------------------------------------
            //this->buildCustomMenu(index);

            //! ---------------------------------------
            //! switch the panel of the central widget
            //! ---------------------------------------
            switch(theNodeType)
            {
            case SimulationNodeClass::nodeType_StructuralAnalysisSolution:
            case SimulationNodeClass::nodeType_thermalAnalysisSolution:
            case SimulationNodeClass::nodeType_combinedAnalysis:
            {
                emit requestSetActiveCentralTab("maingwindow");
                emit requestHideAllResults();
            }
                break;

            case SimulationNodeClass::nodeType_StructuralAnalysisSolutionInformation:
            case SimulationNodeClass::nodeType_thermalAnalysisSolutionInformation:
            case SimulationNodeClass::nodeType_combinedAnalysisSolutionInformation:
            {
                emit requestSetActiveCentralTab("worksheetViewer");
            }
                break;

            default:
            {
                emit requestSetActiveCentralTab("maingwindow");
            }
                break;
            }

            //! -----------------
            //! hide all markers
            //! -----------------
            emit requestHideAllMarkers(true);

            //! -------------------
            //! actually highlight
            //! -------------------
            switch(theNodeType)
            {                
            case SimulationNodeClass::nodeType_geometryBody:
            {
                //! ------------------------------------------------------------
                //! handle the working modes
                //! "0" => Mesh; "1" => Contacts; "2" => Model; "3" => Solution
                //! ------------------------------------------------------------
                emit requestSetWorkingMode(2);
                emit requestHideAllResults();
                emit requestHideSlicedMeshes();

                //! --------------------------------------------
                //! prepare for highlighting: this version uses
                //! the "Map index". Can also use "Geometry"
                //! --------------------------------------------
                QList<int> listOfBodies;
                for(QList<QModelIndex>::iterator it = selectedIndexes.begin();  it!=selectedIndexes.end(); ++it)
                {
                    QModelIndex anIndex = *it;
                    SimulationNodeClass *aNode = anIndex.data(Qt::UserRole).value<SimulationNodeClass*>();
                    QStandardItem *itemMapIndex = aNode->getPropertyItem("Map index");
                    if(itemMapIndex!=Q_NULLPTR)
                    {
                        int mapIndex = aNode->getPropertyValue<int>("Map index");
                        listOfBodies<<mapIndex;
                    }
                }
                emit requestHighlightBody(listOfBodies);
            }
                break;

            case SimulationNodeClass::nodeType_StructuralAnalysisSolution:
            case SimulationNodeClass::nodeType_thermalAnalysisSolution:
            case SimulationNodeClass::nodeType_combinedAnalysisSolution:
            {
                emit requestUnhighlightBodies(false);
                emit requestSetWorkingMode(2);
                emit requestHideAllResults();
                emit requestHideSlicedMeshes();
                emit requestClearGraph();
                this->changeColor();
            }
                break;

            case SimulationNodeClass::nodeType_StructuralAnalysisSolutionInformation:
            case SimulationNodeClass::nodeType_thermalAnalysisSolutionInformation:
            case SimulationNodeClass::nodeType_combinedAnalysisSolutionInformation:
            {
                CCXSolverMessage solverOutput = theNode->getPropertyValue<CCXSolverMessage>("Solver output");
                emit requestDisplayTextOnConsole(solverOutput.myText);
                emit requestUnhighlightBodies(false);
                emit requestSetWorkingMode(2);
                emit requestHideAllResults();
                emit requestHideSlicedMeshes();
                emit requestClearGraph();
                this->changeColor();
            }
                break;

            case SimulationNodeClass::nodeType_OpenFoamScalarData:
            {
                emit requestSetActiveCentralTab("maingwindow");
            }
                break;

            case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_ImportedTemperatureDistribution:
            {
                //! ---------------------------------------------------------------------
                //! hide the meshes: keep the bodies in wireframe mode for the selection
                //! ---------------------------------------------------------------------
                emit requestHideMeshes();
                emit requestHideAllResults();
                emit requestHideSlicedMeshes();

                SimulationNodeClass *node = index.data(Qt::UserRole).value<SimulationNodeClass*>();
                cout<<"____name: "<<node->getName().toStdString()<<"____"<<endl;
                QStandardItem *itemTemperature = (QStandardItem*)(node->getPropertyValue<void*>("Imported body temperature"));
                cout<<"____item name: "<<itemTemperature->data(Qt::DisplayRole).toString().toStdString()<<"____"<<endl;
                SimulationNodeClass *nodeTemperature = itemTemperature->data(Qt::UserRole).value<SimulationNodeClass*>();
                cout<<"____contained node: "<<nodeTemperature->getName().toStdString()<<"____"<<endl;
                cout<<"____listing the properties of the contained node____"<<endl;
                for(int i=0; i<nodeTemperature->getPropertyItems().length(); i++)
                {
                    cout<<"____"<<nodeTemperature->getPropertyItems().at(i)->data(Qt::UserRole).value<Property>().getName().toStdString()<<"____"<<endl;
                }
                if(nodeTemperature->getPropertyItem("Post object")!=Q_NULLPTR)
                {
                    postObject aPostObject = nodeTemperature->getPropertyValue<postObject>("Post object");
                    emit requestSetWorkingMode(3);
                    emit requestShowAllBodies();    //! check if wireframe... to do
                    emit requestDisplayResult(aPostObject);
                }
                else
                {
                    cout<<"____no data____"<<endl;
                    requestSetWorkingMode(2);
                }
                this->changeColor();
            }
                break;

            case SimulationNodeClass::nodeType_solutionThermalFlux:
            case SimulationNodeClass::nodeType_solutionThermalTemperature:
            {
                //! ---------------------------------------------------------------------
                //! hide the meshes: keep the bodies in wireframe mode for the selection
                //! ---------------------------------------------------------------------
                emit requestHideMeshes();
                emit requestHideAllResults();
                emit requestHideSlicedMeshes();

                //! -------------------------------------------------------
                //! in the working mode "3" "Solution" the selection modes
                //! and the corresponding toolbar buttons are disabled
                //! -------------------------------------------------------
                postObject aPostObject;
                bool isDone = this->retrieveCurrentItemResult(aPostObject);
                if(isDone)
                {
                    emit requestSetWorkingMode(3);
                    emit requestShowAllBodies();    //! check if wireframe... to do
                    emit requestDisplayResult(aPostObject);
                }
                else requestSetWorkingMode(2);

                //! ---------------
                //! switch the tab
                //! ---------------
                emit requestSetActiveCentralTab("maingwindow");
            }
                break;

            case SimulationNodeClass::nodeType_importedBodyScalar:
            case SimulationNodeClass::nodeType_postObject:
            case SimulationNodeClass::nodeType_solutionStructuralNodalDisplacement:
            case SimulationNodeClass::nodeType_solutionStructuralStress:
            case SimulationNodeClass::nodeType_solutionStructuralTotalStrain:
            case SimulationNodeClass::nodeType_solutionStructuralThermalStrain:
            case SimulationNodeClass::nodeType_solutionStructuralMechanicalStrain:
            case SimulationNodeClass::nodeType_solutionStructuralTemperature:
            case SimulationNodeClass::nodeType_solutionStructuralEquivalentPlasticStrain:
            case SimulationNodeClass::nodeType_solutionStructuralContact:
            case SimulationNodeClass::nodeType_solutionStructuralFatigueTool:
            case SimulationNodeClass::nodeType_solutionStructuralGamma:
            case SimulationNodeClass::nodeType_solutionStructuralNodalForces:
            case SimulationNodeClass::nodeType_solutionStructuralReactionForce:
            {
                //! ---------------------------------------------------------------------
                //! hide the meshes: keep the bodies in wireframe mode for the selection
                //! ---------------------------------------------------------------------
                emit requestHideMeshes();
                emit requestHideAllResults();
                emit requestHideSlicedMeshes();

                //! -------------------------------------------------------
                //! in the working mode "3" "Solution" the selection modes
                //! and the corresponding toolbar buttons are disabled
                //! -------------------------------------------------------
                postObject aPostObject;
                bool isDone = this->retrieveCurrentItemResult(aPostObject);
                if(isDone)
                {
                    emit requestSetWorkingMode(3);
                    emit requestShowAllBodies();    //! check if wireframe... to do
                    emit requestDisplayResult(aPostObject);
                }
                else requestSetWorkingMode(2);

                //! ---------------
                //! switch the tab
                //! ---------------
                emit requestSetActiveCentralTab("maingwindow");
            }
                break;

            case SimulationNodeClass::nodeType_mapper:
            {
                emit requestUnhighlightBodies(true);

                //! hide the meshes: keep the bodies for selection
                emit requestSetWorkingMode(3);
                emit requestShowAllBodies();
                emit requestHideAllResults();
                emit requestHideSlicedMeshes();
                this->changeColor();

                //! switch the tab
                emit requestSetActiveCentralTab("maingwindow");
            }
                break;

            case SimulationNodeClass::nodeType_root:
            {
                emit requestUnhighlightBodies(false);
                emit requestHideMeshes();
                emit requestHideSlicedMeshes();
                emit requestSetWorkingMode(2);
                emit requestHideAllResults();
                this->changeColor();
            }
                break;

            case SimulationNodeClass::nodeType_geometry:
            {
                emit requestUnhighlightBodies(false);
                emit requestHideSlicedMeshes();
                emit requestSetWorkingMode(2);
                emit requestHideAllResults();
                this->changeColor();
            }
                break;

            case SimulationNodeClass::nodeType_coordinateSystems:
            {
                emit requestHideAllResults();
                emit requestUnhighlightBodies(true);
                emit requestHideMeshes();
                emit requestHideSlicedMeshes();
                emit requestSetWorkingMode(2);
                this->changeColor();
            }
                break;

            case SimulationNodeClass::nodeType_coordinateSystem_global:
            {
                emit requestHideAllResults();
                emit requestUnhighlightBodies(Standard_True);
                emit requestHideMeshes();
                emit requestHideSlicedMeshes();
                emit requestSetWorkingMode(2);

                //! -------------------------------------
                //! display the global coordinate system
                //! -------------------------------------
                QVector<double> origin{0,0,0};
                QVector<double> directionalDataX{1,0,0};
                QVector<double> directionalDataY{0,1,0};
                QVector<double> directionalDataZ{0,0,1};

                QVector<QVector<double>> dirData;
                dirData.push_back(directionalDataX);
                dirData.push_back(directionalDataY);
                dirData.push_back(directionalDataZ);
                emit requestDisplayTrihedron(origin,dirData);
            }
                break;

            case SimulationNodeClass::nodeType_coordinateSystem:
            {
                emit requestHideAllResults();
                emit requestUnhighlightBodies(Standard_True);
                emit requestHideMeshes();
                emit requestHideSlicedMeshes();
                emit requestSetWorkingMode(2);
                this->changeColor();

                //! ------------------------------
                //! display the coordinate system
                //! ------------------------------
                SimulationNodeClass *theNode = myTreeView->currentIndex().data(Qt::UserRole).value<SimulationNodeClass*>();
                QVector<double> baseOrigin = theNode->getPropertyItem("Base origin")->data(Qt::UserRole).value<Property>().getData().value<QVector<double>>();
                QVector<QVector<double>> baseDirectionalData = theNode->getPropertyItem("Base directional data")->data(Qt::UserRole).value<Property>().getData().value<QVector<QVector<double>>>();
                //! ---------------------------
                //! show the coordinate system
                //! ---------------------------
                emit requestDisplayTrihedron(baseOrigin,baseDirectionalData);
            }
                break;

            case SimulationNodeClass::nodeType_remotePointRoot:
            {
                emit requestHideAllResults();
                emit requestUnhighlightBodies(true);
                emit requestHideMeshes();
                emit requestHideSlicedMeshes();
                emit requestSetWorkingMode(2);
            }
                break;

            case SimulationNodeClass::nodeType_remotePoint:
            {
                emit requestHideAllResults();
                emit requestUnhighlightBodies(false);
                emit requestHideMeshes();
                emit requestHideSlicedMeshes();
                emit requestSetWorkingMode(2);

                theNode->getModel()->blockSignals(true);    //! avoids calling handleItemChange()
                bool isDone = markerBuilder::addMarker(this->getCurrentNode(), mySimulationDataBase);
                theNode->getModel()->blockSignals(false);   //! reconnect - unblock signals

                if(isDone == true) this->displayMarker();
            }
                break;

            case SimulationNodeClass::nodeType_pointMass:
            {
                emit requestSetWorkingMode(2);
                emit requestHideAllResults();
                emit requestHideSlicedMeshes();
                theNode->getModel()->blockSignals(true);    //! avoids calling handleItemChange()
                bool isDone = markerBuilder::addMarker(this->getCurrentNode(), mySimulationDataBase);
                theNode->getModel()->blockSignals(false);   //! reconnect - unblock signals
                if(isDone == true) this->displayMarker();
            }
                break;

            case SimulationNodeClass::nodeType_meshControl:
            {
                emit requestHideAllResults();
                emit requestUnhighlightBodies(true);
                emit requestSetWorkingMode(0);

                SimulationNodeClass *meshRootNode = this->getTreeItem(SimulationNodeClass::nodeType_meshControl)->data(Qt::UserRole).value<SimulationNodeClass*>();
                bool areMeshNodeVisible = meshRootNode->getPropertyValue<bool>("Show mesh nodes");
                emit requestShowMeshes(areMeshNodeVisible);
                this->changeColor();
                this->requestClearGraph();
            }
                break;

            case SimulationNodeClass::nodeType_meshBodyMeshControl:
            case SimulationNodeClass::nodeType_meshMethod:
            case SimulationNodeClass::nodeType_meshPrismaticLayer:
            case SimulationNodeClass::nodeType_meshVertexSize:
            case SimulationNodeClass::nodeType_meshEdgeSize:
            case SimulationNodeClass::nodeType_meshFaceSize:
            case SimulationNodeClass::nodeType_meshGrading:
            case SimulationNodeClass::nodeType_meshMaxElementSize:
            case SimulationNodeClass::nodeType_meshMinElementSize:
            case SimulationNodeClass::nodeType_meshGenericSizing:
            {
                emit requestHideAllResults();
                emit requestUnhighlightBodies(Standard_True);
                emit requestSetWorkingMode(0);
                emit requestHideMeshes();
                emit requestHideSlicedMeshes();
                this->changeColor();
                emit requestClearGraph();
            }
                break;

            case SimulationNodeClass::nodeType_thermalAnalysis:
            case SimulationNodeClass::nodeType_structuralAnalysis:
            case SimulationNodeClass::nodeType_combinedAnalysis:
            {
                emit requestHideAllResults();
                emit requestUnhighlightBodies(Standard_True);
                emit requestSetWorkingMode(2);
                emit requestHideMeshes();
                emit requestHideSlicedMeshes();
                this->changeColor();
                emit requestClearGraph();
            }
                break;

            case SimulationNodeClass::nodeType_structuralAnalysisSettings:
            case SimulationNodeClass::nodeType_thermalAnalysisSettings:
            case SimulationNodeClass::nodeType_combinedAnalysisSettings:
            {
                cout<<"SimulationManager::highlighter()->____analysis settings selected____"<<endl;
                emit requestHideAllResults();
                emit requestUnhighlightBodies(Standard_True);
                emit requestSetWorkingMode(2);
                emit requestHideMeshes();
                emit requestHideSlicedMeshes();
                this->changeColor();

                //! -------------------------------------------------------------
                //! set the model for the load table: the argument of the SIGNAL
                //! MUST BE the index of an "Analysis settings" item
                //! Here "index" is already the index of "Analysis settings"
                //! -------------------------------------------------------------
                emit requestTabularData(index);

                QList<int> columnsToShow;
                CustomTableModel *tabData = index.data(Qt::UserRole).value<SimulationNodeClass*>()->getTabularDataModel();
                if(tabData!=Q_NULLPTR)
                {
                    //! ----------------------------------
                    //! show all the columns of the table
                    //! ----------------------------------
                    for(int k=0; k<tabData->columnCount(); k++) columnsToShow<<k;

                    //! ------------------------------------------
                    //! show only the first three columns instead
                    //! ------------------------------------------
                    //columnsToShow<<0<<1<<2;

                    //! show the columns of the table
                    emit requestShowColumns(columnsToShow);

                    //! hide the row with Time = 0
                    emit requestHideFirstRow();
                    emit requestClearGraph();
                }
            }
                break;

            case SimulationNodeClass::nodeType_thermalAnalysisAdiabaticWall:
            {
                emit requestHideAllResults();
                emit requestUnhighlightBodies(true);
                emit requestHideMeshes();
                emit requestHideSlicedMeshes();
                emit requestSetWorkingMode(2);
                this->changeColor();
                emit requestClearGraph();

                //! --------------
                //! set the model
                //! --------------
                QModelIndex index_analysisSettings = this->getAnalysisSettingsItemFromCurrentItem()->index();
                emit requestTabularData(index_analysisSettings);
                emit requestHideFirstRow();

                emit requestClearGraph();

                bool isDone = markerBuilder::addMarker(this->getCurrentNode(), mySimulationDataBase);
                if(isDone == true) this->displayMarker();
            }

            case SimulationNodeClass::nodeType_thermalAnalysisTemperature:
            case SimulationNodeClass::nodeType_thermalAnalysisConvection:
            case SimulationNodeClass::nodeType_thermalAnalysisThermalFlux:
            case SimulationNodeClass::nodeType_thermalAnalysisThermalFlow:
            case SimulationNodeClass::nodeType_thermalAnalysisThermalPower:
            {
                emit requestHideAllResults();
                emit requestUnhighlightBodies(true);
                emit requestHideMeshes();
                emit requestHideSlicedMeshes();
                emit requestSetWorkingMode(2);
                this->changeColor();

                //! --------------
                //! set the model
                //! --------------
                QModelIndex index_analysisSettings = this->getAnalysisSettingsItemFromCurrentItem()->index();
                emit requestTabularData(index_analysisSettings);

                //! ---------------------------------
                //! show the first row with Time = 0
                //! ---------------------------------
                emit requestShowFirstRow();

                //! -----------------------------------------------------------
                //! calculate the number of columns to show => in the table <=
                //! -----------------------------------------------------------
                QList<int> columnsToShow;
                columnsToShow << TABULAR_DATA_STEP_NUMBER_COLUMN << TABULAR_DATA_STEP_END_TIME_COLUMN << mainTreeTools::getColumnsToRead(myTreeView);
                if(columnsToShow.length()>=2)
                {
                    emit requestShowColumns(columnsToShow);

                    //! ------------------------------------
                    //! remove the column showing the times
                    //! ------------------------------------
                    columnsToShow.removeFirst();
                    CustomTableModel *tabData = index_analysisSettings.data(Qt::UserRole).value<SimulationNodeClass*>()->getTabularDataModel();
                    emit requestShowGraph(tabData,columnsToShow);
                }
                bool isDone = markerBuilder::addMarker(this->getCurrentNode(), mySimulationDataBase);
                if(isDone == true) this->displayMarker();
            }
                break;

            case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_CompressionOnlySupport:
            case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_CylindricalSupport:
            case SimulationNodeClass::nodeType_structuralAnalysisBoundaryContidion_FixedSupport:
            case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_FrictionlessSupport:
            {
                emit requestHideAllResults();
                emit requestUnhighlightBodies(true);
                emit requestHideMeshes();
                emit requestHideSlicedMeshes();
                emit requestSetWorkingMode(2);
                this->changeColor();
                emit requestClearGraph();

                //! --------------
                //! set the model
                //! --------------
                QModelIndex index_analysisSettings = this->getAnalysisSettingsItemFromCurrentItem()->index();
                emit requestTabularData(index_analysisSettings);
                emit requestHideFirstRow();
                emit requestClearGraph();

                bool isDone = markerBuilder::addMarker(this->getCurrentNode(), mySimulationDataBase);
                if(isDone == true) this->displayMarker();
            }
                break;

            case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Displacement:
            case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteDisplacement:
            case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteRotation:
            case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RotationalVelocity:
            case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Force:
            case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteForce:
            case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Moment:
            case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Pressure:
            case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Acceleration:
            case SimulationNodeClass::nodeType_structuralAnalysisThermalCondition:
            case SimulationNodeClass::nodeType_modelChange:
            {
                emit requestHideAllResults();
                emit requestUnhighlightBodies(true);
                emit requestHideMeshes();
                emit requestHideSlicedMeshes();
                emit requestSetWorkingMode(2);
                this->changeColor();

                //! --------------
                //! set the model
                //! --------------
                QModelIndex index_analysisSettings = mainTreeTools::getAnalysisSettingsItemFromCurrentItem(myTreeView)->index();
                emit requestTabularData(index_analysisSettings);

                //! ---------------------------------------------------------------------
                //! show the first row with Time = 0, apart from the item "Model change"
                //! ---------------------------------------------------------------------
                if(theNodeType==SimulationNodeClass::nodeType_modelChange) emit requestHideFirstRow();
                else emit requestShowFirstRow();

                //! -----------------------------------------------------------
                //! calculate the number of columns to show => in the table <=
                //! -----------------------------------------------------------
                QList<int> columnsToShow;
                columnsToShow << TABULAR_DATA_STEP_NUMBER_COLUMN << TABULAR_DATA_STEP_END_TIME_COLUMN << mainTreeTools::getColumnsToRead(myTreeView);
                if(columnsToShow.length()>=3)
                {
                    emit requestShowColumns(columnsToShow);

                    //! ------------------------------------
                    //! remove the column showing the times
                    //! ------------------------------------
                    columnsToShow.removeFirst();
                    CustomTableModel *tabData = index_analysisSettings.data(Qt::UserRole).value<SimulationNodeClass*>()->getTabularDataModel();
                    //cout<<"____"<<tabData->rowCount()<<", "<<tabData->columnCount()<<"____"<<endl;
                    //for(int k=0; k<columnsToShow.length(); k++) cout<<"____column: "<<columnsToShow[k]<<"____"<<endl;
                    emit requestShowGraph(tabData,columnsToShow);
                }

                bool isDone = markerBuilder::addMarker(this->getCurrentNode(), mySimulationDataBase);
                if(isDone == true) this->displayMarker();
            }
                break;

            case SimulationNodeClass::nodeType_structuralAnalysisBoltPretension:
            {
                emit requestHideAllResults();
                emit requestUnhighlightBodies(Standard_True);
                emit requestHideMeshes();
                emit requestHideSlicedMeshes();
                emit requestSetWorkingMode(2);
                this->changeColor();

                //! set the model
                QModelIndex index_analysisSettings = this->getAnalysisSettingsItemFromCurrentItem()->index();

                emit requestTabularData(index_analysisSettings);

                //! hide the first row with Time = 0
                emit requestHideFirstRow();

                //! -----------------------------------------------------------
                //! calculate the number of columns to show => in the table <=
                //! -----------------------------------------------------------
                QList<int> columnsToShow;
                columnsToShow << TABULAR_DATA_STEP_NUMBER_COLUMN << TABULAR_DATA_STEP_END_TIME_COLUMN << mainTreeTools::getColumnsToRead(myTreeView);
                if(columnsToShow.length()>=3)
                {
                    emit requestShowColumns(columnsToShow);

                    //! ------------------------------------
                    //! remove the column showing the times
                    //! ------------------------------------
                    columnsToShow.removeFirst();
                    CustomTableModel *tabData = index_analysisSettings.data(Qt::UserRole).value<SimulationNodeClass*>()->getTabularDataModel();
                    emit requestShowGraph(tabData,columnsToShow);
                }

                bool isDone = markerBuilder::addMarker(this->getCurrentNode(), mySimulationDataBase);
                if(isDone == true) this->displayMarker();
            }
                break;

            case SimulationNodeClass::nodeType_connection:
            {
                emit requestHideAllResults();
                emit requestUnhighlightBodies(Standard_True);
                emit requestSetWorkingMode(1);
                emit requestHideSlicedMeshes();
                this->changeColor();
                emit requestClearGraph();
            }
                break;

            case SimulationNodeClass::nodeType_connectionGroup:
            {
                emit requestHideAllResults();
                emit requestUnhighlightBodies(Standard_True);
                emit requestHideMeshes();
                emit requestSetWorkingMode(1);
                emit requestHideSlicedMeshes();
                this->changeColor();
                emit requestClearGraph();
            }
                break;

            case SimulationNodeClass::nodeType_connectionPair:
            {
                emit requestHideAllResults();
                emit requestSetWorkingMode(1);
                emit requestHideSlicedMeshes();
                //this->changeColor();
                emit requestClearGraph();
                //! -----------------------------------------
                //! show the double view port - experimental
                //! -----------------------------------------
                emit requestShowDoubleViewPort(true);

                QVector<GeometryTag> vecLocM = theNode->getPropertyItem("Tags master")->data(Qt::UserRole).value<Property>().getData().value<QVector<GeometryTag>>();
                QVector<GeometryTag> vecLocS = theNode->getPropertyItem("Tags slave")->data(Qt::UserRole).value<Property>().getData().value<QVector<GeometryTag>>();

                std::vector<int> indexes_master;
                std::vector<int> indexes_slave;
                for(QVector<GeometryTag>::iterator it = vecLocM.begin(); it!=vecLocM.end(); ++it)
                {
                    GeometryTag loc = *it;
                    indexes_master.push_back(loc.parentShapeNr);
                }
                for(QVector<GeometryTag>::iterator it = vecLocS.begin(); it!=vecLocS.end(); ++it)
                {
                    GeometryTag loc = *it;
                    indexes_slave.push_back(loc.parentShapeNr);
                }

                //! -------------------------------------------
                //! hide all the bodies from the two viewports
                //! -------------------------------------------
                //TopTools_IndexedMapOfShape shapeMap = mySimulationDataBase->bodyMap;
                QMap<int,TopoDS_Shape> shapeMap = mySimulationDataBase->bodyMap;
                TColStd_ListOfInteger allBodies;

                for(QMap<int,TopoDS_Shape>::iterator it = shapeMap.begin(); it!=shapeMap.end(); ++it)
                {
                    int bodyIndex = it.key();
                    allBodies.Append(bodyIndex);
                }
                emit requestHideBodiesFromViewers(allBodies);

                //! -----------------------------------------------------------
                //! indexes of the bodies interested by the contact definition
                //! -----------------------------------------------------------
                indexes_master = tools::clearFromDuplicates(indexes_master);
                indexes_slave = tools::clearFromDuplicates(indexes_slave);

                TColStd_ListOfInteger masterListNbs, slaveListNbs;
                for(std::vector<int>::iterator it = indexes_master.begin(); it!=indexes_master.end(); ++it)
                {
                    //!cout<<"___adding master body number: "<<*it<<endl;
                    masterListNbs.Append(*it);
                }
                for(std::vector<int>::iterator it = indexes_slave.begin(); it!=indexes_slave.end(); ++it)
                {
                    //!cout<<"___adding slave body number: "<<*it<<endl;
                    slaveListNbs.Append(*it);
                }

                emit requestShowBodyOnMasterViewPort(masterListNbs, slaveListNbs);
                emit requestShowBodyOnSlaveViewPort(slaveListNbs, masterListNbs);
                this->changeColor();
            }
                break;

            case SimulationNodeClass::nodeType_namedSelection:
            {
                cout<<"SimulationManager::highlighter()->____named selection root item clicked____"<<endl;

                emit requestHideAllResults();
                emit requestUnhighlightBodies(Standard_True);
                emit requestHideMeshes();
                emit requestHideSlicedMeshes();
                emit requestSetWorkingMode(2);
                this->changeColor();
                emit requestClearGraph();
            }
                break;

            case SimulationNodeClass::nodeType_namedSelectionGeometry:
            //case SimulationNodeClass::nodeType_namedSelectionElement:
            case SimulationNodeClass::nodeType_namedSelectionNodal:
            {
                emit requestHideAllResults();
                emit requestUnhighlightBodies(Standard_True);
                emit requestHideMeshes();
                emit requestHideSlicedMeshes();
                emit requestSetWorkingMode(2);
                this->changeColor();
                emit requestClearGraph();
            }
                break;

            case SimulationNodeClass::nodeType_namedSelectionElement:
            {
                emit requestHideAllResults();
                emit requestUnhighlightBodies(Standard_True);
                emit requestSetWorkingMode(0);  //! "On mesh" working mode
                emit requestHideMeshes();
                emit requestHideSlicedMeshes();
                emit requestClearGraph();

                QStandardItem *item = theNode->getPropertyItem("Selected elements");
                if(item!=Q_NULLPTR)
                {
                    const occHandle(Ng_MeshVS_DataSource3D) &miniMesh = item->data(Qt::UserRole).value<Property>().getData().value<occHandle(Ng_MeshVS_DataSource3D)>();
                    this->displayFaceMesh(miniMesh,Quantity_NOC_AQUAMARINE1,true);
                }
                emit requestClearGraph();
            }
                break;
            }
        }
    }
    else
    {
        //! --------------------------------------
        //! the user has deselected all the items
        //! --------------------------------------
        emit requestUnhighlightBodies(true);
    }
}

//! ---------------------
//! function: setContext
//! details:
//! ---------------------
void SimulationManager::setContext(const occHandle(AIS_InteractiveContext) &aCTX)
{
    if(!aCTX.IsNull())
    {
        myCTX = aCTX;
        theTextWriter = new writeLabelClass(aCTX,this);
    }
    else
    {
        QMessageBox::critical(this, tr("SimulationManager::setContext()"),tr("Error: the context is null"));
    }
}

//! -----------------------------------------
//! function: createContent
//! details:  set the layout and the widgets
//! -----------------------------------------
#include <QScrollBar>
#include "simulationmanagerdelegate.h"
void SimulationManager::createContent()
{
    cout<<"SimulationManager::createContent()->____function called____"<<endl;
    QVBoxLayout *mainLayout = new QVBoxLayout();
    mainLayout->setMargin(0);
    //!mainLayout->setContentsMargins(QMargins());
    //!mainLayout->addWidget(new QSizeGrip(this), 0, Qt::AlignBottom | Qt::AlignRight);

    //! --------------
    //! the tree view
    //! --------------
    myTreeView = new QTreeView(this);

    //! ----------------------
    //! horizontal scrollbars
    //! ----------------------
    myTreeView->horizontalScrollBar()->setEnabled(true);
    myTreeView->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOn);
    myTreeView->header()->setStretchLastSection(false);
    //connect(myTreeView,SIGNAL(viewportEntered()),this,SLOT(onMouseEnter()));
    myTreeView->header()->setSectionResizeMode(QHeaderView::ResizeToContents);

    //! delegate
    //SimulationManagerDelegate *delegate = new SimulationManagerDelegate(this);
    //myTreeView->setItemDelegate(delegate);

    //! --------------------------------------------------------------------
    //! allow multiple selection
    //!
    //! Qt documentation
    //! When the user selects an item in the usual way, the selection
    //! is cleared and the new item selected. However, if the user
    //! presses the Ctrl key when clicking on an item, the clicked
    //! item gets toggled and all other items are left untouched.
    //! If the user presses the Shift key while clicking on an item,
    //! all items between the current item and the clicked item are
    //! selected or unselected, depending on the state of the clicked item.
    //! Multiple items can be selected by dragging the mouse over them.
    //! --------------------------------------------------------------------
    myTreeView->setSelectionMode(QAbstractItemView::ExtendedSelection);

    mySelectionModel = myTreeView->selectionModel();
    myTreeView->header()->hide();

    QString style = "QTreeView::branch:has-siblings:!adjoins-item{border-image:"
                    "url(:/stylesheets/vline.png) 0;}"
                    "QTreeView::bra"
                    "ch:has-siblings:adjoins-item {border-image:"
                    "url(:/stylesheets/branch-more.png) 0;}"
                    "QTreeView::branch:!has-children:!has-siblings:adjoins-item {border-image:"
                    "url(:/stylesheets/branch-end.png) 0;}"
                    "QTreeView::branch:has-children:!has-siblings:closed,"
                    "QTreeView::branch:closed:has-children:has-siblings {border-image: none;image:"
                    "url(:/stylesheets/branch-closed.png);}"
                    "QTreeView::branch:open:has-children:!has-siblings,"
                    "QTreeView::branch:open:has-children:has-siblings {border-image: none;image:"
                    "url(:/stylesheets/branch-open.png);}";
                    //"QTreeView::indicator:checked{image: url(:/stylesheets/icon_checkmark green.png);}"
                    //"QTreeView::indicator:!checked{image: url(:/stylesheets/icon_checkmark yellow.png);}";


    //QTreeView::hover {text-decoration: underline;}
    myTreeView->setStyleSheet(style);
    mainLayout->addWidget(myTreeView);
    setLayout(mainLayout);

    //! ----------------
    //! create the menu
    //! ----------------
    myContextMenu = new QMenu(this);

    //! ----------------------------
    //! set the context menu policy
    //! ----------------------------
    this->setContextMenuPolicy(Qt::CustomContextMenu);
    connect(this,SIGNAL(customContextMenuRequested(const QPoint&)),this,SLOT(showContextMenu(const QPoint&)));
}

//! --------------------------
//! function: buildCustomMenu
//! details:
//! --------------------------
void SimulationManager::buildCustomMenu(const QModelIndex &modelIndex)
{
    //! ------------------------------
    //! check if something is running
    //! ------------------------------
    bool isEnabled = this->isSomethingRunning()? false:true;
    myContextMenu->clear();

    SimulationNodeClass *node = modelIndex.data(Qt::UserRole).value<SimulationNodeClass*>();
    SimulationNodeClass::nodeType nodeType = node->getType();
    if(modelIndex.isValid())
    {
        switch(node->getFamily())
        {
        case SimulationNodeClass::nodeType_root: contextMenuBuilder::buildModelRootContextMenu(myContextMenu); break;
        case SimulationNodeClass::nodeType_import: contextMenuBuilder::buildImportContextMenu(myContextMenu,false,isEnabled); break;
        case SimulationNodeClass::nodeType_geometry:
        {
            contextMenuBuilder::buildGeometryContextMenu(myContextMenu,true,isEnabled);
        }
            break;
        case SimulationNodeClass::nodeType_coordinateSystems: contextMenuBuilder::buildCoordinateSystemMenu(myContextMenu); break;
        case SimulationNodeClass::nodeType_remotePointRoot: contextMenuBuilder::buildRemotePointContextMenu(myContextMenu); break;
        case SimulationNodeClass::nodeType_namedSelection: contextMenuBuilder::buildNamedSelectionContextMenu(myContextMenu); break;

        case SimulationNodeClass::nodeType_connection:
        case SimulationNodeClass::nodeType_connectionGroup:
            contextMenuBuilder::buildContactContextMenu(myContextMenu);
            break;

        case SimulationNodeClass::nodeType_meshControl:
            contextMenuBuilder::buildMeshContextMenu(myContextMenu,true,isEnabled);
            break;

        case SimulationNodeClass::nodeType_structuralAnalysis:
            contextMenuBuilder::buildStructuralAnalysisContextMenu(myContextMenu,true,isEnabled);
            if(nodeType==SimulationNodeClass::nodeType_OpenFoamScalarData) break;
            if(nodeType==SimulationNodeClass::nodeType_mapper) break;
            break;

        case SimulationNodeClass::nodeType_combinedAnalysis:
            contextMenuBuilder::buildStructuralAnalysisContextMenu(myContextMenu,false,isEnabled);
            if(nodeType==SimulationNodeClass::nodeType_OpenFoamScalarData) break;
            if(nodeType==SimulationNodeClass::nodeType_mapper) break;
            contextMenuBuilder::buildThermalAnalysisContextMenu(myContextMenu,false,isEnabled);
            break;

        case SimulationNodeClass::nodeType_particlesInFieldsAnalysis:
            contextMenuBuilder::buildParticlesInFieldsContextMenu(myContextMenu,true,isEnabled);
            break;

        case SimulationNodeClass::nodeType_thermalAnalysis:
            contextMenuBuilder::buildThermalAnalysisContextMenu(myContextMenu,true,isEnabled);
            break;

        case SimulationNodeClass::nodeType_StructuralAnalysisSolution:
            contextMenuBuilder::buildStructuralSolutionContextMenu(myContextMenu);
            break;

        case SimulationNodeClass::nodeType_thermalAnalysisSolution: contextMenuBuilder::buildThermalResultsContextMenu(myContextMenu); break;
        case SimulationNodeClass::nodeType_postObject: contextMenuBuilder::buildPostObjectContextMenu(myContextMenu); break;
        case SimulationNodeClass::nodeType_combinedAnalysisSolution: contextMenuBuilder::buildCombinedAnalysisResultsContextMenu(myContextMenu,true,isEnabled); break;
        }
    }
}

//! --------------------------------
//! function: show the context menu
//! details:
//! --------------------------------
void SimulationManager::showContextMenu(const QPoint &pos)
{    
    const QPoint &globalPos = this->mapToGlobal(pos);

    //! build the context menu
    const QModelIndex &index = myTreeView->indexAt(pos);
    if(index.isValid()) this->buildCustomMenu(index);

    //! handling the selected action
    QAction* selectedAction = myContextMenu->exec(globalPos);

    //! check if the item is valid
    if(selectedAction)  handleItem(selectedAction->data().toInt());
}

//! ------------------------------------------
//! function: deleteItem
//! details:  delete items from the main tree
//! ------------------------------------------
void SimulationManager::deleteItem(QList<QModelIndex> indexesList)
{
    cout<<"SimulationManager::deleteItem()->____function called____"<<endl;
    QList<QModelIndex> itemListToDelete;

    //! -------------------------------------
    //! list of items that cannot be deleted
    //! -------------------------------------
    QList<SimulationNodeClass::nodeType> undeletableItems;
    undeletableItems<<//SimulationNodeClass::nodeType_structuralAnalysis<<
                      SimulationNodeClass::nodeType_structuralAnalysisSettings<<
                      SimulationNodeClass::nodeType_StructuralAnalysisSolution<<
                      SimulationNodeClass::nodeType_StructuralAnalysisSolutionInformation<<
                      SimulationNodeClass::nodeType_thermalAnalysisSettings<<
                      //SimulationNodeClass::nodeType_thermalAnalysis<<
                      SimulationNodeClass::nodeType_thermalAnalysisSolution<<
                      SimulationNodeClass::nodeType_thermalAnalysisSolutionInformation<<
                      SimulationNodeClass::nodeType_geometry<<
                      SimulationNodeClass::nodeType_geometryBody<<
                      SimulationNodeClass::nodeType_connection<<
                      SimulationNodeClass::nodeType_meshControl<<
                      SimulationNodeClass::nodeType_namedSelection<<
                      SimulationNodeClass::nodeType_coordinateSystems<<
                      SimulationNodeClass::nodeType_coordinateSystem_global;

    //! ------------------------------------------------------------------------------
    //! "Displacement" "Remote displacement" "Remote rotation" have the option "free"
    //! ------------------------------------------------------------------------------
    QList<SimulationNodeClass::nodeType> specialItems;
    specialItems<<SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Displacement<<
                  SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteDisplacement<<
                  SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteRotation;

    specialItems<<SimulationNodeClass::nodeType_namedSelectionGeometry<<
                  SimulationNodeClass::nodeType_structuralAnalysisThermalCondition<<
                  SimulationNodeClass::nodeType_structuralAnalysisBoltPretension<<
                  SimulationNodeClass::nodeType_thermalAnalysisConvection<<
                  SimulationNodeClass::nodeType_thermalAnalysisRadiation<<
                  SimulationNodeClass::nodeType_thermalAnalysisTemperature<<
                  SimulationNodeClass::nodeType_thermalAnalysisThermalFlow<<
                  SimulationNodeClass::nodeType_thermalAnalysisThermalFlux<<
                  SimulationNodeClass::nodeType_thermalAnalysisThermalPower;

    //! --------------------------------------------------------------------
    //! list of the selected items
    //! if the list passed via argument is empty, search into the selection
    //! --------------------------------------------------------------------
    if(indexesList.isEmpty()) indexesList<<myTreeView->selectionModel()->selectedIndexes();
    if(indexesList.isEmpty()) return;   //! nothing to delete

    cout<<"SimulationManager::deleteItem()->____number of selected items for deletion: "<<indexesList.length()<<"____"<<endl;

    //! -------------------------------------------------
    //! check if the selection contains an analysis root
    //! -------------------------------------------------
    bool allSelectedAnalysisRoot = false;
    for(int i=0; i<indexesList.length(); i++)
    {
        const QModelIndex &modelIndex = indexesList.at(i);
        SimulationNodeClass *aNode = modelIndex.data(Qt::UserRole).value<SimulationNodeClass*>();
        if(aNode->isAnalysisRoot()==true)
        {
            allSelectedAnalysisRoot = true;
            continue;
        }
        else allSelectedAnalysisRoot = false;
    }
    if(allSelectedAnalysisRoot)
    {
        int val = QMessageBox::warning(this,APPNAME,"This will delete analysis setup and results.\n Are you sure?",QMessageBox::Ok, QMessageBox::Cancel);
        if((QMessageBox::StandardButton)(val)==QMessageBox::Cancel) return;
        cout<<"____deleting____"<<endl;
        return;
    }

    //! -------------------------------------------
    //! retrieve the analysis settings item
    //! (here used in some of the following cases)
    //! -------------------------------------------
    SimulationNodeClass* nodeAnalysisSetting = myTreeView->currentIndex().parent().child(0,0).data(Qt::UserRole).value<SimulationNodeClass*>();
    CustomTableModel *tabData = nodeAnalysisSetting->getTabularDataModel();

    for(int i=0; i<indexesList.length(); i++)
    {
        const QModelIndex &modelIndex = indexesList.at(i);
        SimulationNodeClass *aNode = modelIndex.data(Qt::UserRole).value<SimulationNodeClass*>();
        SimulationNodeClass::nodeType theType = aNode->getType();

        cout<<"SimulationManager::deleteItem()->____trying to delete item: \""<<aNode->type().toStdString()<<"\"____"<<endl;

        //! -----------------------------------------------------------------------------
        //! items defined by a one or three components having the "Define by" components
        //! -----------------------------------------------------------------------------
        if(!undeletableItems.contains(theType) && !specialItems.contains(theType))
        {
            //! ---------------------------------------------------------
            //! add to the list for final removal fro from the tree view
            //! ---------------------------------------------------------
            itemListToDelete<<modelIndex;

            //! --------------------------------------------
            //! items/nodes having the "Define by" property
            //! they require an update of the tabular data
            //! --------------------------------------------
            if(aNode->getPropertyItem("Define by")!=Q_NULLPTR)
            {
                int SC = mainTreeTools::calculateStartColumn(myTreeView);
                Property::defineBy theDefineBy = aNode->getPropertyValue<Property::defineBy>("Define by");
                int count;
                switch(theDefineBy)
                {
                case Property::defineBy_normal:
                case Property::defineBy_vector:
                    count = 1;
                    break;
                case Property::defineBy_components:
                    count = 3;
                    break;
                }
                tabData->removeColumns(SC,count);
            }            
        }

        //! ------------------------------------------------------------------------------------
        //! handle the special case of a "Displacement"/"Remote displacement"/"Remote rotation"
        //! ------------------------------------------------------------------------------------
        if(theType==SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Displacement ||
                theType ==SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteDisplacement ||
                theType ==SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteRotation)
        {
            int SC = mainTreeTools::calculateStartColumn(myTreeView);
            int count = 0;
            Property::defineBy theDefineBy = aNode->getPropertyValue<Property::defineBy>("Define by");
            switch(theDefineBy)
            {
            case Property::defineBy_vector:
            case Property::defineBy_normal:
                count = 1; break;
            case Property::defineBy_components:
            {
                Property::loadDefinition loadDef_X = aNode->getPropertyValue<Property::loadDefinition>("X component");
                Property::loadDefinition loadDef_Y = aNode->getPropertyValue<Property::loadDefinition>("Y component");
                Property::loadDefinition loadDef_Z = aNode->getPropertyValue<Property::loadDefinition>("Z component");
                if(loadDef_X!=Property::loadDefinition_free) count++;
                if(loadDef_Y!=Property::loadDefinition_free) count++;
                if(loadDef_Z!=Property::loadDefinition_free) count++;
            }
                break;
            }
            tabData->removeColumns(SC,count);
            itemListToDelete<<modelIndex;
        }

        //! ------------------------------------------------
        //! handle the special case of a "Named selection"
        //! -----------------------------------------------
        if(theType==SimulationNodeClass::nodeType_namedSelectionGeometry)
        {
            itemListToDelete<<modelIndex;

            //! ---------------------------------
            //! parse the simulation setup items
            //! ---------------------------------
            QStandardItem *modelRoot = Geometry_RootItem->parent();
            for(int n=0; n<modelRoot->rowCount(); n++)
            {
                QStandardItem *analysisRoot = modelRoot->child(n,0);
                SimulationNodeClass *nodeAnalysisRoot = analysisRoot->data(Qt::UserRole).value<SimulationNodeClass*>();
                if(nodeAnalysisRoot->isAnalysisRoot()==false) continue;

                //! --------------------------------------------------------------
                //! scan the simulation setup nodes of tthe current analysis root
                //! --------------------------------------------------------------
                for(int i=1; i<analysisRoot->rowCount()-1; i++)
                {
                    QStandardItem *setUpItem = analysisRoot->child(i,0);
                    SimulationNodeClass *nodeSetUp = setUpItem->data(Qt::UserRole).value<SimulationNodeClass*>();
                    if(nodeSetUp->isSimulationSetUpNode()== false) return;

                    //! --------------------------------------------------------------------
                    //! if a named selection is contained replace with a geometry selection
                    //! and change the scoping method accordingly
                    //! --------------------------------------------------------------------
                    if(nodeSetUp->getPropertyItem("Named selection")==Q_NULLPTR) continue;
                    {
                        nodeSetUp->getModel()->blockSignals(true);

                        QVariant data;
                        data.setValue(QVector<GeometryTag>());
                        nodeSetUp->replaceProperty("Geometry",Property("Tags",data,Property::Property::PropertyGroup_Scope));
                        nodeSetUp->replaceProperty("Tags",Property("Tags",data,Property::Property::PropertyGroup_Scope));
                        data.setValue(Property::ScopingMethod_GeometrySelection);
                        nodeSetUp->replaceProperty("Scoping method",Property("Scoping method",data,Property::PropertyGroup_Scope));
                        nodeSetUp->getModel()->blockSignals(false);
                    }

                }
            }

        }

        //! -------------------------------------------------
        //! handle the special case of a "Thermal condition"
        //! handle the special case of a "Bolt pretension"
        //! -------------------------------------------------
        if(theType == SimulationNodeClass::nodeType_structuralAnalysisThermalCondition)
        {
            int SC = mainTreeTools::calculateStartColumn(myTreeView);
            int count = 1;
            tabData->removeColumns(SC,count);
            itemListToDelete<<modelIndex;
        }
        if(theType == SimulationNodeClass::nodeType_structuralAnalysisBoltPretension)
        {
            int SC = mainTreeTools::calculateStartColumn(myTreeView);
            int count = 3;
            tabData->removeColumns(SC,count);
            itemListToDelete<<modelIndex;
        }

        //! --------------------------------------------------------------
        //! handle the thermal analysis items having tabular data defined
        //! --------------------------------------------------------------
        if(theType == SimulationNodeClass::nodeType_thermalAnalysisThermalFlux ||
                theType == SimulationNodeClass::nodeType_thermalAnalysisRadiation ||
                theType == SimulationNodeClass::nodeType_thermalAnalysisTemperature ||
                theType == SimulationNodeClass::nodeType_thermalAnalysisThermalFlow)
        {
            int SC = mainTreeTools::calculateStartColumn(myTreeView);
            int count = 1;
            tabData->removeColumns(SC,count);
            itemListToDelete<<modelIndex;
        }
        if(theType == SimulationNodeClass::nodeType_thermalAnalysisConvection)
        {
            int SC = mainTreeTools::calculateStartColumn(myTreeView);
            int count = 2;
            tabData->removeColumns(SC,count);
            itemListToDelete<<modelIndex;
        }
    }

    //! ---------------------------------------------------------------------------
    //! delete the selected items from the model
    //! see Ref.
    //! https://www.qtcentre.org/threads/
    //! 32721-How-can-I-remove-a-list-of-selected-items-in-the-QListView-in-QT-4-6
    //! ---------------------------------------------------------------------------
    myTreeView->setUpdatesEnabled(false);
    QList<QModelIndex> indexes; indexes<<itemListToDelete;
    std::sort(indexes.begin(), indexes.end());
    for (int i = indexes.count() - 1; i > -1; --i)
        myModel->removeRow(indexes.at(i).row(),indexes.at(i).parent());
    myTreeView->setUpdatesEnabled(true);
}

//! ---------------------
//! function: handleItem
//! details:
//! ---------------------
#include "maintreetools.h"
#include <exportingtools.h>
#include "bolttool.h"
#include <TopLoc_Datum3D.hxx>
#include "shapecomparison.h"
#include <GProp_GProps.hxx>

void SimulationManager::handleItem(int type)
{
    cout<<"SimulationManager::handleItem()->____function called: action Nr: "<<type<<"_____"<<endl;
    switch(type)
    {
    case 78: this->createSimulationNode(SimulationNodeClass::nodeType_pointMass); break;
    case 85:
    {
        SimulationNodeClass *curNode = myTreeView->currentIndex().data(Qt::UserRole).value<SimulationNodeClass*>();
        const QVector<GeometryTag> tags = curNode->getPropertyValue<QVector<GeometryTag>>("Geometry");
        if(tags.isEmpty() || tags.size()>1) return;
        TopAbs_ShapeEnum shapeType = tags.first().subShapeType;
        if(shapeType!=TopAbs_SOLID) return;

        int bodyIndex = tags.at(0).parentShapeNr;
        const occHandle(Ng_MeshVS_DataSource3D) &volumeMesh = occHandle(Ng_MeshVS_DataSource3D)::DownCast(mySimulationDataBase->ArrayOfMeshDS.value(bodyIndex));
        if(volumeMesh.IsNull()) return;

        //! --------------------------------------
        //! retrieve the "bolt" coordinate system
        //! (origin, direction)
        //! --------------------------------------
        void *p = curNode->getPropertyValue<void*>("Coordinate system");
        QStandardItem *itemCS = static_cast<QStandardItem*>(p);
        SimulationNodeClass *nodeCS = itemCS->data(Qt::UserRole).value<SimulationNodeClass*>();
        QVector<double> boltDirection;
        QVector<double> planeOrigin;
        if(nodeCS->getType()==SimulationNodeClass::nodeType_coordinateSystem_global)
        {
            planeOrigin.push_back(0); planeOrigin.push_back(0); planeOrigin.push_back(0);
            boltDirection = nodeCS->getPropertyValue<QVector<double>>("Z axis data");
        }
        else
        {
            planeOrigin = nodeCS->getPropertyValue<QVector<double>>("Base origin");
            QVector<QVector<double>> baseDirData = nodeCS->getPropertyValue<QVector<QVector<double>>>("Base directional data");
            boltDirection = baseDirData.at(2);
        }

        double xP = planeOrigin[0]; double yP = planeOrigin[1]; double zP = planeOrigin[2];
        double a = boltDirection[0]; double b = boltDirection[1]; double c = boltDirection[2];
        double d = -a*xP-b*yP-c*zP;

        boltTool bt(volumeMesh);
        //std::vector<std::pair<int,std::string>> vecCCXFaceDefs;
        std::vector<std::pair<int,int>> vecCCXFaceDefs;

        occHandle(MeshVS_DataSource) slicedMeshDS;
        bool isDone = bt.sliceMeshWithPlane(a,b,c,d,slicedMeshDS,vecCCXFaceDefs);

        //! -----------
        //! diagnostic
        //! -----------
        const occHandle(Ng_MeshVS_DataSourceFace) &m = occHandle(Ng_MeshVS_DataSourceFace)::DownCast(slicedMeshDS);
        FILE *f_slicedMesh = fopen("D:/slicedMesh.txt","w");
        for(int n=1; n<m->GetAllNodes().Extent(); n++)
        {
            int localNodeID = n;
            int globalNodeID = m->myNodesMap.FindKey(localNodeID);
            const std::vector<double> &p = m->getNodeCoordinates(localNodeID);
            double x = p.data()[0];
            double y = p.data()[1];
            double z = p.data()[2];
            fprintf(f_slicedMesh,"%d\t%d\t%lf\t%lf\t%lf\n",localNodeID,globalNodeID,x,y,z);
        }
        fclose(f_slicedMesh);
        //! ---------------
        //! end diagnostic
        //! ---------------

        if(isDone)
        {
            emit requestHideMeshes();
            displayFaceMesh(slicedMeshDS,Quantity_NOC_AQUAMARINE1,true);
        }
    }
        break;

    case 86: this->replicateBolt(); break;

    case 26:
    {
        //! --------------------------------------------------------------------------------------------------------
        //! "Promote to remote point": active when a "Remote force" or a "Remote displacement" is defined through
        //! a direct geometry selection or a named selection. When the action is triggered the reference point
        //! P(xR, yR, zR) (internally calculated as the centroid of the geometry of or the named selection) is used
        //! for creating a "Remote point" item under the "Remote points" tree branch. This remote point is placed
        //! in P(xR, yR, zR), and it is attached to the "Remote force" or "Remote displacement" geometry selection.
        //! At the same time the Scoping method of the "Remote force" of the "Remote displacement" is changed
        //! to "Remote point", the "Remote point" pointing to the just created one
        //! --------------------------------------------------------------------------------------------------------
        if(this->getTreeItem(SimulationNodeClass::nodeType_remotePointRoot)==Q_NULLPTR)
        {
            mySimulationDataBase->createRemotePointRoot();
            //! hide the dummy "Select from list" item
            myTreeView->setRowHidden(0,this->getTreeItem(SimulationNodeClass::nodeType_remotePointRoot)->index(),true);
        }
        QVariant options;
        QVector<double> refPoint = this->getCurrentNode()->getPropertyValue<QVector<double>>("Reference point");

        cout<<"---------------------------------------------------------"<<endl;
        cout<<" creating a new remote point using as reference point: ("<<refPoint.at(0)<<", "<<refPoint.at(1)<<", "<<refPoint.at(2)<<")____"<<endl;
        options.setValue(refPoint);

        //! -----------------------------------------------------------------------------------------------
        //! note: before creating the new remote point, record the current node, which is a "Remote force"
        //! or a "Remote displacement", since once created, the current item becomes the just created one
        //! (i.e. the "promoted" "Remote point")
        //! -----------------------------------------------------------------------------------------------
        SimulationNodeClass *curNode = this->getCurrentNode();
        curNode->getModel()->blockSignals(true);

        this->createSimulationNode(SimulationNodeClass::nodeType_remotePoint,options);
        QExtendedStandardItem *newItemRP = static_cast<QExtendedStandardItem*>(mySimulationDataBase->getModel()->itemFromIndex(myTreeView->currentIndex()));
        SimulationNodeClass *theNewRPNode = this->getCurrentNode();
        theNewRPNode->getModel()->blockSignals(true);

        QVariant data;

        //! ------------------------------------------------------------------------------
        //! record the "Scoping method", "Geometry"/"Named Selection", "Tags" of the
        //! "Remote force"/"Remote displacement" and replace the corresponding properties
        //! of the just created "Remote point"
        //! ------------------------------------------------------------------------------
        if(curNode->getPropertyItem("Scoping method")!=Q_NULLPTR)
        {
            cout<<" replacing the \"Scoping method\" method of the promoted remote point "<<endl;
            Property::ScopingMethod scopingMethod = curNode->getPropertyValue<Property::ScopingMethod>("Scoping method");
            data.setValue(scopingMethod);
            Property prop("Scoping method",data,Property::PropertyGroup_Scope);
            theNewRPNode->replaceProperty("Scoping method",prop);
            //curNode->removeProperty("Scoping method");
        }
        if(curNode->getPropertyItem("Geometry")!=Q_NULLPTR)
        {
            cout<<" replacing the \"Geometry\" property of the promoted remote point"<<endl;
            QVector<GeometryTag> vecLoc = curNode->getPropertyValue<QVector<GeometryTag>>("Geometry");
            data.setValue(vecLoc);
            Property prop("Geometry",data,Property::PropertyGroup_Scope);
            theNewRPNode->replaceProperty("Geometry",prop);
            //curNode->removeProperty("Geometry");
        }
        if(curNode->getPropertyItem("Named selection")!=Q_NULLPTR)
        {
            cout<<" replacing the \"Named selection\" property of the promoted remote point"<<endl;
            void *p = curNode->getPropertyItem("Named selection")->data(Qt::UserRole).value<Property>().getData().value<void*>();
            data.setValue(p);
            Property prop("Named selection",data,Property::PropertyGroup_Scope);

            if(theNewRPNode->getPropertyItem("Geometry")!=NULL)
            {
                theNewRPNode->removeProperty("Geometry");
                theNewRPNode->addProperty(prop,1);
            }
            if(theNewRPNode->getPropertyItem("Named selection")!=NULL)
            {
                //! use "replace" if the property "Named selection" exists
                theNewRPNode->replaceProperty("Named selection",prop);
            }
            //curNode->removeProperty("Named selection");
        }
        if(curNode->getPropertyItem("Tags")!=Q_NULLPTR)
        {
            cout<<" replacing the \"Tags\" property of the promoted remote point"<<endl;
            QVector<GeometryTag> vecLoc = curNode->getPropertyValue<QVector<GeometryTag>>("Tags");
            data.setValue(vecLoc);
            Property prop("Tags",data,Property::PropertyGroup_Scope);
            theNewRPNode->replaceProperty("Tags",prop);
            //curNode->removeProperty("Tags");
        }

        //! ----------------------------------------------------
        //! force the "new" scoping method through remote point
        //! ----------------------------------------------------
        data.setValue(Property::ScopingMethod_RemotePoint);
        Property prop_scopingMethod("Scoping method",data,Property::PropertyGroup_Scope);
        curNode->replaceProperty("Scoping method",prop_scopingMethod);

        //! -----------------------------------------
        //! force the "new" "Remote points" property
        //! under "Scoping method"
        //! -----------------------------------------
        curNode->removeProperty("Geometry");
        curNode->removeProperty("Named selection");
        void *p = (void*)newItemRP;
        data.setValue(p);
        Property prop_remotePoint("Remote points",data,Property::PropertyGroup_Scope);
        curNode->addProperty(prop_remotePoint,1);

        //! ----------------------------------------------------------------------------
        //! disable the "Advanced" properties of the "Remote force"/"Remote properties"
        //! (the control of the coupling type and of the DOFs coupling is transferred
        //! to the promoted remote point ... to do
        //! ----------------------------------------------------------------------------
        // can be done in "GeneralDelegate"

        //! -------------------------------------
        //! remove all the "Advanced" properties
        //! -------------------------------------
        curNode->getModel()->removeRows(curNode->getPropertyItem("Coupling")->index().row(),curNode->getPropertyItem("Coupling")->parent()->rowCount(),
                                        curNode->getPropertyItem("Coupling")->parent()->index());

        //! ---------------------------------------------------------------------------------
        //! replace the "Advanced" properties of the "Remote force" or "Remote displacement"
        //! with the "Advanced properties" of the remote point
        //! ---------------------------------------------------------------------------------
        SimulationNodeClass *nodeNewRP = newItemRP->data(Qt::UserRole).value<SimulationNodeClass*>();
        QStandardItem *itemSeparatorAdvanced = nodeNewRP->getPropertyItem("Coupling")->parent();
        for(int row = 0; row<itemSeparatorAdvanced->rowCount(); row++)
        {
            QString name = itemSeparatorAdvanced->child(row,0)->data(Qt::DisplayRole).toString();
            QVariant data = itemSeparatorAdvanced->child(row,1)->data(Qt::UserRole).value<Property>().getData();
            Property::PropertyGroup group = itemSeparatorAdvanced->child(row,1)->data(Qt::UserRole).value<Property>().getGroup();
            Property curProp(name,data,group);
            curNode->addProperty(curProp,row);
        }

        curNode->getModel()->blockSignals(false);
        newItemRP->data(Qt::UserRole).value<SimulationNodeClass*>()->getModel()->blockSignals(false);
    }
        break;
    case 27: this->createSimulationNode(SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteRotation); break;
    case 28: this->createSimulationNode(SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteDisplacement); break;
    case 29: this->createSimulationNode(SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteForce); break;
    case 31: this->createSimulationNode(SimulationNodeClass::nodeType_structuralAnalysisBoundaryContidion_FixedSupport); break;
    case 32:
        this->request2DBodySelectionMode(true);
        this->createSimulationNode(SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_CylindricalSupport);
        break;
    case 33:
        this->request2DBodySelectionMode(true);
        this->createSimulationNode(SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_FrictionlessSupport);
        break;
    case 108:
        this->request2DBodySelectionMode(true);
        this->createSimulationNode(SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_CompressionOnlySupport);
        break;
    case 34: this->createSimulationNode(SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Force); break;
    case 35: this->createSimulationNode(SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Moment); break;
    case 36:
        this->request2DBodySelectionMode(true);
        this->createSimulationNode(SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Pressure);
        break;
    case 37: this->createSimulationNode(SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Displacement); break;
    case 38: this->createSimulationNode(SimulationNodeClass::nodeType_structuralAnalysisThermalCondition); break;
    case 8:
        this->request3DBodySelectionMode(true);
        this->createSimulationNode(SimulationNodeClass::nodeType_structuralAnalysisThermalCondition);
        break;
    case 97:
        this->createSimulationNode(SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_ImportedTemperatureDistribution);
        break;
    case 30: this->createSimulationNode(SimulationNodeClass::nodeType_namedSelectionGeometry); break;
    case 44: this->createSimulationNode(SimulationNodeClass::nodeType_connectionGroup); break;
    case 45:
    {
        mySimulationDataBase->createRemotePointRoot();
        //! -----------------------------------------------------------
        //! immediately hide the empty remote point "Select from list"
        //! -----------------------------------------------------------
        //myTreeView->setRowHidden(0,this->getTreeItem(SimulationNodeClass::nodeType_remotePointRoot)->index(),true);
    }
        break;
    case 46: this->createSimulationNode(SimulationNodeClass::nodeType_remotePoint); break;
    case 47: this->createSimulationNode(SimulationNodeClass::nodeType_coordinateSystem); break;
    case 48: this->swapContact(); break;
    case 49:
    {
        //! -----------------------------------------------------
        //! since this action is triggered by the user, use
        //! connectionPairGenerationOptions with isManual = true
        //! -----------------------------------------------------
        connectionPairGenerationOption aConnectionPairGenerationOptions;
        aConnectionPairGenerationOptions.manual = true;
        QVariant addOptions;
        addOptions.setValue(aConnectionPairGenerationOptions);
        //! ------------------------------------
        //! actually create the connection item
        //! ------------------------------------
        this->createSimulationNode(SimulationNodeClass::nodeType_connectionPair,addOptions);
    }
        break;
    case 50:
    {
        QVector<GeometryTag> mergedMasterTags, mergedSlaveTags;
        QList<QModelIndex> insel = myTreeView->selectionModel()->selectedIndexes();

        //! ----------------------------
        //! record the scope properties
        //! ----------------------------
        for(QList<QModelIndex>::iterator it = insel.begin(); it!=insel.end(); it++)
        {
            QModelIndex index = *it;
            SimulationNodeClass *node = index.data(Qt::UserRole).value<SimulationNodeClass*>();
            mergedMasterTags.append(node->getPropertyValue<QVector<GeometryTag>>("Tags master"));
            mergedSlaveTags.append(node->getPropertyValue<QVector<GeometryTag>>("Tags slave"));
        }
        QVariant data;
        data.setValue(mergedMasterTags);
        Property prop_masterScope("Master",data,Property::PropertyGroup_Scope);
        Property prop_masterTags("Tags master",data,Property::PropertyGroup_Scope);

        data.setValue(mergedSlaveTags);
        Property prop_slaveScope("Slave",data,Property::PropertyGroup_Scope);
        Property prop_slaveTags("Tags slave",data,Property::PropertyGroup_Scope);

        this->deleteItem(insel);

        this->createSimulationNode(SimulationNodeClass::nodeType_connectionGroup);
        this->createSimulationNode(SimulationNodeClass::nodeType_connectionPair);

        SimulationNodeClass *node = myTreeView->currentIndex().data(Qt::UserRole).value<SimulationNodeClass*>();
        data.setValue(QString("Multiple to multiple"));
        myModel->itemFromIndex(myTreeView->currentIndex())->setData(data,Qt::DisplayRole);
        node->getModel()->blockSignals(true);
        node->setName(QString("Multiple to multiple"));
        node->replaceProperty("Master",prop_masterScope);
        node->replaceProperty("Tags master",prop_masterTags);
        node->replaceProperty("Slave",prop_slaveScope);
        node->replaceProperty("Tags slave",prop_slaveTags);
        node->getModel()->blockSignals(false);
    }
        break;

    case 51: this->createSimulationNode(SimulationNodeClass::nodeType_meshBodyMeshControl); break;
    case 72: this->createSimulationNode(SimulationNodeClass::nodeType_meshFaceSize); break;
    case 74: this->createSimulationNode(SimulationNodeClass::nodeType_meshMethod); break;
    case 75: this->createSimulationNode(SimulationNodeClass::nodeType_meshPrismaticLayer); break;
    case 79: this->createSimulationNode(SimulationNodeClass::nodeType_meshMeshType); break;
    case 71: this->createSimulationNode(SimulationNodeClass::nodeType_mapper); break;
    case 73: this->createSimulationNode(SimulationNodeClass::nodeType_structuralAnalysisBoltPretension); break;
    case 53: this->createSimulationNode(SimulationNodeClass::nodeType_importedBodyScalar); break;
    case 100: this->createSimulationNode(SimulationNodeClass::nodeType_OpenFoamScalarData); break;
    case 104:
    {
        mySimulationDataBase->createThermalAnalysisRootNode();
        //! ------------------------------------------
        //! make the last inserted branch the current
        //! ------------------------------------------
        QStandardItem *itemRoot = Geometry_RootItem->parent();
        int Nb = itemRoot->rowCount();
        myTreeView->selectionModel()->setCurrentIndex(myModel->indexFromItem(itemRoot->child(Nb-1,0)),QItemSelectionModel::Current);
        this->setTheActiveAnalysisBranch();
    }
        break;
    case 105:
    {
        mySimulationDataBase->createStructuralAnalysisRootNode();
        //! ------------------------------------------
        //! make the last inserted branch the current
        //! ------------------------------------------
        QStandardItem *itemRoot = Geometry_RootItem->parent();
        int Nb = itemRoot->rowCount();
        myTreeView->selectionModel()->setCurrentIndex(myModel->indexFromItem(itemRoot->child(Nb-1,0)),QItemSelectionModel::Current);
        this->setTheActiveAnalysisBranch();
    }
        break;
    case 107:
    {
        mySimulationDataBase->createCombinedAnalysisRootNode();
        //! ------------------------------------------
        //! make the last inserted branch the current
        //! ------------------------------------------
        QStandardItem *itemRoot = Geometry_RootItem->parent();
        int Nb = itemRoot->rowCount();
        myTreeView->selectionModel()->setCurrentIndex(myModel->indexFromItem(itemRoot->child(Nb-1,0)),QItemSelectionModel::Current);
        this->setTheActiveAnalysisBranch();
    }
    case 106:
    {
        //! ------------------------------------------------------
        //! create a named selection from an item/items selection
        //! ------------------------------------------------------
        const QList<QModelIndex> &indexes = myTreeView->selectionModel()->selectedIndexes();
        QVector<GeometryTag> tags;
        for(int n=0; n<indexes.length(); n++)
        {
            SimulationNodeClass *node = indexes[n].data(Qt::UserRole).value<SimulationNodeClass*>();
            if(node->isAnalysisResult()==false && node->isSimulationSetUpNode()==false) continue;
            if(node->getPropertyItem("Tags")==Q_NULLPTR) continue;
            const QVector<GeometryTag> curTags = node->getPropertyValue<QVector<GeometryTag>>("Tags");
            tags.append(curTags);
        }
        if(tags.length()==0) return;
        this->createSimulationNode(SimulationNodeClass::nodeType_namedSelectionGeometry);
        SimulationNodeClass *namedSelectionNode = myTreeView->currentIndex().data(Qt::UserRole).value<SimulationNodeClass*>();
        QVariant data;
        data.setValue(tags);
        namedSelectionNode->replaceProperty("Geometry",Property("Geometry",data,Property::PropertyGroup_Scope));
        namedSelectionNode->replaceProperty("Tags",Property("Tags",data,Property::PropertyGroup_Scope));
    }
        break;
    case 54:
        //! call the mesher - volume mesh
        myIsMeshingRunning = true;
        this->buildMesh(true);
        break;
    case 55:
        //! call the mesher - surface mesh
        myIsMeshingRunning = true;
        this->buildMesh(false);
        break;
    case 56:
        //! clear generated data
        if(QMessageBox::question(this, tr(APPNAME),tr("Do you want to clear all the meshes?"), QMessageBox::Yes | QMessageBox::No)==QMessageBox::Yes)
        {
            emit requestClearMesh();
            updateMeshStatistics();
        }
        break;
    case 57: this->changeNodeSuppressionStatus(Property::SuppressionStatus_Suppressed); break;
    case 58: this->changeNodeSuppressionStatus(Property::SuppressionStatus_Active); break;
    case 68: this->createSimulationNode(SimulationNodeClass::nodeType_meshEdgeSize); break;
    case 69: this->createSimulationNode(SimulationNodeClass::nodeType_meshVertexSize); break;
    case 59: this->invertSuppressionSet(); break;
    case 61: this->unsuppressAllBodies(); break;
    case 62: this->createSimulationNode(SimulationNodeClass::nodeType_repairTool); break;
    case 63:
    {
        //! --------------------------------------
        //! hide bodies selected in the main tree
        //! --------------------------------------
        TColStd_ListOfInteger bodyNumbers;
        QList<QModelIndex> selectedIndexes = mySelectionModel->selectedIndexes();
        for(int i=0; i<selectedIndexes.length(); i++)
        {
            QModelIndex index = selectedIndexes.at(i);
            SimulationNodeClass *theCurNode = index.data(Qt::UserRole).value<SimulationNodeClass*>();
            int bodyIndex = theCurNode->getPropertyValue<int>("Map index");
            bodyNumbers.Append(bodyIndex);

            //! ----------------------------------
            //! function: change the icon opacity
            //! ----------------------------------
            QExtendedStandardItem *item = static_cast<QExtendedStandardItem*>(myModel->itemFromIndex(index));
            tools::changeIconOpacity(item,true);

            //! -----------------------------------
            //! synch the "Visible" property value
            //! -----------------------------------
            QVariant data;
            data.setValue(false);
            Property prop_visible("Visible",data,Property::PropertyGroup_GraphicProperties);
            theCurNode->replaceProperty("Visible",prop_visible);
        }
        emit requestHideBody(bodyNumbers);
    }
        break;
    case 64:
    {
        //! ---------------------------------------------------------------------
        //! hide all other bodies
        //! [1] find the selected body indexes and collect into a list
        //! [2] scan all the body indexes one by one. If the current body index
        //!     is not contained into the previous list, add to "allOtherBodies"
        //! ---------------------------------------------------------------------
        TColStd_ListOfInteger selectedBodiesIndexes;
        QList<QModelIndex> selectedIndexes = mySelectionModel->selectedRows(0);
        for(int i=0; i<selectedIndexes.length(); i++)
        {
            SimulationNodeClass* selectedNode = selectedIndexes.at(i).data(Qt::UserRole).value<SimulationNodeClass*>();
            int selectedBodyIndex = selectedNode->getPropertyItem("Map index")->data(Qt::UserRole).value<Property>().getData().toInt();
            selectedBodiesIndexes.Append(selectedBodyIndex);
        }

        QExtendedStandardItem *geomeytryRootItem = static_cast<QExtendedStandardItem*>(this->getTreeItem(SimulationNodeClass::nodeType_geometry));
        TColStd_ListOfInteger allOtherBodies;
        int NbBodies = geomeytryRootItem->rowCount();
        for(int i=0; i<NbBodies; i++)
        {
            QExtendedStandardItem *theCurItem = static_cast<QExtendedStandardItem*>(geomeytryRootItem->child(i,0));
            SimulationNodeClass *theCurNode = theCurItem->data(Qt::UserRole).value<SimulationNodeClass*>();
            int bodyIndex = theCurNode->getPropertyItem("Map index")->data(Qt::UserRole).value<Property>().getData().toInt();
            if(!selectedBodiesIndexes.Contains(bodyIndex))
            {
                allOtherBodies.Append(bodyIndex);

                //! ---------------------------
                //! synch the "Visible" switch
                //! ---------------------------
                QVariant data;
                data.setValue(false);
                Property prop_visible("Visible",data,Property::PropertyGroup_GraphicProperties);
                theCurNode->replaceProperty("Visible",prop_visible);

                //! -------------------------------------------
                //! change icon opacity => set isOpaque = true
                //! -------------------------------------------
                tools::changeIconOpacity(theCurItem,true);
            }
            else
            {
                //! --------------------------------------------
                //! change icon opacity => set isOpaque = false
                //! --------------------------------------------
                tools::changeIconOpacity(theCurItem,false);
            }

        }
        if(!allOtherBodies.IsEmpty())
        {
            requestHideBody(allOtherBodies);
        }
    }
        break;
    case 65:
    {
        //! show body
        TColStd_ListOfInteger bodiesToShow;
        QList<QModelIndex> selectedIndexes = mySelectionModel->selectedIndexes();
        for(int i=0; i<selectedIndexes.length(); i++)
        {
            SimulationNodeClass *curNode = selectedIndexes.at(i).data(Qt::UserRole).value<SimulationNodeClass*>();
            int bodyIndex = curNode->getPropertyItem("Map index")->data(Qt::UserRole).value<Property>().getData().toInt();

            cout<<"action 65:____body to show: "<<bodyIndex<<"____"<<endl;

            Property::SuppressionStatus isSuppressed = curNode->getPropertyItem("Suppressed")->data(Qt::UserRole).value<Property>().getData().value<Property::SuppressionStatus>();
            if(isSuppressed == Property::SuppressionStatus_Active)
            {
                QVariant data;
                data.setValue(true);
                Property prop_visible("Visible",data,Property::PropertyGroup_GraphicProperties);
                curNode->replaceProperty("Visible",prop_visible);
                bodiesToShow.Append(bodyIndex);

                //! ---------------------------------------------
                //! change icon opacity => set isOpaque = false
                //! ---------------------------------------------
                QExtendedStandardItem *curItem = static_cast<QExtendedStandardItem*>(myModel->itemFromIndex(selectedIndexes.at(i)));
                tools::changeIconOpacity(curItem,false);
            }
        }
        if(!bodiesToShow.IsEmpty()) emit requestShowBody(bodiesToShow);
    }
        break;
    case 66:
    {
        QExtendedStandardItem *itemGeometryRoot = static_cast<QExtendedStandardItem*>(this->getTreeItem(SimulationNodeClass::nodeType_geometry));
        QVariant data;
        data.setValue(true);
        Property prop_visible("Visible",data,Property::PropertyGroup_GraphicProperties);
        int NbBodies = itemGeometryRoot->rowCount();
        for(int i=0; i<NbBodies; i++)
        {
            QExtendedStandardItem *itemBody = static_cast<QExtendedStandardItem*>(itemGeometryRoot->child(i,0));
            SimulationNodeClass *curNode = itemBody->data(Qt::UserRole).value<SimulationNodeClass*>();
            curNode->replaceProperty("Visible",prop_visible);

            //! ------------------------------------------
            //! change icon opacity: set isOpaque == true;
            //! ------------------------------------------
            tools::changeIconOpacity(itemBody,false);
        }
        emit requestShowAllBodies();
    }
        break;

    case 70: this->createSimulationNode(SimulationNodeClass::nodeType_namedSelectionGeometry); break;
    case 87: this->createSimulationNode(SimulationNodeClass::nodeType_namedSelectionElement); break;
    case 88:
    {
        QList<QModelIndex> modelIndexList = myTreeView->selectionModel()->selectedIndexes();
        if(modelIndexList.size()<2) return;
        QVector<GeometryTag> tags;
        for(int i=0; i<modelIndexList.size(); i++)
        {
            const QModelIndex &curIndex = modelIndexList[i];
            SimulationNodeClass *curNode = curIndex.data(Qt::UserRole).value<SimulationNodeClass*>();
            tags.append(curNode->getPropertyValue<QVector<GeometryTag>>("Tags"));
        }
        //! --------------------------
        //! delete the selected items
        //! --------------------------
        //this->deleteItem(modelIndexList);
        this->createSimulationNode(SimulationNodeClass::nodeType_namedSelectionGeometry);

        SimulationNodeClass *node = myTreeView->currentIndex().data(Qt::UserRole).value<SimulationNodeClass*>();
        node->getModel()->blockSignals(true);
        QVariant data;
        data.setValue(tags);
        Property prop_scope("Geometry",data,Property::PropertyGroup_Scope);
        Property prop_tags("Tags",data,Property::PropertyGroup_Scope);
        node->replaceProperty("Geometry",prop_scope);
        node->replaceProperty("Tags",prop_tags);
        node->getModel()->blockSignals(false);
    }
        break;

    case 76: this->exportSTEPFile(); break;
    case 90: this->createSimulationNode(SimulationNodeClass::nodeType_thermalAnalysisTemperature); break;
    case 91: this->createSimulationNode(SimulationNodeClass::nodeType_thermalAnalysisConvection); break;
    case 92: this->createSimulationNode(SimulationNodeClass::nodeType_thermalAnalysisThermalFlux); break;
    case 93: this->createSimulationNode(SimulationNodeClass::nodeType_thermalAnalysisThermalFlow); break;
    case 94: this->createSimulationNode(SimulationNodeClass::nodeType_thermalAnalysisAdiabaticWall); break;
    case 96: this->createSimulationNode(SimulationNodeClass::nodeType_thermalAnalysisThermalPower); break;
    case 98: this->duplicateItem(); break;

    case 99:
    {
        //! ------------------------------------------------------------
        //! removing a mesh control causes the a mesh to become invalid
        //! ------------------------------------------------------------
        SimulationNodeClass *curNode = myTreeView->currentIndex().data(Qt::UserRole).value<SimulationNodeClass*>();
        SimulationNodeClass::nodeType theFamily = curNode->getFamily();
        if(theFamily==SimulationNodeClass::nodeType_meshControl)
        {
            QVector<GeometryTag> vecLocs = curNode->getPropertyValue<QVector<GeometryTag>>("Tags");
            std::vector<int> parentShapes;

            for(QVector<GeometryTag>::iterator it = vecLocs.begin(); it!=vecLocs.end(); ++it)
            {
                const GeometryTag &loc = *it;
                int n = loc.parentShapeNr;
                if(std::find(parentShapes.begin(),parentShapes.end(),n)==parentShapes.end()) parentShapes.push_back(n);
            }
            emit requestMeshInvalidate(parentShapes);
        }
        this->deleteItem();
    }
        break;

    case 200:
    {
        cout<<"SimulationManager::handleItem()->____update post object____"<<endl;
        SimulationNodeClass *curNode = myTreeView->currentIndex().data(Qt::UserRole).value<SimulationNodeClass*>();
        postObject curPostObject = curNode->getPropertyItem("Post object")->data(Qt::UserRole).value<Property>().getData().value<postObject>();
        curPostObject.update(static_cast<meshDataBase*>(mySimulationDataBase));
        //! ----
        //! to be implemented. Not sure if needed ...
    }
        break;
    case 201: this->createSimulationNode(SimulationNodeClass::nodeType_solutionStructuralNodalDisplacement,0); break;
    case 205: this->clearGeneratedData(); break;
    case 206: this->createSimulationNode(SimulationNodeClass::nodeType_solutionStructuralNodalDisplacement,1); break;

    //! --------------------------------
    //! 202 -> equivalent stress
    //! 202 -> stress intensity
    //! 207 -> maximum principal stress
    //! 208 -> middle principal stress
    //! 209 -> minimum principal stress
    //! --------------------------------
    case 202: this->createSimulationNode(SimulationNodeClass::nodeType_solutionStructuralStress,0); break;
    case 210: this->createSimulationNode(SimulationNodeClass::nodeType_solutionStructuralStress,1); break;
    case 207: this->createSimulationNode(SimulationNodeClass::nodeType_solutionStructuralStress,2); break;
    case 208: this->createSimulationNode(SimulationNodeClass::nodeType_solutionStructuralStress,3); break;
    case 209: this->createSimulationNode(SimulationNodeClass::nodeType_solutionStructuralStress,4); break;

    //! -------------------------------------------------
    //! 203 -> insert equivalent strain
    //! 214 -> insert strain intensity
    //! 211 -> insert max principal strain
    //! 212 -> insert middle principal strain
    //! 213 -> insert min principal strain
    //! 215 -> insert equivalent mechanical strain
    //! 216 -> insert mechanical strain intensity
    //! 219 -> insert max principal mechanical strain
    //! 220 -> insert middle principal mechanical strain
    //! -------------------------------------------------
    case 203: this->createSimulationNode(SimulationNodeClass::nodeType_solutionStructuralTotalStrain,0);break;
    case 214: this->createSimulationNode(SimulationNodeClass::nodeType_solutionStructuralTotalStrain,1); break;
    case 211: this->createSimulationNode(SimulationNodeClass::nodeType_solutionStructuralTotalStrain,2); break;
    case 212: this->createSimulationNode(SimulationNodeClass::nodeType_solutionStructuralTotalStrain,3); break;
    case 213: this->createSimulationNode(SimulationNodeClass::nodeType_solutionStructuralTotalStrain,4); break;
    case 215: this->createSimulationNode(SimulationNodeClass::nodeType_solutionStructuralTotalStrain,5); break;
    case 216: this->createSimulationNode(SimulationNodeClass::nodeType_solutionStructuralTotalStrain,6); break;
    case 219: this->createSimulationNode(SimulationNodeClass::nodeType_solutionStructuralTotalStrain,7); break;
    case 220: this->createSimulationNode(SimulationNodeClass::nodeType_solutionStructuralTotalStrain,8); break;

    //! 221 -> insert min principal mechanical strain
    case 221: this->createSimulationNode(SimulationNodeClass::nodeType_solutionStructuralMechanicalStrain,9); break;

    //! ----------------------------------------------
    //! 222 -> insert equivalent thermal strain
    //! 223 -> insert thermal strain intensity
    //! 224 -> insert max principal thermal strain
    //! 225 -> insert middle principal thermal strain
    //! 226 -> insert min principal thermal strain
    //! ----------------------------------------------
    case 222: this->createSimulationNode(SimulationNodeClass::nodeType_solutionStructuralThermalStrain,10); break;
    case 223: this->createSimulationNode(SimulationNodeClass::nodeType_solutionStructuralThermalStrain,11); break;
    case 224: this->createSimulationNode(SimulationNodeClass::nodeType_solutionStructuralThermalStrain,12); break;
    case 225: this->createSimulationNode(SimulationNodeClass::nodeType_solutionStructuralThermalStrain,13);break;
    case 226: this->createSimulationNode(SimulationNodeClass::nodeType_solutionStructuralThermalStrain,14); break;

    case 227: this->createSimulationNode(SimulationNodeClass::nodeType_solutionStructuralFatigueTool); break;
    case 228: this->createSimulationNode(SimulationNodeClass::nodeType_solutionStructuralTemperature); break;
    case 229: this->createSimulationNode(SimulationNodeClass::nodeType_solutionStructuralEquivalentPlasticStrain); break;

    //! ------------------------------------------------------------------
    //! 230 -> insert total "Nodal forces"
    //! directional nodal forces (option "1" => "x" direction by default)
    //! ------------------------------------------------------------------
    case 230: this->createSimulationNode(SimulationNodeClass::nodeType_solutionStructuralNodalForces,0); break;
    case 231: this->createSimulationNode(SimulationNodeClass::nodeType_solutionStructuralNodalForces,1); break;

    //! ----------------------------------------
    //! contact pressure => option "0"
    //! contact frictional stress => option "1"
    //! contact penetration => option "2"
    //! contact sliding => option "3"
    //! ----------------------------------------
    case 232: this->createSimulationNode(SimulationNodeClass::nodeType_solutionStructuralContact,0); break;
    case 233: this->createSimulationNode(SimulationNodeClass::nodeType_solutionStructuralContact,1); break;
    case 234: this->createSimulationNode(SimulationNodeClass::nodeType_solutionStructuralContact,2); break;
    case 235: this->createSimulationNode(SimulationNodeClass::nodeType_solutionStructuralContact,3); break;

    //! ------------------------------------------------------------------
    //! 246 -> insert total "Reaction forces"
    //! directional nodal forces (option "1" => "x" direction by default)
    //! ------------------------------------------------------------------
    case 246: this->createSimulationNode(SimulationNodeClass::nodeType_solutionStructuralReactionForce,0); break;
    case 247: this->createSimulationNode(SimulationNodeClass::nodeType_solutionStructuralReactionForce,1); break;

    //! ------------------------------------------------------------------
    //! 245 -> insert total "Gamma"
    //! ------------------------------------------------------------------
    case 245: this->createSimulationNode(SimulationNodeClass::nodeType_solutionStructuralGamma); break;

    //! -----------------
    //!  results
    //! -----------------
    case 204:
    {
        emit requestHideAllResults();
        this->callPostEngineEvaluateResult();
    }
        break;
    case 236: this->translateOpenFoamScalarData(); break;
    case 237: this->interpolate(); break;
    case 239: this->createSimulationNode(SimulationNodeClass::nodeType_modelChange); break;
    case 240: this->createSimulationNode(SimulationNodeClass::nodeType_solutionThermalTemperature); break;
    case 241: this->createSimulationNode(SimulationNodeClass::nodeType_solutionThermalFlux); break;
    case 40: this->createSimulationNode(SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Acceleration); break;
    case 41: this->createSimulationNode(SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RotationalVelocity); break;

    //! -------
    //! rename
    //! -------
    case 101:
    {
        //! -------------------
        //! start editing name
        //! -------------------
        QStandardItem *curItem = myModel->itemFromIndex(myTreeView->currentIndex());
        curItem->setEditable(true);
        myTreeView->edit(curItem->index());
        QItemDelegate *theDelegate = static_cast<QItemDelegate*>(myTreeView->itemDelegate());
        disconnect(theDelegate,SIGNAL(closeEditor(QWidget*,QAbstractItemDelegate::EndEditHint)),this,SLOT(renameItem()));
        connect(theDelegate,SIGNAL(closeEditor(QWidget*,QAbstractItemDelegate::EndEditHint)),this,SLOT(renameItem()));
    }
        break;

    case 103: this->renameItemBasedOnDefinition(); break;
    case 102: this->deleteAllChildrenItems(); break;
    case 43: this->createAutomaticConnections(); break;
    case 304:
    {
        ;
    }
        break;

    //! --------------------------------------
    //! generate and preview prismatic layers
    //! --------------------------------------
    case 67:
    {
        cout<<"____case 75: begin generating prismatic layers____"<<endl;
        this->previewPrismaticLayer();
    }
        break;
    case 89:
    {
        //! -----------------
        //! mesh metric item
        //! -----------------
        this->createSimulationNode(SimulationNodeClass::nodeType_meshMeshMetric);
    }
        break;
    case 83: mainTreeTools::formNewPart(myTreeView); break;
    case 84:
    {
        SimulationNodeClass *node = myTreeView->currentIndex().data(Qt::UserRole).value<SimulationNodeClass*>();
        int bodyIndex = node->getPropertyValue<int>("Map index");
        int res = QMessageBox::warning(this,"Delete body","Are you sure to delete body?\nThe operation is not recoverable",QMessageBox::Cancel,QMessageBox::Ok);
        if(res == QMessageBox::Ok)
        {
            mySimulationDataBase->removeBodyFromDataBase(bodyIndex);

            //! -----------------------------------
            //! remove the item from the tree view
            //! -----------------------------------
            myTreeView->model()->removeRow(myTreeView->currentIndex().row(),myTreeView->currentIndex().parent());

            //! -------------------
            //! refresh the viewer
            //! -------------------
            emit requestRefreshViewer();
        }
    }
        break;
    case 42: this->resetAndUpdateModel(); break;
    case 238:
    {
        //! -------------------------------------------------------
        //! retrieve the mesh data sources of a boundary condition
        //! -------------------------------------------------------
        SimulationNodeClass* node = myTreeView->currentIndex().data(Qt::UserRole).value<SimulationNodeClass*>();
        if(node->getPropertyItem("Mesh data sources")!=NULL)
        {
            IndexedMapOfMeshDataSources meshDataSources = node->getPropertyValue<IndexedMapOfMeshDataSources>("Mesh data sources");
            emit requestShowMeshDataSources(meshDataSources);
        }
    }
        break;

    case 300: mySimulationDataBase->createParticlesInFieldsRootNode(); break;
    case 301: this->createSimulationNode(SimulationNodeClass::nodeType_electrostaticPotential); break;

#ifdef COSTAMP_VERSION
    case 1000:
        cout<<"____case 1000: start the time step builder tool____"<<endl;
        this->COSTAMP_startTimeStepBuilder();
        break;
#endif
    }

    //for(myCTX->InitSelected();myCTX->MoreSelected();myCTX->NextSelected()) myCTX->ClearSelected(false);
    //myCTX->UpdateCurrentViewer();

    //! expand the tree branch if not expanded
    myTreeView->expand(myTreeView->selectionModel()->currentIndex());
}

//! -----------------------------------------------------------
//! function: showHealingElements
//! details:  slot - display elements used to repair the mesh
//! -----------------------------------------------------------
void SimulationManager::showHealingElements()
{
    cout<<"SimulationManager::showHealingElements()->____function called____"<<endl;
    this->requestHideMeshes();

    std::vector<int> vecBodyIndexes;
    ListOfShape listOfShapes;
    myCTX->InitSelected();
    if(myCTX->MoreSelected())
    {
        for(myCTX->InitSelected();myCTX->MoreSelected();myCTX->NextSelected())
        {
            const TopoDS_Shape &curShape = myCTX->SelectedShape();
            listOfShapes.Append(curShape);
        }
        QVector<GeometryTag> vecLoc = TopologyTools::generateLocationPairs(mySimulationDataBase,listOfShapes);
        for(QVector<GeometryTag>::iterator it = vecLoc.begin(); it!=vecLoc.end(); ++it)
        {
            GeometryTag loc = *it;
            vecBodyIndexes.push_back(loc.parentShapeNr);
        }
        //! eliminate duplicates
        vecBodyIndexes = tools::clearFromDuplicates(vecBodyIndexes);
    }
    else
    {
        for(int bodyIndex=1; bodyIndex<=mySimulationDataBase->bodyMap.size();bodyIndex++)
            vecBodyIndexes.push_back(bodyIndex);
    }

    for(std::vector<int>::iterator it = vecBodyIndexes.begin(); it!=vecBodyIndexes.end(); ++it)
    {
        int bodyIndex = *it;
        occHandle(Ng_MeshVS_DataSourceFace) meshDS = occHandle(Ng_MeshVS_DataSourceFace)::DownCast(mySimulationDataBase->ArrayOfMeshDSOnFaces.getValue(bodyIndex,0));
        if(!meshDS.IsNull())
        {
            cout<<"Body index: "<<bodyIndex<<" Number of elements: "<<meshDS->GetAllElements().Extent()<<endl;

            //! cosmetic for displaying the face mesh
            occHandle(MeshVS_Mesh) aFaceMesh = new MeshVS_Mesh();
            //! set the face mesh data source
            aFaceMesh->SetDataSource(meshDS);

            //! create the presentation builder
            occHandle(MeshVS_MeshPrsBuilder) aBuilder = new MeshVS_MeshPrsBuilder(aFaceMesh);

            //! add the presentation builder
            aFaceMesh->AddBuilder(aBuilder,Standard_False);

            //! the aspect
            Graphic3d_MaterialAspect myAspect(Graphic3d_NOM_GOLD);
            myAspect.SetColor(Quantity_NOC_RED);

            aFaceMesh->GetDrawer()->SetMaterial(MeshVS_DA_FrontMaterial, myAspect);
            aFaceMesh->GetDrawer()->SetBoolean(MeshVS_DA_DisplayNodes, Standard_False);
            aFaceMesh->GetDrawer()->SetBoolean(MeshVS_DA_ShowEdges, Standard_True);
            aFaceMesh->GetDrawer()->SetColor(MeshVS_DA_EdgeColor,Quantity_NOC_BLACK);
            aFaceMesh->SetDisplayMode(MeshVS_DMF_Shading);

            if(aFaceMesh.IsNull()) cout<<"occPreGLWidget::showFaceMesh()->____null face mesh____"<<endl;
            else
            {
                TColStd_ListOfInteger listOfBodies;
                for(int i=1; i<=mySimulationDataBase->bodyMap.size(); i++)
                    listOfBodies.Append(i);
                emit requestHideBody(listOfBodies);
                myCTX->Display(aFaceMesh,Standard_False);
            }
        }
        else
        {
            QMessageBox::information(this,APPNAME,"No surface elements have been added during surface mesh healing",QMessageBox::Ok);
        }
    }
}

//! -----------------------------------------------
//! function: getGeometrySelection
//! details:  check if something has been selected
//! -----------------------------------------------
ListOfShape SimulationManager::getGeometrySelection() const
{
    ListOfShape scope;
    scope.Clear();
    if(!myCTX.IsNull())
    {
        for(myCTX->InitSelected(); myCTX->MoreSelected(); myCTX->NextSelected())
            scope.Append(myCTX->SelectedShape());
        //! 24/12/2018 - commented
        //! myCTX->ClearSelected();
    }
    return scope;
}

//! -----------------------
//! function: addParentTag
//! details:  helper
//! -----------------------
void SimulationManager::addParentTimeTag(SimulationNodeClass *aNode)
{
    //! ---------------------------------------------
    //! aim - retrieve the root of the current item
    //! and generate the parent time tag property
    //! ---------------------------------------------
    if(myTreeView->currentIndex().parent().isValid()==false) return;

    QVariant data;
    SimulationNodeClass *n = myTreeView->currentIndex().parent().data(Qt::UserRole).value<SimulationNodeClass*>();
    //if(n->getType()!=SimulationNodeClass::nodeType_structuralAnalysis && n->getType()!=SimulationNodeClass::nodeType_thermalAnalysis)
    //{
    //    n = myTreeView->currentIndex().data(Qt::UserRole).value<SimulationNodeClass*>();
    //}

    //! -----------------------------
    //! check if the property exists
    //! -----------------------------
    if(n->getPropertyItem("Time tag")==Q_NULLPTR)
    {
        cerr<<"____Time tag property not defined for node: "<<aNode->getName().toStdString()<<"____"<<endl;
        return;
    }

    const QString &timeTag = n->getPropertyValue<QString>("Time tag");
    data.setValue(timeTag);
    aNode->addProperty(Property("Parent time tag",data,Property::PropertyGroup_Identifier));
}
//! -------------------------------
//! function: createSimulationNode
//! details:
//! -------------------------------
void SimulationManager::createSimulationNode(SimulationNodeClass::nodeType type, QVariant addOptions)
{    
    //! ---------
    //! the node
    //! ---------
    SimulationNodeClass *aNode;

    //! -------------
    //! generic data
    //! -------------
    QVariant data;

    //! -----------------------------
    //! check the geometry selection
    //! -----------------------------
    const ListOfShape &scope = this->getGeometrySelection();

    //! -------------
    //! the new item
    //! -------------
    QExtendedStandardItem *item = new QExtendedStandardItem();

    //! -------------
    //! MODEL CHANGE
    //! -------------
    if(type == SimulationNodeClass::nodeType_modelChange)
    {
        aNode = nodeFactory::nodeFromScratch(type);
        aNode->setParent(this);

        item->setData(aNode->getName(),Qt::DisplayRole);
        data.setValue(aNode);
        item->setData(data,Qt::UserRole);
        mainTreeTools::getCurrentSimulationRoot(myTreeView)->insertRow(this->getInsertionRow(),item);

        //! ------------------------------------
        //! retrieve the item Analysis Settings
        //! ------------------------------------
        QStandardItem *curItem = myModel->itemFromIndex(mySelectionModel->currentIndex());
        SimulationNodeClass *curNode = curItem->data(Qt::UserRole).value<SimulationNodeClass*>();
        QStandardItem *itemAnalysisSettings;
        if(curNode->isAnalysisRoot()) itemAnalysisSettings = curItem->child(0,0);
        else itemAnalysisSettings = curItem->parent()->child(0,0);
        SimulationNodeClass *nodeAnalysisSettings = itemAnalysisSettings->data(Qt::UserRole).value<SimulationNodeClass*>();

        //! ---------------------------------------------------------------
        //! add "Time info" - these properties must be added here, since
        //! the nodeFactory does not have access to the simulation manager
        //! ---------------------------------------------------------------
        load aLoad;
        QVariant aLoadData; aLoadData.setValue(Property::modelChangeActivationStatus_Inactive);
        QVector<QVariant> loadData; loadData.push_back(aLoadData);
        aLoad.setType(Property::loadType_modelChange);
        aLoad.setData(loadData);
        nodeAnalysisSettings->getTabularDataModel()->appendColumn(aLoad);

        //! ---------------------------------
        //! add the switch to the node model
        //! ---------------------------------
        data.setValue(Property::modelChangeActivationStatus_Inactive);
        Property prop_activationStatus("Activation status",data,Property::PropertyGroup_Definition);
        aNode->addProperty(prop_activationStatus);
    }
    //! ------------------------------------------
    //! HANDLING STRUCTURAL POST PROCESSING ITEMS
    //! ------------------------------------------
    else if(type ==SimulationNodeClass::nodeType_solutionStructuralNodalDisplacement ||
            type ==SimulationNodeClass::nodeType_solutionStructuralTotalStrain ||
            type ==SimulationNodeClass::nodeType_solutionStructuralMechanicalStrain ||
            type ==SimulationNodeClass::nodeType_solutionStructuralThermalStrain ||
            type ==SimulationNodeClass::nodeType_solutionStructuralEquivalentPlasticStrain ||
            type ==SimulationNodeClass::nodeType_solutionStructuralStress ||
            type ==SimulationNodeClass::nodeType_solutionStructuralFatigueTool ||
            type ==SimulationNodeClass::nodeType_solutionStructuralNodalForces ||
            type ==SimulationNodeClass::nodeType_solutionStructuralContact  ||
            type ==SimulationNodeClass::nodeType_solutionStructuralGamma  ||
            type ==SimulationNodeClass::nodeType_solutionStructuralReactionForce)
    {
        aNode = nodeFactory::nodeFromScratch(type,mySimulationDataBase, myCTX, addOptions);
        aNode->setParent(this);

        item->setData(aNode->getName(),Qt::DisplayRole);
        data.setValue(aNode);
        item->setData(data,Qt::UserRole);

        //! --------------------------------------------------------------
        //! search for the "Solution" item starting from the current item
        //! --------------------------------------------------------------
        QStandardItem *curItem = myModel->itemFromIndex(myTreeView->currentIndex());
        SimulationNodeClass *curNode = curItem->data(Qt::UserRole).value<SimulationNodeClass*>();
        QStandardItem *solutionItem;
        if(curNode->isSolution()) solutionItem = curItem;
        else if(curNode->isSolutionInformation()) solutionItem = curItem->parent();
        if(curNode->isAnalysisResult()) solutionItem = curItem->parent();
        //! ----------------
        //! append the item
        //! ----------------
        solutionItem->appendRow(item);
    }
    else if (type ==SimulationNodeClass::nodeType_solutionStructuralTemperature)
    {
        aNode = nodeFactory::nodeFromScratch(type,mySimulationDataBase);
        aNode->setParent(this);

        item->setData(aNode->getName(),Qt::DisplayRole);
        data.setValue(aNode);
        item->setData(data,Qt::UserRole);
        //! --------------------------------------------------------------
        //! search for the "Solution" item starting from the current item
        //! --------------------------------------------------------------
        QStandardItem *curItem = myModel->itemFromIndex(myTreeView->currentIndex());
        SimulationNodeClass *curNode = curItem->data(Qt::UserRole).value<SimulationNodeClass*>();
        QStandardItem *solutionItem;
        if(curNode->isSolution()) solutionItem = curItem;
        else if(curNode->isSolutionInformation()) solutionItem = curItem->parent();
        if(curNode->isAnalysisResult()) solutionItem = curItem->parent();
        //! ----------------
        //! append the item
        //! ----------------
        solutionItem->appendRow(item);
    }
    //! ---------------------------------------
    //! HANDLING THERMAL POST PROCESSING ITEMS
    //! ---------------------------------------
    else if(type ==SimulationNodeClass::nodeType_solutionThermalTemperature ||
            type ==SimulationNodeClass::nodeType_solutionThermalFlux)
    {
        aNode = nodeFactory::nodeFromScratch(type,mySimulationDataBase);
        aNode->setParent(this);

        item->setData(aNode->getName(),Qt::DisplayRole);
        data.setValue(aNode);
        item->setData(data,Qt::UserRole);
        //! --------------------------------------------------------------
        //! search for the "Solution" item starting from the current item
        //! --------------------------------------------------------------
        QStandardItem *curItem = myModel->itemFromIndex(myTreeView->currentIndex());
        SimulationNodeClass *curNode = curItem->data(Qt::UserRole).value<SimulationNodeClass*>();
        QStandardItem *solutionItem;
        if(curNode->isSolution()) solutionItem = curItem;
        else if(curNode->isSolutionInformation()) solutionItem = curItem->parent();
        if(curNode->isAnalysisResult()) solutionItem = curItem->parent();
        //! ----------------
        //! append the item
        //! ----------------
        solutionItem->appendRow(item);
    }
    //! ------------------------------------------
    //! HANDLING MAPPER ROOT/IMPORTED BODY SCALAR
    //! ------------------------------------------
    else if(type==SimulationNodeClass::nodeType_importedBodyScalar ||
            type==SimulationNodeClass::nodeType_mapper ||
            type==SimulationNodeClass::nodeType_OpenFoamScalarData)
    {
        aNode = nodeFactory::nodeFromScratch(type);
        aNode->setParent(this);
        data.setValue(aNode);
        item->setData(data,Qt::UserRole);
        data.setValue(aNode->getName());
        item->setData(data,Qt::DisplayRole);

        //! -------------------------
        //! get tree insertion point
        //! -------------------------
        switch(type)
        {
        case SimulationNodeClass::nodeType_importedBodyScalar:
        case SimulationNodeClass::nodeType_OpenFoamScalarData:
            this->myModel->itemFromIndex(myTreeView->currentIndex())->appendRow(item);
            break;
        case SimulationNodeClass::nodeType_mapper:
            mainTreeTools::getCurrentSimulationRoot(myTreeView)->insertRow(this->getInsertionRow(),item);
            break;
        }
    }
    //! ---------------------------------
    //! HANDLING COORDINATE SYSTEM ITEMS
    //! ---------------------------------
    else if(type==SimulationNodeClass::nodeType_coordinateSystem)
    {
        //! --------------------------------------------------------------
        //! count the number of coordinate systems for handling the names
        //! --------------------------------------------------------------
        int NS = CoordinateSystems_RootItem->rowCount();
        QString CS_Name = QString("Coordinate system %1").arg(NS);

        aNode = nodeFactory::nodeFromScratch(SimulationNodeClass::nodeType_coordinateSystem,mySimulationDataBase,myCTX);
        aNode->setParent(this);
        aNode->setName(CS_Name);
        data.setValue(aNode);
        item->setData(data, Qt::UserRole);
        data.setValue(aNode->getName());
        item->setData(data, Qt::DisplayRole);
        CoordinateSystems_RootItem->insertRow(NS,item);
    }
    //! -----------------------
    //! HANDLING REMOTE POINTS
    //! -----------------------
    else if(type==SimulationNodeClass::nodeType_remotePoint)
    {
        RemotePoint_RootItem = this->getTreeItem(SimulationNodeClass::nodeType_remotePointRoot);
        if(RemotePoint_RootItem==Q_NULLPTR)
        {
            cerr<<"SimulationManager::createSimulationNode()->___the remote point root is NULL____"<<endl;
            return;
        }

        //! --------------------------------------------------------------
        //! passing a NULL meshDataBase and a NULL AIS_InteractiveContext
        //! since not needed
        //! --------------------------------------------------------------
        aNode = nodeFactory::nodeFromScratch(type,0,0,addOptions);

        //! ----------------------------------------------------------------------
        //! this part is not currently implemented within the "nodeFactory" class
        //! since it does not have access to the SimulationManager
        //! -----------------------------------------------------------------------
        QExtendedStandardItem *itemCSroot = this->getTreeItem(SimulationNodeClass::nodeType_coordinateSystems);
        QExtendedStandardItem *itemGlobalCS = static_cast<QExtendedStandardItem*>(itemCSroot->child(0,0));
        void *itemGlobalCSvoid = (void*)itemGlobalCS;
        data.setValue(itemGlobalCSvoid);
        Property prop_coordinateSystem("Coordinate system",data,Property::PropertyGroup_Scope);
        aNode->addProperty(prop_coordinateSystem,3);
        data.setValue(aNode);
        aNode->setParent(this);

        item->setData("Remote point", Qt::DisplayRole);
        item->setData(data, Qt::UserRole);

        RemotePoint_RootItem->appendRow(item);

        bool isDone = markerBuilder::addMarker(this->getCurrentNode(), mySimulationDataBase);
        if(isDone == true) this->displayMarker();
    }
    //! ------------------------------------------
    //! HANDLING THERMAL BOUNDARY CONDITION ITEMS
    //! ------------------------------------------
    else if(type==SimulationNodeClass::nodeType_thermalAnalysisAdiabaticWall)
    {
        aNode = nodeFactory::nodeFromScratch(type,mySimulationDataBase,myCTX);
        aNode->setParent(this);
        data.setValue(aNode->getName());
        item->setData(data,Qt::DisplayRole);
        data.setValue(aNode);
        item->setData(data,Qt::UserRole);
        mainTreeTools::getCurrentSimulationRoot(myTreeView)->insertRow(this->getInsertionRow(),item);

        emit request2DBodySelectionMode(true);
        emit requestStartEditingScope();
    }
    else if(type==SimulationNodeClass::nodeType_thermalAnalysisRadiation ||
            type==SimulationNodeClass::nodeType_thermalAnalysisTemperature ||
            type==SimulationNodeClass::nodeType_thermalAnalysisThermalFlow ||
            type==SimulationNodeClass::nodeType_thermalAnalysisThermalFlux ||
            type==SimulationNodeClass::nodeType_thermalAnalysisThermalPower)
    {
        aNode = nodeFactory::nodeFromScratch(type,mySimulationDataBase,myCTX);
        aNode->setParent(this);
        data.setValue(aNode->getName());
        item->setData(data,Qt::DisplayRole);
        data.setValue(aNode);
        item->setData(data,Qt::UserRole);
        mainTreeTools::getCurrentSimulationRoot(myTreeView)->insertRow(this->getInsertionRow(),item);

        //! ------------------------------------------------------------------------------
        //! access the "Analysis settings" item: if the current node is a simulation root
        //! the Analysis Settings item is the first child
        //! ------------------------------------------------------------------------------
        SimulationNodeClass *nodeAnalysisSettings = Q_NULLPTR;
        SimulationNodeClass *curNode = myTreeView->currentIndex().data(Qt::UserRole).value<SimulationNodeClass*>();
        if(curNode->isAnalysisRoot())
        {
            nodeAnalysisSettings = myTreeView->currentIndex().child(0,0).data(Qt::UserRole).value<SimulationNodeClass*>();
        }
        else nodeAnalysisSettings = myTreeView->currentIndex().parent().child(0,0).data(Qt::UserRole).value<SimulationNodeClass*>();
        if(nodeAnalysisSettings==Q_NULLPTR) return;

        //! -----------------------------
        //! append a column to the table
        //! -----------------------------
        load aLoad;
        switch(type)
        {
        case SimulationNodeClass::nodeType_thermalAnalysisTemperature: aLoad.setType(Property::loadType_temperatureMagnitude); break;
        case SimulationNodeClass::nodeType_thermalAnalysisThermalFlux: aLoad.setType(Property::loadType_thermalFluxMagnitude); break;
        case SimulationNodeClass::nodeType_thermalAnalysisThermalFlow: aLoad.setType(Property::loadType_thermalFlowMagnitude); break;
        case SimulationNodeClass::nodeType_thermalAnalysisThermalPower: aLoad.setType(Property::loadType_thermalPowerMagnitude); break;
        }
        nodeAnalysisSettings->getTabularDataModel()->appendColumn(aLoad);
    }
    else if(type==SimulationNodeClass::nodeType_thermalAnalysisConvection)
    {
        aNode = nodeFactory::nodeFromScratch(type,mySimulationDataBase,myCTX);
        aNode->setParent(this);
        data.setValue(aNode->getName());
        item->setData(data,Qt::DisplayRole);
        data.setValue(aNode);
        item->setData(data,Qt::UserRole);
        mainTreeTools::getCurrentSimulationRoot(myTreeView)->insertRow(this->getInsertionRow(),item);

        //! ------------------------------------------------------------------------------
        //! access the "Analysis settings" item: if the current node is a simulation root
        //! the Analysis Settings item is the first child
        //! ------------------------------------------------------------------------------
        SimulationNodeClass *nodeAnalysisSettings = Q_NULLPTR;
        SimulationNodeClass *curNode = myTreeView->currentIndex().data(Qt::UserRole).value<SimulationNodeClass*>();
        if(curNode->isAnalysisRoot())
        {
            nodeAnalysisSettings = myTreeView->currentIndex().child(0,0).data(Qt::UserRole).value<SimulationNodeClass*>();
        }
        else nodeAnalysisSettings = myTreeView->currentIndex().parent().child(0,0).data(Qt::UserRole).value<SimulationNodeClass*>();
        if(nodeAnalysisSettings==Q_NULLPTR) return;

        //! ----------------------
        //! film coefficient load
        //! ----------------------
        QVector<QVariant> vecData;
        data.setValue(0.0);
        vecData.push_back(data);
        load aLoad1(vecData,Property::loadType_thermalConvectionFilmCoefficientMagnitude);
        nodeAnalysisSettings->getTabularDataModel()->appendColumn(aLoad1);

        //! ---------------------------
        //! reference temperature load
        //! ---------------------------
        vecData.clear();
        data.setValue(0.0);
        vecData.push_back(data);
        load aLoad2(vecData,Property::loadType_thermalConvectionReferenceTemperatureMagnitude);
        nodeAnalysisSettings->getTabularDataModel()->appendColumn(aLoad2);
    }
    //! ----------------------------------------------
    //! HANDLING STRUCTURAL BOUNDARY CONDITIONS ITEMS
    //! ----------------------------------------------
    else if(type==SimulationNodeClass::nodeType_structuralAnalysisBoltPretension)
    {
        QStandardItem *itemGlobalCS = this->getTreeItem(SimulationNodeClass::nodeType_coordinateSystem_global);
        QVariant addOptions;
        void *p = (void*)(itemGlobalCS);
        addOptions.setValue(p);
        aNode = nodeFactory::nodeFromScratch(type,mySimulationDataBase,myCTX,addOptions);
        aNode->setParent(this);
        data.setValue(aNode);
        item->setData(data,Qt::UserRole);
        data.setValue(aNode->getName());
        item->setData(data,Qt::DisplayRole);

        //! -------------------------------------------------------------------------------
        //! tabular data for "Bolt pretension"
        //! loadType_boltStatusDefinedBy: 0 => Force; 1 => Adjustment; 2 => Open 3 => Lock
        //! loadType_boltForce: the preload force
        //! loadType_boltAdjustment: the adjustment [L]
        //! -------------------------------------------------------------------------------
        cout<<"SimulationManager::createSimulationNode()->____creating tabular data for bolt____"<<endl;

        QVector<QVariant> values;
        data.setValue(Property::boltStatusDefinedBy_load);
        values.push_back(data);
        load loadBoltStatusDefinedBy(values,Property::loadType_boltStatusDefinedBy);

        QVector<QVariant> values1;
        data.setValue(0.0);     //! value of the pretension [double]
        values1.push_back(data);
        load loadForce(values1,Property::loadType_boltForce);
        QVector<QVariant> values2;
        data.setValue(0.0);     //! value of the displacement [double]
        values2.push_back(data);
        load loadAdjustment(values2,Property::loadType_boltAdjustment);

        //! -------------------------------------------
        //! append the three columns to the load table
        //! -------------------------------------------
        SimulationNodeClass *nodeAnalysisSettings = this->getAnalysisSettingsNodeFromCurrentItem();
        nodeAnalysisSettings->getTabularDataModel()->appendColumn(loadBoltStatusDefinedBy);
        nodeAnalysisSettings->getTabularDataModel()->appendColumn(loadForce);
        nodeAnalysisSettings->getTabularDataModel()->appendColumn(loadAdjustment);

        markerBuilder::addMarker(this->getCurrentNode(),mySimulationDataBase);
        mainTreeTools::getCurrentSimulationRoot(myTreeView)->insertRow(this->getInsertionRow(),item);
    }
    else if(type==SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_CylindricalSupport ||
            type==SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_FrictionlessSupport ||
            type==SimulationNodeClass::nodeType_structuralAnalysisBoundaryContidion_FixedSupport ||
            type==SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_CompressionOnlySupport)
    {
        aNode = nodeFactory::nodeFromScratch(type,mySimulationDataBase,myCTX);
        aNode->setParent(this);
        data.setValue(aNode->getName());
        item->setData(data,Qt::DisplayRole);
        data.setValue(aNode);
        item->setData(data,Qt::UserRole);
        mainTreeTools::getCurrentSimulationRoot(myTreeView)->insertRow(this->getInsertionRow(),item);
    }
    else if(type==SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Displacement ||
            type==SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteDisplacement ||
            type==SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteRotation)
    {
        aNode = nodeFactory::nodeFromScratch(type,mySimulationDataBase,myCTX);
        aNode->setParent(this);

        //! -------------------------------------------------------------------
        //! the "Coordinate system" property is not handled by the nodeFactory
        //! since the nodeFactory does not have access to the main tree: the
        //! property is created and added here
        //! -------------------------------------------------------------------
        QExtendedStandardItem *itemCSroot = this->getTreeItem(SimulationNodeClass::nodeType_coordinateSystems);
        QExtendedStandardItem *itemGlobalCS = static_cast<QExtendedStandardItem*>(itemCSroot->child(0,0));
        void *itemGlobalCSvoid = (void*)itemGlobalCS;
        data.setValue(itemGlobalCSvoid);
        Property prop_CS("Coordinate system",data,Property::PropertyGroup_Definition);
        aNode->addProperty(prop_CS);

        data.setValue(aNode);
        item->setData(data,Qt::UserRole);
        data.setValue(aNode->getName());
        item->setData(data,Qt::DisplayRole);
        mainTreeTools::getCurrentSimulationRoot(myTreeView)->insertRow(this->getInsertionRow(),item);

        //! ---------------------------------------------------------------
        //! begin editing X component: don't parse the item immediately...
        //! ---------------------------------------------------------------
        if(scope.Size()==0) emit requestStartEditingScope();
        else emit requestStartEditingXcomponent();
    }
    else if(type==SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Moment ||
            type==SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Force ||
            type==SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteForce ||
            type==SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Acceleration ||
            type==SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RotationalVelocity ||
            type==SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Pressure ||
            type==SimulationNodeClass::nodeType_structuralAnalysisThermalCondition)
    {
        aNode = nodeFactory::nodeFromScratch(type,mySimulationDataBase,myCTX);
        aNode->setParent(this);
        data.setValue(aNode);
        item->setData(data,Qt::UserRole);
        data.setValue(aNode->getName());
        item->setData(data,Qt::DisplayRole);

        //! -------------------------------------------------------------------------------------------
        //! Force, moment, pressure, temperature, bolts forces/displacements values are stored in the
        //! form of tabular data contained in the Analysis settings node. There is one single table
        //! for all the loads, instead of a table for each load.
        //! In the following, the tabular data for the node are created, and added to the tabular data
        //! model in "Analysis settings"
        //! ------------------------------------------------------------------------------------------

        //! ------------------------------------
        //! access the "Analysis settings" item
        //! ------------------------------------
        SimulationNodeClass *nodeAnalysisSettings = this->getAnalysisSettingsNodeFromCurrentItem();

        //! --------------------------------------------------------------------------
        //! the default constructor load::load() inits a load of type "loadType_none"
        //! --------------------------------------------------------------------------
        load aLoad;
        data.setValue(0.0);
        QVector<QVariant> vecData{data};
        aLoad.setData(vecData);
        switch(type)
        {
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Moment:
        {
            aLoad.setType(Property::loadType_momentMagnitude);
            nodeAnalysisSettings->getTabularDataModel()->appendColumn(aLoad);
        }
            break;

        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Force:
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteForce:
        {
            aLoad.setType(Property::loadType_forceMagnitude);
            nodeAnalysisSettings->getTabularDataModel()->appendColumn(aLoad);
        }
            break;

        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Acceleration:
        {
            aLoad.setType(Property::loadType_accelerationMagnitude);
            nodeAnalysisSettings->getTabularDataModel()->appendColumn(aLoad);
        }
            break;

        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RotationalVelocity:
        {
            aLoad.setType(Property::loadType_rotationalVelocityMagnitude);
            nodeAnalysisSettings->getTabularDataModel()->appendColumn(aLoad);
        }
            break;

        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Pressure:
        {
            aLoad.setType(Property::loadType_pressureMagnitude);
            nodeAnalysisSettings->getTabularDataModel()->appendColumn(aLoad);
        }
            break;

        case SimulationNodeClass::nodeType_structuralAnalysisThermalCondition:
        {
            //! --------------------------------------------------------------------
            //! read the default setting for the property "Environment temperature"
            //! --------------------------------------------------------------------
            double Tenv_default = this->getTreeItem(SimulationNodeClass::nodeType_structuralAnalysis)->data(Qt::UserRole).value<SimulationNodeClass*>()
                    ->getPropertyValue<double>("Environment temperature");

            //! --------------------------------------------------------
            //! modify the tabular data in Analysis settings
            //! column of data containing the environment temperature
            //! --------------------------------------------------------
            vecData.clear();
            data.setValue(Tenv_default);
            vecData.push_back(data);
            aLoad.setType(Property::loadType_thermalConditionTemperature);
            aLoad.setData(vecData);
            nodeAnalysisSettings->getTabularDataModel()->appendColumn(aLoad);
        }
            break;
        }
        markerBuilder::addMarker(aNode,mySimulationDataBase);
        mainTreeTools::getCurrentSimulationRoot(myTreeView)->insertRow(this->getInsertionRow(),item);
    }
    //! --------------------------
    //! imported body temperature
    //! --------------------------
    else if(type==SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_ImportedTemperatureDistribution)
    {
        aNode = nodeFactory::nodeFromScratch(type,mySimulationDataBase);
        aNode->setParent(this);

        item->setData(aNode->getName(),Qt::DisplayRole);
        data.setValue(aNode);
        item->setData(data,Qt::UserRole);

        mainTreeTools::getCurrentSimulationRoot(myTreeView)->insertRow(this->getInsertionRow(),item);
    }
    //! -----------
    //! POINT MASS
    //! -----------
    else if(type==SimulationNodeClass::nodeType_pointMass)
    {
        emit request2DBodySelectionMode(true);
        aNode = nodeFactory::nodeFromScratch(type,mySimulationDataBase,myCTX);
        aNode->setParent(this);
        item->setData(aNode->getName(),Qt::DisplayRole);
        QVariant data;
        data.setValue(aNode);
        item->setData(data,Qt::UserRole);
        Geometry_RootItem->appendRow(item);
    }
    //! -----------------------
    //! HANDLING MESH CONTROLS
    //! -----------------------
    else if(type==SimulationNodeClass::nodeType_meshMethod ||
            type==SimulationNodeClass::nodeType_meshBodyMeshMethod ||
            type==SimulationNodeClass::nodeType_meshBodyMeshControl ||
            type==SimulationNodeClass::nodeType_meshMeshType ||
            type==SimulationNodeClass::nodeType_meshMeshMetric)
    {
        emit request3DBodySelectionMode(true);
        data.setValue(scope);
        aNode = nodeFactory::nodeFromScratch(type,mySimulationDataBase,myCTX,data);
        aNode->setParent(this);        
        QVariant data;
        data.setValue(aNode->getName());
        item->setData(data,Qt::DisplayRole);
        data.setValue(aNode);
        item->setData(data,Qt::UserRole);
        Mesh_RootItem->appendRow(item);
    }
    else if(type==SimulationNodeClass::nodeType_meshFaceSize)
    {
        emit request2DBodySelectionMode(true);
        aNode = nodeFactory::nodeFromScratch(type,mySimulationDataBase,myCTX);
        aNode->setParent(this);
        item->setData(aNode->getName(),Qt::DisplayRole);
        QVariant data;
        data.setValue(aNode);
        item->setData(data,Qt::UserRole);
        Mesh_RootItem->appendRow(item);
    }
    else if(type==SimulationNodeClass::nodeType_meshEdgeSize)
    {
        emit request1DBodySelectionMode(true);
        aNode = nodeFactory::nodeFromScratch(type,mySimulationDataBase,myCTX);
        aNode->setParent(this);
        item->setData(aNode->getName(),Qt::DisplayRole);
        QVariant data;
        data.setValue(aNode);
        item->setData(data,Qt::UserRole);
        Mesh_RootItem->appendRow(item);
    }
    else if(type==SimulationNodeClass::nodeType_meshVertexSize)
    {
        emit request0DBodySelectionMode(true);
        aNode = nodeFactory::nodeFromScratch(type,mySimulationDataBase,myCTX);
        aNode->setParent(this);
        item->setData(aNode->getName(),Qt::DisplayRole);
        QVariant data;
        data.setValue(aNode);
        item->setData(data,Qt::UserRole);
        Mesh_RootItem->appendRow(item);
    }
    else if(type==SimulationNodeClass::nodeType_meshPrismaticLayer)
    {
        emit request2DBodySelectionMode(true);
        aNode = nodeFactory::nodeFromScratch(type,mySimulationDataBase,myCTX);
        aNode->setParent(this);
        item->setData(aNode->getName(),Qt::DisplayRole);
        QVariant data;
        data.setValue(aNode);
        item->setData(data,Qt::UserRole);
        Mesh_RootItem->appendRow(item);
    }
    //! -----------------------
    //! HANDLING CONTACT PAIRS
    //! -----------------------
    else if(type==SimulationNodeClass::nodeType_connectionGroup)
    {
        emit request3DBodySelectionMode(true);
        aNode =nodeFactory::nodeFromScratch(SimulationNodeClass::nodeType_connectionGroup,mySimulationDataBase,myCTX);
        aNode->setParent(this);
        QVariant data;
        data.setValue(aNode);
        item->setData(data,Qt::UserRole);
        item->setData(aNode->getName(),Qt::DisplayRole);
        Connections_RootItem->appendRow(item);
    }
    else if(type==SimulationNodeClass::nodeType_connectionPair)
    {
        //! ---------------------------------------------------------------------------
        //! A "Contact pair" is always inserted as child of a "Connection group" item.
        //! If the "Connection group" does not exists, it is created
        //! ---------------------------------------------------------------------------
        //! the current clicked item
        QModelIndex theCurrentIndex = myTreeView->currentIndex();
        QExtendedStandardItem *theCurrentClickedItem = static_cast<QExtendedStandardItem*>(this->myModel->itemFromIndex(theCurrentIndex));
        SimulationNodeClass *theNode = theCurrentIndex.data(Qt::UserRole).value<SimulationNodeClass*>();
        SimulationNodeClass::nodeType theNodeType = theNode->getType();

        QString parentTimeTag;
        QExtendedStandardItem *itemConnectionGroup;
        switch(theNodeType)
        {
        case SimulationNodeClass::nodeType_connectionPair:
        {
            //! [a] retrieve the connection group
            itemConnectionGroup = static_cast<QExtendedStandardItem*>(theCurrentClickedItem->parent());

            //! [b] get the time tag of the "Connection group"
            parentTimeTag = itemConnectionGroup->data(Qt::UserRole).value<SimulationNodeClass*>()->getPropertyValue<QString>("Time tag");
         }
            break;

        case SimulationNodeClass::nodeType_connectionGroup:
        {
            //! [a] retrieve the connection group
            itemConnectionGroup = static_cast<QExtendedStandardItem*>(theCurrentClickedItem);

            //! [b] get the time tag of the "Connection group"
            parentTimeTag = itemConnectionGroup->data(Qt::UserRole).value<SimulationNodeClass*>()->getPropertyValue<QString>("Time tag");
        }
            break;

        case SimulationNodeClass::nodeType_connection:
        {
            //! the "Connection" root has been clicked
            //! [a] create a connection group, which is made "current" immediately
            this->createSimulationNode(SimulationNodeClass::nodeType_connectionGroup);
            itemConnectionGroup = static_cast<QExtendedStandardItem*>(this->myModel->itemFromIndex(myTreeView->currentIndex()));

            //! [b] get the time tag of the "Connection group"
            parentTimeTag = itemConnectionGroup->data(Qt::UserRole).value<SimulationNodeClass*>()->getPropertyValue<QString>("Time tag");
        }
            break;
        }

        //! ----------------------------------------------------------------------------------------
        //! the contact pair is created with the time tag of the Connection group it is attached to
        //! => add the "Time tag" to the contact pair generation options
        //! ----------------------------------------------------------------------------------------
        connectionPairGenerationOption curConnectionGenerationOptions = addOptions.value<connectionPairGenerationOption>();
        curConnectionGenerationOptions.timeTag = parentTimeTag;
        cout<<"____"<<parentTimeTag.toStdString()<<"____"<<endl;
        addOptions.setValue(curConnectionGenerationOptions);

        aNode = nodeFactory::nodeFromScratch(SimulationNodeClass::nodeType_connectionPair,mySimulationDataBase,myCTX,addOptions);
        aNode->setParent(this);
        QVariant data;
        data.setValue(aNode);
        item->setData(data,Qt::UserRole);
        item->setData(aNode->getName(),Qt::DisplayRole);
        itemConnectionGroup->appendRow(item);
        itemConnectionGroup->setEditable(false);
    }
    //! --------------------------
    //! HANDLING NAMED SELECTIONS
    //! --------------------------
    else if(type==SimulationNodeClass::nodeType_namedSelectionGeometry)
    {
        aNode =nodeFactory::nodeFromScratch(SimulationNodeClass::nodeType_namedSelectionGeometry,mySimulationDataBase,myCTX);
        aNode->setParent(this);
        data.setValue(aNode);
        item->setData(data,Qt::UserRole);
        item->setData(aNode->getName(),Qt::DisplayRole);
        NamedSelection_RootItem->appendRow(item);
    }
    else if(type==SimulationNodeClass::nodeType_namedSelectionElement)
    {
        aNode =nodeFactory::nodeFromScratch(SimulationNodeClass::nodeType_namedSelectionElement,mySimulationDataBase,myCTX);
        aNode->setParent(this);
        data.setValue(aNode);
        item->setData(data,Qt::UserRole);
        item->setData(aNode->getName(),Qt::DisplayRole);
        NamedSelection_RootItem->appendRow(item);
    }

    //! -----------------------------
    //! HANDLING PARTICLES IN FIELDS
    //! -----------------------------
    else if(type==SimulationNodeClass::nodeType_electrostaticPotential)
    {
        emit request2DBodySelectionMode(true);
        data.setValue(scope);
        aNode = nodeFactory::nodeFromScratch(type,mySimulationDataBase,myCTX,data);
        aNode->setParent(this);
        QVariant data;
        data.setValue(aNode->getName());
        item->setData(data,Qt::DisplayRole);
        data.setValue(aNode);
        item->setData(data,Qt::UserRole);
        markerBuilder::addMarker(this->getCurrentNode(),mySimulationDataBase);
        mainTreeTools::getCurrentSimulationRoot(myTreeView)->insertRow(this->getInsertionRow(),item);
    }
    //! ------------------------------
    //! make the new item the current
    //! ------------------------------
    myTreeView->setCurrentIndex(item->index());
    this->highlighter();

    //! ---------------------------------
    //! set the new item as not editable
    //! ---------------------------------
    item->setEditable(false);

    //! ------------------------------
    //! parse the new item definition
    //! ------------------------------
    //parser::parseItem(item);

    //! ----------------
    //! parent time tag
    //! ----------------
    if(aNode->getType()!=SimulationNodeClass::nodeType_connectionPair)
        this->addParentTimeTag(aNode);

    connect(aNode->getModel(),SIGNAL(itemChanged(QStandardItem*)),this,SLOT(handleItemChange(QStandardItem*)));

    //! set the current analysis branch
    this->setTheActiveAnalysisBranch();
}

//! -----------------------
//! function: getTreeItem
//! details:
//! -----------------------
QExtendedStandardItem* SimulationManager::getTreeItem(SimulationNodeClass::nodeType theNodeType)
{
    //cout<<"SimulationManager::getTreeItem()->____function called____"<<endl;
    if(myModel==Q_NULLPTR)
    {
        cerr<<"SimulationManager::getTreeItem()->____NULL model____"<<endl;
        return Q_NULLPTR;
    }

    QList<QExtendedStandardItem*> items;
    this->getTreeItemsRecursively(myModel,items);
    for(QList<QExtendedStandardItem*>::iterator it = items.begin(); it!=items.end(); it++)
    {
        QExtendedStandardItem* curItem = *it;
        SimulationNodeClass *curNode = curItem->data(Qt::UserRole).value<SimulationNodeClass*>();

        if(curNode==Q_NULLPTR) return Q_NULLPTR;

        SimulationNodeClass::nodeType curNodeType = curNode->getType();
        if(curNodeType==theNodeType)
        {
            //cout<<"SimulationManager::getTreeItem()->____found item: "<<curItem->data(Qt::UserRole).value<SimulationNodeClass*>()->getName().toStdString()<<"____"<<endl;
            return curItem;
        }
    }
    //cerr<<"SimulationManager::getTreeItem()->____item not found____"<<endl;
    return Q_NULLPTR;
}

//! -------------------------------
//! function: getAllTreeItemOfType
//! details:
//! -------------------------------
QList<QExtendedStandardItem*> SimulationManager::getAllTreeItemOfType(SimulationNodeClass::nodeType theNodeType)
{
    if(myModel==Q_NULLPTR) return QList<QExtendedStandardItem*>();
    QList<QExtendedStandardItem*> items, itemsout;
    this->getTreeItemsRecursively(myModel,items);
    for(QList<QExtendedStandardItem*>::iterator it = items.begin(); it!=items.end(); ++it)
    {
        QExtendedStandardItem* curItem = *it;
        if(curItem->data(Qt::UserRole).value<SimulationNodeClass*>()->getType()==theNodeType) itemsout.append(curItem);
    }
    return itemsout;
}

//! ----------------------
//! function: setDataBase
//! details:  old...
//! ----------------------
void SimulationManager::setDataBase(simulationDataBase *aDB)
{
    mySimulationDataBase = aDB;
    myModel = aDB->getModel();

    Geometry_RootItem = this->getTreeItem(SimulationNodeClass::nodeType_geometry);
    CoordinateSystems_RootItem = this->getTreeItem(SimulationNodeClass::nodeType_coordinateSystems);
    RemotePoint_RootItem = this->getTreeItem(SimulationNodeClass::nodeType_remotePointRoot);
    Mesh_RootItem = this->getTreeItem(SimulationNodeClass::nodeType_meshControl);
    Connections_RootItem = this->getTreeItem(SimulationNodeClass::nodeType_connection);
    NamedSelection_RootItem = this->getTreeItem(SimulationNodeClass::nodeType_namedSelection);

    myTreeView->setModel(myModel);

    //! ------------------------------------------
    //! hide the empty "Named selection"
    //! hide the empty "Remote point"
    //! hide the dummy contact "Select from list"
    //! ------------------------------------------
    if(NamedSelection_RootItem!=Q_NULLPTR) myTreeView->setRowHidden(0, NamedSelection_RootItem->index(), true);
    if(RemotePoint_RootItem!=Q_NULLPTR) myTreeView->setRowHidden(0, RemotePoint_RootItem->index(), true);
    if(Connections_RootItem!=Q_NULLPTR) myTreeView->setRowHidden(0, Connections_RootItem->index(), true);

    myTreeView->expandAll();
}

//! -----------------------
//! function: loadCADModel
//! details:
//! -----------------------
bool SimulationManager::loadCADModel(const QString &fileName,
                                     TopoDS_Compound &shapeFromReader,
                                     QList<QString> &listOfNames,
                                     const occHandle(QOccProgressIndicator) &aProgressIndicator)
{
    cout<<"SimulationManager::loadCADModel()->____function called____"<<endl;

    //! ------------------------------------------
    //! the model loader: init with the file name
    //! and load the geometry
    //! ------------------------------------------
    ModelLoader aModelLoader(fileName);
    bool isDone = aModelLoader.perform(shapeFromReader, listOfNames, aProgressIndicator);

    //! ------------------------------------------------------------------
    //! check the list of names: if a name is empty assign a default name
    //! ------------------------------------------------------------------
    for(int i=0, index = 0; i<listOfNames.length(); i++)
    {
        if(listOfNames.at(i).isEmpty())
        {
            index++;
            listOfNames.replace(i,QString("body %1").arg(index));
        }
    }

    QStandardItem *importItem = this->getTreeItem(SimulationNodeClass::nodeType_import);
    QVariant data;
    data.setValue(fileName);
    Property prop_sourceFilePath("Source file path",data,Property::PropertyGroup_Definition);
    cout<<"____replacing property \"Source file path\"____"<<endl;

    importItem->data(Qt::UserRole).value<SimulationNodeClass*>()->replaceProperty("Source file path",prop_sourceFilePath);

    return isDone;
}

//! ----------------------------------------
//! function: createSimulationDataBaseEmpty
//! details:
//! ----------------------------------------
void SimulationManager::createSimulationDataBaseEmpty()
{
    cout<<"SimulationManager::createSimulationDataBaseEmpty()->____function called____"<<endl;

    //! -------------------------------
    //! create the simulation database
    //! -------------------------------
    TopoDS_Shape aShape;
    mySimulationDataBase = new simulationDataBase(aShape, "", this);

    //! ------------------------------------------------------------------
    //! access the simulation database model: the class member "myModel"
    //! must be initialized here since "SimulationManager::getTreeItem()"
    //! which is used after internally needs it
    //! ------------------------------------------------------------------
    myModel = mySimulationDataBase->getModel();

    //! -----------------------------
    //! attach the view to the model
    //! -----------------------------
    myTreeView->setModel(myModel);

    //! ------------------------
    //! set the selection model
    //! ------------------------
    mySelectionModel = myTreeView->selectionModel();

    //! ------------------------
    //! trigger the highlighter
    //! ------------------------
    disconnect(myTreeView,SIGNAL(clicked(QModelIndex)),this,SLOT(highlighter(QModelIndex)));
    connect(myTreeView,SIGNAL(clicked(QModelIndex)),this,SLOT(highlighter(QModelIndex)));

    //! -----------------------------------
    //! set up the main tree item pointers
    //! -----------------------------------
    QStandardItem *Importer_Item = this->getTreeItem(SimulationNodeClass::nodeType_import);
    SimulationNodeClass *importerNode = Importer_Item->data(Qt::UserRole).value<SimulationNodeClass*>();
    connect(importerNode->getModel(),SIGNAL(itemChanged(QStandardItem*)),this,SLOT(handleItemChange(QStandardItem*)));

    Geometry_RootItem = this->getTreeItem(SimulationNodeClass::nodeType_geometry);
    CoordinateSystems_RootItem = this->getTreeItem(SimulationNodeClass::nodeType_coordinateSystems);
    RemotePoint_RootItem = this->getTreeItem(SimulationNodeClass::nodeType_remotePointRoot);
    Mesh_RootItem = this->getTreeItem(SimulationNodeClass::nodeType_meshControl);
    Connections_RootItem = this->getTreeItem(SimulationNodeClass::nodeType_connection);
    NamedSelection_RootItem = this->getTreeItem(SimulationNodeClass::nodeType_namedSelection);

    //! -------------------------------------------------
    //! hide the dummy named selectio "Select from list"
    //! hide the dummy contact pair "Select from list"
    //! hide the dummy remote point "Select from list"
    //! -------------------------------------------------
    if(NamedSelection_RootItem!=Q_NULLPTR) myTreeView->setRowHidden(0, NamedSelection_RootItem->index(), true);
    if(Connections_RootItem!=Q_NULLPTR) myTreeView->setRowHidden(0, Connections_RootItem->index(), true);
    if(RemotePoint_RootItem!=Q_NULLPTR) myTreeView->setRowHidden(0, RemotePoint_RootItem->index(), true);

    //! ----------------
    //! expand the tree
    //! ----------------
    myTreeView->expandAll();

    //! -----------------------------------------------------
    //! set the simulation data base pointer for the clipper
    //! -----------------------------------------------------
    emit requestSetClipperDataBase(mySimulationDataBase);
}

/*
//! -----------------------------------
//! function: createSimulationDataBase
//! details:  old version
//! -----------------------------------
void SimulationManager::createSimulationDataBase(const TopoDS_Shape &shapeFromReader,
                                                 const QString &fileName,
                                                 const QList<QString> &listOfNames)
{
    cout<<"SimulationManager::createSimulationDataBase()->____function called____"<<endl;

    //! -------------------------------
    //! create the simulation database
    //! -------------------------------
    mySimulationDataBase = new simulationDataBase(shapeFromReader, fileName, this);

    //! ------------------------------------------------------------------
    //! access the simulation database model: the class member "myModel"
    //! must be initialized here since "SimulationManager::getTreeItem()"
    //! which is used after internally needs it
    //! ------------------------------------------------------------------
    myModel = mySimulationDataBase->getModel();

    if(!shapeFromReader.IsNull())
    {
        //! ---------------------------------------------------
        //! fill the array of the names read from the CAD file
        //! ---------------------------------------------------
        int index = 1;
        for(int i=0; i<listOfNames.length(); i++, index++)
        {
            const QString &name = listOfNames.at(i);
            TCollection_AsciiString acname(name.toStdString().c_str());
            mySimulationDataBase->MapOfBodyNames.insert(index,name);
            ccout(QString("SimulationManager::createSimulationDataBase()->____found body with name: ").append(name).append("____"));
        }

        //! -----------------------------------------------
        //! update the names of the items and of the nodes
        //! -----------------------------------------------
        mySimulationDataBase->transferNames();
    }

    //! ----------------------------------
    //! access the roots of the main tree
    //! ----------------------------------
    Geometry_RootItem = this->getTreeItem(SimulationNodeClass::nodeType_geometry);
    CoordinateSystems_RootItem = this->getTreeItem(SimulationNodeClass::nodeType_coordinateSystems);
    RemotePoint_RootItem = this->getTreeItem(SimulationNodeClass::nodeType_remotePointRoot);
    Mesh_RootItem = this->getTreeItem(SimulationNodeClass::nodeType_meshControl);
    Connections_RootItem = this->getTreeItem(SimulationNodeClass::nodeType_connection);
    NamedSelection_RootItem = this->getTreeItem(SimulationNodeClass::nodeType_namedSelection);

    //! ----------------------------------
    //! init the tree view with the model
    //! ----------------------------------
    myTreeView->setModel(myModel);

    //! ---------------------------------------------------------------------------------------
    //! The "Analysis settings" item is internally created by the "SimulationDataBase" class.
    //! It has not the SIGNAL/SLOT connection linking the changes triggered by the user to the
    //! SimulationManager::handleItemChange(). This connection is created here: for other,
    //! dynamically inserted items, it is created at the instant of item insertion.
    //!
    //! The node "Analysis setting" can be of type structural or thermal (mutually exclusive)
    //!
    //! Warning: Remember to do that also at reload time
    //! ---------------------------------------------------------------------------------------
    //SimulationNodeClass *nodeAnalysisSettings;
    //QExtendedStandardItem *analysisSettingsItem = this->getTreeItem(SimulationNodeClass::nodeType_structuralAnalysisSettings);
    //if(analysisSettingsItem ==NULL) analysisSettingsItem = this->getTreeItem(SimulationNodeClass::nodeType_thermalAnalysisSettings);
    //nodeAnalysisSettings = analysisSettingsItem->data(Qt::UserRole).value<SimulationNodeClass*>();
    //connect(nodeAnalysisSettings->getModel(),SIGNAL(itemChanged(QStandardItem*)),this,SLOT(handleItemChange(QStandardItem*)));

    //! -------------------------------------------------
    //! hide the dummy named selectio "Select from list"
    //! hide the dummy contact pair "Select from list"
    //! hide the dummy remote point "Select from list"
    //! -------------------------------------------------
    if(NamedSelection_RootItem!=Q_NULLPTR) myTreeView->setRowHidden(0, NamedSelection_RootItem->index(), true);
    if(Connections_RootItem!=Q_NULLPTR) myTreeView->setRowHidden(0, Connections_RootItem->index(), true);
    if(RemotePoint_RootItem!=Q_NULLPTR) myTreeView->setRowHidden(0, RemotePoint_RootItem->index(), true);

    this->buildMeshIO();

    //! -----------
    //! expand all
    //! -----------
    myTreeView->expandAll();
}
*/

//! -----------------------------------
//! function: createSimulationDataBase
//! details:  new version
//! -----------------------------------
void SimulationManager::createSimulationDataBase(const TopoDS_Shape &shapeFromReader,
                                                 const QString &fileName,
                                                 const QList<QString> &listOfNames)
{
    cout<<"SimulationManager::createSimulationDataBase()->____function called____"<<endl;

    //! -------------------------------
    //! create the simulation database
    //! -------------------------------
    mySimulationDataBase->update(shapeFromReader);

    if(!shapeFromReader.IsNull())
    {
        //! ---------------------------------------------------
        //! fill the array of the names read from the CAD file
        //! ---------------------------------------------------
        int index = 1;
        for(int i=0; i<listOfNames.length(); i++, index++)
        {
            const QString &name = listOfNames.at(i);
            mySimulationDataBase->MapOfBodyNames.insert(index,name);
            cout<<"SimulationManager::createSimulationDataBase()->____found body with name: "<<name.toStdString()<<"____"<<endl;
        }

        //! -----------------------------------------------
        //! update the names of the items and of the nodes
        //! -----------------------------------------------
        mySimulationDataBase->transferNames();
    }

    //! -----------
    //! expand all
    //! -----------
    myTreeView->expandAll();

#ifdef COSTAMP_VERSION
    //this->COSTAMP_addProcessParameters();
#endif
    //ccout("SimulationManager::createSimulationDataBase()->____database created____");
}

//! --------------------------------------------------
//! function: transferMeshNodes
//! details:  this function scan the mesh tree items
//!           and copy the properties (mesh controls)
//!           into the array of parameters
//! --------------------------------------------------
void SimulationManager::transferMeshNodes()
{
    cout<<"SimulationManager::transferMeshNodes()->____function called____"<<endl;

    //! -----------------------------------------------------------------------
    //! first reset the meshing parameters using the default parameters values
    //! -----------------------------------------------------------------------
    mySimulationDataBase->calculateDefaultMeshingParameters();
    mySimulationDataBase->applyDefaultMeshingParameters();

    //! ----------------------------------------------
    //! the "global" mesh controls the mesh root item
    //! ----------------------------------------------

    //! ----------
    //! relevance
    //! ----------
    QExtendedStandardItem *itemMeshRoot = this->getTreeItem(SimulationNodeClass::nodeType_meshControl);

    int r = itemMeshRoot->data(Qt::UserRole).value<SimulationNodeClass*>()->getPropertyItem("Relevance")->data(Qt::UserRole).value<Property>().getData().toInt();
    //! -------------------
    //! function relevance
    //! -------------------
    double R = 2*(1-1/(1+exp(-0.025*double(r))));
    cout<<"SimulationManager::transferMeshNodes()->____relevance: "<<r<<"____"<<endl;

    //! ----------
    //! smoothing
    //! ----------
    int smoothing = itemMeshRoot->data(Qt::UserRole).value<SimulationNodeClass*>()->getPropertyItem("Smoothing")->data(Qt::UserRole).value<Property>().getData().toInt();
    int nsteps;
    switch(smoothing)
    {
    case 0: nsteps = 0; break;
    case 1: nsteps = 2; break;
    case 2: nsteps = 5; break;
    case 3: nsteps = 8; break;
    }

    //! ------------------
    //! initial size seed
    //! ------------------
    int val = itemMeshRoot->data(Qt::UserRole).value<SimulationNodeClass*>()->getPropertyValue<int>("Initial size seed");
    switch(val)
    {
    case 0: case 1:
    {
        //! ----------------
        //! active assembly
        //! ----------------
        double diag,minlen;
        const TopoDS_Shape &wholeShape = mySimulationDataBase->shape();
        mySimulationDataBase->shapeSize(wholeShape,diag,minlen);
        for(int bodyIndex = 1; bodyIndex<=mySimulationDataBase->bodyMap.size(); bodyIndex++)
        {
            mySimulationDataBase->ArrayOfMaxBodyElementSize.insert(bodyIndex,R*diag/(10.0));
        }
    }
        break;

    case 2:
        //! --------
        //! by part
        //! --------
        for(int bodyIndex = 1; bodyIndex<=mySimulationDataBase->bodyMap.size(); bodyIndex++)
        {
            double diag,minlen;
            const TopoDS_Shape &curShape = mySimulationDataBase->bodyMap.value(bodyIndex);
            mySimulationDataBase->shapeSize(curShape,diag,minlen);
            mySimulationDataBase->ArrayOfMaxBodyElementSize.insert(bodyIndex,R*diag/(10.0));
        }
        break;
    }
    val = itemMeshRoot->data(Qt::UserRole).value<SimulationNodeClass*>()->getPropertyValue<int>("Element midside nodes");

    //! ---------------------------
    //! midside nodes (mesh order)
    //! ---------------------------
    Property::meshOrder meshOrder;
    switch (val)
    {
    case 0: meshOrder = Property::meshOrder_First; break;   //! program controlled
    case 1: meshOrder = Property::meshOrder_First; break;   //! dropped
    case 2: meshOrder = Property::meshOrder_Second; break;  //! kept
    }

    //! ------------------------
    //! straight sided elements
    //! ------------------------
    bool isStraight = itemMeshRoot->data(Qt::UserRole).value<SimulationNodeClass*>()->getPropertyValue<bool>("Straight sided elements");

    //! ----------------
    //! smoothing steps
    //! ----------------
    for(QMap<int,TopoDS_Shape>::iterator it=mySimulationDataBase->bodyMap.begin(); it!=mySimulationDataBase->bodyMap.end(); ++it)
    {
        int bodyIndex = it.key();
        mySimulationDataBase->ArrayOfSmoothingSteps.insert(bodyIndex,nsteps);
        mySimulationDataBase->ArrayOfMeshOrder.insert(bodyIndex,meshOrder);
        mySimulationDataBase->mapOfIsElementStraight.insert(bodyIndex,isStraight);
    }

    //! ---------------------------------------------------------
    //! local mesh controls: found within each mesh control item
    //! ---------------------------------------------------------
    if(Mesh_RootItem->hasChildren())
    {
        for(int i=0;i<Mesh_RootItem->rowCount();i++)
        {
            //! -------------------------------
            //! get the i-th mesh control node
            //! -------------------------------
            SimulationNodeClass *meshControlNode = Mesh_RootItem->child(i,0)->data(Qt::UserRole).value<SimulationNodeClass*>();

            switch(meshControlNode->getType())
            {
            case SimulationNodeClass::nodeType_meshBodyMeshControl:
            {
                int Nb = 0;
                double Grading, minElementSize, maxElementSize;
                bool isSuppressed = true;
                bool isScopeOK = false;
                bool isGradingOK = false;
                bool isMinElementSizeOK = false;
                bool isMaxElementSizeOK = false;

                //! retrieve the vector of properties
                QVector<QExtendedStandardItem*> propertyItems = meshControlNode->getPropertyItems();

                QVector<Property> props;

                //! items containing <Property>(ies)
                for(QVector<QExtendedStandardItem*>::iterator it=propertyItems.begin(); it!=propertyItems.end();++it)
                {
                    QExtendedStandardItem *curItem = *it;
                    if(curItem->data(Qt::UserRole).canConvert<Property>())
                        props.append(curItem->data(Qt::UserRole).value<Property>());
                }

                std::vector<int> bodyIndex_;    //! a vector of body indexes
                for(QVector<Property>::iterator it = props.begin(); it!=props.end(); ++it)
                {                    
                    const Property &curProp = *it;
                    if(curProp.getGroup()==Property::PropertyGroup_Scope)
                    {
                        QVector<GeometryTag> vecLocs;
                        if(curProp.getName()=="Tags")
                        {
                            vecLocs = curProp.getData().value<QVector<GeometryTag>>();

                            cout<<"SimulationManager::transferMeshNodes()->____property \"Tags\" found; extent: "<<vecLocs.size()<<"____"<<endl;

                            int j=0;
                            for(QVector<GeometryTag>::iterator vecIt = vecLocs.begin();  vecIt!=vecLocs.end(); ++vecIt)
                            {
                                GeometryTag curLoc = *vecIt;
                                bodyIndex_.push_back(curLoc.parentShapeNr);
                                Nb = j+1;
                                j++;
                            }
                            if(Nb>0)isScopeOK = true;
                            else isScopeOK = false;
                        }
                    }
                    else if(curProp.getGroup() == Property::PropertyGroup_Definition)
                    {
                        if(curProp.getName()=="Suppressed")
                        {
                            Property::SuppressionStatus ss = curProp.getData().value<Property::SuppressionStatus>();
                            if(ss==Property::SuppressionStatus_Active)isSuppressed = false;
                            else isSuppressed = true;
                            //!cout<<"SimulationManager::transferMeshNodes()->____found property suppression status: "<<isSuppressed<<"____"<<endl;
                        }
                        else if(curProp.getName()=="Min element size")
                        {
                            minElementSize = curProp.getData().toDouble();
                            //!cout<<"SimulationManager::transferMeshNodes()->____found property minimum element size: "<<minElementSize<<"____"<<endl;
                            isMinElementSizeOK = true;
                        }
                        else if(curProp.getName()=="Max element size")
                        {
                            maxElementSize = curProp.getData().toDouble();
                            //!cout<<"SimulationManager::transferMeshNodes()->____found property maximum element size: "<<maxElementSize<<"____"<<endl;
                            isMaxElementSizeOK = true;
                        }
                        else if(curProp.getName()=="Grading")
                        {
                            Grading = curProp.getData().toDouble();
                            //!cout<<"SimulationManager::transferMeshNodes()->____found property grading: "<<Grading<<"____"<<endl;
                            isGradingOK = true;
                        }
                    }

                    //! --------------------------------
                    //! copy the values into the arrays
                    //! --------------------------------
                    if(isScopeOK && isGradingOK && isMinElementSizeOK && isMaxElementSizeOK && !isSuppressed)
                    {
                        //!cout<<"SimulationManager::transferMeshNodes()->____copying the values____"<<endl;
                        for(int j=0;j<Nb;j++)
                        {
                            mySimulationDataBase->ArrayOfGradingValue.insert(bodyIndex_.at(j),Grading);
                            mySimulationDataBase->ArrayOfMaxBodyElementSize.insert(bodyIndex_.at(j),maxElementSize);
                            mySimulationDataBase->ArrayOfMinBodyElementSize.insert(bodyIndex_.at(j),minElementSize);
                        }
                    }
                }
            }
                break;

            case SimulationNodeClass::nodeType_meshMeshType:
            {
                int meshType = 0;

                //! -----------------------------------------------------------
                //! check if the control is Active and if the Scope is defined
                //! -----------------------------------------------------------
                bool isActive = false;
                Property::SuppressionStatus ss = meshControlNode->getPropertyValue<Property::SuppressionStatus>("Suppressed");
                if(ss==Property::SuppressionStatus_Active) isActive = true;

                bool isScopeOK = false;
                QVector<GeometryTag> vecLoc = meshControlNode->getPropertyValue<QVector<GeometryTag>>("Tags");
                if(!vecLoc.isEmpty()) isScopeOK = true;

                //! ---------------------------------
                //! "Active" is on and "Scope" is ok
                //! ---------------------------------
                if(isActive == true && isScopeOK == true)
                {
                    meshType = meshControlNode->getPropertyValue<int>("Mesh type");
                }

                //! -----------------------------
                //! actually transfer parameters
                //! -----------------------------
                for(QVector<GeometryTag>::iterator it = vecLoc.begin(); it!=vecLoc.end(); ++it)
                {
                    GeometryTag aLoc = *it;
                    int bodyIndex = aLoc.parentShapeNr;
                    mySimulationDataBase->ArrayOfMeshType.insert(bodyIndex,meshType);
                }
            }
                break;

            case SimulationNodeClass::nodeType_meshFaceSize:
            {
                //!cout<<"SimulationManager::transferMeshNodes()->____face sizing control found____"<<endl;
                double faceMeshSize;
                bool isSuppressed = true;
                bool isFaceMeshSizeOK = false;
                bool isScopeOK = false;

                //! retrieve the vector of properties
                QVector<QExtendedStandardItem*> propertyItems = meshControlNode->getPropertyItems();
                QVector<Property> props;

                for(QVector<QExtendedStandardItem*>::Iterator it = propertyItems.begin();it!=propertyItems.end();++it)
                {
                    QExtendedStandardItem* item = *it;
                    props.push_back(item->data(Qt::UserRole).value<Property>());
                }

                std::vector<pair<int,int>> vecPairs;   //! a vector of pairs (bodyNumber, faceNumber)
                for(int k=0;k<propertyItems.length();k++)
                {
                    if(props.at(k).getGroup()==Property::PropertyGroup_Scope)
                    {
                        QVector<GeometryTag> vecLocs;
                        if(props.at(k).getName()=="Tags")
                        {
                            vecLocs = props.at(k).getData().value<QVector<GeometryTag>>();
                            if(vecLocs.size()>0)
                            {
                                //!cout<<"SimulationManager::transferMeshNodes()->____property scope found; pairs: "<<vecLocs.size()<<"____"<<endl;
                                std::pair<int,int> aPair;
                                for(int i=0; i<vecLocs.size();i++)
                                {
                                    GeometryTag loc = vecLocs.at(i);
                                    //cout<<"---->("<<loc.parentShapeNr<<", "<<loc.subTopNr<<" - "<<loc.isParent<<")"<<endl;
                                    aPair.first = loc.parentShapeNr;
                                    aPair.second = loc.subTopNr;
                                    vecPairs.push_back(aPair);
                                }
                                isScopeOK=true;
                            }
                        }
                    }
                    else if(props.at(k).getGroup()==Property::PropertyGroup_Definition)
                    {
                        if(props.at(k).getName()=="Suppressed")
                        {
                            //! the control is "suppressed", i.e. is not accounted for
                            Property::SuppressionStatus ss = props.at(k).getData().value<Property::SuppressionStatus>();
                            if(ss==Property::SuppressionStatus_Active)isSuppressed = false;
                            else isSuppressed = true;
                            //!cout<<"SimulationManager::transferMeshNodes()->____found property suppression status: "<<isSuppressed<<"____"<<endl;
                        }
                        else if(props.at(k).getName()=="Face sizing")
                        {
                            faceMeshSize = props.at(k).getData().toDouble();
                            isFaceMeshSizeOK=true;
                        }
                    }
                }
                //! --------------------------------
                //! copy the values into the arrays
                //! --------------------------------
                if(!isSuppressed && isFaceMeshSizeOK && isScopeOK)
                {
                    //!cout<<"SimulationManager::transferMeshNodes()->____copying the values____"<<endl;
                    for(int n=0;n<vecPairs.size();n++)
                    {
                        mySimulationDataBase->MapOfElementSizeOnFace.setValue(vecPairs[n].first, vecPairs[n].second, faceMeshSize);
                        mySimulationDataBase->MapOfIsFaceModified.setValue(vecPairs[n].first, vecPairs[n].second, true);
                    }
                }
            }
                break;

            case SimulationNodeClass::nodeType_meshEdgeSize:
            {
                //!cout<<"SimulationManager::transferMeshNodes()->____edge sizing control found____"<<endl;
                double edgeMeshSize;
                int numberOfDivisions;
                int typeOfSizing;

                bool isSuppressed = true;
                bool isEdgeMeshSizeOK = false;
                bool isNumberOfDivisionsOK = false;
                bool isScopeOK = false;

                //! retrieve the vector of properties
                QVector<QExtendedStandardItem*> propertyItems = meshControlNode->getPropertyItems();
                QVector<Property> props;

                for(QVector<QExtendedStandardItem*>::Iterator it = propertyItems.begin();it!=propertyItems.end();++it)
                {
                    QExtendedStandardItem* item = *it;
                    props.push_back(item->data(Qt::UserRole).value<Property>());
                }

                std::vector<pair<int,int>> vecPairs;   //! a vector of pairs (bodyNumber, edgeNumber)
                for(int k=0;k<propertyItems.length();k++)
                {
                    if(props.at(k).getGroup()==Property::PropertyGroup_Scope)
                    {
                        QVector<GeometryTag> vecLocs;
                        if(props.at(k).getName()=="Tags")
                        {
                            vecLocs = props.at(k).getData().value<QVector<GeometryTag>>();
                            if(vecLocs.size()>0)
                            {
                                //!cout<<"SimulationManager::transferMeshNodes()->____property scope found; pairs: "<<vecLocs.size()<<"____"<<endl;
                                std::pair<int,int> aPair;
                                for(int i=0; i<vecLocs.size();i++)
                                {
                                    GeometryTag loc = vecLocs.at(i);
                                    //cout<<"---->("<<loc.parentShapeNr<<", "<<loc.subTopNr<<" - "<<loc.isParent<<")"<<endl;
                                    aPair.first = loc.parentShapeNr;
                                    aPair.second = loc.subTopNr;
                                    vecPairs.push_back(aPair);
                                }
                                isScopeOK=true;
                            }
                        }
                    }
                    else if(props.at(k).getGroup()==Property::PropertyGroup_Definition)
                    {
                        if(props.at(k).getName()=="Suppressed")
                        {
                            Property::SuppressionStatus ss = props.at(k).getData().value<Property::SuppressionStatus>();
                            if(ss==Property::SuppressionStatus_Active)isSuppressed = false;
                            else isSuppressed = true;
                            //!cout<<"SimulationManager::transferMeshNodes()->____found property suppression status: "<<isSuppressed<<"____"<<endl;
                        }
                        else if(props.at(k).getName()=="Sizing type")
                        {
                            typeOfSizing = props.at(k).getData().toInt();
                            //! this is always OK;
                        }
                        //! -------------------------------------------------------------------------------------
                        //! both the properties "Element size" and "Number of divisions" are defined and present
                        //! in the DetailViewer: one of them is hidden, according to the value of "Sizing type"
                        //! -------------------------------------------------------------------------------------
                        else if(props.at(k).getName()=="Element size")
                        {
                            edgeMeshSize = props.at(k).getData().toDouble();
                            isEdgeMeshSizeOK=true;
                        }
                        else if(props.at(k).getName()=="Number of divisions")
                        {
                            numberOfDivisions = props.at(k).getData().toInt();
                            isNumberOfDivisionsOK=true;
                        }
                    }
                }

                //! --------------------------------
                //! copy the values into the arrays
                //! --------------------------------
                if(!isSuppressed && (isEdgeMeshSizeOK || isNumberOfDivisionsOK) && isScopeOK)
                {
                    //!cout<<"SimulationManager::transferMeshNodes()->____copying the values for edges____"<<endl;
                    for(int n=0;n<vecPairs.size();n++)
                    {
                        mySimulationDataBase->MapOfElementSizeOnEdge.setValue(vecPairs[n].first, vecPairs[n].second, edgeMeshSize);
                        mySimulationDataBase->MapOfNumberOfDivisionOnEdge.setValue(vecPairs[n].first, vecPairs[n].second, numberOfDivisions);
                        mySimulationDataBase->MapOfSizingTypeOnEdge.setValue(vecPairs[n].first, vecPairs[n].second,typeOfSizing);
                        mySimulationDataBase->MapOfIsEdgeModified.setValue(vecPairs[n].first, vecPairs[n].second, true);
                    }
                }
            }
                break;

            case SimulationNodeClass::nodeType_meshVertexSize:
            {
                //!cout<<"SimulationManager::transferMeshNodes()->____vertex sizing control found____"<<endl;
                double vertexMeshSize;
                double pinball;
                bool isSuppressed = true;
                bool isVertexMeshSizeOK = false;
                bool isScopeOK = false;
                bool isPinballOK = true;    //! pinball is always ok

                //! retrieve the vector of properties
                QVector<QExtendedStandardItem*> propertyItems = meshControlNode->getPropertyItems();
                QVector<Property> props;

                for(QVector<QExtendedStandardItem*>::Iterator it = propertyItems.begin();it!=propertyItems.end();++it)
                {
                    QExtendedStandardItem* item = *it;
                    props.push_back(item->data(Qt::UserRole).value<Property>());
                }

                //! a vector of pairs (bodyNumber, vertexNumber)
                std::vector<pair<int,int>> vecPairs;
                for(int k=0;k<propertyItems.length();k++)
                {
                    if(props.at(k).getGroup()==Property::PropertyGroup_Scope)
                    {
                        QVector<GeometryTag> vecLocs;
                        if(props.at(k).getName()=="Tags")
                        {
                            vecLocs = props.at(k).getData().value<QVector<GeometryTag>>();
                            if(vecLocs.size()>0)
                            {
                                //!cout<<"SimulationManager::transferMeshNodes()->____property scope found; pairs: "<<vecLocs.size()<<"____"<<endl;
                                std::pair<int,int> aPair;
                                for(int i=0; i<vecLocs.size();i++)
                                {
                                    GeometryTag loc = vecLocs.at(i);
                                    aPair.first = loc.parentShapeNr;
                                    aPair.second = loc.subTopNr;
                                    vecPairs.push_back(aPair);
                                }
                                isScopeOK=true;
                            }
                        }
                    }
                    else if(props.at(k).getGroup()==Property::PropertyGroup_Definition)
                    {
                        if(props.at(k).getName()=="Suppressed")
                        {
                            //! the control is "suppressed", i.e. is not accounted for
                            Property::SuppressionStatus ss = props.at(k).getData().value<Property::SuppressionStatus>();
                            if(ss==Property::SuppressionStatus_Active)isSuppressed = false;
                            else isSuppressed = true;
                            //!cout<<"SimulationManager::transferMeshNodes()->____found property suppression status: "<<isSuppressed<<"____"<<endl;
                        }
                        else if(props.at(k).getName()=="Element size")
                        {
                            vertexMeshSize = props.at(k).getData().toDouble();
                            isVertexMeshSizeOK=true;
                        }
                        else if(props.at(k).getName()=="Pinball")
                        {
                            pinball = props.at(k).getData().toDouble();
                            isPinballOK = true;
                        }
                    }
                }

                //! --------------------------------
                //! copy the values into the arrays
                //! --------------------------------
                if(!isSuppressed && isVertexMeshSizeOK && isScopeOK && isPinballOK)
                {
                    //!cout<<"SimulationManager::transferMeshNodes()->____copying the values for vertexes____"<<endl;
                    for(int n=0;n<vecPairs.size();n++)
                    {
                        mySimulationDataBase->MapOfElementSizeOnVertex.setValue(vecPairs[n].first, vecPairs[n].second, vertexMeshSize);
                        mySimulationDataBase->MapOfIsVertexModified.setValue(vecPairs[n].first, vecPairs[n].second, true);
                        mySimulationDataBase->MapOfVertexPinball.setValue(vecPairs[n].first, vecPairs[n].second, pinball);
                    }
                }

            }
                break;

            case SimulationNodeClass::nodeType_meshMethod:
            {
                Property::meshEngine2D surfaceMesher;
                Property::meshEngine3D volumeMesher;
                bool useBRep;
                Property::meshEngine2D Tessellator;
                double angularDeflection;               //! property of the Tessellator
                double linearDeflection;                //! property of the Tessellator
                double minFaceSize;                     //! property of the Tessellator
                double maxFaceSize;                     //! property of the Tessellator
                bool isDefeaturingOn;
                bool isHealingOn;
                bool isSimplificationOn;
                int meshSimplificationBy;
                double meshSimplificationParameter;
                bool preserveBoundaryConditionsEdges;
                bool projectMeshPointsOntoGeometry;

                //! -------------------------------
                //! check if the control is Active
                //! -------------------------------
                bool isActive = false;
                Property::SuppressionStatus ss = meshControlNode->getPropertyValue<Property::SuppressionStatus>("Suppressed");
                if(ss==Property::SuppressionStatus_Active) isActive = true;

                //! ------------------------------
                //! check if the scope is defined
                //! ------------------------------
                bool isScopeOK = false;
                QVector<GeometryTag> vecLoc = meshControlNode->getPropertyValue<QVector<GeometryTag>>("Tags");
                if(!vecLoc.isEmpty()) isScopeOK = true;

                //! ---------------------------------
                //! "Active" is on and "Scope" is ok
                //! ---------------------------------
                if(isActive == true && isScopeOK == true)
                {
                    useBRep = meshControlNode->getPropertyValue<bool>("Patch conforming");
                    if(useBRep)
                    {
                        surfaceMesher = meshControlNode->getPropertyValue<Property::meshEngine2D>("Surface mesher");
                        volumeMesher = meshControlNode->getPropertyValue<Property::meshEngine3D>("Volume mesher");
                        Tessellator = Property::meshEngine2D_OCC_STL;   //! set but not used
                        angularDeflection = -1;                         //! set but not used
                        linearDeflection = -1;                          //! set but not used
                        minFaceSize = 0.174;                            //! set but not used
                        maxFaceSize = 100;                              //! set but not used
                        isDefeaturingOn = false;                        //! set but not used
                        isHealingOn = false;                            //! set but not used
                        isSimplificationOn = false;                     //! set but not used
                        meshSimplificationBy = -1;                      //! set but not used
                        meshSimplificationParameter = -1;               //! set but not used
                        preserveBoundaryConditionsEdges = false;        //! set but not used
                        projectMeshPointsOntoGeometry = false;          //! set but not used
                    }
                    else
                    {
                        //! ----------------------------------------------------------------------------------------
                        //! when the volume mesher "TetgenBR" is selected, the property "Surface mesher"
                        //! is removed from the interface. Initialize the property with Property::meshEngine2D_NULL
                        //! -----------------------------------------------------------------------------------------
                        if(meshControlNode->getPropertyItem("Surface mesher")!=Q_NULLPTR)
                        {
                            surfaceMesher = meshControlNode->getPropertyValue<Property::meshEngine2D>("Surface mesher");
                        }
                        else
                        {
                            surfaceMesher = Property::meshEngine2D_NULL;
                            //!cout<<"____setting a NULL surface mesher____"<<endl;
                        }

                        volumeMesher = meshControlNode->getPropertyValue<Property::meshEngine3D>("Volume mesher");
                        Tessellator = meshControlNode->getPropertyValue<Property::meshEngine2D>("Tessellator");

                        angularDeflection = meshControlNode->getPropertyValue<double>("Angular deflection");
                        linearDeflection = meshControlNode->getPropertyValue<double>("Linear deflection");

                        //! -----------------------------------------------------------------
                        //! if the volume mesher is "TetWild" the defeaturing controls
                        //! are not active (they are not necessary for that kind of mesher):
                        //! this case must be handled here by checking the existance of the
                        //! "Defeaturing" item
                        //! -----------------------------------------------------------------
                        if(meshControlNode->getPropertyItem("Defeaturing")!=Q_NULLPTR)
                        {
                            isDefeaturingOn = meshControlNode->getPropertyValue<bool>("Defeaturing");
                            if(isDefeaturingOn)
                            {
                                isHealingOn = meshControlNode->getPropertyValue<bool>("Healing");
                                isSimplificationOn = meshControlNode->getPropertyValue<bool>("Simplification");
                                if(isSimplificationOn)
                                {
                                    meshSimplificationBy = meshControlNode->getPropertyValue<int>("By");
                                    if(meshSimplificationBy==0) meshSimplificationParameter = meshControlNode->getPropertyValue<double>("Level");
                                    if(meshSimplificationBy==1) meshSimplificationParameter = meshControlNode->getPropertyValue<double>("Pair distance");
                                }
                                else
                                {
                                    meshSimplificationBy = -1;                      //! set but not used
                                    meshSimplificationParameter = -1;               //! set but not used
                                }
                            }
                            else
                            {
                                isHealingOn = false;                            //! set but not used
                                meshSimplificationBy = -1;                      //! set but not used
                                meshSimplificationParameter = -1;               //! set but not used
                            }

                            if(meshControlNode->getPropertyItem("Preserve boundary conditions edges")!=Q_NULLPTR)
                                preserveBoundaryConditionsEdges = meshControlNode->getPropertyValue<bool>("Preserve boundary conditions edges");
                            if(meshControlNode->getPropertyItem("Project points on geometry")!=Q_NULLPTR)
                                projectMeshPointsOntoGeometry = meshControlNode->getPropertyValue<bool>("Project points on geometry");
                        }
                    }

                    //! --------------------------
                    //! TetWild mesher parameters
                    //! --------------------------
                    int envelopeSizing = 0;     //! relative
                    double relativeEnvelopeSize = 0.0025;
                    double absoluteEnvelopeSize = 0.10;
                    if(meshControlNode->getPropertyItem("Envelope sizing")!=Q_NULLPTR)
                    {
                        envelopeSizing = meshControlNode->getPropertyValue<int>("Envelope sizing");
                        if(envelopeSizing==0)
                        {
                            relativeEnvelopeSize = meshControlNode->getPropertyValue<double>("Relative size");
                        }
                        else
                        {
                            absoluteEnvelopeSize = meshControlNode->getPropertyValue<double>("Absolute size");
                        }
                    }
                    int idealLenghtSizing = 0;              //! relative
                    double relativeIdealLength = 0.10;
                    double absoluteIdealLength = 100.0;
                    if(meshControlNode->getPropertyItem("Ideal length sizing")!=Q_NULLPTR)
                    {
                        idealLenghtSizing = meshControlNode->getPropertyValue<int>("Ideal length sizing");
                        if(idealLenghtSizing==0)
                        {
                            relativeIdealLength = meshControlNode->getPropertyValue<double>("Relative length");
                        }
                        else
                        {
                            absoluteIdealLength = meshControlNode->getPropertyValue<double>("Absolute length");
                        }
                    }

                    Property::meshOrder mo = meshControlNode->getPropertyValue<Property::meshOrder>("Mesh order");

                    //! --------------
                    //! run in memory
                    //! --------------
                    bool runInMemory = false;
                    switch(volumeMesher)
                    {
                    case Property::meshEngine3D_Tetgen:
                    case Property::meshEngine3D_Tetgen_BR:
                    case Property::meshEngine3D_TetWild:
                        runInMemory = meshControlNode->getPropertyValue<bool>("Run in memory");
                        break;
                    default: runInMemory = false; break;
                    }

                    //! -----------------------------
                    //! actually transfer parameters
                    //! -----------------------------
                    for(QVector<GeometryTag>::iterator it = vecLoc.begin(); it!=vecLoc.end(); ++it)
                    {
                        GeometryTag aLoc = *it;
                        int bodyIndex = aLoc.parentShapeNr;

                        mySimulationDataBase->mapOfUseBRep.insert(bodyIndex,useBRep);
                        mySimulationDataBase->ArrayOfMesh2DEngine.insert(bodyIndex,surfaceMesher);
                        mySimulationDataBase->ArrayOfMesh3DEngine.insert(bodyIndex,volumeMesher);
                        mySimulationDataBase->MapOfIsBodyDefeaturingOn.insert(bodyIndex,isDefeaturingOn);
                        mySimulationDataBase->MapOfIsBodyHealingOn.insert(bodyIndex,isHealingOn);
                        mySimulationDataBase->MapOfIsMeshSimplificationOn.insert(bodyIndex,isSimplificationOn);
                        mySimulationDataBase->MapOfMeshSimplificationBy.insert(bodyIndex,meshSimplificationBy);
                        mySimulationDataBase->MapOfBodyDefeaturingParameterValue.insert(bodyIndex,meshSimplificationParameter);
                        mySimulationDataBase->ArrayOfMeshOrder.insert(bodyIndex,mo);

                        mySimulationDataBase->mapOfDiscretizer.insert(bodyIndex,Tessellator);
                        mySimulationDataBase->mapOfAngularDeflection.insert(bodyIndex,angularDeflection);
                        mySimulationDataBase->mapOfLinearDeflection.insert(bodyIndex,linearDeflection);
                        mySimulationDataBase->mapOfMinFaceSize.insert(bodyIndex,minFaceSize);
                        mySimulationDataBase->mapOfMaxFaceSize.insert(bodyIndex,maxFaceSize);
                        mySimulationDataBase->MapOfIsMeshingRunningInMemory.insert(bodyIndex,runInMemory);

                        mySimulationDataBase->mapOfEnvelopeSizingType.insert(bodyIndex,envelopeSizing);
                        mySimulationDataBase->mapOfRelativeEnvelopeSize.insert(bodyIndex,relativeEnvelopeSize);
                        mySimulationDataBase->mapOfAbsoluteEnvelopeSize.insert(bodyIndex,absoluteEnvelopeSize);
                        mySimulationDataBase->mapOfIdealLengthSizingType.insert(bodyIndex,relativeIdealLength);
                        mySimulationDataBase->mapOfIdealLengthRelativeSize.insert(bodyIndex,relativeIdealLength);
                        mySimulationDataBase->mapOfIdealLengthAbsoluteSize.insert(bodyIndex,absoluteIdealLength);

                        mySimulationDataBase->mapOfFeaturePreserving.insert(bodyIndex,preserveBoundaryConditionsEdges);
                        mySimulationDataBase->mapOfGeometryCorrection.insert(bodyIndex,projectMeshPointsOntoGeometry);

                        //! ----------------------------
                        //! diagnostic - can be removed
                        //! ----------------------------
                        cout<<"/----------------------------------------------/"<<endl;
                        cout<<"/ Summary of the trasferred mesh parameters    /"<<endl;
                        cout<<"/----------------------------------------------/"<<endl;
                        cout<<" Use BRep: "<<(useBRep? "TRUE": "FALSE")<<endl;
                        cout<<" Surface mesher: ";
                        switch(surfaceMesher)
                        {
                        case Property::meshEngine2D_Netgen: cout<<"Netgen"<<endl; break;
                        case Property::meshEngine2D_Netgen_STL: cout<<"Netgen STL 2D"<<endl; break;
                        case Property::meshEngine2D_OCC_ExpressMesh: cout<<"Express mesh"<<endl; break;
                        }
                        cout<<" Volume mesher: ";
                        switch(volumeMesher)
                        {
                        case Property::meshEngine3D_Netgen: cout<<"Netgen"<<endl; break;
                        case Property::meshEngine3D_Netgen_STL: cout<<"Netgen STL"<<endl; break;
                        case Property::meshEngine3D_Tetgen: cout<<"Tetgen"<<endl; break;
                        case Property::meshEngine3D_Tetgen_BR: cout<<"Tetgen BR"<<endl; break;
                        case Property::meshEngine3D_TetWild: cout<<"TetWild"<<endl; break;
                        }

                        cout<<" Mesh order: "<<(mo==Property::meshOrder_First? "First":"Second")<<endl;
                        cout<<" Use BRep: "<<(useBRep==true? "TRUE":"FALSE")<<endl;
                        cout<<" Surface Tessellator: "<<(Tessellator==Property::meshEngine2D_OCC_STL? "Standard STL":"Express mesh")<<"____"<<endl;
                        cout<<" Angular deflection: "<<angularDeflection<<endl;
                        cout<<" Linear deflection: "<<linearDeflection<<endl;
                        cout<<" Healing: "<<(isHealingOn==true? "ON":"OFF")<<endl;
                        cout<<" Defeaturing: "<<(isDefeaturingOn==true? "ON":"OFF")<<endl;
                        cout<<" By: ";

                        switch(meshSimplificationBy)
                        {
                        case 0: cout<<"Level: "<<meshSimplificationParameter<<endl; break;
                        case 1: cout<<"pair distance: "<<meshSimplificationParameter<<endl; break;
                        default: cout<<meshSimplificationBy<<endl;
                        }

                        cout<<" Preferred mode for meshing: "<<(runInMemory==true? "in memory":"on disk")<<endl;
                        cout<<"/----------------------------------------------/"<<endl;
                    }
                }
            }
                break;

            case SimulationNodeClass::nodeType_meshPrismaticLayer:
            {
                //! -------------------------------
                //! check if the control is Active
                //! -------------------------------
                bool isActive = false;
                Property::SuppressionStatus ss = meshControlNode->getPropertyValue<Property::SuppressionStatus>("Suppressed");
                if(ss==Property::SuppressionStatus_Active) isActive = true;
                //! ------------------------------
                //! check if the scope is defined
                //! ------------------------------
                bool isScopeOK = false;
                QVector<GeometryTag> vecLoc = meshControlNode->getPropertyValue<QVector<GeometryTag>>("Tags");
                if(!vecLoc.isEmpty()) isScopeOK = true;

                QVector<GeometryTag> boundaryVecLoc = meshControlNode->getPropertyValue<QVector<GeometryTag>>("Boundary tags");
                if(!boundaryVecLoc.isEmpty()) isScopeOK = true;

                if(isActive == true && isScopeOK == true)
                {
                    cout<<"_____transferring prismatic layer parameters to database____"<<endl;

                    //! ----------------------------------------------------
                    //! retrieve the prismatic mesh parameters from the GUI
                    //! ----------------------------------------------------
                    prismaticLayerParameters p;

                    p.typeOfSizing = meshControlNode->getPropertyValue<prismaticLayer_sizing>("Options");
                    switch(p.typeOfSizing)
                    {
                    case prismaticLayer_sizing_FirstLayerThickness: p.firstLayerThickness = meshControlNode->getPropertyValue<double>("First layer height"); break;
                    case prismaticLayer_sizing_TotalThickness: p.totalThickness = meshControlNode->getPropertyValue<double>("Total thickness"); break;
                    }

                    p.NbLayers = meshControlNode->getPropertyValue<int>("Number of layers");
                    p.expRatio = meshControlNode->getPropertyValue<double>("Expansion ratio");
                    p.lockBoundary = meshControlNode->getPropertyValue<bool>("Lock boundary");
                    p.checkSelfIntersections = meshControlNode->getPropertyValue<bool>("Check self intersections");
                    p.checkMutualIntersections = meshControlNode->getPropertyValue<bool>("Check mutual intersections");
                    p.boundaryMeshType = meshControlNode->getPropertyValue<int>("Boundary mesh type");
                    p.generationAlgorithm = meshControlNode->getPropertyValue<int>("Algorithm");
                    p.curvatureSensitivityForShrink = meshControlNode->getPropertyValue<double>("Curvature sensitivity");
                    p.NbGuidingVectorSmoothingSteps = meshControlNode->getPropertyValue<int>("Guiding vectors smoothing steps");
                    p.NbLayerThicknessSmoothingSteps = meshControlNode->getPropertyValue<int>("Thickness smoothing steps");
                    p.curvatureSensitivityForGuidingVectorsSmoothing = meshControlNode->getPropertyValue<double>("Guiding vector smoothing - curvature sensitivity");
                    p.curvatureSensitivityForThicknessSmoothing = meshControlNode->getPropertyValue<double>("Thickness smoothing - curvature sensitivity");

                    //! ------------------------------------------------------------------------
                    //! retrieve the prismatic faces locations: cycle over "Scope" using "Tags"
                    //! ------------------------------------------------------------------------
                    for(int n=0; n<vecLoc.size(); n++)
                    {
                        const GeometryTag &curBodyLoc = vecLoc.at(n);
                        int bodyIndex = curBodyLoc.parentShapeNr;

                        //! ---------------------------------------
                        //! activate inflation on the current boby
                        //! ---------------------------------------
                        mySimulationDataBase->HasPrismaticFaces.insert(bodyIndex,true);

                        //! --------------------------------------------
                        //! cycle over "Boundary" using "Boundary tags"
                        //! --------------------------------------------
                        QList<int> prismaticFacesOnBody;
                        for(int i=0; i<boundaryVecLoc.size(); i++)
                        {
                            const GeometryTag &aLoc = boundaryVecLoc.at(i);
                            int bodyIndex2 = aLoc.parentShapeNr;
                            int faceIndex2 = aLoc.subTopNr;
                            if(bodyIndex==bodyIndex2) prismaticFacesOnBody<<faceIndex2;
                        }
                        mySimulationDataBase->prismaticFaces.insert(bodyIndex,prismaticFacesOnBody);
                        mySimulationDataBase->prismaticMeshParameters.insert(bodyIndex,p);
                    }

                    //! -------------------------------
                    //! begin sending console messages
                    //! -------------------------------
                    cout<<"\n@____TRANSFERRING PRISMATIC MESH PARAMETERS____@"<<endl;
                    for(QMap<int,QList<int>>::iterator it=mySimulationDataBase->prismaticFaces.begin(); it!=mySimulationDataBase->prismaticFaces.end(); ++it)
                    {
                        int bodyIndex = it.key();
                        QString msg("@____prismatic faces on body index: ");

                        const QList<int> &faceNbList = it.value();
                        for(int i=0; i<faceNbList.length()-1; i++)
                        {
                            msg.append(QString("%1 ,").arg(faceNbList.at(i)));
                        }
                        msg.append(QString("%1}____@").arg(faceNbList.last()));

                        const prismaticLayerParameters &pl = mySimulationDataBase->prismaticMeshParameters.value(bodyIndex);
                        cout<<"@____number of layers: "<<pl.NbLayers<<"____@"<<endl;

                        switch(pl.typeOfSizing)
                        {
                        case prismaticLayer_sizing_FirstLayerThickness:
                            cout<<"@____type of sizing: first layer thickness____@"<<endl;
                            cout<<"@____first layer height: "<<pl.firstLayerThickness<<"____@"<<endl;
                            break;
                        case prismaticLayer_sizing_TotalThickness:
                            cout<<"@____type of sizing: total thickness____@"<<endl;
                            cout<<"@____total thickness: "<<pl.totalThickness<<"____@"<<endl;
                            break;
                        }

                        cout<<"@____expansion ratio: "<<pl.expRatio<<"____@"<<endl;
                        cout<<"@____boundary locked: "<<(pl.lockBoundary==true? "YES":"NO")<<"____@"<<endl;
                        cout<<"@____algorithm "<<pl.generationAlgorithm<<"____@"<<endl;
                        cout<<"@____boundary mesh type: "<<(pl.boundaryMeshType==0? "Hybrid":"Tetrahedral")<<"____@"<<endl;
                        //! ------------------------
                        //! end of console messages
                        //! ------------------------
                    }
                    cout<<"@__________END OF TRANSFER PROCESS_____________@\n"<<endl;
                }
            }
                break;
            }
        }
    }
}

//! ------------------------------------------------------------------
//! function: handleMeshItemChange
//! details:  if the content of a mesh control is changed by the user
//!           the mesh(es) to which the control is(are) applied must
//!           become invalid
//! ------------------------------------------------------------------
void SimulationManager::handleMeshItemChange(QStandardItem *item)
{
    cout<<"SimulationManager::handleMeshItemChange()->____function called____"<<endl;
    Q_UNUSED(item)

    QExtendedStandardItem *curItem = static_cast<QExtendedStandardItem*>(myModel->itemFromIndex(myTreeView->currentIndex()));
    SimulationNodeClass *curNode = curItem->data(Qt::UserRole).value<SimulationNodeClass*>();

    //! ----------------------------------------------------
    //! if a change has occurred into a "Mesh metric" item,
    //! the current mesh should not be affected
    //! ----------------------------------------------------
    if(curNode->getType() == SimulationNodeClass::nodeType_meshMeshMetric) return;

    QVector<GeometryTag> vecLocs = curNode->getPropertyValue<QVector<GeometryTag>>("Tags");
    std::vector<int> parentShapes;
    for(int k=0; k<vecLocs.size();k++)
    {
        const GeometryTag &loc = vecLocs.at(k);
        int n = loc.parentShapeNr;
        if(std::find(parentShapes.begin(),parentShapes.end(),n)==parentShapes.end()) parentShapes.push_back(n);
    }
    emit requestMeshInvalidate(parentShapes);
}

//! --------------------------------------------------------------------
//! function: updateRemotePointAbsCoordinates
//! details:  enter the remote points root;
//!           the function updates the absolute coordinates of all
//!           the defined remote points (introduced because a change
//!           in the definition of a coordinate system can in principle
//!           affect the definition of several remote points)
//! --------------------------------------------------------------------
void SimulationManager::updateRemotePointAbsCoordinates()
{
    cout<<"SimulationManager::updateRemotePointAbsCoordinates()->____function called____"<<endl;
    QExtendedStandardItem *itemRemotePointRoot = this->getTreeItem(SimulationNodeClass::nodeType_remotePointRoot);
    if(itemRemotePointRoot==Q_NULLPTR) return;

    //! ------------------------------------------------------------------------------
    //! start from row = 1, not row = 0, since it exists an hidden remote point which
    //! is used as "Select from list" in the interface combobox
    //! ------------------------------------------------------------------------------
    for(int row = 1; row<itemRemotePointRoot->rowCount(); row++)
    {
        cout<<"SimulationManager::updateRemotePointAbsCoordinates()->____updating Remote point root child N: "<<row<<"____"<<endl;
        QStandardItem *curItem = itemRemotePointRoot->child(row,0);
        SimulationNodeClass *curNode = curItem->data(Qt::UserRole).value<SimulationNodeClass*>();

        //! ----------------------------------------------------------------
        //! retrieve the system of reference the remote point is defined in
        //! ----------------------------------------------------------------
        void *p = curNode->getPropertyItem("Coordinate system")->data(Qt::UserRole).value<Property>().getData().value<void*>();
        QExtendedStandardItem *item = (QExtendedStandardItem*)p;
        SimulationNodeClass *CS = item->data(Qt::UserRole).value<SimulationNodeClass*>();

        double xP = curNode->getPropertyValue<double>("X coordinate");
        double yP = curNode->getPropertyValue<double>("Y coordinate");
        double zP = curNode->getPropertyValue<double>("Z coordinate");

        if(CS->getType()!=SimulationNodeClass::nodeType_coordinateSystem_global)
        {
            //! -------------------------------------------------------------------------------
            //! the remote point is defined with respect to a user defined system of reference
            //! -------------------------------------------------------------------------------
            QVector<double> CS_origin = CS->getPropertyValue<QVector<double>>("Base origin");
            QVector<double> X_axisData = CS->getPropertyValue<QVector<double>>("X axis data");
            QVector<double> Y_axisData = CS->getPropertyValue<QVector<double>>("Y axis data");
            QVector<double> Z_axisData = CS->getPropertyValue<QVector<double>>("Z axis data");

            double xO = CS_origin.at(0);
            double yO = CS_origin.at(1);
            double zO = CS_origin.at(2);

            double xnew = xO + xP*X_axisData.at(0) + yP*Y_axisData.at(0) + zP*Z_axisData.at(0);
            double ynew = yO + yP*X_axisData.at(1) + yP*Y_axisData.at(1) + zP*Z_axisData.at(1);
            double znew = zO + zP*X_axisData.at(2) + yP*Y_axisData.at(2) + zP*Z_axisData.at(2);

            QVariant data;
            data.setValue(xnew);
            curNode->replaceProperty("X abs coordinate",Property("X abs coordinate",data,Property::PropertyGroup_Scope));
            data.setValue(ynew);
            curNode->replaceProperty("Y abs coordinate",Property("Y abs coordinate",data,Property::PropertyGroup_Scope));
            data.setValue(znew);
            curNode->replaceProperty("Z abs coordinate",Property("Z abs coordinate",data,Property::PropertyGroup_Scope));
        }
        else
        {
            //! ---------------------------------------------------------------------------
            //! the remote point is defined with respect to the "Global coordinate system"
            //! ---------------------------------------------------------------------------
            QVariant data;
            data.setValue(xP);
            curNode->replaceProperty("X abs coordinate",Property("X abs coordinate",data,Property::PropertyGroup_Scope));
            data.setValue(yP);
            curNode->replaceProperty("Y abs coordinate",Property("Y abs coordinate",data,Property::PropertyGroup_Scope));
            data.setValue(zP);
            curNode->replaceProperty("Z abs coordinate",Property("Z abs coordinate",data,Property::PropertyGroup_Scope));
        }

        //! ------------------------------------------------------------
        //! immediately update the position of the marker in the viewer
        //! ------------------------------------------------------------
        markerBuilder::addMarker(curNode,mySimulationDataBase);
        emit requestHideAllMarkers(false);
        this->displayMarker();
    }
}

//! ------------------------------------------------------------
//! function: handleItemChange
//! details:  the "item" argument is an item of the node model!
//! ------------------------------------------------------------
void SimulationManager::handleItemChange(QStandardItem *item)
{
    cout<<"*------------------------------------------------*"<<endl;
    cout<<"*-  handleItemChange()->____function called____ -*"<<endl;
    cout<<"*------------------------------------------------*"<<endl;

    //! ------------------------------------------------------------------------------------
    //! "Tags", "Master tags", "Slave tags" changes are changed by the widget
    //! "Shape selector", by consequence they must be excluded (otherwise this function
    //! would be called two times both for the change in the "Shape selector", both for the
    //! change in the "Tags"/"Master tags"/"Slave tags" property
    //! ------------------------------------------------------------------------------------
    QString propertyName = item->data(Qt::UserRole).value<Property>().getName();

    cout<<"____changed property: \""<<propertyName.toStdString()<<"\"____"<<endl;
    if(propertyName =="Source file path")
    {
        cout<<"____clear the geometry data base____"<<endl;
        return;
    }

    if(propertyName== "Tags" ||  propertyName== "Master tags" || propertyName =="Slave tags" || propertyName== "Boundary tags")
        return;

    //! -------------------------
    //! simulation item and node
    //! -------------------------
    QExtendedStandardItem *curItem = static_cast<QExtendedStandardItem*>(myModel->itemFromIndex(myTreeView->currentIndex()));
    SimulationNodeClass *curNode = curItem->data(Qt::UserRole).value<SimulationNodeClass*>();
    SimulationNodeClass::nodeType family = curNode->getFamily();
    SimulationNodeClass::nodeType type = curNode->getType();

    //! ------------------------------------------------
    //! parse the item after the change
    //! reconnection is at the end of this function [*]
    //! ------------------------------------------------
    disconnect(curNode->getModel(),SIGNAL(itemChanged(QStandardItem*)),this,SLOT(handleItemChange(QStandardItem*)));
    //curNode->getModel()->blockSignals(true);

    switch(family)
    {
    case SimulationNodeClass::nodeType_geometry:
    {
        if(type==SimulationNodeClass::nodeType_pointMass)
        {
            if(propertyName=="X coordinate" || propertyName=="Y coordinate" || propertyName=="Z coordinate")
            {
                markerBuilder::addMarker(curNode,mySimulationDataBase);
                emit requestHideAllMarkers(false);
                this->displayMarker();
            }
        }
    }
        break;
    case SimulationNodeClass::nodeType_StructuralAnalysisSolution:
    case SimulationNodeClass::nodeType_thermalAnalysisSolution:
    case SimulationNodeClass::nodeType_postObject:
    {
        if(propertyName =="Type ")
        {
            int componentType = curNode->getPropertyValue<int>("Type ");
            QString newName;

            switch(type)
            {
            case SimulationNodeClass::nodeType_solutionStructuralNodalDisplacement:
            {
                switch(componentType)
                {
                case 0: newName = "Total displacement"; break;
                case 1: newName = "Directional displacement X"; break;
                case 2: newName = "Directional displacement Y"; break;
                case 3: newName = "Directional displacement Z"; break;
                }
            }
                break;
            case SimulationNodeClass::nodeType_solutionStructuralContact:
            {
                switch(componentType)
                {
                case 0: newName = "Contact pressure"; break;
                case 1: newName = "Frictional stress"; break;
                case 2: newName = "Penetration"; break;
                case 3: newName = "Sliding"; break;
                }
            }
                break;
            case SimulationNodeClass::nodeType_solutionStructuralNodalForces:
            {
                switch(componentType)
                {
                case 0: newName = "Total"; break;
                case 1: newName = "X direction"; break;
                case 2: newName = "Y direction"; break;
                case 3: newName = "Z direction"; break;
                }
            }
                break;
            case SimulationNodeClass::nodeType_solutionStructuralStress:
            {
                switch(componentType)
                {
                case 0: newName = "Equivalent stress"; break;
                case 1: newName = "Stress intensity"; break;
                case 2: newName = "Maximum principal stress"; break;
                case 3: newName = "Middle principal stress"; break;
                case 4: newName = "Minimum principal stress"; break;
                case 5: newName = "Normal stress X"; break;
                case 6: newName = "Normal stress Y"; break;
                case 7: newName = "Normal stress Z"; break;
                case 8: newName = "Shear stress XY"; break;
                case 9: newName = "Shear stress YZ"; break;
                case 10: newName = "Shear stress XZ"; break;
                }
            }
                break;
            case SimulationNodeClass::nodeType_solutionStructuralTotalStrain:
            case SimulationNodeClass::nodeType_solutionStructuralMechanicalStrain:
            case SimulationNodeClass::nodeType_solutionStructuralThermalStrain:
            {
                switch(componentType)
                {
                case 0: newName = "Equivalent strain"; break;
                case 1: newName = "Strain intensity"; break;
                case 2: newName = "Maximum principal strain"; break;
                case 3: newName = "Middle principal strain"; break;
                case 4: newName = "Minimum principal strain"; break;
                }
            }
                break;
            }

            //! -----------------------------------------------------------------
            //! change the node name and the item display name (Qt::DisplayRole)
            //! -----------------------------------------------------------------
            curNode->setName(newName);
            QVariant data;
            data.setValue(newName);
            static_cast<QStandardItemModel*>(myTreeView->model())->setData(myTreeView->currentIndex(),data);
        }
        if(propertyName =="Scale type" || propertyName == "# intervals" || propertyName =="Min" || propertyName =="Max")
        {
            //! ------------------------
            //! reorganize the switches
            //! ------------------------
            myDetailViewer->handleColorBoxScaleChanged();

            //! ------------------------------------------
            //! get the properties of the color box scale
            //! ------------------------------------------
            int scaleType = curNode->getPropertyValue<int>("Scale type");
            double minValue, maxValue;
            if(scaleType == 1)
            {
                minValue = curNode->getPropertyValue<double>("Min");
                maxValue = curNode->getPropertyValue<double>("Max");
            }
            int NbIntervals = curNode->getPropertyValue<int>("# intervals");

            //! ----------------------------------------
            //! update the post object
            //! (first check if the post object exists)
            //! ----------------------------------------
            QStandardItem *itemPostObject = curNode->getPropertyItem("Post object");
            if(itemPostObject!=Q_NULLPTR)
            {
                postObject thePostObject = itemPostObject->data(Qt::UserRole).value<Property>().getData().value<postObject>();
                emit requestHideSingleResult(thePostObject);
                myPostEngine->updateResultScale(thePostObject,scaleType,minValue,maxValue,NbIntervals);

                //! --------------------------
                //! update the property value
                //! --------------------------
                QVariant data;
                data.setValue(thePostObject);
                Property prop_postObject("Post object",data,Property::PropertyGroup_GraphicObjects);
                curNode->replaceProperty("Post object",prop_postObject);

                //! ----------------
                //! show the result
                //! ----------------
                emit requestDisplayResult(thePostObject);
            }
        }
        if(propertyName =="Mapping")
        {
            int val = curNode->getPropertyValue<int>("Mapping");
            cout<<"____mapping: "<<val<<"____"<<endl;
            if(curNode->getPropertyItem("Post object")!=Q_NULLPTR)
            {
                postObject thePostObject = curNode->getPropertyValue<postObject>("Post object");
                emit requestHideSingleResult(thePostObject);
            }
        }
        if(propertyName =="By")
        {
            //! ------------------------
            //! reorganize the switches
            //! ------------------------
            myDetailViewer->handleByChanged();

            if(curNode->getPropertyItem("Post object")!=Q_NULLPTR)
            {
                const postObject &curPostObject = curNode->getPropertyValue<postObject>("Post object");
                emit requestHideSingleResult(curPostObject);
                curNode->removeProperty("Post object");
                this->callPostEngineEvaluateResult_private(curItem,false);
            }
        }
        if(propertyName =="Display time" || propertyName=="Set number")
        {
            if(curNode->getPropertyItem("Post object")!=Q_NULLPTR)
            {
                const postObject &curPostObject = curNode->getPropertyValue<postObject>("Post object");
                emit requestHideSingleResult(curPostObject);
                curNode->removeProperty("Post object");
                this->callPostEngineEvaluateResult_private(curItem,false);
            }
        }
    }
        break;

    case SimulationNodeClass::nodeType_coordinateSystems:
    {
        if(type==SimulationNodeClass::nodeType_coordinateSystem)
        {
            cout<<"____handling coordinate systems____"<<endl;

            //! ------------------------------------------------------------
            //! update the definition of the current clip plane, if defined
            //! to many calls... correct... to do...
            //! ------------------------------------------------------------
            //curNode->getModel()->blockSignals(true);
            clipTool *theClipTool = static_cast<clipTool*>(tools::getWidgetByName("clipTool"));
            theClipTool->updateCSDataByExternalCSChange(curItem);
            //curNode->getModel()->blockSignals(false);

            //! ----------------------------------------------------------------
            //! update the remote points, which can contain in their definition
            //! one of the defined coordinate systems
            //! ----------------------------------------------------------------
            QStandardItem *itemRemotePointRoot = this->getTreeItem(SimulationNodeClass::nodeType_remotePointRoot);
            if(itemRemotePointRoot!=NULL)
            {
                //! ------------------------------------------------------
                //! disconnect/block signals. Reconnection at the end [*]
                //! starts from row = 1, jumping row = 0, since there
                //! is the hidden "Select from list" face remote point
                //! ------------------------------------------------------
                for(int row = 1; row<itemRemotePointRoot->rowCount(); row++)
                {
                    QStandardItem *curItem = itemRemotePointRoot->child(row,0);
                    SimulationNodeClass *curNode = curItem->data(Qt::UserRole).value<SimulationNodeClass*>();
                    curNode->getModel()->blockSignals(true);
                }

                //! ---------------------------------------------------
                //! updating absolute coordinates of the remote points
                //! ---------------------------------------------------
                this->updateRemotePointAbsCoordinates();

                //! ---------------------
                //! unnblock signals [*]
                //! ---------------------
                for(int row = 1; row<itemRemotePointRoot->rowCount(); row++)
                {
                    QStandardItem *curItem = itemRemotePointRoot->child(row,0);
                    SimulationNodeClass *curNode = curItem->data(Qt::UserRole).value<SimulationNodeClass*>();
                    curNode->getModel()->blockSignals(false);
                }
            }

            //! ---------------------------------------------------------------------------
            //! update also the definition of:
            //! - Remote force[s] and Remote displacement[s] => update the reference point
            //!   used for the definition of the reference node
            //! ---------------------------------------------------------------------------

            //! ------------------------------------------------------------------------------------
            //! block/disconnect from SimulationManager::handleItemChange() all the children of all
            //! the simulation roots, in order to avoid multiple calls of this function, causing
            //! application crash. Then unblock/reconnect/at the end [*]
            //! ------------------------------------------------------------------------------------
            QList<QExtendedStandardItem*> listOfSimulationRootItems =
                    this->getAllTreeItemOfType(SimulationNodeClass::nodeType_structuralAnalysis);

            for(int n = 0; n<listOfSimulationRootItems.length(); n++)
            {
                cerr<<"____simulation root n. "<<n<<"____"<<endl;
                QStandardItem* aStructuralAnalysisRootItem = listOfSimulationRootItems.at(n);

                //! --------------------------------------------------------------
                //! "-1" makes evident the fact that we stop at before "Solution"
                //! starts from row = 1, jumping over "Analysis settings"
                //! --------------------------------------------------------------
                for(int row = 1; row<=aStructuralAnalysisRootItem->rowCount()-1; row++)
                {
                    QStandardItem * curItem = aStructuralAnalysisRootItem->child(row,0);
                    SimulationNodeClass *curNode = curItem->data(Qt::UserRole).value<SimulationNodeClass*>();
                    curNode->getModel()->blockSignals(true);
                }
            }

            //! -----------------
            //! updating routine
            //! -----------------
            //QStandardItem *itemRemotePointRoot = this->getTreeItem(SimulationNodeClass::nodeType_remotePointRoot);
            for(int n = 0; n<listOfSimulationRootItems.length(); n++)
            {
                QStandardItem* aStructuralAnalysisRootItem = listOfSimulationRootItems.at(n);

                //! ----------------------------------------------------
                //! start from row = 1, jumping over "Analysis settings"
                //! also jump over "Solution" item
                //! ----------------------------------------------------
                for(int row = 1; row<=aStructuralAnalysisRootItem->rowCount()-1; row++)
                {
                    QStandardItem * curItem = aStructuralAnalysisRootItem->child(row,0);
                    SimulationNodeClass *curNode = curItem->data(Qt::UserRole).value<SimulationNodeClass*>();
                    SimulationNodeClass::nodeType nodeType = curNode->getType();
                    if(nodeType == SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteForce ||
                            nodeType == SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteDisplacement ||
                            nodeType == SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteRotation)
                    {
                        cout<<"____\"Remote force\" or \"Remote displacement\" found____"<<endl;
                        if(curNode->getPropertyValue<Property::ScopingMethod>("Scoping method")==Property::ScopingMethod_RemotePoint)
                        {
                            cout<<"____\"Scoping method\" defined through \"Remote point\" found____"<<endl;
                            for(int k=0; k<itemRemotePointRoot->rowCount(); k++)
                            {
                                cout<<"____comparing____"<<endl;

                                //! ---------------------------------------------------------------------------------
                                //! compare the "Remote point" withinn the scope definition of the "Remote force" or
                                //! the "Remote displacement" with each remote point present in the "Remote points"
                                //! branch
                                //! ---------------------------------------------------------------------------------
                                QStandardItem *curItemRemotePointInRemotePointBranch = itemRemotePointRoot->child(k,0);
                                void *p = curNode->getPropertyItem("Remote points")->data(Qt::UserRole).value<Property>().getData().value<void*>();
                                QStandardItem *curItemRemotePointInItem = (QStandardItem*) p;
                                if(curItemRemotePointInRemotePointBranch==curItemRemotePointInItem)
                                {
                                    cout<<"____\"Remote point\" belonging to the remote point list found____"<<endl;

                                    //! update the remote force or remote displacement reference point
                                    cout<<"____updating the reference point____"<<endl;
                                    SimulationNodeClass *nodeRP = curItemRemotePointInRemotePointBranch->data(Qt::UserRole).value<SimulationNodeClass*>();
                                    double xRefPoint = nodeRP->getPropertyValue<double>("X abs coordinate");
                                    double yRefPoint = nodeRP->getPropertyValue<double>("Y abs coordinate");
                                    double zRefPoint = nodeRP->getPropertyValue<double>("Z abs coordinate");
                                    QVector<double> refPoint;
                                    refPoint.push_back(xRefPoint);
                                    refPoint.push_back(yRefPoint);
                                    refPoint.push_back(zRefPoint);
                                    QVariant data;
                                    data.setValue(refPoint);
                                    Property prop_refPoint("Reference point",data,Property::PropertyGroup_Hidden);
                                    curNode->replaceProperty("Reference point",prop_refPoint);
                                    break;
                                }
                                else continue;
                            }
                        }
                    }
                }
            }

            //! ----------------------
            //! unblock/reconnect [*]
            //! ----------------------
            for(int n = 0; n<listOfSimulationRootItems.length(); n++)
            {
                QStandardItem* aStructuralAnalysisRootItem = listOfSimulationRootItems.at(n);
                for(int row = 1; row<=aStructuralAnalysisRootItem->rowCount()-1; row++)
                {
                    QStandardItem * curItem = aStructuralAnalysisRootItem->child(row,0);
                    SimulationNodeClass *curNode = curItem->data(Qt::UserRole).value<SimulationNodeClass*>();
                    curNode->getModel()->blockSignals(false);
                }
            }
        }
    }
        break;

    case SimulationNodeClass::nodeType_remotePointRoot:
    {
        cout<<"SimulationManager::handleItemChange()->____handling remote point____"<<endl;

        //! ------------------------------------------------------
        //! block signals in order to avoid multiple calls of
        //! "SimulationManager::handleItemChange()" causing crash
        //! Unlock is at the end [*]
        //! ------------------------------------------------------
        QStandardItem *itemRemotePointRoot = this->getTreeItem(SimulationNodeClass::nodeType_remotePointRoot);
        for(int row = 1; row<itemRemotePointRoot->rowCount(); row++)
        {
            QStandardItem *curItem = itemRemotePointRoot->child(row,0);
            SimulationNodeClass *curNode = curItem->data(Qt::UserRole).value<SimulationNodeClass*>();
            curNode->getModel()->blockSignals(true);
        }

        //! ---------------------------------------------------
        //! updating absolute coordinates of the remote points
        //! ---------------------------------------------------
        this->updateRemotePointAbsCoordinates();

        //! -----------------------------------------
        //! update the position of the sphere marker
        //! -----------------------------------------
        double x = curNode->getPropertyValue<double>("X abs coordinate");
        double y = curNode->getPropertyValue<double>("Y abs coordinate");
        double z = curNode->getPropertyValue<double>("Z abs coordinate");

        gp_Pnt P(x,y,z);
        emit requestDisplaySphericalMarker(P);

        //! -----------------------------------------------------------------------
        //! "Coupling" from "Kinematic"/"Distributed" to "Distributed"/"Kinematic"
        //! Remove/add rotational DOFs
        //! -----------------------------------------------------------------------
        if(propertyName=="Coupling") myDetailViewer->handleCouplingChanged();
        //! -------------------------------------------------------------------------------------
        //! "DOFs selection" from "Manual"/"Program controlled" to "Program controlled"/"Manual"
        //! Remove/add switches
        //! -------------------------------------------------------------------------------------
        if(propertyName=="DOFs selection") myDetailViewer->handleDOFselectorChange();

        //! ---------------------
        //! unnblock signals [*]
        //! ---------------------
        for(int row = 1; row<itemRemotePointRoot->rowCount(); row++)
        {
            QStandardItem *curItem = itemRemotePointRoot->child(row,0);
            SimulationNodeClass *curNode = curItem->data(Qt::UserRole).value<SimulationNodeClass*>();
            curNode->getModel()->blockSignals(false);
        }

        //! --------------------------------------------------------
        //! if a child of the remote point root is changed
        //! all the related boundary conditions (remote forces,
        //! remote displacements) must be updated. The properties
        //! to be update is the "Reference point" and all the group
        //! of the "Advanced" properties, if the "Scoping method"
        //! is defined through the "Remote points" property
        //! --------------------------------------------------------

        //! ------------------------------------------------------------------------------------
        //! block/disconnect from SimulationManager::handleItemChange() all the children of all
        //! the simulation roots, in order to avoid multiple calls of this function, causing
        //! application crash. Then unblock/reconnect/at the end [*]
        //! ------------------------------------------------------------------------------------
        QList<QExtendedStandardItem*> listOfSimulationRootItems =
                this->getAllTreeItemOfType(SimulationNodeClass::nodeType_structuralAnalysis);

        for(int n = 0; n<listOfSimulationRootItems.length(); n++)
        {
            cout<<"SimulationManager::handleItemChange()->____simulation root n. "<<n<<"____"<<endl;
            QStandardItem* aStructuralAnalysisRootItem = listOfSimulationRootItems.at(n);

            //! --------------------------------------------------------------
            //! "-1" makes evident the fact that we stop at before "Solution"
            //! "1" that we jump over "Analysis settings"
            //! --------------------------------------------------------------
            for(int row = 1; row<=aStructuralAnalysisRootItem->rowCount()-1; row++)
            {
                QStandardItem * curItem = aStructuralAnalysisRootItem->child(row,0);
                SimulationNodeClass *curNode = curItem->data(Qt::UserRole).value<SimulationNodeClass*>();
                curNode->getModel()->blockSignals(true);
            }
        }

        //! -----------------
        //! updating routine
        //! -----------------
        //QStandardItem *itemRemotePointRoot = this->getTreeItem(SimulationNodeClass::nodeType_remotePointRoot);
        for(int n = 0; n<listOfSimulationRootItems.length(); n++)
        {
            QStandardItem* aStructuralAnalysisRootItem = listOfSimulationRootItems.at(n);
            for(int row = 0; row<=aStructuralAnalysisRootItem->rowCount()-1; row++)
            {
                QStandardItem * curItem = aStructuralAnalysisRootItem->child(row,0);
                SimulationNodeClass *curNode = curItem->data(Qt::UserRole).value<SimulationNodeClass*>();
                SimulationNodeClass::nodeType nodeType = curNode->getType();
                if(nodeType == SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteForce ||
                        nodeType == SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteDisplacement ||
                        nodeType == SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteRotation)
                {
                    cout<<"SimulationManager::handleItemChange()->____\"Remote force\" or \"Remote displacement\" found____"<<endl;
                    if(curNode->getPropertyValue<Property::ScopingMethod>("Scoping method")==Property::ScopingMethod_RemotePoint)
                    {
                        cout<<"SimulationManager::handleItemChange()->____\"Scoping method\" defined through \"Remote point\" found____"<<endl;
                        for(int k=0; k<itemRemotePointRoot->rowCount(); k++)
                        {
                            cout<<"SimulationManager::handleItemChange()->____comparing cycle k: "<<k<<"____"<<endl;

                            //! ---------------------------------------------------------------------------------
                            //! compare the "Remote point" withinn the scope definition of the "Remote force" or
                            //! the "Remote displacement" with each remote point present in the "Remote points"
                            //! branch
                            //! ---------------------------------------------------------------------------------
                            QStandardItem *curItemRemotePointInRemotePointBranch = itemRemotePointRoot->child(k,0);
                            void *p = curNode->getPropertyItem("Remote points")->data(Qt::UserRole).value<Property>().getData().value<void*>();
                            QStandardItem *curItemRemotePointInItem = (QStandardItem*) p;
                            if(curItemRemotePointInRemotePointBranch==curItemRemotePointInItem)
                            {
                                cout<<"SimulationManager::handleItemChange()->____\"Remote point\" belonging to the remote point list found____"<<endl;

                                //! ---------------------------------------------------------------
                                //! update the remote force or remote displacement reference point
                                //! ---------------------------------------------------------------
                                cout<<"SimulationManager::handleItemChange()->____updating the reference point____"<<endl;
                                SimulationNodeClass *nodeRP = curItemRemotePointInRemotePointBranch->data(Qt::UserRole).value<SimulationNodeClass*>();
                                double xRefPoint = nodeRP->getPropertyValue<double>("X abs coordinate");
                                double yRefPoint = nodeRP->getPropertyValue<double>("Y abs coordinate");
                                double zRefPoint = nodeRP->getPropertyValue<double>("Z abs coordinate");
                                QVector<double> refPoint;
                                refPoint.push_back(xRefPoint);
                                refPoint.push_back(yRefPoint);
                                refPoint.push_back(zRefPoint);
                                QVariant data;
                                data.setValue(refPoint);
                                Property prop_refPoint("Reference point",data,Property::PropertyGroup_Hidden);
                                curNode->replaceProperty("Reference point",prop_refPoint);

                                //! ----------------------------------------------------------------
                                //! remove the "Advanced" properties and replace them with the ones
                                //! from the "Remote point"
                                //! ----------------------------------------------------------------
                                cout<<"SimulationManager::handleItemChange()->____removing \"Advanced properties\"____"<<endl;
                                QStandardItem *itemCoupling = curNode->getPropertyItem("Coupling");
                                QStandardItem *itemSeparatorAdvanced = itemCoupling->parent();
                                curNode->getModel()->removeRows(itemCoupling->row(),itemSeparatorAdvanced->rowCount(),itemSeparatorAdvanced->index());

                                cout<<"SimulationManager::handleItemChange()->____rebuilding \"Advanced properties\"____"<<endl;
                                QStandardItem *separatorAdvancedInRemotePointBranch = nodeRP->getPropertyItem("Coupling")->parent();
                                cout<<"SimulationManager::handleItemChange()->____number of property to read from \"Remote point\": "<<separatorAdvancedInRemotePointBranch->rowCount()<<"____"<<endl;

                                for(int row = 0; row<separatorAdvancedInRemotePointBranch->rowCount(); row++)
                                {
                                    QString name = separatorAdvancedInRemotePointBranch->child(row,0)->data(Qt::DisplayRole).toString();
                                    cout<<"SimulationManager::handleItemChange()->____copying property: "<<name.toStdString()<<"____"<<endl;
                                    QVariant data = separatorAdvancedInRemotePointBranch->child(row,1)->data(Qt::UserRole).value<Property>().getData();
                                    Property::PropertyGroup group = separatorAdvancedInRemotePointBranch->child(row,1)->data(Qt::UserRole).
                                            value<Property>().getGroup();
                                    Property prop(name,data,group);
                                    curNode->addProperty(prop,row);
                                }
                                break;
                            }
                            else continue;
                        }
                    }
                }
            }
        }

        //! ----------------------
        //! unblock/reconnect [*]
        //! ----------------------
        for(int n = 0; n<listOfSimulationRootItems.length(); n++)
        {
            QStandardItem* aStructuralAnalysisRootItem = listOfSimulationRootItems.at(n);
            for(int row = 0; row<=aStructuralAnalysisRootItem->rowCount()-1; row++)
            {
                QStandardItem * curItem = aStructuralAnalysisRootItem->child(row,0);
                SimulationNodeClass *curNode = curItem->data(Qt::UserRole).value<SimulationNodeClass*>();
                curNode->getModel()->blockSignals(false);
            }
        }
    }
        break;

    case SimulationNodeClass::nodeType_meshControl:
    {
        //! this "if" is because the mesh root should not be parsed
        if(type!=SimulationNodeClass::nodeType_meshControl)
        {
            this->handleMeshItemChange(curItem);
            this->changeColor();
        }
        else
        {
            //! at the moment all the changes occurring in the mesh root nodes are handled by
            //! the signal "globalMeshControlChanged()" emitted by the delegate: this in turn
            //! activates the DetailViewer which calls (always through signal/slot connection
            //! this->handleGlobalMeshControlChange(), which invalidates all the meshes.
            //! This is working, but in the future try to bring the simulation manager the control
            //! of the mesh root changes
        }
    }
        break;

    case SimulationNodeClass::nodeType_structuralAnalysis:
    case SimulationNodeClass::nodeType_thermalAnalysis:
    {
        switch(type)
        {
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteForce:
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteDisplacement:
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteRotation:
        {
            cout<<"SimulationManager::handleItemChange()->____handling coordinate \"Remote force\" or \"Remote displacement\" item____"<<endl;
            Property::ScopingMethod prop_ScopingMethod = curNode->getPropertyValue<Property::ScopingMethod>("Scoping method");
            switch(prop_ScopingMethod)
            {
            case Property::ScopingMethod_GeometrySelection:
            case Property::ScopingMethod_NamedSelection:
                break;

            case Property::ScopingMethod_RemotePoint:
            {
                //! ----------------------------------------------------------------
                //! update the reference point by reading the new coordinates from
                //! the remote point ("X abs absolute", ...)
                //! ----------------------------------------------------------------
                if(curNode->getPropertyItem("Remote points")!=NULL)
                {
                    void *p = curNode->getPropertyItem("Remote points")->data(Qt::UserRole).value<Property>().getData().value<void*>();

                    QExtendedStandardItem *itemRP = (QExtendedStandardItem*) p;
                    SimulationNodeClass *nodeRP = itemRP->data(Qt::UserRole).value<SimulationNodeClass*>();

                    double xRP = nodeRP->getPropertyValue<double>("X abs coordinate");
                    double yRP = nodeRP->getPropertyValue<double>("Y abs coordinate");
                    double zRP = nodeRP->getPropertyValue<double>("Z abs coordinate");

                    QVector<double> coords;
                    coords.push_back(xRP);
                    coords.push_back(yRP);
                    coords.push_back(zRP);

                    QVariant data;
                    data.setValue(coords);
                    Property prop_referencePoint("Reference point",data,Property::PropertyGroup_Hidden);
                    curNode->replaceProperty("Reference point",prop_referencePoint);
                }
            }
                break;
            }
        }
            break;
        }

        //! -----------------------------------------------------------------------
        //! "Coupling" from "Kinematic"/"Distributed" to "Distributed"/"Kinematic"
        //! Remove/add rotational DOFs
        //! -----------------------------------------------------------------------
        if(propertyName=="Coupling") myDetailViewer->handleCouplingChanged();
        //! -------------------------------------------------------------------------------------
        //! "DOFs selection" from "Manual"/"Program controlled" to "Program controlled"/"Manual"
        //! Remove/add switches
        //! -------------------------------------------------------------------------------------
        if(propertyName=="DOFs selection") myDetailViewer->handleDOFselectorChange();
    }
        break;

    case SimulationNodeClass::nodeType_connection:
    {
        cout<<"____handling contact____"<<endl;
        cout<<"____"<<myTreeView->selectionModel()->selectedIndexes().length()<<"____"<<endl;

        for(int n=0; n<myTreeView->selectionModel()->selectedIndexes().length(); n++)
        {
            curNode = myTreeView->selectionModel()->selectedIndexes().at(n).data(Qt::UserRole).value<SimulationNodeClass*>();
            //! ----------------------------------------------------------------------
            //! handle the changes in the contact item (add/remove) switches/controls
            //! ----------------------------------------------------------------------
            if(propertyName=="Type")
            {
                Property::contactType theContactType;

                //! ---------------------------------
                //! check if a dummy node is defined
                //! ---------------------------------
                DetailViewer *dw = static_cast<DetailViewer*>(tools::getWidgetByName("detailViewer"));
                if(dw->getCurrentMultipleSelectionNode()!=Q_NULLPTR)
                {
                    theContactType = dw->getCurrentMultipleSelectionNode()->getPropertyValue<Property::contactType>("Type");
                    QVariant data;
                    data.setValue(theContactType);
                    curNode->getModel()->blockSignals(true);
                    curNode->replaceProperty("Type",Property("Type",data,Property::PropertyGroup_Definition));
                    curNode->getModel()->blockSignals(false);
                }
                else
                {
                    item = curNode->getPropertyItem("Type");
                    theContactType = item->data(Qt::UserRole).value<Property>().getData().value<Property::contactType>();
                }
                //item = curNode->getPropertyItem("Type");
                //Property::contactType theContactType = item->data(Qt::UserRole).value<Property>().getData().value<Property::contactType>();

                //! -----------------------------------------------------------------------------------------
                //! handle the "Friction coefficient" control: in case of a "Frictional" contact pair add it
                //! if it was previously removed; in case of a "Frictionless" or "Bonded" contact pair,
                //! remove it, if it was previously added
                //! -----------------------------------------------------------------------------------------
                switch(theContactType)
                {
                case Property::contactType_frictional:
                {
                    cout<<"____handling frictional____"<<endl;
                    //! -----------------------------------
                    //! add the friction if it was removed
                    //! -----------------------------------
                    QStandardItem *itemFriction = curNode->getPropertyItem("Friction coefficient");
                    if(itemFriction==NULL)
                    {
                        QVariant data;
                        data.setValue(0.2);
                        Property prop_friction("Friction coefficient",data,Property::PropertyGroup_Definition);
                        curNode->addProperty(prop_friction,2);
                    }
                    //! ------------------------------------------
                    //! add the "Small sliding" if it was removed
                    //! ------------------------------------------
                    QStandardItem *itemSmallSliding = curNode->getPropertyItem("Small sliding");
                    if(itemSmallSliding==NULL)
                    {
                        QVariant data;

                        //! ----------------------------------------------------------------------------------
                        //! Warning: lock/unlock the "Small sliding" control, according to the "Behavior"
                        //! "Small sliding" = "On" only if the contact pair is node to face (i.e. Asymmetric)
                        //! ----------------------------------------------------------------------------------
                        switch(curNode->getPropertyItem("Behavior")->data(Qt::UserRole).value<Property>().getData().value<Property::contactBehavior>())
                        {
                        case Property::contactBehavior_asymmetric: data.setValue(int(1)); break;
                        case Property::contactBehavior_symmetric: data.setValue(int(0)); break;
                        }
                        Property prop_smallSliding("Small sliding",data,Property::PropertyGroup_Definition);
                        curNode->addProperty(prop_smallSliding,3);
                    }
                    //! -------------------------------------------------------------------------------------------------
                    //! if the "Small sliding" is already present, and the contact pair is face to face (i.e. symmetric)
                    //! force "Small sliding" = "Off" - as before
                    //! -------------------------------------------------------------------------------------------------
                    else
                    {
                        Property::contactBehavior behavior = curNode->getPropertyItem("Behavior")->data(Qt::UserRole).value<Property>().getData().value<Property::contactBehavior>();
                        if(behavior==Property::contactBehavior_symmetric)
                        {
                            QVariant data;
                            data.setValue(int(0));
                            Property prop_smallSliding("Small sliding",data,Property::PropertyGroup_Definition);
                            curNode->replaceProperty("Small sliding",prop_smallSliding);
                        }
                    }
                    //! -------------------------------
                    //! add "Lambda" if it was removed
                    //! -------------------------------
                    QStandardItem *itemLambda = curNode->getPropertyItem("Lambda");
                    if(itemLambda==NULL)
                    {
                        QVariant data;
                        data.setValue(0);
                        Property prop_Lambda("Lambda",data,Property::PropertyGroup_Advanced);
                        curNode->addProperty(prop_Lambda,4);
                    }

                    //! ------------------------------------------------------------------
                    //! re-init the contact with "overpressure linear"
                    //! add "Overpressure", "K", "Sigma infty", "CO" if they were removed
                    //! ------------------------------------------------------------------
                    QVariant data;
                    if(curNode->getPropertyItem("K")==NULL)
                    {
                        data.setValue(0.0);
                        Property prop_overpressure("Overpressure",data,Property::PropertyGroup_Advanced);
                        curNode->addProperty(prop_overpressure,0);
                    }
                    if(curNode->getPropertyItem("K")==NULL)
                    {
                        data.setValue(0.0);
                        Property prop_K("K",data,Property::PropertyGroup_Advanced);
                        curNode->addProperty(prop_K,1);
                    }
                    if(curNode->getPropertyItem("Sigma infty")==NULL)
                    {
                        data.setValue(0.0);
                        Property prop_sigmaInfty("Sigma infty",data,Property::PropertyGroup_Advanced);
                        curNode->addProperty(prop_sigmaInfty,2);
                    }
                    if(curNode->getPropertyItem("C0")==NULL)
                    {
                        data.setValue(0.0);
                        Property prop_C0("C0",data,Property::PropertyGroup_Advanced);
                        curNode->addProperty(prop_C0,3);
                    }
                    //! -------------------------------
                    //! add "Lambda" if it was removed
                    //! -------------------------------
                    if(curNode->getPropertyItem("Lambda")==NULL)
                    {
                        QVariant data;
                        data.setValue(0);
                        Property prop_lambda("Lambda",data,Property::PropertyGroup_Advanced);
                        curNode->addProperty(prop_lambda,4);
                    }
                }
                    break;

                case Property::contactType_frictionless:
                {
                    cout<<"____handling frictionless____"<<endl;
                    //! ------------------------------------
                    //! remove the friction if it was added
                    //! ------------------------------------
                    QStandardItem *itemFriction = curNode->getPropertyItem("Friction coefficient");
                    if(itemFriction!=NULL) curNode->removeProperty("Friction coefficient");
                    //! ------------------------------------------
                    //! add the "Small sliding" if it was removed
                    //! ------------------------------------------
                    QStandardItem *itemSmallSliding = curNode->getPropertyItem("Small sliding");
                    if(itemSmallSliding==NULL)
                    {
                        QVariant data;

                        //! ----------------------------------------------------------------------------------
                        //! Warning: force the "Small sliding" = "0" in case of surface to surface contact
                        //! "Small sliding" = "On" only if the contact pair is node to face (i.e. Asymmetric)
                        //! ----------------------------------------------------------------------------------
                        switch(curNode->getPropertyItem("Behavior")->data(Qt::UserRole).value<Property>().getData().value<Property::contactBehavior>())
                        {
                        case Property::contactBehavior_asymmetric: data.setValue(int(1)); break;
                        case Property::contactBehavior_symmetric: data.setValue(int(0)); break;
                        }
                        Property prop_smallSliding("Small sliding",data,Property::PropertyGroup_Definition);
                        curNode->addProperty(prop_smallSliding,3);
                    }
                    //! -------------------------------------------------------------------------------------------------
                    //! if the "Small sliding" is already present, and the contact pair is face to face (i.e. symmetric)
                    //! force "Small sliding" = "Off" - as before
                    //! -------------------------------------------------------------------------------------------------
                    else
                    {
                        QVariant data;
                        data.setValue(int(0));
                        Property prop_smallSliding("Small sliding",data,Property::PropertyGroup_Definition);
                        curNode->replaceProperty("Small sliding",prop_smallSliding);
                    }

                    //! ------------------------------------------------------------------
                    //! re-init the contact with "overpressure linear"
                    //! add "Overpressure", "K", "Sigma infty", "CO" if they were removed
                    //! ------------------------------------------------------------------
                    QVariant data;
                    if(curNode->getPropertyItem("K")==NULL)
                    {
                        data.setValue(0.0);
                        Property prop_overpressure("Overpressure",data,Property::PropertyGroup_Advanced);
                        curNode->addProperty(prop_overpressure,0);
                    }
                    if(curNode->getPropertyItem("K")==NULL)
                    {
                        data.setValue(0.0);
                        Property prop_K("K",data,Property::PropertyGroup_Advanced);
                        curNode->addProperty(prop_K,1);
                    }
                    if(curNode->getPropertyItem("Sigma infty")==NULL)
                    {
                        data.setValue(0.0);
                        Property prop_sigmaInfty("Sigma infty",data,Property::PropertyGroup_Advanced);
                        curNode->addProperty(prop_sigmaInfty,2);
                    }
                    if(curNode->getPropertyItem("C0")==NULL)
                    {
                        data.setValue(0.0);
                        Property prop_C0("C0",data,Property::PropertyGroup_Advanced);
                        curNode->addProperty(prop_C0,3);
                    }

                    //! ---------------------------
                    //! remove "Lambda" if present
                    //! ---------------------------
                    if(curNode->getPropertyItem("Lambda")!=NULL) curNode->removeProperty("Lambda");
                }
                    break;

                case Property::contactType_bonded:
                {
                    cout<<"____handling bonded____"<<endl;
                    //! ------------------------------------
                    //! remove the friction if it was added
                    //! ------------------------------------
                    curNode->removeProperty("Friction coefficient");

                    //! -------------------------------------------
                    //! remove the "Small sliding" if it was added
                    //! -------------------------------------------
                    curNode->removeProperty("Small sliding");

                    //! ---------------------------------------------------
                    //! remove all the "Advanced" controls, excluding "C0"
                    //! ---------------------------------------------------
                    curNode->removeProperty("Overpressure");
                    curNode->removeProperty("K");
                    curNode->removeProperty("Sigma infty");
                    curNode->removeProperty("Lambda");
                    curNode->removeProperty("P0");
                    curNode->removeProperty("Beta");
                    if(curNode->getPropertyItem("C0")==NULL)
                    {
                        QVariant data;
                        data.setValue(0.0);
                        Property prop_C0("C0",data,Property::PropertyGroup_Advanced);
                        curNode->addProperty(prop_C0);
                    }
                }
                    break;

                case Property::contactType_tied:
                {
                    cout<<"____handing tied____"<<endl;
                    //! ------------------
                    //! removing friction
                    //! ------------------
                    curNode->removeProperty("Friction coefficient");

                    //! -----------------------------------
                    //! force the "Behavior" = "Symmetric"
                    //! -----------------------------------
                    Property::contactBehavior behavior = Property::contactBehavior_symmetric;
                    QVariant data;
                    data.setValue(behavior);
                    Property prop_behavior("Behavior",data,Property::PropertyGroup::PropertyGroup_Definition);
                    curNode->replaceProperty("Behavior",prop_behavior);

                    //! ------------------------------------------
                    //! add the "Small sliding" if it was removed
                    //! force "Small sliding" = "0" (Off)
                    //! ------------------------------------------
                    QStandardItem *itemSmallSliding = curNode->getPropertyItem("Small sliding");
                    if(itemSmallSliding==NULL)
                    {
                        QVariant data;
                        data.setValue(0);
                        Property prop_smallSliding("Small sliding",data,Property::PropertyGroup_Definition);
                        curNode->addProperty(prop_smallSliding,3);
                    }
                    //! ------------------------------------------------------------
                    //! if it was already present force "Small sliding" = "0" (Off)
                    //! ------------------------------------------------------------
                    else
                    {
                        QVariant data;
                        data.setValue(0);
                        Property prop_smallSliding("Small sliding",data,Property::PropertyGroup_Definition);
                        curNode->replaceProperty("Small sliding",prop_smallSliding);
                    }
                    //! ------------------------------------
                    //! force "Overpressure" tied
                    //! ------------------------------------
                    Property::overpressureFunction overpressure = Property::overpressureFunction_tied;
                    data.setValue(overpressure);
                    Property prop_overpressure("Overpressure",data,Property::PropertyGroup_Advanced);
                    curNode->replaceProperty("Overpressure",prop_overpressure);

                    //! ------------------------------------
                    //! remove all the "Advanced"; keep "K"
                    //! ------------------------------------
                    curNode->removeProperty("Sigma infty");
                    curNode->removeProperty("C0");
                    curNode->removeProperty("P0");
                    curNode->removeProperty("Lambda");
                }
                    break;
                }
            }

            //! ------------------------------
            //! handle the "Behavior" control
            //! ------------------------------
            if(propertyName =="Behavior")
            {
                Property::contactBehavior theBehavior = item->data(Qt::UserRole).value<Property>().getData().value<Property::contactBehavior>();
                switch(theBehavior)
                {
                case Property::contactBehavior_symmetric:
                {
                    //! --------------------------------------------------------------------------------------
                    //! If the contact pair is "Frictional" or "Frictionless" it has the "Behavior" property.
                    //! -> Force "Small sliding" = "Off" in case of face to face (i.e) symmetric contact pair
                    //! -> Remove "Sigma infty"
                    //! -> Force "C0" = "0.0"
                    //! --------------------------------------------------------------------------------------
                    QStandardItem *item = curNode->getPropertyItem("Type");
                    Property::contactType type = item->data(Qt::UserRole).value<Property>().getData().value<Property::contactType>();
                    if(type == Property::contactType_frictional || type == Property::contactType_frictionless)
                    {
                        QVariant data;

                        //! ------------------------------
                        //! force "Small sliding" = "Off"
                        //! ------------------------------
                        data.setValue(int(0));
                        Property prop_smallSliding("Small sliding",data,Property::PropertyGroup_Definition);
                        curNode->replaceProperty("Small sliding",prop_smallSliding);

                        //! ---------------------
                        //! remove "Sigma infty"
                        //! ---------------------
                        curNode->removeProperty("Sigma infty");

                        //! -------------------
                        //! force "C0" = "0.0"
                        //! -------------------
                        data.setValue(double(0.0));
                        Property prop_C0("C0",data,Property::PropertyGroup_Advanced);
                        curNode->replaceProperty("C0",prop_C0);
                    }

                    //! case bonded?
                }
                    break;

                case Property::contactBehavior_asymmetric:
                {
                    QStandardItem *item = curNode->getPropertyItem("Type");
                    Property::contactType type = item->data(Qt::UserRole).value<Property>().getData().value<Property::contactType>();
                    if(type == Property::contactType_frictional || type == Property::contactType_frictionless)
                    {
                        //! -----------------------------------------------
                        //! Add "Sigma infty" if it was previously removed
                        //! -----------------------------------------------
                        QStandardItem *itemSigmaInfty = curNode->getPropertyItem("Sigma infty");
                        if(itemSigmaInfty==NULL)
                        {
                            QVariant data;
                            data.setValue(0);
                            Property prop_sigmaInfty("Sigma infty",data,Property::PropertyGroup_Advanced);
                            curNode->addProperty(prop_sigmaInfty,2);
                        }
                        //! ---------------
                        //! "C0" = "0.001"
                        //! ---------------
                        QVariant data;
                        data.setValue(double(0.0));
                        Property prop_C0("C0",data,Property::PropertyGroup_Advanced);
                        curNode->replaceProperty("C0",prop_C0);
                    }
                }
                    break;

                    //! case bonded?
                }
            }

            //! ----------------------------------
            //! handle the "Overpressure" control
            //! ----------------------------------
            if(propertyName == "Overpressure")
            {
                QVariant data;
                Property::overpressureFunction theOverpressure = item->data(Qt::UserRole).value<Property>().getData().value<Property::overpressureFunction>();
                switch(theOverpressure)
                {
                case Property::overpressureFunction_linear:
                {
                    cerr<<"____handling overpressure linear____"<<endl;
                    //! -----------------------------------------------------------------------
                    //! remove all the properties defining an exponential overclosure function
                    //! -----------------------------------------------------------------------
                    item = curNode->getPropertyItem("P0");
                    if(item!=NULL) curNode->removeProperty("P0");
                    item = curNode->getPropertyItem("Beta");
                    if(item!=NULL) curNode->removeProperty("Beta");

                    //! ------------------------------------------------------------------
                    //! add the properties defining a (quasi) linear overclosure function
                    //! ------------------------------------------------------------------
                    item = curNode->getPropertyItem("K");
                    if(item==NULL)
                    {
                        data.setValue(0);
                        Property prop_K("K",data,Property::PropertyGroup_Advanced);
                        curNode->addProperty(prop_K,1);
                    }
                    item = curNode->getPropertyItem("Sigma infty");
                    if(item==NULL)
                    {
                        data.setValue(0);
                        Property prop_sigmaInfty("Sigma infty",data,Property::PropertyGroup_Advanced);
                        curNode->addProperty(prop_sigmaInfty,2);
                    }
                    item = curNode->getPropertyItem("C0");
                    if(item==NULL)
                    {
                        data.setValue(0);
                        Property prop_C0("C0",data,Property::PropertyGroup_Advanced);
                        curNode->addProperty(prop_C0,3);
                    }
                }
                    break;

                case Property::overpressureFunction_exponential:
                {
                    cerr<<"____handling overpressure exponential____"<<endl;
                    //! -----------------------------------------------------------------
                    //! remove all the properties defining a linear overclosure function
                    //! -----------------------------------------------------------------
                    curNode->removeProperty("K");
                    curNode->removeProperty("Sigma infty");

                    //! ----------------------------------------------------------------
                    //! add the properties defining an exponential overclosure function
                    //! ----------------------------------------------------------------
                    item = curNode->getPropertyItem("P0");
                    if(item==NULL)
                    {
                        data.setValue(0);
                        Property prop_P0("P0",data,Property::PropertyGroup_Advanced);
                        curNode->addProperty(prop_P0,1);
                    }
                    item = curNode->getPropertyItem("C0");
                    if(item==NULL)
                    {
                        data.setValue(0);
                        Property prop_C0("C0",data,Property::PropertyGroup_Advanced);
                        curNode->addProperty(prop_C0,2);
                    }
                }
                    break;
                }

                Property::contactType type = curNode->getPropertyItem("Type")->data(Qt::UserRole).value<Property>().getData().value<Property::contactType>();
                switch(type)
                {
                case Property::contactType_frictional:
                    if(curNode->getPropertyItem("Lambda")==NULL)
                    {
                        QVariant data;
                        data.setValue(0);
                        Property prop_lambda("Lambda",data,Property::PropertyGroup_Advanced);
                        curNode->addProperty(prop_lambda,4);
                    }
                    break;

                case Property::contactType_frictionless:
                    if(curNode->getPropertyItem("Lambda")!=NULL) curNode->removeProperty("Lambda");
                    break;
                }
            }
            cout<<"____handling contacts: exiting____"<<endl;
        }
    }
        break;
    }

    //! -----------------
    //! reconnection [*]
    //! -----------------
    connect(curNode->getModel(),SIGNAL(itemChanged(QStandardItem*)),this,SLOT(handleItemChange(QStandardItem*)));
    //curNode->getModel()->blockSignals(false);

    //! ------
    //! parse
    //! ------
    //parser::parseItem(curItem);

    cout<<"*-------------------------------------------------*"<<endl;
    cout<<"*-  handleItemChange()->____exiting function____ -*"<<endl;
    cout<<"*-------------------------------------------------*"<<endl;
}

//! --------------------------------------------------
//! function: prepareForMeshing
//! details:  1. remove the obsolete meshes
//!           2. build the list of the shapes to mesh
//! --------------------------------------------------
TopTools_ListOfShape SimulationManager::prepareForMeshing()
{
    cout<<"SimulationManager::prepareForMeshing()->____function called____"<<endl;

    //! -----------------------------------------------------
    //! request the OCC viewer to remove the obsolete meshes
    //! -----------------------------------------------------
    emit requestRemoveObsoleteMeshes();

    //! list of the TopoDS_Shape to be meshed
    TopTools_ListOfShape shapeToBeMeshed;

    //! --------------------------------------
    //! 1. if the user selection is not empty
    //! --------------------------------------
    myCTX->InitSelected();
    if(myCTX->MoreSelected())
    {
        for(myCTX->InitSelected();myCTX->MoreSelected();myCTX->NextSelected())
        {
            const occHandle(AIS_Shape) &theAIS = occHandle(AIS_Shape)::DownCast(myCTX->SelectedInteractive());
            const TopoDS_Shape &theShape = theAIS->Shape();
            int bodyIndex = mySimulationDataBase->bodyMap.key(theShape);

            //! ---------------------------------------------------------------------------
            //! check if the shape is not suppressed && the shape has not a valid mesh yet
            //! ---------------------------------------------------------------------------
            if(mySimulationDataBase->ArrayOfMeshDS.value(bodyIndex).IsNull() &&
                    mySimulationDataBase->MapOfIsActive.value(bodyIndex)==true)
                shapeToBeMeshed.Append(theShape);
        }
        myCTX->ClearSelected(true);
        cout<<"SimulationManager::prepareForMeshing()->____number of shapes to mesh: "<<shapeToBeMeshed.Extent()<<"____"<<endl;
        return shapeToBeMeshed;
    }
    //! -------------------------------------------------
    //! 2. the selection is empty => mesh all the bodies
    //! -------------------------------------------------
    for(QMap<int,TopoDS_Shape>::iterator it = mySimulationDataBase->bodyMap.begin(); it!=mySimulationDataBase->bodyMap.end(); ++it)
    {
        int bodyIndex = it.key();        
        const TopoDS_Shape &theShape = it.value();

        //! ---------------------------------------------------------------------
        //! check if the shape is not suppressed and it has not a valid mesh yet
        //! ---------------------------------------------------------------------
        if(mySimulationDataBase->ArrayOfMeshDS.value(bodyIndex)==occHandle(MeshVS_DataSource)() &&
                mySimulationDataBase->MapOfIsActive.value(bodyIndex)==true)
        {
            shapeToBeMeshed.Append(theShape);
        }
    }
    return shapeToBeMeshed;
}

//! ----------------------------
//! function: buildMesh
//! details:  generate the mesh
//! ----------------------------
void SimulationManager::buildMesh(bool isVolumeMesh)
{
    cout<<"SimulationManager::buildMesh()->____function called____"<<endl;

    //! ---------------------------
    //! change the status variable
    //! ---------------------------
    myIsMeshingRunning = true;

    //! -------------------------------------------------------------
    //! transfer mesh parameters to the "array" part of the database
    //! -------------------------------------------------------------
    this->transferMeshNodes();

    //! ------------------------------------------------------
    //! prepare for meshing: retrieve the list shapes to mesh
    //! ------------------------------------------------------
    const TopTools_ListOfShape &shapeToBeMeshed = this->prepareForMeshing();
    cout<<"SimulationManager::buildMesh()->____number of shape to mesh: "<<shapeToBeMeshed.Extent()<<"____"<<endl;

    //! -----------------------
    //! list of bodies to mesh
    //! -----------------------
    QList<int> theListOfBodies;
    TopTools_ListIteratorOfListOfShape anIter;
    for(anIter.Initialize(shapeToBeMeshed);anIter.More();anIter.Next())
    {
        const TopoDS_Shape &theShape = anIter.Value();
        int bodyIndex = mySimulationDataBase->bodyMap.key(theShape);
        theListOfBodies.append(bodyIndex);
    }

    //! -----------------------
    //! the progress indicator
    //! -----------------------
    QWidget *progressIndicatorWidget = tools::getWidgetByName("progressIndicator");
    QProgressIndicator *progressIndicator = static_cast<QProgressIndicator*>(progressIndicatorWidget);
    if(progressIndicator == Q_NULLPTR) cout<<"SimulationManager::buildMesh()->_____warning: the progress indicator is NULL____"<<endl;

    //! ----------------------------------------------------------------------------
    //! scan the contactd and the boundary conditions of all the analysis branches
    //! In case of patch independent meshing method, with geometry defeaturing and
    //! simplification, the boundary of the patches of the boundary conditions will
    //! be preserved (if the "Preserve boundary condition edges" selector is ON
    //! ----------------------------------------------------------------------------
    QVector<GeometryTag> vecTags;
    int type = 4;       //! on faces
    mainTreeTools::getAllBoundaryConditionsTags(myTreeView,type,vecTags);

    //! --------------------------------------
    //! immediately transfer to mesh database
    //! --------------------------------------
    mySimulationDataBase->featuredTags = vecTags;

    //! --------------------------------
    //! meshing process on a new thread
    //! --------------------------------
    myMeshingServer->init(mySimulationDataBase,theListOfBodies,isVolumeMesh,progressIndicator);
    disconnect(myMeshingServer,SIGNAL(meshingFinished(bool)),this,SLOT(handleMeshingResults(bool)));
    connect(myMeshingServer,SIGNAL(meshingFinished(bool)),this,SLOT(handleMeshingResults(bool)));

    //! --------------
    //! start meshing
    //! --------------
    emit myMeshingServer->operate();

    //! ----------------------------------------------
    //! single thread (this solution is working fine)
    //! ----------------------------------------------
    //MesherClass *theMesher = new MesherClass(mySimulationDataBase,theListOfBodies,isVolumeMesh,this);
    //theMesher->generateMesh();
    //this->handleMeshingResults();
}

//! -------------------------------
//! function: handleMeshingResults
//! details:
//! -------------------------------
void SimulationManager::handleMeshingResults(bool isMeshingSuccessfull)
{
    if(isMeshingSuccessfull)
    {
        emit requestBuildMeshIOs();

        //! ------------------------------------------------
        //! build mesh IO - experimental
        //! put the mesh interactive objects into the items
        //! ------------------------------------------------
        //this->buildMeshIO();

        emit requestSetWorkingMode(0);
        SimulationNodeClass *nodeMeshRoot = this->getTreeItem(SimulationNodeClass::nodeType_meshControl)->data(Qt::UserRole).value<SimulationNodeClass*>();
        bool showMeshNodes = nodeMeshRoot->getPropertyValue<bool>("Show mesh nodes");

        //! ----------------
        //! show the meshes
        //! ----------------
        emit requestShowMeshes(showMeshNodes);
    }
    //! -----------------------
    //! update mesh statistics
    //! -----------------------
    this->updateMeshStatistics();

    //! ---------------
    //! update display
    //! ---------------
    this->changeColor();

    //! ---------------------------
    //! change the status variable
    //! ---------------------------
    myIsMeshingRunning = false;
}

//! -------------------------------
//! function: updateMeshStatistics
//! details:
//! -------------------------------
void SimulationManager::updateMeshStatistics()
{
    //! ----------------------------------------
    //! update the number of nodes and elements
    //! ----------------------------------------
    int NN3D=0, NE3D=0, NN2D=0, NE2D=0;
    for(int i=1; i<=mySimulationDataBase->ArrayOfMeshDS.size(); i++)
    {
        int bodyIndex = mySimulationDataBase->bodyMap.key(mySimulationDataBase->bodyMap.value(i));
        if(!mySimulationDataBase->ArrayOfMeshDS.value(bodyIndex).IsNull())
        {
            NN3D = NN3D + mySimulationDataBase->ArrayOfMeshDS.value(bodyIndex)->GetAllNodes().Extent();
            NE3D = NE3D + mySimulationDataBase->ArrayOfMeshDS.value(bodyIndex)->GetAllElements().Extent();
        }
    }
    for(int i=1; i<=mySimulationDataBase->ArrayOfMeshDS2D.size(); i++)
    {
        int bodyIndex = mySimulationDataBase->bodyMap.key(mySimulationDataBase->bodyMap.value(i));
        if(!mySimulationDataBase->ArrayOfMeshDS2D.value(bodyIndex).IsNull())
        {
            NN2D = NN2D + mySimulationDataBase->ArrayOfMeshDS2D.value(bodyIndex)->GetAllNodes().Extent();
            NE2D = NE2D + mySimulationDataBase->ArrayOfMeshDS2D.value(bodyIndex)->GetAllElements().Extent();
        }
    }

    int NN, NE;
    if(NN3D != 0) { NN = NN3D; NE = NE3D; }
    else { NN = NN2D; NE = NE2D; }

    //! ---------------------------
    //! update the mesh statistics
    //! ---------------------------
    Property prop_NN("Number of nodes",NN,Property::PropertyGroup_Statistics);
    Property prop_NE("Number of elements",NE,Property::PropertyGroup_Statistics);
    SimulationNodeClass *meshRootNode = Mesh_RootItem->data(Qt::UserRole).value<SimulationNodeClass*>();

    //! -----------------------
    //! block the node signals
    //! -----------------------
    meshRootNode->getModel()->blockSignals(true);

    meshRootNode->replaceProperty("Number of nodes",prop_NN);
    meshRootNode->replaceProperty("Number of elements",prop_NE);

    //! ------------------------
    //! unlock the node signals
    //! ------------------------
    meshRootNode->getModel()->blockSignals(false);
}

//! --------------------------------------------------------
//! function: changeNodeSuppressionStatus
//! details:  also handle visibility. Works on items and on
//!           a list of shapes selected in the OCC viewer
//! --------------------------------------------------------
void SimulationManager::changeNodeSuppressionStatus(Property::SuppressionStatus newSuppressionStatus)
{
    cout<<"SimulationManager::changeNodeSuppressionStatus()->____function called____"<<endl;

    //! ----------------------------
    //! cast the suppression status
    //! ----------------------------
    bool isSuppressed = newSuppressionStatus==Property::SuppressionStatus_Active? false:true;

    QVariant data;
    data.setValue(newSuppressionStatus);
    Property prop_newSuppressionStatus("Suppressed",data,Property::PropertyGroup_Definition);

    //! ---------------------------------------------------------------
    //! body indexes (if the selected items contains bodies):
    //! the list could be empty if nothing is selected from the viewer
    //! or no body has been selected from the main tree
    //! ---------------------------------------------------------------
    TColStd_ListOfInteger ListOfBodyNumbers;

    //! --------------------------------------------------------------------
    //! handle a selection from the OCC viewer: by definition the selection
    //! contains TopoDS_Shape(s) since done within the viewer
    //! --------------------------------------------------------------------
    myCTX->InitSelected();
    if(myCTX->MoreSelected())
    {
        //! --------------------------
        //! retrieve the body indices
        //! --------------------------
        TopTools_ListOfShape aListOfShape;
        for(myCTX->InitSelected(); myCTX->MoreSelected(); myCTX->NextSelected())
        {
            const occHandle(AIS_Shape) &theAISShape = occHandle(AIS_Shape)::DownCast(myCTX->SelectedInteractive());
            const TopoDS_Shape &theShape = theAISShape->Shape();
            aListOfShape.Append(theShape);
            int bodyIndex = mySimulationDataBase->bodyMap.key(theShape);
            ListOfBodyNumbers.Append(bodyIndex);
        }
        //! ------------------------------------------------------------------------
        //! scan the main tree and update the properties "Suppressed" and "Visible"
        //! ------------------------------------------------------------------------
        QExtendedStandardItem *itemGeometry = this->getTreeItem(SimulationNodeClass::nodeType_geometry);
        for(int i=0; i<itemGeometry->rowCount(); i++)
        {
            QExtendedStandardItem *itemBody = static_cast<QExtendedStandardItem*>(itemGeometry->child(i,0));
            SimulationNodeClass *nodeBody = itemBody->data(Qt::UserRole).value<SimulationNodeClass*>();
            int bodyIndex = nodeBody->getPropertyValue<int>("Map index");
            if(ListOfBodyNumbers.Contains(bodyIndex))
            {
                //! -----------------------
                //! block the node signals
                //! -----------------------
                nodeBody->getModel()->blockSignals(true);

                //! ----------------------------------------------------
                //! synch the "Suppressed" property and update the icon
                //! ----------------------------------------------------
                nodeBody->replaceProperty("Suppressed",prop_newSuppressionStatus);
                if(isSuppressed) itemBody->setIcon(QIcon(":/icons/icon_suppress.png"));
                else itemBody->setIcon(itemBody->getIcon(nodeBody->getType()));

                //! -----------------------------------
                //! synch also the visibility property
                //! -----------------------------------
                bool isVisible = newSuppressionStatus==Property::SuppressionStatus_Active? true:false;
                data.setValue(isVisible);
                Property prop_visible("Visible",data,Property::PropertyGroup_GraphicProperties);
                nodeBody->replaceProperty("Visible",prop_visible);

                //! ------------------------
                //! unlock the node signals
                //! ------------------------
                nodeBody->getModel()->blockSignals(false);
            }
        }
    }
    //! --------------------------------------
    //! handle a selection from the main tree
    //! --------------------------------------
    else if(mySelectionModel->hasSelection())
    {
        QList<QStandardItem*> listOfItems;
        for(int i=0; i<mySelectionModel->selectedIndexes().length(); i++)
        {
            QModelIndex index = mySelectionModel->selectedIndexes().at(i);
            QStandardItem *item = myModel->itemFromIndex(index);
            listOfItems.append(item);
        }
        //! ------------------------
        //! scan the selected items
        //! ------------------------
        for(int k=0; k<listOfItems.length(); k++)
        {
            QExtendedStandardItem *theCurItem = static_cast<QExtendedStandardItem*>(listOfItems.at(k));
            SimulationNodeClass *theCurNode = theCurItem->data(Qt::UserRole).value<SimulationNodeClass*>();
            SimulationNodeClass::nodeType theNodeType = theCurNode->getType();

            //! -----------------------
            //! block the node signals
            //! -----------------------
            theCurNode->getModel()->blockSignals(true);

            //! ---------------------------------------------------------------------------
            //! if the item is a geometry item, store the bodyIndex into ListOfBodyNumbers
            //! ---------------------------------------------------------------------------
            if(theNodeType == SimulationNodeClass::nodeType_geometryBody)
            {
                int bodyIndex = theCurNode->getPropertyValue<int>("Map index");
                ListOfBodyNumbers.Append(bodyIndex);

                //! -----------------------------
                //! synch the "Visible" property
                //! -----------------------------
                bool isVisible = newSuppressionStatus==Property::SuppressionStatus_Active? true:false;
                data.setValue(isVisible);
                Property prop_visible("Visible",data,Property::PropertyGroup_GraphicProperties);
                theCurNode->replaceProperty("Visible",prop_visible);
            }

            //! ----------------------------------------------------
            //! change the suppression property and update the icon
            //! ----------------------------------------------------
            theCurNode->replaceProperty("Suppressed", prop_newSuppressionStatus);
            if(isSuppressed) theCurItem->setIcon(QIcon(":/icons/icon_suppress.png"));
            else theCurItem->setIcon(theCurItem->getIcon(theCurNode->getType()));

            //! ------------------------
            //! unlock the node signals
            //! ------------------------
            theCurNode->getModel()->blockSignals(false);
        }
    }

    //! ------------------------------------------
    //! set the suppression flag into the "array"
    //! part of the simulation database
    //! ------------------------------------------
    for(TColStd_ListIteratorOfListOfInteger anIt(ListOfBodyNumbers); anIt.More(); anIt.Next())
    {
        int bodyIndex = anIt.Value();
        mySimulationDataBase->MapOfIsActive.insert(bodyIndex, isSuppressed==true? false: true);
    }

    //! --------------------------------------------------------
    //! hide/display bodies according to the suppression status
    //! --------------------------------------------------------
    if(isSuppressed) emit requestHideBody(ListOfBodyNumbers);
    else emit requestShowBody(ListOfBodyNumbers);

    //! -------------------------------------------
    //! delete the meshes of the suppressed bodies
    //! -------------------------------------------
    std::vector<int> meshOfBodies;
    for(TColStd_ListIteratorOfListOfInteger anIt(ListOfBodyNumbers); anIt.More(); anIt.Next())
    {
        int bodyIndex = anIt.Value();
        meshOfBodies.push_back(bodyIndex);
    }
    this->requestMeshInvalidate(meshOfBodies);
}

//! -------------------------------------------------------
//! function: unsuppressAllBodies
//! details:  available when at least a body is suppressed
//! -------------------------------------------------------
void SimulationManager::unsuppressAllBodies()
{
    TColStd_ListOfInteger ListOfBodyNumbers;
    for(int k=0;k<Geometry_RootItem->rowCount();k++)
    {
        QExtendedStandardItem *item = static_cast<QExtendedStandardItem*>(Geometry_RootItem->child(k,0));
        SimulationNodeClass *theNode = item->data(Qt::UserRole).value<SimulationNodeClass*>();
        Property::SuppressionStatus ss = theNode->getPropertyValue<Property::SuppressionStatus>("Suppressed");

        //! -----------------------
        //! block the node signals
        //! -----------------------
        theNode->getModel()->blockSignals(true);

        //! -------------------------
        //! unsuppress if suppressed
        //! -------------------------
        if(ss==Property::SuppressionStatus_Suppressed)
        {
            QVariant data;
            data.setValue(Property::SuppressionStatus_Active);
            Property prop_suppressed("Suppressed",data,Property::PropertyGroup_Definition);
            theNode->replaceProperty("Suppressed",prop_suppressed);

            //! ------------------------------------------
            //! change the icon and handle the visibility
            //! ------------------------------------------
            item->setIcon(item->getIcon(theNode->getType()));
            if(theNode->getPropertyItem("Visible")!=Q_NULLPTR)
            {
                data.setValue(true);
                Property prop_visible("Visible",data,Property::PropertyGroup_GraphicProperties);
                theNode->replaceProperty("Visible",prop_visible);
            }

            int bodyIndex = theNode->getPropertyValue<int>("Map index");
            ListOfBodyNumbers.Append(bodyIndex);

            //! --------------------------------------------
            //! change the content of the geometry database
            //! --------------------------------------------
            mySimulationDataBase->MapOfIsActive.insert(bodyIndex,true);
        }

        //! ------------------------
        //! unlock the node signals
        //! ------------------------
        theNode->getModel()->blockSignals(false);
    }    
    if(ListOfBodyNumbers.Extent()!=0) emit requestShowBody(ListOfBodyNumbers);
}

//! ----------------------------------------------------------
//! function: handleVisibilityChange
//! details:  the user has switched the combobox to "Visible"
//! ----------------------------------------------------------
void SimulationManager::handleVisibilityChange(bool newIsVisible)
{
    TColStd_ListOfInteger ListOfBodyNumbers;
    for(int i=0; i<mySelectionModel->selectedIndexes().length(); i++)
    {
        QModelIndex index = mySelectionModel->selectedIndexes().at(i);
        SimulationNodeClass *theNode = myModel->itemFromIndex(index)->data(Qt::UserRole).value<SimulationNodeClass*>();
        int bodyIndex = theNode->getPropertyValue<int>("Map index");

        //! ----------------------------
        //! check if the body is active
        //! ----------------------------
        Property::SuppressionStatus ss = theNode->getPropertyValue<Property::SuppressionStatus>("Suppressed");
        if(ss==Property::SuppressionStatus_Suppressed) continue;
        ListOfBodyNumbers.Append(bodyIndex);
    }
    if(newIsVisible) emit requestShowBody(ListOfBodyNumbers);
    else emit requestHideBody(ListOfBodyNumbers);
}

//! --------------------------------------------------------------------------
//! function: itemsFromSelection
//! details:  for a given list of shapes return the list of the corresponging
//!           geometry items
//! --------------------------------------------------------------------------
QList<QStandardItem*> SimulationManager::ItemListFromListOfShape(TopTools_ListOfShape *listOfShapes)
{
    //!cout<<"SimulationManager::ItemListFromListOfShape()->____number of shapes to convert in items: "<<listOfShapes->Extent()<<"____"<<endl;
    QList<QStandardItem*> theListOfItems;
    int N = Geometry_RootItem->rowCount();
    TopTools_ListIteratorOfListOfShape anIter;
    for(anIter.Initialize(*listOfShapes);anIter.More();anIter.Next())
    {
        const TopoDS_Shape &theCurShape = anIter.Value();
        int bodyIndex = mySimulationDataBase->bodyMap.key(theCurShape);
        for(int k=0; k<N; k++)
        {
            QStandardItem *item = Geometry_RootItem->child(k,0);
            int theCurMapIndex = item->data(Qt::UserRole).value<SimulationNodeClass*>()->getPropertyValue<int>("Map index");
            if(bodyIndex == theCurMapIndex)
            {
                theListOfItems.append(item);
                break;
            }
        }
    }
    return theListOfItems;
}

//! ----------------------------------------------------------
//! function: synchVisibility
//! details:  update the "Visible" property of the node model
//!           when "hide"/"show" are called from the viewer
//! ----------------------------------------------------------
void SimulationManager::synchVisibility()
{
    //!cout<<"SimulationManager::synchVisibility()->____function called____"<<endl;
    AIS_ListOfInteractive listOfDisplayedObjects, listOfHiddenObjects;
    TopTools_ListOfShape listOfDisplayedShapes, listOfHiddenShapes;

    myCTX->DisplayedObjects(AIS_KOI_Shape, -1, listOfDisplayedObjects);
    myCTX->ErasedObjects(AIS_KOI_Shape, -1, listOfHiddenObjects);

    AIS_ListIteratorOfListOfInteractive anIter;
    for(anIter.Initialize(listOfDisplayedObjects);anIter.More();anIter.Next())
    {
        const occHandle(AIS_Shape) &theAISShape = occHandle(AIS_Shape)::DownCast(anIter.Value());
        if(theAISShape.IsNull()) continue;
        const TopoDS_Shape &theShape = theAISShape->Shape();
        listOfDisplayedShapes.Append(theShape);
    }

    //!cout<<"SimulationManager::synchVisibility()->____displayed TopoDS_Shape: "<<listOfDisplayedShapes.Extent()<<"____"<<endl;

    QVariant data;
    data.setValue(true);
    Property property_visible("Visible",data,Property::PropertyGroup_GraphicProperties);
    QList<QStandardItem*> theListOfDisplayedItems = this->ItemListFromListOfShape(&listOfDisplayedShapes);
    for(int i=0;i<theListOfDisplayedItems.length();i++)
    {
        QExtendedStandardItem *item = static_cast<QExtendedStandardItem *>(theListOfDisplayedItems.at(i));
        SimulationNodeClass *theNode = item->data(Qt::UserRole).value<SimulationNodeClass*>();
        theNode->getModel()->blockSignals(false);
        theNode->replaceProperty("Visible",property_visible);
        tools::changeIconOpacity(item,false);
        theNode->getModel()->blockSignals(false);
    }

    for(anIter.Initialize(listOfHiddenObjects);anIter.More();anIter.Next())
    {
        const occHandle(AIS_Shape) &theAISShape = occHandle(AIS_Shape)::DownCast(anIter.Value());
        if(theAISShape.IsNull()) continue;
        const TopoDS_Shape &theShape = theAISShape->Shape();
        listOfHiddenShapes.Append(theShape);
    }

    //!cout<<"SimulationManager::synchVisibility()->____hidden TopoDS_Shape: "<<listOfHiddenShapes.Extent()<<"____"<<endl;

    data.setValue(false);
    Property property_hidden("Visible",data,Property::PropertyGroup_GraphicProperties);
    QList<QStandardItem*> theListOfHiddenItems = this->ItemListFromListOfShape(&listOfHiddenShapes);
    for(int i=0;i<theListOfHiddenItems.length();i++)
    {
        QExtendedStandardItem *item = static_cast<QExtendedStandardItem *>(theListOfHiddenItems.at(i));
        SimulationNodeClass *theNode = item->data(Qt::UserRole).value<SimulationNodeClass*>();
        theNode->getModel()->blockSignals(true);
        theNode->replaceProperty("Visible",property_hidden);
        tools::changeIconOpacity(item,true);
        theNode->getModel()->blockSignals(false);
    }
}

//! ----------------------------------------------
//! function: duplicateItem
//! details:  duplicate the current selected item
//!           or an item passed by argument
//! ----------------------------------------------
void SimulationManager::duplicateItem(QExtendedStandardItem *item)
{
    //! ----------------------
    //! the item to be copied
    //! ----------------------
    QStandardItem *theItemToCopy;
    if(item==Q_NULLPTR)
    {
        QModelIndex index = myTreeView->currentIndex();
        if(!index.isValid()) return;
        theItemToCopy= myModel->itemFromIndex(myTreeView->currentIndex());
        if(theItemToCopy==Q_NULLPTR) return;
    }
    else
    {
        theItemToCopy = item;
        if(item==Q_NULLPTR) return;
    }
    SimulationNodeClass *theNode = theItemToCopy->data(Qt::UserRole).value<SimulationNodeClass*>();

    //! --------------------------------------------------
    //! scan the node properties and copy them one by one
    //! put the copied properties into "props"
    //! --------------------------------------------------
    QVector<QExtendedStandardItem*> itemVec = theNode->getPropertyItems();

    QVectorIterator<QExtendedStandardItem*> anIter(itemVec);
    QVector<Property> props;
    while(anIter.hasNext())
    {
        QExtendedStandardItem* theCurPropertyItem = anIter.next();
        Property theCurProperty = theCurPropertyItem->data(Qt::UserRole).value<Property>();
        Property copyProperty(theCurProperty);  // only for testing the copy constructor
        props.push_back(copyProperty);
    }

    //! -----------------------------------------
    //! retrieve the node type and the node name
    //! -----------------------------------------
    SimulationNodeClass::nodeType theType = theNode->getType();
    QString theNodeName = theNode->getName();

    //! ------------------------------------------------------------------------------------
    //! create the node specific property - warning: come non specific node properties (i.e
    //! independent on the node type) are created at node constructio time, by consequence
    //! when duplicating, one of these properties can result present more than one time.
    //! For this reason in "SimulationNodeClass" the patch [*] has been indroduced
    //! [Ref. 1 - in SimulationNodeClass]
    //! ------------------------------------------------------------------------------------
    SimulationNodeClass *theNewNode = new SimulationNodeClass(theNodeName,theType,props,this);

    QVariant data;
    data.setValue(theNewNode);
    QExtendedStandardItem *itemCopy = new QExtendedStandardItem();
    itemCopy->setData(data,Qt::UserRole);

    //! ----------------------------------------------------
    //! get the parent and build a name for the copied item
    //! ----------------------------------------------------
    QStandardItem *theParent = theItemToCopy->parent();
    int N=1;
    for(int k=0; k<theParent->rowCount()-1; k++)
    {
        QStandardItem *curItem = theParent->child(k,0);
        SimulationNodeClass::nodeType curType = curItem->data(Qt::UserRole).value<SimulationNodeClass*>()->getType();
        if(curType==theNode->getType()) N++;
    }
    QString copyName = theNode->getName()+QString(" %1").arg(N);
    itemCopy->setData(copyName,Qt::DisplayRole);

    //! -----------------------------------------
    //! copy the tabular data when present
    //! 1) retrieve the "Analysis settings" item
    //! 2) retrieve the column contaning data
    //! 3) copy it and append to the table
    //! -----------------------------------------
    if(theNode->hasTabularData())
    {
        QStandardItem *itemAnslysisSettings = theItemToCopy->parent()->child(0,0);
        SimulationNodeClass *nodeAnalysisSettings = itemAnslysisSettings->data(Qt::UserRole).value<SimulationNodeClass*>();
        CustomTableModel *theTabularData = nodeAnalysisSettings->getTabularDataModel();

        //! -----------------------------------------------------------------------
        //! determine the number of columns to copy and append to the tabular data
        //! -----------------------------------------------------------------------
        int NbCol = mainTreeTools::getColumnsToRead(myTreeView).length();
        int columnToCopy = mainTreeTools::calculateStartColumn(myTreeView);

        for(int i=0; i<NbCol; i++)
        {
            load loadToCopy = theTabularData->getColumn(columnToCopy+i);
            theTabularData->appendColumn(loadToCopy);
        }
    }

    //! ---------------------
    //! replace the time tag
    //! ---------------------
    disconnect(theNewNode->getModel(),SIGNAL(itemChanged(QStandardItem*)),this,SLOT(handleItemChange(QStandardItem*)));
    QThread::msleep(5);
    QString timeTag = QString::fromStdString(SimulationNodeClass::generateTimeString());
    data.setValue(timeTag);
    Property prop_timeTag("Time tag",timeTag,Property::PropertyGroup_Identifier);
    theNewNode->replaceProperty("Time tag",prop_timeTag);
    connect(theNewNode->getModel(),SIGNAL(itemChanged(QStandardItem*)),this,SLOT(handleItemChange(QStandardItem*)));

    //! ----------------
    //! attach the item
    //! ----------------
    theParent->insertRow(this->getInsertionRow(),itemCopy);
    myTreeView->setCurrentIndex(itemCopy->index());
}

//! -------------------------------
//! function: changeElementControl
//! details:
//! -------------------------------
void SimulationManager::ChangeElementControl()
{
    //!cout<<"SimulationManager::ChangeElementControl()->____change element control____"<<endl;
    Property::elementControl theElementControl = Geometry_RootItem->data(Qt::UserRole).value<SimulationNodeClass*>()
            ->getPropertyValue<Property::elementControl>("Element control");

    int N = this->getAllTreeItemOfType(SimulationNodeClass::nodeType_geometryBody).length();

    //! if the element control is "Program controlled" remove all the items "Integration scheme"
    if(theElementControl==Property::elementControl_programControlled)
    {
        if(Geometry_RootItem->hasChildren())
        {
            for(int k=0; k<N; k++)
            {
                QExtendedStandardItem *bodyItem = static_cast<QExtendedStandardItem*>(Geometry_RootItem->child(k,0));
                SimulationNodeClass* theBodyNode = bodyItem->data(Qt::UserRole).value<SimulationNodeClass*>();
                theBodyNode->getModel()->blockSignals(true);
                theBodyNode->removeProperty("Integration scheme");
                theBodyNode->getModel()->blockSignals(false);
            }
        }
    }

    //! if the element control is "Manual" add "Integration scheme"
    else if(theElementControl==Property::elementControl_manual)
    {
        if(Geometry_RootItem->hasChildren())
        {
            for(int k=0; k<N; k++)
            {
                QExtendedStandardItem *bodyItem = static_cast<QExtendedStandardItem*>(Geometry_RootItem->child(k,0));
                SimulationNodeClass* theBodyNode = bodyItem->data(Qt::UserRole).value<SimulationNodeClass*>();
                theBodyNode->getModel()->blockSignals(true);
                QVariant data;
                data.setValue(Property::integrationScheme_full);
                Property property_integrationScheme("Integration scheme",data,Property::PropertyGroup_Definition);
                theBodyNode->addProperty(property_integrationScheme);
                theBodyNode->getModel()->blockSignals(false);
            }
        }
    }
}

//! ----------------------
//! function: changeColor
//! details:
//! ----------------------
void SimulationManager::changeColor()
{
    QModelIndex index = myTreeView->currentIndex();
    SimulationNodeClass *node = index.data(Qt::UserRole).value<SimulationNodeClass*>();
    SimulationNodeClass::nodeType nodeType = node->getType();
    if(node->isSolution() ||
            node->isSolutionInformation() ||
            node->isAnalysisResult() ||
            node->isAnalysisRoot() ||
            node->isAnalysisSettings() ||
            nodeType == SimulationNodeClass::nodeType_root ||
            nodeType == SimulationNodeClass::nodeType_geometry ||
            nodeType == SimulationNodeClass::nodeType_geometryBody ||
            nodeType == SimulationNodeClass::nodeType_geometryPart ||
            nodeType == SimulationNodeClass::nodeType_meshControl ||
            nodeType == SimulationNodeClass::nodeType_coordinateSystem ||
            nodeType == SimulationNodeClass::nodeType_coordinateSystem_global ||
            nodeType == SimulationNodeClass::nodeType_coordinateSystems ||
            nodeType == SimulationNodeClass::nodeType_remotePointRoot ||
            nodeType == SimulationNodeClass::nodeType_remotePoint ||
            nodeType == SimulationNodeClass::nodeType_namedSelection ||
            nodeType == SimulationNodeClass::nodeType_connection ||
            nodeType == SimulationNodeClass::nodeType_connectionGroup)
    {
        //! ---------------------------
        //! reset all the shape colors
        //! ---------------------------
        emit requestResetCustomColors(true);
        return;
    }
    if(nodeType == SimulationNodeClass::nodeType_connectionPair)
    {
        emit requestResetCustomColors(true);

        const QVector<GeometryTag> &tagsMaster = node->getPropertyValue<QVector<GeometryTag>>("Tags master");
        if(tagsMaster.isEmpty()==false)
        {
            QMap<GeometryTag,TopoDS_Shape> subShapesMap;
            for(int i=0; i<tagsMaster.size(); i++)
            {
                const GeometryTag &curTag = tagsMaster[i];
                if(subShapesMap.contains(curTag)) continue;
                const TopoDS_Shape &aSubShape = this->fromTagToShape(curTag);
                subShapesMap.insert(curTag,aSubShape);
            }
            emit requestApplyCustomColor(subShapesMap,Quantity_NOC_BLUE1,false);
        }

        const QVector<GeometryTag> &tagsSlave = node->getPropertyValue<QVector<GeometryTag>>("Tags slave");
        if(tagsSlave.isEmpty()==false)
        {
            QMap<GeometryTag,TopoDS_Shape> subShapesMap;
            for(int i=0; i<tagsSlave.size(); i++)
            {
                const GeometryTag &curTag = tagsSlave[i];
                if(subShapesMap.contains(curTag)) continue;
                const TopoDS_Shape &aSubShape = this->fromTagToShape(curTag);
                subShapesMap.insert(curTag,aSubShape);
            }
            emit requestApplyCustomColor(subShapesMap,Quantity_NOC_RED,true);
        }
    }
    else
    {
        const QVector<GeometryTag> &tags = node->getPropertyValue<QVector<GeometryTag>>("Tags");
        emit requestResetCustomColors(true);
        QMap<GeometryTag,TopoDS_Shape> subShapesMap;
        for(int i=0; i<tags.size(); i++)
        {
            const GeometryTag &curTag = tags[i];
            if(subShapesMap.contains(curTag)) continue;
            const TopoDS_Shape &aSubShape = this->fromTagToShape(curTag);
            subShapesMap.insert(curTag,aSubShape);
        }
        Quantity_NameOfColor aColor = graphicsTools::getModelFeatureColor(nodeType);
        emit requestApplyCustomColor(subShapesMap,aColor,true);
    }
}

/*
//! ----------------------
//! function: changeColor
//! details:
//! ----------------------
void SimulationManager::changeColor()
{
    //! -------------------------
    //! this version uses "Tags"
    //! -------------------------
    //cout<<"SimulationManager::changeColor()->____function called____"<<endl;
    QModelIndex index = mySelectionModel->currentIndex();
    if(!index.isValid()) return;

    QVariant options;

    SimulationNodeClass *node = myModel->itemFromIndex(index)->data(Qt::UserRole).value<SimulationNodeClass*>();
    SimulationNodeClass::nodeType theFamily = node->getFamily();
    SimulationNodeClass::nodeType theType = node->getType();

    ListOfShape list1, list2;
    Quantity_NameOfColor color1, color2;

    switch(theFamily)
    {
    case SimulationNodeClass::nodeType_connection:
    {
        //! -------------------------------------------------------------
        //! handle the particular case of a contact pairs, which must be
        //! highlighted with two colors, one for the master and one for
        //! the slave face(s)
        //! -------------------------------------------------------------
        if(theType != SimulationNodeClass::nodeType_connection && theType != SimulationNodeClass::nodeType_connectionGroup)
        {
            QExtendedStandardItem* itemMasterTags = node->getPropertyItem("Tags master");
            QExtendedStandardItem* itemSlaveTags = node->getPropertyItem("Tags slave");
            color1 = MASTER_COLOR;
            color2 = SLAVE_COLOR;

            QVector<GeometryTag> vecLocs;
            GeometryTag loc;
            if(itemMasterTags!=NULL)
            {
                vecLocs = itemMasterTags->data(Qt::UserRole).value<Property>().getData().value<QVector<GeometryTag>>();
                std::vector<int> vecParentShapes;
                for(QVector<GeometryTag>::iterator it = vecLocs.begin(); it!=vecLocs.end(); ++it)
                {
                    loc = *it;
                    TopAbs_ShapeEnum type = loc.subShapeType;
                    int parentShapeIndex = loc.parentShapeNr;
                    int childShapeIndex = loc.subTopNr;
                    TopoDS_Shape shape;

                    vecParentShapes.push_back(parentShapeIndex);

                    switch(type)
                    {
                    case TopAbs_SOLID:
                        shape = mySimulationDataBase->bodyMap.value(parentShapeIndex);
                        break;
                    case TopAbs_FACE:
                        shape = mySimulationDataBase->MapOfBodyTopologyMap.value(parentShapeIndex).faceMap.FindKey(childShapeIndex);
                        break;
                    case TopAbs_EDGE:
                        shape = mySimulationDataBase->MapOfBodyTopologyMap.value(parentShapeIndex).edgeMap.FindKey(childShapeIndex);
                        break;
                    case TopAbs_VERTEX:
                        shape = mySimulationDataBase->MapOfBodyTopologyMap.value(parentShapeIndex).vertexMap.FindKey(childShapeIndex);
                        break;
                    }
                    //! -------------------------------------
                    //! check if the parent shape is visible
                    //! -------------------------------------
                    cout<<"____Analyzing master visibility____"<<endl;

                    vecParentShapes = tools::clearFromDuplicates(vecParentShapes);

                    QVector<int> newVecParentShapes = QVector<int>::fromStdVector(vecParentShapes);

                    QExtendedStandardItem *itemGeometryRoot = this->getTreeItem(SimulationNodeClass::nodeType_geometry);
                    QExtendedStandardItem *itemGeometryBody;
                    for(int i=0; i<itemGeometryRoot->rowCount(); i++)
                    {
                        itemGeometryBody = static_cast<QExtendedStandardItem*>(itemGeometryRoot->child(i,0));
                        int mapIndex = itemGeometryBody->data(Qt::UserRole).value<SimulationNodeClass*>()->getPropertyItem("Map index")
                                ->data(Qt::UserRole).value<Property>().getData().toInt();
                        if(newVecParentShapes.contains(mapIndex))
                        {
                            cout<<"____body index: "<<mapIndex<<"____"<<endl;
                            bool isVisible = itemGeometryBody->data(Qt::UserRole).value<SimulationNodeClass*>()->getPropertyItem("Visible")
                                    ->data(Qt::UserRole).value<Property>().getData().toBool();
                            if(isVisible)
                            {
                                list1.Append(shape);
                                break;
                            }
                        }
                    }
                }
                emit requestDisplayShapeCopy1(list1,color1);
            }
            if(itemSlaveTags!=NULL)
            {
                std::vector<int> vecParentShapes;
                vecLocs = itemSlaveTags->data(Qt::UserRole).value<Property>().getData().value<QVector<GeometryTag>>();
                for(QVector<GeometryTag>::iterator it = vecLocs.begin(); it!=vecLocs.end(); ++it)
                {
                    loc = *it;
                    TopAbs_ShapeEnum type = loc.subShapeType;
                    int parentShapeIndex = loc.parentShapeNr;
                    int childShapeIndex = loc.subTopNr;
                    TopoDS_Shape shape;
                    vecParentShapes.push_back(parentShapeIndex);
                    switch(type)
                    {
                    case TopAbs_SOLID:
                        shape = mySimulationDataBase->bodyMap.value(parentShapeIndex);
                        break;
                    case TopAbs_FACE:
                        shape = mySimulationDataBase->MapOfBodyTopologyMap.value(parentShapeIndex).faceMap.FindKey(childShapeIndex);
                        break;
                    case TopAbs_EDGE:
                        shape = mySimulationDataBase->MapOfBodyTopologyMap.value(parentShapeIndex).edgeMap.FindKey(childShapeIndex);
                        break;
                    case TopAbs_VERTEX:
                        shape = mySimulationDataBase->MapOfBodyTopologyMap.value(parentShapeIndex).vertexMap.FindKey(childShapeIndex);
                        break;
                    }
                    //! -------------------------------------
                    //! check if the parent shape is visible
                    //! -------------------------------------
                    vecParentShapes = tools::clearFromDuplicates(vecParentShapes);

                    //cout<<"____vecParentShapes: "<<vecParentShapes.size()<<"____"<<endl;
                    //for(int i=0; i<vecParentShapes.size(); i++)
                    //    cout<<"____vecParentShape["<<i<<"] = "<<vecParentShapes.at(i)<<"____"<<endl;

                    QVector<int> newVecParentShapes = QVector<int>::fromStdVector(vecParentShapes);

                    QExtendedStandardItem *itemGeometryRoot = this->getTreeItem(SimulationNodeClass::nodeType_geometry);
                    QExtendedStandardItem *itemGeometryBody;
                    for(int i=0; i<itemGeometryRoot->rowCount(); i++)
                    {
                        itemGeometryBody = static_cast<QExtendedStandardItem*>(itemGeometryRoot->child(i,0));
                        int mapIndex = itemGeometryBody->data(Qt::UserRole).value<SimulationNodeClass*>()->getPropertyItem("Map index")
                                ->data(Qt::UserRole).value<Property>().getData().toInt();

                        if(newVecParentShapes.contains(mapIndex))
                        {
                            cout<<"____body index: "<<mapIndex<<"____"<<endl;
                            bool isVisible = itemGeometryBody->data(Qt::UserRole).value<SimulationNodeClass*>()->getPropertyItem("Visible")
                                    ->data(Qt::UserRole).value<Property>().getData().toBool();
                            if(isVisible)
                            {
                                list2.Append(shape);
                                break;
                            }
                        }
                    }
                }
            }
            emit requestDisplayShapeCopy2(list2,color2);
        }
    }
        break;

    default:
    {
        //! ---------------------------------------------------------
        //! all the items requiring a color, apart from contact pair
        //! which is handled within the previous case
        //! ---------------------------------------------------------
        if(theType !=SimulationNodeClass::nodeType_root &&
                theType !=SimulationNodeClass::nodeType_geometry &&
                theType !=SimulationNodeClass::nodeType_connection &&
                theType !=SimulationNodeClass::nodeType_meshControl &&
                theType !=SimulationNodeClass::nodeType_coordinateSystem &&
                theType !=SimulationNodeClass::nodeType_coordinateSystems &&
                theType !=SimulationNodeClass::nodeType_coordinateSystem_global &&
                theType !=SimulationNodeClass::nodeType_structuralAnalysis &&
                theType !=SimulationNodeClass::nodeType_thermalAnalysis &&
                theType !=SimulationNodeClass::nodeType_structuralAnalysisSettings &&
                theType !=SimulationNodeClass::nodeType_thermalAnalysisSettings &&
                theType !=SimulationNodeClass::nodeType_StructuralAnalysisSolution &&
                theType !=SimulationNodeClass::nodeType_mapper &&
                theType !=SimulationNodeClass::nodeType_importedBodyScalar &&
                theType !=SimulationNodeClass::nodeType_OpenFoamScalarData &&
                theType !=SimulationNodeClass::nodeType_remotePoint &&
                theType !=SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Acceleration &&
                theType !=SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RotationalVelocity &&
                theType !=SimulationNodeClass::nodeType_solutionStructuralNodalDisplacement &&
                theType !=SimulationNodeClass::nodeType_solutionStructuralStress &&
                theType !=SimulationNodeClass::nodeType_solutionStructuralMechanicalStrain &&
                theType !=SimulationNodeClass::nodeType_solutionStructuralThermalStrain &&
                theType !=SimulationNodeClass::nodeType_solutionStructuralTotalStrain)
        {
            QExtendedStandardItem* itemTags = node->getPropertyItem("Tags");
            color1 = graphicsTools::getModelFeatureColor(theType);
            QVector<GeometryTag> vecLocs;
            if(itemTags!=Q_NULLPTR)
            {
                GeometryTag loc;
                vecLocs = itemTags->data(Qt::UserRole).value<Property>().getData().value<QVector<GeometryTag>>();
                for(QVector<GeometryTag>::iterator it = vecLocs.begin(); it!=vecLocs.end(); ++it)
                {
                    loc = *it;
                    TopAbs_ShapeEnum type = loc.subShapeType;
                    int parentShapeIndex = loc.parentShapeNr;
                    int childShapeIndex = loc.subTopNr;

                    TopoDS_Shape shape;
                    switch(type)
                    {
                    case TopAbs_SOLID:
                    {
                        //cout<<"SimulationManager::changeColor()->____function called for solids____"<<endl;
                        shape = mySimulationDataBase->bodyMap.value(parentShapeIndex);
                    }
                        break;

                    case TopAbs_FACE:
                    {
                        //cout<<"SimulationManager::changeColor()->____function called for faces____"<<endl;
                        shape = mySimulationDataBase->MapOfBodyTopologyMap.value(parentShapeIndex).faceMap.FindKey(childShapeIndex);
                    }
                        break;

                    case TopAbs_EDGE:
                    {
                        //cout<<"SimulationManager::changeColor()->____function called for edges____"<<endl;
                        shape = mySimulationDataBase->MapOfBodyTopologyMap.value(parentShapeIndex).edgeMap.FindKey(childShapeIndex);
                        //! ------------------------------------
                        //! check if the item is a mesh control
                        //! ------------------------------------
                        SimulationNodeClass *curNode = myTreeView->currentIndex().data(Qt::UserRole).value<SimulationNodeClass*>();
                        //SimulationNodeClass *curNode = myModel->itemFromIndex(myTreeView->selectionModel()->currentIndex())->data(Qt::UserRole).value<SimulationNodeClass*>();
                        bool isMeshControl = curNode->getType() == SimulationNodeClass::nodeType_meshEdgeSize? true:false;
                        options.setValue(isMeshControl);
                    }
                        break;

                    case TopAbs_VERTEX:
                    {
                        //cout<<"SimulationManager::changeColor()->____function called for vertexes____"<<endl;
                        SimulationNodeClass *curNode = myTreeView->currentIndex().data(Qt::UserRole).value<SimulationNodeClass*>();
                        //SimulationNodeClass *curNode = myModel->itemFromIndex(myTreeView->selectionModel()->currentIndex())->data(Qt::UserRole).value<SimulationNodeClass*>();
                        shape = mySimulationDataBase->MapOfBodyTopologyMap.value(parentShapeIndex).vertexMap.FindKey(childShapeIndex);
                        //! ------------------------------------------------------------------------------
                        //! check if the item is a mesh control: if yes provide the radius of the pinball
                        //! ------------------------------------------------------------------------------
                        if(curNode->getType() == SimulationNodeClass::nodeType_meshVertexSize)
                        {
                            //cout<<"SimulationManager::changeColor()->____function called for a vertex mesh control____"<<endl;
                            double radius = curNode->getPropertyItem("Pinball")->data(Qt::UserRole).value<Property>().getData().toDouble();
                            options.setValue(radius);
                        }
                    }
                        break;
                    }
                    //! ---------------------------------------------------
                    //! such a shape could be null after a model refresh
                    //! (full change of the geometry, or simply an update)
                    //! ---------------------------------------------------
                    if(!shape.IsNull()) list1.Append(shape);
                }
            }
        }
        else
        {
            list1.Clear();
            list2.Clear();
        }
    }
        break;
    }
    emit requestDisplayShapeCopy(list1,list2,color1,color2,options);
    //cout<<"SimulationManager::changeColor()->____exiting____"<<endl;
}
*/

//! -----------------------------------------------
//! function: swapContact
//! details:  flip the master and the slave scopes
//! -----------------------------------------------
void SimulationManager::swapContact()
{
    cout<<"SimulationManager::swapContact()->____function called____"<<endl;
    QModelIndex index = myTreeView->currentIndex();
    SimulationNodeClass *node = myModel->itemFromIndex(index)->data(Qt::UserRole).value<SimulationNodeClass*>();
    if(node->getType()!=SimulationNodeClass::nodeType_connectionPair) return;

    //! ----------------------
    //! lock the node signals
    //! ----------------------
    node->getModel()->blockSignals(true);

    //! -------------------------
    //! the old master and slave
    //! -------------------------
    const QVector<GeometryTag> &scope1 = node->getPropertyValue<QVector<GeometryTag>>("Master");
    const QVector<GeometryTag> &scope2 = node->getPropertyValue<QVector<GeometryTag>>("Slave");

    QVariant data;

    //! ----------------------------------------------------
    //! change the scoping method into "Geometry selection"
    //! ----------------------------------------------------
    data.setValue(Property::ScopingMethod_GeometrySelection);
    Property property_scopingMethod("Scoping method",data,Property::PropertyGroup_Scope);
    data.setValue(property_scopingMethod);
    node->replaceProperty("Scoping method",property_scopingMethod);

    //! --------------------------------
    //! the old slave is the new master
    //! --------------------------------
    data.setValue(scope2);
    Property property_scopeMaster("Master",data,Property::PropertyGroup_Scope);
    node->replaceProperty("Master",property_scopeMaster);
    Property property_tagsMaster("Tags master",data,Property::PropertyGroup_Scope);
    node->replaceProperty("Tags master",property_tagsMaster);

    //! ---------------------------------
    //! the old masters is the new slave
    //! ---------------------------------
    data.setValue(scope1);
    Property property_scopeSlave("Slave",data,Property::PropertyGroup_Scope);
    node->replaceProperty("Slave",property_scopeSlave);
    Property property_tagsSlave("Tags slave",data,Property::PropertyGroup_Scope);
    node->replaceProperty("Tags slave",property_tagsSlave);

    //! ------------------------
    //! unlock the node signals
    //! ------------------------
    node->getModel()->blockSignals(false);

    //! ------------------
    //! update the colors
    //! ------------------
    ListOfShape shapeScope1,shapeScope2;

    //! QVector<GeometryTag> contains homogeneous shapes (same type)
    if(scope1.first().isParent)
    {
        for(QVector<GeometryTag>::const_iterator it = scope1.cbegin(); it!= scope1.cend(); ++it)
        {
            GeometryTag aLoc = *it;
            shapeScope1.Append(mySimulationDataBase->bodyMap.value(aLoc.parentShapeNr));
        }
    }
    else
    {
        switch(scope1.first().subShapeType)
        {
        case TopAbs_FACE:
            for(QVector<GeometryTag>::const_iterator it = scope1.cbegin(); it!= scope1.cend(); ++it)
            {
                GeometryTag aLoc = *it;
                int bodyIndex = aLoc.parentShapeNr;
                int subShapeNr = aLoc.subTopNr;
                shapeScope1.Append(mySimulationDataBase->MapOfBodyTopologyMap.value(bodyIndex).faceMap.FindKey(subShapeNr));
            }
            break;
        case TopAbs_EDGE:
            for(QVector<GeometryTag>::const_iterator it = scope1.cbegin(); it!= scope1.cend(); ++it)
            {
                GeometryTag aLoc = *it;
                int bodyIndex = aLoc.parentShapeNr;
                int subShapeNr = aLoc.subTopNr;
                shapeScope1.Append(mySimulationDataBase->MapOfBodyTopologyMap.value(bodyIndex).edgeMap.FindKey(subShapeNr));
            }
            break;
        case TopAbs_VERTEX:
            for(QVector<GeometryTag>::const_iterator it = scope1.cbegin(); it!= scope1.cend(); ++it)
            {
                GeometryTag aLoc = *it;
                int bodyIndex = aLoc.parentShapeNr;
                int subShapeNr = aLoc.subTopNr;
                shapeScope1.Append(mySimulationDataBase->MapOfBodyTopologyMap.value(bodyIndex).edgeMap.FindKey(subShapeNr));
            }
            break;
        }
    }

    //! QVector<GeometryTag> contains homogeneous shapes (same type)
    if(scope2.first().isParent)
    {
        for(QVector<GeometryTag>::const_iterator it = scope2.begin(); it!= scope2.end(); ++it)
        {
            GeometryTag aLoc = *it;
            shapeScope2.Append(mySimulationDataBase->bodyMap.value(aLoc.parentShapeNr));
        }
    }
    else
    {
        switch(scope2.first().subShapeType)
        {
        case TopAbs_FACE:
            for(QVector<GeometryTag>::const_iterator it = scope2.cbegin(); it!= scope2.cend(); ++it)
            {
                GeometryTag aLoc = *it;
                int bodyIndex = aLoc.parentShapeNr;
                int subShapeNr = aLoc.subTopNr;
                shapeScope2.Append(mySimulationDataBase->MapOfBodyTopologyMap.value(bodyIndex).faceMap.FindKey(subShapeNr));
            }
            break;
        case TopAbs_EDGE:
            for(QVector<GeometryTag>::const_iterator it = scope2.cbegin(); it!= scope2.cend(); ++it)
            {
                GeometryTag aLoc = *it;
                int bodyIndex = aLoc.parentShapeNr;
                int subShapeNr = aLoc.subTopNr;
                shapeScope2.Append(mySimulationDataBase->MapOfBodyTopologyMap.value(bodyIndex).edgeMap.FindKey(subShapeNr));
            }
            break;
        case TopAbs_VERTEX:
            for(QVector<GeometryTag>::const_iterator it = scope2.cbegin(); it!= scope2.cend(); ++it)
            {
                GeometryTag aLoc = *it;
                int bodyIndex = aLoc.parentShapeNr;
                int subShapeNr = aLoc.subTopNr;
                shapeScope2.Append(mySimulationDataBase->MapOfBodyTopologyMap.value(bodyIndex).edgeMap.FindKey(subShapeNr));
            }
            break;
        }
    }

    emit requestDisplayShapeCopy(shapeScope2,shapeScope1,MASTER_COLOR,SLAVE_COLOR);
    emit requestDisplayShapeCopy1(shapeScope1,MASTER_COLOR);
    emit requestDisplayShapeCopy1(shapeScope2,SLAVE_COLOR);
}

//! ------------------------
//! function: startAnalysis
//! details:
//! ------------------------
bool SimulationManager::startAnalysis(const QString &projectDataDir)
{
    cout<<"SimulationManager::startAnalysis()->____function called: solution data directory: "<<projectDataDir.toStdString()<<"____"<<endl;

    //! ---------------------------------------------------------------
    //! retrieve the current item and check if it is a simulation root
    //! ---------------------------------------------------------------
    QModelIndex curIndex = myTreeView->currentIndex();
    SimulationNodeClass *curAnalysisRootNode = curIndex.data(Qt::UserRole).value<SimulationNodeClass*>();
    if(curAnalysisRootNode->isAnalysisRoot()==false)
    {
        QMessageBox::information(this,APPNAME,tr("Please select an analysis root item"));
        return false;
    }

    //! ---------------------------------------
    //! check if the active bodies have a mesh
    //! ---------------------------------------
    bool isMeshOk = false;
    for(QMap<int,TopoDS_Shape>::const_iterator it = mySimulationDataBase->bodyMap.cbegin(); it!=mySimulationDataBase->bodyMap.cend(); ++it)
    {
        int bodyIndex = it.key();
        if(mySimulationDataBase->MapOfIsActive.value(bodyIndex)==false) continue;
        const occHandle(MeshVS_DataSource) &aMeshDataSource = mySimulationDataBase->ArrayOfMeshDS.value(bodyIndex);
        if(aMeshDataSource.IsNull()==false) isMeshOk = true;
    }
    if(isMeshOk==false)
    {
        QMessageBox::information(this,APPNAME,tr("Not all the bodies have been meshed\nThe simulation cannot be started"));
        return false;
    }

    //! ----------------------------------
    //! switch the tab of the main window
    //! ----------------------------------
    emit requestSetActiveCentralTab("worksheetViewer");

    //! --------------------------------------------------------------------
    //! store the projectDataDir (..\..\<project name>_files) into the tree
    //! (the method callPostEngineEvaluateResult() (action 204) uses it)
    //! --------------------------------------------------------------------
    QVariant data;
    data.setValue(projectDataDir);
    Property prop_projectFileDir("Project files dir",data,Property::PropertyGroup_Information);
    data.setValue(prop_projectFileDir);

    //! -----------------------------
    //! retrieve the "Solution" item
    //! -----------------------------
    QStandardItem *itemSimulationRoot = myModel->itemFromIndex(curIndex);
    QStandardItem *itemAnalysisSettings = itemSimulationRoot->child(0,0);
    QStandardItem *itemSolution = itemSimulationRoot->child(itemSimulationRoot->rowCount()-1);
    QStandardItem *itemSolutionInformation = itemSolution->child(0,0);
    SimulationNodeClass *nodeAnalysisSettigs = itemAnalysisSettings->data(Qt::UserRole).value<SimulationNodeClass*>();
    SimulationNodeClass *nodeSolution = itemSolution->data(Qt::UserRole).value<SimulationNodeClass*>();
    SimulationNodeClass *nodeSolutionInformation = itemSolutionInformation->data(Qt::UserRole).value<SimulationNodeClass*>();

    //! --------------------------------------------------------
    //! insert the "Project files dir" into the "Solution" item
    //! --------------------------------------------------------
    if(nodeSolution->getPropertyItem("Project files dir")!=Q_NULLPTR) nodeSolution->replaceProperty("Project files dir",prop_projectFileDir);
    else nodeSolution->addProperty(prop_projectFileDir,1);

    //! -------------------------------
    //! clear the solution information
    //! -------------------------------
    mainTreeTools::resetSolutionInformation(nodeSolutionInformation);

    //! ---------------------------
    //! read the number of threads
    //! ---------------------------
    int NbThreads = nodeAnalysisSettigs->getPropertyValue<int>("Number of threads");

    //! ---------------------------------------------------------------------------------
    //! create the directory for the solution data
    //! "SolutionData_<time tag of the analysis root>" within ../../<project name>_files
    //! ---------------------------------------------------------------------------------
    QDir curDir(projectDataDir);
    QString solutionDataDir = QString("SolutionData_")+curAnalysisRootNode->getPropertyValue<QString>("Time tag");
    if(!curDir.cd(solutionDataDir))
    {
        bool isDone = curDir.mkdir(solutionDataDir);
        if(isDone == false)
        {
            QMessageBox::critical(this,APPNAME,"Cannot start run: cannot create\n directory for analysis files",QMessageBox::Ok);
            return false;
        }
        curDir.cd(solutionDataDir);
    }
    cout<<"SimulationManager::startAnalysis()->____\"SolutionData\" full path: "<<curDir.absolutePath().toStdString()<<"____"<<endl;

    //! ----------------------------
    //! write the solver input file
    //! ----------------------------
    QString fileName = curDir.absolutePath()+"/input.inp";
    cout<<"SimulationManager::startAnalysis()->____input file full path: "<<fileName.toStdString()<<"____"<<endl;

    //! ---------------------------------------------------
    //! generate the mesh data sources for BC and contacts
    //! ---------------------------------------------------
    bool generateDual = false;
    //if(curAnalysisRootNode->getType()==SimulationNodeClass::nodeType_thermalAnalysis) generateDual = true;
    this->generateBoundaryConditionsMeshDS(generateDual);

    //! ----------------------------------
    //! work on the stack - single thread
    //! ----------------------------------
    writeSolverFileClass inputFileGenerator(mySimulationDataBase,static_cast<QExtendedStandardItem*>(itemSimulationRoot));
    QWidget *piw = tools::getWidgetByName("progressIndicator");
    QProgressIndicator *aProgressIndicator = static_cast<QProgressIndicator*>(piw);
    inputFileGenerator.setProgressIndicator(aProgressIndicator);
    inputFileGenerator.setName(fileName);
    bool isDone = inputFileGenerator.perform();

    if(isDone == false)
    {
        cout<<"SimulationManager::startAnalysis()->____input file not generated____"<<endl;
        return false;
    }

    //! ------------------------------------------
    //! enable the items with boundary conditions
    //! skip the "Solution" item (the last child)
    //! ------------------------------------------
    for(int i=0; i<itemSimulationRoot->rowCount()-1; i++) itemSimulationRoot->child(i,0)->setEnabled(true);

    //! ---------------------------------
    //! retrieve the final analysis time
    //! ---------------------------------
    SimulationNodeClass *nodeAnalysisSettings = itemAnalysisSettings->data(Qt::UserRole).value<SimulationNodeClass*>();
    CustomTableModel *tabData = nodeAnalysisSettings->getTabularDataModel();
    double endTime = tabData->dataRC(tabData->rowCount()-1,1).toDouble();

    cout<<"SimulationManager::startAnalysis()->____required analysis end time: "<<endTime<<"____"<<endl;

    //! ------------------------
    //! update the progress bar
    //! ------------------------
    QProgressEvent *e = new QProgressEvent(QProgressEvent_Init,0,int(100.0*endTime),0,"Running CCX",QProgressEvent_None,-1,-1,-1,"Running CCX");
    QApplication::postEvent(aProgressIndicator,e);
    QApplication::processEvents();

    //! ---------------------------------------------------------------------------------------
    //! create the solver manager - put the progress indicator into the constructor. To do ...
    //! ---------------------------------------------------------------------------------------
    CCXSolverManager1 *CCXsm = new CCXSolverManager1(this);

    CCXsm->setNbProcessors(NbThreads);
    CCXsm->setInputFile(fileName);
    CCXsm->setProgressIndicator(aProgressIndicator);

    disconnect(CCXsm,SIGNAL(CCXRunFinished()),myTimer,SLOT(stop()));
    connect(CCXsm,SIGNAL(CCXRunFinished()),myTimer,SLOT(stop()));

    disconnect(CCXsm,SIGNAL(CCXRunFinished()),this,SLOT(configureAndStartPostEngine()));
    connect(CCXsm,SIGNAL(CCXRunFinished()),this,SLOT(configureAndStartPostEngine()));

    // gives some problems in case of "fast" analysis - hang on CCXConsoleToFile::perform()
    //disconnect(CCXsm,SIGNAL(CCXRunFinished()),this,SLOT(retrieveSolverInfo()));
    //connect(CCXsm,SIGNAL(CCXRunFinished()),this,SLOT(retrieveSolverInfo()));

    disconnect(CCXsm,&CCXSolverManager1::CCXRunFinished,aProgressIndicator,&QProgressIndicator::hide);
    connect(CCXsm,&CCXSolverManager1::CCXRunFinished,aProgressIndicator,&QProgressIndicator::hide);
    CCXsm->start();

    //! ---------------------------------------------------------
    //! configure and start the internal timer
    //! retrieve the update interval from "Solution information"
    //! ---------------------------------------------------------
    double updateInterval = nodeSolutionInformation->getPropertyValue<double>("Update interval");
    int updateInterval_ms = int(updateInterval);
    myTimer->setInterval(updateInterval_ms*1000);
    disconnect(CCXsm,SIGNAL(CCXRunFinished()),myTimer,SLOT(stop()));
    connect(CCXsm,SIGNAL(CCXRunFinished()),myTimer,SLOT(stop()));
    myTimer->start();

    //! ---------------------------------
    //! set the current running analysis
    //! ---------------------------------
    myCurrentRunningAnalysis = itemSimulationRoot;
    Global::status().curRunningAnalysis = itemSimulationRoot;

    //! --------------------------------------
    //! clean the tree from mesh data sources
    //! --------------------------------------
    this->deleteDataSourcesFromModel();
    return true;
}


//! ------------------------------
//! function: configurePostEngine
//! details:
//! ------------------------------
void SimulationManager::configurePostEngine()
{
    cout<<"SimulationManager::configurePostEngine()->____function called____"<<endl;
}

//! --------------------------------
//! function: startPostEngine
//! details:  unused for the moment
//! --------------------------------
void SimulationManager::startPostEngine()
{
    cout<<"SimulationManager::startPostEngine()->____function called____"<<endl;
}

//! --------------------------------------
//! function: configureAndStartPostEngine
//! details:
//! --------------------------------------
void SimulationManager::configureAndStartPostEngine()
{
    cout<<"SimulationManager::configureAndStartPostEngine()->____function called____"<<endl;

    //! ------------------------------------------------------
    //! the run has ended: reset the current running analysis
    //! ------------------------------------------------------
    //myCurrentRunningAnalysis = Q_NULLPTR;
    //Global::status().curRunningAnalysis = Q_NULLPTR;

    //! -----------------
    //! set the database
    //! -----------------
    myPostEngine->setDataBase(mySimulationDataBase);

    //! -------------
    //! set the name
    //! -------------
    QStandardItem *itemSolution = myCurrentRunningAnalysis->child(myCurrentRunningAnalysis->rowCount()-1,0);
    SimulationNodeClass *nodeSolution = itemSolution->data(Qt::UserRole).value<SimulationNodeClass*>();
    QString solutionDirectory = nodeSolution->getPropertyValue<QString>("Project files dir")+"/SolutionData_"+nodeSolution->getPropertyValue<QString>("Parent time tag");

    if(solutionDirectory.isEmpty()) return;

    QString FRDfile = solutionDirectory+QString("/input.frd");

    if(QFile(FRDfile).exists()==false) return;

    myPostEngine->setResultsFile(FRDfile);

    //! --------------------------------------------------------------------------------------
    //! start the post engine: the FrdReader splits the .frd file into several files
    //! this is done up to the final available simulation time (handle unconverged solutions)
    //! --------------------------------------------------------------------------------------
    myPostEngine->perform();

    //! ------------------------------------------------------
    //! the run has ended: reset the current running analysis
    //! ------------------------------------------------------
    //myCurrentRunningAnalysis = Q_NULLPTR;
    //Global::status().curRunningAnalysis = Q_NULLPTR;
}

//! ----------------------------------------------
//! function: handleSolutionComponentChanged
//! details:  immediately update the colored mesh
//! ----------------------------------------------
void SimulationManager::handleSolutionComponentChanged()
{
    cout<<"SimulationManager::handleSolutionComponentChanged()->____function called____"<<endl;

    QModelIndex index = myTreeView->currentIndex();
    QStandardItem *curItem = static_cast<QStandardItemModel*>(myTreeView->model())->itemFromIndex(index);

    SimulationNodeClass *curNode = curItem->data(Qt::UserRole).value<SimulationNodeClass*>();
    QStandardItem *postObjectItem = curNode->getPropertyItem("Post object");
    if(postObjectItem!=Q_NULLPTR)
    {
        cout<<"SimulationManager::handleSolutionComponentChanged()->____start updating post object and viewer____"<<endl;
        postObject curPostObject = postObjectItem->data(Qt::UserRole).value<Property>().getData().value<postObject>();

        //! delete from the viewer the previous result
        emit requestHideSingleResult(curPostObject);

        //! retrieve the new solution component
        //int component = curNode->getPropertyItem("Type ")->data(Qt::UserRole).value<Property>().getData().toInt();

        //! method postObject::update
        //curPostObject.update(mySimulationDataBase, component);

        //! display the new result
        //emit requestDisplayResult(curPostObject);
    }
}

//! -------------------------------
//! function: writeSolverInputFile
//! details:
//! -------------------------------
void SimulationManager::writeSolverInputFile()
{
    QStandardItem *curItem = myModel->itemFromIndex(myTreeView->currentIndex());
    SimulationNodeClass *curNode = curItem->data(Qt::UserRole).value<SimulationNodeClass*>();
    if(curNode->isAnalysisRoot()==false)
    {
        QMessageBox::information(this,APPNAME,tr("Please select an analysis root item"));
        return;
    }
    cout<<"SimulationManager::writeSolverInputFile()->____writing input file for Analysis root: \""<<curNode->getName().toStdString()<<"\"____"<<endl;
    QString selectedFilter;
    QString fileName = QFileDialog::getSaveFileName(this,"Save solver input file",tools::getWorkingDir(),INP_FILES,&selectedFilter,0);
    if(fileName.isEmpty()) return;

    //! ----------------------
    //! write the solver file
    //! ----------------------
    bool generateDual = false;
    if(curNode->getType()==SimulationNodeClass::nodeType_thermalAnalysis) generateDual = true;
    this->generateBoundaryConditionsMeshDS(generateDual);
    //this->createSimulationNode(SimulationNodeClass::nodeType_thermalAnalysisAdiabaticWall);

    //! -----------------------
    //! prepare the parameters
    //! -----------------------
    std::vector<void*> parameters;
    parameters.push_back((void*)(mySimulationDataBase));
    parameters.push_back((void*)(QExtendedStandardItem*)(curItem));
    parameters.push_back((void*)(&fileName));

    myInputFileGenerator->setParameters(parameters);
    myInputFileGenerator->start();
}

//! ---------------------------------------------------
//! function: resizeTabularData
//! details:  add or remove rows from the tabular data
//! ---------------------------------------------------
void SimulationManager::resizeTabularData()
{
    cout<<"SimulationManager::resizeTabularData()->____function called____"<<endl;

    //! access the "Analysis Settings" node/item
    QModelIndex index = myTreeView->currentIndex();
    if(index.data(Qt::UserRole).canConvert<SimulationNodeClass*>())
    {
        SimulationNodeClass *node = index.data(Qt::UserRole).value<SimulationNodeClass*>();

        int newNumberOfSteps = node->getPropertyItem("Number of steps")->data(Qt::UserRole).value<Property>().getData().toInt();
        QExtendedStandardItem *itemSimulationRoot = static_cast<QExtendedStandardItem*>(myModel->itemFromIndex(index.parent()));

        //! update the tabular data contained in the "Analysis settings" item
        QExtendedStandardItem *AnalysisSettingsItem = static_cast<QExtendedStandardItem*>(itemSimulationRoot->child(0,0));
        SimulationNodeClass *AnalysisSettingsNode = AnalysisSettingsItem->data(Qt::UserRole).value<SimulationNodeClass*>();

        //! here "node" and "AnalysisSettingsNode" are the same thing ...
        //! correct. To do ...
        //!
        CustomTableModel *theTabularDataModel = AnalysisSettingsNode->getTabularDataModel();

        //! the current number of rows in the table
        //! 09/01/2018 "-1" has been added since after the creation of the table
        //! an "initial" row with "Time = 0" has been added, so when the number
        //! of row is "2" the Number of steps is "1"
        int curNumberOfSteps = theTabularDataModel->rowCount()-1;

        //! append rows
        if(newNumberOfSteps>curNumberOfSteps)
        {
            int numberOfRowsToInsert = newNumberOfSteps-curNumberOfSteps;
            theTabularDataModel->insertRows(theTabularDataModel->rowCount(),numberOfRowsToInsert);
        }
        else //! remove rows
        {
            //! {0, 1, 2} => remove N=2 rows
            //! last = rowCount-1 = 2;
            //! first = last - N + 1
            int numberOfRowsToRemove = curNumberOfSteps-newNumberOfSteps;
            int last = theTabularDataModel->rowCount()-1;
            int first = last-numberOfRowsToRemove+1;
            theTabularDataModel->removeRows(first,last);
        }
    }
    else
    {
        QMessageBox::critical(this, tr("SimulationManager::resizeTabularData()"),tr("Error: node Analysis Settings not found"));
    }
}

//! -----------------------------------------------------------------------
//! function: HandleTabularData
//! details:  handle the tabular data. This is done here because the class
//!           has direct access to item "Analysis settings". Switches are
//!           removed/added by the DetailViewer
//! -----------------------------------------------------------------------
void SimulationManager::HandleTabularData()
{
    cout<<"SimulationManager::HandleTabularData()->____function called____"<<endl;
    SimulationNodeClass *theCurNode = myTreeView->currentIndex().data(Qt::UserRole).value<SimulationNodeClass*>();
    Property::defineBy theDefineBy = theCurNode->getPropertyValue<Property::defineBy>("Define by");

    //! ------------------------------------
    //! access the "Analysis settings" item
    //! ------------------------------------
    SimulationNodeClass *nodeAnalysisSettings = this->getAnalysisSettingsNodeFromCurrentItem();

    //! -----------------------------------------
    //! here handle the property change
    //! -----------------------------------------
    CustomTableModel *tabData = nodeAnalysisSettings->getTabularDataModel();

    int startColumn = mainTreeTools::calculateStartColumn(myTreeView);

    //! ----------------------------------
    //! check if the data in table exsist
    //! ----------------------------------
    if(theDefineBy==Property::defineBy_components)
    {
        //! -----------------------------------------------------------------------
        //! remove the column for magnitude and add the columns for the components
        //! -----------------------------------------------------------------------
        tabData->removeColumns(startColumn, 1, QModelIndex());

        QVector<QVariant> values;
        QVariant data;
        data.setValue(0.0);
        for(int i=0; i<tabData->rowCount(); i++) values.push_back(data);
        load load_componentX, load_componentY, load_componentZ;
        load_componentX.setData(values);
        load_componentY.setData(values);
        load_componentZ.setData(values);

        //! ---------------------------
        //! establish the type of load
        //! ---------------------------
        SimulationNodeClass::nodeType theNodeType = theCurNode->getType();
        switch(theNodeType)
        {
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Force:
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteForce:
            load_componentX.setType(Property::loadType_forceX);
            load_componentY.setType(Property::loadType_forceY);
            load_componentZ.setType(Property::loadType_forceZ);
            tabData->setLoadToInsert(load_componentX);
            tabData->insertColumns(startColumn,1);
            tabData->setLoadToInsert(load_componentY);
            tabData->insertColumns(startColumn+1,1);
            tabData->setLoadToInsert(load_componentZ);
            tabData->insertColumns(startColumn+2,1);
            break;
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Moment:
            load_componentX.setType(Property::loadType_momentX);
            load_componentY.setType(Property::loadType_momentY);
            load_componentZ.setType(Property::loadType_momentZ);
            tabData->setLoadToInsert(load_componentX);
            tabData->insertColumns(startColumn,1);
            tabData->setLoadToInsert(load_componentY);
            tabData->insertColumns(startColumn+1,1);
            tabData->setLoadToInsert(load_componentZ);
            tabData->insertColumns(startColumn+2,1);
            break;
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Acceleration:
            load_componentX.setType(Property::loadType_accelerationX);
            load_componentY.setType(Property::loadType_accelerationY);
            load_componentZ.setType(Property::loadType_accelerationZ);
            tabData->setLoadToInsert(load_componentX);
            tabData->insertColumns(startColumn,1);
            tabData->setLoadToInsert(load_componentY);
            tabData->insertColumns(startColumn+1,1);
            tabData->setLoadToInsert(load_componentZ);
            tabData->insertColumns(startColumn+2,1);
            break;
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RotationalVelocity:
            load_componentX.setType(Property::loadType_rotationalVelocityX);
            load_componentY.setType(Property::loadType_rotationalVelocityY);
            load_componentZ.setType(Property::loadType_rotationalVelocityZ);
            tabData->setLoadToInsert(load_componentX);
            tabData->insertColumns(startColumn,1);
            tabData->setLoadToInsert(load_componentY);
            tabData->insertColumns(startColumn+1,1);
            tabData->setLoadToInsert(load_componentZ);
            tabData->insertColumns(startColumn+2,1);
            break;
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteRotation:
        {
            bool b0, b1, b2;
            b0 = theCurNode->getOldXLoadDefinition()!=Property::loadDefinition_free? true:false;
            b1 = theCurNode->getOldYLoadDefinition()!=Property::loadDefinition_free? true:false;
            b2 = theCurNode->getOldZLoadDefinition()!=Property::loadDefinition_free? true:false;
            if((b0==true && b1==false && b2==false)||(b0==false && b1==true && b2==false)||(b0==false && b1==false && b2==true))
            {
                //! -----------------------
                //! add only one component
                //! -----------------------
                if(b0==true)    //! add "X component"
                {
                    load_componentX.setType(Property::loadType_remoteRotationX);
                    tabData->setLoadToInsert(load_componentX);
                    tabData->insertColumns(startColumn,1);
                }
                if(b1==true)    //! add "Y component"
                {
                    load_componentY.setType(Property::loadType_remoteRotationY);
                    tabData->setLoadToInsert(load_componentY);
                    tabData->insertColumns(startColumn,1);
                }
                if(b2==true)    //! add "Z component"
                {
                    load_componentZ.setType(Property::loadType_remoteRotationZ);
                    tabData->setLoadToInsert(load_componentZ);
                    tabData->insertColumns(startColumn,1);
                }
            }
            if((b0==true && b1 == true && b2==false)||(b0==false && b1 == true && b2==true)&&(b0==true && b1 == false && b2==true))
            {
                //! -------------------
                //! add two components
                //! -------------------
                if(b0==true && b1 == true)  //! add "X component" and "Y component"
                {
                    load_componentX.setType(Property::loadType_remoteRotationX);
                    tabData->setLoadToInsert(load_componentX);
                    tabData->insertColumns(startColumn,1);
                    load_componentX.setType(Property::loadType_remoteRotationY);
                    tabData->setLoadToInsert(load_componentX);
                    tabData->insertColumns(startColumn+1,1);
                }
                if(b1==true && b2 == true)  //! add "Y component" and "Z component"
                {
                    load_componentY.setType(Property::loadType_remoteRotationY);
                    tabData->setLoadToInsert(load_componentY);
                    tabData->insertColumns(startColumn,1);
                    load_componentZ.setType(Property::loadType_remoteRotationZ);
                    tabData->setLoadToInsert(load_componentZ);
                    tabData->insertColumns(startColumn+1,1);
                }
                if(b0==true && b2 == true)  //! add "X component" and "Z component"
                {
                    load_componentX.setType(Property::loadType_remoteRotationX);
                    tabData->setLoadToInsert(load_componentX);
                    tabData->insertColumns(startColumn,1);
                    load_componentZ.setType(Property::loadType_remoteRotationZ);
                    tabData->setLoadToInsert(load_componentZ);
                    tabData->insertColumns(startColumn+1,1);
                }
            }
            if(b0==true && b1 == true && b2 == true)
            {
                //! ---------------------
                //! add three components
                //! ---------------------
                load_componentX.setType(Property::loadType_remoteRotationX);
                tabData->setLoadToInsert(load_componentX);
                tabData->insertColumns(startColumn,1);
                load_componentY.setType(Property::loadType_remoteRotationY);
                tabData->setLoadToInsert(load_componentY);
                tabData->insertColumns(startColumn,1);
                load_componentZ.setType(Property::loadType_remoteRotationZ);
                tabData->setLoadToInsert(load_componentZ);
                tabData->insertColumns(startColumn+2,1);
            }
        }
            break;
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Displacement:
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteDisplacement:
        {
            bool b0, b1, b2;
            b0 = theCurNode->getOldXLoadDefinition()!=Property::loadDefinition_free? true:false;
            b1 = theCurNode->getOldYLoadDefinition()!=Property::loadDefinition_free? true:false;
            b2 = theCurNode->getOldZLoadDefinition()!=Property::loadDefinition_free? true:false;

            if((b0==true && b1==false && b2==false)||(b0==false && b1==true && b2==false)||(b0==false && b1==false && b2==true))
            {
                //! -----------------------
                //! add only one component
                //! -----------------------
                if(b0==true)    //! add "X component"
                {
                    load_componentX.setType(Property::loadType_displacementX);
                    tabData->setLoadToInsert(load_componentX);
                    tabData->insertColumns(startColumn,1);
                }
                if(b1==true)    //! add "Y component"
                {
                    load_componentY.setType(Property::loadType_displacementY);
                    tabData->setLoadToInsert(load_componentY);
                    tabData->insertColumns(startColumn,1);
                }
                if(b2==true)    //! add "Z component"
                {
                    load_componentZ.setType(Property::loadType_displacementZ);
                    tabData->setLoadToInsert(load_componentZ);
                    tabData->insertColumns(startColumn,1);
                }
            }
            if((b0==true && b1 == true && b2==false)||(b0==false && b1 == true && b2==true)&&(b0==true && b1 == false && b2==true))
            {
                //! -------------------
                //! add two components
                //! -------------------
                if(b0==true && b1 == true)  //! add "X component" and "Y component"
                {
                    load_componentX.setType(Property::loadType_displacementX);
                    tabData->setLoadToInsert(load_componentX);
                    tabData->insertColumns(startColumn,1);
                    load_componentX.setType(Property::loadType_displacementY);
                    tabData->setLoadToInsert(load_componentX);
                    tabData->insertColumns(startColumn+1,1);
                }
                if(b1==true && b2 == true)  //! add "Y component" and "Z component"
                {
                    load_componentY.setType(Property::loadType_displacementY);
                    tabData->setLoadToInsert(load_componentY);
                    tabData->insertColumns(startColumn,1);
                    load_componentZ.setType(Property::loadType_displacementZ);
                    tabData->setLoadToInsert(load_componentZ);
                    tabData->insertColumns(startColumn+1,1);
                }
                if(b0==true && b2 == true)  //! add "X component" and "Z component"
                {
                    load_componentX.setType(Property::loadType_displacementX);
                    tabData->setLoadToInsert(load_componentX);
                    tabData->insertColumns(startColumn,1);
                    load_componentZ.setType(Property::loadType_displacementZ);
                    tabData->setLoadToInsert(load_componentZ);
                    tabData->insertColumns(startColumn+1,1);
                }
            }
            if(b0==true && b1 == true && b2 == true)
            {
                //! ---------------------
                //! add three components
                //! ---------------------
                load_componentX.setType(Property::loadType_displacementX);
                tabData->setLoadToInsert(load_componentX);
                tabData->insertColumns(startColumn,1);
                load_componentY.setType(Property::loadType_displacementY);
                tabData->setLoadToInsert(load_componentY);
                tabData->insertColumns(startColumn,1);
                load_componentZ.setType(Property::loadType_displacementZ);
                tabData->setLoadToInsert(load_componentZ);
                tabData->insertColumns(startColumn+2,1);
            }
        }
            break;
        }
    }
    else
    {
        cout<<"SimulationManager::HandleTabularData()->____switched to vector____"<<endl;
        //! ---------------------------------------------------------------------------
        //! theDefineBy==defineBy_vector
        //! remove the columns for the components and add the column for the magnitude
        //! ---------------------------------------------------------------------------
        int count=0;
        if(theCurNode->getType() == SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Displacement ||
                theCurNode->getType() == SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteDisplacement ||
                theCurNode->getType() == SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteRotation)
        {
            if(theCurNode->getOldXLoadDefinition()!=Property::loadDefinition_free)
            {
                count++;
                cout<<"____found X component to be removed____"<<endl;
            }
            if(theCurNode->getOldYLoadDefinition()!=Property::loadDefinition_free)
            {
                count++;
                cout<<"____found Y component to be removed____"<<endl;
            }
            if(theCurNode->getOldZLoadDefinition()!=Property::loadDefinition_free)
            {
                count++;
                cout<<"____found Z component to be removed____"<<endl;
            }
        }
        else
        {
            count=3;
        }
        cout<<"SimulationManager::HandleTabularData()->____removing: "<<count<<" columns____"<<endl;

        //tabData->removeColumns(startColumn, count, QModelIndex());
        if(count!=0)
        {
            tabData->removeColumns(startColumn, count, QModelIndex());
        }
        //cout<<"SimulationManager::HandleTabularData()->____number of columns after removal: "<<tabData->columnCount()<<"____"<<endl;
        //cout<<"SimulationManager::HandleTabularData()->____start adding the magnitude____"<<endl;

        QVector<QVariant> values;
        load load_magnitude(values,Property::loadType_none);

        //! ---------------------------
        //! establish the type of load
        //! ---------------------------
        SimulationNodeClass::nodeType theNodeType = theCurNode->getType();
        switch(theNodeType)
        {
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Force:
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteForce:
            load_magnitude.setType(Property::loadType_forceMagnitude);
            break;
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Acceleration:
            load_magnitude.setType(Property::loadType_accelerationMagnitude);
            break;
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RotationalVelocity:
            load_magnitude.setType(Property::loadType_rotationalVelocityMagnitude);
            break;
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Moment:
            load_magnitude.setType(Property::loadType_momentMagnitude);
            break;
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Displacement:
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteDisplacement:
            load_magnitude.setType(Property::loadType_displacementMagnitude);
            break;
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteRotation:
            load_magnitude.setType(Property::loadType_remoteRotationMagnitude);
            break;
        }

        tabData->setLoadToInsert(load_magnitude);

        //! ---------------------------------------------------------
        //! "Magnitude" => one single column inserted into the table
        //! ---------------------------------------------------------
        tabData->insertColumns(startColumn,1);
    }

    emit requestTabularData(this->getAnalysisSettingsItemFromCurrentItem()->index());

    QList<int> N1;
    N1 << TABULAR_DATA_STEP_END_TIME_COLUMN << mainTreeTools::getColumnsToRead(myTreeView);

    cout<<"\\--------------------------------------------------------\\"<<endl;
    for(int n=0; n<N1.length(); n++) cout<<"\\ N (by maintreetools) = "<<N1.at(n)<<endl;
    cout<<"\\--------------------------------------------------------\\"<<endl;

    emit requestShowGraph(tabData,N1);
}

//! ------------------------------------------------------------------------
//! function: getInsertionRow
//! details:  when adding a simulation setup item, that item must be placed
//!           between the "Analysis settings" item and the "Solution" item.
//!           The function finds the right row for inserting the new item
//! ------------------------------------------------------------------------
int SimulationManager::getInsertionRow() const
{
    QModelIndex theCurIndex = myTreeView->currentIndex();
    QStandardItem *theCurItem = myModel->itemFromIndex(theCurIndex);
    SimulationNodeClass* theCurNode = theCurItem->data(Qt::UserRole).value<SimulationNodeClass*>();
    int insertionRow;
    if(theCurNode->isAnalysisRoot())
    {
        insertionRow = theCurItem->rowCount()-1;    //! the (-1) inserts before the "Solution" item
    }
    else
    {
        insertionRow = theCurItem->parent()->rowCount();
        SimulationNodeClass *parentNode = theCurItem->parent()->data(Qt::UserRole).value<SimulationNodeClass*>();
        if(parentNode->isAnalysisRoot()) insertionRow--;
    }
    cout<<"____insertion row: "<<insertionRow<<"____"<<endl;
    return insertionRow;
}

//! -------------------------------------------------
//! function: getAnalysisSettingsNodeFromCurrentItem
//! details:
//! -------------------------------------------------
SimulationNodeClass* SimulationManager::getAnalysisSettingsNodeFromCurrentItem() const
{
    QStandardItem *itemAnalysisSettings = this->getAnalysisSettingsItemFromCurrentItem();
    if(itemAnalysisSettings==Q_NULLPTR) return Q_NULLPTR;
    return itemAnalysisSettings->data(Qt::UserRole).value<SimulationNodeClass*>();
}

//! -------------------------------------------------
//! function: getAnalysisSettingsItemFromCurrentItem
//! details:
//! -------------------------------------------------
QExtendedStandardItem* SimulationManager::getAnalysisSettingsItemFromCurrentItem() const
{
    cout<<"SimulationManager::getAnalysisSettingsItemFromCurrentItem()->____function called____"<<endl;

    QModelIndex currentModelIndex = myTreeView->currentIndex();
    QStandardItem *curItem = this->getModel()->itemFromIndex(currentModelIndex);
    SimulationNodeClass *curNode = currentModelIndex.data(Qt::UserRole).value<SimulationNodeClass*>();

    //! ----------------------------------------------
    //! case 1: the current item is a simulation root
    //! ----------------------------------------------
    if(curNode->isAnalysisRoot())
    {
        QStandardItem *item = curItem->child(0,0);
        return static_cast<QExtendedStandardItem*>(item);
    }
    //! ---------------------------------------------------------
    //! case 2: the current item is a child of a simulation root
    //! ---------------------------------------------------------
    if(curNode->isSimulationSetUpNode() || curNode->isAnalysisSettings())
    {
        QStandardItem *item = curItem->parent()->child(0,0);
        return static_cast<QExtendedStandardItem*>(item);
    }

    //! ---------------------------------------------------
    //! case 3: the current item is a post processing item
    //! ---------------------------------------------------
    if(curNode->isAnalysisResult() || curNode->isSolutionInformation())
    {
        QStandardItem *item = curItem->parent()->parent()->child(0,0);
        return static_cast<QExtendedStandardItem*>(item);
    }
    return Q_NULLPTR;
}

//! -----------------------------------------------------
//! function: handleFilmCoefficientLoadDefinitionChanged
//! details:
//! -----------------------------------------------------
void SimulationManager::handleFilmCoefficientLoadDefinitionChanged(const QString &textData)
{
    //! ----------------------------------------------------------------
    //! avoids ping-pong - reconnection at the end of this function [*]
    //! ----------------------------------------------------------------
    TableWidget *tableWidget = static_cast<TableWidget*>(tools::getWidgetByName("messagesAndLoadsWidget"));
    DetailViewer *detailViewer = static_cast<DetailViewer*>(tools::getWidgetByName("detailViewer"));
    disconnect(tableWidget,SIGNAL(requestUpdateDetailViewer(QModelIndex,QModelIndex,QVector<int>)),detailViewer,SLOT(updateDetailViewerFromTabularData(QModelIndex,QModelIndex,QVector<int>)));

    //! --------------------------
    //! retrieve the tabular data
    //! --------------------------
    SimulationNodeClass *nodeAnalysisSettings = this->getAnalysisSettingsNodeFromCurrentItem();
    CustomTableModel *tabData = nodeAnalysisSettings->getTabularDataModel();

    //! ------------------------------------------------------
    //! generate the "constant" and the "tabular" load vector
    //! ------------------------------------------------------
    QVector<QVariant> dataVecConst, dataVecTabular;
    dataVecConst.push_back(QVariant(0.0));
    dataVecTabular.push_back(QVariant(0.0));
    double curValue = textData.toDouble();
    for(int k=1;k<tabData->rowCount();k++)
    {
        dataVecConst.push_back(QVariant(curValue));
        dataVecTabular.push_back(QVariant(0.0));
    }

    //! --------------------------------------------------------------
    //! if the loadDefinition is "Constant" set the following values:
    //! time = 0 => load value = 0
    //! time = 1 => load value = <what written into the editor>
    //! time = 2 => load value = <what written into the editor>
    //! [...]
    //! --------------------------------------------------------------
    SimulationNodeClass *theCurNode = myTreeView->currentIndex().data(Qt::UserRole).value<SimulationNodeClass*>();
    Property::loadDefinition theLoadDefinition = theCurNode->getPropertyValue<Property::loadDefinition>("Film coefficient");

    int startColumn = mainTreeTools::calculateStartColumn(myTreeView);

    if(theLoadDefinition==Property::loadDefinition_constant)
    {
        //! -------------------------------------------------------------------------------------------------------------------
        //! this line should been removed for the following reason:
        //! 1) the current implementation of insertColumns() within the CustomTableModel class is actually is based
        //!    on the QVector::replace() function, by consequence it performs "insert and replace"
        //! 2) since the tabular data model is linked to the QtCharts::QXYModelMapper, if one removes a column
        //!    the warning QtCharts::QXYModelMapperPrivate::initializeXYFromModel "Invalid Y coordinate index in model mapper."
        //!    is generated
        //! -------------------------------------------------------------------------------------------------------------------
        tabData->removeColumns(startColumn,1,QModelIndex());

        //! ------------------------------------
        //! handle a "constant" definition load
        //! ------------------------------------
        load aLoad(dataVecConst,Property::loadType_thermalConvectionFilmCoefficientMagnitude);
        tabData->setLoadToInsert(aLoad);
        tabData->insertColumns(startColumn,1,QModelIndex());
        theCurNode->updateOldMagnitude(Property::loadDefinition_constant);
    }
    else
    {
        tabData->removeColumns(startColumn,1,QModelIndex());

        //! ---------------------------------------------------------
        //! if "tabular" has been selected reset the values (policy)
        //! ---------------------------------------------------------
        load aLoad(dataVecTabular,Property::loadType_thermalConvectionFilmCoefficientMagnitude);
        tabData->setLoadToInsert(aLoad);
        tabData->insertColumns(startColumn, 1, QModelIndex());
        theCurNode->updateOldMagnitude(Property::loadDefinition_tabularData);
    }

    //! ------------------------------
    //! see note [*] at the beginning
    //! ------------------------------
    connect(tableWidget,SIGNAL(requestUpdateDetailViewer(QModelIndex,QModelIndex,QVector<int>)),detailViewer,SLOT(updateDetailViewerFromTabularData(QModelIndex,QModelIndex,QVector<int>)));
}

//! -----------------------------------------------------------
//! function: handleReferenceTemperatureLoadDefinitionChanged
//! details:
//! -----------------------------------------------------------
void SimulationManager::handleReferenceTemperatureLoadDefinitionChanged(const QString &textData)
{
    //! ----------------------------------------------------------------
    //! avoids ping-pong - reconnection at the end of this function [*]
    //! ----------------------------------------------------------------
    TableWidget *tableWidget = static_cast<TableWidget*>(tools::getWidgetByName("messagesAndLoadsWidget"));
    DetailViewer *detailViewer = static_cast<DetailViewer*>(tools::getWidgetByName("detailViewer"));
    disconnect(tableWidget,SIGNAL(requestUpdateDetailViewer(QModelIndex,QModelIndex,QVector<int>)),detailViewer,SLOT(updateDetailViewerFromTabularData(QModelIndex,QModelIndex,QVector<int>)));

    //! --------------------------
    //! retrieve the tabular data
    //! --------------------------
    SimulationNodeClass *nodeAnalysisSettings = this->getAnalysisSettingsNodeFromCurrentItem();
    CustomTableModel *tabData = nodeAnalysisSettings->getTabularDataModel();

    //! ------------------------------------------------------
    //! generate the "constant" and the "tabular" load vector
    //! ------------------------------------------------------
    QVector<QVariant> dataVecConst, dataVecTabular;
    dataVecConst.push_back(QVariant(0.0));
    dataVecTabular.push_back(QVariant(0.0));
    double curValue = textData.toDouble();
    for(int k=1;k<tabData->rowCount();k++)
    {
        dataVecConst.push_back(QVariant(curValue));
        dataVecTabular.push_back(QVariant(0.0));
    }

    //! --------------------------------------------------------------
    //! if the loadDefinition is "Constant" set the following values:
    //! time = 0 => load value = 0
    //! time = 1 => load value = <what written into the editor>
    //! time = 2 => load value = <what written into the editor>
    //! [...]
    //! --------------------------------------------------------------
    SimulationNodeClass *theCurNode = myTreeView->currentIndex().data(Qt::UserRole).value<SimulationNodeClass*>();
    Property::loadDefinition theLoadDefinition = theCurNode->getPropertyValue<Property::loadDefinition>("Reference temperature");

    int startColumn = mainTreeTools::calculateStartColumn(myTreeView);

    if(theLoadDefinition==Property::loadDefinition_constant)
    {
        //! -------------------------------------------------------------------------------------------------------------------
        //! this line should been removed for the following reason:
        //! 1) the current implementation of insertColumns() within the CustomTableModel class is actually is based
        //!    on the QVector::replace() function, by consequence it performs "insert and replace"
        //! 2) since the tabular data model is linked to the QtCharts::QXYModelMapper, if one removes a column
        //!    the warning QtCharts::QXYModelMapperPrivate::initializeXYFromModel "Invalid Y coordinate index in model mapper."
        //!    is generated
        //! -------------------------------------------------------------------------------------------------------------------
        tabData->removeColumns(startColumn+1,1,QModelIndex());

        //! ------------------------------------
        //! handle a "constant" definition load
        //! ------------------------------------
        load aLoad(dataVecConst,Property::loadType_thermalConvectionReferenceTemperatureMagnitude);
        tabData->setLoadToInsert(aLoad);
        tabData->insertColumns(startColumn+1,1,QModelIndex());
        theCurNode->updateOldMagnitude(Property::loadDefinition_constant);
    }
    else
    {
        tabData->removeColumns(startColumn,1,QModelIndex());

        //! ---------------------------------------------------------
        //! if "tabular" has been selected reset the values (policy)
        //! ---------------------------------------------------------
        load aLoad(dataVecTabular,Property::loadType_thermalConvectionReferenceTemperatureMagnitude);
        tabData->setLoadToInsert(aLoad);
        tabData->insertColumns(startColumn, 1, QModelIndex());
        theCurNode->updateOldMagnitude(Property::loadDefinition_tabularData);
    }

    //! ------------------------------
    //! see note [*] at the beginning
    //! ------------------------------
    connect(tableWidget,SIGNAL(requestUpdateDetailViewer(QModelIndex,QModelIndex,QVector<int>)),detailViewer,SLOT(updateDetailViewerFromTabularData(QModelIndex,QModelIndex,QVector<int>)));
}

//! -----------------------------------------------
//! function: handleLoadMagnitudeDefinitionChanged
//! details:
//! -----------------------------------------------
void SimulationManager::handleLoadMagnitudeDefinitionChanged(const QString& textData)
{
    cout<<"SimulationManager::handleLoadMagnitudeDefinitionChanged()->____"<<textData.toStdString()<<"____"<<endl;

    //! ----------------------------------------------------------------
    //! avoids ping-pong - reconnection at the end of this function [*]
    //! ----------------------------------------------------------------
    TableWidget *tableWidget = static_cast<TableWidget*>(tools::getWidgetByName("messagesAndLoadsWidget"));
    DetailViewer *detailViewer = static_cast<DetailViewer*>(tools::getWidgetByName("detailViewer"));
    disconnect(tableWidget,SIGNAL(requestUpdateDetailViewer(QModelIndex,QModelIndex,QVector<int>)),detailViewer,SLOT(updateDetailViewerFromTabularData(QModelIndex,QModelIndex,QVector<int>)));

    //! --------------------------
    //! retrieve the tabular data
    //! --------------------------
    SimulationNodeClass *nodeAnalysisSettings = this->getAnalysisSettingsNodeFromCurrentItem();
    CustomTableModel *tabData = nodeAnalysisSettings->getTabularDataModel();

    //! ------------------------------------------------------
    //! generate the "constant" and the "tabular" load vector
    //! ------------------------------------------------------
    QVector<QVariant> dataVecConst, dataVecTabular;
    dataVecConst.push_back(QVariant(0.0));
    dataVecTabular.push_back(QVariant(0.0));
    double curValue = textData.toDouble();
    for(int k=1;k<tabData->rowCount();k++)
    {
        dataVecConst.push_back(QVariant(curValue));
        dataVecTabular.push_back(QVariant(0.0));
    }

    //! --------------------------------------------------------------
    //! if the loadDefinition is "Constant" set the following values:
    //! time = 0 => load value = 0
    //! time = 1 => load value = <what written into the editor>
    //! time = 2 => load value = <what written into the editor>
    //! [...]
    //! --------------------------------------------------------------
    SimulationNodeClass *theCurNode = myTreeView->currentIndex().data(Qt::UserRole).value<SimulationNodeClass*>();
    Property::loadDefinition theLoadDefinition = theCurNode->getPropertyValue<Property::loadDefinition>("Magnitude");

    Property::loadType aLoadType;
    SimulationNodeClass::nodeType aNodeType = theCurNode->getType();
    switch(aNodeType)
    {
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Acceleration:
        aLoadType = Property::loadType_accelerationMagnitude;
        break;
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RotationalVelocity:
        aLoadType = Property::loadType_rotationalVelocityMagnitude;
        break;
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Force:
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteForce:
        aLoadType = Property::loadType_forceMagnitude;
        break;
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Moment:
        aLoadType = Property::loadType_momentMagnitude;
        break;
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Displacement:
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteDisplacement:
        aLoadType = Property::loadType_displacementMagnitude;
        break;
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteRotation:
        aLoadType = Property::loadType_remoteRotationMagnitude;
        break;
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Pressure:
        aLoadType = Property::loadType_pressureMagnitude;
        break;
    case SimulationNodeClass::nodeType_thermalAnalysisTemperature:
        aLoadType = Property::loadType_temperatureMagnitude;
        break;
    case SimulationNodeClass::nodeType_thermalAnalysisThermalFlow:
        aLoadType = Property::loadType_thermalFlowMagnitude;
        break;
    case SimulationNodeClass::nodeType_thermalAnalysisThermalFlux:
        aLoadType = Property::loadType_thermalFluxMagnitude;
        break;
    }

    int startColumn = mainTreeTools::calculateStartColumn(myTreeView);

    if(theLoadDefinition==Property::loadDefinition_constant)
    {
        //! -------------------------------------------------------------------------------------------------------------------
        //! this line should been removed for the following reason:
        //! 1) the current implementation of insertColumns() within the CustomTableModel class is actually is based
        //!    on the QVector::replace() function, by consequence it performs "insert and replace"
        //! 2) since the tabular data model is linked to the QtCharts::QXYModelMapper, if one removes a column
        //!    the warning QtCharts::QXYModelMapperPrivate::initializeXYFromModel "Invalid Y coordinate index in model mapper."
        //!    is generated
        //! -------------------------------------------------------------------------------------------------------------------
        tabData->removeColumns(startColumn,1,QModelIndex());

        //! ------------------------------------
        //! handle a "constant" definition load
        //! ------------------------------------
        load aLoad(dataVecConst,aLoadType);
        tabData->setLoadToInsert(aLoad);
        tabData->insertColumns(startColumn,1,QModelIndex());
        theCurNode->updateOldMagnitude(Property::loadDefinition_constant);
    }
    else
    {
        tabData->removeColumns(startColumn,1,QModelIndex());

        //! ---------------------------------------------------------
        //! if "tabular" has been selected reset the values (policy)
        //! ---------------------------------------------------------
        load aLoad(dataVecTabular,aLoadType);
        tabData->setLoadToInsert(aLoad);
        tabData->insertColumns(startColumn, 1, QModelIndex());
        theCurNode->updateOldMagnitude(Property::loadDefinition_tabularData);
    }

    QList<int> N;
    N << TABULAR_DATA_STEP_END_TIME_COLUMN << mainTreeTools::getColumnsToRead(myTreeView);

    if(N.length()>=2)
    {
        //! ------------------------------------------------
        //! this means that at least a component is present
        //! ------------------------------------------------
        SimulationNodeClass *nodeAnalysisSettings = this->getAnalysisSettingsNodeFromCurrentItem();
        CustomTableModel *tabData = nodeAnalysisSettings->getTabularDataModel();
        emit requestShowGraph(tabData,N);
    }
    else requestClearGraph();

    //! ------------------------------
    //! see note [*] at the beginning
    //! ------------------------------
    connect(tableWidget,SIGNAL(requestUpdateDetailViewer(QModelIndex,QModelIndex,QVector<int>)),detailViewer,SLOT(updateDetailViewerFromTabularData(QModelIndex,QModelIndex,QVector<int>)));
}

//! ---------------------------------------------------------------------
//! function: handleLoadXDefinitionChanged
//! details:  add or remove the X component of the displacement/rotation
//! ---------------------------------------------------------------------
void SimulationManager::handleLoadXDefinitionChanged(const QString &textData)
{
    cout<<"SimulationManager::handleLoadXDefinitionChanged()->____function called____"<<textData.toStdString()<<"____"<<endl;

    //! ----------------------------------------------------------------
    //! avoids ping-pong - reconnection at the end of this function [*]
    //! ----------------------------------------------------------------
    TableWidget *tableWidget = static_cast<TableWidget*>(tools::getWidgetByName("messagesAndLoadsWidget"));
    DetailViewer *detailViewer = static_cast<DetailViewer*>(tools::getWidgetByName("detailViewer"));
    disconnect(tableWidget,SIGNAL(requestUpdateDetailViewer(QModelIndex,QModelIndex,QVector<int>)),detailViewer,SLOT(updateDetailViewerFromTabularData(QModelIndex,QModelIndex,QVector<int>)));

    //! --------------------------
    //! retrieve the tabular data
    //! --------------------------
    SimulationNodeClass *nodeAnalysisSettings = this->getAnalysisSettingsNodeFromCurrentItem();
    CustomTableModel *tabData = nodeAnalysisSettings->getTabularDataModel();

    //! ------------------------------------------------------
    //! generate the "constant" and the "tabular" load vector
    //! ------------------------------------------------------
    QVector<QVariant> dataVecConst, dataVecTabular;
    dataVecConst.push_back(QVariant(0.0));
    dataVecTabular.push_back(QVariant(0.0));
    double curValue = textData.toDouble();
    for(int k=1;k<tabData->rowCount();k++)
    {
        dataVecConst.push_back(QVariant(curValue));
        dataVecTabular.push_back(QVariant(0.0));
    }

    //! -----------------------
    //! calculate the position
    //! -----------------------
    int startColumn = mainTreeTools::calculateStartColumn(myTreeView);

    SimulationNodeClass *theCurNode = myTreeView->currentIndex().data(Qt::UserRole).value<SimulationNodeClass*>();
    if(theCurNode->getType()==SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Displacement ||
            theCurNode->getType()==SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteDisplacement ||
            theCurNode->getType()==SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteRotation)
    {
        Property::loadDefinition oldXdef =theCurNode->getOldXLoadDefinition();
        Property::loadDefinition loadDefinition_Xcomponent = theCurNode->getPropertyValue<Property::loadDefinition>("X component");

        //! -------------------------------------------
        //! remove the previous X component if present
        //! -------------------------------------------
        if(oldXdef!=Property::loadDefinition_free) tabData->removeColumns(startColumn,1,QModelIndex());

        if(loadDefinition_Xcomponent==Property::loadDefinition_tabularData)
        {
            //! --------------------------------------
            //! handle the "tabular" definition load
            //! --------------------------------------
            load load_componentX;
            load_componentX.setData(dataVecTabular);
            if(theCurNode->getType()==SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteRotation)
            {
                load_componentX.setType(Property::loadType_remoteRotationX);
            }
            else
            {
                load_componentX.setType(Property::loadType_displacementX);
            }
            tabData->setLoadToInsert(load_componentX);
            tabData->insertColumns(startColumn, 1, QModelIndex());
        }
        else if(loadDefinition_Xcomponent==Property::loadDefinition_constant)
        {
            //! --------------------------------------
            //! handle the "constant" definition load
            //! --------------------------------------
            load load_componentX;
            load_componentX.setData(dataVecConst);
            if(theCurNode->getType()==SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteRotation)
            {
                load_componentX.setType(Property::loadType_remoteRotationX);
            }
            else
            {
                load_componentX.setType(Property::loadType_displacementX);
            }
            tabData->setLoadToInsert(load_componentX);
            tabData->insertColumns(startColumn, 1, QModelIndex());
        }
        else if(loadDefinition_Xcomponent==Property::loadDefinition_free)
        {
            if(oldXdef!=Property::loadDefinition_free)
            {
                //! ------------------
                //! remove the column
                //! ------------------
                int startColumn = mainTreeTools::calculateStartColumn(myTreeView);
                tabData->removeColumns(startColumn,1,QModelIndex());
            }
        }
        theCurNode->updateOldLoadDefinition(1,loadDefinition_Xcomponent);
    }
    else
    {
        //! -------------------------------------------------------------
        //! what is not displacement/remote displacement/remote rotation
        //! -------------------------------------------------------------
        Property::loadDefinition loadDefinition_Xcomponent = theCurNode->getPropertyValue<Property::loadDefinition>("X component");

        Property::loadType aLoadType;
        SimulationNodeClass::nodeType aNodeType = theCurNode->getType();
        switch(aNodeType)
        {
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Acceleration:
            aLoadType = Property::loadType_accelerationX;
            break;
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RotationalVelocity:
            aLoadType = Property::loadType_rotationalVelocityX;
            break;
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Force:
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteForce:
            aLoadType = Property::loadType_forceX;
            break;
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Moment:
            aLoadType = Property::loadType_momentX;
            break;
        }

        if(loadDefinition_Xcomponent==Property::loadDefinition_constant)
        {
            tabData->removeColumns(startColumn,1,QModelIndex());
            //! --------------------------------------
            //! handle the "constant" definition load
            //! --------------------------------------
            load aLoad(dataVecConst,aLoadType);
            tabData->setLoadToInsert(aLoad);
            tabData->insertColumns(startColumn, 1, QModelIndex());
        }
        else
        {
            tabData->removeColumns(startColumn,1,QModelIndex());
            //! -------------------------------------
            //! handle the "tabular" definition load
            //! -------------------------------------
            load aLoad(dataVecTabular,aLoadType);
            tabData->setLoadToInsert(aLoad);
            tabData->insertColumns(startColumn, 1, QModelIndex());
        }
    }

    QList<int> N;
    N << TABULAR_DATA_STEP_END_TIME_COLUMN << mainTreeTools::getColumnsToRead(myTreeView);

    if(N.length()>=2)
    {
        //! ------------------------------------------------
        //! this means that at least a component is present
        //! ------------------------------------------------
        SimulationNodeClass *nodeAnalysisSettings = this->getAnalysisSettingsNodeFromCurrentItem();
        CustomTableModel *tabData = nodeAnalysisSettings->getTabularDataModel();
        emit requestShowGraph(tabData,N);
    }
    else requestClearGraph();

    //! ------------------------------
    //! see note [*] at the beginning
    //! ------------------------------
    connect(tableWidget,SIGNAL(requestUpdateDetailViewer(QModelIndex,QModelIndex,QVector<int>)),
            detailViewer,SLOT(updateDetailViewerFromTabularData(QModelIndex,QModelIndex,QVector<int>)));
}

//! -------------------------------------------------------------------------
//! function: handleLoadYDefinitionChanged
//! details:  add or remove the Y component of the displacement on the basis
//!           of what chosen in the detail viewer
//! -------------------------------------------------------------------------
void SimulationManager::handleLoadYDefinitionChanged(const QString &textData)
{
    cout<<"SimulationManager::handleLoadYDefinitionChanged()->____function called____"<<textData.toStdString()<<"____"<<endl;

    //! ----------------------------------------------------------------
    //! avoids ping-pong - reconnection at the end of this function [*]
    //! ----------------------------------------------------------------
    TableWidget *tableWidget = static_cast<TableWidget*>(tools::getWidgetByName("messagesAndLoadsWidget"));
    DetailViewer *detailViewer = static_cast<DetailViewer*>(tools::getWidgetByName("detailViewer"));
    disconnect(tableWidget,SIGNAL(requestUpdateDetailViewer(QModelIndex,QModelIndex,QVector<int>)),
               detailViewer,SLOT(updateDetailViewerFromTabularData(QModelIndex,QModelIndex,QVector<int>)));


    SimulationNodeClass *theCurNode = myTreeView->currentIndex().data(Qt::UserRole).value<SimulationNodeClass*>();

    if(theCurNode->getType()==SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Displacement ||
            theCurNode->getType()==SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteDisplacement ||
            theCurNode->getType()==SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteRotation)
    {
        Property::loadDefinition oldYdef = theCurNode->getOldYLoadDefinition();
        SimulationNodeClass *nodeAnalysisSettings = this->getAnalysisSettingsNodeFromCurrentItem();
        CustomTableModel *tabData = nodeAnalysisSettings->getTabularDataModel();
        Property::loadDefinition loadDefinition_Ycomponent = theCurNode->getPropertyValue<Property::loadDefinition>("Y component");

        if(loadDefinition_Ycomponent==Property::loadDefinition_tabularData)
        {
            //! --------------------------------------
            //! add the column for the Y displacement
            //! calculate the point of insertion
            //! --------------------------------------
            int startColumn = mainTreeTools::calculateStartColumn(myTreeView);

            //! ---------------------------------------------------------------------------------------
            //! check if X component is present (in order to calculate the right column for insertion)
            //! ---------------------------------------------------------------------------------------
            Property::loadDefinition loadDefinition_Xcomponent = theCurNode->getPropertyValue<Property::loadDefinition>("X component");

            if(loadDefinition_Xcomponent!=Property::loadDefinition_free) startColumn++;
            if(oldYdef!=Property::loadDefinition_free) tabData->removeColumns(startColumn,1,QModelIndex());

            //! -------------------------------------
            //! handle the "Tabular" definition load
            //! -------------------------------------
            QVector<QVariant> dataVec;
            QVariant value;
            value.setValue(0.0);
            for(int i=0; i<tabData->rowCount(); i++)
            {
                dataVec.push_back(value);
            }
            load load_componentY;
            load_componentY.setData(dataVec);

            if(theCurNode->getType()==SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteRotation)
            {
                load_componentY.setType(Property::loadType_remoteRotationY);
            }
            else
            {
                load_componentY.setType(Property::loadType_displacementY);
            }
            tabData->setLoadToInsert(load_componentY);
            tabData->insertColumns(startColumn, 1, QModelIndex());
        }
        if(loadDefinition_Ycomponent==Property::loadDefinition_constant)
        {
            //! --------------------------------------
            //! add the column for the Y displacement
            //! calculate the point of insertion
            //! --------------------------------------
            int startColumn = mainTreeTools::calculateStartColumn(myTreeView);

            //! ---------------------------------------------------------------------------------------
            //! check if X component is present (in order to calculate the right column for insertion)
            //! ---------------------------------------------------------------------------------------
            Property::loadDefinition loadDefinition_Xcomponent = theCurNode->getPropertyValue<Property::loadDefinition>("X component");
            if(loadDefinition_Xcomponent!=Property::loadDefinition_free) startColumn++;

            if(oldYdef!=Property::loadDefinition_free) tabData->removeColumns(startColumn,1,QModelIndex());

            //! ----------------------------------------
            //! define the load
            //! handle a "constant" definition load
            //! ----------------------------------------
            QVector<QVariant> dataVec;
            QVariant value;
            value.setValue(0.0);
            dataVec.push_back(value);
            for(int k=1;k<tabData->rowCount();k++)
            {
                value.setValue(textData.toDouble());
                dataVec.append(value);
            }
            load load_componentY;
            load_componentY.setData(dataVec);
            if(theCurNode->getType()==SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteRotation)
            {
                load_componentY.setType(Property::loadType_remoteRotationY);
            }
            else
            {
                load_componentY.setType(Property::loadType_displacementY);
            }
            tabData->setLoadToInsert(load_componentY);
            tabData->insertColumns(startColumn, 1, QModelIndex());
        }
        else if(loadDefinition_Ycomponent==Property::loadDefinition_free && oldYdef!=Property::loadDefinition_free)
        {
            //! -------
            //! remove
            //! -------
            int startColumn = mainTreeTools::calculateStartColumn(myTreeView);
            Property::loadDefinition loadDefinition_Xcomponent = theCurNode->getPropertyValue<Property::loadDefinition>("X component");
            if(loadDefinition_Xcomponent!=Property::loadDefinition_free) startColumn++;
            if(oldYdef!=Property::loadDefinition_free) tabData->removeColumns(startColumn,1,QModelIndex());
        }
        theCurNode->updateOldLoadDefinition(2,loadDefinition_Ycomponent);
    }
    else
    {
        //! -------------------------
        //! what is not displacement
        //! -------------------------
        SimulationNodeClass *nodeAnalysisSettings = this->getAnalysisSettingsNodeFromCurrentItem();
        CustomTableModel *tabData = nodeAnalysisSettings->getTabularDataModel();
        Property::loadDefinition loadDefinition_Ycomponent = theCurNode->getPropertyValue<Property::loadDefinition>("Y component");

        Property::loadType aLoadType;
        SimulationNodeClass::nodeType aNodeType = theCurNode->getType();
        switch(aNodeType)
        {
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Acceleration:
            aLoadType = Property::loadType_accelerationY;
            break;
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RotationalVelocity:
            aLoadType = Property::loadType_rotationalVelocityY;
            break;
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Force:
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteForce:
            aLoadType = Property::loadType_forceY;
            break;
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Moment:
            aLoadType = Property::loadType_momentY;
            break;
        }

        int startColumn = mainTreeTools::calculateStartColumn(myTreeView);

        if(loadDefinition_Ycomponent==Property::loadDefinition_constant)
        {
            tabData->removeColumns(startColumn+1,1,QModelIndex());

            //! --------------------------------------
            //! handle the "constant" definition load
            //! --------------------------------------
            QVector<QVariant> dataVec;
            QVariant value;
            value.setValue(0.0);
            dataVec.push_back(value);
            for(int k=1;k<tabData->rowCount();k++)
            {
                value.setValue(textData.toDouble());
                dataVec.push_back(value);
            }
            load aLoad(dataVec,aLoadType);
            tabData->setLoadToInsert(aLoad);
            tabData->insertColumns(startColumn+1, 1, QModelIndex());
        }
        else
        {
            //! -------------------------------------
            //! handle the "Tabular" definition load
            //! -------------------------------------
            tabData->removeColumns(startColumn+1,1,QModelIndex());
            QVector<QVariant> dataVec;
            QVariant value;
            value.setValue(0.0);
            for(int k=0;k<tabData->rowCount();k++)
            {
                dataVec.push_back(value);
            }
            load aLoad(dataVec,aLoadType);
            tabData->setLoadToInsert(aLoad);
            tabData->insertColumns(startColumn+1, 1, QModelIndex());
        }
    }

    QList<int> N;
    N << TABULAR_DATA_STEP_END_TIME_COLUMN << mainTreeTools::getColumnsToRead(myTreeView);

    if(N.length()>=2)
    {
        //! ------------------------------------------------
        //! this means that at least a component is present
        //! ------------------------------------------------
        SimulationNodeClass *nodeAnalysisSettings = this->getAnalysisSettingsNodeFromCurrentItem();
        CustomTableModel *tabData = nodeAnalysisSettings->getTabularDataModel();
        emit requestShowGraph(tabData,N);
    }
    else emit requestClearGraph();

    //! ------------------------------
    //! see note [*] at the beginning
    //! ------------------------------
    connect(tableWidget,SIGNAL(requestUpdateDetailViewer(QModelIndex,QModelIndex,QVector<int>)),
            detailViewer,SLOT(updateDetailViewerFromTabularData(QModelIndex,QModelIndex,QVector<int>)));
}

//! --------------------------------------------------------------------------
//! function: handleLoadZDefinitionChanged
//! details:  add or remove the Z component of the displacement on the basis
//!           of what chosen in the detail viewer
//! --------------------------------------------------------------------------
void SimulationManager::handleLoadZDefinitionChanged(const QString &textData)
{
    cout<<"SimulationManager::handleLoadZDefinitionChanged()->____function called____"<<textData.toStdString()<<"____"<<endl;

    //! ----------------------------------------------------------------
    //! avoids ping-pong - reconnection at the end of this function [*]
    //! ----------------------------------------------------------------
    TableWidget *tableWidget = static_cast<TableWidget*>(tools::getWidgetByName("messagesAndLoadsWidget"));
    DetailViewer *detailViewer = static_cast<DetailViewer*>(tools::getWidgetByName("detailViewer"));
    disconnect(tableWidget,SIGNAL(requestUpdateDetailViewer(QModelIndex,QModelIndex,QVector<int>)),
               detailViewer,SLOT(updateDetailViewerFromTabularData(QModelIndex,QModelIndex,QVector<int>)));

    SimulationNodeClass *theCurNode = myTreeView->currentIndex().data(Qt::UserRole).value<SimulationNodeClass*>();
    if(theCurNode->getType()==SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Displacement  ||
            theCurNode->getType()==SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteDisplacement ||
            theCurNode->getType()==SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteRotation)
    {
        Property::loadDefinition oldZdef = theCurNode->getOldZLoadDefinition();
        SimulationNodeClass *nodeAnalysisSettings = this->getAnalysisSettingsNodeFromCurrentItem();
        CustomTableModel *tabData = nodeAnalysisSettings->getTabularDataModel();
        Property::loadDefinition loadDefinition_Zcomponent = theCurNode->getPropertyValue<Property::loadDefinition>("Z component");

        cout<<"\\----------------------------------------------------------------------------------------\\"<<endl;
        cout<<"\\ SimulationManager::handleLoadMagnitudeDefinitionChanged()->____LOAD DEFINITION: "<<loadDefinition_Zcomponent<<"____"<<endl;
        cout<<"\\----------------------------------------------------------------------------------------\\"<<endl;

        if(loadDefinition_Zcomponent==Property::loadDefinition_tabularData)
        {
            //! --------------------------------------
            //! add the column for the Z displacement
            //! calculate the point of insertion
            //! --------------------------------------
            int startColumn = mainTreeTools::calculateStartColumn(myTreeView);

            //! --------------------------------
            //! check if X component is present
            //! --------------------------------
            Property::loadDefinition loadDefinition_Xcomponent = theCurNode->getPropertyValue<Property::loadDefinition>("X component");
            if(loadDefinition_Xcomponent!=Property::loadDefinition_free) startColumn++;

            //! --------------------------------
            //! check if Y component is present
            //! --------------------------------
            Property::loadDefinition loadDefinition_Ycomponent = theCurNode->getPropertyValue<Property::loadDefinition>("Y component");
            if(loadDefinition_Ycomponent!=Property::loadDefinition_free) startColumn++;

            if(oldZdef!=Property::loadDefinition_free) tabData->removeColumns(startColumn,1,QModelIndex());

            //! --------------------------
            //! handle the "Tabular" load
            //! --------------------------
            QVector<QVariant> dataVec;
            QVariant value;
            value.setValue(0.0);
            for(int i=0; i<tabData->rowCount(); i++)
            {
                dataVec.push_back(value);
            }
            load load_componentZ;
            load_componentZ.setData(dataVec);
            if(theCurNode->getType()==SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteRotation)
            {
                load_componentZ.setType(Property::loadType_remoteRotationZ);
            }
            else
            {
                load_componentZ.setType(Property::loadType_displacementZ);
            }
            tabData->setLoadToInsert(load_componentZ);
            tabData->insertColumns(startColumn, 1, QModelIndex());
        }
        else if(loadDefinition_Zcomponent==Property::loadDefinition_constant)
        {
            //! --------------------------------------
            //! add the column for the Z displacement
            //! calculate the point of insertion
            //! --------------------------------------
            int startColumn = mainTreeTools::calculateStartColumn(myTreeView);

            //! --------------------------------
            //! check if X component is present
            //! --------------------------------
            Property::loadDefinition loadDefinition_Xcomponent = theCurNode->getPropertyValue<Property::loadDefinition>("X component");
            if(loadDefinition_Xcomponent!=Property::loadDefinition_free) startColumn++;

            //! --------------------------------
            //! check if Y component is present
            //! --------------------------------
            Property::loadDefinition loadDefinition_Ycomponent = theCurNode->getPropertyValue<Property::loadDefinition>("Y component");
            if(loadDefinition_Ycomponent!=Property::loadDefinition_free) startColumn++;

            if(oldZdef!=Property::loadDefinition_free) tabData->removeColumns(startColumn,1,QModelIndex());

            //! --------------------------------------
            //! handle the "constant" definition load
            //! --------------------------------------
            QVector<QVariant> dataVec;
            QVariant value;
            value.setValue(0.0);
            dataVec.push_back(value);
            for(int k=1;k<tabData->rowCount();k++)
            {
                value.setValue(textData.toDouble());
                dataVec.push_back(value);
            }
            load load_componentZ;
            load_componentZ.setData(dataVec);
            if(theCurNode->getType()==SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteRotation)
            {
                load_componentZ.setType(Property::loadType_remoteRotationZ);
            }
            else
            {
                load_componentZ.setType(Property::loadType_displacementZ);
            }
            tabData->setLoadToInsert(load_componentZ);
            tabData->insertColumns(startColumn, 1, QModelIndex());
        }
        else if(loadDefinition_Zcomponent==Property::loadDefinition_free && oldZdef!=Property::loadDefinition_free)
        {
            //! -----------------------------------------
            //! remove the column for the X displacement
            //! calculate the point of removal
            //! -----------------------------------------
            int startColumn = mainTreeTools::calculateStartColumn(myTreeView);

            //! -----------------------------------------------------------------------------
            //! check if X component is present for calculating the right column for removal
            //! -----------------------------------------------------------------------------
            Property::loadDefinition loadDefinition_Xcomponent = theCurNode->getPropertyValue<Property::loadDefinition>("X component");
            if(loadDefinition_Xcomponent!=Property::loadDefinition_free) startColumn++;

            //! -----------------------------------------------------------------------------
            //! check if Y component is present for calculating the right column for removal
            //! -----------------------------------------------------------------------------
            Property::loadDefinition loadDefinition_Ycomponent = theCurNode->getPropertyValue<Property::loadDefinition>("Y component");
            if(loadDefinition_Ycomponent!=Property::loadDefinition_free) startColumn++;

            tabData->removeColumns(startColumn,1,QModelIndex());
        }
        theCurNode->updateOldLoadDefinition(3,loadDefinition_Zcomponent);
    }
    else
    {
        //! -------------------------
        //! what is not displacement
        //! -------------------------
        SimulationNodeClass *nodeAnalysisSettings = this->getAnalysisSettingsNodeFromCurrentItem();
        CustomTableModel *tabData = nodeAnalysisSettings->getTabularDataModel();
        Property::loadDefinition loadDefinition_Zcomponent = theCurNode->getPropertyItem("Z component")->data(Qt::UserRole).value<Property>().getData().value<Property::loadDefinition>();

        Property::loadType aLoadType;
        SimulationNodeClass::nodeType aNodeType = theCurNode->getType();
        switch(aNodeType)
        {
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Acceleration:
            aLoadType = Property::loadType_accelerationZ;
            break;
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RotationalVelocity:
            aLoadType = Property::loadType_rotationalVelocityZ;
            break;
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Force:
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteForce:
            aLoadType = Property::loadType_forceZ;
            break;
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Moment:
            aLoadType = Property::loadType_momentZ;
            break;
        }

        int startColumn = mainTreeTools::calculateStartColumn(myTreeView);

        if(loadDefinition_Zcomponent==Property::loadDefinition_constant)
        {
            tabData->removeColumns(startColumn+2,1,QModelIndex());

            //! --------------------------------------
            //! handle the "constant" definition load
            //! --------------------------------------
            QVector<QVariant> dataVec;
            QVariant value;
            value.setValue(0.0);
            dataVec.push_back(value);
            for(int k=1;k<tabData->rowCount();k++)
            {
                value.setValue(textData.toDouble());
                dataVec.push_back(value);
            }
            load aLoad(dataVec,aLoadType);
            tabData->setLoadToInsert(aLoad);
            tabData->insertColumns(startColumn+2, 1, QModelIndex());
        }
        else
        {
            tabData->removeColumns(startColumn+2,1,QModelIndex());

            //! --------------------------------------
            //! handle the "tabular" definition load
            //! --------------------------------------
            QVector<QVariant> dataVec;
            QVariant value;
            value.setValue(0.0);
            dataVec.push_back(value);
            for(int k=0;k<tabData->rowCount();k++)
            {
                dataVec.push_back(value);
            }
            load aLoad(dataVec,aLoadType);
            tabData->setLoadToInsert(aLoad);
            tabData->insertColumns(startColumn+2, 1, QModelIndex());
        }
    }

    QList<int> N;
    N << TABULAR_DATA_STEP_END_TIME_COLUMN << mainTreeTools::getColumnsToRead(myTreeView);

    if(N.length()>=2)
    {
        //! ------------------------------------------------
        //! this means that at least a component is present
        //! ------------------------------------------------
        SimulationNodeClass *nodeAnalysisSettings = this->getAnalysisSettingsNodeFromCurrentItem();
        CustomTableModel *tabData = nodeAnalysisSettings->getTabularDataModel();
        emit requestShowGraph(tabData,N);
    }
    else requestClearGraph();

    //! ------------------------------
    //! see note [*] at the beginning
    //! ------------------------------
    connect(tableWidget,SIGNAL(requestUpdateDetailViewer(QModelIndex,QModelIndex,QVector<int>)),
            detailViewer,SLOT(updateDetailViewerFromTabularData(QModelIndex,QModelIndex,QVector<int>)));
}

//! -------------------------------
//! function: createItemDescriptor
//! details:  under construction
//! -------------------------------
QString SimulationManager::createItemDescriptor() const
{
    //! init with an empty string
    QString textDescriptor ="";

    SimulationNodeClass *theNode = myTreeView->currentIndex().data(Qt::UserRole).value<SimulationNodeClass*>();

    //! create a text descriptor for boundary conditions items
    if((theNode->getFamily()==SimulationNodeClass::nodeType_structuralAnalysis && theNode->getType()!=SimulationNodeClass::nodeType_structuralAnalysisSettings) &&
            theNode->getType()!=SimulationNodeClass::nodeType_structuralAnalysisSettings)
    {
        //! retrieve the current time step
        SimulationNodeClass *nodeAnalysisSettings = this->getAnalysisSettingsNodeFromCurrentItem();
        int currentStepNumber = nodeAnalysisSettings->getPropertyValue<int>("Current step number");

        //! retrieve the tabular data
        CustomTableModel *tabData = nodeAnalysisSettings->getTabularDataModel();

        Q_UNUSED(currentStepNumber)
        Q_UNUSED(tabData)

        switch (theNode->getType())
        {
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_CylindricalSupport:
            textDescriptor = "Cylindrical support\n";
            break;
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryContidion_FixedSupport:
            textDescriptor = "Fixed support\n";
            break;
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_FrictionlessSupport:
            textDescriptor = "Frictionless support\n";
            break;
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Acceleration:
            textDescriptor = "Acceleration\n";
            break;
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RotationalVelocity:
            textDescriptor = "Rotational velocity\n";
            break;
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Force:
            textDescriptor = "Force\n";
            break;
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteForce:
            textDescriptor = "Remote force\n";
            break;
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Pressure:
            textDescriptor = "Pressure\n";
            break;
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Displacement:
            textDescriptor = "Displacement\n";
            break;
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteDisplacement:
            textDescriptor ="Remote displacement\n";
            break;
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteRotation:
            textDescriptor ="Remote rotation\n";
            break;
        }
    }
    return textDescriptor;
}

//! ------------------------------------
//! function: getTreeItemsRecursively
//! details:  scan recursively the tree
//! ------------------------------------
void SimulationManager::getTreeItemsRecursively(QStandardItemModel* model, QList<QExtendedStandardItem*> &items, QModelIndex parent)
{
    for(int r = 0; r<model->rowCount(parent); ++r)
    {
        QModelIndex index = model->index(r, 0, parent);
        // here is your applicable code
        QExtendedStandardItem *item = static_cast<QExtendedStandardItem*>(model->itemFromIndex(index));
        items.push_back(item);
        if(model->hasChildren(index))
        {
            getTreeItemsRecursively(model, items, index);
        }
    }
}

//! --------------------------------------------------------------------------
//! function: saveSimulationDataBase
//! details:  "savingDir" is the absolute path of the directory in which the
//!           "/SDB" directory for saving main tree items will be created
//!           "fileName" is relative
//! -----------------------------------------------------------------------
void SimulationManager::saveSimulationDataBase(const QString &savingDir, const QString &fileName)
{
    cout<<"SimulationManager::saveSimulationDataBase()->____file name: "<<fileName.toStdString()<<"____"<<endl;
    cout<<"SimulationManager::saveSimulationDataBase()->____saving tree items into folder: "<<savingDir.toStdString()<<"____"<<endl;

    //! ---------------------------
    //! [0] get all the tree items
    //! ---------------------------
    QList<QExtendedStandardItem*> items;
    this->getTreeItemsRecursively(myModel,items);

    //! ----------------------------------------------------------------------
    //! update the "Projct files dir" property for all the "Simulation" items
    //! ----------------------------------------------------------------------
    QString projectFilesDir = savingDir;
    projectFilesDir.chop(4);
    QVariant data;
    data.setValue(projectFilesDir);
    QStandardItem *itemRoot = Geometry_RootItem->parent();
    for(int n=0; n<itemRoot->rowCount(); n++)
    {
        QStandardItem *itemAnalysisRoot = itemRoot->child(n,0);
        SimulationNodeClass *curAnalysisRoot = itemAnalysisRoot->data(Qt::UserRole).value<SimulationNodeClass*>();
        if(curAnalysisRoot->isAnalysisRoot()==false) continue;

        QStandardItem *itemSolution = itemAnalysisRoot->child(itemAnalysisRoot->rowCount()-1,0);
        SimulationNodeClass *nodeSolution = itemSolution->data(Qt::UserRole).value<SimulationNodeClass*>();

        //! sanity check
        if(nodeSolution->isSolution()==false) continue;

        nodeSolution->replaceProperty("Project files dir",Property("Project files dir",data,Property::PropertyGroup_Information));
    }

    //! --------------------------------------------
    //! [1] set the directory for serializing items
    //! --------------------------------------------
    mySerializer->setSavingDirPath(savingDir);

    //! --------------------------------------------------
    //! [2] progress indicator - number of items to write
    //! --------------------------------------------------
    int N = items.length();

    //! ------------------------------------------------------------
    //! calculate the total number of meshVS_dataSource(s) to write
    //! ------------------------------------------------------------
    int NmeshDS = 0;
    for(int bodyIndex=1; bodyIndex<=mySimulationDataBase->bodyMap.size(); bodyIndex++)
    {
        if(!mySimulationDataBase->ArrayOfMeshDS.value(bodyIndex).IsNull())NmeshDS++;
        if(!mySimulationDataBase->ArrayOfMeshDS2D.value(bodyIndex).IsNull())NmeshDS++;

        //! ----------------------------------------------
        //! calculate the number of face mesh datasources
        //! ----------------------------------------------
        TopologyMap topologyMap = mySimulationDataBase->MapOfBodyTopologyMap.value(bodyIndex);
        for(int faceNr=1; faceNr<=topologyMap.faceMap.Extent();faceNr++)
        {
            if(!mySimulationDataBase->ArrayOfMeshDSOnFaces.getValue(bodyIndex,faceNr).IsNull()) NmeshDS++;
        }
    }

    //! -----------------------------
    //! add the final "packing" step
    //! -----------------------------
    N += NmeshDS + 1;

    //! ----------------------------------
    //! initialize the progress indicator
    //! ----------------------------------
    int written =  0;
    QWidget *progressIndicatorWidget = tools::getWidgetByName("progressIndicator");
    QProgressIndicator *pi = static_cast<QProgressIndicator*>(progressIndicatorWidget);

    //! -----------------------------
    //! deactivate the secondary bar
    //! -----------------------------
    pi->setSecondaryBarVisible(false);
    QProgressEvent *e = new QProgressEvent(QProgressEvent_Init,0,N,written,"Saving database");
    QApplication::postEvent(progressIndicatorWidget,e);
    QApplication::processEvents();

    //! ------------------
    //! saving tree items
    //! ------------------
    for(int i=0; i<items.length(); i++)
    {
        QExtendedStandardItem* curItem = items.at(i);
        QString text = QString("Saving tree item: '").append(curItem->data(Qt::DisplayRole).toString()).append("'");

        //! ------------------------------
        //! update the progress indicator
        //! ------------------------------
        e = new QProgressEvent(QProgressEvent_Update,0,N,written,text);
        QApplication::postEvent(progressIndicatorWidget,e);
        QApplication::processEvents();

        cout<<"* SAVING ITEM: "<<curItem->data(Qt::DisplayRole).toString().toStdString()<<endl;
        mySerializer->setItem(curItem);
        mySerializer->serialize(written);

        //! ------------------------------
        //! update the progress indicator
        //! ------------------------------
        written++;
        e = new QProgressEvent(QProgressEvent_Update,0,N,written,text);
        QApplication::postEvent(progressIndicatorWidget,e);
        QApplication::processEvents();
    }

    //! ----------------------------------------
    //! create a folder for storing the meshes
    //! the folder location is created at the
    //! same level of the SDB directory
    //! ----------------------------------------
    QDir curDir(savingDir);             //! QDir is initialized using the path of SDB
    curDir.cdUp();
    curDir.mkdir("MDB");
    curDir.cd("MDB");                   //! now current directory is "<savingDir>/MDB"
    curDir.mkdir("Volume");
    curDir.mkdir("Surface");
    curDir.mkdir("Faces");

    //! ----------------------------
    //! Saving the main body meshes
    //! ----------------------------
    curDir.cd("Volume");                //! current dir: "<savingDir>/MDB/Volume"
    cout<<"SimulationManager::saveSimulationDataBase()->____"<<curDir.absolutePath().toStdString()<<"____"<<endl;

    for(QMap<int,TopoDS_Shape>::iterator it = mySimulationDataBase->bodyMap.begin(); it!= mySimulationDataBase->bodyMap.end(); it++)
    {
        QApplication::processEvents();
        int bodyIndex = it.key();
        cout<<"SimulationManager::saveSimulationDataBase()->____saving volume mesh for body Nr: "<<bodyIndex<<"____"<<endl;
        if(!mySimulationDataBase->ArrayOfMeshDS.value(bodyIndex).IsNull())
        {
            //! ------------------------------
            //! update the progress indicator
            //! ------------------------------
            QString bodyName = mySimulationDataBase->MapOfBodyNames.value(bodyIndex);
            QString text = QString("Saving mesh: '")+bodyName+QString("'");
            e = new QProgressEvent(QProgressEvent_Update,0,N,written,text);
            QApplication::postEvent(progressIndicatorWidget,e);

            //! -------------------------------------------------------------
            //! create the name of the volume mesh file for the current body
            //! -------------------------------------------------------------
            QString meshFileName = curDir.absoluteFilePath(QString("%1").arg(bodyIndex));

            const occHandle(Ng_MeshVS_DataSource3D) &aMeshVS_DS = occHandle(Ng_MeshVS_DataSource3D)::
                    DownCast(mySimulationDataBase->ArrayOfMeshDS.value(bodyIndex));

            cout<<"SimulationManager::saveSimulationDataBase()->____writing volume mesh for body Nr: "<<bodyIndex<<" done____"<<endl;
            aMeshVS_DS->writeMesh(meshFileName,3);

            //! ------------------------------
            //! update the progress indicator
            //! ------------------------------
            written++;
            e = new QProgressEvent(QProgressEvent_Update,0,N,written,text);
            QApplication::postEvent(progressIndicatorWidget,e);
            QApplication::processEvents();
        }
    }

    //! ------------------------------------
    //! Saving the main body surface meshes
    //! ------------------------------------
    curDir.cdUp();
    curDir.cd("Surface");               //! now current directory is "<savingDir>/MDB/Surface"
    cout<<"SimulationManager::saveSimulationDataBase()->____"<<curDir.absolutePath().toStdString()<<"____"<<endl;

    for(QMap<int,TopoDS_Shape>::iterator it = mySimulationDataBase->bodyMap.begin(); it!= mySimulationDataBase->bodyMap.end(); it++)
    {
        QApplication::processEvents();
        int bodyIndex = it.key();
        if(!mySimulationDataBase->ArrayOfMeshDS2D.value(bodyIndex).IsNull())
        {
            ccout(QString("SimulationManager::saveSimulationDataBase()->____saving surface mesh for body Nr: %1____").arg(bodyIndex));
            QString bodyName = mySimulationDataBase->MapOfBodyNames.value(bodyIndex);
            QString text = QString("Saving surface mesh: '").append(bodyName).append("'");
            e = new QProgressEvent(QProgressEvent_Update,0,N,written,text);
            QApplication::postEvent(progressIndicatorWidget,e);

            QString meshFileName = curDir.absoluteFilePath(QString("%1_2D").arg(bodyIndex));
            const occHandle(Ng_MeshVS_DataSource2D) &aMeshVS_DS = occHandle(Ng_MeshVS_DataSource2D)::
                    DownCast(mySimulationDataBase->ArrayOfMeshDS2D.value(bodyIndex));
            aMeshVS_DS->writeMesh(meshFileName,2);

            //! ------------------------------
            //! update the progress indicator
            //! ------------------------------
            written++;
            QProgressEvent *e = new QProgressEvent(QProgressEvent_Update,0,N,written,text);
            QApplication::postEvent(progressIndicatorWidget,e);
            QApplication::processEvents();
        }
    }

    //! -----------------------
    //! Saving the face meshes
    //! -----------------------
    curDir.cdUp();
    curDir.cd("Faces");                 //! now current directory is "<savingDir>/MDB/Faces"
    cout<<"SimulationManager::saveSimulationDataBase()->____"<<curDir.absolutePath().toStdString()<<"____"<<endl;

    for(QMap<int,TopoDS_Shape>::iterator it = mySimulationDataBase->bodyMap.begin(); it!= mySimulationDataBase->bodyMap.end(); it++)
    {
        int bodyIndex = it.key();
        int Nfaces = mySimulationDataBase->MapOfBodyTopologyMap.value(bodyIndex).faceMap.Extent();
        for(int subTopNr=1; subTopNr<=Nfaces; subTopNr++)
        {
            if(mySimulationDataBase->ArrayOfMeshDSOnFaces.getValue(bodyIndex,subTopNr).IsNull()) continue;

            cout<<"SimulationManager::saveSimulationDataBase()->____saving face mesh Nr: "<<subTopNr<<" of body Nr: "<<bodyIndex<<"____"<<endl;
            QString bodyName = QString("%1_%2").arg(bodyIndex).arg(subTopNr);
            QString text = QString("Saving face mesh: '").append(bodyName).append("'");

            //! ------------------------------
            //! update the progress indicator
            //! ------------------------------
            e = new QProgressEvent(QProgressEvent_Update,0,N,written,text);
            QApplication::postEvent(progressIndicatorWidget,e);
            QApplication::processEvents();

            QString meshFileName = curDir.absoluteFilePath(QString("%1_%2").arg(bodyIndex).arg(subTopNr));
            const occHandle(Ng_MeshVS_DataSourceFace) &aMeshVS_DS = occHandle(Ng_MeshVS_DataSourceFace)::
                    DownCast(mySimulationDataBase->ArrayOfMeshDSOnFaces.getValue(bodyIndex,subTopNr));

            aMeshVS_DS->writeMesh(meshFileName,2);

            //! ------------------------------
            //! update the progress indicator
            //! ------------------------------
            written++;
            QProgressEvent *e = new QProgressEvent(QProgressEvent_Update,0,N,written,text);
            QApplication::postEvent(progressIndicatorWidget,e);
            QApplication::processEvents();
        }
    }

    //! -----------------------------------
    //! change the directory
    //! going up to "<project_files>/MDB/"
    //! -----------------------------------
    curDir.cdUp();
    cout<<"SimulationManager::saveSimulationDataBase()->____"<<curDir.absolutePath().toStdString()<<"____"<<endl;

    //! -------------------------------------
    //! go up to directory "<project_files>"
    //! -------------------------------------
    curDir.cdUp();
    cout<<"SimulationManager::saveSimulationDataBase()->____"<<curDir.absolutePath().toStdString()<<"____"<<endl;

    //! ------------------------------
    //! update the progress indicator
    //! ------------------------------
    e = new QProgressEvent(QProgressEvent_Update,0,N,written,QString("Packing files"));
    QApplication::postEvent(progressIndicatorWidget,e);
    QApplication::processEvents();
    QString _filesDir = curDir.absolutePath();

    //! ----------------------------------
    //! prepare final name of the archive
    //! going up to <myWorkingDir>
    //! ----------------------------------
    curDir.cdUp();
    cout<<"SimulationManager::saveSimulationDataBase()->____"<<curDir.absolutePath().toStdString()<<"____"<<endl;
    QString fullFileName = curDir.absolutePath()+("/")+fileName;

    //! --------------
    //! packing files
    //! --------------
    cout<<"SimulationManager::saveSimulationDataBase()->____packing files in dir: "<<_filesDir.toStdString()<<" into file: "<<fullFileName.toStdString()<<"____"<<endl;
    JlCompress::compressDir(fullFileName,_filesDir,true);

    //! -----------------------
    //! reset the progress bar
    //! -----------------------
    written++;
    e = new QProgressEvent(QProgressEvent_Reset,0,9999,0,"",QProgressEvent_Reset,0,9999,0);
    QApplication::postEvent(progressIndicatorWidget,e);
    QApplication::processEvents();

    //! -----------------------------
    //! reactivate the secondary bar
    //! -----------------------------
    pi->setSecondaryBarVisible(true);
}

//! ----------------------------------------------
//! function: nodeListBuilder
//! details:  rebuild the list of nodes from disk
//! ----------------------------------------------
QList<SimulationNodeClass*> SimulationManager::nodeListBuilder(const QString &savedFilesDir)
{
    cout<<"SimulationManager::nodeListBuilder->____function called____"<<endl;

    QDirIterator *dir1Iterator = new QDirIterator(savedFilesDir,QDirIterator::NoIteratorFlags);
    QList<SimulationNodeClass*> listOfNodes;

    //! ------------------------------
    //! update the progress indicator
    //! ------------------------------
    QWidget* theProgressIndicator = tools::getWidgetByName("progressIndicator");
    QProgressEvent *pgrEvent = new QProgressEvent(QProgressEvent_Init,0,100,0,"Searching the tree items");
    QApplication::postEvent(theProgressIndicator,pgrEvent);
    QApplication::processEvents();

    Sleep(1000);

    //! ---------------------------
    //! count the number of nodes
    //! ---------------------------
    int N=0;
    while(dir1Iterator->hasNext())
    {
        QString nodeFileName = dir1Iterator->next();
        QFileInfo info(nodeFileName);
        if(!info.isDir()) N++;
    }

    //! ------------------------------
    //! update the progress indicator
    //! ------------------------------
    pgrEvent = new QProgressEvent(QProgressEvent_Init,0,N,0,"Begin loading tree items");
    QApplication::postEvent(theProgressIndicator,pgrEvent);
    QApplication::processEvents();

    Sleep(1000);

    //! -----------------------------------
    //! scan the directory with node files
    //! -----------------------------------
    QDirIterator dirIterator(savedFilesDir, QDirIterator::NoIteratorFlags);
    int val = 0;
    while(dirIterator.hasNext())
    {
        QString nodeFileName = dirIterator.next();
        QFileInfo info(nodeFileName);
        if(!info.isDir())
        {
            QString nodeName, nodeTypeName;
            QVector<Property> vecProp = myDeserializer->deserialize(nodeFileName, nodeTypeName, nodeName);

            //! -----------
            //! workaround
            //! -----------
            int n = SimulationNodeClass::typeToInt(nodeTypeName);
            SimulationNodeClass::nodeType theType = static_cast<SimulationNodeClass::nodeType>(n);
            //! ---------------
            //! end workaround
            //! ---------------

            //! ----------------------------------------
            //! rebuild the node and append to the list
            //! ----------------------------------------
            SimulationNodeClass *aNode = new SimulationNodeClass(nodeName,theType,vecProp,this);

            //! ----------------------------------------------------------------------
            //! if the node is "Analysis settings" rebuild the tabular data from file
            //! ----------------------------------------------------------------------
            if(aNode->isAnalysisSettings())
            {
                QString timeTag = aNode->getPropertyValue<QString>("Time tag");
                QString tabularDataFileAbsolutePath = savedFilesDir+"/Tabular data/Tabular data_"+timeTag+".txt";

                cout<<"SimulationManager::nodeListBuilder()->____"<<tabularDataFileAbsolutePath.toStdString()<<"____"<<endl;
                QVector<load> vecLoad = myDeserializer->readTabularDataFromFile(tabularDataFileAbsolutePath);
                aNode->createTabularData(vecLoad, false);
            }
            listOfNodes<<aNode;

            //! ------------------------------
            //! update the progress indicator
            //! ------------------------------
            val++;
            pgrEvent = new QProgressEvent(QProgressEvent_Update,0,N,val,QString("Rebuilding tree: item '").append(nodeName).append("'"));
            QApplication::postEvent(theProgressIndicator,pgrEvent);
            QApplication::processEvents();
        }
    }
    pgrEvent = new QProgressEvent(QProgressEvent_Reset,0,9999,0,QString(""),QProgressEvent_Reset,0,9999,0);
    QApplication::postEvent(theProgressIndicator,pgrEvent);
    QApplication::processEvents();

    return listOfNodes;
}

//! --------------------------------------------------------
//! function: buildDataBaseFromDisk
//! details:  build the analysis database from disk
//!           fileName is the full path of the archive file
//! --------------------------------------------------------
void SimulationManager::buildDataBaseFromDisk(const QString &fileName)
{
    cout<<"SimulationManager::buildDataBaseFromDisk->____function called. File: "<<fileName.toStdString()<<"____"<<endl;

    QString project_filesDir = fileName;
    project_filesDir.chop(4);
    project_filesDir.append("_files");

    QDir curDir(project_filesDir);
    curDir.cd("SDB");

    QString SDB_Dir = curDir.absolutePath();
    cout<<"SimulationManager::buildDataBaseFromDisk->____reading from directory: "<<SDB_Dir.toStdString()<<"____"<<endl;

    //! --------------------------
    //! build the nodes from file
    //! --------------------------
    QList<SimulationNodeClass*> listOfNodes = this->nodeListBuilder(SDB_Dir);
    for(int i=0; i<listOfNodes.length(); i++)
    {
        //! -------------------------------------------------------------------
        //! do not add the "Import" node when the data base is built from disk
        //! -------------------------------------------------------------------
        if(listOfNodes.at(i)->getType()==SimulationNodeClass::nodeType_import) continue;
        cout<<"SimulationManager::buildDataBaseFromDisk->____: "<<listOfNodes.at(i)->getName().toStdString()<<"____"<<endl;
    }

    //! -------------------------------
    //! build the simulation data base
    //! -------------------------------
    cout<<"SimulationManager::buildDataBaseFromDisk()->____rebuilding database from file: "<<fileName.toStdString()<<"____"<<endl;

    //!geometryDataBase *aDB = new geometryDataBase(listOfNodes,this);
    //!meshDataBase *aDB = new meshDataBase(listOfNodes,this);
    simulationDataBase *aDB = new simulationDataBase(listOfNodes,fileName,this);

    cout<<"SimulationManager::buildDataBaseFromDisk->____data base rebuilt____"<<endl;

    //! --------------------------
    //! initialize the tree model
    //! --------------------------
    myModel = aDB->getModel();

    //! --------------------------------
    //! set the model for the QTreeView
    //! --------------------------------
    myTreeView->setModel(myModel);
    myTreeView->expandAll();

    //! -----------------------------
    //! allow for multiple selection
    //! -----------------------------
    myTreeView->setSelectionMode(QAbstractItemView::ExtendedSelection);

    //! -------------------------------------------
    //! retrieve the selection model from the tree
    //! -------------------------------------------
    mySelectionModel = myTreeView->selectionModel();

    //! ------------------------------
    //! initialize the private member
    //! ------------------------------
    mySimulationDataBase = aDB;

    //! -------------------
    //! retrieve the roots
    //! -------------------
    Geometry_RootItem = this->getTreeItem(SimulationNodeClass::nodeType_geometry);
    RemotePoint_RootItem = this->getTreeItem(SimulationNodeClass::nodeType_remotePointRoot);
    CoordinateSystems_RootItem = this->getTreeItem(SimulationNodeClass::nodeType_coordinateSystems);
    NamedSelection_RootItem = this->getTreeItem(SimulationNodeClass::nodeType_namedSelection);
    Connections_RootItem = this->getTreeItem(SimulationNodeClass::nodeType_connection);
    Mesh_RootItem = this->getTreeItem(SimulationNodeClass::nodeType_meshControl);

    //! -------------------------------------------------
    //! hide the dummy named selectio "Select from list"
    //! hide the dummy contact pair "Select from list"
    //! hide the dummy remote point "Select from list"
    //! -------------------------------------------------
    if(NamedSelection_RootItem!=Q_NULLPTR) myTreeView->setRowHidden(0, NamedSelection_RootItem->index(), true);
    if(Connections_RootItem!=Q_NULLPTR) myTreeView->setRowHidden(0, Connections_RootItem->index(), true);
    if(RemotePoint_RootItem!=Q_NULLPTR) myTreeView->setRowHidden(0, RemotePoint_RootItem->index(), true);

    //! ----------------------------------------
    //! update the "Project files dir" property
    //! for all the "Solution" items
    //! ----------------------------------------
    QVariant data;
    data.setValue(project_filesDir);
    QStandardItem *rootItem = Geometry_RootItem->parent();
    for(int n=0; n<rootItem->rowCount(); n++)
    {
        QStandardItem *itemAnalysisRoot = rootItem->child(n,0);
        SimulationNodeClass *nodeAnalysisRoot = itemAnalysisRoot->data(Qt::UserRole).value<SimulationNodeClass*>();
        if(nodeAnalysisRoot->isAnalysisRoot()==false) continue;
        QStandardItem *itemSolution = itemAnalysisRoot->child(itemAnalysisRoot->rowCount()-1,0);
        SimulationNodeClass *nodeSolution = itemSolution->data(Qt::UserRole).value<SimulationNodeClass*>();
        if(nodeSolution->isSolution()==false) continue;
        nodeSolution->replaceProperty("Project files dir",Property("Project files dir",data,Property::PropertyGroup_Information));
    }

    //! ----------------------------------------------------
    //! reconnect signals/slots for handling item changes
    //! ----------------------------------------------------
    QList<QExtendedStandardItem*> allItems;
    this->getTreeItemsRecursively(this->myModel,allItems);    
    for(QList<QExtendedStandardItem*>::iterator it = allItems.begin(); it!=allItems.end(); ++it)
    {
        QExtendedStandardItem *anItem = *it;
        SimulationNodeClass *aNode = anItem->data(Qt::UserRole).value<SimulationNodeClass*>();

        //! ----------------------------------------------------------------
        //! during reload, parse the mesh control items,
        //! Excluding the "Mesh" root
        //! ----------------------------------------------------------------
        //if(aNode->getFamily()==SimulationNodeClass::nodeType_meshControl &&
        //        aNode->getType()!=SimulationNodeClass::nodeType_meshControl)
        //{
        //    parser::parseItem(anItem);
        //}

        //! ----------------------------------------------------------------
        //! during reload, parse the simulation setup items
        //! excluding the "Structural analysis" and the "Analysis settings"
        //! ----------------------------------------------------------------
        //if(aNode->getFamily()==SimulationNodeClass::nodeType_structuralAnalysis &&
        //        aNode->getFamily()!=SimulationNodeClass::nodeType_structuralAnalysis &&
        //        aNode->getFamily()!=SimulationNodeClass::nodeType_structuralAnalysisSettings)
        //{
        //    parser::parseItem(anItem);
        //}
        //! --------------------------------------------------------------
        //! during the reload, parse the post processing items, excluding
        //! the "Solution" root, and the "Solution information"
        //! --------------------------------------------------------------
        //if(aNode->getFamily()==SimulationNodeClass::nodeType_StructuralAnalysisSolution &&
        //        aNode->getType()!=SimulationNodeClass::nodeType_StructuralAnalysisSolution &&
        //        aNode->getType()!=SimulationNodeClass::nodeType_StructuralAnalysisSolutionInformation)
        //{
        //    parser::parseItem(anItem);
        //}

        disconnect(aNode->getModel(),SIGNAL(itemChanged(QStandardItem*)),this,SLOT(handleItemChange(QStandardItem*)));  //! avoids multiple connections
        connect(aNode->getModel(),SIGNAL(itemChanged(QStandardItem*)),this,SLOT(handleItemChange(QStandardItem*)));
    }

    //! ------------------------------------------------------------
    //! connection - triggering "SimulationManager:::highlighter()"
    //! ------------------------------------------------------------
    disconnect(myTreeView,SIGNAL(clicked(QModelIndex)),this,SLOT(highlighter(QModelIndex)));
    connect(myTreeView,SIGNAL(clicked(QModelIndex)),this,SLOT(highlighter(QModelIndex)));

    //! -------------------------------------------
    //! set the mesh data base for the post engine
    //! -------------------------------------------
    myPostEngine->setDataBase(aDB);

    //! ---------------------------------------
    //! aim: rebuild the post processing items
    //! ---------------------------------------
    QStandardItem *modelRootItem = Geometry_RootItem->parent();
    for(int n = 0; n < modelRootItem->rowCount(); n++)
    {
        QStandardItem *curAnalysisRootItem = modelRootItem->child(n,0);
        SimulationNodeClass *curAnalysisRootNode = curAnalysisRootItem->data(Qt::UserRole).value<SimulationNodeClass*>();
        if(curAnalysisRootNode->isAnalysisRoot()==false) continue;

        //! -------------------------------------------------
        //! set the input file path base for the post engine
        //! -------------------------------------------------
        QStandardItem *curSolutionItem = curAnalysisRootItem->child(curAnalysisRootItem->rowCount()-1,0);
        SimulationNodeClass *curSolutionNode = curSolutionItem->data(Qt::UserRole).value<SimulationNodeClass*>();

        if(curSolutionNode->isSolution()==false) continue;
        if(curSolutionNode->getPropertyItem("Project files dir")==Q_NULLPTR) continue;

        //! -------------------------------------------------
        //! set the input file path base for the post engine
        //! -------------------------------------------------
        QString timeTag = curAnalysisRootNode->getPropertyValue<QString>("Time tag");
        QString resultsFilePath = project_filesDir+"/SolutionData_"+timeTag+"/input.frd";
        myPostEngine->setResultsFile(resultsFilePath);

        //! -------------------------------------------------------------------------------
        //! scan the children of "Solution", jumping over the first "Solution information"
        //! -------------------------------------------------------------------------------
        for(int i=1; i<curSolutionItem->rowCount(); i++)
        {
            QStandardItem *curPostProcessingItem = curSolutionItem->child(i,0);
            SimulationNodeClass *curPostProcessingNode = curPostProcessingItem->data(Qt::UserRole).value<SimulationNodeClass*>();

            //! ------------------
            //! working exception
            //! ------------------
            if(curPostProcessingNode->getType()==SimulationNodeClass::nodeType_solutionStructuralFatigueTool) continue;

            bool immediatelyDisplay = false;
            this->callPostEngineEvaluateResult_private(curPostProcessingItem,immediatelyDisplay);
        }

        /*
        //! -----------------------------------------------------------------------
        //! rebuild the output of the simulation monitor - disabled for the moment
        //! since each simulation has its own solver messages
        //! -----------------------------------------------------------------------
        QDirIterator aDirIterator(project_filesDir,QDirIterator::Subdirectories);
        while(aDirIterator.hasNext())
        {
            QString curSubDirPath = aDirIterator.next();
            //cout<<"____"<<curSubDirPath.toStdString()<<"____"<<endl;
            if(curSubDirPath.contains("SolutionData_",Qt::CaseSensitive)==false) continue;

            //! -----------------------------------------
            //! retrieve the time tag from the file name
            //! -----------------------------------------
            QString timeTag = curSubDirPath.split("/").last();
            timeTag = timeTag.split(("_")).last();

            //cout<<"____time tag: "<<timeTag.toStdString()<<"____"<<endl;

            //! ---------------------------------------------------------------
            //! navigate the main tree and retrieve the "Solution information"
            //! ---------------------------------------------------------------
            QStandardItem *itemRoot = Geometry_RootItem->parent();
            for(int n=0; n<itemRoot->rowCount(); n++)
            {
                QStandardItem *itemAnalysisRoot = itemRoot->child(n,0);
                SimulationNodeClass *nodeAnalysisRoot = itemAnalysisRoot->data(Qt::UserRole).value<SimulationNodeClass*>();
                if(nodeAnalysisRoot->isAnalysisRoot()==false) continue;

                //cout<<"____"<<nodeAnalysisRoot->getPropertyValue<QString>("Time tag").toStdString()<<"____"<<endl;

                if(nodeAnalysisRoot->getPropertyValue<QString>("Time tag")!=timeTag) continue;

                QStandardItem *itemSolution = itemAnalysisRoot->child(itemAnalysisRoot->rowCount()-1,0);
                SimulationNodeClass *solutionInformationNode = itemSolution->child(0,0)->data(Qt::UserRole).value<SimulationNodeClass*>();

                //cout<<"____"<<solutionInformationNode->getPropertyValue<QString>("Time tag").toStdString()<<"____"<<endl;

                //! -----------------------------------
                //! read the "RawSolverOutputFile.txt"
                //! -----------------------------------
                QFile solverOutputFile(curSubDirPath+"/RawSolverOutput.txt");
                if(solverOutputFile.exists())
                {
                    QString aString;
                    QString fileName = curSubDirPath+"/RawSolverOutput.txt";

                    //cout<<"____"<<fileName.toStdString()<<"____"<<endl;

                    FILE *fileRawSolverOutput = fopen(fileName.toStdString().c_str(),"r");
                    if(fileRawSolverOutput!=NULL)
                    {
                        for(;feof(fileRawSolverOutput)==0;)
                        {
                            char line [256];
                            fgets(line, sizeof line, fileRawSolverOutput);
                            aString.append(line);
                            cout<<line<<endl;
                        }

                        solutionInformationNode->getModel()->blockSignals(true);

                        QVariant data;
                        data.setValue(aString);
                        Property prop_solverOutput("Solver output",data,Property::PropertyGroup_Hidden);
                        if(solutionInformationNode->getPropertyItem("Solver output")==Q_NULLPTR)
                            solutionInformationNode->addProperty(prop_solverOutput);
                        else solutionInformationNode->replaceProperty("Solver output",prop_solverOutput);

                        solutionInformationNode->getModel()->blockSignals(false);
                    }
                    fclose(fileRawSolverOutput);
                }
            }
        }
        */

    }

    //! -----------------------------
    //! transfer DislayRole to array
    //! -----------------------------
    QStandardItem *itemGeometryRoot = this->getTreeItem(SimulationNodeClass::nodeType_geometry);
    int NbBodies = itemGeometryRoot->rowCount();
    for(int i=0; i<NbBodies; i++)
    {
        QStandardItem *itemBody = itemGeometryRoot->child(i,0);
        QString bodyName = itemBody->data(Qt::DisplayRole).toString();
        if(itemBody->data(Qt::UserRole).value<SimulationNodeClass*>()->getType()==SimulationNodeClass::nodeType_pointMass) continue;
        int mapIndex = itemBody->data(Qt::UserRole).value<SimulationNodeClass*>()->getPropertyValue<int>("Map index");
        mySimulationDataBase->MapOfBodyNames.insert(mapIndex,bodyName);
    }
}

//! -----------------------------
//! function: renameItem
//! details:  rename a tree item
//! -----------------------------
void SimulationManager::renameItem()
{
    cout<<"SimulationManager::renameItem()->____function called____"<<endl;

    QStandardItem *curItem = myModel->itemFromIndex(myTreeView->currentIndex());
    QString newName = curItem->data(Qt::DisplayRole).toString();
    if(!curItem->isEditable()) curItem->setEditable(true);
    if(newName.isEmpty() || newName.isNull())
    {
        cerr<<"SimulationManager::renameItem()->____invalid name: \"default\" used____"<<endl;
        newName = "default";
    }
    curItem->data(Qt::UserRole).value<SimulationNodeClass*>()->setName(newName);

    //! -----------------------------------------------------
    //! update the name of the node and make it not editable
    //! -----------------------------------------------------
    this->updateNodeName();
    curItem->setEditable(false);
}

//! ----------------------------------------------------------------------
//! function: updateNodeName
//! details:  [1] update the name old the node
//!           [2] update the property "Name" under the "Definition" group
//!           [3] update the array of names
//! ----------------------------------------------------------------------
void SimulationManager::updateNodeName()
{
    cout<<"SimulationManager::updateNodeName()->____function called: updating the name of a CAD part____"<<endl;
    QString itemName = myTreeView->currentIndex().data(Qt::DisplayRole).toString();

    //! ---------------------
    //! update the node name
    //! ---------------------
    SimulationNodeClass *theNode = myTreeView->currentIndex().data(Qt::UserRole).value<SimulationNodeClass*>();
    theNode->setName(itemName);

    //! -------------------------------------------------
    //! if the item is a body update the vector of names
    //! -------------------------------------------------
    if(theNode->getType()==SimulationNodeClass::nodeType_geometryBody)
    {
        //! ----> questione da indagare approfonditamente <---- //
        //const TopoDS_Shape &theShape = theNode->getPropertyItem("Geometry")->
        //        data(Qt::UserRole).value<Property>().getData().value<ListOfShape>().First();
        //int mapIndex = mySimulationDataBase->bodyMap.key(theShape);
        int mapIndex = theNode->getPropertyValue<int>("Map index");

        //! ----------------------------------------------
        //! update the property in the "Definition" group
        //! ----------------------------------------------
        QVariant data;
        data.setValue(itemName);
        theNode->replaceProperty("Name",Property("Name",data,Property::PropertyGroup_Definition));

        //! ----------------------------------------------
        //! update the "array-style" part of the database
        //! ----------------------------------------------
        mySimulationDataBase->MapOfBodyNames.insert(mapIndex,QString(itemName.toStdString().c_str()));
    }
}

//! ----------------------------------------------
//! function: interpolate
//! details:  function called by the context menu
//! ----------------------------------------------
void SimulationManager::interpolate()
{
    //! -----------------------------------------------------------------------
    //! this node is contained into the "Mapper" item
    //! (SimulationManager::::interpolate is triggered within the context menu
    //! -----------------------------------------------------------------------
    SimulationNodeClass *node = myTreeView->currentIndex().data(Qt::UserRole).value<SimulationNodeClass*>();

    //! -----------------------------------
    //! retrieve the "Step selection mode"
    //! -----------------------------------
    QExtendedStandardItem *item = node->getPropertyItem("Step selection mode");
    if(item!=NULL)
    {
        int mode = item->data(Qt::UserRole).value<Property>().getData().toInt();
        this->interpolatePrivate(mode);
    }
}

//! -----------------------------
//! function; interpolatePrivate
//! details:
//! -----------------------------
void SimulationManager::interpolatePrivate(int mode)
{
    cout<<"SimulationManager::interpolatePrivate()->____function called____"<<endl;;

    //! ----------------------------------------------------------
    //! store the current model item ("Mapper")
    //! the results of the interpolation will be appended to this
    //! ----------------------------------------------------------
    QStandardItem* curItem = static_cast<QStandardItemModel*>(myTreeView->model())->itemFromIndex(myTreeView->currentIndex());

    //! -----------------
    //! the current node
    //! -----------------
    SimulationNodeClass *node = myTreeView->currentIndex().data(Qt::UserRole).value<SimulationNodeClass*>();
    QString parentTimeTag = node->getPropertyValue<QString>("Time tag");

    //! ---------------------------------
    //! retrieve the progress indicator
    //! --------------------------------
    QProgressIndicator *aProgressIndicator = static_cast<QProgressIndicator*>(tools::getWidgetByName("progressIndicator"));
    QProgressEvent *prgEvent;

    //! -----------------------
    //! reset the progress bar
    //! -----------------------
    if(aProgressIndicator)
    {
        cout<<"____reset progress indicator____"<<endl;
        prgEvent = new QProgressEvent(QProgressEvent_Reset,0,9999,0,"",QProgressEvent_Reset,0,9999,0);
        QApplication::postEvent(aProgressIndicator,prgEvent);
        QApplication::processEvents();
    }

    //! ---------------------------
    //! retrieve the target bodies
    //! ---------------------------
    QVector<GeometryTag> vecLocs = node->getPropertyValue<QVector<GeometryTag>>("Tags");

    //! ---------------------------------------
    //! retrieve the number of remapping steps
    //! ---------------------------------------
    bool remapFlag = node->getPropertyValue<bool>("Remap");
    int NbRemappingSteps = 0;
    if(node->getPropertyItem("Remapping steps")!=NULL)
    {
        NbRemappingSteps = node->getPropertyValue<int>("Remapping steps");
    }

    //! -------------------------------
    //! retrieve the number of buckets
    //! -------------------------------
    int NbBucketsX = node->getPropertyValue<int>("X buckets");
    int NbBucketsY = node->getPropertyValue<int>("Y buckets");
    int NbBucketsZ = node->getPropertyValue<int>("Z buckets");

    //! ------------------------
    //! interpolation algorithm
    //! ------------------------
    int algo = node->getPropertyValue<int>("Algorithm");

    //! ------------------------------------------------------------------------------------
    //! vector of sources nodes, contains a vector of scalar values associated to each node
    //! ------------------------------------------------------------------------------------
    std::vector<Mapper3DClass::nodeSource> vecSourceNodes;

    //! ----------------------------------------
    //! number of times: "1" for case "default"
    //! ----------------------------------------
    int NbStep = 1;

    //! --------------------------
    //! list of times (by string)
    //! --------------------------
    QList<QString> listTimeS;

    //! ---------------
    //! measuring time
    //! ---------------
    std::chrono::time_point<std::chrono::system_clock> start, end;

    switch(mode)
    {
    case 4:
    {
        //! ---------------------------------
        //! the directory of the source data
        //! ---------------------------------
        QString directoryDataPath = node->getPropertyValue<QString>("Source directory");
        if(directoryDataPath=="") return;

        QDir curDir;
        curDir.cd(directoryDataPath);
        curDir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks);

        cout<<"SimulationManager::interpolatePrivate()->____working on directory "<<curDir.absolutePath().toStdString()<<"____"<<endl;
        QDirIterator dirIterator(curDir, QDirIterator::NoIteratorFlags);

        //! ---------------
        //! read the times
        //! ---------------
        int NbFiles = curDir.entryList().length();
        cout<<"SimulationManager::interpolatePrivate()->____number of \"time\" files: "<<NbFiles<<"____"<<endl;

        //! --------------
        //! send progress
        //! --------------
        if(aProgressIndicator!=Q_NULLPTR)
        {
            cout<<"____send first event____"<<endl;
            QString msg("Start reading files");
            prgEvent = new QProgressEvent(QProgressEvent_Init,0,NbFiles-1,0,msg,QProgressEvent_Init,0,NbFiles-1,0);
            QApplication::postEvent(aProgressIndicator,prgEvent);
            QApplication::processEvents();
            Sleep(500);
        }

        int c=0;
        while(dirIterator.hasNext())
        {
            dirIterator.next();
            if(dirIterator.fileInfo().isFile())
            {
                QString filePath = dirIterator.filePath();

                //! --------------
                //! send progress
                //! --------------
                if(aProgressIndicator!=Q_NULLPTR)
                {
                    QString msg = QString("Reading file ")+filePath;
                    prgEvent = new QProgressEvent(QProgressEvent_Init,0,NbFiles-1,c,msg,QProgressEvent_Init,0,NbFiles-1,c);
                    QApplication::postEvent(aProgressIndicator,prgEvent);
                    QApplication::processEvents();
                    Sleep(250);
                }

                //! -------------------------------------------
                //! check if the file name is valid
                //! the file must end with "_<time value>.txt"
                //! -------------------------------------------
                bool isSuffixValid;
                QString timeS = filePath.split("_").last();

                //! -------------------------
                //! remove the ".txt" suffix
                //! -------------------------
                timeS.chop(4);
                timeS.toDouble(&isSuffixValid);
                if(!isSuffixValid) continue;

                listTimeS<<timeS;

                //! -------------------
                //! read a source file
                //! -------------------
                FILE *f = fopen(filePath.toStdString().c_str(),"r");
                if(f==NULL) continue;
                if(c==0)
                {
                    cout<<"SimulationManager::interpolatePrivate()->____reading file: "<<filePath.toStdString()<<"____"<<endl;
                    double x,y,z,aVal;
                    for(;feof(f)==NULL;)
                    {
                        if(4==fscanf(f,"%lf%lf%lf%lf",&x,&y,&z,&aVal))
                        {
                            Mapper3DClass::nodeSource aNodeSource;
                            aNodeSource.x[0] = x;
                            aNodeSource.x[1] = y;
                            aNodeSource.x[2] = z;
                            aNodeSource.vecVal.push_back(aVal);
                            vecSourceNodes.push_back(aNodeSource);
                        }
                    }
                    fclose(f);
                }
                else
                {
                    cout<<"SimulationManager::interpolatePrivate()->____reading file: "<<filePath.toStdString()<<"____"<<endl;
                    double x,y,z,aVal;
                    for(int k=0; feof(f)==0; k++)
                    {
                        if(4==fscanf(f,"%lf%lf%lf%lf",&x,&y,&z,&aVal))
                        {
                            Mapper3DClass::nodeSource aNodeSource = vecSourceNodes[k];
                            aNodeSource.vecVal.push_back(aVal);
                            vecSourceNodes[k] = aNodeSource;
                        }
                    }
                }
                c++;
            }
        }
        NbStep = listTimeS.size();
    }
        break;

    default:
    {
        //! --------------
        //! send progress
        //! --------------
        if(aProgressIndicator!=Q_NULLPTR)
        {
            cout<<"____send first event____"<<endl;
            QString msg("Start reading files");
            prgEvent = new QProgressEvent(QProgressEvent_Init,0,1,0,msg,QProgressEvent_None,0,9999,1);
            QApplication::postEvent(aProgressIndicator,prgEvent);
            QApplication::processEvents();
            Sleep(500);
        }

        //! -------------------------
        //! retrieve the source file
        //! -------------------------
        QString filePath = node->getPropertyValue<QString>("Source file");
        if(filePath=="") return;
        cout<<"SimulationManager::interpolatePrivate()->____reading file: "<<filePath.toStdString()<<"____"<<endl;

        //! --------------------------------------
        //! interpolation on a single source file
        //! a dummy list of times
        //! --------------------------------------
        listTimeS<<"9999";

        //! ------------------------
        //! update the progress bar
        //! -------------------------
        if(aProgressIndicator!=Q_NULLPTR)
        {
            QString msg = QString("Reading the source file");
            int NbBodies = int(vecLocs.size());
            prgEvent = new QProgressEvent(QProgressEvent_Init,0,NbBodies-1,0,msg,QProgressEvent_Init,0,100,0);
            QApplication::postEvent(aProgressIndicator,prgEvent);
            QApplication::processEvents();
            Sleep(2000);
        }

        //! ---------------------
        //! read the source file
        //! ---------------------
        FILE *f = fopen(filePath.toStdString().c_str(),"r");
        if(f!=NULL)
        {
            for(;feof(f)==0;)
            {
                Mapper3DClass::nodeSource aNodeSource;
                double x,y,z,v;
                if(4==fscanf(f,"%lf%lf%lf%lf",&x,&y,&z,&v))
                {
                    aNodeSource.x[0] = x;
                    aNodeSource.x[1] = y;
                    aNodeSource.x[2] = z;
                    aNodeSource.vecVal.push_back(v);
                    vecSourceNodes.push_back(aNodeSource);
                }
            }
            fclose(f);
        }
    }
        break;
    }

    //! -----------------------------
    //! results of the interpolation
    //! -----------------------------
    QList<QMap<int,std::pair<double,double>>> listMapMinMax;
    QList<QMap<GeometryTag,QList<QMap<int,double>>>> listMapOfRes;

    //! --------------------------
    //! start the mapping process
    //! --------------------------
    QMap<GeometryTag,QString> computationTimesMap;
    for(QVector<GeometryTag>::iterator it = vecLocs.begin();it!=vecLocs.end();++it)
    {
        //! -------------
        //! start chrono
        //! -------------
        start = std::chrono::system_clock::now();

        GeometryTag loc = *it;
        int bodyIndex = it->parentShapeNr;

        cout<<"SimulationManager::interpolatePrivate()-> start the mapping process on body Nb "<<bodyIndex<<"______"<<endl;

        const occHandle(MeshVS_DataSource) &theMeshVS_DataSource = mySimulationDataBase->ArrayOfMeshDS.value(bodyIndex);

        if(!theMeshVS_DataSource.IsNull())
        {
            cout<<"SimulationManager::interpolatePrivate()->____start mapper____"<<endl;

            //! ---------------------------
            //! mapping on the target mesh
            //! ---------------------------
            Mapper3DClass mapper(theMeshVS_DataSource);
            mapper.setProgressIndicator(aProgressIndicator);
            mapper.setSource(vecSourceNodes);

            //! --------------------------------------------------------------------------
            //! after removing case 0
            //! merge switch cases into the already prepared mapper::perform(int theAlgo)
            //! --------------------------------------------------------------------------
            //! ----------------
            //! select the algo
            //! ----------------
            switch(algo)
            {
            case 0:
            {
                //! --------------------------------
                //! interpolate using nearest point
                //! --------------------------------
                mapper.setRemap(remapFlag);
                if(remapFlag == true) mapper.setRemappingSteps(NbRemappingSteps);
                double pinball = node->getPropertyValue<double>("Pinball");
                mapper.performNearest(pinball);
            }
                break;

            case 1:
            {
                //! --------------------------------
                //! interpolate using nearest point
                //! --------------------------------
                mapper.setNbBuckets(NbBucketsX,NbBucketsY,NbBucketsZ);
                mapper.splitSourceIntoBuckets();
                mapper.setRemap(remapFlag);
                if(remapFlag == true) mapper.setRemappingSteps(NbRemappingSteps);
                double pinball = node->getPropertyValue<double>("Pinball");
                mapper.performNearestNeighboring(pinball);
            }
                break;

            case 2:
            {
                //! ----------------------------------
                //! interpolate using shape functions
                //! ----------------------------------
                mapper.setNbBuckets(NbBucketsX,NbBucketsY,NbBucketsZ);
                mapper.splitSourceIntoBuckets();
                mapper.setRemap(remapFlag);
                if(remapFlag == true) mapper.setRemappingSteps(NbRemappingSteps);
                mapper.performShapeFunctions();
            }
                break;
            }

            cout<<"SimulationManager::interpolatePrivate()->____mapping finished____"<<endl;

            for(int pos=0;pos<NbStep;pos++)
            {
                QMap<GeometryTag,QList<QMap<int,double>>> mapOfRes;
                QMap<int,std::pair<double,double>> mapMinMax;
                QList<QMap<int,double>> listOfRes;

                cout<<"SimulationManager::interpolatePrivate()->____retrieve list of map of result at step "<<pos<<"____"<<endl;

                //! retrieve the mapper interpolation result at time "pos"
                mapper.retrieveResMap(pos);

                //! record the min-max values for each body undergoing interpolation at time pos
                std::pair<double,double> aPair = mapper.getMinMax();

                //! retrieve the results at time pos
                listOfRes<<mapper.getResults();

                if(it==vecLocs.begin())
                {
                    cout<<"SimulationManager::interpolatePrivate()->inserting results on first body number = "<<bodyIndex<<endl;
                    mapMinMax.insert(bodyIndex,aPair);
                    mapOfRes.insert(loc,listOfRes);

                    //! insert all results type into the "time" list of results
                    listMapMinMax<<mapMinMax;
                    listMapOfRes<<mapOfRes;
                }
                else
                {
                    cout<<"SimulationManager::interpolatePrivate()->inserting results on next body number =  "<<bodyIndex<<endl;
                    mapMinMax=listMapMinMax[pos];
                    mapOfRes=listMapOfRes[pos];

                    mapMinMax.insert(bodyIndex,aPair);
                    mapOfRes.insert(loc,listOfRes);

                    //! --------------------------------------------
                    //! insert all results into the list of results
                    //! --------------------------------------------
                    listMapMinMax.replace(pos,mapMinMax);
                    listMapOfRes.replace(pos,mapOfRes);
                }
            }
            //cout<<"SimulationManager::interpolatePrivate()->____exit if cycle______"<<endl;
        }
        //cout<<"SimulationManager::interpolatePrivate()->____exit iteration over bodies______"<<endl;

        //! ------------
        //! stop chrono
        //! ------------
        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;
        computationTimesMap.insert(loc,QString("Interpolation time [s]: %1").arg(elapsed_seconds.count()));
    }

    //! -----------------------------
    //! close the progress indicator
    //! -----------------------------
    if(aProgressIndicator!=Q_NULLPTR)
    {
        prgEvent = new QProgressEvent(QProgressEvent_Reset,0,100,0,"",QProgressEvent_Reset,0,100,0);
        QApplication::postEvent(aProgressIndicator,prgEvent);
        QApplication::processEvents();
    }
    cout<<"SimulationManager::interpolatePrivate()->____mapping process finished______"<<endl;

    //! ---------------------
    //! retrieve the results
    //! ---------------------
    int t=0;
    for(QList<QString>::iterator itt = listTimeS.begin();itt!=listTimeS.end();++itt)
    {
        //! ----------------------------
        //! source time value by string
        //! ----------------------------
        QString timeS = *itt;

        //! -----------------------------------
        //! hide all the bodies and the meshes
        //! -----------------------------------
        TColStd_ListOfInteger l;
        for(int i=1; i<=mySimulationDataBase->bodyMap.size(); i++) l.Append(i);
        emit requestHideBody(l);
        emit requestHideMeshes();

        //! ---------------------------------------------------
        //! create a post object (colored mesh with color box)
        //! ---------------------------------------------------
        QString postObjectName = QString("Interpolation result on");

        //! -----------------------------------------
        //! prepare the name(s) of the postObject(s)
        //! -----------------------------------------
        mapOfMeshDataSources meshMap;
        for(QVector<GeometryTag>::const_iterator it = vecLocs.cbegin();it!=vecLocs.cend();++it)
        {
            GeometryTag loc = *it;
            int bodyIndex=loc.parentShapeNr;

            //! ------------------------------
            //! pile up the mesh data sources
            //! ------------------------------
            const occHandle(MeshVS_DataSource) &theMeshVS_DataSource = mySimulationDataBase->ArrayOfMeshDS.value(bodyIndex);
            meshMap.insert(loc, theMeshVS_DataSource);

            //! ---------------------------
            //! use the name of the bodies
            //! ---------------------------
            QString name = mySimulationDataBase->MapOfBodyNames.value(bodyIndex);
            postObjectName.append(" ").append(name).append("\n").append(computationTimesMap.value(loc));

            cout<<"Simulation manager::interpolatePrivate()->____the interpolation post object has name: \""<<postObjectName.toStdString()<<"\"____"<<endl;
        }

        //! ---------------------------------------------
        //! append the timestamp to the color box legend
        //! ---------------------------------------------
        postObjectName.append("\n").append(tools::timeStamp()).append("\n");

        QMap<GeometryTag,QList<QMap<int,double>>> mapOfRes = listMapOfRes.at(t);

        //! -----------------------
        //! create the post object
        //! -----------------------
        postObject aPostObject(mapOfRes,vecLocs,postObjectName);

        //! ------------
        //! 1-st column
        //! ------------
        int component = 0;

        //! -----------------------
        //! generica data continer
        //! -----------------------
        QVariant data;

        bool isAutoscale = true;
        aPostObject.buildMeshIO(meshMap,-1e80,1e80,10,isAutoscale,component);
        data.setValue(aPostObject);
        Property prop_postObject("Post object",data,Property::PropertyGroup_GraphicObjects);
        data.setValue(prop_postObject);

        SimulationNodeClass *nodeResult = nodeFactory::nodeFromScratch(SimulationNodeClass::nodeType_postObject);

        //! ----------------------------------
        //! add the post object as a property
        //! ----------------------------------
        nodeResult->addProperty(prop_postObject);

        //! -----------------------------------------------------
        //! add the timeTag copied from the parent item "Mapper"
        //! -----------------------------------------------------
        data.setValue(parentTimeTag);
        Property prop_timeTag("Time tag",data,Property::PropertyGroup_Identifier);
        nodeResult->addProperty(prop_timeTag);

        //! ----------------------------------------------------------
        //! add the "Source time" and "Analysis time" properties only
        //! for interpolation over time
        //! ----------------------------------------------------------
        switch(mode)
        {
        case 4:
        {
            data.setValue(timeS.toDouble());
            Property prop_sourceTime("Source time",data,Property::PropertyGroup_Definition);
            nodeResult->addProperty(prop_sourceTime);

            //! at the beginning the source time == analysis time
            data.setValue(timeS.toDouble());
            Property prop_analysisTime("Analysis time",data,Property::PropertyGroup_Definition);
            nodeResult->addProperty(prop_analysisTime);
        }
            break;

        default:
        {
            data.setValue(1.0);
            Property prop_analysisTime("Analysis time",data,Property::PropertyGroup_Definition);
            nodeResult->addProperty(prop_analysisTime);
        }
            break;
        }

        //! ------------------------------------------------------------------
        //! replace since "tags" have been already created by the nodeFactory
        //! ------------------------------------------------------------------
        data.setValue(vecLocs);
        Property prop_tags("Tags",data,Property::PropertyGroup_Scope);
        data.setValue(prop_tags);
        nodeResult->replaceProperty("Tags",prop_tags);

        //! ----------------
        //! create the item
        //! ----------------
        data.setValue(nodeResult);
        QExtendedStandardItem *itemResult = new QExtendedStandardItem();
        QString itemName;
        switch(mode)
        {
        case 4:
        {
            itemName = QString("Interpolation result T = %1").arg(timeS.toDouble());
            itemResult->setIcon(QIcon(":/icons/icon_clock.png"));
        }
            break;
        default:
        {
            itemName = QString("Interpolation result");
            itemResult->setIcon(QIcon(":/icons/icon_locator.png"));
        }
            break;
        }
        itemResult->setData(itemName,Qt::DisplayRole);
        itemResult->setData(data,Qt::UserRole);

        //! ----------------------------
        //! append to the "Mapper" item
        //! ----------------------------
        curItem->appendRow(itemResult);

        connect(nodeResult->getModel(),SIGNAL(itemChanged(QStandardItem*)),this,SLOT(handleItemChange(QStandardItem*)));

        t++;
    }
}

//! -----------------------------------------------------------------
//! function: retrieveCurrentItemResult
//! details:  retrieve a result in the form of a MeshVS_Mesh object
//!           from the current item
//! -----------------------------------------------------------------
bool SimulationManager::retrieveCurrentItemResult(postObject &aPostObject)
{
    ccout("SimulationManager::retrieveItemMesh()->____function called____");
    SimulationNodeClass* curNode = myTreeView->currentIndex().data(Qt::UserRole).value<SimulationNodeClass*>();
    QExtendedStandardItem* item = curNode->getPropertyItem("Post object");
    if(item!=Q_NULLPTR)
    {
        aPostObject = item->data(Qt::UserRole).value<Property>().getData().value<postObject>();
        return true;
    }
    return false;
}

//! -------------------------------------------------------------
//! function: retrieveAllResults
//! details:  retrieve the results of an analysis run and/or the
//!           result of an inerpolation
//! -------------------------------------------------------------
QList<postObject> SimulationManager::retrieveAllResults()
{
    cout<<"SimulationManager::retrieveAllResults()->____function called____"<<endl;
    QList<postObject> results;

    //! -----------------
    //! the current node
    //! -----------------
    QStandardItem *curItem = myModel->itemFromIndex(myTreeView->currentIndex());
    SimulationNodeClass *curNode = myTreeView->currentIndex().data(Qt::UserRole).value<SimulationNodeClass*>();

    if(curNode->isAnalysisResult()) return results;

    //! ----------------------------
    //! the current "Solution" item
    //! ----------------------------
    QStandardItem *curSolutionItem = curItem->parent();
    for(int i=1; i<curSolutionItem->rowCount(); i++)
    {
        QStandardItem *curResultItem = curSolutionItem->child(i,0);
        SimulationNodeClass *curResultNode = curResultItem->data(Qt::UserRole).value<SimulationNodeClass*>();

        QExtendedStandardItem *itemPostObject = curResultNode->getPropertyItem("Post object");
        if(itemPostObject==Q_NULLPTR) continue;

        const postObject &aPostObject = itemPostObject->data(Qt::UserRole).value<Property>().getData().value<postObject>();
        results<<aPostObject;
    }

    return results;
}

//! ----------------------------------
//! function: handleGlobalMeshControl
//! details:  make all meshes invalid
//! ----------------------------------
void SimulationManager::handleGlobalMeshControlChange()
{
    std::vector<int> l;
    for(QMap<int,TopoDS_Shape>::iterator it = mySimulationDataBase->bodyMap.begin(); it!=mySimulationDataBase->bodyMap.end(); it++)
    {
        int bodyIndex = it.key();
        l.push_back(bodyIndex);
    }
    this->requestMeshInvalidate(l);
}

//! ---------------------------------------------------------------------
//! function: translateOpenFoamScalarData
//! details:  start the translation process - suitable for scalar values
//! ---------------------------------------------------------------------
#include "openfoamcontroller.h"
bool SimulationManager::translateOpenFoamScalarData()
{
    cout<<"SimulationManager::translateOpenFoamScalarData->____function called____"<<endl;

    //! ---------------------------------------------
    //! retrieve the settings from the detail viewer
    //! ---------------------------------------------
    SimulationNodeClass *theCurNode = myTreeView->currentIndex().data(Qt::UserRole).value<SimulationNodeClass*>();
    QExtendedStandardItem *theItem;

    //! -------------------------
    //! get the source directory
    //! -------------------------
    theItem = theCurNode->getPropertyItem("Source directory");
    const QString &sourceDirectory = theItem->data(Qt::UserRole).value<Property>().getData().toString();
    cout<<"SimulationManager::translateOpenFoamScalarData->____source directory: "<<sourceDirectory.toStdString()<<"____"<<endl;

    //! -------------------------
    //! get the target directory
    //! -------------------------
    theItem = theCurNode->getPropertyItem("Target directory");
    const QString &targetDirectory = theItem->data(Qt::UserRole).value<Property>().getData().toString();

    //! ----------------------------------------------------
    //! get the file mode 0 - one file 1 - split into files
    //! ----------------------------------------------------
    theItem = theCurNode->getPropertyItem("Split data");
    int fileMode = theItem->data(Qt::UserRole).value<Property>().getData().toInt();

    if(sourceDirectory.isEmpty() || targetDirectory.isEmpty()) return false;

    //! ---------------------------
    //! get the progress indicator
    //! ---------------------------
    QProgressIndicator *aProgressIndicator = static_cast<QProgressIndicator*>(tools::getWidgetByName("progressIndicator"));

    //! ---------------
    //! another thread
    //! ---------------
    openFoamController *anOpenFoamController = new openFoamController(sourceDirectory,targetDirectory,fileMode,aProgressIndicator,this);

#ifdef COSTAMP_VERSION
    anOpenFoamController->setTimeFolders(tSbList);
#endif

    //! --------------------------------------------------------------------------
    //! start the thread - this will also lock the items within the detail viewer
    //! --------------------------------------------------------------------------
    anOpenFoamController->operate(theCurNode);


    return true;
}

//! ---------------------------------------
//! function: callPostEngineEvaluateResult
//! details:  evaluate a single result
//! ---------------------------------------
void SimulationManager::callPostEngineEvaluateResult()
{
    cout<<"SimulationManager::callPostEngineEvaluateResult()->____function called____"<<endl;

    //! ----------------------
    //! the current item/node
    //! ----------------------
    QModelIndex index = myTreeView->currentIndex();
    QExtendedStandardItem *curItem = static_cast<QExtendedStandardItem*>(myModel->itemFromIndex(index));
    SimulationNodeClass *curNode = curItem->data(Qt::UserRole).value<SimulationNodeClass*>();

    if(curNode->isSolutionInformation())
    {
        cout<<"SimulationManager::callPostEngineEvaluateResult()->____solutionInformation Item____"<<endl;

        //! -------------------------------------------
        //! the current node is "Solution information"
        //! jump over the first child (i=0)
        //! -------------------------------------------
        QExtendedStandardItem *itemSolution = static_cast<QExtendedStandardItem*>(curItem->parent());
        for(int i=1; i<itemSolution->rowCount(); i++)
        {
            QExtendedStandardItem *postProcessingItem = static_cast<QExtendedStandardItem*>(itemSolution->child(i,0));
            this->callPostEngineEvaluateResult_private(postProcessingItem,false);
        }
    }
    if(curNode->isSolution())
    {
        cout<<"SimulationManager::callPostEngineEvaluateResult()->____solution Item____"<<endl;
        //! --------------------------------
        //! the current node is "Solution"
        //! jump over the first child (i=0)
        //! --------------------------------
        QExtendedStandardItem *itemSolution = curItem;
        for(int i=1; i<itemSolution->rowCount(); i++)
        {
            QExtendedStandardItem *postProcessingItem = static_cast<QExtendedStandardItem*>(itemSolution->child(i,0));
            this->callPostEngineEvaluateResult_private(postProcessingItem, false);
        }
    }
    if(curNode->isAnalysisResult())
    {
        //! ------------------------------------------
        //! the current item is a postprocessing item
        //! ------------------------------------------
        bool immediatelyDisplay = true;
        cout<<"SimulationManager::callPostEngineEvaluateResult()->____postProcessing Item____"<<endl;
        this->callPostEngineEvaluateResult_private(curItem, immediatelyDisplay);
    }

    //! parse the item
    //parser::PostProcessingItem(treeItem);
}

//! ----------------------------------------------------------------------
//! function: callPostEngineEvaluateResult_private
//! details:  retrieve scope and time info, then call the "postEngine"
//!           which creates the nodes. If "immediately display" is "true"
//!           the results are immediately displayed on the screen
//! ----------------------------------------------------------------------
void SimulationManager::callPostEngineEvaluateResult_private(QStandardItem *curItem, bool immediatelyDisplay)
{
    cout<<"SimulationManager::callPostEngineEvaluateResult_private()->____function called____"<<endl;

    //! --------------
    //! the node type
    //! --------------
    SimulationNodeClass *curNode = curItem->data(Qt::UserRole).value<SimulationNodeClass*>();
    SimulationNodeClass::nodeType type = curNode->getType();

    //! -----------------------------------------------------------
    //! single post processing items coming from an interpolation
    //! should be excluded from this
    //! -----------------------------------------------------------
    if(type==SimulationNodeClass::nodeType_postObject) return;

    //! -------------------------------------
    //! retrieve the location of the results
    //! -------------------------------------
    QStandardItem *itemSolution = curItem->parent();
    SimulationNodeClass *nodeSolution = itemSolution->data(Qt::UserRole).value<SimulationNodeClass*>();
    QString projectFilesDir = nodeSolution->getPropertyValue<QString>("Project files dir");
    QString timeTag = nodeSolution->getPropertyValue<QString>("Parent time tag");
    QString resultsFilePath = projectFilesDir+"/SolutionData_"+timeTag+"/input.frd";

    myPostEngine->setResultsFile(resultsFilePath);
    QStandardItem *itemSolutionInformation = itemSolution->child(0,0);
    SimulationNodeClass *nodeSolutionInformation = itemSolutionInformation->data(Qt::UserRole).value<SimulationNodeClass*>();
    QMap<double,QVector<int>> dtm = nodeSolutionInformation->getPropertyValue<QMap<double,QVector<int>>>("Discrete time map");
    myPostEngine->setDiscreteTimeMap(dtm);

    //! ----------------------
    //! retrieve the location
    //! ----------------------
    QVector<GeometryTag> vecLoc = curNode->getPropertyValue<QVector<GeometryTag>>("Tags");

    //! ----------------------------------------------------------------
    //! check if a mesh for each location exists
    //! This avoids application crash when calling "Evaluate result(s)"
    //! and the seleted geometry has not a mesh yet
    //! ----------------------------------------------------------------
    bool isMeshOK = true;
    for(QVector<GeometryTag>::iterator it = vecLoc.begin(); it!=vecLoc.end(); ++it)
    {
        const GeometryTag &loc = *it;
        if(loc.isParent)
        {
            if(mySimulationDataBase->ArrayOfMeshDS.value(loc.parentShapeNr).IsNull())
            {
                isMeshOK=false;
                cout<<"SimulationManager::callPostEngineEvaluateResult_private()->____mesh is null____"<<endl;
                break;
            }
        }
        else
        {
            int bodyIndex = loc.parentShapeNr;
            int subShapeIndex = loc.subTopNr;
            TopAbs_ShapeEnum type = loc.subShapeType;
            if(type == TopAbs_FACE)
            {
                if(mySimulationDataBase->ArrayOfMeshDSOnFaces.getValue(bodyIndex,subShapeIndex).IsNull())
                {
                    isMeshOK = false;
                    cout<<"____mesh is null____"<<endl;
                    break;
                }
            }
            if(type == TopAbs_EDGE)
            {
                if(mySimulationDataBase->ArrayOfMeshDSOnEdges.getValue(bodyIndex,subShapeIndex).IsNull())
                {
                    cout<<"SimulationManager::callPostEngineEvaluateResult_private()->____mesh is null____"<<endl;
                    isMeshOK = false;
                    break;
                }
            }
            if(type == TopAbs_VERTEX)
            {
                //! to do...
                cout<<"SimulationManager::callPostEngineEvaluateResult_private()->____mesh is null____"<<endl;
                isMeshOK = false;
                break;
            }
        }
    }
    if(isMeshOK==false) return;

    switch (type)
    {
    case SimulationNodeClass::nodeType_solutionStructuralFatigueTool:
    {
        QList<double> timeList;
        int component = curNode->getPropertyValue<int>("Component");
        int NbCycles = curNode->getPropertyValue<int>("Number of cycles");
        //cout<<"____number of cycles: "<<NbCycles<<"____"<<endl;
        //exit(1);
        postObject aPostObject;

        //! --------------------------------------------
        //! a results is already present into the item
        //! --------------------------------------------
        if(curNode->getPropertyItem("Post object")!=Q_NULLPTR)
        {
            aPostObject = curNode->getPropertyItem("Post object")->data(Qt::UserRole).value<Property>().getData().value<postObject>();
            aPostObject.update(static_cast<meshDataBase*>(mySimulationDataBase), component);
        }
        else
        {
            // left here for documentation
            //QMap<double,QVector<int>> dTm = myPostEngine->getDTM();

            //! ---------------------------------------
            //! retrieve the solution information item
            //! ---------------------------------------
            QStandardItem *itemSolutionInformation = curItem->parent()->child(0,0);
            SimulationNodeClass *nodeSolutionInformation = itemSolutionInformation->data(Qt::UserRole).value<SimulationNodeClass*>();
            QMap<double,QVector<int>> dTm = nodeSolutionInformation->getPropertyValue<QMap<double,QVector<int>>>("Discrete time map");

            //! ------------------------------------------------------------------
            //! execute only if the discrete time map is not empty
            //! (that is the .frd file exixts and it has been scanned)
            //! This avoids application crash when requiring "Evaluate result(s)"
            //! and no result exixts
            //! ------------------------------------------------------------------
            if(dTm.isEmpty())
            {
                cout<<"----------------------------------------------------------------------------------------------------------------------"<<endl;
                cout<<"SimulationManager::callPostEngineEvaluateResult_private()->____cannot evaluate results: the discrete time is empty____"<<endl;
                cout<<"----------------------------------------------------------------------------------------------------------------------"<<endl;
                QMessageBox::warning(this,APPNAME,"No result available", QMessageBox::Ok);
                return;
            }

            //! ---------------------
            //! set the fatigue algo
            //! ---------------------
            int fatigueAlgo = curNode->getPropertyValue<int>("Fatigue algo");
            myPostEngine->setFatigueModel(fatigueAlgo);

            CustomTableModel *tabData =  this->getAnalysisSettingsNodeFromCurrentItem()->getTabularDataModel();

            QStandardItem *theGeometryRoot=this->getTreeItem(SimulationNodeClass::nodeType_geometry);
            QMap<int,int> materialBodyMap;
            for(int k=0; k<theGeometryRoot->rowCount();k++)
            {
                QStandardItem *theGeometryItem = theGeometryRoot->child(k,0);
                SimulationNodeClass *theCurNode = theGeometryItem->data(Qt::UserRole).value<SimulationNodeClass*>();
                Property::SuppressionStatus theNodeSS = theCurNode->getPropertyValue<Property::SuppressionStatus>("Suppressed");
                if(theNodeSS==Property::SuppressionStatus_Active)
                {
                    int loc = theCurNode->getPropertyValue<int>("Map index");
                    int matNumber = theCurNode->getPropertyValue<int>("Assignment");
                    materialBodyMap.insert(loc,matNumber);
                }
            }
            for(int i=0;i<tabData->rowCount();i++) timeList<<tabData->dataRC(i,1).toDouble();

            //! -----------------------------------------------------------------------------
            //! create the postObject
            //! the post object retrieves the mesh data sources from the simulation database
            //! and internally builds its own interactive mesh objects
            //! -----------------------------------------------------------------------------

            aPostObject = myPostEngine->evaluateFatigueResults(component,vecLoc,timeList,materialBodyMap,NbCycles);
            aPostObject.update(static_cast<meshDataBase*>(mySimulationDataBase), component);
        }
        //! ------------------------------
        //! set the title of the colorbox
        //! ------------------------------
        QVariant data;
        data.setValue(aPostObject);

        Property prop_postObject("Post object",data,Property::PropertyGroup_GraphicObjects);
        curNode->removeProperty("Post object");
        curNode->addProperty(prop_postObject);

        //! ---------------------------------------------------------------------
        //! color box controls: synchronize the post object min, man, scale type
        //! with the color box properties
        //! ---------------------------------------------------------------------
        if(curNode->getPropertyItem("Scale type")==NULL) cerr<<"____NULL property____"<<endl;

        if(curNode->getPropertyValue<int>("Scale type") == 1)
        {
            disconnect(curNode->getModel(),SIGNAL(itemChanged(QStandardItem*)),this,SLOT(handleItemChange(QStandardItem*)));
            double minValue = aPostObject.getMin();
            double maxValue = aPostObject.getMax();
            int NbLevels = aPostObject.getNbLevels();

            data.setValue(minValue);
            Property prop_min("Min",data,Property::PropertyGroup_ColorBox);
            data.setValue(maxValue);
            Property prop_max("Max",data,Property::PropertyGroup_ColorBox);
            data.setValue(NbLevels);
            Property prop_intervals("# intervals",data,Property::PropertyGroup_ColorBox);
            curNode->replaceProperty("Min",prop_min);
            curNode->replaceProperty("Max",prop_max);
            curNode->replaceProperty("# intervals",prop_intervals);
            connect(curNode->getModel(),SIGNAL(itemChanged(QStandardItem*)),this,SLOT(handleItemChange(QStandardItem*)));
        }

        if(immediatelyDisplay == true)
        {
            emit requestHideMeshes();
            emit requestSetWorkingMode(3);
            emit requestDisplayResult(aPostObject);
        }
    }
        break;
    default:
    {
        int component = curNode->getPropertyValue<int>("Type ");
        int mode = curNode->getPropertyValue<int>("Mode number");
        //! ---------------------------------------------
        //! a results is already present into the item
        //! build the mesh object from the internal data
        //! ---------------------------------------------
        postObject aPostObject;
        if(curNode->getPropertyItem("Post object")!=Q_NULLPTR)
        {
            aPostObject = curNode->getPropertyItem("Post object")->data(Qt::UserRole).value<Property>().getData().value<postObject>();
            aPostObject.update(static_cast<meshDataBase*>(mySimulationDataBase), component);
        }
        else
        {
            //! ------------------------------------------------------------
            //! retrieve the time info:
            //! a result is always retrieved using the pair (step, substep)
            //! ------------------------------------------------------------
            int subStepNb, stepNb;

            //! if "By" is "Time" - "Analysis time" is read from the GUI
            double analysisTime;

            //! if "By" is "Set" - "Set number" is read from the GUI
            int setNumber;
            int rmode = curNode->getPropertyValue<int>("By");

            // left here for documentation
            //QMap<double,QVector<int>> dTm = myPostEngine->getDTM();

            //! ---------------------------------------
            //! retrieve the solution information item
            //! ---------------------------------------
            QStandardItem *itemSolutionInformation = curItem->parent()->child(0,0);
            SimulationNodeClass *nodeSolutionInformation = itemSolutionInformation->data(Qt::UserRole).value<SimulationNodeClass*>();
            QMap<double,QVector<int>> dTm = nodeSolutionInformation->getPropertyValue<QMap<double,QVector<int>>>("Discrete time map");


            for(QMap<double,QVector<int>>::iterator it = dTm.begin(); it!=dTm.end(); ++it)
            {
                cout<<"@-------------------------"<<endl;
                cout<<"@ time: "<<it.key()<<endl;
                for(int i=0; i<it.value().size(); i++)
                {
                    cout<<"@ set: "<<it.value().at(0)<<endl;
                    cout<<"@ step: "<<it.value().at(1)<<endl;
                    cout<<"@ substep: "<<it.value().at(2)<<endl;
                }
                cout<<"@-------------------------"<<endl;
            }

            //! ------------------------------------------------------------------
            //! execute only if the discrete time map is not empty
            //! (that is the .frd file exixts and it has been scanned)
            //! This avoids application crash when requiring "Evaluate result(s)"
            //! and no result exixts
            //! ------------------------------------------------------------------
            if(dTm.isEmpty())
            {
                cout<<"----------------------------------------------------------------------------------------------------------------------"<<endl;
                cout<<"SimulationManager::callPostEngineEvaluateResult_private()->____cannot evaluate results: the discrete time is empty____"<<endl;
                cout<<"----------------------------------------------------------------------------------------------------------------------"<<endl;
                QMessageBox::warning(this,APPNAME,"No result available", QMessageBox::Ok);
                return;
            }
            switch(rmode)
            {
            case 0:
            {
                //! --------------------------------------------
                //! from "Analysis time" to ("Step", "Substep")
                //! --------------------------------------------
                cout<<"SimulationManager::callPostEngineEvaluateResult_private()->____retriving data using \"Display time\": "<<endl;
                analysisTime = curNode->getPropertyValue<double>("Display time");
                postTools::getStepSubStepByTimeDTM(dTm,analysisTime, stepNb, subStepNb);
                cout<<"\"StepNb\" = "<<stepNb<<endl;
                cout<<"\"SubStepNb\" = "<<subStepNb<<endl;
            }
                break;

            case 1:
            {
                //! -----------------------------------------
                //! from "Set" number to ("Step", "Substep")
                //! "Analysis time" is returned also
                //! -----------------------------------------
                cout<<"SimulationManager::callPostEngineEvaluateResult_private()->____retriving data using \"Set number\": "<<endl;
                setNumber = curNode->getPropertyValue<int>("Set number");
                postTools::getStepSubStepBySetDTM(dTm, setNumber, analysisTime, stepNb, subStepNb);
            }
                break;
            }

            //! -------------------------------------------------------
            //! the required "Set"/"Display time" are always converted
            //! into the pair ("Step number", "Substep number")
            //! -------------------------------------------------------
            QString keyName;
            switch(type)
            {
            case SimulationNodeClass::nodeType_solutionThermalTemperature: keyName ="NDTEMP"; break;
            case SimulationNodeClass::nodeType_solutionThermalFlux: keyName ="FLUX"; break;
            case SimulationNodeClass::nodeType_solutionStructuralNodalDisplacement: keyName = "DISP"; break;
            case SimulationNodeClass::nodeType_solutionStructuralStress: keyName ="STRESS"; break;
            case SimulationNodeClass::nodeType_solutionStructuralTotalStrain: keyName ="TOSTRAIN"; break;
            case SimulationNodeClass::nodeType_solutionStructuralMechanicalStrain: keyName ="MESTRAIN"; break;
            case SimulationNodeClass::nodeType_solutionStructuralEquivalentPlasticStrain: keyName ="PE"; break;
            case SimulationNodeClass::nodeType_solutionStructuralNodalForces: keyName ="FORC"; break;
            case SimulationNodeClass::nodeType_solutionStructuralReactionForce: keyName = "FORC"; break;
            case SimulationNodeClass::nodeType_solutionStructuralTemperature: keyName ="NDTEMP"; break;
            case SimulationNodeClass::nodeType_solutionStructuralContact: keyName = "CONTACT"; break;
            }

            //! -----------------------------------------------
            //! launch the post engine on the <keyName> result
            //! -----------------------------------------------
            cout<<"@____launching the post engine____"<<endl;
            cout<<"@____keyName: "<<keyName.toStdString()<<"____"<<endl;
            cout<<"@____component: "<<component<<"____"<<endl;
            cout<<"@____step: "<<stepNb<<"____"<<endl;
            cout<<"@____sub step: "<<subStepNb<<"____\n"<<endl;

            //! -----------------------------------------------------------------------------
            //! create the postObject
            //! the post object retrieves the mesh data sources from the simulation database
            //! and internally builds its own interactive mesh objects
            //! -----------------------------------------------------------------------------
            aPostObject = myPostEngine->buildPostObject(keyName,component,subStepNb,stepNb,mode,vecLoc);
            aPostObject.update(static_cast<meshDataBase*>(mySimulationDataBase), component);
        }

        QVariant data;
        data.setValue(aPostObject);
        Property prop_postObject("Post object",data,Property::PropertyGroup_GraphicObjects);

        curNode->removeProperty("Post object");
        curNode->addProperty(prop_postObject);

        //! ---------------------------------------------------------------------
        //! color box controls: synchronize the post object min, man, scale type
        //! with the color box properties
        //! ---------------------------------------------------------------------
        if(curNode->getPropertyItem("Scale type")==Q_NULLPTR) cerr<<"____NULL property____"<<endl;

        if(curNode->getPropertyValue<int>("Scale type") == 1)
        {
            disconnect(curNode->getModel(),SIGNAL(itemChanged(QStandardItem*)),this,SLOT(handleItemChange(QStandardItem*)));
            double minValue = aPostObject.getMin();
            double maxValue = aPostObject.getMax();
            int NbLevels = aPostObject.getNbLevels();

            data.setValue(minValue);
            Property prop_min("Min",data,Property::PropertyGroup_ColorBox);
            data.setValue(maxValue);
            Property prop_max("Max",data,Property::PropertyGroup_ColorBox);
            data.setValue(NbLevels);
            Property prop_intervals("# intervals",data,Property::PropertyGroup_ColorBox);
            curNode->replaceProperty("Min",prop_min);
            curNode->replaceProperty("Max",prop_max);
            curNode->replaceProperty("# intervals",prop_intervals);
            connect(curNode->getModel(),SIGNAL(itemChanged(QStandardItem*)),this,SLOT(handleItemChange(QStandardItem*)));
        }
    }
        break;
    }
}

//! ---------------------------------------------------------------
//! function: evaluateAllResults
//! details:  evaluate all the results under the "Solution" branch
//! ---------------------------------------------------------------
void SimulationManager::evaluateAllResults()
{
    cout<<"SimulationManager::evaluateAllResults()->____function called____"<<endl;

    //! -------------------------------------
    //! retrieve the current "Solution" item
    //! -------------------------------------
    QStandardItem *curItemSolution;
    QStandardItem *curItem = myModel->itemFromIndex(myTreeView->currentIndex());
    SimulationNodeClass *curNode = curItem->data(Qt::UserRole).value<SimulationNodeClass*>();
    if(curNode->isSolution()) curItemSolution = curItem;
    if(curNode->isAnalysisRoot()) curItemSolution = curItem->child(curItem->rowCount()-1,0);

    if(curItemSolution->hasChildren()==false)
    {
        QMessageBox::information(this,APPNAME,"No result has been defined yet",QMessageBox::Ok);
        return;
    }

    for(int row = 1; row<curItemSolution->rowCount(); row++)
    {
        QStandardItem *itemResult = curItemSolution->child(row,0);
        SimulationNodeClass *nodeResult = itemResult->data(Qt::UserRole).value<SimulationNodeClass*>();
        this->callPostEngineEvaluateResult_private(itemResult,false);
    }
}

//! ----------------------
//! function: eventFilter
//! details:
//! ----------------------
bool SimulationManager::eventFilter(QObject *object, QEvent *event)
{
    if(event->type()==QSimulationStatusEvent::type() || event->type()==QCCXSolverMessageEvent::type())
    {
        //! -------------
        //! generic data
        //! -------------
        QVariant data;

        //! -----------------------------------------
        //! retrieve the "Solution information" item
        //! -----------------------------------------
        QStandardItem *itemSolution = myCurrentRunningAnalysis->child(myCurrentRunningAnalysis->rowCount()-1,0);
        SimulationNodeClass *nodeSolution = itemSolution->data(Qt::UserRole).value<SimulationNodeClass*>();
        QStandardItem *itemSolutionInformation = myCurrentRunningAnalysis->child(myCurrentRunningAnalysis->rowCount()-1,0)->child(0,0);
        SimulationNodeClass *nodeSolutionInformation = itemSolutionInformation->data(Qt::UserRole).value<SimulationNodeClass*>();

        //if(event->type()==QSimulationStatusEvent::type())
        //{
        //    cout<<"SimulationManager::eventFilter()->____simulation status event received____"<<endl;
        //}
        //else if(event->type()==QCCXSolverMessageEvent::type())
        if(event->type()==QCCXSolverMessageEvent::type())
        {
             cout<<"SimulationManager::eventFilter()->____CCX solver message event received____"<<endl;
             QCCXSolverMessageEvent *e = static_cast<QCCXSolverMessageEvent*>(event);
             CCXSolverMessage ccxmsg = e->getMessage();
             data.setValue(ccxmsg);
             nodeSolutionInformation->replaceProperty("Solver output",Property("Solver output",data,Property::PropertyGroup_Hidden));

             //! -------------------
             //! read the .sta file
             //! -------------------
             QString projectFilesDir = nodeSolution->getPropertyValue<QString>("Project files dir");
             QString timeTag = nodeSolution->getPropertyValue<QString>("Parent time tag");
             QString stafile = projectFilesDir+"/SolutionData_"+timeTag+"/input.sta";
             QFile f(stafile);
             if(f.exists())
             {
                 cout<<"____.STA FILE FOUND: \""<<stafile.toStdString()<<"\"____"<<endl;
                 QMap<double,QVector<int>> timeinfo;
                 bool isDone = CCXTools::readsta(stafile,timeinfo);
                 if(isDone)
                 {
                     data.setValue(timeinfo);
                     cout<<timeinfo.firstKey()<<endl;
                     nodeSolutionInformation->replaceProperty("Discrete time map",Property("Discrete time map",data,Property::PropertyGroup_Hidden));
                 }
                 else
                 {
                     QVector<int> init;
                     double ini=0.0;
                     init<<0.0<<0.0<<0.0;
                     timeinfo.insert(ini,init);
                     data.setValue(timeinfo);
                     nodeSolutionInformation->replaceProperty("Discrete time map",Property("Discrete time map",data,Property::PropertyGroup_Hidden));
                 }
             }
        }

        //! --------------------------
        //! the event is stopped here
        //! --------------------------
        return true;
    }    
    else
    {
        return QObject::eventFilter(object, event);
    }
}

//! -------------------------
//! function: getCurrentNode
//! details:  helper
//! -------------------------
SimulationNodeClass* SimulationManager::getCurrentNode()
{
    return (myModel->itemFromIndex(myTreeView->currentIndex())->data(Qt::UserRole).value<SimulationNodeClass*>());
}

//! -------------------------
//! function: fromTagToShape
//! details:  helper
//! --------------------------
TopoDS_Shape SimulationManager::fromTagToShape(const GeometryTag &aTag)
{
    int parentShapeIndex = aTag.parentShapeNr;
    int subShapeIndex = aTag.subTopNr;
    TopAbs_ShapeEnum subShapeType = aTag.subShapeType;
    TopoDS_Shape aShape;
    if(aTag.isParent) aShape = mySimulationDataBase->bodyMap.value(parentShapeIndex);
    else
    {
        switch(subShapeType)
        {
        case TopAbs_FACE: aShape = mySimulationDataBase->MapOfBodyTopologyMap.value(parentShapeIndex).faceMap.FindKey(subShapeIndex); break;
        case TopAbs_EDGE: aShape = mySimulationDataBase->MapOfBodyTopologyMap.value(parentShapeIndex).edgeMap.FindKey(subShapeIndex); break;
        case TopAbs_VERTEX: aShape = mySimulationDataBase->MapOfBodyTopologyMap.value(parentShapeIndex).vertexMap.FindKey(subShapeIndex); break;            break;
        }
    }
    return aShape;
}

//! -------------------------
//! function: fromTagToShape
//! details:  helper
//! --------------------------
TopTools_ListOfShape SimulationManager::fromTagToShape(const QVector<GeometryTag> &vecLoc)
{
    TopTools_ListOfShape lshapes;
    for(QVector<GeometryTag>::const_iterator it = vecLoc.cbegin(); it!=vecLoc.cend(); ++it)
    {
        GeometryTag loc = *it;
        int parentShapeIndex = loc.parentShapeNr;
        int subShapeIndex = loc.subTopNr;
        TopAbs_ShapeEnum subShapeType = loc.subShapeType;
        if(loc.isParent)
        {
            lshapes.Append(mySimulationDataBase->bodyMap.value(parentShapeIndex));
        }
        else
        {
            switch(subShapeType)
            {
            case TopAbs_FACE:
                lshapes.Append(mySimulationDataBase->MapOfBodyTopologyMap.value(parentShapeIndex).faceMap.FindKey(subShapeIndex));
                break;
            case TopAbs_EDGE:
                lshapes.Append(mySimulationDataBase->MapOfBodyTopologyMap.value(parentShapeIndex).edgeMap.FindKey(subShapeIndex));
                break;
            case TopAbs_VERTEX:
                lshapes.Append(mySimulationDataBase->MapOfBodyTopologyMap.value(parentShapeIndex).vertexMap.FindKey(subShapeIndex));
                break;
            }
        }
    }
    return lshapes;
}

//! -----------------------------
//! function: retrieveSolverInfo
//! details:
//! -----------------------------
void SimulationManager::retrieveSolverInfo()
{
    static int msTime;
    cout<<"SimulationManager::retrieveSolverInfo()->____time: "<<++msTime<<"____"<<endl;

    if(myCurrentRunningAnalysis==Q_NULLPTR) return;

    int analysisType = -1;
    SimulationNodeClass::nodeType analysisNodeType = myCurrentRunningAnalysis->data(Qt::UserRole).value<SimulationNodeClass*>()->getType();
    switch(analysisNodeType)
    {
    case SimulationNodeClass::nodeType_structuralAnalysis: analysisType = 0; break;
    case SimulationNodeClass::nodeType_thermalAnalysis: analysisType = 1; break;
    case SimulationNodeClass::nodeType_combinedAnalysis: analysisType = 2; break;
    }

    //! ------------------------------------------------------
    //! "Time tag" identificator of the current analysis root
    //! ------------------------------------------------------
    QString timeTag = myCurrentRunningAnalysis->data(Qt::UserRole).value<SimulationNodeClass*>()->getPropertyValue<QString>("Time tag");

    QStandardItem *itemSolution = myCurrentRunningAnalysis->child(myCurrentRunningAnalysis->rowCount()-1);
    SimulationNodeClass *nodeSolution = itemSolution->data(Qt::UserRole).value<SimulationNodeClass*>();
    QString projectFilesDir = nodeSolution->getPropertyValue<QString>("Project files dir");
    QString targetFileName = projectFilesDir+"/SolutionData_"+timeTag+"/SolverOutput.txt";
    QString sourceFileName = projectFilesDir+"/SolutionData_"+timeTag+"/RawSolverOutput_copy.txt";
    QString lockFileName = projectFilesDir+"/SolutionData_"+timeTag+"/lock";

    QStandardItem *itemSolutionInformation = itemSolution->child(0,0);
    SimulationNodeClass *nodeSolutionInformation = itemSolutionInformation->data(Qt::UserRole).value<SimulationNodeClass*>();

    //! ---------------------------------
    //! retrieve the final analysis time
    //! ---------------------------------
    //QStandardItem *itemAnalysisSettings = myCurrentRunningAnalysis->child(0,0);
    //SimulationNodeClass *nodeAnalysisSettings = itemAnalysisSettings->data(Qt::UserRole).value<SimulationNodeClass*>();
    //CustomTableModel *tabData = nodeAnalysisSettings->getTabularDataModel();
    //double endTime = tabData->dataRC(tabData->rowCount()-1,1).toDouble();

    //! --------------------------------------------------------------
    //! check if the source file "RawSolverOutput_copy.txt" is in use.
    //! If a lock file exists, the file is in use; if not the file
    //! can be opened for being processed.
    //! Before processing "RawSolverOutput_copy.txt" deny the access,
    //! by creating a new lock file. Remove the lock file at the end
    //! of the processing.
    //! --------------------------------------------------------------
    if(!QFile::exists(lockFileName.toStdString().c_str()))
    {
        FILE *lockFile = fopen(lockFileName.toStdString().c_str(),"w");
        cout<<"SimulationManager::retrieveSolverInfo()->____lock file created____"<<endl;

        QList<solutionInfo> solInfoList;
        bool simulationError;
        runTerminationData rtd = CCXconsoleToFile::perform(targetFileName,sourceFileName,analysisType,solInfoList,simulationError);

        cout<<"@ -----------------------------------------"<<endl;
        cout<<"@ last av step: "<<rtd.lastAvailableStep<<endl;
        cout<<"@ last av substep: "<<rtd.lastAvailableSubStep<<endl;
        cout<<"@ last av time: "<<rtd.lastAvailableTime<<endl;
        cout<<"@ -----------------------------------------"<<endl;

        //! ------------------------
        //! update the progress bar
        //! ------------------------
        QProgressIndicator *aProgressIndicator = static_cast<QProgressIndicator*>(tools::getWidgetByName("progressIndicator"));
        QProgressEvent *e = new QProgressEvent(QProgressEvent_Update,-1,-1,100.0*rtd.lastAvailableTime,"Running CCX",QProgressEvent_None,-1,-1,-1,"Running CCX");
        QApplication::postEvent(aProgressIndicator,e);
        QApplication::processEvents();

        //! ---------------------------------------
        //! put the convergence data into the item
        //! ---------------------------------------
        QVariant data;
        data.setValue(solInfoList);
        nodeSolutionInformation->replaceProperty("Convergence data",Property("Convergence data",data,Property::PropertyGroup_Hidden));

        //! ---------------------------------------
        //! request convergence data viewer update
        //! ---------------------------------------
        emit requestUpdateConvergenceViewer(solInfoList);
        fclose(lockFile);
        bool lockRemoved = QFile::remove(lockFileName);
        if(lockRemoved) cout<<"SimulationManager::retrieveSolverInfo()->____lock file deleted____"<<endl;
    }
    else
    {
        cout<<"SimulationManager::retrieveSolverInfo()->____already opened: cannot retrieve solution information____"<<endl;
    }
}

//! ----------------------
//! function: updateTimer
//! details:
//! ----------------------
void SimulationManager::updateTimer()
{
    cout<<"SimulationManager::updateTimer()->____function called____"<<endl;

    //! ---------------------------------
    //! minimum update period in seconds
    //! ---------------------------------
    const double min_update_period = 1.5;

    //! -----------------------------------------------
    //! retrieve the current solution information item
    //! -----------------------------------------------
    SimulationNodeClass *nodeSolutionInformation = myTreeView->currentIndex().data(Qt::UserRole).value<SimulationNodeClass*>();
    if(nodeSolutionInformation->isSolutionInformation()==false) return;

    double updateInterval = nodeSolutionInformation->getPropertyValue<double>("Update interval");
    if(updateInterval<=min_update_period) updateInterval = min_update_period;
    cout<<"SimulationManager::updateTimer()->____new timer interval: "<<int(updateInterval*1000)<<" ms____"<<endl;

    //! ---------------------------------------------------------------------------
    //! Starts or restarts the timer with a timeout interval of msec milliseconds.
    //! If the timer is already running, it will be stopped and restarted
    //! ---------------------------------------------------------------------------
    myTimer->setInterval(int(updateInterval*1000));
}

//! -----------------------------
//! function: clearGeneratedData
//! details:
//! -----------------------------
void SimulationManager::clearGeneratedData()
{
    QModelIndex index = myTreeView->currentIndex();
    QStandardItem *curItem = myModel->itemFromIndex(index);
    SimulationNodeClass *curNode = curItem->data(Qt::UserRole).value<SimulationNodeClass*>();
    SimulationNodeClass::nodeType type = curNode->getType();

    //! ------------------
    //! Mesh root clicked
    //! ------------------
    if(type==SimulationNodeClass::nodeType_meshControl)
    {
        for(int bodyIndex = 1; bodyIndex<=mySimulationDataBase->bodyMap.size(); bodyIndex++)
        {
            //! ----------------------------
            //! clear the mesh data sources
            //! ----------------------------
            mySimulationDataBase->ArrayOfMeshDS.insert(bodyIndex, occHandle(MeshVS_DataSource)());
            mySimulationDataBase->ArrayOfMeshDS2D.insert(bodyIndex, occHandle(MeshVS_DataSource)());

            //! ------------------------------------------
            //! delete the face and edge mesh datasources
            //! ------------------------------------------
            for(int faceNr=1; faceNr<=mySimulationDataBase->MapOfBodyTopologyMap.value(bodyIndex).faceMap.Extent(); faceNr++)
                mySimulationDataBase->ArrayOfMeshDSOnFaces.setValue(bodyIndex,faceNr, occHandle(MeshVS_DataSource)());

            for(int edgeNr=1; edgeNr<=mySimulationDataBase->MapOfBodyTopologyMap.value(bodyIndex).edgeMap.Extent(); edgeNr++)
                mySimulationDataBase->ArrayOfMeshDSOnEdges.setValue(bodyIndex,edgeNr, occHandle(MeshVS_DataSource)());

            //! ---------------------------
            //! set the flag to be updated
            //! ---------------------------
            mySimulationDataBase->ArrayOfMeshIsToBeUdpdated.insert(bodyIndex,true);
        }
        //! --------------------------------------------------------------------------
        //! refresh the mesh interactive objects stored within the QStandardItemModel
        //! --------------------------------------------------------------------------
        //this->buildMeshIO();
    }

    //! ---------------------------------------
    //! "Solution"/"Solution information" item
    //! ---------------------------------------
    if(curNode->isSolution() || curNode->isSolutionInformation())
    {
        QStandardItem *itemSolution;
        SimulationNodeClass *nodeSolution;
        SimulationNodeClass *nodeSolutionInformation;
        if(curNode->isSolution())
        {
            itemSolution = curItem;
            nodeSolution = curNode;
            nodeSolutionInformation = curItem->child(0,0)->data(Qt::UserRole).value<SimulationNodeClass*>();
        }
        if(curNode->isSolutionInformation())
        {
            itemSolution = curItem->parent();
            nodeSolution = curItem->parent()->data(Qt::UserRole).value<SimulationNodeClass*>();
            nodeSolutionInformation = curNode;
        }

        //! ------------------------------------------
        //! retrieve the project files directory path
        //! ------------------------------------------
        QString projectFiledDir = nodeSolution->getPropertyValue<QString>("Project files dir");
        if(projectFiledDir.isEmpty()) return;

        //! ----------------------------------------------------------------------
        //! if the command is issued from "Solution" delete also the result files
        //! ----------------------------------------------------------------------
        if(curNode->isSolution())
        {
            int button = QMessageBox::warning(this,APPNAME,"Warning: The results files will be deleted",QMessageBox::Ok,QMessageBox::Cancel);

            if(button==int(QMessageBox::Cancel)) return;

            //! --------------------------------------
            //! Retrieve the time tag of the analysis
            //! --------------------------------------
            QString timeTag = nodeSolution->getPropertyValue<QString>("Parent time tag");
            QString directoryToDeleteFullPath = projectFiledDir+"/SolutionData_"+timeTag;
            cout<<"____directory to delete: "<<directoryToDeleteFullPath.toStdString()<<"____"<<endl;

            QDir curDir(directoryToDeleteFullPath);
            bool isDone = curDir.removeRecursively();
            if(isDone) cout<<"____directory deleted____"<<endl;
            else cout<<"____directory not deleted____"<<endl;

            //! -------------------------------
            //! reset the solution information
            //! -------------------------------
            mainTreeTools::resetSolutionInformation(nodeSolutionInformation);
        }

        //! -----------------------------------------------------------
        //! scan the post processing items and delete the post objects
        //! -----------------------------------------------------------
        for(int n=1; n<itemSolution->rowCount(); n++)
        {
            SimulationNodeClass *nodePostProcessing = itemSolution->child(n,0)->data(Qt::UserRole).value<SimulationNodeClass*>();
            QStandardItem *itemPostObject = nodePostProcessing->getPropertyItem("Post object");
            if(itemPostObject==Q_NULLPTR) continue;
            postObject curPostObject = itemPostObject->data(Qt::UserRole).value<Property>().getData().value<postObject>();
            emit requestHideSingleResult(curPostObject);
            nodePostProcessing->removeProperty("Post object");
        }
    }

    //! ------------------------------------------
    //! post processing item - an analysis result
    //! ------------------------------------------
    if(curNode->isAnalysisResult())
    {
        QStandardItem *itemPostObject = curNode->getPropertyItem("Post object");
        if(itemPostObject==Q_NULLPTR) return;
        postObject curPostObject = itemPostObject->data(Qt::UserRole).value<Property>().getData().value<postObject>();
        emit requestHideSingleResult(curPostObject);
        curNode->removeProperty("Post object");
    }

}

//! ----------------------------------------------------
//! function: handleBoltControls
//! details:  this function enables/disables the "Load"
//!           and "Adjustment" controls, according to
//!           the bolt status definition "Define by"
//! ----------------------------------------------------
void SimulationManager::handleBoltControls()
{
    cout<<"SimulationManager::handleBoltControls()->____function called____"<<endl;

    QWidget *w = tools::getWidgetByName("messagesAndLoadsWidget");
    TableWidget *tableWidget = static_cast<TableWidget*>(w);
    CustomTableModel *tabularDataModel = static_cast<CustomTableModel*>(tableWidget->getTableView()->model());

    QExtendedStandardItem *itemBolt = static_cast<QExtendedStandardItem*>(myModel->itemFromIndex(myTreeView->currentIndex()));
    SimulationNodeClass *nodeBolt = itemBolt->data(Qt::UserRole).value<SimulationNodeClass*>();
    QExtendedStandardItem *itemBoltStatus = nodeBolt->getPropertyItem("Define by");

    //! ------------------------------------------------------------------------
    //! act on the Current step number: modify the right column in tabular data
    //! ------------------------------------------------------------------------
    int currentRow = this->getAnalysisSettingsNodeFromCurrentItem()->getPropertyValue<int>("Current step number");

    //int SC =this->calculateStartColumn();
    int SC = mainTreeTools::calculateStartColumn(myTreeView);

    QModelIndex indexBoltStatusDefinedBy = tabularDataModel->makeIndex(currentRow,SC);
    QVariant data;
    //Property::boltStatusDefinedBy boltStatusDefineBy = itemBoltStatus->data(Qt::UserRole).value<Property>().getData().value<Property::boltStatusDefinedBy>();
    Property::defineBy boltStatusDefineBy = itemBoltStatus->data(Qt::UserRole).value<Property>().getData().value<Property::defineBy>();
    data.setValue(boltStatusDefineBy);
    tabularDataModel->setData(indexBoltStatusDefinedBy,data,Qt::EditRole);
}

//! -----------------------
//! function: showElements
//! details:
//! -----------------------
void SimulationManager::showElements()
{
    cout<<"SimulationManager::showElements()->____function called____"<<endl;

    //! ---------------------
    //! retrieve all results
    //! ---------------------
    QList<postObject> postObjectList= this->retrieveAllResults();

    //! ----------------------------------------------------------------
    //! iterate over the results in order to find the displayed shapes:
    //! the shapes are shown in wireframe mode, and here must be hidden
    //! Moreover scanning the list of post object remove the mesh view
    //! ----------------------------------------------------------------
    std::vector<int> parentShapeIndexes;
    for(QList<postObject>::iterator it = postObjectList.begin(); it!=postObjectList.end(); ++it)
    {
        postObject aPostObject = *it;

        QMap<GeometryTag,QList<QMap<int,double>>> theData = aPostObject.getData();
        for(QMap<GeometryTag,QList<QMap<int,double>>>::iterator it = theData.begin(); it!= theData.end(); ++it)
        {
            GeometryTag aLoc = it.key();
            int bodyIndex = aLoc.parentShapeNr;
            parentShapeIndexes.push_back(bodyIndex);
        }
        bool showElements = true;
        aPostObject.updateView(showElements);
    }

    //! ---------------------------------------------------------
    //! clean from duplicated values the vector of parent shapes
    //! ---------------------------------------------------------
    parentShapeIndexes = tools::clearFromDuplicates(parentShapeIndexes);

    //! -------------------------------------------------------------------
    //! hide all the bodies (which are shown, by default, using wireframe)
    //! -------------------------------------------------------------------
    TColStd_ListOfInteger listOfBodies;
    for(std::vector<int>::iterator it = parentShapeIndexes.begin(); it!=parentShapeIndexes.end(); ++it) listOfBodies.Append(*it);
    emit requestHideBody(listOfBodies);

    //! ------------------------
    //! update the mesh context
    //! ------------------------
    emit requestUpdateMeshView();
}

//! ----------------------------------
//! function: showUndeformedWireframe
//! details:
//! ----------------------------------
void SimulationManager::showUndeformedWireframe()
{
    cerr<<"SimulationManager::showUndeformedWireframe()->____function called____"<<endl;

    //! ---------------------
    //! retrieve all results
    //! ---------------------
    QList<postObject> postObjectList= this->retrieveAllResults();
    std::vector<int> parentShapeIndexes;
    for(QList<postObject>::iterator it = postObjectList.begin(); it!=postObjectList.end(); ++it)
    {
        postObject aPostObject = *it;

        QMap<GeometryTag,QList<QMap<int,double>>> theData = aPostObject.getData();
        for(QMap<GeometryTag,QList<QMap<int,double>>>::iterator it = theData.begin(); it!= theData.end(); ++it)
        {
            GeometryTag aLoc = it.key();
            int bodyIndex = aLoc.parentShapeNr;
            parentShapeIndexes.push_back(bodyIndex);
        }
        bool showElements = false;
        aPostObject.updateView(showElements);
    }

    //! clean the vector of parent shapes
    parentShapeIndexes = tools::clearFromDuplicates(parentShapeIndexes);

    TColStd_ListOfInteger listOfBodies;
    for(std::vector<int>::iterator it = parentShapeIndexes.begin(); it!=parentShapeIndexes.end(); ++it) listOfBodies.Append(*it);
    emit requestShowBody(listOfBodies);

    emit requestUpdateMeshView();
}

//! ------------------------------
//! function: showUndeformedModel
//! details:
//! ------------------------------
void SimulationManager::showUndeformedModel()
{
    cout<<"SimulationManager::showUndeformedModel()->____function called____"<<endl;
    //! to do
}

//! ----------------------
//! function: noWireframe
//! details:
//! ----------------------
void SimulationManager::noWireframe()
{
    cout<<"SimulationManager::noWireframe()->____function called____"<<endl;

    //! ---------------------
    //! retrieve all results
    //! ---------------------
    QList<postObject> postObjectList= this->retrieveAllResults();

    //! ----------------------------------------------------------------
    //! iterate over the results in order to find the displayed shapes:
    //! the shapes are shown in wireframe mode, and here must be hidden
    //! Moreover scanning the list of post object remove the mesh view
    //! ----------------------------------------------------------------
    QList<int> parentShapeIndexes;
    for(QList<postObject>::iterator it = postObjectList.begin(); it!=postObjectList.end(); ++it)
    {
        postObject aPostObject = *it;
        QMap<GeometryTag,QList<QMap<int,double>>> theData = aPostObject.getData();
        for(QMap<GeometryTag,QList<QMap<int,double>>>::iterator it = theData.begin(); it!= theData.end(); ++it)
        {
            const GeometryTag &aLoc = it.key();
            int bodyIndex = aLoc.parentShapeNr;
            if(!parentShapeIndexes.contains(bodyIndex))parentShapeIndexes<<bodyIndex;
        }
        bool showElements = false;
        aPostObject.updateView(showElements);
    }

    //! -------------------------------------------------------------------
    //! hide all the bodies (which are shown, by default, using wireframe)
    //! -------------------------------------------------------------------
    TColStd_ListOfInteger listOfBodies;
    for(QList<int>::iterator it = parentShapeIndexes.begin(); it!=parentShapeIndexes.end(); ++it) listOfBodies.Append(*it);
    emit requestHideBody(listOfBodies);

    //! ------------------------
    //! update the mesh context
    //! ------------------------
    emit requestUpdateMeshView();
}

//! --------------------------------
//! function: upadtePostObjectScale
//! details:  experimental
//! --------------------------------
void SimulationManager::updatePostObjectScale(double scale)
{
    cout<<"SimulationManager::updatePostObjectScale()->____function called. Scale: "<<scale<<"____"<<endl;

    //! -------------------------------------------------
    //! get the current item and the current post object
    //! -------------------------------------------------
    SimulationNodeClass *nodePost = this->getCurrentNode();
    QExtendedStandardItem *itemPostObject = nodePost->getPropertyItem("Post object");
    if(itemPostObject!=NULL)
    {
        cout<<"SimulationManager::updatePostObjectScale()->____post object found____"<<endl;
        postObject *aPostObject = &itemPostObject->data(Qt::UserRole).value<Property>().getData().value<postObject>();
        if(!aPostObject->isEmpty())
        {
            cout<<"SimulationManager::updatePostObjectScale()->____the post object contains data____"<<endl;
            aPostObject->setScale(scale);
            aPostObject->updateScaledView();
        }
        emit requestUpdateMeshView();
    }
}

//! ------------------------------
//! function: readResultsFile
//! details:  read a results file
//! ------------------------------
void SimulationManager::readResultsFile(const QString &fileName, const QString &solutionDataDir)
{
    cout<<"SimulationManager::readResultsFile()->____function called____"<<endl;

    //! -----------------------------------------
    //! retrieve the "Solution information" item
    //! -----------------------------------------
    QStandardItem *itemCurrentRoot =myModel->itemFromIndex(myTreeView->currentIndex());
    QStandardItem *itemSolution = itemCurrentRoot->child(itemCurrentRoot->rowCount()-1,0);
    SimulationNodeClass *nodeSolution = itemSolution->data(Qt::UserRole).value<SimulationNodeClass*>();
    QStandardItem *itemSolutionInformation = itemCurrentRoot->child(itemCurrentRoot->rowCount()-1,0)->child(0,0);
    SimulationNodeClass *nodeSolutionInformation = itemSolutionInformation->data(Qt::UserRole).value<SimulationNodeClass*>();
    QVariant data;

    //! -------------------------------------------
    //! check if the model mesh has been generated
    //! (check only volume mesh)
    //! -------------------------------------------
    bool meshOK = true;
    int NbBodies = mySimulationDataBase->bodyMap.size();
    for(int i=1; i<=NbBodies; i++)
    {
        if(mySimulationDataBase->ArrayOfMeshDS.value(i).IsNull())
        {
            meshOK = false;
            break;
        }
    }
    if(meshOK == true)
    {
        //cout<<"SimulationManager::readResultsFile()->____the mesh is OK____"<<endl;
        //! -------------------
        //! read the .sta file
        //! -------------------
        QString projectFilesDir = nodeSolution->getPropertyValue<QString>("Project files dir");
        QString timeTag = nodeSolution->getPropertyValue<QString>("Parent time tag");
        QString stafile = projectFilesDir+"/SolutionData_"+timeTag+"/input.sta";
        QFile f(stafile);
        if(f.exists())
        {
            //cout<<"____.STA FILE FOUND: \""<<stafile.toStdString()<<"\"____"<<endl;
            QMap<double,QVector<int>> timeinfo;
            bool isDone = CCXTools::readsta(stafile,timeinfo);
            if(isDone)
            {
                data.setValue(timeinfo);
                nodeSolutionInformation->replaceProperty("Discrete time map",Property("Discrete time map",data,Property::PropertyGroup_Hidden));
            }
            else
            {
                QVector<int> init;
                double ini=0.0;
                init<<0.0<<0.0<<0.0;
                timeinfo.insert(ini,init);
                data.setValue(timeinfo);
                nodeSolutionInformation->replaceProperty("Discrete time map",Property("Discrete time map",data,Property::PropertyGroup_Hidden));
            }
        }

        //! ------------------------------------------------------------------------
        //! copy the .frd file (fileName) into the "<project name>_files" directory
        //! ------------------------------------------------------------------------
        QFile frdFile(fileName);
        QString copyOfFrdFile = solutionDataDir+"/input.frd";
        frdFile.copy(copyOfFrdFile);

        //! -------------
        //! set the name
        //! -------------
        myPostEngine->setResultsFile(copyOfFrdFile);

        //! -------------------------------------------------------------
        //! the method "perform" splits the .frd file into several files
        //! (myPostEngine interally calls FrdReader::split()
        //! -------------------------------------------------------------
        bool isDone = myPostEngine->perform();
        if(isDone==false) QMessageBox::critical(this,APPNAME,"Error in reading the results file",QMessageBox::Ok);
        this->evaluateAllResults();
    }
}

//! --------------------------------
//! function: previewPrismaticLayer
//! details:  experimental
//! --------------------------------
#include <prismaticlayer.h>
bool SimulationManager::previewPrismaticLayer()
{
    cout<<"SimulationManager::previewPrismaticLayer()->____function called____"<<endl;
    /*
    SimulationNodeClass *curNode = this->getCurrentNode();

    //! -------------------------------------------------------------------------------
    //! retrieve the prismatic layer parameters
    //! (generation (first layerheight/total thickness), Expansion ratio, number of layers
    //! -------------------------------------------------------------------------------
    int options = curNode->getPropertyValue<int>("Options");
    prismaticLayerParameters p;
    switch(options)
    {
    case 0:
    {
        prismaticLayer_sizing prismSizing = prismaticLayer_sizing_FirstLayerThickness;
        int NbLayers = curNode->getPropertyValue<int>("Number of layers");
        double firstLayerHeight = curNode->getPropertyValue<double>("First layer height");
        double expRatio = curNode->getPropertyValue<double>("Expansion ratio");

        p.typeOfSizing = prismSizing;
        p.NbLayers = NbLayers;
        p.firstLayerThickness = firstLayerHeight;
        p.expRatio = expRatio;
        p.totalThickness = -1;          //! the value will be changed internally
    }
        break;

    case 1:
    {
        prismaticLayer_sizing prismSizing = prismaticLayer_sizing_TotalThickness;
        int NbLayers = curNode->getPropertyValue<int>("Number of layers");
        double totalThickness = curNode->getPropertyValue<double>("Total thickness");
        double expRatio = curNode->getPropertyValue<double>("Expansion ratio");

        p.typeOfSizing = prismSizing;
        p.NbLayers = NbLayers;
        p.firstLayerThickness = -1;     //! the value will be changed internally
        p.expRatio = expRatio;
        p.totalThickness = totalThickness;
    }
        break;
    }

    //! -------------------------------------------
    //! retrieve the other parameters from the GUI
    //! -------------------------------------------
    p.numberOfModulationDiffusionSteps = curNode->getPropertyValue<int>("Modulation diffusion steps");
    p.modulationDiffusionCutoff = curNode->getPropertyValue<double>("Modulation diffusion cutoff");
    p.modulationCoefficientTransferPercentage = curNode->getPropertyValue<double>("Modulation coefficient transfer percentage");
    p.lockBoundary = curNode->getPropertyValue<double>("Lock boundary");
    if(p.lockBoundary==true) cout<<"____BOUNDARY LOCKED____"<<endl; else cout<<"____BOUNDARY FREE____"<<endl;

    //! --------------------------------------------------
    //! "Tags" of the "Scope" property (should be bodies)
    //! --------------------------------------------------
    const QVector<GeometryTag> &vecLoc = this->getCurrentNode()->getPropertyValue<QVector<GeometryTag>>("Tags");

    //! -------------------------------------------------------------
    //! "Boundary tags" of the "Boundary" property (should be faces)
    //! -------------------------------------------------------------
    const QVector<GeometryTag> &boundaryVecLoc = this->getCurrentNode()->getPropertyValue<QVector<GeometryTag>>("Boundary tags");

    //! ----------------------------------------------
    //! the prismatic faces:
    //! key => bodyIndex
    //! value => list of prismatic faces on bodyIndex
    //! ----------------------------------------------
    QMap<int,QList<int>> prismaticFaces;

    //! --------------------------------
    //! cycle over "Scope" using "Tags"
    //! --------------------------------
    for(int n=0; n<vecLoc.size(); n++)
    {
        const GeometryTag &curBodyLoc = vecLoc.at(n);
        int bodyIndex = curBodyLoc.parentShapeNr;

        //! --------------------------------------------
        //! cycle over "Boundary" using "Boundary tags"
        //! --------------------------------------------
        QList<int> prismaticFacesOnBody;
        for(int i=0; i<boundaryVecLoc.size(); i++)
        {
            const GeometryTag &aLoc = boundaryVecLoc.at(i);
            int bodyIndex2 = aLoc.parentShapeNr;
            int faceIndex2 = aLoc.subTopNr;
            if(bodyIndex==bodyIndex2) prismaticFacesOnBody<<faceIndex2;
        }
        prismaticFaces.insert(bodyIndex,prismaticFacesOnBody);
    }

    //! ----------------------
    //! debug window messages
    //! ----------------------
    for(QMap<int,QList<int>>::iterator it = prismaticFaces.begin(); it!=prismaticFaces.cend(); ++it)
    {
        int bodyIndex = it.key();
        QList<int> listOfFaceNumbers = it.value();
        ccout(QString("____prismatic layers. List of faces for body index: %1____").arg(bodyIndex));
        for(int n=0; n<listOfFaceNumbers.length(); n++) ccout(QString("____face nr: %1____").arg(listOfFaceNumbers.at(n)));
    }
    ccout("____end of the list____");

    //! -----------------------------------
    //! create the prismatic layer builder
    //! -----------------------------------
    prismaticLayer plbuilder(mySimulationDataBase);
    plbuilder.setParameters(p);

    //! ----------------------
    //! cycle over the bodies
    //! ----------------------
    QList<occHandle(Ng_MeshVS_DataSourceFace)> meshes;
    for(QMap<int,QList<int>>::iterator it = prismaticFaces.begin(); it!=prismaticFaces.end(); ++it)
    {
        int bodyIndex = it.key();
        plbuilder.setBody(bodyIndex);
        QList<int> prismaticFacesOnBody = prismaticFaces.value(bodyIndex);
        plbuilder.setPrismaticFaces(prismaticFacesOnBody);

        //! --------------------------------
        //! displace the previous displaced
        //! --------------------------------
        plbuilder.inflateMesh(meshes);
    }
    //! ------------------------------------------------------
    //! test volume mesh starting from the last inflated mesh
    //! ------------------------------------------------------
    //MeshTools::surfaceMeshDSToNetgenMesh(meshes.last());

    //! -----------------------------------
    //! test prismatic layer mesh building
    //! -----------------------------------
    occHandle(Ng_MeshVS_DataSource3D) prismaticMeshDS;

    prismaticLayer::buildPrismaticElements(meshes,prismaticMeshDS);
    //plbuilder.buildPrismaticElements(meshes,prismaticMeshDS);
    prismaticLayer::displayMesh(prismaticMeshDS);
    */
    return true;
}

//! ------------------------
//! function: displayMarker
//! details:
//! ------------------------
void SimulationManager::displayMarker()
{
    //cout<<"SimulationManager::displayMarker()->____function called____"<<endl;
    SimulationNodeClass *curNode = this->getCurrentNode();
    if(curNode->getPropertyItem("Graphic object")!=NULL)
    {
        occHandle(AIS_Shape) theMarker;
        switch(curNode->getType())
        {
        case SimulationNodeClass::nodeType_structuralAnalysisBoltPretension:
            theMarker = curNode->getPropertyValue<AIS_DoubleArrowMarker_handle_reg>("Graphic object");
            break;

        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Acceleration:
            theMarker = curNode->getPropertyValue<AIS_ArrowMarker_handle_reg>("Graphic object");
            break;

        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Moment:
            theMarker = curNode->getPropertyValue<AIS_CurvedArrowMarker_handle_reg>("Graphic object");
            break;

        case SimulationNodeClass::nodeType_remotePoint:
        case SimulationNodeClass::nodeType_pointMass:
            theMarker = curNode->getPropertyValue<AIS_SphereMarker_handle_reg>("Graphic object");
            break;
        }
        theMarker->SetZLayer(Graphic3d_ZLayerId_TopOSD);
        myCTX->Display(theMarker,1,-1,true,false);
    }
}

//! ----------------------
//! function: buildMeshIO
//! details:
//! ----------------------
void SimulationManager::buildMeshIO()
{
    for(QMap<int,TopoDS_Shape>::iterator it=mySimulationDataBase->bodyMap.begin(); it!=mySimulationDataBase->bodyMap.end(); it++)
    {
        int bodyIndex = it.key();
        for(int k=0; k<Geometry_RootItem->rowCount(); k++)
        {
            SimulationNodeClass *curNode = Geometry_RootItem->child(k)->data(Qt::UserRole).value<SimulationNodeClass*>();
            if(curNode->getType()==SimulationNodeClass::nodeType_pointMass) continue;
            int mapIndex = curNode->getPropertyValue<int>("Map index");
            if(mapIndex==bodyIndex)
            {
                occHandle(MeshVS_DataSource) curDataSource3D = mySimulationDataBase->ArrayOfMeshDS.value(bodyIndex);
                occHandle(MeshVS_DataSource) curDataSource2D = mySimulationDataBase->ArrayOfMeshDS2D.value(bodyIndex);

                QVariant data;
                if(!curDataSource3D.IsNull())
                {
                    MeshVS_Mesh_handle_reg aMeshVS_Mesh3D = new MeshVS_Mesh;
                    aMeshVS_Mesh3D->SetDataSource(curDataSource3D);
                    data.setValue(aMeshVS_Mesh3D);
                    Property prop_meshVS3D("3D mesh",data,Property::PropertyGroup_MeshDataSources);
                    if(curNode->getPropertyItem("3D mesh")==NULL) curNode->addProperty(prop_meshVS3D);
                    else curNode->replaceProperty("3D mesh",prop_meshVS3D);
                }
                else
                {
                    curNode->removeProperty("3D mesh");
                }
                if(!curDataSource2D.IsNull())
                {
                    MeshVS_Mesh_handle_reg aMeshVS_Mesh2D = new MeshVS_Mesh;
                    aMeshVS_Mesh2D->SetDataSource(curDataSource2D);
                    data.setValue(aMeshVS_Mesh2D);
                    Property prop_meshVS2D("2D mesh",data,Property::PropertyGroup_MeshDataSources);
                    if(curNode->getPropertyItem("2D mesh")==NULL) curNode->addProperty(prop_meshVS2D);
                    else curNode->replaceProperty("2D mesh",prop_meshVS2D);
                }
                else
                {
                    curNode->removeProperty("2D mesh");
                }
            }
        }
    }
}

//! -----------------------------------------------------------
//! function: renameItemBasedOnDefinition
//! details:  change the display name of the item on the basis
//!           of the node type and locatione
//! -----------------------------------------------------------
void SimulationManager::renameItemBasedOnDefinition()
{
    SimulationNodeClass *node = myTreeView->currentIndex().data(Qt::UserRole).value<SimulationNodeClass*>();

    if(node->getType() == SimulationNodeClass::nodeType_connectionPair)
    {
        //! ------------------------------------------------
        //! build a name for a contact pair
        //! do not rename if master or slave tags are empty
        //! ------------------------------------------------
        const QVector<GeometryTag> &vecLocMasterTags = node->getPropertyValue<QVector<GeometryTag>>("Tags master");
        if(vecLocMasterTags.size()==0) return;
        const QVector<GeometryTag> &vecLocSlaveTags = node->getPropertyValue<QVector<GeometryTag>>("Tags slave");
        if(vecLocSlaveTags.size()==0) return;

        //! -------------
        //! build a name
        //! -------------
        QString locationNames = "";
        if(vecLocSlaveTags.size()>1) locationNames.append("Multiple to ");
        else
        {
            int bodyIndex = vecLocSlaveTags.at(0).parentShapeNr;
            locationNames.append(mySimulationDataBase->MapOfBodyNames.value(bodyIndex)).append(" to ");
        }
        if(vecLocMasterTags.size()>1) locationNames.append("multiple");
        else
        {
            int bodyIndex = vecLocMasterTags.at(0).parentShapeNr;
            locationNames.append(mySimulationDataBase->MapOfBodyNames.value(bodyIndex));
        }

        QString itemName = locationNames;
        QStandardItem *curItem = myModel->itemFromIndex(myTreeView->currentIndex());
        curItem->setEditable(true);
        QVariant data;
        data.setValue(itemName);
        curItem->setData(data,Qt::DisplayRole);
        node->setName(itemName);
        curItem->setEditable(false);
    }
    else
    {
        QString controlName;

        switch(node->getType())
        {
        case SimulationNodeClass::nodeType_meshBodyMeshControl: controlName = QString("Body sizing"); break;
        case SimulationNodeClass::nodeType_meshMethod: controlName = QString("Meshing method"); break;
        case SimulationNodeClass::nodeType_meshFaceSize: controlName = QString("Face sizing"); break;
        case SimulationNodeClass::nodeType_meshEdgeSize: controlName = QString("Edge sizing"); break;
        case SimulationNodeClass::nodeType_meshVertexSize: controlName = QString("Vertex sizing"); break;
        case SimulationNodeClass::nodeType_meshPrismaticLayer: controlName = QString("Prismatic layer meshing"); break;
        case SimulationNodeClass::nodeType_connectionPair: controlName = QString("Contact pair"); break;
        case SimulationNodeClass::nodeType_solutionThermalTemperature: controlName = QString("Temperature distribution"); break;
        case SimulationNodeClass::nodeType_solutionThermalFlux: controlName = QString("Thermal flux"); break;
        case SimulationNodeClass::nodeType_thermalAnalysisAdiabaticWall: controlName = QString("Adiabatic wall"); break;
        case SimulationNodeClass::nodeType_thermalAnalysisRadiation: controlName = QString("Thermal radiation"); break;
        case SimulationNodeClass::nodeType_thermalAnalysisTemperature: controlName = QString("Temperature"); break;
        case SimulationNodeClass::nodeType_thermalAnalysisThermalFlow: controlName = QString("Thermal flow"); break;
        case SimulationNodeClass::nodeType_thermalAnalysisThermalFlux: controlName = QString("Thermal flux"); break;
        case SimulationNodeClass::nodeType_thermalAnalysisConvection: controlName = QString("Convection"); break;
        case SimulationNodeClass::nodeType_structuralAnalysisBoltPretension: controlName = QString("Bolt preload"); break;
        case SimulationNodeClass::nodeType_coordinateSystem: controlName = QString("Coordinate system"); break;
        }

        //! -------------
        //! build a name
        //! -------------
        QString locationNames;
        int bodyIndex;
        const QVector<GeometryTag> &locs = node->getPropertyValue<QVector<GeometryTag>>("Tags");

        //! ----------------------------------
        //! do not remane if there is not tag
        //! ----------------------------------
        if(locs.size()==0) return;

        int i;
        for(i=0; i<locs.length()-1; i++)
        {
            const GeometryTag &aLoc = locs.at(i);
            bodyIndex = aLoc.parentShapeNr;
            locationNames.append(mySimulationDataBase->MapOfBodyNames.value(bodyIndex)).append(", ");
        }
        const GeometryTag &aLoc = locs.at(i);
        bodyIndex = aLoc.parentShapeNr;
        locationNames.append(mySimulationDataBase->MapOfBodyNames.value(bodyIndex));

        QString itemName = controlName + " on " + locationNames;
        QStandardItem *curItem = myModel->itemFromIndex(myTreeView->currentIndex());
        curItem->setEditable(true);
        QVariant data;
        data.setValue(itemName);
        curItem->setData(data,Qt::DisplayRole);
        node->setName(itemName);
        curItem->setEditable(false);
    }
}

#ifdef COSTAMP_VERSION
//! -------------------------------
//! function: startTimeStepBuilder
//! details:
//! -------------------------------
void SimulationManager::COSTAMP_startTimeStepBuilder()
{
    cout<<"SimulationManager::startTimeStepBuilder()->____function called____"<<endl;
    SimulationNodeClass *curNode = myTreeView->currentIndex().data(Qt::UserRole).value<SimulationNodeClass*>();
    const QString &timeHistoryFileLoc = curNode->getPropertyValue<QString>("Time history file");    
    QString program = QString("D:/Work/Qt/build_pro26.0_OCC7.3.0/release/TimeStepBuilder.exe");
    QStringList arguments;
    arguments<<myCurrentProjectDir<<timeHistoryFileLoc;
    QProcess *tsbProcess = new QProcess(this);
    tsbProcess->start(program,arguments);
    tsbProcess->waitForFinished(-1);

    //! -------------------------------------------
    //! remove previous items if present. To do...
    //! -------------------------------------------

    bool returnCode = this->COSTAMP_addProcessParameters();
    if(returnCode == false) cout<<"____addProcessParameters error___"<<endl;
}

void SimulationManager::COSTAMP_changeTimeHistoryFile(const QString &fileLocation)
{
    SimulationNodeClass *curNode = myTreeView->currentIndex().data(Qt::UserRole).value<SimulationNodeClass*>();
    QVariant data;
    data.setValue(fileLocation);
    Property prop_fileLoc("Time history file",data,Property::PropertyGroup_Definition);
    curNode->replaceProperty("Time history file",prop_fileLoc);
}

bool SimulationManager::COSTAMP_addProcessParameters()
{
    cout<<"SimulationManager::COSTAMP_addProcessParameters()->____function called____"<<endl;
    //tSbList.clear();

    //! Path of the configuration file
    QString dirPath = myCurrentProjectDir;
    cout<<"SimulationManager::COSTAMP_addProcessParameters()->____dirPath "<<myCurrentProjectDir.toStdString()<<endl;
    QString tsbFile= dirPath+"/timepoints.out";
    //! Path of the OF folder
    SimulationNodeClass *tsbNode = myTreeView->currentIndex().data(Qt::UserRole).value<SimulationNodeClass*>();
    QString &timeHistoryFileLoc = tsbNode->getPropertyValue<QString>("Time history file");
    timeHistoryFileLoc.chop(timeHistoryFileLoc.split("/").last().length());

    //! Path of the OF mapped data
    QDir dir;
    dir.current();
    dir.cd(timeHistoryFileLoc);
    dir.mkdir("Mapped");
    const QString mappedFilePath = timeHistoryFileLoc+"Mapped";
    cout<<"SimulationManager::COSTAMP_addProcessParameters()->____config file "<<tsbFile.toStdString()<<endl;
    //! ---------------------------------
    //! read the configuration file
    //! detils: config file must exist
    //! ---------------------------------
    ifstream is;
    QFile file(tsbFile);
    if(!file.exists())
        return false;
    else
    {
        is.open(tsbFile.toStdString());
        std::string val;
        //! timeStepType
        //!  0: ClosedAssemblyWithoutPressure,
        //!  1: ClosedAssemblyWithPressure,
        //!  2: OpenAssembly,
        std::vector<int> timeStepNr,type;
        std::vector<double>  prevTime,curTime;
        int n=0;
        if(is.is_open())
            while(!is.eof())
            {
                std::getline(is,val);
                int tStepNr,tsType;
                double pTime,cTime;
                if(4 == sscanf(val.c_str(),"%d%d%lf%lf",&tStepNr,&tsType, &pTime,&cTime))
                {
                    //! insert the current row into the vectors
                    timeStepNr.push_back(tStepNr);
                    type.push_back(tsType);
                    prevTime.push_back(pTime);
                    curTime.push_back(cTime);
                    n++;
                }
            }
        is.close();
        tSbList = curTime;
        QVariant data;
        int closureIndex, prexIndex, modelChangeIndex,tSbIndex,mapperIndex;
        closureIndex = -1;
        prexIndex = -1;
        modelChangeIndex = -1;

        //! ------------------------------------------------------------
        //! for "createSimulationNode()" which needs the "current" item
        //! ------------------------------------------------------------
        myTreeView->setCurrentIndex(StaticAnalysis_RootItem->index().child(0,0));
        int curRow = 1;
        tSbIndex = curRow;
        cout<<"curRow= "<<tSbIndex<<endl;
        SimulationNodeClass *nodeAnalysisSettings = StaticAnalysis_RootItem->child(0,0)->data(Qt::UserRole).value<SimulationNodeClass*>();
        CustomTableModel *tabData = nodeAnalysisSettings->getTabularDataModel();

        //! ------------------------------------------
        //! set the number of step and timestep policy
        //! detils:
        //! -------------------------------------------
        nodeAnalysisSettings->getModel()->blockSignals(true);
        int NbTstep = int(timeStepNr.size());
        Property property_numberOfSteps("Number of steps",data,Property::PropertyGroup_StepControls);
        nodeAnalysisSettings->replaceProperty("Number of steps",property_numberOfSteps);
        this->resizeTabularData();
        nodeAnalysisSettings->getModel()->blockSignals(false);

        for(int i=0; i<NbTstep;i++)
        {
            tabData->setDataRC(curTime.at(i),i+1,1,Qt::EditRole);
            data.setValue(NbTstep);           //! the default Number of steps
        }
        curRow++;
        cout<<"curRow= "<<curRow<<endl;

        //! -----------------------
        //! create a fixed support
        //! -----------------------
        this->createSimulationNode(SimulationNodeClass::nodeType_structuralAnalysisBoundaryContidion_FixedSupport);
        curRow++;

        //! -------------------------------
        //! create mapper and OFtranslator
        //! -------------------------------
        this->createSimulationNode(SimulationNodeClass::nodeType_mapper);
        mapperIndex = curRow;
        cout<<"curRow= "<<mapperIndex<<endl;
        curRow++;
        myTreeView->setCurrentIndex(StaticAnalysis_RootItem->index().child(mapperIndex,0));
        SimulationNodeClass *mapperNode = myTreeView->currentIndex().data(Qt::UserRole).value<SimulationNodeClass*>();
        QExtendedStandardItem *mapperItem = this->getTreeItem(mapperNode->getType());
        this->createSimulationNode(SimulationNodeClass::nodeType_OpenFoamScalarData);
        myTreeView->setCurrentIndex(mapperItem->index().child(0,0));
        SimulationNodeClass *ofNode = myTreeView->currentIndex().data(Qt::UserRole).value<SimulationNodeClass*>();
        ofNode->getModel()->blockSignals(true);
        data.setValue(mappedFilePath);       //! target directory
        Property property_targetDir("Target directory",data,Property::PropertyGroup_Definition);
        ofNode->replaceProperty("Target directory",property_targetDir);
        data.setValue(timeHistoryFileLoc);       //! source directory
        Property property_sourceDir("Source directory",data,Property::PropertyGroup_Definition);
        ofNode->replaceProperty("Source directory",property_sourceDir);
        data.setValue(0);       //! split in single file
        Property property_split("Split data",data,Property::PropertyGroup_OutputSettings);
        ofNode->replaceProperty("Split data",property_split);
        ofNode->getModel()->blockSignals(false);

        myTreeView->setCurrentIndex(StaticAnalysis_RootItem->index().child(mapperIndex,0));
        this->createSimulationNode(SimulationNodeClass::nodeType_importedBodyScalar);
        myTreeView->setCurrentIndex(mapperItem->index().child(1,0));
        SimulationNodeClass *importedBSNode = myTreeView->currentIndex().data(Qt::UserRole).value<SimulationNodeClass*>();
        importedBSNode->getModel()->blockSignals(true);
        //! Automatic time Hystory
        int stepSM = 4;
        data.setValue(stepSM);
        Property property_sSm("Step selection mode",data,Property::PropertyGroup_Definition);
        importedBSNode->replaceProperty("Step selection mode",property_sSm);
        //! Mapping Algo
        int algo = 1;
        data.setValue(algo);
        Property property_algo("Algorithm",data,Property::PropertyGroup_Advanced);
        importedBSNode->replaceProperty("Algorithm",property_algo);
        //! Time hystory path
        importedBSNode->removeProperty("Source file");
        if(importedBSNode->getPropertyItem("Source directory")==NULL)
        {
            data.setValue(mappedFilePath);
            Property prop_sourceDir("Source directory",data,Property::PropertyGroup_Definition);
            importedBSNode->addProperty(prop_sourceDir,2);
        }
        //! Pinball
        int pin = 50;
        data.setValue(pin);       //! Automatic time Hystory
        Property property_pin("Pinball",data,Property::PropertyGroup_Advanced);
        importedBSNode->replaceProperty("Pinball",property_pin);
        importedBSNode->getModel()->blockSignals(false);

        //! -------------------------------
        //! create load boundary condition
        //! -------------------------------
        int nBclosure = 0;
        int nBpressure = 0;
        int nBopen = 0;
        myTreeView->setCurrentIndex(StaticAnalysis_RootItem->index().child(0,0));
        for(int i=0; i<NbTstep;i++)
        {
            //! ------------------------------
            //! create the force closure node
            //! ------------------------------
            if(type.at(i) == 0 && nBclosure == 0)
            {
                this->createSimulationNode(SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteForce);
                closureIndex = curRow;
                cout<<"curRow= "<<closureIndex<<endl;

                curRow++;
                QStandardItem *curItem =StaticAnalysis_RootItem->child(closureIndex,0);
                SimulationNodeClass *curNode = curItem->data(Qt::UserRole).value<SimulationNodeClass*>();
                curNode->getModel()->blockSignals(true);
                QString newName1="Closure Force";
                curNode->setName(newName1);
                data.setValue(newName1);
                curItem->setData(data,Qt::DisplayRole);
                data.setValue(Property::loadDefinition_tabularData);
                Property prop_loadMagnitude("Magnitude",data,Property::PropertyGroup_Definition);
                curNode->replaceProperty("Magnitude",prop_loadMagnitude);
                nBclosure++;
                curNode->getModel()->blockSignals(false);
            }
            //! -------------------------------
            //! create the inner pressure node
            //! -------------------------------
            if(type.at(i)==1 && nBpressure == 0)
            {
                this->createSimulationNode(SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Pressure);
                prexIndex = curRow;
                cout<<"curRow= "<<prexIndex<<endl;

                curRow++;
                QStandardItem *curItem =StaticAnalysis_RootItem->child(prexIndex,0);
                SimulationNodeClass *curNode = curItem->data(Qt::UserRole).value<SimulationNodeClass*>();
                curNode->getModel()->blockSignals(true);
                QString newName2="Inner Pressure";
                curNode->setName(newName2);
                data.setValue(newName2);
                curItem->setData(data,Qt::DisplayRole);
                data.setValue(Property::loadDefinition_tabularData);
                Property prop_loadMagnitude("Magnitude",data,Property::PropertyGroup_Definition);
                curNode->replaceProperty("Magnitude",prop_loadMagnitude);
                curNode->getModel()->blockSignals(false);
                nBpressure++;
                curNode->getModel()->blockSignals(false);
            }
            //! -----------------------
            //! create modelChange
            //! -----------------------
            if(type.at(i)==2 && nBopen == 0)
            {
                this->createSimulationNode(SimulationNodeClass::nodeType_modelChange);
                modelChangeIndex = curRow;
                cout<<"curRow= "<<modelChangeIndex<<endl;

                curRow++;
                QStandardItem *curItem =StaticAnalysis_RootItem->child(modelChangeIndex,0);
                SimulationNodeClass *curNode = curItem->data(Qt::UserRole).value<SimulationNodeClass*>();
                curNode->getModel()->blockSignals(true);
                data.setValue(1);   //! contact
                Property prop_itemType("Item type",data,Property::PropertyGroup_Definition);
                curNode->replaceProperty("Item type",prop_itemType);
                curNode->removeProperty("Scoping method");
                curNode->removeProperty("Tags");
                data.setValue(0);
                Property prop_contact("Contact",data,Property::PropertyGroup_Scope);
                curNode->replaceProperty("Geometry",prop_contact);
                nBopen++;
                curNode->getModel()->blockSignals(false);
            }
        }
        for(int i=0; i<NbTstep;i++)
        {
            int stepNb = timeStepNr.at(i);
            if(type.at(i)!=2)
            {
                if(closureIndex!=-1)
                {
                    myTreeView->setCurrentIndex(StaticAnalysis_RootItem->index().child(closureIndex,0));
                    double force = 18000000.0;
                    QList<int> columns = mainTreeTools::getColumnsToRead(myTreeView);
                    tabData->setDataRC(force,stepNb,columns.at(0),Qt::EditRole);
                    cout<<"closureIndex "<<closureIndex<<" column n "<<columns.at(0)<<endl;
                }
                if(type.at(i)==0 && prexIndex!=-1)
                {
                    myTreeView->setCurrentIndex(StaticAnalysis_RootItem->index().child(prexIndex,0));
                    double prex = 0;
                    QList<int> columns = mainTreeTools::getColumnsToRead(myTreeView);
                    tabData->setDataRC(prex,stepNb,columns.at(0),Qt::EditRole);
                }
                if(type.at(i)==1)
                {
                    myTreeView->setCurrentIndex(StaticAnalysis_RootItem->index().child(prexIndex,0));
                    double prex = 60;
                    QList<int> columns = mainTreeTools::getColumnsToRead(myTreeView);
                    tabData->setDataRC(prex,stepNb,columns.at(0),Qt::EditRole);
                }
            }
            else if(type.at(i)==2)
            {
                myTreeView->setCurrentIndex(StaticAnalysis_RootItem->index().child(modelChangeIndex,0));
                QList<int> columns = mainTreeTools::getColumnsToRead(myTreeView);
                int mChangeValue=-1;
                tabData->setDataRC(mChangeValue,stepNb,columns.at(0),Qt::EditRole);
                if(prexIndex!=-1)
                {
                    myTreeView->setCurrentIndex(StaticAnalysis_RootItem->index().child(prexIndex,0));
                    double prex = 0;
                    QList<int> columns = mainTreeTools::getColumnsToRead(myTreeView);
                    tabData->setDataRC(prex,stepNb,columns.at(0),Qt::EditRole);
                }
                if(closureIndex!=-1)
                {
                    myTreeView->setCurrentIndex(StaticAnalysis_RootItem->index().child(closureIndex,0));
                    double load = 0;
                    QList<int> columns = mainTreeTools::getColumnsToRead(myTreeView);
                    tabData->setDataRC(load,stepNb,columns.at(0),Qt::EditRole);
                }
            }
        }
        myTreeView->setCurrentIndex(StaticAnalysis_RootItem->index().child(0,0));
        return true;
    }
}
#endif

//! ------------------------------
//! function: resetAndUpdateModel
//! details:
//! ------------------------------
void SimulationManager::resetAndUpdateModel()
{
    cout<<"SimulationManager::updateModel()->____function called____"<<endl;

    //! ----------------------------------------------------------
    //! reset the database
    //!
    //! - reset the geometry data base: it removes the geometry
    //!   items and clear the data base inner maps
    //! - reset the map of mesh controls: it clears the data base
    //!   inner maps; it does not remove the mesh items
    //! ----------------------------------------------------------
    cout<<"SimulationManager::updateModel()->____resetting database____"<<endl;
    ccout("SimulationManager::updateModel()->____resetting database____");
    mySimulationDataBase->resetDataBase();

    //! -------------------------
    //! reload the geometry file
    //! -------------------------
    QStandardItem *item = getTreeItem(SimulationNodeClass::nodeType_import);
    SimulationNodeClass *nodeImport = item->data(Qt::UserRole).value<SimulationNodeClass*>();
    QString filePath = nodeImport->getPropertyValue<QString>("Source file path");

    Handle (QOccProgressIndicator) aProgressIndicator = new QOccProgressIndicator(0,100,this,0);

    //! --------------------
    //! call the CAD reader
    //! --------------------
    TopoDS_Compound shapeFromReader;
    QList<QString> listOfNames;
    loadCADModel(filePath,shapeFromReader,listOfNames,aProgressIndicator);

    //! ------------------------------------------------------------------
    //! check the list of names: if a name is empty assign a default name
    //! ------------------------------------------------------------------
    for(int i=0, index = 0; i<listOfNames.length(); i++)
    {
        if(listOfNames.at(i).isEmpty())
        {
            index++;
            listOfNames.replace(i,QString("body %1").arg(index));
        }
    }

    //! ---------------------------------------------------
    //! update the geometry data base: set the main shape,
    //! rebuild the geometry items
    //! ---------------------------------------------------
    mySimulationDataBase->update(shapeFromReader);

    if(!shapeFromReader.IsNull())
    {
        //! ---------------------------------------------------
        //! fill the array of the names read from the CAD file
        //! ---------------------------------------------------
        int index = 1;
        for(int i=0; i<listOfNames.length(); i++, index++)
        {
            const QString &name = listOfNames.at(i);
            mySimulationDataBase->MapOfBodyNames.insert(index,name);
            ccout(QString("SimulationManager::createSimulationDataBase()->____found body with name: ").append(name).append("____"));
        }
        //! -----------------------------------------------
        //! update the names of the items and of the nodes
        //! -----------------------------------------------
        mySimulationDataBase->transferNames();
    }

    //! -----------------
    //! update viewports
    //! -----------------
    emit requestUpdateViewport();
}

//! ------------------------------------
//! function: customMesherBuildFaceMesh
//! details:
//! ------------------------------------
#include "customMesher/custommesher.h"
void SimulationManager::customMesherBuildFaceMesh(const TopoDS_Face &aFace)
{
    cout<<"SimulationManager::customMesherBuildFaceMesh()->____function called____"<<endl;
    CustomMesher aCustomMesher(aFace,this->getDataBase());
    QMap<GeometryTag,occHandle(Ng_MeshVS_DataSourceFace)> faceMeshDSs = aCustomMesher.meshAllFaces();
    this->displayFaceMesh(faceMeshDSs.first());
    myCTX->ClearSelected(true);
}

//! --------------------------
//! function: displayFaceMesh
//! details:  utility
//! --------------------------
void SimulationManager::displayFaceMesh(const occHandle(MeshVS_DataSource) &aMeshDS,
                                        Quantity_NameOfColor aColorName,
                                        bool showMeshEdges)
{
    occHandle(MeshVS_Mesh) meshIO = new MeshVS_Mesh();
    meshIO->SetDataSource(aMeshDS);
    occHandle(MeshVS_MeshPrsBuilder) aBuilder = new MeshVS_MeshPrsBuilder(meshIO);
    meshIO->AddBuilder(aBuilder,false);
    Graphic3d_MaterialAspect aspect(Graphic3d_NOM_GOLD);
    aspect.SetColor(aColorName);

    meshIO->GetDrawer()->SetMaterial(MeshVS_DA_FrontMaterial,aspect);
    meshIO->GetDrawer()->SetBoolean(MeshVS_DA_DisplayNodes,false);
    meshIO->GetDrawer()->SetBoolean(MeshVS_DA_ShowEdges,showMeshEdges);
    meshIO->GetDrawer()->SetColor(MeshVS_DA_EdgeColor,Quantity_NOC_BLACK);
    meshIO->SetDisplayMode(MeshVS_DMF_Shading);
    myCTX->Display(meshIO,true);
}

//! -------------------------------------
//! function: createAutomaticConnections
//! details:
//! -------------------------------------
#include <contactFinder.h>
#include "connectionpairgenerationoptions.h"
void SimulationManager::createAutomaticConnections()
{
    cout<<"SimulationManager::createAutomaticConnections()->____function called____"<<endl;

    SimulationNodeClass *curNode = myTreeView->currentIndex().data(Qt::UserRole).value<SimulationNodeClass*>();
    QExtendedStandardItem *itemTags = curNode->getPropertyItem("Tags");
    QVector<GeometryTag> vecLoc = itemTags->data(Qt::UserRole).value<Property>().getData().value<QVector<GeometryTag>>();

    if(vecLoc.isEmpty()) return;

    //! --------------------------------------
    //! tolerance for contact pairs detection
    //! --------------------------------------
    double tolerance = curNode->getPropertyValue<double>("Tolerance");

    //! ---------------------------------
    //! time tag of the connection group
    //! ---------------------------------
    QString parentTimeTag = curNode->getPropertyValue<QString>("Time tag");

    //! ---------------------------------------------------
    //! grouping options - compatibility with old versions
    //! default "By master body"
    //! ---------------------------------------------------
    int grouping = 1;
    if(curNode->getPropertyItem("Grouping")!=Q_NULLPTR) grouping = curNode->getPropertyValue<int>("Grouping");

    //! ----------------------------------------------------
    //! angular criterion - compatibility with old versions
    //! default angle 30 degrees
    //! ----------------------------------------------------
    double angularCriterion = 30;
    if(curNode->getPropertyItem("Angular criterion")!=Q_NULLPTR) curNode->getPropertyValue<double>("Angular criterion");

    //! -----------
    //! make pairs
    //! -----------
    std::vector<std::pair<GeometryTag,GeometryTag>> vectorOfTagPairs;
    size_t NbBodies = vecLoc.size();
    for(int i=1; i<NbBodies; i++)
    {
        const GeometryTag &firstTag = vecLoc.at(i-1);
        for(int j=i+1; j<=NbBodies; j++)
        {
            //cout<<"____(i,j) = ("<<i<<", "<<j<<")____"<<endl;
            const GeometryTag &secondTag = vecLoc.at(j-1);
            std::pair<GeometryTag,GeometryTag> aTagPair;
            aTagPair.first = firstTag; aTagPair.second = secondTag;
            vectorOfTagPairs.push_back(aTagPair);
        }
    }

    //! --------------------------------------------------------------------
    //! sanity check - this could be removed when extending to self-contact
    //! detection, indeed for that case a single body is enough
    //! --------------------------------------------------------------------
    if(NbBodies<=1)
    {
        QMessageBox::information(this,"Contact finder","Insert at least two bodies in selector",QMessageBox::Ok);
        return;
    }

    //! ---------------------------------------------
    //! define the result: a vector of mesh pairs
    //! indexed as the input vector of geometry tags
    //! ---------------------------------------------
    std::vector<std::pair<QVector<GeometryTag>,QVector<GeometryTag>>> allContactPairs;

    //! -------------------------------------------------------------
    //! create an instance of contactFinder
    //! Note: using the default constructor the angular criterion is
    //! inizialied as angle = 20.0
    //! -------------------------------------------------------------
    contactFinder aContactFinder;
    aContactFinder.setBodyPairs(vectorOfTagPairs);
    aContactFinder.setDataBase(this->getDataBase());
    aContactFinder.setAngularCriterion(angularCriterion);

    //! -------------------------------------------------------
    //! set the progress indicator for the contactFinder tool
    //! -------------------------------------------------------
    QProgressIndicator *myProgressIndicator = static_cast<QProgressIndicator*>(tools::getWidgetByName("progressIndicator"));
    if(myProgressIndicator!=Q_NULLPTR)
    {
        myProgressIndicator->setSecondaryBarVisible(true);
        aContactFinder.setProgressIndicator(myProgressIndicator);
    }

    //! --------------------------------------------------------------------------
    //! perform: if the process has been intentionally stopped by the user return
    //! grouping:
    //! "0" => by bodies
    //! "1" => by master face
    //! "2" => by slave face
    //! "3" => by master body
    //! "4" => by slave body
    //! "5" => none - ungrouped
    //! --------------------------------------------------------------------------
    bool stopped = aContactFinder.perform(vectorOfTagPairs,allContactPairs,tolerance,grouping);
    if(stopped) return;

    //! -----------------------
    //! create the model items
    //! -----------------------
    for(std::vector<std::pair<QVector<GeometryTag>,QVector<GeometryTag>>>::iterator it = allContactPairs.begin(); it!=allContactPairs.end(); it++)
    {
        const std::pair<QVector<GeometryTag>,QVector<GeometryTag>> &curPair = *it;
        int masterBodyIndex = curPair.first[0].parentShapeNr;
        int slaveBodyIndex = curPair.second[0].parentShapeNr;

        //! -------------------------------
        //! a name for the connection item
        //! -------------------------------
        QString masterName = this->getDataBase()->MapOfBodyNames.value(masterBodyIndex);
        QString slaveName = this->getDataBase()->MapOfBodyNames.value(slaveBodyIndex);
        QString connectionName = slaveName + " to " + masterName;

        //cout<<"____creating item: (masterBodyIndex, slaveBodyIndex) = ("<<masterBodyIndex<<", "<<slaveBodyIndex<<")____"<<endl;

        //! -----------------------------------------------------
        //! connection pair generation option: "isManual==false"
        //! -----------------------------------------------------
        QVariant data;
        connectionPairGenerationOption options;
        options.manual = false;
        options.connectionPairName = connectionName;
        data.setValue(options);
        SimulationNodeClass *aContactNode = nodeFactory::nodeFromScratch(SimulationNodeClass::nodeType_connectionPair,0,0,data);

        aContactNode->getModel()->blockSignals(true);
        data.setValue(curPair.first);
        Property prop_master("Master",data,Property::PropertyGroup_Scope);
        Property prop_tagsMaster("Tags master",data,Property::PropertyGroup_Scope);
        aContactNode->replaceProperty("Master",prop_master);
        aContactNode->replaceProperty("Tags master",prop_tagsMaster);
        data.setValue(curPair.second);
        Property prop_slave("Slave",data,Property::PropertyGroup_Scope);
        Property prop_tagsSlave("Tags slave",data,Property::PropertyGroup_Scope);
        aContactNode->replaceProperty("Slave",prop_slave);
        aContactNode->replaceProperty("Tags slave",prop_tagsSlave);
        data.setValue(parentTimeTag);
        aContactNode->replaceProperty("Parent time tag",Property("Parent time tag",data,Property::PropertyGroup_Identifier));
        aContactNode->getModel()->blockSignals(false);

        connect(aContactNode->getModel(),SIGNAL(itemChanged(QStandardItem*)),this,SLOT(handleItemChange(QStandardItem*)));

        //! --------------
        //! item creation
        //! --------------
        QExtendedStandardItem *itemConnection = new QExtendedStandardItem();
        data.setValue(connectionName);
        itemConnection->setData(data,Qt::DisplayRole);
        data.setValue(aContactNode);
        itemConnection->setData(data,Qt::UserRole);

        //! -------------------
        //! append to the tree
        //! -------------------
        this->myModel->itemFromIndex(myTreeView->currentIndex())->appendRow(itemConnection);
    }
}

//! -------------------------------------
//! function: setTheActiveAnalysisBranch
//! details:
//! -------------------------------------
void SimulationManager::setTheActiveAnalysisBranch()
{
    static QStandardItem *theActiveAnalysis_old;

    //! -----------------------------------
    //! pointer to the standard item model
    //! -----------------------------------
    QModelIndex theModelIndex = myTreeView->currentIndex();
    //QStandardItemModel *theModel = static_cast<QStandardItemModel*>(myTreeView->model());

    QStandardItem *theCurrentItem = myModel->itemFromIndex(theModelIndex);
    SimulationNodeClass *theCurrentNode = theCurrentItem->data(Qt::UserRole).value<SimulationNodeClass*>();

    //! ----------------------
    //! Simulation setup item
    //! "Solution" item
    //! "Analysis settings"
    //! ----------------------
    if(theCurrentNode->isSimulationSetUpNode() || theCurrentNode->isSolution() || theCurrentNode->isAnalysisSettings())
    {
        myActiveAnalysisBranch = theCurrentItem->parent();
        theActiveAnalysis_old = myActiveAnalysisBranch;
        cout<<"@ ------------------------------------------------------------@"<<endl;
        cout<<"@ the current analysis branch is: "<<myActiveAnalysisBranch->data(Qt::UserRole).value<SimulationNodeClass*>()->getName().toStdString()<<endl;
        cout<<"@ ------------------------------------------------------------@"<<endl;
        return;
    }
    //! --------------------------
    //! the item is Analysis root
    //! --------------------------
    if(theCurrentNode->isAnalysisRoot())
    {
        myActiveAnalysisBranch = theCurrentItem;
        theActiveAnalysis_old = myActiveAnalysisBranch;
        cout<<"@ ------------------------------------------------------------@"<<endl;
        cout<<"@ the current analysis branch is: "<<myActiveAnalysisBranch->data(Qt::UserRole).value<SimulationNodeClass*>()->getName().toStdString()<<endl;
        cout<<"@ ------------------------------------------------------------@"<<endl;
        return;
    }
    //! -----------------------
    //! "Solution information"
    //! a post processing item
    //! -----------------------
    if(theCurrentNode->isSolutionInformation() || theCurrentNode->isAnalysisResult())
    {
        myActiveAnalysisBranch = theCurrentItem->parent()->parent();
        theActiveAnalysis_old = myActiveAnalysisBranch;
        cout<<"@ ------------------------------------------------------------@"<<endl;
        cout<<"@ the current analysis branch is: "<<myActiveAnalysisBranch->data(Qt::UserRole).value<SimulationNodeClass*>()->getName().toStdString()<<endl;
        cout<<"@ ------------------------------------------------------------@"<<endl;
        return;
    }

    //! ------------------------------------------------------------
    //! the current item is can not be retrieved from current item:
    //! the previously value is used (which can be NULL)
    //! ------------------------------------------------------------
    myActiveAnalysisBranch = theActiveAnalysis_old;

    if(myActiveAnalysisBranch==NULL)
    {
        cout<<"@ ------------------------------------@"<<endl;
        cout<<"@ the current analysis branch is NULL @"<<endl;
        cout<<"@ ------------------------------------@"<<endl;
    }
    else
    {
        cout<<"@ ------------------------------------------------------------@"<<endl;
        cout<<"@ the current analysis branch is: "<<myActiveAnalysisBranch->data(Qt::UserRole).value<SimulationNodeClass*>()->getName().toStdString()<<endl;
        cout<<"@ ------------------------------------------------------------@"<<endl;
    }
}

//! ----------------------------------
//! function: getActiveAnalysisBranch
//! details:
//! ----------------------------------
QStandardItem* SimulationManager::getActiveAnalysisBranch()
{
    return myActiveAnalysisBranch;
}

//! -------------------------
//! function: exportSTEPFile
//! details:
//! -------------------------
void SimulationManager::exportSTEPFile()
{
    cout<<"SimulationManager::exportSTEPFile()->____function called____"<<endl;

    //! --------------------------------------
    //! prepare stuff for building a compound
    //! --------------------------------------
    BRep_Builder aBuilder;
    TopoDS_Compound aComp;
    aBuilder.MakeCompound(aComp);

    //! --------------------------------------------------------
    //! check if something has been selected within the context
    //! --------------------------------------------------------
    if(!myCTX.IsNull())
    {
        for(myCTX->InitSelected();myCTX->MoreSelected();myCTX->NextSelected())
        {
            TopoDS_Shape aShape = myCTX->SelectedShape();
            if(!aShape.IsNull()) aBuilder.Add(aComp,aShape);
        }
        exportingTools::exportSTEP(aComp);
    }
}

//! -------------------------
//! function: exportSTEPFile
//! details:
//! -------------------------
void SimulationManager::exportBREPFile()
{
    cout<<"SimulationManager::exportBREPFile()->____function called____"<<endl;
}

//! --------------------------------
//! function: deleteAllChildrenItem
//! details:
//! --------------------------------
void SimulationManager::deleteAllChildrenItems()
{
    //! -------------------------
    //! delete all children items
    //! -------------------------
    QStandardItem *itemParent = myModel->itemFromIndex(myTreeView->currentIndex());
    SimulationNodeClass *nodeRoot = itemParent->data(Qt::UserRole).value<SimulationNodeClass*>();
    if(QMessageBox::information(this,APPNAME,"All children items will be removed.\nAre you sure?",QMessageBox::Ok,QMessageBox::Cancel)==QMessageBox::Ok)
    {
        switch(nodeRoot->getType())
        {
        case SimulationNodeClass::nodeType_meshControl:
            myModel->removeRows(0,itemParent->rowCount(),itemParent->index());
            emit requestClearMesh();
            break;

        case SimulationNodeClass::nodeType_coordinateSystems:
        case SimulationNodeClass::nodeType_structuralAnalysis:
        case SimulationNodeClass::nodeType_thermalAnalysis:
            //! -------------------------------------------------------------------
            //! start removing from the end up to row number = 1
            //! (do not remove the "Global coordinate system", "Analysis settins")
            //! Please do not use the built-in QStandardItemModel::removeRows(...)
            //! myModel->removeRows(1,itemParent->rowCount(),itemParent->index());
            //! since not working (Qt bug?)
            //! -------------------------------------------------------------------
            for(int i=itemParent->rowCount(); i!=0; i--) myModel->removeRow(i,itemParent->index());
            break;

        case SimulationNodeClass::nodeType_remotePointRoot:
        case SimulationNodeClass::nodeType_namedSelection:
        case SimulationNodeClass::nodeType_connection:
            //! ----------------------------------------------------------
            //! also in this case stop at the first child, the dummy item
            //! "Select from list", which is hidden
            //! ----------------------------------------------------------
            for(int i=itemParent->rowCount(); i!=0; i--) myModel->removeRow(i,itemParent->index());
            break;

        case SimulationNodeClass::nodeType_connectionGroup:
            myModel->removeRows(0,itemParent->rowCount(),itemParent->index());
            break;

        case SimulationNodeClass::nodeType_importedBodyScalar:
            for(int i = itemParent->rowCount(); i>=0; i--) myModel->removeRow(i,itemParent->index());
            break;
        }
    }
}

//! -------------------------------------------
//! function: generateBoundaryConditionsMeshDS
//! details:
//! -------------------------------------------
void SimulationManager::generateBoundaryConditionsMeshDS(bool computeDual)
{
    cout<<"SimulationManager::generateBoundaryConditionsMeshDS()->____function called____"<<endl;

    //! ----------------------------------------------------------
    //! initialize the following map
    //! key =>   bodyIndex
    //! value => the sub mesh data sources are exact: true/false
    //!
    //! the map value depends on the meshing method applied to
    //! a body: if patch conforming each face/edge/vertes has its
    //! own submesh; if patch independend each face/edge/vertex
    //! should be generated using the proximity method
    //! ----------------------------------------------------------
    QMap<int,bool> mapOfIsMeshDSExact;
    for(QMap<int,TopoDS_Shape>::iterator it = mySimulationDataBase->bodyMap.begin(); it!= mySimulationDataBase->bodyMap.end(); ++it)
    {
        int bodyIndex = it.key();
        mapOfIsMeshDSExact.insert(bodyIndex,true);
    }
    int NbMeshControls = Mesh_RootItem->rowCount();
    for(int k=0; k<NbMeshControls; k++)
    {
        QStandardItem *meshControl = Mesh_RootItem->child(k,0);
        SimulationNodeClass *meshNode = meshControl->data(Qt::UserRole).value<SimulationNodeClass*>();
        if(meshNode->getType()==SimulationNodeClass::nodeType_meshMethod)
        {
            Property::meshEngine2D meshEngine = meshNode->getPropertyValue<Property::meshEngine2D>("Surface mesher");
            switch(meshEngine)
            {
            case Property::meshEngine2D_Netgen:
            case Property::meshEngine2D_Netgen_STL:
            case Property::meshEngine2D_OCC_ExpressMesh:
            {
                const QVector<GeometryTag> &vecLoc = meshNode->getPropertyValue<QVector<GeometryTag>>("Tags");
                for(int k=0; k<vecLoc.size(); k++)
                {
                    int bodyIndex = vecLoc.at(k).parentShapeNr;
                    mapOfIsMeshDSExact.insert(bodyIndex,true);
                }
            }
                break;
            case Property::meshEngine2D_NULL:
            {
                Property::meshEngine3D meshEngine3D = meshNode->getPropertyValue<Property::meshEngine3D>("Volume mesher");
                if(meshEngine3D == Property::meshEngine3D_Tetgen_BR)
                {
                    //! Tetgen with boundary recovery
                    const QVector<GeometryTag> &vecLoc = meshNode->getPropertyValue<QVector<GeometryTag>>("Tags");
                    for(int k=0; k<vecLoc.size(); k++)
                    {
                        int bodyIndex = vecLoc.at(k).parentShapeNr;
                        mapOfIsMeshDSExact.insert(bodyIndex,true);
                    }
                }
                else    // this "else" is unessential, since the map values are initialized as "false"
                {
                    //! TetWild - no boundary recovery. The face mesh datasources
                    //! must be rebuilt
                    const QVector<GeometryTag> &vecLoc = meshNode->getPropertyValue<QVector<GeometryTag>>("Tags");
                    for(int k=0; k<vecLoc.size(); k++)
                    {
                        int bodyIndex = vecLoc.at(k).parentShapeNr;
                        mapOfIsMeshDSExact.insert(bodyIndex,false);
                    }
                }
            }
                break;
            }
        }
    }

    //! ----------------------------
    //! diagnostic - can be removed
    //! ----------------------------
    cout<<"\n\\-----------------------------------------------------\\"<<endl;
    for(QMap<int,bool>::iterator it= mapOfIsMeshDSExact.begin(); it!=mapOfIsMeshDSExact.end(); it++)
        cout<<"\\ body nr: "<< it.key()<<" mesh: "<<(it.value()==true? "Exact":"Approximated")<<endl;
    cout<<"\\-----------------------------------------------------\\\n"<<endl;

    //! --------------------------------------
    //! create a data source builder and init
    //! --------------------------------------
    faceDataSourceBuilder aBuilder;
    aBuilder.setDataBase(mySimulationDataBase);
    aBuilder.setMapOfIsMeshExact(mapOfIsMeshDSExact);

    IndexedMapOfMeshDataSources quelCheResta;
    if(computeDual==true)
    {
        //! ---------------
        //! quel che resta
        //! ---------------
        //IndexedMapOfMeshDataSources quelCheResta;
        for(QMap<int,occHandle(MeshVS_DataSource)>::iterator it = mySimulationDataBase->ArrayOfMeshDS2D.begin();
            it != mySimulationDataBase->ArrayOfMeshDS2D.end(); it++)
        {
            int bodyIndex = it.key();
            const occHandle(Ng_MeshVS_DataSource2D) &surfaceMesh = occHandle(Ng_MeshVS_DataSource2D)::DownCast(it.value());
            const occHandle(Ng_MeshVS_DataSourceFace) &faceMesh = new Ng_MeshVS_DataSourceFace(surfaceMesh);
            quelCheResta.insert(bodyIndex,faceMesh);
        }
    }

    //! --------------------------------------------
    //! scan the rows of the current analysis setup
    //! --------------------------------------------
    QStandardItem *curAnalysisRootItem = myModel->itemFromIndex(myTreeView->currentIndex());
    int NbRows = curAnalysisRootItem->rowCount();

    for(int n=1; n<NbRows-1; n++) //skip the analysis settings item
    {
        QVector<GeometryTag> patchConformingTags;
        QVector<GeometryTag> nonPatchConformingTags;
        //! -------------------
        //! working on an item
        //! -------------------
        QStandardItem *curItem = curAnalysisRootItem->child(n,0);
        SimulationNodeClass *curNode = curItem->data(Qt::UserRole).value<SimulationNodeClass*>();
        cout<<"____"<<curNode->getName().toStdString()<<"____"<<endl;

        SimulationNodeClass::nodeType nodeType = curNode->getType();
        Property::SuppressionStatus isSuppressed = curNode->getPropertyValue<Property::SuppressionStatus>("Suppressed");
        if(isSuppressed == Property::SuppressionStatus_Active)
        {
            if(nodeType == SimulationNodeClass::nodeType_mapper ||
                    nodeType == SimulationNodeClass::nodeType_modelChange ||
                    nodeType == SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_ImportedTemperatureDistribution ||
                    nodeType == SimulationNodeClass::nodeType_structuralAnalysisBoltPretension ||
                    nodeType == SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Acceleration
        #ifdef COSTAMP_VERSION
                    || nodeType == SimulationNodeClass::nodeType_timeStepBuilder
        #endif
                    )
            {
                continue;
            }
            else
            {
                cout<<"____valid BC detected____"<<endl;
                const QVector<GeometryTag> &vecLoc = curNode->getPropertyValue<QVector<GeometryTag>>("Tags");
                for(int i=0; i<vecLoc.size(); i++)
                {
                    int bodyIndex = vecLoc.at(i).parentShapeNr;
                    bool isMeshDSExactOnBody = mapOfIsMeshDSExact.value(bodyIndex);
                    if(isMeshDSExactOnBody) patchConformingTags<<vecLoc.at(i);
                    else nonPatchConformingTags<<vecLoc.at(i);
                }

                //! --------------------------
                //! work on exact datasources
                //! --------------------------
                aBuilder.setFaces(patchConformingTags);
                IndexedMapOfMeshDataSources exactMeshDS;
                aBuilder.perform2(exactMeshDS,true);

                //! ----------------------------
                //! work on inexact datasources
                //! ----------------------------
                aBuilder.setFaces(nonPatchConformingTags);
                IndexedMapOfMeshDataSources inexactMeshDS;
                aBuilder.perform2(inexactMeshDS,false);

                //! --------------------------------------------------------
                //! merge into a single map
                //! details: in case of non patch conforming the DSbuilder
                //! calcualte the exact DS on STL mesh,
                //! KEEP the for on inexactMeshDS after the exact One
                //! --------------------------------------------------------
                IndexedMapOfMeshDataSources finalMapOfMeshDS;
                for(IndexedMapOfMeshDataSources::iterator it = exactMeshDS.begin(); it!=exactMeshDS.end(); it++)
                {
                    int bodyIndex = it.key();
                    occHandle(MeshVS_DataSource) aMeshDS = it.value();
                    finalMapOfMeshDS.insert(bodyIndex,aMeshDS);
                }
                for(IndexedMapOfMeshDataSources::iterator it = inexactMeshDS.begin(); it!=inexactMeshDS.end(); it++)
                {
                    int bodyIndex = it.key();
                    occHandle(MeshVS_DataSource) aMeshDS = it.value();
                    finalMapOfMeshDS.insert(bodyIndex,aMeshDS);
                }

                //! ---------------------------------
                //! put into the simulation database
                //! substitute nrow ... to do ...
                //! ---------------------------------
                QVariant data;
                data.setValue(finalMapOfMeshDS);
                Property prop_meshDataSources("Mesh data sources",data,Property::PropertyGroup_MeshDataSources);
                curNode->removeProperty("Mesh data sources");
                curNode->addProperty(prop_meshDataSources);
                //mySimulationDataBase->MapItemBoudaryConditionDS.insert(n,finalMapOfMeshDS);

                //! ------------
                //! subtraction
                //! ------------
                if(computeDual==true)
                {
                    for(IndexedMapOfMeshDataSources::iterator it=finalMapOfMeshDS.begin(); it!=finalMapOfMeshDS.end(); it++)
                    {
                        int bodyIndex = it.key();
                        occHandle(Ng_MeshVS_DataSourceFace) A = occHandle(Ng_MeshVS_DataSourceFace)::DownCast(quelCheResta.value(bodyIndex));
                        occHandle(Ng_MeshVS_DataSourceFace) B = occHandle(Ng_MeshVS_DataSourceFace)::DownCast(it.value());
                        occHandle(Ng_MeshVS_DataSourceFace) C = new Ng_MeshVS_DataSourceFace(A,B);
                        quelCheResta.insert(bodyIndex,C);
                    }
                }
            }
        }
    }

    cout<<"SimulationManager::generateBoundaryConditionsMeshDS()->____data source created____"<<endl;

    //! ------------------------------------
    //! scan the rows of the contact groups
    //! ------------------------------------
    if(Connections_RootItem->hasChildren())
    {
        int NbContactGroup = Connections_RootItem->rowCount();
        cout<<"SimulationManager::generateBoundaryConditionsMeshDS()->____number of contact groups: "<<NbContactGroup<<"____"<<endl;

        for(int h=0;h<NbContactGroup;h++)
        {
            //! the current connection group
            QStandardItem *itemConnectionGroup = Connections_RootItem->child(h,0);
            //! number of contacts under the current connection group

            if(itemConnectionGroup->hasChildren())
            {
                int NbContactPairs = itemConnectionGroup->rowCount();
                for(int i=0; i<2; i++)
                {
                    for(int n=0; n<NbContactPairs; n++)
                    {
                        QVector<GeometryTag> patchConformingTags;
                        QVector<GeometryTag> nonPatchConformingTags;
                        //! -------------------
                        //! working on an item
                        //! -------------------
                        QStandardItem *curItem = itemConnectionGroup->child(n,0);
                        SimulationNodeClass *curNode = curItem->data(Qt::UserRole).value<SimulationNodeClass*>();

                        Property::SuppressionStatus isSuppressed = curNode->getPropertyItem("Suppressed")->data(Qt::UserRole).value<Property>().getData().value<Property::SuppressionStatus>();
                        if(isSuppressed == Property::SuppressionStatus_Active)
                        {
                            QVector<GeometryTag> vecLoc;
                            if(i==0) //Master
                                vecLoc = curNode->getPropertyValue<QVector<GeometryTag>>("Tags master");
                            else //Slave
                                vecLoc = curNode->getPropertyValue<QVector<GeometryTag>>("Tags slave");

                            for(int ii=0; ii<vecLoc.size(); ii++)
                            {
                                int bodyIndex = vecLoc.at(ii).parentShapeNr;
                                bool isMeshDSExactOnBody = mapOfIsMeshDSExact.value(bodyIndex);
                                if(isMeshDSExactOnBody) patchConformingTags<<vecLoc.at(ii);
                                else nonPatchConformingTags<<vecLoc.at(ii);
                            }
                            //! --------------------------
                            //! work on exact datasources
                            //! --------------------------
                            aBuilder.setFaces(patchConformingTags);
                            IndexedMapOfMeshDataSources exactMeshDS;
                            aBuilder.perform2(exactMeshDS,true);
                            //! ----------------------------
                            //! work on inexact datasources
                            //! ----------------------------
                            aBuilder.setFaces(nonPatchConformingTags);
                            IndexedMapOfMeshDataSources inexactMeshDS;
                            aBuilder.perform2(inexactMeshDS,false);

                            //! ------------------------
                            //! merge into a single map   // check if is it a real merge or a replacement
                            //! ------------------------
                            IndexedMapOfMeshDataSources finalMapOfMeshDS;
                            for(IndexedMapOfMeshDataSources::iterator it = exactMeshDS.begin(); it!=exactMeshDS.end(); it++)
                            {
                                int bodyIndex = it.key();
                                occHandle(MeshVS_DataSource) aMeshDS = it.value();
                                finalMapOfMeshDS.insert(bodyIndex,aMeshDS);
                            }
                            for(IndexedMapOfMeshDataSources::iterator it = inexactMeshDS.begin(); it!=inexactMeshDS.end(); it++)
                            {
                                int bodyIndex = it.key();
                                occHandle(MeshVS_DataSource) aMeshDS = it.value();
                                finalMapOfMeshDS.insert(bodyIndex,aMeshDS);
                            }

                            //! ---------------------------------
                            //! put into the simulation database
                            //! substitute nrow ... to do ...
                            //! ---------------------------------
                            QVariant data;
                            data.setValue(finalMapOfMeshDS);
                            if(i==0)
                            {
                                Property prop_masterMeshDataSources("Master mesh data source",data,Property::PropertyGroup_MeshDataSources);
                                curNode->removeProperty("Master mesh data source");
                                curNode->addProperty(prop_masterMeshDataSources);
                            }
                            else
                            {
                                Property prop_slaveMeshDataSources("Slave mesh data source",data,Property::PropertyGroup_MeshDataSources);
                                curNode->removeProperty("Slave mesh data source");
                                curNode->addProperty(prop_slaveMeshDataSources);
                            }

                            if(computeDual==true)
                            {
                                //! ------------
                                //! subtraction
                                //! ------------
                                for(IndexedMapOfMeshDataSources::iterator it=finalMapOfMeshDS.begin(); it!=finalMapOfMeshDS.end(); it++)
                                {
                                    int bodyIndex = it.key();
                                    occHandle(Ng_MeshVS_DataSourceFace) A = occHandle(Ng_MeshVS_DataSourceFace)::DownCast(quelCheResta.value(bodyIndex));
                                    occHandle(Ng_MeshVS_DataSourceFace) B = occHandle(Ng_MeshVS_DataSourceFace)::DownCast(it.value());
                                    occHandle(Ng_MeshVS_DataSourceFace) C = new Ng_MeshVS_DataSourceFace(A,B);
                                    quelCheResta.insert(bodyIndex,C);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    //! -----------------------------------
    //! scan the rows of the geometry root
    //! -----------------------------------
    for(int n=0; n<Geometry_RootItem->rowCount();n++)
    {
        QVector<GeometryTag> patchConformingTags;
        QVector<GeometryTag> nonPatchConformingTags;
        //! -------------------
        //! working on an item
        //! -------------------
        QStandardItem *curItem = Geometry_RootItem->child(n,0);
        SimulationNodeClass *curNode = curItem->data(Qt::UserRole).value<SimulationNodeClass*>();
        SimulationNodeClass::nodeType nodeType = curNode->getType();
        Property::SuppressionStatus isSuppressed = curNode->getPropertyValue<Property::SuppressionStatus>("Suppressed");
        if(isSuppressed == Property::SuppressionStatus_Active)
        {
            if(nodeType == SimulationNodeClass::nodeType_pointMass)
            {
                const QVector<GeometryTag> &vecLoc = curNode->getPropertyValue<QVector<GeometryTag>>("Tags");
                for(int i=0; i<vecLoc.size(); i++)
                {
                    int bodyIndex = vecLoc.at(i).parentShapeNr;
                    bool isMeshDSExactOnBody = mapOfIsMeshDSExact.value(bodyIndex);
                    if(isMeshDSExactOnBody) patchConformingTags<<vecLoc.at(i);
                    else nonPatchConformingTags<<vecLoc.at(i);
                }
                //! --------------------------
                //! work on exact datasources
                //! --------------------------
                aBuilder.setFaces(patchConformingTags);
                IndexedMapOfMeshDataSources exactMeshDS;
                aBuilder.perform2(exactMeshDS,true);

                //! ----------------------------
                //! work on inexact datasources
                //! ----------------------------
                aBuilder.setFaces(nonPatchConformingTags);
                IndexedMapOfMeshDataSources inexactMeshDS;
                aBuilder.perform2(inexactMeshDS,false);

                //! --------------------------------------------------------
                //! merge into a single map
                //! details: in case of non patch conforming the DSbuilder
                //! calcualte the exact DS on STL mesh,
                //! KEEP the for on inexactMeshDS after the exact One
                //! --------------------------------------------------------
                IndexedMapOfMeshDataSources finalMapOfMeshDS;
                for(IndexedMapOfMeshDataSources::iterator it = exactMeshDS.begin(); it!=exactMeshDS.end(); it++)
                {
                    int bodyIndex = it.key();
                    occHandle(MeshVS_DataSource) aMeshDS = it.value();
                    finalMapOfMeshDS.insert(bodyIndex,aMeshDS);
                }
                for(IndexedMapOfMeshDataSources::iterator it = inexactMeshDS.begin(); it!=inexactMeshDS.end(); it++)
                {
                    int bodyIndex = it.key();
                    occHandle(MeshVS_DataSource) aMeshDS = it.value();
                    finalMapOfMeshDS.insert(bodyIndex,aMeshDS);
                }

                //! ---------------------------------
                //! put into the simulation database
                //! substitute nrow ... to do ...
                //! ---------------------------------
                QVariant data;
                data.setValue(finalMapOfMeshDS);
                Property prop_meshDataSources("Mesh data sources",data,Property::PropertyGroup_MeshDataSources);
                curNode->removeProperty("Mesh data sources");
                curNode->addProperty(prop_meshDataSources);

                //! ------------
                //! subtraction
                //! ------------
                if(computeDual==true)
                {
                    for(IndexedMapOfMeshDataSources::iterator it=finalMapOfMeshDS.begin(); it!=finalMapOfMeshDS.end(); it++)
                    {
                        int bodyIndex = it.key();
                        occHandle(Ng_MeshVS_DataSourceFace) A = occHandle(Ng_MeshVS_DataSourceFace)::DownCast(quelCheResta.value(bodyIndex));
                        occHandle(Ng_MeshVS_DataSourceFace) B = occHandle(Ng_MeshVS_DataSourceFace)::DownCast(it.value());
                        occHandle(Ng_MeshVS_DataSourceFace) C = new Ng_MeshVS_DataSourceFace(A,B);
                        quelCheResta.insert(bodyIndex,C);
                    }
                }
            }
        }
    }
}

//! ---------------------------------------------------------
//! function: replicateBolt
//! details:  this code will be transferred to mainTreeTools
//! ---------------------------------------------------------
#include <BRepGProp.hxx>
#include <GProp_PrincipalProps.hxx>
void SimulationManager::replicateBolt()
{
    FILE *diagFile = fopen("D:\\test.txt","w");

    //! -----------------------
    //! the current model item
    //! -----------------------
    QStandardItemModel *model = static_cast<QStandardItemModel*>(myTreeView->model());
    QExtendedStandardItem *curItem = static_cast<QExtendedStandardItem*>(model->itemFromIndex(myTreeView->currentIndex()));

    //! ----------------------------
    //! the current simulation node
    //! ----------------------------
    SimulationNodeClass *curNode = myTreeView->currentIndex().data(Qt::UserRole).value<SimulationNodeClass*>();
    const QVector<GeometryTag> tags = curNode->getPropertyValue<QVector<GeometryTag>>("Geometry");
    if(tags.isEmpty() || tags.size()>1) return;
    TopAbs_ShapeEnum shapeType = tags.first().subShapeType;
    if(shapeType!=TopAbs_SOLID) return;

    //! ------------------------------------------------
    //! the "reference" solid and its coordinate system
    //! ------------------------------------------------
    int bodyIndex = tags.at(0).parentShapeNr;
    TopoDS_Shape boltShape = mySimulationDataBase->bodyMap.value(bodyIndex);

    //! ---------------------------
    //! compute the center of mass
    //! ---------------------------
    GProp_GProps prop;
    BRepGProp::VolumeProperties(boltShape,prop);
    GProp_PrincipalProps pprop = prop.PrincipalProperties();
    gp_Pnt O = prop.CentreOfMass();
    double xO = O.X(); double yO = O.Y(); double zO = O.Z();

    cout<<"____center of mass ("<<xO<<", "<<yO<<", "<<zO<<")____"<<endl;

    gp_Vec a1 = pprop.FirstAxisOfInertia();
    gp_Vec a2 = pprop.SecondAxisOfInertia();
    gp_Vec a3 = pprop.ThirdAxisOfInertia();
    a1.Normalize(); a2.Normalize(); a3.Normalize();

    double uXx = a1.X(); double uXy = a1.Y(); double uXz= a1.Z();
    double uYx = a2.X(); double uYy = a2.Y(); double uYz= a2.Z();
    double uZx = a3.X(); double uZy = a3.Y(); double uZz= a3.Z();

    //! --------------------------------------
    //! retrieve the "bolt" coordinate system
    //! (reference point, direction)
    //! --------------------------------------
    void *p = curNode->getPropertyValue<void*>("Coordinate system");
    QStandardItem *itemCS = static_cast<QStandardItem*>(p);
    SimulationNodeClass *nodeCS = itemCS->data(Qt::UserRole).value<SimulationNodeClass*>();
    QVector<double> boltDirection;
    QVector<double> planeOrigin;
    if(nodeCS->getType()==SimulationNodeClass::nodeType_coordinateSystem_global)
    {
        planeOrigin.push_back(0); planeOrigin.push_back(0); planeOrigin.push_back(0);
        boltDirection = nodeCS->getPropertyValue<QVector<double>>("Z axis data");
    }
    else
    {
        planeOrigin = nodeCS->getPropertyValue<QVector<double>>("Base origin");
        QVector<QVector<double>> baseDirData = nodeCS->getPropertyValue<QVector<QVector<double>>>("Base directional data");
        boltDirection = baseDirData.at(2);
    }
    //! -------------------------------------------
    //! coordinates of the reference point as read
    //! from the interface
    //! -------------------------------------------
    double xRP = planeOrigin[0];
    double yRP = planeOrigin[1];
    double zRP = planeOrigin[2];

    //! -----------------------------------
    //! position of the reference point
    //! into the local system of reference
    //! -----------------------------------
    double XRP = (xRP-xO)*uXx+(yRP-yO)*uXy+(zRP-zO)*uXz;
    double YRP = (xRP-xO)*uYx+(yRP-yO)*uYy+(zRP-zO)*uYz;
    double ZRP = (xRP-xO)*uZx+(yRP-yO)*uZy+(zRP-zO)*uZz;

    fprintf(diagFile,"%lf\t%lf\t%lf\n",XRP,YRP,ZRP);

    for(myCTX->InitSelected();myCTX->MoreSelected(); myCTX->NextSelected())
    {
        TopoDS_Shape curSelectedShape = myCTX->SelectedShape();
        TopLoc_Location curLoc = curSelectedShape.Location();
        TopoDS_Shape aShape = boltShape.Located(curLoc);
        if(!aShape.IsEqual(boltShape)) continue;
        cout<<"____found identical____"<<endl;

        //! ---------------------------
        //! compute the center of mass
        //! ---------------------------
        GProp_GProps prop_;
        BRepGProp::VolumeProperties(curSelectedShape,prop_);
        gp_Pnt O_ = prop_.CentreOfMass();
        GProp_PrincipalProps pprop_ = prop_.PrincipalProperties();
        double xO_ = O_.X(); double yO_ = O_.Y(); double zO_ = O_.Z();
        cout<<"____center of mass ("<<xO_<<", "<<yO_<<", "<<zO_<<")____"<<endl;

        gp_Vec a1_ = pprop_.FirstAxisOfInertia();
        gp_Vec a2_ = pprop_.SecondAxisOfInertia();
        gp_Vec a3_ = pprop_.ThirdAxisOfInertia();
        a1_.Normalize(); a2_.Normalize(); a3_.Normalize();

        double uXx_ = a1_.X(); double uXy_ = a1_.Y(); double uXz_= a1_.Z();
        double uYx_ = a2_.X(); double uYy_ = a2_.Y(); double uYz_= a2_.Z();
        double uZx_ = a3_.X(); double uZy_ = a3_.Y(); double uZz_= a3_.Z();

        //! ---------------------------------------
        //! parameters for a new coordinate system
        //! "Base origin", "Base directional data"
        //! ---------------------------------------
        double xRP_ = xO_+XRP*uXx_+YRP*uYx_+ZRP*uZx_;
        double yRP_ = yO_+XRP*uXy_+YRP*uYy_+ZRP*uZy_;
        double zRP_ = zO_+XRP*uXz_+YRP*uYz_+ZRP*uZz_;

        fprintf(diagFile,"%lf\t%lf\t%lf\n",xRP_,yRP_,zRP_);

        QVector<double> origin_{xRP_,yRP_,zRP_};
        QVariant data;
        data.setValue(origin_);
        Property prop_baseOrigin("Base origin",data,Property::PropertyGroup_Transformations);

        QVector<double> xAxisData{uXx_,uXy_,uXz_};
        QVector<double> yAxisData{uYx_,uYy_,uYz_};
        QVector<double> zAxisData{uZx_,uZy_,uZz_};
        QVector<QVector<double>> axisData{xAxisData,yAxisData,zAxisData};
        data.setValue(axisData);
        Property prop_axisData("Base directional data",data,Property::PropertyGroup_Transformations);

        //! -------------------------------
        //! build a coordinate system item
        //! -------------------------------
        QExtendedStandardItem *itemCSRoot = this->getTreeItem(SimulationNodeClass::nodeType_coordinateSystems);
        int NbCoordinateSystems = itemCSRoot->rowCount()-1;
        QString coordinateSystemNodeName = QString("Coordinate system %1").arg(++NbCoordinateSystems);
        SimulationNodeClass *nodeCS = nodeFactory::nodeFromScratch(SimulationNodeClass::nodeType_coordinateSystem);

        //! -----------------------
        //! block new node signals
        //! -----------------------
        //nodeCS->getModel()->blockSignals(true);

        nodeCS->setName(coordinateSystemNodeName);
        nodeCS->replaceProperty("Base origin",prop_baseOrigin);
        nodeCS->replaceProperty("Base directional data",prop_axisData);

        QExtendedStandardItem *itemCS = new QExtendedStandardItem();
        data.setValue(nodeCS);
        itemCS->setData(data,Qt::UserRole);
        data.setValue(coordinateSystemNodeName);
        itemCS->setData(data,Qt::DisplayRole);

        //! -------------------------
        //! unblock new node signals
        //! -------------------------
        //nodeCS->getModel()->blockSignals(false);

        //! -----------------------------
        //! append new coordinate system
        //! -----------------------------
        QList<QStandardItem*> itemList;
        itemList<<itemCS;
        itemCSRoot->appendRow(itemList);
        //myTreeView->setCurrentIndex(itemCSRoot->index());

        //! -------------------------------------------------------
        //! programmatically add the new bolt: use "duplicateItem"
        //! -------------------------------------------------------
        this->duplicateItem(curItem);
        SimulationNodeClass *boltNode = curItem->data(Qt::UserRole).value<SimulationNodeClass*>();

        //! -----------------------
        //! block new node signals
        //! -----------------------
        //boltNode->getModel()->blockSignals(true);

        //! -----------------------------------
        //! replace the default scope and tags
        //! -----------------------------------
        ListOfShape listOfShape; listOfShape.Append(curSelectedShape);
        QVector<GeometryTag> tags = TopologyTools::generateLocationPairs(mySimulationDataBase,listOfShape);
        data.setValue(tags);
        Property prop_scope("Geometry",data,Property::PropertyGroup_Scope);
        Property prop_tags("Tags",data,Property::PropertyGroup_Scope);
        boltNode->replaceProperty("Geometry",prop_scope);
        boltNode->replaceProperty("Tags",prop_tags);

        //! --------------------------------------
        //! replace the default coordinate system
        //! --------------------------------------
        void *p = (void*)itemCS;
        data.setValue(p);
        Property prop_coordinateSystem("Coordinate system",data,Property::PropertyGroup_Definition);
        boltNode->replaceProperty("Coordinate system",prop_coordinateSystem);

        //! ------------------------------
        //! replace the "Reference point"
        //! ------------------------------
        data.setValue(origin_);
        Property prop_referencePoint("Reference point",data,Property::PropertyGroup_Hidden);
        boltNode->replaceProperty("Reference point",prop_referencePoint);

        //! -------------------------
        //! unblock new node signals
        //! -------------------------
        //boltNode->getModel()->blockSignals(true);
    }
    fclose(diagFile);
}

//! -------------------------------
//! function: invertSuppressionSet
//! details:
//! -------------------------------
void SimulationManager::invertSuppressionSet()
{
    cout<<"SimulationManager::invertSuppressionSet()->____function called____"<<endl;

    //! --------------------------------------------------------------------
    //! handle a selection from the OCC viewer: by definition the selection
    //! contains "active" (since visible) TopoDS_Shape(s)
    //! --------------------------------------------------------------------
    myCTX->InitSelected();
    if(myCTX->MoreSelected())
    {
        //! ---------------------
        //! list of body numbers
        //! ---------------------
        TColStd_ListOfInteger ListOfActiveBodyNumbers;
        TColStd_ListOfInteger ListOfSuppressedBodyNumbers;

        //! ----------------------------------------------------------------
        //! the list of all the AIS_Shapes currently present in the context
        //! ----------------------------------------------------------------
        AIS_ListOfInteractive listOfAISShapes;
        myCTX->ObjectsInside(listOfAISShapes,AIS_KOI_Shape,0);

        //! --------------------------------------------------------
        //! retrieve the body indices of the active, visible bodies
        //! meanwhile collect the indices of the suppressed bodies
        //! --------------------------------------------------------
        TopTools_ListOfShape aListOfShape;
        for(myCTX->InitSelected(); myCTX->MoreSelected(); myCTX->NextSelected())
        {
            const occHandle(AIS_Shape) &theAISShape = occHandle(AIS_Shape)::DownCast(myCTX->SelectedInteractive());
            const TopoDS_Shape &theShape = theAISShape->Shape();
            aListOfShape.Append(theShape);
            int bodyIndex = mySimulationDataBase->bodyMap.key(theShape);
            ListOfActiveBodyNumbers.Append(bodyIndex);

            //! ---------------------------------------------------------
            //! remove the visible shape from the list of all the shapes
            //! at the end the list of all the shapes will contains also
            //! suppressed non visible shapes
            //! ---------------------------------------------------------
            listOfAISShapes.Remove(theAISShape);
        }

        for(AIS_ListIteratorOfListOfInteractive it(listOfAISShapes); it.More(); it.Next())
        {
            const occHandle(AIS_Shape) &theAISShape = occHandle(AIS_Shape)::DownCast(it.Value());
            const TopoDS_Shape &theShape = theAISShape->Shape();
            int suppressedBodyIndex = mySimulationDataBase->bodyMap.key(theShape);
            ListOfSuppressedBodyNumbers.Append(suppressedBodyIndex);
        }

        //! -----------------
        //! test the results
        //! -----------------
        cout<<"----CURRENT ACTIVE BODIES----"<<endl;
        for(TColStd_ListOfInteger::iterator it = ListOfActiveBodyNumbers.begin(); it!=ListOfActiveBodyNumbers.end(); ++it)
        {
            cout<<"____currently active: "<<*it<<"____"<<endl;
        }
        cout<<"----CURRENT SUPPRESSED BODIES----"<<endl;
        for(TColStd_ListOfInteger::iterator it = ListOfSuppressedBodyNumbers.begin(); it!=ListOfSuppressedBodyNumbers.end(); ++it)
        {
            cout<<"____currently suppressed: "<<*it<<"____"<<endl;
        }
    }
    /*
        //! ------------------------------------------------------------------------
        //! scan the main tree and update the properties "Suppressed" and "Visible"
        //! meanwhile fill the list of the currently suppressed bodies
        //! ------------------------------------------------------------------------
        QExtendedStandardItem *itemGeometry = this->getTreeItem(SimulationNodeClass::nodeType_geometry);
        for(int i=0; i<itemGeometry->rowCount(); i++)
        {
            QExtendedStandardItem *itemBody = static_cast<QExtendedStandardItem*>(itemGeometry->child(i,0));
            SimulationNodeClass *nodeBody = itemBody->data(Qt::UserRole).value<SimulationNodeClass*>();

            //! -----------------------
            //! block the node signals
            //! -----------------------
            nodeBody->getModel()->blockSignals(true);

            int bodyIndex = nodeBody->getPropertyValue<int>("Map index");
            if(ListOfBodyNumbers.Contains(bodyIndex))
            {
                //! ----------------------------------------------------
                //! synch the "Suppressed" property and update the icon
                //! ----------------------------------------------------
                nodeBody->replaceProperty("Suppressed",prop_newSuppressionStatus);
                if(isSuppressed) itemBody->setIcon(QIcon(":/icons/icon_suppress.png"));
                else itemBody->setIcon(itemBody->getIcon(nodeBody->getType()));

                //! -----------------------------------
                //! synch also the visibility property
                //! -----------------------------------
                bool isVisible = newSuppressionStatus==Property::SuppressionStatus_Active? true:false;
                data.setValue(isVisible);
                Property prop_visible("Visible",data,Property::PropertyGroup_GraphicProperties);
                nodeBody->replaceProperty("Visible",prop_visible);
            }
            else
            {

            }

            //! ------------------------
            //! unlock the node signals
            //! ------------------------
            nodeBody->getModel()->blockSignals(false);
        }
    }
*/
}

//! --------------------------------------
//! function: computeAndDisplayMeshMetric
//! details:
//! --------------------------------------
#include <qhistogramdata.h>
#include <qhistogram.h>
#include <tetqualityclass.h>
#include <elementtypes.h>
void SimulationManager::computeAndDisplayMeshMetric()
{
    cout<<"SimulationManager::computeAndDisplayMeshMetric()->____function called____"<<endl;

    //! -------------------------
    //! preliminary sanity check
    //! -------------------------
    QModelIndex modelIndex = myTreeView->currentIndex();
    if(!modelIndex.isValid()) return;

    SimulationNodeClass *node = modelIndex.data(Qt::UserRole).value<SimulationNodeClass*>();
    if(node->getType()!=SimulationNodeClass::nodeType_meshMeshMetric) return;

    QVector<GeometryTag> vecLoc = node->getPropertyValue<QVector<GeometryTag>>("Tags");
    if(vecLoc.isEmpty()) return;
    if(vecLoc.at(0).subShapeType!=TopAbs_SOLID) return;

    //! ---------------------------------------
    //! build the bins and initialize the data
    //! ---------------------------------------
    histogramData hData;    // empty data container
    QVariant data;

    //! ----------------------------------
    //! compute all the mesh quality data
    //! ----------------------------------
    int metricType = node->getPropertyValue<int>("Metric type");

    //! ---------------------------------
    //! vector of {q0,q1,q2,q3,V} values
    //! ---------------------------------
    std::vector<std::vector<double>> qualityValues;

    //! -----------------------
    //! scan the volume meshes
    //! -----------------------
    //TetQualityClass aTetQualityChecker;
    for(QVector<GeometryTag>::iterator it = vecLoc.begin(); it!=vecLoc.end(); it++)
    {
        int bodyIndex = it->parentShapeNr;

        const occHandle(Ng_MeshVS_DataSource3D) &aVolumeMesh = occHandle(Ng_MeshVS_DataSource3D)::DownCast(mySimulationDataBase->ArrayOfMeshDS.value(bodyIndex));
        if(aVolumeMesh.IsNull()) continue;

        //! ---------------------------------------------
        //! iterate over the volume elements of the mesh
        //! ---------------------------------------------
        for(TColStd_MapIteratorOfPackedMapOfInteger eIt(aVolumeMesh->GetAllElements()); eIt.More(); eIt.Next())
        {
            int NbNodes, buf[4];
            TColStd_Array1OfInteger nodeIDs(*buf,1,4);
            int globalElementID = eIt.Key();

            //! -------------
            //! limit to TET
            //! -------------
            ElemType eType;
            if(!aVolumeMesh->GetElementType(eType,globalElementID,false)) continue;
            if(eType!=TET) continue;

            aVolumeMesh->GetNodesByElement(globalElementID,nodeIDs,NbNodes);
            std::vector<mesh::meshPoint> aTet;
            for(int n=1; n<=NbNodes; n++)
            {
                int localNodeID = aVolumeMesh->myNodesMap.FindIndex(nodeIDs(n));
                const std::vector<double> &c = aVolumeMesh->getNodeCoordinates(localNodeID);
                aTet.push_back(mesh::meshPoint(c[0],c[1],c[2]));
            }

            //! --------------------------------
            //! do quality computation on a tet
            //! --------------------------------
            TetQualityClass aTetQualityChecker(aTet);
            //aTetQualityChecker.setPoints(aTet);
            double q0,q1,q2,V;
            aTetQualityChecker.getQualityMeasure(q0,q1,q2,V);

            std::vector<double> qualityParameters;
            qualityParameters.push_back(q0);
            qualityParameters.push_back(q1);
            qualityParameters.push_back(q2);
            qualityParameters.push_back(V);

            //cout<<"____("<<q0<<", "<<q1<<", "<<q2<<", "<<V<<")____"<<endl;
            qualityValues.push_back(qualityParameters);
        }
    }

    cout<<"____quality computation for the mesh done____"<<endl;

    //! ----------------------------
    //! array of values into vector
    //! ----------------------------
    cout<<"____metric type: "<<metricType<<"____"<<endl;
    if(metricType == 3)
    {
        data.setValue(hData);
        Property prop_meshMetricData("Metric data",data,Property::PropertyGroup_Definition);
        node->replaceProperty("Metric data",prop_meshMetricData);
        emit requestDisplayMeshMetricHystogram(hData);
        return;
    }

    //! ----------------------------------------------------
    //! values for a type of metric among the available one
    //! ----------------------------------------------------
    std::vector<double> valuesForMetric;
    for(std::vector<std::vector<double>>::iterator it = qualityValues.begin(); it!= qualityValues.end(); it++)
    {
        std::vector<double> metricValues = *it;
        double val = metricValues[metricType];
        valuesForMetric.push_back(val);
    }

    //! ---------------------------------------------------------------------
    //! min & max of the quality values for the specific metric "metricType"
    //! ---------------------------------------------------------------------
    double xMin = *(std::min_element(valuesForMetric.begin(),valuesForMetric.end()));
    double xMax = *(std::max_element(valuesForMetric.begin(),valuesForMetric.end()));
    //cout<<"____("<<xMin<<", "<<xMax<<")____"<<endl;

    const int NBINS = 25;
    double dx = (xMax-xMin)/double(NBINS);

    //! -----------------------
    //! invalid hystogram data
    //! -----------------------
    if(dx == 0)
    {
        data.setValue(hData);
        Property prop_meshMetricData("Metric data",data,Property::PropertyGroup_Definition);
        node->replaceProperty("Metric data",prop_meshMetricData);
        emit requestDisplayMeshMetricHystogram(hData);
        return;
    }

    //! -----------------
    //! pile-up the bins
    //! -----------------
    std::vector<hbin> vecBin;
    for(int n=0; n<NBINS; n++)
    {
        double x_left = xMin + n*dx;
        double x_right = x_left+dx;
        hbin aBin(x_left,x_right);
        vecBin.push_back(aBin);
        std::pair<hbin,double> element;
        element.first = aBin;
        element.second = 0;
        hData.insert(aBin,0);
    }

    //! -------------------------
    //! build the hystogram data
    //! -------------------------
    std::map<hbin,double> hysData;
    hData.getData(hysData);

    for(std::vector<double>::iterator it = valuesForMetric.begin(); it!=valuesForMetric.end(); it++)
    {
        double curMetricVal = *it;
        int binNb = floor((curMetricVal-xMin)/dx);
        if(binNb==NBINS) binNb--;
        const hbin &curBin = vecBin[binNb];
        std::map<hbin,double>::iterator it_ = hysData.find(curBin);
        (*it_).second++;
    }

    hData.setData(hysData);

    //! --------------------------------------------------
    //! replace property and display the mesh metric data
    //! --------------------------------------------------
    data.setValue(hysData);
    Property prop_meshMetricData("Metric data",data,Property::PropertyGroup_Definition);
    node->replaceProperty("Metric data",prop_meshMetricData);
    emit requestDisplayMeshMetricHystogram(hysData);
}

//! -------------------------------------
//! function: deleteDataSourcesFromModel
//! details:
//! -------------------------------------
void SimulationManager::deleteDataSourcesFromModel()
{
    cout<<"SimulationManager::deleteDataSourcesFromModel()->____function called____"<<endl;

    //! ------------------------------------------
    //! delete the data sources from the treeview
    //! ------------------------------------------
    QList<QExtendedStandardItem*> items;
    this->getTreeItemsRecursively(myModel,items);
    int NbDeleted = 0;
    int NbItems = items.length();
    for(int n=0; n<NbItems; n++)
    {
        SimulationNodeClass *curNode = items[n]->data(Qt::UserRole).value<SimulationNodeClass*>();
        if(curNode->isSimulationSetUpNode())
        {
            bool isDone = curNode->removeProperty("Mesh data sources");
            if(isDone)
            {
                NbDeleted++;
            }
        }
        if(curNode->getType()==SimulationNodeClass::nodeType_connectionPair)
        {
            bool isDone = curNode->removeProperty("Master mesh data source");
            if(isDone) NbDeleted++;
            isDone = curNode->removeProperty("Slave mesh data source");
            NbDeleted++;
        }
    }
    cout<<"SimulationManager::deleteMeshDataSourcesFromTree()->____"<<NbDeleted<<" mesh data sources____"<<endl;
}

