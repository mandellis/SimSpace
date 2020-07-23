//! ----------------
//! custom includes
//! ----------------
#include "contextmenubuilder.h"
#include "simulationnodeclass.h"
#include "qextendedstandarditem.h"
#include "postobject.h"
#include "ccout.h"
#include "tools.h"
#include "mainwindow.h"
#include "simulationmanager.h"

//! ---
//! Qt
//! ---
#include <QIcon>
#include <QAction>

//! ----
//! OCC
//! ----
#include <AIS_Shape.hxx>
#include <TopAbs_ShapeEnum.hxx>

//! -----------------------------
//! function: buildCommonActions
//! details:
//! -----------------------------
void contextMenuBuilder::buildCommonActions(QMenu *contextMenu, bool isEnabled)
{
    //! ---------------------------------------------
    //! preliminary: retrieve the simulation manager
    //! ---------------------------------------------
    SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
    QModelIndex curIndex = sm->myTreeView->currentIndex();
    SimulationNodeClass *node = curIndex.data(Qt::UserRole).value<SimulationNodeClass*>();

    //! -----------------
    //! action duplicate
    //! -----------------
    contextMenuBuilder::addActionDuplicate(contextMenu);

    //! add separator
    contextMenu->addSeparator();

    //! --------------------------
    //! suppression/unsuppression
    //! --------------------------
    if(node->getPropertyItem("Suppressed")!=Q_NULLPTR)
    {
        QList<QModelIndex> selectedIndexes = sm->myTreeView->selectionModel()->selectedIndexes();
        if(selectedIndexes.length()==1)
        {
            Property::SuppressionStatus ss = node->getPropertyValue<Property::SuppressionStatus>("Suppressed");
            if(ss==Property::SuppressionStatus_Suppressed) contextMenuBuilder::addActionUnsuppress(contextMenu);
            else contextMenuBuilder::addActionSuppress(contextMenu);
        }
        else
        {
            contextMenuBuilder::addActionSuppress(contextMenu);
            contextMenuBuilder::addActionUnsuppress(contextMenu);
        }
    }

    //! add separator
    contextMenu->addSeparator();

    //! --------------
    //! delete action
    //! --------------
    contextMenuBuilder::addActionDelete(contextMenu);
    contextMenuBuilder::addActionDeleteAllChildrenItems(contextMenu);

    //! add a separator
    contextMenu->addSeparator();

    //! ---------------------------------------------
    //! add action rename/rename based on definition
    //! ---------------------------------------------
    contextMenuBuilder::addActionRename(contextMenu);
    if(node->getType()!=node->getFamily()) contextMenuBuilder::addActionRenameBasedOnDefinition(contextMenu);
}

//! -------------------------------
//! function: buildMeshContextMenu
//! details:
//! -------------------------------
void contextMenuBuilder::buildMeshContextMenu(QMenu *contextMenu, bool addCommonActions, bool isEnabled)
{    
    if(!isEnabled) return;

    //! ---------------------------------------------
    //! preliminary: retrieve the simulation manager
    //! ---------------------------------------------
    SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
    SimulationNodeClass *node = sm->myTreeView->currentIndex().data(Qt::UserRole).value<SimulationNodeClass*>();

    //! ---------------
    //! submenu insert
    //! ---------------
    QMenu *menuInsert = contextMenu->addMenu("Insert");
    menuInsert->setIcon(QIcon(":/icons/icon_insert.png"));

    //! ------------
    //! mesh method
    //! ------------
    QAction *ActionInsertMeshMethod = menuInsert->addAction("Method");
    ActionInsertMeshMethod->setIcon(QIcon(":/icons/icon_mesh method.png"));
    ActionInsertMeshMethod->setData(74);

    //! add separator
    menuInsert->addSeparator();

    //! -------------
    //! Type of mesh
    //! -------------
    QAction *ActionInsertTypeOfMesh = menuInsert->addAction("Type of mesh");
    ActionInsertTypeOfMesh->setIcon(QIcon(":/icons/icon_type of mesh.png"));
    ActionInsertTypeOfMesh->setData(79);

    //! add separator
    menuInsert->addSeparator();

    //! ------------
    //! body sizing
    //! ------------
    QAction *ActionInsertMeshBodyControl = menuInsert->addAction("Body sizing");
    ActionInsertMeshBodyControl->setIcon(QIcon(":/icons/icon_volume mesh.png"));
    ActionInsertMeshBodyControl->setData(51);

    //! ------------
    //! face sizing
    //! ------------
    QAction *ActionInsertFaceSizing = menuInsert->addAction("Face sizing");
    ActionInsertFaceSizing->setIcon(QIcon(":/icons/icon_mesh face sizing.png"));
    ActionInsertFaceSizing->setData(72);

    //! ------------
    //! edge sizing
    //! ------------
    QAction *ActionInsertEdgeSizing = menuInsert->addAction("Edge sizing");
    ActionInsertEdgeSizing->setIcon(QIcon(":/icons/icon_edge sizing.png"));
    ActionInsertEdgeSizing->setData(68);

    //! --------------
    //! vertex sizing
    //! --------------
    QAction *ActionInsertVertexSizing = menuInsert->addAction("Vertex sizing");
    ActionInsertVertexSizing->setIcon(QIcon(":/icons/icon_point.png"));
    ActionInsertVertexSizing->setData(69);

    //! add separator
    menuInsert->addSeparator();

    //! ----------------
    //! prismatic layer
    //! ----------------
    QAction *ActionInsertPrismaticLayer = menuInsert->addAction("Prismatic layer");
    ActionInsertPrismaticLayer->setIcon(QIcon(":/icons/icon_prismatic layer.png"));
    ActionInsertPrismaticLayer->setData(75);

    //! add a separator
    menuInsert->addSeparator();

    //! ------------
    //! mesh metric
    //! ------------
    QAction* ActionInsertMeshMetric = menuInsert->addAction("Mesh metric");
    ActionInsertMeshMetric->setIcon(QIcon(":/icons/icon_metric.png"));
    ActionInsertMeshMetric->setData(89);

    //! add a separator
    contextMenu->addSeparator();

    //! -------------------------------------------------------------
    //! generate the surface mesh - this option sould be deactivated
    //! when the volume mesher is set "to TetgenBR"
    //! -------------------------------------------------------------
    QAction *ActionGenerateSurfaceMesh;
    ActionGenerateSurfaceMesh = contextMenu->addAction("Preview mesh");
    ActionGenerateSurfaceMesh->setIcon(QIcon(":/icons/icon_solve.png"));
    ActionGenerateSurfaceMesh->setData(55);

    //! ------------------
    //! generate the mesh
    //! ------------------
    QAction *ActionGenerateVolumeMesh;
    ActionGenerateVolumeMesh = contextMenu->addAction("Generate mesh");
    ActionGenerateVolumeMesh->setIcon(QIcon(":/icons/icon_solve.png"));
    ActionGenerateVolumeMesh->setData(54);

    //! ------------------
    //! preview inflation
    //! ------------------
    if(sm->getCurrentNode()->getType()==SimulationNodeClass::nodeType_meshPrismaticLayer)
    {
        QAction *ActionPreviewInflation = contextMenu->addAction("Preview inflation");
        ActionPreviewInflation->setIcon(QIcon(":/icons/icon_prismatic layer.png"));
        ActionPreviewInflation->setData(67);
    }

    if(addCommonActions==true)
    {
        //! add a separator
        contextMenu->addSeparator();

        //! a mesh control is selected
        if(node->getType()!=node->getFamily()) contextMenuBuilder::buildCommonActions(contextMenu);
    }

    //! add a separator
    contextMenu->addSeparator();

    //! ------------------------------------------------------
    //! clear the mesh - intentionally at the end of the menu
    //! ------------------------------------------------------
    QAction *ActionClearGeneratedData = contextMenu->addAction("Clear all meshes");
    ActionClearGeneratedData->setIcon(QIcon(":/icons/icon_clear data.png"));
    ActionClearGeneratedData->setData(56);
}

//! ----------------------------------
//! function: buildContactContextMenu
//! details:
//! ----------------------------------
void contextMenuBuilder::buildContactContextMenu(QMenu *contextMenu, bool addCommonActions, bool isEnabled)
{
    if(!isEnabled) return;

    //! ---------------------------------------------
    //! preliminaty: retrieve the simulation manager
    //! ---------------------------------------------
    SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
    SimulationNodeClass *node = sm->myTreeView->currentIndex().data(Qt::UserRole).value<SimulationNodeClass*>();
    SimulationNodeClass::nodeType nodeType = node->getType();

    //! ---------------
    //! submenu insert
    //! ---------------
    QMenu *menuInsert = contextMenu->addMenu("Insert");
    menuInsert->setIcon(QIcon(":/icons/icon_insert.png"));

    //! ----------------------
    //! action insert contact
    //! ----------------------
    QAction* ActionInsertContact = menuInsert->addAction("Manual contact");
    ActionInsertContact->setIcon(QIcon(":/icons/icon_add contact.png"));
    ActionInsertContact->setData(49);

    //! ----------------------------
    //! action insert contact group
    //! ----------------------------
    QAction* ActionInsertContactGroup = menuInsert->addAction("Contact group");
    ActionInsertContactGroup->setIcon(QIcon(":/icons/icon_folder.png"));
    ActionInsertContactGroup->setData(44);

    //! add separator
    contextMenu->addSeparator();

    //! --------------------------------------
    //! action generate automatic connections
    //! --------------------------------------
    if(nodeType==SimulationNodeClass::nodeType_connectionGroup)
    {
        QAction *ActionGenerateAutomaticContact = contextMenu->addAction("Create automatic connections");
        ActionGenerateAutomaticContact->setIcon(QIcon(":/icons/icon_solve.png"));
        ActionGenerateAutomaticContact->setData(43);
    }

    //! -------------------------------------------------------
    //! a contact pair or a connection group has been selected
    //! -------------------------------------------------------
    if(nodeType == SimulationNodeClass::nodeType_connectionPair || nodeType == SimulationNodeClass::nodeType_connectionGroup)
    {
        //! --------------------
        //! action flip contact
        //! --------------------
        if(node->getType()==SimulationNodeClass::nodeType_connectionPair)
        {
            QAction *ActionFlipContact = contextMenu->addAction("Swap master slave");
            ActionFlipContact->setIcon(QIcon(":/icons/icon_flip contact.png"));
            ActionFlipContact->setData(48);
        }

        //! add separator
        contextMenu->addSeparator();

        //! -----------------------------------------------------
        //! check if all the selected items are connection pairs
        //! -----------------------------------------------------
        bool allItemsContactPairs = true;
        for(int i=0; i<sm->myTreeView->selectionModel()->selectedIndexes().size(); i++)
        {
            if(sm->myTreeView->selectionModel()->selectedIndexes().at(i).data(Qt::UserRole).value<SimulationNodeClass*>()->getType()!=
                    SimulationNodeClass::nodeType_connectionPair)
            {
                allItemsContactPairs = false;
                break;
            }
        }
        if(allItemsContactPairs && sm->myTreeView->selectionModel()->selectedIndexes().length()>1)
        {
            QAction *ActionMergeSelectedContacts = contextMenu->addAction("Merge selected contacts");
            ActionMergeSelectedContacts->setIcon(QIcon(":/icons/icon_merge.png"));
            ActionMergeSelectedContacts->setData(50);
        }

        //! add separator
        contextMenu->addSeparator();

        //! -----------------------
        //! add the common actions
        //! -----------------------
        if(addCommonActions==true) contextMenuBuilder::buildCommonActions(contextMenu,isEnabled);
    }
}

//! ------------------------------------
//! function: buildCoordinateSystemMenu
//! details:
//! ------------------------------------
void contextMenuBuilder::buildCoordinateSystemMenu(QMenu *contextMenu, bool addCommonActions, bool isEnabled)
{
    if(!isEnabled) return;

    //! ---------------------------------------------
    //! preliminary: retrieve the simulation manager
    //! ---------------------------------------------
    SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
    SimulationNodeClass *node = sm->myTreeView->currentIndex().data(Qt::UserRole).value<SimulationNodeClass*>();
    SimulationNodeClass::nodeType nodeType = node->getType();

    //! ------------
    //! menu insert
    //! ------------
    QMenu *menuInsert = contextMenu->addMenu("Insert");
    menuInsert->setIcon(QIcon(":/icons/icon_insert.png"));

    //! --------------------------------
    //! action insert coordinate system
    //! --------------------------------
    QAction *ActionInsertCoordinateSystem = menuInsert->addAction("Coordinate system");
    ActionInsertCoordinateSystem->setIcon(QIcon(":/icons/icon_system of reference.png"));
    ActionInsertCoordinateSystem->setData(47);

    //! add separator
    contextMenu->addSeparator();

    //! ---------------------------------
    //! action delete all children items
    //! ---------------------------------
    if(node->getType()==SimulationNodeClass::nodeType_coordinateSystems)
    {
        contextMenuBuilder::addActionDeleteAllChildrenItems(contextMenu);
    }

    //! --------------------------------------------------------------
    //! the global coordinate system cannot be removed nor duplicated
    //! --------------------------------------------------------------
    if(nodeType==SimulationNodeClass::nodeType_coordinateSystem)
    {
        if(addCommonActions==true)
        {
            contextMenu->addSeparator();
            contextMenuBuilder::buildCommonActions(contextMenu,isEnabled);
        }
    }
}

//! ----------------------------------------------------------------
//! function: buildRemotePointContextMenu
//! details:  once inserted the Remote point root cannot be removed
//! ----------------------------------------------------------------
void contextMenuBuilder::buildRemotePointContextMenu(QMenu *contextMenu, bool addCommonActions, bool isEnabled)
{
    if(isEnabled==false) return;

    //! ---------------------------------------------
    //! preliminary: retrieve the simulation manager
    //! ---------------------------------------------
    SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
    SimulationNodeClass *node = sm->myTreeView->currentIndex().data(Qt::UserRole).value<SimulationNodeClass*>();

    //! ---------------
    //! submenu insert
    //! ---------------
    QMenu *menuInsert = contextMenu->addMenu("Insert");
    menuInsert->setIcon(QIcon(":/icons/icon_insert.png"));

    //! ---------------------------
    //! action insert remote point
    //! ---------------------------
    QAction* ActionInsertRemotePoint = menuInsert->addAction("Insert remote point");
    ActionInsertRemotePoint->setIcon(QIcon(":/icons/icon_remote point.png"));
    ActionInsertRemotePoint->setData(46);

    if(node->getType()!=SimulationNodeClass::nodeType_remotePointRoot)
    {
        if(addCommonActions==true)
        {
            contextMenu->addSeparator();
            contextMenuBuilder::buildCommonActions(contextMenu,isEnabled);
        }
    }
}

//! -----------------------------------------
//! function: buildNamedSelectionContextMenu
//! details:
//! -----------------------------------------
void contextMenuBuilder::buildNamedSelectionContextMenu(QMenu *contextMenu, bool addCommonActions, bool isEnabled)
{
    if(isEnabled==false) return;

    //! ---------------------------------------------
    //! preliminary: retrieve the simulation manager
    //! ---------------------------------------------
    SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
    SimulationNodeClass *node = sm->myTreeView->currentIndex().data(Qt::UserRole).value<SimulationNodeClass*>();
    SimulationNodeClass::nodeType aType = node->getType();

    //! ---------------
    //! submenu insert
    //! ---------------
    QMenu *menuInsert = contextMenu->addMenu("Insert");
    menuInsert->setIcon(QIcon(":/icons/icon_insert.png"));

    //! ------------------------------
    //! action insert named selection
    //! ------------------------------
    QAction* ActionInsertNamedSelection = menuInsert->addAction("Insert named selection");
    ActionInsertNamedSelection->setIcon(QIcon(":/icons/icon_named selection geometry.png"));
    ActionInsertNamedSelection->setData(70);

    //! add separator
    contextMenu->addSeparator();

    //! -------------------------------------
    //! action insert mesh element selection
    //! -------------------------------------
    QAction* ActionInsertMeshElementSelection = menuInsert->addAction("Insert mesh element selection");
    ActionInsertMeshElementSelection->setIcon(QIcon(":/icons/icon_volume mesh.png"));
    ActionInsertMeshElementSelection->setData(87);

    //! a named selection is selected
    if(aType!=SimulationNodeClass::nodeType_namedSelection)
    {
        //! ----------------------
        //! action merge selected
        //! ----------------------
        QList<QModelIndex> modelIndexList = sm->myTreeView->selectionModel()->selectedIndexes();

        //! check if the selected items are all named selections
        bool allNamedSelectionGeometry = true;
        for(int i=0; i<modelIndexList.size(); i++)
        {
            QModelIndex curIndex = modelIndexList[i];
            SimulationNodeClass *curNode = curIndex.data(Qt::UserRole).value<SimulationNodeClass*>();
            if(curNode->getType()!=SimulationNodeClass::nodeType_namedSelectionGeometry)
            {
                allNamedSelectionGeometry = false;
                break;
            }
        }
        if(allNamedSelectionGeometry && modelIndexList.size()>1)
        {
            QAction* ActionMergeNamedSelection = contextMenu->addAction("Merge named selections");
            ActionMergeNamedSelection->setIcon(QIcon(":/icons/icon_merge.png"));
            ActionMergeNamedSelection->setData(88);

            //! add separator
            contextMenu->addSeparator();
        }

        //! -------------------
        //! add common actions
        //! -------------------
        if(addCommonActions==true)
        {
            contextMenu->addSeparator();
            contextMenuBuilder::buildCommonActions(contextMenu,isEnabled);

            //! -------------------------------------------------------------
            //! delete all children items: activate this menu item only
            //! if at least one child exists; actually at least two children
            //! if the current item has an empty/dummy/hidden item
            //! (typically the "Select from list" item
            //! -------------------------------------------------------------
            QStandardItem *curItem = sm->getModel()->itemFromIndex(sm->myTreeView->currentIndex());
            SimulationNodeClass *curNode = curItem->data(Qt::UserRole).value<SimulationNodeClass*>();
            SimulationNodeClass::nodeType type = curNode->getType();
            int NbHiddenItems = 0;
            if(type==SimulationNodeClass::nodeType_connection ||
                    type==SimulationNodeClass::nodeType_remotePointRoot ||
                    type==SimulationNodeClass::nodeType_namedSelection) NbHiddenItems = 1;
            if(curItem->rowCount()>NbHiddenItems)
            {
                QAction *ActionDeleteAllChildrenItems = contextMenu->addAction("Delete all children items");
                ActionDeleteAllChildrenItems->setIcon(QIcon(":/icons/icon_delete.png"));
                ActionDeleteAllChildrenItems->setData(102);
            }
        }
    }
}

//! -----------------------------------
//! function: buildGeometryContextMenu
//! details:
//! -----------------------------------
void contextMenuBuilder::buildGeometryContextMenu(QMenu *contextMenu, bool addCommonActions, bool isEnabled)
{
    if(isEnabled==false) return;
    SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
    SimulationNodeClass *node = sm->myTreeView->currentIndex().data(Qt::UserRole).value<SimulationNodeClass*>();

    //! -------------------------------------------------
    //! an item under the geometry root has been clicked
    //! -------------------------------------------------
    if(node->getType()!=node->getFamily())
    {
        QMenu *menuInsert = contextMenu->addMenu("Insert");
        menuInsert->setIcon(QIcon(":/icons/icon_insert.png"));

        //! -----------
        //! Point mass
        //! -----------
        QAction *ActionInsertPointMass = menuInsert->addAction("Point mass");
        ActionInsertPointMass->setIcon(QIcon(":/icons/icon_point mass.png"));
        ActionInsertPointMass->setData(78);

        //! add a separator
        contextMenu->addSeparator();

        if(node->getType()==SimulationNodeClass::nodeType_geometryBody)
        {
            //! --------------------------
            //! action insert repair tool
            //! --------------------------
            QAction *ActionInsertRepairTool = menuInsert->addAction("Repair tool");
            ActionInsertRepairTool->setIcon(QIcon(":/icons/icon_repair tool.png"));
            ActionInsertRepairTool->setData(62);

            //! add separator
            contextMenu->addSeparator();

            QMenu *menuTools = contextMenu->addMenu("Export");
            menuTools->setIcon(QIcon(":/icons/icon_exporting tools.png"));

            //! ------------------------
            //! action export step file
            //! ------------------------
            QAction *ActionExportStepFile = menuTools->addAction("Export STEP file");
            ActionExportStepFile->setIcon(QIcon(":/icons/icon_export step.png"));
            ActionExportStepFile->setData(76);

            //! ------------------------
            //! action export BREP file
            //! ------------------------
            QAction *ActionExportBREPFile = menuTools->addAction("Export BREP file");
            ActionExportBREPFile->setIcon(QIcon(":/icons/icon_export BREP.png"));
            ActionExportBREPFile->setData(77);
        }

        if(node->getType()==SimulationNodeClass::nodeType_pointMass) contextMenuBuilder::addActionDelete(contextMenu);

        //! add separator
        contextMenu->addSeparator();

        //! ------------------------------------------------------------
        //! the common actions are added also when something is running
        //! ------------------------------------------------------------
        if(addCommonActions)
        {
            if(node->getType()!=SimulationNodeClass::nodeType_pointMass)
            {
                bool isVisible = node->getPropertyValue<bool>("Visible");
                if(isVisible)
                {
                    //! -----
                    //! hide
                    //! -----
                    QAction *ActionHide = contextMenu->addAction("Hide");
                    ActionHide->setIcon(QIcon(":/icons/icon_lamp OFF.png"));
                    ActionHide->setData(63);

                    QAction *ActionHideAllOtherBodies = contextMenu->addAction("Hide all other bodies");
                    ActionHideAllOtherBodies->setIcon(QIcon(":/icons/icon_lamp OFF.png"));
                    ActionHideAllOtherBodies->setData(64);

                    //! ---------------------------------------------------------------------
                    //! a visible body has been clicked: if some body is not visible add the
                    //! action "Show all bodies
                    //! ---------------------------------------------------------------------
                    //! check if some body is hidden
                    QExtendedStandardItem *itemGeometryRoot = static_cast<QExtendedStandardItem*>(static_cast<QStandardItemModel*>(sm->myTreeView->model())
                                                                                                  ->itemFromIndex(sm->myTreeView->currentIndex().parent()));
                    int NbPointMasses = 0;
                    for(int i=0; i<itemGeometryRoot->rowCount(); i++)
                    {
                        QStandardItem *item = itemGeometryRoot->child(i,0);
                        SimulationNodeClass *aNode = item->data(Qt::UserRole).value<SimulationNodeClass*>();
                        if(aNode->getType()==SimulationNodeClass::nodeType_pointMass) NbPointMasses++;
                    }
                    int NbBodies = itemGeometryRoot->rowCount()-NbPointMasses;
                    bool bodyHidden = false;
                    for(int i=0; i<NbBodies; i++)
                    {
                        QExtendedStandardItem *theCurItem = static_cast<QExtendedStandardItem*>(itemGeometryRoot->child(i,0));
                        SimulationNodeClass *theCurNode = theCurItem->data(Qt::UserRole).value<SimulationNodeClass*>();
                        bool isVisible = theCurNode->getPropertyValue<bool>("Visible");
                        if(!isVisible)
                        {
                            bodyHidden = true;
                            break;
                        }
                    }
                    if(bodyHidden)
                    {
                        QAction *ActionShowAllBodies = contextMenu->addAction("Show all bodies");
                        ActionShowAllBodies->setIcon(QIcon(":/icons/icon_lamp ON.png"));
                        ActionShowAllBodies->setData(66);
                    }
                }
                else
                {
                    //! -----
                    //! show
                    //! -----
                    QAction *ActionShow = contextMenu->addAction("Show body");
                    ActionShow->setIcon(QIcon(":/icons/icon_lamp ON.png"));
                    ActionShow->setData(65);

                    QAction *ActionShowAllBodies = contextMenu->addAction("Show all bodies");
                    ActionShowAllBodies->setIcon(QIcon(":/icons/icon_lamp ON.png"));
                    ActionShowAllBodies->setData(66);
                }
            }

            //! add separator
            contextMenu->addSeparator();

            //! ----------------------------------------------
            //! check if more than one item has been selected
            //! ----------------------------------------------
            int NbSelectedItems = sm->myTreeView->selectionModel()->selectedIndexes().length();

            if(NbSelectedItems>1)
            {
                //! --------------
                //! form new part
                //! --------------
                QAction *ActionFormNewPart = contextMenu->addAction("Form new part");
                ActionFormNewPart->setIcon(QIcon(":/icons/icon_form new part.png"));
                ActionFormNewPart->setData(83);

                //! add separator
                contextMenu->addSeparator();
            }

            //! -----------------------------
            //!  check if the item is a part
            //! -----------------------------
            if(node->getType()==SimulationNodeClass::nodeType_geometryPart && isEnabled==true)
            {
                //! -------------
                //! explode part
                //! -------------
                QAction *ActionExplodePart = contextMenu->addAction("Explode part");
                ActionExplodePart->setIcon(QIcon(":/icons/icon_explode part.png"));
                ActionExplodePart->setData(82);

                //! add separator
                contextMenu->addSeparator();
            }

            //! ------------
            //! suppression
            //! ------------
            if(isEnabled==true)
            {
                Property::SuppressionStatus ss = node->getPropertyValue<Property::SuppressionStatus>("Suppressed");
                if(ss==Property::SuppressionStatus_Active)
                {
                    contextMenuBuilder::addActionSuppress(contextMenu);
                    contextMenuBuilder::addActionSuppressAllOther(contextMenu);

                    //! check if some body is suppressed
                    QExtendedStandardItem *itemGeometryRoot = static_cast<QExtendedStandardItem*>(static_cast<QStandardItemModel*>(sm->myTreeView->model())
                                                                                                  ->itemFromIndex(sm->myTreeView->currentIndex().parent()));
                    int NbBodies = itemGeometryRoot->rowCount();
                    bool bodySuppressed = false;
                    for(int i=0; i<NbBodies; i++)
                    {
                        QExtendedStandardItem *theCurItem = static_cast<QExtendedStandardItem*>(itemGeometryRoot->child(i,0));
                        SimulationNodeClass *theCurNode = theCurItem->data(Qt::UserRole).value<SimulationNodeClass*>();
                        Property::SuppressionStatus suppressed = theCurNode->getPropertyValue<Property::SuppressionStatus>("Suppressed");
                        if(suppressed == Property::SuppressionStatus_Suppressed)
                        {
                            bodySuppressed = true;
                            break;
                        }
                    }
                    if(bodySuppressed) contextMenuBuilder::addActionUnsuppressAllBodies(contextMenu);
                }
                else
                {
                    contextMenuBuilder::addActionUnsuppress(contextMenu);
                    contextMenuBuilder::addActionUnsuppressAllBodies(contextMenu);
                    contextMenuBuilder::addActionInvertSuppressionSet(contextMenu);
                }
            }

            //! add separator
            contextMenu->addSeparator();
        }
    }

    //! ---------------------------------------------
    //! action rename - also the root can be renamed
    //! ---------------------------------------------
    contextMenuBuilder::addActionRename(contextMenu);
}

//! ------------------------------------
//! function: buildModelRootContextMenu
//! details:
//! ------------------------------------
void contextMenuBuilder::buildModelRootContextMenu(QMenu *contextMenu, bool addCommonActions, bool isEnabled)
{
    Q_UNUSED(addCommonActions)
    if(isEnabled==false) return;

    //! ------------
    //! menu insert
    //! ------------
    QMenu *menuInsert = contextMenu->addMenu("Insert");
    menuInsert->setIcon(QIcon(":/icons/icon_insert.png"));

    //! ---------------------
    //! add thermal analysis
    //! ---------------------
    QAction *ActionAddStructuralAnalysisBranch = menuInsert->addAction("Structural analysis");
    ActionAddStructuralAnalysisBranch->setIcon(QIcon(":/icons/icon_static structural.png"));
    ActionAddStructuralAnalysisBranch->setData(105);

    menuInsert->addSeparator();

    //! ---------------------
    //! add thermal analysis
    //! ---------------------
    QAction *ActionAddThermalAnalysisBranch = menuInsert->addAction("Thermal analysis");
    ActionAddThermalAnalysisBranch->setIcon(QIcon(":/icons/icon_thermal analysis.png"));
    ActionAddThermalAnalysisBranch->setData(104);

    menuInsert->addSeparator();

    //! ----------------------
    //! add combined analysis
    //! ----------------------
    QAction *ActionAddCombinedAnalysisBranch = menuInsert->addAction("Combined analysis");
    ActionAddCombinedAnalysisBranch->setIcon(QIcon(":/icons/icon_combined analysis.png"));
    ActionAddCombinedAnalysisBranch->setData(107);

    menuInsert->addSeparator();

    //! ---------------------------------
    //! add particles in fields analysis
    //! ---------------------------------
    QAction *ActionAddParcitlesInFieldsAnalysisBranch = menuInsert->addAction("Particles in fields");
    ActionAddParcitlesInFieldsAnalysisBranch->setIcon(QIcon(":/icons/icon_a field.png"));
    ActionAddParcitlesInFieldsAnalysisBranch->setData(300);

    menuInsert->addSeparator();

    //! --------------------------------------------------
    //! add remote points root - only one root is allowed
    //! --------------------------------------------------
    SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
    if(sm->getTreeItem(SimulationNodeClass::nodeType_remotePointRoot)==Q_NULLPTR)
    {
        QAction *ActionInsertRemotePointRoot = menuInsert->addAction("Insert remote point");
        ActionInsertRemotePointRoot->setIcon(QIcon(":/icons/icon_remote point.png"));
        ActionInsertRemotePointRoot->setData(45);
        contextMenu->addSeparator();
    }
}

//! ---------------------------------------------
//! function: buildStructuralAnalysisContextMenu
//! details:
//! ---------------------------------------------
void contextMenuBuilder::buildStructuralAnalysisContextMenu(QMenu *contextMenu, bool addCommonActions, bool isEnabled)
{
    if(isEnabled == false) return;

    //! ---------------------------------------
    //! preliminary: retrieve the current node
    //! ---------------------------------------
    SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
    SimulationNodeClass *node = sm->myTreeView->currentIndex().data(Qt::UserRole).value<SimulationNodeClass*>();
    SimulationNodeClass::nodeType nodeType = node->getType();

#ifdef COSTAMP_VERSION
    if(nodeType == SimulationNodeClass::nodeType_timeStepBuilder)
    {
        //! --------------------------
        //! add run time step builder
        //! --------------------------
        QAction *ActionRunTimeStepBuilder = contextMenu->addAction("Run");
        ActionRunTimeStepBuilder->setIcon(QIcon(":/icons/icon_solve.png"));
        ActionRunTimeStepBuilder->setData(1000);
        contextMenu->addAction(ActionRunTimeStepBuilder);

        //! add separator
        contextMenu->addSeparator();
    }
#endif

    //! --------------------------------------------
    //! add the button Run for the Open Foam reader
    //! --------------------------------------------
    if(nodeType==SimulationNodeClass::nodeType_OpenFoamScalarData)
    {
        QAction *actionRun = contextMenu->addAction("Run");
        actionRun->setIcon(QIcon(":/icons/icon_solve.png"));
        actionRun->setData(236);

        //! add a separator
        contextMenu->addSeparator();

        //! --------------------
        //! acition "Duplicate"
        //! --------------------
        contextMenuBuilder::addActionDuplicate(contextMenu);

        //! add a separator
        contextMenu->addSeparator();

        //! ------------------
        //! action "Suppress"
        //! ------------------
        contextMenuBuilder::addActionSuppress(contextMenu);

        //! add a separator
        contextMenu->addSeparator();

        //! ----------------
        //! action "Delete"
        //! ----------------
        contextMenuBuilder::addActionDelete(contextMenu);

        //! add a separator
        contextMenu->addSeparator();

        //! ----------------
        //! action "Rename"
        //! ----------------
        contextMenuBuilder::addActionRename(contextMenu);

        //! add a separator
        contextMenu->addSeparator();

        return;
    }
    //! ----------------
    //! menu for mapper
    //! ----------------
    if(nodeType==SimulationNodeClass::nodeType_mapper)
    {
        //! -----------------------------------------
        //! initial temperature, initial stress, ...
        //! -----------------------------------------
        QMenu *subMenu = contextMenu->addMenu("Select type");
        subMenu->setIcon(QIcon(":/icons/icon_select type.png"));

        //! ------------------------------
        //! body temperature distribution
        //! ------------------------------
        QAction *ActionInsertBodyTemperatureDist = subMenu->addAction("Body temperature distribution");
        ActionInsertBodyTemperatureDist->setIcon(QIcon(":/icons/icon_insert body temperature dist.png"));
        ActionInsertBodyTemperatureDist->setData(53);

        //! -----------------------------------------------------
        //! insert initial stress - to be implemented. To do ...
        //! -----------------------------------------------------
        QAction *ActionInsertInitialStress = subMenu->addAction("Initial stress");
        ActionInsertInitialStress->setIcon(QIcon(":/icons/icon_insert initial stress.png"));
        ActionInsertInitialStress->setData(80);

        //! -----------------------------------------------------
        //! insert initial strain - to be implemented. To do ...
        //! -----------------------------------------------------
        QAction *ActionInsertInitialStrain = subMenu->addAction("Initial strain");
        ActionInsertInitialStrain->setIcon(QIcon(":/icons/icon_insert initial strain.png"));
        ActionInsertInitialStrain->setData(81);

        //! ------------------------------------------------------------
        //! insert here the tools for various types of data translation
        //! ------------------------------------------------------------
        QMenu *subMenu1 = contextMenu->addMenu("Tools");
        subMenu1->setIcon(QIcon(":/icons/icon_tools.png"));

        //! -----------------
        //! open foam reader
        //! -----------------
        QAction *ActionInsertOpenFoamScalarFieldTranslator = subMenu1->addAction("OpenFoam scalar field translator");
        ActionInsertOpenFoamScalarFieldTranslator->setIcon(QIcon(":/icons/icon_openfoam.png"));
        ActionInsertOpenFoamScalarFieldTranslator->setData(100);

        //! add separator
        contextMenu->addSeparator();

        //! --------------------
        //! acition "Duplicate"
        //! --------------------
        contextMenuBuilder::addActionDuplicate(contextMenu);

        //! add a separator
        contextMenu->addSeparator();

        //! ------------------
        //! action "Suppress"
        //! ------------------
        contextMenuBuilder::addActionSuppress(contextMenu);

        //! add a separator
        contextMenu->addSeparator();

        //! ----------------
        //! action "Delete"
        //! ----------------
        contextMenuBuilder::addActionDelete(contextMenu);

        //! add a separator
        contextMenu->addSeparator();

        //! ----------------
        //! action "Rename"
        //! ----------------
        contextMenuBuilder::addActionRename(contextMenu);

        //! add a separator
        contextMenu->addSeparator();

        return;
    }

    //! ----------------------------------------------------------------------------
    //! don't insert boundary conditions if an imported body scalar, a post object,
    //! ( = result of a time sequence interpolation) has been selected
    //! ----------------------------------------------------------------------------
    if(nodeType==SimulationNodeClass::nodeType_mapper ||
            nodeType==SimulationNodeClass::nodeType_importedBodyScalar ||
            nodeType==SimulationNodeClass::nodeType_postObject) return;

    //! -----------------------------------------------
    //! sub menu insert structural boundary conditions
    //! -----------------------------------------------
    QMenu *menuInsertStructural = contextMenu->addMenu("Structural");
    menuInsertStructural->setIcon(QIcon(":/icons/icon_insert.png"));

    //! -----------------------------------------
    //! insert standard earth gravity. To do ...
    //! -----------------------------------------

    //! --------------------
    //! insert acceleration
    //! --------------------
    QAction *ActionInsertAcceleration = menuInsertStructural->addAction("Acceleration");
    ActionInsertAcceleration->setIcon(QIcon(":/icons/icon_acceleration.png"));
    ActionInsertAcceleration->setData(40);

    //! ---------------------------
    //! insert rotational velocity
    //! ---------------------------
    QAction *ActionRotationalVelocity = menuInsertStructural->addAction("Rotational velocity");
    ActionRotationalVelocity->setIcon(QIcon(":/icons/icon_rotational velocity.png"));
    ActionRotationalVelocity->setData(41);

    //! add a separator
    menuInsertStructural->addSeparator();

    //! ---------------------
    //! insert fixed support
    //! ---------------------
    QAction *ActionInsertFixedSupport = menuInsertStructural->addAction("Fixed support");
    ActionInsertFixedSupport->setIcon(QIcon(":/icons/icon_BC fixed support.png"));
    ActionInsertFixedSupport->setData(31);

    //! ---------------------------
    //! insert cylindrical support
    //! ---------------------------
    QAction *ActionInsertCylindricalSupport = menuInsertStructural->addAction("Cylindrical support");
    ActionInsertCylindricalSupport->setIcon(QIcon(":/icons/icon_BC cylindrical support.png"));
    ActionInsertCylindricalSupport->setData(32);

    //! ----------------------------
    //! insert frictionless support
    //! ----------------------------
    QAction *ActionInsertFrictionlessSupport = menuInsertStructural->addAction("Frictionless support");
    ActionInsertFrictionlessSupport->setIcon(QIcon(":/icons/icon_BC frictionless support.png"));
    ActionInsertFrictionlessSupport->setData(33);

    //! -------------------------
    //! compression only support
    //! -------------------------
    QAction *ActionInsertCompressionOnlySupport = menuInsertStructural->addAction("Compression only support");
    ActionInsertCompressionOnlySupport->setIcon(QIcon(":/icons/icon_compression only support.png"));
    ActionInsertCompressionOnlySupport->setData(108);

    //! add a separator
    menuInsertStructural->addSeparator();

    //! -------------------------
    //! insert thermal condition
    //! -------------------------
    QAction *ActionInsertThermalCondition = menuInsertStructural->addAction("Thermal condition");
    ActionInsertThermalCondition->setIcon(QIcon(":/icons/icon_thermal analysis.png"));
    ActionInsertThermalCondition->setData(38);

    //! add a separator
    menuInsertStructural->addSeparator();

    //! --------------------
    //! insert displacement
    //! --------------------
    QAction *ActionInsertDisplacement = menuInsertStructural->addAction("Displacement");
    ActionInsertDisplacement->setIcon(QIcon(":/icons/icon_BC displacement.png"));
    ActionInsertDisplacement->setData(37);

    //! -------------
    //! insert force
    //! -------------
    QAction *ActionInsertForce = menuInsertStructural->addAction("Force");
    ActionInsertForce->setIcon(QIcon(":/icons/icon_BC force.png"));
    ActionInsertForce->setData(34);

    //! --------------
    //! insert moment
    //! --------------
    QAction *ActionInsertMoment = menuInsertStructural->addAction("Moment");
    ActionInsertMoment->setIcon(QIcon(":/icons/icon_BC moment.png"));
    ActionInsertMoment->setData(35);

    //! ----------------
    //! insert pressure
    //! ----------------
    QAction *ActionInsertPressure = menuInsertStructural->addAction("Pressure");
    ActionInsertPressure->setIcon(QIcon(":/icons/icon_BC pressure.png"));
    ActionInsertPressure->setData(36);

    //! add a separator
    menuInsertStructural->addSeparator();

    //! --------------------
    //! insert remote force
    //! --------------------
    QAction *ActionInsertRemoteForce = menuInsertStructural->addAction("Remote force");
    ActionInsertRemoteForce->setIcon(QIcon(":/icons/icon_remote force.png"));
    ActionInsertRemoteForce->setData(29);

    //! ---------------------------
    //! insert remote displacement
    //! ---------------------------
    QAction *ActionInsertRemoteDisplacement = menuInsertStructural->addAction("Remote displacement");
    ActionInsertRemoteDisplacement->setIcon(QIcon(":/icons/icon_remote displacement.png"));
    ActionInsertRemoteDisplacement->setData(28);

    //! -----------------------
    //! insert remote rotation
    //! -----------------------
    QAction *ActionInsertRemoteRotation = menuInsertStructural->addAction("Remote rotation");
    ActionInsertRemoteRotation->setIcon(QIcon(":/icons/icon_rotation.png"));
    ActionInsertRemoteRotation->setData(27);

    //! add a separator
    menuInsertStructural->addSeparator();

    //! -----------------------
    //! insert bolt pretension
    //! -----------------------
    QAction *ActionInsertBoltPretension = menuInsertStructural->addAction("Bolt pretension");
    ActionInsertBoltPretension->setIcon(QIcon(":/icons/icon_bolt.png"));
    ActionInsertBoltPretension->setData(73);

    //! add a separator
    menuInsertStructural->addSeparator();

    //! ----------------------------
    //! insert imported temperature
    //! ----------------------------
    QAction *ActionInsertImportedBodyTemperature = menuInsertStructural->addAction("Imported body temperature");
    ActionInsertImportedBodyTemperature->setIcon(QIcon(":/icons/icon_insert body temperature dist.png"));
    ActionInsertImportedBodyTemperature->setData(97);

    //! add a separator
    menuInsertStructural->addSeparator();

    //! --------------
    //! external data
    //! --------------
    QAction *ActionInsertExternalData = menuInsertStructural->addAction("External data");
    ActionInsertExternalData->setIcon(QIcon(":/icons/icon_mapping.png"));
    ActionInsertExternalData->setData(71);

    //! add a separator
    menuInsertStructural->addSeparator();

    //! -------------
    //! model change
    //! -------------
    QAction *ActionInsertModelChange = menuInsertStructural->addAction("Model change");
    ActionInsertModelChange->setIcon(QIcon(":/icons/icon_balordo.png"));
    ActionInsertModelChange->setData(239);

    //! add a separator
    contextMenu->addSeparator();

    if(nodeType==SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteForce ||
            nodeType==SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteDisplacement ||
            nodeType==SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteRotation ||
            nodeType==SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Moment)
    {
        if(node->getPropertyValue<Property::ScopingMethod>("Scoping method")!=Property::ScopingMethod_RemotePoint)
        {
            //! ------------------------------------------------
            //! insert promote to "Remote point" (only if
            //! the current scoping method is defined through
            //! a geometry or a named selection, i.e. there is
            //! something to "promote")
            //! ------------------------------------------------
            if(!node->getPropertyValue<QVector<GeometryTag>>("Tags").isEmpty())
            {
                QAction *ActionPromoteToRemotePoint = contextMenu->addAction("Promote to remote point");
                ActionPromoteToRemotePoint->setIcon(QIcon(":/icons/icon_pin.png"));
                ActionPromoteToRemotePoint->setData(26);
            }
        }

        if(node->getPropertyValue<Property::ScopingMethod>("Scoping method")==Property::ScopingMethod_NamedSelection)
        {
            //! -----------------------------------------------
            //! insert promote to "Promote to named selection"
            //! (only if the current scoping method is through
            //! a geometry selection, i.e. there is something
            //! to "promote")
            //! -----------------------------------------------
            if(!node->getPropertyValue<QVector<GeometryTag>>("Tags").isEmpty())
            {
                QAction *ActionPromoteToNamedSelection = contextMenu->addAction("Promote to named selection");
                ActionPromoteToNamedSelection->setIcon(QIcon(":/icons/icon_pin.png"));
                ActionPromoteToNamedSelection->setData(27);
            }
        }
    }

    if(nodeType==SimulationNodeClass::nodeType_structuralAnalysisBoltPretension)
    {
        //! add a separator
        contextMenu->addSeparator();

        //! ---------------------------
        //! action display sliced mesh
        //! ---------------------------
        QAction *ActionDisplaySlicedMesh = contextMenu->addAction("Show pre-tensioning elements");
        ActionDisplaySlicedMesh->setIcon(QIcon(":/icons/icon_bolt.png"));
        ActionDisplaySlicedMesh->setData(85);

        //! add a separator
        contextMenu->addSeparator();

        //! ----------------------------------
        //! replicate bolt on twin geometries
        //! ----------------------------------
        QAction *ActionReplicateOnTwinGeometries = contextMenu->addAction("Replicate on twin geometries");
        ActionReplicateOnTwinGeometries->setIcon(QIcon(":/icons/icon_replicate on twins.png"));
        ActionReplicateOnTwinGeometries->setData(86);

        //! add a separator
        contextMenu->addSeparator();
    }

    //! -----------------------
    //! add the common actions
    //! -----------------------
    if(addCommonActions==true && isEnabled==true)
    {
        if(nodeType!=SimulationNodeClass::nodeType_thermalAnalysis &&
                nodeType!=SimulationNodeClass::nodeType_thermalAnalysisSettings &&
                nodeType!=SimulationNodeClass::nodeType_structuralAnalysis &&
                nodeType!=SimulationNodeClass::nodeType_structuralAnalysisSettings)
        {
            contextMenu->addSeparator();
            contextMenuBuilder::buildCommonActions(contextMenu,isEnabled);
        }

        //! add separator
        contextMenu->addSeparator();

        contextMenuBuilder::addActionCreateNamedSelection(contextMenu);
    }

    //! ---------------------------------------
    //! add the action "Solve" if the selected
    //! item is an Analysis root
    //! ---------------------------------------
    if(node->isAnalysisRoot() && addCommonActions == true)
    {
        contextMenuBuilder::addActionRename(contextMenu);
        contextMenu->addSeparator();
        contextMenuBuilder::addActionDelete(contextMenu);
    }
}

//! ---------------------------------
//! function: buildImportContextMenu
//! details:
//! ---------------------------------
void contextMenuBuilder::buildImportContextMenu(QMenu *contextMenu, bool addCommonActions, bool isEnabled)
{
    Q_UNUSED(addCommonActions)
    if(isEnabled==false) return;

    SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
    SimulationNodeClass *node = sm->myTreeView->currentIndex().data(Qt::UserRole).value<SimulationNodeClass*>();

    QString sourceFileLocation = node->getPropertyValue<QString>("Source file path");
    if(!sourceFileLocation.isEmpty())
    {
        //! ----------------
        //! Action Generate
        //! ----------------
        QAction *actionGenerate = contextMenu->addAction("Generate");
        actionGenerate->setIcon(QIcon(":/icons/icon_solve.png"));
        actionGenerate->setData(42);
    }
}

//! ------------------------------------------
//! function: buildThermalAnalysisContextMenu
//! details:
//! ------------------------------------------
void contextMenuBuilder::buildThermalAnalysisContextMenu(QMenu *contextMenu, bool addCommonActions, bool isEnabled)
{
    if(isEnabled == false) return;

    //! ------------
    //! preliminary
    //! ------------
    SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
    SimulationNodeClass *node = sm->myTreeView->currentIndex().data(Qt::UserRole).value<SimulationNodeClass*>();
    SimulationNodeClass::nodeType nodeType = node->getType();

    //! --------------------------------------------
    //! sub menu insert thermal boundary conditions
    //! --------------------------------------------
    QMenu *menuInsertThermal = contextMenu->addMenu("Thermal");
    menuInsertThermal->setIcon(QIcon(":/icons/icon_insert.png"));

    contextMenu->addMenu(menuInsertThermal);

    //! -------------------
    //! insert temperature
    //! -------------------
    QAction *actionInsertTemperature = menuInsertThermal->addAction("Temperature");
    actionInsertTemperature->setIcon(QIcon(":/icons/icon_temperature.png"));
    actionInsertTemperature->setData(90);

    //! add a separator
    menuInsertThermal->addSeparator();

    //! ------------------
    //! insert convection
    //! ------------------
    QAction *actionInsertConvection = menuInsertThermal->addAction("Convection");
    actionInsertConvection->setIcon(QIcon(":/icons/icon_convection.png"));
    actionInsertConvection->setData(91);

    //! add a separator
    menuInsertThermal->addSeparator();

    //! ----------------------
    //! insert adiabatic wall
    //! ----------------------
    QAction *actionInsertAdiabaticWall = menuInsertThermal->addAction("Adiabatic wall");
    actionInsertAdiabaticWall->setIcon(QIcon(":/icons/icon_adiabatic wall.png"));
    actionInsertAdiabaticWall->setData(94);

    //! --------------------
    //! insert thermal flux
    //! --------------------
    QAction *actionInsertThermalFlux = menuInsertThermal->addAction("Thermal flux");
    actionInsertThermalFlux->setIcon(QIcon(":/icons/icon_thermal flux.png"));
    actionInsertThermalFlux->setData(92);

    //! --------------------
    //! insert thermal flow
    //! --------------------
    QAction *actionInsertThermalFlow = menuInsertThermal->addAction("Thermal flow");
    actionInsertThermalFlow->setIcon(QIcon(":/icons/icon_thermal flow.png"));
    actionInsertThermalFlow->setData(93);

    //! add a separator
    menuInsertThermal->addSeparator();

    //! -----------------
    //! insert radiation
    //! -----------------
    QAction *actionInsertRadiation = menuInsertThermal->addAction("Radiation");
    actionInsertRadiation->setIcon(QIcon(":/icons/icon_radiation.png"));
    actionInsertRadiation->setData(95);

    //! add a separator
    menuInsertThermal->addSeparator();

    //! ---------------------
    //! insert thermal power
    //! ---------------------
    QAction *actionInsertThermalPower = menuInsertThermal->addAction("Thermal power");
    actionInsertThermalPower->setIcon(QIcon(":/icons/icon_thermal power.png"));
    actionInsertThermalPower->setData(96);

    //! add a separator
    contextMenu->addSeparator();

    //! -----------------------
    //! add the common actions
    //! -----------------------
    if(addCommonActions==true && isEnabled==true)
    {
        if(nodeType!=SimulationNodeClass::nodeType_thermalAnalysis &&
                nodeType!=SimulationNodeClass::nodeType_thermalAnalysisSettings &&
                nodeType!=SimulationNodeClass::nodeType_structuralAnalysis &&
                nodeType!=SimulationNodeClass::nodeType_structuralAnalysisSettings)
        {
            contextMenu->addSeparator();
            contextMenuBuilder::buildCommonActions(contextMenu,isEnabled);
        }

        //! add separator
        contextMenu->addSeparator();

        contextMenuBuilder::addActionCreateNamedSelection(contextMenu);
    }

    if(node->isAnalysisRoot() && addCommonActions == true)
    {
        contextMenuBuilder::addActionRename(contextMenu);
        contextMenu->addSeparator();
        contextMenuBuilder::addActionDelete(contextMenu);
    }
}

//! --------------------------------------------
//! function: buildParticlesInFieldsContextMenu
//! details:
//! --------------------------------------------
void contextMenuBuilder::buildParticlesInFieldsContextMenu(QMenu *contextMenu, bool addCommonActions, bool isEnabled)
{
    if(isEnabled==false) return;

    //! ----------------
    //! sub menu insert
    //! ----------------
    QMenu *menuInsert = contextMenu->addMenu("Insert");
    menuInsert->setIcon(QIcon(":/icons/icon_insert.png"));

    //! --------------------------
    //! insert electrostatic wall
    //! --------------------------
    QAction *actionInsertElectrostaticWall = menuInsert->addAction("Electrostatic wall");
    actionInsertElectrostaticWall->setIcon(QIcon(":/icons/icon_electrostatic potential.png"));
    actionInsertElectrostaticWall->setData(301);
}

//! -------------------------------------
//! function: buildPostObjectContextMenu
//! details:
//! -------------------------------------
void contextMenuBuilder::buildPostObjectContextMenu(QMenu *contextMenu)
{
    contextMenuBuilder::addActionDuplicate(contextMenu);

    //! add separator
    contextMenu->addSeparator();

    //! ------------------
    //! add action delete
    //! ------------------
    contextMenuBuilder::addActionDelete(contextMenu);

    //! add separator
    contextMenu->addSeparator();

    //! ------------------
    //! add action rename
    //! ------------------
    contextMenuBuilder::addActionRename(contextMenu);

}

//! ---------------------------------------------
//! function: buildStructuralSolutionContextMenu
//! details:
//! ---------------------------------------------
void contextMenuBuilder::buildStructuralSolutionContextMenu(QMenu *contextMenu, bool addCommonActions, bool isEnabled)
{
    //! ---------------------------------------------
    //! preliminary: retrieve the simulation manager
    //! ---------------------------------------------
    SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
    SimulationNodeClass *node = sm->myTreeView->currentIndex().data(Qt::UserRole).value<SimulationNodeClass*>();

    QMenu *menuInsert = contextMenu->addMenu("Insert");
    menuInsert->setIcon(QIcon(":/icons/icon_insert.png"));

    QMenu *subMenuNodalDisplacement = menuInsert->addMenu("Deformation");
    QMenu *subMenuStress = menuInsert->addMenu("Stress");

    //! add a separator
    menuInsert->addSeparator();

    //! ------------------------------------------------------------------------------
    //! total strain - mechanical strain - thermal strain - equivalent plastic strain
    //! ------------------------------------------------------------------------------
    QMenu *subMenuTotalStrain = menuInsert->addMenu("Total strain");
    QMenu *subMenuMechanicalStrain = menuInsert->addMenu("Mechanical strain");
    QMenu *subMenuThermalStrain = menuInsert->addMenu("Thermal strain");
    QMenu *subMenuEquivalentPlasticStrain = menuInsert->addMenu("Equivalent plastic strain");

    //! ----------------
    //! add a separator
    //! ----------------
    menuInsert->addSeparator();

    //! -----------------
    //! sub-menu contact
    //! -----------------
    QMenu *subMenuContact = menuInsert->addMenu("Contact");
    subMenuContact->setIcon(QIcon(":/icons/icon_insert contact.png"));

    //! add a separator
    menuInsert->addSeparator();

    //! ---------------------
    //! sub menu temperature
    //! ---------------------
    QMenu *subMenuTemperature = menuInsert->addMenu("Temperature");

    //! add a separator
    menuInsert->addSeparator();

    //! ----------------------
    //! sub menu nodal forces
    //! ----------------------
    QMenu *subMenuNodalForces =menuInsert->addMenu("Nodal forces");
    subMenuNodalForces->setIcon(QIcon(":/icons/icon_nodal force.png"));

    //! ----------------------------
    //! sub menu nodal displacement
    //! ----------------------------
    QAction *ActionInsertTotalDisplacement = subMenuNodalDisplacement->addAction("Total");
    ActionInsertTotalDisplacement->setIcon(QIcon(":/icons/icon_deformation.png"));
    ActionInsertTotalDisplacement->setData(201);

    QAction *ActionInsertDirectionalDisplacement = subMenuNodalDisplacement->addAction("Directional");
    ActionInsertDirectionalDisplacement->setIcon(QIcon(":/icons/icon_deformation.png"));
    ActionInsertDirectionalDisplacement->setData(206);

    //! ----------------
    //! sub menu stress
    //! ----------------
    QAction *ActionInsertEquivalentStress = subMenuStress->addAction("Stress equivalent");
    ActionInsertEquivalentStress->setIcon(QIcon(":/icons/icon_spring.png"));
    ActionInsertEquivalentStress->setData(202);

    QAction *ActionInsertStressIntensity = subMenuStress->addAction("Stress intensity");
    ActionInsertStressIntensity->setIcon(QIcon(":/icons/icon_spring.png"));
    ActionInsertStressIntensity->setData(210);

    QAction *ActionInsertMaximumPrincipal = subMenuStress->addAction("Maximum principal");
    ActionInsertMaximumPrincipal->setIcon(QIcon(":/icons/icon_spring.png"));
    ActionInsertMaximumPrincipal->setData(207);

    QAction *ActionInsertMiddlePrincipal = subMenuStress->addAction("Middle principal");
    ActionInsertMiddlePrincipal->setIcon(QIcon(":/icons/icon_spring.png"));
    ActionInsertMiddlePrincipal->setData(208);

    QAction *ActionInsertMinimumPrincipal = subMenuStress->addAction("Minimum principal");
    ActionInsertMinimumPrincipal->setIcon(QIcon(":/icons/icon_spring.png"));
    ActionInsertMinimumPrincipal->setData(209);

    {
        //! ----------------------
        //! sub menu total strain
        //! ----------------------
        QAction *ActionInsertEquivalentStrain = subMenuTotalStrain->addAction("Equivalent total strain");
        ActionInsertEquivalentStrain->setIcon(QIcon(":/icons/icon_spring.png"));
        ActionInsertEquivalentStrain->setData(203);

        QAction *ActionInsertStrainIntensity = subMenuTotalStrain->addAction("Total strain intensity");
        ActionInsertStrainIntensity->setIcon(QIcon(":/icons/icon_spring.png"));
        ActionInsertStrainIntensity->setData(214);

        QAction *ActionInsertMaxPrincipalStrain = subMenuTotalStrain->addAction("Total maximum principal");
        ActionInsertMaxPrincipalStrain->setIcon(QIcon(":/icons/icon_spring.png"));
        ActionInsertMaxPrincipalStrain->setData(211);

        QAction *ActionInsertMiddlePrincipalStrain = subMenuTotalStrain->addAction("Total middle principal");
        ActionInsertMiddlePrincipalStrain->setIcon(QIcon(":/icons/icon_spring.png"));
        ActionInsertMiddlePrincipalStrain->setData(212);

        QAction *ActionInsertMinPrincipalStrain = subMenuTotalStrain->addAction("Total minimum principal");
        ActionInsertMinPrincipalStrain->setIcon(QIcon(":/icons/icon_spring.png"));
        ActionInsertMinPrincipalStrain->setData(213);
    }

    {
        //! ---------------------------
        //! sub menu mechanical strain
        //! ---------------------------
        QAction *ActionInsertMechanicalEquivalentStrain = subMenuMechanicalStrain->addAction("Equivalent mechanical strain");
        ActionInsertMechanicalEquivalentStrain->setIcon(QIcon(":/icons/icon_spring.png"));
        ActionInsertMechanicalEquivalentStrain->setData(215);

        QAction *ActionInsertMechanicalStrainIntensity = subMenuMechanicalStrain->addAction("Mechanical strain intensity");
        ActionInsertMechanicalStrainIntensity->setIcon(QIcon(":/icons/icon_spring.png"));
        ActionInsertMechanicalStrainIntensity->setData(216);

        QAction *ActionInsertMechanicalMaxPrincipalStrain = subMenuMechanicalStrain->addAction("Mechanical maximum principal");
        ActionInsertMechanicalMaxPrincipalStrain->setIcon(QIcon(":/icons/icon_spring.png"));
        ActionInsertMechanicalMaxPrincipalStrain->setData(219);

        QAction *ActionInsertMechanicalMiddlePrincipalStrain = subMenuMechanicalStrain->addAction("Mechanical middle principal");
        ActionInsertMechanicalMiddlePrincipalStrain->setIcon(QIcon(":/icons/icon_spring.png"));
        ActionInsertMechanicalMiddlePrincipalStrain->setData(220);

        QAction *ActionInsertMechanicalMinPrincipalStrain = subMenuMechanicalStrain->addAction("Mechanical minimum principal");
        ActionInsertMechanicalMinPrincipalStrain->setIcon(QIcon(":/icons/icon_spring.png"));
        ActionInsertMechanicalMinPrincipalStrain->setData(221);
    }

    {
        //! ------------------------
        //! sub menu thermal strain
        //! ------------------------
        QAction *ActionInsertThermalEquivalentStrain = subMenuThermalStrain->addAction("Equivalent thermal strain");
        ActionInsertThermalEquivalentStrain->setIcon(QIcon(":/icons/icon_spring.png"));
        ActionInsertThermalEquivalentStrain->setData(222);

        QAction *ActionInsertThermalStrainIntensity = subMenuThermalStrain->addAction("Thermal strain intensity");
        ActionInsertThermalStrainIntensity->setIcon(QIcon(":/icons/icon_spring.png"));
        ActionInsertThermalStrainIntensity->setData(223);

        QAction *ActionInsertThermalMaxPrincipalStrain = subMenuThermalStrain->addAction("Thermal maximum principal");
        ActionInsertThermalMaxPrincipalStrain->setIcon(QIcon(":/icons/icon_spring.png"));
        ActionInsertThermalMaxPrincipalStrain->setData(224);

        QAction *ActionInsertThermalMiddlePrincipalStrain = subMenuThermalStrain->addAction("Thermal middle principal");
        ActionInsertThermalMiddlePrincipalStrain->setIcon(QIcon(":/icons/icon_spring.png"));
        ActionInsertThermalMiddlePrincipalStrain->setData(225);

        QAction *ActionInsertThermalMinPrincipalStrain = subMenuThermalStrain->addAction("Thermal minimum principal");
        ActionInsertThermalMinPrincipalStrain->setIcon(QIcon(":/icons/icon_spring.png"));
        ActionInsertThermalMinPrincipalStrain->setData(226);
    }

    //! ----------------------
    //! sub menu nodal forces
    //! ----------------------
    QAction *ActionInsertTotalForce = subMenuNodalForces->addAction("Total nodal force");
    ActionInsertTotalForce->setIcon(QIcon(":/icons/icon_nodal force.png"));
    ActionInsertTotalForce->setData(230);

    QAction *ActionInsertDirectionalForces = subMenuNodalForces->addAction("Directional nodal force");
    ActionInsertDirectionalForces->setIcon(QIcon(":/icons/icon_nodal force.png"));
    ActionInsertDirectionalForces->setData(231);

    //! -----------------
    //! sub menu contact
    //! -----------------
    QAction *ActionInsertContactPressure = subMenuContact->addAction("Contact pressure");
    ActionInsertContactPressure->setIcon(QIcon(":/icons/icon_contact pressure.png"));
    ActionInsertContactPressure->setData(232);

    QAction *ActionInsertContactFrictionalStress = subMenuContact->addAction("Frictional stress");
    ActionInsertContactFrictionalStress->setIcon(QIcon(":/icons/icon_contact friction.png"));
    ActionInsertContactFrictionalStress->setData(233);

    QAction *ActionInsertContactPenetration = subMenuContact->addAction("Penetration");
    ActionInsertContactPenetration->setIcon(QIcon(":/icons/icon_contact penetration.png"));
    ActionInsertContactPenetration->setData(234);

    QAction *ActionInsertContactSliding = subMenuContact->addAction("Sliding");
    ActionInsertContactSliding->setIcon(QIcon(":/icons/icon_contact sliding.png"));
    ActionInsertContactSliding->setData(235);

    //! ---------------------
    //! sub menu temperature
    //! ---------------------
    QAction *ActionInsertTemperature = subMenuTemperature->addAction("Temperature");
    ActionInsertTemperature->setIcon(QIcon(":/icons/icon_insert body temperature dist.png"));
    ActionInsertTemperature->setData(228);

    //! ------------------------------------
    //! sub menu equivalent plastic strain
    //! ------------------------------------
    QAction *ActionInsertEPS = subMenuEquivalentPlasticStrain->addAction("Equivalent plastic strain");
    ActionInsertEPS->setIcon(QIcon(":/icons/icon_deformation.png"));            //! change item icon
    ActionInsertEPS->setData(229);

    //! add separator
    contextMenu->addSeparator();

    //! -------------
    //! fatigue tool
    //! -------------
    QAction *ActionInsertFatigueTool = menuInsert->addAction("Fatigue tool");
    ActionInsertFatigueTool->setIcon(QIcon(":/icons/icon_crack.png"));
    ActionInsertFatigueTool->setData(227);

    //! add separator
    contextMenu->addSeparator();

    if(node->getFamily() == SimulationNodeClass::nodeType_StructuralAnalysisSolution
            && node->getType() != SimulationNodeClass::nodeType_StructuralAnalysisSolutionInformation
            && node->getType() != SimulationNodeClass::nodeType_StructuralAnalysisSolution)
    {

        QAction *ActionEvaluateResult = contextMenu->addAction("Evaluate result");
        ActionEvaluateResult->setIcon(QIcon(":/icons/icon_solve.png"));
        ActionEvaluateResult->setData(204);
    }

    //! --------------------------------------------------
    //! generated on "Solution" or "Solution information"
    //! --------------------------------------------------
    if(node->getType() !=SimulationNodeClass::nodeType_combinedAnalysisSolution && node->getType()!=SimulationNodeClass::nodeType_combinedAnalysisSolutionInformation)
    {
        if(node->getType() == SimulationNodeClass::nodeType_StructuralAnalysisSolution ||
                node->getType() == SimulationNodeClass::nodeType_StructuralAnalysisSolutionInformation)
        {
            QAction *ActionEvaluateResult = contextMenu->addAction("Evaluate all results");
            ActionEvaluateResult->setIcon(QIcon(":/icons/icon_solve.png"));
            ActionEvaluateResult->setData(204);
        }

        //! add separator
        contextMenu->addSeparator();

        //! ---------------------
        //! Clear generated data
        //! ---------------------
        QAction *ActionClearGeneratedData = contextMenu->addAction("Clear generated data");
        ActionClearGeneratedData->setIcon(QIcon(":/icons/icon_clear data.png"));
        ActionClearGeneratedData->setData(205);
    }

    //! add separator
    contextMenu->addSeparator();

    if(addCommonActions)
    {
        //! ----------------------
        //! insert common actions
        //! ----------------------
        if(node->getType()!=SimulationNodeClass::nodeType_StructuralAnalysisSolutionInformation &&
                node->getType()!=SimulationNodeClass::nodeType_StructuralAnalysisSolution)
        {
            contextMenu->addSeparator();
            contextMenuBuilder::buildCommonActions(contextMenu,isEnabled);
        }
    }

    //! add separator
    contextMenu->addSeparator();

    QAction *ActionExport = contextMenu->addAction("Export");
    ActionExport->setIcon(QIcon(":/icons/icon_exporting tools.png"));
    ActionExport->setData(109);

    //! add separator
    contextMenu->addSeparator();

    contextMenuBuilder::addActionCreateNamedSelection(contextMenu);
}

//! -----------------------------------------
//! function: buildThermalResultsContextMenu
//! details:
//! -----------------------------------------
void contextMenuBuilder::buildThermalResultsContextMenu(QMenu *contextMenu, bool addCommonActions, bool isEnabled)
{
    if(isEnabled == false) return;

    //! ---------------------------------------------
    //! preliminary: retrieve the simulation manager
    //! ---------------------------------------------
    SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
    SimulationNodeClass *node = sm->myTreeView->currentIndex().data(Qt::UserRole).value<SimulationNodeClass*>();

    QMenu *menuInsert = contextMenu->addMenu("Insert");
    menuInsert->setIcon(QIcon(":/icons/icon_insert.png"));

    //! --------------------------
    //! action insert temperature
    //! --------------------------
    QAction *ActionInsertTemperature = menuInsert->addAction("Temperature");
    ActionInsertTemperature->setIcon(QIcon(":/icons/icon_temperature.png"));
    ActionInsertTemperature->setData(240);

    //! -------------------
    //! action insert flux
    //! -------------------
    QAction *ActionInsertFlux = menuInsert->addAction("Flux");
    ActionInsertFlux->setIcon(QIcon(":/icons/icon_thermal flux.png"));
    ActionInsertFlux->setData(241);

    //! add separator
    contextMenu->addSeparator();

    if(node->isAnalysisResult())
    {
        QAction *ActionEvaluateResult = contextMenu->addAction("Evaluate result");
        ActionEvaluateResult->setIcon(QIcon(":/icons/icon_solve.png"));
        ActionEvaluateResult->setData(204);
    }

    //! --------------------------------------------------
    //! generated on "Solution" or "Solution information"
    //! --------------------------------------------------
    if(node->getType() !=SimulationNodeClass::nodeType_combinedAnalysisSolution && node->getType()!=SimulationNodeClass::nodeType_combinedAnalysisSolutionInformation)
    {
        if(node->getType() == SimulationNodeClass::nodeType_thermalAnalysisSolution ||
                node->getType() == SimulationNodeClass::nodeType_thermalAnalysisSolutionInformation)
        {
            QAction *ActionEvaluateResult = contextMenu->addAction("Evaluate all results");
            ActionEvaluateResult->setIcon(QIcon(":/icons/icon_solve.png"));
            ActionEvaluateResult->setData(204);
        }

        //! add separator
        contextMenu->addSeparator();

        //! ---------------------
        //! Clear generated data
        //! ---------------------
        QAction *ActionClearGeneratedData = contextMenu->addAction("Clear generated data");
        ActionClearGeneratedData->setIcon(QIcon(":/icons/icon_clear data.png"));
        ActionClearGeneratedData->setData(205);
    }

    //! add separator
    contextMenu->addSeparator();

    if(addCommonActions)
    {
        //! ----------------------
        //! insert common actions
        //! ----------------------
        if(node->getType()!=SimulationNodeClass::nodeType_thermalAnalysisSolutionInformation &&
                node->getType()!=SimulationNodeClass::nodeType_thermalAnalysisSolution)
        {
            contextMenu->addSeparator();
            contextMenuBuilder::buildCommonActions(contextMenu,isEnabled);
        }
    }

    //! add separator
    contextMenu->addSeparator();

    QAction *ActionExport = contextMenu->addAction("Export");
    ActionExport->setIcon(QIcon(":/icons/icon_exporting tools.png"));
    ActionExport->setData(109);

    //! add separator
    contextMenuBuilder::addActionCreateNamedSelection(contextMenu);

    /*
    if(node->getType()!=SimulationNodeClass::nodeType_thermalAnalysisSolutionInformation)
    {
        //! add separator
        contextMenu->addSeparator();

        contextMenuBuilder::addActionDuplicate(contextMenu);

        //! add separator
        contextMenu->addSeparator();

        contextMenuBuilder::addActionDelete(contextMenu);

        //! add separator
        contextMenu->addSeparator();

        contextMenuBuilder::addActionRename(contextMenu);
        contextMenuBuilder::addActionRenameBasedOnDefinition(contextMenu);

        //! add separator
        contextMenu->addSeparator();

        contextMenuBuilder::addActionCreateNamedSelection(contextMenu);
    }
    */
}

//! -------------------------------------------
//! function: buildCombinedAnalysisContextMenu
//! details:
//! -------------------------------------------
void contextMenuBuilder::buildCombinedAnalysisResultsContextMenu(QMenu* contextMenu, bool addCommonActions, bool isEnabled)
{
    if(isEnabled==false) return;

    //! ---------------------------------------------
    //! preliminary: retrieve the simulation manager
    //! ---------------------------------------------
    SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
    SimulationNodeClass *node = sm->myTreeView->currentIndex().data(Qt::UserRole).value<SimulationNodeClass*>();

    //! -------------------
    //! structural results
    //! -------------------
    QMenu *structuralSolutionMenu = contextMenu->addMenu("Structural");
    structuralSolutionMenu->setIcon(QIcon(":/icons/icon_static structural.png"));
    contextMenuBuilder::buildStructuralSolutionContextMenu(structuralSolutionMenu,false,isEnabled);

    contextMenu->addSeparator();

    //! ----------------
    //! thermal results
    //! ----------------
    QMenu *thermalSolutionMenu = contextMenu->addMenu("Thermal");
    thermalSolutionMenu->setIcon(QIcon(":/icons/icon_thermal analysis.png"));
    contextMenuBuilder::buildThermalResultsContextMenu(thermalSolutionMenu,false,isEnabled);

    if(addCommonActions)
    {
        //! ----------------------
        //! insert common actions
        //! ----------------------
        if(node->getType()!=SimulationNodeClass::nodeType_combinedAnalysisSolutionInformation &&
                node->getType()!=SimulationNodeClass::nodeType_combinedAnalysisSolution)
        {
            contextMenu->addSeparator();
            contextMenuBuilder::buildCommonActions(contextMenu,isEnabled);
        }
    }

    //! add separator
    contextMenu->addSeparator();

    if(node->getType() == SimulationNodeClass::nodeType_combinedAnalysisSolution ||
      node->getType() == SimulationNodeClass::nodeType_combinedAnalysisSolutionInformation)
    {
        QAction *ActionEvaluateResult = contextMenu->addAction("Evaluate all results");
        ActionEvaluateResult->setIcon(QIcon(":/icons/icon_solve.png"));
        ActionEvaluateResult->setData(204);

        //! add separator
        contextMenu->addSeparator();

        //! ---------------------
        //! Clear generated data
        //! ---------------------
        QAction *ActionClearGeneratedData = contextMenu->addAction("Clear generated data");
        ActionClearGeneratedData->setIcon(QIcon(":/icons/icon_clear data.png"));
        ActionClearGeneratedData->setData(205);
    }

    //! add separator
    contextMenu->addSeparator();

    contextMenuBuilder::addActionCreateNamedSelection(contextMenu);
}

//! ------------------------------------------------
//! function: -
//! details:  functions for adding "common" actions
//! ------------------------------------------------
void contextMenuBuilder::addActionDelete(QMenu *contextMenu)
{
    QAction *ActionDelete = contextMenu->addAction("Delete");
    ActionDelete->setIcon(QIcon(":/icons/icon_delete.png"));
    ActionDelete->setData(99);
}
void contextMenuBuilder::addActionDuplicate(QMenu *contextMenu)
{
    QAction *ActionDuplicate = contextMenu->addAction("Duplicate");
    ActionDuplicate->setIcon(QIcon(":/icons/icon_duplicate.png"));
    ActionDuplicate->setData(98);
}
void contextMenuBuilder::addActionSuppress(QMenu *contextMenu)
{
    QAction *ActionSuppress = contextMenu->addAction("Suppress");
    ActionSuppress->setIcon(QIcon(":/icons/icon_suppress.png"));
    ActionSuppress->setData(57);
}
void contextMenuBuilder::addActionSuppressAllOther(QMenu *contextMenu)
{
    QAction *ActionSuppressAllOther = contextMenu->addAction("Suppress all other bodies");
    ActionSuppressAllOther->setIcon(QIcon(":/icons/icon_suppress.png"));
    ActionSuppressAllOther->setData(60);
}
void contextMenuBuilder::addActionUnsuppress(QMenu *contextMenu)
{
    QAction *ActionUnsuppress = contextMenu->addAction("Unsuppress");
    ActionUnsuppress->setIcon(QIcon(":/icons/icon_unsuppress.png"));
    ActionUnsuppress->setData(58);
}
void contextMenuBuilder::addActionInvertSuppressionSet(QMenu *contextMenu)
{
    QAction *ActionInvertSuppressionSet = contextMenu->addAction("Invert suppression set");
    ActionInvertSuppressionSet->setIcon(QIcon(":/icons/icon_invert suppression.png"));
    ActionInvertSuppressionSet->setData(59);
}
void contextMenuBuilder::addActionRename(QMenu *contextMenu)
{
    QAction *ActionRename = contextMenu->addAction("Rename");
    ActionRename->setIcon(QIcon(":/icons/icon_rename.png"));
    ActionRename->setData(101);
}
void contextMenuBuilder::addActionRenameBasedOnDefinition(QMenu *contextMenu)
{
    QAction *ActionRename = contextMenu->addAction("Rename based on definition");
    ActionRename->setIcon(QIcon(":/icons/icon_rename.png"));
    ActionRename->setData(103);
}
void contextMenuBuilder::addActionUnsuppressAllBodies(QMenu *contextMenu)
{
    QAction *ActionUnsuppressAllBodies = contextMenu->addAction("Unsuppress all bodies");
    ActionUnsuppressAllBodies->setIcon(QIcon(":/icons/icon_invert suppression.png"));
    ActionUnsuppressAllBodies->setData(61);
}
void contextMenuBuilder::addActionDeleteAllChildrenItems(QMenu *contextMenu)
{
    //! ----------------------------------------------------------
    //! delete all children items
    //! activate this menu item only if at least one child exists
    //! ----------------------------------------------------------
    SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
    if(sm->getModel()->itemFromIndex(sm->myTreeView->currentIndex())->rowCount()>0)
    {
        //! add separator
        contextMenu->addSeparator();

        QAction *ActionDeleteAllChildrenItems = contextMenu->addAction("Delete all children items");
        ActionDeleteAllChildrenItems->setIcon(QIcon(":/icons/icon_delete.png"));
        ActionDeleteAllChildrenItems->setData(102);
    }
}
void contextMenuBuilder::addActionCreateNamedSelection(QMenu *contextMenu)
{
    SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
    const QList<QModelIndex> &selectedIndexes = sm->myTreeView->selectionModel()->selectedIndexes();
    QVector<GeometryTag> geometryTags;
    for(int n=0; n<selectedIndexes.length(); n++)
    {
        SimulationNodeClass *node = selectedIndexes[n].data(Qt::UserRole).value<SimulationNodeClass*>();
        if(node->isSimulationSetUpNode()==false && node->isAnalysisResult()==false) return;
        if(node->getPropertyItem("Tags")!=Q_NULLPTR)
        {
            QVector<GeometryTag> curTags = node->getPropertyValue<QVector<GeometryTag>>("Tags");
            geometryTags.append(curTags);
        }
    }
    if(geometryTags.length()==0) return;
    QAction *ActionCreateNamedSelection = contextMenu->addAction("Create named selection");
    ActionCreateNamedSelection->setIcon(QIcon(":/icons/icon_named selection geometry.png"));
    ActionCreateNamedSelection->setData(106);
}
//! -----------------------------------------------
//! remote displacement                         28
//! remote force                                29
//! remote rotation                             27
//! fixed support                               31
//! cylindrical support                         32
//! compression only support                   108
//! frictionless support                        33
//! force                                       34
//! moment                                      35
//! pressure                                    36
//! displacement                                37
//! thermal condition                           38
//! standard earth gravity                      39
//! acceleration                                40
//! rotational velocity                         41
//! imported body temperature                   97
//! update geometry source file                 42
//! hide                                        63
//! hide all other bodies                       64
//! show body                                   65
//! show all bodies                             66
//! suppress                                    57
//! unsuppress                                  58
//! invert suppression set                      59
//! suppress all other bodies                   60
//! unsuppress all bodies                       61
//! insert repair tool                          62
//! update geometry from source                 73
//! form new part                               83
//! delete body                                 84
//!
//! insert remote point root                    45
//! insert remote point                         46
//! insert coordinate system                    47
//!
//! flip contact                                48
//! insert contact                              49
//! merge selected contacts                     50
//! insert contact group                        44
//! generate automatic connections              43
//!
//! insert mesh method                          74
//! insert body mesh control                    51
//! insert edge sizing                          68
//! insert vertex sizing
//! generate volume mesh                        54
//! generate surface mesh                       55
//! clear all meshes                            56
//! preview inflation                           67
//! insert type of mesh                         79
//! inset mesh metric                           89
//! insert vertex sizing                        69
//! insert named selection                      70
//! merge named selections                      88
//! insert mesh element selection               87
//! insert face sizing                          72
//! inset prismatic layer                       75
//!
//! insert imported body temperature dist       53
//! insert initial stress                       80
//! insert initial strain                       81
//! insert external data                        71
//! insert bolt pretension                      73
//! insert point mass                           78
//! insert openfoam scalar field translator    100
//! rename                                     101
//! rename based on definition                 103
//! duplicate                                   98
//! delete                                      99
//! delete all children items                  102
//! create named selection                     106
//!
//! add structural analysis                    105
//! add thermal analysis                       104
//! add combined analysis                      107
//!
//! update post object                         200
//! insert total displacement                  201
//! insert directional displacement            206
//! insert stress equivalent                   202
//! insert maximum principal                   207
//! insert middle principal                    208
//! insert minimum principal                   209
//!
//! insert equivalent strain                   203
//! insert strain intensity                    214
//! insert max principal strain                211
//! insert middle principal strain             212
//! insert minimum principal strain            213
//!
//! insert equivalent mechanical strain        215
//! insert strain mechanical intensity         216
//! insert max principal mechanical strain     219
//! insert middle principal mechanical strain  220
//! insert minimum principal mechanical strain 221
//!
//! insert equivalent mechanical strain        222
//! insert strain mechanical intensity         223
//! insert max principal mechanical strain     224
//! insert middle principal mechanical strain  225
//! insert minimum principal mechanical strain 226
//! insert fatigue tool                        227
//! insert temperature                         228
//! insert equivalent plastic strain           229
//! insert total nodal force                   230
//! inesrt directional nodal force             231
//! insert contact pressure                    232
//! insert contact frictional stress           233
//! insert contact penetration                 234
//! insert contact sliding                     235
//!
//! insert temperature (thermal solution)      240
//! inset thermal flux (thermal solution)      241
//!
//! evaluate result                            204
//! clear generated data                       205
//!
//! run openfoam reader                        236
//! run mapper                                 237
//! show underlying mesh                       238
//! insert model change                        239
//! run time step builder                     1000
//!
//! export step file                            76
//! export BRep file                            77
//! export (result)                            109
//!
//! display sliced volume mesh                  85
//! replicate on twin geometries                86
//!
//! insert temperature                          90
//! insert convection                           91
//! insert thermal flux                         92
//! insert thermal flow                         93
//! insert thermal power                        96
//! insert adiabatic walls                      94
//! insert radiation                            95
//!
//! insert particles in fields                 300
//! insert electrostatic wall                  301
//! -----------------------------------------------
