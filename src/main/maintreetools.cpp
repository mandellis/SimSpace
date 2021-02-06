//! ----------------
//! custom includes
//! ----------------
#include "maintreetools.h"
#include "simulationnodeclass.h"
#include "src/gui/tabularData/tabulardatacolumns.h"
#include "qextendedstandarditem.h"
#include "ccxsolvermessage.h"
#include "solutioninfo.h"

//! ---
//! Qt
//! ---
#include <QStandardItemModel>

//! ----
//! C++
//! ----
#include <memory>

//! ---------------------------
//! function: getColumnsToRead
//! details:
//! ---------------------------
QList<int> mainTreeTools::getColumnsToRead(QTreeView *tree)
{
    //cout<<"mainTreeTools::getColumnsToRead()->____function called____"<<endl;

    QModelIndex currentIndex=tree->currentIndex();
    QStandardItem *anItem = static_cast<QStandardItemModel*>(tree->model())->itemFromIndex(currentIndex);
    const QList<int> &theColumnsToShow = mainTreeTools::getColumnsToRead(anItem);
    return theColumnsToShow;
}

//! ---------------------------
//! function: getColumnsToRead
//! details:
//! ---------------------------
QList<int> mainTreeTools::getColumnsToRead(QStandardItem *anItem)
{
    //cout<<"mainTreeTools::getColumnsToRead()->____function called____"<<endl;

    QList<int> theColumnsToShow;
    int SC = mainTreeTools::calculateStartColumn(anItem);
    SimulationNodeClass *aNode = anItem->data(Qt::UserRole).value<SimulationNodeClass*>();
    SimulationNodeClass::nodeType theType = aNode->getType();

    //! -----------------
    //! always 3 columns
    //! -----------------
    QList<SimulationNodeClass::nodeType> threeColumns;
    threeColumns<<SimulationNodeClass::nodeType_structuralAnalysisBoltPretension;

    //! -----------------
    //! always 2 columns
    //! -----------------
    QList<SimulationNodeClass::nodeType> twoColumns;
    twoColumns<<SimulationNodeClass::nodeType_thermalAnalysisConvection;

    //! ----------------
    //! always 1 column
    //! ----------------
    QList<SimulationNodeClass::nodeType> oneColumn;
    oneColumn<<SimulationNodeClass::nodeType_modelChange<<
               SimulationNodeClass::nodeType_structuralAnalysisThermalCondition<<
               SimulationNodeClass::nodeType_thermalAnalysisTemperature<<
               SimulationNodeClass::nodeType_thermalAnalysisThermalFlow<<
               SimulationNodeClass::nodeType_thermalAnalysisThermalFlux<<
               SimulationNodeClass::nodeType_thermalAnalysisThermalPower<<
               SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Pressure<<
               SimulationNodeClass::nodeType_modelChange<<
               SimulationNodeClass::nodeType_electrostaticPotential;

    //! -----------------
    //! always 0 columns
    //! -----------------
    QList<SimulationNodeClass::nodeType> zeroColumns;
    zeroColumns<<SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_CylindricalSupport<<
                 SimulationNodeClass::nodeType_structuralAnalysisBoundaryContidion_FixedSupport<<
                 SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_FrictionlessSupport<<
                 SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_CompressionOnlySupport<<
                 SimulationNodeClass::nodeType_thermalAnalysisAdiabaticWall<<
                 SimulationNodeClass::nodeType_particlesInFieldsParticlePack;

    //! ------------
    //! 1/3 columns
    //! ------------
    QList<SimulationNodeClass::nodeType> oneThreeColumns;
    oneThreeColumns<<SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Force<<
                     SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteForce<<
                     SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Moment<<
                     SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Acceleration<<
                     SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RotationalVelocity<<
                     SimulationNodeClass::nodeType_magneticField;

    //! ----------------
    //! 0/1/2/3 columns
    //! ----------------
    QList<SimulationNodeClass::nodeType> zeroOneTwoThreeColumns;
    zeroOneTwoThreeColumns<<SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Displacement<<
                             SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteDisplacement<<
                             SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteRotation;


    if(threeColumns.contains(theType)) theColumnsToShow<<SC<<SC+1<<SC+2;
    if(twoColumns.contains(theType)) theColumnsToShow<<SC<<SC+1;
    if(oneColumn.contains(theType)) theColumnsToShow<<SC;
    if(oneThreeColumns.contains(theType))
    {
        Property::defineBy theDefineBy = aNode->getPropertyValue<Property::defineBy>("Define by");
        if(theDefineBy==Property::defineBy_components) theColumnsToShow<<SC<<SC+1<<SC+2;
        else theColumnsToShow<<SC;
    }
    if(zeroOneTwoThreeColumns.contains(theType))
    {
        Property::defineBy theDefineBy = aNode->getPropertyValue<Property::defineBy>("Define by");
        if(theDefineBy==Property::defineBy_components)
        {
            //! ---------------------------------------------
            //! show only the active displacement components
            //! ---------------------------------------------
            Property::loadDefinition theLoadDefinitionXcomponent = aNode->getPropertyValue<Property::loadDefinition>("X component");
            Property::loadDefinition theLoadDefinitionYcomponent = aNode->getPropertyValue<Property::loadDefinition>("Y component");
            Property::loadDefinition theLoadDefinitionZcomponent = aNode->getPropertyValue<Property::loadDefinition>("Z component");
            bool b0,b1,b2;

            if(theLoadDefinitionXcomponent!=Property::loadDefinition_free) b0=true; else b0=false;
            if(theLoadDefinitionYcomponent!=Property::loadDefinition_free) b1=true; else b1=false;
            if(theLoadDefinitionZcomponent!=Property::loadDefinition_free) b2=true; else b2=false;

            if((b0 == true && b1 == false && b2 == false) ||
                    (b0 == false && b1 == true && b2 == false) ||
                    (b0 == false && b1 == false && b2 == true))
            {
                theColumnsToShow<<SC;
            }
            if((b0 == true && b1 == true && b2 == false) ||
                    (b0 == true && b1 == false && b2 == true) ||
                    (b0 == false && b1 == true && b2 == true))
            {
                theColumnsToShow<<SC<<SC+1;
            }
            if(b0 == true && b1 == true && b2 == true)
            {
                theColumnsToShow<<SC<<SC+1<<SC+2;
            }
            if(b0 == false && b1 == false && b2 == false)
            {
                //! all components are free => no column to read
            }
        }
        else theColumnsToShow<<SC;
    }

    //! ----------------------------
    //! diagnostic - can be removed
    //! ----------------------------
    //cout<<"mainTreeTools::getColumnsToRead()->____columns to read: {";
    //int i; for(i=0;i<theColumnsToShow.length()-1;i++) cout<<theColumnsToShow.at(i)<<",";
    //cout<<theColumnsToShow.at(i)<<"}"<<endl;
    //! ---------------
    //! end diagnostic
    //! ---------------
    return theColumnsToShow;
}

//! -------------------------------
//! function: calculateStartColumn
//! details:
//! -------------------------------
int mainTreeTools::calculateStartColumn(QStandardItem *anItem)
{
    //cout<<"mainTreeTools::calculateStartColumn()->____function called____"<<endl;

    //! -----------------
    //! the start column
    //! -----------------
    int startColumn = 0;

    //! --------------------------------
    //! the row of the item in the tree
    //! --------------------------------
    QModelIndex currentIndex=anItem->index();
    int rrow = currentIndex.row();

    //! ----------------------------------
    //! retrieve the simulation root item
    //! ----------------------------------
    QStandardItem* itemSimulationRoot = anItem->parent();

    //! --------------------------------------------------------------------------------
    //! the for cycle starts from the second item (the previous is "Analysis settings")
    //! --------------------------------------------------------------------------------
    int offset = 0;
    for(int k=1; k<rrow; k++)
    {
        QStandardItem *curItem = itemSimulationRoot->child(k,0);
        SimulationNodeClass *curNode = curItem->data(Qt::UserRole).value<SimulationNodeClass*>();
        SimulationNodeClass::nodeType curType = curNode->getType();
        Property::defineBy theDefineBy;
        int delta = 0;

        switch (curType)
        {
            //! ---------------------
            //! "0" columns in table
            //! ---------------------
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_ImportedTemperatureDistribution:
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_CompressionOnlySupport:
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_CylindricalSupport:
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_FrictionlessSupport:
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryContidion_FixedSupport:
        case SimulationNodeClass::nodeType_thermalAnalysisAdiabaticWall:
        case SimulationNodeClass::nodeType_particlesInFieldsParticlePack:
#ifdef COSTAMP_VERSION
        case SimulationNodeClass::nodeType_timeStepBuilder:
#endif
        case SimulationNodeClass::nodeType_mapper:
            delta = -1;
            break;

            //! --------------------
            //! "1" column in table
            //! --------------------
        case SimulationNodeClass::nodeType_modelChange:
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Pressure:
        case SimulationNodeClass::nodeType_structuralAnalysisThermalCondition:
        case SimulationNodeClass::nodeType_thermalAnalysisRadiation:
        case SimulationNodeClass::nodeType_thermalAnalysisTemperature:
        case SimulationNodeClass::nodeType_thermalAnalysisThermalFlow:
        case SimulationNodeClass::nodeType_thermalAnalysisThermalPower:
        case SimulationNodeClass::nodeType_thermalAnalysisThermalFlux:
        case SimulationNodeClass::nodeType_electrostaticPotential:
            delta = 0;
            break;

            //! ------------
            //! "2" columns
            //! ------------
        case SimulationNodeClass::nodeType_thermalAnalysisConvection:
            delta = 1;
            break;

            //! -------------------
            //! "1" or "3" columns
            //! -------------------
        case SimulationNodeClass::nodeType_magneticField:
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Force:
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Moment:
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteForce:
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Acceleration:
            theDefineBy = curNode->getPropertyValue<Property::defineBy>("Define by");
            if(theDefineBy==Property::defineBy_vector) delta = 0;
            else delta = 2;
            break;

            //! -----------------------------------------------
            //! always "3" columns in table: "bolt pretension"
            //! -----------------------------------------------
        case SimulationNodeClass::nodeType_structuralAnalysisBoltPretension:
            delta = 2;
            break;

            //! ----------------------------------------------------------------
            //! the "Displacement" and "Remote displacement" are special cases,
            //! since their components could have the option "Free"
            //! ----------------------------------------------------------------
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Displacement:
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteDisplacement:
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteRotation:
            theDefineBy = curNode->getPropertyValue<Property::defineBy>("Define by");
            if(theDefineBy==Property::defineBy_vector) delta=0;
            else
            {
                delta = 0;
                int subDelta = 0;
                Property::loadDefinition loadDefinitionXcomponent = curNode->getPropertyValue<Property::loadDefinition>("X component");
                Property::loadDefinition loadDefinitionYcomponent = curNode->getPropertyValue<Property::loadDefinition>("Y component");
                Property::loadDefinition loadDefinitionZcomponent = curNode->getPropertyValue<Property::loadDefinition>("Z component");
                if(loadDefinitionXcomponent!=Property::loadDefinition_free)subDelta++;
                if(loadDefinitionYcomponent!=Property::loadDefinition_free)subDelta++;
                if(loadDefinitionZcomponent!=Property::loadDefinition_free)subDelta++;
                delta = subDelta-1;
            }
            break;
        }
        offset = offset + delta;
    }

    //! --------------------------------------------------
    //! number of columns defining the "Analysis setting"
    //! --------------------------------------------------
    int initNumberOfColumns = NUMBER_OF_COLUMNS_BEFORE_BC_DATA;
    startColumn = (rrow+initNumberOfColumns)+offset;
    //cout<<"mainTreeTools::calculateStartColumn()->____exiting function: "<<startColumn<<"____"<<endl;
    return startColumn;
}

//! -------------------------------
//! function: calculateStartColumn
//! details:
//! -------------------------------
int mainTreeTools::calculateStartColumn(QTreeView *tree)
{
    //! -----------------
    //! the start column
    //! -----------------
    int startColumn = 0;

    QModelIndex currentIndex=tree->currentIndex();
    QStandardItem *currentItem = static_cast<QStandardItemModel*>(tree->model())->itemFromIndex(currentIndex);
    startColumn = mainTreeTools::calculateStartColumn(static_cast<QExtendedStandardItem*>(currentItem));
    return startColumn;
}

//! -------------------------------
//! function: getAllTreeItemOfType
//! details:
//! -------------------------------
QList<QStandardItem*> mainTreeTools::getAllTreeItemOfType(QTreeView *tree, SimulationNodeClass::nodeType theNodeType)
{
    cout<<"mainTreeTools::getAllTreeItemOfType()->____function called____"<<endl;
    QStandardItemModel *theModel = static_cast<QStandardItemModel*>(tree->model());
    if(theModel==Q_NULLPTR)
    {
        cout<<"SimulationManager::getAllTreeItemOfType()->____the tree has not a model____"<<endl;
        return QList<QStandardItem*>();
    }
    QList<QStandardItem*> items, itemsout;
    mainTreeTools::getTreeItemsRecursively(theModel,items);
    for(QList<QStandardItem*>::iterator it = items.begin(); it!=items.end(); ++it)
    {
        QStandardItem* curItem = *it;
        if(curItem->data(Qt::UserRole).value<SimulationNodeClass*>()->getType()==theNodeType) itemsout.append(curItem);
    }
    return itemsout;
}

//! ----------------------------------
//! function: getTreeItemsRecursively
//! details:
//! ----------------------------------
void mainTreeTools::getTreeItemsRecursively(QStandardItemModel* model, QList<QStandardItem *> &items, QModelIndex parent)
{
    for(int r = 0; r<model->rowCount(parent); ++r)
    {
        QModelIndex index = model->index(r, 0, parent);
        QExtendedStandardItem *item = static_cast<QExtendedStandardItem*>(model->itemFromIndex(index));
        items.push_back(item);
        if(model->hasChildren(index))
        {
            getTreeItemsRecursively(model, items, index);
        }
    }
}

//! --------------------------------------------------------------------------------
//! function: getTreeItemRecursively
//! details:  return the first occurrence of item containing a node of type "aType"
//! --------------------------------------------------------------------------------
QStandardItem* mainTreeTools::getFirstTreeItemOfType(SimulationNodeClass::nodeType aType, QStandardItemModel* model)
{
    QStandardItem *curItem = Q_NULLPTR;
    QList<QStandardItem*> items;
    mainTreeTools::getTreeItemsRecursively(model,items);
    for(QList<QStandardItem*>::iterator it = items.begin(); it!= items.end(); it++)
    {
        curItem = *it;
        SimulationNodeClass *curNode = curItem->data(Qt::UserRole).value<SimulationNodeClass*>();
        if(curNode->getType()==aType) break;
    }
    return curItem;
}

//! ---------------------------------
//! function: getTreeItemsFromShapes
//! details:
//! ---------------------------------
bool mainTreeTools::getTreeItemsFromShapes(QTreeView *tree, const std::vector<TopoDS_Shape> &vecShapes, std::vector<QStandardItem*> &vecItems)
{
    QStandardItemModel *model = (QStandardItemModel*)tree->model();
    QStandardItem *item = mainTreeTools::getFirstTreeItemOfType(SimulationNodeClass::nodeType_geometry,model);
    if(item==Q_NULLPTR)
    {
        cout<<"mainTreeTools::getTreeItemsFromShapes()->____item geometry not found____"<<endl;
        return false;
    }
    for(int i=0; i<item->rowCount(); i++)               // scan the children, jumping over point mass items
    {
        QStandardItem *childItem = item->child(i,0);
        SimulationNodeClass *node = childItem->data(Qt::UserRole).value<SimulationNodeClass*>();
        if(node->getType()==SimulationNodeClass::nodeType_pointMass) continue;
        TopoDS_Shape shapeInItem = node->getPropertyValue<TopoDS_Shape>("Shape");
        if(std::find(vecShapes.begin(), vecShapes.end(), shapeInItem) == vecShapes.end()) continue;
        vecItems.push_back(childItem);
    }
    if(vecItems.size()==0)
    {
        cout<<"mainTreeTools::getTreeItemsFromShapes()->____no items found in tree____"<<endl;
        return false;
    }
    else
    {
        cout<<"mainTreeTools::getTreeItemsFromShapes()->____nr. "<<vecItems.size()<<" items found____"<<endl;
        return true;
    }
}

//! -------------------------------------
//! function: formNewPart
//! details:  experimental. Not complete
//! -------------------------------------
#include "src/registeredMetatypes/topods_shape_reg.h"
#include "src/utils/topologytools.h"
void mainTreeTools::formNewPart(QTreeView *tree)
{
    cout<<"mainTreeTools::formNewPart()->____function called____"<<endl;

    QExtendedStandardItem *itemPart = new QExtendedStandardItem();
    QVariant data;
    data.setValue(QString("Part"));
    itemPart->setData(data,Qt::DisplayRole);
    int currentRow = tree->currentIndex().row();
    QModelIndex parent = tree->currentIndex().parent();

    QStandardItemModel *model = static_cast<QStandardItemModel*>(tree->model());
    model->insertRow(currentRow-1,parent);

    TopTools_ListOfShape shapeList;

    QList<QModelIndex> modelIndexes = tree->selectionModel()->selectedIndexes();
    for(QList<QModelIndex>::iterator it = modelIndexes.begin(); it!=modelIndexes.end(); ++it)
    {
        cout<<"____adding shape____"<<endl;
        QModelIndex curIndex = *it;
        SimulationNodeClass *curNode = curIndex.data(Qt::UserRole).value<SimulationNodeClass*>();
        TopoDS_Shape curShape = curNode->getPropertyValue<TopoDS_Shape_Reg>("Shape");
        shapeList.Append(curShape);
    }

}
//! -----------------------------------
//! function: getCurrentSimulationRoot
//! details:
//! -----------------------------------
QStandardItem* mainTreeTools::getCurrentSimulationRoot(QTreeView *treeView)
{
    QModelIndex curIndex = treeView->currentIndex();
    QStandardItemModel *treeModel = static_cast<QStandardItemModel*>(treeView->model());
    QStandardItem *curItem = treeModel->itemFromIndex(curIndex);
    SimulationNodeClass* node = curIndex.data(Qt::UserRole).value<SimulationNodeClass*>();

    if(node->isAnalysisRoot()) return curItem;
    if(node->isAnalysisSettings() || node->isSolution() || node->isSimulationSetUpNode()) return curItem->parent();
    if(node->isAnalysisResult() || node->isChildSimulationSetUpNode()) return curItem->parent()->parent();
    if(node->isNephewSimulationSetUpNode()) return curItem->parent()->parent()->parent();
    return Q_NULLPTR;
}

//! -----------------------------------------
//! function: resetSolutionInformation
//! details:  reset the solution information
//!           reset the discrete time map
//!           reset the solver output
//! -----------------------------------------
void mainTreeTools::resetSolutionInformation(SimulationNodeClass *nodeSolutionInformation)
{
    nodeSolutionInformation->getModel()->blockSignals(true);
    QVariant data;
    data.setValue(CCXSolverMessage());
    nodeSolutionInformation->replaceProperty("Solver output",Property("Solver output",data,Property::PropertyGroup_Hidden));
    QMap<double,QVector<int>> amap;
    data.setValue(amap);
    nodeSolutionInformation->replaceProperty("Discrete time map",Property("Discrete time map",data,Property::PropertyGroup_Hidden));

    QList<solutionInfo> solInfoList;
    data.setValue(solInfoList);
    nodeSolutionInformation->replaceProperty("Convergence data",Property("Convergence data",data,Property::PropertyGroup_Hidden));
    nodeSolutionInformation->getModel()->blockSignals(false);
}

//! ----------------------
//! function: addSolution
//! details:
//! ----------------------
void mainTreeTools::addSolution(QStandardItem *analysisRootItem)
{
    QVector<Property> props;
    QVariant data;

    //! -------------------------
    //! create the "Status" item
    //! -------------------------
    Property::solutionInformation theSolutionStatus = Property::solutionInformation_solveRequired;
    data.setValue(theSolutionStatus);
    props.push_back(Property("Status",data,Property::PropertyGroup_Information));

    //! --------------------
    //! "Project files dir"
    //! --------------------
    data.setValue(QString("undefined"));
    props.push_back(Property("Project files dir",data,Property::PropertyGroup_Information));

    SimulationNodeClass *nodeAnalysisRoot = analysisRootItem->data(Qt::UserRole).value<SimulationNodeClass*>();
    SimulationNodeClass::nodeType typeOfSolution = SimulationNodeClass::nodeType_NULL;
    switch(nodeAnalysisRoot->getType())
    {
    case SimulationNodeClass::nodeType_structuralAnalysis: typeOfSolution = SimulationNodeClass::nodeType_StructuralAnalysisSolution; break;
    case SimulationNodeClass::nodeType_thermalAnalysis: typeOfSolution = SimulationNodeClass::nodeType_thermalAnalysisSolution; break;
    case SimulationNodeClass::nodeType_combinedAnalysis: typeOfSolution = SimulationNodeClass::nodeType_combinedAnalysisSolution; break;
    case SimulationNodeClass::nodeType_CFDAnalysis: typeOfSolution = SimulationNodeClass::nodeType_CFDAnalysisSolution; break;
    case SimulationNodeClass::nodeType_particlesInFieldsAnalysis: typeOfSolution = SimulationNodeClass::nodeType_particlesInFieldsSolution; break;
    }

    //! ----------------
    //! create the node
    //! ----------------
    SimulationNodeClass *nodeSolution = new SimulationNodeClass("Solution",typeOfSolution,props);

    //! -----------------------------
    //! time tag and parent time tag
    //! -----------------------------
    nodeSolution->addTimeTag();
    QString parentTimeTag = analysisRootItem->data(Qt::UserRole).value<SimulationNodeClass*>()->getPropertyValue<QString>("Time tag");
    data.setValue(parentTimeTag);
    nodeSolution->addProperty(Property("Parent time tag",data,Property::PropertyGroup_Identifier));

    //! ----------------
    //! create the item
    //! ----------------
    QExtendedStandardItem *SolutionItem = new QExtendedStandardItem();
    data.setValue(nodeSolution);
    SolutionItem->setData(data,Qt::UserRole);
    SolutionItem->setData(nodeSolution->getName(),Qt::DisplayRole);

    analysisRootItem->appendRow(SolutionItem);
}

//! ---------------------------------
//! function: addSolutionInformation
//! details:
//! ---------------------------------
void mainTreeTools::addSolutionInformation(QStandardItem* solutionItem)
{
    QVariant data;
    QVector<Property> props;

    //! ------------------------------
    //! 0 => solver txt messages
    //! 1 => force convergence
    //! 2 => displacement convergence
    //! 3 => line search
    //! 4 => time step size
    //! ------------------------------
    data.setValue(int(0));
    Property prop_solutionInfoType("Solution information",data,Property::PropertyGroup_SolutionInfo);

    double updateInterval = 2.5;
    data.setValue(updateInterval);
    Property prop_updateInterval("Update interval",data,Property::PropertyGroup_SolutionInfo);

    props.push_back(prop_solutionInfoType);
    props.push_back(prop_updateInterval);

    //! ---------------------------
    //! hidden: solver text output
    //! ---------------------------
    CCXSolverMessage msg;
    data.setValue(msg);
    Property prop_solverOutput("Solver output",data,Property::PropertyGroup_Hidden);
    props.push_back(prop_solverOutput);

    //! ------------------
    //! discrete time map
    //! ------------------
    QMap<double,QVector<int>> dtm;
    data.setValue(dtm);
    Property prop_discreteTimeMap("Discrete time map",data,Property::PropertyGroup_Hidden);
    props.push_back(prop_discreteTimeMap);

    //! -----------------
    //! convergence data
    //! -----------------
    QList<solutionInfo> solInfoList;
    data.setValue(solInfoList);
    Property prop_solInfoList("Convergence data",data,Property::PropertyGroup_Hidden);
    props.push_back(prop_solInfoList);

    SimulationNodeClass *analysisRootNode = solutionItem->parent()->data(Qt::UserRole).value<SimulationNodeClass*>();
    SimulationNodeClass::nodeType analysisRootType = analysisRootNode->getType();
    SimulationNodeClass::nodeType typeOfSolutionInformation = SimulationNodeClass::nodeType_NULL;

    switch(analysisRootType)
    {
    case SimulationNodeClass::nodeType_structuralAnalysis: typeOfSolutionInformation = SimulationNodeClass::nodeType_StructuralAnalysisSolutionInformation; break;
    case SimulationNodeClass::nodeType_thermalAnalysis: typeOfSolutionInformation = SimulationNodeClass::nodeType_thermalAnalysisSolutionInformation; break;
    case SimulationNodeClass::nodeType_combinedAnalysis: typeOfSolutionInformation = SimulationNodeClass::nodeType_combinedAnalysisSolutionInformation; break;
    case SimulationNodeClass::nodeType_particlesInFieldsAnalysis: typeOfSolutionInformation = SimulationNodeClass::nodeType_particlesInFieldsSolutionInformation; break;
    }

    //! ----------------
    //! create the node
    //! ----------------
    SimulationNodeClass *nodeSolutionInformation = new SimulationNodeClass("Solution information",typeOfSolutionInformation,props);

    //! -----------------------------
    //! time tag and parent time tag
    //! -----------------------------
    nodeSolutionInformation->addTimeTag();
    SimulationNodeClass *nodeSolution = solutionItem->data(Qt::UserRole).value<SimulationNodeClass*>();
    QString parentTimeTag = nodeSolution->getPropertyValue<QString>("Time tag");
    data.setValue(parentTimeTag);
    nodeSolutionInformation->addProperty(Property("Parent time tag",data,Property::PropertyGroup_Identifier));

    //! -------------------------------------
    //! create the solution information item
    //! -------------------------------------
    QExtendedStandardItem *itemSolutionInformation = new QExtendedStandardItem();

    data.setValue(nodeSolutionInformation);
    itemSolutionInformation->setData(data,Qt::UserRole);
    data.setValue(nodeSolutionInformation->getName());
    itemSolutionInformation->setData(data,Qt::DisplayRole);
    solutionItem->appendRow(itemSolutionInformation);
}


//! ---------------------------------------
//! function: getAllBoundaryConditionsTags
//! details:
//! ---------------------------------------
void mainTreeTools::getAllBoundaryConditionsTags(QTreeView *tree, int type, std::vector<GeometryTag> &vecTags)
{
    cout<<"mainTreeTools::getAllBoundaryConditionsTags()->____function called____"<<endl;

    //! -----------
    //! model root
    //! -----------
    QStandardItemModel *model = (QStandardItemModel*)(tree->model());
    QStandardItem *RootItem = model->invisibleRootItem()->child(0,0);

    if(RootItem==Q_NULLPTR) return;
    for(int i=0; i<RootItem->rowCount();i++)
    {
        QStandardItem *item = RootItem->child(i,0);
        SimulationNodeClass *node = item->data(Qt::UserRole).value<SimulationNodeClass*>();
        //! -----------------------
        //! connections root found
        //! -----------------------
        if(node->getType()==SimulationNodeClass::nodeType_connection)
        {
            QStandardItem *Connections_RootItem = item;
            for(int i=0; i<Connections_RootItem->rowCount(); i++)
            {
                QStandardItem *connectionGroup = Connections_RootItem->child(i,0);
                for(int k=0; k<connectionGroup->rowCount(); k++)
                {
                    QStandardItem *connectionItem = connectionGroup->child(k,0);
                    SimulationNodeClass *node = connectionItem->data(Qt::UserRole).value<SimulationNodeClass*>();
                    std::vector<GeometryTag> masterTags = node->getPropertyValue<std::vector<GeometryTag>>("Tags master");
                    std::vector<GeometryTag> slaveTags = node->getPropertyValue<std::vector<GeometryTag>>("Tags slave");
                    for(int m = 0; m<masterTags.size(); m++)
                    {
                        int shapeType = int(masterTags[m].subShapeType);
                        if(shapeType!=type) continue;
                        if(std::find(masterTags.begin(),masterTags.end(),masterTags[m])==vecTags.end()) vecTags.push_back(masterTags[m]);
                    }
                    for(int m=0; m<slaveTags.size(); m++)
                    {
                        int shapeType = int(slaveTags[m].subShapeType);
                        if(shapeType!=type) continue;
                        if(std::find(vecTags.begin(),vecTags.end(),slaveTags[m])==vecTags.end()) vecTags.push_back(slaveTags[m]);
                    }
                }
            }
        }

        //! --------------------
        //! analysis root found
        //! --------------------
        if(node->isAnalysisRoot())
        {
            for(int j=1; j<item->rowCount();j++)
            {
                QStandardItem *itemBC = item->child(j,0);
                SimulationNodeClass *nodeBC = itemBC->data(Qt::UserRole).value<SimulationNodeClass*>();
                if(nodeBC->isSimulationSetUpNode()==false || nodeBC->isChildSimulationSetUpNode() || nodeBC->isNephewSimulationSetUpNode()) continue;
                if(nodeBC->getPropertyItem("Tags")==Q_NULLPTR) continue;
                std::vector<GeometryTag> BCTags = nodeBC->getPropertyValue<std::vector<GeometryTag>>("Tags");
                for(int m=0; m<BCTags.size(); m++)
                {
                    if(std::find(vecTags.begin(),vecTags.end(),BCTags[m])==vecTags.end())
                    {
                        int shapeType = int(BCTags[m].subShapeType);
                        if(shapeType!=type) continue;
                        vecTags.push_back(BCTags[m]);
                    }
                }
            }
        }

        if(!node->isAnalysisRoot() || node->getType()==SimulationNodeClass::nodeType_connection) continue;
    /* NO need to keep mesh control boundary
        if(node->getType()==SimulationNodeClass::nodeType_meshControl)
        {
            cout<<"____MESH CONTROL____"<<endl;
            for(int j=0; j<item->rowCount();j++)
            {
                QStandardItem *itemMeshControl = item->child(j,0);
                SimulationNodeClass *nodeMeshControl = itemMeshControl->data(Qt::UserRole).value<SimulationNodeClass*>();

                if(nodeMeshControl->getType()!=SimulationNodeClass::nodeType_meshFaceSize)
                    continue;
                std::vector<GeometryTag> faceSizingTags = nodeMeshControl->getPropertyValue<std::vector<GeometryTag>>("Tags");
                for(int m=0; m<faceSizingTags.size(); m++)
                {
                    if(std::find(vecTags.begin(),vecTags.end(),faceSizingTags[m])==vecTags.end())
                    {
                        int shapeType = int(faceSizingTags[m].subShapeType);
                        if(shapeType!=type) continue;
                        vecTags.push_back(faceSizingTags[m]);
                    }
                }
            }
        }
    */
    }
    cout<<"mainTreeTools::getAllBoundaryConditionsTags()->____function exiting____"<<endl;
}

//! -----------------------------------
//! function: getAnalysisNodeFromIndex
//! details:
//! -----------------------------------
SimulationNodeClass* mainTreeTools::getAnalysisSettingsNodeFromIndex(QModelIndex curIndex)
{
    if(curIndex.isValid()==false) return Q_NULLPTR;
    SimulationNodeClass *nodeAnalysisSettings = Q_NULLPTR;
    SimulationNodeClass *curNode = curIndex.data(Qt::UserRole).value<SimulationNodeClass*>();
    if(curNode->isAnalysisRoot()) nodeAnalysisSettings = curIndex.child(0,0).data(Qt::UserRole).value<SimulationNodeClass*>();
    if(curNode->isAnalysisSettings()) nodeAnalysisSettings = curIndex.data(Qt::UserRole).value<SimulationNodeClass*>();
    if(curNode->isSimulationSetUpNode()) nodeAnalysisSettings = curIndex.parent().child(0,0).data(Qt::UserRole).value<SimulationNodeClass*>();
    if(curNode->isSolution()) nodeAnalysisSettings = curIndex.parent().child(0,0).data(Qt::UserRole).value<SimulationNodeClass*>();
    if(curNode->isSolutionInformation()) nodeAnalysisSettings = curIndex.parent().parent().child(0,0).data(Qt::UserRole).value<SimulationNodeClass*>();
    if(curNode->isAnalysisResult()) nodeAnalysisSettings = curIndex.parent().parent().child(0,0).data(Qt::UserRole).value<SimulationNodeClass*>();
    if(curNode->isChildSimulationSetUpNode()) nodeAnalysisSettings = curIndex.parent().parent().child(0,0).data(Qt::UserRole).value<SimulationNodeClass*>();
    if(curNode->isNephewSimulationSetUpNode()) nodeAnalysisSettings = curIndex.parent().parent().parent().child(0,0).data(Qt::UserRole).value<SimulationNodeClass*>();

    return nodeAnalysisSettings;
}

//! -------------------------------------------------
//! function: getAnalysisSettingsNodeFromCurrentItem
//! details:
//! -------------------------------------------------
QStandardItem* mainTreeTools::getAnalysisSettingsItemFromCurrentItem(QTreeView *treeView)
{
    cout<<"mainTreeTools::getAnalysisSettingsNodeFromCurrentItem()->____function called____"<<endl;

    QModelIndex currentModelIndex = treeView->currentIndex();
    QStandardItem *curItem = static_cast<QStandardItemModel*>(treeView->model())->itemFromIndex(currentModelIndex);
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

    //! -----------------------------------------------------------------------------------
    //! case 3: the current item is a post processing item or a child of a simulation node
    //! -----------------------------------------------------------------------------------
    if(curNode->isAnalysisResult() || curNode->isSolutionInformation() || curNode->isChildSimulationSetUpNode())
    {
        QStandardItem *item = curItem->parent()->parent()->child(0,0);
        return static_cast<QExtendedStandardItem*>(item);
    }

    //! ----------------------------------------------------------
    //! case 4: the current item is a nephew of a simulation node
    //! ----------------------------------------------------------
    if(curNode->isNephewSimulationSetUpNode())
    {
        QStandardItem *item = curItem->parent()->parent()->parent()->child(0,0);
        return static_cast<QExtendedStandardItem*>(item);
    }
    return Q_NULLPTR;
}

//! -------------------------------------------------
//! function: getAnalysisSettingsItemFromCurrentItem
//! details:
//! -------------------------------------------------
SimulationNodeClass* mainTreeTools::getAnalysisSettingsNodeFromCurrentItem(QTreeView *treeView)
{
    QStandardItem *item = mainTreeTools::getAnalysisSettingsItemFromCurrentItem(treeView);
    SimulationNodeClass *node = item->data(Qt::UserRole).value<SimulationNodeClass*>();
    return node;
}

//! --------------------------------------------
//! function: getTreeItem
//! details:  already used in SimulationManager
//! --------------------------------------------
QExtendedStandardItem* mainTreeTools::getTreeItem(QStandardItemModel *model, SimulationNodeClass::nodeType theNodeType)
{
    if(model==Q_NULLPTR)
    {
        cerr<<"SimulationManager::getTreeItem()->____NULL model____"<<endl;
        return Q_NULLPTR;
    }

    //! retrieve the root nodes
    QList<QStandardItem*> items;
    mainTreeTools::getTreeItemsRecursively(model,items);
    for(QList<QStandardItem*>::iterator it = items.begin(); it!=items.end(); it++)
    {
        QStandardItem* curItem = *it;
        SimulationNodeClass *curNode = curItem->data(Qt::UserRole).value<SimulationNodeClass*>();

        if(curNode==Q_NULLPTR) return Q_NULLPTR;

        SimulationNodeClass::nodeType curNodeType = curNode->getType();
        if(curNodeType==theNodeType)
        {
            return static_cast<QExtendedStandardItem*>(curItem);
        }
    }
    return Q_NULLPTR;
}


//! ----------------------------------------------------------------
//! function: ItemFromScope
//! details:  for a given shape in a geometry item, return the item
//! ----------------------------------------------------------------
QExtendedStandardItem* mainTreeTools::ItemFromScope(QStandardItemModel *model,const TopoDS_Shape &aShape)
{
    QStandardItem *theGeometryRoot=mainTreeTools::getTreeItem(model,SimulationNodeClass::nodeType_geometry);
    int N = theGeometryRoot->rowCount();
    for(int k=0; k<N;k++)
    {
        QStandardItem *aGeometryItem = theGeometryRoot->child(k,0);
        SimulationNodeClass *aNode = aGeometryItem->data(Qt::UserRole).value<SimulationNodeClass*>();
        if(aNode->getType()==SimulationNodeClass::nodeType_pointMass) continue;
        //int mapIndex = aNode->getPropertyValue<int>("Map index");
        //TopoDS_Shape theShape = myDB->bodyMap.value(mapIndex);
        //if(theShape==aShape) return static_cast<QExtendedStandardItem*>(aGeometryItem);
        TopoDS_Shape shapeInItem = aNode->getPropertyValue<TopoDS_Shape>("Shape");
        if(shapeInItem == aShape)
        return static_cast<QExtendedStandardItem*>(aGeometryItem);

    }
    return Q_NULLPTR;
}

//! -------------------------------
//! function: getAllTreeItemOfType
//! details:
//! -------------------------------
QList<QStandardItem*> mainTreeTools::getAllTreeItemOfType(QStandardItemModel *model, SimulationNodeClass::nodeType theNodeType)
{
    if(model==Q_NULLPTR) return QList<QStandardItem*>();
    QList<QStandardItem*> items, itemsout;
    mainTreeTools::getTreeItemsRecursively(model,items);
    for(QList<QStandardItem*>::iterator it = items.begin(); it!=items.end(); ++it)
    {
        QStandardItem* curItem = *it;
        if(curItem->data(Qt::UserRole).value<SimulationNodeClass*>()->getType()==theNodeType) itemsout.append(curItem);
    }
    return itemsout;
}

//! ------------------------------------------------------------------------
//! function: getInsertionRow
//! details:  when adding a simulation setup item, that item must be placed
//!           between the "Analysis settings" item and the "Solution" item.
//!           The function finds the right row for inserting the new item
//! ------------------------------------------------------------------------
const int mainTreeTools::getInsertionRow(QTreeView *tree)
{
    QModelIndex theCurIndex = tree->currentIndex();
    QStandardItem *theCurItem = static_cast<QStandardItemModel*>(tree->model())->itemFromIndex(theCurIndex);
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
    //cout<<"____insertion row: "<<insertionRow<<"____"<<endl;
    return insertionRow;
}
