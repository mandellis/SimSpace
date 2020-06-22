//! ----------------
//! custom includes
//! ----------------
#include "maintreetools.h"
#include "simulationnodeclass.h"
#include "tabulardatacolumns.h"
#include "qextendedstandarditem.h"
#include "ccxsolvermessage.h"
#include "solutioninfo.h"

//! ---
//! Qt
//! ---
#include <QStandardItemModel>

#include <memory>

//! --------------------------
//! function: getColumnToRead
//! details:
//! --------------------------
QList<int> mainTreeTools::getColumnsToRead(QTreeView *tree)
{
    cout<<"mainTreeTools::getColumnsToRead()->____function called____"<<endl;

    int SC = mainTreeTools::calculateStartColumn(tree);
    QModelIndex currentIndex=tree->currentIndex();

    QList<int> theColumnsToShow;
    SimulationNodeClass *aNode = currentIndex.data(Qt::UserRole).value<SimulationNodeClass*>();
    SimulationNodeClass::nodeType theType = aNode->getType();

    if(aNode->getPropertyItem("Define by")!=Q_NULLPTR)
    {
        //! ------------------------------------
        //! items with the "Define by" property
        //! ------------------------------------
        if(theType == SimulationNodeClass::nodeType_structuralAnalysisBoltPretension)
        {
            theColumnsToShow<<SC<<SC+1<<SC+2;
        }
        if(theType!= SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_CylindricalSupport ||
                theType!= SimulationNodeClass::nodeType_structuralAnalysisBoundaryContidion_FixedSupport ||
                theType!= SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_FrictionlessSupport)
        {
            Property::defineBy theDefineBy = aNode->getPropertyValue<Property::defineBy>("Define by");

            if(theDefineBy==Property::defineBy_components)
            {
                SimulationNodeClass::nodeType theType = aNode->getType();

                //! -------------------------------------------------------------------
                //! the displacement is a special case, since it has the "free" option
                //! -------------------------------------------------------------------
                if(theType==SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Displacement ||
                        theType==SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteDisplacement ||
                        theType==SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteRotation)
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
                else
                {
                    //! ---------------------------
                    //! read all the three columns
                    //! ---------------------------
                    theColumnsToShow<<SC<<SC+1<<SC+2;
                }
            }
            else
            {
                //! ------------------------------------------
                //! definition through a scalar (1 component)
                //! ------------------------------------------
                theColumnsToShow<<SC;
            }
        }
        else
        {
            //! no column to read for supports
        }
    }
    else
    {
        //! ------------------------------------------------------------
        //! items having tabular data, without the "Define by" property
        //! ------------------------------------------------------------
        switch(theType)
        {
        case SimulationNodeClass::nodeType_modelChange:
        case SimulationNodeClass::nodeType_structuralAnalysisThermalCondition:
        case SimulationNodeClass::nodeType_thermalAnalysisThermalFlow:
        case SimulationNodeClass::nodeType_thermalAnalysisThermalFlux:
        case SimulationNodeClass::nodeType_thermalAnalysisThermalPower:

        {
            theColumnsToShow<<SC;
        }
            break;
        case SimulationNodeClass::nodeType_thermalAnalysisConvection:
        {
            theColumnsToShow<<SC<<SC+1;
        }
            break;
        //! ------------------------------------------------------------
        //! items without tabular data and without "Define by" property
        //! default case left for documentation
        //! ------------------------------------------------------------
        default:
        {
            ;
        }
            break;
        }
    }

    //! ----------------------------
    //! diagnostic - can be removed
    //! ----------------------------
    cout<<"mainTreeTools::getColumnsToRead()->____columns to read: {";
    int i; for(i=0;i<theColumnsToShow.length()-1;i++) cout<<theColumnsToShow.at(i)<<",";
    cout<<theColumnsToShow.at(i)<<"}"<<endl;
    //! ---------------
    //! end diagnostic
    //! ---------------
    return theColumnsToShow;
}

//! --------------------------
//! function: getColumnToRead
//! details:
//! --------------------------
QList<int> mainTreeTools::getColumnsToRead(QExtendedStandardItem *anItem)
{
    cout<<"mainTreeTools::getColumnsToRead()->____function called____"<<endl;

    int SC = mainTreeTools::calculateStartColumn(anItem);
    QModelIndex currentIndex=anItem->index();

    QList<int> theColumnsToShow;
    SimulationNodeClass *aNode = currentIndex.data(Qt::UserRole).value<SimulationNodeClass*>();
    SimulationNodeClass::nodeType theType = aNode->getType();

    //! ------------------------------------------------
    //! these nodes do not have the "defineBy" property
    //! ------------------------------------------------
    if(aNode->getPropertyItem("Define by")!=Q_NULLPTR)
    {
        if(theType == SimulationNodeClass::nodeType_structuralAnalysisBoltPretension)
        {
            theColumnsToShow<<SC<<SC+1<<SC+2;
        }
        if(theType!= SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_CylindricalSupport ||
                theType!= SimulationNodeClass::nodeType_structuralAnalysisBoundaryContidion_FixedSupport ||
                theType!= SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_FrictionlessSupport)
        {
            Property::defineBy theDefineBy = aNode->getPropertyValue<Property::defineBy>("Define by");
            if(theDefineBy==Property::defineBy_components)
            {
                SimulationNodeClass::nodeType theType = aNode->getType();

                //! -------------------------------------------------------------------
                //! the displacement is a special case, since it has the "free" option
                //! -------------------------------------------------------------------
                if(theType==SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Displacement ||
                        theType==SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteDisplacement ||
                        theType==SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteRotation)
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
                else
                {
                    //! ---------------------------
                    //! read all the three columns
                    //! ---------------------------
                    theColumnsToShow<<SC<<SC+1<<SC+2;
                }
            }
            else
            {
                //! ------------------------------------------
                //! definition through a scalar (1 component)
                //! ------------------------------------------
                theColumnsToShow<<SC;
            }
        }
        else
        {
            //! no column to read for supports
        }
    }
    else
    {
        switch(theType)
        {
        case SimulationNodeClass::nodeType_modelChange:
        case SimulationNodeClass::nodeType_structuralAnalysisThermalCondition:
        case SimulationNodeClass::nodeType_thermalAnalysisTemperature:
        case SimulationNodeClass::nodeType_thermalAnalysisThermalFlow:
        case SimulationNodeClass::nodeType_thermalAnalysisThermalFlux:
        case SimulationNodeClass::nodeType_thermalAnalysisThermalPower:

        {
            theColumnsToShow<<SC;
        }
            break;
        case SimulationNodeClass::nodeType_thermalAnalysisConvection:
        {
            theColumnsToShow<<SC<<SC+1;
        }
            break;

        //! ------------------------------------------------------------
        //! items without tabular data and without "Define by" property
        //! default case left here for documentation
        //! ------------------------------------------------------------
        default:
        {
            ;
        }
            break;
        }
    }

    //! ----------------------------
    //! diagnostic - can be removed
    //! ----------------------------
    cout<<"mainTreeTools::getColumnsToRead()->____columns to read: {";
    int i; for(i=0;i<theColumnsToShow.length()-1;i++) cout<<theColumnsToShow.at(i)<<",";
    cout<<theColumnsToShow.at(i)<<"}"<<endl;
    //! ---------------
    //! end diagnostic
    //! ---------------
    return theColumnsToShow;
}

//! -------------------------------
//! function: calculateStartColumn
//! details:
//! -------------------------------
int mainTreeTools::calculateStartColumn(QExtendedStandardItem *anItem)
{
    //cout<<"mainTreeTools::calculateStartColumn()->____function called____"<<endl;
    //! -----------------
    //! the start column
    //! -----------------
    int startColumn = 0;

    QModelIndex currentIndex=anItem->index();
    //! --------------------------------
    //! the row of the item in the tree
    //! --------------------------------
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
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_ImportedTemperatureDistribution:
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_CompressionOnlySupport:
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_CylindricalSupport:
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_FrictionlessSupport:
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryContidion_FixedSupport:
        case SimulationNodeClass::nodeType_thermalAnalysisAdiabaticWall:
#ifdef COSTAMP_VERSION
        case SimulationNodeClass::nodeType_timeStepBuilder:
#endif
        case SimulationNodeClass::nodeType_mapper:
            delta = -1;
            break;

        case SimulationNodeClass::nodeType_modelChange:
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Pressure:
        case SimulationNodeClass::nodeType_structuralAnalysisThermalCondition:
        case SimulationNodeClass::nodeType_thermalAnalysisRadiation:
        case SimulationNodeClass::nodeType_thermalAnalysisTemperature:
        case SimulationNodeClass::nodeType_thermalAnalysisThermalFlow:
        case SimulationNodeClass::nodeType_thermalAnalysisThermalPower:
        case SimulationNodeClass::nodeType_thermalAnalysisThermalFlux:
            delta = 0;
            break;
        case SimulationNodeClass::nodeType_thermalAnalysisConvection:
            delta = 1;
            break;
        case SimulationNodeClass::nodeType_magneticField:
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Force:
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Moment:
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteForce:
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Acceleration:
            //! [...] add here other items supporting a vectorial definition
            theDefineBy = curNode->getPropertyValue<Property::defineBy>("Define by");
            if(theDefineBy==Property::defineBy_vector) delta = 0;
            else delta = 2;
            break;

            //! ----------------
            //! bolt pretension
            //! ----------------
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

    //! --------------------------------
    //! the row of the item in the tree
    //! --------------------------------
    int rrow = currentIndex.row();

    //! ----------------------------------
    //! retrieve the simulation root item
    //! ----------------------------------
    QModelIndex index = currentIndex.parent();
    QStandardItemModel *treeModel = static_cast<QStandardItemModel*>(tree->model());
    QStandardItem* itemSimulationRoot = treeModel->itemFromIndex(index);

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
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_ImportedTemperatureDistribution:
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_CylindricalSupport:
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_CompressionOnlySupport:
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_FrictionlessSupport:
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryContidion_FixedSupport:
        case SimulationNodeClass::nodeType_thermalAnalysisAdiabaticWall:
#ifdef COSTAMP_VERSION
        case SimulationNodeClass::nodeType_timeStepBuilder:
#endif
        case SimulationNodeClass::nodeType_mapper:
            delta = -1;
            break;

        case SimulationNodeClass::nodeType_modelChange:
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Pressure:
        case SimulationNodeClass::nodeType_structuralAnalysisThermalCondition:
        case SimulationNodeClass::nodeType_thermalAnalysisRadiation:
        case SimulationNodeClass::nodeType_thermalAnalysisTemperature:
        case SimulationNodeClass::nodeType_thermalAnalysisThermalFlow:
        case SimulationNodeClass::nodeType_thermalAnalysisThermalPower:
        case SimulationNodeClass::nodeType_thermalAnalysisThermalFlux:
            delta = 0;
            break;
        case SimulationNodeClass::nodeType_thermalAnalysisConvection:
            delta = 1;
            break;
        case SimulationNodeClass::nodeType_magneticField:
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Force:
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Moment:
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteForce:
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Acceleration:
            //! [...] add here other items supporting a vectorial definition
            theDefineBy = curNode->getPropertyValue<Property::defineBy>("Define by");
            if(theDefineBy==Property::defineBy_vector) delta = 0;
            else delta = 2;
            break;

            //! ----------------
            //! bolt pretension
            //! ----------------
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
        // here is your applicable code
        QExtendedStandardItem *item = static_cast<QExtendedStandardItem*>(model->itemFromIndex(index));
        items.push_back(item);
        if(model->hasChildren(index))
        {
            getTreeItemsRecursively(model, items, index);
        }
    }
}

//! -------------------------------------
//! function: formNewPart
//! details:  experimental. Not complete
//! -------------------------------------
#include "topods_shape_reg.h"
#include "topologytools.h"
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
    if(node->isAnalysisResult()) return curItem->parent()->parent();
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
void mainTreeTools::getAllBoundaryConditionsTags(QTreeView *tree, int type, QVector<GeometryTag> &vecTags)
{
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
                    QVector<GeometryTag> masterTags = node->getPropertyValue<QVector<GeometryTag>>("Tags master");
                    QVector<GeometryTag> slaveTags = node->getPropertyValue<QVector<GeometryTag>>("Tags slave");
                    for(int m = 0; m<masterTags.length(); m++)
                    {
                        int shapeType = int(masterTags[m].subShapeType);
                        if(shapeType!=type) continue;
                        if(vecTags.contains(masterTags[m])==false) vecTags.push_back(masterTags[m]);
                    }
                    for(int m=0; m<slaveTags.length(); m++)
                    {
                        int shapeType = int(slaveTags[m].subShapeType);
                        if(shapeType!=type) continue;
                        if(vecTags.contains(slaveTags[m])==false) vecTags.push_back(slaveTags[m]);
                    }
                }
            }
        }

        //! --------------------
        //! analysis root found
        //! --------------------
        if(node->isAnalysisRoot())
        {
            for(int j=0; j<item->rowCount();j++)
            {
                QStandardItem *itemBC = item->child(j,0);
                SimulationNodeClass *nodeBC = itemBC->data(Qt::UserRole).value<SimulationNodeClass*>();
                if(nodeBC->isSimulationSetUpNode()==false) continue;
                QVector<GeometryTag> BCTags = nodeBC->getPropertyValue<QVector<GeometryTag>>("Tags");
                for(int m=0; m<BCTags.length(); m++)
                {
                    if(vecTags.contains(BCTags[m])==false)
                    {
                        int shapeType = int(BCTags[m].subShapeType);
                        if(shapeType!=type) continue;
                        vecTags.push_back(BCTags[m]);
                    }
                }
            }
        }

        if(node->getType()==SimulationNodeClass::nodeType_meshControl)
        {
            for(int j=0; j<item->rowCount();j++)
            {
                QStandardItem *itemMeshControl = item->child(j,0);
                SimulationNodeClass *nodeFaceMeshControl = itemMeshControl->data(Qt::UserRole).value<SimulationNodeClass*>();
                if(nodeFaceMeshControl->getType()!=SimulationNodeClass::nodeType_meshFaceSize) continue;
                QVector<GeometryTag> faceSizingTags = nodeFaceMeshControl->getPropertyValue<QVector<GeometryTag>>("Tags");
                for(int m=0; m<faceSizingTags.length(); m++)
                {
                    if(vecTags.contains(faceSizingTags[m])==false)
                    {
                        int shapeType = int(faceSizingTags[m].subShapeType);
                        if(shapeType!=type) continue;
                        vecTags.push_back(faceSizingTags[m]);
                    }
                }
            }
        }
    }
}
