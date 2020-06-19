//! ----------------
//! custom includes
//! ----------------
#include "maintreetools.h"
#include "simulationnodeclass.h"
#include "tabulardatacolumns.h"
#include "qextendedstandarditem.h"
#include "topods_shape_reg.h"
#include "topologytools.h"

//! ---
//! Qt
//! ---
#include <QStandardItemModel>

//! ---------------------------
//! function: getColumnsToRead
//! details:
//! ---------------------------
QList<int> mainTreeTools::getColumnsToRead(QTreeView *tree)
{
    QList<int> listOfColumns;

    //! the starting column
    int startColumn = 0;

    //! the row of the item in the tree
    int rrow = tree->currentIndex().row();

    //! --------------------------------------
    //! retrieve the "Static Structural" item
    //! --------------------------------------
    QModelIndex index = tree->currentIndex().parent();
    QStandardItemModel *treeModel = static_cast<QStandardItemModel*>(tree->model());
    QStandardItem* itemSimulationRoot = treeModel->itemFromIndex(index);

    //! --------------------------
    //! retrieve the current item
    //! --------------------------
    QStandardItem *curItem = itemSimulationRoot->child(rrow,0);
    SimulationNodeClass *curNode = curItem->data(Qt::UserRole).value<SimulationNodeClass*>();

    //! the for cycle starts from the second item (the previous is "Analysis settings")
    int offset = 0;
    int k=1;
    do
    {
        QStandardItem *curItem = itemSimulationRoot->child(k,0);
        SimulationNodeClass *curNode = curItem->data(Qt::UserRole).value<SimulationNodeClass*>();
        SimulationNodeClass::nodeType curType = curNode->getType();
        Property::defineBy theDefineBy;

        int delta = 0;
        switch (curType)
        {
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_CylindricalSupport:
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_FrictionlessSupport:
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryContidion_FixedSupport:
        case SimulationNodeClass::nodeType_mapper:
        case SimulationNodeClass::nodeType_thermalAnalysisAdiabaticWall:
#ifdef COSTAMP_VERSION
        case SimulationNodeClass::nodeType_timeStepBuilder:
            #endif
                delta = -1;
            break;

        case SimulationNodeClass::nodeType_modelChange:
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Pressure:
        case SimulationNodeClass::nodeType_structuralAnalysisThermalCondition:
        case SimulationNodeClass::nodeType_thermalAnalysisTemperature:
        case SimulationNodeClass::nodeType_thermalAnalysisThermalFlux:
        case SimulationNodeClass::nodeType_thermalAnalysisThermalFlow:
        case SimulationNodeClass::nodeType_thermalAnalysisConvection:
        case SimulationNodeClass::nodeType_thermalAnalysisThermalPower:
            delta = 0;
            break;

        case SimulationNodeClass::nodeType_structuralAnalysisBoltPretension:
            delta = 2;
            break;

        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Force:
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteForce:
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Moment:
            theDefineBy = curNode->getPropertyValue<Property::defineBy>("Define by");
            if(theDefineBy==Property::defineBy_vector) delta = 0;
            else delta = 2;
            break;

            //! ---------------------------------------------------------------------------------------
            //! the "Displacement" is a special case, since its component could have the option "Free"
            //! ---------------------------------------------------------------------------------------
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Displacement:
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteDisplacement:
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteRotation:
            theDefineBy = curNode->getPropertyValue<Property::defineBy>("Define by");
            if(theDefineBy==Property::defineBy_vector || theDefineBy==Property::defineBy_normal) delta=0;
            else
            {
                delta = 0;
                int subDelta = 0;
                Property::loadDefinition loadDefinitionYcomponent = curNode->getPropertyValue<Property::loadDefinition>("X component");
                Property::loadDefinition loadDefinitionXcomponent = curNode->getPropertyValue<Property::loadDefinition>("Y component");
                Property::loadDefinition loadDefinitionZcomponent = curNode->getPropertyValue<Property::loadDefinition>("Z component");
                if(loadDefinitionXcomponent!=Property::loadDefinition_free)subDelta++;
                if(loadDefinitionYcomponent!=Property::loadDefinition_free)subDelta++;
                if(loadDefinitionZcomponent!=Property::loadDefinition_free)subDelta++;
                delta = subDelta-1;
            }
            break;
        }
        offset = offset + delta;
        k++;
    }
    while(k<rrow-1);

    //! --------------------------------------------------
    //! number of columns defining the "Analysis setting"
    //! --------------------------------------------------
    int initNumberOfColumns = NUMBER_OF_COLUMNS_BEFORE_BC_DATA;
    startColumn = (rrow+initNumberOfColumns)+offset;

    //SimulationNodeClass* theNode = anItem->data(Qt::UserRole).value<SimulationNodeClass*>();
    SimulationNodeClass::nodeType theNodeType = curNode->getType();
    switch(theNodeType)
    {
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Displacement:
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteDisplacement:
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteRotation:
    {
        Property::defineBy theDefineBy = curNode->getPropertyItem("Define by")->data(Qt::UserRole).value<Property>().getData().value<Property::defineBy>();
        if(theDefineBy==Property::defineBy_components)
        {
            QExtendedStandardItem* Xcomp = curNode->getPropertyItem("X component");
            QExtendedStandardItem* Ycomp = curNode->getPropertyItem("Y component");
            QExtendedStandardItem* Zcomp = curNode->getPropertyItem("Z component");
            Property::loadDefinition displXdefinition = Xcomp->data(Qt::UserRole).value<Property>().getData().value<Property::loadDefinition>();
            Property::loadDefinition displYdefinition = Ycomp->data(Qt::UserRole).value<Property>().getData().value<Property::loadDefinition>();
            Property::loadDefinition displZdefinition = Zcomp->data(Qt::UserRole).value<Property>().getData().value<Property::loadDefinition>();
            bool b0 = displXdefinition!=Property::loadDefinition_free? true:false;
            bool b1 = displYdefinition!=Property::loadDefinition_free? true:false;
            bool b2 = displZdefinition!=Property::loadDefinition_free? true:false;

            if((b0==true && b1==false && b2==false)||
                    (b0==false && b1==true && b2==false)||
                    (b0==false && b1==false && b2==true))
            {
                listOfColumns<<startColumn;
            }
            else if((b0==true && b1 == true && b2 == false) ||
                    (b0==true && b1 == false && b2 == true) ||
                    (b0==false && b1 == true && b2 == true))
            {
                listOfColumns<<startColumn<<startColumn+1;
            }
            else if(b0==true && b1==true && b2==true)
            {
                Property::defineBy theDefineBy = curNode->getPropertyValue<Property::defineBy>("Define by");
                if(theDefineBy==Property::defineBy_components) listOfColumns<<startColumn<<startColumn+1<<startColumn+2;
                else listOfColumns<<startColumn;
            }
        }
        else listOfColumns<<startColumn;
    }
        break;
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Force:
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteForce:
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Moment:
    {
        Property::defineBy theDefineBy = curNode->getPropertyValue<Property::defineBy>("Define by");
        if(theDefineBy==Property::defineBy_components)
        {
            listOfColumns<<startColumn<<startColumn+1<<startColumn+2;
        }
        else listOfColumns<<startColumn;
    }
        break;
    case SimulationNodeClass::nodeType_structuralAnalysisBoltPretension:
        listOfColumns<<startColumn<<startColumn+1<<startColumn+2;
        break;
    default:
        listOfColumns<<startColumn;
        break;
    }
    return listOfColumns;
}

//! --------------------------
//! function: getColumnToRead
//! details:
//! --------------------------
QList<int> mainTreeTools::getColumnsToRead1(QTreeView *tree)
{
    cout<<"mainTreeTools::getColumnsToRead1()->____function called____"<<endl;

    int SC = mainTreeTools::calculateStartColumn(tree);
    cout<<"mainTreeTools::getColumnsToRead1()->____start column: "<<SC<<"____"<<endl;
    QList<int> theColumnsToShow;
    SimulationNodeClass *aNode = tree->currentIndex().data(Qt::UserRole).value<SimulationNodeClass*>();
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
            if(aNode->getPropertyItem("Define by")==Q_NULLPTR)
            {
                cerr<<"____Error: defineBy not existing => now crashing____"<<endl;
            }
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
        case SimulationNodeClass::nodeType_thermalAnalysisThermalFlux:
        case SimulationNodeClass::nodeType_thermalAnalysisThermalFlow:
        case SimulationNodeClass::nodeType_thermalAnalysisConvection:
        case SimulationNodeClass::nodeType_thermalAnalysisThermalPower:
        {
            cout<<"\\---------------------------------------------------------------------------\\"<<endl;
            cout<<"\\- mainTreeTools1::getColumnsToRead1()-> calling the function for tree item \\"<<endl;
            cout<<"\\- without \"Define by\" and with tabular data (i.e. Thermal condition)     \\"<<endl;
            cout<<"\\---------------------------------------------------------------------------\\"<<endl;
            theColumnsToShow<<SC;
        }
            break;

        default:
        {
            //! no tabular data => no columns to read
            cout<<"\\---------------------------------------------------------------------------\\"<<endl;
            cout<<"\\- mainTreeTools1::getColumnsToRead1()-> no tabular data: no column to read \\"<<endl;
            cout<<"\\---------------------------------------------------------------------------\\"<<endl;
        }
            break;
        }
    }

    //! ----------------------------
    //! diagnostic - can be removed
    //! ----------------------------
    cout<<"mainTreeTools::getColumnsToRead1()->____columns to read: {";
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
int mainTreeTools::calculateStartColumn(QTreeView *tree)
{
    //cout<<"mainTreeTools::calculateStartColumn()->____function called____"<<endl;

    //! -----------------
    //! the start column
    //! -----------------
    int startColumn = 0;

    //! --------------------------------
    //! the row of the item in the tree
    //! --------------------------------
    int rrow = tree->currentIndex().row();

    //! ----------------------------------
    //! retrieve the simulation root item
    //! ----------------------------------
    QModelIndex index = tree->currentIndex().parent();
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
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_CylindricalSupport:
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_FrictionlessSupport:
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryContidion_FixedSupport:
#ifdef COSTAMP_VERSION
        case SimulationNodeClass::nodeType_timeStepBuilder:
#endif
        case SimulationNodeClass::nodeType_mapper:
        case SimulationNodeClass::nodeType_thermalAnalysisAdiabaticWall:
            delta = -1;
            break;

        case SimulationNodeClass::nodeType_modelChange:
        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Pressure:
        case SimulationNodeClass::nodeType_structuralAnalysisThermalCondition:
        case SimulationNodeClass::nodeType_thermalAnalysisTemperature:
        case SimulationNodeClass::nodeType_thermalAnalysisThermalFlux:
        case SimulationNodeClass::nodeType_thermalAnalysisThermalFlow:
        case SimulationNodeClass::nodeType_thermalAnalysisConvection:
        case SimulationNodeClass::nodeType_thermalAnalysisThermalPower:
            delta = 0;
            break;

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
    //cout<<"mainTreeTools::calculateStartColumn()->____exiting function____"<<endl;
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

//! ------------------------------------
//! function: getTreeItemsRecursively
//! details:  scan recursively the tree
//! ------------------------------------
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
    SimulationNodeClass* node = curIndex.data(Qt::UserRole).value<SimulationNodeClass*>();

    //! -------------------------------------------------
    //! if the selected item contains an simulation root
    //! node, immediately return
    //! -------------------------------------------------
    if(node->isAnalysisRoot()) return treeModel->itemFromIndex(curIndex);

    //! ------------------------------------------------------
    //! else get the parent, and then the child: if the child
    //! is an analysis settings item return the parent
    //! ------------------------------------------------------
    QModelIndex parentIndex = curIndex.parent();
    QModelIndex firstChild = parentIndex.child(0,0);
    SimulationNodeClass *nodeChild = firstChild.data(Qt::UserRole).value<SimulationNodeClass*>();
    if(nodeChild->isAnalysisSettings())
    {
        return treeModel->itemFromIndex(parentIndex);
    }
    return Q_NULLPTR;
}
