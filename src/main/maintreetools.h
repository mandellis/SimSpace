#ifndef MAINTREETOOLS_H
#define MAINTREETOOLS_H

//! ---
//! Qt
//! ---
#include <QList>
#include <QTreeView>
#include <QStandardItem>
#include <QModelIndex>

//! ----------------
//! custom includes
//! ----------------
#include "geometrytag.h"
#include "simulationnodeclass.h"

//! ----
//! C++
//! ----
#include <set>

class mainTreeTools
{
public:

    mainTreeTools(){;}

    static int calculateStartColumn(QTreeView *tree, int columnsBeforeBC);
    static int calculateStartColumn(QStandardItem *anItem, int columnsBeforeBC);

    static QList<int> getColumnsToRead(QTreeView *tree, int columnsBeforeBC);
    static QList<int> getColumnsToRead(QStandardItem *anItem, int columnsBeforeBC);

    static QStandardItem* getCurrentSimulationRoot(QTreeView *treeView);

    static void formNewPart(QTreeView *tree);

    static QList<QStandardItem*> getAllTreeItemOfType(QTreeView *tree,SimulationNodeClass::nodeType theNodeType);
    static void getTreeItemsRecursively(QStandardItemModel* model, QList<QStandardItem*> &items, QModelIndex parent = QModelIndex());

    static void addSolution(QStandardItem *analysisRootItem);
    static void addSolutionInformation(QStandardItem* solutionItem);
    static void resetSolutionInformation(SimulationNodeClass *nodeSolutionInformation);

    //! ----------------------------------------------------------
    //! given a simulation tree retrieve all the tags of boundary
    //! conditions, including contacts
    //! type: 2 solid, 3 shell, 4 face, 5 wire, 6 edge, 7 vertex
    //! ----------------------------------------------------------
    static void getAllBoundaryConditionsTags(QTreeView *tree, int type=0, std::vector<GeometryTag> &tags =std::vector<GeometryTag>());

    static SimulationNodeClass* getAnalysisSettingsNodeFromIndex(QModelIndex curIndex);
    static QStandardItem* getAnalysisSettingsItemFromCurrentItem(QTreeView *treeView);
    static SimulationNodeClass* getAnalysisSettingsNodeFromCurrentItem(QTreeView *treeView);
    static QStandardItem* getFirstTreeItemOfType(SimulationNodeClass::nodeType aType, QStandardItemModel* model);
    static bool getTreeItemsFromShapes(QTreeView *tree, const std::vector<TopoDS_Shape> &vecShapes, std::vector<QStandardItem*> &vecItems);
    static QExtendedStandardItem* getTreeItem(QStandardItemModel* model, SimulationNodeClass::nodeType theNodeType);
    static QExtendedStandardItem* ItemFromScope(QStandardItemModel* model, const TopoDS_Shape &aShape);
    static QList<QStandardItem*> getAllTreeItemOfType(QStandardItemModel* model, SimulationNodeClass::nodeType theNodeType);
    static const int mainTreeTools::getInsertionRow(QTreeView *tree);
    static QStandardItem* getSolutionItemFromCurrentItem(QTreeView *treeView);

};

#endif // MAINTREETOOLS_H
