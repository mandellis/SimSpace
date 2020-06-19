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

class mainTreeTools
{
public:

    mainTreeTools(){;}

    //static QList<int> getColumnsToRead1(QTreeView *tree);
    static int calculateStartColumn(QTreeView *tree);
    static int calculateStartColumn(QExtendedStandardItem *anItem);

    static QList<int> getColumnsToRead(QTreeView *tree);
    static QList<int> getColumnsToRead(QExtendedStandardItem *anItem);

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
    static void getAllBoundaryConditionsTags(QTreeView *tree, int type=0, QVector<GeometryTag> &tags =QVector<GeometryTag>());
};

#endif // MAINTREETOOLS_H
