#ifndef MAINTREETOOLS_H
#define MAINTREETOOLS_H

#include <QList>
#include <QTreeView>
#include <QStandardItem>
#include <QModelIndex>
#include "simulationnodeclass.h"

class mainTreeTools
{
public:

    mainTreeTools(){;}

    static QList<int> getColumnsToRead(QTreeView *tree);
    static QList<int> getColumnsToRead1(QTreeView *tree);
    static int calculateStartColumn(QTreeView *tree);
    static void formNewPart(QTreeView *tree);
    static QList<QStandardItem*> getAllTreeItemOfType(QTreeView *tree,SimulationNodeClass::nodeType theNodeType);
    static void getTreeItemsRecursively(QStandardItemModel* model, QList<QStandardItem*> &items, QModelIndex parent = QModelIndex());
    static QStandardItem* getCurrentSimulationRoot(QTreeView *treeView);
};

#endif // MAINTREETOOLS_H
