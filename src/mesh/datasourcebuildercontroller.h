#ifndef DATASOURCEBUILDERCONTROLLER_H
#define DATASOURCEBUILDERCONTROLLER_H

//! ---
//! Qt
//! ---
#include <QObject>
#include <QThread>
#include <QStandardItem>

//! ----------------
//! custom includes
//! ----------------
#include <simulationdatabase.h>

class dataSourceBuilderController: public QObject
{
    Q_OBJECT

public:

    dataSourceBuilderController();

    void setDataBase(simulationDataBase *sDB) { myDB = sDB; }
    void setItem(QStandardItem *anItem) { myItem = anItem; }

private:

    QThread workerThread;
    QStandardItem *myItem;
    simulationDataBase *myDB;
};

#endif // DATASOURCEBUILDERCONTROLLER_H
