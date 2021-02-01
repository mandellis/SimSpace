#ifndef TABLEWIDGET_H
#define TABLEWIDGET_H

//! ----------------
//! custom includes
//! ----------------
#include "usermessagesmodel.h"

//! ---
//! Qt
//! ---
#include <QtCharts/QChart>
#include <QtCharts/QChartView>
#include <QtCharts/QLineSeries>
#include <QtCharts/QVXYModelMapper>
#include <QtCharts/QValueAxis>
#include <QHeaderView>
#include <QTabWidget>
#include <QPoint>

//! ----
//! C++
//! ----
#include <iostream>

using namespace QtCharts;

class QTableView;

//class TableWidget: public QWidget
class TableWidget: public QTabWidget
{
    Q_OBJECT

public:

    TableWidget(QWidget *parent = 0);
    TableWidget::~TableWidget()
    {
        std::cout<<"TableWidget::~TableWidget()->____DESTRUCTOR CALLED____"<<std::endl;
    }

private:

    QTableView *tableView;
    QTableView *messageTable;

public:

    QTableView *getTableView() { return tableView; }

public slots:

    void setTheModel(QModelIndex index);
    void showColumns(QList<int> columnsToShow);
    void hideAllColumnsAndHeaders();
    void hideFirstRow();
    void showFirstRow();
    void synchBoltPretension();
    void setMessagesModel(userMessagesModel *messageModel);


private slots:

    void handleModelDataChanged(QModelIndex topLeftIndex, QModelIndex bottomRightIndex, QVector<int> roles);
    void changeWindowTitle(int tabNumber);
    void showContextMenu(QPoint aPoint);

signals:

    void requestUpdateDetailViewer(QModelIndex tl,QModelIndex br,QVector<int> roles);
    void tabularDataChanged();
};

#endif // TABLEWIDGET_H
