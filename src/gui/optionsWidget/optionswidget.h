#ifndef OPTIONSWIDGET_H
#define OPTIONSWIDGET_H

#include <QWidget>
#include <QTreeView>
#include <QTableView>
#include <QStandardItem>
#include <QStandardItemModel>
#include <QItemSelectionModel>
#include <QDialog>

class optionsWidget : public QWidget
{
    Q_OBJECT

public:

    optionsWidget(QWidget *parent = 0);
    ~optionsWidget();

private:

    QTreeView *myTreeView;
    QTableView *myTableView;
    QStandardItemModel *myModel;
    QItemSelectionModel *mySelectionModel;
    QStandardItem *myRoot;

private:

    void createItems();

public slots:

    void showContent(QModelIndex newIndex, QModelIndex oldIndex);

};

#endif // OPTIONSWIDGET_H
