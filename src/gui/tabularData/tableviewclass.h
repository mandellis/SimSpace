#ifndef TABLEVIEWCLASS_H
#define TABLEVIEWCLASS_H

#include <QWidget>
#include <QTableView>
#include <QStyledItemDelegate>

class tableViewClass: public QTableView
{
    Q_OBJECT

public:

    //! constructor
    tableViewClass(QWidget *parent=0);

public slots:

    //! set the model
    void setTheModel(QModelIndex index);
};

#endif // TABLEVIEWCLASS_H
