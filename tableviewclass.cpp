//! -----------------
//! custrom includes
//! -----------------
#include "tableviewclass.h"
#include "simulationnodeclass.h"
#include "customtablemodel.h"
#include "tableviewclassitemdelegate.h"

//! ---
//! Qt
//! ---
#include <QHeaderView>

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
tableViewClass::tableViewClass(QWidget *parent):QTableView(parent)
{
    cout<<"tableViewClass::tableViewClass()->____constructor called____"<<endl;
    tableViewClassItemDelegate *delegate = new tableViewClassItemDelegate();
    this->setItemDelegate(delegate);
    this->setStyleSheet("QHeaderView::section { background-color: lightgray } "
                        "QTableCornerButton::section { background: lightgray; border: 2px outset lightgray;}");
}

//! ----------------------
//! function: setTheModel
//! details:
//! ----------------------
void tableViewClass::setTheModel(QModelIndex index)
{
    if(index.isValid())
    {
        SimulationNodeClass *node = index.data(Qt::UserRole).value<SimulationNodeClass*>();
        cout<<"tableViewClass::setTheModel()->____setTheModel() on node: "<<node->getName().toStdString()<<"____"<<endl;
        if(node!=NULL && node->getTabularDataModel()!=NULL)
        {
            CustomTableModel *theModel = node->getTabularDataModel();
            this->setModel(theModel);
        }
    }
    this->resizeColumnsToContents();
    this->resizeRowsToContents();
}
