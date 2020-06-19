//! ----------------
//! custom includes
//! ----------------
#include "tablewidget.h"
#include "customtablemodel.h"
#include "tableviewclass.h"
#include "simulationnodeclass.h"
#include "tableviewclassitemdelegate.h"
#include "customtablemodel.h"

//! ---
//! Qt
//! ---
#include <QHBoxLayout>
#include <QTableView>
//#include <QValueAxis>
//#include <QXYLegendMarker>
#include <QGraphicsLayout>
#include <QScrollBar>
#include <QMenu>
#include <QAction>
#include <QIcon>

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
TableWidget::TableWidget(QWidget *parent): QTabWidget(parent)
{
    std::cout<<"TableWidget::TableWidget()->____CONSTRUCTOR CALLED____"<<std::endl;

    this->setObjectName("messagesAndLoadsWidget");

    //! -------------------------------
    //! the table containing the loads
    //! -------------------------------
    //QHBoxLayout *mainLayout = new QHBoxLayout();
    tableView = new QTableView(this);
    tableView->horizontalHeader()->setDefaultAlignment(Qt::AlignHCenter);
    tableView->verticalHeader()->setDefaultAlignment(Qt::AlignVCenter);
    tableView->setContentsMargins(0,0,0,0);

    tableViewClassItemDelegate *delegate = new tableViewClassItemDelegate(this);
    tableView->setItemDelegate(delegate);

    tableView->setStyleSheet("QHeaderView::section { background-color: lightgray } "
                        "QTableCornerButton::section { background: lightgray; border: 2px outset lightgray;}");

    connect(delegate,SIGNAL(boltStatusDefineByChanged()),this,SLOT(synchBoltPretension()));

    //! ------------
    //! add the tab
    //! ------------
    this->addTab(tableView,"Tabular data");
    this->setTabPosition(QTabWidget::South);
    this->setTabShape(QTabWidget::Triangular);

    //! ---------------------------
    //! the table for the messages
    //! ---------------------------
    messageTable = new QTableView(this);
    messageTable->setContentsMargins(0,0,0,0);
    messageTable->setAlternatingRowColors(true);
    messageTable->verticalScrollBar()->setVisible(true);

    //! ---------------------------------------------------
    //! with this line the table fits the container widget
    //! ---------------------------------------------------
    messageTable->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
    messageTable->resizeRowsToContents();

    //! -------------------
    //! enable custom menu
    //! -------------------
    messageTable->setContextMenuPolicy(Qt::CustomContextMenu);
    connect(messageTable,SIGNAL(customContextMenuRequested(QPoint)),this,SLOT(showContextMenu(QPoint)));

    //! ------------
    //! add the tab
    //! ------------
    this->addTab(messageTable,"Messages");

    //! ----------------------------------------------------
    //! update the window title according to the active tab
    //! ----------------------------------------------------
    connect(this,SIGNAL(tabBarClicked(int)),this,SLOT(changeWindowTitle(int)));
}

//! ---------------------------
//! function: setMessagesModel
//! details:
//! ---------------------------
void TableWidget::setMessagesModel(userMessagesModel *messageModel)
{
    messageTable->setModel(messageModel);
}

//! ----------------------------
//! function: changeWindowTitle
//! details:
//! ----------------------------
void TableWidget::changeWindowTitle(int tabNumber)
{
    QString title;
    switch(tabNumber)
    {
    case 0: title = QString("Tabular data"); break;
    case 1: title = QString("Messages"); break;
    }
    this->parentWidget()->setWindowTitle(title);
}

//! ------------------------------
//! function: synchBoltPretension
//! details:
//! ------------------------------
void TableWidget::synchBoltPretension()
{
    //cout<<"TableWidget::synchBoltPretension()->____synch bolt pretension item____"<<endl;
    QModelIndex index = tableView->currentIndex();
    cout<<"TableWidget::synchBoltPretension()->____synch bolt pretension item: row: "<<index.row()<<" column: "<<index.column()<<"____"<<endl;
    CustomTableModel *tabularDataModel = static_cast<CustomTableModel*>(tableView->model());
    tabularDataModel->autoUpdateTable(index,index);
}

//! -----------------------------------
//! function: hideAllColumnsAndHeaders
//! details:
//! -----------------------------------
void TableWidget::hideAllColumnsAndHeaders()
{
    tableView->verticalHeader()->setHidden(true);
    int Ncol = tableView->model()->columnCount();
    for(int k=0; k<Ncol; k++)
    {
        tableView->setColumnHidden(k,true);
    }
}

//! ---------------------------------
//! function: handleModelDataChanged
//! details:
//! ---------------------------------
void TableWidget::handleModelDataChanged(QModelIndex topLeftIndex, QModelIndex bottomRightIndex, QVector<int> roles)
{
    cout<<"TableWidget::handleModelDataChanged()->____function called. Emitting requestUpdateDetailViewer()____"<<endl;
    emit requestUpdateDetailViewer(topLeftIndex,bottomRightIndex,roles);
}

//! ----------------------
//! function: setTheModel
//! details:
//! ----------------------
void TableWidget::setTheModel(QModelIndex index)
{
    if(index.isValid())
    {
        SimulationNodeClass *node = index.data(Qt::UserRole).value<SimulationNodeClass*>();
        if(node!=NULL && node->getTabularDataModel()!=NULL)
        {
            CustomTableModel *theModel = node->getTabularDataModel();
            tableView->setModel(theModel);

            disconnect(theModel,SIGNAL(dataChanged(QModelIndex,QModelIndex,QVector<int>)),this,SLOT(handleModelDataChanged(QModelIndex,QModelIndex,QVector<int>)));
            connect(theModel,SIGNAL(dataChanged(QModelIndex,QModelIndex,QVector<int>)),this,SLOT(handleModelDataChanged(QModelIndex,QModelIndex,QVector<int>)));

            tableView->resizeColumnsToContents();
            tableView->resizeRowsToContents();

            //! -------------------------------------------------------
            //! if the model is empty hide the headers, else show them
            //! -------------------------------------------------------
            if(theModel->rowCount()==0)tableView->horizontalHeader()->hide();
            else tableView->horizontalHeader()->show();
        }
    }
}

//! ----------------------------------------------------
//! function: showColumns
//! details:  show specific columns of the tabular data
//! ----------------------------------------------------
void TableWidget::showColumns(QList<int> columnsToShow)
{
    //cout<<"TableWidget::showColumns()->____function called____"<<endl;
    int Ncolumns = tableView->model()->columnCount();
    for(int i=0; i<Ncolumns; i++)
    {
        if(columnsToShow.contains(i))
        {
            //cout<<"TableWidget::showingColumns()->____showing columns n. "<<i<<"____"<<endl;
            tableView->setColumnHidden(i,false);
        }
        else tableView->setColumnHidden(i,true);
    }
}

//! ----------------------------------------------------------
//! function: hideFirstRow
//! details:  used when clicking an "Analysis settings" items
//! ----------------------------------------------------------
void TableWidget::hideFirstRow()
{
    tableView->setRowHidden(0,true);
}

//! -----------------------------------------------------
//! function: showFirstRow
//! details:  used when clicking a simulation setup item
//! -----------------------------------------------------
void TableWidget::showFirstRow()
{
    tableView->setRowHidden(0,false);
}

//! --------------------------
//! function: ShowContextMenu
//! details:
//! --------------------------
void TableWidget::showContextMenu(QPoint aPoint)
{
    cout<<"TableWidget::ShowContextMenu()->____context menu activated____"<<endl;
    QPoint pos = this->mapToGlobal(aPoint);
    QMenu *ctxMenu= new QMenu(this);

    //! ----------------------------------------------------
    //! add the option "Remove" only when items are present
    //! ----------------------------------------------------
    if(messageTable->model()->rowCount(QModelIndex())>0)
    {
        QAction *actionRemove = ctxMenu->addAction("Remove");
        actionRemove->setIcon(QIcon(":/icons/icon_remove.png"));
        actionRemove->setData(0);

        ctxMenu->addSeparator();

        QAction *actionClear = ctxMenu->addAction("Clear");
        actionClear->setIcon(QIcon(":/icons/icon_clear.png"));
        actionClear->setData(1);
    }

    QAction* selectedItem = ctxMenu->exec(pos);
    if(selectedItem)
    {
        switch(selectedItem->data().toInt())
        {
        case 0:
        {
            QModelIndex clickedIndex = messageTable->indexAt(aPoint);
            int row = clickedIndex.row();
            ((userMessagesModel*)(messageTable->model()))->removeMessage(row);
        }
            break;

        case 1:
        {
            ((userMessagesModel*)(messageTable->model()))->clear();
        }
            break;
        }
    }
}
