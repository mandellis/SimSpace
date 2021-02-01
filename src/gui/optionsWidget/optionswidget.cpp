//! custom includes
#include "optionswidget.h"
#include "colorselector.h"
#include "celldata.h"
#include "tabledelegate.h"
#include "qtablestandarditem.h"

//! Qt
#include <QHBoxLayout>
#include <QSplitter>
#include <QVariant>
#include <QStandardItem>

//! C++
#include <iostream>

using namespace std;

//! -----------------------------------------------------
//! function: constructor
//! details:
//! -----------------------------------------------------
optionsWidget::optionsWidget(QWidget *parent): QWidget(parent)
{
    //! create the content

    //! [1] layout
    QHBoxLayout *hLayout = new QHBoxLayout(this);
    hLayout->setMargin(0);
    hLayout->setContentsMargins(0,0,0,0);
    hLayout->setSpacing(0);

    //! [2] splitter
    QSplitter *split = new QSplitter(this);
    split->setOrientation(Qt::Horizontal);
    hLayout->addWidget(split);

    //! [2] the tree view
    myTreeView = new QTreeView(this);
    split->addWidget(myTreeView);

    //! [3] create the table view
    myTableView = new QTableView(this);
    tableDelegate *tabDelegate = new tableDelegate(this);
    myTableView->setItemDelegate(tabDelegate);
    split->addWidget(myTableView);

    //! [4] create the model
    myModel = new QStandardItemModel(this);
    //myModel->setColumnCount(2);

    //! [5] set the model
    myTreeView->setModel(myModel);

    //! [6] create the selection model
    mySelectionModel = myTreeView->selectionModel();
    connect(mySelectionModel,SIGNAL(currentChanged(QModelIndex,QModelIndex)),this,SLOT(showContent(QModelIndex,QModelIndex)));

    //! [7] get the root
    myRoot = myModel->invisibleRootItem();

    //! [8] create the items
    this->createItems();
}

//! ----------------------
//! function: createItems
//! details:
//! ----------------------
#include "viewoptions.h"
void optionsWidget::createItems()
{
    QVariant data;

    //! ----------------------
    //! [0] project management
    //! ----------------------
    QStandardItem* item_projectManagement = new QStandardItem();
    item_projectManagement->setText("Project management");

    //! -----------------------
    //! [1] graphic appearence
    //! -----------------------
    QStandardItem* item_appearence = new QStandardItem();
    item_appearence->setText("Appearence");

    //! [1.1] cell type of gradient
    //! 0 - no gradient 1 - vertical  - 2 horizontal
    cell TypeOfGradient;
    viewOptions tog = viewOptions::typeOfGradient_uniform;
    TypeOfGradient.name = QString("Type of gradient");
    TypeOfGradient.data.setValue(tog);

    //! [1.2] cell first color: initialize with white
    cell firstColor;
    firstColor.name = QString("First color");
    QColor aColor(Qt::white);
    firstColor.data.setValue(aColor);

    //! [1.3] cell second color: initialize with white
    cell secondColor;
    secondColor.name = QString("Second color");
    secondColor.data.setValue(aColor);

    //! insert the cells into the item
    QList<cell> itemContent;
    itemContent<<TypeOfGradient<<firstColor<<secondColor;
    data.setValue(itemContent);
    item_appearence->setData(data,Qt::UserRole);
    //! --------------------------------------------------------------------------
    //! end of graphic appearence
    //! --------------------------------------------------------------------------

    //! --------------------------------------------------------------------------
    //! [2] graphic interaction
    //! --------------------------------------------------------------------------
    QStandardItem* item_graphicInteraction = new QStandardItem();
    item_graphicInteraction->setText("Graphic interaction");
    //! --------------------------------------------------------------------------
    //! end of graphic interaction
    //! --------------------------------------------------------------------------

    myRoot->appendRow(item_projectManagement);
    myRoot->appendRow(item_appearence);
    myRoot->appendRow(item_graphicInteraction);
}

//! -----------------------------------------------------
//! function: showContent
//! details:  create a local model and show it in a table
//! -----------------------------------------------------
void optionsWidget::showContent(QModelIndex newIndex, QModelIndex oldIndex)
{
    Q_UNUSED(oldIndex)

    cout<<"optionsWidget::showContent()->____function called____"<<endl;

    //! title
    QString windowTitle = newIndex.data(Qt::DisplayRole).toString();
    this->setWindowTitle(windowTitle);

    //! build the local model with 2 rows
    QStandardItemModel *localModel = new QStandardItemModel(this);
    localModel->setColumnCount(2);

    cout<<"optionsWidget::showContent()->____item clicked: "<<newIndex.data(Qt::DisplayRole).value<QString>().toStdString()<<"____"<<endl;

    //! the list of cell
    QList<cell> listOfCell = newIndex.data(Qt::UserRole).value<QList<cell>>();
    cout<<"optionsWidget::showContent()->____number of cells: "<<listOfCell.length()<<"____"<<endl;

    int row = 0;
    for(QList<cell>::iterator it = listOfCell.begin(); it!= listOfCell.end(); ++it, ++row)
    {
        const cell &aCell = *it;
        QString name = aCell.name;
        cout<<"optionsWidget::showContent()->____creating cell: "<<name.toStdString()<<"____"<<endl;

        QVariant data;
        QVariant content = aCell.data;

        //! item name
        QTableStandardItem *itemName = new QTableStandardItem;
        data.setValue(name);
        itemName->setText(name);
        itemName->setData(data,Qt::UserRole);

        //! item cell content
        QTableStandardItem *itemData = new QTableStandardItem;
        data.setValue(content);
        itemData->setData(data,Qt::UserRole);

        localModel->setItem(row,0,itemName);
        localModel->setItem(row,1,itemData);
    }
    //! set and show the model
    myTableView->setModel(localModel);
}

//! -----------------------------------------------------
//! function: destructor
//! details:
//! -----------------------------------------------------
optionsWidget::~optionsWidget()
{
    cout<<"optionsWidget::~optionsWidget()->____DESTRUCTOR CALLED____"<<endl;
}
