//! custom includes
#include "boundaryvaluemanager.h"

//! Qt
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QPushButton>
#include <QMenu>
#include <QMessageBox>
#include <QAction>

//! C++
#include <iostream>

using namespace std;

void BoundaryValueManager::buildContextMenu()
{
    //! create menu actions
    actionConstant = new QAction("Constant",this);
    actionTabular = new QAction("Tabular",this);
    actionFunction = new QAction("Function",this);
    actionFree = new QAction("Free",this);

    actionConstant->setCheckable(true);
    actionTabular->setCheckable(true);
    actionFunction->setCheckable(true);
    actionFree->setCheckable(true);

    actionConstant->setData(1);
    actionTabular->setData(2);
    actionFunction->setData(3);
    actionFree->setData(4);

    QActionGroup *group = new QActionGroup(this);
    group->addAction(actionConstant);
    group->addAction(actionTabular);
    group->addAction(actionFunction);
    group->addAction(actionFree);

    myContextMenu = new QMenu(this);

    myContextMenu->addAction(actionConstant);
    myContextMenu->addAction(actionTabular);
    myContextMenu->addAction(actionFunction);
    myContextMenu->addAction(actionFree);
}

//! -----------------------------------------------------------------//
//! function: constructor I                                          //
//! details:                                                         //
//! -----------------------------------------------------------------//
BoundaryValueManager::BoundaryValueManager(QWidget *parent): QWidget(parent), myTypeOfValue(Property::typeOfValue_constant)
{
    //! build the context menu content
    buildContextMenu();

    QHBoxLayout *lay = new QHBoxLayout(this);
    lineEdit = new ShapeSelectorBox();

    QString name;
    switch(myTypeOfValue)
    {
    case Property::typeOfValue_constant:
        lineEdit->setReadOnly(false);
        name="Constant";
        break;
    case Property::typeOfValue_tabular:
        lineEdit->setReadOnly(true);
        name="Tabular";
        break;
    case Property::typeOfValue_function:
        lineEdit->setReadOnly(false);
        name="Function";
        break;
    case Property::typeOfValue_free:
        lineEdit->setReadOnly(true);
        name="Free";
        break;
    }
    lineEdit->setText(name);

    but = new QPushButton();
    but->setIcon(QIcon(":/icons/icon_right arrow.png"));
    this->setContextMenuPolicy(Qt::CustomContextMenu);

    lay->addWidget(lineEdit);
    lay->addWidget(but);
    lay->setMargin(0);
    //lay->setContentsMargins(0,0,0,0);
    //but->hide();

    //connect(lineEdit,SIGNAL(activateSelector()),but,SLOT(show()));
    connect(but,SIGNAL(released()),this,SLOT(showContextMenu()));
}

//! -----------------------------------------------------------------//
//! function: show context menu                                      //
//! details:                                                         //
//! -----------------------------------------------------------------//
void BoundaryValueManager::setValue(Property::typeOfValue theTypeOfValue)
{
    myTypeOfValue = theTypeOfValue;

    QString name;
    switch(myTypeOfValue)
    {
    case Property::typeOfValue_constant:
        actionConstant->setChecked(true);
        lineEdit->setReadOnly(false);
        name="Constant (ramped)";
        emit typeOfValue_constant_selected();
        break;
    case Property::typeOfValue_tabular:
        actionTabular->setChecked(true);
        lineEdit->setReadOnly(true);
        name="Tabular";
        emit typeOfValue_tabular_selected();
        break;
    case Property::typeOfValue_function:
        actionFunction->setChecked(true);
        lineEdit->setReadOnly(false);
        cout<<"BoundaryValueManager::setValue()->____function behavior selected____"<<endl;
        emit typeOfValue_function_selected();
        name="Function";
        break;
    case Property::typeOfValue_free:
        actionFree->setChecked(true);
        lineEdit->setReadOnly(true);
        emit typeOfValue_free_selected();
        name="Free";
        break;
    }
    lineEdit->setText(name);
}

//! -----------------------------------------------------------------//
//! function: show context menu                                      //
//! details:                                                         //
//! -----------------------------------------------------------------//
void BoundaryValueManager::showContextMenu()
{
    cout<<"showing context menu"<<endl;
    //QPoint globalPos = this->mapToGlobal(pos);

    QAction* selectedItem = myContextMenu->exec(this->mapToGlobal(but->geometry().bottomRight()));

    // check if the selection is valid
    if(selectedItem)
    {
        switch(selectedItem->data().toInt())
        {
        case 1:
            //! constant
            lineEdit->setText("Constant (ramped)");
            lineEdit->setReadOnly(false);
            myTypeOfValue=Property::typeOfValue_constant;
            emit typeOfValue_constant_selected();
            break;
        case 2:
            //! tabular
            lineEdit->setText("Tabular");
            lineEdit->setReadOnly(false);
            myTypeOfValue=Property::typeOfValue_tabular;
            emit typeOfValue_tabular_selected();
            break;
        case 3:
            //! function
            lineEdit->setText("Function");
            lineEdit->setReadOnly(false);
            myTypeOfValue=Property::typeOfValue_function;
            emit typeOfValue_function_selected();
            break;
        case 4:
            //! free
            lineEdit->setText("Free");
            lineEdit->setReadOnly(true);
            myTypeOfValue=Property::typeOfValue_free;
            emit typeOfValue_free_selected();
            break;
        }
        QPalette palette;
        palette.setColor(QPalette::Base,Qt::white);
        lineEdit->setPalette(palette);

        emit editingFinished();
    }
}


