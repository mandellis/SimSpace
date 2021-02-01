#include "lineedit.h"
#include <QToolButton>
#include <QStyle>
#include <QAction>
#include <QMenu>
#include <QDoubleValidator>

#include <iostream>
using namespace std;

//! ------------------------
//! function: createContent
//! details:
//! ------------------------
void LineEdit::createContent(bool hasFreeOption)
{
    //!cout<<"LineEdit::createContent()->____function called____"<<endl;
    QDoubleValidator *doubleValidator = new QDoubleValidator();
    this->setValidator(doubleValidator);

    Button = new QToolButton(this);
    Button->setIcon(QIcon(":/icons/icon_right arrow.png"));
    Button->setCursor(Qt::ArrowCursor);
    Button->setStyleSheet("QToolButton { border: none; padding: 0px; }");
    Button->hide();

    //! create the menu actions
    actionConstant = new QAction("Constant",this);
    actionConstant->setCheckable(true);
    actionConstant->setData(1);

    actionTabular = new QAction("Tabular",this);
    actionTabular->setCheckable(true);
    actionTabular->setData(2);

    QActionGroup *group = new QActionGroup(this);
    group->addAction(actionConstant);
    group->addAction(actionTabular);

    myContextMenu = new QMenu(this);

    myContextMenu->addAction(actionConstant);
    myContextMenu->addAction(actionTabular);

    if(hasFreeOption==true)
    {
        actionFree = new QAction("Free",this);
        actionFree->setCheckable(true);
        actionFree->setData(3);
        group->addAction(actionFree);
        myContextMenu->addAction(actionFree);
    }

    int frameWidth = style()->pixelMetric(QStyle::PM_DefaultFrameWidth);
    setStyleSheet(QString("QLineEdit { padding-right: %1px; } ").arg(Button->sizeHint().width() + frameWidth + 1));
    QSize msz = minimumSizeHint();
    this->setMinimumSize(qMax(msz.width(), Button->sizeHint().height() + frameWidth * 2 + 2),
                         qMax(msz.height(), Button->sizeHint().height() + frameWidth * 2 + 2));
}

//! ---------------------------------
//! function: setData
//! details:  called by the delegate
//! ---------------------------------
void LineEdit::setData(Property::loadDefinition theLoadDefinition)
{
    //cout<<"LineEdit::setData()->____function called____"<<endl;
    myLoadDefinition = theLoadDefinition;
    switch(myLoadDefinition)
    {
    case Property::loadDefinition_constant:
        actionConstant->setChecked(true);
        this->setReadOnly(false);
        break;
    case Property::loadDefinition_tabularData:
        actionTabular->setChecked(true);
        this->setText("Tabular data");
        this->setReadOnly(true);
        break;
    case Property::loadDefinition_free:
        actionFree->setChecked(true);
        this->setText("Free");
        this->setReadOnly(true);
        break;
    }
}

//! --------------------------
//! function: showContextMenu
//! details:
//! --------------------------
void LineEdit::showContextMenu()
{
    //cout<<"LineEdit::showContextMenu()->____Showing context menu____"<<endl;
    QAction* selectedItem = myContextMenu->exec(this->mapToGlobal(Button->geometry().bottomRight()));
    if(selectedItem)
    {
        //cerr<<"____selectedItem: "<<selectedItem->data().toInt()<<"____"<<endl;
        switch(selectedItem->data().toInt())
        {
        case 1:
            myLoadDefinition = Property::loadDefinition_constant;
            break;
        case 2:
            myLoadDefinition = Property::loadDefinition_tabularData;
            break;
        case 3:
            myLoadDefinition = Property::loadDefinition_free;
            break;
        }
        //cout<<"LineEdit::showContextMenu()->____selection done____"<<endl;
        //cout<<"LineEdit::showContextMenu()->____load definition: "<<myLoadDefinition<<"____"<<endl;
        //cout<<"LineEdit::showContextMenu()->____emitting signal editing finished____"<<endl;
        //emit editingFinished();
        this->emitEditingFinished();
    }
}

//! ------------------------------
//! function: emitEditingFinished
//! details:
//! ------------------------------
void LineEdit::emitEditingFinished()
{
    //cout<<"LineEdit::emitEditingFinished()->____emit editingFinished()____"<<endl;
    myTextData="";
    emit editingFinished(myTextData);
}

//! --------------------------------------------------------
//! function: emitEditingFinishedAfterTextChange
//! details:  a text change automatically sets the
//!           the loadDefinition to loadDefinition_constant
//! --------------------------------------------------------
void LineEdit::emitEditingFinishedAfterTextChange()
{
    //cout<<"LineEdit::emitEditingFinishedAfterTextChange()->____emit editingFinished()____"<<endl;
    //cout<<"LineEdit::emitEditingFinishedAfterTextChange()->____final text data: "<<myTextData.toStdString()<<"____"<<endl;
    myLoadDefinition = Property::loadDefinition_constant;
    emit editingFinished(myTextData);
}

//! -----------------------
//! function: changeMyText
//! details:
//! -----------------------
void LineEdit::changeMyText()
{
    myTextData = this->text();
    cout<<"text data: "<<myTextData.toStdString()<<endl;
}

//! ------------------
//! function: getText
//! details:
//! ------------------
QString LineEdit::getText() const
{
    return myTextData;
}

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
LineEdit::LineEdit(bool hasFreeOption, QWidget *parent): QLineEdit(parent)
{
    //!cout<<"LineEdit::LineEdit()->____constructor called____"<<endl;
    this->createContent(hasFreeOption);
    this->setContextMenuPolicy(Qt::CustomContextMenu);

    //! -------------------------------------------------
    //! the default status of the button at the creation
    //! -------------------------------------------------
    Button->show();
    if(hasFreeOption)actionFree->setChecked(true);
    else actionConstant->setChecked(true);

    //! ------------
    //! connections
    //! ------------
    connect(Button,SIGNAL(released()),this,SLOT(showContextMenu()));
    connect(this,SIGNAL(returnPressed()),this,SLOT(emitEditingFinishedAfterTextChange()));
    connect(this,SIGNAL(textEdited(QString)),this,SLOT(changeMyText()));
    connect(this,SIGNAL(selectionChanged()),this,SLOT(unselect()));
}

//! ----------------------
//! function: resizeEvent
//! details:
//! ----------------------
void LineEdit::resizeEvent(QResizeEvent *)
{
    QSize sz = Button->sizeHint();
    int frameWidth = style()->pixelMetric(QStyle::PM_DefaultFrameWidth);
    Button->move(rect().right() - frameWidth - sz.width(),(rect().bottom() + 1 - sz.height())/2);
}

//! ----------------------------
//! function: getLoadDefinition
//! details:
//! ----------------------------
Property::loadDefinition LineEdit::getLoadDefinition() const
{
    return myLoadDefinition;
}
