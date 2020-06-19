#include "itemselector.h"

#include <QHBoxLayout>
#include <QGridLayout>

#include <iostream>
using namespace std;

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
itemSelector::itemSelector(QWidget *parent):QLineEditExtended(parent)
{
    this->createContent();
    this->createConnections();
    this->hidePushButtons();
}

//! ------------------------
//! function: constructor I
//! details:
//! ------------------------
itemSelector::itemSelector(QTreeView *aTreeView, QWidget *parent):QLineEditExtended(parent),treeView(aTreeView)
{
    this->createContent();
    this->createConnections();
    this->hidePushButtons();
}

//! ----------------------------
//! function: createConnections
//! details:
//! ----------------------------
void itemSelector::createConnections()
{
    connect(this,SIGNAL(activateSelector()),this,SLOT(showPushButtons()));
    connect(bap,SIGNAL(clicked()),this,SLOT(setAccepted()));
    connect(bcn,SIGNAL(clicked()),this,SLOT(setRejected()));
    connect(this,SIGNAL(activateSelector()),this,SLOT(setStrongFocus()));
}

void itemSelector::setStrongFocus()
{
    this->setFocusPolicy(Qt::StrongFocus);
}

//! ------------------------
//! function: createContent
//! details:
//! ------------------------
void itemSelector::createContent()
{
    this->updateText();

    //! -----------------------------------
    //! hide the line edit blinking cursor
    //! -----------------------------------
    this->setReadOnly(true);

    QHBoxLayout *h = new QHBoxLayout(this);
    h->setMargin(0);
    h->setContentsMargins(0,0,0,0);
    h->setSpacing(0);

    QGridLayout *g = new QGridLayout();
    g->setHorizontalSpacing(0);
    g->setMargin(0);
    g->setContentsMargins(0,0,0,0);

    h->addLayout(g);

    bap = new QPushButton("apply");
    bcn = new QPushButton("cancel");
    g->addWidget(bap,0,0);
    g->addWidget(bcn,0,1);
}

//! -------------------------
//! function: showPushButton
//! details:
//! -------------------------
void itemSelector::showPushButtons()
{
    //if(!(bap->isVisible() && bcn->isVisible()))
    {
        bap->show();
        bcn->show();
    }
}

//! -------------------------
//! function: hidePushButton
//! details:
//! -------------------------
void itemSelector::hidePushButtons()
{
    cout<<"itemSelector::hidePushButtons()->____function called____"<<endl;
    //if((bap->isVisible() && bcn->isVisible()))
    {
        cout<<"itemSelector::hidePushButtons()->____function called____"<<endl;
        bap->hide();
        bcn->hide();
    }
}

//! ----------------------
//! function: setAccepted
//! details:
//! ----------------------
void itemSelector::setAccepted()
{
    hidePushButtons();
    updateText();
    emit editingFinished();
}

//! ----------------------
//! function: setRejected
//! details:
//! ----------------------
void itemSelector::setRejected()
{
    hidePushButtons();
    updateText();
    emit editingFinished();
}

//! ---------------------
//! function: updateText
//! details:
//! ---------------------
void itemSelector::updateText()
{
    //! -----------------------------
    //! set text (number of objects)
    //! -----------------------------
    int NbObjects = objects.length();
    this->setText(QString("%1 objects").arg(NbObjects));

    // change the backgound color: 0 Objects => yellow
    // >0 Objects => white ... to do
}

//! ---------------------
//! function: setObjects
//! details:
//! ---------------------
void itemSelector::setObjects(const QList<void*> objectList)
{
    int Nb = objectList.length();
}
