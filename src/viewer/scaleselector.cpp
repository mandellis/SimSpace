//! ----------------
//! custom includes
//! ----------------
#include "scaleselector.h"
//#include "global.h"

//! ---
//! Qt
//! ---
#include <QDoubleValidator>
#include <QLineEdit>

//! ----
//! C++
//! ----
#include <iostream>
using namespace std;

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
ScaleSelector::ScaleSelector(QWidget *parent):QComboBox(parent)
{
    //cout<<"ScaleSelector::ScaleSelector()->____CONSTRUCTOR CALLED____"<<endl;
    this->addItem("Undeformed",0.0);            //! Nr. 0
    this->addItem("True scale",1.0);            //! Nr. 1
    this->addItem("x 2",2.0);                   //! Nr. 2
    this->addItem("x 5",5.0);                   //! Nr. 3
    this->addItem("x 10",10.0);                 //! Nr. 4
    this->addItem("x 20",20.0);                 //! Nr. 5
    this->addItem("x 50",50.0);                 //! Nr. 6
    this->addItem("x 100",100.0);               //! Nr. 7
    this->addItem("x 200",200.0);               //! Nr. 8
    this->addItem("x 500",500.0);               //! Nr. 9
    this->addItem("x 1000",1000.0);             //! Nr. 10

    connect(this,SIGNAL(currentIndexChanged(int)),this,SLOT(emitScaleChanged()));
}

//! ---------------------------
//! function: emitScaleChanged
//! details:
//! ---------------------------
void ScaleSelector::emitScaleChanged()
{
    double scale = this->currentData().toDouble();
    //Global::status().myResultPresentation.theScale = scale;
    emit scaleChanged(scale);
}

//! ---------------------
//! function: destructor
//! details:
//! ---------------------
ScaleSelector::~ScaleSelector()
{
    cout<<"ScaleSelector::~ScaleSelector()->____DESTRUCTOR CALLED____"<<endl;
}
