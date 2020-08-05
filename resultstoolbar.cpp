//! ----------------
//! custom includes
//! ----------------
#include "resultstoolbar.h"
#include <scaleselector.h>

//! ---
//! Qt
//! ---
#include <QMenu>
#include <QAction>

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
ResultsToolBar::ResultsToolBar(const QString &title, QWidget *parent):QToolBar(title,parent)
{
    //! ----------------
    //! set object name
    //! ----------------
    this->setObjectName("resultsToolBar");

    //! ------------------------------
    //! add the scale selector widget
    //! ------------------------------
    ScaleSelector *scaleSelector = new ScaleSelector();

    //! ----------------------
    //! set the initial scale
    //! ----------------------
    this->addWidget(scaleSelector);
    connect(scaleSelector,SIGNAL(scaleChanged(double)),this,SLOT(emitRequestUpdatePostObjectScale(double)));
    myResultPresentation.isDeformedView = false;
    myResultPresentation.theScale = 1.0;

    //! -----------------------
    //! create a button widget
    //! -----------------------
    combinedViewSelectorButton = new QPushButtonExtended(this);
    combinedViewSelectorButton->setIcon(QIcon(":/icons/icon_no wireframe.png"));
    myResultPresentation.theCombinedView = resultPresentation::combinedView_resultOnly;

    //! -----------------------------
    //! create a menu for the button
    //! -----------------------------
    combinedViewSelectorMenu = new QMenu(this);

    //! ----------------------
    //! action "No wireframe"
    //! ----------------------
    actionNoWireframe = combinedViewSelectorMenu->addAction(QIcon(":/icons/icon_no wireframe.png"),"No wireframe");
    actionNoWireframe->setToolTip("No wireframe");
    connect(actionNoWireframe,SIGNAL(triggered(bool)),this,SLOT(emitRequestNoWireframe()));
    combinedViewSelectorMenu->addAction(actionNoWireframe);

    //! -----------------------------------
    //! action "Show undeformed wireframe"
    //! -----------------------------------
    actionShowUndeformedWireframe = combinedViewSelectorMenu->addAction(QIcon(":/icons/icon_wireframe.png"),"Show undeformed wireframe");
    actionShowUndeformedWireframe->setToolTip("Show undeformed wireframe");
    connect(actionShowUndeformedWireframe,SIGNAL(triggered(bool)),this,SLOT(emitRequestShowUndeformedWireframe()));
    combinedViewSelectorMenu->addAction(actionShowUndeformedWireframe);

    //! -------------------------------
    //! action "Show undeformed model"
    //! -------------------------------
    actionShowUndeformedModel = combinedViewSelectorMenu->addAction(QIcon(":/icons/icon_undeformed model.png"),"Show undeformed model");
    actionShowUndeformedModel->setToolTip("Show undeformed model");
    connect(actionShowUndeformedModel,SIGNAL(triggered(bool)),this,SLOT(emitRequestShowUndeformedModel()));
    combinedViewSelectorMenu->addAction(actionShowUndeformedModel);

    //! -----------------------
    //! action "Show elements"
    //! -----------------------
    actionShowElements = combinedViewSelectorMenu->addAction(QIcon(":/icons/icon_show elements.png"),"Show elements");
    actionShowElements->setToolTip("Show elements");
    connect(actionShowElements,SIGNAL(triggered(bool)),SLOT(emitRequestShowMeshElements()));
    combinedViewSelectorMenu->addAction(actionShowElements);

    //! --------------------
    //! set the button menu
    //! --------------------
    combinedViewSelectorButton->setMenu(combinedViewSelectorMenu);

    //! -------------------------------------
    //! add the widget button to the toolbar
    //! -------------------------------------
    this->addWidget(combinedViewSelectorButton);
}

//! ---------------------------------
//! function: emitRequestNoWireframe
//! details:
//! ---------------------------------
void ResultsToolBar::emitRequestNoWireframe()
{
    cout<<"ResultsToolBar()->____emit request \"no wireframe\"____"<<endl;
    combinedViewSelectorButton->setIcon(QIcon(":/icons/icon_no wireframe.png"));
    myResultPresentation.theCombinedView = resultPresentation::combinedView_resultOnly;
    emit requestSetResultPresentation(myResultPresentation);
}

//! -----------------------------------------
//! function: emitRequestShowUndeformedModel
//! details:
//! -----------------------------------------
void ResultsToolBar::emitRequestShowUndeformedModel()
{
    cout<<"ResultsToolBar()->____emit request show undeformed model____"<<endl;
    combinedViewSelectorButton->setIcon(QIcon(":/icons/icon_undeformed model.png"));
    myResultPresentation.theCombinedView = resultPresentation::combinedView_undeformedModel;
    emit requestSetResultPresentation(myResultPresentation);
}

//! ---------------------------------------------
//! function: emitRequestShowUndeformedWireframe
//! details:
//! ---------------------------------------------
void ResultsToolBar::emitRequestShowUndeformedWireframe()
{
    cout<<"ResultsToolBar()->____emit request show undeformed wireframe____"<<endl;
    combinedViewSelectorButton->setIcon(QIcon(":/icons/icon_wireframe.png"));
    myResultPresentation.theCombinedView = resultPresentation::combinedView_undeformedWireFrame;
    emit requestSetResultPresentation(myResultPresentation);
}

//! --------------------------------------
//! function: emitRequestShowMeshElements
//! details:
//! --------------------------------------
void ResultsToolBar::emitRequestShowMeshElements()
{
    cout<<"ResultsToolBar()->____emit request show mesh elements____"<<endl;
    combinedViewSelectorButton->setIcon(QIcon(":/icons/icon_show elements.png"));
    myResultPresentation.theCombinedView = resultPresentation::combinedView_meshVisible;
    emit requestSetResultPresentation(myResultPresentation);
}

//! ---------------------
//! function: updateIcon
//! details:
//! ---------------------
void ResultsToolBar::updateIcon(QAction *action)
{
    QIcon icon = action->icon();
    QWidget *widget = this->widgetForAction(action);
    QPushButtonExtended *button = static_cast<QPushButtonExtended*>(widget);
    button->setIcon(icon);
}

//! -------------------------------------------
//! function: emitRequestUpdatePostObjectScale
//! details:
//! -------------------------------------------
void ResultsToolBar::emitRequestUpdatePostObjectScale(double scale)
{
    cout<<"ResultsToolBar::emitRequestUpdatePostObjectScale()->____function called. Scale: "<<scale<<"____"<<endl;
    myResultPresentation.theScale = scale;
    emit requestSetResultPresentation(myResultPresentation);
    //emit requestUpdatePostObjectScale(scale);
}

//! ---------------------
//! function: saveStatus
//! details:
//! ---------------------
bool ResultsToolBar::saveStatus(const std::string &fileName)
{
    FILE *f = fopen(fileName.c_str(),"w");
    if(f==NULL) return false;
    fprintf(f,"%d\t%d\t%lf\n",
            myResultPresentation.theCombinedView,
            myResultPresentation.isDeformedView,
            myResultPresentation.theScale);
    return true;
}

//! --------------------
//! function: setStatus
//! details:
//! --------------------
void ResultsToolBar::setStatus(const resultPresentation &aResultPresentation)
{
    //actionNoWireframe->blockSignals(true);
    //actionShowElements->blockSignals(true);
    //actionShowUndeformedModel->blockSignals(true);
    //actionShowUndeformedWireframe->blockSignals(true);
    switch(aResultPresentation.theCombinedView)
    {
    case resultPresentation::combinedView_resultOnly:
        combinedViewSelectorMenu->setActiveAction(actionNoWireframe);
        break;
    case resultPresentation::combinedView_meshVisible:
        combinedViewSelectorMenu->setActiveAction(actionShowElements);
        break;
    case resultPresentation::combinedView_undeformedModel:
        combinedViewSelectorMenu->setActiveAction(actionShowUndeformedModel);
        break;
    case resultPresentation::combinedView_undeformedWireFrame:
        combinedViewSelectorMenu->setActiveAction(actionShowUndeformedWireframe);
        break;
    }
}
