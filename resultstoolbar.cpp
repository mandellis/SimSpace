//! ----------------
//! custom includes
//! ----------------
#include "resultstoolbar.h"
#include <scaleselector.h>
#include "global.h"

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

    this->setContentsMargins(0,0,0,0);

    //! ---------------------
    //! type of presentation
    //! ---------------------
    typeOfPresentationButton = new QPushButtonExtended(this);
    typeOfPresentationButton->setIcon(QIcon(":/icons/icon_isostrips.png"));

    //! -----------------------------
    //! create a menu for the button
    //! -----------------------------
    typeOfPresentationMenu = new QMenu(this);

    actionUseIsoStrips = typeOfPresentationMenu->addAction("Countouring");
    actionUseIsoStrips->setIcon(QIcon(":/icons/icon_isostrips.png"));
    actionUseIsoStrips->setToolTip("Iso-strips");
    connect(actionUseIsoStrips,SIGNAL(triggered(bool)),this,SLOT(emitRequstUseIsoStrips()));

    //! icon stolen from FREECAD - please change it
    actionUseIsoSurfaces = typeOfPresentationMenu->addAction("Isosurfaces");
    actionUseIsoSurfaces->setIcon(QIcon(":/icons/icon_isosurfaces.svg"));
    actionUseIsoSurfaces->setToolTip("Iso-surfaces");
    connect(actionUseIsoSurfaces,SIGNAL(triggered(bool)),this,SLOT(emitRequestUseIsoSurface()));

    actionUseSmoothNodal = typeOfPresentationMenu->addAction("Smooth nodal");
    actionUseSmoothNodal->setIcon(QIcon(":/icons/icon_color gradient 1.jpg"));
    actionUseSmoothNodal->setToolTip("Use nodal results");
    connect(actionUseSmoothNodal,SIGNAL(triggered(bool)),this,SLOT(emitRequestUseSmoothNodal()));

    actionUseIsoLines = typeOfPresentationMenu->addAction("Isolines");
    actionUseIsoLines->setIcon(QIcon(":/icons/icon_isolines.png"));
    actionUseIsoLines->setToolTip("Iso-lines");
    connect(actionUseIsoLines,SIGNAL(triggered(bool)),this,SLOT(emitRequestUseIsoLines()));

    typeOfPresentationButton->setMenu(typeOfPresentationMenu);

    this->addWidget(typeOfPresentationButton);

    //! ------------------------------
    //! add the scale selector widget
    //! ------------------------------
    ScaleSelector *scaleSelector = new ScaleSelector();

    //! ----------------------
    //! set the initial scale
    //! ----------------------
    this->addWidget(scaleSelector);
    connect(scaleSelector,SIGNAL(scaleChanged(double)),this,SLOT(emitRequestUpdatePostObjectScale(double)));

    //! -----------------------
    //! create a button widget
    //! -----------------------
    combinedViewSelectorButton = new QPushButtonExtended(this);
    combinedViewSelectorButton->setIcon(QIcon(":/icons/icon_no wireframe.png"));

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
    Global::status().myResultPresentation.theCombinedView = resultPresentation::combinedView_resultOnly;
    emit requestUpdateViewerStatus();
}

//! -----------------------------------------
//! function: emitRequestShowUndeformedModel
//! details:
//! -----------------------------------------
void ResultsToolBar::emitRequestShowUndeformedModel()
{
    cout<<"ResultsToolBar()->____emit request show undeformed model____"<<endl;
    combinedViewSelectorButton->setIcon(QIcon(":/icons/icon_undeformed model.png"));
    Global::status().myResultPresentation.theCombinedView = resultPresentation::combinedView_undeformedModel;
    emit requestUpdateViewerStatus();
}

//! ---------------------------------------------
//! function: emitRequestShowUndeformedWireframe
//! details:
//! ---------------------------------------------
void ResultsToolBar::emitRequestShowUndeformedWireframe()
{
    cout<<"ResultsToolBar()->____emit request show undeformed wireframe____"<<endl;
    combinedViewSelectorButton->setIcon(QIcon(":/icons/icon_wireframe.png"));
    Global::status().myResultPresentation.theCombinedView = resultPresentation::combinedView_undeformedWireFrame;
    emit requestUpdateViewerStatus();
}

//! --------------------------------------
//! function: emitRequestShowMeshElements
//! details:
//! --------------------------------------
void ResultsToolBar::emitRequestShowMeshElements()
{
    cout<<"ResultsToolBar()->____emit request show mesh elements____"<<endl;
    combinedViewSelectorButton->setIcon(QIcon(":/icons/icon_show elements.png"));
    Global::status().myResultPresentation.theCombinedView = resultPresentation::combinedView_meshVisible;
    emit requestUpdateViewerStatus();
}

//! -------------------------------------------
//! function: emitRequestUpdatePostObjectScale
//! details:
//! -------------------------------------------
void ResultsToolBar::emitRequestUpdatePostObjectScale(double scale)
{
    cout<<"ResultsToolBar::emitRequestUpdatePostObjectScale()->____function called. Scale: "<<scale<<"____"<<endl;
    Global::status().myResultPresentation.theScale = scale;
    emit requestUpdateViewerStatus();
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

//! ---------------------------------
//! function: emitRequstUseIsoStrips
//! details:
//! ---------------------------------
void ResultsToolBar::emitRequstUseIsoStrips()
{
    typeOfPresentationButton->setIcon(QIcon(":/icons/icon_isostrips.png"));
    Global::status().myResultPresentation.theTypeOfPresentation = resultPresentation::typeOfPresentation_isostrips;
    emit requestUpdateViewerStatus();
}

//! ----------------------------------
//! function: emitRequstUseIsoSurface
//! details:
//! ----------------------------------
void ResultsToolBar::emitRequestUseIsoSurface()
{
    typeOfPresentationButton->setIcon(QIcon(":/icons/icon_isosurfaces.svg"));
    Global::status().myResultPresentation.theTypeOfPresentation = resultPresentation::typeOfPresentation_isosurfaces;
    emit requestUpdateViewerStatus();
}

//! -----------------------------------
//! function: emitRequstUseSmoothNodal
//! details:
//! -----------------------------------
void ResultsToolBar::emitRequestUseSmoothNodal()
{
    typeOfPresentationButton->setIcon(QIcon(":/icons/icon_color gradient 1.jpg"));
    Global::status().myResultPresentation.theTypeOfPresentation = resultPresentation::typeOfPresentation_nodalresults;
    emit requestUpdateViewerStatus();
}

//! --------------------------------
//! function: emitRequstUseIsoLines
//! details:
//! --------------------------------
void ResultsToolBar::emitRequestUseIsoLines()
{
    typeOfPresentationButton->setIcon(QIcon(":/icons/icon_isolines.png"));
    Global::status().myResultPresentation.theTypeOfPresentation = resultPresentation::typeOfPresentation_isolines;
    emit requestUpdateViewerStatus();
}

//! --------------------
//! function: setStatus
//! details:
//! --------------------
void ResultsToolBar::setStatus(const resultPresentation &aResultPresentation)
{
    switch(aResultPresentation.theCombinedView)
    {
    case resultPresentation::combinedView_resultOnly: combinedViewSelectorMenu->setActiveAction(actionNoWireframe); break;
    case resultPresentation::combinedView_meshVisible: combinedViewSelectorMenu->setActiveAction(actionShowElements); break;
    case resultPresentation::combinedView_undeformedModel: combinedViewSelectorMenu->setActiveAction(actionShowUndeformedModel); break;
    case resultPresentation::combinedView_undeformedWireFrame: combinedViewSelectorMenu->setActiveAction(actionShowUndeformedWireframe); break;
    }
}
