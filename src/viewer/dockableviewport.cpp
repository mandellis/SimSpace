//! ----------------
//! custom includes
//! ----------------
#include "dockableViewPort.h"

//! ---
//! Qt
//! ---
#include <QVBoxLayout>

//! ---------------------
//! function: constuctor
//! details:
//! ---------------------
dockableViewPort::dockableViewPort(QWidget *parent):QDockWidget(parent)
{
    //! --------------------------
    //! create the widget content
    //! --------------------------
    this->createContent();
}

//! ------------------------
//! function: createContent
//! details:
//! ------------------------
void dockableViewPort::createContent()
{
    //! ------------------
    //! widget and layout
    //! ------------------
    myWidget = new QWidget();
    myWidget->setContentsMargins(0,0,0,0);

    QVBoxLayout *vl = new QVBoxLayout(myWidget);
    vl->setContentsMargins(0,0,0,0);
    vl->setSpacing(0);

    //! -------------------
    //! viewport instances
    //! -------------------
    myViewPort = new occPreGLWidget();
    myViewPort->setWindowTitle("Master");
    myViewPort->setContentsMargins(0,0,0,0);

    //! ----------------
    //! add viewe ports
    //! ----------------
    vl->addWidget(myViewPort);

    //! ----------------------------
    //! set the widget for the dock
    //! ----------------------------
    this->setWidget(myWidget);

    //! -------------
    //! focus policy
    //! -------------
    this->setFocusPolicy(Qt::ClickFocus);
}

//! -----------------------
//! function: focusInEvent
//! details:
//! -----------------------
void dockableViewPort::focusInEvent(QFocusEvent *focusEvent)
{
    Q_UNUSED(focusEvent)

    cerr<<"___focus on widget: "<<this->objectName().toStdString()<<"____"<<endl;
    cerr<<"___"<<focusEvent->MouseButtonPress<<"____"<<endl;
    emit requestChangeDelegateContext(this->myViewPort->getContext());
}

//! -------------------------------
//! function: setAction3D_Rotation
//! details:
//! -------------------------------
void dockableViewPort::setAction3D_Rotation()
{
    myViewPort->setAction3D_Rotation();
}

//! --------------------------
//! function: setAction3D_Pan
//! details:
//! --------------------------
void dockableViewPort::setAction3D_Pan()
{
    myViewPort->setAction3D_Pan();
}

//! ------------------------------------
//! function: setAction3D_WindowZooming
//! details:
//! ------------------------------------
void dockableViewPort::setAction3D_WindowZooming()
{
    myViewPort->setAction3D_WindowZooming();
}

//! -----------------
//! function: FitAll
//! details:
//! -----------------
void dockableViewPort::FitAll()
{
    myViewPort->FitAll();
}

//! ---------------------------
//! function: setSelectionMode
//! details:
//! ---------------------------
void dockableViewPort::setSelectionMode(CurSelectionMode selectionMode)
{
    myViewPort->setSelectionMode(selectionMode);
}

//! ---------------------
//! function: displayCAD
//! details:
//! ---------------------
void dockableViewPort::displayCAD()
{
    myViewPort->displayCAD();
}

//! -------------------
//! function: hideBody
//! details:
//! -------------------
void dockableViewPort::hideBody(const TColStd_ListOfInteger &bodyListNbs)
{
    myViewPort->hideBody(bodyListNbs);
}

//! -------------------
//! function: showBody
//! details:
//! -------------------
void dockableViewPort::showBody(const TColStd_ListOfInteger &bodyListNbs)
{
    myViewPort->showBody(bodyListNbs);
}

//! ---------------------
//! function: getContext
//! details:
//! ---------------------
const occHandle(AIS_InteractiveContext)& dockableViewPort::getContext() const
{
    return myViewPort->getContext();
}

//! ---------------------------
//! function: displayShapeCopy
//! details:
//! ---------------------------
void dockableViewPort::displayShapeCopy(const TopTools_ListOfShape &list1,
                      const TopTools_ListOfShape &list2,
                      Quantity_NameOfColor color1,
                      Quantity_NameOfColor color2,
                      QVariant options)
{
    myViewPort->displayShapeCopy(list1,list2,color1,color2,options);
}

//! ----------------------------
//! function: displayShapeCopy1
//! details:
//! ----------------------------
void dockableViewPort::displayShapeCopy1(const TopTools_ListOfShape &listShapes, Quantity_NameOfColor color)
{
    myViewPort->displayShapeCopy1(listShapes,color);
}

void dockableViewPort::displayShapeCopy2(const TopTools_ListOfShape &listShapes, Quantity_NameOfColor color)
{
    myViewPort->displayShapeCopy2(listShapes,color);
}
