#include "dockableViewPort.h"

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

void dockableViewPort::setAction3D_Rotation()
{
    myViewPort->setAction3D_Rotation();
}

void dockableViewPort::setAction3D_Pan()
{
    myViewPort->setAction3D_Pan();
}

void dockableViewPort::setAction3D_WindowZooming()
{
    myViewPort->setAction3D_WindowZooming();
}

void dockableViewPort::FitAll()
{
    myViewPort->FitAll();
}

void dockableViewPort::setSelectionMode(CurSelectionMode selectionMode)
{
    myViewPort->setSelectionMode(selectionMode);
}
