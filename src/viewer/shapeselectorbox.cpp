//! ----------------
//! custom includes
//! ----------------
#include "shapeselectorbox.h"

//! ---
//! Qt
//! ---
#include <QEvent>

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
ShapeSelectorBox::ShapeSelectorBox(QWidget *parent): QLineEdit(parent)
{
    //! install an event filter for the mouse click
    installEventFilter(this);
}

//! -----------------------------------------------
//! function: install event filter for mouse click
//! details:
//! -----------------------------------------------
bool ShapeSelectorBox::eventFilter(QObject* object, QEvent* event)
{
    //! left mouse button pressed
    if(object == this && event->type() == QEvent::MouseButtonPress)
    {
        emit activateSelector();
        return false; // lets the event continue to the edit
    }
    return QObject::event(event);
}
