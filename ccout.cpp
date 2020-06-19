//! custom includes
#include "ccout.h"
#include <qconsoleevent.h>
#include <iostream>
using namespace std;

//! Qt
#include <QWidget>
#include <QList>
#include <QApplication>

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
ccout::ccout(const QString &aString, const QString &theTargetWidgetName)
{
    //! search for the target widget
    QWidget *widget;
    QList<QWidget*> widgets = QApplication::allWidgets();
    for(int i=0; i<widgets.length(); i++)
    {
        if(widgets.at(i)->objectName()==theTargetWidgetName) widget = widgets.at(i);
    }
    //! post an event containing the aString message
    if(widget!=NULL)
    {
        QConsoleEvent *e = new QConsoleEvent(aString);
        QApplication::postEvent(widget,e);
    }
}
