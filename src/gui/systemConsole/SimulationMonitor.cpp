//! ----------------
//! custom includes
//! ----------------
#include "SimulationMonitor.h"
#include "qprogressevent.h"
#include "qconsoleevent.h"
#include "src/utils/tools.h"
//#include "ccxconsoletofile.h"
#include "src/utils/ccout.h"

//! ---
//! Qt
//! ---
#include <QApplication>

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
SimulationMonitor::SimulationMonitor(bool isCreatedAsRunning,
                                   bool isMenuBarVisible,
                                   bool isContinuosLog,
                                   const QString &myLogFile,
                                   QWidget *parent):
    systemConsole(isCreatedAsRunning, isMenuBarVisible, isContinuosLog, myLogFile, parent)
{
    cout<<"SimulationMonitor::SimulationMonitor()->____CONSTRUCTOR CALLED____"<<endl;

    //! install event filter
    this->installEventFilter(this);
}

//! ---------------------
//! function: destructor
//! details:
//! ---------------------
SimulationMonitor::~SimulationMonitor()
{
    cout<<"SimulationMonitor::~SimulationMonitor()->____DESTRUCTOR CALLED____"<<endl;
    emit SimulationMonitorClosed();
}

//! ----------------------
//! function: eventFilter
//! details:
//! ----------------------
bool SimulationMonitor::eventFilter(QObject* object, QEvent* event)
{
    if(object == this && event->type() == QConsoleEvent::type())
    {
        if(this->isRunning)
        {
            //cout<<"SimulationMonitor::eventFilter()->____event received____"<<endl;            
            //! ------------------------------
            //! append the text to the viewer
            //! ------------------------------
            QString text = static_cast<QConsoleEvent*>(event)->getMessage();
            myTextEdit->appendPlainText(text);
        }

        //! the event is not stopped here
        return false;
    }
    else
    {
        //! standard event processing
        return QObject::eventFilter(object, event);
    }
}
