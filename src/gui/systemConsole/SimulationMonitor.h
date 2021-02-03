#ifndef CCXSYSTEMCONSOLE_H
#define CCXSYSTEMCONSOLE_H

//! custom includes
#include "systemconsole.h"
#include "solutioninfo.h"

//! Qt
#include <QWidget>

#include <iostream>

class SimulationMonitor: public systemConsole
{
    Q_OBJECT

public:

    //! constructor
    explicit SimulationMonitor(bool isCreatedAsRunning = true,
                              bool isMenuBarVisible = true,
                              bool isContinuosLog = false,
                              const QString &myLogFile = QString(),
                              QWidget *parent = 0);

    //! destructor
    virtual ~SimulationMonitor();

    //! set log file name
    virtual void setLogFile(const QString& logFilePathAbsolute)
    {
        myLogFile = logFilePathAbsolute;
    }

protected:

    bool eventFilter(QObject* object, QEvent* event);

private:

signals:

    void SimulationMonitorClosed();
};

#endif // CCXSYSTEMCONSOLE_H
