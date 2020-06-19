#ifndef MEMORYPROFILER_H
#define MEMORYPROFILER_H

//! Qt
#include <QDockWidget>

//! ---------------------------------------------------------------------------------------------------------
//! experimental
//! https://docs.microsoft.com/it-it/windows/desktop/psapi/collecting-memory-usage-information-for-a-process
//! ---------------------------------------------------------------------------------------------------------
#include <windows.h>
#include <stdio.h>
#include <psapi.h>

class QTimer;
class QProgressBar;

class memoryProfiler : public QDockWidget
{
    Q_OBJECT

public:

    //! constructor
    memoryProfiler(int msUpdate = 2500, QWidget *parent = 0);

    //! destructor
    ~memoryProfiler();

private:

    //std::vector<double> myDummyVector;

    int myAvailablePhysicalMemory;

    QTimer *myTimer;
    QProgressBar *myProgressBar;

    SIZE_T getProcessMemoryInfo();

    //! get the available physical system memory
    int getSystemMemory();

    //! create the widgets
    void createContent();

private slots:

    void updateProgressBar();

    void handleTimer(bool isWidgetVisible);
};

#endif // MEMORYPROFILER_H
