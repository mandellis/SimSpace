#define MAX_NUMBER_OF_GRAPH_POINT 20

//! custom includes
#include "memoryprofiler.h"

//! Qt
#include <QVBoxLayout>
#include <QProgressBar>
#include <QTimer>
#include <QApplication>
#include <QSizeGrip>

//! C++
#include <iostream>
using namespace std;

//! ------------------------
//! function: createContent
//! details:
//! ------------------------
void memoryProfiler::createContent()
{
    QWidget *container = new QWidget(this);
    container->setContentsMargins(0,0,0,0);

    QVBoxLayout *vl = new QVBoxLayout();
    vl->setContentsMargins(0,0,0,0);
    vl->setMargin(0);

    myProgressBar = new QProgressBar(this);
    myProgressBar->setContentsMargins(0,0,0,0);
    myProgressBar->setRange(0,myAvailablePhysicalMemory);
    myProgressBar->setAlignment(Qt::AlignCenter);
    myProgressBar->setFormat("Memory %v [MB]");

    vl->addWidget(myProgressBar);

    QSizeGrip *sizeGrip = new QSizeGrip(this);
    //vl->addWidget(sizeGrip);
    vl->addWidget(sizeGrip, 0, Qt::AlignBottom | Qt::AlignRight);

    container->setLayout(vl);
    this->setWidget(container);
}

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
memoryProfiler::memoryProfiler(int msUpdate, QWidget *parent):QDockWidget(parent)
{
    myAvailablePhysicalMemory = this->getSystemMemory();

    this->createContent();
    myTimer = new QTimer(this);
    myTimer->setInterval(msUpdate);
    //myTimer->start();

    connect(myTimer,SIGNAL(timeout()),this,SLOT(updateProgressBar()));
    connect(this,SIGNAL(visibilityChanged(bool)),this,SLOT(handleTimer(bool)));
}

//! ----------------------------
//! function: updateProgressBar
//! details:
//! ----------------------------
void memoryProfiler::updateProgressBar()
{
    // for testing purposes
    // for(int k=0; k<10000000; k++) myDummyVector.push_back(0.0);

    SIZE_T memory = getProcessMemoryInfo();

    //! ----------------------------------------------------
    //! recompute the free system memory every 5 iterations
    //! ----------------------------------------------------
    int inUse = memory/(1024*1024);
    myProgressBar->setValue(inUse);

    //cout<<"____memory [MB]: "<<memory<<"____"<<endl;
}

//! -------------------------------------------------------------------------------------------------------------------
//! function: getProcessMemoryInfo
//! details:  https://docs.microsoft.com/it-it/windows/desktop/psapi/collecting-memory-usage-information-for-a-process
//! -------------------------------------------------------------------------------------------------------------------
SIZE_T memoryProfiler::getProcessMemoryInfo()
{
    //cout<<"____APPLICATION PID: "<<QApplication::applicationPid()<<endl;
    DWORD processID = QApplication::applicationPid();
    HANDLE hProcess;
    PROCESS_MEMORY_COUNTERS pmc;
    hProcess = OpenProcess(  PROCESS_QUERY_INFORMATION |
                             PROCESS_VM_READ,
                             FALSE, processID );
    hProcess = OpenProcess(  PROCESS_QUERY_INFORMATION |
                                        PROCESS_VM_READ,
                                        FALSE, processID );
    if (NULL == hProcess)
        return 0;

    if (GetProcessMemoryInfo(hProcess,&pmc,sizeof(pmc)))
    {
        /*
        printf( "\tPageFaultCount: 0x%08X\n", pmc.PageFaultCount );
        printf( "\tPeakWorkingSetSize: 0x%08X\n",
                pmc.PeakWorkingSetSize );
        printf( "\tWorkingSetSize: 0x%08X\n", pmc.WorkingSetSize );
        printf( "\tQuotaPeakPagedPoolUsage: 0x%08X\n",
                pmc.QuotaPeakPagedPoolUsage );
        printf( "\tQuotaPagedPoolUsage: 0x%08X\n",
                pmc.QuotaPagedPoolUsage );
        printf( "\tQuotaPeakNonPagedPoolUsage: 0x%08X\n",
                pmc.QuotaPeakNonPagedPoolUsage );
        printf( "\tQuotaNonPagedPoolUsage: 0x%08X\n",
                pmc.QuotaNonPagedPoolUsage );
        printf( "\tPagefileUsage: 0x%08X\n", pmc.PagefileUsage );
        printf( "\tPeakPagefileUsage: 0x%08X\n",
                pmc.PeakPagefileUsage );
                */
    }
    CloseHandle(hProcess);
    return pmc.WorkingSetSize;
}

//! ----------------------------------
//! function: getSystemMemory
//! details:  get the physical memory
//! ----------------------------------
int memoryProfiler::getSystemMemory()
{
    MEMORYSTATUSEX memory_status;
    ZeroMemory(&memory_status, sizeof(MEMORYSTATUSEX));
    memory_status.dwLength = sizeof(MEMORYSTATUSEX);
    int availablePhysicalSystemMemory = 0;
    if (GlobalMemoryStatusEx(&memory_status))
    {
        availablePhysicalSystemMemory = memory_status.ullAvailPhys/(1024*1024);
        printf("____available physical system memory %d [MB]____\n",availablePhysicalSystemMemory);
    }
    else
    {
        printf("____Unknown RAM____\n");
    }
    return availablePhysicalSystemMemory;
}

//! ----------------------
//! function: destructor
//! details:
//! ----------------------
memoryProfiler::~memoryProfiler()
{
    cout<<"memoryProfiler::~memoryProfiler()->____destructor called____"<<endl;
}

//! ----------------------
//! function: handleTimer
//! details:
//! ----------------------
void memoryProfiler::handleTimer(bool isWidgetVisible)
{
    if(isWidgetVisible)
    {
        if(!myTimer->isActive()) myTimer->start();
    }
    else
    {
        if(myTimer->isActive()) myTimer->stop();
    }
}
