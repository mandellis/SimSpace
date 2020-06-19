#ifndef CCXSOLVERMANAGER_H
#define CCXSOLVERMANAGER_H

//! ---
//! Qt
//! ---
#include <QWidget>
#include <QDialog>
#include <QThread>
#include <QEvent>
#include <QProcess>

//! ----------------
//! custom includes
//! ----------------
#include "qprogressevent.h"

//! ----
//! C++
//! ----
#include <iostream>
using namespace std;


class QProgressBar;
class QPushButton;

class CCXSolverManager: public QDialog
{
    Q_OBJECT

private:

    QThread workerThread;
    QWidget *myTargetWidget;
    QString myInputFile;
    QProgressBar *prgbar;
    QPushButton *myStopButton;
    int myNbProcs=1;

private:

    void createContent();

public:

    //! constructor
    CCXSolverManager(const QString &inputFileName, QWidget *parent=0);

    //! destructor
    ~CCXSolverManager();

    //! set the number of processors
    void setNbProcessors(int NbProcs) { myNbProcs = NbProcs; }

protected:

    virtual bool eventFilter(QObject *object, QEvent *event);

signals:

    void operate(int NbProcs);
    //void solutionReady();

    void solutionReady(const QString &whereIsSolution=QString());

private slots:

    void terminateRun();
    void closeManager(int code);

public slots:

    //! set the input file
    void setInputFile(const QString &inputFile);

    //! function that emits a signal if the solution process has been successfull
    void emitSolutionProcessFinished(int code);

signals:

    //! emit stop timer
    void stopTimer();
};

#endif // CCXSOLVERMANAGER_H
