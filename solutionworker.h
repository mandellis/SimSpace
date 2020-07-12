#ifndef SOLUTIONWORKER_H
#define SOLUTIONWORKER_H

//! ---
//! Qt
//! ---
#include <QObject>
#include <QProcess>
#include <QString>

//! ----
//! C++
//! ----
#include <iostream>
#include <fstream>
using namespace std;

//! --
//! C
//! --
#include <stdio.h>

class QTimer;

class SolutionWorker : public QObject
{
    Q_OBJECT

public:

    explicit SolutionWorker(const QString &inputFileName, QObject *parent = 0);
    ~SolutionWorker();

public slots:

    void perform(int NbProcs);
    void emitResultReady(int code);
    void retrieveSolverData();

public:

    QString messages;

signals:

    void resultReady(int code);
    void resultReady1(int code, QProcess::ExitStatus es);

private:

    QProcess *shellProcess;
    //QTimer *myTimer;
    QString myInputFileName;
    QWidget *mySolutionMonitorPanel;
    QWidget *mySimulationManager;
    QString CCX_PATH;

private:

private slots:

    //void sendEnter();

public slots:

};

#endif // SOLUTIONWORKER_H
