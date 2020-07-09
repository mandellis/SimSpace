#ifndef CCXSOLVERMANAGER1_H
#define CCXSOLVERMANAGER1_H

#include <QWidget>
#include <QThread>
#include <QString>

class QProcess;
class QProgressIndicator;

class CCXSolverManager1: public QThread
{
    Q_OBJECT

private:

    void run() Q_DECL_OVERRIDE;
    void copy_file(const char* srce_file, const char* dest_file);

public:

    CCXSolverManager1(QObject *parent = 0);
    void setNbProcessors(int NbProcessors=2);
    void setInputFile(const QString &inputFile);
    void setProgressIndicator(QProgressIndicator *aProgressIndicator);

private:

    QProcess *myProcess;
    QString myInputFileName;
    QProgressIndicator *myProgressIndicator;
    QString myCCXLog;
    QWidget *mySolutionMonitorPanel;
    QWidget *mySimulationManager;
    int myNbProcessors;

private slots:

    void retrieveSolverData();

signals:

    void CCXRunFinished();
};

#endif // CCXSOLVERMANAGER1_H
