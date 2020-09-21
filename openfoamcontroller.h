#ifndef OPENFOAMCONTROLLER_H
#define OPENFOAMCONTROLLER_H

//! Qt
#include <QObject>
#include <QThread>
#include <customtablemodel.h>

#include "simulationnodeclass.h"

class QProgressIndicator;

class openFoamController: public QObject
{
    Q_OBJECT

private:

    QThread workerThread;
    QProgressIndicator *myProgressIndicator;
    SimulationNodeClass *myNode;
    QString myTargetDirectory;
#ifdef COSTAMP_VERSION
    std::vector<double> myTimeFolders;
#endif
public:

    //! constructor
    openFoamController(const QString &sourceDirPath,
                       const QString &targetDirPath,
                       int fileMode,
                       QProgressIndicator *aProgressIndicator,
                       QObject *parent);

    //! set progress indicator
    void setProgressIndicator(QProgressIndicator *aProgressIndicator) { myProgressIndicator = aProgressIndicator; }
#ifdef COSTAMP_VERSION
    void setTimeFolders(const std::vector<double> timeFolders)
    {myTimeFolders = timeFolders;}
#endif
private slots:

    //! for stopping the thread
    void stopThread();

    //! lock node
    void lockNode(SimulationNodeClass *node);



signals:

    //! for starting the thread
    void operate(SimulationNodeClass *node);

    //! task finished or stopped - experimental
    void taskFinishedOrStopped();
};

#endif // OPENFOAMCONTROLLER_H
