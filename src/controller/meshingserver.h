#ifndef MESHINGSERVER_H
#define MESHINGSERVER_H

//! ---
//! Qt
//! ---
#include <QObject>
#include <QThread>

//! ----------------
//! custom includes
//! ----------------
#include <meshdatabase.h>
#include <userMessage.h>

//! ----
//! C++
//! ----
#include <vector>

class QTimer;
class QProgressIndicator;

class MeshingServer: public QObject
{

    Q_OBJECT

private:

    QThread workerThread;
    QProgressIndicator *myProgressIndicator;
    bool myIsMeshingRunning;

    //! ------------------------------------------------------------
    //! this timer is used to trigger the communication with netgen
    //! and retrieving its progress information
    //! ------------------------------------------------------------
    QTimer *myNetgenInquireTimer;

public:

    //! ------------
    //! constructor
    //! ------------
    MeshingServer(QProgressIndicator *aProgressIndicator, QObject *parent=0);

    //! ------------
    //! constructor
    //! ------------
    MeshingServer(meshDataBase *mDB, const QList<int> &listOfBodies, bool isVolume, QProgressIndicator *aProgressIndicator, QObject *parent=0);

    //! -----
    //! init
    //! -----
    void init(meshDataBase *mDB, const QList<int> &listOfBodies, bool isVolume,
              QProgressIndicator *aProgressIndicator = Q_NULLPTR);

    //! -----------
    //! destructor
    //! -----------
    ~MeshingServer()
    {
        cout<<"MeshingServer::~MeshingServer()->____DESTRUCTOR CALLED____"<<endl;
        workerThread.wait();
    }

    //! -----------
    //! get status
    //! -----------
    bool isMeshingRunning();

signals:

    void operate();
    void meshingFinished(bool isSuccessfull);

public slots:

    void handleResult();
    void stopThread();

    //! ---------------------
    //! show netgen progress
    //! ---------------------
    void showNetgenProgress();

    //! ---------------------
    //! netgen enquire timer
    //! ---------------------
    void stopNetgenEnquireTime();
    void startNetgenEnquireTime();

    //! ------
    //! abort
    //! ------
    void abortSTLNetgen();
};

#endif // MESHINGSERVER_H
