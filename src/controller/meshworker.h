#ifndef MESHWORKER_H
#define MESHWORKER_H

//! ---
//! Qt
//! ---
#include <QObject>
#include <QList>

//! ----------------
//! custom includes
//! ----------------
#include <mesherclass.h>
#include <meshdatabase.h>

class QProgressIndicator;

class meshWorker: public QObject
{
    Q_OBJECT

public:

    //! constructor
    meshWorker(meshDataBase *mDB, QList<int> bodies, bool isVolume, QProgressIndicator *aProgressIndicator=0, QObject *parent=0);

    //! destructor
    meshWorker::~meshWorker()
    {
        cout<<"meshWorker::~meshWorker()->____destructor called____"<<endl;
    }

    //! set progress indicator
    void setProgressIndicator(QProgressIndicator *aProgressIndicator);

    //! get progress indicator
    QProgressIndicator* getProgressIndicator() { return myProgressIndicator; }

private:

    MesherClass *theMesher;
    QList<int> myListOfBodies;
    bool myIsVolume;
    QProgressIndicator *myProgressIndicator;

private slots:

    void update();

    //
    void emitRequestStartingNetgenEnquireTimer();
    void emitRequestStoppingNetgenEnquireTimer();

public slots:

    void perform();
    void emitResultReady();

signals:

    void resultReady();

    //!
    void requestStartingNetgenEnquireTimer();
    void requestStoppingNetgenEnquireTimer();
};

#endif // MESHWORKER_H
