#ifndef QPROGRESSEVENT_H
#define QPROGRESSEVENT_H

#include <QEvent>
#include <QString>

class QThread;

enum QProgressEventAction
{
    QProgressEvent_Init,
    QProgressEvent_Reset,
    QProgressEvent_Update,
    QProgressEvent_EnableStop,
    QProgressEvent_DisableStop,
    QProgressEvent_None
};

class QProgressEvent : public QEvent
{

public:

    //! ------------
    //! constructor
    //! ------------
    QProgressEvent(QProgressEventAction action = QProgressEvent_Reset, int min=0, int max=100, int val=0, QString aMessage="",
                   QProgressEventAction action1 = QProgressEvent_Reset, int min1=0, int max1=100, int val1=0, QString currentRunningtask ="",
                   QObject *object =0):
        QEvent(static_cast<QEvent::Type>(2000)),
        theAction(action),
        theMin(min),
        theMax(max),
        theVal(val),
        theAction1(action1),
        theMin1(min1),
        theMax1(max1),
        theVal1(val1),
        myMessage(aMessage),
        myTask(currentRunningtask),
        myObject(object)
    {
        ;
    }

    //! -----------
    //! destructor
    //! -----------
    virtual ~QProgressEvent() { ; }

    QEvent::Type QProgressEvent::customEventType = static_cast<QEvent::Type>(QEvent::registerEventType());

    static QEvent::Type type()
    {
        return static_cast<QEvent::Type>(2000);
    }

public:

    QProgressEventAction action() { return theAction; }
    QProgressEventAction action1() { return theAction1; }

    inline void setAction(QProgressEventAction anAction) { theAction = anAction; }
    inline void setAction1(QProgressEventAction anAction1) { theAction1 = anAction1; }

    inline void setRange(int min, int max) { theMin = min; theMax = max;}
    inline void setRange1(int min1, int max1) { theMin1 = min1; theMax1 = max1;}

    inline int min() { return theMin; }
    inline int max() { return theMax; }
    inline int val() { return theVal; }

    inline int min1() { return theMin1; }
    inline int max1() { return theMax1; }
    inline int val1() { return theVal1; }

    void setVal (int val)
    {
        theAction = QProgressEventAction::QProgressEvent_Update;
        theVal = val;
    }
    void setVal1 (int val1)
    {
        theAction1 = QProgressEventAction::QProgressEvent_Update;
        theVal1 = val1;
    }

    inline void setMessage(const QString &message) { myMessage = message; }

    inline QString getMessage() { return myMessage; }
    inline QString getTask() { return myTask; }
    inline QObject *getObject() { return myObject; }

private:

    QProgressEventAction theAction;
    int theMin;
    int theMax;
    int theVal;

    QProgressEventAction theAction1;
    int theMin1;
    int theMax1;
    int theVal1;

    QString myMessage;
    QString myTask;
    QObject *myObject;
};

#endif // QPROGRESSEVENT_H
