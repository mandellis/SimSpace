#ifndef QCONSOLEEVENT_H
#define QCONSOLEEVENT_H

#include <QEvent>
#include <QString>

class QConsoleEvent : public QEvent
{

public:

    //! constructor
    QConsoleEvent(QString aMessage=""): QEvent(static_cast<QEvent::Type>(4000)),message(aMessage)
    {
        ;
    }

    //! destructor
    virtual ~QConsoleEvent() { ; }

    QEvent::Type QConsoleEvent::customEventType = static_cast<QEvent::Type>(QEvent::registerEventType());

    static QEvent::Type type()
    {
        return static_cast<QEvent::Type>(4000);
    }

public:

    QString getMessage() {return message;}

private:

    QString message;
};

#endif // QCONSOLEEVENT_H
