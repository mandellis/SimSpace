#ifndef QCCXSOLVERMESSAGEEVENT_H
#define QCCXSOLVERMESSAGEEVENT_H

//! ---
//! Qt
//! ---
#include <QEvent>
#include <QString>

//! ----------------
//! custom includes
//! ----------------
#include "ccxsolvermessage.h"

class QCCXSolverMessageEvent: public QEvent
{

public:

    //! ------------
    //! constructor
    //! ------------
    QCCXSolverMessageEvent(CCXSolverMessage message):
        QEvent(static_cast<QEvent::Type>(2002))
    {
        myCCXSolverMessage = message;
    }

    //! -----------
    //! destructor
    //! -----------
    virtual ~QCCXSolverMessageEvent() { ; }

    QEvent::Type QCCXSolverMessageEvent::customEventType = static_cast<QEvent::Type>(QEvent::registerEventType());

    static QEvent::Type type()
    {
        return static_cast<QEvent::Type>(2002);
    }

public:

    CCXSolverMessage getMessage() {return myCCXSolverMessage;}

private:

    CCXSolverMessage myCCXSolverMessage;
};

#endif // QCCXSOLVERMESSAGEEVENT_H
