#ifndef QBACKGROUNDEVENT_H
#define QBACKGROUNDEVENT_H

#include <QEvent>
#include <QColor>

class QBackgroundEvent : public QEvent
{

public:

    //! constructor
    QBackgroundEvent(int aGradient, QColor aFirstColor, QColor aSecondColor):
        QEvent(static_cast<QEvent::Type>(3000)),
        gradient(aGradient),
        firstColor(aFirstColor),
        secondColor(aSecondColor)
    {
        ;
    }

    //! destructor
    virtual ~QBackgroundEvent() { ; }

    QEvent::Type QBackgroundEvent::customEventType = static_cast<QEvent::Type>(QEvent::registerEventType());

    static QEvent::Type type()
    {
        return static_cast<QEvent::Type>(3000);
    }

public:

    int getGradient() {return gradient; }
    QColor getFirstColor() {return firstColor; }
    QColor getSecondColor() {return secondColor; }

    int R1() {return firstColor.red(); }
    int G1() {return firstColor.green(); }
    int B1() {return firstColor.blue(); }

    int R2() {return secondColor.red(); }
    int G2() {return secondColor.green(); }
    int B2() {return secondColor.blue(); }

private:

    int gradient;
    QColor firstColor;
    QColor secondColor;
};

#endif // QBACKGROUNDEVENT_H
