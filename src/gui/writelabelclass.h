#ifndef WRITELABELCLASS_H
#define WRITELABELCLASS_H

//! redefinition of opencascade Handle() macro
#include "occhandle.h"

#include <AIS_InteractiveContext.hxx>
#include <AIS_TextLabel.hxx>
#include <QString>
#include <QObject>

class writeLabelClass: public QObject
{
    Q_OBJECT

private:

    occHandle(AIS_InteractiveContext) myCTX;
    occHandle(AIS_TextLabel) myText;

public:

    ~writeLabelClass();
    writeLabelClass(QObject *parent=0);
    writeLabelClass(const occHandle(AIS_InteractiveContext) &aCTX, QObject *parent=0);
    void setContext(const occHandle(AIS_InteractiveContext) &aCTX);

public slots:

    void setText(const QString &aText);
};

#endif // WRITELABELCLASS_H
