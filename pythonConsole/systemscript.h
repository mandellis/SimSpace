#ifndef SYSTEMSCRIPT_H
#define SYSTEMSCRIPT_H

#include <QObject>
#include "mainwindow.h"

class SystemScript : public QObject
{

    Q_OBJECT

public:

    explicit SystemScript(QObject *parent = nullptr);

public slots:

    void importFile(QString address);

private:

    MainWindow *mw;

};

#endif // SYSTEM_H
