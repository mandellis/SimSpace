#include "systemscript.h"
#include <QObject>
#include <QString>

SystemScript::SystemScript(QObject *parent) : QObject(parent)
{
    mw = static_cast<MainWindow*>(this->parent());
}

void SystemScript::importFile(QString address)
{
    mw->importFile(address);
}
