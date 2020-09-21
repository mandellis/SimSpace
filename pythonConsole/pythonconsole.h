#ifndef PYTHONCONSOLE_H
#define PYTHONCONSOLE_H

//! Qt
#include <QObject>
#include <QMainWindow>

//! PythonQt
#include "gui/PythonQtScriptingConsole.h"

using namespace std;

class PythonConsole : public QObject
{
    friend class MeshScript;

    friend class SystemScript;

    Q_OBJECT

public:

    //! constructor
    explicit PythonConsole(QObject *parent);

    PythonQtScriptingConsole *getConsole();

private:

    PythonQtScriptingConsole *pyQtScrCons;

};

#endif // PYTHONCONSOLE_H
