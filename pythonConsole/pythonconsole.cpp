#include "pythonconsole.h"
#include "PythonQt.h"
#include "PythonQt_QtAll.h"
#include "gui/PythonQtScriptingConsole.h"
#include "systemscript.h"
#include "meshscript.h"

//! ---------------------------------------
//! function: constructor
//! details:
//! ---------------------------------------

PythonConsole::PythonConsole(QObject *parent) : QObject(parent)
{
    QMainWindow *mw = static_cast<QMainWindow*>(this->parent());
    SystemScript *system = new SystemScript(mw);
    MeshScript *mesh = new MeshScript(mw);

    PythonQt::init(PythonQt::IgnoreSiteModule | PythonQt::RedirectStdOut);
    PythonQt_QtAll::init();
    PythonQtObjectPtr pyContext = PythonQt::self()->getMainModule();
    pyQtScrCons = new PythonQtScriptingConsole(mw, pyContext);

    pyContext.addObject("system", system);
    pyContext.addObject("mesh", mesh);
    pyQtScrCons->show();
}

PythonQtScriptingConsole *PythonConsole::getConsole()
{
    return pyQtScrCons;
}


