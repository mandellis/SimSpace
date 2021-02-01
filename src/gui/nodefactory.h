#ifndef NODEFACTORY_H
#define NODEFACTORY_H

//! ----------------
//! custom includes
//! ----------------
#include "simulationnodeclass.h"
#include <meshdatabase.h>

//! ----
//! OCC
//! ----
#include <AIS_InteractiveContext.hxx>

//! ----
//! Qt
//! ----
#include <QString>
#include <QVariant>

//! ----
//! C++
//! ----
#include <iostream>
using namespace std;

class nodeFactory
{
public:

    nodeFactory();

    //! ----------------
    //! create the node
    //! ----------------
    static SimulationNodeClass* nodeFromScratch(SimulationNodeClass::nodeType type,
                                                meshDataBase *mDB=NULL,
                                                const occHandle(AIS_InteractiveContext)& CTX=NULL,
                                                QVariant addOptions = QVariant());

    //! unused
    static SimulationNodeClass* nodeFromFile(const QString &fileName);

private:

};

#endif // NODEFACTORY_H
