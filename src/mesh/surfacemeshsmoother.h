#ifndef SURFACEMESHSMOOTHER_H
#define SURFACEMESHSMOOTHER_H

//! --------------------------------------
//! replacement for opencascade::handle<>
//! --------------------------------------
#include "occhandle.h"

//! ----------------
//! custom includes
//! ----------------
#include "meshdatabase.h"

//! ----
//! OCC
//! ----
#include <TopoDS_Shape.hxx>
#include "qprogressindicator.h"

//! ---
//! Qt
//! ---
#include <QObject>

class surfaceMeshSmoother: public QObject
{
    Q_OBJECT

public:

    //! constructor
    surfaceMeshSmoother(QProgressIndicator *aProgressIndicator, QObject *parent=0);

    //! perform
    bool perform(meshDataBase *mDB, int bodyIndex, bool rebuildVolumeMesh=true);

    //! set progress indicator
    void setProgressIndicator(QProgressIndicator *aProgressIndicator);

private:

    //! progress indicator
    QProgressIndicator *myProgressIndicator;

signals:

    void requestEnableStop();
    void requestDisableStop();
};

#endif // SURFACEMESHSMOOTHER_H
