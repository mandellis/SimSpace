#ifndef TETHEX_H
#define TETHEX_H

//! ----------------
//! custom includes
//! ----------------
#include "occhandle.h"
#include "meshdatabase.h"
#include <ng_meshvs_datasource3d.h>
#include <ng_meshvs_datasource2d.h>
#include "qprogressindicator.h"
#include <userMessage.h>

//! ---
//! Qt
//! ---
#include <QString>

class TetHex
{
public:

    //! ----------------------
    //! constructor - default
    //! ----------------------
    TetHex()
    {
        myTask = QString("TetHex conversion");
    }

    //! ------------
    //! constructor
    //! ------------
    TetHex(meshDataBase *aMeshDB):myMDB(aMeshDB)
    {
        myTask = QString("TetHex conversion");
    }

    //! ------------
    //! setDataBase
    //! ------------
    void setDataBase(meshDataBase *aMDB) { myMDB = aMDB; }

    //! --------
    //! perform
    //! --------
    userMessage perform(int bodyIndex);

    //! ---------
    //! makeHexa
    //! ---------
    userMessage makeHexa(const occHandle(Ng_MeshVS_DataSource3D) &aTetMesh, occHandle(Ng_MeshVS_DataSource3D) &HexaMeshDS);

    //! -----------------------
    //! set progress indicator
    //! -----------------------
    void setProgressIndicator(QProgressIndicator *aProgressIndicator)
    {
        if(aProgressIndicator!=Q_NULLPTR) myProgressIndicator = aProgressIndicator;
    }

private:

    //! ----------
    //! task name
    //! ----------
    QString myTask;

    //! ---------------
    //! mesh data base
    //! ---------------
    meshDataBase *myMDB;

    //! ---------------------------------------------------------
    //! the progress indicator - this class is not a QObject:
    //! if you want to disable the "Stop" button of the progress
    //! indicator an event is posted
    //! ---------------------------------------------------------
    QProgressIndicator *myProgressIndicator;

private:

    //! --------------------------------------------------
    //! enable/disable the progress indicator stop button
    //! --------------------------------------------------
    void setStopButtonEnabled(bool isEnabled);
};

#endif // TETHEX_H
