//! -----------------------------------------------------------------------------------------------------
//! reference
//! https://stackoverflow.com/questions/45724741/qt-c-create-a-global-variable-accessible-to-all-classes
//! -----------------------------------------------------------------------------------------------------
#ifndef GLOBAL_H
#define GLOBAL_H

//! ---
//! Qt
//! ---
#include <QString>

//! ----
//! C++
//! ----
#include <vector>

//! ----------------
//! custom includes
//! ----------------
#include <userMessage.h>
#include "usermessagesmodel.h"

//! ---
//! Qt
//! ---
#include <QStandardItem>

//! ---------------------------
//! codes: "0" stopped or stop
//!        "1" run or running
//! ---------------------------
class Global
{
public:

    //! -----------------------------
    //! "0" stopped or to be stopped
    //! "1" run or to be started
    //! -----------------------------
    int code = 1;

    //! ----------------------------------------------------
    //! percent of an execution
    //! widely used when retrieving netgen meshing progress
    //! ----------------------------------------------------
    double percent = 0.0;

    //! -------------------
    //! key name of a task
    //! -------------------
    QString task ="";

    //! -------------------------------------------------------
    //! this flags indicates if when showing results
    //! the full volume mesh is drawn or only the surface mesh
    //! -------------------------------------------------------
    bool isVolumeMeshShownAsSurface = false;

    //! -------------------------------------------------
    //! the user message model storing the user messages
    //! -------------------------------------------------
    userMessagesModel *myMessages;

    //! -------------------------
    //! current running analysis
    //! -------------------------
    QStandardItem *curRunningAnalysis;

public:

    Global() = default;
    Global(const Global&) = delete;
    Global(Global&&) = delete;

    static Global& status()
    {
        static Global global;
        return global;
    }
};

#endif // GLOBAL_H
