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
#include "resultpresentation.h"

//! ---
//! Qt
//! ---
#include <QStandardItem>

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

    //! --------------------
    //! result presentation
    //! --------------------
    resultPresentation myResultPresentation;

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
