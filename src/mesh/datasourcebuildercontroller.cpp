//! ----------------
//! custom includes
//! ----------------
#include "datasourcebuildercontroller.h"
#include <datasourcebuilder.h>

//! ---
//! Qt
//! ---
#include <QThread>

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
dataSourceBuilderController::dataSourceBuilderController()
{    
    //! --------------
    //! Worker thread
    //! --------------
    dataSourceBuilder *worker = new dataSourceBuilder();
    worker->moveToThread(&workerThread);

}
