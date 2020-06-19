//! ----------------
//! custom includes
//! ----------------
#include "datasourcebuildercontroller.h"
#include <facedatasourcebuilder.h>

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
    faceDataSourceBuilder *worker = new faceDataSourceBuilder();
    worker->moveToThread(&workerThread);

}
