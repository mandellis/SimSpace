//! Qt
#include <QMetaType>

//! custom includes
#include "modelloadercontroller.h"
#include "modelloader.h"

//! OCC
#include <TopoDS_Shape.hxx>

ModelLoaderController::ModelLoaderController(const QString &fileName,
                                             TopoDS_Shape &shapeFromReader,
                                             TColStd_ListOfAsciiString &listOfNames,
                                             const opencascade::handle<QOccProgressIndicator> aPi,
                                             QObject *parent):
    QObject(parent),
    myShapeFromReader(shapeFromReader),
    myListOfNames(listOfNames),
    myProgressIndicator(aPi)
{
    //! the model loader has no parent, since it has to be moved
    //! on a different thread
    ModelLoader *aModelLoader = new ModelLoader(fileName);
    aModelLoader->moveToThread(&modelLoaderThread);

    connect(&modelLoaderThread,SIGNAL(finished()),aModelLoader,SLOT(deleteLater()));
    connect(this,SIGNAL(operate()),aModelLoader,SLOT(perform(shapeFromReader,listOfNames,aPi)));
    connect(aModelLoader,SIGNAL(resultReady(bool isDone)),this,SLOT(handleResults(bool isDone)));

    //! start the event loop of the thread
    modelLoaderThread.start();
}
