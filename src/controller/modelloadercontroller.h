#ifndef MODELLOADERCONTROLLER_H
#define MODELLOADERCONTROLLER_H

//! Qt
#include <QThread>
#include <QObject>
#include <QString>

//! OCC
#include <TopoDS_Shape.hxx>
#include <TColStd_ListOfAsciiString.hxx>

//! custom includes
#include "src/utils/modelloader.h"
#include "qoccprogressindicator.h"

class ModelLoaderController: public QObject
{
    Q_OBJECT

private:

    QThread modelLoaderThread;
    TopoDS_Shape &myShapeFromReader;
    TColStd_ListOfAsciiString &myListOfNames;
    occHandle(QOccProgressIndicator) myProgressIndicator;

public:

    //! constructor
    ModelLoaderController(const QString &fileName,
                          TopoDS_Shape &shapeFromReader,
                          TColStd_ListOfAsciiString &listOfNames,
                          const occHandle(QOccProgressIndicator) aPi,
                          QObject *parent=0);

    //! destructor
    ~ModelLoaderController()
    {
        modelLoaderThread.quit();
        modelLoaderThread.wait();
    }

public slots:

    void handleResults(bool isDone)
    {
        bool a = isDone;
    }

signals:

    void operate();
};

#endif // MODELLOADERCONTROLLER_H
