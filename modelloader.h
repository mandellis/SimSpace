//! -------------------------------------------
//! An interface for loading importing tools
//! 03/08/2018 - support for STEP files
//! ...
//! -------------------------------------------
#ifndef MODELLOADER_H
#define MODELLOADER_H

//! Qt
#include <QObject>
#include <QList>
#include <QString>

//! OCC
#include <TopoDS_Shape.hxx>
#include <TopoDS_Compound.hxx>
#include <TColStd_ListOfAsciiString.hxx>

//! custom includes
#include "qoccprogressindicator.h"

class STEPimporter;

class ModelLoader: public QObject
{
    Q_OBJECT

private:

    QString myFileName;
    STEPimporter *myStepImporter;

    enum fileType
    {
        fileType_STEP,
        fileType_IGES,
        fileType_STL,
        fileType_Parasolid,
        fileType_NONE
    }
    myFileType;

public:

    ModelLoader(const QString &fileName, QObject *parent=0);

public slots:

    bool perform(TopoDS_Compound &shapeFromReader, QList<QString> &listOfNames,
                 const occHandle(QOccProgressIndicator) &aProgressIndicator);


signals:

      void resultReady(bool isDone);
};

#endif // MODELLOADER_H
