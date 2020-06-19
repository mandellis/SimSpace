#ifndef STEPIMPORTER_H
#define STEPIMPORTER_H

//! redefinition of opencascade Handle() macro
#include "occhandle.h"

//! Qt
#include <QObject>
#include <QString>

//! OCC
#include <TopoDS_Shape.hxx>
#include <TopoDS_Compound.hxx>
#include <TColStd_ListOfAsciiString.hxx>
#include <STEPCAFControl_Reader.hxx>

//! Qt
#include <QList>
#include <QString>

//!C++
#include <stdio.h>
#include <windows.h>
#include <eh.h>
#include <string>
using namespace std;

class QOccProgressIndicator;

class STEPimporter: public QObject
{
    Q_OBJECT

public:

    STEPimporter(QObject *parent=0);

    //! import a STEP file
    bool import(const QString &fileName, TopoDS_Compound &Comp, QList<QString> &listOfNames);

    //! return the list of the names
    void STEP_GetEntityName(const TopoDS_Shape & theShape, STEPCAFControl_Reader* aReader, std::string &acName);

    void SEfunc(bool &isDone, opencascade::handle<TDocStd_Document> &step_doc);

    //! progress indicator
    void setProgressIndicator(const occHandle(QOccProgressIndicator) &aProgressIndicator);

private:

    STEPCAFControl_Reader reader;

    occHandle(QOccProgressIndicator) myProgress;
};

#endif // STEPIMPORTER_H
