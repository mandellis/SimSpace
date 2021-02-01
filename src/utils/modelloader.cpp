//! custom includes
#include "modelloader.h"
#include "stepimporter.h"
#include <geometrydatabase.h>
#include "src/main/mydefines.h"

//! OCC
#include <TopoDS_Shape.hxx>
#include <TColStd_ListOfAsciiString.hxx>
#include <IGESControl_Reader.hxx>
#include <Interface_Static.hxx>

//! Qt
#include <QMessageBox>

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
ModelLoader::ModelLoader(const QString &fileName, QObject *parent):
    QObject(parent),
    myFileName(fileName)
{
    cout<<"ModelLoader::ModelLoader()->____constructor called____"<<endl;

    if(!myFileName.isEmpty() && !myFileName.isNull())
    {
        if(myFileName.endsWith(".stp")|myFileName.endsWith(".step")|myFileName.endsWith(".STEP"))
        {
            //! -------------------
            //! set the file type
            //! -------------------
            myFileType = fileType_STEP;

            /*
            //! ---------------------------
            //! configure the STEP reader
            //! read the old parameters
            //! ---------------------------
            Standard_Real ival = Interface_Static::IVal("read.precision.mode");
            Standard_Real rp = Interface_Static::RVal("read.precision.val");
            cout<<"ModelLoader::ModelLoader()->____mode: "<<ival<<"____"<<endl;
            cout<<"ModelLoader::ModelLoader()->____precision: "<<rp<<"____"<<endl;

            //! ----------------------------------
            //! set the new parameters and reread
            //! ----------------------------------
            Interface_Static::SetIVal("read.precision.mode",1);
            Interface_Static::SetRVal("read.precision.val",0.1);
            ival = Interface_Static::IVal("read.precision.mode");
            rp = Interface_Static::RVal("read.precision.val");
            cout<<"ModelLoader::ModelLoader()->____mode: "<<ival<<"____"<<endl;
            cout<<"ModelLoader::ModelLoader()->____precision: "<<rp<<"____"<<endl;
            */

            //! --------------------
            //! create the importer
            //! --------------------
            myStepImporter = new STEPimporter(this);
        }
    }
    else
    {
        //! set the NONE type
        myFileType = fileType_NONE;
    }
}

//! ------------------
//! function: perform
//! details:
//! ------------------
bool ModelLoader::perform(TopoDS_Compound &shapeFromReader,
                          QList<QString> &listOfNames,
                          const occHandle(QOccProgressIndicator) &aProgressIndicator)
{
    cout<<"ModelLoader::perform()->____function called____"<<endl;

    if(myFileType!=fileType_NONE)
    {
        bool isDone = false;
        switch(myFileType)
        {
        case fileType_STEP:
        {
            //! --------------------------------------------------------
            //! install the progress indicator
            //! the STEPimporter class is able to internally handle the
            //! case of a NULL QOccProgressIndicator handle
            //! --------------------------------------------------------
            cout<<"ModelLoader::perform()->____reading a .STEP file____"<<endl;
            if(!aProgressIndicator.IsNull()) myStepImporter->setProgressIndicator(aProgressIndicator);
            isDone = myStepImporter->import(myFileName, shapeFromReader, listOfNames);
        }
            break;

        default:
            //! ---------------------------
            //! file types not handled yet
            //! ---------------------------
            isDone = false;
            QMessageBox::warning(0,QString(APPNAME),QString("Cannot load geometry: format not supported"),QMessageBox::Ok);
            break;
        }
        emit resultReady(isDone);
        return isDone;
    }
    return false;
}
