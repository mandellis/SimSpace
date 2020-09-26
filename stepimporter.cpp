#define STEP_IMPORTER_CODE 1

//! custom includes
#include "stepimporter.h"
#include "se_exception.h"
#include "mydefines.h"
#include "qoccprogressindicator.h"
#include "steptools.h"

//! OCC
#include <XCAFApp_Application.hxx>
#include <XCAFDoc_ColorTool.hxx>
#include <XCAFDoc_ShapeTool.hxx>
#include <XCAFDoc_DocumentTool.hxx>
#include <TDocStd_Document.hxx>
#include <TDF_LabelSequence.hxx>
#include <STEPCAFControl_Reader.hxx>
#include <Quantity_Color.hxx>
#include <BRep_Builder.hxx>
#include <TopoDS_Compound.hxx>
#include <XSControl_WorkSession.hxx>
#include <XSControl_TransferReader.hxx>
#include <StepRepr_Representation.hxx>
#include <TCollection_HAsciiString.hxx>
#include <TColStd_ListOfAsciiString.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Shape.hxx>
#include <Transfer_TransientProcess.hxx>

#include <QApplication>

//! Qt
#include <QMessageBox>

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
STEPimporter::STEPimporter(QObject *parent):QObject(parent)
{
    ;
}

//! ----------------------
//! function: setProgress
//! details:
//! ----------------------
void STEPimporter::setProgressIndicator(const occHandle(QOccProgressIndicator) &aProgressIndicator)
{
    myProgress = aProgressIndicator;
}

//! -------------------------
//! function: trans_funct_N
//! details:
//! -------------------------
void trans_func(unsigned int u, EXCEPTION_POINTERS* pExp)
{
    Q_UNUSED(pExp)
    Q_UNUSED(u)
    printf( "In trans_func.\n" );
    throw SE_Exception();
}

//! ------------------------------------
//! function: SEfunct
//! details:  system exception function
//! ------------------------------------
void STEPimporter::SEfunc(bool &isDone, occHandle(TDocStd_Document) &step_doc)
{
    isDone = reader.Transfer(step_doc);
}

//! -----------------------
//! function: STEPimporter
//! details:
//! -----------------------
#include <Interface_Static.hxx>
bool STEPimporter::import(const QString &fileName, TopoDS_Compound &Comp, QList<QString> &listOfNames)
{
    cout<<"STEPimporter::import()->____importing file: "<<fileName.toStdString()<<"____"<<endl;

#if STEP_IMPORTER_CODE == 1
    //! Initiate an XCAF Application to handle the STEP XCAF Document
    static occHandle(XCAFApp_Application) dummy_app = XCAFApp_Application::GetApplication();

    //! Create an XCAF Document containing the STEP file
    occHandle(TDocStd_Document) step_doc;

    //! Check if a STEP File is already open the handle
    //! If so, close it to prevent Segmentation Faults when trying to create a new document
    if(dummy_app->NbDocuments()>0)
    {
       dummy_app->GetDocument(1,step_doc);
       dummy_app->Close(step_doc);
    }
    dummy_app->NewDocument("STEP-XCAF",step_doc);   // the "step_doc" is created here    

    //! ----------
    //! configure
    //! ----------
    reader.SetColorMode(Standard_True);
    reader.SetNameMode(Standard_True);
    reader.SetGDTMode(Standard_False);
    reader.SetLayerMode(Standard_False);
    reader.SetPropsMode(Standard_False);
    reader.SetSHUOMode(Standard_False);

    //! ----------------------
    //! reading the STEP file
    //! ----------------------
    if(myProgress)
    {
        if(myProgress->UserBreak()) return false;
        myProgress->NewScope(20, "Loading file");
        myProgress->Show();
    }

    Standard_Integer stat = reader.ReadFile(fileName.toStdString().c_str());
    if(myProgress) myProgress->EndScope();

    if(stat!=IFSelect_RetDone)
    {
       return false;
    }

    //! ------------------------------------------
    //! error handling for the "Transfer" process
    //! ------------------------------------------
    bool isDone;
    try
    {
        if(myProgress)
        {
            if(myProgress->UserBreak()) return false;
            myProgress->NewScope(60, "Transferring data");
            myProgress->Show();
        }

        _set_se_translator(trans_func);
        this->SEfunc(isDone, step_doc);

        if(myProgress) myProgress->EndScope();
    }
    catch(SE_Exception)
    {
        QMessageBox::critical(0,QString(APPNAME),QString("Unhandled exception during STEP model import"),QMessageBox::Ok,QMessageBox::NoButton);
        cerr<<"STEPimporter::import()->____unhandled error during STEP model import____"<<endl;

        if(myProgress) myProgress->EndScope();
        return false;
    }
    catch(...)
    {
        QMessageBox::critical(0,QString(APPNAME),QString("OCC error during STEP model import"),QMessageBox::Ok,QMessageBox::NoButton);
        cerr<<"STEPimporter::import()->____error in transfer function____"<<endl;

        if(myProgress) myProgress->EndScope();
        return false;
    }

    //! ----------------------
    //! "Transfer" process OK
    //! ----------------------
    if(isDone==true)
    {
        TopoDS_Builder builder;
        builder.MakeCompound(Comp);

        //! ----------------------------------------------------------------
        //! List out the available colours in the STEP File as Colour Names
        //! ----------------------------------------------------------------
        occHandle(XCAFDoc_ColorTool) step_colour_content = XCAFDoc_DocumentTool::ColorTool(step_doc->Main());
        TDF_LabelSequence all_colours;
        step_colour_content->GetColors(all_colours);
        for(int i=1; i<=all_colours.Length(); i++)
        {
            Quantity_Color col;
            step_colour_content->GetColor(all_colours.Value(i),col);
        }

        //! -----------------------------------------------
        //! List out the available shapes in the STEP File
        //! -----------------------------------------------
        occHandle(XCAFDoc_ShapeTool) step_shape_content = XCAFDoc_DocumentTool::ShapeTool(step_doc->Main());
        TDF_LabelSequence allShapes;

        step_shape_content->GetShapes(allShapes);
        //step_shape_content->GetFreeShapes(allShapes);
        cout<<"STEPimporter::import()->____Number of shapes: "<<allShapes.Length()<<"____"<<endl;
        for(int i=1; i<=allShapes.Length();i++) cout<<"____shape type: "<<allShapes.Value(i)<<"____"<<endl;

        for(int i=1; i<=allShapes.Length(); i++)
        {
            TopoDS_Shape shape1;
            std::string acName;

            step_shape_content->GetShape(allShapes.Value(i), shape1);
            switch(shape1.ShapeType())
            {
            case TopAbs_COMPOUND:
            {
                cout<<"STEPimporter::import()->____found TopoDS_COMPOUND____"<<endl;
            }
                break;

            case TopAbs_COMPSOLID:
            {
                cout<<"STEPimporter::import()->____found TopoDS_COMPSOLID____"<<endl;
                builder.Add(Comp,shape1);
                STEPimporter::STEP_GetEntityName(shape1,&reader,acName);
                if(acName=="") acName = "3D part";
                cout<<"STEPimporter::import()->____name: "<<acName<<"____"<<endl;
                listOfNames<<(QString::fromLatin1(acName.c_str()));
            }
                break;

            case TopAbs_SOLID:
            {
                cout<<"STEPimporter::import()->____found TopoDS_SOLID____"<<endl;
                builder.Add(Comp,shape1);
                STEPimporter::STEP_GetEntityName(shape1,&reader,acName);
                if(acName=="") acName = "3D body";
                cout<<"STEPimporter::import()->____name: "<<acName<<"____"<<endl;
                listOfNames<<(QString::fromLatin1(acName.c_str()));
            }
                break;

            case TopAbs_SHELL:
            {
                cout<<"STEPimporter::import()->____found TopoDS_SHELL____"<<endl;
                builder.Add(Comp,shape1);
                STEPimporter::STEP_GetEntityName(shape1,&reader,acName);
                if(acName=="") acName = "Surface body";
                cout<<"STEPimporter::import()->____name: "<<acName<<"____"<<endl;
                listOfNames<<(QString::fromLatin1(acName.c_str()));
            }
                break;

            case TopAbs_FACE:
            {
                cout<<"STEPimporter::import()->____found TopoDS_FACE____"<<endl;
                builder.Add(Comp,shape1);
                STEPimporter::STEP_GetEntityName(shape1,&reader,acName);
                if(acName=="") acName = "Surface body";
                cout<<"STEPimporter::import()->____name: "<<acName<<"____"<<endl;
                listOfNames<<(QString::fromLatin1(acName.c_str()));
            }
                break;

            case TopAbs_WIRE:
            {
                cout<<"STEPimporter::import()->____found TopoDS_WIRE____"<<endl;
                builder.Add(Comp,shape1);
                STEPimporter::STEP_GetEntityName(shape1,&reader,acName);
                if(acName=="") acName = "Line body";
                cout<<"STEPimporter::import()->____name: "<<acName<<"____"<<endl;
                listOfNames<<(QString::fromLatin1(acName.c_str()));
            }
                break;

            case TopAbs_EDGE:
            {
                cout<<"STEPimporter::import()->____found TopoDS_EDGE____"<<endl;
                builder.Add(Comp,shape1);
                STEPimporter::STEP_GetEntityName(shape1,&reader,acName);
                if(acName=="") acName = "Line body";
                cout<<"STEPimporter::import()->____name: "<<acName<<"____"<<endl;
                listOfNames<<(QString::fromLatin1(acName.c_str()));
            }
                break;
            }
        }
        return true;
    }
    else
    {
        return false;
    }
#elif STEP_IMPORTER_CODE ==0
    //! ----------------------------------------------
    //! a minimal version of the step importer
    //! this will not activate the progress indicator
    //! ----------------------------------------------
    STEPControl_Reader aReader;

    IFSelect_ReturnStatus stat = aReader.ReadFile(fileName.toStdString().c_str());
    if(stat!=IFSelect_RetDone)
    {
        cerr<<"STEPimporter::import()->____the STEP file cannot be imported____"<<endl;
        return false;
    }

    Interface_Static::SetIVal("read.step.assembly.level",1);
    aReader.TransferRoots();
    TopoDS_Builder compBuilder;
    compBuilder.MakeCompound(Comp);
    TopoDS_Shape shapeFromReader = aReader.OneShape();
    for(TopExp_Explorer anExp(shapeFromReader,TopAbs_SOLID); anExp.More(); anExp.Next())
    {
        TopoDS_Solid aSolid = TopoDS::Solid(anExp.Current());
        if(aSolid.IsNull()) continue;
        compBuilder.Add(Comp,aSolid);
    }

    TColStd_ListOfAsciiString names;
    STEPTools::dumpSTEPLabels(aReader,names);
    for(TColStd_ListIteratorOfListOfAsciiString anIt(names); anIt.More(); anIt.Next())
    {
        TCollection_AsciiString cname = anIt.Value();
        cout<<"STEPimporter::import()->____name: "<<cname.ToCString()<<"____"<<endl;
        QString bodyName = QString::fromLatin1(cname.ToCString());
        listOfNames<<bodyName;
    }
    return true;
#endif
}

//! -----------------------------
//! function: STEP_GetEntityName
//! details:
//! -----------------------------
void STEPimporter::STEP_GetEntityName(const TopoDS_Shape &theShape, STEPCAFControl_Reader *aReader, std::string &acName)
{
    const occHandle(XSControl_WorkSession)& theSession = aReader->Reader().WS();
    const occHandle(XSControl_TransferReader)& aTransferReader = theSession->TransferReader();

    occHandle(Standard_Transient) anEntity = aTransferReader->EntityFromShapeResult(theShape, 1);

    if(anEntity.IsNull())
    {
        // as just mapped
        anEntity = aTransferReader->EntityFromShapeResult (theShape,-1);
    }
    if(anEntity.IsNull())
    {
        // as anything
        anEntity = aTransferReader->EntityFromShapeResult (theShape,4);
    }
    if(anEntity.IsNull())
    {
        cout<<"Warning: XSInterVertex_STEPReader::ReadAttributes() entity not found"<<endl;
    }
    else
    {
        occHandle(StepRepr_RepresentationItem) aReprItem;
        aReprItem = occHandle(StepRepr_RepresentationItem)::DownCast(anEntity);

        if (aReprItem.IsNull())
        {
            cout<<"Error: STEPReader::ReadAttributes(): StepRepr_RepresentationItem Is NULL"<<endl;
        }
        else acName = aReprItem->Name()->ToCString();
    }
}

bool STEPimporter::import1(const QString &fileName, TopoDS_Compound &Comp, QList<QString> &listOfNames)
{
    BRep_Builder builder;
    builder.MakeCompound(Comp);

    STEPControl_Reader aReader;
    TopoDS_Shape aShape;
    if(fileName.isNull()) return false;
    if(fileName.isEmpty()) return false;

    Standard_CString name = fileName.toStdString().c_str();
    if (aReader.ReadFile((Standard_CString)name)!=IFSelect_RetDone) return false;

    //occHandle(Message_ProgressIndicator) pi = new ProgressIndicator(100);
    //aReader.WS()->MapReader()->SetProgress(pi);
    //pi->NewScope(100, "Reading STEP file...");
    //pi->Show();

    // Root transfers
    int nbr = aReader.NbRootsForTransfer();
    //aReader.PrintCheckTransfer (failsonly, IFSelect_ItemsByEntity);
    for (int n = 1; n<= nbr; n++)
    {
        cout<<"____transferring root "<<n<<"____"<<endl;
        aReader.TransferRoot(n);
    }
    //pi->EndScope();

    // Collecting resulting entities
    int nbs = aReader.NbShapes();
    if (nbs == 0)
    {
        cout<<"____no shapes in file____"<<endl;
        return false;
    }
    else
    {
        //Handle(StepData_StepModel) Model = aReader.StepModel();
        //Handle(XSControl_WorkSession) ws = aReader.WS();
        //Handle(XSControl_TransferReader) tr = ws->TransferReader();

        //std::map<int, Quantity_Color> hash_col;
        //ReadColors(aReader.WS(), hash_col);
        //ReadNames(aReader.WS());

        for (int i=1; i<=nbs; i++)
        {
            cout<<"____transferring shape: "<<i<<"____"<<endl;
            aShape = aReader.Shape(i);

            // load each solid as an own object
            TopExp_Explorer ex;
            for (ex.Init(aShape, TopAbs_SOLID); ex.More(); ex.Next())
            {
                // get the shape
                const TopoDS_Solid& aSolid = TopoDS::Solid(ex.Current());
            }

            // load all non-solids now
            for (ex.Init(aShape, TopAbs_SHELL, TopAbs_SOLID); ex.More(); ex.Next())
            {
                // get the shape
                const TopoDS_Shell& aShell = TopoDS::Shell(ex.Current());
            }

            // put all other free-flying shapes into a single compound

            for (ex.Init(aShape, TopAbs_FACE, TopAbs_SHELL); ex.More(); ex.Next()) {
                if (!ex.Current().IsNull()) {
                    builder.Add(Comp, ex.Current());
                }
            }
            for (ex.Init(aShape, TopAbs_WIRE, TopAbs_FACE); ex.More(); ex.Next()) {
                if (!ex.Current().IsNull()) {
                    builder.Add(Comp, ex.Current());
                }
            }
            for (ex.Init(aShape, TopAbs_EDGE, TopAbs_WIRE); ex.More(); ex.Next()) {
                if (!ex.Current().IsNull()) {
                    builder.Add(Comp, ex.Current());
                }
            }
            for (ex.Init(aShape, TopAbs_VERTEX, TopAbs_EDGE); ex.More(); ex.Next()) {
                if (!ex.Current().IsNull()) {
                    builder.Add(Comp, ex.Current());
                }
            }

            //if (!emptyComp) {
            //    std::string name = fi.fileNamePure();
            //    Part::Feature *pcFeature = static_cast<Part::Feature*>(pcDoc->addObject
            //        ("Part::Feature", name.c_str()));
            //    pcFeature->Shape.setValue(comp);
            //}
        }
    }
    return true;
}
