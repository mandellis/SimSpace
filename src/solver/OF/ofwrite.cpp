//! ----------------
//! custom includes
//! ----------------
#include "ofwrite.h"
#include "simulationdatabase.h"
#include "src/main/simulationmanager.h"
#include "src/utils/vectortool.h"
#include "src/main/mydefines.h"
#include "src/utils/geomtoolsclass.h"
#include "src/gui/tabularData/customtablemodel.h"
#include "src/utils/tools.h"
#include "ng_meshvs_datasourceface.h"
#include "qprogressevent.h"
#include "src/utils/ccout.h"
#include "postobject.h"
#include "occmeshtoccxmesh.h"
#include "src/gui/tabularData/tabulardatacolumns.h"
#include "src/connections/contactparameters.h"
#include "src/utils/geomtoolsclass.h"
#include "src/main/maintreetools.h"
#include <datasourcebuilder.h>
#include <src/utils/feaTool/bolttool.h>

//! -------
//! global
//! -------
#include "src/utils/global.h"

//! ---
//! Qt
//! ---
#include <QString>
#include <QTime>
#include <QDir>
#include <QApplication>
#include <QMessageBox>
#include <QThread>

//! ----
//! OCC
//! ----
#include <TColStd_PackedMapOfInteger.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopTools_ListIteratorOfListOfShape.hxx>
#include <TColStd_MapIteratorOfPackedMapOfInteger.hxx>
#include <MeshVS_DataSource.hxx>
#include <MeshVS_EntityType.hxx>
#include <TColgp_Array1OfPnt.hxx>
#include <GProp_PGProps.hxx>
#include <BRepGProp.hxx>
#include <gp_Trsf.hxx>
#include <gp_Dir.hxx>
#include <gp_TrsfForm.hxx>

//! ----------------------------
//! function: constructor I
//! details:  requires perform
//! ----------------------------
ofwrite::ofwrite(simulationDataBase *aDB, QExtendedStandardItem* aSimulationRoot, QObject *parent):
    QObject(parent)
{
    //cout<<"writeSolverFileClass::writeSolverFileClass()->____constructor I called____"<<endl;
    //! the simulation data-base
    myDB = aDB;

    //! init the progress indicator
    myProgressIndicator = Q_NULLPTR;

    //! the simulation item root
    mySimulationRoot = aSimulationRoot;

    //! set format
    myInputFile.setf(ios::scientific);
    myInputFile.precision(EXPFORMAT_PRECISION);

    vecMatNames.push_back("Water");
    vecMatNames.push_back("Air");
}


//! -------------------------------
//! function: setProgressIndicator
//! details:
//! -------------------------------
void ofwrite::setProgressIndicator(QProgressIndicator *aProgressIndicator)
{
    myProgressIndicator = aProgressIndicator;
}

//! ------------------
//! function: perform
//! details:
//! ------------------
bool ofwrite::perform()
{

    //! -----------------------------------
    //! reset the running status of global
    //! -----------------------------------
    Global::status().code = 1;

    //! ----------------------
    //! init the progress bar
    //! ----------------------
    int done = 0;
    int Nevents = 7;
    if(myProgressIndicator!=Q_NULLPTR)
    {
        //! ---------------------------------
        //! hide the additional progress bar
        //! ---------------------------------
        myProgressIndicator->setSecondaryBarVisible(false);

        QProgressEvent *e = new QProgressEvent(QProgressEvent_Init,0,Nevents,0,"Writing solver input file",
                                               QProgressEvent_None,-1,-1,-1,"Writing CCX input file");
        QApplication::postEvent(myProgressIndicator,e);
        QApplication::processEvents();
        QThread::msleep(1000);
    }

    //! --------------------
    //! update the progress
    //! --------------------
    if(myProgressIndicator!=Q_NULLPTR)
    {
        done++;
        QProgressEvent *e = new QProgressEvent(QProgressEvent_Update,0,Nevents-1,done,"Connectivity maps generated",
                                               QProgressEvent_None,-1,-1,-1,"Writing CCX solver input file");
        QApplication::postEvent(myProgressIndicator,e);
        QApplication::processEvents();
        QThread::msleep(250);

        if(Global::status().code==0)
        {
            cout<<"writeSolverFileClass::perform()->____process stopped____"<<endl;
            return false;
        }
    }

    //! ----------------------------------
    //! a default name for the input dir
    //! ----------------------------------
    if(myFileDir=="") myFileDir ="input.inp";  //qui ci va nome salvataggio

    //! open the file and set current directory
    //myInputFile.open(myFileName.toStdString());
    //QString inputName = myFileName;

    //! number of items within the tree
    int N = mySimulationRoot->rowCount();

    //! read the "Analysis settings" item
    SimulationNodeClass *nodeAnalysisSettings = mySimulationRoot->child(0,0)->data(Qt::UserRole).value<SimulationNodeClass*>();
    CustomTableModel *tabData = nodeAnalysisSettings->getTabularDataModel();

    return true;
}

void ofwrite::clock()
{
    QTime time = QTime::currentTime();
    QString text = time.toString("mm:ss:zzz");
    cout<<text.toStdString()<<endl;
}

