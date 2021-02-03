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

    vecMatNames.push_back("Structural_steel");
    vecMatNames.push_back("Bilinear_steel");
    vecMatNames.push_back("H11_fatigue");
    vecMatNames.push_back("F22_fatigue");
    vecMatNames.push_back("B16_fatigue");
    vecMatNames.push_back("F6NM_fatigue");
    vecMatNames.push_back("F92_fatigue");
    vecMatNames.push_back("A479_fatigue");
    vecMatNames.push_back("SA479_XM19_fatigue");
    vecMatNames.push_back("SA182-B8M_CL2");
    vecMatNames.push_back("SA182-F316");
    vecMatNames.push_back("SA352-LCB");
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
bool ofwrite::perform() {return true;}

void ofwrite::clock()
{
    QTime time = QTime::currentTime();
    QString text = time.toString("mm:ss:zzz");
    cout<<text.toStdString()<<endl;
}

