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
#include <cfdTool/occtoof.h>
#include "src/utils/cfdTool/occtoof.h"

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
    //myInputDir.setf(ios::scientific);
    //myInputDir.precision(EXPFORMAT_PRECISION);

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
                                               QProgressEvent_None,-1,-1,-1,"Writing OF input file");
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
                                               QProgressEvent_None,-1,-1,-1,"Writing OF solver input file");
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
    if(myFileDir=="") myFileDir ="D://test";  //qui ci va nome salvataggio

    QDir curDir;
    curDir.cd(myFileDir);
    curDir.mkdir("0");
    curDir.mkdir("system");
    curDir.mkdir("constant");

    //! open the file and set current directory
    //myInputFile.open(myFileName.toStdString());
    //QString inputName = myFileName;

    //! number of items within the tree
    int N = mySimulationRoot->rowCount();

    //! read the "Analysis settings" item
    SimulationNodeClass *nodeAnalysisSettings = mySimulationRoot->child(0,0)->data(Qt::UserRole).value<SimulationNodeClass*>();
    CustomTableModel *tabData = nodeAnalysisSettings->getTabularDataModel();

    std::string path = myFileDir.toStdString() + "/0";
    myU.open((path+"/U").c_str());
    myP.open((path+"/p").c_str());
    myK.open((path+"/k").c_str());
    myEPS.open((path+"/epsilon").c_str());
    myW.open((path+"/omega").c_str());
    myNUT.open((path+"/nuTilda").c_str());
    myNu.open((path+"/nut").c_str());
/*
    std::vector<ofstream> allBCfiles;
    allBCfiles.push_back(myU);
    allBCfiles.push_back(myP);
    allBCfiles.push_back(myK);
    allBCfiles.push_back(myEPS);
    allBCfiles.push_back(myW);
    allBCfiles.push_back(myNUT);
    allBCfiles.push_back(myNu);
*/
    of::printMark(myU);
    of::printHeading(myU,std::string("volVectorField"),std::string("U"),std::vector<int>());
    myU<<"dimensions      [0 1 -1 0 0 0 0];"<<endl;
    myU<<"internalField   uniform (0 0 0);"<<endl;
    myU<<"boundaryField"<<endl;
    myU<<"{"<<endl;

    of::printMark(myP);
    of::printHeading(myP,std::string("volScalarField"),std::string("p"),std::vector<int>());
    myP<<"dimensions      [0 2 -2 0 0 0 0];"<<endl;
    myP<<"internalField   uniform 0;"<<endl;
    myP<<"boundaryField"<<endl;
    myP<<"{"<<endl;

    of::printMark(myK);
    of::printHeading(myK,std::string("volScalarField"),std::string("k"),std::vector<int>());
    myK<<"dimensions      [0 2 -2 0 0 0 0];"<<endl;
    myK<<"internalField   uniform 1e-3;"<<endl;
    myK<<"boundaryField"<<endl;
    myK<<"{"<<endl;


    of::printHeading(myEPS,std::string("volScalarField"),std::string("epsilon"),std::vector<int>());
    myEPS<<"dimensions      [0 2 -3 0 0 0 0];"<<endl;
    myEPS<<"internalField   uniform 10;"<<endl;
    myEPS<<"boundaryField"<<endl;
    myEPS<<"{"<<endl;

    of::printHeading(myW,std::string("volScalarField"),std::string("omega"),std::vector<int>());
    myW<<"dimensions      [0 0 -1 0 0 0 0];"<<endl;
    myW<<"internalField   uniform 100000;"<<endl;
    myW<<"boundaryField"<<endl;
    myW<<"{"<<endl;

    of::printHeading(myNUT,std::string("volScalarField"),std::string("nuTilda"),std::vector<int>());
    myNUT<<"dimensions      [0 2 -1 0 0 0 0];"<<endl;
    myNUT<<"internalField   uniform 0;"<<endl;
    myNUT<<"boundaryField"<<endl;
    myNUT<<"{"<<endl;

    of::printHeading(myNu,std::string("volScalarField"),std::string("nut"),std::vector<int>());
    myNu<<"dimensions      [0 2 -1 0 0 0 0];"<<endl;
    myNu<<"internalField   uniform 0;"<<endl;
    myNu<<"boundaryField"<<endl;
    myNu<<"{"<<endl;

    std::map<std::string,std::map<int,occHandle(Ng_MeshVS_DataSourceFace)>> oFMap;
    for(int k=1; k<N-1; k++)
    {
        QString itemName = itemNameClearSpaces(mySimulationRoot->child(k,0)->data(Qt::DisplayRole).toString());
        cout<<"ofWriteClass::perform()->____found Item of type____"<<itemName.toStdString()<<"___"<<endl;

        SimulationNodeClass *theNode = mySimulationRoot->child(k,0)->data(Qt::UserRole).value<SimulationNodeClass*>();
        QStandardItem *theItem = mySimulationRoot->child(k,0)->data(Qt::UserRole).value<QStandardItem*>();

        Property::SuppressionStatus theNodeSS = theNode->getPropertyValue<Property::SuppressionStatus>("Suppressed");
        SimulationNodeClass::nodeType theNodeType = theNode->getType();

        if(theNodeSS==Property::SuppressionStatus_Active)
        {
            QString SetName = itemName.append("_").append(QString("%1").arg(k));
            IndexedMapOfMeshDataSources anIndexedMapOfFaceMeshDS;
            anIndexedMapOfFaceMeshDS = theNode->getPropertyValue<IndexedMapOfMeshDataSources>("Mesh data sources");
            std::map<int,occHandle(Ng_MeshVS_DataSourceFace)> tempMap;
            for (IndexedMapOfMeshDataSources::iterator iter = anIndexedMapOfFaceMeshDS.begin();
                 iter!=anIndexedMapOfFaceMeshDS.end(); ++iter)
            {
                int key = iter.key();
                const occHandle(MeshVS_DataSource) &valore = iter.value();
                const occHandle(Ng_MeshVS_DataSourceFace) &corrVal =
                        occHandle(Ng_MeshVS_DataSourceFace)::DownCast(valore);
                tempMap.insert(std::make_pair(key, corrVal));
            }

            //! possible BC value
            double p,U,k,w,eps,T,rho,v;

            QList<int> columnList = mainTreeTools::getColumnsToRead(theItem,tabData->getColumnBeforeBC());
            v = tabData->dataRC(1,columnList.at(0)).value<double>();
            std::vector<double> vec;
            vec.push_back(v);

            switch(theNodeType)
            {
            case SimulationNodeClass::nodeType_CFDAnalysisBoundaryConditionVelocity:
            {
                std::vector<double> p;
                SetName.prepend("inlet_");
                of::printBoundary(myU,std::string("fixedValue"),vec);
                of::printBoundary(myP,std::string("zeroGradient"),p);
                of::printBoundary(myK,std::string("fixedValue"),std::string("$internalField"));
                of::printBoundary(myEPS,std::string("fixedValue"),std::string("$internalField"));
                of::printBoundary(myW,std::string("fixedValue"),std::string("$internalField"));
                of::printBoundary(myNu,std::string("calculated"),std::string("$internalField"));
                of::printBoundary(myNUT,std::string("calculated"),std::string("$internalField"));
            }
                break;

            case SimulationNodeClass::nodeType_CFDAnalysisBoundaryConditionPressure:
            {
                SetName.prepend("pressure_");
            }
                break;

            case SimulationNodeClass::nodeType_CFDAnalysisBoundaryConditionWall:
            {
                SetName.prepend("wall_");
            }
                break;
            }

            oFMap.insert(std::make_pair(SetName.toStdString(), tempMap));



        }
    }


    const occHandle(MeshVS_DataSource) &bodyMesh =  myDB->ArrayOfMeshDS.value(1); //change in vector/array of MeshVS_DS
    bool writeMesh = of::occToOF(bodyMesh, oFMap, myFileDir.toStdString());


    return true;
}

void ofwrite::clock()
{
    QTime time = QTime::currentTime();
    QString text = time.toString("mm:ss:zzz");
    cout<<text.toStdString()<<endl;
}

void ofwrite::writeBC(std::ofstream of, std::string stype, double value)
{

}
