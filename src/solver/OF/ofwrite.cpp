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
    int Nevents = 1;
    if(myProgressIndicator!=Q_NULLPTR)
    {
        //! ---------------------------------
        //! hide the additional progress bar
        //! ---------------------------------
        myProgressIndicator->setSecondaryBarVisible(false);

        QProgressEvent *e = new QProgressEvent(QProgressEvent_Init,0,Nevents,0,"Writing solver input file",
                                               QProgressEvent_None,-1,-1,-1,"Writing OF solver input file");
        QApplication::postEvent(myProgressIndicator,e);
        QApplication::processEvents();
        QThread::msleep(1000);
    }

    //! ----------------------------------
    //! a default name for the input dir
    //! ----------------------------------
    if(myFileDir=="") myFileDir ="OFinput.inp";

    QDir curDir;
    curDir.cd(myFileDir);
    curDir.mkdir("0");
    curDir.mkdir("system");
    curDir.mkdir("constant");

    //! number of items within the tree
    int N = mySimulationRoot->rowCount();

    //! read the "Analysis settings" item
    SimulationNodeClass *nodeAnalysisSettings = mySimulationRoot->child(0,0)->data(Qt::UserRole).value<SimulationNodeClass*>();
    CustomTableModel *tabData = nodeAnalysisSettings->getTabularDataModel();

    std::string path = myFileDir.toStdString() + "/0";
    myU.open((path+"/U").c_str());
    myP.open((path+"/p").c_str());
    myK.open((path+"/k").c_str());
    myW.open((path+"/omega").c_str());
    myNu.open((path+"/nut").c_str());

    std::vector<double> vect{0, 0, 0};

    of::printMark(myU);
    of::printHeading(myU,std::string("volVectorField"),std::string("U"),std::vector<int>());
    of::printBoundaryIntro(myU, 1, -1, vect);

    of::printMark(myP);
    of::printHeading(myP,std::string("volScalarField"),std::string("p"),std::vector<int>());
    of::printBoundaryIntro(myP, 2, -2, 0);

    of::printMark(myK);
    of::printHeading(myK,std::string("volScalarField"),std::string("k"),std::vector<int>());
    of::printBoundaryIntro(myK, 2, -2, 0.003);

    of::printMark(myW);
    of::printHeading(myW,std::string("volScalarField"),std::string("omega"),std::vector<int>());
    of::printBoundaryIntro(myW, 0, -1, 100000);

    of::printMark(myNu);
    of::printHeading(myNu,std::string("volScalarField"),std::string("nut"),std::vector<int>());
    of::printBoundaryIntro(myW, 2, -1, 0);

    std::map<std::string,std::map<int,occHandle(Ng_MeshVS_DataSourceFace)>> oFMap;
    for(int k=1; k<N-1; k++)
    {
        QString itemName = itemNameClearSpaces(mySimulationRoot->child(k,0)->data(Qt::DisplayRole).toString());
        cout<<"ofWriteClass::perform()->____found Item of type____"<<itemName.toStdString()<<"___"<<endl;

        SimulationNodeClass *theNode = mySimulationRoot->child(k,0)->data(Qt::UserRole).value<SimulationNodeClass*>();
        QStandardItem *theItem = static_cast<QStandardItem*>(mySimulationRoot->child(k,0));
        Property::SuppressionStatus theNodeSS = theNode->getPropertyValue<Property::SuppressionStatus>("Suppressed");
        SimulationNodeClass::nodeType theNodeType = theNode->getType();

        if(theNodeSS==Property::SuppressionStatus_Active)
        {
            QString setName = itemName.append("_").append(QString("%1").arg(k));
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

            QList<int> columnList;
            std::vector<double> vec;

            if (theNodeType != SimulationNodeClass::nodeType_CFDAnalysisBoundaryConditionWall)
            {
                columnList = mainTreeTools::getColumnsToRead(theItem,tabData->getColumnBeforeBC());
                double val;
                for (int i = 0; i < columnList.size(); i++)
                {
                    val = tabData->dataRC(1,columnList.at(i)).value<double>();
                    vec.push_back(val);
                }
            }

            std::string fullName;

            switch(theNodeType)
            {
            case SimulationNodeClass::nodeType_CFDAnalysisBoundaryConditionVelocity:
            {
                fullName = setName.prepend("velocity_").toStdString();

                of::printBoundary(myU,fullName,std::string("fixedValue"),vec);
                of::printBoundary(myP,fullName,std::string("zeroGradient"),std::vector<double>());
                of::printBoundary(myK,fullName,std::string("fixedValue"),std::string("$internalField"));
                of::printBoundary(myW,fullName,std::string("fixedValue"),std::string("$internalField"));
                of::printBoundary(myNu,fullName,std::string("calculated"),std::string("uniform 0"));
            }
                break;

            case SimulationNodeClass::nodeType_CFDAnalysisBoundaryConditionPressure:
            {
                fullName = setName.prepend("pressure_").toStdString();

                of::printBoundary(myU,fullName,std::string("zeroGradient"),std::vector<double>());
                of::printBoundary(myP,fullName,std::string("fixedValue"),vec);
                of::printBoundary(myK,fullName,std::string("zeroGradient"),std::vector<double>());
                of::printBoundary(myW,fullName,std::string("zeroGradient"),std::vector<double>());
                of::printBoundary(myNu,fullName,std::string("calculated"),std::string("uniform 0"));
            }
                break;

            case SimulationNodeClass::nodeType_CFDAnalysisBoundaryConditionWall:
            {
                fullName = setName.prepend("wall_").toStdString();

                of::printBoundary(myU,fullName,std::string("noSlip"),std::vector<double>());
                of::printBoundary(myP,fullName,std::string("zeroGradient"),std::vector<double>());
                of::printBoundary(myK,fullName,std::string("kqRWallFunction"),std::string("$internalField"));
                of::printBoundary(myW,fullName,std::string("omegaWallFunction"),std::string("$internalField"));
                of::printBoundary(myNu,fullName,std::string("nutkWallFunction"),std::string("uniform 0"));
            }
                break;
            }

            oFMap.insert(std::make_pair(fullName, tempMap));
        }
    }

    of::closeBoundary(myU);
    of::closeBoundary(myP);
    of::closeBoundary(myK);
    of::closeBoundary(myW);
    of::closeBoundary(myNu);

  /*  myU.close();
    myP.close();
    myK.close();
    myEPS.close();
    myW.close();
    myNUT.close();
    myNu.close();
*/

    if(myProgressIndicator!=Q_NULLPTR)
    {
        done++;
        QProgressEvent *e = new QProgressEvent(QProgressEvent_Update,0,Nevents-1,done,"Step 0 files generated",
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

    const occHandle(MeshVS_DataSource) &bodyMesh =  myDB->ArrayOfMeshDS.value(1); //change in vector/array of MeshVS_DS
    bool writeMesh = of::occToOF(bodyMesh, oFMap, myFileDir.toStdString());
    cout<<myFileDir.toStdString()<<endl;

    if(myProgressIndicator!=Q_NULLPTR)
    {
        done++;
        QProgressEvent *e = new QProgressEvent(QProgressEvent_Reset,-1,-1,-1,"",
                                               QProgressEvent_Reset,-1,-1,-1,"");
        QApplication::postEvent(myProgressIndicator,e);
        QApplication::processEvents();
        QThread::msleep(250);
    }

    if (writeMesh) return true;
    return false;
}

void ofwrite::clock()
{
    QTime time = QTime::currentTime();
    QString text = time.toString("mm:ss:zzz");
    cout<<text.toStdString()<<endl;
}

//void ofwrite::writeControlDict(std::ofstream of, std::string stype, double value)
//{

//}
