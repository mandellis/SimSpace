//! ----------------
//! custom includes
//! ----------------
#include "writesolverfileclass.h"
#include "simulationdatabase.h"
#include "simulationmanager.h"
#include "vectortool.h"
#include "mydefines.h"
#include "geomtoolsclass.h"
#include "customtablemodel.h"
#include "tools.h"
#include "ng_meshvs_datasourceface.h"
#include "qprogressevent.h"
#include "ccout.h"
#include "postobject.h"
#include "occmeshtoccxmesh.h"
#include "tabulardatacolumns.h"
#include "contactparameters.h"
#include "geomtoolsclass.h"
#include "maintreetools.h"
#include <facedatasourcebuilder.h>
#include <bolttool.h>

//! -------
//! global
//! -------
#include "global.h"
#include "map"

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
writeSolverFileClass::writeSolverFileClass(simulationDataBase *aDB, QExtendedStandardItem* aSimulationRoot, QObject *parent):
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
void writeSolverFileClass::setProgressIndicator(QProgressIndicator *aProgressIndicator)
{
    myProgressIndicator = aProgressIndicator;
}

//! ------------------
//! function: perform
//! details:
//! ------------------
bool writeSolverFileClass::perform()
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

    //! ----------------------------
    //! build the connectivity maps
    //! ----------------------------
    for(QMap<int,TopoDS_Shape>::iterator it = myDB->bodyMap.begin(); it!=myDB->bodyMap.end(); it++)
    {
        if(Global::status().code==0)
        {
            cout<<"writeSolverFileClass::perform()->____process stopped____"<<endl;
            return false;
        }

        int bodyIndex = it.key();
        occHandle(Ng_MeshVS_DataSource3D) curVolumeMesh = occHandle(Ng_MeshVS_DataSource3D)::DownCast(myDB->ArrayOfMeshDS.value(bodyIndex));

        std::map<meshElement2D,std::vector<std::pair<int,int>>> facesToElements;
        curVolumeMesh->buildCCXFaceToElementConnectivity(facesToElements);

        //! ------------------------------------------------------------
        //! create the face to elements connectivity map, for each body
        //! ------------------------------------------------------------
        std::pair<int, std::map<meshElement2D,std::vector<std::pair<int,int>>>> p;
        p.first = bodyIndex;
        p.second = facesToElements;
        bigMap.insert(p);
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
        QThread::msleep(500);

        if(Global::status().code==0)
        {
            cout<<"writeSolverFileClass::perform()->____process stopped____"<<endl;
            return false;
        }
    }

    //! ----------------------------------
    //! a default name for the input file
    //! ----------------------------------
    if(myFileName=="") myFileName ="input.inp";

    //! open the file and set current directory
    myInputFile.open(myFileName.toStdString());
    QString inputName = myFileName;

    //! number of items within the tree
    int N = mySimulationRoot->rowCount();

    //! read the "Analysis settings" item
    SimulationNodeClass *nodeAnalysisSettings = mySimulationRoot->child(0,0)->data(Qt::UserRole).value<SimulationNodeClass*>();
    CustomTableModel *tabData = nodeAnalysisSettings->getTabularDataModel();

    //! -------------------------
    //! write nodes and elements
    //! -------------------------
    QMap<int,QList<int>> nodeListByBody;
    this->writeNodesAndElements(inputName,nodeListByBody);

    //! --------------------
    //! update the progress
    //! --------------------
    if(myProgressIndicator!=Q_NULLPTR)
    {
        done++;
        QProgressEvent *e = new QProgressEvent(QProgressEvent_Update,0,Nevents-1,done,"Sending nodes and elements",
                                               QProgressEvent_None,-1,-1,-1,"Writing CCX solver input file");
        QApplication::postEvent(myProgressIndicator,e);
        QApplication::processEvents();
        QThread::msleep(500);

        if(Global::status().code==0)
        {
            cout<<"writeSolverFileClass::perform()->____process stopped____"<<endl;
            return false;
        }
    }

    //! retrieve the type of simulation => unused for the moment <=
    SimulationNodeClass::nodeType theSimulationType = mySimulationRoot->data(Qt::UserRole).value<SimulationNodeClass*>()->getType();
    Q_UNUSED(theSimulationType);

    //! -------------------------
    //! scan the simulation tree
    //! .  contacts root
    //! .  named selections root
    //! .  simulation root
    //! -------------------------

    //! --------------------------------------------
    //! [1] write element/node sets: read the setup
    //! --------------------------------------------
    for(int k=1; k<N-1; k++)
    {
        if(Global::status().code==0)
        {
            cout<<"writeSolverFileClass::perform()->____process stopped____"<<endl;
            return false;
        }

        QString itemName = itemNameClearSpaces(mySimulationRoot->child(k,0)->data(Qt::DisplayRole).toString());
        cout<<"writeSolverFileClass::perform()->____found Item of type____"<<itemName.toStdString()<<"___"<<endl;

        SimulationNodeClass *theItemNode = mySimulationRoot->child(k,0)->data(Qt::UserRole).value<SimulationNodeClass*>();
        SimulationNodeClass::nodeType theNodeType = theItemNode->getType();
        if(theNodeType==SimulationNodeClass::nodeType_mapper
                || theNodeType==SimulationNodeClass::nodeType_modelChange
                || theNodeType==SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_ImportedTemperatureDistribution
                || theNodeType==SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Acceleration
        #ifdef COSTAMP_VERSION
                || theNodeType==SimulationNodeClass::nodeType_timeStepBuilder
        #endif
                )
            continue;
        Property::SuppressionStatus theNodeSS = theItemNode->getPropertyValue<Property::SuppressionStatus>("Suppressed");

        if(theNodeSS==Property::SuppressionStatus_Active)
        {
            //! ---------------
            //! find the scope
            //! ---------------
            //! append to the item name the number of the row
            QString SetName = itemName.append("_").append(QString("%1").arg(k));

            //! ------------------------------------------------------------------------------
            //! retrive the IndexedMapOfMeshDS of the BC - they do not exist in case of bolts
            //! ------------------------------------------------------------------------------
            IndexedMapOfMeshDataSources anIndexedMapOfFaceMeshDS;
            if(theNodeType!=SimulationNodeClass::nodeType_structuralAnalysisBoltPretension)
                anIndexedMapOfFaceMeshDS = theItemNode->getPropertyValue<IndexedMapOfMeshDataSources>("Mesh data sources");
            switch(theNodeType)
            {
            case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_CompressionOnlySupport:
            {
                double K = theItemNode->getPropertyValue<double>("K");
                double F = theItemNode->getPropertyValue<double>("Sigma infinity");

                this->writeGapElement(anIndexedMapOfFaceMeshDS,SetName,K,F);
                myInputFile<<"*INCLUDE, INPUT="<<SetName.toStdString()<<".gap"<<endl;
            }
                break;
            case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Acceleration:
            {
                cout<<"Acceleration"<<endl
                ;
            }
                break;
            case SimulationNodeClass::nodeType_structuralAnalysisBoltPretension:
            {
                QString refNodeName = itemName+"_RN";
                //! ---------------------------------
                //! retrieve the axis of the bolt
                //! and the reference point location
                //! ---------------------------------
                void *cs = theItemNode->getPropertyValue<void*>("Coordinate system");
                QStandardItem *itemCS = static_cast<QStandardItem*>(cs);
                SimulationNodeClass *nodeCS = itemCS->data(Qt::UserRole).value<SimulationNodeClass*>();

                QVector<double> origin;     //! this is used for the reference node
                QVector<double> boltAxis;

                if(nodeCS->getType()==SimulationNodeClass::nodeType_coordinateSystem_global)
                {
                    for(int i=0; i<3; i++) origin.push_back(0);
                    boltAxis = nodeCS->getPropertyValue<QVector<double>>("Z axis data");
                }
                else
                {
                    QVector<double> baseOrigin = nodeCS->getPropertyValue<QVector<double>>("Base origin");
                    for(int i=0; i<3; i++) origin.push_back(baseOrigin[i]);
                    QVector<QVector<double>> baseDirectionalData = nodeCS->getPropertyValue<QVector<QVector<double>>>("Base directional data");
                    boltAxis = baseDirectionalData.at(2);
                }
                myInputFile<<"*NODE, NSET = "<<refNodeName.toStdString()<<endl;
                myInputFile<<++totalNumberOfNodes<<", "<<origin[0]<<", "<<origin[1]<<", "<<origin[2]<<endl;

                //! ----------
                //! eccezione
                //! ----------
                std::vector<GeometryTag> tags = theItemNode->getPropertyValue<std::vector<GeometryTag>>("Tags");
                int boltBodyIndex = tags.at(0).parentShapeNr;
                double d = -boltAxis[0]*origin[0]-boltAxis[1]*origin[1]-boltAxis[2]*origin[2];
                this->writeElementSurfaceBolt(SetName,boltBodyIndex,boltAxis[0],boltAxis[1],boltAxis[2],d);
                myInputFile<<"*PRE-TENSION SECTION, SURFACE = "<<SetName.toStdString()<<", NODE ="<<totalNumberOfNodes<<endl;
                myInputFile<<boltAxis[0]<<", "<<boltAxis[1]<<","<<boltAxis[2]<<endl;
            }
                break;
            case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Force:
            case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Moment:
            case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteForce:
            case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteDisplacement:
            {
                this->writeElementSurface(SetName,anIndexedMapOfFaceMeshDS);
                //! ------------------------------------------------
                //! retrieve the coordinates of the reference point
                //! ------------------------------------------------
                QVector<double> refPoint = theItemNode->getPropertyValue<QVector<double>>("Reference point");
                double rPx,rPy,rPz;
                rPx = refPoint.at(0);
                rPy = refPoint.at(1);
                rPz = refPoint.at(2);

#ifdef COSTAMP_VERSION
                Property::defineBy theDefineBy = theItemNode->getPropertyValue<Property::defineBy>("Define by");
                if(theDefineBy == Property::defineBy_vector)
                {
                    QVector<double> dirData = theItemNode->getPropertyValue<QVector<double>>("Direction");
                    double a0 = dirData.at(3);
                    double a1 = dirData.at(4);
                    double a2 = dirData.at(5);
                    int h=500; //to handle
                    rPx = refPoint.at(0)+h*a0;
                    rPy = refPoint.at(1)+h*a1;
                    rPz = refPoint.at(2)+h*a2;
                }
#endif

                //! ------------------------------
                //! retrieve the type of coupling
                //! ------------------------------
                int TOC = theItemNode->getPropertyValue<int>("Coupling");
                bool isXcoupled;

                myInputFile<<"*NODE,NSET="<<(QString("CM_")+SetName).toStdString()<<endl;
                myInputFile<<++totalNumberOfNodes<<","<<rPx<<","<<rPy<<","<<rPz<<endl;
                myInputFile<<"*COUPLING,REF NODE="<<totalNumberOfNodes<<",SURFACE="<<SetName.toStdString()<<",CONSTRAINT NAME="<<SetName.toStdString()<<endl;

                switch(TOC)
                {
                case 0:
                {
                    //! -----------------
                    //! case "Kinematic"
                    //! -----------------
                    myInputFile<<"*KINEMATIC"<<endl;
                    //! -------------------------
                    //! retrieve the DOFs status
                    //! -------------------------
                    bool isXcoupled = theItemNode->getPropertyValue<int>("X component ");
                    bool isYcoupled = theItemNode->getPropertyValue<int>("Y component ");
                    bool isZcoupled = theItemNode->getPropertyValue<int>("Z component ");
                    if(isXcoupled) myInputFile<<"1,1"<<endl;
                    if(isYcoupled) myInputFile<<"2,2"<<endl;
                    if(isZcoupled) myInputFile<<"3,3"<<endl;
                }
                    break;
                case 1:
                {
                    //! -------------------
                    //! case "Distributed"
                    //! -------------------
                    myInputFile<<"*DISTRIBUTING"<<endl;
                    //! -------------------------
                    //! retrieve the DOFs status
                    //! -------------------------
                    isXcoupled = theItemNode->getPropertyValue<int>("X component ");
                    bool isYcoupled = theItemNode->getPropertyValue<int>("Y component ");
                    bool isZcoupled = theItemNode->getPropertyValue<int>("Z component ");
                    bool isROTXcoupled = theItemNode->getPropertyValue<int>("X rotation");
                    bool isROTYcoupled = theItemNode->getPropertyValue<int>("Y rotation");
                    bool isROTZcoupled = theItemNode->getPropertyValue<int>("Z rotation");
                    if(isXcoupled) myInputFile<<"1,1"<<endl;
                    if(isYcoupled) myInputFile<<"2,2"<<endl;
                    if(isZcoupled) myInputFile<<"3,3"<<endl;
                    if(isROTXcoupled) myInputFile<<"4,4"<<endl;
                    if(isROTYcoupled) myInputFile<<"5,5"<<endl;
                    if(isROTZcoupled) myInputFile<<"6,6"<<endl;
                }
                    break;
                }
            }
                break;
            case SimulationNodeClass::nodeType_structuralAnalysisBoundaryContidion_FixedSupport:
            case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Displacement:
            case SimulationNodeClass::nodeType_structuralAnalysisThermalCondition:
            case SimulationNodeClass::nodeType_thermalAnalysisTemperature:
            case SimulationNodeClass::nodeType_thermalAnalysisThermalPower:
            {
                this->writeNodalSet(SetName,anIndexedMapOfFaceMeshDS);
            }
                break;
            case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_FrictionlessSupport:
            case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_CylindricalSupport:
            {
                QString extension=".nam";
                QString nsetName=SetName+extension;

                QString absFileName = myFileName.split("/").last();
                QString dirName = myFileName;
                dirName.chop(absFileName.size());
                nsetName.prepend(dirName);
                ofstream myNSet;
                myNSet.setf(ios::scientific);
                myNSet.precision(EXPFORMAT_PRECISION);
                myNSet.open(nsetName.toStdString());
                for(QMap<int,opencascade::handle<MeshVS_DataSource>>::const_iterator it = anIndexedMapOfFaceMeshDS.cbegin(); it!= anIndexedMapOfFaceMeshDS.cend(); ++it)
                {
                    occHandle(Ng_MeshVS_DataSourceFace) &faceMeshDS = occHandle(Ng_MeshVS_DataSourceFace)::DownCast(it.value());
                    faceMeshDS->computeNormalAtNodes();
                    QMap<int,QList<double>> nodeNormals = faceMeshDS->myNodeNormals;

                    int bodyIndex = it.key();
                    int offset = 0;
                    for(int k=1; k<bodyIndex; k++)
                    {
                        if(myDB->MapOfIsActive.value(k)==true)
                            offset = offset+myDB->ArrayOfMeshDS.value(k)->GetAllNodes().Extent();
                    }
                    double buf[3];
                    TColStd_Array1OfReal aCoords(*buf,1,3);
                    Standard_Integer nbNodes;
                    MeshVS_EntityType aType;

                    for(TColStd_MapIteratorOfPackedMapOfInteger anIter(faceMeshDS->GetAllNodes()); anIter.More(); anIter.Next())
                    {
                        Standard_Integer globalNodeID = anIter.Key();
                        QList<double> nodeNormal = nodeNormals.value(globalNodeID);

                        if(!faceMeshDS->GetGeom(globalNodeID,Standard_False,aCoords,nbNodes,aType)) continue;
                        Standard_Real x = aCoords(1);
                        Standard_Real y = aCoords(2);
                        Standard_Real z = aCoords(3);

                        //! tangent plane equation ax+by+cz+d=0   z=-(d+ax+by)/c
                        double a = nodeNormal.at(0);
                        double b = nodeNormal.at(1);
                        double c = nodeNormal.at(2);
                        double d = -(a*x+b*y+c*z);
                        double delta = 1.0;
                        //delta=1;

                        double p1_x,p1_y,p1_z,p2_x,p2_y,p2_z;

                        if(a!=0.0)
                        {
                            p1_z = z+delta;
                            p1_y = y+delta;
                            p1_x = -(d+c*p1_z+b*p1_y)/a;

                            p2_z = z;
                            p2_y = y;
                            p2_x = -(d+c*p2_z+b*p2_y)/a;
                        }
                        else if(b!=0.0)
                        {
                            p1_x = x+delta;
                            p1_z = z+delta;
                            p1_y = -(d+c*p1_z+a*p1_x)/b;

                            p2_x = x;
                            p2_z = z;
                            p2_y = -(d+c*p2_z+a*p2_x)/b;
                        }
                        else if(c!=0)
                        {
                            p1_x = x+delta;
                            p1_y = y+delta;
                            p1_z = -(d+a*p1_x+b*p1_y)/c;

                            p2_x = x;
                            p2_y = y;
                            p2_z = -(d+a*p2_x+b*p2_y)/c;
                        }

                        QString name = SetName+QString("_")+(QString("%1").arg(globalNodeID+offset));
                        myNSet<<"*NSET, NSET = "<<name.toStdString()<<endl;
                        myNSet<<globalNodeID+offset<<endl;
                        myNSet<<"*TRANSFORM, NSET = "<<name.toStdString()<<" , TYPE = R" <<endl;
                        myNSet<<p1_x<<","<<p1_y<<","<<p1_z<<","<<p2_x<<","<<p2_y<<","<<p2_z<<endl;
                        myNSet<<"*BOUNDARY"<<endl;
                        myNSet<<name.toStdString()<<","<<"3,3"<<endl;
                    }
                }
                myNSet<<endl;
                myNSet.close();
                myInputFile<<"*INCLUDE, INPUT="<<nsetName.split("/").last().toStdString()<<endl;
            }
                break;
            case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Pressure:
            case SimulationNodeClass::nodeType_thermalAnalysisConvection:
            case SimulationNodeClass::nodeType_thermalAnalysisAdiabaticWall:
            case SimulationNodeClass::nodeType_thermalAnalysisThermalFlow:
            case SimulationNodeClass::nodeType_thermalAnalysisThermalFlux:
            {
                this->writeElementSurface(SetName,anIndexedMapOfFaceMeshDS);
            }
                break;
            default:
                break;
            }
        }
    }

    //! -----------------
    //! write point mass
    //! -----------------
    QStandardItem *theGeometryRoot=this->getTreeItem(SimulationNodeClass::nodeType_geometry);
    cout<<"writeSolverFileClass::perform()->____writing Point Mass___"<<endl;
    for(int k=0; k<theGeometryRoot->rowCount();k++)
    {
        if(Global::status().code==0)
        {
            cout<<"writeSolverFileClass::perform()->____process stopped____"<<endl;
            return false;
        }

        QStandardItem *theGeometryItem = theGeometryRoot->child(k,0);
        SimulationNodeClass *theCurNode = theGeometryItem->data(Qt::UserRole).value<SimulationNodeClass*>();
        Property::SuppressionStatus theNodeSS = theCurNode->getPropertyValue<Property::SuppressionStatus>("Suppressed");
        if(theNodeSS==Property::SuppressionStatus_Active &&
                theCurNode->getType()==SimulationNodeClass::nodeType_pointMass)
        {
            QString itemName = itemNameClearSpaces(theGeometryRoot->child(k,0)->data(Qt::DisplayRole).toString());
            itemName.append("_").append(QString("%1").arg(k));
            //! retrive the IndexedMapOfMeshDS of the BC
            IndexedMapOfMeshDataSources anIndexedMapOfFaceMeshDS;
            anIndexedMapOfFaceMeshDS = theCurNode->getPropertyValue<IndexedMapOfMeshDataSources>("Mesh data sources");

            this->writeElementSurface(itemName,anIndexedMapOfFaceMeshDS);    //TO DO DS is missing for point mass

            //std::vector<GeometryTag> scope = theCurNode->getPropertyValue<std::vector<GeometryTag>>("Tags");
            double mass = theCurNode->getPropertyValue<double>("Mass");
            double Jx = theCurNode->getPropertyValue<double>("Jx");
            double Jy = theCurNode->getPropertyValue<double>("Jy");
            double Jz = theCurNode->getPropertyValue<double>("Jz");

            Q_UNUSED (Jx)
            Q_UNUSED (Jy)
            Q_UNUSED (Jz)

            double x = theCurNode->getPropertyValue<double>("X coordinate");
            double y = theCurNode->getPropertyValue<double>("Y coordinate");
            double z = theCurNode->getPropertyValue<double>("Z coordinate");

            myInputFile<<"*NODE,NSET="<<(QString("CM_")+itemName).toStdString()<<endl;
            myInputFile<<++totalNumberOfNodes<<","<<x<<","<<y<<","<<z<<endl;
            myInputFile<<"*ELEMENT, ELSET ="<<(QString("E")+itemName).toStdString()<<", TYPE=MASS"<<endl;
            myInputFile<<++totalNumberOfElements<<","<<totalNumberOfNodes<<endl;
            myInputFile<<"*MASS, ELSET ="<<(QString("E"
                                                    "")+itemName).toStdString()<<endl;
            myInputFile<<mass<<endl;
            myInputFile<<"*COUPLING,REF NODE="<<totalNumberOfNodes<<",SURFACE="<<itemName.toStdString()<<",CONSTRAINT NAME="<<itemName.toStdString()<<endl;
            myInputFile<<"*KINEMATIC"<<endl;
            myInputFile<<"1,1"<<endl;
            myInputFile<<"2,2"<<endl;
            myInputFile<<"3,3"<<endl;
        }
    }

    //! --------------------
    //! update the progress
    //! --------------------
    if(myProgressIndicator!=Q_NULLPTR)
    {
        done++;
        QProgressEvent *e = new QProgressEvent(QProgressEvent_Update,0,Nevents-1,done,"Sending boundary conditions",
                                               QProgressEvent_None,-1,-1,-1,"Writing CCX solver input file");
        QApplication::postEvent(myProgressIndicator,e);
        QApplication::processEvents();
        QThread::msleep(500);

        if(Global::status().code==0)
        {
            cout<<"writeSolverFileClass::perform()->____process stopped____"<<endl;
            return false;
        }
    }

    //! -------------------------------
    //! [2] read the connections group
    //! -------------------------------
    //! map for storage of master and slave name
    QMap<QString,pair<QString,QString>> contactMapName;

    QExtendedStandardItem *theConnectionItem = this->getTreeItem(SimulationNodeClass::nodeType_connection);
    for(int n=0; n<theConnectionItem->rowCount(); n++)
    {
        if(Global::status().code==0)
        {
            cout<<"writeSolverFileClass::perform()->____process stopped____"<<endl;
            return false;
        }

        //! the current connection group
        QStandardItem *itemConnectionGroup = theConnectionItem->child(n,0);

        //! number of contacts under the current connection group
        int NbContactPairs = itemConnectionGroup->rowCount();

        for(int k=0; k<NbContactPairs; k++)
        {
            QStandardItem *item = itemConnectionGroup->child(k,0);
            SimulationNodeClass *node = item->data(Qt::UserRole).value<SimulationNodeClass*>();

            Property::SuppressionStatus ss = node->getPropertyValue<Property::SuppressionStatus>("Suppressed");
            if(ss==Property::SuppressionStatus_Active)
            {
                QString timeTag = node->getPropertyValue<QString>("Time tag");
                std::pair<QString ,QString> masterSlaveName;

                IndexedMapOfMeshDataSources anIndexedMapOfFaceMeshDS_Slave = node->getPropertyValue<IndexedMapOfMeshDataSources>("Slave mesh data source");
                IndexedMapOfMeshDataSources anIndexedMapOfFaceMeshDS_Master = node->getPropertyValue<IndexedMapOfMeshDataSources>("Master mesh data source");

                Property::contactBehavior theContactBehavior = node->getPropertyValue<Property::contactBehavior>("Behavior");
                switch(theContactBehavior)
                {
                case Property::contactBehavior_asymmetric:
                {
                    //! -------------------------------------------------------
                    //! find the scope of the slave, the write a NODAL surface
                    //! -------------------------------------------------------
                    QString name = itemNameClearSpaces(item->data(Qt::DisplayRole).toString().append("_%1").arg(n).append("%1").arg(k+1).append("_NODAL_SLAVE"));

                    this->writeNodalSet(name,anIndexedMapOfFaceMeshDS_Slave);
                    myInputFile<<(QString("*SURFACE,NAME=S_").append(name).append(",TYPE=NODE\n").append(name).append("\n")).toStdString();

                    //! ------------------------------------------------------------
                    //! find the scope of the master, then write an ELEMENT surface
                    //! ------------------------------------------------------------
                    QString mName=itemNameClearSpaces(item->data(Qt::DisplayRole).toString().append("_%1").arg(n).append("%1").arg(k+1).append("_ELEMENT_MASTER"));
                    this->writeElementSurface(mName,anIndexedMapOfFaceMeshDS_Master);

                    masterSlaveName.first = QString("S_").append(name);
                    masterSlaveName.second= mName;
                }
                    break;

                case Property::contactBehavior_symmetric:
                {
                    //! -----------------------------------------------------------
                    //! find the scope of the slave, then write an ELEMENT surface
                    //! -----------------------------------------------------------
                    QString sName=itemNameClearSpaces(item->data(Qt::DisplayRole).toString().append("_%1").arg(n).append("%1").arg(k+1).append("_ELEMENT_SLAVE"));
                    this->writeElementSurface(sName,anIndexedMapOfFaceMeshDS_Slave);

                    //! -----------------------------------------------------------
                    //! find the scope of the master, the write an ELEMENT surface
                    //! -----------------------------------------------------------
                    QString mName=itemNameClearSpaces(item->data(Qt::DisplayRole).toString().append("_%1").arg(n).append("%1").arg(k+1).append("_ELEMENT_MASTER"));
                    this->writeElementSurface(mName,anIndexedMapOfFaceMeshDS_Master);
                    masterSlaveName.first=sName;
                    masterSlaveName.second=mName;
                }
                    break;
                }
                contactMapName.insert(timeTag,masterSlaveName);
            }
        }
    }

    //! ------------------------------------------------------
    //! [3] write the "contact pair" headers: rescan the tree
    //! ------------------------------------------------------
    int NtotCP = 0;     // total number of contact pair
    for(int n=0; n<theConnectionItem->rowCount(); n++)
    {
        if(Global::status().code==0)
        {
            cout<<"writeSolverFileClass::perform()->____process stopped____"<<endl;
            return false;
        }

        //! the current connection group
        QStandardItem *itemConnectionGroup = theConnectionItem->child(n,0);

        //! number of contacts under the current connection group
        int NbContactPairs = itemConnectionGroup->rowCount();

        //! update the total number of contact pairs
        NtotCP = NtotCP + NbContactPairs;

        for(int k=0; k<NbContactPairs; k++)
        {
            QStandardItem *item = itemConnectionGroup->child(k,0);
            SimulationNodeClass *node = item->data(Qt::UserRole).value<SimulationNodeClass*>();
            Property::SuppressionStatus ss = node->getPropertyValue<Property::SuppressionStatus>("Suppressed");
            if(ss==Property::SuppressionStatus_Active)
            {
                QString timeTag = node->getPropertyValue<QString>("Time tag");
                Property::contactType theContactType = node->getPropertyValue<Property::contactType>("Type");
                Property::contactBehavior theContactBehavior = node->getPropertyValue<Property::contactBehavior>("Behavior");
                Property::contactFormulation theContactFormulation = node->getPropertyValue<Property::contactFormulation>("Formulation");
                Property::overpressureFunction theOverPressure = node->getPropertyValue<Property::overpressureFunction>("Overpressure");

                //! Normal stiffness
                double K,KF,KN;
                KF = node->getPropertyValue<double>("K");
                if(KF == 0) KF=1;

                std::vector<GeometryTag> tagsMaster = node->getPropertyValue<std::vector<GeometryTag>>("Tags master");
                std::vector<GeometryTag> tagsSlave = node->getPropertyValue<std::vector<GeometryTag>>("Tags slave");
                QList<occHandle(Ng_MeshVS_DataSourceFace)> masterFaces,slaveFaces;

                for(std::vector<GeometryTag>::iterator it = tagsMaster.begin(); it!= tagsMaster.end(); ++it)
                {
                    GeometryTag aLoc = *it;
                    int bodyIndex = aLoc.parentShapeNr;
                    int faceIndex = aLoc.subTopNr;

                    const occHandle(Ng_MeshVS_DataSourceFace) &faceMesh = occHandle(Ng_MeshVS_DataSourceFace)::DownCast(myDB->ArrayOfMeshDSOnFaces.getValue(bodyIndex,faceIndex));
                    masterFaces<<faceMesh;
                }
                for(std::vector<GeometryTag>::iterator itt = tagsSlave.begin(); itt!= tagsSlave.end(); ++itt)
                {
                    GeometryTag aLoc = *itt;
                    int bodyIndex = aLoc.parentShapeNr;
                    int faceIndex = aLoc.subTopNr;

                    const occHandle(Ng_MeshVS_DataSourceFace) &faceMesh = occHandle(Ng_MeshVS_DataSourceFace)::DownCast(myDB->ArrayOfMeshDSOnFaces.getValue(bodyIndex,faceIndex));
                    slaveFaces<<faceMesh;
                }
                K = contactParameters::calc_K(masterFaces,slaveFaces);
                KN = K*KF;

                //! C0
                double C0 = node->getPropertyValue<double>("C0");

                //! sigmaInfinity
                double sigmaInfty = node->getPropertyValue<double>("Sigma infinity");

                //! Small sliding
                int smallSliding = node->getPropertyValue<int>("Small sliding");;

                //! fiction coefficient
                frictionCoefficient = node->getPropertyValue<double>("Friction coefficient");

                //! lambda
                lambda = node->getPropertyValue<double>("Lambda");
                if(lambda == 0)
                {
                    lambda = KN/20.0;
                }

                //! gap conductance
                gapConductance = node->getPropertyValue<double>("Gap conductance");

                switch(theContactType)
                {
                case Property::contactType_bonded:
                {
                        //! 1st) *TIE,<TOLERANCE: C0>,<NAME OF THE CONNECTION>
                        //! 2st) <SLAVE SURF NAME>,<MASTER SURF NAME>
                        //! ---------------------------------------------------
                        QString slaveName = itemNameClearSpaces(item->data(Qt::DisplayRole).toString().append("_%1").arg(n).append("%1").arg(k+1).append("_NODAL_SLAVE"));
                        QString masterName = itemNameClearSpaces(item->data(Qt::DisplayRole).toString().append("_%1").arg(n).append("%1").arg(k+1).append("_ELEMENT_MASTER"));

                        double C0 = node->getPropertyValue<double>("C0");
                        if(C0==0.0) C0=1.0;
                        //! ---------
                        //! 1st) row
                        //! ---------
                        //myInputFile<<"*TIE, ADJUST = NO, POSITION TOLERANCE="<<C0<<", NAME="<<
                        myInputFile<<"*TIE, POSITION TOLERANCE="<<C0<<", NAME="<<
                                     itemNameClearSpaces((item->data(Qt::DisplayRole).toString().append("_%1").arg(n).append("%1").arg(k+1))).toStdString()<<endl;
                        //! ---------
                        //! 2nd) row
                        //! ---------
                        myInputFile<<(QString("S_").append(slaveName)).toStdString()<<","<<masterName.toStdString()<<"\n";
                    }
                        break;

                    case Property::contactType_frictional:
                    case Property::contactType_frictionless:
                    {
                        //! ----------------------------------------------------------------------------------------------------------
                        //! asymmetric - frictional or frictionless contact
                        //!
                        //! 1st) *CONTACT PAIR,INTERACTION=<name of the contact pair>,TYPE=NODE TO SURFACE,SMALL SLIDING <if defined>
                        //! 2nd) <slave node set>, <master face set>
                        //! 3rd) *SURFACE INTERACTION, NAME = <name of the contact pair>
                        //! 4th) *SURFACE BEHAVIOR, PRESSURE-OVERCLOSURE = <LINEAR, EXPONENTIAL, TABULAR>
                        //! 5th) *FRICTION
                        //! 6th) <friction coefficient>, <lambda>
                        //! ----------------------------------------------------------------------------------------------------------
                        QString name=itemNameClearSpaces(item->data(Qt::DisplayRole).toString().append("_%1").arg(n).append("%1").arg(k+1));
                        Property::overpressureFunction theOverPressure = node->getPropertyValue<Property::overpressureFunction>("Overpressure");
                        int smallSliding = node->getPropertyValue<int>("Small sliding");;

                        //! ---------
                        //! 1st) row
                        //! ---------
                        myInputFile<<"*CONTACT PAIR, INTERACTION = "<<name.toStdString()<<", TYPE = NODE TO SURFACE";
                        if(smallSliding==1) myInputFile<<", SMALL SLIDING"<<endl; else myInputFile<<endl;

                        //! ---------
                        //! 2nd) row
                        //! ---------
                        myInputFile<<(QString("S_").append(itemNameClearSpaces(item->data(Qt::DisplayRole).toString().append("_%1").arg(n).append("%1").arg(k+1))).
                                      append("_NODAL_SLAVE")).toStdString()<<", "<<
                                     (itemNameClearSpaces(item->data(Qt::DisplayRole).toString()).append("_%1").arg(n).append("%1").arg(k+1).append("_ELEMENT_MASTER")).toStdString()<<"\n";

                        //! ---------
                        //! 3rd) row
                        //! ---------
                        myInputFile<<"*SURFACE INTERACTION, NAME = "<<(itemNameClearSpaces(item->data(Qt::DisplayRole).toString()).append("_%1").arg(n).append("%1").arg(k+1)).toStdString()<<"\n";

                        //! ---------
                        //! 4th) row
                        //! ---------
                        myInputFile<<"*SURFACE BEHAVIOR, PRESSURE-OVERCLOSURE=";
                        switch(theOverPressure)
                        {
                        case Property::overpressureFunction_linear:
                        {
                            double C0 = node->getPropertyValue<double>("C0");
                            KF = node->getPropertyValue<double>("K");
                            double sigmaInfty = node->getPropertyValue<double>("Sigma infinity");

                            myInputFile<<"LINEAR"<<endl;
                            if(KF == 0) KF=1;
                            //{
                            std::vector<GeometryTag> tagsMaster = node->getPropertyValue<std::vector<GeometryTag>>("Tags master");
                            std::vector<GeometryTag> tagsSlave = node->getPropertyValue<std::vector<GeometryTag>>("Tags slave");
                            QList<occHandle(Ng_MeshVS_DataSourceFace)> masterFaces,slaveFaces;

                            for(std::vector<GeometryTag>::iterator it = tagsMaster.begin(); it!= tagsMaster.end(); ++it)
                            {
                                GeometryTag aLoc = *it;
                                int bodyIndex = aLoc.parentShapeNr;
                                int faceIndex = aLoc.subTopNr;

                                const occHandle(Ng_MeshVS_DataSourceFace) &faceMesh = occHandle(Ng_MeshVS_DataSourceFace)::DownCast(myDB->ArrayOfMeshDSOnFaces.getValue(bodyIndex,faceIndex));
                                masterFaces<<faceMesh;
                            }
                            for(std::vector<GeometryTag>::iterator itt = tagsSlave.begin(); itt!= tagsSlave.end(); ++itt)
                            {
                                GeometryTag aLoc = *itt;
                                int bodyIndex = aLoc.parentShapeNr;
                                int faceIndex = aLoc.subTopNr;

                                const occHandle(Ng_MeshVS_DataSourceFace) &faceMesh = occHandle(Ng_MeshVS_DataSourceFace)::DownCast(myDB->ArrayOfMeshDSOnFaces.getValue(bodyIndex,faceIndex));
                                slaveFaces<<faceMesh;
                            }
                            K = contactParameters::calc_K(masterFaces,slaveFaces);
                            //}
                            if(sigmaInfty == 0.0)
                            {
                                //! calculix suggest 0.25% of the maximum stress expected
                                //! we use 0.25% of the tensile yield strenght // TO DO....
                                sigmaInfty = 0.25*750;
                            }
                            if(C0 == 0.0)
                            {
                                //! calculix default 10e-3
                                C0 = 10.0e-3;
                            }
                            myInputFile<<K*KF<<", "<<sigmaInfty<<", "<<C0<<endl;
                        }
                            break;
                        case Property::overpressureFunction_exponential:
                        {
                            double P0 = node->getPropertyValue<double>("P0");
                            double C0 = node->getPropertyValue<double>("C0");
                            myInputFile<<"EXPONENTIAL"<<endl;
                            myInputFile<<C0<<", "<<P0<<endl;
                        }
                            break;
                        }

                        if(theContactType == Property::contactType_frictional)
                        {
                            double frictionCoefficient;
                            double lambda;
                            frictionCoefficient = node->getPropertyValue<double>("Friction coefficient");
                            lambda = node->getPropertyValue<double>("Lambda");
                            if(lambda == 0)
                            {
                                lambda = K*KF/20.0;
                            }
                            //! ---------
                            //! 5th) row
                            //! ---------
                            myInputFile<<"*FRICTION"<<endl;

                            //! ---------
                            //! 6th) row
                            //! ---------
                            myInputFile<<frictionCoefficient<<", "<<lambda<<endl;
                        }
                        //! ----------------
                        //! GAP CONDUCTANCE
                        //! ----------------
                        myInputFile<<"*GAP CONDUCTANCE "<<endl;
                        //! default parameters conductance of 100 for all contact pressure and all temeprature
                        myInputFile<<"1e9,,273 "<<endl;
                    }
                        break;
                    }
                }
                    break;
                case Property::contactBehavior_symmetric:
                {
                    switch(theContactType)
                    {
                    case Property::contactType_bonded:
                    {
                        //! ---------------------------------------------------
                        //! symmetric - bonded contact
                        //!
                        //! 1st) *TIE,<TOLERANCE: C0>,<NAME OF THE CONNECTION>
                        //! 2st) <SLAVE SURF NAME>,<MASTER SURF NAME>
                        //! ---------------------------------------------------
                        QString slaveName = itemNameClearSpaces(item->data(Qt::DisplayRole).toString().append("_%1").arg(n).append("%1").arg(k+1).append("_ELEMENT_SLAVE"));
                        QString masterName = itemNameClearSpaces(item->data(Qt::DisplayRole).toString().append("_%1").arg(n).append("%1").arg(k+1).append("_ELEMENT_MASTER"));
                        double C0 = node->getPropertyValue<double>("C0");
                        if (C0==0.0) C0=10.0;

                        //! ---------
                        //! 1st) row
                        //! ---------
                        //myInputFile<<"*TIE, ADJUST=NO, POSITION TOLERANCE="<<C0<<", NAME="<<
                        myInputFile<<"*TIE, POSITION TOLERANCE="<<C0<<", NAME="<<
                                     itemNameClearSpaces((item->data(Qt::DisplayRole).toString().append("_%1").arg(n).append("%1").arg(k+1))).toStdString()<<endl;

                        //! ---------
                        //! 2nd) row
                        //! ---------
                        myInputFile<<slaveName.toStdString()<<","<<masterName.toStdString()<<"\n";
                    }
                        break;
                    case Property::contactType_frictional:
                    case Property::contactType_frictionless:
                    {
                        //! -------------------------------------------------------------------------------------------------------------
                        //! symmetric - frictional or frictionless contact
                        //!
                        //! 1st) *CONTACT PAIR,INTERACTION=<name of the contact pair>,TYPE=SURFACE TO SURFACE,SMALL SLIDING <if defined>
                        //! 2nd) <slave face set>, <master face set>
                        //! 3rd) *SURFACE INTERACTION, NAME = <name of the contact pair>
                        //! 4th) *SURFACE BEHAVIOR, PRESSURE-OVERCLOSURE = <LINEAR, EXPONENTIAL, TABULAR>
                        //! 5th) *FRICTION
                        //! 6th) <friction coefficient>, <lambda>
                        //! -------------------------------------------------------------------------------------------------------------
                        QString slaveName = itemNameClearSpaces(item->data(Qt::DisplayRole).toString().append("_%1").arg(n).append("%1").arg(k+1).append("_ELEMENT_SLAVE"));
                        QString masterName = itemNameClearSpaces(item->data(Qt::DisplayRole).toString().append("_%1").arg(n).append("%1").arg(k+1).append("_ELEMENT_MASTER"));

                        Property::overpressureFunction theOverPressure = node->getPropertyValue<Property::overpressureFunction>("Overpressure");;

                        //! -------------------------------------------------------------------------
                        //! 1st) row - Warning: here the "SMALL SLIDING" parameter cannot be defined
                        //! -------------------------------------------------------------------------
                        myInputFile<<"*CONTACT PAIR, INTERACTION = "<<(itemNameClearSpaces(item->data(Qt::DisplayRole).toString()).append("_%1").arg(n).append("%1").arg(k+1)).toStdString()<<
                                     ", TYPE = SURFACE TO SURFACE\n";

                        //! ---------
                        //! 2nd) row
                        //! ---------
                        myInputFile<<slaveName.toStdString()<<", "<< masterName.toStdString()<<"\n";

                        //! ---------
                        //! 3rd) row
                        //! ---------
                        myInputFile<<"*SURFACE INTERACTION, NAME = "<<(itemNameClearSpaces(item->data(Qt::DisplayRole).toString()).append("_%1").arg(n).append("%1").arg(k+1)).toStdString()<<"\n";

                        //! ---------
                        //! 4th) row
                        //! ---------
                        myInputFile<<"*SURFACE BEHAVIOR, PRESSURE-OVERCLOSURE=";
                        switch(theOverPressure)
                        {
                        case Property::overpressureFunction_linear:
                        {
                            myInputFile<<"LINEAR"<<endl;

                            //! in case of face to face behavior "Sigma infty" is irrelevant
                            //! in case of face to face behavior "C0 is irrelevant

                            KF = node->getPropertyValue<double>("K");
                            if(KF==0) KF=1.0;
                            //{
                            std::vector<GeometryTag> tagsMaster = node->getPropertyValue<std::vector<GeometryTag>>("Tags master");
                            std::vector<GeometryTag> tagsSlave = node->getPropertyValue<std::vector<GeometryTag>>("Tags slave");

                            QList<occHandle(Ng_MeshVS_DataSourceFace)> masterFaces,slaveFaces;

                                for(std::vector<GeometryTag>::iterator it = tagsMaster.begin(); it!= tagsMaster.end(); ++it)
                                {
                                    GeometryTag aLoc = *it;
                                    int bodyIndex = aLoc.parentShapeNr;
                                    int faceIndex = aLoc.subTopNr;

                                    const occHandle(Ng_MeshVS_DataSourceFace) &faceMesh = occHandle(Ng_MeshVS_DataSourceFace)::DownCast(myDB->ArrayOfMeshDSOnFaces.getValue(bodyIndex,faceIndex));
                                    masterFaces<<faceMesh;
                                }
                                for(std::vector<GeometryTag>::iterator itt = tagsSlave.begin(); itt!= tagsSlave.end(); ++itt)
                                {
                                    GeometryTag aLoc = *itt;
                                    int bodyIndex = aLoc.parentShapeNr;
                                    int faceIndex = aLoc.subTopNr;

                                    const occHandle(Ng_MeshVS_DataSourceFace) &faceMesh = occHandle(Ng_MeshVS_DataSourceFace)::DownCast(myDB->ArrayOfMeshDSOnFaces.getValue(bodyIndex,faceIndex));
                                    slaveFaces<<faceMesh;
                                }
                                K = contactParameters::calc_K(masterFaces,slaveFaces);
                            //}
                            myInputFile<<K*KF<<endl;
                        }
                            break;

                        case Property::overpressureFunction_exponential:
                        {
                            myInputFile<<"EXPONENTIAL"<<endl;
                            double C0 = node->getPropertyValue<double>("C0");
                            double P0 = node->getPropertyValue<double>("P0");
                            myInputFile<<C0<<", "<<P0<<endl;
                        }
                            break;
                        }

                        if(theContactType == Property::contactType_frictional)
                        {
                            double frictionCoefficient;
                            double lambda;
                            frictionCoefficient = node->getPropertyValue<double>("Friction coefficient");;
                            lambda = node->getPropertyValue<double>("Lambda");

                            if(lambda==0)
                            {
                                lambda = K*KF/20.0;
                            }
                            //! ---------
                            //! 5th) row
                            //! ---------
                            myInputFile<<"*FRICTION"<<endl;

                            //! ---------
                            //! 6th) row
                            //! ---------
                            myInputFile<<frictionCoefficient<<", "<<lambda<<endl;
                        }
                        //! ----------------
                        //! GAP CONDUCTANCE
                        //! ----------------
                        myInputFile<<"*GAP CONDUCTANCE "<<endl;
                        //! default parameters conductance of 1e9 for all contact pressure and all temeprature
                        myInputFile<<"1e9,,273 "<<endl;
                    }
                        break;
                    case Property::contactType_noSeparation:
                    {
                        QString slaveName = itemNameClearSpaces(item->data(Qt::DisplayRole).toString().append("_%1").arg(n).append("%1").arg(k+1).append("_ELEMENT_SLAVE"));
                        QString masterName = itemNameClearSpaces(item->data(Qt::DisplayRole).toString().append("_%1").arg(n).append("%1").arg(k+1).append("_ELEMENT_MASTER"));

                        //! -------------------------------------------------------------------------
                        //! 1st) row - Warning: here the "SMALL SLIDING" parameter cannot be defined
                        //! -------------------------------------------------------------------------
                        myInputFile<<"*CONTACT PAIR, INTERACTION = "<<(itemNameClearSpaces(item->data(Qt::DisplayRole).toString()).append("_%1").arg(n).append("%1").arg(k+1)).toStdString()<<
                                     ", TYPE = SURFACE TO SURFACE\n";

                        //! ---------
                        //! 2nd) row
                        //! ---------
                        myInputFile<<slaveName.toStdString()<<", "<< masterName.toStdString()<<"\n";

                        //! ---------
                        //! 3rd) row
                        //! ---------
                        myInputFile<<"*SURFACE INTERACTION, NAME = "<<(itemNameClearSpaces(item->data(Qt::DisplayRole).toString()).append("_%1").arg(n).append("%1").arg(k+1)).toStdString()<<"\n";

                        //! ---------
                        //! 4th) row
                        //! ---------
                        myInputFile<<"*SURFACE BEHAVIOR, PRESSURE-OVERCLOSURE=";
                        myInputFile<<"TIED"<<endl;

                        //! in case of face to face behavior "Sigma infty" is irrelevant
                        //! in case of face to face behavior "C0 is irrelevant

                        KF = node->getPropertyValue<double>("K");
                        if(K==0) KF=1.0;
                        //{
                            std::vector<GeometryTag> tagsMaster = node->getPropertyValue<std::vector<GeometryTag>>("Tags master");
                            std::vector<GeometryTag> tagsSlave = node->getPropertyValue<std::vector<GeometryTag>>("Tags slave");

                            QList<occHandle(Ng_MeshVS_DataSourceFace)> masterFaces,slaveFaces;

                            for(std::vector<GeometryTag>::iterator it = tagsMaster.begin(); it!= tagsMaster.end(); ++it)
                            {
                                GeometryTag aLoc = *it;
                                int bodyIndex = aLoc.parentShapeNr;
                                int faceIndex = aLoc.subTopNr;

                                const occHandle(Ng_MeshVS_DataSourceFace) &faceMesh = occHandle(Ng_MeshVS_DataSourceFace)::DownCast(myDB->ArrayOfMeshDSOnFaces.getValue(bodyIndex,faceIndex));
                                masterFaces<<faceMesh;
                            }
                            for(std::vector<GeometryTag>::iterator itt = tagsSlave.begin(); itt!= tagsSlave.end(); ++itt)
                            {
                                GeometryTag aLoc = *itt;
                                int bodyIndex = aLoc.parentShapeNr;
                                int faceIndex = aLoc.subTopNr;

                                const occHandle(Ng_MeshVS_DataSourceFace) &faceMesh = occHandle(Ng_MeshVS_DataSourceFace)::DownCast(myDB->ArrayOfMeshDSOnFaces.getValue(bodyIndex,faceIndex));
                                slaveFaces<<faceMesh;
                            }
                            K = contactParameters::calc_K(masterFaces,slaveFaces);
                        //}
                        myInputFile<<K*KF<<endl;

                        //! ----------------
                        //! GAP CONDUCTANCE
                        //! ----------------
                        myInputFile<<"*GAP CONDUCTANCE "<<endl;
                        //! default parameters conductance of 1e9 for all contact pressure and all temeprature
                        myInputFile<<"1000000000,,273 "<<endl;

                        //! ---------
                        //! 5th) row
                        //! ---------
                        myInputFile<<"*FRICTION"<<endl;

                        double lambda = K*KF/20.0;
                        //! ---------
                        //! 6th) row
                        //! ---------
                        myInputFile<<1<<", "<<lambda<<endl;
                    }
                        break;
                    }
                }
                    break;
                }
            }
        }
    }

    //! --------------------
    //! update the progress
    //! --------------------
    if(myProgressIndicator!=Q_NULLPTR)
    {
        done++;
        QProgressEvent *e = new QProgressEvent(QProgressEvent_Update,0,Nevents-1,done,"Sending contacts",
                                               QProgressEvent_None,-1,-1,-1,"Writing CCX solver input file");
        QApplication::postEvent(myProgressIndicator,e);
        QApplication::processEvents();
        QThread::msleep(500);

        if(Global::status().code==0)
        {
            cout<<"writeSolverFileClass::perform()->____process stopped____"<<endl;
            return false;
        }
    }

    //! ------------------------------------------------------------------------------------
    //! [4] write the "boundary" headers: rescan the tree starting from "Static structural"
    //! ------------------------------------------------------------------------------------
    for(int k=1; k<mySimulationRoot->rowCount()-1; k++)
    {
        if(Global::status().code==0)
        {
            cout<<"writeSolverFileClass::perform()->____process stopped____"<<endl;
            return false;
        }

        //if(mySimulationRoot==Q_NULLPTR)
        //{
        //    cerr<<"writeSolverFileClass::perform()->____the simulation root is NULL____"<<endl;
        //    exit(100);
        //}
        QString itemName = itemNameClearSpaces(mySimulationRoot->child(k,0)->data(Qt::DisplayRole).toString());

        SimulationNodeClass *theCurNode = mySimulationRoot->child(k,0)->data(Qt::UserRole).value<SimulationNodeClass*>();
        SimulationNodeClass::nodeType theNodeType = theCurNode->getType();
        Property::SuppressionStatus theNodeSS = theCurNode->getPropertyValue<Property::SuppressionStatus>("Suppressed");

        if(theNodeSS==Property::SuppressionStatus_Active)
        {
            QString SetName;
            switch(theNodeType)
            {
            case SimulationNodeClass::nodeType_structuralAnalysisBoundaryContidion_FixedSupport:
                SetName = itemName.append("_").append(QString("%1").arg(k));
                myInputFile<<"*BOUNDARY"<<endl;
                myInputFile<<SetName.toStdString()<<",1,3"<<endl;
                break;
            case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_CylindricalSupport:
                SetName = itemName.append("_").append(QString("%1").arg(k));
            {
                /*
                gp_Ax1 theCylAxis;
                gp_Pnt P0, P1;

                //! retrieve the scope through Tags
                std::vector<GeometryTag> vecLoc = theCurNode->getPropertyItem("Tags")->data(Qt::UserRole).value<Property>().getData().value<std::vector<GeometryTag>>();

                int i=0;
                for(std::vector<GeometryTag>::iterator it = vecLoc.begin(); it!=vecLoc.end(); ++it)
                {
                    ++i;
                    GeometryTag loc = *it;
                    int bodyIndex = loc.parentShapeNr;
                    int childShape = loc.subTopNr;
                    const TopoDS_Shape &scope = myDB->MapOfBodyTopologyMap.value(bodyIndex).faceMap.FindKey(childShape);
                    QString SetName_it=SetName+QString("%1").arg(i);

                    TopoDS_Face theCylFace = TopoDS::Face(scope);

                    GeomToolsClass::getCylindricalFaceInfo(theCylFace,theCylAxis,P0);
                    x0 = P0.X(); y0 = P0.Y(); z0 = P0.Z();
                    gp_Vec vec = gp_Vec(theCylAxis.Direction()).Multiplied(1.0);
                    P1 = P0;
                    P1.Translate(vec);
                    x1 = P1.X(); y1 = P1.Y(); z1 = P1.Z();
                    myInputFile<<"*TRANSFORM, NSET="<<SetName_it.toStdString()<<", TYPE=C"<<endl;
                    myInputFile<<x0<<","<<y0<<","<<z0<<","<<x1<<","<<y1<<","<<z1<<endl;

                    myInputFile<<"*BOUNDARY"<<endl;

                    Property::DOFfreedom theDOFfreedom_radial = theCurNode->getPropertyItem("Radial")->data(Qt::UserRole).value<Property>().getData().value<Property::DOFfreedom>();
                    Property::DOFfreedom theDOFfreedom_axial = theCurNode->getPropertyItem("Axial")->data(Qt::UserRole).value<Property>().getData().value<Property::DOFfreedom>();
                    Property::DOFfreedom theDOFfreedom_tangential = theCurNode->getPropertyItem("Tangential")->data(Qt::UserRole).value<Property>().getData().value<Property::DOFfreedom>();
                    bool br_isFixed = theDOFfreedom_radial==Property::DOFfreedom_fixed? true:false;
                    bool ba_isFixed = theDOFfreedom_axial==Property::DOFfreedom_fixed? true:false;
                    bool bt_isFixed = theDOFfreedom_tangential==Property::DOFfreedom_fixed? true:false;
                    if(br_isFixed)myInputFile<<SetName_it.toStdString()<<",1,1"<<endl;
                    if(ba_isFixed)myInputFile<<SetName_it.toStdString()<<",3,3"<<endl;
                    if(bt_isFixed)myInputFile<<SetName_it.toStdString()<<",2,2"<<endl;

                }*/
            }
                break;
            case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_CompressionOnlySupport:
            {
                ;
            }
                break;
            //! boundary condition for now written in the writing set section
            case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_FrictionlessSupport:
            {
                //! the set name
                //SetName = itemName+QString("_%1").arg(k);

                //IndexedMapOfMeshDataSources anIndexedMapOfFaceMeshDS;
                //anIndexedMapOfFaceMeshDS = theItemNode->getPropertyValue<IndexedMapOfMeshDataSources>("Mesh data sources");

                //SetName = itemName+append("_").append(QString("%1").arg(k));
                /*
                //! retrieve the the scope using Tags
                std::vector<GeometryTag> vecLoc = theCurNode->getPropertyValue<std::vector<GeometryTag>>("Tags");
                int i=0;
                ListOfShape scope;
                for(std::vector<GeometryTag>::iterator it = vecLoc.begin(); it!=vecLoc.end(); ++it)
                {
                    ++i;
                    QString SName=SetName+QString("%1").arg(i);
                    GeometryTag loc = *it;
                    int bodyIndex = loc.parentShapeNr;
                    int childShape = loc.subTopNr;
                    TopoDS_Shape aShape = myDB->MapOfBodyTopologyMap.value(bodyIndex).faceMap.FindKey(childShape);

                    if(aShape.IsNull())
                    {
                        cout<<"____NULL shape____"<<endl;
                        exit(1);
                    }
                    scope.Append(aShape);

                    TopoDS_Face theFace = TopoDS::Face(scope.First());
                    if(theFace.IsNull())
                    {
                        cout<<"____NULL face____"<<endl;
                        exit(1);
                    }

                    //! find a vertex of the face
                    TopExp_Explorer anExp(theFace,TopAbs_VERTEX);
                    TopoDS_Vertex aVertex = TopoDS::Vertex(anExp.Current());
                    if(aVertex.IsNull())
                    {
                        cout<<"____NULL vertex____"<<endl;
                        exit(1);
                    }
                    gp_Pnt C1 = BRep_Tool::Pnt(aVertex);

                    GeomAbs_SurfaceType theSurfaceType;
                    GeomToolsClass::getFaceType(theFace,theSurfaceType);
                    if(theSurfaceType==GeomAbs_Plane)
                    {
                        //! handle a planar face
                        //! "theAxis" is normal to the planar face
                        gp_Pnt C0;                                  //! centroid
                        gp_Ax1 theAxis;                             //! the Axis
                        GeomToolsClass::getPlanarFaceInfo(theFace,theAxis,C0);

                        double x0,y0,z0,x1,y1,z1;
                        x0 = C0.X(); y0 = C0.Y(); z0 = C0.Z();

                        //! create the point (a): start from global origin and
                        //! translate of 1.0 along the normal
                        gp_Pnt a(0,0,0);  //! origin
                        gp_Vec unitVector(theAxis.Direction());
                        a.Translate(unitVector);

                        //! preliminary to the creation of the point (b)
                        //! create a vector from the face centroid to one vertex
                        gp_Vec C0C1(C0,C1);

                        cout<<" Centroid (x0,y0,z0) = ("<<x0<<", "<<y0<<", "<<z0<<")"<<endl;
                        x1 = C1.X(); y1 = C1.Y(); z1 = C1.Z();
                        cout<<"A vertex (x1,y1,z1) = ("<<x1<<", "<<y1<<", "<<z1<<")"<<endl;

                        cout<<" point (a) = ("<<a.X()<<", "<<a.Y()<<", "<<a.Z()<<")"<<endl;

                        //! create the point (b)
                        gp_Pnt b(0,0,0);    //! origin
                        gp_Dir dirb(C0C1);
                        b.Translate(dirb);

                        cout<<" point (b) = ("<<b.X()<<", "<<b.Y()<<", "<<b.Z()<<")"<<endl;

                        myInputFile<<"*TRANSFORM, NSET="<<SName.toStdString()<<", TYPE=R"<<endl;
                        myInputFile<<a.X()<<","<<a.Y()<<","<<a.Z()<<","<<b.X()<<","<<b.Y()<<","<<b.Z()<<endl;

                        myInputFile<<"*BOUNDARY"<<endl;
                        //! lock local X DOF since the local CS is defined using the surface normal
                        //! vector as X local axis
                        myInputFile<<SName.toStdString()<<",1"<<endl;

                    }
                    else if(theSurfaceType==GeomAbs_Cylinder)
                    {
                        //! handle a cylindrical face
                        //! what follows is exactly equivalent to a "Cylindrical support"
                        //! with radial DOF locked
                        gp_Pnt C0;                                  //! origin of the symmetry axis
                        gp_Ax1 theAxis;                             //! the Axis
                        GeomToolsClass::getCylindricalFaceInfo(theFace,theAxis,C0);

                        cout<<"Origin of the cylindrical axis: ("<<C0.X()<<", "<<C0.Y()<<", "<<C0.Z()<<")"<<endl;
                        gp_Pnt a = C0;
                        gp_Pnt b = a;
                        b.Translate(theAxis.Direction());

                        cout<<"point (a) = ("<<a.X()<<", "<<a.Y()<<", "<<a.Z()<<")"<<endl;
                        cout<<"point (b) = ("<<b.X()<<", "<<b.Y()<<", "<<b.Z()<<")"<<endl;

                        myInputFile<<"*TRANSFORM, NSET="<<SName.toStdString()<<", TYPE=R"<<endl;
                        myInputFile<<a.X()<<","<<a.Y()<<","<<a.Z()<<","<<b.X()<<","<<b.Y()<<","<<b.Z()<<endl;

                        myInputFile<<"*BOUNDARY"<<endl;
                        //! lock local X DOF which correspond to the cyl radial direction
                        myInputFile<<SName.toStdString()<<",1"<<endl;

                    }
                    else if(theSurfaceType==GeomAbs_Cone || theSurfaceType==GeomAbs_Sphere ||
                            theSurfaceType==GeomAbs_BezierSurface || theSurfaceType ==GeomAbs_BSplineSurface)
                    {
                        //! These cases cannot be handled using an unique local CS defined
                        //! for the whole NSET. The nodes of the surface must be scanned one by one,
                        //! and foreach node the TRANSFORM card, followed by BOUNDARY must be used
                        //! to do ...
                    }
                }
            */
            }
                break;
            }
        }
    }

    //! --------------------
    //! update the progress
    //! --------------------
    if(myProgressIndicator!=Q_NULLPTR)
    {
        done++;
        QProgressEvent *e = new QProgressEvent(QProgressEvent_Update,0,Nevents-1,done,"Sending boundary conditions",
                                               QProgressEvent_None,-1,-1,-1,"Writing CCX solver input file");
        QApplication::postEvent(myProgressIndicator,e);
        QApplication::processEvents();
        QThread::msleep(500);

        if(Global::status().code==0)
        {
            cout<<"writeSolverFileClass::perform()->____process stopped____"<<endl;
            return false;
        }
    }

    /*
    //! ----------------------------------
    //! [5] write the material properties
    //! ----------------------------------
    QStandardItem *theGeometryRoot=this->getTreeItem(SimulationNodeClass::nodeType_geometry);
    for(int k=0; k<theGeometryRoot->rowCount();k++)
    {
        QStandardItem *theGeometryItem = theGeometryRoot->child(k,0);
        SimulationNodeClass *theCurNode = theGeometryItem->data(Qt::UserRole).value<SimulationNodeClass*>();
        Property::SuppressionStatus theNodeSS = theCurNode->getPropertyValue<Property::SuppressionStatus>("Suppressed");

        if(theNodeSS==Property::SuppressionStatus_Active)
        {
            int matNumber = theCurNode->getPropertyValue<int>("Assignment");
            std::string matName = vecMatNames.at(matNumber);

            QString materialDBPath = QCoreApplication::applicationDirPath();
            QString matIncludeFileAbsPosition = materialDBPath+"/material/"+QString::fromStdString(matName)+".mat";
            std::string bodyName =theGeometryItem->data(Qt::DisplayRole).toString().toStdString();
            myInputFile<<"*MATERIAL, Name="<<matName<<endl;
            myInputFile<<"*INCLUDE, INPUT="<<matIncludeFileAbsPosition.toStdString()<<endl;
            myInputFile<<"*SOLID SECTION, ELSET= E"<<bodyName<<", MATERIAL="<<matName<<endl;
        }
    }
    */
    //! ----------------------------------
    //! [5] write the material properties
    //! ----------------------------------
    for(int p=0;p<vecMatNames.size();p++)
    {
        std::string matName = vecMatNames.at(p);
        QString materialDBPath = QCoreApplication::applicationDirPath();
        QString matIncludeFileAbsPosition = materialDBPath+"/material/"+QString::fromStdString(matName)+".mat";
        myInputFile<<"*MATERIAL, Name="<<matName<<endl;
        myInputFile<<"*INCLUDE, INPUT="<<matIncludeFileAbsPosition.toStdString()<<endl;
    }

    for(int k=0; k<theGeometryRoot->rowCount();k++)
    {
        if(Global::status().code==0)
        {
            cout<<"writeSolverFileClass::perform()->____process stopped____"<<endl;
            return false;
        }

        QStandardItem *theGeometryItem = theGeometryRoot->child(k,0);
        SimulationNodeClass *theCurNode = theGeometryItem->data(Qt::UserRole).value<SimulationNodeClass*>();
        Property::SuppressionStatus theNodeSS = theCurNode->getPropertyValue<Property::SuppressionStatus>("Suppressed");

        if(theNodeSS==Property::SuppressionStatus_Active && theCurNode->getType()!=SimulationNodeClass::nodeType_pointMass)
        {
            int matNumber = theCurNode->getPropertyValue<int>("Assignment");
            std::string matName = vecMatNames.at(matNumber);
            std::string bodyName =theGeometryItem->data(Qt::DisplayRole).toString().toStdString();
            myInputFile<<"*SOLID SECTION, ELSET= E"<<bodyName<<", MATERIAL="<<matName<<endl;
        }
    }

    //! --------------------
    //! update the progress
    //! --------------------
    if(myProgressIndicator!=Q_NULLPTR)
    {
        done++;
        QProgressEvent *e = new QProgressEvent(QProgressEvent_Update,0,Nevents-1,done,"Sending materials",
                                               QProgressEvent_None,-1,-1,-1,"Writing CCX solver input file");
        QApplication::postEvent(myProgressIndicator,e);
        QApplication::processEvents();
        QThread::msleep(500);

        if(Global::status().code==0)
        {
            cout<<"writeSolverFileClass::perform()->____process stopped____"<<endl;
            return false;
        }
    }

    //! ----------------------------
    //! [6] write initial condition
    //! ----------------------------
#ifdef COSTAMP_VERSION
    bool initialTempDistr = false;
    for(int k=1; k<mySimulationRoot->rowCount()-1; k++)
    {
        code = Global::status().code;
        if(code==0) return;

        cout<<" - writing BC "<<mySimulationRoot->child(k,0)->data(Qt::DisplayRole).toString().toStdString()<<endl;
        //QString itemName = itemNameClearSpaces(mySimulationRoot->child(k,0)->data(Qt::DisplayRole).toString());

        SimulationNodeClass *theCurNode = mySimulationRoot->child(k,0)->data(Qt::UserRole).value<SimulationNodeClass*>();
        SimulationNodeClass::nodeType theNodeType= theCurNode->getType();
        QExtendedStandardItem *theCurItem = static_cast<QExtendedStandardItem*>(mySimulationRoot->child(k,0));

        Property::SuppressionStatus theNodeSS = theCurNode->getPropertyItem("Suppressed")->data(Qt::UserRole).value<Property>().getData().value<Property::SuppressionStatus>();

        if(theNodeSS==Property::SuppressionStatus_Active &&
                theNodeType==SimulationNodeClass::nodeType_mapper )
        {
            //! scan the imported body scalar under the "Mapper" item
            int n=0;
            QExtendedStandardItem *itemBodyScalar = static_cast<QExtendedStandardItem*>(theCurItem->child(0,0));
            SimulationNodeClass *ImportedBodyScalarNode = theCurItem->child(0,0)->data(Qt::UserRole).value<SimulationNodeClass*>();
            QString itemName = itemNameClearSpaces(mySimulationRoot->child(k,0)->data(Qt::DisplayRole).toString());
            SimulationNodeClass::nodeType ImportedBodyScalarType= ImportedBodyScalarNode->getType();
            if(ImportedBodyScalarType!=SimulationNodeClass::nodeType_OpenFoamScalarData)
            {
                QString extension=".t";
                QString index=QString::number(n);
                //QString time=QString::number(i);
                //QString tname = "Initial_"+SetName+"_"+index+"_"+time+extension;
                QString tname = "Initial_"+itemName+"_"+index+"_"+extension;
                QString absFileName = myFileName.split("/").last();
                QString dirName = myFileName;
                dirName.chop(absFileName.size());
                //tname.prepend(dirName);
                postObject pObject;

                //! scan the interpolation results at different times
                SimulationNodeClass *node = itemBodyScalar->child(n,0)->data(Qt::UserRole).value<SimulationNodeClass*>();
                Property::SuppressionStatus ss = node->getPropertyValue<Property::SuppressionStatus>("Suppressed");
                if(ss==Property::SuppressionStatus_Active)
                {
                    //!search the Graphic Object property and extract temperature distribution data map
                    pObject =  node->getPropertyItem("Post object")->data(Qt::UserRole).value<Property>().getData().value<postObject>();
                    this->writeTemperatureHistory(pObject, tname);
                }
                myInputFile<<"*INITIAL CONDITIONS, TYPE=TEMPERATURE"<<endl;
                myInputFile<<"*INCLUDE, INPUT = "<<tname.toStdString()<<endl;
                initialTempDistr = true;
            }
        }
    }
    if(initialTempDistr==false)
    {
        double theEnvTemperature = mySimulationRoot->data(Qt::UserRole).value<SimulationNodeClass*>()->getPropertyValue<double>("Environment temperature");
        myInputFile<<"*INITIAL CONDITIONS, TYPE=TEMPERATURE"<<endl;
        myInputFile<<"NALL, "<<theEnvTemperature<<endl;
    }
#endif
#ifndef COSTAMP_VERSION

    bool initialTempDistr = false;
    for(int k=1; k<mySimulationRoot->rowCount()-1; k++)
    {
        if(Global::status().code==0)
        {
            cout<<"writeSolverFileClass::perform()->____process stopped____"<<endl;
            return false;
        }

        SimulationNodeClass *theCurNode = mySimulationRoot->child(k,0)->data(Qt::UserRole).value<SimulationNodeClass*>();
        SimulationNodeClass::nodeType theNodeType= theCurNode->getType();

        Property::SuppressionStatus theNodeSS = theCurNode->getPropertyValue<Property::SuppressionStatus>("Suppressed");
        if(theNodeSS == Property::SuppressionStatus_Suppressed) continue;

        if(theNodeType!=SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_ImportedTemperatureDistribution) continue;
        int couplingTime = theCurNode->getPropertyValue<int>("Structural time step");
        if(couplingTime==0)
        {
            void *p =theCurNode->getPropertyValue<void*>("Imported body temperature");
            QStandardItem* itemTemperatureDist = (QStandardItem*)(p);
            SimulationNodeClass *nodeTemperatureDist = itemTemperatureDist->data(Qt::UserRole).value<SimulationNodeClass*>();
            QStandardItem *itemPostObject = nodeTemperatureDist->getPropertyItem("Post object");
            if(itemPostObject!=Q_NULLPTR)
            {
                QString itemName = itemNameClearSpaces(mySimulationRoot->child(k,0)->data(Qt::DisplayRole).toString());
                QString extension=".t";

                QString tname = "Initial_"+itemName+extension;
                QString absFileName = myFileName.split("/").last();
                QString dirName = myFileName;
                dirName.chop(absFileName.size());
                tname.prepend(dirName);

                myInputFile<<"*INITIAL CONDITIONS, TYPE=TEMPERATURE"<<endl;
                postObject pObject = itemPostObject->data(Qt::UserRole).value<Property>().getData().value<postObject>();
                this->writeTemperatureHistory(pObject, tname);
            }
            initialTempDistr = true;
            break;
        }
    }
    if(initialTempDistr==false)
    {
        double theEnvTemperature = mySimulationRoot->data(Qt::UserRole).value<SimulationNodeClass*>()->getPropertyValue<double>("Environment temperature");
        myInputFile<<"*INITIAL CONDITIONS, TYPE=TEMPERATURE"<<endl;
        myInputFile<<"NALL, "<<theEnvTemperature<<endl;
    }
#endif

    //! --------------------
    //! update the progress
    //! --------------------
    if(myProgressIndicator!=Q_NULLPTR)
    {
        done++;
        QProgressEvent *e = new QProgressEvent(QProgressEvent_Update,0,Nevents-1,done,"Sending temperatures",
                                               QProgressEvent_None,-1,-1,-1,"Writing CCX solver input file");
        QApplication::postEvent(myProgressIndicator,e);
        QApplication::processEvents();
        QThread::msleep(500);

        if(Global::status().code==0)
        {
            cout<<"writeSolverFileClass::perform()->____process stopped____"<<endl;
            return false;
        }
    }


    //! ---------------
    //! write the STEP
    //! ---------------
    // consider using the "Number of Step" property instead of tabData->rowCount()
    bool NLgeom;
    if(mySimulationRoot->data(Qt::UserRole).value<SimulationNodeClass*>()->getType() != SimulationNodeClass::nodeType_thermalAnalysis)
        NLgeom = nodeAnalysisSettings->getPropertyValue<bool>("Large deflection");

    int NbSteps = tabData->rowCount()-1;
    for(int i=1;i<=NbSteps;i++)
    {
        if(Global::status().code==0)
        {
            cout<<"writeSolverFileClass::perform()->____process stopped____"<<endl;
            return false;
        }

        //Property::analysisType analysisType = nodeAnalysisSettings->getPropertyValue<Property::analysisType>("Analysis type");
        //Property::timeIntegration timeIntegration = nodeAnalysisSettings->getPropertyValue<Property::timeIntegration>("Static/Transient");
        Property::analysisType analysisType = tabData->dataRC(i,TABULAR_DATA_ANALYSIS_TYPE_COLUMN,Qt::EditRole).value<Property::analysisType>();
        Property::timeIntegration timeIntegration = tabData->dataRC(i,TABULAR_DATA_TIME_INTEGRATION_COLUMN,Qt::EditRole).value<Property::timeIntegration>();
        switch(analysisType)
        {
        case Property::analysisType_structural:
        {
            if(NLgeom) myInputFile<<"*STEP, NLGEOM"<<endl;
            else myInputFile<<"*STEP, NLGEOM = NO"<<endl;
            if(timeIntegration==Property::timeIntegration_steadyState) myInputFile<<"*STATIC, ";
            else myInputFile<<"*DYNAMIC, ";
        }
            break;
        case Property::analysisType_thermal:
        {
            myInputFile<<"*STEP, INC=10000"<<endl;
            if(timeIntegration==Property::timeIntegration_steadyState) myInputFile<<"*HEAT TRANSFER, STEADY STATE,";
            else myInputFile<<"*HEAT TRANSFER,";
        }
            break;
        case Property::analysisType_modal:
        {
            if(i==1)
                myInputFile<<"*STEP"<<endl;
            else
                //! Perturbation has to be active is a previous static step is present
                myInputFile<<"*STEP, PERTURBATION"<<endl;
            myInputFile<<"*FREQUENCY,";
        }
            break;
        case Property::analysisType_frequencyResponse:
        {
            myInputFile<<"*STEP, "<<endl;
            myInputFile<<"*MODAL DYNAMIC,";
        }
            break;
        case Property::analysisType_uncoupledTemperatureDisplacement:
        {
            myInputFile<<"*STEP, INC=10000"<<endl;
            if(timeIntegration==Property::timeIntegration_steadyState) myInputFile<<"*UNCOUPLED TEMPERATURE-DISPLACEMENT, STEADY STATE,";
            else myInputFile<<"*UNCOUPLED TEMPERATURE-DISPLACEMENT,";
        }
            break;
        case Property::analysisType_coupledTemperatureDisplacement:
        {
            myInputFile<<"*STEP, INC=10000"<<endl;
            if(timeIntegration==Property::timeIntegration_steadyState) myInputFile<<"*COUPLED TEMPERATURE-DISPLACEMENT, STEADY STATE,";
            else myInputFile<<"*COUPLED TEMPERATURE-DISPLACEMENT,";
        }
            break;
        }

        //! ------------
        //! solver type
        //! ------------
        Property::solverType theSolverType = tabData->dataRC(i,TABULAR_DATA_SOLVER_TYPE_COLUMN).value<Property::solverType>();
        if(theSolverType == Property::solverType_direct || theSolverType == Property::solverType_programControlled) {; }
        else myInputFile<<"SOLVER = ITERATIVE CHOLESKY, ";

        double TimeWidth = tabData->dataRC(i,1).toDouble()-tabData->dataRC(i-1,1).toDouble();
        QVector<int> N = tabData->dataRC(i,3).value<QVector<int>>();

        switch(analysisType)
        {
        case Property::analysisType_structural:
        case Property::analysisType_thermal:
        case Property::analysisType_coupledTemperatureDisplacement:
        case Property::analysisType_uncoupledTemperatureDisplacement:
        {
            int theAutoTimeStepping = N.at(3);
            if(theAutoTimeStepping==2)
            {
                double TincIni= TimeWidth/N[0];
                myInputFile<<"DIRECT"<<endl;
                myInputFile<<QString("%1").arg(TincIni).toStdString()<<", "<<QString("%1").arg(TimeWidth).toStdString()<<endl;
            }
            else if (theAutoTimeStepping==1)
            {
                double TincMax= TimeWidth/N[1];
                double TincMin= TimeWidth/N[2];
                double TincIni= TimeWidth/N[0];
                myInputFile<<endl;
                myInputFile<<QString("%1 , %2, %3, %4").arg(TincIni).arg(TimeWidth).arg(TincMin).arg(TincMax).toStdString()<<endl;
            }
            else if (theAutoTimeStepping==0)
            {
                double TincMax= TimeWidth/N[1];
                double TincMin= TimeWidth/N[2];
                double TincIni= TimeWidth/N[0];
                myInputFile<<endl;
                myInputFile<<QString("%1 , %2, %3, %4").arg(TincIni).arg(TimeWidth).arg(TincMin).arg(TincMax).toStdString()<<endl;
            }
        }
            break;
        case Property::analysisType_modal:
        {
            //! type of properties: SOLVER, STORAGE, GLOBAL and CYCMPC
            myInputFile<<endl;
            //! Number of eigenFrequency to compute
            myInputFile<<"2"<<endl;
            //myInputFile<<"1"<<endl;
        }
            break;
        case Property::analysisType_frequencyResponse:
        {
            //! type of properties: SOLVER, STORAGE, GLOBAL and CYCMPC
            myInputFile<<endl;

        }
            break;
        }
        myInputFile<<"*RESTART, WRITE, FREQUENCY=1"<<endl;

        //! -----------------------------------------------------
        //! control parameters: time incrementation - first line
        //! -----------------------------------------------------
        //! Time incrementation policy
        QVector<int> timeIncrementationParameters = tabData->dataRC(i,TABULAR_DATA_TIME_INCREMENTATION_PARAMETERS_COLUMN,Qt::EditRole).value<QVector<int>>();

        int I_0 = timeIncrementationParameters.at(0);
        int I_R = timeIncrementationParameters.at(1);
        int I_P = timeIncrementationParameters.at(2);
        int I_C = timeIncrementationParameters.at(3);
        int I_L = timeIncrementationParameters.at(5);
        int I_G = timeIncrementationParameters.at(4);
        int I_S = timeIncrementationParameters.at(6);   //! currently unused
        int I_A = timeIncrementationParameters.at(7);
        int I_J = timeIncrementationParameters.at(8);   //! currently unused
        int I_T = timeIncrementationParameters.at(9);   //! currently unused

        myInputFile<<"*CONTROLS,PARAMETERS=TIME INCREMENTATION"<<endl;
        myInputFile<<I_0<<",";
        myInputFile<<I_R<<",";
        myInputFile<<I_P<<",";
        myInputFile<<I_C<<",";
        myInputFile<<I_L<<",";
        myInputFile<<I_G<<",";
        if(I_S==-1) myInputFile<<","; else myInputFile<<I_S<<",";
        myInputFile<<I_A<<",";
        if(I_J==-1) myInputFile<<","; else myInputFile<<I_J<<",";
        if(I_T==-1) myInputFile<<","<<endl; else myInputFile<<I_T<<","<<endl;

        //! ------------------------------------------------------
        //! control parameters: time incrementation - second line
        //! ------------------------------------------------------
        QVector<double> cutBackFactors = tabData->dataRC(i,TABULAR_DATA_CUTBACK_PARAMETERS_COLUMN,Qt::EditRole).value<QVector<double>>();

        double D_f = cutBackFactors.at(0);
        double D_C = cutBackFactors.at(1);
        double D_B = cutBackFactors.at(2);
        double D_A = cutBackFactors.at(3);
        double D_S = cutBackFactors.at(4);
        double D_H = cutBackFactors.at(5);
        double D_D = cutBackFactors.at(6);
        double W_G = cutBackFactors.at(7);

        myInputFile<<D_f<<",";
        myInputFile<<D_C<<",";
        myInputFile<<D_B<<",";
        myInputFile<<D_A<<",";
        if(D_S==-1.0) myInputFile<<",";
        else myInputFile<<D_S<<",";
        if(D_H==-1.0) myInputFile<<",";
        else myInputFile<<D_H<<",";
        myInputFile<<D_D<<",";
        if(W_G==-1.0) myInputFile<<","<<endl;
        else myInputFile<<W_G<<endl;

        //! -----------------------------------------
        //! Convergence parameters: field parameters
        //! -----------------------------------------
        QVector<double> fieldParameters = tabData->dataRC(i,TABULAR_DATA_FIELD_PARAMETERS_COLUMN,Qt::EditRole).value<QVector<double>>();

        double R_alpha_n = fieldParameters.at(0);
        double C_alpha_n = fieldParameters.at(1);
        double q_alpha_0 = fieldParameters.at(2);
        double q_alpha_u = fieldParameters.at(3);
        double R_alpha_P = fieldParameters.at(4);
        double epsilon_alpha = fieldParameters.at(5);
        double C_alpha_epsilon = fieldParameters.at(6);
        double R_alpha_l = fieldParameters.at(7);

        myInputFile<<"*CONTROLS,PARAMETERS=FIELD"<<endl;
        myInputFile<<R_alpha_n<<",";
        myInputFile<<C_alpha_n<<",";
        if(q_alpha_0==0.0) myInputFile<<","; else myInputFile<<q_alpha_0<<",";
        if(q_alpha_u==0.0) myInputFile<<","; else myInputFile<<q_alpha_u<<",";
        myInputFile<<R_alpha_P<<",";
        myInputFile<<epsilon_alpha<<",";
        myInputFile<<C_alpha_epsilon<<",";
        myInputFile<<R_alpha_l<<endl;

        //! ------------
        //! line search
        //! ------------
        QVector<double> lineSearch = tabData->dataRC(i,TABULAR_DATA_LINE_SEARCH_PARAMETERS_COLUMN,Qt::EditRole).value<QVector<double>>();
        double minLineSearch = lineSearch.at(0);
        double maxLineSearch = lineSearch.at(1);

        myInputFile<<"*CONTROLS,PARAMETERS=LINE SEARCH"<<endl;
        myInputFile<<",";
        myInputFile<<maxLineSearch<<",";
        myInputFile<<minLineSearch<<",";
        myInputFile<<","<<endl;

        myInputFile<<"*CONTROLS,PARAMETERS=CONTACT"<<endl;
        myInputFile<<"0.001,";
        myInputFile<<"0.1,";
        myInputFile<<"100,";
        myInputFile<<"25,"<<endl;

        //! --------------------
        //! boundary conditions
        //! --------------------
        for(int k=1; k<mySimulationRoot->rowCount()-1; k++)
        {
            if(Global::status().code==0)
            {
                cout<<"writeSolverFileClass::perform()->____process stopped____"<<endl;
                return false;
            }

            cout<<"writeSolverFileClass::perform()->_____writing "<<mySimulationRoot->child(k,0)->data(Qt::DisplayRole).toString().toStdString()<<endl;
            QString itemName = itemNameClearSpaces(mySimulationRoot->child(k,0)->data(Qt::DisplayRole).toString());

            SimulationNodeClass *theCurNode = mySimulationRoot->child(k,0)->data(Qt::UserRole).value<SimulationNodeClass*>();
            SimulationNodeClass::nodeType theNodeType= theCurNode->getType();

            Property::SuppressionStatus theNodeSS = theCurNode->getPropertyValue<Property::SuppressionStatus>("Suppressed");

            if(theNodeSS==Property::SuppressionStatus_Active &&
                    theNodeType!=SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_CylindricalSupport &&
                    theNodeType!=SimulationNodeClass::nodeType_structuralAnalysisBoundaryContidion_FixedSupport &&
                    theNodeType!=SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_FrictionlessSupport &&
                    theNodeType!=SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_CompressionOnlySupport &&
                    analysisType!=Property::analysisType_modal
       #ifdef COSTAMP_VERSION
               && theNodeType!=SimulationNodeClass::nodeType_timeStepBuilder
       #endif
               )
            {
                QExtendedStandardItem *theCurItem = static_cast<QExtendedStandardItem*>(mySimulationRoot->child(k,0));
                QString SetName = itemName.append("_").append(QString("%1").arg(k));
                switch(theNodeType)
                {
                case SimulationNodeClass::nodeType_structuralAnalysisBoltPretension:
                {
                    QList<int> ColumnList = mainTreeTools::getColumnsToRead(theCurItem);

                    QString RN = SetName+"_RN";

                    Property::boltStatusDefinedBy theBoltDefineBy = tabData->dataRC(i,ColumnList.at(0)).value<Property::boltStatusDefinedBy>();
                    switch(theBoltDefineBy)
                    {
                    case Property::boltStatusDefinedBy_load:
                    {
                        double loadValue = tabData->dataRC(i,ColumnList.at(1)).toDouble();
                        myInputFile<<"*CLOAD"<<endl;
                        myInputFile<<RN.toStdString()<<","<<1<<","<<loadValue<<endl;
                    }
                        break;
                    case Property::boltStatusDefinedBy_adjustment:
                    {
                        double loadValue = tabData->dataRC(i,ColumnList.at(2)).toDouble();
                        myInputFile<<"*BOUNDARY"<<endl;
                        myInputFile<<RN.toStdString()<<","<<1<<","<<1<<","<<loadValue<<endl;
                    }
                        break;
                    case Property::boltStatusDefinedBy_lock:
                    {
                        double loadValue = 0.0;
                        myInputFile<<"*BOUNDARY,FIXED"<<endl;
                        myInputFile<<RN.toStdString()<<","<<1<<","<<3<<","<<loadValue<<endl;
                    }
                        break;
                    }
                }
                    break;
                case SimulationNodeClass::nodeType_structuralAnalysisThermalCondition:
                {
                    QList<int> ColumnList = mainTreeTools::getColumnsToRead(theCurItem);
                    for(int n=0;n<ColumnList.length();n++)
                    {
                        //! -----------------------------------
                        //! the loadValue is the "Temperature"
                        //! -----------------------------------
                        double loadValue = tabData->dataRC(i,ColumnList[n]).toDouble(); //if it is a temperature distribution add an exception with an include to an external file
                        myInputFile<<"*TEMPERATURE"<<endl;
                        myInputFile<<SetName.toStdString()<<", "<<loadValue<<endl;
                    }
                }
                    break;
                case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_ImportedTemperatureDistribution:
                {
                    int couplingTime = theCurNode->getPropertyValue<int>("Structural time step");
                    if(i==couplingTime)
                    {
                        //void *p =theCurNode->getPropertyValue<void*>("Imported body temperature");
                        void *p =theCurNode->getPropertyItem("Imported body temperature")->data(Qt::UserRole).value<Property>().getData().value<void*>();
                        QExtendedStandardItem* itemTemperatureDist = (QExtendedStandardItem*)(p);
                        SimulationNodeClass *nodeTemperatureDist = itemTemperatureDist->data(Qt::UserRole).value<SimulationNodeClass*>();
                        QStandardItem *itemPostObject = nodeTemperatureDist->getPropertyItem("Post object");
                        if(itemPostObject!=Q_NULLPTR)
                        {
                            QString extension=".t";
                            QString tname = itemName+QString("_%1_").arg(i)+extension;
                            QString absFileName = myFileName.split("/").last();
                            QString dirName = myFileName;
                            dirName.chop(absFileName.size());
                            tname.prepend(dirName);
                            postObject pObject = itemPostObject->data(Qt::UserRole).value<Property>().getData().value<postObject>();
                            myInputFile<<"*TEMPERATURE"<<endl;
                            this->writeTemperatureHistory(pObject, tname);
                        }
                        else
                        {
                            int button = QMessageBox::warning(0,APPNAME,"Problem in writing the imported body temperature.\n"
                                                                        "Do you wish to continue?",QMessageBox::Ok,QMessageBox::Abort);
                            if((QMessageBox::StandardButton)(button)==QMessageBox::Abort) return false;
                        }
                    }
                }
                    break;
                case SimulationNodeClass::nodeType_mapper:
                {
                    //! scan the imported body scalar under the "Mapper" item
                    for(int n=0; n<theCurItem->rowCount(); n++)
                    {
                        //! the current mapper
                        QExtendedStandardItem *itemBodyScalar = static_cast<QExtendedStandardItem*>(theCurItem->child(n,0));
                        SimulationNodeClass *ImportedBodyScalarNode = theCurItem->child(n,0)->data(Qt::UserRole).value<SimulationNodeClass*>();

                        SimulationNodeClass::nodeType ImportedBodyScalarType= ImportedBodyScalarNode->getType();
                        if(ImportedBodyScalarType!=SimulationNodeClass::nodeType_OpenFoamScalarData)
                        {
                            QString extension=".t";
                            QString index=QString::number(n);
                            QString time=QString::number(i);
                            QString tname = SetName+"_"+index+"_"+time+extension;
                            QString absFileName = myFileName.split("/").last();
                            QString dirName = myFileName;
                            dirName.chop(absFileName.size());
                            tname.prepend(dirName);

                            int stepSelection = ImportedBodyScalarNode->getPropertyValue<int>("Step selection mode");
                            postObject pObject;
                            myInputFile<<"*TEMPERATURE"<<endl;

                            switch(stepSelection)
                            {
                            case 0:     //All time steps
                            {
                                //! scan the interpolation results at different times
                                SimulationNodeClass *node = itemBodyScalar->child(i-1,0)->data(Qt::UserRole).value<SimulationNodeClass*>();
                                Property::SuppressionStatus ss = node->getPropertyValue<Property::SuppressionStatus>("Suppressed");
                                if(ss==Property::SuppressionStatus_Active)
                                {
                                    //!search the Graphic Object property and extract temperature distribution data map
                                    pObject =  node->getPropertyItem("Post object")->data(Qt::UserRole).value<Property>().getData().value<postObject>();
                                    this->writeTemperatureHistory(pObject, tname);
                                }
                            }
                                break;

                            case 1:     //First
                            {
                                if (i==1)
                                {
                                    //! scan the interpolation results at different times
                                    SimulationNodeClass *node = itemBodyScalar->child(i-1,0)->data(Qt::UserRole).value<SimulationNodeClass*>();
                                    Property::SuppressionStatus ss = node->getPropertyValue<Property::SuppressionStatus>("Suppressed");
                                    if(ss==Property::SuppressionStatus_Active)
                                    {
                                        //!search the Graphic Object property and extract temperature distribution data map
                                        pObject =  node->getPropertyItem("Post object")->data(Qt::UserRole).value<Property>().getData().value<postObject>();
                                        this->writeTemperatureHistory(pObject, tname);
                                    }
                                }
                            }
                                break;

                            case 2:     //Last
                            {
                                cout<<i<<","<<NbSteps<<endl;
                                if (i==NbSteps)
                                {
                                    //! scan the interpolation results at different times
                                    SimulationNodeClass *node = itemBodyScalar->child(i-1,0)->data(Qt::UserRole).value<SimulationNodeClass*>();
                                    Property::SuppressionStatus ss = node->getPropertyValue<Property::SuppressionStatus>("Suppressed");
                                    if(ss==Property::SuppressionStatus_Active)
                                    {
                                        //!search the Graphic Object property and extract temperature distribution data map
                                        pObject =  node->getPropertyItem("Post object")->data(Qt::UserRole).value<Property>().getData().value<postObject>();
                                        this->writeTemperatureHistory(pObject, tname);
                                    }
                                }
                            }
                                break;

                            case 3:     //By number
                            {
                                int tNumber =  ImportedBodyScalarNode->getPropertyValue<int>("Step number");
                                if (i==tNumber)
                                {
                                    //! scan the interpolation results at different times
                                    SimulationNodeClass *node = itemBodyScalar->child(i-1,0)->data(Qt::UserRole).value<SimulationNodeClass*>();
                                    Property::SuppressionStatus ss = node->getPropertyValue<Property::SuppressionStatus>("Suppressed");
                                    if(ss==Property::SuppressionStatus_Active)
                                    {
                                        //!search the Graphic Object property and extract temperature distribution data map
                                        pObject =  node->getPropertyItem("Post object")->data(Qt::UserRole).value<Property>().getData().value<postObject>();
                                        this->writeTemperatureHistory(pObject, tname);
                                    }
                                }
                            }
                                break;

                            case 4:     //Automatic Time Stepping
                            {
                                //! scan the interpolation results at different times
                                SimulationNodeClass *node = itemBodyScalar->child(i-1,0)->data(Qt::UserRole).value<SimulationNodeClass*>();
                                Property::SuppressionStatus ss = node->getPropertyValue<Property::SuppressionStatus>("Suppressed");
                                if(ss==Property::SuppressionStatus_Active)
                                {
                                    //!search the Graphic Object property and extract temperature distribution data map
                                    pObject =  node->getPropertyItem("Post object")->data(Qt::UserRole).value<Property>().getData().value<postObject>();
                                    this->writeTemperatureHistory(pObject, tname);
                                }

                            }
                                break;
                            }
                        }
                    }
                }
                    break;
                case SimulationNodeClass::nodeType_modelChange:
                {
                    int itemType = theCurNode->getPropertyValue<int>("Item type");
                    QString type;

                    QList<int> ColumnList = mainTreeTools::getColumnsToRead(theCurItem);
                    int loadValue = tabData->dataRC(i,ColumnList.at(0)).toInt();
                    if(loadValue!=0)
                    {
                        QString addORemove;
                        if(loadValue==1) addORemove = "ADD";
                        if(loadValue==-1) addORemove = "REMOVE";
                        //! find the scope based on type
                        switch(itemType)
                        {
                        case 0:
                        {
                            type="ELEMENT";
                            std::vector<GeometryTag> scope = theCurNode->getPropertyValue<std::vector<GeometryTag>>("Tags");

                            myInputFile<<"*MODEL CHANGE, TYPE = "<<type.toStdString()<<", "<<addORemove.toStdString()<<endl;
                            for(int i=0; i<scope.size();i++)
                            {
                                GeometryTag aLoc = scope.at(i);
                                int bodyIndex = aLoc.parentShapeNr;
                                //! retrieve the name of the body from the data base
                                std::string bodyName = myDB->MapOfBodyNames.value(bodyIndex).toStdString();
                                myInputFile<<"E"<<bodyName<<",";
                            }

                            myInputFile<<endl;
                        }
                            break;
                        case 1:
                        {
                            type="CONTACT PAIR";
                            void* p = theCurNode->getPropertyItem("Contact pair")->data(Qt::UserRole).value<Property>().getData().value<void*>();
                            QExtendedStandardItem *theCPItem = static_cast<QExtendedStandardItem*>(p);
                            SimulationNodeClass *nodeCP = theCPItem->data(Qt::UserRole).value<SimulationNodeClass*>();
                            QString timeTag = nodeCP->getPropertyValue<QString>("Time tag");
                            pair<QString,QString> masterSlave = contactMapName.value(timeTag);
                            myInputFile<<"*MODEL CHANGE, TYPE = "<<type.toStdString()<<", "<<addORemove.toStdString()<<endl;
                            myInputFile<<masterSlave.first.toStdString()<<","<<masterSlave.second.toStdString()<<endl;
                        }
                            break;
                        }
                    }
                }
                    break;
                case SimulationNodeClass::nodeType_thermalAnalysisTemperature:
                {
                    QList<int> columnList = mainTreeTools::getColumnsToRead(theCurItem);
                    double loadValue = tabData->dataRC(i,columnList.at(0)).toDouble();
                    if(loadValue!=0)
                    {
                        myInputFile<<"*BOUNDARY"<<endl;
                        myInputFile<<SetName.toStdString()<<", 11,11, "<<loadValue<<endl;
                    }
                }
                    break;
                case SimulationNodeClass::nodeType_thermalAnalysisConvection:
                {
                    QList<int> ColumnList = mainTreeTools::getColumnsToRead(theCurItem);                            ;
                    double loadValue = tabData->dataRC(i,ColumnList.at(0)).toDouble();
                    //double refTemperature = theCurNode->getPropertyValue<double>("Reference temperature");
                    double refTemperature = tabData->dataRC(i,ColumnList.at(1)).toDouble();
                    if(loadValue!=0)
                    this->writeFilm(loadValue,SetName,refTemperature);
                }
                    break;
                case SimulationNodeClass::nodeType_thermalAnalysisAdiabaticWall:
                {
                    double loadValue = 0;
                    this->writeDflux(loadValue,SetName);
                }
                    break;
                case SimulationNodeClass::nodeType_thermalAnalysisThermalFlow:
                case SimulationNodeClass::nodeType_thermalAnalysisThermalFlux:
                {
                    QList<int> ColumnList = mainTreeTools::getColumnsToRead(theCurItem);                            ;
                    double loadValue = tabData->dataRC(i,ColumnList.at(0)).toDouble();
                    this->writeDflux(loadValue,SetName);
                }
                    break;
                case SimulationNodeClass::nodeType_thermalAnalysisThermalPower:
                {
                    QList<int> ColumnList = mainTreeTools::getColumnsToRead(theCurItem);
                    //! -----------------------------------
                    //! the loadValue is the "Temperature"
                    //! -----------------------------------
                    double loadValue = tabData->dataRC(i,ColumnList[0]).toDouble(); //if it is a temperature distribution add an exception with an include to an external file
                    myInputFile<<"*DFLUX"<<endl;
                    myInputFile<<SetName.toStdString()<<",BF, "<<loadValue<<endl;
                }
                    break;
                default:
                {
                    QList<int> ColumnList = mainTreeTools::getColumnsToRead(theCurItem);
                    Property::defineBy theDefineBy = theCurNode->getPropertyValue<Property::defineBy>("Define by");
                    Property::ScopingMethod scopingMethod = theCurNode->getPropertyValue<Property::ScopingMethod>("Scoping method");

                    if(theDefineBy==Property::defineBy_components)
                    {
                        Property::loadDefinition loadDefinitionXcomponent,loadDefinitionYcomponent,loadDefinitionZcomponent;
                        loadDefinitionXcomponent = theCurNode->getPropertyValue<Property::loadDefinition>("X component");
                        loadDefinitionYcomponent = theCurNode->getPropertyValue<Property::loadDefinition>("Y component");
                        loadDefinitionZcomponent = theCurNode->getPropertyValue<Property::loadDefinition>("Z component");

                        QList<int> listOfLabel;
                        if(loadDefinitionXcomponent!=Property::loadDefinition_free)listOfLabel<<1;
                        if(loadDefinitionYcomponent!=Property::loadDefinition_free)listOfLabel<<2;
                        if(loadDefinitionZcomponent!=Property::loadDefinition_free)listOfLabel<<3;

                        //! -----------------------------------------------------------------------
                        //! now it is very far from this point, but "i" is the current step number
                        //! -----------------------------------------------------------------------
                        double loadValueX=0.0;
                        double loadValueY=0.0;
                        double loadValueZ=0.0;
                        bool skipX = false;
                        bool skipY = false;
                        bool skipZ = false;
                        for(int n=0;n<ColumnList.length();n++)
                        {
                            if(loadDefinitionXcomponent!=Property::loadDefinition_free && skipX==false)
                            {
                                loadValueX = tabData->dataRC(i,ColumnList.at(n)).toDouble();
                                skipX = true;
                                continue;
                            }
                            if(loadDefinitionYcomponent!=Property::loadDefinition_free && skipY==false)
                            {
                                loadValueY = tabData->dataRC(i,ColumnList.at(n)).toDouble();
                                skipY = true;
                                continue;
                            }
                            if(loadDefinitionZcomponent!=Property::loadDefinition_free && skipZ==false)
                            {
                                loadValueZ = tabData->dataRC(i,ColumnList.at(n)).toDouble();
                                skipZ = true;
                                continue;
                            }
                        }

                        //! ---------------------------------
                        //! retrieve the system of reference
                        //! ---------------------------------
                        void* p = theCurNode->getPropertyItem("Coordinate system")->data(Qt::UserRole).value<Property>().getData().value<void*>();
                        QExtendedStandardItem *theCSItem = static_cast<QExtendedStandardItem*>(p);
                        SimulationNodeClass *nodeCS = theCSItem->data(Qt::UserRole).value<SimulationNodeClass*>();

                        //! ---------------------------------------------------------------------------------
                        //! init the directional data with [1;0;0][0;1;0][0;0;1] => global coordinate system
                        //! ---------------------------------------------------------------------------------
                        QVector<double> xAxisData {1,0,0};
                        QVector<double> yAxisData {0,1,0};
                        QVector<double> zAxisData {0,0,1};

                        if(nodeCS->getType()!=SimulationNodeClass::nodeType_coordinateSystem_global)
                        {
                            xAxisData = nodeCS->getPropertyValue<QVector<double>>("X axis data");
                            yAxisData = nodeCS->getPropertyValue<QVector<double>>("Y axis data");
                            zAxisData = nodeCS->getPropertyValue<QVector<double>>("Z axis data");
                        }

                        double loadX_global = loadValueX*xAxisData.at(0)+loadValueY*yAxisData.at(0)+loadValueZ*zAxisData.at(0);
                        double loadY_global = loadValueX*xAxisData.at(1)+loadValueY*yAxisData.at(1)+loadValueZ*zAxisData.at(1);
                        double loadZ_global = loadValueX*xAxisData.at(2)+loadValueY*yAxisData.at(2)+loadValueZ*zAxisData.at(2);

                        //! -----------------------------------------------------
                        //! the component of the load evaluated in the Global CS
                        //! -----------------------------------------------------
                        switch(theNodeType)
                        {
                        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Acceleration:
                        {
                            //! -----------------------------
                            //! Acceleration
                            //! details:
                            //! treat acceleration as gravity
                            //! ------------------------------

                            //std::vector<GeometryTag> vecLoc = theCurNode->getPropertyValue<std::vector<GeometryTag>>("Tags");
                            myInputFile<<"*DLOAD"<<endl;

                            double loadValue = pow((pow(loadX_global,2)+pow(loadY_global,2)+pow(loadZ_global,2)),0.5);
                            if(loadValue!=0.0)
                            {
                                myInputFile<<"*DLOAD"<<endl;

                                double x,y,z;
                                x = loadX_global/loadValue;
                                y = loadY_global/loadValue;
                                z = loadZ_global/loadValue;
                                /*
                            for(int i=0; i<vecLoc.size();i++)
                            {
                                GeometryTag aLoc = vecLoc.at(i);
                                int bodyIndex = aLoc.parentShapeNr;

                                //! retrieve the name of the body from the data base
                                std::string bodyName = myDB->MapOfBodyNames.value(bodyIndex).toStdString();
                                myInputFile<<"E"<<bodyName<<", GRAV, "<<loadValue<<" ,"<<x<<", "<<y<<", "<<z<<endl;
                            }
                            */
                                for(int i=0; i<theGeometryRoot->rowCount();i++)
                                {
                                    std::string bodyName;
                                    QStandardItem *aGeometryItem = theGeometryRoot->child(i,0);
                                    SimulationNodeClass *aNode = aGeometryItem->data(Qt::UserRole).value<SimulationNodeClass*>();
                                    Property::SuppressionStatus aNodeSS = aNode->getPropertyValue<Property::SuppressionStatus>("Suppressed");

                                    if(aNodeSS==Property::SuppressionStatus_Active)
                                    {
                                        if(aNode->getType()==SimulationNodeClass::nodeType_pointMass)
                                        {
                                            QString bodyNameP = itemNameClearSpaces(theGeometryRoot->child(i,0)->data(Qt::DisplayRole).toString());
                                            bodyNameP.append("_").append(QString("%1").arg(i));
                                            bodyName = bodyNameP.toStdString();
                                        }
                                        else
                                        {
                                            //! retrieve the name of the body from the data base
                                            int mapIndex = aNode->getPropertyValue<int>("Map index");
                                            bodyName = myDB->MapOfBodyNames.value(mapIndex).toStdString();
                                        }
                                        myInputFile<<"E"<<bodyName<<", GRAV, "<<loadValue<<" ,"<<x<<", "<<y<<", "<<z<<endl;
                                    }
                                }
                            }
                        }
                            break;

                        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteDisplacement:
                        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Displacement:
                        {
                            if( scopingMethod == Property::ScopingMethod_RemotePoint )
                            {
                               QString front = "CM_";
                               SetName = front.append(SetName);
                            }
                            //! -------------
                            //! Displacement
                            //! -------------
                            myInputFile<<"*BOUNDARY"<<endl;
                            if(loadDefinitionXcomponent!=Property::loadDefinition_free)
                            {
                                myInputFile<<SetName.toStdString()<<", 1, 1, "<<loadX_global<<endl;
                            }
                            if(loadDefinitionYcomponent!=Property::loadDefinition_free)
                            {
                                myInputFile<<SetName.toStdString()<<", 2, 2, "<<loadY_global<<endl;
                            }
                            if(loadDefinitionZcomponent!=Property::loadDefinition_free)
                            {
                                myInputFile<<SetName.toStdString()<<", 3, 3, "<<loadZ_global<<endl;
                            }
                        }
                            break;

                        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Force:
                        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteForce:
                        {
                            //! ------
                            //! Force
                            //! ------
                            myInputFile<<"*CLOAD"<<endl;
                            myInputFile<<(QString("CM_")+SetName).toStdString()<<", 1, "<<loadX_global<<endl;
                            myInputFile<<(QString("CM_")+SetName).toStdString()<<", 2, "<<loadY_global<<endl;
                            myInputFile<<(QString("CM_")+SetName).toStdString()<<", 3, "<<loadZ_global<<endl;
                        }
                            break;

                        case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Moment:
                        {
                            //! -------
                            //! Moment
                            //! -------
                            myInputFile<<"*CLOAD"<<endl;
                            myInputFile<<(QString("CM_")+SetName).toStdString()<<", 4, "<<loadX_global<<endl;
                            myInputFile<<"*CLOAD"<<endl;
                            myInputFile<<(QString("CM_")+SetName).toStdString()<<", 5, "<<loadY_global<<endl;
                            myInputFile<<"*CLOAD"<<endl;
                            myInputFile<<(QString("CM_")+SetName).toStdString()<<", 6, "<<loadZ_global<<endl;
                            cout<<" - writing moment (Mx, My, Mz) = ("<<loadX_global<<", "<<loadY_global<<", "<<loadZ_global<<")"<<endl;
                        }
                            break;
                        }
                    }
                    if(theDefineBy==Property::defineBy_normal)
                    {
                        for(int p=0;p<ColumnList.length();p++)
                        {
                            double loadValue = tabData->dataRC(i,ColumnList.at(p)).toDouble();
                            QString aName = SetName;
                            aName.append(QString("_%1").arg(i));
                            switch(theNodeType)
                            {
                            //! ---------
                            //! Pressure
                            //! ---------
                            case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Pressure:
                                this->writeDload(loadValue,aName);
                                break;
                            case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Displacement:
                            case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteDisplacement:
                            case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Force:
                            case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteForce:
                            case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Moment:
                            case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Acceleration:
                            {
                                //! do nothing: here only for documentation
                            }
                                break;
                            }
                        }
                    }
                    if(theDefineBy==Property::defineBy_vector)
                    {
                        //! -------------------------------------------------------------
                        //! in this case, by definition the loadValue is the "Amplitude"
                        //! -------------------------------------------------------------
                        QVector<double> dirData1 = theCurNode->getPropertyValue<QVector<double>>("Direction");
                        double a0 = dirData1.at(3);
                        double a1 = dirData1.at(4);
                        double a2 = dirData1.at(5);

                        QVector<double> dirData;
                        dirData.push_back(a0);
                        dirData.push_back(a1);
                        dirData.push_back(a2);

                        int n;
                        for(n=0;n<ColumnList.length();n++)
                        {
                            double loadValue = tabData->dataRC(i,ColumnList[n]).toDouble();
                            //! --------------------
                            //! find the components
                            //! --------------------
                            QVector<double> comps = vectorTool::getComponents(loadValue, dirData);
                            double Xcomp = comps.at(0);
                            double Ycomp = comps.at(1);
                            double Zcomp = comps.at(2);

                            switch(theNodeType)
                            {
                            case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Displacement:
                            case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteDisplacement:
                            {
                                if(scopingMethod == Property::ScopingMethod_RemotePoint)
                                    SetName = ("CM_")+SetName;

                                myInputFile<<"*BOUNDARY"<<endl;
                                myInputFile<<SetName.toStdString()<<", 1, 1, "<<Xcomp<<endl;
                                myInputFile<<SetName.toStdString()<<", 2, 2, "<<Ycomp<<endl;
                                myInputFile<<SetName.toStdString()<<", 3, 3, "<<Zcomp<<endl;
                            }
                                break;
                            case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Force:
                            case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteForce:
                            {
#ifdef COSTAMP_VERSION
                                if(loadValue ==0.0)
                                {
                                    myInputFile<<"*BOUNDARY"<<endl;
                                    myInputFile<<(QString("CM_")+SetName).toStdString()<<", 1,3,0 "<<endl;
                                }
                                else
                                {
#endif
                                    myInputFile<<"*CLOAD"<<endl;
                                    myInputFile<<(QString("CM_")+SetName).toStdString()<<", 1, "<<Xcomp<<endl;
                                    myInputFile<<(QString("CM_")+SetName).toStdString()<<", 2, "<<Ycomp<<endl;
                                    myInputFile<<(QString("CM_")+SetName).toStdString()<<", 3, "<<Zcomp<<endl;
#ifdef COSTAMP_VERSION
                                }
#endif
                            }
                                break;
                            case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Moment:
                            {
                                myInputFile<<"*CLOAD"<<endl;
                                myInputFile<<(QString("CM_")+SetName).toStdString()<<", 4, "<<Xcomp<<endl;
                                myInputFile<<(QString("CM_")+SetName).toStdString()<<", 5, "<<Ycomp<<endl;
                                myInputFile<<(QString("CM_")+SetName).toStdString()<<", 6, "<<Zcomp<<endl;
                            }
                                break;
                            case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Acceleration:
                            {
                                //std::vector<GeometryTag> vecLoc = theCurNode->getPropertyValue<std::vector<GeometryTag>>("Tags");
                                myInputFile<<"*DLOAD"<<endl;
                                /*
                                for(int i=0; i<vecLoc.size();i++)
                                {
                                    GeometryTag aLoc = vecLoc.at(i);
                                    int bodyIndex = aLoc.parentShapeNr;

                                    //! retrieve the name of the body from the data base
                                    std::string bodyName = myDB->MapOfBodyNames.value(bodyIndex).toStdString();
                                    myInputFile<<"E"<<bodyName<<", GRAV, "<<loadValue<<" ,"<<Xcomp<<", "<<Ycomp<<", "<<Zcomp<<endl;
                                }
                                */
                                for(int i=0; i<theGeometryRoot->rowCount();i++)
                                {
                                    std::string bodyName;
                                    QStandardItem *aGeometryItem = theGeometryRoot->child(i,0);
                                    SimulationNodeClass *aNode = aGeometryItem->data(Qt::UserRole).value<SimulationNodeClass*>();
                                    if(aNode->getType()==SimulationNodeClass::nodeType_pointMass)
                                    {
                                        QString bodyNameP = itemNameClearSpaces(theGeometryRoot->child(i,0)->data(Qt::DisplayRole).toString());
                                        bodyNameP.append("_").append(QString("%1").arg(i));
                                        bodyName = bodyNameP.toStdString();
                                    }
                                    else
                                    {
                                        int mapIndex = aNode->getPropertyValue<int>("Map index");
                                        //! retrieve the name of the body from the data base
                                        bodyName = myDB->MapOfBodyNames.value(mapIndex).toStdString();
                                    }
                                    myInputFile<<"E"<<bodyName<<", GRAV, "<<loadValue<<" ,"<<dirData.at(0)<<", "<<dirData.at(1)<<", "<<dirData.at(2)<<endl;
                                }
                            }
                                break;
                            }
                        }
                    }
                }
                    break;
                }
            }
        }

        //! ------------------------------------------------------------------------
        //! output file settings:
        //!
        //! - nodal: Always stored: nodal displacements, temperature.
        //!          Reaction forces if required
        //!
        //! - element: Always stored: energy, error, herror.
        //!            Stress tensor (average values) if required
        //!
        //!            Equivalent creep strain && total lagrangian/eulerian strain
        //!            && mechanical lagrangian/eulerian strain if required
        //!
        //! - contact: relative displacement && contact stress && contact energy &&
        //!            PCON if required
        //! -------------------------------------------------------------------------
        QVector<int> outputControls = tabData->dataRC(i,TABULAR_DATA_OUPUT_CONTROLS_COLUMN,Qt::EditRole).value<QVector<int>>();
        int reactionForcesFlag = outputControls.at(2);
        int stressFlag = outputControls.at(0);
        int strainFlag = outputControls.at(1);
        int contactDataFlag = outputControls.at(3);
        int error_energyFlag = 0;

        //! --------------------------------------
        //! retrieve the results frequency output
        //! --------------------------------------
        QVector<int> frequencyOutput = tabData->dataRC(i,TABULAR_DATA_STORE_RESULTS_AT_COLUMN,Qt::EditRole).value<QVector<int>>();
        int frequency = frequencyOutput.at(1);

        //! -----------------------------------------------
        //! Node results - displacement are always written
        //! -----------------------------------------------
        myInputFile<<"*NODE FILE";
        myInputFile<<", LAST ITERATIONS";
        myInputFile<<", FREQUENCY ="<<frequency<<endl;
        myInputFile<<"U,NT,PU";
        if(reactionForcesFlag==1) myInputFile<<",RF"<<endl;
        else myInputFile<<endl;

        //! --------------------------------------------------------------
        //! element results: energy, error (structural), herror (thermal)
        //! aren't written insert corresponding Flag in interface
        //! --------------------------------------------------------------
        myInputFile<<"*EL FILE";
        myInputFile<<", LAST ITERATIONS";
        myInputFile<<", FREQUENCY ="<<frequency<<endl;
        if(stressFlag==1) myInputFile<<"S,";
        if(error_energyFlag==1) myInputFile<<"ENER,ERR,HER";
        if(strainFlag==1) myInputFile<<",CEEQ,E,ME,PEEQ";
        myInputFile<<",HFL"; //ThermalFlux
        myInputFile<<endl;
        if(contactDataFlag==1)
        {
            myInputFile<<"*CONTACT FILE";
            myInputFile<<", FREQUENCY ="<<frequency<<endl;
            myInputFile<<"CDIS,CSTR,CELS,PCON"<<endl;
        }

        //! -------------------------
        //! end output file settings
        //! -------------------------
        myInputFile<<"*END STEP"<<endl;
    }

    cout<<"writeSolverFileClass::perform()->____INPUT FILE SUCCESSFULLY WRITTEN____"<<endl;

    //! --------------------
    //! update the progress
    //! --------------------
    if(myProgressIndicator!=Q_NULLPTR)
    {
        done++;
        QProgressEvent *e = new QProgressEvent(QProgressEvent_Reset,-1,-1,-1,"",
                                               QProgressEvent_Reset,-1,-1,-1,"");
        QApplication::postEvent(myProgressIndicator,e);
        QApplication::processEvents();
        QThread::msleep(500);
    }

    myInputFile.close();
    return true;
}

//! ------------------------------
//! function: writeElementSurface
//! details:  write SURFACE block
//! ------------------------------
void writeSolverFileClass::writeElementSurface(QString SetName,
                                               const IndexedMapOfMeshDataSources &anIndexedMapOfFaceMeshDS)
{
    QString extension=".surf";

    //! this is the absolute path on the disk
    QString eSurfName=SetName+extension;
    QString absFileName = myFileName.split("/").last();
    QString dirName = myFileName;
    dirName.chop(absFileName.size());
    eSurfName.prepend(dirName);

    ofstream myESurf;
    myESurf.setf(ios::scientific);
    myESurf.precision(EXPFORMAT_PRECISION);
    myESurf.open(eSurfName.toStdString());
    std::vector<int> theElementIDs;
    std::vector<int> theFaceNumbers;

    //! -----------------------------------------------------------------------
    //! searching the 3D elements containing, as faces, the selected triangles
    //! -----------------------------------------------------------------------
    this->createElementSurface(theElementIDs,theFaceNumbers,anIndexedMapOfFaceMeshDS);

    myESurf<<"*SURFACE, NAME = "<<SetName.toStdString()<<endl;
    for(int k=0; k<theElementIDs.size(); k++)
    {
        myESurf<<theElementIDs[k]<<", S"<<theFaceNumbers[k]<<endl;
    }
    myESurf.close();
    myInputFile<<"*INCLUDE, INPUT="<<eSurfName.split("/").last().toStdString()<<endl;
}

//! ----------------------------------
//! function: writeElementSurfaceBolt
//! details:
//! ----------------------------------
void writeSolverFileClass::writeElementSurfaceBolt(QString SetName, int boltBodyIndex, double a, double b, double c, double d)
{
    //! --------------------------------------
    //! this is the absolute path on the disk
    //! --------------------------------------
    QString eSurfName=SetName+".surf";
    QString absFileName = myFileName.split("/").last();
    QString dirName = myFileName;
    dirName.chop(absFileName.size());
    eSurfName.prepend(dirName);

    ofstream myESurf;
    myESurf.setf(ios::scientific);
    myESurf.precision(EXPFORMAT_PRECISION);
    myESurf.open(eSurfName.toStdString());

    //! ---------------------------------
    //! write the header for the SURFACE
    //! ---------------------------------
    myESurf<<"*SURFACE, NAME = "<<SetName.toStdString()<<endl;

    //! --------------------------
    //! create the element offset
    //! --------------------------
    int offset = 0;
    int K = 0;
    for(int m=1; m<boltBodyIndex; m++)
    {
        if (myDB->ArrayOfMeshDS.value(m).IsNull()) K=0;
        else K=myDB->ArrayOfMeshDS.value(m)->GetAllElements().Extent();
        offset = offset+K;
    }

    //! ---------------------
    //! create the bolt tool
    //! ---------------------
    boltTool abt(occHandle(Ng_MeshVS_DataSource3D)::DownCast(myDB->ArrayOfMeshDS.value(boltBodyIndex)));

    occHandle(MeshVS_DataSource) slicedMeshDS;
    std::vector<std::pair<int,int>> vecCCXFaceDefs;
    bool isDone = abt.sliceMeshWithPlane(a,b,c,d,slicedMeshDS,vecCCXFaceDefs);
    if(isDone)
    {
        for(int n=0; n<vecCCXFaceDefs.size(); n++)
        {
            const std::pair<int,int> &aCCXFaceDef = vecCCXFaceDefs[n];
            int globalElementID = aCCXFaceDef.first;
            int CCXnumber = aCCXFaceDef.second;
            myESurf<<globalElementID+offset<<", S"<<CCXnumber<<endl;
        }
    }

    myESurf.close();
    myInputFile<<"*INCLUDE, INPUT="<<eSurfName.split("/").last().toStdString()<<endl;
}

//! -----------------------------------
//! function: write nodes and elements
//! details:
//! -----------------------------------
void writeSolverFileClass::writeNodesAndElements(QString aName,QMap<int,QList<int>> &nodeListByBody)
{
    QString MeshName=aName.split("/").last();
    MeshName=MeshName.split(".").first();
    MeshName=MeshName.append(".msh");

    //! this is the absolute path on the disk
    QString absFileName = myFileName.split("/").last();
    QString dirName = myFileName;
    dirName.chop(absFileName.size());
    MeshName.prepend(dirName);

    myMesh.setf(ios::scientific);
    myMesh.precision(EXPFORMAT_PRECISION);
    myMesh.open(MeshName.toStdString());

    //! number of bodies: main shapes
    int Nb = myDB->bodyMap.size();

    //! -------------------------------------------------
    //! write the NODE section - use 3D mesh datasources
    //! -------------------------------------------------

    //! increment in the number of the node
    int anIncrement = 0;
    for(int bodyIndex = 1; bodyIndex<=Nb; bodyIndex++)
    {
        const occHandle(MeshVS_DataSource) &aMeshVS_DataSource =  myDB->ArrayOfMeshDS.value(bodyIndex);
        if(!aMeshVS_DataSource.IsNull())
        {
            QList<int> listOfNodes;

            //! retrieve and write the nodes
            const TColStd_PackedMapOfInteger& aNodes = aMeshVS_DataSource->GetAllNodes();
            if(aNodes.Extent()!=0)
            {
                //! retrieve the name of the body from the data base
                std::string bodyName = myDB->MapOfBodyNames.value(bodyIndex).toStdString();
                myMesh<<"*NODE, NSET= N"<<bodyName<<endl;

                double buf[3];
                TColStd_Array1OfReal coords(*buf,1,3);
                int nbNodes;
                MeshVS_EntityType aType;

                //! ----------------------------
                //! write the nodes coordinates
                //! ----------------------------
                for (TColStd_MapIteratorOfPackedMapOfInteger anIter(aNodes); anIter.More();anIter.Next())
                {
                    int globalNodeID = anIter.Key();
                    if (!aMeshVS_DataSource->GetGeom(globalNodeID,false,coords,nbNodes,aType)) continue;
                    myMesh<<globalNodeID+anIncrement<<","<<coords(1)<<","<<coords(2)<<","<<coords(3)<<endl;
                    listOfNodes<<globalNodeID+anIncrement;
                }
                nodeListByBody.insert(bodyIndex,listOfNodes);
            }
            anIncrement = anIncrement+aNodes.Extent();
            totalNumberOfNodes=anIncrement;
        }
    }

    //! ---------------------------------------------
    //! write the ELEMENT section - all the elements
    //! ---------------------------------------------
    int anElementIncrement = 0;
    int aNodeIncrement = 0;

    //! --------------------------------
    //! retrieve the integration scheme
    //! --------------------------------
    QExtendedStandardItem *theGeometryRoot=static_cast<QExtendedStandardItem*>(this->getTreeItem(SimulationNodeClass::nodeType_geometry));
    Property::elementControl theElementControl =  theGeometryRoot->data(Qt::UserRole).value<SimulationNodeClass*>()
            ->getPropertyValue<Property::elementControl>("Element control");

    //! if flag==true perform the check of the "Integration scheme" for each body
    bool flag = theElementControl==Property::elementControl_programControlled? false:true;

    for(int bodyIndex=1; bodyIndex<=Nb; bodyIndex++)
    {
        //! retrieve the integration scheme
        Property::integrationScheme theIntegrationScheme = Property::integrationScheme_full;
        if(flag==true)
        {
            const TopoDS_Shape &theBody = myDB->bodyMap.value(bodyIndex);
            QExtendedStandardItem *theGeometryItem = this->ItemFromScope(theBody);
            SimulationNodeClass *theNode = theGeometryItem->data(Qt::UserRole).value<SimulationNodeClass*>();
            theIntegrationScheme = theNode->getPropertyValue<Property::integrationScheme>("Integration scheme");
        }

        int Ntet4 = 0;
        int Ntet10 = 0;
        int Nhexa8 = 0;
        int Nhexa20 = 0;
        int Npyr5 = 0;
        int Npyr13 = 0;
        int Nprism6 = 0;
        int Nprism15 = 0;

        //! -----------------------------------------------------------------------
        //! use the 3D mesh datasources (Ng_MeshVS_DataSource3D), as for the nodes
        //! -----------------------------------------------------------------------
        const occHandle(MeshVS_DataSource) &aMeshDS =  myDB->ArrayOfMeshDS.value(bodyIndex);
        if(!aMeshDS.IsNull())
        {
            TColStd_PackedMapOfInteger mapOfElements = aMeshDS->GetAllElements();
            TColStd_MapIteratorOfPackedMapOfInteger anIter(mapOfElements);
            for(;anIter.More();anIter.Next())
            {
                int globalElementID = anIter.Key();
                int NbNodes, aNodesBuf[20];
                TColStd_Array1OfInteger nodeIDs(*aNodesBuf,1,20);
                aMeshDS->GetNodesByElement(globalElementID,nodeIDs,NbNodes);

                //! get the element type
                switch(NbNodes)
                {
                case 4:     Ntet4++; break;
                case 10:    Ntet10++; break;
                case 8:     Nhexa8++; break;
                case 20:    Nhexa20++; break;
                case 5:     Npyr5++; break;
                case 13:    Npyr13++; break;
                case 6:     Nprism6++; break;
                case 15:    Nprism15++; break;
                }
            }

            //! ---------------------------------------------------------------------
            //! allocate the arrays for the element grouping
            //! array<element-type>[N] = [elementID, nodeID1, nodeID2, ..., nodeIDn]
            //! ---------------------------------------------------------------------
            int **arrayTet4 = new int*[Ntet4];
            for(int i=0; i<Ntet4; i++)arrayTet4[i] = new int[5];

            int **arrayTet10 = new int*[Ntet10];
            for(int i=0; i<Ntet10; i++)arrayTet10[i] = new int[11];

            int **arrayHexa8 = new int*[Nhexa8];
            for(int i=0; i<Nhexa8; i++)arrayHexa8[i] = new int[9];

            int **arrayHexa20 = new int*[Nhexa20];
            for(int i=0; i<Nhexa20; i++)arrayHexa20[i] = new int[21];

            int **arrayPyr5 = new int*[Npyr5];
            for(int i=0; i<Npyr5; i++)arrayPyr5[i] = new int[7];

            int **arrayPyr13 = new int*[Npyr13];
            for(int i=0; i<Npyr13; i++)arrayPyr13[i] = new int[14];

            int **arrayPrism6 = new int*[Nprism6];
            for(int i=0; i<Nprism6; i++)arrayPrism6[i] = new int[7];

            int **arrayPrism15 = new int*[Nprism15];
            for(int i=0; i<Nprism15; i++)arrayPrism15[i] = new int[16];

            anIter.Initialize(mapOfElements);
            int ntet4 = 0;
            int ntet10 = 0;
            int nprism6 = 0;
            int npyr5 = 0;
            int nhexa8 = 0;

            for(int j=0;anIter.More();anIter.Next(),j++)
            {
                if(j==mapOfElements.Extent()) break;
                int aKey = anIter.Key();
                int NbNodes, aNodesBuf[20];                         //! a 20 points buffer
                TColStd_Array1OfInteger nodeIDs(*aNodesBuf,1,20);   //! this holds up to second order hexa
                aMeshDS->GetNodesByElement(aKey,nodeIDs,NbNodes);

                //! get the element type using the number of nodes
                switch(NbNodes)
                {
                case 4:
                {
                    //! a TET4 element has been found
                    arrayTet4[j-nprism6-npyr5][0]=aKey+anElementIncrement;
                    for(int i=1; i<=4; i++) arrayTet4[j-nprism6-npyr5-nhexa8-ntet10][i]=nodeIDs.Value(i)+aNodeIncrement;
                    ntet4++;
                }
                    break;

                case 10:
                {
                    //! a TET10 element has been found
                    arrayTet10[j][0]=aKey+anElementIncrement;
                    //for(int i=1; i<=10; i++) arrayTet10[j][i]=nodeIDs.Value(i)+aNodeIncrement;
                    for(int i=1; i<=10; i++) arrayTet10[j-ntet4-npyr5-nprism6-nhexa8][i]=nodeIDs.Value(i)+aNodeIncrement;
                    ntet10++;
                }
                    break;

                case 5:
                {
                    //! a PYR5 element has been found
                    arrayPyr5[j-ntet4-nprism6][0]=aKey+anElementIncrement;
                    //for(int i=1; i<=6; i++) arrayPyr5[j-ntet4-nprism6-nhexa8][i]=nodeIDs.Value(i)+aNodeIncrement;
                    for(int i=1; i<=6; i++) arrayPyr5[j-ntet4-nprism6-nhexa8-ntet10][i]=nodeIDs.Value(i)+aNodeIncrement;
                    npyr5++;
                }
                    break;

                case 13:
                {
                    //! a PYR13 element has been found
                    arrayPyr13[j][0]=aKey+anElementIncrement;
                    for(int i=1; i<=13; i++) arrayPyr13[j][i]=nodeIDs.Value(i)+aNodeIncrement;
                }
                    break;

                case 8:
                {
                    //! a HEXA8 element has been found
                    arrayHexa8[j][0]=aKey+anElementIncrement;
                    //for(int i=1; i<=8; i++) arrayHexa8[j-ntet4-nprism6-npyr5][i]=nodeIDs.Value(i)+aNodeIncrement;
                    for(int i=1; i<=8; i++) arrayHexa8[j-ntet4-npyr5-nprism6-ntet10][i]=nodeIDs.Value(i)+aNodeIncrement;
                    nhexa8++;
                }
                    break;

                case 20:
                {
                    //! a HEXA20 element has been found
                    arrayHexa20[j][0]=aKey+anElementIncrement;
                    for(int i=1; i<=20; i++) arrayHexa20[j][i]=nodeIDs.Value(i)+aNodeIncrement;
                }
                    break;

                case 6:
                {
                    //! a PRISM6 element has been found
                    arrayPrism6[j-ntet4-npyr5][0]=aKey+anElementIncrement;
                    //for(int i=1; i<=6; i++) arrayPrism6[j-ntet4-npyr5-nhexa8][i]=nodeIDs.Value(i)+aNodeIncrement;
                    for(int i=1; i<=6; i++) arrayPrism6[j-ntet4-npyr5-nhexa8-ntet10][i]=nodeIDs.Value(i)+aNodeIncrement;
                    nprism6++;
                }
                    break;

                case 15:
                {
                    //! a PRISM15 element has been found
                    arrayPrism15[j][0]=aKey+anElementIncrement;
                    for(int i=1; i<=15; i++) arrayPrism15[j][i]=nodeIDs.Value(i)+aNodeIncrement;
                }
                    break;
                }
            }
            //! -------------------------------------------------
            //! retrieve the name of the body from the data base
            //! -------------------------------------------------
            std::string bodyName = myDB->MapOfBodyNames.value(bodyIndex).toStdString();
            //! ------------------------------------------------------------------------------------
            //! When writing element, add them to a multi-map
            //! The first index of the multimap (the Key) is the element number (elementID)
            //! The second element of the multimap (the value) is a triangular (quadrangular)
            //! face (see the struct definition at the beginning). This allow to access the
            //! element number (element ID) using a triad {V1,V2,V3} (or a quaternion {V1,V2,V3,V4})
            //! -------------------------------------------------------------------------------------

            //! ---------------------
            //! write the tet4 group
            //! ---------------------
            if(Ntet4>0)
            {
                switch(theIntegrationScheme)
                {
                case Property::integrationScheme_full: myMesh<<"*ELEMENT,TYPE = C3D4"<<endl; break;
                case Property::integrationScheme_reduced: myMesh<<"*ELEMENT,TYPE = C3R4"<<endl; break;
                }

                for(int k=0; k<Ntet4; k++)
                {
                    myMesh<<arrayTet4[k][0]<<","<<arrayTet4[k][1]<<
                                             ","<<arrayTet4[k][2]<<
                                             ","<<arrayTet4[k][4]<<
                                             ","<<arrayTet4[k][3]<<endl;
                }
            }
            for(int i=0; i<Ntet4; i++) delete[] arrayTet4[i]; delete [] arrayTet4;

            //! ----------------------
            //! write the tet10 group
            //! ----------------------
            if(Ntet10>0)
            {
                switch(theIntegrationScheme)
                {
                case Property::integrationScheme_full: myMesh<<"*ELEMENT,TYPE = C3D10"<<endl; break;
                case Property::integrationScheme_reduced: myMesh<<"*ELEMENT,TYPE = C3R10"<<endl; break;
                }

                for(int k=0; k<Ntet10; k++)
                {
                    myMesh<<arrayTet10[k][0]<<","<<arrayTet10[k][1]<<
                                              ","<<arrayTet10[k][2]<<
                                              ","<<arrayTet10[k][4]<<
                                              ","<<arrayTet10[k][3]<<
                                              ","<<arrayTet10[k][5]<<
                                              ","<<arrayTet10[k][9]<<
                                              ","<<arrayTet10[k][10]<<
                                              ","<<arrayTet10[k][7]<<
                                              ","<<arrayTet10[k][6]<<
                                              ","<<arrayTet10[k][8]<<endl;
                }
            }
            for(int i=0; i<Ntet10; i++) delete[] arrayTet10[i]; delete[] arrayTet10;

            //! -----------------------
            //! write the prism6 group
            //! -----------------------
            if(Nprism6>0)
            {
                switch(theIntegrationScheme)
                {
                case Property::integrationScheme_full:
                case Property::integrationScheme_reduced:
                    myMesh<<"*ELEMENT,TYPE = C3D6"<<endl; break;
                }

                for(int k=0; k<Nprism6; k++)
                {
                    myMesh<<arrayPrism6[k][0]<<","<<arrayPrism6[k][1]<<
                                               ","<<arrayPrism6[k][2]<<
                                               ","<<arrayPrism6[k][3]<<
                                               ","<<arrayPrism6[k][4]<<
                                               ","<<arrayPrism6[k][5]<<
                                               ","<<arrayPrism6[k][6]<<endl;

                }
            }
            for(int i=0; i<Nprism6; i++) delete[] arrayPrism6[i]; delete[] arrayPrism6;

            //! ------------------------
            //! write the prism15 group
            //! ------------------------
            if(Nprism15>0)
            {
                switch(theIntegrationScheme)
                {
                case Property::integrationScheme_full:
                case Property::integrationScheme_reduced:
                    myMesh<<"*ELEMENT,TYPE = C3D15"<<endl; break;
                }

                for(int k=0; k<Nprism15; k++)
                {
                    myMesh<<arrayPrism15[k][0]<<
                                                 ","<<arrayPrism15[k][1]<<
                                                 ","<<arrayPrism15[k][2]<<
                                                 ","<<arrayPrism15[k][3]<<
                                                 ","<<arrayPrism15[k][4]<<
                                                 ","<<arrayPrism15[k][5]<<
                                                 ","<<arrayPrism15[k][6]<<
                                                 ","<<arrayPrism15[k][7]<<
                                                 ","<<arrayPrism15[k][8]<<
                                                 ","<<arrayPrism15[k][9]<<
                                                 ","<<arrayPrism15[k][10]<<
                                                 ","<<arrayPrism15[k][11]<<
                                                 ","<<arrayPrism15[k][12]<<
                                                 ","<<arrayPrism15[k][13]<<
                                                 ","<<arrayPrism15[k][14]<<
                                                 ","<<arrayPrism15[k][15]<<endl;

                    int offset = 0;
                    int J=0;

                    for(int j=1; j<bodyIndex; j++)
                    {
                        if (myDB->ArrayOfMeshDS.value(j).IsNull()) J=0;
                        else J=myDB->ArrayOfMeshDS.value(j)->GetAllElements().Extent();
                        offset = offset+J;
                    }
                }
            }
            for(int i=0; i<Nprism15; i++) delete[] arrayPrism15[i]; delete[] arrayPrism15;

            //! ----------------------------------------------------------
            //! write the pyr5 group
            //! details: in calculix they are a degenerate form of a C3D6
            //! ----------------------------------------------------------
            if(Npyr5>0)
            {
                switch(theIntegrationScheme)
                {
                case Property::integrationScheme_full: myMesh<<"*ELEMENT,TYPE = C3D8"<<endl; break;
                case Property::integrationScheme_reduced: myMesh<<"*ELEMENT,TYPE = C3D8R"<<endl; break;
                }

                for(int k=0; k<Npyr5; k++)
                {
                    myMesh<<arrayPyr5[k][0]<<
                                              ","<<arrayPyr5[k][3]<<
                                              ","<<arrayPyr5[k][2]<<
                                              ","<<arrayPyr5[k][5]<<
                                              ","<<arrayPyr5[k][4]<<
                                              ","<<arrayPyr5[k][1]<<
                                              ","<<arrayPyr5[k][1]<<
                                              ","<<arrayPyr5[k][1]<<
                                              ","<<arrayPyr5[k][1]<<endl;

                }
            }
            for(int i=0; i<Npyr5; i++) delete[] arrayPyr5[i]; delete[] arrayPyr5;

            //! ----------------------
            //! write the hexa8 group
            //! ----------------------
            if(Nhexa8>0)
            {
                switch(theIntegrationScheme)
                {
                case Property::integrationScheme_full: myMesh<<"*ELEMENT,TYPE = C3D8"<<endl; break;
                case Property::integrationScheme_reduced: myMesh<<"*ELEMENT,TYPE = C3R8"<<endl; break;
                }

                for(int k=0; k<Nhexa8; k++)
                {
                    myMesh<<arrayHexa8[k][0]<<","<<arrayHexa8[k][4]<<
                                              ","<<arrayHexa8[k][3]<<
                                              ","<<arrayHexa8[k][2]<<
                                              ","<<arrayHexa8[k][1]<<
                                              ","<<arrayHexa8[k][8]<<
                                              ","<<arrayHexa8[k][7]<<
                                              ","<<arrayHexa8[k][6]<<
                                              ","<<arrayHexa8[k][5]<<
                                              endl;
                }
            }
            for(int i=0; i<Nhexa8; i++) delete[] arrayHexa8[i]; delete [] arrayHexa8;

            myMesh<<"*ELSET, ELSET = E"<<bodyName<<", GENERATE"<<endl;
            myMesh<<anElementIncrement+1<<", ";
            //! ---------------------------------------------
            //! increment the element key (assembly support)
            //! ---------------------------------------------
            anElementIncrement = anElementIncrement + Ntet4 + Ntet10 + Nhexa8 + Nhexa20 + Npyr5 + Npyr13 + Nprism6 + Nprism15;
            myMesh<<anElementIncrement<<endl;
            totalNumberOfElements+=anElementIncrement;

            //! ---------------------------------------------
            //! increment the node number (assembly support)
            //! ---------------------------------------------
            aNodeIncrement = aNodeIncrement + aMeshDS->GetAllNodes().Extent();
        }
    }
    ++totalNumberOfElements;
    myMesh<<"*NSET, NSET=NALL, GENERATE"<<endl;
    myMesh<<"1,"<<totalNumberOfNodes<<endl;
    myMesh.close();
    myInputFile<<"*INCLUDE, INPUT="<<MeshName.split("/").last().toStdString()<<endl;
}

//! -------------------------------
//! function: writeNodalSet
//! details:  write the NSET block
//! -------------------------------
void writeSolverFileClass::writeNodalSet(QString SetName,
                                         const IndexedMapOfMeshDataSources &anIndexedMapOfFaceMeshDS)
{
    QString extension=".nam";
    //! this is the absolute path on the disk
    QString nsetName=SetName+extension;
    QString absFileName = myFileName.split("/").last();
    QString dirName = myFileName;
    dirName.chop(absFileName.size());
    nsetName.prepend(dirName);

    ofstream myNSet;
    myNSet.setf(ios::scientific);
    myNSet.precision(EXPFORMAT_PRECISION);
    myNSet.open(nsetName.toStdString());
    myNSet<<"*NSET, NSET = "<<SetName.toStdString()<<endl;

    std::vector<int> theNodeIDs;
    this->createNodalSet(anIndexedMapOfFaceMeshDS,theNodeIDs);

    for(int i=0, block =1; i<theNodeIDs.size(); i++, block++)
    {
        int nodeID = theNodeIDs[i];
        myNSet<<nodeID<<",";
        //! a maximum of 12 entries (nodes) per line (16 allowed)
        if(block%12==0)myNSet<<endl;
    }
    myNSet<<endl;
    myNSet.close();
    myInputFile<<"*INCLUDE, INPUT="<<nsetName.split("/").last().toStdString()<<endl;
}

//! --------------------------
//! function: writeElementSet
//! details:
//! --------------------------
void writeSolverFileClass::writeElementSet(std::vector<GeometryTag> vecLoc,QList<QString> &bodyNameList)
{
    for(int i=0; i<vecLoc.size();i++)
    {
        GeometryTag aLoc = vecLoc.at(i);
        int bodyIndex = aLoc.parentShapeNr;

        //! retrieve the name of the body from the data base
        //! Example: *DLOAD Eall,GRAV,9810.,0.,0.,-1.
        QString bodyName = myDB->MapOfBodyNames.value(bodyIndex);
        bodyNameList<<bodyName;
    }
}
//! ----------------------------------
//! function: clock
//! details:  for diagnostic purposes
//! ----------------------------------
void writeSolverFileClass::clock()
{
    QTime time = QTime::currentTime();
    QString text = time.toString("mm:ss:zzz");
    cout<<text.toStdString()<<endl;
}

//!------------------------------------------------
//! function: create a nodal set
//! details:  create a calculix-style node section
//! -----------------------------------------------
void writeSolverFileClass::createNodalSet(const IndexedMapOfMeshDataSources &anIndexedMapOfFaceMeshDS,
                                          std::vector<int> &theNodeIDs)
{
    for(QMap<int,opencascade::handle<MeshVS_DataSource>>::const_iterator it = anIndexedMapOfFaceMeshDS.cbegin(); it!= anIndexedMapOfFaceMeshDS.cend(); ++it)
    {
        int bodyIndex = it.key();
        int offset = 0;
        for(int k=1; k<bodyIndex; k++)
        {
            if(myDB->MapOfIsActive.value(k)==true)
                offset = offset+myDB->ArrayOfMeshDS.value(k)->GetAllNodes().Extent();
        }

        //! current mesh datasource
        occHandle(MeshVS_DataSource) theCurMeshDS = it.value();
        for(TColStd_MapIteratorOfPackedMapOfInteger anIter(theCurMeshDS->GetAllNodes()); anIter.More(); anIter.Next())
        {
            int nodeID = anIter.Key()+offset;
            theNodeIDs.push_back(nodeID);
        }
    }
}

//! -------------------------------
//! function: createElementSurface
//! details:
//! -------------------------------
void writeSolverFileClass::createElementSurface(std::vector<int> &theElementIDs,
                                                std::vector<int> &theFaceNumbers,
                                                const IndexedMapOfMeshDataSources &anIndexedMapOfFaceMeshDS)
{
    for(QMap<int,opencascade::handle<MeshVS_DataSource>>::const_iterator it = anIndexedMapOfFaceMeshDS.cbegin(); it!= anIndexedMapOfFaceMeshDS.cend(); ++it)
    {
        int bodyIndex = it.key();
        int offset = 0;
        int K = 0;
        for(int m=1; m<bodyIndex; m++)
        {
            if (myDB->ArrayOfMeshDS.value(m).IsNull()) K=0;
            else K=myDB->ArrayOfMeshDS.value(m)->GetAllElements().Extent();
            offset = offset+K;
        }
        std::map<meshElement2D,std::vector<std::pair<int,int>>> facesToElements = bigMap.at(bodyIndex);

        opencascade::handle<MeshVS_DataSource> faceMeshDS = anIndexedMapOfFaceMeshDS.value(bodyIndex);
        for(TColStd_MapIteratorOfPackedMapOfInteger it(faceMeshDS->GetAllElements()); it.More(); it.Next())
        {
            int surfaceElementID = it.Key();

            meshElement2D aMeshElement2D;
            int NbNodes, nbuf[8];
            TColStd_Array1OfInteger nodeIDs(*nbuf,1,8);
            faceMeshDS->GetNodesByElement(surfaceElementID,nodeIDs,NbNodes);
            for(int i=1; i<=NbNodes; i++)
            {
                aMeshElement2D.nodeIDs<<nodeIDs(i);
            }
            switch(NbNodes)
            {
            case 3: aMeshElement2D.type = TRIG; break;
            case 4: aMeshElement2D.type = QUAD; break;
            }
            ////std::map<meshElement2D,std::vector<std::pair<int,int>>> facesToElements = bigMap.at(bodyIndex);
            std::map<meshElement2D,std::vector<std::pair<int,int>>>::iterator itt = facesToElements.find(aMeshElement2D);

            if(itt!=facesToElements.end())
            {
                std::pair<meshElement2D,std::vector<std::pair<int,int>>> apair = *itt;
                std::vector<std::pair<int,int>> vecFaceElements = apair.second;
                std::pair<int,int> p = vecFaceElements.at(0);

                int globalElementID = p.first;
                theElementIDs.push_back(globalElementID+offset);

                //! ------------------------------------------
                //! workaround for converting "S<n>" into "n"
                //! ------------------------------------------
                int CCXnumber = p.second;
                theFaceNumbers.push_back(CCXnumber);
            }
        }
     }
}

//! --------------------------------------------
//! function: getTreeItem
//! details:  already used in SimulationManager
//! --------------------------------------------
QExtendedStandardItem* writeSolverFileClass::getTreeItem(SimulationNodeClass::nodeType theNodeType)
{
    //! retrieve the root nodes
    QStandardItemModel *myModel = myDB->getModel();
    for(int k=0; k<myModel->invisibleRootItem()->child(0,0)->rowCount(); k++)
    {
        QExtendedStandardItem *item = static_cast<QExtendedStandardItem*>(myModel->invisibleRootItem()->child(0,0)->child(k,0));
        if(item->data(Qt::UserRole).value<SimulationNodeClass*>()->getType()==theNodeType)
        {
            return item;
        }
    }
    return NULL;
}

//! ----------------------------------------------------------------
//! function: ItemFromScope
//! details:  for a given shape in a geometry item, return the item
//! ----------------------------------------------------------------
QExtendedStandardItem* writeSolverFileClass::ItemFromScope(const TopoDS_Shape &aShape)
{
    QStandardItem *theGeometryRoot=this->getTreeItem(SimulationNodeClass::nodeType_geometry);
    int N = theGeometryRoot->rowCount();
    for(int k=0; k<N;k++)
    {
        QStandardItem *aGeometryItem = theGeometryRoot->child(k,0);
        SimulationNodeClass *aNode = aGeometryItem->data(Qt::UserRole).value<SimulationNodeClass*>();
        if(aNode->getType()==SimulationNodeClass::nodeType_pointMass) continue;
        int mapIndex = aNode->getPropertyValue<int>("Map index");
        TopoDS_Shape theShape = myDB->bodyMap.value(mapIndex);
        if(theShape==aShape) return static_cast<QExtendedStandardItem*>(aGeometryItem);
    }
    return Q_NULLPTR;
}

//! ------------------------------------------
//! function: writeDload
//! details:  suitable for applyng a pressure
//! ------------------------------------------
void writeSolverFileClass::writeDload(double aLoad, QString aName)
{
    ofstream myDload;
    myDload.setf(ios::scientific);
    myDload.precision(EXPFORMAT_PRECISION);

    ifstream mySet;

    QString extension=".dlo";
    QString extension1=".surf";
    //! this is the absolute path on the disk
    QString name = aName+extension;
    QString setName = aName;
    setName.chop(2);
    setName=setName+extension1;

    QString absFileName = myFileName.split("/").last();
    QString dirName = myFileName;
    dirName.chop(absFileName.size());
    name.prepend(dirName);
    setName.prepend(dirName);

    mySet.open(setName.toStdString());

    myDload.open(name.toStdString());

    //! write the header for the DLOAD
    myDload<<"*DLOAD"<<endl;
    //! assign LoadValue to each surface element
    //! In case of pressure a negative value means traction
    int theElementIDs;
    int theFaceNumbers;
    std::string s;
    getline(mySet,s);
    while(!mySet.eof())
    {
        getline(mySet,s);
        if(s.size()==0) break;
        char c;
        sscanf(s.c_str(),"%d, %c%d",&theElementIDs,&c,&theFaceNumbers);
        myDload<<theElementIDs<<",P"<<theFaceNumbers<<", "<<aLoad<<endl;
    }
    myDload.close();
    myInputFile<<"*INCLUDE, INPUT= "<<name.split("/").last().toStdString()<<endl;
}

//! ----------------------------------------------
//! function: writeDflux
//! details:  suitable for applyng a thermal flux
//! ----------------------------------------------
void writeSolverFileClass::writeDflux(double aLoad, QString aName)
{
    ofstream myDload;
    myDload.setf(ios::scientific);
    myDload.precision(EXPFORMAT_PRECISION);

    ifstream mySet;

    QString extension=".dlo";
    QString extension1=".surf";

    //! this is the absolute path on the disk
    QString name = aName+extension;
    QString setName = aName+extension1;

    QString absFileName = myFileName.split("/").last();
    QString dirName = myFileName;
    dirName.chop(absFileName.size());
    name.prepend(dirName);
    setName.prepend(dirName);

    mySet.open(setName.toStdString());
    myDload.open(name.toStdString());

    //! write the header for the DLOAD
    myDload<<"*DFLUX"<<endl;
    //! assign LoadValue to each surface element
    //! In case of pressure a negative value means traction
    int theElementIDs;
    int theFaceNumbers;
    std::string s;
    getline(mySet,s);
    while(!mySet.eof())
    {
        getline(mySet,s);
        if(s.size()==0) break;
        char c;
        sscanf(s.c_str(),"%d, %c%d",&theElementIDs,&c,&theFaceNumbers);
        myDload<<theElementIDs<<",S"<<theFaceNumbers<<", "<<aLoad<<endl;
    }
    myDload.close();
    myInputFile<<"*INCLUDE, INPUT= "<<name.split("/").last().toStdString()<<endl;
}

//! --------------------------------------------------
//! function: writeFilm
//! details:  suitable for applyng a film coefficient
//! --------------------------------------------------
void writeSolverFileClass::writeFilm(double aLoad, QString aName,double refTemperature)
{
    ofstream myDload;
    myDload.setf(ios::scientific);
    myDload.precision(EXPFORMAT_PRECISION);

    ifstream mySet;

    QString extension=".dlo";
    QString extension1=".surf";

    //! this is the absolute path on the disk
    QString name = aName+extension;
    QString setName = aName+extension1;

    QString absFileName = myFileName.split("/").last();
    QString dirName = myFileName;
    dirName.chop(absFileName.size());
    name.prepend(dirName);
    setName.prepend(dirName);

    mySet.open(setName.toStdString());

    myDload.open(name.toStdString());

    //! write the header for the DLOAD
    myDload<<"*FILM"<<endl;
    //! assign LoadValue to each surface element
    //! In case of pressure a negative value means traction
    int theElementIDs;
    int theFaceNumbers;
    std::string s;
    getline(mySet,s);
    while(!mySet.eof())
    {
        getline(mySet,s);
        if(s.size()==0) break;
        char c;
        sscanf(s.c_str(),"%d, %c%d",&theElementIDs,&c,&theFaceNumbers);
        myDload<<theElementIDs<<",F"<<theFaceNumbers<<", "<<refTemperature<<","<<aLoad<<endl;
    }
    myDload.close();
    myInputFile<<"*INCLUDE, INPUT= "<<name.split("/").last().toStdString()<<endl;
}

//! ----------------------------------
//! function: writeTemperatureHistory
//! details:
//! ----------------------------------
void writeSolverFileClass::writeTemperatureHistory(postObject pObject, QString tName)
{
    ofstream myTemperature;
    myTemperature.setf(ios::scientific);
    myTemperature.precision(EXPFORMAT_PRECISION);
    myTemperature.open(tName.toStdString());

    std::map<GeometryTag,std::vector<std::map<int,double>>> Tdata = pObject.getData();
    for(std::map<GeometryTag,std::vector<std::map<int,double>>>::iterator mapIt = Tdata.begin(); mapIt!=Tdata.end(); ++mapIt)
    {
        GeometryTag aloc = mapIt->first;

        std::vector<std::map<int,double>> lres = mapIt->second;     //extract the list of results
        std::map<int,double> res = lres.at(0);                      //extract the result from the list
        std::map<int,int> trans = OCCMeshToCCXmesh::perform(aloc,myDB);

        int l=1;
        for(std::map<int,double>::const_iterator it =res.cbegin(); it!=res.cend(); it++)
        {
            int keyOfl = trans.find(l)->first;
            myTemperature<<keyOfl<<", "<<it->second<<endl;
            l++;
        }
    }
    myInputFile<<"*INCLUDE, INPUT="<<tName.split("/").last().toStdString()<<endl;
}

//! --------------------------
//! function: writeGapElement
//! details:
//! --------------------------
void writeSolverFileClass::writeGapElement(const IndexedMapOfMeshDataSources &anIndexedMapOfFaceMeshDS,
                                           QString setName, double K, double F)
{
    QString gapuniName = setName + ".gap";
    //! this is the absolute path on the disk
    QString absFileName = myFileName.split("/").last();
    QString dirName = myFileName;
    dirName.chop(absFileName.size());
    gapuniName.prepend(dirName);

    ofstream myGapuni;
    myGapuni.setf(ios::scientific);
    myGapuni.precision(EXPFORMAT_PRECISION);
    myGapuni.open(gapuniName.toStdString());

    //! ----------------------------
    //! write the copied face mesh
    //! ----------------------------
    myGapuni<<"*NODE, NSET= N"<<setName.toStdString()<<endl;

    //QMap<std::pair<int,int>,QList<double>> gapInfo;
    std::map<std::pair<int,int>,QList<double>> gapInfo;

    int offset=0;
    for(QMap<int,occHandle(MeshVS_DataSource)>::const_iterator it = anIndexedMapOfFaceMeshDS.cbegin(); it!= anIndexedMapOfFaceMeshDS.cend(); ++it)
    {
        int bodyIndex = it.key();
        for(int k=1; k<bodyIndex; k++)
        {
            if(myDB->MapOfIsActive.value(k)==true)
                offset = offset+myDB->ArrayOfMeshDS.value(k)->GetAllNodes().Extent();
        }

        //! current mesh datasource
        occHandle(MeshVS_DataSource) theCurMeshDS = it.value();

        occHandle(Ng_MeshVS_DataSourceFace) &faceMeshDS = occHandle(Ng_MeshVS_DataSourceFace)::DownCast(theCurMeshDS);
        if(faceMeshDS->myNodeNormals.isEmpty()) faceMeshDS->computeNormalAtNodes();
        const QMap<int,QList<double>> &nodeNormals = faceMeshDS->myNodeNormals;

        for(TColStd_MapIteratorOfPackedMapOfInteger anIter(theCurMeshDS->GetAllNodes()); anIter.More(); anIter.Next())
        {
            totalNumberOfNodes++;
            double aCoordsBuf[3];
            TColStd_Array1OfReal aCoords(*aCoordsBuf,1,3);
            int nbNodes;
            MeshVS_EntityType aType;

            std::pair<int,int> aNodePair;
            QList<double> curNormal = nodeNormals.value(anIter.Key());
            aNodePair.first = anIter.Key()+offset;
            aNodePair.second=totalNumberOfNodes;

            theCurMeshDS->GetGeom(anIter.Key(),false,aCoords,nbNodes,aType);
            myGapuni<<aNodePair.second<<","<<aCoords(1)<<","<<aCoords(2)<<","<<aCoords(3)<<endl;

            //gapInfo.insert(aNodePair,curNormal);
            gapInfo.insert(std::make_pair(aNodePair,curNormal));
        }
    }

    //! ----------------------------
    //! lock the copied face's nodes
    //! ----------------------------
    myGapuni<<"*BOUNDARY"<<endl;
    myGapuni<<"N"<<setName.toStdString()<<","<<"1,"<<"3"<<endl;

    for(std::map<std::pair<int,int>,QList<double>>::iterator it=gapInfo.begin();it!=gapInfo.end();it++)
    //for(QMap<std::pair<int,int>,QList<double>>::iterator it=gapInfo.begin();it!=gapInfo.end();it++)
    {
        totalNumberOfElements++;
        //std::pair<int,int> aPair = it.key();
        //QList<double> nodeNormal = it.value();
        std::pair<int,int> aPair = it->first;
        QList<double> nodeNormal = it->second;

        //! write the gap element (connection btw faceMeshDS and the copied one)
        myGapuni<<"*ELEMENT,TYPE = GAPUNI, ELSET= G"<<totalNumberOfElements<<endl;
        myGapuni<<totalNumberOfElements<<", "<<aPair.first<<", "<<aPair.second<<endl;

        if(K==0.0) K=10e12;
        if(F==0.0) F=10e-5;

        //! write the 0 gap condition
        myGapuni<<"*GAP, ELSET= G"<<totalNumberOfElements<<endl;
        myGapuni<<"0.0, "<<nodeNormal.at(0)<<", "<<nodeNormal.at(1)<<", "<<nodeNormal.at(2)<<",,"<<K<<","<<F<<endl;
    }
}
