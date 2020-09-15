//! ----------------
//! custom includes
//! ----------------
#include "simulationdatabase.h"
#include "property.h"
#include "load.h"
#include "nodefactory.h"
#include "postobject.h"
#include "nodefactory.h"
#include <customtablemodel.h>
#include "ccxsolvermessage.h"
#include "maintreetools.h"
#include "solutioninfo.h"

//! ----
//! OCC
//! ----
#include <MeshVS_DataSource.hxx>
#include <TColStd_IndexedMapOfInteger.hxx>
#include <TColStd_MapIteratorOfPackedMapOfInteger.hxx>

//! ----
//! C++
//! ----
#include <iostream>

/*
//!---------------------------------
//! function: constructor - default
//! details:
//! --------------------------------
simulationDataBase::simulationDataBase(QObject *parent): meshDataBase(parent)
{
    cerr<<"simulationDataBase::simulationDataBase()->____default constructor called____"<<endl;

    //! create a connections root node
    this->createConnectionsRootNode();

    //! create a named selection node
    this->createNamedSelectionRootNode();

    //! create a analysis root node
    //! this->createThermalAnalysisRootNode();

    //! create a structural analysis root node
    this->createStructuralAnalysisRootNode();
    cerr<<"simulationDataBase::simulationDataBase()->____default constructor exiting____"<<endl;
}
*/

//!--------------------------------
//! function: constructor III
//! details:  from a list of nodes
//! -------------------------------
simulationDataBase::simulationDataBase(const QList<SimulationNodeClass*> listOfNodes,
                                       const QString &archiveFileName,
                                       QObject *parent):
    meshDataBase(listOfNodes, archiveFileName, parent)
{
    cout<<"simulationDataBase::simulationDataBase()->____CONSTRUCTOR III (FROM FILE) CALLED____"<<endl;

    QList<SimulationNodeClass*> listOfNodes_ = listOfNodes;
    for(int n=myRootItem->rowCount()-1; n>=0; n--)
    {
        QStandardItem *item = myRootItem->child(n,0);
        SimulationNodeClass *node = item->data(Qt::UserRole).value<SimulationNodeClass*>();
        if(node->isAnalysisRoot()) myRootItem->removeRow(n);
    }

    //! -------------
    //! generic data
    //! -------------
    QVariant data;

    //! ----------------
    //! Connection root
    //! ----------------
    for(QList<SimulationNodeClass*>::iterator it=listOfNodes_.begin(); it!=listOfNodes_.end();it++)
    {
        SimulationNodeClass *curNode = *it;
        SimulationNodeClass::nodeType curNodeType = curNode->getType();
        if(curNodeType != SimulationNodeClass::nodeType_connection) continue;

        //! ---------------
        //! child of model
        //! ---------------
        cout<<"____loading \"Connection\" root____"<<endl;
        QExtendedStandardItem *item = new QExtendedStandardItem();
        data.setValue(curNode);
        item->setData(data, Qt::UserRole);
        item->setData(curNode->getName(),Qt::DisplayRole);
        myRootItem->appendRow(item);
        ConnectionsRootItem = item;

        listOfNodes_.erase(it);
        break;
    }

    //! ---------------------------------
    //! dummy contact "Select from list"
    //! ---------------------------------
    for(QList<SimulationNodeClass*>::iterator it=listOfNodes_.begin(); it!=listOfNodes_.end();it++)
    {
        SimulationNodeClass *curNode = *it;
        SimulationNodeClass::nodeType curNodeType = curNode->getType();
        if(curNodeType != SimulationNodeClass::nodeType_connectionPair) continue;
        if(curNode->getName()!="Select from list") continue;

        cout<<"____loading dummy contact \"Select from list\"____"<<endl;
        QExtendedStandardItem *item = new QExtendedStandardItem();
        data.setValue(curNode);
        item->setData(data, Qt::UserRole);
        item->setData(curNode->getName(),Qt::DisplayRole);
        ConnectionsRootItem->appendRow(item);

        listOfNodes_.erase(it);
        //cout<<"____residual items: "<<listOfNodes_.length()<<"____"<<endl;
        break;
    }

    //! ------------------
    //! Connection groups
    //! ------------------
    for(QList<SimulationNodeClass*>::iterator it=listOfNodes_.begin(); it!=listOfNodes_.end();)
    {
        SimulationNodeClass *curNode = *it;
        SimulationNodeClass::nodeType curNodeType = curNode->getType();
        if(curNodeType != SimulationNodeClass::nodeType_connectionGroup)
        {
            it++;
            continue;
        }

        cout<<"____loading \"Connection group\"____"<<endl;
        QExtendedStandardItem *item = new QExtendedStandardItem();
        data.setValue(curNode);
        item->setData(data, Qt::UserRole);
        item->setData(curNode->getName(),Qt::DisplayRole);
        ConnectionsRootItem->appendRow(item);

        it = listOfNodes_.erase(it);
        //cout<<"____residual items: "<<listOfNodes_.length()<<"____"<<endl;
    }

    //! -----------------
    //! Connection pairs
    //! -----------------
    for(QList<SimulationNodeClass*>::iterator it=listOfNodes_.begin(); it!=listOfNodes_.end();)
    {
        SimulationNodeClass *curNode = *it;
        SimulationNodeClass::nodeType curNodeType = curNode->getType();
        if(curNodeType != SimulationNodeClass::nodeType_connectionPair)
        {
            it++;
            continue;
        }

        //! -----------------------------------------------------
        //! jump over the dummy contact pair "Select from list"
        //! this should never be reached since the node has been
        //! already erased from the list before
        //! -----------------------------------------------------
        if(curNode->getName() == "Select from list")
        {
            cout<<"____loading \"Connection pair\": \"Select from list found\". Jumping over____"<<endl;
            it++;
            continue;
        }

        //cout<<"____"<<curNode->getName().toStdString()<<"____"<<endl;
        //cout<<"____"<<curNode->getPropertyValue<QString>("Time tag").toStdString()<<"____"<<endl;

        if(curNode->getPropertyItem("Parent time tag")==Q_NULLPTR) exit(3);
        QString parentTimeTag = curNode->getPropertyValue<QString>("Parent time tag");

        //! ---------------------------
        //! scan the Connection groups
        //! ---------------------------
        for(int n=0; n<ConnectionsRootItem->rowCount(); n++)
        {
            QStandardItem *itemConnectionGroup = ConnectionsRootItem->child(n,0);
            SimulationNodeClass *curConnectionGroup = itemConnectionGroup->data(Qt::UserRole).value<SimulationNodeClass*>();

            if(curConnectionGroup->getType()==SimulationNodeClass::nodeType_connectionPair) continue;

            QString timeTag = curConnectionGroup->getPropertyValue<QString>("Time tag");
            if(timeTag == parentTimeTag)
            {
                cout<<"____loading \"Connection pair\": \""<<curNode->getName().toStdString()<<"\"____"<<endl;
                QExtendedStandardItem *item = new QExtendedStandardItem();
                data.setValue(curNode);
                item->setData(data, Qt::UserRole);
                item->setData(curNode->getName(),Qt::DisplayRole);
                itemConnectionGroup->appendRow(item);

                it = listOfNodes_.erase(it);
            }
        }
    }

    //! ---------------------
    //! Named selection root
    //! ---------------------
    for(QList<SimulationNodeClass*>::iterator it=listOfNodes_.begin(); it!=listOfNodes_.end();it++)
    {
        SimulationNodeClass *curNode = *it;
        if(curNode->getType()!=SimulationNodeClass::nodeType_namedSelection) continue;

        cout<<"____loading \"Named selection\" root____"<<endl;
        QExtendedStandardItem *item = new QExtendedStandardItem();
        data.setValue(curNode);
        item->setData(data, Qt::UserRole);
        item->setData(curNode->getName(),Qt::DisplayRole);
        myRootItem->appendRow(item);
        NamedSelectionRootItem = item;

        listOfNodes_.erase(it);
        break;
    }

    //cout<<"____residual items: "<<listOfNodes_.length()<<"____"<<endl;

    //! ------------------------
    //! dummy "Named selection"
    //! ------------------------
    for(QList<SimulationNodeClass*>::iterator it=listOfNodes_.begin(); it!=listOfNodes_.end();it++)
    {
        SimulationNodeClass *curNode = *it;
        SimulationNodeClass::nodeType curNodeType = curNode->getType();
        QString name = curNode->getName();
        if(curNodeType!=SimulationNodeClass::nodeType_namedSelectionGeometry) continue;
        if(name!="Select from list") continue;

        QExtendedStandardItem *item = new QExtendedStandardItem();
        data.setValue(curNode);
        item->setData(data, Qt::UserRole);
        item->setData(curNode->getName(),Qt::DisplayRole);
        cout<<"____loading dummy named selection \"Select from list\" item____"<<endl;
        NamedSelectionRootItem->appendRow(item);

        listOfNodes_.erase(it);
        break;
    }

    cout<<"____dummy \"Named selection\" added____"<<endl;

    //! -----------------------
    //! Named selections items
    //! -----------------------
    for(QList<SimulationNodeClass*>::iterator it=listOfNodes_.begin(); it!=listOfNodes_.end();)
    {
        SimulationNodeClass *curNode = *it;
        SimulationNodeClass::nodeType curNodeType = curNode->getType();
        if(curNodeType!=SimulationNodeClass::nodeType_namedSelectionGeometry)
        {
            it++;
            continue;
        }

        QExtendedStandardItem *item = new QExtendedStandardItem();
        data.setValue(curNode);
        item->setData(data, Qt::UserRole);
        item->setData(curNode->getName(),Qt::DisplayRole);
        cout<<"____loading named selection: \""<<curNode->getName().toStdString()<<"\"____"<<endl;
        NamedSelectionRootItem->appendRow(item);

        it = listOfNodes_.erase(it);
    }

    cout<<"____\"Named selection\" loaded____"<<endl;

    //! -----------------------------
    //! 1 - add the simulation roots
    //! -----------------------------
    for(QList<SimulationNodeClass*>::iterator it=listOfNodes_.begin(); it!=listOfNodes_.end();)
    {
        SimulationNodeClass *curNode = *it;
        if(curNode->isAnalysisRoot()==false)
        {
            it++;
            continue;
        }

        //! -------------------------
        //! create a simulation root
        //! -------------------------
        cout<<"____loading simulation root \""<<curNode->getName().toStdString()<<"\" item____"<<endl;
        QExtendedStandardItem *item = new QExtendedStandardItem();
        data.setValue(curNode);
        item->setData(data, Qt::UserRole);
        item->setData(curNode->getName(),Qt::DisplayRole);
        myRootItem->appendRow(item);

        it = listOfNodes_.erase(it);
    }

    cout<<"____Analysis root items____"<<endl;

    //! ------------------------------------
    //! 2 - add the analysis settings items
    //! ------------------------------------
    for(QList<SimulationNodeClass*>::iterator it=listOfNodes_.begin(); it!=listOfNodes_.end();)
    {
        SimulationNodeClass *curNode = *it;
        if(curNode->isAnalysisSettings()==false)
        {
            it++;
            continue;
        }
        QString curParentTimeTag = curNode->getPropertyValue<QString>("Parent time tag");

        for(int n=0; n<myRootItem->rowCount(); n++)
        {
            QStandardItem *curSimulationRootItem = myRootItem->child(n,0);
            SimulationNodeClass *curSimulationNodeRoot = curSimulationRootItem->data(Qt::UserRole).value<SimulationNodeClass*>();
            if(curSimulationNodeRoot->isAnalysisRoot()==false) continue;
            QString timeTag = curSimulationNodeRoot->getPropertyValue<QString>("Time tag");
            if(curParentTimeTag == timeTag)
            {
                //! ---------------------------------------------
                //! create and append the Analysis settings item
                //! ---------------------------------------------
                cout<<"____loading analysis settings: \""<<curNode->getName().toStdString()<<"\" item____"<<endl;
                QExtendedStandardItem *item = new QExtendedStandardItem();
                data.setValue(curNode);
                item->setData(data, Qt::UserRole);
                item->setData(curNode->getName(),Qt::DisplayRole);
                curSimulationRootItem->appendRow(item);

                it = listOfNodes_.erase(it);
            }
        }
    }

    cout<<"____\"Analysis settings\" items added____"<<endl;

    //! -----------------------------------
    //! 3 - add the simulation setup nodes
    //! -----------------------------------
    //! -----------------------------------------------------
    //! scan the main tree and find the analysis setup nodes
    //! -----------------------------------------------------
    for(int n=0; n<myRootItem->rowCount(); n++)
    {
        QStandardItem *curSimulationRootItem = myRootItem->child(n,0);
        SimulationNodeClass *curSimulationNodeRoot = curSimulationRootItem->data(Qt::UserRole).value<SimulationNodeClass*>();
        if(curSimulationNodeRoot->isAnalysisRoot()==false) continue;

        //! ------------------------------------------------
        //! an analysis root has been found within the tree
        //! ------------------------------------------------
        QString analysisRootTimeTag = curSimulationNodeRoot->getPropertyValue<QString>("Time tag");
        cout<<"____ANALYSIS ROOT: "<<analysisRootTimeTag.toStdString()<<"____"<<endl;

        //! ---------------------------------------------------
        //! retrieve the nodes belonging to the current branch
        //! ---------------------------------------------------
        std::map<QString,SimulationNodeClass*> timeTagToNodeMap;
        std::vector<unsigned long long int> vecKeys;
        for(QList<SimulationNodeClass*>::iterator it=listOfNodes_.begin(); it!=listOfNodes_.end();)
        {
            SimulationNodeClass *curNode = *it;
            bool isSetUpNode = curNode->isSimulationSetUpNode();
            bool isChildSetUpNode = curNode->isChildSimulationSetUpNode();

            QString curParentTimeTag = curNode->getPropertyValue<QString>("Parent time tag");
            QString curTimeTag = curNode->getPropertyValue<QString>("Time tag");

            if(isSetUpNode==false)
            {
                it++;
                continue;
            }

            if(curParentTimeTag!=analysisRootTimeTag)
            {
                it++;
                continue;
            }
            //cout<<"____NODE TO ATTACH FOUND____"<<endl;
            //cout<<"____PARENT TIME TAG: "<<curParentTimeTag.toStdString()<<"____"<<endl;
            //cout<<"____TIME TAG: "<<curTimeTag.toStdString()<<"____"<<endl;

            //! ------------------------------
            //! remove the node from the list
            //! ------------------------------
            it = listOfNodes_.erase(it);

            //! ------------------------------------
            //! store the current node into the map
            //! ------------------------------------
            std::pair<QString,SimulationNodeClass*> element;
            element.first = curTimeTag;
            element.second = curNode;

            timeTagToNodeMap.insert(element);
            vecKeys.push_back(curTimeTag.toULongLong());
            cout<<"____a check =>"<<curTimeTag.toULongLong()<<"____"<<endl;
        }
        //! ---------------------------
        //! sort the map using the key
        //! ---------------------------
        std::sort(vecKeys.begin(),vecKeys.end());


        //! ---------------------------------------------------------
        //! retrieve the child nodes belonging to the current branch
        //! ---------------------------------------------------------
        std::map<QString,std::pair<QString,SimulationNodeClass*>> timeTagToChildNodeMap;
        for(QList<SimulationNodeClass*>::iterator it=listOfNodes_.begin(); it!=listOfNodes_.end();)
        {
            SimulationNodeClass *curNode = *it;
            //bool isSetUpNode = curNode->isSimulationSetUpNode();
            bool isChildSetUpNode = curNode->isChildSimulationSetUpNode();

            if(isChildSetUpNode == false)
            {
                continue;
            }

            QString curChildParentTimeTag = curNode->getPropertyValue<QString>("Parent time tag");
            QString curChildTimeTag = curNode->getPropertyValue<QString>("Time tag");

            if(curChildParentTimeTag!=analysisRootTimeTag)
            {
                it++;
                continue;
            }
            //cout<<"____NODE TO ATTACH FOUND____"<<endl;
            //cout<<"____PARENT TIME TAG: "<<curParentTimeTag.toStdString()<<"____"<<endl;
            //cout<<"____TIME TAG: "<<curTimeTag.toStdString()<<"____"<<endl;

            for(int n=0; n<vecKeys.size(); n++)
            {
                //QString key = QString("%1").arg(vecKeys[n]);
                unsigned long long int parentKey = vecKeys[n];
                if(parentKey == curChildParentTimeTag.toLongLong())
                {
                    //! ------------------------------------
                    //! store the current node into the map
                    //! ------------------------------------
                    std::pair<QString,SimulationNodeClass*> element;
                    element.first = curChildTimeTag;
                    element.second = curNode;

                    timeTagToNodeMap.insert(curChildParentTimeTag,element);
                }
            }

            //! ------------------------------
            //! remove the node from the list
            //! ------------------------------
            it = listOfNodes_.erase(it);
            //cout<<"____a check =>"<<curTimeTag.toULongLong()<<"____"<<endl;
        }

        //! -----------------------------
        //! append in the original order
        //! -----------------------------
        for(int n=0; n<vecKeys.size(); n++)
        {
            QString key = QString("%1").arg(vecKeys[n]);

            //cout<<"____Attaching node with time tag: "<<key.toStdString()<<"____"<<endl;

            std::map<QString,SimulationNodeClass*>::iterator mapit = timeTagToNodeMap.find(key);
            if(mapit!=timeTagToNodeMap.end())
            {
                std::pair<QString,SimulationNodeClass*> element = *mapit;
                SimulationNodeClass* aNode = element.second;

                QExtendedStandardItem *item = new QExtendedStandardItem();
                data.setValue(aNode);
                item->setData(data, Qt::UserRole);
                item->setData(aNode->getName(),Qt::DisplayRole);
                curSimulationRootItem->appendRow(item);


            }
        }
    }
    cout<<"____simulation set up nodes added____"<<endl;

    //! ---------------------------
    //! 4 - add the solution items
    //! ---------------------------
    for(QList<SimulationNodeClass*>::iterator it=listOfNodes_.begin(); it!=listOfNodes_.end();)
    {
        SimulationNodeClass *curNode = *it;
        if(curNode->isSolution()==false)
        {
            it++;
            continue;
        }
        QString curParentTimeTag = curNode->getPropertyValue<QString>("Parent time tag");

        for(int n=0; n<myRootItem->rowCount(); n++)
        {
            QStandardItem *curRootItem = myRootItem->child(n,0);
            SimulationNodeClass *curNodeRoot = curRootItem->data(Qt::UserRole).value<SimulationNodeClass*>();
            if(!curNodeRoot->isAnalysisRoot()) continue;

            QString curRootTimeTag = curNodeRoot->getPropertyValue<QString>("Time tag");
            if(curParentTimeTag == curRootTimeTag)
            {
                //! ------------------------------------
                //! create and append the Solution item
                //! ------------------------------------
                cout<<"____loading Solution item____"<<endl;
                QExtendedStandardItem *item = new QExtendedStandardItem();
                data.setValue(curNode);
                item->setData(data, Qt::UserRole);
                item->setData(curNode->getName(),Qt::DisplayRole);
                curRootItem->appendRow(item);

                it = listOfNodes_.erase(it);
            }
        }
    }

    cout<<"____\"Solution\" items added____"<<endl;

    //! -----------------------------------------
    //! 5 - add the "Solution information" items
    //! -----------------------------------------
    for(QList<SimulationNodeClass*>::iterator it=listOfNodes_.begin(); it!=listOfNodes_.end();)
    {
        SimulationNodeClass *curNode = *it;
        cout<<"____"<<curNode->getName().toStdString()<<"____"<<endl;
        //if(curNode->isSolutionInformation()==true) exit(19);

        if(curNode->isSolutionInformation()==false)
        {
            it++;
            continue;
        }
        cout<<"____Solution information found____"<<endl;
        QString curParentTimeTag = curNode->getPropertyValue<QString>("Parent time tag");
        //exit(20);

        //! ---------------------------------
        //! search for an Analysis root item
        //! ---------------------------------
        for(int n=0; n<myRootItem->rowCount(); n++)
        {
            QStandardItem *curSimulationRootItem = myRootItem->child(n,0);
            SimulationNodeClass *curSimulationNodeRoot = curSimulationRootItem->data(Qt::UserRole).value<SimulationNodeClass*>();
            if(curSimulationNodeRoot->isAnalysisRoot()==false) continue;

            //! ----------------------------------------------------------------
            //! Retrieve the item "Solution" within the current analysis branch
            //! ----------------------------------------------------------------
            QStandardItem *curSolutionItem = curSimulationRootItem->child(curSimulationRootItem->rowCount()-1,0);
            SimulationNodeClass *curSolutionNode = curSolutionItem->data(Qt::UserRole).value<SimulationNodeClass*>();
            QString timeTag = curSolutionNode->getPropertyValue<QString>("Time tag");
            if(timeTag == curParentTimeTag)
            {
                //! ------------------------------------------------
                //! create and append the Solution information item
                //! ------------------------------------------------
                cout<<"____loading \"Solution information\"____"<<endl;

                //! --------------------------------
                //! compatibility with old versions
                //! --------------------------------
                if(curNode->getPropertyItem("Discrete time map")==Q_NULLPTR)
                {
                    QMap<double,QVector<int>> dtm;
                    QVariant data;
                    data.setValue(dtm);
                    curNode->addProperty(Property("Discrete time map",data,Property::PropertyGroup_Hidden),0);
                }
                if(curNode->getPropertyItem("Solver output")==Q_NULLPTR)
                {
                    CCXSolverMessage msg;
                    QVariant data;
                    data.setValue(msg);
                    curNode->addProperty(Property("Solver output",data,Property::PropertyGroup_Hidden),1);
                }

                QExtendedStandardItem *item = new QExtendedStandardItem();
                data.setValue(curNode);
                item->setData(data, Qt::UserRole);
                item->setData(curNode->getName(),Qt::DisplayRole);
                curSolutionItem->appendRow(item);

                it = listOfNodes_.erase(it);
                cout<<"____residual items: "<<listOfNodes_.length()<<"____"<<endl;
            }
        }
    }

    cout<<"____\"Solution information\" items added____"<<endl;

    //! ----------------------------------
    //! 6 - add the post processing items
    //! ----------------------------------
    for(QList<SimulationNodeClass*>::iterator it=listOfNodes_.begin(); it!=listOfNodes_.end();)
    {
        SimulationNodeClass *curNode = *it;
        if(curNode->isAnalysisResult()==false)
        {
            it++;
            continue;
        }
        cout<<"____analysis result: \""<<curNode->getName().toStdString()<<"\" found____"<<endl;
        QString curParentTimeTag = curNode->getPropertyValue<QString>("Parent time tag");

        //! ---------------------------------
        //! search for an Analysis root item
        //! ---------------------------------
        for(int n=0; n<myRootItem->rowCount(); n++)
        {
            QStandardItem *curSimulationRootItem = myRootItem->child(n,0);
            SimulationNodeClass *curSimulationNodeRoot = curSimulationRootItem->data(Qt::UserRole).value<SimulationNodeClass*>();
            if(curSimulationNodeRoot->isAnalysisRoot()==false) continue;



            //! -----------------------------
            //! retrieve the "Solution" item
            //! -----------------------------
            QStandardItem *curSolutionItem = curSimulationRootItem->child(curSimulationRootItem->rowCount()-1,0);
            SimulationNodeClass *curSolutionNode = curSolutionItem->data(Qt::UserRole).value<SimulationNodeClass*>();
            QString timeTag = curSolutionNode->getPropertyValue<QString>("Time tag");

            if(timeTag == curParentTimeTag)
            {
                QExtendedStandardItem *item = new QExtendedStandardItem();
                data.setValue(curNode);
                item->setData(data, Qt::UserRole);
                item->setData(curNode->getName(),Qt::DisplayRole);
                curSolutionItem->appendRow(item);

                it = listOfNodes_.erase(it);
                cout<<"____residual items: "<<listOfNodes_.length()<<"____"<<endl;
            }
        }
    }
}

//!--------------------------
//! function: constructor II
//! details:
//! -------------------------
simulationDataBase::simulationDataBase(const TopoDS_Shape &shape, const QString& theFilePath, QObject *parent):
    meshDataBase(shape, theFilePath, parent)
{
    cout<<"simulationDataBase::simulationDataBase()->____CONSTRUCTOR II CALLED____"<<endl;

    //! create a connections root node
    this->createConnectionsRootNode();

    //! create a named selection node
    this->createNamedSelectionRootNode();

    //! create a structural analysis root node
    //this->createStructuralAnalysisRootNode();
}

//!----------------------
//! function: destructor
//! details:
//! ---------------------
simulationDataBase::~simulationDataBase()
{
    cout<<"simulationDataBase::~simulationDataBase()->____DESTRUCTOR CALLED____"<<endl;
}

//! ------------------------------
//! function: createStandardModel
//! details:
//! ------------------------------
void simulationDataBase::createStandardModel()
{
    cout<<"simulationDataBase::createStandardModel()->____function called____"<<endl;
    meshDataBase::createStandardModel();
    this->createStructuralAnalysisRootNode();
}

//! ------------------------------------
//! function: createConnectionsRootNode
//! details:
//! ------------------------------------
void simulationDataBase::createConnectionsRootNode()
{
    cout<<"simulationDataBase::createConnectionsRootNode()->____function called____"<<endl;

    QVector<Property> props;

    Property property_autoDetection("Auto detection",false,Property::PropertyGroup_AutoDetection);
    Property property_transparency("Transparency",0.7,Property::PropertyGroup_Transparency);
    props.push_back(property_autoDetection);
    props.push_back(property_transparency);

    SimulationNodeClass *nodeConnectionsRoot = new SimulationNodeClass("Connections",SimulationNodeClass::nodeType_connection,props,this);

    //! --------------------------------
    //! create the connection root item
    //! --------------------------------
    ConnectionsRootItem = new QExtendedStandardItem();
    ConnectionsRootItem->setData("Connections",Qt::DisplayRole);
    QVariant data;
    data.setValue(nodeConnectionsRoot);
    ConnectionsRootItem->setData(data, Qt::UserRole);

    //! -----------------------------
    //! time tag and parent time tag
    //! -----------------------------
    nodeConnectionsRoot->addTimeTag();
    QString timeTag = myRootItem->data(Qt::UserRole).value<SimulationNodeClass*>()->getPropertyValue<QString>("Time tag");
    data.setValue(timeTag);
    nodeConnectionsRoot->addProperty(Property("Parent time tag",data,Property::PropertyGroup_Identifier));

    //! ----------------------------------------------
    //! append the connections item to the model root
    //! ----------------------------------------------
    myRootItem->appendRow(ConnectionsRootItem);

    //! ---------------------------------------
    //! create a dummy item "Select from list"
    //! ---------------------------------------
    SimulationNodeClass *nodeEmptyContact = nodeFactory::nodeFromScratch(SimulationNodeClass::nodeType_connectionPair);
    nodeEmptyContact->setName(QString("Select from list"));

    //! -------------------------------------------------------------------------
    //! the time tag property is created by the node factory: do not add it here
    //! the Parent time tag should be replaced
    //! -------------------------------------------------------------------------
    //nodeEmptyContact->addTimeTag();
    timeTag = nodeConnectionsRoot->getPropertyValue<QString>("Time tag");
    data.setValue(timeTag);
    nodeEmptyContact->replaceProperty("Parent time tag",Property("Parent time tag",data,Property::PropertyGroup_Identifier));

    data.setValue(nodeEmptyContact);
    QExtendedStandardItem *itemnodeEmptyContact = new QExtendedStandardItem();
    itemnodeEmptyContact->setData(data,Qt::UserRole);
    data.setValue(QString("Select from list"));
    itemnodeEmptyContact->setData(data,Qt::DisplayRole);
    ConnectionsRootItem->appendRow(itemnodeEmptyContact);

    cout<<"simulationDataBase::createConnectionsRootNode()->____function exiting____"<<endl;
}

//! ---------------------------------------
//! function: createNamedSelectionRootNode
//! details:
//! ---------------------------------------
void simulationDataBase::createNamedSelectionRootNode()
{
    cout<<"simulationDataBase::createNamedSelectionRootNode()->____function called____"<<endl;
    QVariant data;
    QVector<Property> props;
    Property property_showAnnotation("Show annotations",false,Property::PropertyGroup_Display);
    props.push_back(property_showAnnotation);

    SimulationNodeClass *nodeNamedSelectionRoot = new SimulationNodeClass("Named selections",SimulationNodeClass::nodeType_namedSelection,props,this);

    //! -----------------------------
    //! time tag and parent time tag
    //! -----------------------------
    nodeNamedSelectionRoot->addTimeTag();
    QString timeTag = myRootItem->data(Qt::UserRole).value<SimulationNodeClass*>()->getPropertyValue<QString>("Time tag");
    data.setValue(timeTag);
    nodeNamedSelectionRoot->addProperty(Property("Parent time tag",data,Property::PropertyGroup_Identifier));
    NamedSelectionRootItem = new QExtendedStandardItem();

    NamedSelectionRootItem->setData("Named selections",Qt::DisplayRole);
    data.setValue(nodeNamedSelectionRoot);
    NamedSelectionRootItem->setData(data, Qt::UserRole);

    //! --------------------------------------------------
    //! an empty named selection item: this will be the
    //! first combo box item when activating the editor
    //! For this item "Geometry" is empty and also "Tags"
    //! --------------------------------------------------
    QVector<Property> propsEmptyNS;
    std::vector<GeometryTag> vecLoc;
    data.setValue(vecLoc);
    Property property_scope("Geometry",data,Property::PropertyGroup_Scope);
    Property property_tags("Tags",data,Property::PropertyGroup_Scope);
    propsEmptyNS.push_back(property_scope);
    propsEmptyNS.push_back(property_tags);

    SimulationNodeClass *nodeEmptyNS = new SimulationNodeClass("Select from list",SimulationNodeClass::nodeType_namedSelectionGeometry,propsEmptyNS,this);

    //! -----------------------------
    //! time tag and parent time tag
    //! -----------------------------
    nodeEmptyNS->addTimeTag();
    timeTag = myRootItem->data(Qt::UserRole).value<SimulationNodeClass*>()->getPropertyValue<QString>("Time tag");
    data.setValue(timeTag);
    Property("Parent time tag",data,Property::PropertyGroup_Identifier);

    QExtendedStandardItem *itemNS_empty = new QExtendedStandardItem();
    itemNS_empty->setData(nodeEmptyNS->getName(),Qt::DisplayRole);
    data.setValue(nodeEmptyNS);
    itemNS_empty->setData(data,Qt::UserRole);
    NamedSelectionRootItem->appendRow(itemNS_empty);

    //! ------------------------------
    //! append to the model root item
    //! ------------------------------
    myRootItem->appendRow(NamedSelectionRootItem);
    //cout<<"simulationDataBase::createNamedSelectionRootNode()->____function exiting____"<<endl;
}

//! -----------------------------------------
//! function: createCombinedAnalysisRootNode
//! details:
//! -----------------------------------------
void simulationDataBase::createCombinedAnalysisRootNode()
{
    cout<<"simulationDataBase::createStructuralAnalysisRootNode()->____creating a structural analysis branch____"<<endl;

    QVariant data;
    QString name = "Combined analysis";

    //! ---------------
    //! the properties
    //! ---------------
    QVector<Property> props;

    //! ------------------------------------------------
    //! a default value for the environment temperature
    //! ------------------------------------------------
    double Tenv = 295;
    Property prop_envTemp("Environment temperature",QVariant(Tenv),Property::PropertyGroup_EnvironmentTemperature);
    props.push_back(prop_envTemp);

    //! --------------
    //! node creation
    //! --------------
    SimulationNodeClass *nodeCombinedAnalysis = new SimulationNodeClass(name,SimulationNodeClass::nodeType_combinedAnalysis,props,this);

    //! -----------------------------
    //! time tag and parent time tag
    //! -----------------------------
    nodeCombinedAnalysis->addTimeTag();
    QString timeTag = myRootItem->data(Qt::UserRole).value<SimulationNodeClass*>()->getPropertyValue<QString>("Time tag");
    data.setValue(timeTag);
    nodeCombinedAnalysis->addProperty(Property("Parent time tag",data,Property::PropertyGroup_Identifier));

    //! --------------------------------------------------
    //! create the item and append to the model root item
    //! --------------------------------------------------
    QExtendedStandardItem *CombinedAnalysisItem = new QExtendedStandardItem();
    CombinedAnalysisItem->setData("Combined analysis",Qt::DisplayRole);
    data.setValue(nodeCombinedAnalysis);
    CombinedAnalysisItem->setData(data, Qt::UserRole);
    myRootItem->appendRow(CombinedAnalysisItem);

    //! ---------------------------------------------
    //! create the item combined "Analysis settings"
    //! ---------------------------------------------
    QExtendedStandardItem *CombinedAnalysisSettingsItem = new QExtendedStandardItem();

    //! ----------------------------
    //! the default Number of steps
    //! ----------------------------
    data.setValue(1);
    Property property_numberOfSteps("Number of steps",data,Property::PropertyGroup_StepControls);

    //! --------------------------
    //! the "Current step number"
    //! --------------------------
    data.setValue(1);
    Property property_currentStepNumber("Current step number",data,Property::PropertyGroup_StepControls);

    //! ------------------------------
    //! the default analysis end time
    //! ------------------------------
    double endTime =1.0;
    data.setValue(endTime);
    Property property_stepEndTime("Step end time",data,Property::PropertyGroup_StepControls);

    //! --------------
    //! analysis type
    //! --------------
    Property::analysisType analysisType = Property::analysisType_structural;
    data.setValue(analysisType);
    Property property_analysisType("Analysis type",data,Property::PropertyGroup_StepControls);

    //! -----------------------
    //! steady state/transient
    //! -----------------------
    Property::timeIntegration timeIntegration = Property::timeIntegration_steadyState;
    data.setValue(timeIntegration);
    Property property_timeIntegration("Static/Transient",data,Property::PropertyGroup_StepControls);

    //! --------------
    //! time stepping
    //! --------------
    Property::autoTimeStepping theTimeStepping = Property::autoTimeStepping_ProgramControlled;
    data.setValue(theTimeStepping);
    Property property_autoTimeStepping("Auto time stepping",data,Property::PropertyGroup_StepControls);

    //! ------------
    //! Solver type
    //! ------------
    Property::solverType theSolverType = Property::solverType_programControlled;
    data.setValue(theSolverType);
    Property property_solverType("Solver type",data,Property::PropertyGroup_SolverControls);

    //! -----------------------------
    //! non linear geometry
    //! 0 -> NL "Off" 1 -> NL "On"
    //! -----------------------------
    Property property_NLgeometry("Large deflection",int(0),Property::PropertyGroup_SolverControls);

    //! ------------------
    //! number of threads
    //! ------------------
    int numberOfThreads = 4;
    data.setValue(numberOfThreads);
    Property property_numberOfThreads("Number of threads",data,Property::PropertyGroup_SolverControls);

    props.push_back(property_numberOfSteps);
    props.push_back(property_currentStepNumber);
    props.push_back(property_stepEndTime);
    props.push_back(property_analysisType);
    props.push_back(property_timeIntegration);
    props.push_back(property_autoTimeStepping);
    props.push_back(property_solverType);
    props.push_back(property_numberOfThreads);
    props.push_back(property_NLgeometry);

    //! -------------------------------------
    //! Time incrementation
    //! 0 = > Custom 1 => program controlled
    //! -------------------------------------
    int timeIncrementation = 0;
    data.setValue(timeIncrementation);
    Property property_timeIncrementation("Time incrementation",data,Property::PropertyGroup_TimeIncrementation);
    props.push_back(property_timeIncrementation);

    data.setValue(timeIncrementation);
    QVector<QVariant> values4;
    values4.push_back(data);
    load loadTimeIncrementationType(values4,Property::loadType_timeIncrementation);

    int I_0 = 4;
    int I_R = 8;
    int I_P = 9;
    int I_C = 16;
    int I_L = 10;
    int I_G = 4;
    int I_S = -1;       //! "-1" means "unused"
    int I_A = 5;
    int I_J = -1;       //! "-1" means "unused"
    int I_T = -1;       //! "-1" means "unused"

    data.setValue(I_0);
    Property property_I_0("I_0",data,Property::PropertyGroup_TimeIncrementation);
    props.push_back(property_I_0);
    data.setValue(I_R);
    Property property_I_R("I_R",data,Property::PropertyGroup_TimeIncrementation);
    props.push_back(property_I_R);
    data.setValue(I_P);
    Property property_I_P("I_P",data,Property::PropertyGroup_TimeIncrementation);
    props.push_back(property_I_P);
    data.setValue(I_C);
    Property property_I_C("I_C",data,Property::PropertyGroup_TimeIncrementation);
    props.push_back(property_I_C);
    data.setValue(I_L);
    Property property_I_L("I_L",data,Property::PropertyGroup_TimeIncrementation);
    props.push_back(property_I_L);
    data.setValue(I_G);
    Property property_I_G("I_G",data,Property::PropertyGroup_TimeIncrementation);
    props.push_back(property_I_G);
    data.setValue(I_S);
    Property property_I_S("I_S",data,Property::PropertyGroup_TimeIncrementation);
    props.push_back(property_I_S);
    data.setValue(I_A);
    Property property_I_A("I_A",data,Property::PropertyGroup_TimeIncrementation);
    props.push_back(property_I_A);
    data.setValue(I_J);
    Property property_I_J("I_J",data,Property::PropertyGroup_TimeIncrementation);
    props.push_back(property_I_J);
    data.setValue(I_T);
    Property property_I_T("I_T",data,Property::PropertyGroup_TimeIncrementation);
    props.push_back(property_I_T);

    QVector<int> timeIncrementationParameters;
    timeIncrementationParameters.push_back(I_0);
    timeIncrementationParameters.push_back(I_R);
    timeIncrementationParameters.push_back(I_P);
    timeIncrementationParameters.push_back(I_C);
    timeIncrementationParameters.push_back(I_L);
    timeIncrementationParameters.push_back(I_G);
    timeIncrementationParameters.push_back(I_S);
    timeIncrementationParameters.push_back(I_A);
    timeIncrementationParameters.push_back(I_J);
    timeIncrementationParameters.push_back(I_T);

    data.setValue(timeIncrementationParameters);
    QVector<QVariant> values3;
    values3.append(data);
    load loadTimeIncrementation(values3,Property::loadType_timeIncrementationParameters);

    //! ----------------------------------------------------------
    //! cutback type: "0" => "Program controlled" "1" => "Custom"
    //! ----------------------------------------------------------
    int typeOfCutback = 0;
    data.setValue(typeOfCutback);
    Property property_cutBackType("Cutback factors",data,Property::PropertyGroup_CutBack);
    props.push_back(property_cutBackType);

    data.setValue(typeOfCutback);
    QVector<QVariant> values5;
    values5.push_back(data);
    load loadCutBackType(values5,Property::loadType_cutBack);

    //! -------------------------------------
    //! cutback factors: CCX solver defaluts
    //! -------------------------------------
    double D_f = 0.25;
    double D_C = 0.5;
    double D_B = 0.75;
    double D_A = 0.85;
    double D_S = -1.0;      //! "-1.0" means "currently unused"
    double D_H = -1.0;      //! "-1.0" means "currently unused"
    double D_D = 1.5;
    double W_G = -1.0;      //! "-1.0" means "currently unused"

    data.setValue(D_f);
    Property property_D_f("D_f",data,Property::PropertyGroup_CutBack);
    data.setValue(D_C);
    Property property_D_C("D_C",data,Property::PropertyGroup_CutBack);
    data.setValue(D_B);
    Property property_D_B("D_B",data,Property::PropertyGroup_CutBack);
    data.setValue(D_A);
    Property property_D_A("D_A",data,Property::PropertyGroup_CutBack);
    data.setValue(D_S);
    Property property_D_S("D_S",data,Property::PropertyGroup_CutBack);
    data.setValue(D_H);
    Property property_D_H("D_H",data,Property::PropertyGroup_CutBack);
    data.setValue(D_D);
    Property property_D_D("D_D",data,Property::PropertyGroup_CutBack);
    data.setValue(W_G);
    Property property_W_G("W_G",data,Property::PropertyGroup_CutBack);

    props.push_back(property_D_f);
    props.push_back(property_D_C);
    props.push_back(property_D_B);
    props.push_back(property_D_A);
    props.push_back(property_D_S);
    props.push_back(property_D_H);
    props.push_back(property_D_D);
    props.push_back(property_W_G);

    QVector<double> cutBackParameters;
    cutBackParameters.push_back(D_f);
    cutBackParameters.push_back(D_C);
    cutBackParameters.push_back(D_B);
    cutBackParameters.push_back(D_A);
    cutBackParameters.push_back(D_S);
    cutBackParameters.push_back(D_H);
    cutBackParameters.push_back(D_D);
    cutBackParameters.push_back(W_G);
    data.setValue(cutBackParameters);
    QVector<QVariant> values6;
    values6.push_back(data);
    load loadCutBackParameters(values6,Property::loadType_cutBackParameters);

    //! ---------------------------------------------------------
    //! Convergence criteria - "Flux convergence"
    //! 0 = > Deactivated 1 => On 2 => Program controlled
    //! The item is created with the option "2" using default
    //! values for the parameters
    //! ---------------------------------------------------------
    int fluxConvergenceType = 2;
    data.setValue(fluxConvergenceType);
    Property property_fluxConvergence("Flux convergence",data,Property::PropertyGroup_ConvergenceCriteria);
    props.append(property_fluxConvergence);

    //! q_alpha_u   => "0" is Program controlled
    double q_alpha_u = 0.0;
    data.setValue(q_alpha_u);
    Property property_value("--Value",data,Property::PropertyGroup_ConvergenceCriteria);
    props.append(property_value);

    //! q_alpha_0   => "0" is Program controlled
    double q_alpha_0 = 0.0;
    data.setValue(q_alpha_0);
    Property property_q_alpha_0("--q_alpha_0",data,Property::PropertyGroup_ConvergenceCriteria);
    props.append(property_q_alpha_0);

    //! First tolerance in case of non zero flux
    double R_alpha_n = 0.001;
    data.setValue(R_alpha_n);
    Property property_R_alpha_n("--R_alpha_n",data,Property::PropertyGroup_ConvergenceCriteria);
    props.append(property_R_alpha_n);

    //! Second tolerance in case of non zero flux
    double R_alpha_P = 0.02;
    data.setValue(R_alpha_P);
    Property property_R_alpha_P("--R_alpha_P",data,Property::PropertyGroup_ConvergenceCriteria);
    props.append(property_R_alpha_P);

    //! Third tolerance parameter: R_alpha_l
    double R_alpha_l = 1e-8;
    data.setValue(R_alpha_l);
    Property property_R_alpha_l("--R_alpha_l",data,Property::PropertyGroup_ConvergenceCriteria);
    props.append(property_R_alpha_l);

    //! zero flux parameter
    double epsilon_alpha = 1e-5;
    data.setValue(epsilon_alpha);
    Property property_epsilon_n("--epsilon_alpha",data,Property::PropertyGroup_ConvergenceCriteria);
    props.append(property_epsilon_n);

    //! --------------------------------------------------
    //! Convergence criteria: "Solution convergence"
    //! 0 = > Deactivated 1 => On 2 => Program controlled
    //! --------------------------------------------------
    int solutionConvergenceType = 2;
    data.setValue(solutionConvergenceType);
    Property property_solutionConvergence("Solution convergence",data,Property::PropertyGroup_ConvergenceCriteria);
    props.append(property_solutionConvergence);

    //! First tolerance in case of non zero flux
    double C_alpha_n = 0.02;
    data.setValue(C_alpha_n);
    Property property_C_alpha_n("--C_alpha_n",data,Property::PropertyGroup_ConvergenceCriteria);
    props.append(property_C_alpha_n);

    //! Second tolerance in case of zero flux
    double C_alpha_epsilon = 0.001;
    data.setValue(C_alpha_epsilon);
    Property property_C_alpha_epsilon("--C_alpha_epsilon",data,Property::PropertyGroup_ConvergenceCriteria);
    props.append(property_C_alpha_epsilon);

    //! --------------------------------------------------------------------------------------------
    //! field parameters
    //! (R_alpha_n,C_alpha_n,q_alpha_0,q_alpha_u,R_alpha_P,epsilon_alpha,C_alpha_epsilon,R_alpha_l)
    //! --------------------------------------------------------------------------------------------
    QVector<double> fieldParameters;

    fieldParameters.append(R_alpha_n);
    fieldParameters.append(C_alpha_n);
    fieldParameters.append(q_alpha_0);
    fieldParameters.append(q_alpha_u);
    fieldParameters.append(R_alpha_P);
    fieldParameters.append(epsilon_alpha);
    fieldParameters.append(C_alpha_epsilon);
    fieldParameters.append(R_alpha_l);

    data.setValue(fieldParameters);
    QVector<QVariant> values;
    values.append(data);

    load loadFieldParameters(values,Property::loadType_fieldParameters);

    data.setValue(fluxConvergenceType);
    QVector<QVariant> values1;
    values1.append(data);
    load loadFluxConvergence(values1,Property::loadType_fluxConvergence);

    data.setValue(solutionConvergenceType);
    QVector<QVariant> values2;
    values2.append(data);
    load loadSolutionConvergence(values2,Property::loadType_solutionConvergence);

    //! --------------------------------------------
    //! line search
    //! "0" => "Program controlled" "1" => "Custom"
    //! --------------------------------------------
    int lineSearch = 0;
    data.setValue(lineSearch);
    Property property_lineSearch("Line search",data,Property::PropertyGroup_LineSearch);

    double minLineSearch = 0.25;
    data.setValue(minLineSearch);
    Property property_minLineSearch("Min value",data,Property::PropertyGroup_LineSearch);

    double maxLineSearch = 1.01;
    data.setValue(maxLineSearch);
    Property property_maxLineSearch("Max value",data,Property::PropertyGroup_LineSearch);

    props.push_back(property_lineSearch);
    props.push_back(property_minLineSearch);
    props.push_back(property_maxLineSearch);

    QVector<QVariant> values7;
    data.setValue(lineSearch);
    values7.push_back(data);
    load loadLineSearch(values7,Property::loadType_lineSearch);

    QVector<QVariant> values8;
    QVector<double> lineSearchParameters;
    lineSearchParameters.push_back(minLineSearch);
    lineSearchParameters.push_back(maxLineSearch);
    data.setValue(lineSearchParameters);
    values8.push_back(data);
    load loadLineSearchParameters(values8,Property::loadType_lineSearchParameters);

    //! -------------------------------
    //! "Output settings"
    //! "0" => do not save "1" => save
    //! -------------------------------
    QVector<int> outputSettings;

    int saveStress = 1;
    data.setValue(saveStress);
    outputSettings.push_back(saveStress);
    Property property_stress("Stress",data,Property::PropertyGroup_OutputSettings);
    props.push_back(property_stress);

    int saveStrain = 1;
    data.setValue(saveStrain);
    outputSettings.push_back(saveStrain);
    Property property_strain("Strain",data,Property::PropertyGroup_OutputSettings);
    props.push_back(property_strain);

    int saveReactionForces = 1;
    data.setValue(saveReactionForces);
    outputSettings.push_back(saveReactionForces);
    Property property_RF("Reaction forces",data,Property::PropertyGroup_OutputSettings);
    props.push_back(property_RF);

    int saveContactData = 1;
    data.setValue(saveContactData);
    outputSettings.push_back(saveContactData);
    Property property_contactData("Contact data",data,Property::PropertyGroup_OutputSettings);
    props.push_back(property_contactData);

    data.setValue(outputSettings);
    QVector<QVariant> values9; values9.push_back(data);
    load loadOuputSettings(values9,Property::loadType_outputSettings);

    //! -------------------------------------
    //! options for storing analysis results
    //! "0" => "All time points" (defaults)
    //! "1" => "Last time point"
    //! "2" => "Specified recurrence rate"
    //! -------------------------------------
    int storeResultsAt = 0;
    data.setValue(storeResultsAt);
    Property property_storeResultsAt("Store results at",data,Property::PropertyGroup_OutputSettings);
    props.push_back(property_storeResultsAt);

    int FREQUENCY = 1;
    QVector<int> storeResultsAtFlags;
    storeResultsAtFlags.push_back(storeResultsAt);
    storeResultsAtFlags.push_back(FREQUENCY);
    data.setValue(storeResultsAtFlags);
    QVector<QVariant> values10; values10.push_back(data);
    load loadStoreResultsAt(values10,Property::loadType_storeResultsAt);

    //! ------------------------------------
    //! create the "Analysis settings" node
    //! ------------------------------------
    SimulationNodeClass *nodeCombinedAnalysisSettings =
            new SimulationNodeClass("Analysis settings",SimulationNodeClass::nodeType_combinedAnalysisSettings,props,this);

    //! -------------------------------------------------
    //! create the tabular data for the time step policy
    //! -------------------------------------------------
    QVector<QVariant> stepNumbers;
    QVector<QVariant> stepEndTimes;
    QVector<QVariant> solverTypes;
    QVector<QVariant> vecSubsteps;

    //! time step number - a default value
    data.setValue(1);
    stepNumbers.push_back(data);

    //! Step end time - a default value
    double aStepEndTime = 1.0;
    data.setValue(aStepEndTime);
    stepEndTimes.push_back(data);

    //! Solver type - a default value
    data.setValue(Property::solverType_programControlled);
    solverTypes.push_back(data);

    //! ----------------------------------------------------------------------
    //! time step divisions
    //! the last value is for the time step policy (Auto = 0, ON = 1, OFF =2)
    //! ----------------------------------------------------------------------
    QVector<int> T;
    T.append(5); T.append(2); T.append(10);
    //! the time step policy
    T.append(0);

    data.setValue(T);
    vecSubsteps.append(data);

    load aLoad1(stepNumbers,Property::loadType_stepNumber);
    load aLoad2(stepEndTimes,Property::loadType_stepEndTime);
    load aLoad3(solverTypes,Property::loadType_solverType);
    load aLoad4(vecSubsteps,Property::loadType_autoTimeStepping);

    //! --------------
    //! analysis type
    //! --------------
    data.setValue(analysisType);
    QVector<QVariant> vecAnalysisType;
    vecAnalysisType.push_back(data);
    load loadAnalysisType(vecAnalysisType,Property::loadType_analysisType);

    //! --------------------------
    //! static/transient in table
    //! --------------------------
    data.setValue(timeIntegration);
    QVector<QVariant> vecTimeIntegration;
    vecTimeIntegration.push_back(data);
    load loadTimeIntegration(vecTimeIntegration,Property::loadType_timeIntegration);

    //! -------------------------
    //! build the internal table
    //! -------------------------
    QVector<load> vecLoad;

    //! 1-st time step numner (column # 0)
    vecLoad.push_back(aLoad1);
    //! 2-nd step end time (column # 1)
    vecLoad.push_back(aLoad2);
    //! 3-rd column solver type (column # 2)
    vecLoad.push_back(aLoad3);
    //! 4-th column time step policy (column # 3)
    vecLoad.push_back(aLoad4);
    //! 5-th column (column # 4)
    vecLoad.push_back(loadFieldParameters);
    //! 6-th column (column # 5)
    vecLoad.push_back(loadFluxConvergence);
    //! 7-th column (column # 6)
    vecLoad.push_back(loadSolutionConvergence);
    //! 8-th column (column # 7)
    vecLoad.push_back(loadTimeIncrementation);
    //! 9-th column (column # 8)
    vecLoad.push_back(loadTimeIncrementationType);
    //! 10-th column (column # 9)
    vecLoad.push_back(loadCutBackType);
    //! 11-th column (column #10)
    vecLoad.push_back(loadCutBackParameters);
    //! 12-th column (column #11)
    vecLoad.push_back(loadLineSearch);
    //! 13-th column (column #12)
    vecLoad.push_back(loadLineSearchParameters);
    //! 14-th column (column #13)
    vecLoad.push_back(loadOuputSettings);
    //! 15-th column (column #14)
    vecLoad.push_back(loadStoreResultsAt);
    //! 16-th column (column #15)
    vecLoad.push_back(loadAnalysisType);
    //! 17-th column (column #16)
    vecLoad.push_back(loadTimeIntegration);

    //! ------------------------
    //! create the tabular data
    //! ------------------------
    nodeCombinedAnalysisSettings->createTabularData(vecLoad);

    //! -----------------------------
    //! time tag and parent time tag
    //! -----------------------------
    nodeCombinedAnalysisSettings->addTimeTag();
    QString parentTimeTag = CombinedAnalysisItem->data(Qt::UserRole).value<SimulationNodeClass*>()->getPropertyValue<QString>("Time tag");
    data.setValue(parentTimeTag);
    nodeCombinedAnalysisSettings->addProperty(Property("Parent time tag",data,Property::PropertyGroup_Identifier));

    //! ------------------------------------------------------
    //! create the item and append to the "Combined analysis"
    //! ------------------------------------------------------
    data.setValue(nodeCombinedAnalysisSettings);
    CombinedAnalysisSettingsItem->setData(data,Qt::UserRole);
    CombinedAnalysisSettingsItem->setData("Analysis settings", Qt::DisplayRole);
    CombinedAnalysisItem->appendRow(CombinedAnalysisSettingsItem);

    mainTreeTools::addSolution(CombinedAnalysisItem);
    mainTreeTools::addSolutionInformation(CombinedAnalysisItem->child(CombinedAnalysisItem->rowCount()-1,0));
}

//! --------------------------------
//! function: createRemotePointRoot
//! details:
//! --------------------------------
void simulationDataBase::createRemotePointRoot()
{
    cout<<"simulationDataBase::createRemotePointRoot()->____function called____"<<endl;
    QVariant data;
    QVector<Property> props;
    SimulationNodeClass *node = new SimulationNodeClass("Remote point", SimulationNodeClass::nodeType_remotePointRoot, props, this);

    //! -----------------------------
    //! time tag and parent time tag
    //! -----------------------------
    node->addTimeTag();
    QString timeTag = myRootItem->data(Qt::UserRole).value<SimulationNodeClass*>()->getPropertyValue<QString>("Time tag");
    data.setValue(timeTag);
    node->addProperty(Property("Parent time tag",data,Property::PropertyGroup_Identifier));

    RemotePointsRootItem = new QExtendedStandardItem();
    data.setValue(node);
    RemotePointsRootItem->setData(data,Qt::UserRole);
    RemotePointsRootItem->setData("Remote points", Qt::DisplayRole);

    //! -----------------------------------------------
    //! create a dummy remote point "Select from list"
    //! -----------------------------------------------
    SimulationNodeClass *emptyRemotePointNode = nodeFactory::nodeFromScratch(SimulationNodeClass::nodeType_remotePoint,this,0,QVariant());
    emptyRemotePointNode->setName(QString("Select from list"));

    //! -----------------------------
    //! time tag and parent time tag
    //! -----------------------------
    emptyRemotePointNode->addTimeTag();
    timeTag = node->getPropertyValue<QString>("Time tag");
    data.setValue(timeTag);
    emptyRemotePointNode->addProperty(Property("Parent time tag",data,Property::PropertyGroup_Identifier));

    //! ---------------------------
    //! create and append the item
    //! ---------------------------
    QExtendedStandardItem *itemEmptyRemotePoint = new QExtendedStandardItem();
    itemEmptyRemotePoint->setData(QString("Select from list"),Qt::DisplayRole);
    data.setValue(emptyRemotePointNode);
    itemEmptyRemotePoint->setData(data,Qt::UserRole);
    RemotePointsRootItem->appendRow(itemEmptyRemotePoint);

    //! --------------------------------------------------------
    //! always insert after "Geometry" and "Coordinate systems"
    //! --------------------------------------------------------
    myRootItem->insertRow(GeometryRootItem->row()+1,RemotePointsRootItem);
}

//! -------------------------------------------
//! function: createStructuralAnalysisRootNode
//! details:
//! -------------------------------------------
void simulationDataBase::createStructuralAnalysisRootNode()
{
    cout<<"simulationDataBase::createStructuralAnalysisRootNode()->____creating a structural analysis branch____"<<endl;

    QVariant data;
    QString name = "Static structural";

    //! ---------------
    //! the properties
    //! ---------------
    QVector<Property> props;

    //! ------------------------------------------------
    //! a default value for the environment temperature
    //! ------------------------------------------------
    double Tenv = 295;
    Property prop_envTemp("Environment temperature",QVariant(Tenv),Property::PropertyGroup_EnvironmentTemperature);
    props.push_back(prop_envTemp);

    //! --------------
    //! node creation
    //! --------------
    SimulationNodeClass *nodeStructuralAnalysis = new SimulationNodeClass(name,SimulationNodeClass::nodeType_structuralAnalysis,props,this);

    //! -----------------------------
    //! time tag and parent time tag
    //! -----------------------------
    nodeStructuralAnalysis->addTimeTag();
    QString timeTag = myRootItem->data(Qt::UserRole).value<SimulationNodeClass*>()->getPropertyValue<QString>("Time tag");
    data.setValue(timeTag);
    nodeStructuralAnalysis->addProperty(Property("Parent time tag",data,Property::PropertyGroup_Identifier));

    //! --------------------------------------------------
    //! create the item and append to the model root item
    //! --------------------------------------------------
    QExtendedStandardItem *StructuralAnalysisItem = new QExtendedStandardItem();
    StructuralAnalysisItem->setData("Static structural",Qt::DisplayRole);
    data.setValue(nodeStructuralAnalysis);
    StructuralAnalysisItem->setData(data, Qt::UserRole);
    myRootItem->appendRow(StructuralAnalysisItem);

    //! -----------------------------------------------
    //! create the item structural "Analysis settings"
    //! -----------------------------------------------
    QExtendedStandardItem *StructuralAnalysisSettingsItem = new QExtendedStandardItem();

    //! ------------------
    //! "Number of steps"
    //! ------------------
    data.setValue(1);
    Property property_numberOfSteps("Number of steps",data,Property::PropertyGroup_StepControls);

    //! ----------------------
    //! "Current step number"
    //! ----------------------
    data.setValue(1);
    Property property_currentStepNumber("Current step number",data,Property::PropertyGroup_StepControls);

    //! ------------------------------
    //! the default analysis end time
    //! ------------------------------
    double endTime =1.0;
    data.setValue(endTime);
    Property property_stepEndTime("Step end time",data,Property::PropertyGroup_StepControls);

    //! --------------
    //! analysis type
    //! --------------
    Property::analysisType analysisType = Property::analysisType_structural;
    data.setValue(analysisType);
    Property property_analysisType("Analysis type",data,Property::PropertyGroup_StepControls);

    //! -----------------------
    //! steady state/transient
    //! -----------------------
    Property::timeIntegration timeIntegration = Property::timeIntegration_steadyState;
    data.setValue(timeIntegration);
    Property property_timeIntegration("Static/Transient",data,Property::PropertyGroup_StepControls);

    //! --------------
    //! time stepping
    //! --------------
    Property::autoTimeStepping theTimeStepping = Property::autoTimeStepping_ProgramControlled;
    data.setValue(theTimeStepping);
    Property property_autoTimeStepping("Auto time stepping",data,Property::PropertyGroup_StepControls);

    //! --------------
    //! "Solver type"
    //! --------------
    Property::solverType theSolverType = Property::solverType_programControlled;
    data.setValue(theSolverType);
    Property property_solverType("Solver type",data,Property::PropertyGroup_SolverControls);

    //! ---------------------------
    //! non linear geometry
    //! 0 -> NL "Off" 1 -> NL "On"
    //! ---------------------------
    Property property_NLgeometry("Large deflection",int(0),Property::PropertyGroup_SolverControls);

    //! ------------------
    //! number of threads
    //! ------------------
    int numberOfThreads = 4;
    data.setValue(numberOfThreads);
    Property property_numberOfThreads("Number of threads",data,Property::PropertyGroup_SolverControls);

    props.push_back(property_numberOfSteps);
    props.push_back(property_currentStepNumber);
    props.push_back(property_stepEndTime);
    props.push_back(property_analysisType);
    props.push_back(property_timeIntegration);
    props.push_back(property_autoTimeStepping);
    props.push_back(property_solverType);
    props.push_back(property_numberOfThreads);
    props.push_back(property_NLgeometry);

    //! -------------------------------------
    //! Time incrementation
    //! 0 = > Custom 1 => program controlled
    //! -------------------------------------
    int timeIncrementation = 0;
    data.setValue(timeIncrementation);
    Property property_timeIncrementation("Time incrementation",data,Property::PropertyGroup_TimeIncrementation);
    props.push_back(property_timeIncrementation);

    data.setValue(timeIncrementation);
    QVector<QVariant> values4;
    values4.push_back(data);
    load loadTimeIncrementationType(values4,Property::loadType_timeIncrementation);

    int I_0 = 4;
    int I_R = 8;
    int I_P = 9;
    int I_C = 16;
    int I_L = 10;
    int I_G = 4;
    int I_S = -1;       //! "-1" means "unused"
    int I_A = 5;
    int I_J = -1;       //! "-1" means "unused"
    int I_T = -1;       //! "-1" means "unused"

    data.setValue(I_0);
    Property property_I_0("I_0",data,Property::PropertyGroup_TimeIncrementation);
    props.push_back(property_I_0);
    data.setValue(I_R);
    Property property_I_R("I_R",data,Property::PropertyGroup_TimeIncrementation);
    props.push_back(property_I_R);
    data.setValue(I_P);
    Property property_I_P("I_P",data,Property::PropertyGroup_TimeIncrementation);
    props.push_back(property_I_P);
    data.setValue(I_C);
    Property property_I_C("I_C",data,Property::PropertyGroup_TimeIncrementation);
    props.push_back(property_I_C);
    data.setValue(I_L);
    Property property_I_L("I_L",data,Property::PropertyGroup_TimeIncrementation);
    props.push_back(property_I_L);
    data.setValue(I_G);
    Property property_I_G("I_G",data,Property::PropertyGroup_TimeIncrementation);
    props.push_back(property_I_G);
    data.setValue(I_S);
    Property property_I_S("I_S",data,Property::PropertyGroup_TimeIncrementation);
    props.push_back(property_I_S);
    data.setValue(I_A);
    Property property_I_A("I_A",data,Property::PropertyGroup_TimeIncrementation);
    props.push_back(property_I_A);
    data.setValue(I_J);
    Property property_I_J("I_J",data,Property::PropertyGroup_TimeIncrementation);
    props.push_back(property_I_J);
    data.setValue(I_T);
    Property property_I_T("I_T",data,Property::PropertyGroup_TimeIncrementation);
    props.push_back(property_I_T);

    QVector<int> timeIncrementationParameters;
    timeIncrementationParameters.push_back(I_0);
    timeIncrementationParameters.push_back(I_R);
    timeIncrementationParameters.push_back(I_P);
    timeIncrementationParameters.push_back(I_C);
    timeIncrementationParameters.push_back(I_L);
    timeIncrementationParameters.push_back(I_G);
    timeIncrementationParameters.push_back(I_S);
    timeIncrementationParameters.push_back(I_A);
    timeIncrementationParameters.push_back(I_J);
    timeIncrementationParameters.push_back(I_T);

    data.setValue(timeIncrementationParameters);
    QVector<QVariant> values3;
    values3.append(data);
    load loadTimeIncrementation(values3,Property::loadType_timeIncrementationParameters);

    //! ----------------------------------------------------------
    //! cutback type: "0" => "Program controlled" "1" => "Custom"
    //! ----------------------------------------------------------
    int typeOfCutback = 0;
    data.setValue(typeOfCutback);
    Property property_cutBackType("Cutback factors",data,Property::PropertyGroup_CutBack);
    props.push_back(property_cutBackType);

    data.setValue(typeOfCutback);
    QVector<QVariant> values5;
    values5.push_back(data);
    load loadCutBackType(values5,Property::loadType_cutBack);

    //! -------------------------------------
    //! cutback factors: CCX solver defaluts
    //! -------------------------------------
    double D_f = 0.25;
    double D_C = 0.5;
    double D_B = 0.75;
    double D_A = 0.85;
    double D_S = -1.0;      //! "-1.0" means "currently unused"
    double D_H = -1.0;      //! "-1.0" means "currently unused"
    double D_D = 1.5;
    double W_G = -1.0;      //! "-1.0" means "currently unused"

    data.setValue(D_f);
    Property property_D_f("D_f",data,Property::PropertyGroup_CutBack);
    data.setValue(D_C);
    Property property_D_C("D_C",data,Property::PropertyGroup_CutBack);
    data.setValue(D_B);
    Property property_D_B("D_B",data,Property::PropertyGroup_CutBack);
    data.setValue(D_A);
    Property property_D_A("D_A",data,Property::PropertyGroup_CutBack);
    data.setValue(D_S);
    Property property_D_S("D_S",data,Property::PropertyGroup_CutBack);
    data.setValue(D_H);
    Property property_D_H("D_H",data,Property::PropertyGroup_CutBack);
    data.setValue(D_D);
    Property property_D_D("D_D",data,Property::PropertyGroup_CutBack);
    data.setValue(W_G);
    Property property_W_G("W_G",data,Property::PropertyGroup_CutBack);

    props.push_back(property_D_f);
    props.push_back(property_D_C);
    props.push_back(property_D_B);
    props.push_back(property_D_A);
    props.push_back(property_D_S);
    props.push_back(property_D_H);
    props.push_back(property_D_D);
    props.push_back(property_W_G);

    QVector<double> cutBackParameters;
    cutBackParameters.push_back(D_f);
    cutBackParameters.push_back(D_C);
    cutBackParameters.push_back(D_B);
    cutBackParameters.push_back(D_A);
    cutBackParameters.push_back(D_S);
    cutBackParameters.push_back(D_H);
    cutBackParameters.push_back(D_D);
    cutBackParameters.push_back(W_G);
    data.setValue(cutBackParameters);
    QVector<QVariant> values6;
    values6.push_back(data);
    load loadCutBackParameters(values6,Property::loadType_cutBackParameters);

    //! ---------------------------------------------------------
    //! Convergence criteria - "Flux convergence"
    //! 0 = > Deactivated 1 => On 2 => Program controlled
    //! The item is created with the option "2" using default
    //! values for the parameters
    //! ---------------------------------------------------------
    int fluxConvergenceType = 2;
    data.setValue(fluxConvergenceType);
    Property property_fluxConvergence("Flux convergence",data,Property::PropertyGroup_ConvergenceCriteria);
    props.append(property_fluxConvergence);

    //! q_alpha_u   => "0" is Program controlled
    double q_alpha_u = 0.0;
    data.setValue(q_alpha_u);
    Property property_value("--Value",data,Property::PropertyGroup_ConvergenceCriteria);
    props.append(property_value);

    //! q_alpha_0   => "0" is Program controlled
    double q_alpha_0 = 0.0;
    data.setValue(q_alpha_0);
    Property property_q_alpha_0("--q_alpha_0",data,Property::PropertyGroup_ConvergenceCriteria);
    props.append(property_q_alpha_0);

    //! First tolerance in case of non zero flux
    double R_alpha_n = 0.001;
    data.setValue(R_alpha_n);
    Property property_R_alpha_n("--R_alpha_n",data,Property::PropertyGroup_ConvergenceCriteria);
    props.append(property_R_alpha_n);

    //! Second tolerance in case of non zero flux
    double R_alpha_P = 0.02;
    data.setValue(R_alpha_P);
    Property property_R_alpha_P("--R_alpha_P",data,Property::PropertyGroup_ConvergenceCriteria);
    props.append(property_R_alpha_P);

    //! Third tolerance parameter: R_alpha_l
    double R_alpha_l = 1e-8;
    data.setValue(R_alpha_l);
    Property property_R_alpha_l("--R_alpha_l",data,Property::PropertyGroup_ConvergenceCriteria);
    props.append(property_R_alpha_l);

    //! zero flux parameter
    double epsilon_alpha = 1e-5;
    data.setValue(epsilon_alpha);
    Property property_epsilon_n("--epsilon_alpha",data,Property::PropertyGroup_ConvergenceCriteria);
    props.append(property_epsilon_n);

    //! --------------------------------------------------
    //! Convergence criteria: "Solution convergence"
    //! 0 = > Deactivated 1 => On 2 => Program controlled
    //! --------------------------------------------------
    int solutionConvergenceType = 2;
    data.setValue(solutionConvergenceType);
    Property property_solutionConvergence("Solution convergence",data,Property::PropertyGroup_ConvergenceCriteria);
    props.append(property_solutionConvergence);

    //! First tolerance in case of non zero flux
    double C_alpha_n = 0.02;
    data.setValue(C_alpha_n);
    Property property_C_alpha_n("--C_alpha_n",data,Property::PropertyGroup_ConvergenceCriteria);
    props.append(property_C_alpha_n);

    //! Second tolerance in case of zero flux
    double C_alpha_epsilon = 0.001;
    data.setValue(C_alpha_epsilon);
    Property property_C_alpha_epsilon("--C_alpha_epsilon",data,Property::PropertyGroup_ConvergenceCriteria);
    props.append(property_C_alpha_epsilon);

    //! --------------------------------------------------------------------------------------------
    //! field parameters
    //! (R_alpha_n,C_alpha_n,q_alpha_0,q_alpha_u,R_alpha_P,epsilon_alpha,C_alpha_epsilon,R_alpha_l)
    //! --------------------------------------------------------------------------------------------
    QVector<double> fieldParameters;

    fieldParameters.append(R_alpha_n);
    fieldParameters.append(C_alpha_n);
    fieldParameters.append(q_alpha_0);
    fieldParameters.append(q_alpha_u);
    fieldParameters.append(R_alpha_P);
    fieldParameters.append(epsilon_alpha);
    fieldParameters.append(C_alpha_epsilon);
    fieldParameters.append(R_alpha_l);

    data.setValue(fieldParameters);
    QVector<QVariant> values;
    values.append(data);

    load loadFieldParameters(values,Property::loadType_fieldParameters);

    data.setValue(fluxConvergenceType);
    QVector<QVariant> values1;
    values1.append(data);
    load loadFluxConvergence(values1,Property::loadType_fluxConvergence);

    data.setValue(solutionConvergenceType);
    QVector<QVariant> values2;
    values2.append(data);
    load loadSolutionConvergence(values2,Property::loadType_solutionConvergence);

    //! --------------------------------------------
    //! line search
    //! "0" => "Program controlled" "1" => "Custom"
    //! --------------------------------------------
    int lineSearch = 0;
    data.setValue(lineSearch);
    Property property_lineSearch("Line search",data,Property::PropertyGroup_LineSearch);

    double minLineSearch = 0.25;
    data.setValue(minLineSearch);
    Property property_minLineSearch("Min value",data,Property::PropertyGroup_LineSearch);

    double maxLineSearch = 1.01;
    data.setValue(maxLineSearch);
    Property property_maxLineSearch("Max value",data,Property::PropertyGroup_LineSearch);

    props.push_back(property_lineSearch);
    props.push_back(property_minLineSearch);
    props.push_back(property_maxLineSearch);

    QVector<QVariant> values7;
    data.setValue(lineSearch);
    values7.push_back(data);
    load loadLineSearch(values7,Property::loadType_lineSearch);

    QVector<QVariant> values8;
    QVector<double> lineSearchParameters;
    lineSearchParameters.push_back(minLineSearch);
    lineSearchParameters.push_back(maxLineSearch);
    data.setValue(lineSearchParameters);
    values8.push_back(data);
    load loadLineSearchParameters(values8,Property::loadType_lineSearchParameters);

    //! -------------------------------
    //! "Output settings"
    //! "0" => do not save "1" => save
    //! -------------------------------
    QVector<int> outputSettings;

    int saveStress = 1;
    data.setValue(saveStress);
    outputSettings.push_back(saveStress);
    Property property_stress("Stress",data,Property::PropertyGroup_OutputSettings);
    props.push_back(property_stress);

    int saveStrain = 1;
    data.setValue(saveStrain);
    outputSettings.push_back(saveStrain);
    Property property_strain("Strain",data,Property::PropertyGroup_OutputSettings);
    props.push_back(property_strain);

    int saveReactionForces = 1;
    data.setValue(saveReactionForces);
    outputSettings.push_back(saveReactionForces);
    Property property_RF("Reaction forces",data,Property::PropertyGroup_OutputSettings);
    props.push_back(property_RF);

    int saveContactData = 1;
    data.setValue(saveContactData);
    outputSettings.push_back(saveContactData);
    Property property_contactData("Contact data",data,Property::PropertyGroup_OutputSettings);
    props.push_back(property_contactData);

    data.setValue(outputSettings);
    QVector<QVariant> values9; values9.push_back(data);
    load loadOuputSettings(values9,Property::loadType_outputSettings);

    //! -------------------------------------
    //! options for storing analysis results
    //! "0" => "All time points" (defaults)
    //! "1" => "Last time point"
    //! "2" => "Specified recurrence rate"
    //! -------------------------------------
    int storeResultsAt = 0;
    data.setValue(storeResultsAt);
    Property property_storeResultsAt("Store results at",data,Property::PropertyGroup_OutputSettings);
    props.push_back(property_storeResultsAt);

    int FREQUENCY = 1;
    QVector<int> storeResultsAtFlags;
    storeResultsAtFlags.push_back(storeResultsAt);
    storeResultsAtFlags.push_back(FREQUENCY);
    data.setValue(storeResultsAtFlags);
    QVector<QVariant> values10; values10.push_back(data);
    load loadStoreResultsAt(values10,Property::loadType_storeResultsAt);

    //! ------------------------------------
    //! create the "Analysis settings" node
    //! ------------------------------------
    SimulationNodeClass *nodeStructuralAnalysisSettings =
            new SimulationNodeClass("Analysis settings",SimulationNodeClass::nodeType_structuralAnalysisSettings,props,this);

    //! -------------------------------------------------
    //! create the tabular data for the time step policy
    //! -------------------------------------------------
    QVector<QVariant> stepNumbers;
    QVector<QVariant> stepEndTimes;
    QVector<QVariant> solverTypes;
    QVector<QVariant> vecSubsteps;

    //! time step number - a default value
    data.setValue(1);
    stepNumbers.push_back(data);

    //! Step end time - a default value
    double aStepEndTime = 1.0;
    data.setValue(aStepEndTime);
    stepEndTimes.push_back(data);

    //! Solver type - a default value
    data.setValue(Property::solverType_programControlled);
    solverTypes.push_back(data);

    //! ----------------------------------------------------------------------
    //! time step divisions
    //! the last value is for the time step policy (Auto = 0, ON = 1, OFF =2)
    //! ----------------------------------------------------------------------
    QVector<int> T;
    T.append(5); T.append(2); T.append(10);
    //! the time step policy
    T.append(0);

    data.setValue(T);
    vecSubsteps.append(data);

    load aLoad1(stepNumbers,Property::loadType_stepNumber);
    load aLoad2(stepEndTimes,Property::loadType_stepEndTime);
    load aLoad3(solverTypes,Property::loadType_solverType);
    load aLoad4(vecSubsteps,Property::loadType_autoTimeStepping);

    //! --------------
    //! analysis type
    //! --------------
    data.setValue(analysisType);
    QVector<QVariant> vecAnalysisType;
    vecAnalysisType.push_back(data);
    load loadAnalysisType(vecAnalysisType,Property::loadType_analysisType);

    //! --------------------------
    //! static/transient in table
    //! --------------------------
    data.setValue(timeIntegration);
    QVector<QVariant> vecTimeIntegration;
    vecTimeIntegration.push_back(data);
    load loadTimeIntegration(vecTimeIntegration,Property::loadType_timeIntegration);

    //! -------------------------
    //! build the internal table
    //! -------------------------
    QVector<load> vecLoad;

    //! 1-st time step numner (column # 0)
    vecLoad.push_back(aLoad1);
    //! 2-nd step end time (column # 1)
    vecLoad.push_back(aLoad2);
    //! 3-rd column solver type (column # 2)
    vecLoad.push_back(aLoad3);
    //! 4-th column time step policy (column # 3)
    vecLoad.push_back(aLoad4);
    //! 5-th column (column # 4)
    vecLoad.push_back(loadFieldParameters);
    //! 6-th column (column # 5)
    vecLoad.push_back(loadFluxConvergence);
    //! 7-th column (column # 6)
    vecLoad.push_back(loadSolutionConvergence);
    //! 8-th column (column # 7)
    vecLoad.push_back(loadTimeIncrementation);
    //! 9-th column (column # 8)
    vecLoad.push_back(loadTimeIncrementationType);
    //! 10-th column (column # 9)
    vecLoad.push_back(loadCutBackType);
    //! 11-th column (column #10)
    vecLoad.push_back(loadCutBackParameters);
    //! 12-th column (column #11)
    vecLoad.push_back(loadLineSearch);
    //! 13-th column (column #12)
    vecLoad.push_back(loadLineSearchParameters);
    //! 14-th column (column #13)
    vecLoad.push_back(loadOuputSettings);
    //! 15-th column (column #14)
    vecLoad.push_back(loadStoreResultsAt);
    //! 16-th column (column #15)
    vecLoad.push_back(loadAnalysisType);
    //! 17-th column (column #16)
    vecLoad.push_back(loadTimeIntegration);

    //! ------------------------
    //! create the tabular data
    //! ------------------------
    nodeStructuralAnalysisSettings->createTabularData(vecLoad);

    //! -----------------------------
    //! time tag and parent time tag
    //! -----------------------------
    nodeStructuralAnalysisSettings->addTimeTag();
    QString parentTimeTag = StructuralAnalysisItem->data(Qt::UserRole).value<SimulationNodeClass*>()->getPropertyValue<QString>("Time tag");
    data.setValue(parentTimeTag);
    nodeStructuralAnalysisSettings->addProperty(Property("Parent time tag",data,Property::PropertyGroup_Identifier));

    //! --------------------------------------------------------
    //! create the item and append to the "Structural analysis"
    //! --------------------------------------------------------
    data.setValue(nodeStructuralAnalysisSettings);
    StructuralAnalysisSettingsItem->setData(data,Qt::UserRole);
    StructuralAnalysisSettingsItem->setData("Analysis settings", Qt::DisplayRole);
    StructuralAnalysisItem->appendRow(StructuralAnalysisSettingsItem);

#ifdef COSTAMP_VERSION

    //! --------------------
    //! "Time step builder"
    //! --------------------
    props.clear();

    name = "Time step builder";

    //! props suppression status
    data.setValue(Property::SuppressionStatus_Active);
    Property prop_ss("Suppressed",data,Property::PropertyGroup_Definition);
    props.push_back(prop_ss);

    //! props time hystory file
    data.setValue(QString(""));
    Property prop_timeHystory("Time history file",data,Property::PropertyGroup_Definition);
    props.push_back(prop_timeHystory);

    SimulationNodeClass *node = new SimulationNodeClass(name,SimulationNodeClass::nodeType_timeStepBuilder,props);

    //! -----------------------------
    //! time tag and parent time tag
    //! -----------------------------
    node->addTimeTag();
    parentTimeTag = StructuralAnalysisItem->data(Qt::UserRole).value<SimulationNodeClass*>()->getPropertyValue<QString>("Time tag");
    data.setValue(parentTimeTag);
    node->addProperty(Property("Parent time tag",data,Property::PropertyGroup_Identifier));

    //! ------------------------------------
    //! create the "Time step builder" item
    //! ------------------------------------
    QExtendedStandardItem *item_timeStepBuilder = new QExtendedStandardItem();
    item_timeStepBuilder->setData(QString("Time step builder"),Qt::DisplayRole);
    item_timeStepBuilder->setEditable(false);
    data.setValue(node);
    item_timeStepBuilder->setData(data,Qt::UserRole);
    StructuralAnalysisItem->appendRow(item_timeStepBuilder);

#endif

    //! ---------------------------
    //! create the "Status" item
    //! ---------------------------
    props.clear();
    Property::solutionInformation theSolutionStatus = Property::solutionInformation_solveRequired;
    data.setValue(theSolutionStatus);
    Property prop_status("Status",data,Property::PropertyGroup_Information);
    props.push_back(prop_status);

    //! --------------------
    //! "Project files dir"
    //! --------------------
    data.setValue(QString("undefined"));
    Property prop_solutionFileDir("Project files dir",data,Property::PropertyGroup_Information);
    props.push_back(prop_solutionFileDir);
    SimulationNodeClass *nodeSolution = new SimulationNodeClass("Solution",SimulationNodeClass::nodeType_StructuralAnalysisSolution,props,this);
    data.setValue(nodeSolution);

    //! -----------------------------
    //! time tag and parent time tag
    //! -----------------------------
    nodeSolution->addTimeTag();
    parentTimeTag = StructuralAnalysisItem->data(Qt::UserRole).value<SimulationNodeClass*>()->getPropertyValue<QString>("Time tag");
    data.setValue(parentTimeTag);
    nodeSolution->addProperty(Property("Parent time tag",data,Property::PropertyGroup_Identifier));

    //! ----------------
    //! create the item
    //! ----------------
    QExtendedStandardItem *SolutionItem = new QExtendedStandardItem();
    data.setValue(nodeSolution);
    SolutionItem->setData(data,Qt::UserRole);
    SolutionItem->setData("Solution",Qt::DisplayRole);

    StructuralAnalysisItem->appendRow(SolutionItem);

    //! ---------------------------------------
    //! create the item "Solution information"
    //! ---------------------------------------
    props.clear();

    //! ------------------------------
    //! 0 => solver txt messages
    //! 1 => force convergence
    //! 2 => displacement convergence
    //! 3 => line search
    //! 4 => time step size
    //! ------------------------------
    data.setValue(int(0));
    Property prop_solutionInfoType("Solution information",data,Property::PropertyGroup_SolutionInfo);
    double updateInterval = 2.5;
    data.setValue(updateInterval);
    Property prop_updateInterval("Update interval",data,Property::PropertyGroup_SolutionInfo);

    props.push_back(prop_solutionInfoType);
    props.push_back(prop_updateInterval);

    //! ---------------------------
    //! hidden: solver text output
    //! ---------------------------
    CCXSolverMessage msg;
    data.setValue(msg);
    Property prop_solverOutput("Solver output",data,Property::PropertyGroup_Hidden);
    props.push_back(prop_solverOutput);

    //! ------------------
    //! discrete time map
    //! ------------------
    QMap<double,QVector<int>> dtm;
    data.setValue(dtm);
    Property prop_discreteTimeMap("Discrete time map",data,Property::PropertyGroup_Hidden);
    props.push_back(prop_discreteTimeMap);

    //! -----------------
    //! convergence data
    //! -----------------
    QList<solutionInfo> solInfoList;
    data.setValue(solInfoList);
    Property prop_solInfoList("Convergence data",data,Property::PropertyGroup_Hidden);
    props.push_back(prop_solInfoList);

    SimulationNodeClass *solInfoNode =
            new SimulationNodeClass("Solution information",SimulationNodeClass::nodeType_StructuralAnalysisSolutionInformation,props,this);

    //! -----------------------------------------------------------
    //! time tag and parent time tag - the time tag property
    //! has been already created by the node factory for this node
    //! -----------------------------------------------------------
    solInfoNode->addTimeTag();
    parentTimeTag = nodeSolution->getPropertyValue<QString>("Time tag");
    data.setValue(parentTimeTag);
    solInfoNode->addProperty(Property("Parent time tag",data,Property::PropertyGroup_Identifier));

    //! -------------------------------------
    //! create the solution information item
    //! -------------------------------------
    QExtendedStandardItem *itemSolutionInformation = new QExtendedStandardItem;

    data.setValue(solInfoNode);
    itemSolutionInformation->setData(data,Qt::UserRole);
    data.setValue(solInfoNode->getName());
    itemSolutionInformation->setData(data,Qt::DisplayRole);
    SolutionItem->appendRow(itemSolutionInformation);

    cout<<"simulationDataBase::createStructuralAnalysisRootNode()->____creating a static structural: done____"<<endl;
}

//! ----------------------------------------
//! function: createThermalAnalysisRootNode
//! details:
//! ----------------------------------------
void simulationDataBase::createThermalAnalysisRootNode()
{
    cout<<"simulationDataBase::createThermalAnalysisRootNode()->____creating a thermal analysis branch____"<<endl;

    QVariant data;
    QString name = "Thermal analysis";

    //! ---------------
    //! the properties
    //! ---------------
    QVector<Property> props;

    //! ------------------------------------------------
    //! a default value for the environment temperature
    //! ------------------------------------------------
    double Tenv = 295;
    Property prop_envTemp("Environment temperature",QVariant(Tenv),Property::PropertyGroup_EnvironmentTemperature);
    props.push_back(prop_envTemp);

    //! --------------
    //! node creation
    //! --------------
    SimulationNodeClass *nodeThermalAnalysis = new SimulationNodeClass(name,SimulationNodeClass::nodeType_thermalAnalysis,props,this);

    //! -----------------------------
    //! time tag and parent time tag
    //! -----------------------------
    nodeThermalAnalysis->addTimeTag();
    QString timeTag = myRootItem->data(Qt::UserRole).value<SimulationNodeClass*>()->getPropertyValue<QString>("Time tag");
    data.setValue(timeTag);
    nodeThermalAnalysis->addProperty(Property("Parent time tag",data,Property::PropertyGroup_Identifier));

    //! --------------------------------------------------
    //! create the item and append to the model root item
    //! --------------------------------------------------
    QExtendedStandardItem *ThermalAnalysisItem = new QExtendedStandardItem();
    ThermalAnalysisItem->setData("Thermal analysis",Qt::DisplayRole);
    data.setValue(nodeThermalAnalysis);
    ThermalAnalysisItem->setData(data, Qt::UserRole);
    myRootItem->appendRow(ThermalAnalysisItem);

    //! ------------------------------------
    //! create the item "Analysis settings"
    //! ------------------------------------
    QExtendedStandardItem *ThermalAnalysisSettingsItem = new QExtendedStandardItem();

    //! ------------------
    //! "Number of steps"
    //! ------------------
    data.setValue(1);
    Property property_numberOfSteps("Number of steps",data,Property::PropertyGroup_StepControls);

    //! ----------------------
    //! "Current step number"
    //! ----------------------
    data.setValue(1);
    Property property_currentStepNumber("Current step number",data,Property::PropertyGroup_StepControls);

    //! ----------------
    //! "Step end time"
    //! ----------------
    double endTime =1.0;
    data.setValue(endTime);
    Property property_stepEndTime("Step end time",data,Property::PropertyGroup_StepControls);

    //! --------------
    //! analysis type
    //! --------------
    Property::analysisType analysisType = Property::analysisType_thermal;
    data.setValue(analysisType);
    Property property_analysisType("Analysis type",data,Property::PropertyGroup_StepControls);

    //! -----------------------
    //! steady state/transient
    //! -----------------------
    Property::timeIntegration timeIntegration = Property::timeIntegration_steadyState;
    data.setValue(timeIntegration);
    Property property_timeIntegration("Static/Transient",data,Property::PropertyGroup_StepControls);

    //! -------------
    //! time steping
    //! -------------
    Property::autoTimeStepping theTimeStepping = Property::autoTimeStepping_ProgramControlled;
    data.setValue(theTimeStepping);
    Property property_autoTimeStepping("Auto time stepping",data,Property::PropertyGroup_StepControls);

    //! ------------
    //! Solver type
    //! ------------
    Property::solverType theSolverType = Property::solverType_programControlled;
    data.setValue(theSolverType);
    Property property_solverType("Solver type",data,Property::PropertyGroup_SolverControls);

    //! -----------------------------
    //! non linear geometry
    //! 0 -> NL "Off" 1 -> NL "On"
    //! -----------------------------
    //Property property_NLgeometry("Large deflection",int(0),Property::PropertyGroup_SolverControls);

    //! ------------------
    //! number of threads
    //! ------------------
    int numberOfThreads = 4;
    data.setValue(numberOfThreads);
    Property property_numberOfThreads("Number of threads",data,Property::PropertyGroup_SolverControls);

    props.push_back(property_numberOfSteps);
    props.push_back(property_currentStepNumber);
    props.push_back(property_stepEndTime);
    props.push_back(property_analysisType);
    props.push_back(property_timeIntegration);
    props.push_back(property_autoTimeStepping);
    props.push_back(property_solverType);
    props.push_back(property_numberOfThreads);

    //! -------------------------------------
    //! Time incrementation
    //! 0 = > Custom 1 => program controlled
    //! -------------------------------------
    int timeIncrementation = 0;
    data.setValue(timeIncrementation);
    Property property_timeIncrementation("Time incrementation",data,Property::PropertyGroup_TimeIncrementation);
    props.push_back(property_timeIncrementation);

    data.setValue(timeIncrementation);
    QVector<QVariant> values4;
    values4.push_back(data);
    load loadTimeIncrementationType(values4,Property::loadType_timeIncrementation);

    int I_0 = 4;
    int I_R = 8;
    int I_P = 9;
    int I_C = 16;
    int I_L = 10;
    int I_G = 4;
    int I_S = -1;       //! "-1" means "unused"
    int I_A = 5;
    int I_J = -1;       //! "-1" means "unused"
    int I_T = -1;       //! "-1" means "unused"

    data.setValue(I_0);
    Property property_I_0("I_0",data,Property::PropertyGroup_TimeIncrementation);
    props.push_back(property_I_0);
    data.setValue(I_R);
    Property property_I_R("I_R",data,Property::PropertyGroup_TimeIncrementation);
    props.push_back(property_I_R);
    data.setValue(I_P);
    Property property_I_P("I_P",data,Property::PropertyGroup_TimeIncrementation);
    props.push_back(property_I_P);
    data.setValue(I_C);
    Property property_I_C("I_C",data,Property::PropertyGroup_TimeIncrementation);
    props.push_back(property_I_C);
    data.setValue(I_L);
    Property property_I_L("I_L",data,Property::PropertyGroup_TimeIncrementation);
    props.push_back(property_I_L);
    data.setValue(I_G);
    Property property_I_G("I_G",data,Property::PropertyGroup_TimeIncrementation);
    props.push_back(property_I_G);
    data.setValue(I_S);
    Property property_I_S("I_S",data,Property::PropertyGroup_TimeIncrementation);
    props.push_back(property_I_S);
    data.setValue(I_A);
    Property property_I_A("I_A",data,Property::PropertyGroup_TimeIncrementation);
    props.push_back(property_I_A);
    data.setValue(I_J);
    Property property_I_J("I_J",data,Property::PropertyGroup_TimeIncrementation);
    props.push_back(property_I_J);
    data.setValue(I_T);
    Property property_I_T("I_T",data,Property::PropertyGroup_TimeIncrementation);
    props.push_back(property_I_T);

    QVector<int> timeIncrementationParameters;
    timeIncrementationParameters.push_back(I_0);
    timeIncrementationParameters.push_back(I_R);
    timeIncrementationParameters.push_back(I_P);
    timeIncrementationParameters.push_back(I_C);
    timeIncrementationParameters.push_back(I_L);
    timeIncrementationParameters.push_back(I_G);
    timeIncrementationParameters.push_back(I_S);
    timeIncrementationParameters.push_back(I_A);
    timeIncrementationParameters.push_back(I_J);
    timeIncrementationParameters.push_back(I_T);

    data.setValue(timeIncrementationParameters);
    QVector<QVariant> values3;
    values3.append(data);
    load loadTimeIncrementation(values3,Property::loadType_timeIncrementationParameters);

    //! ----------------------------------------------------------
    //! cutback type: "0" => "Program controlled" "1" => "Custom"
    //! ----------------------------------------------------------
    int typeOfCutback = 0;
    data.setValue(typeOfCutback);
    Property property_cutBackType("Cutback factors",data,Property::PropertyGroup_CutBack);
    props.push_back(property_cutBackType);

    data.setValue(typeOfCutback);
    QVector<QVariant> values5;
    values5.push_back(data);
    load loadCutBackType(values5,Property::loadType_cutBack);

    //! -------------------------------------
    //! cutback factors: CCX solver defaluts
    //! -------------------------------------
    double D_f = 0.25;
    double D_C = 0.5;
    double D_B = 0.75;
    double D_A = 0.85;
    double D_S = -1.0;      //! "-1.0" means "currently unused"
    double D_H = -1.0;      //! "-1.0" means "currently unused"
    double D_D = 1.5;
    double W_G = -1.0;      //! "-1.0" means "currently unused"

    data.setValue(D_f);
    Property property_D_f("D_f",data,Property::PropertyGroup_CutBack);
    data.setValue(D_C);
    Property property_D_C("D_C",data,Property::PropertyGroup_CutBack);
    data.setValue(D_B);
    Property property_D_B("D_B",data,Property::PropertyGroup_CutBack);
    data.setValue(D_A);
    Property property_D_A("D_A",data,Property::PropertyGroup_CutBack);
    data.setValue(D_S);
    Property property_D_S("D_S",data,Property::PropertyGroup_CutBack);
    data.setValue(D_H);
    Property property_D_H("D_H",data,Property::PropertyGroup_CutBack);
    data.setValue(D_D);
    Property property_D_D("D_D",data,Property::PropertyGroup_CutBack);
    data.setValue(W_G);
    Property property_W_G("W_G",data,Property::PropertyGroup_CutBack);

    props.push_back(property_D_f);
    props.push_back(property_D_C);
    props.push_back(property_D_B);
    props.push_back(property_D_A);
    props.push_back(property_D_S);
    props.push_back(property_D_H);
    props.push_back(property_D_D);
    props.push_back(property_W_G);

    QVector<double> cutBackParameters;
    cutBackParameters.push_back(D_f);
    cutBackParameters.push_back(D_C);
    cutBackParameters.push_back(D_B);
    cutBackParameters.push_back(D_A);
    cutBackParameters.push_back(D_S);
    cutBackParameters.push_back(D_H);
    cutBackParameters.push_back(D_D);
    cutBackParameters.push_back(W_G);
    data.setValue(cutBackParameters);
    QVector<QVariant> values6;
    values6.push_back(data);
    load loadCutBackParameters(values6,Property::loadType_cutBackParameters);

    //! ---------------------------------------------------------
    //! Convergence criteria - "Flux convergence"
    //! 0 = > Deactivated 1 => On 2 => Program controlled
    //! The item is created with the option "2" using default
    //! values for the parameters
    //! ---------------------------------------------------------
    int fluxConvergenceType = 2;
    data.setValue(fluxConvergenceType);
    Property property_fluxConvergence("Flux convergence",data,Property::PropertyGroup_ConvergenceCriteria);
    props.append(property_fluxConvergence);

    //! q_alpha_u   => "0" is Program controlled
    double q_alpha_u = 0.0;
    data.setValue(q_alpha_u);
    Property property_value("--Value",data,Property::PropertyGroup_ConvergenceCriteria);
    props.append(property_value);

    //! q_alpha_0   => "0" is Program controlled
    double q_alpha_0 = 0.0;
    data.setValue(q_alpha_0);
    Property property_q_alpha_0("--q_alpha_0",data,Property::PropertyGroup_ConvergenceCriteria);
    props.append(property_q_alpha_0);

    //! First tolerance in case of non zero flux
    double R_alpha_n = 0.005;
    data.setValue(R_alpha_n);
    Property property_R_alpha_n("--R_alpha_n",data,Property::PropertyGroup_ConvergenceCriteria);
    props.append(property_R_alpha_n);

    //! Second tolerance in case of non zero flux
    double R_alpha_P = 0.02;
    data.setValue(R_alpha_P);
    Property property_R_alpha_P("--R_alpha_P",data,Property::PropertyGroup_ConvergenceCriteria);
    props.append(property_R_alpha_P);

    //! Third tolerance parameter: R_alpha_l
    double R_alpha_l = 1e-8;
    data.setValue(R_alpha_l);
    Property property_R_alpha_l("--R_alpha_l",data,Property::PropertyGroup_ConvergenceCriteria);
    props.append(property_R_alpha_l);

    //! zero flux parameter
    double epsilon_alpha = 1e-5;
    data.setValue(epsilon_alpha);
    Property property_epsilon_n("--epsilon_alpha",data,Property::PropertyGroup_ConvergenceCriteria);
    props.append(property_epsilon_n);

    //! --------------------------------------------------
    //! Convergence criteria: "Solution convergence"
    //! 0 = > Deactivated 1 => On 2 => Program controlled
    //! --------------------------------------------------
    int solutionConvergenceType = 2;
    data.setValue(solutionConvergenceType);
    Property property_solutionConvergence("Solution convergence",data,Property::PropertyGroup_ConvergenceCriteria);
    props.append(property_solutionConvergence);

    //! First tolerance in case of non zero flux
    double C_alpha_n = 0.02;
    data.setValue(C_alpha_n);
    Property property_C_alpha_n("--C_alpha_n",data,Property::PropertyGroup_ConvergenceCriteria);
    props.append(property_C_alpha_n);

    //! Second tolerance in case of zero flux
    double C_alpha_epsilon = 0.001;
    data.setValue(C_alpha_epsilon);
    Property property_C_alpha_epsilon("--C_alpha_epsilon",data,Property::PropertyGroup_ConvergenceCriteria);
    props.append(property_C_alpha_epsilon);

    //! --------------------------------------------------------------------------------------------
    //! field parameters
    //! (R_alpha_n,C_alpha_n,q_alpha_0,q_alpha_u,R_alpha_P,epsilon_alpha,C_alpha_epsilon,R_alpha_l)
    //! --------------------------------------------------------------------------------------------
    QVector<double> fieldParameters;

    fieldParameters.append(R_alpha_n);
    fieldParameters.append(C_alpha_n);
    fieldParameters.append(q_alpha_0);
    fieldParameters.append(q_alpha_u);
    fieldParameters.append(R_alpha_P);
    fieldParameters.append(epsilon_alpha);
    fieldParameters.append(C_alpha_epsilon);
    fieldParameters.append(R_alpha_l);

    data.setValue(fieldParameters);
    QVector<QVariant> values;
    values.append(data);

    load loadFieldParameters(values,Property::loadType_fieldParameters);

    data.setValue(fluxConvergenceType);
    QVector<QVariant> values1;
    values1.append(data);
    load loadFluxConvergence(values1,Property::loadType_fluxConvergence);

    data.setValue(solutionConvergenceType);
    QVector<QVariant> values2;
    values2.append(data);
    load loadSolutionConvergence(values2,Property::loadType_solutionConvergence);

    //! --------------------------------------------
    //! line search
    //! "0" => "Program controlled" "1" => "Custom"
    //! --------------------------------------------
    int lineSearch = 0;
    data.setValue(lineSearch);
    Property property_lineSearch("Line search",data,Property::PropertyGroup_LineSearch);

    double minLineSearch = 0.25;
    data.setValue(minLineSearch);
    Property property_minLineSearch("Min value",data,Property::PropertyGroup_LineSearch);

    double maxLineSearch = 1.01;
    data.setValue(maxLineSearch);
    Property property_maxLineSearch("Max value",data,Property::PropertyGroup_LineSearch);

    props.push_back(property_lineSearch);
    props.push_back(property_minLineSearch);
    props.push_back(property_maxLineSearch);

    QVector<QVariant> values7;
    data.setValue(lineSearch);
    values7.push_back(data);
    load loadLineSearch(values7,Property::loadType_lineSearch);

    QVector<QVariant> values8;
    QVector<double> lineSearchParameters;
    lineSearchParameters.push_back(minLineSearch);
    lineSearchParameters.push_back(maxLineSearch);
    data.setValue(lineSearchParameters);
    values8.push_back(data);
    load loadLineSearchParameters(values8,Property::loadType_lineSearchParameters);

    //! -------------------------------
    //! "Output settings"
    //! "0" => do not save "1" => save
    //! -------------------------------
    QVector<int> outputSettings;

    int saveStress = 1;
    data.setValue(saveStress);
    outputSettings.push_back(saveStress);
    Property property_stress("Stress",data,Property::PropertyGroup_OutputSettings);
    props.push_back(property_stress);

    int saveStrain = 1;
    data.setValue(saveStrain);
    outputSettings.push_back(saveStrain);
    Property property_strain("Strain",data,Property::PropertyGroup_OutputSettings);
    props.push_back(property_strain);

    int saveReactionForces = 1;
    data.setValue(saveReactionForces);
    outputSettings.push_back(saveReactionForces);
    Property property_RF("Reaction forces",data,Property::PropertyGroup_OutputSettings);
    props.push_back(property_RF);

    int saveContactData = 1;
    data.setValue(saveContactData);
    outputSettings.push_back(saveContactData);
    Property property_contactData("Contact data",data,Property::PropertyGroup_OutputSettings);
    props.push_back(property_contactData);

    data.setValue(outputSettings);
    QVector<QVariant> values9; values9.push_back(data);
    load loadOuputSettings(values9,Property::loadType_outputSettings);

    //! -------------------------------------
    //! options for storing analysis results
    //! "0" => "All time points" (defaults)
    //! "1" => "Last time point"
    //! "2" => "Specified recurrence rate"
    //! -------------------------------------
    int storeResultsAt = 0;
    data.setValue(storeResultsAt);
    Property property_storeResultsAt("Store results at",data,Property::PropertyGroup_OutputSettings);
    props.push_back(property_storeResultsAt);

    int FREQUENCY = 1;
    QVector<int> storeResultsAtFlags;
    storeResultsAtFlags.push_back(storeResultsAt);
    storeResultsAtFlags.push_back(FREQUENCY);
    data.setValue(storeResultsAtFlags);
    QVector<QVariant> values10; values10.push_back(data);
    load loadStoreResultsAt(values10,Property::loadType_storeResultsAt);

    //! ------------------------------------
    //! create the "Analysis settings" node
    //! ------------------------------------
    SimulationNodeClass *nodeThermalAnalysisSettings =
            new SimulationNodeClass("Analysis settings",SimulationNodeClass::nodeType_thermalAnalysisSettings,props,this);

    //! -------------------------------------------------
    //! create the tabular data for the time step policy
    //! -------------------------------------------------
    QVector<QVariant> stepNumbers;
    QVector<QVariant> stepEndTimes;
    QVector<QVariant> solverTypes;
    QVector<QVariant> vecSubsteps;

    //! time step number - a default value
    data.setValue(1);
    stepNumbers.push_back(data);

    //! Step end time - a default value
    double aStepEndTime = 1.0;
    data.setValue(aStepEndTime);
    stepEndTimes.push_back(data);

    //! Solver type - a default value
    data.setValue(Property::solverType_programControlled);
    solverTypes.push_back(data);

    //! ----------------------------------------------------------------------
    //! time step divisions
    //! the last value is for the time step policy (Auto = 0, ON = 1, OFF =2)
    //! ----------------------------------------------------------------------
    QVector<int> T;
    T.append(5); T.append(2); T.append(10);
    //! the time step policy
    T.append(0);

    data.setValue(T);
    vecSubsteps.append(data);

    load aLoad1(stepNumbers,Property::loadType_stepNumber);
    load aLoad2(stepEndTimes,Property::loadType_stepEndTime);
    load aLoad3(solverTypes,Property::loadType_solverType);
    load aLoad4(vecSubsteps,Property::loadType_autoTimeStepping);

    //! --------------
    //! analysis type
    //! --------------
    data.setValue(analysisType);
    QVector<QVariant> vecAnalysisType;
    vecAnalysisType.push_back(data);
    load loadAnalysisType(vecAnalysisType,Property::loadType_analysisType);

    //! --------------------------
    //! static/transient in table
    //! --------------------------
    data.setValue(timeIntegration);
    QVector<QVariant> vecTimeIntegration;
    vecTimeIntegration.push_back(data);
    load loadTimeIntegration(vecTimeIntegration,Property::loadType_timeIntegration);

    //! -------------------------
    //! build the internal table
    //! -------------------------
    QVector<load> vecLoad;

    //! 1-st time step numner (column # 0)
    vecLoad.push_back(aLoad1);
    //! 2-nd step end time (column # 1)
    vecLoad.push_back(aLoad2);
    //! 3-rd column solver type (column # 2)
    vecLoad.push_back(aLoad3);
    //! 4-th column time step policy (column # 3)
    vecLoad.push_back(aLoad4);
    //! 5-th column (column # 4)
    vecLoad.push_back(loadFieldParameters);
    //! 6-th column (column # 5)
    vecLoad.push_back(loadFluxConvergence);
    //! 7-th column (column # 6)
    vecLoad.push_back(loadSolutionConvergence);
    //! 8-th column (column # 7)
    vecLoad.push_back(loadTimeIncrementation);
    //! 9-th column (column # 8)
    vecLoad.push_back(loadTimeIncrementationType);
    //! 10-th column (column # 9)
    vecLoad.push_back(loadCutBackType);
    //! 11-th column (column #10)
    vecLoad.push_back(loadCutBackParameters);
    //! 12-th column (column #11)
    vecLoad.push_back(loadLineSearch);
    //! 13-th column (column #12)
    vecLoad.push_back(loadLineSearchParameters);
    //! 14-th column (column #13)
    vecLoad.push_back(loadOuputSettings);
    //! 15-th column (column #14)
    vecLoad.push_back(loadStoreResultsAt);
    //! 16-th column (column #15)
    vecLoad.push_back(loadAnalysisType);
    //! 17-th column (column #16)
    vecLoad.push_back(loadTimeIntegration);

    //! ------------------------
    //! create the tabular data
    //! ------------------------
    nodeThermalAnalysisSettings->createTabularData(vecLoad);

    //! -----------------------------
    //! time tag and parent time tag
    //! -----------------------------
    nodeThermalAnalysisSettings->addTimeTag();
    QString parentTimeTag = nodeThermalAnalysis->getPropertyValue<QString>("Time tag");
    data.setValue(parentTimeTag);
    nodeThermalAnalysisSettings->addProperty(Property("Parent time tag",data,Property::PropertyGroup_Identifier));

    //! -----------------------------------------------------
    //! create the item and append to the "Thermal analysis"
    //! -----------------------------------------------------
    data.setValue(nodeThermalAnalysisSettings);
    ThermalAnalysisSettingsItem->setData(data,Qt::UserRole);
    ThermalAnalysisSettingsItem->setData("Analysis settings", Qt::DisplayRole);
    ThermalAnalysisItem->appendRow(ThermalAnalysisSettingsItem);

    //! ---------------------------
    //! create the "Solution" item
    //! ---------------------------
    props.clear();
    Property::solutionInformation theSolutionInformation = Property::solutionInformation_solveRequired;
    data.setValue(theSolutionInformation);
    Property prop_status("Status",data,Property::PropertyGroup_Information);
    props.push_back(prop_status);

    //! --------------------
    //! "Project files dir"
    //! --------------------
    data.setValue(QString("undefined"));
    Property prop_solutionFileDir("Project files dir",data,Property::PropertyGroup_Information);
    props.push_back(prop_solutionFileDir);

    SimulationNodeClass *nodeSolution = new SimulationNodeClass("Solution",SimulationNodeClass::nodeType_thermalAnalysisSolution,props,this);
    data.setValue(nodeSolution);

    //! -----------------------------
    //! time tag and parent time tag
    //! -----------------------------
    nodeSolution->addTimeTag();
    parentTimeTag = ThermalAnalysisItem->data(Qt::UserRole).value<SimulationNodeClass*>()->getPropertyValue<QString>("Time tag");
    data.setValue(parentTimeTag);
    nodeSolution->addProperty(Property("Parent time tag",data,Property::PropertyGroup_Identifier));

    //! ----------------
    //! create the item
    //! ----------------
    QExtendedStandardItem *SolutionItem = new QExtendedStandardItem();
    data.setValue(nodeSolution);
    SolutionItem->setData(data,Qt::UserRole);
    SolutionItem->setData("Solution",Qt::DisplayRole);

    ThermalAnalysisItem->appendRow(SolutionItem);

    //! ---------------------------------------
    //! create the item "Solution information"
    //! ---------------------------------------
    props.clear();

    //! ------------------------------
    //! 0 => solver txt messages
    //! 1 => force convergence
    //! 2 => displacement convergence
    //! 3 => line search
    //! 4 => time step size
    //! ------------------------------
    data.setValue(int(0));
    Property prop_solutionInfoType("Solution information",data,Property::PropertyGroup_SolutionInfo);
    double updateInterval = 2.5;
    data.setValue(updateInterval);
    Property prop_updateInterval("Update interval",data,Property::PropertyGroup_SolutionInfo);

    props.push_back(prop_solutionInfoType);
    props.push_back(prop_updateInterval);

    //! ---------------------------
    //! hidden: solver text output
    //! ---------------------------
    CCXSolverMessage msg;
    data.setValue(msg);
    Property prop_solverOutput("Solver output",data,Property::PropertyGroup_Hidden);
    props.push_back(prop_solverOutput);

    //! ------------------
    //! discrete time map
    //! ------------------
    QMap<double,QVector<int>> dtm;
    data.setValue(dtm);
    Property prop_discreteTimeMap("Discrete time map",data,Property::PropertyGroup_Hidden);
    props.push_back(prop_discreteTimeMap);

    //! -----------------
    //! convergence data
    //! -----------------
    QList<solutionInfo> solInfoList;
    data.setValue(solInfoList);
    Property prop_solInfoList("Convergence data",data,Property::PropertyGroup_Hidden);
    props.push_back(prop_solInfoList);

    SimulationNodeClass *solInfoNode =
            new SimulationNodeClass("Solution information",SimulationNodeClass::nodeType_thermalAnalysisSolutionInformation,props,this);

    //! -----------------------------------------------------------
    //! time tag and parent time tag - the time tag property
    //! has been already created by the node factory for this node
    //! -----------------------------------------------------------
    solInfoNode->addTimeTag();
    parentTimeTag = nodeSolution->getPropertyValue<QString>("Time tag");
    data.setValue(parentTimeTag);
    solInfoNode->addProperty(Property("Parent time tag",data,Property::PropertyGroup_Identifier));

    //! -------------------------------------
    //! create the solution information item
    //! -------------------------------------
    QExtendedStandardItem *itemSolutionInformation = new QExtendedStandardItem;

    data.setValue(solInfoNode);
    itemSolutionInformation->setData(data,Qt::UserRole);
    data.setValue(solInfoNode->getName());
    itemSolutionInformation->setData(data,Qt::DisplayRole);
    SolutionItem->appendRow(itemSolutionInformation);

    cout<<"simulationDataBase::createThermalAnalysisRootNode()->____creating a thermal analysis branch: done____"<<endl;
}


//! ------------------------------------------
//! function: createParticlesInFieldsRootNode
//! details:
//! ------------------------------------------
void simulationDataBase::createParticlesInFieldsRootNode()
{
    QVariant data;
    QString name = "Particles in fields";

    //! ---------------
    //! the properties
    //! ---------------
    QVector<Property> props;

    //! --------------
    //! node creation
    //! --------------
    SimulationNodeClass *nodeParticlesInFieldsAnalysis = new SimulationNodeClass(name,SimulationNodeClass::nodeType_particlesInFieldsAnalysis,props,this);

    //! -----------------------------
    //! time tag and parent time tag
    //! -----------------------------
    nodeParticlesInFieldsAnalysis->addTimeTag();
    QString timeTag = myRootItem->data(Qt::UserRole).value<SimulationNodeClass*>()->getPropertyValue<QString>("Time tag");
    data.setValue(timeTag);
    nodeParticlesInFieldsAnalysis->addProperty(Property("Parent time tag",data,Property::PropertyGroup_Identifier));

    //! --------------------------------------------------
    //! create the item and append to the model root item
    //! --------------------------------------------------
    QExtendedStandardItem *ParticlesInFieldsAnalysisItem = new QExtendedStandardItem();
    ParticlesInFieldsAnalysisItem->setData("Particles in fields",Qt::DisplayRole);
    data.setValue(nodeParticlesInFieldsAnalysis);
    ParticlesInFieldsAnalysisItem->setData(data, Qt::UserRole);
    myRootItem->appendRow(ParticlesInFieldsAnalysisItem);

    //! -----------------------------------------------
    //! create the item structural "Analysis settings"
    //! -----------------------------------------------
    QExtendedStandardItem *ParticlesInFieldsAnalysisSettingsItem = new QExtendedStandardItem();

    //! ------------------
    //! "Number of steps"
    //! ------------------
    data.setValue(1);
    Property property_numberOfSteps("Number of steps",data,Property::PropertyGroup_StepControls);
    props.push_back(property_numberOfSteps);

    //! ----------------------
    //! "Current step number"
    //! ----------------------
    data.setValue(1);
    Property property_currentStepNumber("Current step number",data,Property::PropertyGroup_StepControls);
    props.push_back(property_currentStepNumber);

    //! ------------------------------
    //! the default analysis end time
    //! ------------------------------
    double endTime =1.0;
    data.setValue(endTime);
    Property property_stepEndTime("Step end time",data,Property::PropertyGroup_StepControls);
    props.push_back(property_stepEndTime);

    //! ----------
    //! time step
    //! ----------
    double timeStepSize = endTime/1000;
    data.setValue(timeStepSize);
    Property property_autoTimeStepSize("Time step size",data,Property::PropertyGroup_StepControls);
    props.push_back(property_autoTimeStepSize);

    //! ------------------------------------
    //! create the "Analysis settings" node
    //! ------------------------------------
    SimulationNodeClass *nodeParticlesInFieldsAnalysisSettings =
            new SimulationNodeClass("Analysis settings",SimulationNodeClass::nodeType_particlesInFieldsAnalysisSettings,props,this);

    //! -----------------------------
    //! time tag and parent time tag
    //! -----------------------------
    nodeParticlesInFieldsAnalysisSettings->addTimeTag();
    QString parentTimeTag = nodeParticlesInFieldsAnalysis->getPropertyValue<QString>("Time tag");
    data.setValue(parentTimeTag);
    nodeParticlesInFieldsAnalysisSettings->addProperty(Property("Parent time tag",data,Property::PropertyGroup_Identifier));

    //! ------------------------------------
    //! create the "Analisys settings" item
    //! ------------------------------------
    data.setValue(nodeParticlesInFieldsAnalysisSettings);
    ParticlesInFieldsAnalysisSettingsItem->setData(data,Qt::UserRole);
    data.setValue(nodeParticlesInFieldsAnalysisSettings->getName());
    ParticlesInFieldsAnalysisSettingsItem->setData(data,Qt::DisplayRole);
    ParticlesInFieldsAnalysisItem->appendRow(ParticlesInFieldsAnalysisSettingsItem);

    //! ---------------------------
    //! create the "Solution" item
    //! ---------------------------
    //props.clear();
    //Property::solutionInformation theSolutionInformation = Property::solutionInformation_solveRequired;
    //data.setValue(theSolutionInformation);
    //Property prop_status("Status",data,Property::PropertyGroup_Information);
    //props.push_back(prop_status);

    //! --------------------
    //! "Project files dir"
    //! --------------------
    data.setValue(QString("undefined"));
    Property prop_solutionFileDir("Project files dir",data,Property::PropertyGroup_Information);
    props.push_back(prop_solutionFileDir);

    SimulationNodeClass *nodeSolution = new SimulationNodeClass("Solution",SimulationNodeClass::nodeType_particlesInFieldsSolution,props,this);
    data.setValue(nodeSolution);

    //! -----------------------------
    //! time tag and parent time tag
    //! -----------------------------
    nodeSolution->addTimeTag();
    parentTimeTag = ParticlesInFieldsAnalysisItem->data(Qt::UserRole).value<SimulationNodeClass*>()->getPropertyValue<QString>("Time tag");
    data.setValue(parentTimeTag);
    nodeSolution->addProperty(Property("Parent time tag",data,Property::PropertyGroup_Identifier));

    //! ----------------
    //! create the item
    //! ----------------
    QExtendedStandardItem *SolutionItem = new QExtendedStandardItem();
    data.setValue(nodeSolution);
    SolutionItem->setData(data,Qt::UserRole);
    SolutionItem->setData("Solution",Qt::DisplayRole);

    ParticlesInFieldsAnalysisItem->appendRow(SolutionItem);
}
