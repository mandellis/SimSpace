//! ----------------
//! custom includes
//! ----------------
#include "deserializerclass.h"
#include "qextendedstandarditem.h"
#include "src/utils/tools.h"
#include "load.h"
#include "src/main/mydefines.h"

//! ----
//! C++
//! ----
#include <iostream>

//! ---
//! Qt
//! ---
#include <QDir>

using namespace std;

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
deserializerClass::deserializerClass(QObject *parent): QObject(parent)
{
    ;
}

//! --------------------------------------------------
//! function: deserialize
//! details:  [A] read the name (QString)
//!           [B] read the type (enum)
//!           [C] read the number of properties (int)
//!           [D] read the properties
//! --------------------------------------------------
QVector<Property> deserializerClass::deserialize(const QString &nodeFileName,
                                                 QString &nodeType,
                                                 QString &nodeName) const
{
    QVector<Property> vecProp;

    //! read a node file
    ifstream nodeFile;
    nodeFile.open(nodeFileName.toStdString());

    //! [A] read the node name
    nodeName = tools::readQVariant(nodeFile).toString();

    cout<<"\n\n------------------------------------------------------"<<endl;
    cout<<"- READING NODE FILE: "<<nodeName.toStdString()<<endl;

    //! [B] read the node type
    nodeType = tools::readQVariant(nodeFile).toString();
    cout<<"- NODE TYPE: "<<nodeType.toStdString()<<endl;

    //! [C] read the number of properties
    int nprop = tools::readQVariant(nodeFile).toInt();

    cout<<"- NUMBER OF PROPERTIES TO READ: "<<nprop<<endl;
    cout<<"------------------------------------------------------"<<endl;

    //! [D] read the properties
    for(int k=0; k<nprop; k++)
    {
        Property theProp;
        Property::readProperty(nodeFile,theProp);
        vecProp.push_back(theProp);
    }
    nodeFile.close();
    return vecProp;
}

//! -------------------------------
//! function: nodeBuilder
//! details:  build node from file
//! -------------------------------
SimulationNodeClass* deserializerClass::nodeBuilder(const QString &nodeFileName)
{
    QString nodeTypeName, nodeName;
    QVector<Property> vecProp = this->deserialize(nodeFileName, nodeTypeName, nodeName);
    SimulationNodeClass *node = new SimulationNodeClass(nodeName, SimulationNodeClass::nodeType_NULL, vecProp);
    node->setType1(nodeTypeName);
    return node;
}

//! -----------------------------
//! function: getNodeName
//! details:  check if unused...
//! -----------------------------
QString deserializerClass::getNodeName(const QString &nodeFileName)
{
    ifstream nodeFile;
    nodeFile.open(nodeFileName.toStdString());
    QString nodeName = tools::readQVariant(nodeFile).toString();
    nodeFile.close();
    return nodeName;
}

//! ----------------------------------
//! function: readTabularDataFromFile
//! details:
//! ----------------------------------
QVector<load> deserializerClass::readTabularDataFromFile(const QString &tabularDataFileAbsolutePath) const
{
    cout<<"deserializerClass::readTabularDataFromFile()->____function called____"<<endl;
    ifstream tabularDataFile;
    QString tabularDataFileName = tabularDataFileAbsolutePath;
    cout<<"- READING TABULAR DATA FILE: "<<tabularDataFileName.toStdString()<<endl;

    tabularDataFile.open(tabularDataFileName.toStdString());
    int Nloads;
    QVariant data = tools::readQVariant(tabularDataFile);
    Nloads = data.toInt();
    cout<<"- NUMBER OF LOADS: "<<Nloads<<endl;

    QVector<load> vecLoads;
    for(int i=0; i<Nloads; i++)
    {
        cout<<"deserializerClass::readTabularDataFromFile()->____reading load "<<i<<"____"<<endl;

        //! use static function load::readLoad()
        const load &aLoad = load::readLoad(tabularDataFile);

        //load aLoad;
        //load::readLoad(aLoad, tabularDataFile);

        vecLoads.push_back(aLoad);
        cout<<"deserializerClass::readTabularDataFromFile()->____load "<<i<<" read____"<<endl;
    }
    tabularDataFile.close();
    return vecLoads;
}

//! --------------------------
//! function: getNodeType
//! details:  check if unused
//! --------------------------
QString deserializerClass::getNodeType(const QString &nodeFileName)
{
    ifstream nodeFile;
    nodeFile.open(nodeFileName.toStdString());
    tools::readQVariant(nodeFile).toString();
    QString nodeTypeName = tools::readQVariant(nodeFile).toString();
    nodeFile.close();
    return nodeTypeName;
}
