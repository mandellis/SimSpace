//! ----------------
//! custom includes
//! ----------------
#include "property.h"
#include "tools.h"
#include "qextendedstandarditem.h"
#include "postobject.h"
#include "deserializerclass.h"
#include <ng_meshvs_datasource2d.h>
#include <ng_meshvs_datasource3d.h>
#include "topods_shape_reg.h"
#include <indexedmapofmeshdatasources.h>
#include <qhistogram.h>
#include "ccxsolvermessage.h"
#include "ccout.h"
#include "solutioninfo.h"

//! ----
//! OCC
//! ----
#include <BinTools.hxx>
#include <BRep_Builder.hxx>
#include <BRepTools.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Solid.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_CompSolid.hxx>

//! ----------------------------------------------------------
//! function: constructor
//! details:  create an empty property. Set name and set data
//!           should be called
//! ----------------------------------------------------------
Property::Property()
{
    propertyName = "default";
    propertyValue = QVariant();
    propertyGroup = PropertyGroup_NULL;
    this->init(QVariant());
}

//! -----------------------------------------------
//! function: constructor I
//! details:  create a property with name and data
//! -----------------------------------------------
Property::Property(QString name, QVariant data, PropertyGroup group):
    propertyName(name),
    propertyValue(data),
    propertyGroup(group)
{
    this->init(data);
}

//! ---------------------------
//! function: copy constructor
//! details:
//! ---------------------------
Property::Property(const Property &other)
{
    this->setName(other.propertyName);
    this->setData(other.propertyValue);
    this->setGroup(other.propertyGroup);
    //this->init(other.propertyValue);
}

//! ---------------
//! function: init
//! details:
//! ---------------
void Property::init(const QVariant &data)
{
    if(data.canConvert<elementControl>()) ec = data.value<elementControl>();
    if(data.canConvert<meshEngine2D>()) me2d = data.value<meshEngine2D>();
    if(data.canConvert<meshEngine3D>()) me3d = data.value<meshEngine3D>();
    if(data.canConvert<meshOrder>()) mor = data.value<meshOrder>();
    if(data.canConvert<meshType_Surface>()) mts = data.value<meshType_Surface>();
    if(data.canConvert<meshType_Volume>()) mtv = data.value<meshType_Volume>();
    if(data.canConvert<meshMethod>()) mm = data.value<meshMethod>();
    if(data.canConvert<ScopingMethod>()) sm = data.value<ScopingMethod>();
    if(data.canConvert<SuppressionStatus>()) ss = data.value<SuppressionStatus>();
    if(data.canConvert<autoTimeStepping>()) ats = data.value<autoTimeStepping>();
    if(data.canConvert<solverType>()) st = data.value<solverType>();
    if(data.canConvert<solutionInformation>()) sinf = data.value<solutionInformation>();
    if(data.canConvert<integrationScheme>()) is = data.value<integrationScheme>();
    if(data.canConvert<DOFfreedom>()) df = data.value<DOFfreedom>();
    if(data.canConvert<typeOfValue>()) tv = data.value<typeOfValue>();
    if(data.canConvert<defineBy>()) db = data.value<defineBy>();
    if(data.canConvert<solverEngine>()) se = data.value<solverEngine>();
    if(data.canConvert<contactType>()) ct = data.value<contactType>();
    if(data.canConvert<contactBehavior>()) cb = data.value<contactBehavior>();
    if(data.canConvert<contactFormulation>()) cf = data.value<contactFormulation>();
    if(data.canConvert<overpressureFunction>()) of = data.value<overpressureFunction>();
    if(data.canConvert<loadType>()) lt = data.value<loadType>();
    if(data.canConvert<loadDefinition>()) ld = data.value<loadDefinition>();
    if(data.canConvert<typeOfTransformation>()) tt = data.value<typeOfTransformation>();
    if(data.canConvert<modelChangeActivationStatus>()) mcas = data.value<modelChangeActivationStatus>();
    if(data.canConvert<analysisType>()) ant = data.value<analysisType>();
    if(data.canConvert<timeIntegration>()) ti = data.value<timeIntegration>();
}

//! -------------------
//! function: set data
//! details:
//! -------------------
void Property::setData(const QVariant &data)
{
    propertyValue = data;
    this->init(data);
}

//! -------------------------------------------------------
//! function: getPropertyGroup
//! details:  convert the PropertyGroup enum into a string
//! -------------------------------------------------------
const QString Property::getPropertyGroup() const
{
    const QMetaObject &mo = Property::staticMetaObject;
    int index = mo.indexOfEnumerator("PropertyGroup");
    QMetaEnum metaEnum = mo.enumerator(index);
    return metaEnum.valueToKey(propertyGroup);
}

//! -------------------------------
//! function: getPropertyGroupEnum
//! details:
//! -------------------------------
const Property::PropertyGroup Property::getPropertyGroupEnum() const
{
    return this->propertyGroup;
}

//! -------------------------------------------------------
//! function: getPropertyValue
//! details:  convert the PropertyValue enum into a string
//! -------------------------------------------------------
const QString Property::getPropertyValue(const QString &propName) const
{
    const QMetaObject &mo = Property::staticMetaObject;
    QString s = propName;
    s.remove(0,10);
    int index = mo.indexOfEnumerator(qPrintable(s));

    QMetaEnum metaEnum = mo.enumerator(index);
    QString RV;

    if(QString(metaEnum.name())=="elementControl") RV = metaEnum.valueToKey(ec);
    if(QString(metaEnum.name())=="meshEngine2D") RV =  metaEnum.valueToKey(me2d);
    if(QString(metaEnum.name())=="meshEngine3D") RV =  metaEnum.valueToKey(me3d);
    if(QString(metaEnum.name())=="meshType_Surface") RV =  metaEnum.valueToKey(mts);
    if(QString(metaEnum.name())=="meshType_Volume") RV =  metaEnum.valueToKey(mtv);
    if(QString(metaEnum.name())=="meshOrder") RV =  metaEnum.valueToKey(mor);
    if(QString(metaEnum.name())=="meshMethod") RV = metaEnum.valueToKey(mm);
    if(QString(metaEnum.name())=="ScopingMethod") RV =  metaEnum.valueToKey(sm);
    if(QString(metaEnum.name())=="SuppressionStatus") RV =  metaEnum.valueToKey(ss);
    if(QString(metaEnum.name())=="autoTimeStepping") RV =  metaEnum.valueToKey(ats);
    if(QString(metaEnum.name())=="solverType") RV =  metaEnum.valueToKey(st);
    if(QString(metaEnum.name())=="solutionInformation") RV =  metaEnum.valueToKey(sinf);
    if(QString(metaEnum.name())=="integrationScheme") RV =  metaEnum.valueToKey(is);
    if(QString(metaEnum.name())=="DOFfreedom") RV =  metaEnum.valueToKey(df);
    if(QString(metaEnum.name())=="typeOfValue") RV =  metaEnum.valueToKey(tv);
    if(QString(metaEnum.name())=="defineBy") RV =  metaEnum.valueToKey(db);
    if(QString(metaEnum.name())=="solverEngine") RV =  metaEnum.valueToKey(se);
    if(QString(metaEnum.name())=="contactType") RV =  metaEnum.valueToKey(ct);
    if(QString(metaEnum.name())=="contactBehavior") RV =  metaEnum.valueToKey(cb);
    if(QString(metaEnum.name())=="contactFormulation") RV =  metaEnum.valueToKey(cf);
    if(QString(metaEnum.name())=="overpressureFunction") RV =  metaEnum.valueToKey(of);
    if(QString(metaEnum.name())=="typeOfTransformation") RV = metaEnum.valueToKey(tt);
    if(QString(metaEnum.name())=="loadType") RV = metaEnum.valueToKey(lt);
    if(QString(metaEnum.name())=="loadDefinition") RV = metaEnum.valueToKey(ld);
    if(QString(metaEnum.name())=="boltStatusDefinedBy") RV = metaEnum.valueToKey(bs);
    if(QString(metaEnum.name())=="modelChangeActivationStatus") RV = metaEnum.valueToKey(mcas);
    if(QString(metaEnum.name())=="analysisType") RV = metaEnum.valueToKey(ant);
    if(QString(metaEnum.name())=="timeIntegration") RV = metaEnum.valueToKey(ti);

    return RV;
}

//! ----------------------------------------------
//! function: setPropertyGroup
//! details:  set the member value using a string
//! ----------------------------------------------
void Property::setPropertyGroup(const QString &thePropertyGroup)
{
    const QMetaObject &mo = Property::staticMetaObject;
    int index = mo.indexOfEnumerator("PropertyGroup");
    QMetaEnum metaEnum = mo.enumerator(index);
    int value = metaEnum.keyToValue(thePropertyGroup.toLocal8Bit().constData());
    propertyGroup = static_cast<PropertyGroup>(value);
}

//! --------------------------------------------------------------
//! function: setPropertyValue
//! details:  set the enum value of a property using string input
//! --------------------------------------------------------------
void Property::setPropertyValue(const QString &enumName, const QString &thePropertyValue)
{
    const QMetaObject &mo = Property::staticMetaObject;
    QString enumName1 = enumName;
    enumName1.remove(0,10);
    int index = mo.indexOfEnumerator(qPrintable(enumName1));

    QMetaEnum metaEnum = mo.enumerator(index);
    int value = metaEnum.keyToValue(qPrintable(thePropertyValue));

    QVariant data;

    if(QString(metaEnum.name())=="boltStatusDefinedBy")
    {
        bs = static_cast<boltStatusDefinedBy>(value);
        data.setValue(bs);
    }
    if(QString(metaEnum.name())=="meshEngine2D")
    {
        me2d = static_cast<meshEngine2D>(value);
        data.setValue(me2d);
    }
    if(QString(metaEnum.name())=="meshEngine3D")
    {
        me3d = static_cast<meshEngine3D>(value);
        data.setValue(me3d);
    }
    if(QString(metaEnum.name())=="meshType_Surface")
    {
        mts = static_cast<meshType_Surface>(value);
        data.setValue(mts);
    }
    if(QString(metaEnum.name())=="meshType_Volume")
    {
        mtv = static_cast<meshType_Volume>(value);
        data.setValue(mtv);
    }
    if(QString(metaEnum.name())=="meshOrder")
    {
        mor = static_cast<meshOrder>(value);
        data.setValue(mor);
    }
    if(QString(metaEnum.name())=="meshMethod")
    {
        mm = static_cast<meshMethod>(value);
        data.setValue(mm);
    }
    if(QString(metaEnum.name())=="ScopingMethod")
    {
        sm = static_cast<ScopingMethod>(value);
        data.setValue(sm);
    }
    if(QString(metaEnum.name())=="SuppressionStatus")
    {
        ss = static_cast<SuppressionStatus>(value);
        data.setValue(ss);
    }
    if(QString(metaEnum.name())=="autoTimeStepping")
    {
        ats = static_cast<autoTimeStepping>(value);
        data.setValue(ats);
    }
    if(QString(metaEnum.name())=="integrationScheme")
    {
        is = static_cast<integrationScheme>(value);
        data.setValue(is);
    }
    if(QString(metaEnum.name())=="solverType")
    {
        st = static_cast<solverType>(value);
        data.setValue(st);
    }
    if(QString(metaEnum.name())=="elementControl")
    {
        ec = static_cast<elementControl>(value);
        data.setValue(ec);
    }
    if(QString(metaEnum.name())=="DOFfreedom")
    {
        df = static_cast<DOFfreedom>(value);
        data.setValue(df);
    }
    if(QString(metaEnum.name())=="typeOfValue")
    {
        tv = static_cast<typeOfValue>(value);
        data.setValue(tv);
    }
    if(QString(metaEnum.name())=="defineBy")
    {
        db = static_cast<defineBy>(value);
        data.setValue(db);
    }
    if(QString(metaEnum.name())=="solverEngine")
    {
        se = static_cast<solverEngine>(value);
        data.setValue(se);
    }
    if(QString(metaEnum.name())=="contactType")
    {
        ct = static_cast<contactType>(value);
        data.setValue(ct);
    }
    if(QString(metaEnum.name())=="contactBehavior")
    {
        cb = static_cast<contactBehavior>(value);
        data.setValue(cb);
    }
    if(QString(metaEnum.name())=="contactFormulation")
    {
        cf = static_cast<contactFormulation>(value);
        data.setValue(cf);
    }
    if(QString(metaEnum.name())=="overpressureFunction")
    {
        of = static_cast<overpressureFunction>(value);
        data.setValue(of);
    }
    if(QString(metaEnum.name())=="typeOfTransformation")
    {
        tt = static_cast<typeOfTransformation>(value);
        data.setValue(tt);
    }
    if(QString(metaEnum.name())=="solutionInformation")
    {
        sinf = static_cast<solutionInformation>(value);
        data.setValue(sinf);
    }
    if(QString(metaEnum.name())=="analysisType")
    {
        ant = static_cast<analysisType>(value);
        data.setValue(ant);
    }
    this->setData(data);
}

//! ----------------------------------------------------------
//! function: writeProperty
//! details:  write a "Property" on a file
//!           [1] write the key name
//!           [2] write the value. When writing a scope write
//!               first the number of shapes
//!           [3] write the PropertyGroup enum
//! -----------------------------------------------------------
void Property::writeProperty(ofstream& out, const Property &prop)
{
    cout<<"Property::writeProperty->____writing property: "<<prop.getName().toStdString()<<"____"<<endl;
    const QMap<QString,QString> &thePropertyMap = Property::propertyMap();

    //! -------------------
    //! write the key name
    //! -------------------
    QString keyName = prop.getName();
    tools::writeQVariant(QVariant(keyName),out);

    //! -----------------------------------------------------
    //! write the PropertyGroup - PropertyGroup is an "enum"
    //! -----------------------------------------------------
    QString propGroup = prop.getPropertyGroup();
    tools::writeQVariant(QVariant(propGroup),out);

    //! ----------------
    //! write the value
    //! ----------------
    QString enumName = thePropertyMap.value(keyName);

    if(!enumName.isEmpty())
    {
        cout<<"* PROPERTY DEFINED THROUGH \"ENUMERATOR\""<<endl;
        //! -------------------------------------------------------------
        //! the property with "keyName" is defined through an enumerator
        //! -------------------------------------------------------------
        QVariant data;

        //! ---------------------
        //! write the enum value
        //! ---------------------
        QString theEnumValue = prop.getPropertyValue(enumName);

        cout<<"* "<<theEnumValue.toStdString()<<endl;

        data.setValue(theEnumValue);
        tools::writeQVariant(data,out);
    }
    else
    {
        //! ------------------
        //! write other types
        //! ------------------
        QVariant data = prop.getData();
        if(data.canConvert<TopoDS_Shape_Reg>())
        {
            const TopoDS_Shape &theCurShape = data.value<TopoDS_Shape_Reg>();
            //! ---------------------------------------
            //! clean the triangulation and save space
            //! ---------------------------------------
            BRepTools::Clean(theCurShape);
            BRepTools::Write(theCurShape,out);
            out<<endl;
        }
        else if(prop.getData().canConvert<QList<solutionInfo>>())
        {
            const QList<solutionInfo> solInfoList = prop.getData().value<QList<solutionInfo>>();
            int Nb = solInfoList.length();
            out<<Nb<<endl;
            for(int i=0; i<Nb; i++)
            {
                const solutionInfo &aSolInfo = solInfoList[i];
                aSolInfo.write(out);
            }
        }
        else if(prop.getData().canConvert<CCXSolverMessage>())
        {
            prop.getData().value<CCXSolverMessage>().write(out);
        }
        else if(prop.getData().canConvert<QMap<double,QVector<int>>>())
        {
            QMap<double,QVector<int>> amap = prop.getData().value<QMap<double,QVector<int>>>();
            int NbKeys = amap.keys().length();
            out<<NbKeys<<endl;
            for(QMap<double,QVector<int>>::iterator it = amap.begin(); it!=amap.end(); it++)
            {
                double key = it.key();
                out<<key<<endl;

                int NbValues = it.value().length();
                out<<NbValues<<endl;

                for(int n=0; n<NbValues; n++)
                {
                    out<<it.value().at(n)<<endl;
                }
            }
        }
        else if(prop.getData().canConvert<QVector<int>>() ||
                prop.getData().canConvert<QVector<float>>() ||
                prop.getData().canConvert<QVector<bool>>() ||
                prop.getData().canConvert<QVector<double>>())
        {
            if(prop.getData().canConvert<QVector<int>>())
            {
                QVector<int> vec = prop.getData().value<QVector<int>>();
                tools::writeQVector<QVector<int>>(vec,out);
            }
            if(prop.getData().canConvert<QVector<bool>>())
            {
                QVector<bool> vec = prop.getData().value<QVector<bool>>();
                tools::writeQVector<QVector<bool>>(vec,out);
            }
            if(prop.getData().canConvert<QVector<float>>())
            {
                QVector<float> vec = prop.getData().value<QVector<float>>();
                tools::writeQVector<QVector<float>>(vec,out);
            }
            if(prop.getData().canConvert<QVector<double>>())
            {
                QVector<double> vec = prop.getData().value<QVector<double>>();
                tools::writeQVector<QVector<double>>(vec,out);
            }
        }
        else if(prop.getData().canConvert<QVector<QVector<double>>>())
        {
            cout<<"* PROPERTY DEFINED THROUGH \"DOUBLE VECTOR\""<<endl;
            QVector<QVector<double>> tensor2 = prop.getData().value<QVector<QVector<double>>>();
            tools::writeTensor2<QVector<QVector<double>>>(tensor2,out);
        }
        else if(prop.getData().canConvert<QVector<GeometryTag>>())
        {
            cout<<"* PROPERTY DEFINED THROUGH \"QVector<GeometryTag>\""<<endl;
            QVector<GeometryTag> vecLocs = prop.getData().value<QVector<GeometryTag>>();
            //out<<vecLocs.size()<<endl;
            //for(QVector<GeometryTag>::iterator it = vecLocs.begin(); it!=vecLocs.end(); ++it)
            //{
            //    const GeometryTag &aGT = *it;
            //    aGT.write(out);
            //}
            tools::writeVectorOfLocations(vecLocs,out);
        }
        else if(prop.getData().canConvert<void*>())
        {
            //! -------------------------------------------------------
            //! do nothing: void* is written at serializer class level
            //! -------------------------------------------------------
            cout<<"* PROPERTY DEFINED THROUGH A VOID POINTER"<<endl;
            void *p = prop.getData().value<void*>();
            QExtendedStandardItem *item = static_cast<QExtendedStandardItem*>(p);
            SimulationNodeClass *node = item->data(Qt::UserRole).value<SimulationNodeClass*>();
            cout<<"* THE VOID POINTER CONTAINS: "<<node->getName().toStdString()<<endl;
            Property::writeVoid(out,p);
        }
        else if(prop.getData().canConvert<postObject>())
        {
            cout<<"* PROPERTY POST OBJECT"<<endl;
            postObject aPostObject = prop.getData().value<postObject>();
            cout<<"____writing post object: \""<<aPostObject.getName().toStdString()<<"\"____"<<endl;
            aPostObject.write(out);
        }
        else if(prop.getData().canConvert<histogramData>())
        {
            prop.getData().value<histogramData>().writeHistogramData(out);
        }
        //! ---------------------
        //! experimental section
        //! ---------------------
        else if(prop.getData().canConvert<occHandle(Ng_MeshVS_DataSource2D)>())
        {
            cout<<"* PROPERTY SURFACE MESH"<<endl;
            //const occHandle(Ng_MeshVS_DataSource2D) &theMeshDS = prop.getData().value<occHandle(Ng_MeshVS_DataSource2D)>();
        }
        else if(prop.getData().canConvert<occHandle(Ng_MeshVS_DataSource3D)>())
        {
            cout<<"* PROPERTY VOLUME MESH"<<endl;
            //const occHandle(Ng_MeshVS_DataSource3D) &theMeshDS = prop.getData().value<occHandle(Ng_MeshVS_DataSource3D)>();
        }
        //! -------------------------
        //! end experimental section
        //! -------------------------
        else if(prop.getData().canConvert<IndexedMapOfMeshDataSources>())
        {
            cout<<"* PROPERTY INDEXED MAP OF MESH DATASOURCES"<<endl;
            //! ---------------------------------
            //! write the number of data sources
            //! ---------------------------------
            const IndexedMapOfMeshDataSources &im = prop.getData().value<IndexedMapOfMeshDataSources>();
            int NbMeshes = im.keys().length();
            out<<NbMeshes<<endl;
            //! -----------------
            //! write the meshes
            //! -----------------
            for(IndexedMapOfMeshDataSources::const_iterator it= im.cbegin(); it!= im.cend(); it++)
            {
                //! --------------
                //! write the key
                //! --------------
                int key = it.key();
                out<<key<<endl;
                //! ---------------
                //! write the mesh
                //! ---------------
                const occHandle(MeshVS_DataSource) &aDS = it.value();
                const occHandle(Ng_MeshVS_DataSourceFace) &faceDS = occHandle(Ng_MeshVS_DataSourceFace)::DownCast(aDS);
                int meshDimension = 2;
                faceDS->write(out,meshDimension);
            }
        }
        else
        {
            cout<<"* PROPERTY DEFINED THROUGH STANDARD TYPES"<<endl;
            tools::writeQVariant(data,out);
        }
    }
}

//! ----------------------------------------------------------------
//! function: readProperty
//! details:  read a "Property" from ifstream
//!           [1] read the key name
//!           [2] read the value in native format, not QVariant
//!                 When writing a scope first write the number of
//!                 shapes contained in
//!           [3] read the PropertyGroup enum
//! ----------------------------------------------------------------
void Property::readProperty(ifstream &in, Property &prop)
{
    //! ------------------
    //! read the key name
    //! ------------------
    QVariant data = tools::readQVariant(in);
    QString propKeyName = data.toString();
    prop.setName(data.toString());
    cout<<"\n- keyName: ---->"<<prop.getName().toStdString()<<"<----"<<endl;

    //! ---------------------------
    //! read the ProperyGroup enum
    //! ---------------------------
    const QString &propGroup = tools::readQVariant(in).toString();
    prop.setPropertyGroup(propGroup);
    cout<<"- property group: ---->"<<prop.getPropertyGroup().toStdString()<<"<----"<<endl;

    //! ---------------
    //! read the value
    //! ---------------
    if(propKeyName=="Source geometry")
    {
        BRep_Builder B;
        TopoDS_Shape aShape;
        BRepTools::Read(aShape,in,B);
        QVariant data;
        data.setValue(aShape);
        prop.setData(data);
    }
    else if(propKeyName=="Shape")
    {
        QVariant data;
        BRep_Builder B;
        TopoDS_Shape aShape;
        BRepTools::Read(aShape,in,B);
        data.setValue(aShape);
        prop.setData(data);
    }
    else if(propKeyName=="Geometry" || propKeyName=="Master"
            || propKeyName=="Slave" || propKeyName=="Boundary") //! "Boundary" identify the prismatic faces
    {
        QVector<GeometryTag> vecLocs;
        vecLocs = tools::readVectorOfLocations(in);
        data.setValue(vecLocs);
        prop.setData(data);
    }
    //! -----------------------------------------
    //! property defined through a void* pointer
    //! -----------------------------------------
    else if(propKeyName=="Named selection" ||
            propKeyName=="Coordinate system" ||
            propKeyName=="Remote point" ||
            propKeyName=="Boundary named selection" ||
            propKeyName=="Analysis" ||
            propKeyName=="Imported body temperature")
    {
        //! -----------------------
        //! [A] read the node name
        //! -----------------------
        QString nodeName = tools::readQVariant(in).toString();
        cout<<"------ READING NODE: "<<nodeName.toStdString()<<endl;

        //! -----------------------
        //! [B] read the node type
        //! -----------------------
        QString nodeType = tools::readQVariant(in).toString();
        cout<<"------ NODE TYPE: "<<nodeType.toStdString()<<endl;

        //! ----------------------------------
        //! [C] read the number of properties
        //! ----------------------------------
        int nprop = tools::readQVariant(in).toInt();
        cout<<"------ NUMBER OF PROPERTIES TO READ: "<<nprop<<endl;

        //! ------------------------
        //! [D] read the properties
        //! ------------------------
        QVector<Property> vecProp;
        for(int k=0; k<nprop; k++)
        {
            Property theProp;
            Property::readProperty(in,theProp);
            vecProp.push_back(theProp);
        }

        SimulationNodeClass *node = new SimulationNodeClass(nodeName, SimulationNodeClass::nodeType_NULL, vecProp);

        if(propKeyName=="Named selection") node->setType(SimulationNodeClass::nodeType_namedSelectionGeometry);
        else if(propKeyName=="Coordinate system")
        {
            node->setType(SimulationNodeClass::nodeType_coordinateSystem);
        }
        else if(propKeyName=="Remote point") node->setType(SimulationNodeClass::nodeType_remotePoint);
        else if(propKeyName=="Boundary named selection") node->setType(SimulationNodeClass::nodeType_meshPrismaticLayer);
        else if(propKeyName=="Analysis")
        {
            node->setType(SimulationNodeClass::nodeType_thermalAnalysis);
            if(vecProp.size()==0) node->setType(SimulationNodeClass::nodeType_NULL);
        }
        else if(propKeyName=="Imported body temperature")
        {
            node->setType(SimulationNodeClass::nodeType_solutionThermalTemperature);
            if(vecProp.size()==0) node->setType(SimulationNodeClass::nodeType_NULL);
        }

        QExtendedStandardItem *item = new QExtendedStandardItem();
        QVariant data;
        data.setValue(node);
        item->setData(data,Qt::UserRole);
        item->setData(node->getName(),Qt::DisplayRole);
        void *p = (void*)(item);
        data.setValue(p);
        prop.setData(data);
    }
    //! -------------------------------------------------------------
    //! property defined through an indexed map of mesh data sources
    //! -------------------------------------------------------------
    else if(propKeyName=="Master mesh data source" ||
            propKeyName=="Slave mesh data source" ||
            propKeyName =="Mesh data sources")
    {
        //! -------------------------------------
        //! read an indexedMapOfMeshDataSoources
        //! -------------------------------------
        IndexedMapOfMeshDataSources im;
        int NbKeys;
        in>>NbKeys;
        //cout<<"____number of map keys: "<<NbKeys<<"____"<<endl;
        for(int i=0; i<NbKeys; i++)
        {
            //! -------------
            //! read the key
            //! -------------
            int key;
            in>>key;
            //cout<<"____reading mesh with key: "<<key<<"____"<<endl;
            //! --------------
            //! read the mesh
            //! --------------
            occHandle(Ng_MeshVS_DataSourceFace) aMesh = new Ng_MeshVS_DataSourceFace(in);

            im.insert(key,aMesh);
        }
        data.setValue(im);
        prop.setData(data);
    }
    //! -------------------
    //! "Convergence data"
    //! -------------------
    else if(propKeyName =="Convergence data")
    {
        QList<solutionInfo> solInfoList;
        int Nb;
        in>>Nb;
        for(int i=0; i<Nb; i++)
        {
            solutionInfo aSolInfo;
            aSolInfo.read(in);
            solInfoList<<aSolInfo;
        }
        data.setValue(solInfoList);
        prop.setData(data);
    }
    //! --------------------------------
    //! discrete time map - solver info
    //! --------------------------------
    else if(propKeyName == "Discrete time map")
    {
        QMap<double,QVector<int>> amap;
        int NbKeys;
        in>>NbKeys;
        for(int n=0; n<NbKeys; n++)
        {
            double akey;
            in>>akey;
            cout<<"____key: "<<akey<<"____"<<endl;

            int NbValues;
            in>>NbValues;

            QVector<int> avalue;
            for(int i=0; i<NbValues; i++)
            {
                int val;
                in>>val;
                avalue.push_back(val);
                cout<<"____val "<<i<<": "<<val<<"____"<<endl;
            }
            amap.insert(akey,avalue);
        }
        data.setValue(amap);
        prop.setData(data);
    }
    //! ----------------------------------------------
    //! handle properties defined through enumerators
    //! ----------------------------------------------
    else if(propKeyName == "Mesh engine 2D" ||
            propKeyName == "Mesh engine 3D" ||
            propKeyName == "Surface mesher" ||
            propKeyName == "Volume mesher" ||
            propKeyName == "Mesh order" ||
            propKeyName == "Surface mesh type" ||
            propKeyName == "Surface volume type" ||
            propKeyName == "Scoping method" ||
            propKeyName == "Boundary scoping method" ||
            propKeyName == "Suppressed" ||
            propKeyName == "Auto time stepping" ||
            propKeyName == "Integration scheme" ||
            propKeyName == "Solver type" ||
            propKeyName == "Element control" ||
            propKeyName == "Radial" ||
            propKeyName == "Axial" ||
            propKeyName == "Tangential" ||
            propKeyName == "Define by" ||
            propKeyName == "Define by " ||
            propKeyName == "Type" ||
            propKeyName == "Behavior" ||
            propKeyName == "Formulation" ||
            propKeyName == "Overpressure" ||
            propKeyName == "Magnitude" ||
            propKeyName == "X component" ||
            propKeyName == "Y component" ||
            propKeyName == "Z component" ||
            propKeyName == "Status" ||
            propKeyName == "Activation status" ||
            propKeyName == "Analysis type" ||
            propKeyName == "Time integration")
    {
        QString enumName = Property::propertyMap().value(propKeyName);
        cout<<"- enumName: "<<enumName.toStdString()<<endl;
        QString thePropertyValue = tools::readQVariant(in).toString();
        prop.setPropertyValue(enumName,thePropertyValue);
        cout<<"- keyName: "<<propKeyName.toStdString()<<endl;
        cout<<"- value: "<<prop.getPropertyValue(enumName).toStdString()<<endl;
    }
    else if(propKeyName == "X axis data" || propKeyName == "Y axis data" || propKeyName == "Z axis data" ||
            propKeyName == "Base origin" || propKeyName =="Direction" || propKeyName =="Reference point")
    {
        QVector<double> vec = tools::readQVector<double>(in);
        data.setValue(vec);
        prop.setData(data);
    }
    else if(propKeyName =="First color" || propKeyName =="Second color")
    {
        QVector<int> vec = tools::readQVector<int>(in);
        data.setValue(vec);
        prop.setData(data);
    }
    else if(propKeyName =="Base directional data")
    {
        QVector<QVector<double>> tensor2 = tools::readTensor2<double>(in);
        data.setValue(tensor2);
        prop.setData(data);
    }
    else if(propKeyName=="Tags" || propKeyName =="Tags slave"
            || propKeyName =="Tags master" || propKeyName =="Boundary tags")
    {
        QVector<GeometryTag> vecLocs = tools::readVectorOfLocations(in);
        data.setValue(vecLocs);
        prop.setData(data);
    }
    else if(propKeyName=="Post object")
    {
        //! ----------------------------------------------------
        //! this constructor reads the name, the mesh, the data
        //! but needs to update the MeshVS_Mesh
        //! ----------------------------------------------------
        postObject aPostObject(in);
        data.setValue(aPostObject);
        prop.setData(data);
    }
    else if(propKeyName =="Metric data")
    {
        histogramData hisData;
        hisData.readHistogramData(in);
        data.setValue(hisData);
        prop.setData(data);
    }
    else if(propKeyName =="Solver output")
    {
        CCXSolverMessage msg;
        msg.read(in);
        data.setValue(msg);
        prop.setData(data);
    }
    //! ----------------------------------------
    //! char, int, float, double, bool, QString
    //! ----------------------------------------
    else
    {
        QVariant data = tools::readQVariant(in);
        prop.setData(data);

        //! diagnostic - can be removed
        switch(QMetaType::type(data.typeName()))
        {
        case 1: prop.getData().toBool()==true? cout<<"- value: true"<<endl : cout<<"- value: false"<<endl; break;
        case 2: cout<<"- value: "<<prop.getData().toInt()<<endl; break;
        case 6: cout<<"- value: "<<prop.getData().toDouble()<<endl; break;
        case 10: cout<<"- value: "<<prop.getData().toString().toStdString()<<endl; break;
        case 34: cout<<"- value: "<<prop.getData().toChar().toLatin1()<<endl; break;
        case 38: cout<<"- value: "<<prop.getData().toFloat()<<endl; break;
        }
        //! end diagnostic - can be removed
    }
}

//! ---------------------------------------------------
//! function: createPropertyMap
//! details:  key name => enum name - UNUSED? - un cas
//! ---------------------------------------------------
QMap<QString,QString> Property::propertyMap()
{
    QMap<QString,QString> myPropertyMap;

    myPropertyMap.insert("Activation status","Property::modelChangeActivationStatus");
    myPropertyMap.insert("Surface mesher","Property::meshEngine2D");
    myPropertyMap.insert("Tessellator","Property::meshEngine2D");
    myPropertyMap.insert("Volume mesher","Property::meshEngine3D");
    myPropertyMap.insert("Mesh engine 2D","Property::meshEngine2D");
    myPropertyMap.insert("Mesh engine 3D","Property::meshEngine3D");
    myPropertyMap.insert("Mesh order","Property::meshOrder");
    myPropertyMap.insert("Surface mesh type","Property::meshType_Surface");
    myPropertyMap.insert("Volume mesh type","Property::meshType_Volume");
    myPropertyMap.insert("Scoping method","Property::ScopingMethod");
    myPropertyMap.insert("Boundary scoping method","Property::ScopingMethod");
    myPropertyMap.insert("Suppressed","Property::SuppressionStatus");
    myPropertyMap.insert("Auto time stepping","Property::autoTimeStepping");
    myPropertyMap.insert("Solver type","Property::solverType");
    myPropertyMap.insert("Status","Property::solutionInformation");
    myPropertyMap.insert("Integration scheme","Property::integrationScheme");
    myPropertyMap.insert("Element control","Property::elementControl");
    myPropertyMap.insert("Radial","Property::DOFfreedom");
    myPropertyMap.insert("Axial","Property::DOFfreedom");
    myPropertyMap.insert("Tangential","Property::DOFfreedom");
    myPropertyMap.insert("Define by","Property::defineBy");
    myPropertyMap.insert("Define by ","Property::defineBy");
    myPropertyMap.insert("Type","Property::contactType");
    myPropertyMap.insert("Behavior","Property::contactBehavior");
    myPropertyMap.insert("Formulation","Property::contactFormulation");     //! not used for the moment
    myPropertyMap.insert("Overpressure","Property::overpressureFunction");
    myPropertyMap.insert("Magnitude","Property::loadDefinition");
    myPropertyMap.insert("X component","Property::loadDefinition");
    myPropertyMap.insert("Y component","Property::loadDefinition");
    myPropertyMap.insert("Z component","Property::loadDefinition");
    myPropertyMap.insert("Film coefficient","Property::loadDefinition");
    myPropertyMap.insert("Reference temperature","Property::loadDefinition");
    //myPropertyMap.insert("Offset X","Property::typeOfTransformation");
    //myPropertyMap.insert("Offset Y","Property::typeOfTransformation");
    //myPropertyMap.insert("Offset Z","Property::typeOfTransformation");
    //myPropertyMap.insert("Rotation X","Property::typeOfTransformation");
    //myPropertyMap.insert("Rotation Y","Property::typeOfTransformation");
    //myPropertyMap.insert("Rotation Z","Property::typeOfTransformation");
    myPropertyMap.insert("Analysis type","Property::analysisType");
    myPropertyMap.insert("Static/Transient","Property::timeIntegration");
    return myPropertyMap;
}

//! --------------------
//! function: writeVoid
//! details:
//! --------------------
void Property::writeVoid(ofstream& outFile, void *p)
{
    QExtendedStandardItem *aTreeItem = static_cast<QExtendedStandardItem*>(p);
    SimulationNodeClass *aSimNode = aTreeItem->data(Qt::UserRole).value<SimulationNodeClass*>();
    SimulationNodeClass::nodeType theNode = aSimNode->getType();

    if(theNode==SimulationNodeClass::nodeType_namedSelectionGeometry ||
            theNode==SimulationNodeClass::nodeType_coordinateSystem ||
            theNode==SimulationNodeClass::nodeType_remotePoint ||
            theNode==SimulationNodeClass::nodeType_importedBodyScalar ||
            theNode==SimulationNodeClass::nodeType_thermalAnalysis ||
            theNode==SimulationNodeClass::nodeType_solutionThermalTemperature ||
            theNode==SimulationNodeClass::nodeType_meshPrismaticLayer)
    {
        //! ---------------------------
        //! write the name of the node
        //! ---------------------------
        QVariant data;
        data.setValue(aSimNode->getName());
        tools::writeQVariant(data,outFile);

        //! -------------------
        //! write the nodeType
        //! -------------------
        data.setValue(aSimNode->type());
        tools::writeQVariant(data,outFile);

        //! ------------------------------------------------------------
        //! retrieve the number of properties within the node container
        //! ------------------------------------------------------------
        QVector<QExtendedStandardItem*> vecItem = aSimNode->getPropertyItems();
        int nprop = vecItem.size();

        //! -----------------------------------------------------------------
        //! write the number of properties (drive for the reading operation)
        //! -----------------------------------------------------------------
        data.setValue(nprop);
        tools::writeQVariant(data,outFile);

        for(int k=0; k<vecItem.size(); k++)
        {
            const Property &curProperty = vecItem.at(k)->data(Qt::UserRole).value<Property>();
            Property::writeProperty(outFile,curProperty);
        }
    }
    else
    {
        cout<<"____WRITING A NULL NODE____"<<endl;
        //! ----------
        //! NULL node
        //! ----------
        //! ---------------------------
        //! write the name of the node
        //! ---------------------------
        QVariant data;
        data.setValue(aSimNode->getName());
        tools::writeQVariant(data,outFile);

        //! -------------------
        //! write the nodeType
        //! -------------------
        data.setValue(aSimNode->type());
        tools::writeQVariant(data,outFile);

        //! ------------------------------------------------------------
        //! retrieve the number of properties within the node container
        //! ------------------------------------------------------------
        QVector<QExtendedStandardItem*> vecItem = aSimNode->getPropertyItems();
        int nprop = vecItem.size();

        //! -----------------------------------------------------------------
        //! write the number of properties (drive for the reading operation)
        //! -----------------------------------------------------------------
        data.setValue(nprop);
        tools::writeQVariant(data,outFile);

        for(int k=0; k<vecItem.size(); k++)
        {
            const Property &curProperty = vecItem.at(k)->data(Qt::UserRole).value<Property>();
            Property::writeProperty(outFile,curProperty);
        }
    }
}

//! -------------------
//! function: readVoid
//! details:
//! -------------------
void Property::readVoid(ifstream &in, QStandardItem *item)
{
    cout<<"Property::readVoid()->____READING A QSTANDARDITEM____"<<endl;

    //! ----------------------------------------------------------
    //! read the node name (dual of "write the name of the node")
    //! ----------------------------------------------------------
    QString nodeName = tools::readQVariant(in).toString();
    cout<<"Property::readVoid()->____READING NODE: "<<nodeName.toStdString()<<"____"<<endl;

    //! --------------------------------------------------
    //! read the node type (dual of "write the nodeType")
    //! --------------------------------------------------
    QString nodeType = tools::readQVariant(in).toString();
    cout<<"Property::readVoid()->____NODE TYPE: "<<nodeType.toStdString()<<"____"<<endl;

    //! -------------------------------------------------------------------------
    //! read the number of properties (dual of "write the number of properties")
    //! -------------------------------------------------------------------------
    int nprop = tools::readQVariant(in).toInt();
    cout<<"Property::readVoid()->____NUMBER OF PROPERTIES TO READ: "<<nprop<<"____"<<endl;

    //! --------------------
    //! read the properties
    //! --------------------
    QVector<Property> vecProp;
    for(int k=0; k<nprop; k++)
    {
        Property theProp;
        Property::readProperty(in,theProp);
        vecProp.push_back(theProp);
    }

    SimulationNodeClass *aNode = new SimulationNodeClass(nodeName,SimulationNodeClass::nodeType_NULL,vecProp);
    cerr<<nodeType.toStdString()<<endl;
    aNode->setType1(nodeType);

    QVariant data;
    data.setValue(aNode);
    item->setData(data,Qt::UserRole);
    data.setValue(nodeName);
    item->setData(data,Qt::DisplayRole);
}
