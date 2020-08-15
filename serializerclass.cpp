//! ----------------
//! custom includes
//! ----------------
#include "serializerclass.h"
#include "simulationnodeclass.h"
#include "qextendedstandarditem.h"
#include "tools.h"
#include "customtablemodel.h"
#include "postobject.h"

//! ----
//! OCC
//! ----
#include <TopTools_ListIteratorOfListOfShape.hxx>

//! ---
//! Qt
//! ---
#include <QFile>
#include <QDir>

//! ----
//! C++
//! ----
#include <ostream>
#include <Standard_OStream.hxx>

//! --------------------
//! function: setItem()
//! details:
//! --------------------
void serializerClass::setItem(QExtendedStandardItem *anItem)
{
    myItem = anItem;
}

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
serializerClass::serializerClass(QObject* parent): QObject(parent)
{
    ;
}

//! -------------------
//! function: savePath
//! details:
//! -------------------
void serializerClass::setSavingDirPath(const QString &path)
{
    myPath = path;
}

//! ---------------------------------------------------
//! function: serialize
//! details:  [A] write the name (QString)
//!           [B] write the type (enum) - as QString
//!           [C] write the number of properties (int)
//!           [D] write the properties
//! ---------------------------------------------------
void serializerClass::serialize(int n)
{
    static int N;
    int k;

    if(n==-1)
    {
        N++;
        k = N;
    }
    else k=n;

    //! --------------------
    //! build the node name
    //! --------------------
    QString itemName = myItem->data(Qt::DisplayRole).toString();
    QString nodeName = myItem->data(Qt::UserRole).value<SimulationNodeClass*>()->getName();
    QString nodeType = myItem->data(Qt::UserRole).value<SimulationNodeClass*>()->type();
    QString nodeFileName = (QString("%1").arg(k))+"_"+itemName+"_"+nodeName+"_"+nodeType;
    QString fullNodeFileName = myPath+QString("\\")+nodeFileName;

    //! --------------------------
    //! open the file for writing
    //! --------------------------
    ofstream outFile;
    outFile.open(fullNodeFileName.toStdString());

    //! --------------------------------
    //! retrieve the node from the item
    //! --------------------------------
    SimulationNodeClass *node = myItem->data(Qt::UserRole).value<SimulationNodeClass*>();

    //! ---------------------------
    //! write the name of the node
    //! ---------------------------
    QVariant data;
    data.setValue(node->getName());
    tools::writeQVariant(data,outFile);

    //! -------------------
    //! write the nodeType
    //! -------------------
    data.setValue(node->type());
    tools::writeQVariant(data,outFile);

    //! ------------------------------------------------------------
    //! retrieve the number of properties within the node container
    //! ------------------------------------------------------------
    QVector<QExtendedStandardItem*> vecItem = node->getPropertyItems();
    int nprop = vecItem.size();
    for(int k=0; k<vecItem.size(); k++)
    {
        const Property &theCurProperty = vecItem.at(k)->data(Qt::UserRole).value<Property>();
        if(theCurProperty.getName()=="Graphic object") nprop--;
    }
    //! -----------------------------------------------------------------
    //! write the number of properties (drive for the reading operation)
    //! -----------------------------------------------------------------
    data.setValue(nprop);
    tools::writeQVariant(data,outFile);

    //! ------------------------------------------------------
    //! scan the list of properties within the node container
    //! ------------------------------------------------------
    for(int k=0; k<vecItem.size(); k++)
    {
        const Property &theCurProperty = vecItem.at(k)->data(Qt::UserRole).value<Property>();
        //! ------------------------------------------------------
        //! case of a single post object present in the main tree
        //! ------------------------------------------------------
        //if(theCurProperty.getData().canConvert<postObject>())
        //{
        //    cout<<"serializerClass::serialize()->____post object found____"<<endl;
        //    theCurProperty.getData().value<postObject>().write(outFile);
        //    cout<<"serializerClass::serialize()->____post object written____"<<endl;
        //    continue;
        //}
        /*
        //! -----------------------------------------------------------------------------------
        //! if the "Property" contains a void pointer, its "writeProperty" method does nothing
        //! the content of the void pointer is serialized here
        //! -----------------------------------------------------------------------------------
        if(theCurProperty.getData().canConvert<void*>())
        {
            //! ------------------------------------------------------------------------
            //! it happens when the simulation node contains a "Coordinate system" or
            //! a "Named selection", or a "Remote point" through an indirect reference,
            //! i.e. through a void pointer
            //! ------------------------------------------------------------------------
            void *p = theCurProperty.getData().value<void*>();
            QExtendedStandardItem *aTreeItem = static_cast<QExtendedStandardItem*>(p);
            SimulationNodeClass *aSimNode = aTreeItem->data(Qt::UserRole).value<SimulationNodeClass*>();
            SimulationNodeClass::nodeType theNode = aSimNode->getType();

            if(theNode==SimulationNodeClass::nodeType_namedSelectionGeometry ||
                    theNode==SimulationNodeClass::nodeType_coordinateSystem ||
                    theNode==SimulationNodeClass::nodeType_remotePoint ||
                    theNode==SimulationNodeClass::nodeType_importedBodyScalar ||
                    theNode==SimulationNodeClass::nodeType_thermalAnalysisSettings ||
                    theNode==SimulationNodeClass::nodeType_meshPrismaticLayer)
            {
                //! write the name of the node
                QVariant data;
                data.setValue(aSimNode->getName());
                tools::writeQVariant(data,outFile);

                //! write the nodeType
                data.setValue(aSimNode->type());
                tools::writeQVariant(data,outFile);

                //! retrieve the number of properties within the node container
                QVector<QExtendedStandardItem*> vecItem = aSimNode->getPropertyItems();
                int nprop = vecItem.size();

                //! write the number of properties (drive for the reading operation)
                data.setValue(nprop);
                tools::writeQVariant(data,outFile);

                for(int k=0; k<vecItem.size(); k++)
                {
                    const Property &curProperty = vecItem.at(k)->data(Qt::UserRole).value<Property>();
                    if(curProperty.getData().canConvert<postObject>())
                    {
                        curProperty.getData().value<postObject>().write(outFile);
                        continue;
                    }
                    Property::writeProperty(outFile,curProperty);
                }
            }
            continue;
        }
        */

        //! ---------------------------------------------------------
        //! "Graphic object" should not be persistent (i.e. markers)
        //! ---------------------------------------------------------
        if(theCurProperty.getName()=="Graphic object") continue;
        Property::writeProperty(outFile,theCurProperty);
    }

    //! ---------------------------------------------------------------
    //! if the item is "Analysis settings" write also the tabular data
    //! ---------------------------------------------------------------
    if(node->isAnalysisSettings())
    {
        cout<<"\n\n* WRITING THE TABULAR DATA IN A SEPARATE FILE"<<endl;
        QDir curDir;
        curDir.cd(myPath);
        curDir.mkdir(QString("Tabular data"));
        curDir.cd("Tabular data");

        ofstream outFile_TabularData;
        QString timeTag = node->getPropertyValue<QString>("Time tag");
        QString tabularDataFileName = curDir.absoluteFilePath(QString("Tabular data_")+timeTag+".txt");
        outFile_TabularData.open(tabularDataFileName.toStdString());

        //! retrieve the tabular data
        CustomTableModel *tabularDataModel = node->getTabularDataModel();

        //! write the number of loads
        int Ncol = tabularDataModel->columnCount();
        QVariant n;
        n.setValue(Ncol);
        tools::writeQVariant(n, outFile_TabularData);

        //! ------------------------
        //! write the loads on disk
        //! ------------------------
        cout<<"\n* ------------------------------------------ *"<<endl;
        cout<<"* NUMBER OF COLUMNS: "<<Ncol<<endl;
        for(int i=0; i<Ncol; i++)
        {
            const load &theLoad = tabularDataModel->getColumn(i);
            cout<<"* WRITING LOAD ->"<<theLoad.getLoadType().toStdString()<<endl;
            //!load::writeLoad(tabularDataModel->getColumn(i),outFile_TabularData);
            theLoad.write(outFile_TabularData);
        }
        cout<<"* ------------------------------------------ *"<<endl<<endl;
        outFile_TabularData.close();
    }
    outFile.close();
}
