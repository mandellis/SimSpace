#ifndef WRITESOLVERFILECLASS_H
#define WRITESOLVERFILECLASS_H

//! ----------------
//! custom includes
//! ----------------
#include "mydefines.h"
#include "qextendedstandarditem.h"
#include "postobject.h"
#include <indexedmapofmeshdatasources.h>
#include <qprogressindicator.h>
#include <meshelement2d.h>

//! ----
//! C++
//! ----
#include <iostream>
#include <fstream>
#include <map>

//! ----
//! OCC
//! ----
#include <TopTools_ListOfShape.hxx>

//! ---
//! Qt
//! ---
#include <QObject>
#include <QTreeView>

class QProgressDialog;
class simulationDataBase;

class writeSolverFileClass: public QObject
{

private:

    simulationDataBase *myDB;
    QExtendedStandardItem *mySimulationRoot;
    QString myFileName;
    ofstream myInputFile;
    ofstream myMesh;

    QProgressDialog *myProgressDialog;
    QProgressIndicator *myProgressIndicator;
    QTreeView *myMainTree;

    //std::map<int, std::map<meshElement2D,std::vector<std::pair<int,std::string>>>> bigMap;
    std::map<int, std::map<meshElement2D,std::vector<std::pair<int,int>>>> bigMap;


private:

    void writeGapElement(const IndexedMapOfMeshDataSources &anIndexedMapOfFaceMeshDS,
                        QString SetName,double K,double F);

    void createNodalSet(const IndexedMapOfMeshDataSources &anIndexedMapOfFaceMeshDS,
                        std::vector<int> &theNodeIDs);

    void writeNodalSet(QString SetName,
                       const IndexedMapOfMeshDataSources &anIndexedMapOfFaceMeshDS);

    void writeElementSurface(QString SetName,
                             const IndexedMapOfMeshDataSources &anIndexedMapOfMeshDataSources);

    void writeElementSet(QVector<GeometryTag> vecLoc,QList<QString> &bodyNameList);
    void writeNodesAndElements(QString aName,QMap<int,QList<int>> &nodeListByBody);

    QExtendedStandardItem *ItemFromScope(const TopoDS_Shape &aShape);
    QExtendedStandardItem* getTreeItem(SimulationNodeClass::nodeType theNodeType);

    //! -----------------------------------------------------------
    //! create element surface (uses face to element connectivity)
    //! -----------------------------------------------------------
    void createElementSurface(std::vector<int> &theElementIDs,
                              std::vector<int> &theFaceNumbers,
                              const IndexedMapOfMeshDataSources &anIndexedMapOfFaceMeshDS);

    //! ---------------------------------------------
    //! create the element surface for the bolt case
    //! ---------------------------------------------
    void writeElementSurfaceBolt(QString SetName, int boltBodyIndex,
                                 double a, double b, double c, double d);

    void writeDload(double aLoad, QString aName);
    void writeFilm(double aLoad, QString aName,double refTemperature);
    void writeDflux(double aLoad, QString aName);

    void writeTemperatureHistory(postObject pObject, QString tName);

    //! progress bar info
    void progressDialogInfo();

    //! for testing purposes
    void clock();
    inline QString currentTime() { return QTime::currentTime().toString("mm:ss:zzz"); }

    //! -------------------------------------
    //! function: itemNameClearSpaces
    //! details:  remove spaces from QString
    //! -------------------------------------
    QString itemNameClearSpaces(const QString &aName)
    {
        QString ss=aName;
        QStringList s=ss.split(" ");
        ss.clear();
        for (int i=0;i<s.length();i++)
        {
            ss=ss+s.at(i);
        }
        return ss;
    }

public:

    //! constructor
    writeSolverFileClass(simulationDataBase *aDB, QExtendedStandardItem* aSimulationRoot, QObject *parent=0);

    //! set file name
    void setName(const QString &aName) { myFileName = aName; }

    //! set main tree
    void setSimTree(QTreeView *aTreeView) { myMainTree = aTreeView; }

    //! perform
    bool perform();

    //! set progress indicator
    void setProgressIndicator(QProgressIndicator *aProgressIndicator);

private:

    std::vector<std::string> vecMatNames;
};

#endif // WRITESOLVERFILECLASS_H
