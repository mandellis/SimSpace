#ifndef OFWRITE_H
#define OFWRITE_H

//! ----------------
//! custom includes
//! ----------------
#include "src/main/mydefines.h"
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

class ofwrite: public QObject
{

private:

    simulationDataBase *myDB;
    QExtendedStandardItem *mySimulationRoot;
    QString myFileDir;
    //ofstream myInputFile;
    //ofstream myMesh;

    QProgressDialog *myProgressDialog;
    QProgressIndicator *myProgressIndicator;
    QTreeView *myMainTree;

    //std::map<int, std::map<meshElement2D,std::vector<std::pair<int,std::string>>>> bigMap;
    //std::map<int, std::map<meshElement2D,std::vector<std::pair<int,int>>>> bigMap;

    //int totalNumberOfNodes;
    //int totalNumberOfElements;

private:

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

    void makeHeader(ifstream &is);
    void makeSeparator(ifstream &is);
    void makeEndFile(ifstream &is);

public:

    //! constructor
    ofwrite(simulationDataBase *aDB, QExtendedStandardItem* aSimulationRoot, QObject *parent=0);

    //! set file name
    void setName(const QString &aName) { myFileDir = aName; }

    //! set main tree
    void setSimTree(QTreeView *aTreeView) { myMainTree = aTreeView; }

    //! perform
    bool perform();

    //! set progress indicator
    void setProgressIndicator(QProgressIndicator *aProgressIndicator);

private:

    std::vector<std::string> vecMatNames;
};

#endif // OFWRITE_H
