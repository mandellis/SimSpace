#ifndef FRDREADER_H
#define FRDREADER_H

//! Qt
#include <QObject>
#include <QDir>
#include <QString>
#include <QMap>
#include <QMetaType>

//! C++
#include <vector>
#include <iostream>
#include <fstream>

//! custom includes
#include "mydefines.h"
#include <meshdatabase.h>
using namespace std;

class FrdReader: public QObject
{
public:

    FrdReader(const QString &FRDfileName, QObject *parent=0);
    FrdReader(QObject *parent=0);

    virtual ~FrdReader()
    {
        cout<<"FrdReader::~FrdReader()->____DESTRUCTOR CALLED____"<<endl;
    }

private:

    QString myPath;
    QString myFRDfileName;
    ofstream myResultFile;

    //! the mesh database
    //meshDataBase *myDB;

private:

    bool readResults(std::ifstream &is,QString pathifstream, int &sb,int &st,double &t);
public:

    //void setMeshDataBase(meshDataBase *mDB) { myDB = mDB; }
    void setResultFile (const QString &FRDfileName) { myFRDfileName = FRDfileName; }
    bool splitFrdFile(QMap<double,QVector<int>> &discreteTimeMap);
};

#endif // FRDREADER_H
