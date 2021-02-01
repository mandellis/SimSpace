#ifndef FRDREADER_H
#define FRDREADER_H

//! ---
//! Qt
//! ---
#include <QObject>
#include <QDir>
#include <QString>
#include <QMap>
#include <QMetaType>

//! ----
//! C++
//! ----
#include <vector>
#include <iostream>
#include <fstream>

//! ----------------
//! custom includes
//! ----------------
#include "src/main/mydefines.h"
#include "meshdatabase.h"

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

private:

    bool readResults(ifstream &is, QString pathifstream);

public:

    void setResultFile (const QString &FRDfileName) { myFRDfileName = FRDfileName; }
    bool perform(/*QMap<double,QVector<int>> &discreteTimeMap*/);
};

#endif // FRDREADER_H
