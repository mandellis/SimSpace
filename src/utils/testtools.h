#ifndef TESTTOOLS_H
#define TESTTOOLS_H

//! custom includes
#include "src/main/mydefines.h"

//! Qt
#include <QString>

//! OCC
#include <MeshVS_DataSource.hxx>

//! C++
#include <iostream>

using namespace std;

class testTools
{

public:

    testTools();
    static void fakeScalarData(const occHandle(MeshVS_DataSource) &aDS, const QString &fileName, int Npoints=1000);
    static void fakeScalarData1(const occHandle(MeshVS_DataSource) &aDS, const QString &fileName);

};

#endif // TESTTOOLS_H
