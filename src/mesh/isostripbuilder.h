#ifndef ISOSTRIPBUILDER_H
#define ISOSTRIPBUILDER_H

//! ---
//! Qt
//! ---
#include <QObject>
#include <QMap>

//! ----------------
//! custom includes
//! ----------------
#include "occhandle.h"
#include <isostrip.h>

//! ----
//! C++
//! ----
#include <vector>

//! ----
//! OCC
//! ----
#include <TColStd_MapIteratorOfPackedMapOfInteger.hxx>

typedef TColStd_MapIteratorOfPackedMapOfInteger MIPMI;

class MeshVS_DataSource;

class isoStripBuilder : public QObject
{
    Q_OBJECT

public:

    explicit isoStripBuilder(const occHandle(MeshVS_DataSource) &aMeshDS = occHandle(MeshVS_DataSource)(),
                             int NbStrips = 5,
                             QObject *parent = 0);
    void setMeshDataSource(const occHandle(MeshVS_DataSource) &aMeshDS);
    void setNbStrip(int NbStrips);
    void setValues(QMap<int,double> *values);
    void setIsoStrips(std::vector<isoStrip> *theIsoStrips);
    bool perform();

private:

    occHandle(MeshVS_DataSource) myMeshDS;
    int myNbStrips;
    QMap<int,double> *myValues;
    std::vector<isoStrip> *myIsoStrips;

    struct point
    {
        double x,y,z;
        point(double ax=0, double ay=0, double az=0) { x = ax; y = ay; z = az; }
        point(const point &aP) { x = aP.x; y = aP.y; z = aP.z; }
        point operator =(const point &aP) { x = aP.x; y = aP.y; z = aP.z; return *this;}
        point sum(const point &P1, const point &P2)
        {
            point P;
            P.x = P1.x + P2.x; P.y = P1.y + P2.y; P.z = P1.z + P2.z; return P;
        }
        point multiply(const point &P1, double k)
        {
            point P;
            P.x = P1.x*k; P.y = P1.y*k; P = P1.z*k;
            return P;
        }
    };

private:

    void classifyNodes(QMap<int, int> &nodeToIsoStrip);

signals:

public slots:


};

#endif // ISOSTRIPBUILDER_H
