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
#include "hash_c.h"
#include <isostrip.h>
#include <meshelementbycoords.h>
#include <mesh.h>

//! ----
//! C++
//! ----
#include <vector>

//! ----
//! OCC
//! ----
#include <TColStd_MapIteratorOfPackedMapOfInteger.hxx>

//! ------
//! point
//! ------
struct point
{
    double x,y,z;
    double val;
    point(double ax=0, double ay=0, double az=0, double aVal = 0) { x = ax; y = ay; z = az; val = aVal; }
    point(const point &aP) { x = aP.x; y = aP.y; z = aP.z; val = aP.val; }
    point operator =(const point &aP) { x = aP.x; y = aP.y; z = aP.z; val = aP.val; return *this;}
    bool operator == (const point &aP) { if(x == aP.x && y == aP.y && z == aP.z) return true; return false; }
    bool operator < (const point &aP)
    {
        std::size_t seed1, seed2; seed1=seed2=0;
        hash_c<double>(seed1,x); hash_c<double>(seed1,y); hash_c<double>(seed1,z);
        hash_c<double>(seed2,aP.x); hash_c<double>(seed2,aP.y); hash_c<double>(seed2,aP.z);
        if(seed1<seed2) return true;
        else return false;
    }
};


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

private:

    void classifyNodes(QMap<int, int> &pointToIsoStrip);
    void pointCoord(double *c, int globalNodeID);

signals:

public slots:

};

//! ---------------------------------------
//! class: faceTable
//!                      D(300)     C(50)
//!  (0) [0-100]          ____________
//!  (1) [100-200]       |            |
//!  (2) [200-300]       |            |
//!  (3) [300-400]       |            |
//!                       ------------
//!   | 0 | 1 | 2 | 3 |  A(50)      B(100)
//!   -----------------
//!   | A | B | D | D |
//!   | B |   |   |   |
//!   | C |   |   |   |
//!   |   |   |   |   |
//!
//! ---------------------------------------
class faceTable: public std::vector<std::vector<point>>
{

public:

    //! --------------
    //! function: col
    //! details:
    //! --------------
    std::vector<point> col(int aCol) { return this->at(aCol); }

    //! -----------------------
    //! function: insertBefore
    //! details:
    //! -----------------------
    void insertBefore(int aCol, const point &thePoint, const point &theValue)
    {
        std::vector<point> *points = &(this->at(aCol));
        std::vector<point>::iterator pos = std::find(points->begin(), points->end(), thePoint);
        if(pos == points->end()) return;
        points->insert(pos,theValue);
    }

    //! ----------------------
    //! function: insertAfter
    //! details:
    //! ----------------------
    void insertAfter(int aCol, const point &thePoint, const point &theValue)
    {
        std::vector<point> *points = &(this->at(aCol));
        std::vector<point>::iterator pos = std::find(points->begin(), points->end(), thePoint);
        if(pos == points->end()) return;
        points->insert(pos+1,theValue);
    }

    //! ----------------------
    //! function: getElements
    //! details:
    //! ----------------------
    void getElements(std::vector<meshElementByCoords> theElements) const
    {
        for(int coln = 0; coln <this->size(); coln++)
        {
            meshElementByCoords anElement;
            const std::vector<point> &faceElement = this->at(coln);
            switch(faceElement.size())
            {
            case 3: anElement.type = TRIG; break;
            case 4: anElement.type = QUAD; break;
            case 6: anElement.type = TRIG6; break;
            //case 6: anElement.type = HEXAG; break;    // require the definition of another type of element
            }
            for(int i = 0; i<faceElement.size(); i++)
            {
                const point &aPoint = faceElement.at(i);
                anElement.pointList<<mesh::meshPoint(aPoint.x,aPoint.y,aPoint.z);
            }
            theElements.push_back(anElement);
        }
    }
};

#endif // ISOSTRIPBUILDER_H
