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
#include <ng_meshvs_datasourceface.h>

//! ----
//! C++
//! ----
#include <vector>

//! ----
//! OCC
//! ----
#include <TColStd_MapIteratorOfPackedMapOfInteger.hxx>

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
class faceTable: public std::vector<std::vector<isoStripPoint>>
{

public:

    //! --------------
    //! function: col
    //! details:
    //! --------------
    std::vector<isoStripPoint> col(int aCol) { return this->at(aCol); }

    //! -----------------------
    //! function: insertBefore
    //! details:
    //! -----------------------
    void insertBefore(int aCol, const isoStripPoint &thePoint, const isoStripPoint &theValue)
    {
        std::vector<isoStripPoint> *points = &(this->at(aCol));
        std::vector<isoStripPoint>::iterator pos = std::find(points->begin(), points->end(), thePoint);
        if(pos == points->end()) return;
        points->insert(pos,theValue);
    }

    //! ----------------------
    //! function: insertAfter
    //! details:
    //! ----------------------
    void insertAfter(int aCol, const isoStripPoint &thePoint, const isoStripPoint &theValue)
    {
        std::vector<isoStripPoint> *points = &(this->at(aCol));
        std::vector<isoStripPoint>::iterator pos = std::find(points->begin(), points->end(), thePoint);
        if(pos == points->end()) return;
        points->insert(pos+1,theValue);
    }

    //! ----------------------
    //! function: getElements
    //! details:
    //! ----------------------
    void getElements(std::vector<meshElementByCoords> &theElements) const
    {
        for(int coln = 0; coln <this->size(); coln++)
        {
            meshElementByCoords anElement;
            const std::vector<isoStripPoint> &faceElement = this->at(coln);
            switch(faceElement.size())
            {
            case 3:
            {
                cout<<"____INSERTING TRIG____"<<endl;
                anElement.type = TRIG;
            }
                break;
            case 4:
            {
                cout<<"____INSERTING QUAD____"<<endl;
                anElement.type = QUAD;
            }
                break;
            case 5:
            {
                cout<<"____INSERTING PENTA____"<<endl;
                anElement.type = PENTA;
            }
                break;
            case 6:
            {
                cout<<"____INSERTING TRIG6____"<<endl;
                anElement.type = TRIG6;
            }
                break;
            case 7:
            {
                cout<<"____INSERTING EPTA____"<<endl;
                anElement.type = EPTA;
            }
                break;
            case 8:
            {
                cout<<"____INSERTING QUAD8____"<<endl;
                anElement.type = QUAD8;
            }
                break;
            default:
            {
                cout<<"____INSERTING NON STANDARD: NUMBER OF NODES: "<<faceElement.size()<<"____"<<endl;
                exit(99999999);
            }
                break;
            }
            for(int i = 0; i<faceElement.size(); i++)
            {
                const isoStripPoint &aPoint = faceElement[i];
                anElement.pointList<<mesh::meshPoint(aPoint.x,aPoint.y,aPoint.z);
            }
            theElements.push_back(anElement);
        }
    }

    //! ------------------------------
    //! function: getElementsOfColumn
    //! details:
    //! ------------------------------
    void getElementOfColumn(int col, meshElementByCoords &theElement) const
    {
        if(col<0 || col>this->size()-1) return;

        const std::vector<isoStripPoint> &faceElement = this->at(col);
        switch(faceElement.size())
        {
        case 3: theElement.type = TRIG; break;
        case 4: theElement.type = QUAD; break;
        case 6: theElement.type = TRIG6; break;
        //case 6: theElement.type = HEXAG; break;    // require the definition of another type of element
        }
        for(int i = 0; i<faceElement.size(); i++)
        {
            const isoStripPoint &aPoint = faceElement[i];
            theElement.pointList<<mesh::meshPoint(aPoint.x,aPoint.y,aPoint.z);
        }
    }
};



class MeshVS_DataSource;
class isoStripBuilder : public QObject
{
    Q_OBJECT

public:

    explicit isoStripBuilder(const occHandle(MeshVS_DataSource) &aMeshDS = occHandle(MeshVS_DataSource)(),
                             QObject *parent = 0);
    void setMeshDataSource(const occHandle(MeshVS_DataSource) &aMeshDS);
    void setNbStrip(int NbStrips);
    void setValues(const QMap<int,double> &values);
    void setIsoStrips(const std::vector<isoStrip> &theIsoStrips);
    bool perform(std::vector<meshElementByCoords> &vecMeshElements);
    void getAllElements(const std::vector<faceTable> &vecFaces, std::vector<meshElementByCoords> &vecMeshElements);


private:

    occHandle(MeshVS_DataSource) myMeshDS;
    QMap<int,double> myValues;
    std::vector<isoStrip> myIsoStrips;

private:

    void classifyNodes(QMap<int, int> &pointToIsoStrip);
    void pointCoord(double *c, int globalNodeID);

private:

};

#endif // ISOSTRIPBUILDER_H
