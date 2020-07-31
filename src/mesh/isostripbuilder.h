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
#include <algorithm>

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
class faceTable
{

private:

    std::map<int,std::vector<isoStripPoint>> mdata;

public:

    //! --------------
    //! function: col
    //! details:
    //! --------------
    //std::vector<isoStripPoint> col(int aCol) { return this->value(aCol,std::vector<isoStripPoint>()); }

    //! -----------------------------------------------
    //! append a value at the end of a specific column
    //! -----------------------------------------------
    void appendAtCol(int aCol, const isoStripPoint &theValue)
    {
        std::map<int,std::vector<isoStripPoint>>::iterator it = mdata.find(aCol);
        if(it==mdata.end())
        {
            std::vector<isoStripPoint> v{theValue};
            std::pair<int,std::vector<isoStripPoint>> apair;
            apair.first = aCol;
            apair.second = v;
            mdata.insert(apair);
        }
        else
        {
            (*it).second.push_back(theValue);
        }
    }

    //! -----------------------
    //! function: insertBefore
    //! details:
    //! -----------------------
    bool insertBefore(int aCol, const isoStripPoint &thePoint, const isoStripPoint &theValue)
    {
        std::map<int,std::vector<isoStripPoint>>::iterator it = mdata.find(aCol);
        if(it==mdata.end()) return false;

        std::vector<isoStripPoint>::iterator it_ = std::find((*it).second.begin(),(*it).second.end(),thePoint);

        /*
        //! ----------------------------
        //! diagnostic - can be removed
        //! ----------------------------
        cout<<"faceTable::insertBefore()->____searching for point ("<<thePoint.x<<", "<<thePoint.y<<", "<<thePoint.z<<")____"<<endl;
        for(std::vector<isoStripPoint>::iterator it1= (*it).second.begin(); it1 != (*it).second.end(); it1++)
        {
            const isoStripPoint &aP = *it1;
            cout<<"faceTable::insertBefore()->____available point ("<<aP.x<<", "<<aP.y<<", "<<aP.z<<")____"<<endl;
        }
        //! ---------------
        //! end diagnostic
        //! ---------------
        */

        if(it_==(*it).second.end())
        {
            return false;
        }
        else
        {
            (*it).second.insert(it_,theValue);
            return true;
        }
    }

    //! ----------------------
    //! function: insertAfter
    //! details:
    //! ----------------------
    bool insertAfter(int aCol, const isoStripPoint &thePoint, const isoStripPoint &theValue)
    {
        std::map<int,std::vector<isoStripPoint>>::iterator it = mdata.find(aCol);
        if(it==mdata.end()) return false;

        std::vector<isoStripPoint>::iterator it_ = std::find((*it).second.begin(),(*it).second.end(),thePoint);

        //! ----------------------------
        //! diagnostic - can be removed
        //! ----------------------------
        //cout<<"faceTable::insertAfter()->____searching for point ("<<thePoint.x<<", "<<thePoint.y<<", "<<thePoint.z<<")____"<<endl;
        //for(std::vector<isoStripPoint>::iterator it1= (*it).second.begin(); it1 != (*it).second.end(); it1++)
        //{
            //const isoStripPoint &aP = *it1;
            //cout<<"faceTable::insertAfter()->____available point ("<<aP.x<<", "<<aP.y<<", "<<aP.z<<")____"<<endl;
        //}
        //! ---------------
        //! end diagnostic
        //! ---------------

        if(it_==(*it).second.end())
        {
            return false;
        }
        else
        {
            (*it).second.insert(it_+1,theValue);
            return true;
        }
    }

    //! ----------------------------------------------------------------------------------------------
    //! function: getElements
    //! details:  when building the mesh data sources the automatic renumbering of nodes is activated
    //!           so the field meshElementByCoords.ID is left free for storing the isostrip number
    //! ----------------------------------------------------------------------------------------------
    void getElements(std::vector<meshElementByCoords> &theElements) const
    {
        for(std::map<int,std::vector<isoStripPoint>>::const_iterator it = mdata.cbegin(); it!=mdata.cend(); it++)
        {
            std::pair<int,std::vector<isoStripPoint>> apair = *it;
            int col = apair.first;
            const std::vector<isoStripPoint> &points = apair.second;
            meshElementByCoords anElement;
            switch(points.size())
            {
            case 3:
            {
                //cout<<"#____INSERTING TRIG____"<<endl; anElement.type = TRIG;
                anElement.type = TRIG;
                anElement.ID = col;
                //for(int i=0;i<points.size(); i++)
                //{
                //    const isoStripPoint &P = points[i];
                //    cout<<P.x<<"\t"<<P.y<<"\t"<<P.z<<endl;
                //}
                /*continue; */
            }
                break;
            case 4:
            {
                //cout<<"#____INSERTING QUAD____"<<endl;
                //for(int i=0;i<points.size(); i++)
                //{
                //    const isoStripPoint &P = points[i];
                //    cout<<P.x<<"\t"<<P.y<<"\t"<<P.z<<endl;
                //}
                anElement.type = QUAD;
                anElement.ID = col;
                /*continue; */
            }
                break;
            case 5: { /*cout<<"____INSERTING PENTA____"<<endl; */ anElement.type = PENTA; anElement.ID = col; /*continue; */} break;
            case 6: { /* cout<<"____INSERTING TRIG6____"<<endl; */ anElement.type = TRIG6; anElement.ID = col; /*continue; */} break;
            case 7: { /* cout<<"____INSERTING EPTA____"<<endl; */ anElement.type = EPTA; anElement.ID = col; /*continue; */} break;
            case 8: { /* cout<<"____INSERTING QUAD8____"<<endl; */ anElement.type = QUAD8; anElement.ID = col; /*continue; */} break;
            default: { cout<<"____INSERTING NON STANDARD: NUMBER OF NODES: "<<points.size()<<"____"<<endl; continue; exit(99999999); } break;
            }
            for(int i=0;i<points.size(); i++)
            {
                const isoStripPoint &P = points[i];
                anElement.pointList<<mesh::meshPoint(P.x,P.y,P.z);
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
        std::map<int,std::vector<isoStripPoint>>::const_iterator it = mdata.find(col);
        if(it==mdata.end()) return;
        const std::vector<isoStripPoint> &points = (*it).second;
        switch(points.size())
        {
        case 3: theElement.type = TRIG; break;
        case 4: theElement.type = QUAD; break;
        case 5: theElement.type = PENTA; break;
        case 6: theElement.type = TRIG6; break;
        case 7: theElement.type = EPTA; break;
        case 8: theElement.type = QUAD8; break;
        }
        for(int i = 0; i<points.size(); i++)
        {
            const isoStripPoint &aPoint = points[i];
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
    void getAllElements(const std::vector<faceTable> &vecFaceTables, std::vector<meshElementByCoords> &vecMeshElements);
    void getIsoStripElements(const std::vector<faceTable> &vecFaceTables, std::multimap<int, meshElementByCoords> &meshElementsByIsoStripNb);

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
