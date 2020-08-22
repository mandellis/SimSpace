#ifndef ISOSTRIPBUILDER_H
#define ISOSTRIPBUILDER_H

#define USE_STDLIBRARY

#ifdef USE_STDLIBRARY
#define myMultiMap multimapex
#define isoStripList std::vector<int>
#endif
#ifndef USE_STDLIBRARY
#define myMultiMap QMap
#define isoStripList QList<int>
#endif

//! ---
//! Qt
//! ---
#include <QObject>

//! ----------------
//! custom includes
//! ----------------
#include "occhandle.h"
#include <isostrip.h>
#include <mesh.h>
#include <meshelementbycoords.h>
#include <ng_meshvs_datasourceface.h>

//! ----
//! C++
//! ----
#include <vector>
#include <algorithm>
#include <map>
#include <unordered_map>

//! ----
//! OCC
//! ----
#include <TColStd_MapIteratorOfPackedMapOfInteger.hxx>

template<class K, class V>
class multimapex: public std::unordered_multimap<K,V>
{
public:

    //! ------------
    //! insertMulti
    //! ------------
    inline void insertMulti(const K &key, const V &value) { this->insert(std::make_pair(key,value)); }

    //! -------------------------------------
    //! retrieve the list of values of a key
    //! -------------------------------------
    std::vector<V> values(const K &key)
    {
        /* linear search O(N)
        std::vector<V> values;
        auto itr1 = this->lower_bound(key);
        auto itr2 = this->upper_bound(key);
        while (itr1 != itr2)
        {
            if (itr1 -> first == key) values.push_back(itr1->second);
            itr1++;
        }
        return values;
        */
        //! ---------------------------------------------------------------------------------------------
        //! https://thispointer.com/finding-all-values-for-a-key-in-multimap-using-equals_range-example/
        //! ---------------------------------------------------------------------------------------------
        std::vector<V> values;
        std::pair<std::unordered_multimap<K,V>::iterator, std::unordered_multimap<K,V>::iterator> r = this->equal_range(key);
        for (std::unordered_multimap<K,V>::iterator it = r.first; it != r.second; it++) values.push_back(it->second);
        return values;
    }
};

//! -----------------
//! class: faceTable
//! -----------------
class faceTable
{

private:

    std::map<int,std::vector<isoStripPoint>> mdata;

public:

    //! ----------
    //! getColumn
    //! ----------
    std::vector<isoStripPoint> getColumn(int col) const
    {
        std::vector<isoStripPoint> vecPoints;
        std::map<int,std::vector<isoStripPoint>>::const_iterator it = mdata.find(col);
        if(it!=mdata.cend())
        {
            const std::vector<isoStripPoint> &vp = it->second;
            for(int k=0; k<vp.size(); k++) vecPoints.push_back(vp[k]);
        }
        return vecPoints;
    }

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
            case 3: { anElement.type = TRIG; anElement.ID = col; /*continue; */ } break;
            case 4: { anElement.type = QUAD; anElement.ID = col; /*continue; */ } break;
            case 5: { anElement.type = PENTA; anElement.ID = col; /*continue; */} break;
            case 6: { anElement.type = TRIG6; anElement.ID = col; /*continue; */} break;
            case 7: { anElement.type = EPTA; anElement.ID = col; /*continue; */} break;
            case 8: { anElement.type = QUAD8; anElement.ID = col; /*continue; */} break;
            default: { /*cout<<"____INSERTING NON STANDARD: NUMBER OF NODES: "<<points.size()<<"____"<<endl; */ continue; } break;
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

    void setValues(const std::map<int,double> &values);
    void setIsoStrips(const std::vector<isoStrip> &theIsoStrips);

    bool perform(std::vector<meshElementByCoords> &vecMeshElements);
    bool perform1(std::vector<meshElementByCoords> &vecMeshElements);   // a - for the moment - failed attempt to optimize

    void getIsoStripElements(const std::vector<faceTable> &vecFaceTables, std::multimap<int, meshElementByCoords> &meshElementsByIsoStripNb);
    bool performIsoSurface(int NbLevels, std::vector<meshElementByCoords> &vecMeshElements, std::map<int,int> &mapElementLevel, int position = -1);

private:

    occHandle(MeshVS_DataSource) myMeshDS;
    std::map<int,double> myValues;
    std::vector<isoStrip> myIsoStrips;
    std::map<int,std::vector<faceTable>> myMapElementFaceTables;
    std::vector<faceTable> myVecFaceTables;

private:

    void classifyNodes(myMultiMap<int, int> &pointToIsoStrip);
    void pointCoord(double *c, int globalNodeID);
    void getAllElements(const std::vector<faceTable> &vecFaceTables, std::vector<meshElementByCoords> &vecMeshElements);
    bool computeFaceTables();

private:

};

#endif // ISOSTRIPBUILDER_H
