#ifndef MESHELEMENT2D_H
#define MESHELEMENT2D_H

//! ----------------
//! custom includes
//! ----------------
#include <elementtypes.h>

//! ---
//! Qt
//! ---
#include <QList>

//! ----
//! C++
//! ----
#include <algorithm>
#include "src/utils/hash_c.h"

//! -------------
//! mesh element
//! -------------
struct meshElement2D
{
    int ID;
    QList<int> nodeIDs;
    ElemType type;

    //! ------------
    //! constructor
    //! ------------
    meshElement2D(const QList<int> theNodeIDs = QList<int>(), int theID = -1, ElemType aType = NULL_ELEMENT):
        ID(theID),nodeIDs(theNodeIDs),type(aType)
    {
        ;
    }

    //! -----------------
    //! copy constructor
    //! -----------------
    meshElement2D(const meshElement2D &rhs)
    {
        ID = rhs.ID;
        nodeIDs.clear();
        nodeIDs<<rhs.nodeIDs;
        type=rhs.type;
    }

    //! ----------------------
    //! operator = assignment
    //! ----------------------
    inline meshElement2D operator = (const meshElement2D &rhs)
    {
        ID = rhs.ID;
        nodeIDs.clear();
        nodeIDs<<rhs.nodeIDs;
        type = rhs.type;
        return *this;
    }

    //! --------------------------------
    //! operator ==
    //! identity, not aware of ordering
    //! --------------------------------
    bool operator == (const meshElement2D &rhs) const
    {
        //! size check
        if(nodeIDs.size()!=rhs.nodeIDs.size()) return false;

        QList<int> tmp(nodeIDs);
        QList<int> tmp1(rhs.nodeIDs);
        std::sort(tmp.begin(), tmp.end());
        std::sort(tmp1.begin(), tmp1.end());
        int Nb = tmp.length();
        for(int i=0; i<Nb; i++) if(tmp.at(i)!=tmp1.at(i)) return false;
        return true;
    }

    //! ------------------------------------------
    //! sort()
    //! example: {10,3,5,6,1,8} => {1,8,10,3,5,6}
    //! ------------------------------------------
    inline void sort()
    {
        QList<int> reorderedNodeIDs;
        int minNodeID = 1e10;
        int NbNodes = nodeIDs.length();
        for(int i=0; i<NbNodes; i++)
        {
            int curNodeID = nodeIDs.at(i);
            if(curNodeID<minNodeID) minNodeID = curNodeID;
        }
        int minIndex = nodeIDs.indexOf(minNodeID);

        int offset=0;
        for(int i=0; i<NbNodes; i++)
        {
            int k=i+minIndex;
            int nodeID;
            if(k<NbNodes)
            {
                offset++;
                nodeID = nodeIDs.at(k);
            }
            else
            {
                nodeID = nodeIDs.at(i-offset);
            }
            reorderedNodeIDs<<nodeID;
        }
        nodeIDs.clear();
        nodeIDs<<reorderedNodeIDs;
    }

    //! ---------------------------
    //! operator <
    //! not aware of node ordering
    //! ---------------------------
    bool operator < (const meshElement2D &rhs) const
    {
        size_t seed1, seed2;
        seed1 = seed2 = 0;
        int Nb = nodeIDs.size();
        QList<int> tmp(nodeIDs);
        QList<int> tmp1(rhs.nodeIDs);
        std::sort(tmp.begin(), tmp.end());
        std::sort(tmp1.begin(), tmp1.end());
        for(int i=0; i<Nb; i++)
        {
            hash_c<int>(seed1,tmp[i]);
            hash_c<int>(seed2,tmp1[i]);
        }
        if(seed1<seed2) return true; return false;
    }
};

//! ---------------
//! qHash function
//! ---------------
inline uint qHash(const meshElement2D &me, uint seed = 0)
{
    size_t seed1 = seed;
    for(int i=0; i<me.nodeIDs.length(); i++) hash_c<int>(seed1,me.nodeIDs.at(i));
    return uint(seed1);
}

#endif // MESHELEMENT2D_H
