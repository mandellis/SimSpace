#ifndef MESH_H
#define MESH_H

//! Qt
#include <QList>

//! C++
#include <algorithm>

//! custom includes
#include "hash_c.h"
#include <elementtypes.h>
#include <iostream>

namespace mesh
{
//! ---------------
//! tolerant point
//! ---------------
struct tolerantPoint
{
    double x,y,z;
    double eps = 1e-3;
    tolerantPoint(double aX=0, double aY=0, double aZ=0, double aEps = 1e-3):x(aX),y(aY),z(aZ),eps(aEps){;}
    tolerantPoint(const tolerantPoint &rhs) { x = rhs.x; y = rhs.y; z = rhs.z; eps = rhs.eps; }

    bool operator == (const tolerantPoint &aP) const
    {
        double d = sqrt(pow(x-aP.x,2)+pow(y-aP.y,2)+pow(z-aP.z,2));
        if(d<=eps) return true;
        return false;
    }
};

//! -----------
//! mesh point
//! -----------
struct meshPoint
{
    double x,y,z;
    int ID;

    //! -------------
    //! constructors
    //! -------------
    meshPoint(){ x=y=z=0; ID = -1; }
    meshPoint(double xP, double yP, double zP, int pointID = -1):x(xP),y(yP),z(zP),ID(pointID){;}

    //! ------------
    //! operator ==
    //! ------------
    inline bool operator == (const meshPoint &rhs) const
    {
        //if(x==rhs.x && y==rhs.y && z == rhs.z && ID == rhs.ID) return true;
        if(x==rhs.x && y==rhs.y && z == rhs.z) return true;
        return false;
    }

    //! -----------
    //! operator =
    //! -----------
    inline meshPoint operator = (const meshPoint &rhs)
    {
        x = rhs.x; y = rhs.y; z = rhs.z;
        ID = rhs.ID;
        return *this;
    }

    //! -----------
    //! operator <
    //! -----------
    bool operator < (const meshPoint &rhs) const
    {
        std::size_t seed1, seed2; seed1=seed2=0;
        hash_c<double>(seed1,x); hash_c<double>(seed1,y); hash_c<double>(seed1,z);
        hash_c<double>(seed2,rhs.x); hash_c<double>(seed2,rhs.y); hash_c<double>(seed2,rhs.z);
        if(seed1<seed2) return true;
        else return false;
    }
};

//! -------------
//! mesh segment
//! -------------
struct meshSegment
{
    QList<int> nodeIDs;

    //! -----
    //! sort
    //! -----
    inline void sort() { std::sort(nodeIDs.begin(), nodeIDs.end()); }

    //! -----------
    //! operator =
    //! -----------
    inline meshSegment operator = (const meshSegment &rhs)
    {
        for(int i=0; i<rhs.nodeIDs.length(); i++) nodeIDs<<rhs.nodeIDs[i];
        return *this;
    }

    //! ------------
    //! operator ==
    //! ------------
    bool operator == (const meshSegment &rhs) const
    {
        if(nodeIDs.length()!=rhs.nodeIDs.length()) return false;
        int NbNodes = nodeIDs.length();
        QList<int> tmp = nodeIDs;
        QList<int> tmp1 = rhs.nodeIDs;
        std::sort(tmp.begin(),tmp.end());
        std::sort(tmp1.begin(),tmp1.end());
        for(int i=0; i<NbNodes; i++) if(tmp[i]!=tmp1[i]) return false;
        return true;

        /*
        if(nodeIDs.length() == 2)
        {
            //! -------------------------
            //! first order mesh segment
            //! [a1,b1] == [a2,b2]
            //! -------------------------
            int a1 = nodeIDs.at(0);
            int b1 = nodeIDs.at(1);
            int a2 = rhs.nodeIDs.at(0);
            int b2 = rhs.nodeIDs.at(1);
            if((a1 == a2 && b1 == b2) || (a1 == b2 && b1 == a2)) return true;
            return false;
        }
        else if(nodeIDs.length() == 3)
        {
            //! --------------------------
            //! second order mesh segment
            //! [a1,b1,c1] == [a2,b2,c2]
            //! --------------------------
            int a1 = nodeIDs.at(0); int b1 = nodeIDs.at(1); int c1 = nodeIDs.at(2);
            int a2 = rhs.nodeIDs.at(0); int b2 = rhs.nodeIDs.at(1); int c2 = rhs.nodeIDs.at(2);
            if((a1 == a2 && b1 == b2 && c1 == c2) || (a1 == c2 && b1 == b2 && c1 == a2))
                return true;
            return false;
        }
        */
    }

    //! -----------
    //! operator <
    //! -----------
    bool operator <(const meshSegment &rhs) const
    {
        size_t seed1 = 0, seed2 = 0;
        for(int i=0; i<nodeIDs.size(); i++) hash_c<int>(seed1,nodeIDs[i]);
        for(int i=0; i<rhs.nodeIDs.size(); i++) hash_c<int>(seed2,rhs.nodeIDs[i]);
        if(seed1<seed2) return true;
        return false;
    }
};

//! -------------
//! mesh element
//! -------------
struct meshElement
{
    ElemType type;
    int ID;
    std::vector<int> theNodeIDs;

    //! --------------------
    //! default constructor
    //! --------------------
    meshElement(){;}

    //! ------------
    //! constructor
    //! ------------
    meshElement(ElemType anElemType, int anID, const std::vector<int> &vecNodeIDs)
    {
        type = anElemType;
        ID = anID;
        for(std::vector<int>::const_iterator it = vecNodeIDs.cbegin(); it!=vecNodeIDs.cend(); ++it)
            theNodeIDs.push_back(*it);
    }

    //! -----------------
    //! copy constructor
    //! -----------------
    meshElement(const meshElement& rhs)
    {
        type = rhs.type;
        ID = rhs.ID;
        for(std::vector<int>::const_iterator it = rhs.theNodeIDs.cbegin(); it!=rhs.theNodeIDs.cend(); ++it)
            theNodeIDs.push_back(*it);
    }

    //! -----------
    //! operator =
    //! -----------
    meshElement operator = (const meshElement &rhs)
    {
        type = rhs.type;
        ID = rhs.ID;
        for(std::vector<int>::const_iterator it = rhs.theNodeIDs.cbegin(); it!=rhs.theNodeIDs.cend(); ++it)
            theNodeIDs.push_back(*it);
        return *this;
    }

    //! -----------------------------
    //! operator <
    //! QMap::contains, QMap::insert
    //! -----------------------------
    bool operator < (const meshElement &other) const
    {
        std::vector<int> vec(theNodeIDs.begin(), theNodeIDs.end());
        std::vector<int> vec1(other.theNodeIDs.begin(), other.theNodeIDs.end());
        std::sort(vec.begin(), vec.end());
        std::sort(vec1.begin(), vec1.end());

        size_t seed,seed1;
        seed = seed1 = 0;
        int NbNodes = int(theNodeIDs.size());
        int NbNodes1 = int(other.theNodeIDs.size());

        for(int i=0;i<NbNodes;i++) hash_c<int>(seed,vec.at(i));
        for(int i=0;i<NbNodes1;i++) hash_c<int>(seed1,vec1.at(i));

        if(seed<seed1) return true;
        return false;
    }

    //! ---------------------------------------
    //! operator ==
    //! QList::contains, QList::insert, ...
    //! this checks identity in a cyclic sense
    //! ---------------------------------------
    bool operator == (const meshElement &other)
    {        
        int NbNodes = int(this->theNodeIDs.size());

        //! -----------------
        //! check the length
        //! -----------------
        if(NbNodes!=other.theNodeIDs.size()) return false;

        bool found = false;
        int firstIndex, val1, val2;
        for(firstIndex = 0; firstIndex<NbNodes; firstIndex++)
        {
            val1 = theNodeIDs[firstIndex];
            //val2 = other.theNodeIDs.at(0);
            val2 = other.theNodeIDs[0];
            if(val1==val2)
            {
                found = true;
                break;
            }
        }
        if(!found) return false;
        for(int i=1; i<NbNodes; i++)
        {
            //val1 = theNodeIDs.at((i+firstIndex)%NbNodes);
            //val2 = other.theNodeIDs.at(i);
            val1 = theNodeIDs[(i+firstIndex)%NbNodes];
            val2 = other.theNodeIDs[i];
            if(val1!=val2) return false;
        }
        return true;
    }
};

//! ------------------------------------
//! element normal
//! the normal vector of a face element
//! ------------------------------------
struct elementNormal
{
    double nx,ny,nz;
    double tolerance;

    //elementNormal(){;}
    elementNormal(double a=0.0, double b=0.0, double c=0.0, double aTolerance = 0.0174):
        nx(a),ny(b),nz(c),tolerance(aTolerance){;}

    //! -------------------
    //! operator ==
    //! with tolerance eps
    //! -------------------
    bool operator == (const elementNormal &rhs) const
    {
        //const double eps = 1e-6;
        double dot = nx*rhs.nx+ny*rhs.ny+nz*rhs.nz;
        double l1 = sqrt(nx*nx+ny*ny+nz*nz);
        double l2 = sqrt(rhs.nx*rhs.nx+rhs.ny*rhs.ny+rhs.nz*rhs.nz);
        double dot_ = dot/(l1*l2);
        if(dot_<-1) dot_ = -1;
        if(dot>1) dot_ = 1;
        double angle = acos(dot);
        double eps = tolerance<=rhs.tolerance? tolerance:rhs.tolerance;
        if(fabs(angle)<eps) return true;
        return false;
    }

    //! ------------------------------
    //! operator +
    //! details: n1 + n2 => (n1 + n2)
    //! ------------------------------
    inline elementNormal operator + (const elementNormal &normal1)
    {
        nx += normal1.nx;
        ny += normal1.ny;
        nz += normal1.nz;
        double L1 = sqrt(nx*nx+ny*ny+nz*nz);
        nx = nx/L1;
        ny = ny/L1;
        nz = nz/L1;
        return *this;
    }

    //! ---------------
    //! times a factor
    //! ---------------
    inline elementNormal operator *(double f)
    {
        nx = nx*f;
        ny = ny*f;
        nz = nz*f;
        return *this;
    }

    //! -----------
    //! operator =
    //! -----------
    inline elementNormal operator = (const elementNormal &rhs)
    {
        nx = rhs.nx;
        ny = rhs.ny;
        nz = rhs.nz;
        return *this;
    }

    //! -----------
    //! operator <
    //! -----------
    bool operator <(const mesh::elementNormal &other) const
    {
        size_t seed = 0;
        size_t seed1 = 0;
        hash_c<double>(seed,nx);
        hash_c<double>(seed,ny);
        hash_c<double>(seed,nz);
        hash_c<double>(seed1,other.nx);
        hash_c<double>(seed1,other.ny);
        hash_c<double>(seed1,other.nz);
        if(seed<seed1) return true;
        return false;
    }
};



//! end namespace
}
#endif // MESH_H
