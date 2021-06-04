#ifndef GEOMETRYTAG_H
#define GEOMETRYTAG_H

//! ----------------
//! custom includes
//! ----------------
#include <src/utils/hash_c.h>

//! ---
//! Qt
//! ---
#include <QMetaType>

//! ----
//! OCC
//! ----
#include <TopAbs_ShapeEnum.hxx>

//! ----
//! C++
//! ----
#include <fstream>

//! ---------------------
//! struct:  geometryTag
//! details:
//! ---------------------
struct GeometryTag
{
    int parentShapeNr;
    int subTopNr;
    bool isParent;
    TopAbs_ShapeEnum subShapeType;

    //! ------------
    //! constructor
    //! ------------
    GeometryTag(int aParentShapeNr = -1, int aSubShapeNr = -1, bool anIsParent = false, TopAbs_ShapeEnum aType = TopAbs_SHAPE):
        parentShapeNr(aParentShapeNr),
        subTopNr(aSubShapeNr),
        isParent(anIsParent),
        subShapeType(aType)
    {
        ;
    }

    //! -----------------
    //! copy constructor
    //! -----------------
    GeometryTag(const GeometryTag &rhs)
    {
        parentShapeNr = rhs.parentShapeNr;
        subTopNr = rhs.subTopNr;
        isParent = rhs.isParent;
        subShapeType = rhs.subShapeType;
    }

    //! -----------
    //! operator <
    //! -----------
    inline bool operator <(const GeometryTag &rhs) const
    {
        size_t seed1 = 0;
        size_t seed2 = 0;
        hash_c<int>(seed1,parentShapeNr);
        hash_c<int>(seed1,subTopNr);
        hash_c<int>(seed1,static_cast<int>(isParent));
        hash_c<int>(seed1,static_cast<int>(subShapeType));

        hash_c<int>(seed2,rhs.parentShapeNr);
        hash_c<int>(seed2,rhs.subTopNr);
        hash_c<int>(seed2,static_cast<int>(rhs.isParent));
        hash_c<int>(seed2,static_cast<int>(rhs.subShapeType));

        if(seed1<seed2) return true; return false;
    }

    //! -----------
    //! operator =
    //! -----------
    inline GeometryTag operator=(const GeometryTag &other)
    {
        parentShapeNr = other.parentShapeNr;
        subTopNr = other.subTopNr;
        isParent = other.isParent;
        subShapeType = other.subShapeType;
        return *this;
    }

    //! ------------
    //! operator ==
    //! ------------
    inline bool operator ==(const GeometryTag &rhs) const
    {
        if(parentShapeNr == rhs.parentShapeNr && subTopNr == rhs.subTopNr
                && isParent == rhs.isParent && subShapeType == rhs.subShapeType)
            return true;
        return false;
    }

    //! ------
    //! write
    //! ------
    void write(std::ofstream &os) const
    {
        if(os.is_open())
        {
            os<<parentShapeNr<<std::endl;
            os<<subTopNr<<std::endl;
            os<<isParent<<std::endl;
            os<<subShapeType<<std::endl;
        }
    }

    //! -----
    //! read
    //! -----
    GeometryTag read(std::ifstream &is)
    {
        if(is.is_open())
        {
            is>>parentShapeNr;
            is>>parentShapeNr;
            is>>subTopNr;
            is>>isParent;
            int type;
            is>>type;
            subShapeType = static_cast<TopAbs_ShapeEnum>(type);
        }
        return *this;
    }

    //! ----------------
    //! contains parent
    //! ----------------
    bool containsParent(int aParentShapeNumber)
    {
        if(parentShapeNr == aParentShapeNumber) return true;
        return false;
    }
};
Q_DECLARE_METATYPE(GeometryTag)

#endif // GEOMETRYTAG_H
