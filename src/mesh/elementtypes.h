#ifndef ELEMENTTYPES
#define ELEMENTTYPES

#include <vector>
#include <QMetaType>
#include <QList>

enum ElemType
{
    //! 3D elements
    TET,
    TET10,
    PYRAM,
    PYRAM13,
    PRISM,
    PRISM15,
    HEXA,
    HEXA20,

    //! 2D elements
    TRIG,
    TRIG6,
    QUAD,
    QUAD6,
    QUAD8,

    //! null element
    NULL_ELEMENT
};
Q_DECLARE_METATYPE(ElemType)

/*
//! ------------
//! struct node
//! ------------
struct node
{    
    int ID;
    double x,y,z;

    //! operator ==
    inline bool operator ==(const node &rhs)
    {
        return ID==rhs.ID && x==rhs.x && y ==rhs.y && z ==rhs.z;
    }

    //! operator !=
    inline bool operator !=(const node &rhs)
    {
        return ID!=rhs.ID || x!=rhs.x || y !=rhs.y || z !=rhs.z;
    }

    //! default constructor
    node(){;}

    //! constructor
    node(int anID, double aX, double aY, double aZ): ID(anID), x(aX), y(aY), z(aZ){;}

    //! copy constructor
    node(const node &other)
    {
        ID = other.ID;
        x = other.x;
        y = other.y;
        z = other.z;
    }
};

Q_DECLARE_METATYPE(node)
*/
#endif // ELEMENTTYPES
