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
    PENTA,
    EPTA,

    //! null element
    NULL_ELEMENT
};
Q_DECLARE_METATYPE(ElemType)

#endif // ELEMENTTYPES
