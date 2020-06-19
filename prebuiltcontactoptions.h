#ifndef PREBUILTCONTACTOPTIONS_H
#define PREBUILTCONTACTOPTIONS_H

//! custom includes
#include "mydefines.h"

//! Qt
#include <QMetaType>
#include <QVector>

struct preBuiltContactOptions
{
    bool isAutomatic;
    QVector<GeometryTag> vecLocMaster;
    QVector<GeometryTag> vecLocSlave;
};

Q_DECLARE_METATYPE(preBuiltContactOptions)

#endif // PREBUILDCONTACTOPTIONS_H
