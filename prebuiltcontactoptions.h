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
    std::vector<GeometryTag> vecLocMaster;
    std::vector<GeometryTag> vecLocSlave;
};

Q_DECLARE_METATYPE(preBuiltContactOptions)

#endif // PREBUILDCONTACTOPTIONS_H
