#ifndef PARSER_H
#define PARSER_H

#include <geometrydatabase.h>

#include "simulationnodeclass.h"
#include "mydefines.h"

#include <QVector>
#include <QList>
#include <QModelIndex>

class parser
{
public:

    parser();

    //! parse an item
    static bool parseItem(QExtendedStandardItem* item);

    //! parse cylindrical support
    static bool CylindricalSupport(QExtendedStandardItem *item);

    //! parse a fixed support
    static bool FixedSupport(QExtendedStandardItem *item);

    //! parse a frictionless support
    static bool FrictionlessSupport(QExtendedStandardItem *item);

    //! parse a force
    static bool Force(QExtendedStandardItem *item);

    //! parse a force
    static bool Moment(QExtendedStandardItem *item);

    //! parse a pressure
    static bool Pressure(QExtendedStandardItem *item);

    //! thermal condition
    static bool ThermalCondition(QExtendedStandardItem *item);

    //! parse a mesh control
    static bool MeshControl(QExtendedStandardItem *item);

    //! parse a named selection
    static bool NamedSelection(QExtendedStandardItem *item);

    //! parse a contact pair
    static bool ContactPair(QExtendedStandardItem *item);

    //! parse a displacement
    static bool Displacement(QExtendedStandardItem *item);

    //! parse a body load
    static bool BodyLoads(QExtendedStandardItem *item);

    //! parse a post-processing item
    static bool PostProcessingItem(QExtendedStandardItem *item);

    //! -----------------------------------------------------------------
    //! background: no background color => OK; yellow => invalid setting
    //! -----------------------------------------------------------------
    static void setItemBackground(SimulationNodeClass *node, QList<QModelIndex> invalidEntries);
};

#endif // PARSER_H
