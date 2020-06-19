#ifndef RUNTERMINATIONDATA_H
#define RUNTERMINATIONDATA_H

#include <QMetaType>

struct runTerminationData
{
    double lastAvailableTime;
    int lastAvailableStep;
    int lastAvailableSubStep;
};

Q_DECLARE_METATYPE(runTerminationData)

#endif // RUNTERMINATIONDATA_H
