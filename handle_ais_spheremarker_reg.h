#ifndef HANDLE_AIS_SPHEREMARKER_REG_H
#define HANDLE_AIS_SPHEREMARKER_REG_H

//! redefinition of opencascade Handle() macro
#include "occhandle.h"

#include "ais_spheremarker.h"
#include <QMetaType>

typedef occHandle(AIS_SphereMarker) AIS_SphereMarker_handle_reg;
Q_DECLARE_METATYPE(AIS_SphereMarker_handle_reg)

#endif // HANDLE_AIS_SPHEREMARKER_REG_H
