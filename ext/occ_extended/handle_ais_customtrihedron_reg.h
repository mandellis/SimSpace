#ifndef HANDLE_AIS_CUSTOMTRIHEDRON_H
#define HANDLE_AIS_CUSTOMTRIHEDRON_H

//! redefinition of opencascade Handle() macro
#include "occhandle.h"

#include "ais_customtrihedron.h"
#include <QMetaType>

typedef occHandle(AIS_CustomTrihedron) AIS_CustomTrihedron_handle_reg;
Q_DECLARE_METATYPE(AIS_CustomTrihedron_handle_reg)

#endif // HANDLE_AIS_CUSTOMTRIHEDRON_H
