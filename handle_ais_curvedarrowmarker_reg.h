#ifndef HANDLE_AIS_CURVEDARROWMARKER_REG_H
#define HANDLE_AIS_CURVEDARROWMARKER_REG_H

//! redefinition of opencascade Handle() macro
#include "occhandle.h"

#include <QMetaType>
#include <ais_curvedarrowmarker.h>

typedef occHandle(AIS_CurvedArrowMarker) AIS_CurvedArrowMarker_handle_reg;
Q_DECLARE_METATYPE(AIS_CurvedArrowMarker_handle_reg)

#endif // HANDLE_AIS_CURVEDARROWMARKER_REG_H
