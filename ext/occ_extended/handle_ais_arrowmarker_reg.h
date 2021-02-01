#ifndef HANDLE_AIS_ARROWMARKER_REG_H
#define HANDLE_AIS_ARROWMARKER_REG_H

//! redefinition of opencascade Handle() macro
#include "occhandle.h"

#include <QMetaType>
#include "ais_arrowmarker.h"

typedef occHandle(AIS_ArrowMarker) AIS_ArrowMarker_handle_reg;
Q_DECLARE_METATYPE(AIS_ArrowMarker_handle_reg)

#endif // HANDLE_AIS_ARROWMARKER_REG_H
