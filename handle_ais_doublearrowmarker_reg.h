#ifndef HANDLE_AIS_DOUBLEARROWMARKER_H
#define HANDLE_AIS_DOUBLEARROWMARKER_H

//! redefinition of opencascade Handle() macro
#include "occhandle.h"

#include "ais_doublearrowmarker.h"
#include <QMetaType>

typedef occHandle(AIS_DoubleArrowMarker) AIS_DoubleArrowMarker_handle_reg;
Q_DECLARE_METATYPE(AIS_DoubleArrowMarker_handle_reg)

#endif // HANDLE_AIS_DOUBLEARROWMARKER_H
