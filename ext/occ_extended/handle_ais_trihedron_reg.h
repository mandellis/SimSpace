#ifndef HANDLE_AIS_TRIHEDRON_REG_H
#define HANDLE_AIS_TRIHEDRON_REG_H

//! redefinition of opencascade Handle() macro
#include "occhandle.h"

#include <QMetaType>
#include <AIS_Trihedron.hxx>

typedef occHandle(AIS_Trihedron) AIS_Trihedron_handle_reg;
Q_DECLARE_METATYPE(AIS_Trihedron_handle_reg)

#endif // HANDLE_AIS_TRIHEDRON_REG_H
