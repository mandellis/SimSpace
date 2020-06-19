#ifndef HANDLE_AIS_COLOREDSHAPE_REG_H
#define HANDLE_AIS_COLOREDSHAPE_REG_H

#include <QMetaType>
#include <AIS_ColoredShape.hxx>

typedef occHandle(AIS_ColoredShape) AIS_ColoredShape_handle_reg;
Q_DECLARE_METATYPE(AIS_ColoredShape_handle_reg)

#endif // HANDLE_AIS_COLOREDSHAPE_REG_H
