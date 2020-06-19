#ifndef CELLDATA_H
#define CELLDATA_H

#include <QMetaType>

struct cell
{
    QString name;
    QVariant data;
};

Q_DECLARE_METATYPE(cell)

#endif // CELLDATA_H
