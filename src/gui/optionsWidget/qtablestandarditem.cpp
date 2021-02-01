//! custom includes
#include "qtablestandarditem.h"

//! Qt
#include <QColor>

void QTableStandardItem::setData(const QVariant &value, int role)
{
    if(role==Qt::UserRole)
    {
        QStandardItem::setData(value,role);
    }
}

#include "viewoptions.h"
QVariant QTableStandardItem::data(int role) const
{
    if(role==Qt::DisplayRole)
    {
        QVariant data;
        QVariant itemData = QStandardItem::data(Qt::UserRole);
        if(itemData.canConvert<QString>())
        {
            data.setValue(itemData.toString());
        }
        else if(itemData.canConvert<QColor>())
        {
            QString colorName = itemData.value<QColor>().name();
            data.setValue(colorName);
        }
        else if(itemData.canConvert<viewOptions>())
        {
            switch(itemData.value<viewOptions>())
            {
            case typeOfGradient_uniform: data.setValue(QString("Uniform")); break;
            case typeOfGradient_horizontal: data.setValue(QString("Horizontal")); break;
            case typeOfGradient_vertical: data.setValue(QString("Vertical")); break;
            }
        }
        return data;
    }
    else
    {
        return QStandardItem::data(role);
    }
}
