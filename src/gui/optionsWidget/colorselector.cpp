#include "colorselector.h"

//! ----------------------------------------------------------------
//! function: constructor
//! details:
//! ----------------------------------------------------------------
ColorSelector::ColorSelector(QWidget *parent):QComboBox(parent)
{
    this->populateList();
}

//! ----------------------------------------------------------------
//! function: color
//! details:
//! ----------------------------------------------------------------
QColor ColorSelector::color() const
{
    return qvariant_cast<QColor>(itemData(currentIndex(), Qt::DecorationRole));
}

//! ----------------------------------------------------------------
//! function: setColor
//! details:
//! ----------------------------------------------------------------
void ColorSelector::setColor(QColor color)
{
    QVariant a;
    setCurrentIndex(findData(color, int(Qt::DecorationRole)));
}

//! ----------------------------------------------------------------
//! function: populateList
//! details:
//! ----------------------------------------------------------------
void ColorSelector::populateList()
{
    QStringList colorNames = QColor::colorNames();
    for (int i = 0; i < colorNames.size(); ++i)
    {
        QColor color(colorNames[i]);
        insertItem(i, colorNames[i]);
        setItemData(i, color, Qt::DecorationRole);
        setItemData(i, color, Qt::BackgroundRole);
    }
}

Q_DECLARE_OPAQUE_POINTER(ColorSelector*)
