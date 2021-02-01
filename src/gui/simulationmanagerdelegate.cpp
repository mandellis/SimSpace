#include "simulationmanagerdelegate.h"

#include <QPainter>
#include <QFontMetrics>
#include <QRect>
#include <iostream>
using namespace std;

SimulationManagerDelegate::SimulationManagerDelegate(QWidget *parent):QStyledItemDelegate(parent)
{
    ;
}
#include <QLineEdit>
#include <QLabel>
#include "qextendedstandarditem.h"
void SimulationManagerDelegate::paint(QPainter *painter,
                                      const QStyleOptionViewItem &option,
                                      const QModelIndex &index) const
{
    QStyleOptionViewItem opt = option;

    if(opt.state & QStyle::State_HasFocus)
    {
        painter->save();
        QStyledItemDelegate::paint(painter,option,index);
        opt.showDecorationSelected = false;
        painter->restore();
    }
    else
    {
        QStyledItemDelegate::paint(painter,option,index);
    }
}
