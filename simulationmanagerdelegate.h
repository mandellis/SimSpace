#ifndef SIMULATIONMANAGERDELEGATE_H
#define SIMULATIONMANAGERDELEGATE_H

#include <QStyledItemDelegate>
#include <QStyleOptionViewItem>
#include <QModelIndex>

class SimulationManagerDelegate: public QStyledItemDelegate
{
public:

    SimulationManagerDelegate(QWidget *parent=0);

    void paint(QPainter *painter,
               const QStyleOptionViewItem &option,
               const QModelIndex &index) const;
};

#endif // SIMULATIONMANAGERDELEGATE_H
