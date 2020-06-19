#ifndef QTABLESTANDARDITEM_H
#define QTABLESTANDARDITEM_H

#include <QStandardItem>

class QTableStandardItem: public QStandardItem
{
public:

    QTableStandardItem(){;}
    ~QTableStandardItem(){}

    virtual QVariant data(int role) const override;
    virtual void setData(const QVariant &value, int role) override;
};

#endif // QTABLESTANDARDITEM_H


