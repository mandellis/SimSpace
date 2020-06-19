#ifndef QEXTENDEDSTANDARDITEM_H
#define QEXTENDEDSTANDARDITEM_H

//! ----------------
//! custom includes
//! ----------------
#include "property.h"
#include "detailviewer.h"
#include "tools.h"
#include <qhistogram.h>

//! ---
//! Qt
//! ---
#include <QStandardItem>
#include <QIcon>
#include <QObject>

class QExtendedStandardItem : public QStandardItem
{
    Q_GADGET

public:

    QExtendedStandardItem(){}
    ~QExtendedStandardItem(){}

    virtual QVariant data(int role) const override;
    virtual void setData(const QVariant &value, int role) override;

    QIcon getIcon(SimulationNodeClass::nodeType theNodeType) const;

private:

    //! getCurrentNode() - helper function
    inline SimulationNodeClass* getCurrentNode() const
    {
        DetailViewer *dt = static_cast<DetailViewer*>(tools::getWidgetByName("detailViewer"));
        return dt->getNode();
    }

public slots:

    void updateIcon(){;}
};

#endif // QEXTENDEDSTANDARDITEM_H
