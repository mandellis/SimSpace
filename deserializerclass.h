#ifndef DESERIALIZERCLASS_H
#define DESERIALIZERCLASS_H

//! Qt
#include <QObject>
#include <QString>

//! custom includes
#include "simulationnodeclass.h"

class QExtendedStandardItem;

class deserializerClass: public QObject
{
    Q_OBJECT

public:

    deserializerClass(QObject *parent=0);
    QVector<Property> deserialize(const QString &nodeFileName, QString &nodeType, QString &nodeName) const;
    QVector<load> readTabularDataFromFile(const QString &tabularDataFileAbsolutePath) const;

    SimulationNodeClass* nodeBuilder(const QString &nodeFileName);

    QString getNodeName(const QString &nodeFileName);
    QString getNodeType(const QString &nodeFileName);
};

#endif // DESERIALIZERCLASS_H
