#ifndef SERIALIZERCLASS_H
#define SERIALIZERCLASS_H

#include <QObject>
#include "listofshape.h"

class QExtendedStandardItem;

class serializerClass: public QObject
{
    Q_OBJECT

private:

    QExtendedStandardItem* myItem;
    QString myFileFullName;
    QString myPath;

public:

    serializerClass(QObject* parent=0);
    void setItem(QExtendedStandardItem *anItem);
    void setSavingDirPath(const QString &path);
    void serialize(int n=-1);       //! if n=-1 => automatic numbering
};

#endif // SERIALIZERCLASS_H
