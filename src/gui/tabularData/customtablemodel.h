#ifndef CUSTOMTABLEMODEL_H
#define CUSTOMTABLEMODEL_H

//! ---
//! Qt
//! ---
#include <QAbstractTableModel>
#include <QHash>
#include <QRect>
#include <QVector>

//! ----------------
//! custom includes
//! ----------------
#include "load.h"

class CustomTableModel : public QAbstractTableModel
{
    Q_OBJECT

public:

    explicit CustomTableModel(const QVector<load> &vecLoads, bool addFirstRow=true, QObject *parent = 0);

    int rowCount(const QModelIndex &parent = QModelIndex()) const;
    int columnCount(const QModelIndex &parent = QModelIndex()) const;
    QVariant headerData(int section, Qt::Orientation orientation, int role = Qt::DisplayRole) const;
    QVariant data(const QModelIndex &index, int role = Qt::DisplayRole) const;

    //! get data through (row, column)
    QVariant dataRC(int row, int column, int role = Qt::EditRole) const;

    bool setData(const QModelIndex &index, const QVariant &value, int role = Qt::EditRole);

    //! set data through (row, column)
    bool setDataRC(const QVariant &value, int row, int column, int role=Qt::EditRole);

    Qt::ItemFlags flags(const QModelIndex &index) const;
    bool insertRows(int row, int count, const QModelIndex &parent = QModelIndex());
    bool removeRows(int first, int last, const QModelIndex &parent = QModelIndex());
    bool insertColumns(int column, int count, const QModelIndex &parent = QModelIndex());
    bool removeColumns(int column, int count, const QModelIndex &parent = QModelIndex());

    bool appendColumn(const load &aLoad, const QModelIndex &parent = QModelIndex());
    load getColumn(int col) const;
    void addMapping(QString color, QRect area);

    void clearMapping()
    {
        m_mapping.clear();
    }

    int getColumnBeforeBC() {return columBeforeBC;}

private:

    QList<QVector<QVariant>*> m_data;
    QHash<QString, QRect> m_mapping;
    int m_columnCount;
    int m_rowCount;
    QList<Property::loadType> m_loadTypes;
    //QString getHeaderString(int section) const;
    load loadToAppend;

    int columBeforeBC;

public:

    double maxVal() const;
    double minVal() const;
    double minXVal() const;
    double maxXVal() const;
    bool appendRow(int count, const QModelIndex &parent/* = QModelIndex()*/);
    void setLoadToInsert(load aload);
    QModelIndex makeIndex(int row, int column) {return this->createIndex(row,column); }
    QString getHeaderString(int section) const;

public slots:

    void autoUpdateTable(QModelIndex topLeftIndex, QModelIndex);

signals:

    void boltStatusDefineByChangedThroughTabularView();
};

#endif // CUSTOMTABLEMODEL_H
