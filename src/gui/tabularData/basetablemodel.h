#ifndef BASETABLEMODEL_H
#define BASETABLEMODEL_H

//! ---
//! Qt
//! ---
#include <QAbstractTableModel>
#include <QObject>

//! --------
//! STL C++
//! --------
#include <vector>
#include <string>
Q_DECLARE_METATYPE(std::vector<int>)
Q_DECLARE_METATYPE(std::vector<std::vector<int>>)
Q_DECLARE_METATYPE(std::vector<bool>)
Q_DECLARE_METATYPE(std::vector<std::vector<bool>>)
Q_DECLARE_METATYPE(std::vector<float>)
Q_DECLARE_METATYPE(std::vector<std::vector<float>>)
Q_DECLARE_METATYPE(std::vector<double>)
Q_DECLARE_METATYPE(std::vector<std::vector<double>>)
Q_DECLARE_METATYPE(std::vector<std::string>)
Q_DECLARE_METATYPE(std::vector<std::vector<std::string>>)
Q_DECLARE_METATYPE(std::string)

//! ----
//! OCC
//! ----
#include <TopoDS_Shape.hxx>
Q_DECLARE_METATYPE(TopoDS_Shape)

//! -------------
//! STL wrappers
//! -------------
#include "stlwrappers.h"

template <typename T>
class Matrix
{
    SSVector<SSVector<T>> inner_;
    unsigned int dimx_, dimy_;

public:

    //! constructor
    Matrix (unsigned int dimx, unsigned int dimy):
        dimx_(dimx),dimy_(dimy)
    {
        inner_.resize(dimx_*dimy_);
        for(unsigned int row=0; row<dimx_; row++)
            for(unsigned int col=0; col<dimy_; col++)
                inner_.push_back(T());
    }

    //! copy constructor
    Matrix (const Matrix &other)
    {
        inner_.resize(other.dimx_*other.dimy_);
        for(unsigned int row=0; row<other.dimx_; row++)
            for(unsigned int col=0; col<other.dimy_; col++)
                inner_.push_back(other.inner_[dimy_*row + col]);
    }

    T& operator()(unsigned int row, unsigned int col)
    {
        if (row >= dimx_ || col>= dimy_) throw std::out_of_range("matrix indices out of range");
        return inner_[dimy_*row + col];
    }

    T operator()(unsigned int row, unsigned int col) const
    {
        if (row >= dimx_ || col>= dimy_) throw std::out_of_range("matrix indices out of range");
        return inner_[dimy_*row + col];
    }

    size_t size() const { return dimx_*dimy_; }
    void clear() { inner_.clear(); }

    //! insert
    void insert(unsigned int row, unsigned int col, const T &value)
    {
        if (row >= dimx_ || col>= dimy_) throw std::out_of_range("matrix indices out of range");
        inner_[row*dimy_+col] = value;
    }

    //! number of rows and columns
    int NbRows() const { return dimx_; }
    int NbCols() const { return dimy_; }

    //! print table - diagnostic
    void print() const
    {
        cout<<endl;
        for(unsigned int row=0; row<dimx_; row++)
        {
            unsigned int col = 0;
            for(; col<dimy_-1; col++) cout<<inner_[row*dimy_+col]<<"\t";
            cout<<inner_[row*dimy_+dimx_-1]<<endl;
        }
        cout<<endl;
    }
};

//! --------------------------------------------------------
//! class BaseTableModel
//! a basic table in which inner data are organized per row
//! --------------------------------------------------------
class BaseTableModel: public QAbstractTableModel
{
    Q_OBJECT

public:

    //! ----------------------
    //! function: constructor
    //! details:  default
    //! ----------------------
    BaseTableModel(QObject* parent=nullptr):
        QAbstractTableModel(parent)
    {
        m_data.clear();
        m_rowCount = 0;
        m_columnCount = 0;
    }

    //! ----------------------
    //! function: constructor
    //! ----------------------
    BaseTableModel(const SSVector<SSVector<QVariant>>& table, bool byRows = true, QObject* parent=nullptr):
        QAbstractTableModel(parent)
    {
        if(byRows)
        {
            try { m_columnCount = table.at(0).size(); }
            catch(...) { return; }
            m_rowCount = table.size();
            for(int row = 0; row<m_rowCount; row++)
            {
                const SSVector<QVariant> &aRow = table.at(row);
                m_data.push_back(aRow);
            }
        }
        else
        {
            // to do
        }
    }

    //! ----------------------
    //! function: constructor
    //! ----------------------
    BaseTableModel(const std::vector<std::vector<QVariant>>& table, bool byRows = true, QObject* parent=nullptr):
        QAbstractTableModel(parent)
    {
        if(byRows)
        {
            try{ m_columnCount = table.at(0).size(); }
            catch(...) { return; }
            m_rowCount = table.size();
            for(int row = 0; row<m_rowCount; row++)
            {
                const std::vector<QVariant> &aRow = table.at(row);
                m_data.push_back(aRow);
            }
        }
        else
        {
            // to do
        }
    }

public:

    //! ---------
    //! rowCount
    //! ---------
    virtual int columnCount(const QModelIndex &parent = QModelIndex()) const override
    {
        Q_UNUSED(parent);
        try { return (int)m_data.at(0).size(); }
        catch(...) { return 0; }
    }
    //! ------------
    //! columnCount
    //! ------------
    virtual int rowCount(const QModelIndex &parent = QModelIndex()) const override
    {
        Q_UNUSED(parent);
        return int(m_data.size());
    }
    //! -----
    //! data
    //! -----
    virtual QVariant data(const QModelIndex &index, int role = Qt::DisplayRole) const override
    {
        if(!index.isValid()) return QVariant();
        int row = index.row();
        int col = index.column();
        QVariant data;
        switch(role)
        {
        case Qt::DisplayRole:
        {
            data = m_data.at(row).at(col);
            if(data.type() == QMetaType::Int) { data.setValue(QString("%1").arg(data.toInt())); }
            if(data.type() == QMetaType::Float) { data.setValue(QString("%1").arg(data.toFloat())); }
            if(data.type() == QMetaType::Double) { data.setValue(QString("%1").arg(data.toDouble())); }
            if(data.type() == QMetaType::QString) { data.setValue(data.toString()); }
            if(data.type() == QMetaType::type("std::string")) { data.setValue(QString::fromStdString(data.value<std::string>())); }
        }
            break;
        case Qt::EditRole: return m_data.at(index.row()).at(index.column()); break;
        case Qt::TextAlignmentRole: return Qt::AlignCenter; break;
        }
        return data;
    }
    //! ------------------------------------------------------------------------------------
    //! headerData
    //! "For horizontal headers, the section number corresponds to the column number.
    //! Similarly, for vertical headers, the section number corresponds to the row number."
    //! ------------------------------------------------------------------------------------
    virtual QVariant headerData(int section, Qt::Orientation orientation, int role = Qt::DisplayRole) const override
    {
        if(role!=Qt::DisplayRole) return QVariant();
        QVariant data;
        switch(orientation)
        {
        case Qt::Horizontal:
        {
            try { data.setValue(m_vlabels.at(section)); }
            catch(...) { data.setValue(QString("Column %1").arg(section)); }    // a default label
        }
            break;
        case Qt::Vertical:
        {
            try { data.setValue(m_hlabels.at(section)); }
            catch(...) { data.setValue(QString("Row %1").arg(section)); }       // a default label
        }
            break;
        }
        return data;
    }
    //! --------
    //! setData
    //! --------
    virtual bool setData(const QModelIndex &index, const QVariant &value, int role = Qt::EditRole) override
    {
        if(!index.isValid()) return false;
        int row = index.row();
        int col = index.column();

        switch (role)
        {
        case Qt::EditRole:
        {
            QVector<int> roles { role };
            m_data[row].replace(col,value);
            emit dataChanged(index, index, roles);
            //cout<<"____tag00____"<<endl;
            return true;
        }
            break;
        }
        return false;
    }
    //! ------
    //! flags
    //! ------
    virtual Qt::ItemFlags flags(const QModelIndex &index) const override
    {
        return QAbstractItemModel::flags(index) | Qt::ItemIsEditable;
    }
    //! ----------
    //! insertRow
    //! ----------
    virtual bool insertRow(int row, const QModelIndex &parent = QModelIndex())
    {
        Q_UNUSED(parent);
        if(row>=this->rowCount()) return false;

        emit beginInsertRows(QModelIndex(),row,row);
        SSVector<QVariant> addRow = m_data.at(row);
        m_data.insert(row,addRow);
        //! update the number of rows
        m_rowCount = m_data.size();
        emit endInsertRows();

        //! inform the viewer about the position of the new data
        QVector<int> roles { Qt::EditRole };
        QModelIndex topLeftIndex = this->createIndex(row,0);
        QModelIndex bottomRightIndex = this->createIndex(this->rowCount()-1,this->columnCount()-1);
        emit dataChanged(topLeftIndex,bottomRightIndex,roles);
        return true;
    }
    //! -----------
    //! insertRows
    //! -----------
    virtual bool insertRows(int row, int count, const QModelIndex &parent = QModelIndex()) override
    {
        Q_UNUSED(parent);

        emit beginInsertRows(QModelIndex(),row,row+count-1);
        SSVector<QVariant> vec = m_data.at(row);
        for(int c=row; c<row+count; c++) m_data.insert(c,vec);

        //! update the number of rows
        m_rowCount = m_data.size();
        m_columnCount = m_data.at(0).size();
        emit endInsertRows();

        //! inform the viewer about the position of the new data
        QVector<int> roles { Qt::EditRole };
        QModelIndex topLeftIndex = this->createIndex(row,0);
        QModelIndex bottomRightIndex = this->createIndex(this->rowCount()-1,this->columnCount()-1);
        emit dataChanged(topLeftIndex,bottomRightIndex,roles);
        return true;
    }
    //! -----------
    //! removeRows
    //! -----------
    virtual bool removeRows(int first, int last, const QModelIndex &parent = QModelIndex())
    {
        Q_UNUSED(parent)

        emit beginRemoveRows(QModelIndex(),first,last);
        for(int k=last; k>=first; k--) m_data.remove(k);

        //! update the number of rows
        m_rowCount = m_data.size();
        emit endRemoveRows();

        //! inform the viewer about the position of the new data
        QVector<int> roles { Qt::EditRole };
        QModelIndex topLeftIndex = this->createIndex(first,0);
        QModelIndex bottomRightIndex = this->createIndex(this->rowCount()-1,this->columnCount()-1);
        emit dataChanged(topLeftIndex,bottomRightIndex,roles);
        return true;
    }
    //! ----------
    //! removeRow
    //! ----------
    virtual bool removeRow(int row, const QModelIndex &parent = QModelIndex())
    {
        Q_UNUSED(parent)
        if(row<0) return false;
        if(row>=this->rowCount()) return false;

        emit beginRemoveRows(QModelIndex(),row,row);
        m_data.remove(row);

        //! update the number of rows
        m_rowCount = m_data.size();
        emit endRemoveRows();

        //! inform the viewer about the position of the new data
        QVector<int> roles { Qt::EditRole };
        QModelIndex topLeftIndex = this->createIndex(row,0);
        QModelIndex bottomRightIndex = this->createIndex(this->rowCount()-1,this->columnCount()-1);
        emit dataChanged(topLeftIndex,bottomRightIndex,roles);
        return true;
    }
    //! --------------
    //! insertColumns
    //! --------------
    virtual bool insertColumns(int column, int count, const QModelIndex &parent = QModelIndex()) override
    {
        Q_UNUSED(parent)

        int Nrows = this->rowCount();   // number of elements
        if(Nrows==0) return false;      // zero rows

        int first = column;
        int last = column+count-1;
        this->beginInsertColumns(QModelIndex(),first,last);
        for(int row=0; row<Nrows; row++)
        {
            SSVector<QVariant> vec = m_data.at(row);
            for(int col = first; col<count; col++)
            {
                QVariant val = m_data.at(row).at(col);
                vec.insert(col,val);
                m_data.insert(row,vec);
            }
        }
        //! update the number of columns
        m_columnCount = m_data.at(0).size();
        this->endInsertColumns();

        //! inform the viewer about the position of the new data
        QVector<int> roles { Qt::EditRole };
        QModelIndex topLeftIndex = this->createIndex(0,first);
        QModelIndex bottomRightIndex = this->createIndex(this->rowCount()-1,last);
        emit dataChanged(topLeftIndex,bottomRightIndex,roles);
        return true;
    }
    //! -------------
    //! insertColumn
    //! -------------
    virtual bool insertColumn(int column, const QModelIndex &parent = QModelIndex())
    {
        Q_UNUSED(parent)
        if(this->columnCount()==0) return false;
        if(column<0 || column>this->columnCount()-1) return false;
        int Nrows = this->rowCount();

        emit beginInsertColumns(QModelIndex(),column,column);
        for(int row=0; row<Nrows; row++)
        {
            QVariant data = m_data.at(row).at(column);
            m_data.at(row).insert(column,data);
        }
        //! update the number of columns
        m_columnCount = m_data.at(0).size();
        emit endInsertColumns();

        //! inform the viewer about the position of the new data
        QVector<int> roles { Qt::EditRole };
        QModelIndex topLeftIndex = this->createIndex(0,column);
        QModelIndex bottomRightIndex = this->createIndex(this->rowCount()-1,this->columnCount()-1);
        emit dataChanged(topLeftIndex,bottomRightIndex,roles);
        return true;
    }
    //! -------------
    //! removeColumn
    //! -------------
    virtual bool removeColumn(int column, const QModelIndex &parent = QModelIndex())
    {
        Q_UNUSED(parent)

        if(column<0 || column>=this->columnCount()-1) return false;
        int first = column;
        emit beginRemoveColumns(QModelIndex(),column,column);
        for(int row = 0; row<this->rowCount(); row++)
        {
            SSVector<QVariant> vec = m_data.at(row);
            vec.remove(column);
            m_data.replace(row,vec);
        }
        //! update the number of columns
        m_columnCount = m_data.at(0).size();            // check error .at() to do ...
        emit endRemoveColumns();

        //! inform the viewer about the position of the new data
        QVector<int> roles { Qt::EditRole };
        QModelIndex topLeftIndex = this->createIndex(0,first);
        QModelIndex bottomRightIndex = this->createIndex(this->rowCount()-1,this->columnCount()-1);
        emit dataChanged(topLeftIndex,bottomRightIndex,roles);
        return true;
    }
    //! --------------
    //! removeColumns
    //! --------------
    virtual bool removeColumns(int first, int count, const QModelIndex &parent = QModelIndex())
    {
        Q_UNUSED(parent)

        emit beginRemoveColumns(QModelIndex(),first,first+count-1);
        for(int i=0; i<this->rowCount(); i++)
        {
            SSVector<QVariant> dataVec = m_data.at(i);
            for(int k=first+count-1; k>=first; k--) dataVec.remove(k);
            m_data.replace(i,dataVec);
        }

        //! update the number of columns
        m_columnCount = m_data.at(0).size();
        emit endRemoveColumns();

        //! inform the viewer about the position of the new data
        QVector<int> roles { Qt::EditRole };
        QModelIndex topLeftIndex = this->createIndex(0,first);
        QModelIndex bottomRightIndex = this->createIndex(this->rowCount()-1,this->columnCount()-1);
        emit dataChanged(topLeftIndex,bottomRightIndex,roles);
        return true;
    }
    //! --------------
    //! setHeaderData
    //! --------------
    virtual bool setHeaderData(int section, Qt::Orientation orientation, const QVariant &value, int role) override
    {
        if(role != Qt::DisplayRole) return false;     // for the moment only for view
        switch(orientation)
        {
        case Qt::Horizontal:
        {
            std::map<int,QString>::iterator it = m_vlabels.find(section);
            if(it!=m_vlabels.end()) it->second = value.toString();
            else m_vlabels.insert(std::make_pair(section,value.toString()));
        }
            break;
        case Qt::Vertical:
        {
            std::map<int,QString>::iterator it = m_hlabels.find(section);
            if(it!=m_hlabels.end()) it->second = value.toString();
            else m_hlabels.insert(std::make_pair(section,value.toString()));
        }
            break;
        }
        emit headerDataChanged(orientation,section,section);
        return true;
    }
    //! --------
    //! helpers
    //! --------
    virtual QVariant dataByRC(int row, int col, int role = Qt::EditRole) const
    {
        if(row>=this->rowCount() || col>=this->columnCount()) return QVariant();
        QModelIndex index = this->createIndex(row,col);
        return this->data(index,role);
    }
    virtual bool setDataByRC(int row, int col, QVariant data, int role = Qt::EditRole)
    {
        if(row>=this->rowCount() || col>=this->columnCount()) return false;
        QModelIndex index = this->createIndex(row,col);
        this->setData(index,data,role);     // it already emits dataChanged()
        return true;
    }
    virtual bool getColumn(int col, QVariant* column)
    {
        if(col>=this->columnCount()) return false;
        for(int row=0; row<this->rowCount(); row++) column[row] = m_data[row][col];
        return true;
    }
    virtual bool getRow(int row, QVariant* column)
    {
        if(row>=this->rowCount()) return false;
        for(int col=0; col<this->columnCount(); col++) column[col] = m_data[row][col];
        return true;
    }
    virtual bool setColumnHeader(int col, const std::string &label)
    {
        if(col<0 || col>this->columnCount()-1) return false;
        QVariant value;
        value.setValue(QString::fromStdString(label));
        return this->setHeaderData(col,Qt::Horizontal,value,Qt::DisplayRole);
    }
    virtual bool setColumnHeader(int col, const QString &label)
    {
        return this->setColumnHeader(col,label.toStdString());
    }
    virtual bool setRowHeader(int row, const std::string &label)
    {
        if(row<0 || row>this->rowCount()-1) return false;
        QVariant value;
        value.setValue(QString::fromStdString(label));
        return this->setHeaderData(row,Qt::Vertical,value,Qt::DisplayRole);
    }
    virtual bool setRowHeader(int row, const QString &label)
    {
        return this->setRowHeader(row,label.toStdString());
    }
    //! -------------------------------------
    //! overloads of insertColumn, insertRow
    //! -------------------------------------
    virtual bool insertColumn(int pos, std::vector<QVariant> &vecData)
    {
        if(vecData.empty()) return false;
        if(pos<0 || pos>this->columnCount()-1) return false;
        emit beginInsertColumns(QModelIndex(),pos,pos);
        for(int row=0; row<this->rowCount(); row++)
        {
            SSVector<QVariant> aRow = m_data.at(row);
            aRow.insert(pos,vecData.at(row));
            m_data.replace(row,aRow);
        }
        m_columnCount = m_data.at(0).size();
        emit endInsertColumns();

        //! inform the viewer about the position of the new data
        QVector<int> roles { Qt::EditRole };
        QModelIndex topLeftIndex = this->createIndex(0,pos);
        QModelIndex bottomRightIndex = this->createIndex(this->rowCount()-1,this->columnCount()-1);
        emit dataChanged(topLeftIndex,bottomRightIndex,roles);
        return true;
    }
    virtual bool insertRow(int pos, std::vector<QVariant> &vecData)
    {
        if(vecData.empty()) return false;
        if(pos<0 || pos>this->rowCount()-1) return false;
        emit beginInsertRows(QModelIndex(),pos,pos);
        SSVector<QVariant> aRow;
        for(int i=0; i<vecData.size(); i++) aRow.push_back(vecData[i]);
        m_data.insert(pos,aRow);
        m_rowCount = m_data.size();
        emit endInsertRows();

        //! inform the viewer about the position of the new data
        QVector<int> roles { Qt::EditRole };
        QModelIndex topLeftIndex = this->createIndex(pos,0);
        QModelIndex bottomRightIndex = this->createIndex(this->rowCount()-1,this->columnCount()-1);
        emit dataChanged(topLeftIndex,bottomRightIndex,roles);
        return true;
    }
    virtual bool appendRow(std::vector<QVariant> &vecData)
    {
        if(vecData.empty()) return false;
        int pos = this->rowCount();                     // use new row index
        emit beginInsertRows(QModelIndex(),pos,pos);
        SSVector<QVariant> aRow;
        for(int i=0; i<vecData.size(); i++) aRow.push_back(vecData[i]);
        m_data.push_back(aRow);
        m_rowCount = m_data.size();
        emit endInsertRows();

        //! inform the viewer about the position of the new data
        QVector<int> roles { Qt::EditRole };
        QModelIndex topLeftIndex = this->createIndex(pos,0);
        QModelIndex bottomRightIndex = this->createIndex(this->rowCount()-1,this->columnCount()-1);
        emit dataChanged(topLeftIndex,bottomRightIndex,roles);
        return true;
    }
    virtual bool appendColumn(std::vector<QVariant> &vecData)
    {
        if(vecData.empty()) return false;
        int pos = this->columnCount();                     // use new column index
        emit beginInsertColumns(QModelIndex(),pos,pos);
        for(int row=0; row<this->rowCount(); row++)
        {
            //cout<<"____adding element at row: "<<row;
            //cout<<" value: "<<vecData.at(row).toDouble()<<"____"<<endl;
            SSVector<QVariant> aRow = m_data.at(row);
            aRow.push_back(vecData.at(row));
            m_data.replace(row,aRow);
        }
        m_columnCount = m_data.at(0).size();
        emit endInsertColumns();

        //! inform the viewer about the position of the new data
        QVector<int> roles { Qt::EditRole };
        QModelIndex topLeftIndex = this->createIndex(0,pos);
        QModelIndex bottomRightIndex = this->createIndex(this->rowCount()-1,this->columnCount()-1);
        emit dataChanged(topLeftIndex,bottomRightIndex,roles);
        return true;
    }

protected:  // class members

    SSVector<SSVector<QVariant>> m_data;            // inner SSVector => row
    size_t m_columnCount;                           // number of columns
    size_t m_rowCount;                              // number of rows
    std::map<int,QString> m_hlabels;                // horizontal labels for header
    std::map<int,QString> m_vlabels;                // vertical labels for header
};

#endif // BASETABLEMODEL_H
