#ifndef LOAD_H
#define LOAD_H

//! ---
//! Qt
//! ---
#include <QMetaType>
#include <QVector>
#include <QVariant>

//! ----------------
//! custom includes
//! ----------------
#include "myenumvariables.h"
#include "property.h"

class load
{
    Q_GADGET

public:

    //! --------------------
    //! default constructor
    //! --------------------
    load();

    //! -----------------
    //! copy constructor
    //! -----------------
    load(const load &other)
    {
        myLoadType = other.myLoadType;
        for(int i=0; i<other.myValues.size(); i++)
        {
            myValues.push_back(other.myValues.at(i));
        }
    }

    //! constructor I
    load(QVector<QVariant> values, Property::loadType type = Property::loadType_none);

    //! return the number of times (number of values)
    int NbTimes() const;

    //! return the vector of values
    QVector<QVariant> values() const;

    //! return the type
    Property::loadType type() const;

    //! return the type ("loadType_<>") as a QString
    const QString getLoadType() const;

    //! set the load type
    void setLoadType(const QString &theLoadType);

private:

    Property::loadType myLoadType;
    QVector<QVariant> myValues;

public:

    void setData(const QVector<QVariant> &values);
    void setType(Property::loadType theType);
    void write(std::ofstream &out) const;

    //! write load
    //static void writeLoad(const load &value, std::ofstream &out);

    //! read load
    static load readLoad(std::ifstream &in);
};

Q_DECLARE_METATYPE(load)

#endif // LOAD_H
