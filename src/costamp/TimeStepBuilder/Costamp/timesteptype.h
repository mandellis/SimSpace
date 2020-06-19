#ifndef TIMESTEPTYPE_H
#define TIMESTEPTYPE_H

#include <QMetaType>

enum TimeStepType
{
    ClosedAssemblyWithoutPressure,
    ClosedAssemblyWithPressure,
    OpenAssembly,
    ToBeRemoved
};

Q_DECLARE_METATYPE(TimeStepType)

struct TimeStepMarker
{
    double time;
    double value;
    TimeStepType tsType;
    bool confirmed;

    //! ------------
    //! constructor
    //! ------------
    TimeStepMarker(double aTime= -1, double aValue = -1, TimeStepType aType = ToBeRemoved, bool isConfirmed = false)
    {
        time = aTime;
        value = aValue;
        tsType = aType;
        confirmed = isConfirmed;
    }

    //! -----------
    //! isSameTime
    //! -----------
    inline bool isSameTime(const TimeStepMarker &otherTimeStep) const
    {
        if(this->time == otherTimeStep.time) return true;
        return false;
    }

    //! -----------------
    //! copy constructor
    //! -----------------
    TimeStepMarker::TimeStepMarker(const TimeStepMarker &other)
    {
        time = other.time;
        value = other.value;
        tsType = other.tsType;
        confirmed = other.confirmed;
    }

    //! -----------
    //! operator =
    //! -----------
    TimeStepMarker operator = (const TimeStepMarker &other)
    {
        time = other.time;
        value = other.value;
        tsType = other.tsType;
        confirmed = other.confirmed;
        return *this;
    }

    //! ------------------------------
    //! operator ==
    //! defined for QList::contains()
    //! ------------------------------
    inline bool operator ==(const TimeStepMarker &other) const
    {
        if(time==other.time && value==other.value && tsType==other.tsType)
            return true;
        return false;
    }

    //! ---------------------------
    //! operator <
    //! defined for QMap::insert()
    //! ---------------------------
    inline bool operator < (const TimeStepMarker &other) const
    {
        if(time<other.time) return true;
        else return false;
    }
};

Q_DECLARE_METATYPE(TimeStepMarker)

#endif // TIMESTEPTYPE_H
