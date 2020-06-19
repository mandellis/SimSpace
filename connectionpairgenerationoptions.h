#ifndef CONNECTIONPAIRGENERATIONOPTIONS_H
#define CONNECTIONPAIRGENERATIONOPTIONS_H

#include <indexedmapofmeshdatasources.h>

//! ---
//! Qt
//! ---
#include <QMap>
#include <QString>
#include <QMetaType>

//! ----
//! OCC
//! ----
#include "occhandle.h"
#include <MeshVS_DataSource.hxx>
#include <ng_meshvs_datasourceface.h>

struct connectionPairGenerationOption
{
    bool manual;
    QString timeTag;
    QString connectionPairName;

    //! ----------------------
    //! constructor - default
    //! ----------------------
    connectionPairGenerationOption()
    {
        manual = false;
        timeTag = "";
        connectionPairName ="";
    }

    //! ------------
    //! constructor
    //! ------------
    connectionPairGenerationOption(bool isManual,
                                   const QString aConnectionPair = "",
                                   const QString aTimeTag = ""):
        manual(isManual),
        connectionPairName(aConnectionPair),
        timeTag(aTimeTag)
    {
        ;
    }

    //! -----------------
    //! copy constructor
    //! -----------------
    connectionPairGenerationOption(const connectionPairGenerationOption &other)
    {
        manual = other.manual;
        timeTag = other.timeTag;
        connectionPairName = other.connectionPairName;
    }

    //! -----------
    //! operator =
    //! -----------
    connectionPairGenerationOption operator = (const connectionPairGenerationOption &other)
    {
        manual = other.manual;
        timeTag = other.timeTag;
        connectionPairName = other.connectionPairName;
        return *this;
    }
};

Q_DECLARE_METATYPE(connectionPairGenerationOption)

#endif // CONNECTIONPAIRGENERATIONOPTIONS_H
