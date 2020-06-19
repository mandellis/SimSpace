#ifndef MAPOFMESHDATASOURCES_H
#define MAPOFMESHDATASOURCES_H

#include "geometrytag.h"
#include "occhandle.h"
#include <MeshVS_DataSource.hxx>
#include <QMap>
#include <QMetaType>

typedef QMap<GeometryTag, occHandle(MeshVS_DataSource)> mapOfMeshDataSources;

Q_DECLARE_METATYPE(mapOfMeshDataSources)

#endif // MAPOFMESHDATASOURCES_H


