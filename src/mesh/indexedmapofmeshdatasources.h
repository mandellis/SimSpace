#ifndef INDEXEDMAPOFMESHDATASOURCES_H
#define INDEXEDMAPOFMESHDATASOURCES_H

#include <QMap>
#include <QMetaType>

#include "occhandle.h"
#include <MeshVS_DataSource.hxx>

typedef QMap<int,occHandle(MeshVS_DataSource)> IndexedMapOfMeshDataSources;

Q_DECLARE_METATYPE(IndexedMapOfMeshDataSources)

#endif // INDEXEDMAPOFMESHDATASOURCES_H


// perform(MeshVS_DataSource, std::map<std::map<std::string,IndexedMapOfMeshDataSources>>)
