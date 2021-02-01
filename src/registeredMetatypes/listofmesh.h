#ifndef LISTOFMESH_H
#define LISTOFMESH_H

#include <MeshVS_Mesh.hxx>
#include <NCollection_List.hxx>
#include <QMetaType>

typedef NCollection_List<occHandle(MeshVS_Mesh)> ListOfMesh;

Q_DECLARE_METATYPE(ListOfMesh)

#endif // LISTOFMESH_H
