#ifndef FACEDATASOURCEBUILDER_H
#define FACEDATASOURCEBUILDER_H

//! ------------------------------------
//! replacement for opencascade::handle
//! ------------------------------------
#include "occhandle.h"

//! ----------------
//! custom includes
//! ----------------
#include <indexedmapofmeshdatasources.h> //! typedef for QMap<int,occHandle(MeshVS_DataSource)>
#include <meshdatabase.h>
//#include <ng_meshvs_datasourceface.h>
#include "geometrytag.h"

//! ----
//! OCC
//! ----
#include <TopoDS_Face.hxx>

//! ---
//! Qt
//! ---
#include <QList>
#include <QMap>
#include <QObject>

class faceDataSourceBuilder: public QObject
{
    Q_OBJECT

public:

    faceDataSourceBuilder(QObject *parent = 0);
    faceDataSourceBuilder(const QList<TopoDS_Face> &aFaceList, meshDataBase *mDB, QObject *parent = 0);

    void setFaces(const QList<TopoDS_Face> &faceList);
    void setFaces(const QVector<GeometryTag> &vecLoc);
    void setFaces(const QMap<int, QVector<GeometryTag>> &vecLocMap);
    void setMapOfIsMeshExact(const QMap<int,bool> &aMapOfIsMeshExact);
    void setDataBase(meshDataBase *mDB);

public slots:

    bool perform2(IndexedMapOfMeshDataSources &mapOfFaceDS, bool doExact);

private:

    meshDataBase *myMDB;
    QMap<int,bool> myMapOfIsMeshDSExact;
    QList<TopoDS_Face> myListOfFaces;
    QMap<int,QList<TopoDS_Face>> myMapOfListOfFaces;

private:

    QMap<int,QList<TopoDS_Face>> groupFaces();

signals:

    void taskFinished();
};

#endif // FACEDATASOURCEBUILDER_H
