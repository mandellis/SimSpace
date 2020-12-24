#ifndef DATASOURCEBUILDER_H
#define DATASOURCEBUILDER_H

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
#include <TopoDS_Solid.hxx>

//! ---
//! Qt
//! ---
#include <QList>
#include <QMap>
#include <QObject>

class dataSourceBuilder: public QObject
{
    Q_OBJECT

public:

    dataSourceBuilder(QObject *parent = 0);
    //dataSourceBuilder(const QList<TopoDS_Face> &aFaceList, meshDataBase *mDB, QObject *parent = 0);

    //void setFaces(const QList<TopoDS_Face> &faceList);
    void setFaces(const std::vector<GeometryTag> &vecLoc);
    //void setFaces(const QMap<int, std::vector<GeometryTag>> &vecLocMap);

    void setShapes(const std::vector<GeometryTag> &vecLoc);

    void setMapOfIsMeshExact(const QMap<int,bool> &aMapOfIsMeshExact);
    void setDataBase(meshDataBase *mDB);

public slots:

    //bool perform2(IndexedMapOfMeshDataSources &mapOfFaceDS, bool doExact);
    bool perform(IndexedMapOfMeshDataSources &mapOfDS, bool doExact);

private:

    meshDataBase *myMDB;
    QMap<int,bool> myMapOfIsMeshDSExact;
    //QList<TopoDS_Face> myListOfFaces;
    //QList<TopoDS_Solid> myListOfBodies;

    std::vector<TopoDS_Shape> myListOfShape;

    //QMap<int,QList<TopoDS_Face>> myMapOfListOfFaces;
    std::map<int,std::vector<TopoDS_Shape>> myMapOfListOfShape;

private:

    //QMap<int,QList<TopoDS_Face>> groupFaces();
    std::map<int,std::vector<TopoDS_Shape>> groupShapes();


signals:

    void taskFinished();
};

#endif // DATASOURCEBUILDER_H
