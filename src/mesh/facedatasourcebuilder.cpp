//#define USE_MESH_CUTTER
//#define TEST_SURFACE_SMOOTHER

//! ----------------
//! custom includes
//! ----------------
#include <surfacemeshsmoother.h>
#include <facedatasourcebuilder.h>
#include <ng_meshvs_datasourceface.h>
#include <ng_meshvs_datasource2d.h>
#include <ng_meshvs_datasource3d.h>
#include "hash_c.h"

#ifdef USE_MESH_CUTTER
#include <surfacemeshcutter.h>
#endif

//! ----
//! OCC
//! ----
#include <TopoDS.hxx>

//! ----
//! C++
//! ----
#include <map>
#include <iostream>
using namespace std;


//! --------------------------------
//! function: faceDataSourceBuilder
//! details:
//! --------------------------------
faceDataSourceBuilder::faceDataSourceBuilder(QObject *parent): QObject(parent)
{
    ;
}

//! --------------------------------
//! function: faceDataSourceBuilder
//! details:
//! --------------------------------
faceDataSourceBuilder::faceDataSourceBuilder(const QList<TopoDS_Face> &aFaceList, meshDataBase *mDB, QObject *parent):
    myListOfFaces(aFaceList),
    myMDB(mDB),
    QObject(parent)
{
    ;
}

//! ----------------------
//! function: setDataBase
//! details:
//! ----------------------
void faceDataSourceBuilder::setDataBase(meshDataBase *mDB)
{
    myMDB = mDB;
}

//! -------------------
//! function: setFaces
//! details:
//! -------------------
void faceDataSourceBuilder::setFaces(const QList<TopoDS_Face> &faceList)
{
    myListOfFaces.clear();
    for(QList<TopoDS_Face>::const_iterator it = faceList.cbegin(); it!=faceList.cend(); ++it)
    {
        TopoDS_Face aFace = *it;
        myListOfFaces<<aFace;
    }
}

//! ------------------------------
//! function: setMapOfIsMeshExact
//! details:
//! ------------------------------
void faceDataSourceBuilder::setMapOfIsMeshExact(const QMap<int,bool> &aMapOfIsMeshExact)
{
    myMapOfIsMeshDSExact = aMapOfIsMeshExact;
}

//! ------------------------------------------
//! function: setFaces
//! details:  through a vector of GeometryTag
//! ------------------------------------------
void faceDataSourceBuilder::setFaces(const QVector<GeometryTag> &vecLoc)
{
    myListOfFaces.clear();
    for(QVector<GeometryTag>::const_iterator it = vecLoc.cbegin(); it!= vecLoc.cend(); ++it)
    {
        const GeometryTag &loc = *it;
        int bodyIndex = loc.parentShapeNr;
        int faceNr = loc.subTopNr;
        TopoDS_Face aFace = TopoDS::Face(myMDB->MapOfBodyTopologyMap.value(bodyIndex).faceMap.FindKey(faceNr));
        if(!aFace.IsNull()) myListOfFaces<<aFace;
    }
}

//! -------------------------------------------------
//! function: setFaces
//! details:  through a map of vector of GeometryTag
//! -------------------------------------------------
void faceDataSourceBuilder::setFaces(const QMap<int,QVector<GeometryTag>> &vecLocMap)
{
    cout<<"faceDataSourceBuilder::setFaces()->____setting up the face list from input geometry tags____"<<endl;
    for(QMap<int,QVector<GeometryTag>>::const_iterator it = vecLocMap.begin(); it!=vecLocMap.end(); ++it)
    {
        int bcNumber = it.key();
        QVector<GeometryTag> vecLoc = it.value();
        QList<TopoDS_Face> listOfFaces;
        for(QVector<GeometryTag>::const_iterator it = vecLoc.cbegin(); it!= vecLoc.cend(); ++it)
        {
            const GeometryTag &loc = *it;
            int bodyIndex = loc.parentShapeNr;
            int faceNr = loc.subTopNr;
            //cout<<"faceDataSourceBuilder::setFaces()->____setting face with tag ("<<bodyIndex<<", "<<faceNr<<")____"<<endl;
            TopoDS_Face aFace = TopoDS::Face(myMDB->MapOfBodyTopologyMap.value(bodyIndex).faceMap.FindKey(faceNr));

            if(!aFace.IsNull())
            {
                cout<<"faceDataSourceBuilder::setFaces()->____the face with tag("<<bodyIndex<<", "<<faceNr<<") is null____"<<endl;
                continue;
            }
            listOfFaces<<aFace;
        }
        myMapOfListOfFaces.insert(bcNumber,listOfFaces);
    }
}

//! --------------------------------------------------------------------------------
//! function: perform2
//! details:  given the list of tags the class is initialized with, group the faces
//!           according to the parent body; the following map is generated:
//!           key   => body index
//!           value => a mesh data source
//!           If "doExact" is true, the face mesh datasources are taken from the
//!           mesh database, and simpy merged into one mesh datasource; otherwhise,
//!           an approximate mes
//! --------------------------------------------------------------------------------
bool faceDataSourceBuilder::perform2(IndexedMapOfMeshDataSources &mapOfFaceDS, bool doExact)
{
    cout<<"faceDataSourceBuilder::perform()->____function called____"<<endl;
    if(myListOfFaces.isEmpty()) return false;

    //! ------------------------------------------------------
    //! Note: the faces which have been put into the input
    //! variable "faceDSList" can belong to different bodies:
    //! group them according to the parent
    //! ------------------------------------------------------
    cout<<"faceDataSourceBuilder::perform()->____begin grouping faces____"<<endl;
    QMap<int,QList<TopoDS_Face>> bodyFacesMap  = this->groupFaces();

    //! -------------
    //! sanity check
    //! -------------
    if(bodyFacesMap.isEmpty())
    {
        emit taskFinished();
        return false;
    }

    if(doExact)
    {

        for(QMap<int,QList<TopoDS_Face>>::iterator it = bodyFacesMap.begin(); it!=bodyFacesMap.end(); ++it)
        {
            int bodyIndex = it.key();
            if(!bodyFacesMap.contains(bodyIndex)) continue;

            cout<<"faceDataSourceBuilder::perform()->____doing exact for body nr: "<<bodyIndex<<"____"<<endl;

            QList<occHandle(Ng_MeshVS_DataSourceFace)> listOfFaceMeshDS;
            TopTools_IndexedMapOfShape faceMap = myMDB->MapOfBodyTopologyMap.value(bodyIndex).faceMap;

            if(faceMap.IsEmpty()) cout<<"____empty face map____"<<endl;

            const QList<TopoDS_Face> &faces = it.value();
            for(int n=0; n<faces.length(); n++)
            {
                TopoDS_Face aFace = faces.at(n);

                if(aFace.IsNull()) cout<<"____NULL face____"<<endl;

                int faceNr = faceMap.FindIndex(aFace);

                cout<<"____faceNr: "<<faceNr<<"____"<<endl;

                occHandle(Ng_MeshVS_DataSourceFace) aFaceDS = occHandle(Ng_MeshVS_DataSourceFace)::DownCast(myMDB->ArrayOfMeshDSOnFaces.getValue(bodyIndex,faceNr));
                if(aFaceDS.IsNull()) cout<<"____face mesh ds nr: "<<faceNr<<" is null____"<<endl;
                cout<<"____face mesh ds: "<<faceNr<<" OK____"<<endl;

                listOfFaceMeshDS<<aFaceDS;
            }

            cout<<"____tag00____"<<endl;

            if(listOfFaceMeshDS.isEmpty())
            {
                cout<<"____list of face mesh ds empty____"<<endl;
                continue;
            }
            cout<<"____summing "<<listOfFaceMeshDS.length()<<" meshes____"<<endl;

            occHandle(Ng_MeshVS_DataSourceFace) finalFaceDS = new Ng_MeshVS_DataSourceFace(listOfFaceMeshDS);

            cout<<"____tag01____"<<endl;

            if(finalFaceDS.IsNull())
            {
                cout<<"____the final mesh is null____"<<endl;
                continue;
            }
            mapOfFaceDS.insert(bodyIndex,finalFaceDS);
            cout<<"____tag02____"<<endl;
        }
        if(mapOfFaceDS.isEmpty())
        {
            cout<<"____empty mapOfFaces____"<<endl;
            emit taskFinished();
            return false;
        }
        emit taskFinished();
        return true;
    }
    else
    {
        //! ------------------------------------------------------------------------------------------
        //! ALGO
        //! for each body
        //!  for each face of the body
        //!   rebuild the mesh of the current face and put it into a list
        //!   sum the face mesh data sources contained into the list
        //!  insert the resulting mesh into QMap<int,occHandle(Ng_MeshVS_DataSourceFace)> mapOfFaceDS
        //! ------------------------------------------------------------------------------------------
        for(QMap<int,QList<TopoDS_Face>>::iterator it = bodyFacesMap.begin(); it!=bodyFacesMap.end(); ++it)
        {
            int bodyIndex = it.key();
            QList<TopoDS_Face> listOfFaces = it.value();

            //! -----------------------------
            //! the surface mesh of the body
            //! -----------------------------
            occHandle(Ng_MeshVS_DataSource2D) surfaceMesh;

            bool mesh2Dtopological = true;
            if(!mesh2Dtopological)
            {
                //! ----------------------------------------------
                //! data source retrieved from the mesh data base
                //! ----------------------------------------------
                surfaceMesh = occHandle(Ng_MeshVS_DataSource2D)::DownCast(myMDB->ArrayOfMeshDS2D.value(bodyIndex));
            }
            else
            {
                //! -----------------------
                //! topological - in place
                //! -----------------------
                const occHandle(Ng_MeshVS_DataSource3D) &volumeMeshDS =
                        occHandle(Ng_MeshVS_DataSource3D)::DownCast(myMDB->ArrayOfMeshDS.value(bodyIndex));
                if(volumeMeshDS.IsNull())
                {
                    cout<<"faceDataSourceBuilder::perform()->____cannot generate the surface mesh from the volume mesh: "<<endl;
                    return false;
                }
                if(volumeMeshDS->myFaceToElements.isEmpty())
                {
                    cout<<"faceDataSourceBuilder::perform()->____the face to element connectivity info have been not generated yet. Building____"<<endl;
                    volumeMeshDS->buildFaceToElementConnectivity();
                }
                surfaceMesh = new Ng_MeshVS_DataSource2D(volumeMeshDS);

                //! replace the mesh 2D - test
                //! myMDB->ArrayOfMeshDS2D.insert(bodyIndex,surfaceMesh);
            }
            if(surfaceMesh.IsNull()) continue;

    #ifdef USE_MESH_CUTTER
            //! ---------------------
            //! cut the surface mesh
            //! ---------------------
            occHandle(MeshVS_DataSource) cutMesh;
            QList<TopoDS_Shape> ls;
            for(int i=0; i<listOfFaces.length(); i++) ls<<listOfFaces.at(i);
            surfaceMeshCutter::cutSurfaceMesh(surfaceMesh,ls,cutMesh);
    #endif

            //! --------------------------------------------
            //! rebuild the "approximate" MeshVS_DataSource
            //! and put it into a list "listOfMeshDS"
            //! --------------------------------------------
            QList<occHandle(Ng_MeshVS_DataSourceFace)> listOfFaceMeshDS;
            for(int i=0; i<listOfFaces.length(); i++)
            {
                TopoDS_Face curFace = listOfFaces.at(i);
                bool autocomputeProximity = false;
                double proximity = 0.1;
    #ifndef USE_MESH_CUTTER
                occHandle(Ng_MeshVS_DataSourceFace) curFaceRebuiltDS = new Ng_MeshVS_DataSourceFace(curFace,surfaceMesh,autocomputeProximity,proximity);
    #endif
    #ifdef USE_MESH_CUTTER
                occHandle(Ng_MeshVS_DataSourceFace) curFaceRebuiltDS = new Ng_MeshVS_DataSourceFace(curFace,cutMesh,autocomputeProximity,proximity);
    #endif
                //! -------------
                //! sanity check
                //! -------------
                if(curFaceRebuiltDS.IsNull()) continue;
                listOfFaceMeshDS<<curFaceRebuiltDS;
            }

            //! -------------
            //! sanity check
            //! -------------
            if(listOfFaces.isEmpty())
            {
                emit taskFinished();
                return false;
            }

            //! -----------------------------------------------------------
            //! sum the face mesh computed on the current body "bodyIndex"
            //! When retrieving face mesh data sources from the map use
            //! Ng_MeshVS_DataSourceFace::DownCast
            //! -----------------------------------------------------------
            occHandle(Ng_MeshVS_DataSourceFace) finalFaceDS = new Ng_MeshVS_DataSourceFace(listOfFaceMeshDS);

            mapOfFaceDS.insert(bodyIndex, finalFaceDS);
        }
        //cout<<"faceDataSourceBuilder::perform()->____data source created____"<<endl;
        emit taskFinished();
        return true;
    }

    cout<<"____exiting____"<<endl;
}

//! -------------------------------------------------------------
//! function: groupFaces
//! details:  group a list of faces according to the parent body
//! -------------------------------------------------------------
QMap<int, QList<TopoDS_Face> > faceDataSourceBuilder::groupFaces()
{
    cout<<"faceDataSourceBuilder::groupFaces()->____grouping faces____"<<endl;
    QMap<int,QList<TopoDS_Face>> bodyFacesMap;
    for(int i=0; i<myListOfFaces.length(); i++)
    {
        const TopoDS_Face &curFace = myListOfFaces.at(i);
        if(curFace.IsNull()) exit(999999);

        //! --------------------------
        //! find the body of the face
        //! --------------------------
        int bodyIndex;
        bool bodyOfFaceFound = false;
        for(QMap<int,TopoDS_Shape>::iterator it = myMDB->bodyMap.begin(); it!=myMDB->bodyMap.end(); ++it)
        {
            bodyIndex = it.key();
            const TopTools_IndexedMapOfShape &faceMap = myMDB->MapOfBodyTopologyMap.value(bodyIndex).faceMap;
            if(faceMap.Contains(curFace))
            {
                bodyOfFaceFound = true;
                break;
            }
        }
        if(bodyOfFaceFound == false) continue;

        //! --------------------------------------------------
        //! put the face into the map at position "bodyIndex"
        //! --------------------------------------------------
        if(bodyFacesMap.contains(bodyIndex))
        {
            QList<TopoDS_Face> listOfFaces = bodyFacesMap.value(bodyIndex);
            listOfFaces<<curFace;
            bodyFacesMap.insert(bodyIndex,listOfFaces);
        }
        else
        {
            QList<TopoDS_Face> listOfFaces;
            listOfFaces<<curFace;
            bodyFacesMap.insert(bodyIndex,listOfFaces);
        }
    }
    cout<<"faceDataSourceBuilder::groupFaces()->____grouping faces: exiting____"<<endl;
    return bodyFacesMap;
}
