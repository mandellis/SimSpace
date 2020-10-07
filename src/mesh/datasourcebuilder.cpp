//#define USE_MESH_CUTTER
//#define TEST_SURFACE_SMOOTHER

//! ----------------
//! custom includes
//! ----------------
#include <surfacemeshsmoother.h>
#include <datasourcebuilder.h>
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
//! function: dataSourceBuilder
//! details:
//! --------------------------------
dataSourceBuilder::dataSourceBuilder(QObject *parent): QObject(parent)
{
    ;
}
/*
//! --------------------------------
//! function: dataSourceBuilder
//! details:
//! --------------------------------
dataSourceBuilder::dataSourceBuilder(const QList<TopoDS_Face> &aFaceList, meshDataBase *mDB, QObject *parent):
    myListOfFaces(aFaceList),
    myMDB(mDB),
    QObject(parent)
{
    ;
}
*/
//! ----------------------
//! function: setDataBase
//! details:
//! ----------------------
void dataSourceBuilder::setDataBase(meshDataBase *mDB)
{
    myMDB = mDB;
}
/*
//! -------------------
//! function: setFaces
//! details:
//! -------------------
void dataSourceBuilder::setFaces(const QList<TopoDS_Face> &faceList)
{
    myListOfFaces.clear();
    for(QList<TopoDS_Face>::const_iterator it = faceList.cbegin(); it!=faceList.cend(); ++it)
    {
        TopoDS_Face aFace = *it;
        myListOfFaces<<aFace;
    }
}
*/

//! ------------------------------
//! function: setMapOfIsMeshExact
//! details:
//! ------------------------------
void dataSourceBuilder::setMapOfIsMeshExact(const QMap<int,bool> &aMapOfIsMeshExact)
{
    myMapOfIsMeshDSExact = aMapOfIsMeshExact;
}
/*
//! ------------------------------------------
//! function: setFaces
//! details:  through a vector of GeometryTag
//! ------------------------------------------
void dataSourceBuilder::setFaces(const std::vector<GeometryTag> &vecLoc)
{
    myListOfFaces.clear();
    for(std::vector<GeometryTag>::const_iterator it = vecLoc.cbegin(); it!= vecLoc.cend(); ++it)
    {
        const GeometryTag &loc = *it;
        int bodyIndex = loc.parentShapeNr;
        int faceNr = loc.subTopNr;
        TopoDS_Face aFace = TopoDS::Face(myMDB->MapOfBodyTopologyMap.value(bodyIndex).faceMap.FindKey(faceNr));
        if(!aFace.IsNull()) myListOfFaces<<aFace;
    }
}
*/
//! ------------------------------------------
//! function: setScopes
//! details:  through a vector of GeometryTag
//! ------------------------------------------
void dataSourceBuilder::setShapes(const std::vector<GeometryTag> &vecLoc)
{
    cout<<"faceDataSourceBuilder::setShapes()->____function called___"<<endl;

    myListOfShape.clear();
    for(std::vector<GeometryTag>::const_iterator it = vecLoc.cbegin(); it!= vecLoc.cend(); ++it)
    {
        cout<<"faceDataSourceBuilder::setShapes()->____tag00___"<<endl;

        const GeometryTag &loc = *it;
        int bodyIndex = loc.parentShapeNr;
        int faceNr = loc.subTopNr;
        TopAbs_ShapeEnum shapeType = loc.subShapeType;
        switch(shapeType)
        {
        case TopAbs_SOLID:
        {
            cout<<"faceDataSourceBuilder::setShapes()->____tag01___"<<endl;
            TopoDS_Solid aSolid = TopoDS::Solid(myMDB->bodyMap.value(bodyIndex));
            if(!aSolid.IsNull()) myListOfShape.push_back(aSolid);
            cout<<"faceDataSourceBuilder::setShapes()->____tag01___"<<endl;

        }
            break;

        case TopAbs_FACE:
        {
            TopoDS_Face aFace = TopoDS::Face(myMDB->MapOfBodyTopologyMap.value(bodyIndex).faceMap.FindKey(faceNr));
            if(!aFace.IsNull()) myListOfShape.push_back(aFace);
        }
            break;

        default:
        {

        }
            break;

        }

    }
    cout<<"faceDataSourceBuilder::setShapes()->____function exiting___"<<endl;
}


//! -------------------------------------------------------------
//! function: groupShapes
//! details:  group a list of shape according to the parent body
//! -------------------------------------------------------------
std::map<int, std::vector<TopoDS_Shape>> dataSourceBuilder::groupShapes()
{
    cout<<"faceDataSourceBuilder::groupShapes()->____grouping shapes___"<<endl;
    std::map<int,std::vector<TopoDS_Shape>> bodyShapesMap;

    for(std::vector<TopoDS_Shape>::const_iterator i=myListOfShape.cbegin(); i!=myListOfShape.cend();i++)
    {
        const TopoDS_Shape &curShape = *i;
        if(curShape.IsNull()) exit(999999);
        TopAbs_ShapeEnum sType = curShape.ShapeType();
        switch (sType)
        {
        case TopAbs_SOLID:
        {
            //! ----------------
            //! find the body
            //! ---------------
            int bodyIndex = myMDB->bodyMap.key(curShape);

            //! --------------------------------------------------
            //! put the body into the map at position "bodyIndex"
            //! --------------------------------------------------
            std::vector<TopoDS_Shape> listOfShape;
            listOfShape.push_back(curShape);
            bodyShapesMap.insert(std::make_pair(bodyIndex,listOfShape));
        }
            break;
        case TopAbs_FACE:
        {
            //! --------------------------
            //! find the body of the face
            //! --------------------------
            int bodyIndex;
            bool bodyOfFaceFound = false;
            for(QMap<int,TopoDS_Shape>::iterator it = myMDB->bodyMap.begin(); it!=myMDB->bodyMap.end(); ++it)
            {
                TopoDS_Face curFace = TopoDS::Face(curShape);
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
            std::map<int,std::vector<TopoDS_Shape>>::iterator bodyIt =bodyShapesMap.find(bodyIndex);
            if(bodyIt!=bodyShapesMap.end())
            {
                std::vector<TopoDS_Shape> listOfShape = bodyIt->second;
                listOfShape.push_back(curShape);
                bodyShapesMap.insert(std::make_pair(bodyIndex,listOfShape));
            }
            else
            {
                std::vector<TopoDS_Shape> listOfShape;
                listOfShape.push_back(curShape);
                bodyShapesMap.insert(std::make_pair(bodyIndex,listOfShape));
            }
        }
            break;
        default:
            break;
        }
    }
    cout<<"faceDataSourceBuilder::groupShapes()->____grouping shapes: exiting____"<<endl;
    return bodyShapesMap;
}


//! --------------------------------------------------------------------------------
//! function: perform
//! details:  given the list of tags the class is initialized with, group the shape
//!           according to the parent body; the following map is generated:
//!           key   => body index
//!           value => a mesh data source
//!           If "doExact" is true, the mesh datasources are taken from the
//!           mesh database, and simpy merged into one mesh datasource; otherwhise,
//!           an approximate meshDS is created
//! --------------------------------------------------------------------------------
bool dataSourceBuilder::perform(IndexedMapOfMeshDataSources &mapOfDS, bool doExact)
{
    cout<<"faceDataSourceBuilder::perform()->____function called____"<<endl;
    if(myListOfShape.empty()) return false;

    //! ------------------------------------------------------
    //! Note: the faces which have been put into the input
    //! variable "faceDSList" can belong to different bodies:
    //! group them according to the parent
    //! ------------------------------------------------------
    cout<<"faceDataSourceBuilder::perform()->____begin grouping shapes____"<<endl;
    std::map<int,std::vector<TopoDS_Shape>> bodyShapesMap  = this->groupShapes();

    //! -------------
    //! sanity check
    //! -------------
    if(bodyShapesMap.empty())
    {
        emit taskFinished();
        cout<<"faceDataSourceBuilder::perform()->____bodyShapesMap finished____"<<endl;
        return false;
    }

    if(doExact)
    {
        for(std::map<int,std::vector<TopoDS_Shape>>::iterator it = bodyShapesMap.begin(); it!=bodyShapesMap.end(); ++it)
        {
            int bodyIndex = it->first;
            if(bodyShapesMap.count(bodyIndex)==0) continue;
            cout<<"faceDataSourceBuilder::perform()->____doing exact for body nr: "<<bodyIndex<<"____"<<endl;

            const std::vector<TopoDS_Shape> &shapes = it->second;
            TopoDS_Shape firstShape = shapes.at(0);
            switch(firstShape.ShapeType())
            {
            case TopAbs_SOLID:
            {
                for(int n=0; n<shapes.size(); n++)
                {
                    TopoDS_Shape aShape = shapes[n];
                    if(aShape.IsNull()) cout<<"faceDataSourceBuilder::perform1()->____NULL body____"<<endl;
                    int bodyIndex = myMDB->bodyMap.key(aShape);
                    cout<<"faceDataSourceBuilder::perform1()->____working on body nr: "<<bodyIndex<<"____"<<endl;

                    occHandle(Ng_MeshVS_DataSource3D) aBodyDS = occHandle(Ng_MeshVS_DataSource3D)::DownCast(myMDB->ArrayOfMeshDS.value(bodyIndex));
                    if(aBodyDS.IsNull())
                    {
                        continue;
                    }
                    mapOfDS.insert(bodyIndex,aBodyDS);
                }
                if(mapOfDS.isEmpty())
                {
                    cout<<"faceDataSourceBuilder::perform()->____mapOfFaceDS empty____"<<endl;
                    emit taskFinished();
                    return false;
                }
                cout<<"faceDataSourceBuilder::perform()->____function exiting____"<<endl;
                emit taskFinished();
                return true;
            }
                break;
            case TopAbs_FACE:
            {
                QList<occHandle(Ng_MeshVS_DataSourceFace)> listOfFaceMeshDS;
                TopTools_IndexedMapOfShape faceMap = myMDB->MapOfBodyTopologyMap.value(bodyIndex).faceMap;
                if(faceMap.IsEmpty()) cerr<<"faceDataSourceBuilder::perform1()->____strange error: empty face map in geometry____"<<endl;

                for(int n=0; n<shapes.size(); n++)
                {
                    TopoDS_Shape aShape = shapes[n];
                    if(aShape.IsNull()) cout<<"faceDataSourceBuilder::perform1()->____NULL face____"<<endl;
                    int faceNr = faceMap.FindIndex(aShape);

                    cout<<"faceDataSourceBuilder::perform1()->____working on face nr: "<<faceNr<<"____"<<endl;

                    occHandle(Ng_MeshVS_DataSourceFace) aFaceDS = occHandle(Ng_MeshVS_DataSourceFace)::DownCast(myMDB->ArrayOfMeshDSOnFaces.getValue(bodyIndex,faceNr));
                    if(aFaceDS.IsNull())
                    {
                        cout<<"faceDataSourceBuilder::perform1()->____the mesh data source of face nr: "<<faceNr<<" is null: jumping over it____"<<endl;
                        continue;
                    }
                    listOfFaceMeshDS<<aFaceDS;
                }

                if(listOfFaceMeshDS.isEmpty())
                {
                    cout<<"faceDataSourceBuilder::perform1()->____list of face mesh ds empty____"<<endl;
                    continue;
                }
                cout<<"faceDataSourceBuilder::perform1()->____summing "<<listOfFaceMeshDS.size()<<" meshes____"<<endl;
                occHandle(Ng_MeshVS_DataSourceFace) finalFaceDS = new Ng_MeshVS_DataSourceFace(listOfFaceMeshDS);
                cout<<"faceDataSourceBuilder::perform1()->____summation done____"<<endl;

                if(finalFaceDS.IsNull())
                {
                    cout<<"____the final mesh is null____"<<endl;
                    continue;
                }
                mapOfDS.insert(bodyIndex,finalFaceDS);
            }
            if(mapOfDS.isEmpty())
            {
                cout<<"faceDataSourceBuilder::perform1()->____mapOfFaceDS empty____"<<endl;
                emit taskFinished();
                return false;
            }
            cout<<"faceDataSourceBuilder::perform1()->____function exiting____"<<endl;
            emit taskFinished();
            return true;
            }
                break;
            }

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
        for(std::map<int,std::vector<TopoDS_Shape>>::iterator it = bodyShapesMap.begin(); it!=bodyShapesMap.end(); ++it)
        {
            int bodyIndex = it->first;
            std::vector<TopoDS_Shape> listOfShape = it->second;

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
            for(int i=0; i<listOfShape.size(); i++)
            {
                TopoDS_Face curFace = TopoDS::Face(listOfShape.at(i));
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
            if(listOfShape.empty())
            {
                cout<<"faceDataSourceBuilder::perform()->____function exiting____"<<endl;
                emit taskFinished();
                return false;
            }

            //! -----------------------------------------------------------
            //! sum the face mesh computed on the current body "bodyIndex"
            //! When retrieving face mesh data sources from the map use
            //! Ng_MeshVS_DataSourceFace::DownCast
            //! -----------------------------------------------------------
            occHandle(Ng_MeshVS_DataSourceFace) finalFaceDS = new Ng_MeshVS_DataSourceFace(listOfFaceMeshDS);

            mapOfDS.insert(bodyIndex, finalFaceDS);
        }
        cout<<"faceDataSourceBuilder::perform()->____function exiting____"<<endl;
        emit taskFinished();
        return true;
    }
}

/*
//! -------------------------------------------------
//! function: setFaces
//! details:  through a map of vector of GeometryTag
//! -------------------------------------------------
void dataSourceBuilder::setFaces(const QMap<int,std::vector<GeometryTag>> &vecLocMap)
{
    cout<<"faceDataSourceBuilder::setFaces()->____setting up the face list from input geometry tags____"<<endl;
    for(QMap<int,std::vector<GeometryTag>>::const_iterator it = vecLocMap.begin(); it!=vecLocMap.end(); ++it)
    {
        int bcNumber = it.key();
        std::vector<GeometryTag> vecLoc = it.value();
        QList<TopoDS_Face> listOfFaces;
        for(std::vector<GeometryTag>::const_iterator it = vecLoc.cbegin(); it!= vecLoc.cend(); ++it)
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
*/
/*
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

            if(faceMap.IsEmpty()) cerr<<"faceDataSourceBuilder::perform()->____strange error: empty face map in geometry____"<<endl;

            const QList<TopoDS_Face> &faces = it.value();
            for(int n=0; n<faces.length(); n++)
            {
                TopoDS_Face aFace = faces.at(n);

                if(aFace.IsNull()) cout<<"faceDataSourceBuilder::perform()->____NULL face____"<<endl;

                int faceNr = faceMap.FindIndex(aFace);

                cout<<"faceDataSourceBuilder::perform()->____working on face nr: "<<faceNr<<"____"<<endl;

                occHandle(Ng_MeshVS_DataSourceFace) aFaceDS = occHandle(Ng_MeshVS_DataSourceFace)::DownCast(myMDB->ArrayOfMeshDSOnFaces.getValue(bodyIndex,faceNr));
                if(aFaceDS.IsNull())
                {
                    cout<<"faceDataSourceBuilder::perform()->____the mesh data source of face nr: "<<faceNr<<" is null: jumping over it____"<<endl;
                    continue;
                }
                listOfFaceMeshDS<<aFaceDS;
            }

            if(listOfFaceMeshDS.isEmpty())
            {
                cout<<"faceDataSourceBuilder::perform()->____list of face mesh ds empty____"<<endl;
                continue;
            }
            cout<<"faceDataSourceBuilder::perform()->____summing "<<listOfFaceMeshDS.length()<<" meshes____"<<endl;
            occHandle(Ng_MeshVS_DataSourceFace) finalFaceDS = new Ng_MeshVS_DataSourceFace(listOfFaceMeshDS);
            cout<<"faceDataSourceBuilder::perform()->____summation done____"<<endl;


            if(finalFaceDS.IsNull())
            {
                cout<<"____the final mesh is null____"<<endl;
                continue;
            }
            mapOfFaceDS.insert(bodyIndex,finalFaceDS);
        }
        if(mapOfFaceDS.isEmpty())
        {
            cout<<"faceDataSourceBuilder::perform()->____function exiting____"<<endl;
            emit taskFinished();
            return false;
        }
        cout<<"faceDataSourceBuilder::perform()->____function exiting____"<<endl;
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
                cout<<"faceDataSourceBuilder::perform()->____function exiting____"<<endl;
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
        cout<<"faceDataSourceBuilder::perform()->____function exiting____"<<endl;
        emit taskFinished();
        return true;
    }
}

//! -------------------------------------------------------------
//! function: groupFaces
//! details:  group a list of faces according to the parent body
//! -------------------------------------------------------------
QMap<int, QList<TopoDS_Face>> faceDataSourceBuilder::groupFaces()
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
*/
