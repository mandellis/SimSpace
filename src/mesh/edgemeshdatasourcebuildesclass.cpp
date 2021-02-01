//! custom includes
#include <edgemeshdatasourcebuildesclass.h>
#include <ng_meshvs_datasource1d.h>
#include <ng_meshvs_datasource3d.h>
#include <ng_meshvs_datasourceface.h>
#include <src/utils/meshtools.h>

//! OCC
#include <TopExp.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopTools_IndexedDataMapOfShapeShape.hxx>
#include <TopTools_ListIteratorOfListOfShape.hxx>
#include <TColStd_MapIteratorOfPackedMapOfInteger.hxx>
#include <TopTools_MapIteratorOfMapOfShape.hxx>
#include <TColStd_PackedMapOfInteger.hxx>
#include <TopoDS.hxx>
#include <MeshVS_DataSource.hxx>
#include <BRep_Tool.hxx>

//! EXPRESS MESH
#include <OMFVS_DataSource.hxx>

//! C++
#include <vector>

#define VERBOSE
//#define NEW_EDGE_DS_BUILDER

#ifndef NEW_EDGE_DS_BUILDER
//! ------------------------------------------------------------------------
//! An edge is defined as the intersection of two faces
//! Here the structure "shapePair" contains:
//!
//! 1) edgeNr = the number of the current edge (number in the body edgeMap)
//! 2) eln1 = the first face number (number in the body faceMap)
//! 3) eln2 = the second face number
//!
//! ------------------------------------------------------------------------
struct shapePair
{
    int edgeNr;
    int eln1;
    int eln2;
};

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
edgeMeshDataSourceBuildesClass::edgeMeshDataSourceBuildesClass(meshDataBase *theMeshDataBase, QObject *parent):
    QObject(parent), myDB(theMeshDataBase)
{
    cout<<"edgeMeshDataSourceBuildesClass::edgeMeshDataSourceBuildesClass()->____CONSTRUCTOR CALLED____"<<endl;
}

//! ------------------------------------------------------ //
//! function: perform                                      //
//! details:                                               //
//! -------------------------------------------------------//
void edgeMeshDataSourceBuildesClass::perform()
{
    cout<<"edgeMeshDataSourceBuildesClass::perform()->____function called____"<<endl;
    for(int bodyIndex = 1; bodyIndex <= myDB->bodyMap.size(); bodyIndex++)
    {
        const TopoDS_Shape &theShape = myDB->bodyMap.value(bodyIndex);
        //! number of edges for of the current body
        int NEdges = myDB->MapOfBodyTopologyMap.value(bodyIndex).edgeMap.Extent();

        NCollection_Array1<occHandle(Ng_MeshVS_DataSource1D)> edgeMeshArray(1,NEdges);
        if(this->performOne(theShape, edgeMeshArray))
        {
            //! put the edge mesh in the simulation database
            for(int edgeNr = 1; edgeNr<=edgeMeshArray.Length(); edgeNr++)
                myDB->ArrayOfMeshDSOnEdges.setValue(bodyIndex, edgeNr, edgeMeshArray.Value(edgeNr));
        }
    }
}

//! -----------------------------------
//! function: perform on a single body
//! details:
//! -----------------------------------
Standard_Boolean edgeMeshDataSourceBuildesClass::performOne(const TopoDS_Shape &theShape,
                                                            NCollection_Array1<opencascade::handle<Ng_MeshVS_DataSource1D>> &edgeMeshArray)
{
    int bodyIndex = myDB->bodyMap.key(theShape);
    if(!myDB->ArrayOfMeshDS2D.value(bodyIndex).IsNull())
    {
#ifdef VERBOSE
        cout<<"edgeMeshDataSourceBuildesClass::performOnBody->____working on body N. "<<myDB->bodyMap.key(theShape)<<"____"<<endl;
#endif

        std::vector<shapePair> aVectorOfPair;
        TopTools_IndexedDataMapOfShapeListOfShape edgeFaceMap;
        TopExp::MapShapesAndAncestors(theShape,TopAbs_EDGE,TopAbs_FACE,edgeFaceMap);

        //! obtain the face and edge maps
        TopTools_IndexedMapOfShape edgeMap = myDB->MapOfBodyTopologyMap.value(bodyIndex).edgeMap;
        TopTools_IndexedMapOfShape faceMap = myDB->MapOfBodyTopologyMap.value(bodyIndex).faceMap;

        //! scan all the faces
        int theTopLevelFaceIndex;
        int theOtherFaceIndex;

        TopTools_MapIteratorOfMapOfShape faceMapIt;
        for(faceMapIt.Initialize(faceMap);faceMapIt.More();faceMapIt.Next())
        {
            const TopoDS_Shape &theTopLevelFace = faceMapIt.Value();

            //! scan all the edges
            for(int edgeNr=1; edgeNr<=edgeMap.Extent(); edgeNr++)
            {
                //! the current edge
                const TopoDS_Shape &theCurEdge = edgeMap.FindKey(edgeNr);

                bool isSeam = BRep_Tool::IsClosed(TopoDS::Edge(theCurEdge),TopoDS::Face(theTopLevelFace));
                bool isDegenerated = BRep_Tool::Degenerated(TopoDS::Edge(theCurEdge));
                if(!isSeam && !isDegenerated)
                {
                    //! the edge is defined by two faces
                    TopTools_ListOfShape theTwoFaces;

                    //! find the two faces defining the edge through intersection
                    edgeFaceMap.FindFromKey(theCurEdge, theTwoFaces);

                    if(theTwoFaces.Contains(theTopLevelFace))
                    {
                        TopTools_ListIteratorOfListOfShape anIter;
                        anIter.Initialize(theTwoFaces);
                        for(;anIter.More();anIter.Next())
                        {
                            //if(TopoDS::Face(anIter.Value())!=TopoDS::Face(theTopLevelFace))
                            if(anIter.Value()!=theTopLevelFace) theOtherFaceIndex = faceMap.FindIndex(anIter.Value());
                            else theTopLevelFaceIndex = faceMap.FindIndex(anIter.Value());
                        }
                        //! cout<<"edge N. "<<edgeNr<<" defined by faces ("<<theTopLevelFaceIndex<<", "<<theOtherFaceIndex<<")"<<endl;
                        shapePair aPair;
                        aPair.edgeNr = edgeNr;
                        aPair.eln1 = theTopLevelFaceIndex;
                        aPair.eln2 = theOtherFaceIndex;

                        //! append to the vector
                        aVectorOfPair.push_back(aPair);
                    }
                }
                else
                {
                    if(isSeam)
                    {
                        cout<<"____the edge nr: "<<edgeNr<<" is seam____"<<endl;
                        edgeMeshArray.SetValue(edgeNr,occHandle(Ng_MeshVS_DataSource1D)());
                        continue;
                    }
                    if(isDegenerated)
                    {
                        cout<<"____the edge nr: "<<edgeNr<<" is degenerate____"<<endl;
                        edgeMeshArray.SetValue(edgeNr,occHandle(Ng_MeshVS_DataSource1D)());
                        continue;
                    }
                }
            }
        }

        //! -------------------------
        //! eliminate the duplicates
        //! -------------------------
        if(aVectorOfPair.size()!=0)
        {
            //cout<<"edgeMeshDataSourceBuildesClass::performOnBody()->____begin eliminating dulicated values____"<<endl;
            std::vector<shapePair> aCleanPairVector;
            aCleanPairVector.push_back(aVectorOfPair[0]);

            for(int n=0; n<aVectorOfPair.size(); n++)
            {
                //! this eliminates the non unique values: it is O(N^2).
                //! This is not a problem, however use std::set
                int d=0;
                for(int k=0; k<aCleanPairVector.size(); k++)
                {
                    //! comparison
                    bool eq;
                    int edgeNr1 = aVectorOfPair[n].edgeNr;
                    int eln1 = aVectorOfPair[n].eln1;
                    int eln2 = aVectorOfPair[n].eln2;
                    int edgeNr1_ = aCleanPairVector[k].edgeNr;
                    int eln1_ = aCleanPairVector[k].eln1;
                    int eln2_ = aCleanPairVector[k].eln2;

                    (eln1 == eln1_ && eln2 == eln2_ && edgeNr1 == edgeNr1_) || (eln1 == eln2_ && eln2 == eln1_ && edgeNr1 == edgeNr1_)?
                                eq = true : eq = false;
                    //! end comparison
                    if(eq == false)d++;
                }
                if(d==aCleanPairVector.size())
                {
                    aCleanPairVector.push_back(aVectorOfPair[n]);
                }
            }

            //!cout<<"edgeMeshDataSourceBuildesClass::performOnBody()->____duplicate cleared____"<<endl;

            //! ------------------------------------------------------------
            //! scan all the edges using {number of the edge, face1, face2}
            //! ------------------------------------------------------------
            std::vector<shapePair>::iterator it;
            for(it = aCleanPairVector.begin(); it!=aCleanPairVector.end(); ++it)
            {                                
                //! the number of the current edge whose nodes are searched for
                const shapePair &p = *it;
                int edgeNr = p.edgeNr;
#ifdef VERBOSE
                cout<<"____working on edgeNr: "<<edgeNr<<"____"<<endl;
#endif

                //! the indexes of the two face MeshVS_DataSource in the simulation database
                int face1Index = p.eln1;
                int face2Index = p.eln2;

                TColStd_PackedMapOfInteger theMap1, theMap2;

                const occHandle(Ng_MeshVS_DataSourceFace) &faceMeshDS1 = occHandle(Ng_MeshVS_DataSourceFace)::DownCast(myDB->ArrayOfMeshDSOnFaces.getValue(bodyIndex,face1Index));
                const occHandle(Ng_MeshVS_DataSourceFace) &faceMeshDS2 = occHandle(Ng_MeshVS_DataSourceFace)::DownCast(myDB->ArrayOfMeshDSOnFaces.getValue(bodyIndex,face2Index));

#ifdef VERBOSE
                if(!faceMeshDS1.IsNull())cout<<"+++++"<<faceMeshDS1->GetAllNodes().Extent()<<endl;
                if(!faceMeshDS2.IsNull())cout<<"-----"<<faceMeshDS2->GetAllNodes().Extent()<<endl;
#endif

                if(!faceMeshDS1.IsNull() && !faceMeshDS2.IsNull())
                {
                    theMap1 =faceMeshDS1->GetAllNodes();
                    theMap2 =faceMeshDS2->GetAllNodes();
                }
                else
                {
                    cout<<"____error in retrieving the face mesh datasources____"<<endl;
                    edgeMeshArray.SetValue(edgeNr,occHandle(Ng_MeshVS_DataSource1D)());
                    continue;
                }

                //! this map contains the list of common nodes => the nodes of the edge
                TColStd_PackedMapOfInteger theIntersectionMap;
                theIntersectionMap.Intersection(theMap1,theMap2);

                if(theIntersectionMap.Extent()<2)
                {
                    cout<<"____empty or too small (<2) intersection map____"<<endl;
                    edgeMeshArray.SetValue(edgeNr,occHandle(Ng_MeshVS_DataSource1D)());
                    continue;
                }
                else
                {
                    cout<<"____map intersection performed. NN: "<<theIntersectionMap.Extent()<<"____"<<endl;
                }

                //! extract the node coordinates from one of the two face mesh data sources
                cout<<"edgeMeshDataSourceBuildesClass::performOnBody()->____map intersection performed. NN: "<<theIntersectionMap.Extent()<<"____"<<endl;

                const occHandle(TColStd_HArray2OfReal) &theArrayOfCoords = new TColStd_HArray2OfReal(1,theIntersectionMap.Extent(),1,3);

                Standard_Real aCoordsBuf[3];
                TColStd_Array1OfReal theCoords(*aCoordsBuf,1,3);
                Standard_Integer NbNodes;
                MeshVS_EntityType theMeshEntityType;

                TColStd_MapIteratorOfPackedMapOfInteger anIter;
                anIter.Initialize(theIntersectionMap);

                QList<mesh::meshPoint> unsortedEdgePointList;
                QList<mesh::meshPoint> sortedEdgePointList;

                const occHandle(Ng_MeshVS_DataSource2D) &theMainMesh = occHandle(Ng_MeshVS_DataSource2D)::DownCast(myDB->ArrayOfMeshDS2D.value(bodyIndex));
                for(int i=1; anIter.More(); anIter.Next(), i++)
                {
                    int aKey = anIter.Key();

                    //! use the global surface mesh for retrieving the nodes of the edge
                    theMainMesh->GetGeom(aKey,Standard_False,theCoords,NbNodes,theMeshEntityType);

                    Standard_Real x = theCoords(1);
                    Standard_Real y = theCoords(2);
                    Standard_Real z = theCoords(3);

#ifdef VERBOSE
                    //!cout<<"____i: "<<i<<" key: "<<aKey<<" x= "<<x<<" y= "<<y<<" z= "<<z<<"____"<<endl;
#endif
                    theArrayOfCoords->SetValue(i,1,x);
                    theArrayOfCoords->SetValue(i,2,y);
                    theArrayOfCoords->SetValue(i,3,z);

                    mesh::meshPoint P(x,y,z,aKey);
                    unsortedEdgePointList<<P;
                }

                const TopoDS_Edge &theCurEdge = TopoDS::Edge(myDB->MapOfBodyTopologyMap.value(bodyIndex).edgeMap.FindKey(edgeNr));

                MeshTools::sortPointsOnEdge(unsortedEdgePointList,theCurEdge,sortedEdgePointList);

                //! at this stage also "theArrayOfCoords" can be ordered ... to do

                edgeNr_to_listOfPoints.insert(edgeNr,sortedEdgePointList);

                int theMeshOrder = myDB->ArrayOfMeshOrder.value(bodyIndex) == Property::meshOrder_First? 1:2;
                occHandle(Ng_MeshVS_DataSource1D) edgeMeshDS = new Ng_MeshVS_DataSource1D(theIntersectionMap,theArrayOfCoords,theMeshOrder);
                edgeMeshArray.SetValue(edgeNr,edgeMeshDS);
            }
        }
        return true;
    }
    else
    {
        cout<<"edgeMeshDataSourceBuildesClass::performOne->____The data source 2D for the body nr: "<<
              bodyIndex<<" has not been generated yet____"<<endl;
        return false;
    }
}
# endif

#ifdef NEW_EDGE_DS_BUILDER
//! ------------------------------------------------------------------------
//! An edge is defined as the intersection of two faces
//! Here the structure "shapePair" contains:
//!
//! 1) edgeNr = the number of the current edge (number in the body edgeMap)
//! 2) eln1 = the first face number (number in the body faceMap)
//! 3) eln2 = the second face number
//!
//! ------------------------------------------------------------------------
struct shapePair
{
    int edgeNr;
    int eln1;
    int eln2;
};

//! ------------------------------------------------------ //
//! function: constructor                                  //
//! details:                                               //
//! -------------------------------------------------------//
edgeMeshDataSourceBuildesClass::edgeMeshDataSourceBuildesClass(meshDataBase *theMeshDataBase, QObject *parent):
    QObject(parent), myDB(theMeshDataBase)
{
    cout<<"edgeMeshDataSourceBuildesClass::edgeMeshDataSourceBuildesClass()->____CONSTRUCTOR CALLED____"<<endl;
}

//! ------------------------------------------------------ //
//! function: perform                                      //
//! details:                                               //
//! -------------------------------------------------------//
void edgeMeshDataSourceBuildesClass::perform()
{
    cout<<"edgeMeshDataSourceBuildesClass::perform()->____function called____"<<endl;
    for(int bodyIndex = 1; bodyIndex <= myDB->bodyMap.size(); bodyIndex++)
    {
        const TopoDS_Shape &theShape = myDB->bodyMap.value(bodyIndex);
        int NbTopologyEdges = myDB->MapOfBodyTopologyMap.value(bodyIndex).edgeMap.Extent();

        NCollection_Array1<occHandle(Ng_MeshVS_DataSource1D)> edgeMeshArray(1,NbTopologyEdges);
        bool isDone = performOne(theShape, edgeMeshArray);

        if(isDone)
        {
            //! ---------------------------------------------
            //! put the edge mesh in the simulation database
            //! ---------------------------------------------
            for(int edgeNr = 1; edgeNr<=NbTopologyEdges; edgeNr++)
                myDB->ArrayOfMeshDSOnEdges.setValue(bodyIndex, edgeNr, edgeMeshArray.Value(edgeNr));
        }
    }
}

//! ------------------------------------------------------ //
//! function: perform on a single body                     //
//! details:                                               //
//! -------------------------------------------------------//
bool edgeMeshDataSourceBuildesClass::performOne(const TopoDS_Shape &theShape,
                                                NCollection_Array1<opencascade::handle<Ng_MeshVS_DataSource1D>> &edgeMeshArray)
{
    int bodyIndex = myDB->bodyMap.key(theShape);
    cout<<"edgeMeshDataSourceBuildesClass::performOne()->____working on body N. "<<bodyIndex<<"____"<<endl;

    const occHandle(Ng_MeshVS_DataSource2D) &surfaceMesh = occHandle(Ng_MeshVS_DataSource2D)::DownCast(myDB->ArrayOfMeshDS2D.getValue(bodyIndex));
    if(!surfaceMesh.IsNull())
    {
        TopTools_IndexedDataMapOfShapeListOfShape edgeFaceMap;
        TopExp::MapShapesAndAncestors(theShape,TopAbs_EDGE,TopAbs_FACE,edgeFaceMap);

        //! ------------------------------
        //! obtain the face and edge maps
        //! ------------------------------
        TopTools_IndexedMapOfShape edgeMap = myDB->MapOfBodyTopologyMap.value(bodyIndex).edgeMap;
        TopTools_IndexedMapOfShape faceMap = myDB->MapOfBodyTopologyMap.value(bodyIndex).faceMap;

        //! -------------------
        //! scan all the edges
        //! -------------------
        for(int edgeNr=1; edgeNr<=edgeMap.Extent(); edgeNr++)
        {
            //! -----------------
            //! the current edge
            //! -----------------
            cout<<"____working on edgeNr: "<<edgeNr<<"____"<<endl;
            const TopoDS_Edge &theCurEdge = TopoDS::Edge(edgeMap.FindKey(edgeNr));

            //! ---------------------------------
            //! check if the edge is degenerated
            //! ---------------------------------
            bool isDegenerated = BRep_Tool::Degenerated(theCurEdge);
            if(!isDegenerated)
            {
                //! ----------------------------------------------------------
                //! find the two faces defining the edge through intersection
                //! the edge is defined by two faces
                //! ----------------------------------------------------------
                TopTools_ListOfShape theTwoFaces;
                edgeFaceMap.FindFromKey(theCurEdge, theTwoFaces);

                const TopoDS_Face &face1 = TopoDS::Face(theTwoFaces.First());
                const TopoDS_Face &face2 = TopoDS::Face(theTwoFaces.Last());

                bool isSeam = BRep_Tool::IsClosed(face1);   //! can use also face2

                if(!isSeam)
                {
                    int faceIndex1 = faceMap.FindIndex(face1);
                    int faceIndex2 = faceMap.FindIndex(face2);

                    const occHandle(Ng_MeshVS_DataSourceFace) &faceMeshDS1 = occHandle(Ng_MeshVS_DataSourceFace)::DownCast(myDB->ArrayOfMeshDSOnFaces.value((bodyIndex,faceIndex1));
                    const occHandle(Ng_MeshVS_DataSourceFace) &faceMeshDS2 = occHandle(Ng_MeshVS_DataSourceFace)::DownCast(myDB->ArrayOfMeshDSOnFaces.value((bodyIndex,faceIndex2));

                    if(!faceMeshDS1.IsNull() && !faceMeshDS2.IsNull())
                    {
                        TColStd_PackedMapOfInteger theMap1 = faceMeshDS1->GetAllNodes();
                        TColStd_PackedMapOfInteger theMap2 = faceMeshDS2->GetAllNodes();
                        TColStd_PackedMapOfInteger theIntersectionMap;
                        theIntersectionMap.Intersection(theMap1,theMap2);

                        if(!faceMeshDS1.IsNull())cout<<"+++++"<<theMap1.Extent()<<endl;
                        if(!faceMeshDS2.IsNull())cout<<"-----"<<theMap2.Extent()<<endl;

                        if(theIntersectionMap.Extent()>2)
                        {
                            //! ------------------------------------------------------------------------
                            //! extract the node coordinates from one of the two face mesh data sources
                            //! ------------------------------------------------------------------------
                            cout<<"____map intersection done. Map size: "<<theIntersectionMap.Extent()<<"____"<<endl;
                            const occHandle(TColStd_HArray2OfReal) &theArrayOfCoords = new TColStd_HArray2OfReal(1,theIntersectionMap.Extent(),1,3);

                            Standard_Real aCoordsBuf[3];
                            TColStd_Array1OfReal theCoords(*aCoordsBuf,1,3);
                            Standard_Integer NbNodes;
                            MeshVS_EntityType theMeshEntityType;

                            //! -------------------
                            //! the list of points
                            //! -------------------
                            QList<QList<double>> pointList;

                            //! ---------------------------------------------------------------
                            //! use one of the two faces for retrieving the points of the edge
                            //! ---------------------------------------------------------------
                            int i=0;
                            for(TColStd_MapIteratorOfPackedMapOfInteger it(theIntersectionMap); it.More(); it.Next())
                            {
                                int aKey = it.Key();
                                surfaceMesh->GetGeom(aKey,Standard_False,theCoords,NbNodes,theMeshEntityType);

                                Standard_Real x = theCoords(1);
                                Standard_Real y = theCoords(2);
                                Standard_Real z = theCoords(3);

                                cout<<"____i: "<<++i<<" key: "<<aKey<<" x= "<<x<<" y= "<<y<<" z= "<<z<<"____"<<endl;

                                QList<double> P;
                                P<<x<<y<<z;
                                pointList<<P;
                            }

                            //! --------------------------------------------------------------
                            //! sort the points along the edge using the curvilinear abscissa
                            //! --------------------------------------------------------------
                            QList<QList<double>> pointList_sorted = MeshTools::sortPointsOnEdge(pointList,theCurEdge);
                            edgeNr_to_listOfPoints.insert(edgeNr,pointList_sorted);

                            for(int i=0; i<pointList_sorted.length(); i++)
                            {
                                double x = pointList_sorted.at(i).at(0);
                                double y = pointList_sorted.at(i).at(1);
                                double z = pointList_sorted.at(i).at(2);

                                cout<<"____i: "<<i+1<<" x= "<<x<<" y= "<<y<<" z= "<<z<<"____"<<endl;
                            }

                            //! --------------------------------
                            //! create the edge mesh datasource
                            //! --------------------------------
                            int meshOrder = myDB->ArrayOfMeshOrder.value(bodyIndex) == Property::meshOrder_First? 1:2;
                            occHandle(Ng_MeshVS_DataSource1D) edgeMeshDS = new Ng_MeshVS_DataSource1D(theIntersectionMap,theArrayOfCoords,meshOrder);
                            edgeMeshArray.SetValue(edgeNr,edgeMeshDS);
                        }
                        else
                        {
                            //! ---------------------------
                            //! intersection map too small
                            //! insert a null mesh Handle
                            //! ---------------------------
                            cout<<"____intersection map too small for edge nr: "<<edgeNr<<"____"<<endl;
                            edgeMeshArray.SetValue(edgeNr,occHandle(Ng_MeshVS_DataSource1D)());
                        }
                    }
                    else
                    {
                        //! -------------------------------------------------------------
                        //! one of the two face mesh data sources (or both) is(are) null
                        //! insert a null mesh Handle
                        //! -------------------------------------------------------------
                        cout<<"____cannot generate the mesh data source for edge nr: "<<edgeNr<<"____"<<endl;
                        edgeMeshArray.SetValue(edgeNr,occHandle(Ng_MeshVS_DataSource1D)());
                    }
                }
                else
                {
                    //! -------------------
                    //! handle a seam edge
                    //! -------------------
                    cout<<"____The edge nr: "<<edgeNr<<" is seam____"<<endl;
                    edgeMeshArray.SetValue(edgeNr,occHandle(Ng_MeshVS_DataSource1D)());
                }
            }
            else
            {
                //! --------------------------
                //! handle a degenerate edge
                //! insert a null mesh Handle
                //! --------------------------
                cout<<"____The edge nr: "<<edgeNr<<" is degenerate____"<<endl;
                edgeMeshArray.SetValue(edgeNr,occHandle(Ng_MeshVS_DataSource1D)());
            }
        }
        return true;
    }
    else
    {
        cout<<"____the edge data sources cannot be generated____"<<endl;
        return true;
    }
}
#endif
