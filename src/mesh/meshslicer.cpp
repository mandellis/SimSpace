//! ----------------
//! custom includes
//! ----------------
#include <meshslicer.h>
#include <polygon.h>
#include <meshelementbycoords.h>
#include <mesh.h>

//! ----
//! OCC
//! ----
#include <TColStd_MapIteratorOfPackedMapOfInteger.hxx>
#include <TColStd_Array1OfInteger.hxx>
#include <MeshVS_EntityType.hxx>

//! ----
//! C++
//! ----
#include <vector>
#include <iostream>
using namespace std;

//! ---
//! Qt
//! ---
#include <QList>

#ifdef TEST_SLICER
//! -----------------------------------------------------------------------
//! igl - do not move this header this header from here, neither as a joke
//! -----------------------------------------------------------------------
#include <libigl/include/igl/copyleft/cgal/intersect_with_half_space.h>
#endif

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
meshSlicer::meshSlicer(const occHandle(Ng_MeshVS_DataSource3D) &aVolumeMeshDS)
{
    if(aVolumeMeshDS.IsNull()) return;
    this->setMeshDataSource(aVolumeMeshDS);
}

//! ----------------------------
//! function: setMeshDataSource
//! details:
//! ----------------------------
void meshSlicer::setMeshDataSource(const occHandle(Ng_MeshVS_DataSource3D) &aVolumeMeshDS)
{
    myVolumeMeshDS = aVolumeMeshDS;

#ifdef TEST_SLICER
    //! ------------
    //! mesh points
    //! ------------
    int NbPoints = aVolumeMeshDS->GetAllNodes().Extent();
    V.resize(NbPoints,3);
    for(int row=1; row<=NbPoints; row++)
    {
        const std::vector<double> &P = aVolumeMeshDS->getNodeCoordinates(row);
        for(int col=1; col<=3; col++) V(row-1,col-1) = P[col-1];
    }

    //! --------------
    //! mesh elements
    //! --------------
    int NbElements = aVolumeMeshDS->GetAllElements().Extent();
    T.resize(NbElements,4);
    int NbNodes = -1, buf[4];
    TColStd_Array1OfInteger nodeIDs(*buf,1,4);
    TColStd_IndexedMapOfInteger nodesMap = aVolumeMeshDS->myNodesMap;
    for(int i=1; i<=NbElements; i++)
    {
        int globalElementID = aVolumeMeshDS->myElementsMap.FindIndex(i);
        aVolumeMeshDS->GetNodesByElement(globalElementID,nodeIDs,NbNodes);

        int row = i-1;
        T(row,0) = nodesMap.FindIndex(nodeIDs(1))-1;
        T(row,1) = nodesMap.FindIndex(nodeIDs(2))-1;
        T(row,2) = nodesMap.FindIndex(nodeIDs(3))-1;
        T(row,3) = nodesMap.FindIndex(nodeIDs(4))-1;

        std::vector<double> P0 = aVolumeMeshDS->getNodeCoordinates(nodeIDs(1));
        std::vector<double> P1 = aVolumeMeshDS->getNodeCoordinates(nodeIDs(2));
        std::vector<double> P2 = aVolumeMeshDS->getNodeCoordinates(nodeIDs(3));
        std::vector<double> P3 = aVolumeMeshDS->getNodeCoordinates(nodeIDs(4));
    }
#endif
}

#ifndef TEST_SLICER
//! -------------------
//! function: perform
//! details:
//! -------------------
bool meshSlicer::perform(double a, double b, double c, double d, opencascade::handle<Ng_MeshVS_DataSource3D> &slicedVolumeMeshDS)
{
    cout<<"meshSlicer::perform()->____function called____"<<endl;
    if(myVolumeMeshDS.IsNull()) return false;

    QList<meshElementByCoords> interceptedVolumeElements;

    int NbNodes, buf[20];
    TColStd_Array1OfInteger nodeIDs(*buf,1,20);
    for(TColStd_MapIteratorOfPackedMapOfInteger it(myVolumeMeshDS->GetAllElements()); it.More(); it.Next())
    {
        //! -------------------------------------------
        //! the candidate mesh points for building the
        //! volume mesh element
        //! -------------------------------------------
        std::vector<mesh::meshPoint> candidateMeshPoints;

        meshElementByCoords aVolumeMeshElement;
        std::vector<polygon::Point> aCloudOfPoints;

        int globalElementID = it.Key();
        myVolumeMeshDS->GetNodesByElement(globalElementID,nodeIDs,NbNodes);
        for(int i=1; i<=NbNodes; i++)
        {
            int globalNodeID = nodeIDs(i);
            int localNodeID = myVolumeMeshDS->myNodesMap.FindIndex(globalNodeID);
            const std::vector<double> &pc = myVolumeMeshDS->getNodeCoordinates(localNodeID);
            candidateMeshPoints.push_back(mesh::meshPoint(pc[0],pc[1],pc[2],globalNodeID));
            aCloudOfPoints.push_back(polygon::Point(pc[0],pc[1],pc[2]));
        }
        bool intersect = polygon::testPolygonPlaneIntersection(aCloudOfPoints,a,b,c,d);
        if(!intersect) continue;

        //! ------------------------
        //! build a mesh element 3D
        //! ------------------------
        aVolumeMeshElement.ID = globalElementID;
        for(int i=0; i<NbNodes; i++) aVolumeMeshElement.pointList<<candidateMeshPoints[i];
        switch(NbNodes)
        {
        case 4: aVolumeMeshElement.type = TET; break;
        case 5: aVolumeMeshElement.type = PYRAM; break;
        case 6: aVolumeMeshElement.type = PRISM; break;
        case 8: aVolumeMeshElement.type = HEXA; break;
        }
        interceptedVolumeElements<<aVolumeMeshElement;
    }

    //! ---------------------------------
    //! call the volume mesh constructor
    //! ---------------------------------
    slicedVolumeMeshDS = new Ng_MeshVS_DataSource3D(interceptedVolumeElements);
    if(slicedVolumeMeshDS.IsNull()) return false;
    //cout<<"meshSlicer::perform()->____number of nodes: "<<slicedVolumeMeshDS->GetAllNodes().Extent()<<"____"<<endl;
    //cout<<"meshSlicer::perform()->____number of elements: "<<slicedVolumeMeshDS->GetAllElements().Extent()<<"____"<<endl;

    return true;
}

#else
void half_space_slice(Eigen::MatrixXd &vertices, Eigen::MatrixXi &faces, Eigen::Vector3d &point, Eigen::Vector3d &direction)
{
    Eigen::MatrixXd J;
    igl::copyleft::cgal::intersect_with_half_space(vertices, faces, point, direction, vertices, faces, J);
}

//! -------------------
//! function: perform
//! details:
//! -------------------
bool meshSlicer::perform(double a, double b, double c, double d, occHandle(Ng_MeshVS_DataSource3D) &slicedVolumeMeshDS)
{
    Eigen::Vector3d point(0,0,100);
    Eigen::Vector3d direction(a,b,c);
    Eigen::MatrixXd Vs = V;
    Eigen::MatrixXi Ts = T;
    half_space_slice(Vs,Ts,point,direction);

    slicedVolumeMeshDS = new Ng_MeshVS_DataSource3D(Vs,Ts);
    cout<<"meshSlicer::perform()->____nodes: "<<slicedVolumeMeshDS->GetAllNodes().Extent()<<"____"<<endl;
    cout<<"meshSlicer::perform()->____elements: "<<slicedVolumeMeshDS->GetAllElements().Extent()<<"____"<<endl;
    return true;
}
#endif
