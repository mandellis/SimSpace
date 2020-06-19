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

//! -------------------
//! function: perform
//! details:
//! -------------------
bool meshSlicer::perform(double a, double b, double c, double d, opencascade::handle<Ng_MeshVS_DataSource3D> &slicedVolumeMeshDS)
{
    cout<<"meshSlicer::perform()->____function called____"<<endl;
    if(myVolumeMeshDS.IsNull()) return false;

    QList<meshElementByCoords> interceptedVolumeElements;

    cout<<"meshSlicer::perform()->____iterating over the elements____"<<endl;

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
        //int n;
        bool intersect = polygon::testPolygonPlaneIntersection(aCloudOfPoints,a,b,c,d/*,n*/);
        //cout<<"meshSlicer::perform()->________element discarded____"<<endl;
        if(!intersect) continue;

        //cout<<"meshSlicer::perform()->____adding element: "<<globalElementID<<"____"<<endl;

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
    //cout<<"meshSlicer::perform()->____calling constructor____"<<endl;

    slicedVolumeMeshDS = new Ng_MeshVS_DataSource3D(interceptedVolumeElements);
    if(slicedVolumeMeshDS.IsNull()) return false;
    //cout<<"meshSlicer::perform()->____number of nodes: "<<slicedVolumeMeshDS->GetAllNodes().Extent()<<"____"<<endl;
    //cout<<"meshSlicer::perform()->____number of elements: "<<slicedVolumeMeshDS->GetAllElements().Extent()<<"____"<<endl;

    return true;
}
