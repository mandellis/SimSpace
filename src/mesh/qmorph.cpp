#include "qmorph.h"

#include <QList>

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
QMorph::QMorph(const TopoDS_Face &aFace, const occHandle(Ng_MeshVS_DataSourceFace) &trigMesh):
    myFace(aFace), myTrigMesh(trigMesh)
{
    if(myTrigMesh->myBoundarySegments.isEmpty()) myTrigMesh->computeFreeMeshSegments();
    if(myTrigMesh->myNodeToSegmentConnectivity.isEmpty()) myTrigMesh->computeNodeToSegmentConnectivity();

    //! ----------------------
    //! set the initial front
    //! ----------------------
    myFront = myTrigMesh->myBoundarySegments;
}

//! ----------
//! function:
//! details:
//! ----------
void QMorph::setMeshEdgeState()
{
    int NbFrontEdges = myFront.length();
    for(int n=0; n<NbFrontEdges; n++)
    {
        //Ng_MeshVS_DataSourceFace::meshSegment aMeshSegment = myFront.at(n);
        mesh::meshSegment aMeshSegment = myFront.at(n);
        QList<int> nodeIDs = aMeshSegment.nodeIDs;

        int leftNode = nodeIDs.at(0);

        int rightNode = nodeIDs.at(1);
    }
}
