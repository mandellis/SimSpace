//! ----------------
//! custom includes
//! ----------------
#include "surfacemeshcutter.h"
#include <polygon.h>
#include <mesh.h>
#include <meshelementbycoords.h>

//! ----
//! OCC
//! ----
#include <Bnd_Box.hxx>
#include <BRepBndLib.hxx>
#include <BRep_Builder.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Shape.hxx>

surfaceMeshCutter::surfaceMeshCutter()
{
    ;
}

//! -------------------------
//! function: cutSurfaceMesh
//! details:
//! -------------------------
bool surfaceMeshCutter::cutSurfaceMesh(const occHandle(MeshVS_DataSource) &inputMeshDS,
                                       const QList<TopoDS_Shape> &listOfShapes,
                                       occHandle(MeshVS_DataSource) &cutMesh)
{
    //! ------------------------------------------------
    //! build a compound using the input shapes (faces)
    //! ------------------------------------------------
    BRep_Builder aBuilder;
    TopoDS_Compound aComp;
    aBuilder.MakeCompound(aComp);
    for(int i=0; i<listOfShapes.length(); i++) aBuilder.Add(aComp,listOfShapes.at(i));

    //! ---------------------------------------------
    //! compute the bounding box of the input shapes
    //! ---------------------------------------------
    Bnd_Box aBBX;
    BRepBndLib::Add(aComp, aBBX);
    double xmin, ymin, zmin, xmax, ymax, zmax;
    aBBX.Get(xmin,ymin,zmin,xmax,ymax,zmax);
    double dx = (xmax-xmin)*0.1;
    double dy = (ymax-ymin)*0.1;
    double dz = (zmax-zmin)*0.1;

    //! ---------------------------
    //! thin BBX along a direction
    //! ---------------------------
    if(dx<0.1) dx = 10.0;
    if(dy<0.1) dy = 10.0;
    if(dz<0.1) dz = 10.0;

    bool isDone = surfaceMeshCutter::cutSurfaceMesh(xmin-dx,ymin-dy,zmin-dz,xmax+dx,ymax+dy,zmax+dz,inputMeshDS,cutMesh);
    return isDone;
}

//! -------------------------
//! function: cutSurfaceMesh
//! details:
//! -------------------------
bool surfaceMeshCutter::cutSurfaceMesh(double xmin, double ymin, double zmin, double xmax, double ymax, double zmax,
                                       const occHandle(MeshVS_DataSource) &inputMeshDS,
                                       occHandle(MeshVS_DataSource) &cutMesh)
{
    //! -------------
    //! sanity check
    //! -------------
    if(inputMeshDS.IsNull()) return false;
    if(inputMeshDS->GetAllElements().Extent()<1) return false;
    if(inputMeshDS->GetAllNodes().Extent()<3) return false;

    //! -----------------------------------------
    //! the mesh elements resulting from the cut
    //! -----------------------------------------
    std::vector<meshElementByCoords> vecElements;

    //! --------------------------------------------
    //! iterate over the elements of the input mesh
    //! --------------------------------------------
    for(TColStd_MapIteratorOfPackedMapOfInteger it(inputMeshDS->GetAllElements()); it.More(); it.Next())
    {
        int globalElementID = it.Key();
        int nbuf[8], NbNodes;
        TColStd_Array1OfInteger nodeIDs(*nbuf,1,8);
        bool isDone = inputMeshDS->GetNodesByElement(globalElementID,nodeIDs,NbNodes);
        if(!isDone) continue;

        //! ---------------
        //! a mesh element
        //! ---------------
        meshElementByCoords aMeshElement;
        aMeshElement.ID = globalElementID;

        //! ------------------------------------------
        //! compute the center of the current element
        //! ------------------------------------------
        std::vector<polygon::Point> aPolygon;
        for(int k=1; k<=NbNodes; k++)
        {
            int globalNodeID = nodeIDs(k);
            MeshVS_EntityType aType;
            int NbNodes1;
            double buf[3];
            TColStd_Array1OfReal nodeCoords(*buf,1,3);
            inputMeshDS->GetGeom(globalNodeID,false,nodeCoords,NbNodes1,aType);
            double x = nodeCoords(1);
            double y = nodeCoords(2);
            double z = nodeCoords(3);
            polygon::Point P(x,y,z);
            aPolygon.push_back(P);
            mesh::meshPoint aMeshPoint(x,y,z,globalNodeID);
            aMeshElement.pointList<<aMeshPoint;
        }

        //! ----------------------
        //! center of the polygon
        //! ----------------------
        polygon::Point center;
        double area;
        isDone = polygon::getPolygonCenter(aPolygon,center,area);
        if(!isDone) continue;

        //! -------------------------------------------------------------
        //! check if the polygon center is inside the input bounding box
        //! -------------------------------------------------------------
        double xc = center.x; double yc = center.y; double zc = center.z;
        if(xc>=xmin && xc<=xmax && yc>=ymin && yc<=ymax && zc>=zmin && zc<=zmax)
        {
            switch(NbNodes)
            {
            case 3: aMeshElement.type=TRIG; break;
            case 4: aMeshElement.type=QUAD; break;
            case 6: aMeshElement.type=TRIG6; break;
            case 8: aMeshElement.type=QUAD8; break;
            }
            vecElements.push_back(aMeshElement);
        }
    }

    //! -------------
    //! sanity check
    //! -------------
    cutMesh = new Ng_MeshVS_DataSourceFace(vecElements,false,false);
    if(cutMesh.IsNull()) return false;
    return true;
}
