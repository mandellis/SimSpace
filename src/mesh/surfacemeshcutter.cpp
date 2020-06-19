//! ----------------
//! custom includes
//! ----------------
#include "surfacemeshcutter.h"
#include <polygon.h>
#include <mesh.h>

//! ---
//! Qt
//! ---
#include <QList>


surfaceMeshCutter::surfaceMeshCutter()
{
    ;
}

//! -------------------------
//! function: cutSurfaceMesh
//! details:
//! -------------------------
#include <Bnd_Box.hxx>
#include <BRepBndLib.hxx>
#include <BRep_Builder.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Shape.hxx>
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
    cout<<"surfaceMeshCutter::cutSurfaceMesh()->____function called____"<<endl;
    cout<<"surfaceMeshCutter::cutSurfaceMesh()->____("<<xmin<<", "<<xmax<<")____"<<endl;
    cout<<"surfaceMeshCutter::cutSurfaceMesh()->____("<<ymin<<", "<<ymax<<")____"<<endl;
    cout<<"surfaceMeshCutter::cutSurfaceMesh()->____("<<zmin<<", "<<zmax<<")____"<<endl;

    //! -------------
    //! sanity check
    //! -------------
    if(inputMeshDS.IsNull()) return false;
    if(inputMeshDS->GetAllElements().Extent()<1) return false;
    if(inputMeshDS->GetAllNodes().Extent()<3) return false;

    //! ---------------------------------------------------
    //! the points and the elements resulting from the cut
    //! ---------------------------------------------------
    QList<mesh::meshPoint> points;
    QList<mesh::meshElement> elements;

    //! --------------------------------------------
    //! iterate over the elements of the input mesh
    //! --------------------------------------------
    for(TColStd_MapIteratorOfPackedMapOfInteger eMap(inputMeshDS->GetAllElements()); eMap.More(); eMap.Next())
    {
        int globalElementID = eMap.Key();
        int NbNodes;
        int nbuf[8];
        TColStd_Array1OfInteger nodeIDs(*nbuf,1,8);
        bool isDone = inputMeshDS->GetNodesByElement(globalElementID,nodeIDs,NbNodes);
        if(!isDone) continue;       // the nodes of the element cannot be retrieved

        //cout<<"\\-------------start triangle----------\\"<<endl;
        //cout<<"____(";
        //for(int i=1; i<NbNodes; i++) cout<<nodeIDs(i)<<", ";
        //cout<<nodeIDs(NbNodes)<<")____"<<endl;

        //! ----------------
        //! build a polygon
        //! ----------------
        std::vector<mesh::meshPoint> element;   // it will be used further on
        std::vector<polygon::Point> aPolygon;
        for(int k=1; k<=NbNodes; k++)
        {
            int globalNodeID = nodeIDs(k);
            MeshVS_EntityType type;
            int NbNodes1;
            double buf[3];
            TColStd_Array1OfReal nodeCoords(*buf,1,3);
            inputMeshDS->GetGeom(globalNodeID,false,nodeCoords,NbNodes1,type);
            double x = nodeCoords(1);
            double y = nodeCoords(2);
            double z = nodeCoords(3);
            //cout<<"surfaceMeshCutter::cutSurfaceMesh()->____("<<x<<", "<<y<<", "<<z<<")____"<<endl;
            polygon::Point P(x,y,z);
            aPolygon.push_back(P);

            mesh::meshPoint Pm(x,y,z,globalNodeID);     // it will be used further on
            element.push_back(Pm);                      // it will be used further on
        }
        //cout<<"\\-------------end triangle-----------\\"<<endl;
        //cout<<"\\-------------calculating center-----\\"<<endl;
        //! ----------------------
        //! center of the polygon
        //! ----------------------
        polygon::Point center;
        double area;
        isDone = polygon::getPolygonCenter(aPolygon,center,area);
        if(!isDone)
        {
            cout<<"surfaceMeshCutter::cutSurfaceMesh()->____the center of an element has not been found____"<<endl;
            continue;
        }

        //cout<<"____C("<<center.x<<", "<<center.y<<", "<<center.z<<")____"<<endl;
        //cout<<"\\-------------center calculated------\\"<<endl;

        //! -------------------------------------------------------------
        //! check if the polygon center is inside the input bounding box
        //! -------------------------------------------------------------
        double xc = center.x; double yc = center.y; double zc = center.z;
        if(xc>=xmin && xc<=xmax && yc>=ymin && yc<=ymax && zc>=zmin && zc<=zmax)
        {
            //! the type of the element
            ElemType eType;
            switch(NbNodes)
            {
            case 3: eType=TRIG; break;
            case 4: eType=QUAD; break;
            case 6: eType=TRIG6; break;
            case 8: eType=QUAD8; break;
            }
            std::vector<int> IDs;
            for(int j=0;j<NbNodes;j++) IDs.push_back(nodeIDs(j+1));

            mesh::meshElement aMeshElement(eType,globalElementID,IDs);
            if(!elements.contains(aMeshElement)) elements<<aMeshElement;

            for(int l=0; l<element.size(); l++)
            {
                const mesh::meshPoint &elementPoint = element.at(l);
                if(!points.contains(elementPoint)) points<<elementPoint;
            }
        }
    }

    //! -------------
    //! sanity check
    //! -------------
    int NbPoints = points.length();
    int NbElements = elements.length();
    if(NbPoints<3 || NbElements<1) return false;

    cutMesh = new Ng_MeshVS_DataSourceFace(points,elements);

    //! ----------------------------------------------
    //! handle internal errors within the constructor
    //! ----------------------------------------------
    if(cutMesh.IsNull()) return false;
    return true;
}
