#ifndef CUSTOMMESHER_H
#define CUSTOMMESHER_H

//! path to triangle.exe
#define TRIANGLE_PROGRAM_PATH "D:/Work/CustomMesher/FaceMesher/triangle.exe"

//! custom includes
#include "occhandle.h"
#include <meshdatabase.h>
#include <geometrytag.h>

//! OCC
#include <TopoDS_Solid.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Edge.hxx>

//! std
#include <vector>

//! Qt
#include <QObject>
#include <QList>

//! custom
#include "hash_c.h"
#include <mesh.h>

class QProcess;
class Geom_Surface;
class gp_Pnt;

//! -----------------------------
//! the class of the face mesher
//! -----------------------------
class CustomMesher: public QObject
{
    Q_OBJECT

public:

    //! ---------------------------
    //! definition of a mesh point
    //! ---------------------------
    struct meshPoint
    {
        double x,y,z;
        int ID;

        //! -------------
        //! constructors
        //! -------------
        meshPoint(){ x=y=z=-1; ID = -1; }
        meshPoint(double x0,double y0,double z0):x(x0),y(y0),z(z0),ID(-1){;}
        meshPoint(double x0,double y0,double z0,int ID0):x(x0),y(y0),z(z0),ID(ID0){;}

        //! -----------------
        //! copy constructor
        //! -----------------
        meshPoint(const meshPoint &aMeshPoint)
        {
            x = aMeshPoint.x;
            y = aMeshPoint.y;
            z = aMeshPoint.z;
            ID = aMeshPoint.ID;
        }

        //! -----------
        //! operator <
        //! -----------
        inline bool operator <(const meshPoint &rhs)
        {
            std::size_t seed_1 = 0, seed_2 = 0;
            hash_c<double>(seed_1,x); hash_c<double>(seed_1,y); hash_c<double>(seed_1,z);
            hash_c<double>(seed_2,rhs.x); hash_c<double>(seed_2,rhs.y); hash_c<double>(seed_2,rhs.z);
            if(seed_1<seed_2) return true; return false;
        }

        //! -----------
        //! operator =
        //! -----------
        inline meshPoint* operator = (const meshPoint &rhs)
        {
            x = rhs.x; y = rhs.y; z = rhs.z; ID = rhs.ID;
            return this;
        }

        //! ---------------------------------------
        //! operator ==
        //! acting over coordinates with tolerance
        //! ---------------------------------------
        inline bool operator == (const meshPoint &rhs)
        {
            const double tolerance = 1e-6;
            if(fabs(x-rhs.x)<=tolerance && fabs(y-rhs.y)<=tolerance && fabs(z-rhs.z)<=tolerance) return true;
            return false;
        }
    };

public:

    CustomMesher(const TopoDS_Shape &aShape, meshDataBase *mDB);

private:

    //! --------------------------------
    //! get the nodeID from coordinates
    //! --------------------------------
    bool getNodeID(const CustomMesher::meshPoint &aPoint, const QList<CustomMesher::meshPoint> &allMeshPoints, int &nodeID)
    {
        if(allMeshPoints.isEmpty()) return false;
        if(!allMeshPoints.contains(aPoint))
        {
            nodeID = -1;
            return false;
        }
        int index = allMeshPoints.indexOf(aPoint);
        const CustomMesher::meshPoint &P = allMeshPoints.at(index);
        nodeID = P.ID;
        return true;
    }

    //! -----------------------------------------------------------------------------------
    //! given a mesh point check if a list of mesh points contains a node with the same ID
    //! if yes return the first occurrence. To be optimized ... to dop ...
    //! -----------------------------------------------------------------------------------
    bool containsNodeID(const CustomMesher::meshPoint &aPoint, const QList<CustomMesher::meshPoint> &allMeshPoints, CustomMesher::meshPoint &firstOccurrencePoint)
    {
        if(allMeshPoints.isEmpty()) return false;
        int theNodeID = aPoint.ID;
        for(int i=0; i<allMeshPoints.length(); i++)
        {
            const CustomMesher::meshPoint &curMeshPoint = allMeshPoints.at(i);
            int currentNodeID = curMeshPoint.ID;
            //cout<<"containsNodeID()->____comparing: "<<currentNodeID<<" with "<<theNodeID<<"____"<<endl;
            if(theNodeID == currentNodeID)
            {
                firstOccurrencePoint = curMeshPoint;
                return true;
            }
        }
        return false;
    }

private:

    //! -------------------
    //! the mesh data base
    //! -------------------
    meshDataBase *myMeshDB;

    //! -----------------------
    //! the shape to be meshed
    //! -----------------------
    TopoDS_Shape myShape;

    //! ---------
    //! edgeMesh
    //! ---------
    QMap<GeometryTag, std::vector<CustomMesher::meshPoint>> myEdgeMesh;

    //! -----------------------------
    //! the current geometry surface
    //! -----------------------------
    occHandle(Geom_Surface) myCurGeomSurface;

    //! ----------------------------------------
    //! the points describing the face boundary
    //! ----------------------------------------
    std::vector<CustomMesher::meshPoint> myBoundaryPoints;

    //! ----------------------------
    //! QProcess for "Triangle.exe"
    //! ----------------------------
    QProcess *myTriangleMesher;

    //! ------------------
    //! working directory
    //! ------------------
    QString myWorkingDir;

    //! ---------------------------------------------------------
    //! project the points of the boundary onto the (u, v) plane
    //! ---------------------------------------------------------
    QList<CustomMesher::meshPoint> getFaceUVBoundary(const TopoDS_Face &aFace, QMap<GeometryTag, std::vector<meshPoint> > &mapOfEdgeMeshUV);

    //! -----------------
    //! helper functions
    //! -----------------
    double faceArea(const TopoDS_Face &aFace);

    //! -----------------------------------------------------------
    //! get the face geometry tag. It acts onto the private member
    //! -----------------------------------------------------------
    GeometryTag faceTag(const TopoDS_Face &aFace);

    //! --------------------------
    //! get the edge geometry tag
    //! --------------------------
    GeometryTag edgeTag(const TopoDS_Edge &anEdge);

public:

    //! ----------------------
    //! get working directory
    //! ----------------------
    const QString getWorkingDir() const { return myWorkingDir; }

    //! -------------------------------------------------
    //! programmatically set the boundary through points
    //! -------------------------------------------------
    void setBoundary(std::vector<CustomMesher::meshPoint> &boundaryPoints) { myBoundaryPoints = boundaryPoints; }

    //! ----------------------
    //! build the planar mesh
    //! ----------------------
    occHandle(Ng_MeshVS_DataSourceFace) buildFaceMesh(const TopoDS_Face &aFace);

    //! ---------------
    //! mesh all edges
    //! ---------------
    QMap<GeometryTag, std::vector<CustomMesher::meshPoint> > meshAllEdges();

    //! -------------------
    //! mesh all the faces
    //! -------------------
    QMap<GeometryTag, occHandle(Ng_MeshVS_DataSourceFace)> meshAllFaces();


private slots:

    //! ------------------------
    //! QProcess shell messages
    //! ------------------------
    void redirectOutput();
};

#endif // FACEMESHER_H
