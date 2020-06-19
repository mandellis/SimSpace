#include "custommesher.h"

//! Qt
#include <QDir>
#include <QProcess>
#include <QList>

//! OCC
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <gp_Pnt.hxx>
#include <Geom_Surface.hxx>
#include <GeomAdaptor_Surface.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Shape.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopExp.hxx>
#include <meshdatabase.h>
#include <BRep_Tool.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <GeomAdaptor_Curve.hxx>
#include <CPnts_AbscissaPoint.hxx>
#include <TopExp_Explorer.hxx>
#include <TopAbs.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Wire.hxx>
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>

//! Eigen
#include <Eigen/Dense>

using namespace std;

//! ---------------------------------------------------------
//! function: constructor
//! details: initialize with the shape and the mesh database
//! ---------------------------------------------------------
CustomMesher::CustomMesher(const TopoDS_Shape &aShape, meshDataBase *mDB):
    myShape(aShape),
    myMeshDB(mDB)
{
    cout<<"CustomMesher::CustomMesher()->____function called____"<<endl;

    //! -------------------
    //! create the process
    //! -------------------
    myTriangleMesher = new QProcess(this);

    //! -------------------
    //! define the program
    //! -------------------
    QString program = QString(TRIANGLE_PROGRAM_PATH);
    myTriangleMesher->setProgram(program);

    //QString currentPath = QDir::currentPath();
    //QString tmpDir = currentPath+"/tmp";
    //cout<<"____"<<tmpDir.toStdString()<<"____"<<endl;
}

//! --------------------------------------------------------------
//! function: meshAllEdges
//! details:  the discretization of edges is done in the 3D space
//! --------------------------------------------------------------
QMap<GeometryTag,std::vector<CustomMesher::meshPoint>> CustomMesher::meshAllEdges()
{
    cout<<"CustomMesher::meshAllEdges()->____function called____"<<endl;

    //!----------------
    //! the edges mesh
    //! ---------------
    QMap<GeometryTag,std::vector<CustomMesher::meshPoint>> mapOfEdgeMesh;

    //! ---------------
    //! all the points
    //! ---------------
    QList<CustomMesher::meshPoint> allEdgeMeshPoints;

    //! ----------------------------
    //! scan the edges of the shape
    //! ----------------------------
    int nodeID = 0;
    for(TopExp_Explorer edgeExp(myShape,TopAbs_EDGE); edgeExp.More(); edgeExp.Next())
    {
        //! -----------------
        //! the current edge
        //! -----------------
        const TopoDS_Edge &edge = TopoDS::Edge(edgeExp.Current());
        if(BRep_Tool::Degenerated(edge))
        {
            cout<<"CustomMesher::meshEdges()->____a degenerate edge has been found: jumping over it____"<<endl;
            continue;
        }

        //! ----------------------------
        //! build the edge geometry tag
        //! ----------------------------
        const GeometryTag &edgeTag = this->edgeTag(edge);

        cout<<"CustomMesher::meshEdges()->____meshing the edge with tag: ("<<edgeTag.parentShapeNr<<","<<edgeTag.subTopNr<<")____"<<endl;

        BRepAdaptor_Curve BRepAdaptor(edge);
        GeomAdaptor_Curve curve = BRepAdaptor.Curve();

        //! -------------------
        //! length of the edge
        //! -------------------
        CPnts_AbscissaPoint CP;
        CP.Init(curve);
        double edgeLength = CP.Length(curve);

        //! --------------------
        //! temporary criterion
        //! --------------------
        const int NbDivisions = 5;
        double deltaL = edgeLength/NbDivisions;

        //! ---------------------
        //! the mesh of the edge
        //! ---------------------
        std::vector<CustomMesher::meshPoint> points;

        //! -------------------------
        //! retrieve the first point
        //! -------------------------
        double s_first = BRepAdaptor.FirstParameter();
        gp_Pnt P_in = BRepAdaptor.Value(s_first);
        CustomMesher::meshPoint P_first;
        P_first.x = P_in.X();
        P_first.y = P_in.Y();
        P_first.z = P_in.Z();

        int nodeID_old = -1;
        bool hasID = this->getNodeID(P_first,allEdgeMeshPoints,nodeID_old);
        if(hasID)
        {
            //cout<<"____change node ID____"<<endl;
            P_first.ID = nodeID_old;
        }
        else
        {
            nodeID++;
            P_first.ID = nodeID;
        }
        allEdgeMeshPoints<<P_first;
        points.push_back(P_first);
        cout<<P_first.ID<<"\t"<<P_first.x<<"\t"<<P_first.y<<"\t"<<P_first.z<<endl;

        //! --------------------------
        //! retrieve the inner points
        //! --------------------------
        double s_old = s_first;
        const double tol = 1e-6;
        for(int n=1; n<=NbDivisions-1; n++)
        {
            CP.Perform(deltaL,s_old,tol);
            double s = CP.Parameter();
            gp_Pnt P_s = BRepAdaptor.Value(s);
            CustomMesher::meshPoint P;
            P.x = P_s.X(); P.y = P_s.Y(); P.z = P_s.Z();

            hasID = this->getNodeID(P,allEdgeMeshPoints,nodeID_old);
            if(hasID)
            {
                //cout<<"____change node ID____"<<endl;
                P.ID = nodeID_old;
            }
            else
            {
                nodeID++;
                P.ID = nodeID;
            }
            allEdgeMeshPoints<<P;
            points.push_back(P);
            s_old = s;
            cout<<P.ID<<"\t"<<P.x<<"\t"<<P.y<<"\t"<<P.z<<endl;
        }

        //! ------------------------
        //! retrieve the last point
        //! ------------------------
        double s_last = BRepAdaptor.LastParameter();
        gp_Pnt P_fin = BRepAdaptor.Value(s_last);
        CustomMesher::meshPoint P_last;
        P_last.x = P_fin.X(); P_last.y = P_fin.Y(); P_last.z = P_fin.Z();

        hasID = this->getNodeID(P_last,allEdgeMeshPoints,nodeID_old);
        if(hasID)
        {
            //cout<<"____change node ID____"<<endl;
            P_last.ID = nodeID_old;
        }
        else
        {
            nodeID++;
            P_last.ID = nodeID;
        }
        allEdgeMeshPoints<<P_last;
        points.push_back(P_last);

        cout<<P_last.ID<<"\t"<<P_last.x<<"\t"<<P_last.y<<"\t"<<P_last.z<<endl;

        mapOfEdgeMesh.insert(edgeTag,points);
    }

    //! -----------
    //! diagnostic
    //! -----------
    FILE *f = fopen("D:/edgeMeshPoints.txt","w");
    fprintf(f,"#bodyIndex\tedgeIndex\tID\tx\ty\tz\n");
    for(QMap<GeometryTag,std::vector<CustomMesher::meshPoint>>::iterator it= mapOfEdgeMesh.begin(); it!=mapOfEdgeMesh.end(); ++it)
    {
        const GeometryTag &aTag = it.key();
        const std::vector<CustomMesher::meshPoint> &anEdgeMesh = it.value();
        for(int i=0; i<anEdgeMesh.size(); i++)
        {
            const CustomMesher::meshPoint &aPoint = anEdgeMesh.at(i);
            fprintf(f,"%d\t%d\t%d\t%lf\t%lf\t%lf\n",aTag.parentShapeNr,aTag.subTopNr,aPoint.ID,aPoint.x,aPoint.y,aPoint.z);
        }
    }
    fclose(f);
    //! ---------------
    //! end diagnostic
    //! ---------------

    myEdgeMesh = mapOfEdgeMesh;
    return mapOfEdgeMesh;
}

//! -----------------------
//! function: meshAllFaces
//! details:
//! -----------------------
QMap<GeometryTag,occHandle(Ng_MeshVS_DataSourceFace)> CustomMesher::meshAllFaces()
{
    cout<<"CustomMesher::meshAllFaces()->____function called____"<<endl;

    if(myEdgeMesh.isEmpty())
    {
        cout<<"CustomMesher::meshAllFaces()->____meshing edges____"<<endl;
        this->meshAllEdges();
    }

    //! --------------------------------
    //! scan all the faces of the shape
    //! --------------------------------
    QMap<GeometryTag,occHandle(Ng_MeshVS_DataSourceFace)> allFacesMesh;
    for(TopExp_Explorer anExp(myShape,TopAbs_FACE); anExp.More(); anExp.Next())
    {
        //! -------------------
        //! build the face tag
        //! -------------------
        const TopoDS_Face &aFace = TopoDS::Face(anExp.Current());
        const GeometryTag &aFaceTag = faceTag(aFace);

        cout<<"CustomMesher::meshAllFaces()->____working on face with tag ("<<aFaceTag.parentShapeNr<<", "<<aFaceTag.subTopNr<<")____"<<endl;

        //! --------------------------------------------------------
        //! actually mesh the current face: call a private function
        //! --------------------------------------------------------
        occHandle(Ng_MeshVS_DataSourceFace) curFaceMesh = this->buildFaceMesh(aFace);
        allFacesMesh.insert(aFaceTag,curFaceMesh);
    }
    return allFacesMesh;
}

//! --------------------------------------------------------------------------------------
//! function: getFaceUVBoundary
//! details:  perform the projection of each face boundary point onto the (u,v) plane
//!           Algo: scan the edges od the face. Retrieve the edge tag. Enter the
//!           map QMap<int,std::vector<CustomMesher::meshPoints> myEdgesMesh and retrieve
//!           the mesh points of the current edge. Project each point on to the surface
//!           and retrieve the (u,v) coordinates of the point. Pile up the (u,v) pair.
//!           Note: we are taking a point "P" of a previously generated edge mesh, and
//!           projecting it onto the surface of the face, getting "P'".
//!           "P" and "P'" points could have slightly different (x,y,z), because of the
//!           projection operation, so it is better to update, after projection, the
//!           coordinates of "P" using the coordinates of "P'"
//! --------------------------------------------------------------------------------------
QList<CustomMesher::meshPoint> CustomMesher::getFaceUVBoundary(const TopoDS_Face &aFace,
                                                               QMap<GeometryTag,std::vector<meshPoint>> &mapOfEdgeMeshUV)
{
    cout<<"CustomMesher::getFaceUVBoundary()->____function called____"<<endl;

    //! ---------
    //! adaptors
    //! ---------
    BRepAdaptor_Surface adaptor(aFace,true);
    const GeomAdaptor_Surface &s_adaptor = adaptor.Surface();

    //! --------------------
    //! surface of the face
    //! --------------------
    myCurGeomSurface = s_adaptor.Surface();

    //! --------------------------------------------
    //! number of points defining the face boundary
    //! --------------------------------------------
    int NbPoints = 0;

    //! -------------------------------------------------------------
    //! points of the face wire (face boundary - remember: no holes)
    //! -------------------------------------------------------------
    QList<CustomMesher::meshPoint> uniqueBoundaryMeshPoints;

    QList<CustomMesher::meshPoint> allEdgeMeshPointsUV;
    for(TopExp_Explorer edgeExp(aFace,TopAbs_EDGE); edgeExp.More(); edgeExp.Next())
    {
        const TopoDS_Edge &curEdge = TopoDS::Edge(edgeExp.Current());
        if(BRep_Tool::Degenerated(curEdge)) continue;

        const GeometryTag &curEdgeTag = this->edgeTag(curEdge);

        //cout<<"CustomMesher::getFaceUVBoundary()->____working on edge tag("<<curEdgeTag.parentShapeNr<<", "<<curEdgeTag.subTopNr<<")____"<<endl;

        //! -------------------------------------------------------------------
        //! the mesh of the edges in the 3D space (should be generated before)
        //! -------------------------------------------------------------------
        std::vector<CustomMesher::meshPoint> curEdgeMesh = myEdgeMesh.value(curEdgeTag);

        //! ----------------------------------------
        //! the mesh of the edge in the (u,v) space
        //! ----------------------------------------
        std::vector<CustomMesher::meshPoint> curEdgeMeshUV;

        //! ----------------------------------------------------------
        //! this snippet retrieve the mesh points of the current face
        //! boundary and project each point onto the (u,v) plane
        //! ----------------------------------------------------------
        for(int i=0; i<curEdgeMesh.size(); i++)
        {
            //! -----------------------------------
            //! the current point of the edge mesh
            //! -----------------------------------
            const CustomMesher::meshPoint &curPoint = curEdgeMesh.at(i);
            int index = curPoint.ID;

            //! ---------------------------------------------------------------
            //! define the OCC gp_Pnt and project onto the surface of the face
            //! ---------------------------------------------------------------
            gp_Pnt P(curPoint.x,curPoint.y,curPoint.z);

            GeomAPI_ProjectPointOnSurf aProjector;
            aProjector.Init(P, myCurGeomSurface);

            //! ----------------------
            //! get the closest point
            //! ----------------------
            double u,v;
            aProjector.LowerDistanceParameters(u,v);

            //! ----------------------------------------------
            //! the projected point, defined with the same ID
            //! ----------------------------------------------
            meshPoint aMeshPointUV(u,v,0,index);
            meshPoint firstOccurrencePoint;

            //cout<<"____index: "<<index<<"____"<<endl;

            //! -----------------------------------------------
            //! this will (probably) merge the common vertices
            //! -----------------------------------------------
            bool contains = this->containsNodeID(aMeshPointUV,allEdgeMeshPointsUV,firstOccurrencePoint);
            if(contains)
            {
                aMeshPointUV.x = firstOccurrencePoint.x;
                aMeshPointUV.y = firstOccurrencePoint.y;
                aMeshPointUV.z = firstOccurrencePoint.z; // "0" by definition
            }
            else
            {
                NbPoints++;
            }
            allEdgeMeshPointsUV<<aMeshPointUV;
            curEdgeMeshUV.push_back(aMeshPointUV);

            CustomMesher::meshPoint dummy;
            if(!this->containsNodeID(aMeshPointUV,uniqueBoundaryMeshPoints,dummy))
            {
                uniqueBoundaryMeshPoints<<aMeshPointUV;
            }

            //! ---------------------------------------------------
            //! the nearest (hopefully the unique) projected point
            //! ---------------------------------------------------
            const gp_Pnt &theProjectedPoint = aProjector.NearestPoint();

            //! ----------------------------------------------------
            //! update the edge point using the coordinates
            //! of the projection. Let the ID unchanged.
            //! The point ID is the same with slightly shifted
            //! coordinated, due to the projection onto the surface
            //! ----------------------------------------------------
            CustomMesher::meshPoint updatedPoint(curPoint);
            updatedPoint.x = theProjectedPoint.X();
            updatedPoint.y = theProjectedPoint.Y();
            updatedPoint.z = theProjectedPoint.Z();

            //! ---------------------------------------------
            //! actually change the coordinates of the point
            //! ---------------------------------------------
            curEdgeMesh.at(i) = updatedPoint;
        }
        myEdgeMesh.insert(curEdgeTag,curEdgeMesh);
        mapOfEdgeMeshUV.insert(curEdgeTag,curEdgeMeshUV);
    }
    cout<<"\\--------------------------------------"<<endl;
    cout<<"\\ Nb points: "<<NbPoints<<endl;
    cout<<"\\ unique points Nb: "<<uniqueBoundaryMeshPoints.length()<<endl;
    cout<<"\\--------------------------------------"<<endl;

    return uniqueBoundaryMeshPoints;
}

//! --------------------------------------------------------------
//! function: buildFaceMesh
//! details:  Planar Straight Line Graph (PSLG) and mesh the face
//!           The operation is performed on disk, using Triangle
//! --------------------------------------------------------------
occHandle(Ng_MeshVS_DataSourceFace) CustomMesher::buildFaceMesh(const TopoDS_Face &aFace)
{
    cout<<"CustomMesher::buildFaceMesh()->____function called____"<<endl;

    //! -----------------
    //! analyze the face
    //! -----------------
    //QList<TopoDS_Wire> wireList;
    //for(TopExp_Explorer anExp(aFace,TopAbs_WIRE); anExp.More(); anExp.Next())
    //{
    //    const TopoDS_Wire &aWire = TopoDS::Wire(anExp.Current());
    //    wireList<<aWire;
    //}
    //int NbHoles = wireList.length()-1;
    //cout<<"CustomMesher::buildPlanarMesh()->____number of holes: "<<NbHoles<<"____"<<endl;

    //! ---------------------------
    //! create a support directory
    //! ---------------------------
    QString currentPath = QDir::currentPath();
    QString tmpDir = currentPath+"/tmp";

    //! --------------------------------------
    //! check if the directory already exists
    //! --------------------------------------
    QDir aDir;
    aDir.setPath(tmpDir);
    if(!aDir.exists()) aDir.mkdir(tmpDir);

    myWorkingDir = tmpDir;
    cout<<"CustomMesher::buildFaceMesh()->____working directory: "<<myWorkingDir.toStdString()<<"____"<<endl;

    //! ------------------
    //! build the contour
    //! ------------------
    char polyFileName[256];
    QString QPolyFileName = tmpDir+"/"+"aFace.poly";
    sprintf(polyFileName,"%s",QPolyFileName.toStdString().c_str());

    //! --------------------
    //! face in (u,v) plane
    //! --------------------
    cout<<"CustomMesher::buildFaceMesh()->____retrieving the boundary____"<<endl;

    QMap<GeometryTag,std::vector<meshPoint>> mapOfEdgeMeshUV;
    const QList<CustomMesher::meshPoint> &uniqueBoundaryPoints = this->getFaceUVBoundary(aFace,mapOfEdgeMeshUV);
    int NbPoints = uniqueBoundaryPoints.length();

    cout<<"____number of face wire points NbPoints: "<<NbPoints<<"____"<<endl;

    //! -------------------------------------------------------------------------------------------------------
    //! First line: <# of vertices> <dimension (must be 2)> <# of attributes> <# of boundary markers (0 or 1)>
    //! Following lines: <vertex #> <x> <y> [attributes] [boundary marker]
    //! One line: <# of segments> <# of boundary markers (0 or 1)>
    //! Following lines: <segment #> <endpoint> <endpoint> [boundary marker]
    //! One line: <# of holes>
    //! Following lines: <hole #> <x> <y>
    //! Optional line: <# of regional attributes and/or area constraints>
    //! Optional following lines: <region #> <x> <y> <attribute> <maximum area>
    //! -------------------------------------------------------------------------------------------------------
    FILE *polyFile = fopen(polyFileName,"w");

    //! -----------------------
    //! write the node section
    //! -----------------------    
    fprintf(polyFile,"#node section\n");
    fprintf(polyFile,"%d\t2\t0\t0\n",NbPoints);
    int i=0;
    for(QList<CustomMesher::meshPoint>::const_iterator it = uniqueBoundaryPoints.begin(); it!=uniqueBoundaryPoints.cend(); ++it)
    {
        const CustomMesher::meshPoint &P = *it;
        double u = P.x;
        double v = P.y;
        fprintf(polyFile,"%d\t%.9e\t%.9e\n",i+1,u,v);
        i++;
    }
    cout<<"____node section written____"<<endl;

    //! -------------------------------------
    //! write the segment section
    //! for the moment holes are not handled
    //! -------------------------------------
    int NbSegments = NbPoints;
    fprintf(polyFile,"#segment section\n");
    fprintf(polyFile,"%d\t%d\n",NbSegments,0);
    for(QMap<GeometryTag,std::vector<CustomMesher::meshPoint>>::iterator it = mapOfEdgeMeshUV.begin(); it!= mapOfEdgeMeshUV.end(); ++it)
    {
        //const GeometryTag &tag = it.key();
        const std::vector<CustomMesher::meshPoint> &edgeMesh = it.value();
        for(int i=0; i<edgeMesh.size()-1; i++)
        {
            const CustomMesher::meshPoint &firstPoint = edgeMesh.at(i);
            const CustomMesher::meshPoint &secondPoint = edgeMesh.at(i+1);
            int first = firstPoint.ID;
            int second = secondPoint.ID;
            fprintf(polyFile,"%d\t%d\t%d\n",i+1,first,second);
        }
    }
    cout<<"____segment section written____"<<endl;

    //! ----------------------------------------------------
    //! write the holes section: not handled for the moment
    //! ----------------------------------------------------
    fprintf(polyFile,"#holes section\n");
    fprintf(polyFile,"%d\n",0);

    fclose(polyFile);

    //! --------------------------------------------------
    //! set the working directory for the triangle mesher
    //! --------------------------------------------------
    myTriangleMesher->setWorkingDirectory(tmpDir);

    //! ----------
    //! arguments
    //! ----------
    QStringList arguments;
    //QString meshParametersString = QString("-pYq20.0a%1").arg(this->faceArea(aFace)/200.0);
    //arguments<<meshParametersString<<QPolyFileName;
    arguments<<QString("-pYq12.5a10.0")<<QPolyFileName;
    myTriangleMesher->setArguments(arguments);

    //! -----------------------------------
    //! setup a connection before starting
    //! -----------------------------------
    disconnect(myTriangleMesher,SIGNAL(readyReadStandardOutput()),this,SLOT(redirectOutput()));
    connect(myTriangleMesher,SIGNAL(readyReadStandardOutput()),this,SLOT(redirectOutput()));

    //! -----------------
    //! start the mesher
    //! -----------------
    myTriangleMesher->start();
    myTriangleMesher->waitForFinished(-1);

    int exitCode = myTriangleMesher->exitCode();
    cout<<"____\"Triangle\" exit code: "<<exitCode<<"____"<<endl;

    //! ------------------------
    //! read the generated mesh
    //! ------------------------
    occHandle(Ng_MeshVS_DataSourceFace) theFaceMeshDS;
    if(exitCode!=0) return theFaceMeshDS;       // "Triangle" mesher failure

    char tmpName[256];
    QString nodeFileName = tmpDir+"/"+"aFace.1.node";
    sprintf(tmpName,"%s",nodeFileName.toStdString().c_str());
    FILE *nodeFile = fopen(tmpName,"r");
    cout<<"opening the node file: "<<tmpName<<"____"<<endl;
    if(nodeFile==NULL) cout<<"____the file is NULL____"<<endl;

    QString faceFileName = tmpDir+"/"+"aFace.1.ele";
    sprintf(tmpName,"%s",faceFileName.toStdString().c_str());
    FILE *faceFile = fopen(tmpName,"r");
    cout<<"opening the face file: "<<tmpName<<"____"<<endl;
    if(faceFile==NULL) cout<<"____the file is NULL____"<<endl;

    //! --------------------
    //! read the mesh nodes
    //! --------------------
    Eigen::MatrixXd V;

    int nodeID,NbMeshNodes,dim,NbAttributes,bm;
    double u,v;
    fscanf(nodeFile,"%d%d%d%d",&NbMeshNodes,&dim,&NbAttributes,&bm);
    V.resize(NbMeshNodes,2);

    for(int i=1; i<=NbMeshNodes; i++)
    {
        fscanf(nodeFile,"%d%lf%lf%d",&nodeID,&u,&v,&bm);
        V(i-1,0) = u;
        V(i-1,1) = v;
        //cout<<"____(u, v) = ("<<u<<", "<<v<<")____"<<endl;
    }
    fclose(nodeFile);

    //! ---------------------
    //! reading the elements
    //! ---------------------
    Eigen::MatrixXi F;

    int elementID,NbTriangles,a,i1,i2,i3;
    fscanf(faceFile,"%d\t%d\t%d\n",&NbTriangles,&a,&NbAttributes);
    F.resize(NbTriangles,3);

    for(int i=1; i<=NbTriangles; i++)
    {
        fscanf(faceFile,"%d%d%d%d",&elementID,&i1,&i2,&i3);
        F(i-1,0) = i1; F(i-1,1) = i2; F(i-1,2) = i3;
        //cout<<"____(i1, i2, i3) = ("<<i1<<", "<<i2<<", "<<i3<<")____"<<endl;
    }
    fclose(faceFile);

    //! -------------------
    //! inverse projection
    //! -------------------
    FILE *f = fopen("D:/FaceMesh.txt","w");
    Eigen::MatrixXd Vt;
    Vt.resize(NbMeshNodes,3);
    for(int i=1; i<=NbMeshNodes; i++)
    {
        double u = V(i-1,0);
        double v = V(i-1,1);
        const gp_Pnt &Pt = myCurGeomSurface->Value(u,v);

        Vt(i-1,0)= Pt.X();
        Vt(i-1,1)= Pt.Y();
        Vt(i-1,2)= Pt.Z();

        cout<<"____("<<Vt(i-1,0)<<", "<<Vt(i-1,1)<<", "<<Vt(i-1,2)<<")____"<<endl;
        fprintf(f,"%lf\t%lf\t%lf\n",Vt(i-1,0),Vt(i-1,1),Vt(i-1,2));
        //fprintf(f,"%lf\t%lf\n",u,v);
    }
    fclose(f);

    cout<<"CustomMesher::buildFaceMesh()->____start building the face mesh in physical space____"<<endl;
    theFaceMeshDS = new Ng_MeshVS_DataSourceFace(Vt,F);
    cout<<"CustomMesher::buildFaceMesh()->____face mesh built____"<<endl;

    return theFaceMeshDS;
}

//! --------------------------------
//! function: faceTag
//! details:  retrieve the face tag
//! --------------------------------
GeometryTag CustomMesher::faceTag(const TopoDS_Face &aFace)
{
    int bodyIndex = -1;
    bool bodyIndexFound = false;
    for(QMap<int,TopoDS_Shape>::iterator it = myMeshDB->bodyMap.begin(); it!=myMeshDB->bodyMap.end(); ++it)
    {
        bodyIndex = it.key();
        const TopTools_IndexedMapOfShape &faceMap = myMeshDB->MapOfBodyTopologyMap.value(bodyIndex).faceMap;
        if(faceMap.Contains(aFace))
        {
            bodyIndexFound = true;
            break;
        }
    }
    GeometryTag aFaceTag;
    aFaceTag.parentShapeNr = bodyIndex;
    aFaceTag.subShapeType = TopAbs_FACE;
    aFaceTag.subTopNr = myMeshDB->MapOfBodyTopologyMap.value(bodyIndex).faceMap.FindIndex(aFace);
    aFaceTag.isParent = false;

    return aFaceTag;
}

//! --------------------------------
//! function: edgeTag
//! details:  retrieve the edge tag
//! --------------------------------
GeometryTag CustomMesher::edgeTag(const TopoDS_Edge &anEdge)
{
    int bodyIndex = -1;
    bool bodyIndexFound = false;
    for(QMap<int,TopoDS_Shape>::iterator it = myMeshDB->bodyMap.begin(); it!=myMeshDB->bodyMap.end(); ++it)
    {
        bodyIndex = it.key();
        const TopTools_IndexedMapOfShape &edgeMap = myMeshDB->MapOfBodyTopologyMap.value(bodyIndex).edgeMap;
        if(edgeMap.Contains(anEdge))
        {
            bodyIndexFound = true;
            break;
        }
    }
    GeometryTag edgeTag;
    edgeTag.parentShapeNr = bodyIndex;
    edgeTag.subShapeType = TopAbs_EDGE;
    edgeTag.subTopNr = myMeshDB->MapOfBodyTopologyMap.value(bodyIndex).edgeMap.FindIndex(anEdge);
    edgeTag.isParent = false;

    return edgeTag;
}

//! -------------------------
//! function: redirectOutput
//! details:
//! -------------------------
void CustomMesher::redirectOutput()
{
    const std::string &message = myTriangleMesher->readAllStandardOutput().toStdString();
    cout<<message<<endl;
}

//! -----------------------------------------------------------
//! function: faceArea
//! details:  compute the area of the face in the (u, v) plane
//!           for automatic sizing purposes
//! -----------------------------------------------------------
double CustomMesher::faceArea(const TopoDS_Face &aFace)
{
    GProp_GProps SProps;
    BRepGProp::SurfaceProperties(aFace,SProps);
    return SProps.Mass();
}
