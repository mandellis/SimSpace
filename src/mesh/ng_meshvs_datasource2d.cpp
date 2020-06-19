//! ----------------
//! custom includes
//! ----------------
#include "ng_meshvs_datasource2d.h"
#include "meshtools.h"
#include <polygon.h>
#include <mesh.h>

//! ----
//! C++
//! ----
#include <set>
#include <map>

IMPLEMENT_STANDARD_HANDLE(Ng_MeshVS_DataSource2D,MeshVS_DataSource)
IMPLEMENT_STANDARD_RTTIEXT(Ng_MeshVS_DataSource2D,MeshVS_DataSource)

//! max number of nodes
const int NMAX = 8;

//! 2D elements - nodal ordering
const int MTRIG3[3]={1,2,3};                // OK
const int MQUAD4[4]={1,2,3,4};              // OK
const int MTRIG6[6]={1,6,2,4,3,5};          // OK
const int MQUAD8[8]={1,5,2,8,3,6,4,7};      // OK


const int MSURFELEM[5][10]={{1,2,3,0,0,0,0,0,0,0},      // TRIG     OK!
                            {1,2,3,4,0,0,0,0,0,0},      // QUAD     OK!
                            {1,6,2,4,3,5,0,0,0,0},      // TRIG6    OK!
                            {1,2,3,4,5,6,0,0,0,0},      // QUAD6    ?
                            {1,5,2,8,3,6,4,7,0,0}};     // QUAD8    OK!


//! -----------------------------------------------
//! function: constructor
//! details:  build the mesh data source from file
//! -----------------------------------------------
Ng_MeshVS_DataSource2D::Ng_MeshVS_DataSource2D(const std::string &fileName)
{
    //cout<<"Ng_MeshVS_DataSource2D::Ng_MeshVS_DataSource2D()->____constructor from file called____"<<endl;

    std::vector<mesh::meshElement> elements;
    std::vector<mesh::meshPoint> nodes;

    bool isDone = this->readMesh(QString::fromStdString(fileName),nodes,elements);
    if(!isDone)
    {
        cout<<"____error in reading the surface mesh____"<<endl;
        return;
    }

    myNumberOfElements=int(elements.size());
    myNumberOfNodes = int(nodes.size());

    //! --------------
    //! further check
    //! --------------
    if(elements.size()<1) return;
    if(nodes.size()<3) return;

    myElemNodes = new TColStd_HArray2OfInteger(1,myNumberOfElements,1,8);
    myElemNormals = new TColStd_HArray2OfReal(1,myNumberOfElements,1,3);
    myElemType = new TColStd_HArray1OfInteger(1,myNumberOfElements);
    myNodeCoords = new TColStd_HArray2OfReal(1,myNumberOfNodes,1,3);

    //! ------------------------------------------
    //! fill the node map and add the coordinates
    //! ------------------------------------------
    for (int i=1;i<=myNumberOfNodes;i++)
    {
        mesh::meshPoint theCurNode = nodes.at(i-1);
        int globalNodeID = theCurNode.ID;
        myNodes.Add(globalNodeID);
        myNodesMap.Add(globalNodeID);
        myNodeCoords->SetValue(i,1,theCurNode.x);
        myNodeCoords->SetValue(i,2,theCurNode.y);
        myNodeCoords->SetValue(i,3,theCurNode.z);
    }

    //! --------------------------
    //! scan the surface elements
    //! --------------------------
    for(int i=1; i<=myNumberOfElements; i++)
    {
        mesh::meshElement anElement = elements.at(i-1);
        int globalElementID = anElement.ID;

        //! fill the map
        myElements.Add(globalElementID);
        myElementsMap.Add(globalElementID);

        //! set the type
        ElemType type = anElement.type;
        myElemType->SetValue(i,type);

        //! set the nodes of the element
        int NbNodes = int(anElement.theNodeIDs.size());
        switch(type)
        {
        case TRIG:
        {
            for(int k=0; k<NbNodes; k++)
            {
                int globalNodeID = anElement.theNodeIDs.at(k);
                myElemNodes->SetValue(i,k+1,globalNodeID);
            }
        }
            break;

        case QUAD:
        {
            for(int k=0; k<NbNodes; k++)
            {
                int globalNodeID = anElement.theNodeIDs.at(k);
                myElemNodes->SetValue(i,k+1,globalNodeID);
            }
        }
            break;

        case TRIG6:
        {
            for(int k=0; k<NbNodes; k++)
            {
                int globalNodeID = anElement.theNodeIDs.at(k);
                myElemNodes->SetValue(i,k+1,globalNodeID);
            }
        }
            break;
        }
    }

    //! ------------------------
    //! compute normal elements
    //! ------------------------
    this->computeNormalAtElements();
}

//! -------------------------------------------------
//! function: constructor
//! details:  build the mesh data source from Netgen
//! -------------------------------------------------
Ng_MeshVS_DataSource2D::Ng_MeshVS_DataSource2D(Ng_Mesh *aMesh)
{
    //! cout<<"Ng_MeshVS_DataSource2D::Ng_MeshVS_DataSource2D->____constructor from Ng_Mesh pointer called____"<<endl;
    if(aMesh==NULL) return;

    myNumberOfNodes = Ng_GetNP(aMesh);
    myNumberOfElements= Ng_GetNSE(aMesh);

    if(myNumberOfNodes<3)
    {
        cout<<"Ng_MeshVS_DataSource2D::Ng_MeshVS_DataSource2D->____number of nodes: "<<myNumberOfNodes<<". Exiting____"<<endl;
        return;
    }
    if(myNumberOfElements<1)
    {
        cout<<"Ng_MeshVS_DataSource2D::Ng_MeshVS_DataSource2D->____number of elements: "<<myNumberOfElements<<". Exiting____"<<endl;
        return;
    }

    //! --------------------
    //! allocate the arrays
    //! --------------------
    myElemNodes = new TColStd_HArray2OfInteger(1,myNumberOfElements,1,NMAX);
    myElemNormals = new TColStd_HArray2OfReal(1,myNumberOfElements,1,3);
    myElemType = new TColStd_HArray1OfInteger(1,myNumberOfElements);
    myNodeCoords = new TColStd_HArray2OfReal(1,myNumberOfNodes,1,3);

    int p_I[NMAX];

    //! -----------------------------------------
    //! add the elements - set the elements type
    //! -----------------------------------------
    for (int i=1;i<=myNumberOfElements;i++)
    {
        myElements.Add(i);
        myElementsMap.Add(i);
        switch(Ng_GetSurfaceElementType(aMesh,i))
        {
        case NG_TRIG: myElemType->SetValue(i,TRIG); break;
        case NG_TRIG6: myElemType->SetValue(i,TRIG6); break;
        case NG_QUAD: myElemType->SetValue(i,QUAD); break;
        case NG_QUAD6: myElemType->SetValue(i,QUAD6); break;
        case NG_QUAD8: myElemType->SetValue(i,QUAD8); break;
        }

        int page;
        Ng_Surface_Element_Type SE_Type = Ng_GetSurfaceElement(aMesh,i,p_I);
        switch(SE_Type)
        {
        case(NG_TRIG): page=0; break;
        case(NG_QUAD): page=1; break;
        case(NG_TRIG6): page=2; break;
        case(NG_QUAD6): page=3; break;
        case(NG_QUAD8): page=4; break;
        }

        for(int j=0;j<NMAX;j++)
        {
            int n=MSURFELEM[page][j]-1;
            myElemNodes->SetValue(i,j+1,p_I[n]);
        }
    }

    //! --------------------------------------------
    //! add the node list and the nodes coordinates
    //! Netgen generates the surface mesh first, so
    //! the globalNodeIDs are {1,2,3,...}, as the
    //! local ones
    //! --------------------------------------------
    double x_I[3];
    for (int i=1;i<=myNumberOfNodes;i++)
    {
        myNodes.Add(i);
        myNodesMap.Add(i);
        Ng_GetPoint(aMesh,i,x_I);
        myNodeCoords->SetValue(i,1,x_I[0]);
        myNodeCoords->SetValue(i,2,x_I[1]);
        myNodeCoords->SetValue(i,3,x_I[2]);
    }

    //! -------------------------------
    //! compute the normal at elements
    //! -------------------------------
    this->computeNormalAtElements();
}

//! ------------------------------------
//! function: constructor
//! details:  init the mesh with Tetgen
//! ------------------------------------
Ng_MeshVS_DataSource2D::Ng_MeshVS_DataSource2D(const tetgenio &aMesh)
{
    cout<<"Ng_MeshVS_DataSource2D::Ng_MeshVS_DataSource2D()->____function called____"<<endl;

    //! --------------
    //! add the nodes
    //! --------------
    std::vector<std::vector<double>> vectorOfPoints;
    int N_Points = aMesh.numberofpoints;

    //! --------------------------------------------------------
    //! extract the mesh points belonging to the surface mesh
    //! if a point has a marker !=0 the point is on the surface
    //! the position "0" within the list corresponds to the
    //! Tetgen node ID "1"
    //! --------------------------------------------------------
    for(int i=0; i<N_Points; i++)
    {
        int pointMarker = aMesh.pointmarkerlist[i];
        if(pointMarker!=0)
        {
            int globalNodeID = i+1;
            myNodes.Add(globalNodeID);       //! global node ID
            myNodesMap.Add(globalNodeID);    //! global node ID

            std::vector<double> point;
            for(int j=0;j<3;j++)
            {
                double coord = aMesh.pointlist[i*3+j];
                point.push_back(coord);
            }
            vectorOfPoints.push_back(point);
        }
    }
    if(vectorOfPoints.size()<3)
    {
        cout<<"Ng_MeshVS_DataSource2D::Ng_MeshVS_DataSource2D()->____no surface mesh points. Exiting____"<<endl;
        return;
    }

    //! ---------------------------------------------------
    //! allocate the space for nodes and nodal coordinates
    //! ---------------------------------------------------
    myNumberOfNodes = int(vectorOfPoints.size());
    myNodeCoords = new TColStd_HArray2OfReal(1,myNumberOfNodes,1,3);

    cout<<"Ng_MeshVS_DataSource2D::Ng_MeshVS_DataSource2D()->____total number of surface mesh points: "<<myNumberOfNodes<<"____"<<endl;

    //! -------------------------------------------------
    //! vectorOfPoints contains the nodes on the surface
    //! -------------------------------------------------
    int localNodeID=1;
    for(std::vector<std::vector<double>>::iterator itPoints= vectorOfPoints.begin() ;itPoints!= vectorOfPoints.end(); ++itPoints)
    {
        std::vector<double> point = *itPoints;
        for(int m=0; m<3; m++)
        {
            double coord = point.at(m);
            myNodeCoords->SetValue(localNodeID,m+1,coord);
        }
        localNodeID++;
    }
    cout<<"Ng_MeshVS_DataSource2D::Ng_MeshVS_DataSource2D()->____nodes coordinates added____"<<endl;

    //! -------------------------
    //! add the surface elements
    //! "surfElement" is a triangle
    //! ----------------------------
    std::vector<std::vector<int>> surfElements;

    //! -----------------------------------------
    //! all the triangles in the overall 3D mesh
    //! -----------------------------------------
    int N_TriFaces = aMesh.numberoftrifaces;

    cout<<"Ng_MeshVS_DataSource2D::Ng_MeshVS_DataSource2D()->____total number of mesh triangles of the volume mesh: "<<N_TriFaces<<"____"<<endl;

    for(int i=0; i<N_TriFaces; i++)
    {
        //! ------------------------------------
        //! the triangle belongs to the surface
        //! ------------------------------------
        int faceMarker = aMesh.trifacemarkerlist[i];
        if(faceMarker!=0)
        {
            int globalElementID = i+1;
            myElements.Add(globalElementID);          //! global element ID
            myElementsMap.Add(globalElementID);       //! global element ID

            //cout<<"/---------------------------------/"<<endl;
            //cout<<"/   defining element "<<globalElementID<<endl;
            //cout<<"/---------------------------------/"<<endl;

            std::vector<int> triad;
            for(int j=0;j<3;j++)
            {
                //! -------------------------------------
                //! "pos" is in the Tetgen list of nodes
                //! -------------------------------------
                int pos = aMesh.trifacelist[i*3+j];
                //cout<<"____pos: "<<pos<<"____"<<endl;
                triad.push_back(pos);
            }
            surfElements.push_back(triad);
        }
    }
    myNumberOfElements = int(surfElements.size());
    if(myNumberOfElements<1)
    {
        cout<<"Ng_MeshVS_DataSource2D::Ng_MeshVS_DataSource2D()->____the surface mesh has no elements. Exiting____"<<endl;
        return;
    }
    cout<<"Ng_MeshVS_DataSource2D::Ng_MeshVS_DataSource2D()->____total number of surface elements: "<<myNumberOfElements<<"____"<<endl;

    //! ------------------------------------
    //! first order triangle for the moment
    //! ------------------------------------
    myElemNodes = new TColStd_HArray2OfInteger(1,myNumberOfElements,1,3);
    myElemType = new TColStd_HArray1OfInteger(1,myNumberOfElements);
    myElemNormals = new TColStd_HArray2OfReal(1,myNumberOfElements,1,3);

    //! ------------------------------------------
    //! iterate over the triangles of the surface
    //! ------------------------------------------
    int localElementID = 0;
    for(int i=0; i<myNumberOfElements; i++)
    {
        localElementID++;
        const std::vector<int> &aSurfaceTriangle = surfElements.at(i);

        //int k=1;
        for(int j=0; j<aSurfaceTriangle.size(); j++)
        {
            int curNode = aSurfaceTriangle.at(j);
            int globalNodeID = curNode;

            //int pos = C[k-1];
            //k++;
            myElemNodes->SetValue(localElementID,j+1,globalNodeID);
        }

        //! -----------------------------------------------------------
        //! set the element type - only linear triangle for the moment
        //! -----------------------------------------------------------
        myElemType->SetValue(localElementID,TRIG);
    }
    //! --------------------------
    //! compute normal at elements
    //! ---------------------------
    this->computeNormalAtElements();
}

//! ---------------------------------------------------------------
//! function: constructor from Tetgen files
//! details:  the element ordering is arbitrary, since the surface
//!           mesh is obtained from a volume mesh
//! ---------------------------------------------------------------
Ng_MeshVS_DataSource2D::Ng_MeshVS_DataSource2D(const QString &faceFileName, const QString &nodeFileName)
{
    cout<<"Ng_MeshVS_DataSource2D::Ng_MeshVS_DataSource2D()->____constructor from Tetgen files called____"<<endl;

    QList<mesh::meshPoint> meshPointList;

    //! ----------------------------------------------------------
    //! the .node file contains the node of the whole volume mesh
    //! ----------------------------------------------------------
    myNumberOfNodes = 0;
    FILE *nodeFile = fopen(nodeFileName.toStdString().c_str(),"r");
    if(nodeFile==NULL) return;

    //! ----------------------------------------------------------------------------------------
    //! First line: <# of points> <dimension (3)> <# of attributes> <boundary markers (0 or 1)>
    //! ----------------------------------------------------------------------------------------
    cout<<"____Start reading: "<<nodeFileName.toStdString()<<"____"<<endl;
    int NbPoints;
    int dimension;
    int NbAttributes;
    int boundaryMarkerTag;
    int pointNumber;
    double x,y,z;
    int tag;
    fscanf(nodeFile,"%d%d%d%d",&NbPoints,&dimension,&NbAttributes,&boundaryMarkerTag);

    cout<<"____Start reading mesh points___"<<endl;
    int n=0;
    for(int line=1; line<=NbPoints; line++)
    {
        fscanf(nodeFile,"%d%lf%lf%lf%d",&pointNumber,&x,&y,&z,&tag);
        //if(tag!=0)
        {
            n++;
            mesh::meshPoint P;
            P.ID = pointNumber;
            P.x = x;
            P.y = y;
            P.z = z;
            meshPointList<<P;

            myNumberOfNodes++;
            myNodes.Add(pointNumber);
            myNodesMap.Add(pointNumber);
        }
    }
    cout<<"____mesh points read___"<<endl;
    cout<<"____closing the .node file____"<<endl;
    fclose(nodeFile);

    if(myNumberOfNodes<3)
    {
        cout<<"Ng_MeshVS_DataSource2D::Ng_MeshVS_DataSource2D()->____no points on the surface mesh. Exiting____"<<endl;
        return;
    }

    myNodeCoords = new TColStd_HArray2OfReal(1,myNumberOfNodes,1,3);
    for(int i=1; i<=meshPointList.length(); i++)
    {
        const mesh::meshPoint &P = meshPointList.at(i-1);
        myNodeCoords->SetValue(i,1,P.x);
        myNodeCoords->SetValue(i,2,P.y);
        myNodeCoords->SetValue(i,3,P.z);
    }

    //! ----------------------
    //! reading the face file
    //! ----------------------
    FILE *faceFile = fopen(faceFileName.toStdString().c_str(),"r");
    if(faceFile==NULL) return;

    //! --------------------
    //! read the first line
    //! --------------------
    cout<<"____Start reading: "<<faceFileName.toStdString()<<"____"<<endl;
    int NbElements;
    int boundaryTag;
    fscanf(faceFile,"%d%d",&NbElements,&boundaryTag);
    int number;
    int n1,n2,n3,n4,n5,n6;
    int faceTag;
    int localElementID = 0;

    cout<<"____Number of surface elements: "<<NbElements<<"____"<<endl;

    myNumberOfElements = NbElements;
    myElemNodes = new TColStd_HArray2OfInteger(1,NbElements,1,6);
    myElemType = new TColStd_HArray1OfInteger(1,NbElements);
    myElemNormals = new TColStd_HArray2OfReal(1,NbElements,1,3);

    for(int line=1; line<=NbElements; line++)
    {
        if(5==fscanf(faceFile,"%d%d%d%d%d",&number,&n1,&n2,&n3,&faceTag))
        {
            localElementID++;
            myElemType->SetValue(localElementID,TRIG);

            myElementsMap.Add(localElementID);
            myElements.Add(localElementID);

            int node1 = myNodesMap.FindIndex(n1);
            int node2 = myNodesMap.FindIndex(n2);
            int node3 = myNodesMap.FindIndex(n3);

            //! --------------------------------------
            //! check Tetgen nodal ordering ... to do
            //! --------------------------------------
            myElemNodes->SetValue(localElementID,1,node1);
            myElemNodes->SetValue(localElementID,2,node3);
            myElemNodes->SetValue(localElementID,3,node2);
        }
        else if(8==fscanf(faceFile,"%d%d%d%d%d%d%d%d",&number,&n1,&n2,&n3,&n4,&n5,&n6,&faceTag))
        {
            localElementID++;
            myElemType->SetValue(localElementID,TRIG);

            myElementsMap.Add(localElementID);
            myElements.Add(localElementID);

            int node1 = myNodesMap.FindIndex(n1);
            int node2 = myNodesMap.FindIndex(n2);
            int node3 = myNodesMap.FindIndex(n3);
            int node4 = myNodesMap.FindIndex(n4);
            int node5 = myNodesMap.FindIndex(n5);
            int node6 = myNodesMap.FindIndex(n6);

            //! --------------------------------------
            //! check Tetgen nodal ordering ... to do
            //! --------------------------------------
            myElemNodes->SetValue(localElementID,1,node1);
            myElemNodes->SetValue(localElementID,2,node2);
            myElemNodes->SetValue(localElementID,3,node3);
            myElemNodes->SetValue(localElementID,4,node4);
            myElemNodes->SetValue(localElementID,5,node5);
            myElemNodes->SetValue(localElementID,6,node6);
        }
    }

    //! ---------------
    //! close the file
    //! ---------------
    fclose(faceFile);

    //! ---------------------------
    //! compute normal at elements
    //! ---------------------------
    this->computeNormalAtElements();
}

//! ---------------------------------------------------------------------
//! function: constructor
//! details:  build a surface mesh data source using using a volume mesh
//!           data source using the face to element connectivity info
//! ---------------------------------------------------------------------
Ng_MeshVS_DataSource2D::Ng_MeshVS_DataSource2D(const opencascade::handle<Ng_MeshVS_DataSource3D> &a3DMeshVS_DataSource)
{
    //!cout<<"Ng_MeshVS_DataSource2D::Ng_MeshVS_DataSource2D()->____constructor from a volume mesh data source called____"<<endl;
    if(a3DMeshVS_DataSource->GetAllNodes().Extent()<3)
    {
        cout<<"Ng_MeshVS_DataSource2D::Ng_MeshVS_DataSource2D()->____constructor from a volume mesh: wrong number of nodes____"<<endl;
        return;
    }
    if(a3DMeshVS_DataSource->GetAllElements().Extent()<1)
    {
        cout<<"Ng_MeshVS_DataSource2D::Ng_MeshVS_DataSource2D()->____constructor from a volume mesh: no element____"<<endl;
        return;
    }

    const std::vector<meshElement2D> &surfaceElements = a3DMeshVS_DataSource->mySurfaceElements;
    myNumberOfElements = int(surfaceElements.size());
    if(myNumberOfElements<1) return;

    myElemType = new TColStd_HArray1OfInteger(1,myNumberOfElements);
    myElemNodes = new TColStd_HArray2OfInteger(1,myNumberOfElements,1,8);
    myElemType = new TColStd_HArray1OfInteger(1,myNumberOfElements);
    myElemNormals = new TColStd_HArray2OfReal(1,myNumberOfElements,1,3);

    int local2DElementID = 0;
    for(std::vector<meshElement2D>::const_iterator it = surfaceElements.cbegin(); it!=surfaceElements.cend(); ++it)
    {
        local2DElementID++;
        const meshElement2D &aMeshElement = *it;

        //! --------------------------------------------------
        //! element type determined using the number of nodes
        //! --------------------------------------------------
        int NbNodes = aMeshElement.nodeIDs.length();
        switch(NbNodes)
        {
        case 3: myElemType->SetValue(local2DElementID,TRIG); break;
        case 4: myElemType->SetValue(local2DElementID,QUAD); break;
        case 6: myElemType->SetValue(local2DElementID,TRIG6); break;
        case 8: myElemType->SetValue(local2DElementID,QUAD8); break;
        }
        //! -----------------------------------------------------------------------
        //! Note: this constructor generates the surface mesh from a pre-generated
        //! volume mesh, by consequence the elements are stored into the map using
        //! the arbitrary global numbering {1,2,...,NSE}
        //! -----------------------------------------------------------------------
        myElements.Add(local2DElementID);
        myElementsMap.Add(local2DElementID);

        for(int i=0; i<NbNodes; i++)
        {
            //! ---------------------------------------------------
            //! globalNodeID is the node number within the 3D mesh
            //! ---------------------------------------------------
            int globalNodeID = aMeshElement.nodeIDs.at(i);
            myElemNodes->SetValue(local2DElementID,i+1,globalNodeID);

            if(!myNodes.Contains(globalNodeID))
            {
                myNodes.Add(globalNodeID);
                myNodesMap.Add(globalNodeID);
            }
        }
    }

    myNumberOfNodes = myNodes.Extent();
    if(myNumberOfNodes<3) return;
    myNodeCoords = new TColStd_HArray2OfReal(1,myNumberOfNodes,1,3);

    for(int localNodeID=1; localNodeID<=myNumberOfNodes; localNodeID++)
    {
        int globalNodeID = myNodesMap.FindKey(localNodeID);

        //! -----------------------------------------------------------
        //! retrieve the node coordinates from the 3D mesh data source
        //! -----------------------------------------------------------
        int NbNodes;
        MeshVS_EntityType type;
        double buf[3];
        TColStd_Array1OfReal coords(*buf,1,3);
        a3DMeshVS_DataSource->GetGeom(globalNodeID,false,coords,NbNodes,type);

        myNodeCoords->SetValue(localNodeID,1,coords(1));
        myNodeCoords->SetValue(localNodeID,2,coords(2));
        myNodeCoords->SetValue(localNodeID,3,coords(3));
    }

    //! --------------------------------
    //! compute the normals at elements
    //! --------------------------------
    this->computeNormalAtElements();
}

//! ----------------------------
//! function: constructor
//! details:  from StlMesh_Mesh
//! ----------------------------
Ng_MeshVS_DataSource2D::Ng_MeshVS_DataSource2D(const occHandle(StlMesh_Mesh) &aMesh)
{
    //!cout<<"Ng_MeshVS_DataSource2D::Ng_MeshVS_DataSource2D()->____constructor called: surface mesh from StlMesh_Mesh____"<<endl;

    if(aMesh.IsNull()) return;
    if(aMesh->NbTriangles()<1 || aMesh->NbVertices()<3) return;

    //int NbDomains = aMesh->NbDomains();

    myNumberOfNodes = aMesh->NbVertices();
    myNumberOfElements = aMesh->NbTriangles();
    const TColgp_SequenceOfXYZ& aCoords = aMesh->Vertices();
    Standard_Integer len = aCoords.Length(),j;
    myNodeCoords = new TColStd_HArray2OfReal(1,len,1,3);

    gp_XYZ xyz;

    for(int i=1;i<=len;i++)
    {
        myNodes.Add(i);
        myNodesMap.Add(i);

        xyz = aCoords(i);
        myNodeCoords->SetValue(i, 1, xyz.X());
        myNodeCoords->SetValue(i, 2, xyz.Y());
        myNodeCoords->SetValue(i, 3, xyz.Z());
    }

    const StlMesh_SequenceOfMeshTriangle& aSeq = aMesh->Triangles();
    len = aSeq.Length();
    myElemNormals = new TColStd_HArray2OfReal(1,len,1,3);
    myElemNodes = new TColStd_HArray2OfInteger(1,len,1,3);
    myElemType = new TColStd_HArray1OfInteger(1,len);

    for(int i=1;i<=len;i++)
    {
        myElements.Add(i);
        myElementsMap.Add(i);
        myElemType->SetValue(i,TRIG);

        occHandle(StlMesh_MeshTriangle) aTriangle = aSeq.Value(i);
        Standard_Integer V[3];
        Standard_Real nx, ny, nz;

        aTriangle->GetVertexAndOrientation(V[0],V[1],V[2],nx,ny,nz);

        for(j=0;j<3;j++) myElemNodes->SetValue(i,j+1,V[j]);

        myElemNormals->SetValue(i,1,nx);
        myElemNormals->SetValue(i,2,ny);
        myElemNormals->SetValue(i,3,nz);
    }
}

//! ----------------------------------------------
//! function: constructor
//! details:  from a vector of elements and nodes
//! ----------------------------------------------
Ng_MeshVS_DataSource2D::Ng_MeshVS_DataSource2D(const std::vector<mesh::meshElement> &elements, std::vector<mesh::meshPoint> &nodes)
{
    //cout<<"Ng_MeshVS_DataSource2D::Ng_MeshVS_DataSource2D()->____constructor called____"<<endl;
    myNumberOfElements = int(elements.size());
    myNumberOfNodes = int(nodes.size());

    if(myNumberOfElements<1)
    {
        cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____no elements: "<<myNumberOfElements<<"____"<<endl;
        return;
    }
    if(myNumberOfNodes<3)
    {
        cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____wrong number of nodes: "<<myNumberOfNodes<<"____"<<endl;
        return;
    }

    myElemNodes = new TColStd_HArray2OfInteger(1,myNumberOfElements,1,8);
    myElemNormals = new TColStd_HArray2OfReal(1,myNumberOfElements,1,3);
    myElemType = new TColStd_HArray1OfInteger(1,myNumberOfElements);
    myNodeCoords = new TColStd_HArray2OfReal(1,myNumberOfNodes,1,3);

    //! ------------------------------------------
    //! fill the node map and add the coordinates
    //! ------------------------------------------
    //cout<<"____adding nodes____"<<endl;
    for(int i=1;i<=myNumberOfNodes;i++)
    {
        mesh::meshPoint theCurNode = nodes.at(i-1);
        int globalNodeID = theCurNode.ID;
        myNodes.Add(globalNodeID);
        myNodesMap.Add(globalNodeID);
        myNodeCoords->SetValue(i,1,theCurNode.x);
        myNodeCoords->SetValue(i,2,theCurNode.y);
        myNodeCoords->SetValue(i,3,theCurNode.z);
    }

    //! ----------------------------------------------
    //! fill the element maps and define the elements
    //! ----------------------------------------------
    //cout<<"____adding elements____"<<endl;
    int i=1;
    for(std::vector<mesh::meshElement>::const_iterator it = elements.cbegin(); it!=elements.cend(); ++it, i++)
    {
        mesh::meshElement anElement = *it;
        int globalElementID = anElement.ID;
        //cout<<"____global element ID: "<<globalElementID<<"____"<<endl;

        //! -------------------------
        //! fill the map of elements
        //! -------------------------
        myElements.Add(globalElementID);
        myElementsMap.Add(globalElementID);
        //cout<<"____global element ID: "<<globalElementID<<" added____"<<endl;

        //! ----------------------------
        //! set the type of the element
        //! ----------------------------
        ElemType type = anElement.type;
        myElemType->SetValue(i,type);
        //if(type==TRIG) cout<<"____element type \"TRIG\" set____"<<endl;

        //! -----------------------------------
        //! set the inverse transformation [?]
        //! -----------------------------------
        //int MTRIG3_I[3]={1,2,3};
        //int MTRIG6_I[6]={1,3,5,4,6,2};

        //! -----------------------------
        //! set the nodes of the element
        //! -----------------------------
        switch(type)
        {
        case TRIG:
            for(int k=0; k<anElement.theNodeIDs.size(); k++)
            {
                int globalNodeID = anElement.theNodeIDs.at(k);
                //int pos = MTRIG3_I[k];
                myElemNodes->SetValue(i,k+1,globalNodeID);
            }
            break;
        }
    }
    //! ----------------------------
    //! compute normals at elements
    //! ----------------------------
    this->computeNormalAtElements();
}


//! ------------------------------------------------------------
//! function: GetGeom
//! details:  ID is number into the map myNodeMap, myElementMap
//! ------------------------------------------------------------
Standard_Boolean Ng_MeshVS_DataSource2D::GetGeom(const Standard_Integer ID,
                                                 const Standard_Boolean IsElement,
                                                 TColStd_Array1OfReal& Coords,
                                                 Standard_Integer& NbNodes,
                                                 MeshVS_EntityType& Type) const
{
    if(IsElement)
    {
        //cout<<"Ng_MeshVS_DataSource2D::GetGeom()->____function called for element ID: "<<ID<<"____"<<endl;
        int localElementID = myElementsMap.FindIndex(ID);        
        if(localElementID>=1 && localElementID<=myNumberOfElements)
        {
            Type = MeshVS_ET_Face;
            switch(myElemType->Value(localElementID))
            {
            case(TRIG): NbNodes=3; break;
            case(QUAD): NbNodes=4; break;
            case(TRIG6): NbNodes=6; break;
            case(QUAD8): NbNodes=8; break;
            }

            for(int i=1, k=1; i<=NbNodes; i++)
            {
                int globalNodeID = myElemNodes->Value(localElementID,i);
                int localNodeID = myNodesMap.FindIndex(globalNodeID);
                for(int j=1; j<=3; j++)
                {
                    Coords(k)=myNodeCoords->Value(localNodeID,j);
                    k++;
                }
            }
            return Standard_True;
        }
        else
        {
            cerr<<"Ng_MeshVS_DataSourceFace::GetGeom()->____error: element ID: "<<ID<<" out of range____"<<endl;
            return Standard_False;
        }
    }
    else
    {
        //cout<<"Ng_MeshVS_DataSource2D::GetGeom()->____function called for node ID: "<<ID<<"____"<<endl;
        int localNodeID = myNodesMap.FindIndex(ID);
        if(localNodeID>=1 && localNodeID<=myNumberOfNodes)
        {
            Type = MeshVS_ET_Node;
            NbNodes=1;
            Coords(1) = myNodeCoords->Value(localNodeID,1);
            Coords(2) = myNodeCoords->Value(localNodeID,2);
            Coords(3) = myNodeCoords->Value(localNodeID,3);
            return Standard_True;
        }
        else
        {
            cerr<<"Ng_MeshVS_DataSourceFace::GetGeom()->____error: node "<<ID<<" out of range____"<<endl;
            return Standard_False;
        }
    }
}

//! ----------------------
//! function: GetGeomType
//! details:
//! ----------------------
Standard_Boolean Ng_MeshVS_DataSource2D::GetGeomType(const Standard_Integer ID,
                                                    const Standard_Boolean IsElement,
                                                    MeshVS_EntityType& Type) const
{
    //cout<<"Ng_MeshVS_DataSource2D::GetGeomType()->____function called for ID: "<<ID<<"____"<<endl;

    if(IsElement)
    {
        int localElementID = myElementsMap.FindIndex(ID);
        if(localElementID>=1 && localElementID<=myNumberOfElements)
        {
            Type = MeshVS_ET_Face;
            return Standard_True;
        }
        else
        {
            cerr<<"Ng_MeshVS_DataSource2D::GetGeomType()->____element ID: "<<ID<<"out of range____"<<endl;
            return Standard_False;
        }
    }
    else
    {
        int localNodeID = myNodesMap.FindIndex(ID);
        if(localNodeID>=1 && localNodeID<=myNumberOfNodes)
        {
            Type = MeshVS_ET_Node;
            return Standard_True;
        }
        else
        {
            cerr<<"Ng_MeshVS_DataSource2D::GetGeomType()->____node ID: "<<ID<<" out of range____"<<endl;
            return Standard_False;
        }
    }
}

//! ------------------
//! function: GetAddr
//! details:
//! ------------------
Standard_Address Ng_MeshVS_DataSource2D::GetAddr(const Standard_Integer,
                                                const Standard_Boolean) const
{
    return NULL;
}

//! ----------------------------
//! function: GetNodesByElement
//! details:
//! ----------------------------
Standard_Boolean Ng_MeshVS_DataSource2D::GetNodesByElement(const Standard_Integer ID,
                                                          TColStd_Array1OfInteger& theNodeIDs,
                                                          Standard_Integer& theNbNodes) const
{
    //cout<<"Ng_MeshVS_DataSource2D::GetNodesByElement->____function called for ID: "<<ID<<"____"<<endl;
    int localElementID = myElementsMap.FindIndex(ID);
    if(localElementID<1 && localElementID>myNumberOfElements)
    {
        cerr<<"Ng_MeshVS_DataSource2D::GetNodesByElement()->____element ID: "<<ID<<" out of range____"<<endl;
        return Standard_False;
    }
    switch(myElemType->Value(localElementID))
    {
    case(TRIG):
        theNbNodes=3;
        theNodeIDs(1)= myElemNodes->Value(localElementID,1);
        theNodeIDs(2)= myElemNodes->Value(localElementID,2);
        theNodeIDs(3)= myElemNodes->Value(localElementID,3);
        break;
    case(QUAD):
        theNbNodes=4;
        theNodeIDs(1)= myElemNodes->Value(localElementID,1);
        theNodeIDs(2)= myElemNodes->Value(localElementID,2);
        theNodeIDs(3)= myElemNodes->Value(localElementID,3);
        theNodeIDs(4)= myElemNodes->Value(localElementID,4);
        break;
    case(TRIG6):
        theNbNodes=6;
        theNodeIDs(1)= myElemNodes->Value(localElementID,1);
        theNodeIDs(2)= myElemNodes->Value(localElementID,2);
        theNodeIDs(3)= myElemNodes->Value(localElementID,3);
        theNodeIDs(4)= myElemNodes->Value(localElementID,4);
        theNodeIDs(5)= myElemNodes->Value(localElementID,5);
        theNodeIDs(6)= myElemNodes->Value(localElementID,6);
        break;
    case(QUAD8):
        theNbNodes=8;
        theNodeIDs(1)= myElemNodes->Value(localElementID,1);
        theNodeIDs(2)= myElemNodes->Value(localElementID,2);
        theNodeIDs(3)= myElemNodes->Value(localElementID,3);
        theNodeIDs(4)= myElemNodes->Value(localElementID,4);
        theNodeIDs(5)= myElemNodes->Value(localElementID,5);
        theNodeIDs(6)= myElemNodes->Value(localElementID,6);
        theNodeIDs(7)= myElemNodes->Value(localElementID,7);
        theNodeIDs(8)= myElemNodes->Value(localElementID,8);

        /*
            theNodeIDs(1)= myElemNodes->Value(localElementID,1);
            theNodeIDs(2)= myElemNodes->Value(localElementID,5);
            theNodeIDs(3)= myElemNodes->Value(localElementID,2);
            theNodeIDs(4)= myElemNodes->Value(localElementID,8);
            theNodeIDs(5)= myElemNodes->Value(localElementID,3);
            theNodeIDs(6)= myElemNodes->Value(localElementID,6);
            theNodeIDs(7)= myElemNodes->Value(localElementID,4);
            theNodeIDs(8)= myElemNodes->Value(localElementID,7);
            */
        break;
    }
    return Standard_True;
}

//! ----------------------
//! function: GetAllNodes
//! details:
//! ----------------------
const TColStd_PackedMapOfInteger& Ng_MeshVS_DataSource2D::GetAllNodes() const
{
    return myNodes;
}

//! -------------------------
//! function: GetAllElements
//! details:
//! -------------------------
const TColStd_PackedMapOfInteger& Ng_MeshVS_DataSource2D::GetAllElements() const
{
    return myElements;
}

//! -------------------------
//! function: GetElementType
//! details:
//! -------------------------
const bool Ng_MeshVS_DataSource2D::GetElementType(ElemType &eType, int elementID, bool isLocal) const
{
    int ID = -1;
    if(isLocal) ID = elementID;
    else ID = myElementsMap.FindIndex(elementID);

    if(ID<0 || ID >myNumberOfElements)
    {
        eType = NULL_ELEMENT;
        return false;
    }
    eType = static_cast<ElemType>(myElemType->Value(ID));
    return true;
}

//! --------------------
//! function: GetNormal
//! details:
//! --------------------
Standard_Boolean Ng_MeshVS_DataSource2D::GetNormal(const Standard_Integer Id,
                                                   const Standard_Integer Max,
                                                   Standard_Real& nx,
                                                   Standard_Real& ny,
                                                   Standard_Real& nz) const
{
    //cout<<"Ng_MeshVS_DataSource2D::GetNormal()->____function called for ID: "<<Id<<"____"<<endl;
    Q_UNUSED(Max)

    int localElementID = myElementsMap.FindIndex(Id);
    if(localElementID>=1 && localElementID<=myNumberOfElements)
    {
        nx = myElemNormals->Value(localElementID,1);
        ny = myElemNormals->Value(localElementID,2);
        nz = myElemNormals->Value(localElementID,3);
        return Standard_True;
    }
    else
    {
        cerr<<"Ng_MeshVS_DataSource2D::GetNormal()->____error: element ID: "<<Id<<"out of range____"<<endl;
        return Standard_False;
    }
}

//! --------------------------------------------------------------
//! function: writeMesh
//! details:  "dim" mesh dimension value is redundant, since
//!           GetGeomType() is called foreach element in the mesh
//! --------------------------------------------------------------
bool Ng_MeshVS_DataSource2D::writeMesh(const QString &meshFileName, int dim)
{
    cout<<"Ng_MeshVS_DataSource2D::writeMesh()->____function called. Writing mesh file: "<<meshFileName.toStdString()<<"____"<<endl;

    char name[128];
    sprintf(name,"%s",meshFileName.toStdString().c_str());
    FILE *fileout = fopen(name,"w");
    if(fileout==NULL)
    {
        cout<<"Ng_MeshVS_DataSource2D::writeMesh()->____cannot open file for writing mesh____"<<endl;
        return false;
    }

    //! ----------------------
    //! writing the dimension
    //! ----------------------
    fprintf(fileout,"%d\n",dim);

    //! ----------------
    //! write the nodes
    //! ----------------
    const TColStd_PackedMapOfInteger& aNodes = this->GetAllNodes();
    fprintf(fileout,"%d\n",aNodes.Extent());

    double buf[3];
    TColStd_Array1OfReal aCoords(*buf,1,3);
    int nbNodes;
    MeshVS_EntityType aType;
    TColStd_MapIteratorOfPackedMapOfInteger anIter;

    for(anIter.Initialize(aNodes);anIter.More();anIter.Next())
    {
        int aKey = anIter.Key();
        if(this->GetGeom(aKey,false,aCoords,nbNodes,aType)==true)
        {
            //! write the nodeID (global) and the coordinates
            fprintf(fileout,"%d\t%.9e\t%.9e\t%.9e\n",aKey,aCoords.Value(1),aCoords.Value(2),aCoords.Value(3));
        }
    }

    //! -------------------
    //! write the elements
    //! -------------------
    TColStd_PackedMapOfInteger mapOfElements = this->GetAllElements();

    //! -----------------------------
    //! write the number of elements
    //! -----------------------------
    fprintf(fileout,"%d\n",mapOfElements.Extent());

    for(anIter.Initialize(mapOfElements);anIter.More();anIter.Next())
    {
        int aKey = anIter.Key();
        int NbNodes;
        int aNodesBuf[20];
        TColStd_Array1OfInteger nodeIDs(*aNodesBuf,1,20);
        this->GetNodesByElement(aKey,nodeIDs,NbNodes);

        MeshVS_EntityType type;
        this->GetGeomType(aKey,true,type);

        //! -----------------------
        //! write the element type
        //! -----------------------
        char header[16];
        if(type==MeshVS_ET_Face)
        {
            switch(NbNodes)
            {
            case 3: sprintf(header,"TRIG"); break;
            case 4: sprintf(header,"QUAD"); break;
            case 6: sprintf(header,"TRIG6"); break;
            case 8: sprintf(header,"QUAD8"); break;
            }
        }
        else if(type==MeshVS_ET_Volume)
        {
            switch(NbNodes)
            {
            case 4: sprintf(header,"TET"); break;
            case 8: sprintf(header,"HEXA"); break;
            case 10: sprintf(header,"TET10"); break;
            case 20: sprintf(header,"HEXA20"); break;
            }
        }

        //! ----------------------------------------
        //! write the element key, the element type
        //! ----------------------------------------
        fprintf(fileout,"%d\n",aKey);
        fprintf(fileout,"%s\n",header);

        //! ------------------
        //! write the nodeIDs
        //! ------------------
        for(int i=1; i<NbNodes; i++)
        {
            fprintf(fileout,"%d\t",nodeIDs(i));
        }
        fprintf(fileout,"%d\n",nodeIDs(NbNodes));
    }

    //! ---------------
    //! close the file
    //! ---------------
    fclose(fileout);
    cout<<"Ng_MeshVS_DataSource2D::writeMesh()->____exiting____"<<endl;
    return true;
}

//! -------------------
//! function: readMesh
//! details:
//! -------------------
bool Ng_MeshVS_DataSource2D::readMesh(const QString &meshFileName,
                                      std::vector<mesh::meshPoint> &nodes,
                                      std::vector<mesh::meshElement> &elements)
{
    cout<<"Ng_MeshVS_DataSource2D::readMesh()->____function called____"<<endl;
    FILE *filein = fopen(meshFileName.toStdString().c_str(),"r");
    if(filein!=NULL)
    {
        //! -------------------
        //! read the dimension
        //! -------------------
        int dim;
        fscanf(filein,"%d",&dim);

        //! -------------------------
        //! read the number of nodes
        //! -------------------------
        int NbNodes;
        fscanf(filein,"%d",&NbNodes);

        //! ---------------
        //! read the nodes
        //! ---------------
        for(int i=0; i<NbNodes;i++)
        {
            mesh::meshPoint aNode;
            fscanf(filein,"%d%lf%lf%lf",&aNode.ID,&aNode.x,&aNode.y,&aNode.z);
            nodes.push_back(aNode);
        }

        //! ----------------------------
        //! read the number of elements
        //! ----------------------------
        int NbElements;
        fscanf(filein,"%d",&NbElements);

        //! ------------------
        //! read the elements
        //! ------------------
        for(int i=0; i<NbElements; i++)
        {
            mesh::meshElement anElement;
            char eType[16];

            fscanf(filein,"%d",&anElement.ID);
            fscanf(filein,"%s",eType);

            if(strcmp(eType,"TRIG")==0)
            {
                anElement.type= TRIG;
                int V1,V2,V3;
                if(3!=fscanf(filein,"%d%d%d",&V1,&V2,&V3))
                {
                    cout<<"____error in reading the mesh data source____"<<endl;
                    return false;
                }
                anElement.theNodeIDs.push_back(V1);
                anElement.theNodeIDs.push_back(V2);
                anElement.theNodeIDs.push_back(V3);
            }
            if(strcmp(eType,"QUAD")==0)
            {
                anElement.type= QUAD;
                int V1,V2,V3,V4;
                if(4!=fscanf(filein,"%d%d%d%d",&V1,&V2,&V3,&V4))
                {
                    cout<<"____error in reading the mesh data source____"<<endl;
                    return false;
                }
                anElement.theNodeIDs.push_back(V1);
                anElement.theNodeIDs.push_back(V2);
                anElement.theNodeIDs.push_back(V3);
                anElement.theNodeIDs.push_back(V4);
            }

            if(strcmp(eType,"TRIG6")==0)
            {
                anElement.type= TRIG6;
                int V1,V2,V3,V4,V5,V6;
                if(6!=fscanf(filein,"%d%d%d%d%d%d",&V1,&V2,&V3,&V4,&V5,&V6))
                {
                    cout<<"____error in reading the mesh data source____"<<endl;
                    return false;
                }
                anElement.theNodeIDs.push_back(V1);
                anElement.theNodeIDs.push_back(V2);
                anElement.theNodeIDs.push_back(V3);
                anElement.theNodeIDs.push_back(V4);
                anElement.theNodeIDs.push_back(V5);
                anElement.theNodeIDs.push_back(V6);
            }
            else if(strcmp(eType,"TET")==0)
            {
                anElement.type= TET;
                int V1,V2,V3,V4;
                if(4!=fscanf(filein,"%d%d%d%d",&V1,&V2,&V3,&V4))
                {
                    cout<<"____error in reading the mesh data source____"<<endl;
                    return false;
                }
                anElement.theNodeIDs.push_back(V1);
                anElement.theNodeIDs.push_back(V2);
                anElement.theNodeIDs.push_back(V3);
                anElement.theNodeIDs.push_back(V4);
            }
            else if(strcmp(eType,"HEXA")==0)
            {
                anElement.type= HEXA;
                int V1,V2,V3,V4,V5,V6,V7,V8;
                if(8!=fscanf(filein,"%d%d%d%d%d%d%d%d",&V1,&V2,&V3,&V4,&V5,&V6,&V7,&V8))
                {
                    cout<<"____error in reading the mesh data source____"<<endl;
                    return false;
                }
                anElement.theNodeIDs.push_back(V1);
                anElement.theNodeIDs.push_back(V2);
                anElement.theNodeIDs.push_back(V3);
                anElement.theNodeIDs.push_back(V4);
                anElement.theNodeIDs.push_back(V5);
                anElement.theNodeIDs.push_back(V6);
                anElement.theNodeIDs.push_back(V7);
                anElement.theNodeIDs.push_back(V8);
            }
            else if(strcmp(eType,"TET10")==0)
            {
                anElement.type= TET10;
                int V1,V2,V3,V4,V5,V6,V7,V8,V9,V10;
                if(10!=fscanf(filein,"%d%d%d%d%d%d%d%d%d%d",&V1,&V2,&V3,&V4,&V5,&V6,&V7,&V8,&V9,&V10))
                {
                    cout<<"____error in reading the mesh data source____"<<endl;
                    return false;
                }
                anElement.theNodeIDs.push_back(V1);
                anElement.theNodeIDs.push_back(V2);
                anElement.theNodeIDs.push_back(V3);
                anElement.theNodeIDs.push_back(V4);
                anElement.theNodeIDs.push_back(V5);
                anElement.theNodeIDs.push_back(V6);
                anElement.theNodeIDs.push_back(V7);
                anElement.theNodeIDs.push_back(V8);
                anElement.theNodeIDs.push_back(V9);
                anElement.theNodeIDs.push_back(V10);
            }
            elements.push_back(anElement);
        }

        fclose(filein);
        return true;
    }
    else return false;
}

//! ---------------------------------------------
//! function: compute normal at nodes
//! details:  compute an average normal at nodes
//! ---------------------------------------------
void Ng_MeshVS_DataSource2D::computeNormalAtNodes()
{
    //cout<<"Ng_MeshVS_DataSource2D::computeNormalAtNodes()->____function called____"<<endl;
    QMultiMap<int,int> nodeToElements;

    //! -----------------------------------------
    //! key => node number
    //! value => normal of the attached elements
    //! -----------------------------------------
    QMultiMap<int,QList<double>> nodeNormals;

    TColStd_MapIteratorOfPackedMapOfInteger eIt;
    for(eIt.Initialize(myElements);eIt.More();eIt.Next())
    {
        int localElementID = myElementsMap.FindIndex(eIt.Key());

        //! ---------------------------------------------------------------------
        //! the element normal
        //! Note: when using GetNormal, use theglobal element ID:
        //! as input, since this is converted internally to the local element ID
        //! ---------------------------------------------------------------------
        double ne_x = myElemNormals->Value(localElementID,1);
        double ne_y = myElemNormals->Value(localElementID,2);
        double ne_z = myElemNormals->Value(localElementID,3);

        QList<double> anElementNormal;
        anElementNormal<<ne_x<<ne_y<<ne_z;

        switch(myElemType->Value(localElementID))
        {
        case TRIG:
        {
            for(int col=1; col<=3; col++)
            {
                int globalNodeID = myElemNodes->Value(localElementID,col);
                int localNodeID = myNodesMap.FindIndex(globalNodeID);

                //! -------------------------------------------------
                //! build the QMultimap using the global node number
                //! -------------------------------------------------
                nodeNormals.insert(globalNodeID,anElementNormal);

                //! ---------------------------------------------------
                //! a choice for defining node to element connectivity
                //! using local numbering or global numbering
                //! ---------------------------------------------------
                nodeToElements.insert(localNodeID,localElementID);
            }
        }
            break;

            //! other cases ... to do
        }
    }

    //! -------------------------------------------
    //! organize node to element connectivity info
    //! the "nodeNumber" can be local or global
    //! according to the previous choice
    //! -------------------------------------------
    QList<int> keys = nodeToElements.keys();
    int curKey_old = -1;
    for(int i=0; i<keys.length(); i++)
    {
        int curKey = keys.at(i);
        if(curKey==curKey_old) continue;
        curKey_old = curKey;

        QList<int> elementNumbers = nodeToElements.values(curKey);
        myNodeToElements.insert(curKey,elementNumbers);
    }

    //! -------------------------------------------------------
    //! average the normal @ nodes
    //! "nodeNumber" can be global or global, depending on how
    //! the multimap has been generated
    //! -------------------------------------------------------
    QList<int> listOfKeys = nodeNormals.keys();

    int nodeNumber_old = -1;
    for(int i=0; i<listOfKeys.length(); i++)
    {
        int key = listOfKeys.at(i);
        int nodeNumber = key;
        if(nodeNumber==nodeNumber_old) continue;

        nodeNumber_old = nodeNumber;

        QList<QList<double>> normals = nodeNormals.values(nodeNumber);
        int NbNormalsPerNode = normals.length();
        //cout<<"____normal per node nr: "<<NbNormalsPerNode<<"____"<<endl;

        //! -----------------------------------------------------------
        //! in general there is more than one normal @ node
        //! first comb the normals (purge the almost identical,
        //! using the definition of elementeNormal::operator ==
        //! defined for struct Ng_MeshVS_DataSourceFace::elementNormal
        //! -----------------------------------------------------------
        double nn_x, nn_y, nn_z;
        nn_x = nn_y = nn_z = 0;

        QList<QList<double>> listOfElementNormals;

        //! --------------------
        //! inserting the first
        //! --------------------
        listOfElementNormals<<normals.at(0);
        for(int i=1; i<NbNormalsPerNode; i++)
        {
            QList<double> curNormal = normals.at(i);
            bool isDifferent = false;
            for(int k=0; k<listOfElementNormals.length(); k++)
            {
                QList<double> aPreviousNormal = listOfElementNormals.at(k);
                double cos_angle = curNormal.at(0)*aPreviousNormal.at(0)+curNormal.at(1)*aPreviousNormal.at(1)+curNormal.at(2)*aPreviousNormal.at(2);
                if(cos_angle<0.95) isDifferent = true;
                else isDifferent = false;
            }
            if(isDifferent == true) listOfElementNormals<<curNormal;
        }
        for(int i=0; i<listOfElementNormals.length(); i++)
        {
            nn_x = nn_x + listOfElementNormals.at(i).at(0);
            nn_y = nn_y + listOfElementNormals.at(i).at(1);
            nn_z = nn_z + listOfElementNormals.at(i).at(2);
        }

        //! ------------------------------------
        //! one way for calculating the normal
        //! check if suitable ... to do ...
        //! ------------------------------------
        double L = sqrt(pow(nn_x,2)+pow(nn_y,2)+pow(nn_z,2));
        if(L>1e-6) { nn_x = nn_x/L; nn_y = nn_y/L; nn_z = nn_z/L; }
        else { nn_x = nn_y = nn_z = 0.0; }

        QList<double> aveNormal;
        aveNormal<<nn_x<<nn_y<<nn_z;
        myNodeNormals.insert(nodeNumber, aveNormal);
    }
}

//! ----------------------------------
//! function: computeNormalAtElements
//! details:
//! ----------------------------------
void Ng_MeshVS_DataSource2D::computeNormalAtElements()
{
    //cout<<"Ng_MeshVS_DataSource2D::computeNormalAtElements()->____function called____"<<endl;
    for(TColStd_MapIteratorOfPackedMapOfInteger it(myElements); it.More(); it.Next())
    {
        std::vector<polygon::Point> points;

        int globalElementID = it.Key();
        int localElementID = myElementsMap.FindIndex(globalElementID);
        int NbNodes;
        double buf[24];
        TColStd_Array1OfReal coords(*buf,1,24);
        MeshVS_EntityType type;
        this->GetGeom(globalElementID,true,coords,NbNodes,type);
        for(int i=0; i<NbNodes; i++)
        {
            int s = 3*i;
            double x = coords(s+1);
            double y = coords(s+2);
            double z = coords(s+3);
            polygon::Point aPoint(x,y,z);
            points.push_back(aPoint);
        }
        const std::vector<double> &n = polygon::getNormal(points);
        myElemNormals->SetValue(localElementID,1,n[0]);
        myElemNormals->SetValue(localElementID,2,n[1]);
        myElemNormals->SetValue(localElementID,3,n[2]);
    }
}

//! -----------------------------
//! function: getNodeCoordinates
//! details:  helper
//! -----------------------------
std::vector<double> Ng_MeshVS_DataSource2D::getNodeCoordinates(int localNodeID)
{
    std::vector<double> coords;
    if(localNodeID<1 || localNodeID>myNumberOfNodes)
    {
        cerr<<"Ng_MeshVS_DataSource2D::getNodeCoordinates()->____localNodeID: "<<localNodeID<<" is out of range____"<<endl;
        return coords;
    }
    coords.push_back(myNodeCoords->Value(localNodeID,1));
    coords.push_back(myNodeCoords->Value(localNodeID,2));
    coords.push_back(myNodeCoords->Value(localNodeID,3));
    return coords;
}

//! ------------------------------------------------
//! function: computeNodeToElementsConnectivity
//! details:  compute node to elementd connectivity
//! ------------------------------------------------
void Ng_MeshVS_DataSource2D::computeNodeToElementsConnectivity()
{
    cout<<"Ng_MeshVS_DataSource2D::computeNodeToElementsConnectivity()->____function called____"<<endl;

    //! -------------------------------
    //! iterate over the mesh elements
    //! -------------------------------
    for(int localElementID=1; localElementID<=myNumberOfElements; localElementID++)
    {
        cout<<"@____element: "<<localElementID<<" of "<<myNumberOfElements<<endl;

        //! -------------------------------------
        //! iterate over the nodes of an element
        //! -------------------------------------
        int NbNodes;
        int eType = myElemType->Value(localElementID);
        cout<<"++++"<<endl;
        switch(eType)
        {
        case TRIG: NbNodes = 3; break;
        case QUAD: NbNodes = 4; break;
        case TRIG6: NbNodes = 6; break;
        case QUAD8: NbNodes = 8; break;
        }

        for(int k=1; k<=NbNodes; k++)
        {
            int globalNodeID = myElemNodes->Value(localElementID,k);
            cout<<" ("<<globalNodeID;

            int localNodeID = myNodesMap.FindIndex(globalNodeID);
            cout<<", "<<localNodeID<<") "<<endl;

            if(!myNodeToElements.contains(localNodeID))
            {
                QList<int> localElementIDs;
                localElementIDs<<localElementID;
                myNodeToElements.insert(localNodeID,localElementIDs);
            }
            else
            {
                QList<int> localElementIDs = myNodeToElements.value(localNodeID);
                localElementIDs<<localElementID;
                myNodeToElements.insert(localNodeID,localElementIDs);
            }
        }
    }
    cout<<"Ng_MeshVS_DataSource2D::computeNodeToElementsConnectivity()->____exiting function____"<<endl;

    /*
    //! ---------
    //! old code
    //! ---------
    QMultiMap<int,int> nodeToElements;
    for(TColStd_MapIteratorOfPackedMapOfInteger eIt(myElements);eIt.More();eIt.Next())
    {
        int localElementID = myElementsMap.FindIndex(eIt.Key());

        switch(myElemType->Value(localElementID))
        {
        case TRIG:
        {
            for(int col=1; col<=3; col++)
            {
                int globalNodeID = myElemNodes->Value(localElementID,col);
                int localNodeID = myNodesMap.FindIndex(globalNodeID);
                nodeToElements.insert(localNodeID,localElementID);
            }
        }
            break;
            //! other cases ... to do
        }
    }
    QList<int> keys = nodeToElements.keys();
    int curKey_old = -1;
    for(int i=0; i<keys.length(); i++)
    {
        int curKey = keys.at(i);
        if(curKey==curKey_old) continue;
        curKey_old = curKey;

        //! ---------------------------------------------------------
        //! retrieve from the QMultiMap the list of surface elements
        //! associated to the current key (current node)
        //! ---------------------------------------------------------
        QList<int> elementNumbers = nodeToElements.values(curKey);
        myNodeToElements.insert(curKey,elementNumbers);
    }
    */
    //! -----------
    //! diagnostic
    //! -----------
    /*
    FILE *f =fopen("D:/connectivityTable.txt","w");
    for(int i=1; i<=myNumberOfNodes; i++)
    {
        int localNodeID = i;
        QList<int> listOfElements = myNodeToElements.value(localNodeID);
        fprintf(f,"%d\t",localNodeID);
        for(int k=0; k<listOfElements.length()-2; k++)
            fprintf(f,"%d\t",listOfElements.at(k));
        fprintf(f,"%d\n",listOfElements.last());
    }
    fclose(f);
    */
}

//! ---------------------------
//! function: changeNodeCoords
//! details:
//! ---------------------------
bool Ng_MeshVS_DataSource2D::changeNodeCoords(int localNodeID, const std::vector<double> &coords)
{
    //! sanity check
    if(localNodeID<1 || localNodeID>myNumberOfNodes)
    {
        cerr<<"Ng_MeshVS_DataSource2D::changeNodeCoords()->____local node ID out of range when trying to change coordinates____"<<endl;
        return false;
    }
    if(coords.size()==0) return false;

    myNodeCoords->ChangeValue(localNodeID,1) = coords.at(0);
    myNodeCoords->ChangeValue(localNodeID,2) = coords.at(1);
    myNodeCoords->ChangeValue(localNodeID,3) = coords.at(2);
    return true;
}

//! -------------------------------------------------------------------
//! function: displaceMySelf
//! details:  displace the nodes using the "displacementField"
//!           the key of the input map should be the global element ID
//!           update normals at elements, at nodes, and 1D boundary
//! -------------------------------------------------------------------
void Ng_MeshVS_DataSource2D::displaceMySelf(const QMap<int, gp_Vec> &displacementField)
{
    for(QMap<int, gp_Vec>::const_iterator it = displacementField.cbegin(); it!= displacementField.cend(); ++it)
    {
        int globalNodeID = it.key();
        int localNodeID = this->myNodesMap.FindIndex(globalNodeID);
        const gp_Vec &curVec = it.value();
        double dX = curVec.X();
        double dY = curVec.Y();
        double dZ = curVec.Z();
        myNodeCoords->ChangeValue(localNodeID,1) = myNodeCoords->Value(localNodeID,1)+dX;
        myNodeCoords->ChangeValue(localNodeID,2) = myNodeCoords->Value(localNodeID,2)+dY;
        myNodeCoords->ChangeValue(localNodeID,3) = myNodeCoords->Value(localNodeID,3)+dZ;
    }
    if(myElemNormals->Length()==0) this->computeNormalAtElements();
    if(myNodeNormals.size()==0) this->computeNormalAtNodes();
    //if(myBoundarySegments.size()==0) this->computeFreeMeshSegments();
}

//! ------------------------------------------------------------------
//! function: buildTolerantPoints
//! details:  return the vector of all the mesh points with tolerance
//! ------------------------------------------------------------------
std::vector<mesh::tolerantPoint> Ng_MeshVS_DataSource2D::buildTolerantPoints(double tolerance)
{
    std::vector<mesh::tolerantPoint> points;
    for(int i=1; i<=myNumberOfNodes; i++)
    {
        double x = myNodeCoords->Value(i,1);
        double y = myNodeCoords->Value(i,2);
        double z = myNodeCoords->Value(i,3);
        mesh::tolerantPoint aP(x,y,z,tolerance);
        points.push_back((aP));
    }
    return points;
}

//! -------------------------------
//! function: getSegmentsOfElement
//! details:
//! -------------------------------
std::vector<mesh::meshSegment> Ng_MeshVS_DataSource2D::getSegmentsOfElement(int elementID, bool isLocal)
{
    int localElementID, globalElementID;
    if(isLocal)
    {
        localElementID = elementID;
        globalElementID = myElementsMap.FindIndex(elementID);
    }
    else
    {
        localElementID = myElementsMap.FindIndex(elementID);
        globalElementID = elementID;
    }

    std::vector<mesh::meshSegment> elementSegments;
    int NbNodes, buf[8];
    TColStd_Array1OfInteger nodeIDs(*buf,1,8);
    this->GetNodesByElement(globalElementID,nodeIDs,NbNodes);
    for(int i=0; i<NbNodes; i++)
    {
        int f = i%NbNodes+1;
        int s = (i+1)%NbNodes+1;
        mesh::meshSegment aSegment;
        unsigned int n1 = nodeIDs(f);
        unsigned int n2 = nodeIDs(s);
        if(n1<n2) aSegment.nodeIDs<<n1<<n2;
        else aSegment.nodeIDs<<n2<<n1;
        elementSegments.push_back(aSegment);
    }
    return elementSegments;
}

//! --------------------------------------------------------------
//! function: constructor
//! details:  convert a first order mesh into n-th order elements
//! --------------------------------------------------------------
Ng_MeshVS_DataSource2D::Ng_MeshVS_DataSource2D(const occHandle(Ng_MeshVS_DataSource2D) &aMesh, int order)
{
    if(order <=1 || order >2) return;

    //! ------------------------------------
    //! copy the nodes and the elements map
    //! ------------------------------------
    myNodes = aMesh->GetAllNodes();
    myNodesMap = aMesh->myNodesMap;
    myElements = aMesh->GetAllElements();
    myElementsMap = aMesh->myElementsMap;

    //! -----------------------------------
    //! the number of elements is the same
    //! -----------------------------------
    myNumberOfElements = aMesh->GetAllElements().Extent();
    myElemType = new TColStd_HArray1OfInteger(1,myNumberOfElements);
    myElemNodes = new TColStd_HArray2OfInteger(1,myNumberOfElements,1,8);
    myElemNormals = new TColStd_HArray2OfReal(1,myNumberOfElements,1,3);

    //! -----------------------------
    //! copy the element definitions
    //! -----------------------------
    for(TColStd_MapIteratorOfPackedMapOfInteger it(aMesh->GetAllElements()); it.More(); it.Next())
    {
        int globalElementID = it.Key();
        int localElementID = aMesh->myElementsMap.FindIndex(globalElementID);
        int NbNodes, buf[8];
        TColStd_Array1OfInteger nodeIDs(*buf,1,8);
        aMesh->GetNodesByElement(globalElementID,nodeIDs,NbNodes);
        for(int i=0; i<NbNodes; i++)
        {
            int j = i*2+1;
            myElemNodes->SetValue(localElementID,j,nodeIDs(i+1));
        }
    }

    //! ------------------------------
    //! meshSegment <=> point between
    //! ------------------------------
    int additionalNodes = 0;
    int globalNodeID_newNodes = -1;
    for(TColStd_MapIteratorOfPackedMapOfInteger it=aMesh->GetAllNodes(); it.More(); it.Next())
    {
        int akey = it.Key();
        if(akey>globalNodeID_newNodes) globalNodeID_newNodes=akey;
    }
    // do not fall into the temptation or replacing the previous code with
    // int globalNodeID_newNodes = aMesh->GetAllNodes().Extent();

    //! --------------------------
    //! iterate over the elements
    //! --------------------------
    myMidSideNodes = std::make_shared<std::vector<int>>();
    std::map<mesh::meshSegment,mesh::meshPoint> mapSegmentPointBetween;
    for(TColStd_MapIteratorOfPackedMapOfInteger it(aMesh->GetAllElements()); it.More(); it.Next())
    {
        int globalElementID = it.Key();
        const std::vector<mesh::meshSegment> &elementSegments = aMesh->getSegmentsOfElement(globalElementID,false);

        int NbSegments = int(elementSegments.size());
        for(int i=0; i<NbSegments; i++)
        {
            const mesh::meshSegment &curSegment = elementSegments[i];
            //! -------------------------------------------------
            //! compute the mid point of a segment only one time
            //! -------------------------------------------------
            if(mapSegmentPointBetween.count(curSegment)==0)
            {
                additionalNodes++;
                globalNodeID_newNodes++;

                //! --------------------------
                //! compute the point between
                //! --------------------------
                int globalNodeID_1 = curSegment.nodeIDs[0];
                int globalNodeID_2 = curSegment.nodeIDs[1];
                int localNodeID_1 = myNodesMap.FindIndex(globalNodeID_1);
                int localNodeID_2 = myNodesMap.FindIndex(globalNodeID_2);

                const std::vector<double> &P1 = aMesh->getNodeCoordinates(localNodeID_1);
                const std::vector<double> &P2 = aMesh->getNodeCoordinates(localNodeID_2);

                mesh::meshPoint mp(0.5*(P1[0]+P2[0]),0.5*(P1[1]+P2[1]),0.5*(P1[2]+P2[2]),globalNodeID_newNodes);
                std::pair<mesh::meshSegment,mesh::meshPoint> apair;
                apair.first = curSegment;
                apair.second = mp;
                mapSegmentPointBetween.insert(apair);

                //! ----------------------
                //! augment the node maps
                //! ----------------------
                myNodes.Add(globalNodeID_newNodes);
                myNodesMap.Add(globalNodeID_newNodes);
                myMidSideNodes->push_back(globalNodeID_newNodes);
            }
        }
    }

    //! -------------------------------------
    //! number of nodes and node coordinates
    //! -------------------------------------
    myNumberOfNodes = myNodes.Extent(); //aMesh->GetAllNodes().Extent() + additionalNodes;
    myNodeCoords = new TColStd_HArray2OfReal(1,myNumberOfNodes,1,3);

    //! --------------------------
    //! iterate over the elements
    //! --------------------------
    for(TColStd_MapIteratorOfPackedMapOfInteger it(aMesh->GetAllElements()); it.More(); it.Next())
    {
        int globalElementID = it.Key();
        const std::vector<mesh::meshSegment> &elementSegments = aMesh->getSegmentsOfElement(globalElementID,false);
        int NbSegments = int(elementSegments.size());
        for(int i=0; i<NbSegments; i++)
        {
            const mesh::meshSegment &curSegment = elementSegments[i];
            //! -------------------------------------------------
            //! compute the mid point of a segment only one time
            //! -------------------------------------------------
            if(mapSegmentPointBetween.count(curSegment)==0)
            {
                additionalNodes++;
                globalNodeID_newNodes++;

                //! --------------------------
                //! compute the point between
                //! --------------------------
                int globalNodeID_1 = curSegment.nodeIDs[0];
                int globalNodeID_2 = curSegment.nodeIDs[1];
                int localNodeID_1 = myNodesMap.FindIndex(globalNodeID_1);
                int localNodeID_2 = myNodesMap.FindIndex(globalNodeID_2);

                const std::vector<double> &P1 = aMesh->getNodeCoordinates(localNodeID_1);
                const std::vector<double> &P2 = aMesh->getNodeCoordinates(localNodeID_2);

                mesh::meshPoint mp(0.5*(P1[0]+P2[0]),0.5*(P1[1]+P2[1]),0.5*(P1[2]+P2[2]),globalNodeID_newNodes);
                std::pair<mesh::meshSegment,mesh::meshPoint> apair;
                apair.first = curSegment;
                apair.second = mp;
                mapSegmentPointBetween.insert(apair);

                //! ----------------------
                //! augment the node maps
                //! ----------------------
                myNodes.Add(globalNodeID_newNodes);
                myNodesMap.Add(globalNodeID_newNodes);
            }
        }
    }

    //! --------------------------
    //! copy the node coordinates
    //! --------------------------
    for(TColStd_MapIteratorOfPackedMapOfInteger it(aMesh->GetAllNodes());it.More();it.Next())
    {
        int globalNodeID = it.Key();
        int localNodeID = aMesh->myNodesMap.FindIndex(globalNodeID);
        const std::vector<double> &P = aMesh->getNodeCoordinates(localNodeID);
        myNodeCoords->SetValue(localNodeID,1,P[0]);
        myNodeCoords->SetValue(localNodeID,2,P[1]);
        myNodeCoords->SetValue(localNodeID,3,P[2]);
    }

    //! --------------------------
    //! iterate over the elements
    //! --------------------------
    for(TColStd_MapIteratorOfPackedMapOfInteger it(aMesh->GetAllElements()); it.More(); it.Next())
    {
        int globalElementID = it.Key();
        int localElementID = myElementsMap.FindIndex(globalElementID);

        //! ---------------------
        //! set the element type
        //! ---------------------
        ElemType eType;
        aMesh->GetElementType(eType,globalElementID,false);
        switch(eType)
        {
        case TRIG: myElemType->SetValue(localElementID,TRIG6); break;
        case QUAD: myElemType->SetValue(localElementID,QUAD8); break;
        }

        //! ------------------------------------------------------------
        //! add the node coordinates of the additional nodes
        //! complete the definition of each element with the mid points
        //! ------------------------------------------------------------
        const std::vector<mesh::meshSegment> &segmentsOfElement = aMesh->getSegmentsOfElement(globalElementID,false);
        for(int j=0;j<segmentsOfElement.size();j++)
        {
            mesh::meshSegment aseg = segmentsOfElement[j];
            std::map<mesh::meshSegment,mesh::meshPoint>::iterator si = mapSegmentPointBetween.find(aseg);
            std::pair<mesh::meshSegment,mesh::meshPoint> apair = *si;
            mesh::meshPoint aP = apair.second;
            int globalNodeID = aP.ID;
            int localNodeID = myNodesMap.FindIndex(globalNodeID);
            myNodeCoords->SetValue(localNodeID,1,aP.x);
            myNodeCoords->SetValue(localNodeID,2,aP.y);
            myNodeCoords->SetValue(localNodeID,3,aP.z);

            int k = 2*j+2;
            myElemNodes->SetValue(localElementID,k,globalNodeID);
        }
    }

    //! ----------------------------
    //! compute normals at elements
    //! ----------------------------
    this->computeNormalAtElements();
}

//! ------------------------
//! function: renumberNodes
//! details:
//! ------------------------
void Ng_MeshVS_DataSource2D::renumberNodes(const std::map<size_t,int> &mapCoordNodeNumbers,
                                           const std::vector<mesh::tolerantPoint> &vecTolPoints,
                                           const std::vector<int> &GlobalNodeIDs)
{
    int retrievedByProximity = 0;
    int retrievedByHash = 0;

    //! ------------------------------------------------
    //! fill a "parallel" node map: key = local node ID
    //! value = new global node ID
    //! ------------------------------------------------
    std::map<int,int> newGlobalNodeIDs;
    for(int localNodeID = 1; localNodeID<=myNumberOfNodes; localNodeID++)
    {
        double x = myNodeCoords->Value(localNodeID,1);
        double y = myNodeCoords->Value(localNodeID,2);
        double z = myNodeCoords->Value(localNodeID,3);

        //! ------------------------------------------
        //! trying to find the node using hash number
        //! ------------------------------------------
        std::size_t seed = 0;
        hash_c<double>(seed,x);
        hash_c<double>(seed,y);
        hash_c<double>(seed,z);
        std::map<size_t,int>::const_iterator it = mapCoordNodeNumbers.find(seed);
        if(it!=mapCoordNodeNumbers.cend())
        {
            std::pair<size_t,int> aPair = *it;
            int newGlobalNodeID = aPair.second;
            newGlobalNodeIDs.insert(std::make_pair(localNodeID,newGlobalNodeID));
            retrievedByHash++;
        }
        else
        {
            //! ----------------------------------------
            //! trying to find the node using proximity
            //! ----------------------------------------
            std::vector<mesh::tolerantPoint>::const_iterator it1 = std::find(vecTolPoints.cbegin(),vecTolPoints.cend(),mesh::tolerantPoint(x,y,z));
            if(it1!=vecTolPoints.cend())
            {
                int indexOf = std::distance(vecTolPoints.cbegin(),it1);
                int newGlobalNodeID = GlobalNodeIDs[indexOf];
                newGlobalNodeIDs.insert(std::make_pair(localNodeID,newGlobalNodeID));
                retrievedByProximity++;
            }
            else
            {
                //! ------------------------------------------
                //! this part of code should never be reached
                //! ------------------------------------------
                cerr<<"//---------------------------------------------------------------//"<<endl;
                cerr<<"// Ng_MeshVS_DataSource2D::renumberNodes()                       //"<<endl;
                cerr<<"// Cannot found ("<<x<<", "<<y<<", "<<z<<") hash: "<<seed<<endl;
                cerr<<"//---------------------------------------------------------------//"<<endl;
                exit(100);
            }
        }
    }

    cout<<"____retrieved by proximity: "<<double(retrievedByProximity)/double(myNumberOfNodes)*100.0<<" %____"<<endl;
    cout<<"____retrieved by hash: "<<double(retrievedByHash)/double(myNumberOfNodes)*100.0<<" %____"<<endl;
    if(retrievedByHash+retrievedByProximity==myNumberOfNodes) cout<<"____all nodes have been renumbered____"<<endl;

    //! --------------------------------------
    //! update the definition of the elements
    //! --------------------------------------
    for(int localElementID = 1; localElementID<=myNumberOfElements; localElementID++)
    {
        int NbNodes;
        switch(myElemType->Value(localElementID))
        {
        case TRIG: NbNodes = 3; break;
        case TRIG6: NbNodes = 6; break;
        case QUAD6: NbNodes = 6; break;
        case QUAD: NbNodes = 4; break;
        case QUAD8: NbNodes = 8; break;
        }
        for(int n=1; n<=NbNodes; n++)
        {
            int oldGlobalNodeID = myElemNodes->Value(localElementID,n);
            int localNodeID = myNodesMap.FindIndex(oldGlobalNodeID);

            std::map<int,int>::iterator it2 = newGlobalNodeIDs.find(localNodeID);
            int newGlobalNodeID = (*it2).second;
            myElemNodes->ChangeValue(localElementID,n) = newGlobalNodeID;
        }
    }

    //! ---------------------
    //! update the node maps
    //! ---------------------
    myNodesMap.Clear();
    myNodes.Clear();
    for(int localNodeID=1; localNodeID<=myNumberOfNodes; localNodeID++)
    {
        std::map<int,int>::iterator it = newGlobalNodeIDs.find(localNodeID);
        if(it!=newGlobalNodeIDs.end())
        {
            std::pair<int,int> apair = *it;
            int newGlobalNodeID = apair.second;
            myNodesMap.Add(newGlobalNodeID);
            myNodes.Add(newGlobalNodeID);
        }
        else
        {
            //! ------------------------------------------
            //! this part of code should never be reached
            //! ------------------------------------------
            cerr<<"//--------------------------------------------------------------//"<<endl;
            cerr<<"// Ng_MeshVS_DataSource2D::renumberNodes(): this part of code   //"<<endl;
            cerr<<"// shoud never be reached.                                      //"<<endl;
            cerr<<"//--------------------------------------------------------------//"<<endl;
            exit(102);
        }
    }
}
