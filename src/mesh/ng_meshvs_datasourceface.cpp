//#define LIBIGLNORMAL

//! -------------------------------------
//! macro for opencascade::handle<class>
//! -------------------------------------
#include "occhandle.h"

//! ----------------
//! custom includes
//! ----------------
#include <ng_meshvs_datasourceface.h>
#include <triangleaccessor.h>
#include <elementtypes.h>
#include <polygon.h>
#include "meshtools.h"
#include "tools.h"
#include "ccout.h"
#include <polygon.h>
#include <igtools.h>

//! ----
//! OCC
//! ----
#include <MeshVS_EntityType.hxx>
#include <TColStd_Array1OfInteger.hxx>
#include <TColStd_PackedMapOfInteger.hxx>
#include <TColStd_IndexedMapOfInteger.hxx>
#include <TColStd_MapIteratorOfPackedMapOfInteger.hxx>
#include <TColStd_HArray1OfInteger.hxx>
#include <TColStd_HArray2OfInteger.hxx>
#include <TColStd_HArray2OfReal.hxx>
#include <TColStd_MapIteratorOfMapOfInteger.hxx>

#include <TopTools_IndexedMapOfShape.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS.hxx>
#include <TopExp.hxx>
#include <TopAbs_ShapeEnum.hxx>
#include <gp_Vec.hxx>
#include <gp_Pnt.hxx>

//! ----
//! C++
//! ----
#include <iostream>
#include <vector>
#include <set>
#include <algorithm>

//! ---
//! Qt
//! ---
#include <QList>
#include <QMultiMap>

//! ------------
//! igl & Eigen
//! ------------
#include <libigl/include/igl/per_vertex_normals.h>
#include <libigl/include/igl/per_corner_normals.h>
#include <libigl/include/igl/per_face_normals.h>
#include <libigl/include/igl/gaussian_curvature.h>
#include <libigl/include/igl/invert_diag.h>
#include <libigl/include/igl/massmatrix.h>
#include <libigl/include/igl/grad.h>
#include <libigl/include/igl/cotmatrix.h>
#include <libigl/include/igl/principal_curvature.h>
#include <libigl/include/igl/ARAPEnergyType.h>
#include <libigl/include/igl/arap_linear_block.h>
#include <libigl/include/igl/arap.h>
#include <libigl/include/igl/harmonic.h>


//! 2D elements - nodal ordering
const int MTRIG3[3]={1,2,3};                // OK
const int MQUAD4[4]={1,2,3,4};              // OK
const int MTRIG6[6]={1,6,2,4,3,5};          // OK
const int MQUAD6[6]={1,2,3,4,5,6};          // da controllare
const int MQUAD8[8]={1,5,2,8,3,6,4,7};      // OK

const int MSURFELEM[5][10]={{1,2,3,0,0,0,0,0},      // TRIG
                            {1,2,3,4,0,0,0,0},      // QUAD
                            {1,6,2,4,3,5,0,0},     // TRIG6
                            {1,2,3,4,5,6,0,0},      // QUAD6    ?
                            {1,5,2,8,3,6,4,7}};     // QUAD8

IMPLEMENT_STANDARD_HANDLE(Ng_MeshVS_DataSourceFace,MeshVS_DataSource)
IMPLEMENT_STANDARD_RTTIEXT(Ng_MeshVS_DataSourceFace,MeshVS_DataSource)

//! ----------------------
//! function: constructor
//! details:  default
//! ----------------------
Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()
{
    ;
}

//! ------------------------------------------------------
//! function: constructor
//! details:  clone a mesh data source - copy constructor
//! ------------------------------------------------------
Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace(const occHandle(Ng_MeshVS_DataSourceFace) &aFaceMesh)
{
    //cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____cloning constructor called____"<<endl;

    //! -------------
    //! sanity check
    //! -------------
    if(aFaceMesh.IsNull()) return;

    //! -------------------------------------
    //! clone the maps of nodes and elements
    //! -------------------------------------
    myElements = aFaceMesh->GetAllElements();
    myNodes = aFaceMesh->GetAllNodes();
    myElementsMap = aFaceMesh->myElementsMap;
    myNodesMap = aFaceMesh->myNodesMap;

    myNumberOfElements = myElements.Extent();
    myNumberOfNodes = myNodes.Extent();
    if(myNumberOfElements==0 || myNumberOfNodes==0) return;

    //! ------------------
    //! set up the arrays
    //! ------------------
    myNodeCoords = new TColStd_HArray2OfReal(1,myNumberOfNodes,1,3);
    myElemNodes = new TColStd_HArray2OfInteger(1,myNumberOfElements,1,8);
    myElemType = new TColStd_HArray1OfInteger(1,myNumberOfElements);
    myElemNormals = new TColStd_HArray2OfReal(1,myNumberOfElements,1,3);

    TColStd_MapIteratorOfPackedMapOfInteger it;
    for(it.Initialize(myElements);it.More();it.Next())
    {
        int globalElementID = it.Key();
        int NbNodes;
        int bufn[8];
        TColStd_Array1OfInteger nodeIDs(*bufn,1,8);
        aFaceMesh->GetNodesByElement(globalElementID,nodeIDs,NbNodes);
        int localElementID = aFaceMesh->myElementsMap.FindIndex(globalElementID);

        for(int n=1; n<=NbNodes; n++)
        {
            myElemNodes->SetValue(localElementID,n,nodeIDs.Value(n));
        }

        switch(NbNodes)
        {
        case 3: myElemType->SetValue(localElementID,TRIG); break;
        case 4: myElemType->SetValue(localElementID,QUAD); break;
        case 6: myElemType->SetValue(localElementID,TRIG6); break;
        case 8: myElemType->SetValue(localElementID,QUAD8); break;
        default: break;
        }
    }

    for(it.Initialize(myNodes);it.More();it.Next())
    {
        int NbNodes;
        double buf[3];
        MeshVS_EntityType type;
        TColStd_Array1OfReal coords(*buf,1,3);

        int globalNodeID = it.Key();
        int localNodeID = aFaceMesh->myNodesMap.FindIndex(globalNodeID);
        aFaceMesh->GetGeom(globalNodeID,false,coords,NbNodes,type);
        myNodeCoords->SetValue(localNodeID,1,coords(1));
        myNodeCoords->SetValue(localNodeID,2,coords(2));
        myNodeCoords->SetValue(localNodeID,3,coords(3));
    }

    //! -------------------------------
    //! compute the normal at elements
    //! -------------------------------
    this->computeNormalAtElements();
}

//! ----------------------------------------------
//! function: constructor
//! details:  from a vector of elements and nodes
//! ----------------------------------------------
Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace(const std::vector<mesh::meshElement> &elements, std::vector<mesh::meshPoint> &nodes)
{
    myNumberOfElements = int(elements.size());
    myNumberOfNodes = int(nodes.size());

    //cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____source elements = "<<elements.size()<<", source nodes = "<<nodes.size()<<"____"<<endl;

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
    int localElementID = 0;
    for(std::vector<mesh::meshElement>::const_iterator it = elements.cbegin(); it!=elements.cend(); ++it)
    {
        const mesh::meshElement &anElement = *it;
        int globalElementID = anElement.ID;
        localElementID++;

        //! -------------------------
        //! fill the map of elements
        //! -------------------------
        myElements.Add(globalElementID);
        myElementsMap.Add(globalElementID);

        //! ----------------------------
        //! set the type of the element
        //! ----------------------------
        ElemType type = anElement.type;
        myElemType->SetValue(localElementID,type);

        //! -----------------------------
        //! set the nodes of the element
        //! -----------------------------
        switch(type)
        {
        case TRIG:
            for(int k=0; k<anElement.theNodeIDs.size(); k++)
            {
                int globalNodeID = anElement.theNodeIDs.at(k);
                myElemNodes->SetValue(localElementID,k+1,globalNodeID);
            }
            break;
        }
    }

    //cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____elements added____"<<endl;

    //! ----------------------------
    //! compute normals at elements
    //! ----------------------------
    this->computeNormalAtElements();

    //cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____exiting constructor____"<<endl;
}

//! -----------------------------------------------
//! function: constructor
//! details:  build the mesh data source from file
//! -----------------------------------------------
Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace(const std::string &fileName)
{
    //cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____constructor from file: ->"<<fileName<<"<-____"<<endl;
    std::vector<mesh::meshElement> elements;
    std::vector<mesh::meshPoint> nodes;
    bool isDone = this->readMesh(QString::fromStdString(fileName),nodes,elements);
    if(!isDone)
    {
        cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____cannot generate the face mesh for file: "<<fileName<<"____"<<endl;
        return;
    }
    myNumberOfElements = int(elements.size());
    myNumberOfNodes = int(nodes.size());

    if(myNumberOfElements<1)
    {
        cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____"<<fileName<<" no elements defined____"<<endl;
        return;
    }
    if(myNumberOfNodes<3)
    {
        cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____"<<fileName<<" wrong number of nodes____"<<endl;
        return;
    }

    myElemNodes = new TColStd_HArray2OfInteger(1,myNumberOfElements,1,8);
    myElemNormals = new TColStd_HArray2OfReal(1,myNumberOfElements,1,3);
    myElemType = new TColStd_HArray1OfInteger(1,myNumberOfElements);
    myNodeCoords = new TColStd_HArray2OfReal(1,myNumberOfNodes,1,3);

    //! ------------------------------------------
    //! fill the node map and add the coordinates
    //! ------------------------------------------
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
    int i=1;
    for(std::vector<mesh::meshElement>::iterator it = elements.begin(); it!=elements.end(); ++it, i++)
    {
        mesh::meshElement anElement = *it;
        int globalElementID = anElement.ID;

        //! -------------------------
        //! fill the map of elements
        //! -------------------------
        myElements.Add(globalElementID);
        myElementsMap.Add(globalElementID);

        //! ----------------------------
        //! set the type of the element
        //! ----------------------------
        ElemType type = anElement.type;
        myElemType->SetValue(i,type);

        //! -----------------------------
        //! set the nodes of the element
        //! -----------------------------
        switch(type)
        {
        case TRIG:
            for(int k=0; k<anElement.theNodeIDs.size(); k++)
            {
                int globalNodeID = anElement.theNodeIDs.at(k);
                myElemNodes->SetValue(i,k+1,globalNodeID);
            }
            break;

        case QUAD:
            for(int k=0; k<anElement.theNodeIDs.size(); k++)
            {
                int globalNodeID = anElement.theNodeIDs.at(k);
                myElemNodes->SetValue(i,k+1,globalNodeID);
            }
            break;

        case TRIG6:
            for(int k=0; k<anElement.theNodeIDs.size(); k++)
            {
                int globalNodeID = anElement.theNodeIDs.at(k);
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

//! ----------------------------------------------
//! function: constructor
//! details:  from Netgen pointer and face number
//! ----------------------------------------------
Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace(Ng_Mesh *aMesh, int faceNr)
{
    cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____constructor with Netgen pointer (face Nr. "<<faceNr<<") called____"<<endl;
    if(aMesh==NULL)
    {
        cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____NULL Netgen pointer. Exiting____"<<endl;
        return;
    }
    if(Ng_GetNP(aMesh)<3)
    {
        cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____The netgen surface mesh has less than 3 points. Exiting____"<<endl;
        return;
    }
    if(Ng_GetNSE(aMesh)<1)
    {
        cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____The mesh has no elements. Exiting____"<<endl;
        return;
    }

    //! --------------------------------------------------
    //! Ng_OCC_GetNSurfaceElementsByFaceNumber returns -1
    //! if the face has no elements
    //! --------------------------------------------------
    myNumberOfElements = Ng_OCC_GetNSurfaceElementsByFaceNumber(aMesh,faceNr);

    if(myNumberOfElements<1)
    {
        cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____number of elements: "<<myNumberOfElements<<". Exiting____"<<endl;
        return;
    }

    //! -------------------
    //! allocate the array
    //! -------------------
    int *sei = new int[myNumberOfElements];

    //! ---------------------------------------
    //! for a given face retrieve the elements
    //! sei[0], ..., sei[myNumberOfElements-1]
    //! ---------------------------------------
    int res = Ng_OCC_GetSurfaceElementsByFaceNumber(aMesh,faceNr,sei);

    if(res<=0)
    {
        cerr<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____the face has no elements____"<<endl;
        return;
    }

    //! ---------------------
    //! add the element list
    //! ---------------------
    for (int i=1;i<=myNumberOfElements;i++)
    {
        //! -----------------------------------------------------
        //! global element ID: en = sei[], surface element index
        //! -----------------------------------------------------
        int en = sei[i-1]+1;
        myElementsMap.Add(en);          //! (local, global)
        myElements.Add(en);             //! (local, global)
    }

    //! -------------------------------------------------------------------------------------
    //! retrieve the mesh nodes on the face: the resulting vector contains duplicated values
    //! since nodes are retrieved scanning surface elements
    //! -------------------------------------------------------------------------------------
    std::vector<int> nodeVector;
    for(int i=0; i<myNumberOfElements; i++)
    {
        int P[8];
        Ng_Surface_Element_Type type = Ng_GetSurfaceElementType(aMesh,sei[i]+1);
        switch(type)
        {
        case NG_TRIG:
            Ng_GetSurfaceElement(aMesh, sei[i]+1, P);
            for(int k=0;k<3;k++)nodeVector.push_back(P[k]);
            break;
        case NG_QUAD:
            Ng_GetSurfaceElement(aMesh, sei[i]+1, P);
            for(int k=0;k<4;k++)nodeVector.push_back(P[k]);
            break;
        case NG_TRIG6:
            Ng_GetSurfaceElement(aMesh, sei[i]+1, P);
            for(int k=0;k<6;k++)nodeVector.push_back(P[k]);
            break;
        case NG_QUAD8:
            Ng_GetSurfaceElement(aMesh, sei[i]+1, P);
            for(int k=0;k<8;k++)nodeVector.push_back(P[k]);
            break;
        }
    }

    //! -----------------
    //! clean duplicates
    //! -----------------
    std::set<int> B(nodeVector.begin(), nodeVector.end());
    std::vector<int> nodeVector_;
    for(int n=0;n<B.size();n++)
    {
        int x = *std::next(B.begin(), n);
        nodeVector_.push_back(x);
    }

    //! ---------------------------------------
    //! init the number of nodes od the face
    //! nodeVector_[0], ..., nodeVector[NNF-1]
    //! ---------------------------------------
    myNumberOfNodes = (int)nodeVector_.size();
    if(myNumberOfNodes<3)
    {
        cerr<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____number of nodes: "<<myNumberOfNodes<<". Exiting____"<<endl;
        return;
    }

    //! ----------------------------------
    //! allocate the space for the arrays
    //! ----------------------------------
    myElemNodes = new TColStd_HArray2OfInteger(1,myNumberOfElements,1,8);
    myElemNormals = new TColStd_HArray2OfReal(1,myNumberOfElements,1,3);
    myNodeCoords = new TColStd_HArray2OfReal(1,myNumberOfNodes,1,3);
    myElemType = new TColStd_HArray1OfInteger(1,myNumberOfElements);

    //! --------------------------------------------
    //! add the node list and the nodes coordinates
    //! --------------------------------------------
    double x_I[3];
    int i=1;
    for(std::vector<int>::iterator it = nodeVector_.begin(); it!=nodeVector_.end(); ++it)
    {
        //! -------------------------------------------
        //! "nodeID" is the global numbering
        //! "i" is the local numbering
        //! myNodesMap provides the link "ID" <=> "i"
        //! -------------------------------------------
        int nodeID = *it;
        myNodes.Add(nodeID);
        myNodesMap.Add(nodeID);
        Ng_GetPoint(aMesh,nodeID,x_I);
        myNodeCoords->SetValue(i,1,x_I[0]);
        myNodeCoords->SetValue(i,2,x_I[1]);
        myNodeCoords->SetValue(i,3,x_I[2]);
        i++;
    }

    //! ------------------------------------------------
    //! element definition through node numbers
    //! each element is defined through global node IDs
    //! ------------------------------------------------
    int p_I[8];
    for(int i=1; i<=myNumberOfElements;i++)
    {
        int en = sei[i-1]+1;
        Ng_Surface_Element_Type SE_Type = Ng_GetSurfaceElementType(aMesh,en);
        switch(SE_Type)
        {
        case NG_TRIG:
        {
            Ng_GetSurfaceElement(aMesh,en,p_I);
            myElemNodes->SetValue(i,1,p_I[0]);
            myElemNodes->SetValue(i,2,p_I[1]);
            myElemNodes->SetValue(i,3,p_I[2]);
            myElemType->SetValue(i,TRIG);
        }
            break;

        case NG_QUAD:
        {
            Ng_GetSurfaceElement(aMesh,en,p_I);
            myElemNodes->SetValue(i,1,p_I[0]);
            myElemNodes->SetValue(i,2,p_I[1]);
            myElemNodes->SetValue(i,3,p_I[2]);
            myElemNodes->SetValue(i,4,p_I[3]);
            myElemType->SetValue(i,QUAD);
        }
            break;

        case NG_TRIG6:
        {
            Ng_GetSurfaceElement(aMesh,en,p_I);
            myElemNodes->SetValue(i,1,p_I[0]);
            myElemNodes->SetValue(i,2,p_I[5]);
            myElemNodes->SetValue(i,3,p_I[1]);
            myElemNodes->SetValue(i,4,p_I[3]);
            myElemNodes->SetValue(i,5,p_I[2]);
            myElemNodes->SetValue(i,6,p_I[4]);
            myElemType->SetValue(i,TRIG6);
        }
            break;

        case NG_QUAD8:
        {
            Ng_GetSurfaceElement(aMesh,en,p_I);
            myElemNodes->SetValue(i,1,p_I[0]);
            myElemNodes->SetValue(i,2,p_I[4]);
            myElemNodes->SetValue(i,3,p_I[1]);
            myElemNodes->SetValue(i,4,p_I[7]);
            myElemNodes->SetValue(i,5,p_I[2]);
            myElemNodes->SetValue(i,6,p_I[5]);
            myElemNodes->SetValue(i,7,p_I[3]);
            myElemNodes->SetValue(i,8,p_I[6]);
            myElemType->SetValue(i,QUAD8);
        }
            break;
        }
    }

    //! ----------------------------
    //! compute normals at elements
    //! ----------------------------
    this->computeNormalAtElements();
}

//! ------------------------------------------
//! function: constructor
//! details:  from Tetgen i/o and face number
//! ------------------------------------------
Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace(tetgenio *aMesh,int faceNr)
{
    if(aMesh==NULL)
    {
        cerr<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____Tetgen i/o for face Nr. "<<faceNr<<" is null____"<<endl;
        return;
    }
    cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____Tetgen constructor called for face Nr. "<<faceNr<<"____"<<endl;

    //! ---------------------------------------
    //! number of triangles of the Tetgen mesh
    //! ---------------------------------------
    int N_TriFaces = aMesh->numberoftrifaces;
    std::vector<int> nodeVector;

    //! -----------------------
    //! "i" global element ID
    //! -----------------------
    for(int i=0; i<N_TriFaces; i++)
    {
        int faceMarker = aMesh->trifacemarkerlist[i];
        //! ---------------------------------
        //! the triangle belongs to the face
        //! ---------------------------------
        if(faceMarker==faceNr)
        {
            myElements.Add(i+1);      //! "i" global element ID
            myElementsMap.Add(i+1);   //! "i" global element ID
            for(int j=0;j<3;j++)
            {
                //! -------------------------------
                //! "nodeID" is the Tetgen node ID
                //! -------------------------------
                int nodeID = aMesh->trifacelist[i*3+j];
                nodeVector.push_back(nodeID);
            }
        }
    }
    myNumberOfElements = int(myElements.Extent());

    //! ------------------------------------------------
    //! eliminate the duplicates from the list of nodes
    //! ------------------------------------------------
    std::set<int> B(nodeVector.begin(), nodeVector.end());
    std::vector<int> nodeVector_;
    for(int n=0;n<B.size();n++)
    {
        int x = *std::next(B.begin(), n);
        nodeVector_.push_back(x);
    }

    if(nodeVector_.size()<3)
    {
        cerr<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->___the current face mesh: "<<
              faceNr<<" has less than 3 nodes. Exiting____"<<endl;
        return;
    }

    myNumberOfNodes = int(nodeVector_.size());
    myNodeCoords = new TColStd_HArray2OfReal(1,myNumberOfNodes,1,3);
    myElemType = new TColStd_HArray1OfInteger(1,myNumberOfElements);
    myElemNodes = new TColStd_HArray2OfInteger(1,myNumberOfElements,1,3);   // only linear triangles for the moment
    myElemNormals = new TColStd_HArray2OfReal(1,myNumberOfElements,1,3);

    //! ------------------------------------------------
    //! nodeVector_ is clean from the duplicated values
    //! nodeID is the Tetgen nodeID
    //! ------------------------------------------------
    for(int k=0; k<nodeVector_.size(); k++)
    {
        int globalNodeID = nodeVector_.at(k);

        //! ----------------------------------------
        //! as for constructor from Ng_Mesh pointer
        //! ----------------------------------------
        myNodes.Add(globalNodeID);
        myNodesMap.Add(globalNodeID);

        //! ----------------------------------------------------
        //! indexed map of integer: nodeID is the Tetgen nodeID
        //! ----------------------------------------------------
        for(int j=0;j<3;j++)
        {
            double coord = aMesh->pointlist[3*(globalNodeID-1)+j]; // -1 ???? // risposta: giusto: si veda la costruzione con "+1"
            myNodeCoords->SetValue(k+1,j+1,coord);
        }
    }

    //! ------------------------------------------------------------
    //! elements nodes and normals
    //! the face elements are defined through global mesh numbering
    //! "i" here is a local numbering
    //! ------------------------------------------------------------
    //const int C[3]{0,2,1};
    for(int i=1; i<=myNumberOfElements; i++)
    {
        myElemType->SetValue(i,TRIG);

        //! -------------------------------------------------
        //! element definition through global mesh numbering
        //! -------------------------------------------------
        int elementID = myElementsMap.FindKey(i);
        for(int k=0; k<3; k++)
        {
            //int k1 = C[k];
            int globalNodeID = aMesh->trifacelist[3*(elementID-1)+k];
            //myElemNodes->SetValue(i,k1+1,globalNodeID);
            myElemNodes->SetValue(i,k+1,globalNodeID);
        }
    }

    //! ----------------------------
    //! compute normals at elements
    //! ----------------------------
    this->computeNormalAtElements();
}

//! ----------------------------
//! function: constructor
//! details:  from StlMesh_Mesh
//! ----------------------------
Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace(const occHandle(StlMesh_Mesh) &anSTLMesh_Mesh,
                                                   const QList<QMap<int,int>> map_localToGlobal_nodeID,
                                                   const QList<QMap<int,int>> map_localToGlobal_elementID,
                                                   int faceNr)
{
    //cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____constructor from StlMesh_Mesh called for face nr. "<<faceNr<<"____"<<endl;
    if(anSTLMesh_Mesh.IsNull()) return;
    if(anSTLMesh_Mesh->NbTriangles()<1) return;
    if(anSTLMesh_Mesh->NbVertices()<1) return;

    const TColgp_SequenceOfXYZ& aCoords = anSTLMesh_Mesh->Vertices();
    Standard_Integer len = aCoords.Length(), j;
    if(len>0)
    {
        myNodeCoords = new TColStd_HArray2OfReal(1, len, 1, 3);
        myNumberOfNodes = len;

        gp_XYZ xyz;
        for(int i = 1; i <= len; i++)
        {
            myNodes.Add(map_localToGlobal_nodeID.at(faceNr).value(i));
            myNodesMap.Add(map_localToGlobal_nodeID.at(faceNr).value(i));

            xyz = aCoords(i);
            myNodeCoords->SetValue(i, 1, xyz.X());
            myNodeCoords->SetValue(i, 2, xyz.Y());
            myNodeCoords->SetValue(i, 3, xyz.Z());
        }

        const StlMesh_SequenceOfMeshTriangle& aSeq = anSTLMesh_Mesh->Triangles();
        len = aSeq.Length();
        myNumberOfElements = len;
        myElemNormals = new TColStd_HArray2OfReal(1, len, 1, 3);
        myElemNodes = new TColStd_HArray2OfInteger(1, len, 1, 3);
        myElemType = new TColStd_HArray1OfInteger(1,len);

        for(int i = 1; i <= len; i++)
        {
            int globalElementID = map_localToGlobal_elementID.at(faceNr).value(i);
            myElements.Add(globalElementID);
            myElementsMap.Add(globalElementID);

            myElemType->SetValue(i,TRIG);

            occHandle(StlMesh_MeshTriangle) aTriangle = aSeq.Value(i);
            int V[3];
            double nx, ny, nz;

            aTriangle->GetVertexAndOrientation(V[0], V[1], V[2], nx, ny, nz);

            for(j = 0; j < 3; j++)
            {
                int globalNodeID = map_localToGlobal_nodeID.at(faceNr).value(V[j]);
                //int globalNodeID = map_localToGlobal_nodeID.at(faceNr).value(j+1);
                //cout<<"____local node ID: "<<V[j]<<" global Node ID: "<<globalNodeID<<"____"<<endl;
                myElemNodes->SetValue(i, j+1, globalNodeID);
            }

            myElemNormals->SetValue(i,1,nx);
            myElemNormals->SetValue(i,2,ny);
            myElemNormals->SetValue(i,3,nz);
        }
    }
}

//! ----------------------------
//! function: constructor
//! details:  from StlMesh_Mesh
//! ----------------------------
Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace(const occHandle(StlMesh_Mesh) &aMesh)
{
    cout<<"Ng_MeshVS_DataSource2D::Ng_MeshVS_DataSourceFace()->____constructor called: surface mesh from StlMesh_Mesh____"<<endl;

    if(aMesh.IsNull() || aMesh->NbTriangles()<1 || aMesh->NbVertices()<3) return;
    myNumberOfNodes = aMesh->NbVertices();
    myNumberOfElements = aMesh->NbTriangles();
    const TColgp_SequenceOfXYZ& aCoords = aMesh->Vertices();
    int len = aCoords.Length(),i,j;
    myNodeCoords = new TColStd_HArray2OfReal(1,len,1,3);

    gp_XYZ xyz;

    for(i =1;i<=len;i++)
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

    for(i=1;i<=len;i++)
    {
        myElements.Add(i);
        myElementsMap.Add(i);
        myElemType->SetValue(i,TRIG);

        occHandle(StlMesh_MeshTriangle) aTriangle = aSeq.Value(i);
        int V[3];
        double nx, ny, nz;

        aTriangle->GetVertexAndOrientation(V[0],V[1],V[2],nx,ny,nz);

        for(j=0;j<3;j++) myElemNodes->SetValue(i,j+1,V[j]);

        myElemNormals->SetValue(i,1,nx);
        myElemNormals->SetValue(i,2,ny);
        myElemNormals->SetValue(i,3,nz);
    }
}

//! ---------------------------------
//! function: GetGeom
//! details:  access using global ID
//!----------------------------------
Standard_Boolean Ng_MeshVS_DataSourceFace::GetGeom(const Standard_Integer ID,
                                                   const Standard_Boolean IsElement,
                                                   TColStd_Array1OfReal& Coords,
                                                   Standard_Integer& NbNodes,
                                                   MeshVS_EntityType& Type) const
{
    if(IsElement)
    {
        //cout<<"Ng_MeshVS_DataSourceFace::GetGeom()->____function called for element ID: "<<ID<<"____"<<endl;
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
            return true;
        }
        else
        {
            cerr<<"Ng_MeshVS_DataSourceFace::GetGeom()->____error: element ID: "<<ID<<" out of range____"<<endl;
            return false;
        }
    }
    else
    {
        //cout<<"Ng_MeshVS_DataSourceFace::GetGeom()->____function called for node ID: "<<ID<<"____"<<endl;
        int localNodeID = myNodesMap.FindIndex(ID);
        if(localNodeID>=1 && localNodeID<=myNumberOfNodes)
        {
            Type = MeshVS_ET_Node;
            NbNodes=1;
            Coords(1) = myNodeCoords->Value(localNodeID,1);
            Coords(2) = myNodeCoords->Value(localNodeID,2);
            Coords(3) = myNodeCoords->Value(localNodeID,3);
            return true;
        }
        else
        {
            cerr<<"Ng_MeshVS_DataSourceFace::GetGeom()->____error: node "<<ID<<" out of range____"<<endl;
            return false;
        }
    }
}

//! -----------------------
//! function: GetGeomType
//! details:  ID is global
//! -----------------------
Standard_Boolean Ng_MeshVS_DataSourceFace::GetGeomType(const Standard_Integer ID,
                                                       const Standard_Boolean IsElement,
                                                       MeshVS_EntityType& Type) const
{
    if(IsElement)
    {
        int globalElementID = ID;
        int localElementID = myElementsMap.FindIndex(globalElementID);
        if(localElementID>=1 && localElementID<=myNumberOfElements)
        {
            Type = MeshVS_ET_Face;
            return true;
        }
        else
        {
            cerr<<"Ng_MeshVS_DataSourceFace::GetGeomType()->____error: element ID: "<<ID<<" out of range____"<<endl;
            return false;
        }
    }
    else
    {
        if(ID>=1 && ID<=myNumberOfNodes)
        {
            Type = MeshVS_ET_Node;
            return true;
        }
        else
        {
            cerr<<"Ng_MeshVS_DataSourceFace::GetGeomType()->____error: node ID: "<<ID<<" out of range____"<<endl;
            return false;
        }
    }
}

//! ------------------
//! function: GetAddr
//! details:
//! ------------------
Standard_Address Ng_MeshVS_DataSourceFace::GetAddr(const Standard_Integer, const Standard_Boolean) const
{
    return NULL;
}

//! ----------------------------
//! function: GetNodesByElement
//! details:
//! ----------------------------
Standard_Boolean Ng_MeshVS_DataSourceFace::GetNodesByElement(const Standard_Integer ID,
                                                             TColStd_Array1OfInteger& theNodeIDs,
                                                             Standard_Integer& theNbNodes) const
{
    //cout<<"Ng_MeshVS_DataSourceFace::GetNodesByElement()->____function called for element ID: "<<ID<<"____"<<endl;
    int localElementID = myElementsMap.FindIndex(ID);
    if(localElementID>=1 && localElementID<=myNumberOfElements)
    {
        switch(myElemType->Value(localElementID))
        {
        case(TRIG): theNbNodes=3; break;
        case(QUAD): theNbNodes=4; break;
        case(TRIG6): theNbNodes=6; break;
        case(QUAD8): theNbNodes=8; break;
        }
        for(int i=1;i<=theNbNodes;i++) theNodeIDs(i) = myElemNodes->Value(localElementID,i);
        return true;
    }
    cerr<<"Ng_MeshVS_DataSourceFace::GetNodesByElement()->____error: element ID: "<<ID<<" out of range____"<<endl;
    return false;
}

//! ----------------------
//! function: GetAllNodes
//! details:
//! ----------------------
const TColStd_PackedMapOfInteger& Ng_MeshVS_DataSourceFace::GetAllNodes() const
{
    return myNodes;
}

//! -------------------------
//! function: GetAllElements
//! details:
//! -------------------------
const TColStd_PackedMapOfInteger& Ng_MeshVS_DataSourceFace::GetAllElements() const
{
    return myElements;
}

//! -------------------------
//! function: GetElementType
//! details:
//! -------------------------
const bool Ng_MeshVS_DataSourceFace::GetElementType(ElemType &eType, int elementID, bool isLocal) const
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
Standard_Boolean Ng_MeshVS_DataSourceFace::GetNormal(const Standard_Integer Id,
                                                     const Standard_Integer Max,
                                                     Standard_Real& nx,
                                                     Standard_Real& ny,
                                                     Standard_Real& nz) const
{
    //cout<<"Ng_MeshVS_DataSourceFace::GetNormal()->____function called for element ID: "<<Id<<"____"<<endl;
    Q_UNUSED(Max)

    int localElementID = myElementsMap.FindIndex(Id);
    if(localElementID>=1 && localElementID<=myNumberOfElements && Max>=3)
    {
        nx = myElemNormals->Value(localElementID,1);
        ny = myElemNormals->Value(localElementID,2);
        nz = myElemNormals->Value(localElementID,3);
        return Standard_True;
    }
    else
    {
        cerr<<"Ng_MeshVS_DataSourceFace::GetNormal()->____error element ID: "<<Id<<" out of range____"<<endl;
        return Standard_False;
    }
}

//! -------------------
//! function: readMesh
//! details:
//! -------------------
bool Ng_MeshVS_DataSourceFace::readMesh(const QString &meshFileName,
                                        std::vector<mesh::meshPoint> &nodes,
                                        std::vector<mesh::meshElement> &elements)
{
    //cout<<"Ng_MeshVS_DataSourceFace::readMesh()->____"<<meshFileName.toStdString()<<"____"<<endl;
    FILE *filein = fopen(meshFileName.toStdString().c_str(),"r");
    if(filein!=NULL)
    {
        //! --------------------
        //! read mesh dimension
        //! --------------------
        int dim;
        if(1!=fscanf(filein,"%d",&dim))
        {
            cerr<<"Ng_MeshVS_DataSourceFace::readMesh()->____error in reading the mesh dimension for mesh file: "<<meshFileName.toStdString()<<"____"<<endl;
            return false;
        }

        //! -------------------------
        //! read the number of nodes
        //! -------------------------
        int NbNodes;
        if(1!=fscanf(filein,"%d",&NbNodes))
        {
            cout<<"Ng_MeshVS_DataSourceFace::readMesh()->____error in reading the number of nodes for mesh file: "<<meshFileName.toStdString()<<"____"<<endl;
            return false;
        }
        if(NbNodes<3) return false;
        //cout<<"____number of nodes: "<<NbNodes<<"____"<<endl;

        //! ---------------
        //! read the nodes
        //! ---------------
        for(int i=0; i<NbNodes;i++)
        {
            mesh::meshPoint aNode;

            if(4!=fscanf(filein,"%d%lf%lf%lf",&aNode.ID,&aNode.x,&aNode.y,&aNode.z))
            {
                cout<<"Ng_MeshVS_DataSourceFace::readMesh()->____error in reading nodes from mesh file: "<<meshFileName.toStdString()<<"____"<<endl;
                return false;
            }
            nodes.push_back(aNode);
        }
        //! ----------------------------
        //! read the number of elements
        //! ----------------------------
        int NbElements;
        if(1!=fscanf(filein,"%d",&NbElements))
        {
            cout<<"Ng_MeshVS_DataSourceFace::readMesh()->____error in reading the number of elements: "<<NbElements<<" from mesh file____"<<endl;
            return false;
        }
        if(NbElements<1) return false;
        //cout<<"____number of elements: "<<NbElements<<"____"<<endl;

        //! ------------------
        //! read the elements
        //! ------------------
        for(int i=0; i<NbElements; i++)
        {
            mesh::meshElement anElement;
            char eType[128];

            if(1!=fscanf(filein,"%d",&anElement.ID))
            {
                cout<<"Ng_MeshVS_DataSourceFace::readMesh()->____error in reading the element ID from mesh file: "<<meshFileName.toStdString()<<"____"<<endl;
                fclose(filein);
                return false;
            }
            if(1!=fscanf(filein,"%s",eType))
            {
                cout<<"Ng_MeshVS_DataSourceFace::readMesh()->____error in reading the element type from mesh file: "<<meshFileName.toStdString()<<"____"<<endl;
                fclose(filein);
                return false;
            }
            if(strcmp(eType,"TRIG")==0)
            {
                anElement.type= TRIG;
                int V1,V2,V3;
                if(3!=fscanf(filein,"%d%d%d",&V1,&V2,&V3))
                {
                    cout<<"Ng_MeshVS_DataSourceFace::readMesh()->____error in reading the node IDs for mesh file: "<<meshFileName.toStdString()<<"____"<<endl;
                    fclose(filein);
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
                    cout<<"Ng_MeshVS_DataSourceFace::readMesh()->____error in reading the node IDs for mesh file: "<<meshFileName.toStdString()<<"____"<<endl;
                    fclose(filein);
                    return false;
                }
                anElement.theNodeIDs.push_back(V1);
                anElement.theNodeIDs.push_back(V2);
                anElement.theNodeIDs.push_back(V3);
                anElement.theNodeIDs.push_back(V4);
            }
            if(strcmp(eType,"QUAD8")==0)
            {
                anElement.type= QUAD8;
                int V1,V2,V3,V4,V5,V6,V7,V8;
                if(8!=fscanf(filein,"%d%d%d%d%d%d%d%d",&V1,&V2,&V3,&V4,&V5,&V6,&V7,&V8))
                {
                    cout<<"Ng_MeshVS_DataSourceFace::readMesh()->____error in reading the node IDs for mesh file: "<<meshFileName.toStdString()<<"____"<<endl;
                    fclose(filein);
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
            if(strcmp(eType,"TRIG6")==0)
            {
                anElement.type= TRIG6;
                int V1,V2,V3,V4,V5,V6;
                if(6!=fscanf(filein,"%d%d%d%d%d%d",&V1,&V2,&V3,&V4,&V5,&V6))
                {
                    cout<<"Ng_MeshVS_DataSourceFace::readMesh()->____error in reading the node IDs for mesh file: "<<meshFileName.toStdString()<<"____"<<endl;
                    fclose(filein);
                    return false;
                }
                anElement.theNodeIDs.push_back(V1);
                anElement.theNodeIDs.push_back(V2);
                anElement.theNodeIDs.push_back(V3);
                anElement.theNodeIDs.push_back(V4);
                anElement.theNodeIDs.push_back(V5);
                anElement.theNodeIDs.push_back(V6);
            }
            elements.push_back(anElement);
        }
        //! ---------------
        //! close the file
        //! ---------------
        fclose(filein);
        //cout<<"Ng_MeshVS_DataSourceFace::readMesh()->____successfully exiting____"<<endl;
        return true;
    }
    else return false;
}

//!------------------------------------------------------------------------------
//! function: constructor
//! details:  "direct" construction from a list of points and a list of elements
//! -----------------------------------------------------------------------------
Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace(const QList<mesh::meshPoint> &meshPointList,
                                                   const QList<mesh::meshElement> &meshElementList)
{
    //cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____direct constructor called____"<<endl;

    myNumberOfNodes = meshPointList.length();
    myNumberOfElements = meshElementList.length();

    //cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____number of nodes: "<<myNumberOfNodes<<"____"<<endl;
    //cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____number of elements: "<<myNumberOfElements<<"____"<<endl;

    if(myNumberOfNodes<3) return;
    if(myNumberOfElements<1) return;

    myNodeCoords = new TColStd_HArray2OfReal(1,myNumberOfNodes,1,3);
    myElemNodes = new TColStd_HArray2OfInteger(1,myNumberOfElements,1,8);
    myElemNormals = new TColStd_HArray2OfReal(1,myNumberOfElements,1,3);
    myElemType = new TColStd_HArray1OfInteger(1,myNumberOfElements);

    //! -------------
    //! adding nodes
    //! -------------
    for(int i=1; i<=myNumberOfNodes; i++)
    {
        const mesh::meshPoint &curMeshPoint = meshPointList.at(i-1);
        int globalNodeID = curMeshPoint.ID;
        myNodes.Add(globalNodeID);
        myNodesMap.Add(globalNodeID);
        myNodeCoords->SetValue(i,1,curMeshPoint.x);
        myNodeCoords->SetValue(i,2,curMeshPoint.y);
        myNodeCoords->SetValue(i,3,curMeshPoint.z);
    }

    //! -------------
    //! add elements
    //! -------------
    for(int i=1; i<=myNumberOfElements; i++)
    {
        const mesh::meshElement &curMeshElement = meshElementList.at(i-1);
        int globalElementID = curMeshElement.ID;
        const std::vector<int> &nodeList = curMeshElement.theNodeIDs;
        switch(nodeList.size())
        {
        case 3: myElemType->SetValue(i,TRIG); break;
        case 4: myElemType->SetValue(i,QUAD); break;
        case 6: myElemType->SetValue(i,TRIG6); break;
        case 8: myElemType->SetValue(i,QUAD8); break;
        }

        myElements.Add(globalElementID);
        myElementsMap.Add(globalElementID);

        for(int k=0; k<nodeList.size(); k++)
        {
            int globalNodeID = curMeshElement.theNodeIDs.at(k);
            myElemNodes->SetValue(i,k+1,globalNodeID);
        }
    }

    //! ---------------------------
    //! compute normal at elements
    //! ---------------------------
    this->computeNormalAtElements();
}

//! ---------------------------------------------
//! function: compute normal at nodes
//! details:  compute an average normal at nodes
//! ---------------------------------------------
#ifndef LIBIGLNORMAL

void Ng_MeshVS_DataSourceFace::computeNormalAtNodes()
{
    //cout<<"Ng_MeshVS_DataSourceFace::computeNormalAtNodes()->____function called____"<<endl;

    //! ---------------------------------------------
    //! tolerance for normal comparisons: <n> degrees
    //! ---------------------------------------------
    const double tolerance = (10.0/180.0)*3.141592654;

    //! ---------------------------------------
    //! key: node number
    //! value: normal of the attached elements
    //! ---------------------------------------
    QMap<int,QList<mesh::elementNormal>> nodeNormals;
    for(TColStd_MapIteratorOfPackedMapOfInteger eIt(myElements);eIt.More();eIt.Next())
    {
        int globalElementID = eIt.Key();
        int localElementID = myElementsMap.FindIndex(globalElementID);

        //! -------------------
        //! the element normal
        //! -------------------
        double ne_x = myElemNormals->Value(localElementID,1);
        double ne_y = myElemNormals->Value(localElementID,2);
        double ne_z = myElemNormals->Value(localElementID,3);

        mesh::elementNormal anElementNormal(ne_x,ne_y,ne_z,tolerance);
        switch(myElemType->Value(localElementID))
        {
        case TRIG:
        {
            for(int col=1; col<=3; col++)
            {
                int globalNodeID = myElemNodes->Value(localElementID,col);
                int localNodeID = myNodesMap.FindIndex(globalNodeID);

                //! ---------------------------------------------------------
                //! build the map "nodeNormals" using the global node number
                //! ---------------------------------------------------------
                if(!nodeNormals.contains(globalNodeID))
                {
                    QList<mesh::elementNormal> normalList;
                    normalList<<anElementNormal;
                    nodeNormals.insert(globalNodeID,normalList);
                }
                else
                {
                    QList<mesh::elementNormal> normalList = nodeNormals[globalNodeID];
                    normalList<<anElementNormal;
                    nodeNormals.insert(globalNodeID,normalList);
                }
            }
        }
            break;
        }
    }

    //! -----------------------------
    //! average the normals at nodes
    //! -----------------------------
    for(QMap<int,QList<mesh::elementNormal>>::iterator it = nodeNormals.begin(); it!=nodeNormals.end(); it++)
    {
        int globalNodeID = it.key();
        QList<mesh::elementNormal> nodalNormals = it.value();
        std::set<mesh::elementNormal> setOfNodeNormals;

        //! -----------------------------------------------------------------------------
        //! use the insertion into a set for eliminating almost equally directed normals
        //! -----------------------------------------------------------------------------
        for(int j=0; j<nodalNormals.length(); j++)
        {
            setOfNodeNormals.insert(nodalNormals.at(j));
        }
        double ntotx = 0; double ntoty = 0; double ntotz = 0;
        for(std::set<mesh::elementNormal>::iterator it = setOfNodeNormals.begin(); it!=setOfNodeNormals.end(); it++)
        {
            ntotx += (*it).nx;
            ntoty += (*it).ny;
            ntotz += (*it).nz;
            //cout<<"____("<<ntotx<<", "<<ntoty<<", "<<ntotz<<")____"<<endl;
        }
        //! -------------------
        //! average the normal
        //! -------------------
        int NbNormals = int(setOfNodeNormals.size());
        ntotx /= NbNormals;
        ntoty /= NbNormals;
        ntotz /= NbNormals;
        double l = sqrt(ntotx*ntotx+ntoty*ntoty+ntotz*ntotz);
        if(l>1e-9)
        {
            ntotx /= l;
            ntoty /= l;
            ntotz /= l;
        }
        else { ntotx = ntoty = ntotz = 0.0; }
        QList<double> aveNormal {ntotx, ntoty, ntotz};
        myNodeNormals.insert(globalNodeID, aveNormal);
    }

    //! -----------------------------------------------------------------------------------
    //! the old version included the calculation of the node to elements connectivity
    //! for compatibility the new function computeNodeToElementsConnectivity is called here
    //! -----------------------------------------------------------------------------------
    this->computeNodeToElementsConnectivity();
}
#endif

#ifdef LIBIGLNORMAL
void Ng_MeshVS_DataSourceFace::computeNormalAtNodes()
{
    cout<<"Ng_MeshVS_DataSourceFace::computeNormalAtNodes()->____function called with igl library____"<<endl;

    //! -----------------------------
    //! node to element connectivity
    //! -----------------------------
    QMultiMap<int,int> nodeToElements;
    TColStd_MapIteratorOfPackedMapOfInteger eIt;
    for(eIt.Initialize(myElements);eIt.More();eIt.Next())
    {
        int localElementID = myElementsMap.FindIndex(eIt.Key());
        //cout<<"____local el id: "<<localElementID<<endl;
        switch(myElemType->Value(localElementID))
        {
        case TRIG:
        {
            for(int col=1; col<=3; col++)
            {
                int globalNodeID = myElemNodes->Value(localElementID,col);
                int localNodeID = myNodesMap.FindIndex(globalNodeID);
                //! ---------------------------------------------------
                //! a choice for defining node to element connectivity
                //! using local numbering or global numbering
                //! ---------------------------------------------------
                nodeToElements.insert(localNodeID,localElementID);
           }
        }
            break;

        default:
        {
            //! other elements
        }
            break;
        }
    }

    //! -----------------------------------------------
    //! re-organize node to element connectivity info
    //! the "nodeNumber" can be local or global
    //! according to the previous choice
    //! Here the key is the local node ID and the
    //! value is a list of local element IDs
    //! ----------------------------------------------
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

    //! -----------------------------
    //! compute the normals at nodes
    //! -----------------------------
    Eigen::MatrixXd V, N;
    Eigen::MatrixXi F;

    //! -----------------------------------------
    //! convert the current mesh into Eigen form
    //! -----------------------------------------
    iglTools::OCCMeshToIglMesh(this,V,F);
    igl::per_vertex_normals(V,F,N);

    for(int i=0; i<N.rows(); i++)
    {
        double nx = N(i,0);
        double ny = N(i,1);
        double nz = N(i,2);

        QList<double> normal;
        normal<<nx<<ny<<nz;
        int localNodeID = i+1;
        int globalNodeID = this->myNodesMap.FindKey(localNodeID);
        myNodeNormals.insert(globalNodeID,normal);
    }
}
#endif

//! ------------------------------------------------------------------------
//! function: computeFreeMeshSegments
//! details:  retrieve the mesh points on the 1D boundary of the face mesh
//!           The function builds also the segment to element connectivity
//! ------------------------------------------------------------------------
void Ng_MeshVS_DataSourceFace::computeFreeMeshSegments()
{
    //cout<<"Ng_MeshVS_DataSourceFace::computeFreeMeshSegments()->____function called____"<<endl;
    //! ------------------------------------
    //! key => segment; value => element ID
    //! ------------------------------------
    QMultiMap<mesh::meshSegment,int> segmentToElement;

    TColStd_PackedMapOfInteger eMap = this->GetAllElements();
    TColStd_MapIteratorOfPackedMapOfInteger eIt;
    for(eIt.Initialize(eMap);eIt.More();eIt.Next())
    {
        int globalElementID = eIt.Key();
        int localElementID = myElementsMap.FindIndex(globalElementID);
        int eType = myElemType->Value(localElementID);
        int NbNodes;
        switch (eType)
        {
        case TRIG: NbNodes = 3; break;
        case QUAD: NbNodes = 4; break;
        case TRIG6: NbNodes = 6; break;
        case QUAD8: NbNodes = 8; break;
        }

        //! ---------------------------------------------------------
        //! Important note:
        //! the mesh segment will be defined through global node IDs
        //! ---------------------------------------------------------
        QList<int> nodeIDs;
        for(int i=0; i<NbNodes; i++) nodeIDs<<myElemNodes->Value(localElementID,i+1);

        int type = myElemType->Value(localElementID);
        switch(type)
        {
        case TRIG: case QUAD:
        {
            //! ------------------------------------------------
            //! TRIG    {a,b,c} => {a,b}, {b,c}, {c,a}
            //! QUAD    {a,b,c,d} => {a,b}, {b,c}, {c,d}, {d,a}
            //! ------------------------------------------------
            int NbElementSegments = NbNodes;
            for(int i=1; i<=NbElementSegments; i++)
            {            
                mesh::meshSegment aMeshSegment;

                int firstIndex = (i-1)%NbNodes;
                int secondIndex = i%NbNodes;

                aMeshSegment.nodeIDs<<nodeIDs.at(firstIndex);
                aMeshSegment.nodeIDs<<nodeIDs.at(secondIndex);

                //! - panacea -
                aMeshSegment.sort();

                segmentToElement.insert(aMeshSegment,globalElementID);
                //cout<<"element ID: "<<eID<<"____Inserting segment ("<<aMeshSegment.nodeIDs.at(0)<<", "<<aMeshSegment.nodeIDs.at(1)<<")____"<<endl;
            }
        }
            break;

        case TRIG6: case QUAD8: //! MAYBE THIS IS NOT RIGHT... IMPLEMENT THE CORRECT VERSION. TO DO...
        {
            //! ----------------------------------------------------------------
            //! TRIG6   {a,b,c,d,e,f} => {a,b,c}, {c,d,e}, {e,f,a}
            //! QUAD8   {a,b,c,d,e,f,g,h} => {a,b,c}, {c,d,e}, {e,f,g}, {g,h,a}
            //! ----------------------------------------------------------------
            int NbElementSegments = nodeIDs.length()-2;
            for(int i=1; i<=NbElementSegments; i++)
            {
                mesh::meshSegment aMeshSegment;

                int firstIndex = (i-1)%NbNodes;
                int secondIndex = (i)%NbNodes;
                int thirdIndex = (i+1)%NbNodes;

                aMeshSegment.nodeIDs<<nodeIDs.at(firstIndex);
                aMeshSegment.nodeIDs<<nodeIDs.at(secondIndex);
                aMeshSegment.nodeIDs<<nodeIDs.at(thirdIndex);

                //! - panacea -
                aMeshSegment.sort();

                segmentToElement.insert(aMeshSegment,globalElementID);
            }
        }
            break;
        }
    }

    //! -------------------------------------------------------------------
    //! organize connectivity info: convert the QMultimap<meshSegment,int>
    //! into a QMap<meshSegment, QList<int>>
    //! -------------------------------------------------------------------
    //cout<<"Ng_MeshVS_DataSourceFace::computeFreeMeshSegments()->____building segment to element connectivity____"<<endl;
    QList<mesh::meshSegment> listOfKeys = segmentToElement.keys();
    mesh::meshSegment aMeshSegment_old;
    aMeshSegment_old.nodeIDs<<-999<<-998<<-997;

    for(int i=0; i<listOfKeys.length(); i++)
    {
        //Ng_MeshVS_DataSourceFace::meshSegment aMeshSegment = listOfKeys.at(i);
        mesh::meshSegment aMeshSegment = listOfKeys.at(i);
        if(aMeshSegment == aMeshSegment_old) continue;
        aMeshSegment_old = aMeshSegment;

        QList<int> adjacentElements = segmentToElement.values(aMeshSegment);

        //! --------------------------------
        //! segment to element connectivity
        //! --------------------------------
        mySegmentToElement.insert(aMeshSegment,adjacentElements);

        //! --------------------------------------------------------
        //! a segment is a boundary segment if it is connected only
        //! to a surface mesh element
        //! --------------------------------------------------------
        if(adjacentElements.length()==1)
        {
            if(!myBoundarySegments.contains(aMeshSegment))
            {
                //cout<<"BS____"<<aMeshSegment.nodeIDs.at(0)<<", "<<aMeshSegment.nodeIDs.at(1)<<"____"<<endl;
                myBoundarySegments<<aMeshSegment;
            }
        }
    }

    //! -----------------------------------------------
    //! fill the list of nodes on the 1D boundary
    //! Important note:
    //! "myBoundaryPoint" will contain global node IDs
    //! -----------------------------------------------
    for(int i=0; i<myBoundarySegments.length(); i++)
    {
        const mesh::meshSegment &aSegment = myBoundarySegments.at(i);
        for(int k=0; k<aSegment.nodeIDs.length(); k++)
        {
            int curNodeID = aSegment.nodeIDs.at(k);
            if(!myBoundaryPoints.contains(curNodeID)) myBoundaryPoints<<curNodeID;
        }
    }

    //! -----------------
    //! testing purposes
    //! -----------------
    //cout<<"____Number of boundary segments: "<<myBoundarySegments.length()<<"____"<<endl;
    //cout<<"____Number of boundary points: "<<myBoundaryPoints.length()<<"____"<<endl;
}

//! ----------------------------------
//! function: computeNormalAtElements
//! details:
//! ----------------------------------
#define USE_POLYGON
void Ng_MeshVS_DataSourceFace::computeNormalAtElements()
{
    //cout<<"Ng_MeshVS_DataSourceFace::computeNormalAtElements()->____function called____"<<endl;

#ifndef USE_POLYGON
    //! -----------------------------------------------------
    //! calculate the normal vector for each surface element
    //! for the moment this method handle triangles
    //! -----------------------------------------------------
    for(TColStd_MapIteratorOfPackedMapOfInteger it(myElements); it.More(); it.Next())
    {
        int globalElementID = it.Key();
        int localElementID = myElementsMap.FindIndex(globalElementID);

        double x[3],y[3],z[3];
        double nx,ny,nz,L;
        for(int n=1; n<=3; n++)
        {
            int globalNodeID = myElemNodes->Value(localElementID,n);
            int localNodeID = myNodesMap.FindIndex(globalNodeID);
            x[n-1] = myNodeCoords->Value(localNodeID,1);
            y[n-1] = myNodeCoords->Value(localNodeID,2);
            z[n-1] = myNodeCoords->Value(localNodeID,3);
        }

        //! ----------------------
        //!   i       j       k
        //! x1-x0   y1-y0   z1-z0
        //! x2-x0   y2-y0   z2-z0
        //! ----------------------
        nx=(y[1]-y[0])*(z[2]-z[0])-(z[1]-z[0])*(y[2]-y[0]);
        ny=(z[1]-z[0])*(x[2]-x[0])-(x[1]-x[0])*(z[2]-z[0]);
        nz=(x[1]-x[0])*(y[2]-y[0])-(y[1]-y[0])*(x[2]-x[0]);

        L = sqrt(pow(nx,2)+pow(ny,2)+pow(nz,2));
        if(L<1e-20) nx=ny=nz=0.0;
        else { nx=nx/L; ny=ny/L; nz=nz/L; }

        myElemNormals->SetValue(localElementID,1,nx);
        myElemNormals->SetValue(localElementID,2,ny);
        myElemNormals->SetValue(localElementID,3,nz);
    }
#endif
#ifdef USE_POLYGON

    for(TColStd_MapIteratorOfPackedMapOfInteger it(myElements); it.More(); it.Next())
    {
        std::vector<polygon::Point> aPolygon;

        int globalElementID = it.Key();
        int localElementID = myElementsMap.FindIndex(globalElementID);

        //cout<<"@____(global element ID, local element ID) = ("<<globalElementID<<", "<<localElementID<<")____"<<endl;

        int NbNodes;
        double buf[24];
        TColStd_Array1OfReal coords(*buf,1,24);
        MeshVS_EntityType type;
        bool isDone = this->GetGeom(globalElementID,true,coords,NbNodes,type);
        Q_UNUSED(isDone)

        for(int i=0; i<NbNodes; i++)
        {
            int s = 3*i;
            double x = coords(s+1);
            double y = coords(s+2);
            double z = coords(s+3);
            polygon::Point aPoint(x,y,z);
            aPolygon.push_back(aPoint);
            //cout<<"@____("<<x<<", "<<y<<", "<<z<<")____"<<endl;
        }
        std::vector<double> n = polygon::getNormal(aPolygon);
        //cout<<"@____normal ("<<n[0]<<", "<<n[1]<<", "<<n[2]<<")____"<<endl;
        myElemNormals->SetValue(localElementID,1,n[0]);
        myElemNormals->SetValue(localElementID,2,n[1]);
        myElemNormals->SetValue(localElementID,3,n[2]);
    }
#endif
}

//! -------------------------------------------------------------------
//! function: displaceMySelf
//! details:  displace the nodes using the "displacementField"
//!           the key of the input map should be the global element ID
//!           update normals at elements, at nodes, and 1D boundary
//! -------------------------------------------------------------------
void Ng_MeshVS_DataSourceFace::displaceMySelf(const QMap<int, gp_Vec> &displacementField)
{
    for(QMap<int, gp_Vec>::const_iterator it = displacementField.cbegin(); it!= displacementField.cend(); ++it)
    {
        int globalNodeID = it.key();
        int localNodeID = this->myNodesMap.FindIndex(globalNodeID);
        if(localNodeID<1 || localNodeID>myNumberOfNodes) continue;
        const gp_Vec &curVec = it.value();
        double dX = curVec.X();
        double dY = curVec.Y();
        double dZ = curVec.Z();
        myNodeCoords->ChangeValue(localNodeID,1) = myNodeCoords->Value(localNodeID,1)+dX;
        myNodeCoords->ChangeValue(localNodeID,2) = myNodeCoords->Value(localNodeID,2)+dY;
        myNodeCoords->ChangeValue(localNodeID,3) = myNodeCoords->Value(localNodeID,3)+dZ;
    }

    //! -----------------------------
    //! automatismi da rimuovere ...
    //! -----------------------------
    if(myElemNormals->Length()==0) this->computeNormalAtElements();
    if(myNodeNormals.size()==0) this->computeNormalAtNodes();
    if(myBoundarySegments.size()==0) this->computeFreeMeshSegments();
}

/*
//! ---------------------------------------------------------
//! function: constructor
//! details:  construction from tetgen .face and .node files
//! ---------------------------------------------------------
Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace(const QString &faceFileName,
                                                   const QString &nodeFileName,
                                                   int faceNr)
{
    //! --------------------------------------------------------------
    //! the nodes ot the whole 3D mesh are stored into the .node file
    //! --------------------------------------------------------------
    QMap<int,mesh::meshPoint> indexedMapOfAllMeshPoints;

    //! -----------------------------------------------------------
    //! the .node file contains the nodes of the whole volume mesh
    //! -----------------------------------------------------------
    FILE *nodeFile = fopen(nodeFileName.toStdString().c_str(),"r");
    if(nodeFile==NULL)
    {
        cerr<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____the .node file cannot be opened____"<<endl;
        return;
    }

    cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____Start reading node file: "<<nodeFileName.toStdString()<<"____"<<endl;

    //! ----------------------------------------------------------------------------------------
    //! First line: <# of points> <dimension (3)> <# of attributes> <boundary markers (0 or 1)>
    //! ----------------------------------------------------------------------------------------
    int NbPoints;
    int dimension;
    int NbAttributes;
    int boundaryMarkerTag;
    int pointNumber;
    double x,y,z;
    int tag;
    fscanf(nodeFile,"%d%d%d%d",&NbPoints,&dimension,&NbAttributes,&boundaryMarkerTag);

    //! -----------------
    //! console messages
    //! -----------------
    cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____Start reading mesh points___"<<endl;

    //! ------------------------------------
    //! "line" is actually the globalNodeID
    //! ------------------------------------
    for(int line=1; line<=NbPoints; line++)
    {
        fscanf(nodeFile,"%d%lf%lf%lf%d",&pointNumber,&x,&y,&z,&tag);
        mesh::meshPoint aMeshPoint(x,y,z,line);
        indexedMapOfAllMeshPoints.insert(line,aMeshPoint);
    }

    cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____mesh points read____"<<endl;
    fclose(nodeFile);

    //! ----------------------------------
    //! the face file must always be read
    //! ----------------------------------
    FILE *faceFile = fopen(faceFileName.toStdString().c_str(),"r");
    if(faceFile==NULL)
    {
        cerr<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____the .face file cannot be opened____"<<endl;
        return;
    }

    cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____Start reading the .face file: "<<faceFileName.toStdString()<<"____"<<endl;

    //! --------------------
    //! read the first line
    //! --------------------
    int NbElements;
    int boundaryTag;
    fscanf(faceFile,"%d%d",&NbElements,&boundaryTag);

    cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____Overall number of surface elements: "<<NbElements<<"____"<<endl;

    //! -----------------------------------------------------------------------------------------
    //! "line" is actually the globalElementID assigned by Tetgen to a mesh triangle.
    //! Thanks to Tetgen tags only the triangles of the current face are retrieved from the file
    //! -----------------------------------------------------------------------------------------
    int NbFaceElements = 0;
    for(int line=1; line<=NbElements; line++)
    {
        int number;
        int n1,n2,n3,n4,n5,n6;
        int faceTag;
        if(5==fscanf(faceFile,"%d%d%d%d%d",&number,&n1,&n2,&n3,&faceTag))
        {
            if(faceTag==faceNr) NbFaceElements++;
        }
        else if(8==fscanf(faceFile,"%d%d%d%d%d%d%d%d",&number,&n1,&n2,&n3,&n4,&n5,&n6,&faceTag))
        {
            if(faceTag==faceNr) NbFaceElements++;
        }
    }

    cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____face nr: "<<faceNr<<" number of elements: "<<NbFaceElements<<"____"<<endl;

    if(NbFaceElements<1)
    {
        cout<<"____face nr: "<<faceNr<<" has no elements____"<<endl;
        return;
    }

    //! ------------------------------------
    //! rewind to the beginning of the file
    //! and jump to the second line
    //! ------------------------------------
    rewind(faceFile);
    fscanf(faceFile,"%d%d",&NbElements,&boundaryTag);

    myNumberOfElements = NbFaceElements;
    myElemNodes = new TColStd_HArray2OfInteger(1,NbFaceElements,1,6);
    myElemType = new TColStd_HArray1OfInteger(1,NbFaceElements);
    myElemNormals = new TColStd_HArray2OfReal(1,NbFaceElements,1,3);

    QList<int> nodeIDsList;
    int localElementID = 0;

    for(int line=1; line<=NbElements; line++)
    {
        int globalElementID;
        int n1,n2,n3,n4,n5,n6;
        int faceTag;

        if(5==fscanf(faceFile,"%d%d%d%d%d",&globalElementID,&n1,&n2,&n3,&faceTag))
        {
            if(faceTag==faceNr)
            {
                localElementID++;
                myElemType->SetValue(localElementID,TRIG);

                myElementsMap.Add(globalElementID);
                myElements.Add(globalElementID);

                myElemNodes->SetValue(localElementID,1,n1);
                myElemNodes->SetValue(localElementID,2,n3);
                myElemNodes->SetValue(localElementID,3,n2);

                if(!nodeIDsList.contains(n1)) nodeIDsList<<n1;
                if(!nodeIDsList.contains(n2)) nodeIDsList<<n2;
                if(!nodeIDsList.contains(n3)) nodeIDsList<<n3;
            }
        }
        else if(8==fscanf(faceFile,"%d%d%d%d%d%d%d%d",&globalElementID,&n1,&n2,&n3,&n4,&n5,&n6,&faceTag))
        {
            if(faceTag==faceNr)
            {
                localElementID++;
                myElemType->SetValue(localElementID,TRIG6);

                myElementsMap.Add(globalElementID);
                myElements.Add(globalElementID);

                myElemNodes->SetValue(localElementID,1,n1);
                myElemNodes->SetValue(localElementID,2,n3);
                myElemNodes->SetValue(localElementID,3,n2);
                myElemNodes->SetValue(localElementID,4,n3);
                myElemNodes->SetValue(localElementID,5,n4);
                myElemNodes->SetValue(localElementID,6,n5);

                if(!nodeIDsList.contains(n1)) nodeIDsList<<n1;
                if(!nodeIDsList.contains(n2)) nodeIDsList<<n2;
                if(!nodeIDsList.contains(n3)) nodeIDsList<<n3;
                if(!nodeIDsList.contains(n4)) nodeIDsList<<n4;
                if(!nodeIDsList.contains(n5)) nodeIDsList<<n5;
                if(!nodeIDsList.contains(n6)) nodeIDsList<<n6;
            }
        }
    }

    cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____nodes added: "<<nodeIDsList.length()<<"____"<<endl;

    //! ---------------
    //! close the file
    //! ---------------
    fclose(faceFile);

    cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____number of elements for the face: "<<myNumberOfElements<<"____"<<endl;

    if(myNumberOfElements < 1)
    {
        cerr<<"____less than 1 element for face mesh data source: "<<faceNr<<". Exiting____"<<endl;
        return;
    }
    if(nodeIDsList.length()<3)
    {
        cerr<<"____less than 3 nodes for face mesh data source: "<<faceNr<<" (Number of nodes: "<<nodeIDsList.length()<<"). Exiting____"<<endl;
        return;
    }
    myNumberOfNodes = nodeIDsList.length();
    myNodeCoords = new TColStd_HArray2OfReal(1,myNumberOfNodes,1,3);

    for(int i=0; i<nodeIDsList.length(); i++)
    {
        int globalNodeID = nodeIDsList.at(i);
        myNodes.Add(globalNodeID);
        myNodesMap.Add(globalNodeID);
        const mesh::meshPoint &aP = indexedMapOfAllMeshPoints.value(globalNodeID);
        double x = aP.x;
        double y = aP.y;
        double z = aP.z;
        int localNodeID = i+1;
        myNodeCoords->SetValue(localNodeID,1,x);
        myNodeCoords->SetValue(localNodeID,2,y);
        myNodeCoords->SetValue(localNodeID,3,z);
    }
    this->computeNormalAtElements();
    //this->computeNormalAtNodes();         //! correct - not working
}
*/

//! ----------------------------------------------------------------------------
//! function: constructor
//! details:  construction from tetgen .face and indexed map of mesh::meshPoint
//!           this was introduced since when generating several face mesh data
//!           sources the .node file, using the previous constructor, should
//!           always be read, and it is slow in case of large meshes
//! ----------------------------------------------------------------------------
Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace(const QString &faceFileName,
                                                   const QMap<int,mesh::meshPoint> indexedMapOfAllMeshPoints,
                                                   int faceNr)
{
    //! ----------------------------------
    //! the face file must always be read
    //! ----------------------------------
    FILE *faceFile = fopen(faceFileName.toStdString().c_str(),"r");
    if(faceFile==NULL)
    {
        cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____the .face file cannot be opened____"<<endl;
        return;
    }

    cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____Start reading the .face file: "<<faceFileName.toStdString()<<"____"<<endl;

    //! --------------------
    //! read the first line
    //! --------------------
    int NbElements;
    int boundaryTag;
    fscanf(faceFile,"%d%d",&NbElements,&boundaryTag);

    cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____Overall number of surface elements: "<<NbElements<<"____"<<endl;

    //! -----------------------------------------------------------------------------------------
    //! "line" is actually the globalElementID assigned by Tetgen to a mesh triangle.
    //! Thanks to Tetgen tags only the triangles of the current face are retrieved from the file
    //! -----------------------------------------------------------------------------------------
    int NbFaceElements = 0;
    for(int line=1; line<=NbElements; line++)
    {
        int number;
        int n1,n2,n3,n4,n5,n6;
        int faceTag;
        if(5==fscanf(faceFile,"%d%d%d%d%d",&number,&n1,&n2,&n3,&faceTag))
        {
            if(faceTag==faceNr) NbFaceElements++;
        }
        else if(8==fscanf(faceFile,"%d%d%d%d%d%d%d%d",&number,&n1,&n2,&n3,&n4,&n5,&n6,&faceTag))
        {
            if(faceTag==faceNr) NbFaceElements++;
        }
    }

    cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____face nr: "<<faceNr<<" number of elements: "<<NbFaceElements<<"____"<<endl;

    if(NbFaceElements<1)
    {
        cout<<"____face nr: "<<faceNr<<" has no elements____"<<endl;
        return;
    }

    //! ------------------------------------
    //! rewind to the beginning of the file
    //! and jump to the second line
    //! ------------------------------------
    rewind(faceFile);
    fscanf(faceFile,"%d%d",&NbElements,&boundaryTag);

    myNumberOfElements = NbFaceElements;
    myElemNodes = new TColStd_HArray2OfInteger(1,NbFaceElements,1,6);
    myElemType = new TColStd_HArray1OfInteger(1,NbFaceElements);
    myElemNormals = new TColStd_HArray2OfReal(1,NbFaceElements,1,3);

    QList<int> nodeIDsList;
    int localElementID = 0;

    for(int line=1; line<=NbElements; line++)
    {
        int globalElementID;
        int n1,n2,n3,n4,n5,n6;
        int faceTag;

        if(5==fscanf(faceFile,"%d%d%d%d%d",&globalElementID,&n1,&n2,&n3,&faceTag))
        {
            if(faceTag==faceNr)
            {
                localElementID++;
                myElemType->SetValue(localElementID,TRIG);

                myElementsMap.Add(globalElementID);
                myElements.Add(globalElementID);

                myElemNodes->SetValue(localElementID,1,n1);
                myElemNodes->SetValue(localElementID,2,n3);
                myElemNodes->SetValue(localElementID,3,n2);

                if(!nodeIDsList.contains(n1)) nodeIDsList<<n1;
                if(!nodeIDsList.contains(n2)) nodeIDsList<<n2;
                if(!nodeIDsList.contains(n3)) nodeIDsList<<n3;
            }
        }
        else if(8==fscanf(faceFile,"%d%d%d%d%d%d%d%d",&globalElementID,&n1,&n2,&n3,&n4,&n5,&n6,&faceTag))
        {
            if(faceTag==faceNr)
            {
                localElementID++;
                myElemType->SetValue(localElementID,TRIG6);

                myElementsMap.Add(globalElementID);
                myElements.Add(globalElementID);

                myElemNodes->SetValue(localElementID,1,n1);
                myElemNodes->SetValue(localElementID,2,n3);
                myElemNodes->SetValue(localElementID,3,n2);
                myElemNodes->SetValue(localElementID,4,n3);
                myElemNodes->SetValue(localElementID,5,n4);
                myElemNodes->SetValue(localElementID,6,n5);

                if(!nodeIDsList.contains(n1)) nodeIDsList<<n1;
                if(!nodeIDsList.contains(n2)) nodeIDsList<<n2;
                if(!nodeIDsList.contains(n3)) nodeIDsList<<n3;
                if(!nodeIDsList.contains(n4)) nodeIDsList<<n4;
                if(!nodeIDsList.contains(n5)) nodeIDsList<<n5;
                if(!nodeIDsList.contains(n6)) nodeIDsList<<n6;
            }
        }
    }

    cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____nodes added: "<<nodeIDsList.length()<<"____"<<endl;

    //! ---------------
    //! close the file
    //! ---------------
    fclose(faceFile);

    cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____number of elements for the face: "<<myNumberOfElements<<"____"<<endl;

    if(myNumberOfElements < 1)
    {
        cout<<"____less than 1 element for face mesh data source: "<<faceNr<<". Exiting____"<<endl;
        return;
    }
    if(nodeIDsList.length()<3)
    {
        cout<<"____less than 3 nodes for face mesh data source: "<<faceNr<<" (Number of nodes: "<<nodeIDsList.length()<<"). Exiting____"<<endl;
        return;
    }
    myNumberOfNodes = nodeIDsList.length();
    myNodeCoords = new TColStd_HArray2OfReal(1,myNumberOfNodes,1,3);

    for(int i=0; i<nodeIDsList.length(); i++)
    {
        int globalNodeID = nodeIDsList.at(i);
        myNodes.Add(globalNodeID);
        myNodesMap.Add(globalNodeID);
        const mesh::meshPoint &aP = indexedMapOfAllMeshPoints.value(globalNodeID);
        double x = aP.x;
        double y = aP.y;
        double z = aP.z;
        int localNodeID = i+1;
        myNodeCoords->SetValue(localNodeID,1,x);
        myNodeCoords->SetValue(localNodeID,2,y);
        myNodeCoords->SetValue(localNodeID,3,z);
    }
    this->computeNormalAtElements();
    //this->computeNormalAtNodes();         //! correct - not working
}


//! -----------------------------------------------------------------
//! function: constructor
//! details:  an alternative constructor from Tetgen data: the nodes+
//!           and the elements are externally provided by reading
//!           the .node and .face files in advance
//! -----------------------------------------------------------------
Standard_EXPORT Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace(const QList<QList<double>> &listOfAllNodes,
                                                                   const QList<QList<int>> &listOfAllFaces,
                                                                   int faceNr)
{
    cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____constructor from Tetgen files called: faceNr: "<<
          faceNr<<" (using lists)____"<<endl;

    if(listOfAllNodes.length()<3)
    {
        cerr<<"____less than 3 nodes for face mesh data source: "<<faceNr<<"(Number of nodes: "<<listOfAllNodes.length()<<". Exiting____"<<endl;
        return;
    }
    if(listOfAllFaces.length())
    {
        cerr<<"___less than 1 triangle for face mesh data source: "<<faceNr<<". Exiting____"<<endl;
        return;
    }

    cout<<"____number of triangles into the overall mesh: "<<listOfAllFaces.length()<<"____"<<endl;
    cout<<"____number of points into the overall mesh: "<<listOfAllNodes.length()<<"____"<<endl;

    QMap<int,mesh::meshPoint> indexedMapOfAllMeshPoints;
    for(int i=0; i<listOfAllNodes.length(); i++)
    {
        double x = listOfAllNodes.at(i).at(0);
        double y = listOfAllNodes.at(i).at(1);
        double z = listOfAllNodes.at(i).at(2);

        mesh::meshPoint aP(x,y,z,i+1);
        indexedMapOfAllMeshPoints.insert(i+1,aP);
    }

    myNumberOfElements = listOfAllFaces.length();
    myElemNodes = new TColStd_HArray2OfInteger(1,myNumberOfElements,1,6);
    myElemType = new TColStd_HArray1OfInteger(1,myNumberOfElements);
    myElemNormals = new TColStd_HArray2OfReal(1,myNumberOfElements,1,3);

    QList<int> nodeIDsList;
   int localElementID = 0;
    for(int line=1; line<=myNumberOfElements; line++)
    {
        int globalElementID = line;
        int n1,n2,n3,n4,n5,n6;

        QList<int> se = listOfAllFaces.at(line-1);
        if(se.length()==5)
        {
            if(se.at(4)==faceNr)
            {
                localElementID++;

                //! ---------------------------------
                //! these are global node number IDs
                //! ---------------------------------
                n1 = se.at(1);
                n2 = se.at(2);
                n3 = se.at(3);

                myElemType->SetValue(localElementID,TRIG);
                myElementsMap.Add(globalElementID);
                myElements.Add(globalElementID);

                //! --------------------------------------
                //! check Tetgen nodal ordering ... to do
                //! --------------------------------------
                myElemNodes->SetValue(localElementID,1,n1);
                myElemNodes->SetValue(localElementID,2,n3);
                myElemNodes->SetValue(localElementID,3,n2);

                if(!nodeIDsList.contains(n1)) nodeIDsList<<n1;
                if(!nodeIDsList.contains(n2)) nodeIDsList<<n2;
                if(!nodeIDsList.contains(n3)) nodeIDsList<<n3;
            }
        }
        else if(se.length()==8)
        {
            if(se.at(7)==faceNr)
            {
                localElementID++;
                n1 = se.at(1);
                n2 = se.at(2);
                n3 = se.at(3);
                n4 = se.at(4);
                n5 = se.at(5);
                n6 = se.at(6);

                myElemType->SetValue(localElementID,TRIG6);

                myElementsMap.Add(globalElementID);
                myElements.Add(globalElementID);

                //! --------------------------------------
                //! check Tetgen nodal ordering ... to do
                //! --------------------------------------
                myElemNodes->SetValue(localElementID,1,n1);
                myElemNodes->SetValue(localElementID,2,n3);
                myElemNodes->SetValue(localElementID,3,n2);
                myElemNodes->SetValue(localElementID,4,n4);
                myElemNodes->SetValue(localElementID,5,n5);
                myElemNodes->SetValue(localElementID,6,n6);

                if(!nodeIDsList.contains(n1)) nodeIDsList<<n1;
                if(!nodeIDsList.contains(n2)) nodeIDsList<<n2;
                if(!nodeIDsList.contains(n3)) nodeIDsList<<n3;
                if(!nodeIDsList.contains(n4)) nodeIDsList<<n4;
                if(!nodeIDsList.contains(n5)) nodeIDsList<<n5;
                if(!nodeIDsList.contains(n6)) nodeIDsList<<n6;
            }
        }
    }

    if(nodeIDsList.length() < 3)
    {
        cout<<"____less than 3 node for face mesh data source: "<<faceNr<<" (Number of nodes: "<<nodeIDsList.length()<<"). Exiting____"<<endl;
        return;
    }

    myNumberOfNodes = nodeIDsList.length();
    myNodeCoords = new TColStd_HArray2OfReal(1,myNumberOfNodes,1,3);

    cout<<"____face nr: "<<faceNr<<" has "<<nodeIDsList.length()<<" nodes____"<<endl;

    for(int localNodeID=1; localNodeID<=myNumberOfNodes; localNodeID++)
    {
        cout<<"____working on node ID: "<<localNodeID<<"____"<<endl;
        int globalNodeID = nodeIDsList.at(localNodeID-1);

        myNodes.Add(globalNodeID);
        myNodesMap.Add(globalNodeID);

        const QList<double> &aPoint = listOfAllNodes.at(globalNodeID-1);
        cout<<"____OK node ID: "<<localNodeID<<"____"<<endl;

        double x = aPoint.at(0);
        double y = aPoint.at(1);
        double z = aPoint.at(2);

        myNodeCoords->SetValue(localNodeID,1,x);
        myNodeCoords->SetValue(localNodeID,2,y);
        myNodeCoords->SetValue(localNodeID,3,z);
    }

    this->computeNormalAtElements();
    this->computeNormalAtNodes();
    this->computeFreeMeshSegments();
    cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____constructor from Tetgen files (using lists) exiting____"<<endl;
}

//! ---------------------------------------------------------------
//! function: constructor
//! details:  build a mesh data source from a list of data sources
//!           the faces must belong to the same body
//! ---------------------------------------------------------------
Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace(const QList<occHandle(Ng_MeshVS_DataSourceFace)> &faceDSList)
{
    //cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____constructor from a list of "<<faceDSList.length()<<" face mesh data sources____"<<endl;

    //! ------------------
    //! fill the node map
    //! ------------------
    int NbSharedNodes = 0;
    for(int i=0; i<faceDSList.length(); i++)
    {
        const occHandle(Ng_MeshVS_DataSourceFace) &curFaceDS = faceDSList.at(i);
        if(curFaceDS.IsNull()) continue;

        //! -------------
        //! adding nodes
        //! -------------
        TColStd_PackedMapOfInteger curNodeMap = curFaceDS->GetAllNodes();
        for(TColStd_MapIteratorOfPackedMapOfInteger aNodeIt(curNodeMap);aNodeIt.More();aNodeIt.Next())
        {
            int globalNodeID = aNodeIt.Key();
            if(!myNodes.Contains(globalNodeID))
            {
                myNodes.Add(globalNodeID);
                myNodesMap.Add(globalNodeID);
            }
            else NbSharedNodes++;
        }

        //! ----------------
        //! adding elements
        //! ----------------
        TColStd_PackedMapOfInteger curEleMap = curFaceDS->GetAllElements();
        for(TColStd_MapIteratorOfPackedMapOfInteger anElemIt(curEleMap);anElemIt.More();anElemIt.Next())
        {
            int globalElementID = anElemIt.Key();

            //! ---------------------------------------------------------
            //! if the meshing method is "patch conforming" the faces
            //! do not share elements by definition. But if the meshing
            //! method is "patch independent" the faces can share pieces
            //! of the underlying mesh
            //! ---------------------------------------------------------
            if(myElements.Contains(globalElementID)==false)
            {
                myElements.Add(globalElementID);
                myElementsMap.Add(globalElementID);
            }
        }
    }

    myNumberOfNodes = myNodes.Extent();
    if(myNumberOfNodes<3)
    {
        cerr<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____constructor from a list. Wrong number of nodes: "<<myNumberOfNodes<<"____"<<endl;
        return;
    }
    myNumberOfElements = myElements.Extent();
    if(myNumberOfElements<1)
    {
        cerr<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____constructor from a list. Wrong number of elements: "<<myNumberOfElements<<"____"<<endl;
        return;
    }

    //! --------------------
    //! allocate the arrays
    //! --------------------
    myNodeCoords = new TColStd_HArray2OfReal(1,myNumberOfNodes,1,3);
    myElemType = new TColStd_HArray1OfInteger(1,myNumberOfElements);
    myElemNormals = new TColStd_HArray2OfReal(1,myNumberOfElements,1,3);
    myElemNodes = new TColStd_HArray2OfInteger(1,myNumberOfElements,1,8);

    //cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____constructor from a list. Number of shared nodes: "<<NbSharedNodes<<"____"<<endl;
    //cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____constructor from a list. Number of nodes: "<<myNumberOfNodes<<"____"<<endl;
    //cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____constructor from a list. Number of elements: "<<myNumberOfElements<<"____"<<endl;

    //! --------------------
    //! define the elements
    //! --------------------
    //cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____constructor from a list. Start defining the elements____"<<endl;

    //! ----------------------------------------------------------------
    //! note: I can not use QList for defining "alreadyVisitedElements"
    //! since the operator == for QList<mesh::meshElement> check == in
    //! a cyclic sense (is it slow? Is it true?)
    //! ----------------------------------------------------------------
    int NbSharedElements = 0;
    QMap<mesh::meshElement,int> alreadyVisitedElements;
    int localElementID=0;
    for(int i=0; i<faceDSList.length(); i++)
    {
        const occHandle(Ng_MeshVS_DataSourceFace) &curFaceDS = faceDSList.at(i);
        if(curFaceDS.IsNull()) continue;
        TColStd_PackedMapOfInteger curEleMap = curFaceDS->GetAllElements();
        for(TColStd_MapIteratorOfPackedMapOfInteger anElemIt(curEleMap);anElemIt.More();anElemIt.Next())
        {
            //! -----------------------------
            //! set the nodes of the element
            //! -----------------------------
            int globalElementID = anElemIt.Key();

            int NbNodes = 0;
            int buf[8];
            TColStd_Array1OfInteger nodeIDs(*buf,1,8);
            bool isDone = curFaceDS->GetNodesByElement(globalElementID,nodeIDs,NbNodes);
            Q_UNUSED(isDone)

            //! ----------------------------------------------------------
            //! face data sources could share elements (if the ds are not
            //! exact): store elements into a container for preventing
            //! an element from been inserted more than one time
            //! ----------------------------------------------------------
            mesh::meshElement aMeshElement;
            for(int n=1; n<=NbNodes; n++)
            {
                aMeshElement.theNodeIDs.push_back(nodeIDs(n));
            }

            if(alreadyVisitedElements.contains(aMeshElement))
            {
                cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____constructor from a list. A shared element has been found____"<<endl;
                NbSharedElements++;
                continue;
            }

            localElementID++;
            alreadyVisitedElements.insert(aMeshElement,localElementID);

            for(int k=1; k<=NbNodes; k++)
            {
                int globalNodeID = nodeIDs.Value(k);
                //!cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____constructor from a list: globalNodeID: "<<globalNodeID<<"____"<<endl;
                myElemNodes->SetValue(localElementID,k,globalNodeID);
            }

            //! ---------------------
            //! set the element type
            //! ---------------------
            switch(NbNodes)
            {
            case 3: myElemType->SetValue(localElementID,TRIG); break;
            case 4: myElemType->SetValue(localElementID,QUAD); break;
            case 6: myElemType->SetValue(localElementID,TRIG6); break;
            case 8: myElemType->SetValue(localElementID,QUAD8); break;
            }
        }
    }
    //cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____number of shared elements: "<<NbSharedElements<<"____"<<endl;
    //cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____constructor from a list. Elements defined____"<<endl;

    //! -----------------------------------------
    //! insert the coordinates of each mesh node
    //! -----------------------------------------
    //cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____constructor from a list. Start defining nodes____"<<endl;
    int localNodeID = 0;
    TColStd_PackedMapOfInteger alreadyVisitedNodes;
    for(int i=0; i<faceDSList.length(); i++)
    {
        const occHandle(Ng_MeshVS_DataSourceFace) &curFaceDS = faceDSList.at(i);
        if(curFaceDS.IsNull()) continue;

        TColStd_PackedMapOfInteger curNodeMap = curFaceDS->GetAllNodes();
        for(TColStd_MapIteratorOfPackedMapOfInteger aNodeIt(curNodeMap);aNodeIt.More();aNodeIt.Next())
        {
            int NbNodes;
            MeshVS_EntityType type;
            double buf[3];
            TColStd_Array1OfReal coords(*buf,1,3);

            int globalNodeID = aNodeIt.Key();
            if(alreadyVisitedNodes.Contains(globalNodeID)) continue;

            alreadyVisitedNodes.Add(globalNodeID);
            localNodeID++;
            curFaceDS->GetGeom(globalNodeID,false,coords,NbNodes,type);

            myNodeCoords->SetValue(localNodeID,1,coords(1));
            myNodeCoords->SetValue(localNodeID,2,coords(2));
            myNodeCoords->SetValue(localNodeID,3,coords(3));
        }
    }
    //cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____constructor from a list. Nodes defined____"<<endl;

    //! --------------------------------
    //! compute the normals at elements
    //! --------------------------------
    this->computeNormalAtElements();
    //cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____exiting____"<<endl;
}

//! ------------------------------
//! function: computeAnglesAtNode
//! details:
//! ------------------------------
void Ng_MeshVS_DataSourceFace::computeAnglesAtNode()
{
    for(int localNodeID = 1; localNodeID<=myNumberOfNodes; localNodeID++)
    {
        //! -------------------
        //! a node of the mesh
        //! -------------------
        int globalNodeID = myNodesMap.FindKey(localNodeID);

        QList<double> curNodeNormal = myNodeNormals.value(globalNodeID);

        //! -----------
        //! its normal
        //! -----------
        double nnx = curNodeNormal.at(0);
        double nny = curNodeNormal.at(1);
        double nnz = curNodeNormal.at(2);

        QList<double> listOfAnglesAtNode;

        //! --------------------------------------------------------
        //! find the surrounding nodes: use the node to element
        //! connectivity info, which are defines with local numbers
        //! --------------------------------------------------------
        QList<int> localElementIDs = myNodeToElements.value(localNodeID);
        for(int k=0; k<localElementIDs.length(); k++)
        {
            int LocalElementID_1 = localElementIDs.at(k);
            double nx1 = myElemNormals->Value(LocalElementID_1,1);
            double ny1 = myElemNormals->Value(LocalElementID_1,2);
            double nz1 = myElemNormals->Value(LocalElementID_1,3);

            double dot = nx1*nnx+ny1*nny+nz1*nnz;
            if(dot<-1.0) dot = -1.0;
            if(dot>1.0) dot = 1.0;
            double angle = acos(dot);
            listOfAnglesAtNode<<angle;
            //cout<<"____"<<angle*180.0/3.141592654<<"____"<<endl;
            myAnglesAtNode.insert(globalNodeID,listOfAnglesAtNode);
        }
    }
}

//! -------------------------------------------------------------------------------------------------
//! function: computeDiscreteCurvature
//! details:  https://github.com/libigl/libigl/blob/master/tutorial/203_CurvatureDirections/main.cpp
//! -------------------------------------------------------------------------------------------------
void Ng_MeshVS_DataSourceFace::computeDiscreteCurvature(int type)
{
    cout<<"Ng_MeshVS_DataSourceFace::computeDiscreteCurvature()->____function called____"<<endl;
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    iglTools::OCCMeshToIglMesh(this,V,F);

    Eigen::MatrixXd HN;
    Eigen::SparseMatrix<double> L,M,Minv;
    igl::cotmatrix(V,F,L);
    igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_VORONOI,M);
    igl::invert_diag(M,Minv);

    //! -----------------------------
    //! Laplace-Beltrami of position
    //! -----------------------------
    HN = -Minv*(L*V);

    //! ------------------------------------
    //! Extract magnitude as mean curvature
    //! ------------------------------------
    Eigen::VectorXd H = HN.rowwise().norm();
    Eigen::VectorXd K(H.rows());

    //! -------------------------------------------------
    //! Compute curvature directions via quadric fitting
    //! -------------------------------------------------
    Eigen::MatrixXd PD1,PD2;
    Eigen::VectorXd PV1,PV2;
    igl::principal_curvature(V,F,PD1,PD2,PV1,PV2);

    //! ----------------------------------------------
    //! mean curvature "H" and gaussian curvature "K"
    //! ----------------------------------------------
    H = 0.5*(PV1+PV2);
    for(int i=0; i<PV1.rows(); i++) K(i) = PV1(i)*PV2(i);

    //! --------------------------------------
    //! Compute gradient operator: #F*3 by #V
    //! --------------------------------------
    Eigen::SparseMatrix<double> G;
    igl::grad(V,F,G);

    switch(type)
    {
    //! -------------------------------------------
    //! insert the mean curvature and its gradient
    //! -------------------------------------------
    case 0:
    {
        //! -----------------------------------------------
        //! Compute gradient of the curvature scalar field
        //! (gradient magnitude)
        //! -----------------------------------------------
        Eigen::MatrixXd GH = Eigen::Map<const Eigen::MatrixXd>((G*H).eval().data(),F.rows(),3);

        for(int i=0; i<H.rows(); i++)
        {
            int globalNodeID = this->myNodesMap.FindKey(i+1);
            myCurvature.insert(globalNodeID,H(i));
            myCurvatureGradient.insert(globalNodeID,GH(i));
        }
    }
        break;
    //! -----------------------------------------------
    //! insert the gaussian curvature and its gradient
    //! -----------------------------------------------
    case 1:
    {
        //! -----------------------------------------------
        //! Compute gradient of the curvature scalar field
        //! (gradient magnitude)
        //! -----------------------------------------------
        Eigen::MatrixXd GK = Eigen::Map<const Eigen::MatrixXd>((G*K).eval().data(),F.rows(),3);

        for(int i=0; i<K.rows(); i++)
        {
            int globalNodeID = this->myNodesMap.FindKey(i+1);
            myCurvature.insert(globalNodeID,K(i));
            myCurvatureGradient.insert(globalNodeID,GK(i));
        }
    }
        break;
    }
}

//! -------------------------------------
//! function: constructor
//! details:  for the "Triangle" library
//! -------------------------------------
Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace(const QString &nodeFileName, const QString &faceFileName)
{
    cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____constructor for \"triangle\" called____"<<endl;
    char tmpName[256];
    sprintf(tmpName,"%s",nodeFileName.toStdString().c_str());
    FILE *nodeFile = fopen(tmpName,"r");
    if(nodeFile==NULL)
    {
        cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____the node file: "<<
              nodeFileName.toStdString()<<" is null____"<<endl;
        return;
    }
    sprintf(tmpName,"%s",faceFileName.toStdString().c_str());
    FILE *faceFile = fopen(tmpName,"r");
    if(faceFile==NULL)
    {
        cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____the face file is null____"<<endl;
        return;
    }

    //! --------------------
    //! read the .node file
    //! 4  2  0  1
    //! 1    -25  25    1
    //! 2    25  25    1
    //! 3    25  -25    1
    //! 4    -25  -25    1
    //! --------------------
    int nodeID;
    int NbNodes,dim,NbAttributes,bm;
    double x,y;
    fscanf(nodeFile,"%d%d%d%d",&NbNodes,&dim,&NbAttributes,&bm);

    myNumberOfNodes = NbNodes;
    myNodeCoords = new TColStd_HArray2OfReal(1,myNumberOfNodes,1,3);

    for(int i=1; i<=NbNodes; i++)
    {
        fscanf(nodeFile,"%d%lf%lf%d",&nodeID,&x,&y,&bm);
        myNodes.Add(i);
        myNodesMap.Add(i);
        myNodeCoords->SetValue(i,1,x);
        myNodeCoords->SetValue(i,2,y);
        myNodeCoords->SetValue(i,3,0.0);
    }
    fclose(nodeFile);

    //! --------------------------------------------------------------------
    //! read the .face file
    //! First line: <# of triangles> <nodes per triangle> <# of attributes>
    //! Remaining lines: <triangle #> <node> <node> <node> ... [attributes]
    //! 2  3  0
    //! 1       1     4     3
    //! 2       3     2     1
    //! --------------------------------------------------------------------
    int elementID;
    int NbTriangles;
    int a,i1,i2,i3;
    fscanf(faceFile,"%d\t%d\t%d\n",&NbTriangles,&a,&NbAttributes);
    myNumberOfElements = NbTriangles;

    myElemNodes = new TColStd_HArray2OfInteger(1,myNumberOfElements,1,8);
    myElemNormals = new TColStd_HArray2OfReal(1,myNumberOfElements,1,3);
    myElemType = new TColStd_HArray1OfInteger(1,myNumberOfElements);

    for(int i=1; i<=NbTriangles; i++)
    {
        fscanf(faceFile,"%d%d%d%d",&elementID,&i1,&i2,&i3);
        myElements.Add(i);
        myElementsMap.Add(i);
        myElemType->SetValue(i,TRIG);
        myElemNodes->SetValue(i,1,i1);
        myElemNodes->SetValue(i,2,i2);
        myElemNodes->SetValue(i,3,i3);
    }
    fclose(faceFile);

    //! --------------------------------
    //! compute the normals at elements
    //! --------------------------------
    this->computeNormalAtElements();
}

//! ------------------------------
//! function: constructor
//! details:  from Eigen matrices
//! ------------------------------
Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F)
{
    cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____constructor from Eigen matrices called____"<<endl;

    //! -------------
    //! sanity check
    //! -------------
    if(V.rows()<3)
    {
        cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____less than 3 nodes: exiting____"<<endl;
        return;
    }
    if(F.rows()<1)
    {
        cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____zero input triangles: exiting____"<<endl;
        return;
    }

    myNumberOfNodes = V.rows();
    myNodeCoords = new TColStd_HArray2OfReal(1,myNumberOfNodes,1,3);

    myNumberOfElements = F.rows();
    myElemNodes = new TColStd_HArray2OfInteger(1,myNumberOfElements,1,3);
    myElemNormals = new TColStd_HArray2OfReal(1,myNumberOfElements,1,3);
    myElemType = new TColStd_HArray1OfInteger(1,myNumberOfElements);

    for(int i=1; i<=myNumberOfNodes; i++)
    {
        myNodes.Add(i);
        myNodesMap.Add(i);
        int row = i-1;
        myNodeCoords->SetValue(i,1,V(row,0));
        myNodeCoords->SetValue(i,2,V(row,1));
        myNodeCoords->SetValue(i,3,V(row,2));
    }
    for(int i=1; i<=myNumberOfElements; i++)
    {
        myElementsMap.Add(i);
        myElements.Add(i);
        myElemType->SetValue(i,TRIG);
        int row = i-1;
        myElemNodes->SetValue(i,1,F(row,0));
        myElemNodes->SetValue(i,2,F(row,1));
        myElemNodes->SetValue(i,3,F(row,2));
    }

    //! ------------------------------------
    //! compute the normal at each triangle
    //! ------------------------------------
    this->computeNormalAtElements();
}

/*
//! ------------------------------------------------------------------
//! function: constructor
//! details:  input: a list of faces, the overall surfaceMesh
//!           output: the mesh triangles describing the list of faces
//!           this set of triangles is found by intersecting the .stl
//!           mesh of the input faces with the overall surface mesh
//! ------------------------------------------------------------------
#include <StlTransfer.hxx>
#include <igtools.h>
#include <libigl/include/igl/copyleft/cgal/intersect_other.h>
#include <TopoDS_Shell.hxx>
#include "topologytools.h"
#include <meshuvprojection.h>
Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace(const TopoDS_Face &aFace,
                                                   const occHandle(Ng_MeshVS_DataSourceFace) &aFaceMesh)
{
    cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____constructor from geometry face and surface mesh called____"<<endl;

    if(aFace.IsNull()) return;
    if(aFaceMesh.IsNull()) return;
    if(aFaceMesh->GetAllElements().Extent()<1) return;
    if(aFaceMesh->GetAllNodes().Extent()<3) return;

    //! -------------------------------------------
    //! convert the surface mesh into Eigen format
    //! -------------------------------------------
    Eigen::MatrixXd VS;
    Eigen::MatrixXi FS;
    iglTools::OCCMeshToIglMesh(aFaceMesh,VS,FS);

    //! -------------------------------------
    //! retrieve the .stl mesh from the face
    //! -------------------------------------
    occHandle(StlMesh_Mesh) stlMesh = new StlMesh_Mesh();
    StlTransfer::RetrieveMesh(aFace,stlMesh);

    if(stlMesh.IsNull()) return;
    if(stlMesh->NbTriangles()<1 || stlMesh->NbVertices()<3) return;

    //! ------------------------
    //! convert into OCC format
    //! ------------------------
    occHandle(Ng_MeshVS_DataSourceFace) curFaceMeshDS_coarse = new Ng_MeshVS_DataSourceFace(stlMesh);

    if(curFaceMeshDS_coarse.IsNull()) return;
    if(curFaceMeshDS_coarse->GetAllElements().Extent()<1 || curFaceMeshDS_coarse->GetAllNodes().Extent()<3) return;

    //! ------------------------
    //! convert into Eigen form
    //! ------------------------
    Eigen::MatrixXd VSS;
    Eigen::MatrixXi FSS;
    VSS.resize(curFaceMeshDS_coarse->GetAllNodes().Extent(),3);
    FSS.resize(curFaceMeshDS_coarse->GetAllElements().Extent(),3);
    iglTools::OCCMeshToIglMesh(curFaceMeshDS_coarse,VSS,FSS);

    //! --------------------------------------------------------------
    //! intersect the (coarse .stl) summed mesh with the surface mesh
    //! --------------------------------------------------------------
    Eigen::MatrixXi IF;             //! intersecting faces
    bool first_only = false;

    //! ---------------------------------------
    //! intersection within the physical space
    //! ---------------------------------------
    igl::copyleft::cgal::intersect_other(VS,FS,VSS,FSS,first_only,IF);
    int NbIntersectingTrianglePairs = IF.rows();

    //! ----------------------------------------------------------------
    //! intersection performed in the (u,v) space - experimental
    //! 1.  orthogonal projection of each point of the stl mesh onto
    //!     the surface the projection is computed into the (u,v) space
    //! ----------------------------------------------------------------

    //! -------------------------------------------------
    //! convert into Eigen format and perform projection
    //! -------------------------------------------------
    Eigen::MatrixXd V,VProj;
    Eigen::MatrixXi F;
    iglTools::OCCMeshToIglMesh(curFaceMeshDS_coarse,V,F);
    meshUVprojection::projectUV(aFace,V,VProj);

    //! 2.  project the computed mesh onto

    //! -----------------
    //! end experimental
    //! -----------------

    cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____intersection done____"<<endl;
    if(NbIntersectingTrianglePairs==0)
    {
        cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____intersecting triangles not found____"<<endl;
        return;
    }

    cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____number of intersecting pairs: "<<NbIntersectingTrianglePairs<<"____"<<endl;

    //! ------------------------------------
    //! definition of the intersection mesh
    //! ------------------------------------
    QList<int> intersectedTriangles;
    for(int n=0; n<NbIntersectingTrianglePairs; n++)
    {
        //! ----------------------------------------------------------
        //! "0" is the first column => triangles of the mesh (VS, FS)
        //! ----------------------------------------------------------
        int currentTriangle = IF(n,0);
        if(!intersectedTriangles.contains(currentTriangle)) intersectedTriangles<<currentTriangle;
    }

    cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____intersected triangles: "<<intersectedTriangles.length()<<"____"<<endl;

    //! -------------------
    //! number of elements
    //! -------------------
    myNumberOfElements = intersectedTriangles.length();
    myElemType = new TColStd_HArray1OfInteger(1,myNumberOfElements);
    myElemNodes = new TColStd_HArray2OfInteger(1,myNumberOfElements,1,8);     //! up to QUAD8
    myElemNormals = new TColStd_HArray2OfReal(1,myNumberOfElements,1,3);

    //! ------------------
    //! faces (triangles)
    //! ------------------
    QList<int> nodesIndices;
    for(int i=0; i<intersectedTriangles.length(); i++)
    {
        int pos = intersectedTriangles.at(i);
        int localElementID = pos+1;
        int globalElementID = aFaceMesh->myElementsMap.FindKey(localElementID);

        cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____triangle index: "<<pos<<" ";

        myElements.Add(globalElementID);
        myElementsMap.Add(globalElementID);

        int localNodeID1_inSurfaceMesh = FS(pos,0)+1;
        int localNodeID2_inSurfaceMesh = FS(pos,1)+1;
        int localNodeID3_inSurfaceMesh = FS(pos,2)+1;

        cout<<"("<<localNodeID1_inSurfaceMesh<<", "<<localNodeID2_inSurfaceMesh<<", "<<localNodeID3_inSurfaceMesh<<")"<<endl;

        int globalNodeID1_inSurfaceMesh = aFaceMesh->myNodesMap.FindKey(localNodeID1_inSurfaceMesh);
        if(!myNodes.Contains(globalNodeID1_inSurfaceMesh))
        {
            myNodes.Add(globalNodeID1_inSurfaceMesh);
            myNodesMap.Add(globalNodeID1_inSurfaceMesh);
            nodesIndices<<FS(pos,0);
        }
        int globalNodeID2_inSurfaceMesh = aFaceMesh->myNodesMap.FindKey(localNodeID2_inSurfaceMesh);
        if(!myNodes.Contains(globalNodeID2_inSurfaceMesh))
        {
            myNodes.Add(globalNodeID2_inSurfaceMesh);
            myNodesMap.Add(globalNodeID2_inSurfaceMesh);
            nodesIndices<<FS(pos,1);
        }
        int globalNodeID3_inSurfaceMesh = aFaceMesh->myNodesMap.FindKey(localNodeID3_inSurfaceMesh);
        if(!myNodes.Contains(globalNodeID3_inSurfaceMesh))
        {
            myNodes.Add(globalNodeID3_inSurfaceMesh);
            myNodesMap.Add(globalNodeID3_inSurfaceMesh);
            nodesIndices<<FS(pos,2);
        }

        //! -------------
        //! element type
        //! -------------
        myElemType->SetValue(i+1,TRIG);

        //! --------------
        //! element nodes
        //! --------------
        myElemNodes->SetValue(i+1,1,globalNodeID1_inSurfaceMesh);
        myElemNodes->SetValue(i+1,2,globalNodeID2_inSurfaceMesh);
        myElemNodes->SetValue(i+1,3,globalNodeID3_inSurfaceMesh);
    }

    //! ------
    //! nodes
    //! ------
    myNumberOfNodes = nodesIndices.length();
    myNodeCoords = new TColStd_HArray2OfReal(1,myNumberOfNodes,1,3);

    for(int i=0; i<nodesIndices.length(); i++)
    {
        int currentNodeIndex = nodesIndices.at(i);

        double x = VS(currentNodeIndex,0);
        double y = VS(currentNodeIndex,1);
        double z = VS(currentNodeIndex,2);

        myNodeCoords->SetValue(i+1,1,x);
        myNodeCoords->SetValue(i+1,2,y);
        myNodeCoords->SetValue(i+1,3,z);
    }

    //! ---------------------------
    //! compute normal at elements
    //! ---------------------------
    this->computeNormalAtElements();
}
*/

//! ---------------------------------------------------------------------------
//! function: constructor
//! details:  input: a TopoDS_Face face, and the overall surfaceMesh
//!           output: the mesh triangles of the overall surfaceMesh
//!           "overlapping" the input TopoDS_Face
//!
//! algo:    foreach surface element of the input mesh scan the nodes
//!
//!       1) Project the node onto the SURFACE underlying the face and
//!          compute the distance "d". If "d" is smaller than a tolerance "t".
//!          If the point is NOT CLOSE discard the element.
//!       2) if the point is CLOSE, classify the projection. If the projection
//!          does not belong to the face consider the point as NOT CLOSE,
//!          and discard the element
//!
//!          The face mesh data source is built using all the mesh elements
//!          which have not been discarded
//! ---------------------------------------------------------------------------
#include <TopAbs_State.hxx>
#include <BRepClass_FaceClassifier.hxx>
#include <BRep_Tool.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepGProp.hxx>
#include <Geom_Surface.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <GProp_GProps.hxx>
#include <polygon.h>
Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace(const TopoDS_Face &aFace,
                                                   const occHandle(MeshVS_DataSource) &surfaceMesh,
                                                   bool automaticProximity,
                                                   double proximityValue)
{
    //cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____constructor from geometry face and surface mesh called____"<<endl;
    int nn = 0;

    //! ---------------------------------------
    //! sanity check on the input surface mesh
    //! ---------------------------------------
    if(aFace.IsNull()) {
        cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____the input face is NULL___"<<endl;
        return;
    }
    if(surfaceMesh.IsNull()) {
        cerr<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____invalid input mesh____"<<endl;
        return;
    }
    if(surfaceMesh->GetAllElements().Extent()<1 || surfaceMesh->GetAllNodes().Extent()<3) {
        cerr<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____invalid input mesh. Nodes: "<<
              surfaceMesh->GetAllNodes().Extent()<<" Elements: "<<surfaceMesh->GetAllElements().Extent()<<"____"<<endl;
        return;
    }

    //! ----------------------------
    //! proximity value calculation
    //! ----------------------------
    if(automaticProximity==true)
    {
        //! area of the face
        double area = 0;
        GProp_GProps Props;
        try
        {
            BRepGProp::SurfaceProperties(aFace,Props);
            area = Props.Mass();
            proximityValue = 0.01*sqrt(area);
        }
        catch(...)
        {
            //! use the input value in case of an error
            cerr<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____BRepGProp::SurfaceProperties() error when"
                  "calculating the geometry face area____"<<endl;
            cerr<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____Using the proximity value provided"
                  "into the contructor____"<<endl;
        }
    }

    //! ------------------------------------------------
    //! a tool for projecting a point of the input mesh
    //! onto the input TopoDS_Face
    //! ------------------------------------------------
    GeomAPI_ProjectPointOnSurf aProjector;
    BRepAdaptor_Surface adaptor(aFace,true);
    const GeomAdaptor_Surface &s_adaptor = adaptor.Surface();
    opencascade::handle<Geom_Surface> aGeomSurface =  s_adaptor.Surface();

    //! -------------------------------------------------
    //! iterate over the elements of the input face mesh
    //! select which elements should be added to the
    //! face mesh according to the proximity from the
    //! input face
    //! -------------------------------------------------

    //! diagnostic file
    //FILE *f = fopen("D:/centers1.txt","w");

    bool useMeshElementCenter = false;
    for(TColStd_MapIteratorOfPackedMapOfInteger it(surfaceMesh->GetAllElements()); it.More(); it.Next())
    {
        int globalElementID = it.Key();

        int NbNodes;
        MeshVS_EntityType type;
        double buf[24];
        TColStd_Array1OfReal coords(*buf,1,24);
        bool isDone = surfaceMesh->GetGeom(globalElementID,true,coords,NbNodes,type);

        if(!isDone)
        {
            cerr<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____element globalElementID: "<<globalElementID<<" not found____"<<endl;
            continue;
        }
        bool isElementCloseToFace = true;

        //! -------------------------------
        //! use the points of the elements
        //! -------------------------------
        if(useMeshElementCenter==false)
        {
            int closePoints = 0;
            for(int k=0; k<NbNodes; k++)
            {
                //! -------------------------------------------------
                //! check the distance of the points of the element
                //! from the current face "aFace"
                //! -------------------------------------------------
                int s = 3*k;
                double xP = coords(s+1);
                double yP = coords(s+2);
                double zP = coords(s+3);

                aProjector.Init(gp_Pnt(xP,yP,zP),aGeomSurface);
                if(!aProjector.IsDone())
                {
                    //!cout<<"____point discarded because of projection failure____"<<endl;
                    //isElementCloseToFace = false;
                    //closePoints = 0;
                    continue;
                }

                const gp_Pnt &nst = aProjector.NearestPoint();
                double distance = sqrt(pow(xP-nst.X(),2)+pow(yP-nst.Y(),2)+pow(zP-nst.Z(),2));
                if(distance>proximityValue)
                {
                    //!cout<<"____point discarded because of distance criterion____"<<endl;
                    //isElementCloseToFace = false;
                    continue;
                }
                else
                {
                    //isElementCloseToFace = true;
                    //closePoints++;
                }

                //! -------------------------------------------------------------------
                //! https://www.opencascade.com/content/how-find-if-point-belongs-face
                //! -------------------------------------------------------------------
                BRepClass_FaceClassifier aFaceClassifier;
                double faceTolerance = 0.1;
                aFaceClassifier.Perform(aFace,nst,faceTolerance);
                TopAbs_State pointStatus = aFaceClassifier.State();
                if(pointStatus != TopAbs_IN && pointStatus !=TopAbs_ON)
                {
                    //!cout<<"____point discarded by face classifier____"<<endl;
                    //isElementCloseToFace = false;
                    //closePoints--;
                    continue;
                }

                //else
                //{
                    closePoints++;
                    //isElementCloseToFace = true;
                //}
                if(closePoints>=3)
                {
                    nn++;
                    myElements.Add(globalElementID);
                    myElementsMap.Add(globalElementID);
                }
            }
        }
        //! ----------------------------
        //! use only the element center
        //! ----------------------------
        else
        {
            //! ---------------------------------------------
            //! the polygon representing the surface element
            //! ---------------------------------------------
            std::vector<polygon::Point> aPolygon;
            for(int k=0; k<NbNodes; k++)
            {
                //! -------------------------------------------------
                //! check the distance of the points of the element
                //! from the current face "aFace"
                //! -------------------------------------------------
                int s = 3*k;
                double xP = coords(s+1);
                double yP = coords(s+2);
                double zP = coords(s+3);

                polygon::Point P; P.x = xP; P.y =yP; P.z = zP;
                aPolygon.push_back(P);
            }

            //! ------------------------------
            //! get the center of the element
            //! ------------------------------
            polygon::Point C;
            double area;
            polygon::getPolygonCenter(aPolygon,C,area);

            //! diagnostic file-
            //fprintf(f,"%lf\t%lf\t%lf\n",C.x,C.y,C.z);

            //! ------------------------------------------------
            //! project the center of the element onto the face
            //! ------------------------------------------------
            aProjector.Init(gp_Pnt(C.x,C.y,C.z),aGeomSurface);
            if(!aProjector.IsDone())
            {
                cerr<<"____point discarded because of intersection failure____"<<endl;
                isElementCloseToFace = false;
                continue;   //! jump over this surface element
            }

            gp_Pnt nst = aProjector.NearestPoint();
            double distance = sqrt(pow(C.x-nst.X(),2)+pow(C.y-nst.Y(),2)+pow(C.z-nst.Z(),2));
            if(distance>proximityValue)
            {
                cerr<<"____point discarded because of distance criterion. Distance: "<<distance<<"____"<<endl;
                isElementCloseToFace = false;
                continue;   //! jump over this surface element
            }
            else
            {
                isElementCloseToFace = true;
            }

            if(isElementCloseToFace)
            {
                //! -------------------------------------------------------------------
                //! https://www.opencascade.com/content/how-find-if-point-belongs-face
                //! -------------------------------------------------------------------
                BRepClass_FaceClassifier aFaceClassifier;
                double faceTolerance = 0.1;
                aFaceClassifier.Perform(aFace,nst,faceTolerance);
                TopAbs_State pointStatus = aFaceClassifier.State();
                if(pointStatus == TopAbs_IN /*|| pointStatus == TopAbs_ON*/)
                {
                    //cout<<"____point discarded by face classifier____"<<endl;
                    isElementCloseToFace = true;
                }
                else
                {
                    isElementCloseToFace = false;
                    continue;
                }

                if(isElementCloseToFace)
                {
                    nn++;
                    myElements.Add(globalElementID);
                    myElementsMap.Add(globalElementID);
                }
            }
        }
    }

    //! diagnostic file
    //fclose(f);

    myNumberOfElements = myElements.Extent();
    if(myElements.Extent()==0)
    {
        cerr<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____no element found____"<<endl;
        return;
    }
    //cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____number of elements: "<<myElements.Extent()<<"____"<<endl;

    //myNodeCoords = new TColStd_HArray2OfReal(1,myNumberOfNodes,1,3);
    myElemType = new TColStd_HArray1OfInteger(1,myNumberOfElements);
    myElemNodes = new TColStd_HArray2OfInteger(1,myNumberOfElements,1,8);     //! up to QUAD8
    myElemNormals = new TColStd_HArray2OfReal(1,myNumberOfElements,1,3);

    int localElementID = 0;
    for(TColStd_MapIteratorOfPackedMapOfInteger it(myElements);it.More();it.Next())
    {
        localElementID++;
        int globalElementID = it.Key();
        int NbNodes, buf[8];
        TColStd_Array1OfInteger nodeIDs(*buf,1,8);
        surfaceMesh->GetNodesByElement(globalElementID,nodeIDs,NbNodes);

        switch(NbNodes)
        {
        case 3: myElemType->SetValue(localElementID,TRIG); break;
        case 4: myElemType->SetValue(localElementID,QUAD); break;
        case 6: myElemType->SetValue(localElementID,TRIG6); break;
        case 8: myElemType->SetValue(localElementID,QUAD8); break;
        }
        for(int k=1; k<=NbNodes; k++)
        {
            int globalNodeID = nodeIDs.Value(k);
            myElemNodes->SetValue(localElementID,k,globalNodeID);

            if(!myNodes.Contains(globalNodeID))
            {
                myNodes.Add(globalNodeID);
                myNodesMap.Add(globalNodeID);
            }
        }
    }

    myNumberOfNodes = myNodes.Extent();
    if(myNumberOfNodes<3)
    {
        cerr<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____less than three nodes____"<<endl;
        return;
    }
    //cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____number of nodes: "<<myNodes.Extent()<<"____"<<endl;

    myNodeCoords = new TColStd_HArray2OfReal(1,myNumberOfNodes,1,3);

    for(TColStd_MapIteratorOfPackedMapOfInteger it(myNodes); it.More(); it.Next())
    {
        int globalNodeID = it.Key();
        int localNodeID = myNodesMap.FindIndex(globalNodeID);
        int NbNodes;
        double buf[3];
        MeshVS_EntityType type;
        TColStd_Array1OfReal coords(*buf,1,3);
        surfaceMesh->GetGeom(globalNodeID,false,coords,NbNodes,type);

        myNodeCoords->SetValue(localNodeID,1,coords(1));
        myNodeCoords->SetValue(localNodeID,2,coords(2));
        myNodeCoords->SetValue(localNodeID,3,coords(3));
    }

    //! ---------------------------
    //! compute normal at elements
    //! ---------------------------
    this->computeNormalAtElements();
}

//! -----------------------------
//! function: getNodeCoordinates
//! details:  helper
//! -----------------------------
std::vector<double> Ng_MeshVS_DataSourceFace::getNodeCoordinates(int localNodeID)
{
    std::vector<double> coords;
    if(localNodeID<1 || localNodeID>myNumberOfNodes)
    {
        cerr<<"Ng_MeshVS_DataSourceFace::getNodeCoordinates()->____localNodeID: "<<localNodeID<<" is out of range____"<<endl;
        return coords;
    }
    coords.push_back(myNodeCoords->Value(localNodeID,1));
    coords.push_back(myNodeCoords->Value(localNodeID,2));
    coords.push_back(myNodeCoords->Value(localNodeID,3));
    return coords;
}

//! ---------------------------
//! function: changeNodeCoords
//! details:
//! ---------------------------
bool Ng_MeshVS_DataSourceFace::changeNodeCoords(int localNodeID, const std::vector<double> &coords)
{
    //! sanity check
    if(localNodeID<1 || localNodeID>myNumberOfNodes) return false;
    if(coords.size()==0) return false;

    myNodeCoords->ChangeValue(localNodeID,1) = coords.at(0);
    myNodeCoords->ChangeValue(localNodeID,2) = coords.at(1);
    myNodeCoords->ChangeValue(localNodeID,3) = coords.at(2);
    return true;
}

//! -----------------------------------------------------------
//! function: toEigen
//! details:  helper - can be deleted since this functionality
//!           is already present in iglTools
//! -----------------------------------------------------------
void Ng_MeshVS_DataSourceFace::toEigen(Eigen::MatrixXd &V, Eigen::MatrixXi &F)
{
    cout<<"Ng_MeshVS_DataSourceFace::toEigen()->____function called____"<<endl;
    V.resize(myNumberOfNodes,3);
    for(int localNodeID = 1; localNodeID<=myNumberOfNodes; localNodeID++)
    {
        //! to Eigen format
        int pos = localNodeID-1;
        V(pos,0) = myNodeCoords->Value(localNodeID,1);
        V(pos,1) = myNodeCoords->Value(localNodeID,2);
        V(pos,2) = myNodeCoords->Value(localNodeID,3);
    }

    F.resize(myNumberOfElements,3);
    for(int localElementID =1; localElementID<myNumberOfElements; localElementID++)
    {
        int NbNodes = 3;
        int pos = localElementID-1;
        for(int i=1; i<=NbNodes; i++)
        {
            int curLocalNodeID = myElemNodes->Value(localElementID,i);
            //! to Eigen format
            F(pos,i-1) = curLocalNodeID-1;
            //F(pos,i-1) = curLocalNodeID;
        }
    }
    cout<<"Ng_MeshVS_DataSourceFace::toEigen()->____to eigen conversion done____"<<endl;
}

//! -------------------------------------------------------------------------
//! function: constructor
//! details:  questo è un accessorio per ora indispensabile, ma da rimuovere
//! -------------------------------------------------------------------------
Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace(const occHandle(Ng_MeshVS_DataSource2D) &surfaceMesh)
{
    cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____conversion from overall mesh 2D____"<<endl;

    if(surfaceMesh.IsNull()) return;
    if(surfaceMesh->GetAllElements().Extent()<1) return;
    if(surfaceMesh->GetAllNodes().Extent()<3) return;

    //! ------
    //! nodes
    //! ------
    myNumberOfNodes = surfaceMesh->GetAllNodes().Extent();
    myNodeCoords = new TColStd_HArray2OfReal(1,myNumberOfNodes,1,3);

    //! ---------
    //! elements
    //! ---------
    myNumberOfElements = surfaceMesh->GetAllElements().Extent();
    myElemNodes = new TColStd_HArray2OfInteger(1,myNumberOfElements,1,8);
    myElemNormals = new TColStd_HArray2OfReal(1,myNumberOfElements,1,3);
    myElemType = new TColStd_HArray1OfInteger(1,myNumberOfElements);

    cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____arrays OK____"<<endl;

    cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____adding elements____"<<endl;
    TColStd_MapIteratorOfPackedMapOfInteger it;
    int localElementID = 0;
    for(it.Initialize(surfaceMesh->GetAllElements()); it.More(); it.Next())
    {
        localElementID++;
        int globalElementID = it.Key();
        myElements.Add(globalElementID);
        myElementsMap.Add(globalElementID);

        int NbNodes;
        int buf[8];
        TColStd_Array1OfInteger nodeIDs(*buf,1,8);
        surfaceMesh->GetNodesByElement(globalElementID,nodeIDs,NbNodes);

        for(int i=1; i<=NbNodes; i++)
        {
            int globalNodeID = nodeIDs.Value(i);
            myElemNodes->SetValue(localElementID,i,globalNodeID);
            switch(NbNodes)
            {
            case 3: myElemType->SetValue(localElementID,TRIG); break;
            case 4: myElemType->SetValue(localElementID,QUAD); break;
            case 6: myElemType->SetValue(localElementID,TRIG6); break;
            case 8: myElemType->SetValue(localElementID,QUAD8); break;
            }
        }
    }
    cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____elements added____"<<endl;

    cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____adding nodes____"<<endl;
    int localNodeID = 0;
    for(it.Initialize(surfaceMesh->GetAllNodes());it.More(); it.Next())
    {
        localNodeID++;
        int globalNodeID = it.Key();
        myNodes.Add(globalNodeID);
        myNodesMap.Add(globalNodeID);

        int NbNodes;
        MeshVS_EntityType type;
        double buf[3];
        TColStd_Array1OfReal coords(*buf,1,3);
        surfaceMesh->GetGeom(globalNodeID,false,coords,NbNodes,type);
        myNodeCoords->SetValue(localNodeID,1,coords(1));
        myNodeCoords->SetValue(localNodeID,2,coords(2));
        myNodeCoords->SetValue(localNodeID,3,coords(3));
    }
    cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____nodes added____"<<endl;

    //! ---------------------------
    //! compute normal at elements
    //! ---------------------------
    this->computeNormalAtElements();
    cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____constructor done____"<<endl;
}

//! ----------------------------------------------------
//! function: computeNodeToElementsConnectivity
//! details:  uses nodes and elements "local" numbering
//! ----------------------------------------------------
void Ng_MeshVS_DataSourceFace::computeNodeToElementsConnectivity()
{
    myNodeToElements.clear();

    //! -------------------------------
    //! iterate over the mesh elements
    //! -------------------------------
    for(int localElementID=1; localElementID<=myNumberOfElements; localElementID++)
    {
        //! -------------------------------------
        //! iterate over the nodes of an element
        //! -------------------------------------
        int NbNodes;
        switch(myElemType->Value(localElementID))
        {
        case TRIG: NbNodes = 3; break;
        case QUAD: NbNodes = 4; break;
        case TRIG6: NbNodes = 6; break;
        case QUAD8: NbNodes = 8; break;
        }

        for(int k=1; k<=NbNodes; k++)
        {
            int globalNodeID = myElemNodes->Value(localElementID,k);
            int localNodeID = myNodesMap.FindIndex(globalNodeID);
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
}

//! -------------------------------------------
//! function: computeNodeToSegmentConnectivity
//! details:  key: localNodeID
//! -------------------------------------------
void Ng_MeshVS_DataSourceFace::computeNodeToSegmentConnectivity()
{
    cout<<"Ng_MeshVS_DataSourceFace::computeNodeToSegmentConnectivity()->____function called____"<<endl;
    for(int localElementID = 1; localElementID<=myNumberOfNodes; localElementID++)
    {
        int NbNodes = 0;
        switch(myElemType->Value(localElementID))
        {
        case TRIG: NbNodes = 3; break;
        case QUAD: NbNodes = 4; break;
        }

        for(int n = 0; n < NbNodes; n++)
        {
            int globalNodeID = myElemNodes->Value(localElementID,n+1);
            int localNodeID = myNodesMap.FindIndex(globalNodeID);
            //cout<<"Ng_MeshVS_DataSourceFace::computeNodeToSegmentConnectivity()->____("<<localNodeID<<", "<<globalNodeID<<")____"<<endl;

            //! --------------
            //! first segment
            //! --------------
            int first1 = n;
            int second1 = (n+1)%NbNodes;
            mesh::meshSegment f;
            //cout<<"____(first1, second1) = ("<<first1<<", "<<second1<<")____"<<endl;

            f.nodeIDs<<myElemNodes->Value(localElementID,first1+1)<<myElemNodes->Value(localElementID,second1+1);

            //! ---------------
            //! second segment
            //! ---------------
            int first2 = (n+2)%NbNodes;
            int second2 = n;
            mesh::meshSegment s;
            //cout<<"____(first2, second2) = ("<<first2<<", "<<second2<<")____"<<endl;

            s.nodeIDs<<myElemNodes->Value(localElementID,first2+1)<<myElemNodes->Value(localElementID,second2+1);

            //QMap<int,std::set<mesh::meshSegment>> myNodeToSegmentConnectivity;
            if(myNodeToSegmentConnectivity.contains(localNodeID))
            {
                std::set<mesh::meshSegment> listOfSegment = myNodeToSegmentConnectivity.value(localNodeID);
                listOfSegment.insert(f);
                listOfSegment.insert(s);
                myNodeToSegmentConnectivity.insert(localNodeID,listOfSegment);
            }
            else
            {
                std::set<mesh::meshSegment> listOfSegment;
                listOfSegment.insert(f);
                listOfSegment.insert(s);
                myNodeToSegmentConnectivity.insert(localNodeID,listOfSegment);
            }
        }
    }
}

//! -----------------------------------------------------------------
//! class gp_PntExtended: operator < added using hash of coordinates
//! this class is used by the constructor which follows
//! -----------------------------------------------------------------
class gp_PntExtended: public gp_Pnt
{
public:

    bool operator < (const gp_PntExtended &rhs) const
    {
        size_t seed, seed1;
        seed = 0; seed1 = 0;
        hash_c<double>(seed,this->X());
        hash_c<double>(seed,this->Y());
        hash_c<double>(seed,this->Z());
        hash_c<double>(seed1,rhs.X());
        hash_c<double>(seed1,rhs.Y());
        hash_c<double>(seed1,rhs.Z());
        if(seed<seed1) return true;
        return false;
    }
};

//! --------------------------------------------------
//! function: constructor
//! details:  input TopoDS_Shape, face number
//! warining: local and global numbering are the same
//!           both for nodes and both for elements
//! --------------------------------------------------
Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace(const TopoDS_Shape &shape, int faceNr)
{
    TopTools_IndexedMapOfShape faceMap;
    TopExp::MapShapes(shape,TopAbs_FACE,faceMap);
    const TopoDS_Face &face = TopoDS::Face(faceMap.FindKey(faceNr));

    //! sanity check
    if(face.IsNull()) return;

    TriangleAccessor aTool (face);

    myNumberOfElements = aTool.NbTriangles();

    //! sanity check
    if(myNumberOfElements<1) return;

    myElemNodes = new TColStd_HArray2OfInteger(1,myNumberOfElements,1,3);
    myElemType = new TColStd_HArray1OfInteger(1,myNumberOfElements);
    myElemNormals = new TColStd_HArray2OfReal(1,myNumberOfElements,1,3);

    std::set<gp_PntExtended> points;

    for (int iTri = 1; iTri <= aTool.NbTriangles(); iTri++)
    {
        gp_Vec aNorm;
        gp_PntExtended aPnt1, aPnt2, aPnt3;

        aTool.GetTriangle (iTri, aNorm, aPnt1, aPnt2, aPnt3);

        //! ---------------------
        //! set the element type
        //! ---------------------
        myElemType->SetValue(iTri,TRIG);

        //! -----------------------
        //! set the element normal
        //! -----------------------
        myElemNormals->SetValue(iTri,1,aNorm.X());
        myElemNormals->SetValue(iTri,2,aNorm.Y());
        myElemNormals->SetValue(iTri,3,aNorm.Z());

        //! ----------------
        //! map of elements
        //! ----------------
        myElements.Add(iTri);
        myElementsMap.Add(iTri);

        //! -------
        //! points
        //! -------
        points.insert(aPnt1);
        points.insert(aPnt2);
        points.insert(aPnt3);
    }

    myNumberOfNodes = int(points.size());

    //! sanity check
    if(myNumberOfNodes<3) return;

    myNodeCoords = new TColStd_HArray2OfReal(1,myNumberOfNodes,1,3);

    int nodeID = 0; //! dummy initialization
    for(std::set<gp_PntExtended>::iterator it = points.begin(); it!=points.end(); it++)
    {
        nodeID++;
        const gp_PntExtended &P = *it;
        myNodeCoords->SetValue(nodeID,1,P.X());
        myNodeCoords->SetValue(nodeID,2,P.Y());
        myNodeCoords->SetValue(nodeID,3,P.Z());

        myNodes.Add(nodeID);
        myNodesMap.Add(nodeID);
    }
}

//! ------------------------
//! function: renumberNodes
//! details:
//! ------------------------
void Ng_MeshVS_DataSourceFace::renumberNodes(const std::map<size_t,int> &mapCoordNodeNumbers,
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
                cerr<<"// Ng_MeshVS_DataSourceFace::renumberNodes()                     //"<<endl;
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
            cerr<<"// Ng_MeshVS_DataSourceFace::renumberNodes(): this part of code //"<<endl;
            cerr<<"// shoud never be reached.                                      //"<<endl;
            cerr<<"//--------------------------------------------------------------//"<<endl;
            exit(102);
        }
    }
}

//! --------------------------
//! function: setNodeNormal
//! details:  access function
//! --------------------------
void Ng_MeshVS_DataSourceFace::setNodeNormal(int globalNodeID, double *n)
{
    QList<double> normal; normal<<n[0]<<n[1]<<n[2];
    myNodeNormals.insert(globalNodeID,normal);
}


//! ------------------------------------------------------------------
//! function: constructor
//! details:  build the data source from a list of (2D) mesh elements
//! ------------------------------------------------------------------
Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace(const std::vector<meshElementByCoords> &meshElements, bool autoRenumberElements, bool autoRenumberNodes)
{
    //cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____constructor called____"<<endl;

    //! -------------
    //! sanity check
    //! -------------
    if(meshElements.size()==0) return;

    myNumberOfElements = int(meshElements.size());

    myElemType = new TColStd_HArray1OfInteger(1,myNumberOfElements);
    myElemNodes = new TColStd_HArray2OfInteger(1,myNumberOfElements,1,8);
    myElemNormals = new TColStd_HArray2OfReal(1,myNumberOfElements,1,3);

    std::vector<mesh::meshPoint> vecMeshPoints;
    std::vector<int> vecNodeIDs;

    int localElementID = 0;
    int localNodeID = 0;
    for(std::vector<meshElementByCoords>::const_iterator it = meshElements.cbegin(); it!= meshElements.cend(); it++)
    {
        const meshElementByCoords &aMeshElement = *it;

        localElementID++;

        if(autoRenumberElements==true)
        {
            myElements.Add(localElementID);
            myElementsMap.Add(localElementID);
        }
        else
        {
            int globalElementID = aMeshElement.ID;
            myElements.Add(globalElementID);
            myElementsMap.Add(globalElementID);
        }

        int NbPoints = aMeshElement.pointList.length();

        switch(NbPoints)
        {
        case 3: myElemType->SetValue(localElementID,TRIG); break;
        case 6: myElemType->SetValue(localElementID,TRIG6); break;
        case 4: myElemType->SetValue(localElementID,QUAD); break;
        case 8: myElemType->SetValue(localElementID,QUAD8); break;
        }

        //! ---------------------------------------------------------------
        //! this is the mechanism for the automatic numbering of the nodes
        //! ---------------------------------------------------------------
        if(autoRenumberNodes)
        {
            for(int i = 0; i<NbPoints; i++)
            {
                mesh::meshPoint aPoint = aMeshElement.pointList[i];
                std::vector<mesh::meshPoint>::iterator nit = std::find(vecMeshPoints.begin(),vecMeshPoints.end(),aPoint);

                if(nit==vecMeshPoints.end())
                {
                    //! ----------------
                    //! point not found
                    //! ---------------
                    localNodeID++;
                    myNodes.Add(localNodeID);
                    myNodesMap.Add(localNodeID);

                    vecMeshPoints.push_back(aPoint);
                    vecNodeIDs.push_back(localNodeID);

                    myElemNodes->SetValue(localElementID,i+1,localNodeID);
                }
                else
                {
                    //! ----------------------
                    //! point already present
                    //! ----------------------
                    int indexOf = std::distance(vecMeshPoints.begin(),nit);
                    int foundLocalNodeID = vecNodeIDs[indexOf];
                    myElemNodes->SetValue(localElementID,i+1,foundLocalNodeID);
                }
            }
        }
        else
        {
            //! ---------------------------------------
            //! label the nodes using their own labels
            //! ---------------------------------------
            for(int i=0; i<NbPoints; i++)
            {
                mesh::meshPoint aPoint = aMeshElement.pointList[i];
                std::vector<mesh::meshPoint>::iterator nit = std::find(vecMeshPoints.begin(),vecMeshPoints.end(),aPoint);
                if(nit==vecMeshPoints.end())
                {
                    vecMeshPoints.push_back(aPoint);
                }

                if(!myNodes.Contains(aPoint.ID))
                {
                    myNodes.Add(aPoint.ID);
                    myNodesMap.Add(aPoint.ID);
                }
                myElemNodes->SetValue(localElementID,i+1,aPoint.ID);
            }
        }
    }

    myNumberOfNodes = myNodes.Extent();
    myNodeCoords = new TColStd_HArray2OfReal(1,myNumberOfNodes,1,3);

    for(int i=1; i<=myNumberOfNodes; i++)
    {
        mesh::meshPoint aMeshPoint = vecMeshPoints[i-1];
        myNodeCoords->SetValue(i,1,aMeshPoint.x);
        myNodeCoords->SetValue(i,2,aMeshPoint.y);
        myNodeCoords->SetValue(i,3,aMeshPoint.z);
    }

    //! ---------------------------
    //! compute normal at elements
    //! ---------------------------
    this->computeNormalAtElements();
}


//! ----------------------------------
//! function: getSurroundingNodes
//! details:  uses "global" numbering
//! ----------------------------------
bool Ng_MeshVS_DataSourceFace::getSurroundingNodes(int globalNodeID, std::vector<int> &surroundingNodes_global)
{
    //! -------------
    //! sanity check
    //! -------------
    int localNodeID = myNodesMap.FindKey(globalNodeID);
    if(localNodeID<1 || localNodeID>myNumberOfNodes)
    {
        cerr<<"Ng_MeshVS_DataSourceFace::getSurroundingPoints()->____cannot retrieve nodeID: "<<globalNodeID<<". Out of range____"<<endl;
        return false;
    }

    //! ---------------------------------------------------------------
    //! note: the node to element connectivity uses "local" numbering
    //! ---------------------------------------------------------------
    if(myNodeToElements.isEmpty()) this->computeNodeToElementsConnectivity();

    QList<int> attachedElements_localIDs = myNodeToElements.value(localNodeID);
    if(attachedElements_localIDs.isEmpty())
    {
        cerr<<"Ng_MeshVS_DataSourceFace::getSurroundingPoints()->____no attached element for global node ID: "<<globalNodeID<<"____"<<endl;
        return false;
    }
    int NbAttachedElements = attachedElements_localIDs.length();
    for(int n=0; n<NbAttachedElements; n++)
    {
        int localElementID = attachedElements_localIDs[n];
        int globalElementID = myElementsMap.FindKey(localElementID);
        int NbNodes, buf[20];
        TColStd_Array1OfInteger nodeIDs(*buf,1,20);
        this->GetNodesByElement(globalElementID,nodeIDs,NbNodes);
        for(int n=1; n<=NbNodes; n++)
        {
            int globalNodeID_attached = nodeIDs(n);
            if(globalNodeID_attached == globalNodeID) continue; //! jump over the input node
            if(std::find(surroundingNodes_global.begin(), surroundingNodes_global.end(), globalNodeID_attached)!=surroundingNodes_global.end()) continue;
            surroundingNodes_global.push_back(globalNodeID_attached);
        }
    }
    return true;
}


//! ------------------------------------------------------------
//! function: displaceMySelf_asRigidAsPossible
//! details:  the "displacementField" is defined on the handles
//! ------------------------------------------------------------
void Ng_MeshVS_DataSourceFace::displaceMySelf_asRigidAsPossible(const QMap<int, gp_Vec> &displacementField)
{
    cout<<"Ng_MeshVS_DataSourceFace::displaceMySelf_asRigidAsPossible()->____function called____"<<endl;

    //! ---------------------------------------------------------------
    //! convert into Eigen format the mesh: V position of the vertices
    //! ---------------------------------------------------------------
    Eigen::MatrixXd V(myNumberOfNodes,3);
    for(int row=1; row<=myNumberOfNodes; row++)
        for(int col=1; col<=3; col++)
            V(row-1,col-1) = myNodeCoords->Value(row,col);

    //! ---------------------------------------------------------------
    //! convert into Eigen format the mesh: definition of the elements
    //! ---------------------------------------------------------------
    Eigen::MatrixXi T(myNumberOfElements,4);
    for(int i=1; i<=myNumberOfElements; i++)
    {
        int index1 = myElemNodes->Value(i,1);
        int index2 = myElemNodes->Value(i,2);
        int index3 = myElemNodes->Value(i,3);
        int index4 = myElemNodes->Value(i,4);

        int row = i-1;
        T(row,0) = index1-1;
        T(row,1) = index2-1;
        T(row,2) = index3-1;
        T(row,3) = index4-1;
    }

    for(QMap<int,gp_Vec>::const_iterator it = displacementField.cbegin(); it !=displacementField.cend(); it++)
    {
        ;
    }
}

//! --------------------
//! function: writeMesh
//! details:
//! --------------------
bool Ng_MeshVS_DataSourceFace::writeMesh(const QString &meshFileName, int dim)
{
    //cout<<"Ng_MeshVS_DataSourceFace::writeMesh()->____function called. Writing mesh file: "<<meshFileName.toStdString()<<"____"<<endl;

    char name[128];
    sprintf(name,"%s",meshFileName.toStdString().c_str());
    FILE *meshFile = fopen(name,"w");

    if(meshFile==NULL)
    {
        cout<<"Ng_MeshVS_DataSourceFace::writeMesh()->____cannot open file for writing mesh____"<<endl;
        return false;
    }
    //! --------------------
    //! write the dimension
    //! --------------------
    fprintf(meshFile,"%d\n",dim);
    fprintf(meshFile,"%d\n",myNodes.Extent());

    //! ----------------
    //! write the nodes
    //! ----------------
    double aCoordsBuf[3];
    TColStd_Array1OfReal aCoords(*aCoordsBuf,1,3);
    int nbNodes;
    MeshVS_EntityType aType;
    TColStd_MapIteratorOfPackedMapOfInteger anIter;

    for(anIter.Initialize(myNodes);anIter.More();anIter.Next())
    {
        int aKey = anIter.Key();
        if(this->GetGeom(aKey,false,aCoords,nbNodes,aType)==true)
        {
            //! write the nodeID (global) and the coordinates
            fprintf(meshFile,"%d\t%.9e\t%.9e\t%.9e\n",aKey,aCoords.Value(1),aCoords.Value(2),aCoords.Value(3));
        }
    }

    //! -------------------
    //! write the elements
    //! -------------------
    fprintf(meshFile,"%d\n",myElements.Extent());

    for(anIter.Initialize(myElements);anIter.More();anIter.Next())
    {
        int globalElementID = anIter.Key();
        int NbNodes;
        int aNodesBuf[20];
        TColStd_Array1OfInteger nodeIDs(*aNodesBuf,1,20);
        this->GetNodesByElement(globalElementID,nodeIDs,NbNodes);

        int localElementID = myElementsMap.FindIndex(globalElementID);
        ElemType eType = static_cast<ElemType>(myElemType->Value(localElementID));

        char header[16];
        switch(eType)
        {
        case TRIG: sprintf(header,"TRIG"); break;
        case QUAD: sprintf(header,"QUAD"); break;
        case TRIG6: sprintf(header,"TRIG6"); break;
        case QUAD8: sprintf(header,"QUAD8"); break;
        }

        //! ----------------------------------------
        //! write the element key, the element type
        //! ----------------------------------------
        fprintf(meshFile,"%d\n",globalElementID);
        fprintf(meshFile,"%s\n",header);

        //! ------------------
        //! write the nodeIDs
        //! ------------------
        for(int i=1; i<NbNodes; i++) fprintf(meshFile,"%d\t",nodeIDs(i));
        fprintf(meshFile,"%d\n",nodeIDs(NbNodes));
    }
    //! ---------------
    //! close the file
    //! ---------------
    fclose(meshFile);
    //cout<<"Ng_MeshVS_DataSourceFace::writeMesh()->____exiting____"<<endl;
    return true;
}

//! ---------------------------
//! function: write
//! details:  use a C++ stream
//! ---------------------------
bool Ng_MeshVS_DataSourceFace::write(ofstream &os, int dim)
{
    if(!os.is_open()) return false;

    cout<<"Ng_MeshVS_DataSourceFace::write()->____function called____"<<endl;

    //! ----------------
    //! nodes, elements
    //! ----------------
    const TColStd_PackedMapOfInteger &aNodes = this->GetAllNodes();
    const TColStd_PackedMapOfInteger &mapOfElements = this->GetAllElements();

    if(aNodes.Extent()==0)
    {
        cerr<<"Ng_MeshVS_DataSourceFace::write()->____no node to write____"<<endl;
        return false;
    }
    if(mapOfElements.Extent()==0)
    {
        cerr<<"Ng_MeshVS_DataSourceFace::write()->____no element to write____"<<endl;
        return false;
    }

    //! ---------------
    //! mesh dimension
    //! ---------------
    os<<dim<<endl;
    os<<aNodes.Extent()<<endl;

    double buf[3];
    TColStd_Array1OfReal aCoords(*buf,1,3);
    int nbNodes;
    MeshVS_EntityType aType;
    TColStd_MapIteratorOfPackedMapOfInteger anIter;

    //! -------------------------------------------
    //! write the globalNodeID and the coordinates
    //! -------------------------------------------
    for(anIter.Initialize(aNodes);anIter.More();anIter.Next())
    {
        int globalNodeID = anIter.Key();
        if(this->GetGeom(globalNodeID,false,aCoords,nbNodes,aType)==true)
        {
            std::streamsize ss = os.precision();
            os.precision(9);
            os<<globalNodeID<<"\t"<<aCoords(1)<<"\t"<<aCoords(2)<<"\t"<<aCoords(3)<<endl;
            os.precision(ss);
        }
    }

    //! -----------------------------
    //! write the number of elements
    //! -----------------------------
    int NbElements = mapOfElements.Extent();
    os<<NbElements<<endl;

    //! ------------------------------
    //! write the elements definition
    //! ------------------------------
    for(anIter.Initialize(mapOfElements);anIter.More();anIter.Next())
    {
        int aKey = anIter.Key();
        int NbNodes;
        int aNodesBuf[20];
        TColStd_Array1OfInteger nodeIDs(*aNodesBuf,1,20);
        this->GetNodesByElement(aKey,nodeIDs,NbNodes);

        MeshVS_EntityType aType;
        this->GetGeomType(aKey,true,aType);

        char header[16];
        if(aType==MeshVS_ET_Face)
        {
            switch(NbNodes)
            {
            case 3: sprintf(header,"TRIG"); break;
            case 4: sprintf(header,"QUAD"); break;
            case 6: sprintf(header,"TRIG6"); break;
            case 8: sprintf(header,"QUAD8"); break;
            }
        }
        else if(aType==MeshVS_ET_Volume)
        {
            switch(NbNodes)
            {
            case 4: sprintf(header,"TET"); break;
            case 8: sprintf(header,"HEXA"); break;
            case 10: sprintf(header,"TET10"); break;
            case 20: sprintf(header,"HEXA20"); break;
            case 6: sprintf(header,"PRISM"); break;
            case 5: sprintf(header,"PYRAM"); break;
            }
        }

        //! ----------------------------------------
        //! write the element key, the element type
        //! ----------------------------------------
        os<<aKey<<endl;
        os<<header<<endl;

        for(int i=1; i<NbNodes; i++) os<<nodeIDs(i)<<"\t";
        os<<nodeIDs(NbNodes)<<endl;
    }
    return true;
}

//! ----------------------------
//! function: constructor
//! details:  from a C++ stream
//! ----------------------------
Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace(ifstream &stream)
{
    //cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____function called____"<<endl;

    //! -------------------
    //! read the dimension
    //! -------------------
    int meshDimension;
    stream>>meshDimension;

    //! -------------------------
    //! read the number of nodes
    //! -------------------------
    stream>>myNumberOfNodes;

    if(meshDimension ==2 && myNumberOfNodes<3) return;
    if(meshDimension ==3 && myNumberOfNodes<4) return;

    myNodeCoords = new TColStd_HArray2OfReal(1,myNumberOfNodes,1,3);

    //! ---------------
    //! load the nodes
    //! ---------------
    for(int n=1; n<=myNumberOfNodes; n++)
    {
        int globalNodeID;
        double x,y,z;
        stream>>globalNodeID>>x>>y>>z;

        myNodesMap.Add(globalNodeID);
        myNodes.Add(globalNodeID);
        myNodeCoords->SetValue(n,1,x);
        myNodeCoords->SetValue(n,2,y);
        myNodeCoords->SetValue(n,3,z);
    }

    //! ----------------------------
    //! read the number of elements
    //! ----------------------------
    stream>>myNumberOfElements;

    //cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____number of elements: "<<myNumberOfElements<<"____"<<endl;

    if(myNumberOfElements<1) return;
    myElemNodes = new TColStd_HArray2OfInteger(1,myNumberOfElements,1,20);
    myElemType = new TColStd_HArray1OfInteger(1,myNumberOfElements);
    myElemNormals = new TColStd_HArray2OfReal(1,myNumberOfElements,1,3);

    //! ------------------
    //! load the elements
    //! ------------------
    for(int i=1;i<=myNumberOfElements; i++)
    {
        int globalElementID;
        stream>>globalElementID;

        std::string val;
        stream>>val;
        char atypes[256];   // element type cstring
        sscanf(val.c_str(),"%s",atypes);
        if(strcmp(atypes,"TRIG")==0)
        {
            int i1,i2,i3;
            stream>>i1>>i2>>i3;
            myElemNodes->SetValue(i,1,i1);
            myElemNodes->SetValue(i,2,i2);
            myElemNodes->SetValue(i,3,i3);
            myElemType->SetValue(i,TRIG);
        }
        if(strcmp(atypes,"QUAD")==0)
        {
            int i1,i2,i3,i4;
            stream>>i1>>i2>>i3>>i4;
            myElemNodes->SetValue(i,1,i1);
            myElemNodes->SetValue(i,2,i2);
            myElemNodes->SetValue(i,3,i3);
            myElemNodes->SetValue(i,4,i4);
            myElemType->SetValue(i,QUAD);
        }
        if(strcmp(atypes,"TET")==0)
        {
            int i1,i2,i3,i4;
            stream>>i1>>i2>>i3>>i4;
            myElemNodes->SetValue(i,1,i1);
            myElemNodes->SetValue(i,2,i2);
            myElemNodes->SetValue(i,3,i3);
            myElemNodes->SetValue(i,4,i4);
            myElemType->SetValue(i,TET);
        }
        if(strcmp(atypes,"HEXA")==0)
        {
            int i1,i2,i3,i4,i5,i6,i7,i8;
            stream>>i1>>i2>>i3>>i4>>i5>>i6>>i7>>i8;
            myElemNodes->SetValue(i,1,i1);
            myElemNodes->SetValue(i,2,i2);
            myElemNodes->SetValue(i,3,i3);
            myElemNodes->SetValue(i,4,i4);
            myElemNodes->SetValue(i,5,i5);
            myElemNodes->SetValue(i,6,i6);
            myElemNodes->SetValue(i,7,i7);
            myElemNodes->SetValue(i,8,i8);
            myElemType->SetValue(i,HEXA);
        }
        if(strcmp(atypes,"PRISM")==0)
        {
            int i1,i2,i3,i4,i5,i6;
            stream>>i1>>i2>>i3>>i4>>i5>>i6;
            myElemNodes->SetValue(i,1,i1);
            myElemNodes->SetValue(i,2,i2);
            myElemNodes->SetValue(i,3,i3);
            myElemNodes->SetValue(i,4,i4);
            myElemNodes->SetValue(i,5,i5);
            myElemNodes->SetValue(i,6,i6);
            myElemType->SetValue(i,PRISM);
        }
        if(strcmp(atypes,"PYRAM")==0)
        {
            int i1,i2,i3,i4,i5;
            stream>>i1>>i2>>i3>>i4>>i5;
            myElemNodes->SetValue(i,1,i1);
            myElemNodes->SetValue(i,2,i2);
            myElemNodes->SetValue(i,3,i3);
            myElemNodes->SetValue(i,4,i4);
            myElemNodes->SetValue(i,5,i5);
            myElemType->SetValue(i,PYRAM);
        }
        //! note: second order elements not supported. To do ...

        myElementsMap.Add(globalElementID);
        myElements.Add(globalElementID);
    }

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
Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace(const opencascade::handle<Ng_MeshVS_DataSource3D> &aVolumeMeshDS)
{
    //!cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSource2D()->____constructor from a volume mesh data source called____"<<endl;
    if(aVolumeMeshDS->GetAllNodes().Extent()<3) return;
    if(aVolumeMeshDS->GetAllElements().Extent()<1) return;

    const std::vector<meshElement2D> &surfaceElements = aVolumeMeshDS->mySurfaceElements;
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
        aVolumeMeshDS->GetGeom(globalNodeID,false,coords,NbNodes,type);

        myNodeCoords->SetValue(localNodeID,1,coords(1));
        myNodeCoords->SetValue(localNodeID,2,coords(2));
        myNodeCoords->SetValue(localNodeID,3,coords(3));
    }

    //! --------------------------------
    //! compute the normals at elements
    //! --------------------------------
    this->computeNormalAtElements();
}

//! -------------------------------
//! function: getSegmentsOfElement
//! details:
//! -------------------------------
std::vector<mesh::meshSegment> Ng_MeshVS_DataSourceFace::getSegmentsOfElement(int elementID, bool isLocal)
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
Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace(const occHandle(Ng_MeshVS_DataSourceFace) &aMesh, int order)
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
    for(TColStd_MapIteratorOfPackedMapOfInteger eit(aMesh->GetAllElements()); eit.More(); eit.Next())
    {
        int globalElementID = eit.Key();
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
    myNumberOfNodes = aMesh->GetAllNodes().Extent() + additionalNodes;
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


//! -------------------------------------
//! function: constructor
//! details:  meshA - meshB (meshA>meshB)
//! -------------------------------------
Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace(const occHandle(Ng_MeshVS_DataSourceFace) &meshDSA, occHandle(Ng_MeshVS_DataSourceFace) &meshDSB)
{
    cout<<"Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace()->____constructor called____"<<endl;

    TColStd_PackedMapOfInteger elementsMapA = meshDSA->GetAllElements();
    TColStd_PackedMapOfInteger elementsMapB = meshDSB->GetAllElements();
    myElements.Subtraction(elementsMapA,elementsMapB);
    myNumberOfElements = myElements.Extent();

    cout<<"____number of elements of mesh A: "<<elementsMapA.Extent()<<endl;
    cout<<"____number of elements of mesh B: "<<elementsMapB.Extent()<<endl;
    cout<<"____number of elements: "<<myNumberOfElements<<"____"<<endl;

    myElemType = new TColStd_HArray1OfInteger(1,myNumberOfElements);
    myElemNodes = new TColStd_HArray2OfInteger(1,myNumberOfElements,1,8);
    myElemNormals = new TColStd_HArray2OfReal(1,myNumberOfElements,1,3);

    int localElementID = 0;
    for(TColStd_MapIteratorOfPackedMapOfInteger it(myElements); it.More(); it.Next())
    {
        localElementID++;
        int globalElementID = it.Key();
        myElementsMap.Add(globalElementID);
        int NbNodes, buf[8];
        TColStd_Array1OfInteger nodeIDs(*buf,1,8);
        double nx,ny,nz;
        if(elementsMapA.Contains(globalElementID))
        {
            meshDSA->GetNodesByElement(globalElementID,nodeIDs,NbNodes);
            meshDSA->GetNormal(globalElementID,3,nx,ny,nz);
        }
        else
        {
            meshDSB->GetNodesByElement(globalElementID,nodeIDs,NbNodes);
            meshDSB->GetNormal(globalElementID,3,nx,ny,nz);
        }
        for(int n=1; n<=NbNodes; n++)
        {
            int globalNodeID = nodeIDs(n);
            myElemNodes->SetValue(localElementID,n,globalNodeID);
            if(!myNodes.Contains(globalNodeID))
            {
                myNodes.Add(globalNodeID);
                myNodesMap.Add(globalNodeID);
            }
        }
        switch(NbNodes)
        {
        case 3: myElemType->SetValue(localElementID,TRIG); break;
        case 4: myElemType->SetValue(localElementID,QUAD); break;
        case 6: myElemType->SetValue(localElementID,TRIG6); break;
        case 8: myElemType->SetValue(localElementID,QUAD8); break;
        }
        myElemNormals->SetValue(localElementID,1,nx);
        myElemNormals->SetValue(localElementID,2,ny);
        myElemNormals->SetValue(localElementID,3,nz);
    }

    cout<<"____definition of elements OK____"<<endl;

    myNumberOfNodes = myNodes.Extent();
    myNodeCoords = new TColStd_HArray2OfReal(1,myNumberOfNodes,1,3);
    cout<<"____number of nodes of mesh A: "<<meshDSA->GetAllNodes().Extent()<<endl;
    cout<<"____number of nodes of mesh B: "<<meshDSB->GetAllNodes().Extent()<<endl;
    cout<<"____number of nodes: "<<myNumberOfNodes<<"____"<<endl;

    //! ------------------
    //! nodes coordinates
    //! ------------------
    for(int localNodeID = 1; localNodeID <=myNumberOfNodes; localNodeID++)
    {
        int globalNodeID = myNodesMap.FindKey(localNodeID);
        int localNodeID_forSearch;
        std::vector<double> c;
        if(meshDSA->myNodesMap.Contains(globalNodeID))
        {
           localNodeID_forSearch = meshDSA->myNodesMap.FindIndex(globalNodeID);
           c = meshDSA->getNodeCoordinates(localNodeID_forSearch);
        }
        else
        {
            localNodeID_forSearch = meshDSB->myNodesMap.FindIndex(globalNodeID);
            c = meshDSB->getNodeCoordinates(localNodeID_forSearch);
        }
        myNodeCoords->SetValue(localNodeID,1,c[0]);
        myNodeCoords->SetValue(localNodeID,2,c[1]);
        myNodeCoords->SetValue(localNodeID,3,c[2]);
    }

    cout<<"____node coordinates OK____"<<endl;

    //! ---------------------------
    //! compute normal at elements
    //! ---------------------------
    //this->computeNormalAtElements();
}
