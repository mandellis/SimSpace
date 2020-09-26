//! ----------------
//! custom includes
//! ----------------
#include "ng_meshvs_datasource3d.h"
#include "mydefines.h"
#include "ccout.h"
#include <polygon.h>

//! ----
//! OCC
//! ----
#include<MeshVS_HArray1OfSequenceOfInteger.hxx>
#include <MeshVS_EntityType.hxx>
#include <MeshVS_Buffer.hxx>
#include <MeshVS_Tool.hxx>
#include <TColStd_HArray1OfReal.hxx>
#include <TColStd_MapIteratorOfPackedMapOfInteger.hxx>

//! ---
//! Qt
//! ---
#include <QMetaType>
#include <QTime>

//! ----
//! C++
//! ----
#include <fstream>
#include <memory>
#include <set>

//! ----
//! igl
//! ----
#include <libigl/include/igl/ARAPEnergyType.h>
#include <libigl/include/igl/arap_linear_block.h>
#include <libigl/include/igl/arap.h>
#include <libigl/include/igl/harmonic.h>


IMPLEMENT_STANDARD_HANDLE(Ng_MeshVS_DataSource3D,MeshVS_DataSource)
IMPLEMENT_STANDARD_RTTIEXT(Ng_MeshVS_DataSource3D,MeshVS_DataSource)

const int NMAX = 20;
const int MTRIG3[3]={1,2,3};                    //! OK
const int MQUAD4[4]={1,2,3,4};                  //! OK
const int MTRIG6[6]={1,6,2,4,3,5};              //! OK
const int MQUAD8[8]={1,5,2,8,3,6,4,7};          //! OK
const int MTET4[4]={1,2,3,4};                   //! OK
const int MTET10[10]={1,2,3,4,5,6,7,8,9,10};    //! OK

//! --------------------------
//! function: clock
//! details:  time diagnostic
//! --------------------------
QTime Ng_MeshVS_DataSource3D::clock()
{
    QTime time = QTime::currentTime();
    QString text = time.toString("hh::mm:ss:zzz");
    cout<<text.toStdString()<<endl;
    return time;
}

//! -------------------------------------
//! function: constructor
//! details:  from a list of 3D elements
//! -------------------------------------
Ng_MeshVS_DataSource3D::Ng_MeshVS_DataSource3D(const std::vector<meshElementByCoords> &meshElements,
                                               bool autoNumberElements,
                                               bool autoNumberNodes)
{
    cout<<"Ng_MeshVS_DataSource3D::Ng_MeshVS_DataSource3D()->____constructor from a QList meshElementByCoords called____"<<endl;
    myNumberOfElements = (int)meshElements.size();
    if(myNumberOfElements<1)
    {
        cerr<<"Ng_MeshVS_DataSource3D::Ng_MeshVS_DataSource3D()->____the mesh has no element____"<<endl;
        return;
    }

    myElemType = new TColStd_HArray1OfInteger(1,myNumberOfElements);
    myElemNodes = new TColStd_HArray2OfInteger(1,myNumberOfElements,1,20);

    std::map<mesh::meshPoint,int> meshPointMap;

    int n = 0;
    for(int i=1; i<=myNumberOfElements; i++)
    {
        const meshElementByCoords &curMeshElement = meshElements.at(i-1);

        //! -----------------------
        //! fill the elements maps
        //! -----------------------
        int globalElementID = autoNumberElements==true? i:curMeshElement.ID;
        myElements.Add(globalElementID);
        myElementsMap.Add(globalElementID);

        int NbNodes = curMeshElement.pointList.length();

        for(int j=0; j<NbNodes; j++)
        {
            const mesh::meshPoint &curMeshPoint = curMeshElement.pointList.at(j);
            if(meshPointMap.count(curMeshPoint)==false)
            {
                //! ----------------------------------------------------------
                //! if autoNumberNodes == true the ID of the node is the ID
                //! of read from the meshPoint, otherwise it is automatically
                //! incremented
                //! ----------------------------------------------------------
                if(autoNumberNodes)
                {
                    n++;
                    meshPointMap.insert(std::make_pair(curMeshPoint,n));
                }
                else
                {
                    int globalNodeID = curMeshPoint.ID;
                    meshPointMap.insert(std::make_pair(curMeshPoint,globalNodeID));
                }
            }
        }
        //! ---------------------
        //! set the element type
        //! ---------------------
        myElemType->SetValue(i,curMeshElement.type);
    }

    //! ----------------------------------------
    //! add the nodes and the nodes coordinates
    //! ----------------------------------------
    myNumberOfNodes = int(meshPointMap.size());
    if(myNumberOfNodes<3)
    {
        cerr<<"Ng_MeshVS_DataSource3D::Ng_MeshVS_DataSource3D()->____wrong number of nodes (<3)____"<<endl;
        return;
    }
    myNodeCoords = new TColStd_HArray2OfReal(1,myNumberOfNodes,1,3);

    int i=0;
    int localNodeID = 0;
    for(std::map<mesh::meshPoint,int>::iterator it = meshPointMap.begin(); it!=meshPointMap.end(); ++it)
    {
        localNodeID++;
        i = it->second;

        //! --------------------
        //! fill the nodes maps
        //! --------------------
        myNodes.Add(i);
        myNodesMap.Add(i);

        //! ------------------
        //! nodes coordinates
        //! ------------------
        const mesh::meshPoint &aP = it->first;
        myNodeCoords->SetValue(localNodeID,1,aP.x);
        myNodeCoords->SetValue(localNodeID,2,aP.y);
        myNodeCoords->SetValue(localNodeID,3,aP.z);
    }

    //cout<<"Ng_MeshVS_DataSource3D::Ng_MeshVS_DataSource3D()->____Number of points: "<<myNumberOfNodes<<" number of elements: "<<myNumberOfElements<<"____"<<endl;

    //! ---------------------------------------------------
    //! build the definition of the elements through nodes
    //! ---------------------------------------------------
    for(int i=1; i<=myNumberOfElements; i++)
    {
        const meshElementByCoords &curmeshElementByCoords = meshElements.at(i-1);
        int NbNodes = curmeshElementByCoords.pointList.length();
        for(int j=0; j<NbNodes; j++)
        {
            const mesh::meshPoint &curMeshPoint = curmeshElementByCoords.pointList.at(j);
            int indexOfPoint = meshPointMap.at(curMeshPoint);
            myElemNodes->SetValue(i,j+1,indexOfPoint);
        }
    }

    //! ------------------
    //! elements topology
    //! ------------------
    this->buildElementsTopology();

    cout<<"Ng_MeshVS_DataSource3D::Ng_MeshVS_DataSource3D()->____constructor from a QList of meshElementByCoords: OK____"<<endl;
}

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
Ng_MeshVS_DataSource3D::Ng_MeshVS_DataSource3D(const std::vector<mesh::meshElement> &elements, const std::vector<mesh::meshPoint> &nodes)
{
    cout<<"Ng_MeshVS_DataSource3D::Ng_MeshVS_DataSource3D()->____from vector of mesh points and mesh elements____"<<endl;
    myNumberOfElements= int(elements.size());
    myNumberOfNodes = int(nodes.size());

    if(myNumberOfElements<1) return;
    if(myNumberOfNodes<4) return;

    myElemNodes = new TColStd_HArray2OfInteger(1,myNumberOfElements,1,20);
    myElemType = new TColStd_HArray1OfInteger(1,myNumberOfElements);
    myNodeCoords = new TColStd_HArray2OfReal(1,myNumberOfNodes,1,3);

    //FILE *f0 = fopen("D:\\nodes_mesh.txt","w");
    //! -------------------------------------------
    //! fill the node maps and add the coordinates
    //! -------------------------------------------
    for (int i=1;i<=myNumberOfNodes;i++)
    {
        mesh::meshPoint theCurNode = nodes.at(i-1);
        int globalNodeID = theCurNode.ID;
        myNodes.Add(globalNodeID);
        myNodesMap.Add(globalNodeID);
        myNodeCoords->SetValue(i,1,theCurNode.x);
        myNodeCoords->SetValue(i,2,theCurNode.y);
        myNodeCoords->SetValue(i,3,theCurNode.z);
        //fprintf(f0,"%d\t%d\t(%lf\t%lf\t%lf)\n",i,globalNodeID,theCurNode.x,theCurNode.y,theCurNode.z);
    }
    //fclose(f0);

    //FILE *f1 = fopen("D:\\elements_mesh.txt","w");
    //! ----------------------------------------------
    //! fill the element maps and define the elements
    //! ----------------------------------------------
    for (int i=1;i<=myNumberOfElements;i++)
    {
        mesh::meshElement curElement = elements.at(i-1);

        //! -------------------------
        //! fill the map of elements
        //! -------------------------
        int globalElementID = curElement.ID;
        myElements.Add(globalElementID);
        myElementsMap.Add(globalElementID);

        //fprintf(f1,"%d\t%d\t",i,globalElementID);

        //! ----------------------------
        //! set the type of the element
        //! ----------------------------
        ElemType curType = curElement.type;
        myElemType->SetValue(i,curType);

        //! -----------------------------
        //! set the nodes of the element
        //! -----------------------------
        for(int k=0; k<curElement.theNodeIDs.size(); k++)
        {
            int nodeID = curElement.theNodeIDs.at(k);
            myElemNodes->SetValue(i,k+1,nodeID);

            //if(k==curElement.theNodeIDs.size()-1) fprintf(f1,"%d\n",nID);
            //else fprintf(f1,"%d\t",nID);
        }
    }
    //fclose(f1);

    //! ------------------
    //! elements topology
    //! ------------------
    this->buildElementsTopology();

    cout<<"Ng_MeshVS_DataSource3D::Ng_MeshVS_DataSource3D()->____from vector of mesh points and mesh elements: OK____"<<endl;
}

//! -------------------------------------
//! function: constructor
//! details:  build the mesh from a file
//! -------------------------------------
Ng_MeshVS_DataSource3D::Ng_MeshVS_DataSource3D(const std::string &fileName)
{
    cout<<"Ng_MeshVS_DataSource3D::Ng_MeshVS_DataSource3D()->____constructor from file: "<<fileName<<"____"<<endl;

    std::vector<mesh::meshElement> elements;
    std::vector<mesh::meshPoint> nodes;

    cout<<"Ng_MeshVS_DataSource3D::Ng_MeshVS_DataSource3D()->____reading mesh file: ->"<<fileName<<"<-____"<<endl;
    bool isDone = this->readMesh(QString::fromStdString(fileName),nodes,elements);
    if(!isDone)
    {
        cout<<"____error in reading the volume mesh____"<<endl;
        return;
    }

    myNumberOfElements= int(elements.size());
    myNumberOfNodes = int(nodes.size());

    cout<<"Ng_MeshVS_DataSource3D::Ng_MeshVS_DataSource3D()->____number of nodes: "<<myNumberOfNodes<<"____"<<endl;
    cout<<"Ng_MeshVS_DataSource3D::Ng_MeshVS_DataSource3D()->____number of elements: "<<myNumberOfElements<<"____\n"<<endl;

    if(myNumberOfElements<1) return;
    if(myNumberOfNodes<4) return;

    myElemNodes = new TColStd_HArray2OfInteger(1,myNumberOfElements,1,20);
    myElemType = new TColStd_HArray1OfInteger(1,myNumberOfElements);
    myNodeCoords = new TColStd_HArray2OfReal(1,myNumberOfNodes,1,3);

    //! -------------------------------------------
    //! fill the node maps and add the coordinates
    //! -------------------------------------------
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

    //! ----------------------------------------------
    //! fill the element maps and define the elements
    //! ----------------------------------------------
    for (int i=1;i<=myNumberOfElements;i++)
    {
        mesh::meshElement curElement = elements.at(i-1);

        //! -------------------------
        //! fill the map of elements
        //! -------------------------
        int globalElementID = curElement.ID;
        myElements.Add(globalElementID);
        myElementsMap.Add(globalElementID);

        //! ----------------------------
        //! set the type of the element
        //! ----------------------------
        ElemType curType = curElement.type;
        myElemType->SetValue(i,curType);

        //! -----------------------------
        //! set the nodes of the element
        //! -----------------------------
        for(int k=0; k<curElement.theNodeIDs.size(); k++)
        {
            int nID = curElement.theNodeIDs.at(k);
            myElemNodes->SetValue(i,k+1,nID);
        }
    }

    //! ------------------
    //! elements topology
    //! ------------------
    this->buildElementsTopology();

    //! -------------
    //! connectivity
    //! -------------
    this->buildFaceToElementConnectivity();
}

//! -----------------------------------------------
//! function: constructor
//! details:  build the mesh from a Netgen pointer
//! -----------------------------------------------
Ng_MeshVS_DataSource3D::Ng_MeshVS_DataSource3D(Ng_Mesh *aMesh)
{
    cout<<"Ng_MeshVS_DataSource3D::Ng_MeshVS_DataSource3D()->____function called: mesh from Netgen pointer____"<<endl;
    if(aMesh==NULL)
    {
        cout<<"Ng_MeshVS_DataSource3D::Ng_MeshVS_DataSource3D()->____Netgen pointer is NULL. Exiting____"<<endl;
        return;
    }
    if(Ng_GetNP(aMesh)<4)   //! at least 4 points
    {
        cout<<"Ng_MeshVS_DataSource3D::Ng_MeshVS_DataSource3D()->____Netgen volume mesh has no points. Exiting____"<<endl;
        return;
    }
    if(Ng_GetNE(aMesh)<1)   //! at least a tetrahedron
    {
        cout<<"Ng_MeshVS_DataSource3D::Ng_MeshVS_DataSource3D()->____Netgen volume mesh has no element. Exiting____"<<endl;
        return;
    }

    myNumberOfElements = Ng_GetNE(aMesh);
    myNumberOfNodes = Ng_GetNP(aMesh);

    myElemNodes = new TColStd_HArray2OfInteger(1,myNumberOfElements,1,10);
    myElemType = new TColStd_HArray1OfInteger(1,myNumberOfElements);
    myNodeCoords = new TColStd_HArray2OfReal(1,myNumberOfNodes,1,3);

    //! store the elements type
    int P[10];  //! not used - change nglib introducing "Ng_GetElementType(NgMesh*, int)"
    for (int i=1;i<=myNumberOfElements;i++)
    {
        myElements.Add(i);      //! global
        myElementsMap.Add(i);   //! global
        Ng_Volume_Element_Type type = Ng_GetVolumeElement(aMesh,i,P);
        switch(type)
        {
        case NG_TET: myElemType->SetValue(i,TET); break;
        case NG_TET10: myElemType->SetValue(i,TET10); break;
        }
    }

    //! node map and add the coordinates
    //! 3D coordinates of the node
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

    //! scan the elements
    int p_I[10];
    for (int i=1;i<=myNumberOfElements;i++)
    {
        Ng_Volume_Element_Type type = Ng_GetVolumeElement(aMesh,i,p_I);
        switch(type)
        {
        case NG_TET:
        {
            myElemNodes->SetValue(i,1,p_I[0]);
            myElemNodes->SetValue(i,2,p_I[1]);
            myElemNodes->SetValue(i,3,p_I[2]);
            myElemNodes->SetValue(i,4,p_I[3]);
        }
            break;

        case NG_TET10:
        {
            myElemNodes->SetValue(i,1,p_I[0]);
            myElemNodes->SetValue(i,2,p_I[1]);
            myElemNodes->SetValue(i,3,p_I[2]);
            myElemNodes->SetValue(i,4,p_I[3]);
            myElemNodes->SetValue(i,5,p_I[4]);
            myElemNodes->SetValue(i,6,p_I[5]);
            myElemNodes->SetValue(i,7,p_I[6]);
            myElemNodes->SetValue(i,8,p_I[7]);
            myElemNodes->SetValue(i,9,p_I[8]);
            myElemNodes->SetValue(i,10,p_I[9]);
        }
            break;
        }
    }

    //! ------------------
    //! elements topology
    //! ------------------
    this->buildElementsTopology();

    //! -----------------------------
    //! face to element connectivity
    //! -----------------------------
    this->buildFaceToElementConnectivity();
}

//! -------------------------------------------------------
//! function: constructor
//! details:  init with tetgenio
//!           only linear elements supported at the moment
//! -------------------------------------------------------
Ng_MeshVS_DataSource3D::Ng_MeshVS_DataSource3D(const tetgenio &aMesh)
{
    //!cout<<"Ng_MeshVS_DataSource::Ng_MeshVS_DataSource->____function called____"<<endl;
    myNumberOfElements = aMesh.numberoftetrahedra;
    myNumberOfNodes = aMesh.numberofpoints;

    myNodeCoords = new TColStd_HArray2OfReal(1,myNumberOfNodes,1,3);
    for(int i=1, S=0; i<=myNumberOfNodes; i++)
    {
        myNodes.Add(i);     //! global
        myNodesMap.Add(i);  //! global
        for(int k=0;k<3;k++)
        {
            double coord = aMesh.pointlist[S];
            myNodeCoords->SetValue(i,k+1,coord);
            S++;
        }
    }

    myElemType = new TColStd_HArray1OfInteger(1,myNumberOfElements);
    myElemNodes = new TColStd_HArray2OfInteger(1,myNumberOfElements,1,10);

    //! ---------------------------
    //! map of the volume elements
    //! nodes of the element
    //! ---------------------------
    const int C[4]{0,1,3,2};
    for(int S=0, i=1; i<=myNumberOfElements; i++)
    {
        myElements.Add(i);      //! local
        myElementsMap.Add(i);   //! global

        for(int k=0; k<4; k++)
        {
            int k1 = C[k];
            int nodeID = aMesh.tetrahedronlist[S];
            myElemNodes->SetValue(i,k1+1,nodeID);
            S++;
        }
        //! set the element type
        myElemType->SetValue(i,TET);
    }

    //! ------------------
    //! elements topology
    //! ------------------
    this->buildElementsTopology();
}

//! ----------------------------------------------
//! function: constructor
//! details:  constructor from a Tetgen .ele file
//! ----------------------------------------------
Ng_MeshVS_DataSource3D::Ng_MeshVS_DataSource3D(const QString &tetgenEleFileName, const QString &tetgenNodeFileName)
{
    cout<<"Ng_MeshVS_DataSource3D::Ng_MeshVS_DataSource3D()->____constructor from Tetgen files called____"<<endl;

    FILE *tetgenNodeFile = fopen(tetgenNodeFileName.toStdString().c_str(),"r");
    if(tetgenNodeFile==NULL) return;
    FILE *tetgenEleFile = fopen(tetgenEleFileName.toStdString().c_str(),"r");
    if(tetgenEleFile==NULL) return;

    //! ----------------------------------------------------------------------------------------
    //! read the nodes
    //! First line: <# of points> <dimension (3)> <# of attributes> <boundary markers (0 or 1)>
    //! ----------------------------------------------------------------------------------------
    cout<<"____Start reading: "<<tetgenNodeFileName.toStdString()<<"____"<<endl;
    unsigned int NbPoints, dimension, NbAttributes, boundaryMarkerFlag;
    fscanf(tetgenNodeFile,"%d%d%d%d",&NbPoints,&dimension,&NbAttributes,&boundaryMarkerFlag);

    myNumberOfNodes = NbPoints;
    myNodeCoords = new TColStd_HArray2OfReal(1,NbPoints,1,3);

    unsigned int nodeNumber, flag;
    double x,y,z;
    bool isNodeReading = true;
    for(unsigned int line = 1; line<=NbPoints; line++)
    {
        if(5==fscanf(tetgenNodeFile,"%d%lf%lf%lf%d",&nodeNumber,&x,&y,&z,&flag))
        {
            myNodes.Add(nodeNumber);
            myNodesMap.Add(nodeNumber);

            myNodeCoords->SetValue(nodeNumber,1,x);
            myNodeCoords->SetValue(nodeNumber,2,y);
            myNodeCoords->SetValue(nodeNumber,3,z);
        }
        else
        {
            isNodeReading = false;
            cout<<"____error in reading nodes____"<<endl;
            break;
        }
    }

    if(isNodeReading==true)
    {
        cout<<"____Nodes reading OK____"<<endl;
    }

    if(isNodeReading==true)
    {
        cout<<"____Start reading: "<<tetgenEleFileName.toStdString()<<"____"<<endl;
        //!ccout("____Start reading: "+tetgenEleFileName+"____");

        //! ----------------------------------------------------------------------------------------
        //! read the elements
        //! First line: <# of points> <dimension (3)> <# of attributes> <boundary markers (0 or 1)>
        //! ----------------------------------------------------------------------------------------
        unsigned int NbElements;
        unsigned int NbPoints;
        unsigned int region;
        unsigned int eNumber;
        unsigned int n1,n2,n3,n4,n5,n6,n7,n8,n9,n10;
        fscanf(tetgenEleFile,"%d%d%d",&NbElements,&NbPoints,&region);

        myNumberOfElements = NbElements;
        myElemType = new TColStd_HArray1OfInteger(1,NbElements);
        myElemNodes = new TColStd_HArray2OfInteger(1,NbElements,1,10);

        cout<<"____Number of volume elements: "<<NbElements<<"____"<<endl;

        bool isElementReadingOk = true;
        for(unsigned int line = 1; line<=NbElements; line++)
        {
            if(5==fscanf(tetgenEleFile,"%d%d%d%d%d",&eNumber,&n1,&n2,&n3,&n4))
            {
                //! --------------------------
                //! use Tetgen nodal ordering
                //! --------------------------
                myElemType->SetValue(eNumber,TET);
                myElemNodes->SetValue(eNumber,1,n1);
                myElemNodes->SetValue(eNumber,2,n2);
                myElemNodes->SetValue(eNumber,4,n3);
                myElemNodes->SetValue(eNumber,3,n4);

                myElements.Add(eNumber);
                myElementsMap.Add(eNumber);
            }
            else if(11==fscanf(tetgenEleFile,"%d%d%d%d%d%d%d%d%d%d%d",&eNumber,&n1,&n2,&n3,&n4,&n5,&n6,&n7,&n8,&n9,&n10))
            {
                //! --------------------------
                //! use Tetgen nodal ordering
                //! --------------------------
                myElemType->SetValue(eNumber,TET10);
                myElemNodes->SetValue(eNumber,1,n1);
                myElemNodes->SetValue(eNumber,2,n2);
                myElemNodes->SetValue(eNumber,3,n3);
                myElemNodes->SetValue(eNumber,4,n4);
                myElemNodes->SetValue(eNumber,5,n5);
                myElemNodes->SetValue(eNumber,6,n6);
                myElemNodes->SetValue(eNumber,7,n7);
                myElemNodes->SetValue(eNumber,8,n8);
                myElemNodes->SetValue(eNumber,9,n9);
                myElemNodes->SetValue(eNumber,10,n10);
            }
            else
            {
                isElementReadingOk = false;
                cout<<"____error in reading elements____"<<endl;
                break;
            }
        }

        if(isElementReadingOk==true)
        {
            cout<<"____Elements reading OK____"<<endl;
        }

        //! --------------------------------
        //! set the topological information
        //! --------------------------------
        if(isElementReadingOk==true)
        {
            cout<<"____Setting the topological information____"<<endl;
            //! -------------------------------------------
            //! topological information section
            //! sequence of nodes for the 3D visualization
            //! 1st order tet: definition of the faces
            //! Face 1: 1-2-3
            //! Face 2: 1-4-2
            //! Face 3: 2-4-3
            //! Face 4: 3-4-1
            //! -------------------------------------------
            TET4MeshData = new MeshVS_HArray1OfSequenceOfInteger(1,4);
            TColStd_SequenceOfInteger face1, face2, face3, face4;

            face1.Append(0); face1.Append(1); face1.Append(3);
            face2.Append(0); face2.Append(2); face2.Append(1);
            face3.Append(1); face3.Append(2); face3.Append(3);
            face4.Append(3); face4.Append(2); face4.Append(0);

            TET4MeshData->SetValue(1,face1);
            TET4MeshData->SetValue(2,face2);
            TET4MeshData->SetValue(3,face3);
            TET4MeshData->SetValue(4,face4);

            //! -------------------------------------------
            //! sequence of nodes for the 3D visualization
            //! 2nd order tet: definition of the faces
            //! Face 1: 0-6-1-8-3-5
            //! Face 2: 0-9-2-7-4-6
            //! Face 3: 1-7-2-4-3-8
            //! Face 4: 3-4-2-9-0-5
            //! --------------------------------------------
            TET10MeshData = new MeshVS_HArray1OfSequenceOfInteger(1,4);

            face1.Clear(); face2.Clear(); face3.Clear(); face4.Clear();
            face1.Append(0); face1.Append(6); face1.Append(1); face1.Append(8); face1.Append(3); face1.Append(5);
            face2.Append(0); face2.Append(9); face2.Append(2); face2.Append(7); face2.Append(4); face2.Append(6);
            face3.Append(1); face3.Append(7); face3.Append(2); face3.Append(4); face3.Append(3); face3.Append(8);
            face4.Append(3); face4.Append(4); face4.Append(2); face4.Append(9); face4.Append(0); face4.Append(5);

            TET10MeshData->SetValue(1,face1);
            TET10MeshData->SetValue(2,face2);
            TET10MeshData->SetValue(3,face3);
            TET10MeshData->SetValue(4,face4);
        }
        else
        {
            //! ---------------------------------------------
            //! reset the node maps and the node coordinates
            //! ---------------------------------------------
            myNodes.Clear();
            myNodesMap.Clear();
            myNodeCoords.Nullify();
            myElements.Clear();
            myElementsMap.Clear();
            myElemNodes.Nullify();
        }
    }
    cout<<"Ng_MeshVS_DataSource3D::Ng_MeshVS_DataSource3D()->____constructor from Tetgen files OK____"<<endl;
}

//! -----------------------------------------------------------------------------
//! function: constructor
//! details:  merge two mesh datasources: mesh 1 => prismatic mesh
//!                                       mesh 2 => interior, non prismatic mesh
//! -----------------------------------------------------------------------------
Ng_MeshVS_DataSource3D::Ng_MeshVS_DataSource3D(const occHandle(Ng_MeshVS_DataSource3D) &aMesh1,
                                               const occHandle(Ng_MeshVS_DataSource3D) &aMesh2)
{
    clock();

    //! -------------
    //! sanity check
    //! -------------
    if(aMesh1.IsNull() || aMesh2.IsNull()) return;
    if(aMesh1->GetAllElements().Extent()<1 || aMesh2->GetAllElements().Extent()<1)
    {
        cerr<<"Ng_MeshVS_DataSource3D::Ng_MeshVS_DataSource3D()->____less than one element for one mesh____"<<endl;
        return;
    }
    if(aMesh1->GetAllNodes().Extent()<3 || aMesh2->GetAllNodes().Extent()<3)
    {
        cerr<<"Ng_MeshVS_DataSource3D::Ng_MeshVS_DataSource3D()->____wrong number of nodes for one mesh____"<<endl;
        return;
    }

    //! ----------------------------------------------
    //! store the points of the first mesh into a map
    //! ----------------------------------------------
    QMap<mesh::meshPoint,int> nodeMap_mesh1;
    for(TColStd_MapIteratorOfPackedMapOfInteger it = aMesh1->GetAllNodes(); it.More(); it.Next())
    {
        int globalNodeID = it.Key();
        int NbNodes;
        MeshVS_EntityType type;
        double buf[3];
        TColStd_Array1OfReal coords(*buf,1,3);
        if(!aMesh1->GetGeom(globalNodeID,false,coords,NbNodes,type)) continue;
        mesh::meshPoint aP(coords(1),coords(2),coords(3));
        nodeMap_mesh1.insert(aP,globalNodeID);
    }

    //! ---------------------------------------------------------------------
    //! number of common nodes at the prismatic/non-prismatic mesh interface
    //! ---------------------------------------------------------------------
    int NbCommonNodes = 0;

    //! -------------------------------------------------------------------
    //! key:   (global) nodeID of mesh2 which has to be replaced
    //! value: (global) nodeID of mesh1, which must be used as replacement
    //! -------------------------------------------------------------------
    QMap<int,int> mapForRenumbering;

    //! ---------------------------------
    //! scan the nodes of the second map
    //! ---------------------------------
    for(TColStd_MapIteratorOfPackedMapOfInteger it = aMesh2->GetAllNodes(); it.More(); it.Next())
    {
        int globalNodeID = it.Key();
        int NbNodes;
        MeshVS_EntityType type;
        double buf[3];
        TColStd_Array1OfReal coords(*buf,1,3);
        if(!aMesh2->GetGeom(globalNodeID,false,coords,NbNodes,type)) continue;
        mesh::meshPoint aP(coords(1),coords(2),coords(3));

        //! -----------------------------
        //! a common node has been found
        //! -----------------------------
        if(nodeMap_mesh1.contains(aP))
        {
            NbCommonNodes++;
            int nodeIDForReplacement = nodeMap_mesh1.value(aP);
            mapForRenumbering.insert(globalNodeID,nodeIDForReplacement);
        }
    }

    //! ------------------
    //! test - disgnostic
    //! ------------------
    cout<<"____number of common nodes at mesh1/mesh2 interface: "<<NbCommonNodes<<"____"<<endl;
    cout<<"____table of replacements____"<<endl;
    cout<<"____mesh2    mesh1____"<<endl;
    for(QMap<int,int>::iterator it = mapForRenumbering.begin();it!=mapForRenumbering.end();it++)
    {
        cout<<"    "<<it.key()<<"    "<<it.value()<<endl;
    }
    //! ---------------
    //! end diagnostic
    //! ---------------
    clock();

    //! ----------------------------
    //! number of nodes and elements
    //! -----------------------------
    int NbNodes1 = aMesh1->GetAllNodes().Extent();
    int NbNodes2 = aMesh2->GetAllNodes().Extent();
    myNumberOfNodes = NbNodes1+NbNodes2-NbCommonNodes;

    cout<<"Ng_MeshVS_DataSource3D::Ng_MeshVS_DataSource3D()->____mesh1____number of nodes: "<<NbNodes1<<"____"<<endl;
    cout<<"Ng_MeshVS_DataSource3D::Ng_MeshVS_DataSource3D()->____mesh2____number of nodes: "<<NbNodes2<<"____"<<endl;
    cout<<"Ng_MeshVS_DataSource3D::Ng_MeshVS_DataSource3D()->____number of common nodes:   "<<NbCommonNodes<<"____"<<endl;

    //! -------------------
    //! number of elements
    //! -------------------
    int NbElements1 = aMesh1->GetAllElements().Extent();
    int NbElements2 = aMesh2->GetAllElements().Extent();

    cout<<"Ng_MeshVS_DataSource3D::Ng_MeshVS_DataSource3D()->____mesh1____number of elements: "<<NbElements1<<"____"<<endl;
    cout<<"Ng_MeshVS_DataSource3D::Ng_MeshVS_DataSource3D()->____mesh2____number of elements: "<<NbElements2<<"____"<<endl;

    myNumberOfElements = NbElements1+NbElements2;

    //! ------------------
    //! memory allocation
    //! ------------------
    myNodeCoords = new TColStd_HArray2OfReal(1,myNumberOfNodes,1,3);
    myElemType = new TColStd_HArray1OfInteger(1,myNumberOfElements);
    myElemNodes = new TColStd_HArray2OfInteger(1,myNumberOfElements,1,20);  // let the possibility to define an 2-nd order hexahedron

    //! --------------
    //! maps of nodes
    //! --------------
    int localNodeID = 0;
    TColStd_MapIteratorOfPackedMapOfInteger nodeIt;
    for(nodeIt.Initialize(aMesh1->GetAllNodes());nodeIt.More();nodeIt.Next())
    {
        localNodeID++;
        int globalNodeID = nodeIt.Key();
        myNodes.Add(globalNodeID);
        myNodesMap.Add(globalNodeID);

        MeshVS_EntityType type;
        double buf[3];
        TColStd_Array1OfReal coords(*buf,1,3);
        int NbNodes;
        if(aMesh1->GetGeom(globalNodeID,false,coords,NbNodes,type))
        {
            myNodeCoords->SetValue(localNodeID,1,coords(1));
            myNodeCoords->SetValue(localNodeID,2,coords(2));
            myNodeCoords->SetValue(localNodeID,3,coords(3));
        }
        else
        {
            cerr<<"Ng_MeshVS_DataSource3D::Ng_MeshVS_DataSource3D()->____mesh 1: local node ID: "<<localNodeID<<" not found____"<<endl;
            return;
        }
    }

    //! ------------------------------------------------------
    //! offset the global nodeIDs of the second mesh by the
    //! number of nodes of the first mesh minus the number of
    //! common nodes. The local node ID is simply increaded
    //! ------------------------------------------------------
    for(nodeIt.Initialize(aMesh2->GetAllNodes());nodeIt.More();nodeIt.Next())
    {
        int globalNodeID = nodeIt.Key();

        //! ----------------------------------------------------------------------
        //! check if the map for replacements contains the current global node ID
        //! if the current node is contained, it means that it has already been
        //! inserted before
        //! ----------------------------------------------------------------------
        if(mapForRenumbering.contains(globalNodeID)) continue;

        localNodeID++;

        myNodes.Add(globalNodeID+NbNodes1-NbCommonNodes);
        myNodesMap.Add(globalNodeID+NbNodes1-NbCommonNodes);

        MeshVS_EntityType type;
        double buf[3];
        int NbNodes;
        TColStd_Array1OfReal coords(*buf,1,3);

        if(aMesh2->GetGeom(globalNodeID,false,coords,NbNodes,type))
        {
            myNodeCoords->SetValue(localNodeID,1,coords(1));
            myNodeCoords->SetValue(localNodeID,2,coords(2));
            myNodeCoords->SetValue(localNodeID,3,coords(3));
        }
        else
        {
            cerr<<"Ng_MeshVS_DataSource3D::Ng_MeshVS_DataSource3D()->____mesh 2: local node ID: "<<localNodeID<<" not found____"<<endl;
            return;
        }
    }

    //! -----------------------------------
    //! add the elements of the first mesh
    //! -----------------------------------
    int localElementID = 0;
    TColStd_MapIteratorOfPackedMapOfInteger elementIt;
    for(elementIt.Initialize(aMesh1->GetAllElements()); elementIt.More(); elementIt.Next())
    {
        localElementID++;
        int globalElementID = elementIt.Key();
        myElements.Add(globalElementID);
        myElementsMap.Add(globalElementID);

        int NbNodes;
        int buf[20];
        TColStd_Array1OfInteger nodeIDs(*buf,1,20);
        if(aMesh1->GetNodesByElement(globalElementID,nodeIDs,NbNodes))
        {
            //! ---------------------
            //! set the element type
            //! ---------------------
            switch(NbNodes)
            {
            //! first order elements
            case 4: myElemType->SetValue(localElementID,TET); break;
            case 8: myElemType->SetValue(localElementID,HEXA); break;
            case 6: myElemType->SetValue(localElementID,PRISM); break;
            case 5: myElemType->SetValue(localElementID,PYRAM); break;

            //! second order elements
            case 10: myElemType->SetValue(localElementID,TET10); break;
            case 20: myElemType->SetValue(localElementID,HEXA20); break;
            case 15: myElemType->SetValue(localElementID,PRISM15); break;
            case 13: myElemType->SetValue(localElementID,PYRAM13); break;
            }
            //! -----------------------
            //! set the elements nodes
            //! -----------------------
            for(int n=1; n<=NbNodes; n++) myElemNodes->SetValue(localElementID,n,nodeIDs.Value(n));
        }
        else
        {
            cerr<<"Ng_MeshVS_DataSource3D::Ng_MeshVS_DataSource3D()->____mesh 2: global node ID: "<<globalElementID<<" not found____"<<endl;
            return;
        }
    }

    //! -----------
    //! diagnostic
    //! -----------
    QList<std::vector<int>> listOfReplacements;

    //! ------------------------------------
    //! add the elements of the second mesh
    //! ------------------------------------
    for(elementIt.Initialize(aMesh2->GetAllElements()); elementIt.More(); elementIt.Next())
    {
        localElementID++;
        int globalElementID = elementIt.Key();

        myElements.Add(globalElementID+NbElements1);
        myElementsMap.Add(globalElementID+NbElements1);

        int NbNodes;
        int buf[20];
        TColStd_Array1OfInteger nodeIDs(*buf,1,20);
        if(aMesh2->GetNodesByElement(globalElementID,nodeIDs,NbNodes))
        {
            //! ---------------------
            //! set the element type
            //! ---------------------
            switch(NbNodes)
            {
            //! first order elements
            case 4: myElemType->SetValue(localElementID,TET); break;
            case 8: myElemType->SetValue(localElementID,HEXA); break;
            case 6: myElemType->SetValue(localElementID,PRISM); break;
            case 5: myElemType->SetValue(localElementID,PYRAM); break;
            //! second order elements
            case 10: myElemType->SetValue(localElementID,TET10); break;
            case 20: myElemType->SetValue(localElementID,HEXA20); break;
            case 15: myElemType->SetValue(localElementID,PRISM15); break;
            case 13: myElemType->SetValue(localElementID,PYRAM13); break;
            }
            //! ----------------------
            //! set the element nodes
            //! ----------------------
            for(int n=1; n<=NbNodes; n++)
            {
                int globalNodeID = nodeIDs.Value(n);
                if(mapForRenumbering.contains(globalNodeID))
                {
                    //! -----------
                    //! diagnostic
                    //! -----------
                    std::vector<int> replacement;
                    replacement.push_back(globalNodeID);
                    //! ---------------
                    //! end diagnostic
                    //! ---------------

                    //! use the global node ID of the point within the first mesh
                    globalNodeID = mapForRenumbering.value(globalNodeID);

                    //! -----------
                    //! diagnostic
                    //! -----------
                    cout<<" with mesh1 nodeID: "<<globalNodeID<<"____"<<endl;
                    replacement.push_back(globalNodeID);
                    //replacement.push_back(mapForRenumbering.value(globalNodeID));
                    if(!listOfReplacements.contains(replacement))
                    {
                        //globalNodeID = mapForRenumbering.value(globalNodeID);
                        listOfReplacements<<replacement;
                    }
                    //! ---------------
                    //! end diagnostic
                    //! ---------------
                }
                else
                {
                    globalNodeID = nodeIDs.Value(n)+NbNodes1-NbCommonNodes;
                }
                myElemNodes->SetValue(localElementID,n,globalNodeID);
            }
        }
        else
        {
            cerr<<"Ng_MeshVS_DataSource3D::Ng_MeshVS_DataSource3D()->____mesh 2: global node ID: "<<globalElementID<<" not found____"<<endl;
            return;
        }
    }

    //! ------------------
    //! elements topology
    //! ------------------
    this->buildElementsTopology();

    clock();
}

//! ------------------------------
//! function: constructor
//! details:  from Eigen matrices
//! ------------------------------
Ng_MeshVS_DataSource3D::Ng_MeshVS_DataSource3D(const Eigen::MatrixXd &V, const Eigen::MatrixXi &T)
{
    //! -------------
    //! sanity check
    //! -------------
    if(V.rows()<4) return;
    if(T.rows()<1) return;

    myNumberOfNodes = V.rows();
    myNodeCoords = new TColStd_HArray2OfReal(1,myNumberOfNodes,1,3);
    for(int i = 1; i<=myNumberOfNodes; i++)
    {
        int row = i-1;
        myNodes.Add(i);
        myNodesMap.Add(i);
        myNodeCoords->SetValue(i,1,V(row,0));
        myNodeCoords->SetValue(i,2,V(row,1));
        myNodeCoords->SetValue(i,3,V(row,2));
    }

    myNumberOfElements = T.rows();
    myElemNodes = new TColStd_HArray2OfInteger(1,myNumberOfElements,1,4);
    myElemType = new TColStd_HArray1OfInteger(1,myNumberOfElements);
    for(int i=1; i<=myNumberOfElements; i++)
    {
        myElements.Add(i);
        myElementsMap.Add(i);
        myElemType->SetValue(i,TET);
        int row = i-1;
        myElemNodes->SetValue(i,1,T(row,0)+1);
        myElemNodes->SetValue(i,2,T(row,1)+1);
        myElemNodes->SetValue(i,3,T(row,2)+1);
        myElemNodes->SetValue(i,4,T(row,3)+1);
    }

    //! -----------------
    //! element topology
    //! -----------------
    this->buildElementsTopology();
}

//! --------------------------------------------------------------
//! function: Get3DGeom
//! details:  returns information about 3D geometry of an element
//! --------------------------------------------------------------
Standard_Boolean Ng_MeshVS_DataSource3D::Get3DGeom(const Standard_Integer theID,
                                                   Standard_Integer &theNbNodes,
                                                   occHandle(MeshVS_HArray1OfSequenceOfInteger) &theData) const
{
    //!cout<<"Ng_MeshVS_DataSource3D::Get3DGeom()->____Get3DGeom called for ID: "<<theID<<"____"<<endl;
    int localElementID = myElementsMap.FindIndex(theID);
    if(localElementID>=1 && localElementID<=myNumberOfElements)
    {
        switch(myElemType->Value(localElementID))
        {
        //! -----------------
        //! Surface elements
        //! -----------------
        case TRIG: theNbNodes = 3; theData = TRIGMeshData; break;
        case QUAD: theNbNodes = 4;  theData = QUADMeshData; break;
        case TRIG6: theNbNodes = 6; theData = TRIG6MeshData; break;
        case QUAD8: theNbNodes = 8; theData = QUAD8MeshData; break;
        //! ----------------
        //! Volume elements
        //! ----------------
        case TET: theNbNodes = 4; theData = TET4MeshData; break;
        case HEXA: theNbNodes = 8; theData = HEXA8MeshData; break;
        case PRISM: theNbNodes = 6; theData = PRISM6MeshData; break;
        case PYRAM: theNbNodes = 5; theData = PYRAM5MeshData; break;
        case TET10: theNbNodes = 10; theData= TET10MeshData; break;
        case HEXA20: theNbNodes = 20; theData = HEXA20MeshData; break;
        case PRISM15: theNbNodes = 15; theData = PRISM15MeshData; break;
        case PYRAM13: theNbNodes = 13; theData = PYRAM13MeshData; break;
        }
        return true;
    }
    else
    {
        cerr<<"Ng_MeshVS_DataSource3D::Get3DGeom()->____error. Element ID: "<<theID<<" out of range____"<<endl;
        return false;
    }
}

//! ------------------
//! function: GetGeom
//! details:
//! ------------------
Standard_Boolean Ng_MeshVS_DataSource3D::GetGeom(const Standard_Integer ID,
                                                 const Standard_Boolean IsElement,
                                                 TColStd_Array1OfReal& Coords,
                                                 Standard_Integer& NbNodes,
                                                 MeshVS_EntityType& Type) const
{
    if(IsElement)
    {
        //cout<<"____tag000: globalElementID: "<<ID<<"____"<<endl;
        int localElementID = myElementsMap.FindIndex(ID);
        //cout<<"____tag001: localElementID: "<<localElementID<<"____"<<endl;

        if(localElementID>=1 && localElementID<=myNumberOfElements)
        {
            //cout<<"____tag002: in range____"<<endl;

            Type = MeshVS_ET_Volume;
            switch(myElemType->Value(localElementID))
            {
            case TET:
            {
                NbNodes=4;
                //cout<<"____tag003: TET found____"<<endl;
            }
                break;
            case TET10: NbNodes=10; break;
            case HEXA: NbNodes = 8; break;
            case HEXA20: NbNodes = 20; break;
            case PYRAM: NbNodes = 5; break;
            case PYRAM13: NbNodes = 13; break;
            case PRISM: NbNodes = 6; break;
            case PRISM15: NbNodes = 15; break;
            }

            for(int i=1, k=1;i<=NbNodes;i++)
            {
                //cout<<"____tag004: i = "<<i<<"____"<<endl;

                int nodeID = myElemNodes->Value(localElementID,i);

                //cout<<"____tag005: globalNodeID: "<<nodeID<<"____"<<endl;

                for(int j=1; j<=3; j++)
                {
                    //cout<<"____tag006____"<<endl;

                    int localNodeID = myNodesMap.FindIndex(nodeID);
                    /*
                        int globalNodeID = myElemNodes->Value(localElementID,i);
                        int localNodeID = myNodesMap.FindIndex(globalNodeID);
                     */

                    Coords(k) = myNodeCoords->Value(localNodeID,j);
                    //Coords(k) = myNodeCoords->Value(nodeID,j);

                    //cout<<"____tag007____"<<endl;

                    k++;
                }
            }
            return true;
        }
        else
        {
            cerr<<"Ng_MeshVS_DataSource3D::GetGeom()->____error. Element ID: "<<ID<<" is out of range____"<<endl;
            return false;
        }
    }
    else
    {
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
            cerr<<"Ng_MeshVS_DataSource3D::GetGeom()->____error. Node ID: "<<ID<<" is out of range____"<<endl;
            return false;
        }
    }
}


//! -------------------------------------------------------------------
//! function: GetGeomType
//! details:  the input argument "ID" should be considered as "global"
//! -------------------------------------------------------------------
Standard_Boolean Ng_MeshVS_DataSource3D::GetGeomType(const Standard_Integer ID, const Standard_Boolean IsElement, MeshVS_EntityType& Type) const
{
    if(IsElement)
    {
        int localElementID = myElementsMap.FindIndex(ID);
        if(localElementID>=1 && localElementID<=myNumberOfElements)
        {
            Type = MeshVS_ET_Volume;
            return true;
        }
        cerr<<"Ng_MeshVS_DataSource3D::GetGeomType()->____element ID: "<<ID<<"out of range____"<<endl;
        return false;
    }
    else
    {
        int localNodeID = myNodesMap.FindIndex(ID);
        if(localNodeID>=1 && localNodeID<=myNumberOfNodes)
        {
            Type = MeshVS_ET_Node;
            return true;
        }
        cerr<<"Ng_MeshVS_DataSource3D::GetGeomType()->____node ID: "<<ID<<"out of range____"<<endl;
        return false;
    }
}

//! ------------------
//! function: GetAddr
//! details:
//! ------------------
Standard_Address Ng_MeshVS_DataSource3D::GetAddr(const Standard_Integer, const Standard_Boolean) const
{
    return NULL;
}

//! ----------------------------
//! function: GetNodesByElement
//! details:
//! ----------------------------
Standard_Boolean Ng_MeshVS_DataSource3D::GetNodesByElement(const Standard_Integer ID,
                                                           TColStd_Array1OfInteger& theNodeIDs,
                                                           Standard_Integer& theNbNodes) const
{
    //!cout<<"Ng_MeshVS_DataSource3D::GetNodesByElement()->____function called for ID: "<<ID<<"____"<<endl;
    int localElementID = myElementsMap.FindIndex(ID);
    if(localElementID>=1 && localElementID<=myNumberOfElements)
    {
        int eType = myElemType->Value(localElementID);
        switch(eType)
        {
        case TET: theNbNodes=4; break;
        case TET10: theNbNodes = 10; break;
        case HEXA : theNbNodes = 8; break;
        case HEXA20: theNbNodes = 20; break;
        case PRISM: theNbNodes = 6; break;
        case PRISM15: theNbNodes = 15; break;
        case PYRAM: theNbNodes = 5; break;
        case PYRAM13: theNbNodes = 13; break;
        }

        //Standard_Integer aLow = theNodeIDs.Lower();
        //for(int i=0;i<theNbNodes;i++)
        for(int i=1;i<=theNbNodes;i++)
        {
            theNodeIDs(i) = myElemNodes->Value(localElementID,i);
            //theNodeIDs(aLow+i)=myElemNodes->Value(localElementID,i+1);
        }
        return true;
    }
    else
    {
        cerr<<"Ng_MeshVS_DataSource3D::GetNodesByElement()->____function called for ID: "<<ID<<" out of range____"<<endl;
        return false;
    }
}

//! ----------------------
//! function: GetAllNodes
//! details:
//! ----------------------
const TColStd_PackedMapOfInteger& Ng_MeshVS_DataSource3D::GetAllNodes() const
{
    return myNodes;
}

//! -------------------------
//! function: GetAllElements
//! details:
//! -------------------------
const TColStd_PackedMapOfInteger& Ng_MeshVS_DataSource3D::GetAllElements() const
{
    return myElements;
}

//! -------------------------
//! function: GetElementType
//! details:
//! -------------------------
const bool Ng_MeshVS_DataSource3D::GetElementType(ElemType &eType, int elementID, bool isLocal) const
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

//! --------------------------------------------------
//! function: GetNormalsByElement
//! details:  never called by the presentation engine
//! --------------------------------------------------
Standard_Boolean Ng_MeshVS_DataSource3D::GetNormalsByElement (const Standard_Integer Id,
                                                              const Standard_Boolean IsNodal,
                                                              const Standard_Integer MaxNodes,
                                                              occHandle(TColStd_HArray1OfReal)& Normals) const
{
    cout<<"Ng_MeshVS_DataSource3D::GetNormalsByElement()->____function called____"<<endl;
    exit(1000);

    MeshVS_Buffer aCoordsBuf(3*MaxNodes*sizeof(Standard_Real));
    TColStd_Array1OfReal Coords(aCoordsBuf,1,3*MaxNodes);
    Standard_Integer NbNodes;
    MeshVS_EntityType Type;

    Standard_Boolean res = Standard_False;

    //! wrong number of nodes: return false
    if(MaxNodes<= 0) return res;

    //! cannot retrieve the geometry of an element: return false
    if(!GetGeom(Id,Standard_True,Coords,NbNodes,Type)) return res;

    //! a normal for each node
    Standard_Integer aNbNormals = NbNodes;

    occHandle(MeshVS_HArray1OfSequenceOfInteger) aTopo;
    if (Type==MeshVS_ET_Volume)
    {
        //! cannot retrieve the 3D geometry (array of sequence of integers
        //! defining the faces)
        if(!Get3DGeom(Id,NbNodes,aTopo)) return res;

        //! calculate number of normals for faces of volume
        aNbNormals = aTopo->Upper() - aTopo->Lower() + 1;
    }

    //! 3 times (nx,ny,nz) values
    occHandle(TColStd_HArray1OfReal) aNormals = new TColStd_HArray1OfReal(1,3*aNbNormals);

    //! -------------------------------------------------------------
    //! Type == MeshVS_Face && isNodal == true => allNormals = false
    //! Type == MeshVS_Node && isNodal == true => allNormals = true
    //! -------------------------------------------------------------
    Standard_Boolean allNormals = ( Type == MeshVS_ET_Face && IsNodal );

    //! ----------------------------------------
    //! Try returning nodal normals if possible
    //! ----------------------------------------
    for(Standard_Integer k=1; k<=NbNodes && allNormals; k++)
        allNormals = GetNodeNormal(k,
                                   Id,
                                   aNormals->ChangeValue(3*k-2),
                                   aNormals->ChangeValue(3*k-1),
                                   aNormals->ChangeValue(3*k  ));

    //! ------------------------------------------
    //! Nodal normals not available or not needed
    //! ------------------------------------------
    if(!allNormals)
    {
        switch(Type)
        {
        //! ---------------------------------------------------------------
        //! Compute a face normal and duplicate it for all element`s nodes
        //! ---------------------------------------------------------------
        case MeshVS_ET_Face:
        {

            res = GetNormal(Id,
                            MaxNodes,
                            aNormals->ChangeValue(1),
                            aNormals->ChangeValue(2),
                            aNormals->ChangeValue(3));
            if(res)
            {
                for( Standard_Integer k=2; k<=NbNodes; k++ )
                {
                    aNormals->ChangeValue(3 * k - 2) = aNormals->Value(1);
                    aNormals->ChangeValue(3 * k - 1) = aNormals->Value(2);
                    aNormals->ChangeValue(3 * k    ) = aNormals->Value(3);
                }
            }
        }
            break;

            //! ---------------------------------------------------------------
            //! Compute normals for each of volum`s faces - not for each node!
            //! ---------------------------------------------------------------
        case MeshVS_ET_Volume:
        {
            gp_Vec norm;
            Standard_Integer low = Coords.Lower();
            for(Standard_Integer k = aTopo->Lower(), last = aTopo->Upper(), i = 1; k <= last; k++, i++)
            {
                //! the sequence of integers defining the k-th face
                const TColStd_SequenceOfInteger& aSeq = aTopo->Value(k);
                Standard_Integer m = aSeq.Length(), ind;

                norm.SetCoord(0, 0, 0);
                MeshVS_Buffer PolyNodesBuf (3*m*sizeof(Standard_Real));
                TColStd_Array1OfReal PolyNodes(PolyNodesBuf, 0, 3*m);
                PolyNodes.SetValue(0, m);
                for(Standard_Integer j=1; j<=m; j++)
                {
                    ind = aSeq.Value( j );
                    PolyNodes.SetValue( 3*j-2, Coords( low+3*ind   ));
                    PolyNodes.SetValue( 3*j-1, Coords( low+3*ind+1 ));
                    PolyNodes.SetValue( 3*j,   Coords( low+3*ind+2 ));
                }

                //! get the average normal of a polygon
                MeshVS_Tool::GetAverageNormal( PolyNodes, norm );

                aNormals->ChangeValue(i * 3 - 2) = norm.X();
                aNormals->ChangeValue(i * 3 - 1) = norm.Y();
                aNormals->ChangeValue(i * 3    ) = norm.Z();
            }
            res = Standard_True;
        }
            break;

        default:
            return res;
        }
    }
    if (res || allNormals) Normals = aNormals;
    return (res || allNormals);
}

//! --------------------------------------------------
//! function: GetNormal
//! details:  never called by the presentation engine
//! --------------------------------------------------
Standard_Boolean Ng_MeshVS_DataSource3D::GetNormal(const Standard_Integer Id,
                                                   const Standard_Integer Max,
                                                   Standard_Real& nx,
                                                   Standard_Real& ny,
                                                   Standard_Real& nz) const
{
    cout<<"Ng_MeshVS_DataSource3D::GetNormal()->____funcion called____"<<endl;
    exit(1000);
    if(Max<=0) return Standard_False;

    MeshVS_Buffer aCoordsBuf (3*Max*sizeof(Standard_Real));
    TColStd_Array1OfReal Coords (aCoordsBuf,1,3*Max);
    Standard_Integer nbNodes;
    MeshVS_EntityType Type;

    Standard_Boolean res = Standard_False;

    if (!GetGeom (Id, Standard_True, Coords, nbNodes, Type))
        return res;

    if (Type == MeshVS_ET_Face && nbNodes >= 3)
    {
        Standard_Real x1 = Coords(1);
        Standard_Real y1 = Coords(2);
        Standard_Real z1 = Coords(3);
        Standard_Real x2 = Coords(4);
        Standard_Real y2 = Coords(5);
        Standard_Real z2 = Coords(6);
        Standard_Real x3 = Coords((nbNodes-1)*3+1);
        Standard_Real y3 = Coords((nbNodes-1)*3+2);
        Standard_Real z3 = Coords((nbNodes-1)*3+3);
        Standard_Real p1 = x2 - x1, p2 = y2 - y1, p3 = z2 - z1,
                q1 = x3 - x1, q2 = y3 - y1, q3 = z3 - z1;

        nx = p2*q3 - p3*q2;
        ny = p3*q1 - p1*q3;
        nz = p1*q2 - p2*q1;

        Standard_Real len = sqrt (nx*nx + ny*ny + nz*nz);
        if (len <= gp::Resolution())
        {
            nx = ny = nz = 0;
            return res;
        }
        nx /= len; ny /= len; nz /= len;
        res = Standard_True;
    }
    return res;
}

//! --------------------------------------------------------------
//! function: writeMesh
//! details:  "dim" mesh dimension value is redundant, since
//!           GetGeomType() is called foreach element in the mesh
//! --------------------------------------------------------------
bool Ng_MeshVS_DataSource3D::writeMesh(const QString &meshFileName, int dim)
{
    //!cout<<"Ng_MeshVS_DataSource3D::writeMesh()->____function called____"<<endl;
    char name[128];
    sprintf(name,"%s",meshFileName.toStdString().c_str());
    FILE *fileout = fopen(name,"w");

    if(fileout!=NULL)
    {
        //! nodes, elements
        const TColStd_PackedMapOfInteger& aNodes = this->GetAllNodes();
        const TColStd_PackedMapOfInteger &mapOfElements = this->GetAllElements();

        //! write the dimension
        fprintf(fileout,"%d\n",dim);
        fprintf(fileout,"%d\n",aNodes.Extent());

        double aCoordsBuf[3];
        TColStd_Array1OfReal aCoords(*aCoordsBuf,1,3);
        Standard_Integer nbNodes;
        MeshVS_EntityType aType;
        TColStd_MapIteratorOfPackedMapOfInteger anIter;

        //! ------------------------------------
        //! write the nodID and the coordinates
        //! ------------------------------------
        int k=1;
        for(anIter.Initialize(aNodes);anIter.More();anIter.Next())
        {
            Standard_Integer aKey = anIter.Key();
            if(this->GetGeom(k,false,aCoords,nbNodes,aType)==true)
            {
                fprintf(fileout,"%d\t%.9e\t%.9e\t%.9e\n",aKey,aCoords.Value(1),aCoords.Value(2),aCoords.Value(3));
                k++;
            }
        }

        //! -----------------------------
        //! write the number of elements
        //! -----------------------------
        fprintf(fileout,"%d\n",mapOfElements.Extent());

        for(anIter.Initialize(mapOfElements);anIter.More();anIter.Next())
        {
            Standard_Integer aKey = anIter.Key();

            int NbNodes;
            Standard_Integer aNodesBuf[20];
            TColStd_Array1OfInteger nodeIDs(*aNodesBuf,1,20);
            this->GetNodesByElement(aKey,nodeIDs,NbNodes);

            MeshVS_EntityType type;
            this->GetGeomType(aKey,true,type);

            char header[16];
            if(type==MeshVS_ET_Face)
            {
                switch(NbNodes)
                {
                case 3: sprintf(header,"TRIG"); break;
                case 4: sprintf(header,"QUAD"); break;
                case 6: sprintf(header,"TRIG6"); break;
                case 8: sprintf(header,"QUAD8"); break;
                default: break;
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
                case 6: sprintf(header,"PRISM"); break;
                case 5: sprintf(header,"PYRAM"); break;
                default: break;
                }
            }

            //! ----------------------------------------
            //! write the element key, the element type
            //! ----------------------------------------
            fprintf(fileout,"%d\n",aKey);
            fprintf(fileout,"%s\n",header);

            for(int i=1; i<NbNodes; i++) fprintf(fileout,"%d\t",nodeIDs(i));
            fprintf(fileout,"%d\n",nodeIDs(NbNodes));
        }
        fclose(fileout);
        return true;
    }
    else return false;
}

//! ---------------------------
//! function: readFromStream
//! details:  read from stream
//! ---------------------------
bool Ng_MeshVS_DataSource3D::readFromStream(std::ifstream &stream,
                                            std::vector<mesh::meshPoint> &nodes,
                                            std::vector<mesh::meshElement> &elements)
{
    if(!stream.is_open()) return false;

    //! -------------------------------------------
    //! read the dimension and the number of nodes
    //! -------------------------------------------
    int dim, NbNodes;
    stream>>dim;
    stream>>NbNodes;
    if(NbNodes==0) return false;
    if(dim == 3 && NbNodes<4) return false;

    //! ---------------------------
    //! read the nodes coordinates
    //! ---------------------------
    mesh::meshPoint aNode;
    for(int i=0; i<NbNodes;i++)
    {
        stream>>aNode.ID>>aNode.x>>aNode.y>>aNode.z;
        nodes.push_back(aNode);
    }

    //! ----------------------------
    //! read the number of elements
    //! ----------------------------
    int NbElements;
    stream>>NbElements;
    if(NbElements==0) return false;

    //! -----------------------------
    //! read the element definitions
    //! -----------------------------
    for(int i=0; i<NbElements; i++)
    {
        mesh::meshElement anElement;
        char eType[16];

        stream>>anElement.ID;
        stream>>eType;

        if(strcmp(eType,"TRIG")==0)
        {
            anElement.type= TRIG;
            int V1,V2,V3;
            stream>>V1>>V2>>V3;
            anElement.theNodeIDs.push_back(V1);
            anElement.theNodeIDs.push_back(V2);
            anElement.theNodeIDs.push_back(V3);
        }
        if(strcmp(eType,"TRIG6")==0)
        {
            anElement.type= TRIG6;
            int V1,V2,V3,V4,V5,V6;
            stream>>V1>>V2>>V3>>V4>>V5>>V6;
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
            stream>>V1>>V2>>V3>>V4;
            anElement.theNodeIDs.push_back(V1);
            anElement.theNodeIDs.push_back(V2);
            anElement.theNodeIDs.push_back(V3);
            anElement.theNodeIDs.push_back(V4);
        }
        else if(strcmp(eType,"TET10")==0)
        {
            anElement.type= TET10;
            int V1,V2,V3,V4,V5,V6,V7,V8,V9,V10;
            stream>>V1>>V2>>V3>>V4>>V5>>V6>>V7>>V8>>V9>>V10;
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
        else if(strcmp(eType,"PYRAM")==0)
        {
            anElement.type= PYRAM;
            int V1,V2,V3,V4,V5;
            stream>>V1>>V2>>V3>>V4>>V5;
            anElement.theNodeIDs.push_back(V1);
            anElement.theNodeIDs.push_back(V2);
            anElement.theNodeIDs.push_back(V3);
            anElement.theNodeIDs.push_back(V4);
            anElement.theNodeIDs.push_back(V5);
        }
        else if(strcmp(eType,"PRISM")==0)
        {
            anElement.type= PRISM;
            int V1,V2,V3,V4,V5,V6;
            stream>>V1>>V2>>V3>>V4>>V5>>V6;
            anElement.theNodeIDs.push_back(V1);
            anElement.theNodeIDs.push_back(V2);
            anElement.theNodeIDs.push_back(V3);
            anElement.theNodeIDs.push_back(V4);
            anElement.theNodeIDs.push_back(V5);
            anElement.theNodeIDs.push_back(V6);
        }
        elements.push_back(anElement);
    }
    return true;
}

//! ---------------------------
//! function: readMesh
//! details:  read a mesh file
//! ---------------------------
bool Ng_MeshVS_DataSource3D::readMesh(const QString &meshFileName,
                                      std::vector<mesh::meshPoint> &nodes,
                                      std::vector<mesh::meshElement> &elements)
{
    FILE *filein = fopen(meshFileName.toStdString().c_str(),"r");
    if(filein==NULL) return false;

    //! -------------------------------------------
    //! read the dimension and the number of nodes
    //! -------------------------------------------
    int dim, NbNodes;
    fscanf(filein,"%d",&dim);
    fscanf(filein,"%d",&NbNodes);

    //! ---------------------------
    //! read the nodes coordinates
    //! ---------------------------
    mesh::meshPoint aNode;
    for(int i=0; i<NbNodes;i++)
    {
        fscanf(filein,"%d%lf%lf%lf",&aNode.ID,&aNode.x,&aNode.y,&aNode.z);
        nodes.push_back(aNode);
    }

    //! ----------------------------
    //! read the number of elements
    //! ----------------------------
    int NbElements;
    fscanf(filein,"%d",&NbElements);

    //! -----------------------------
    //! read the element definitions
    //! -----------------------------
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
                cout<<"____error in reading the volume mesh____";
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
                cout<<"____error in reading the volume mesh____";
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
                cout<<"____error in reading the volume mesh____";
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
                cout<<"____error in reading the volume mesh____";
                return false;
            }
            anElement.theNodeIDs.push_back(V1);
            anElement.theNodeIDs.push_back(V2);
            anElement.theNodeIDs.push_back(V3);
            anElement.theNodeIDs.push_back(V4);
        }
        else if(strcmp(eType,"HEXA")==0)
        {
            anElement.type = HEXA;
            int V1,V2,V3,V4,V5,V6,V7,V8;
            if(8!=fscanf(filein,"%d%d%d%d%d%d%d%d",&V1,&V2,&V3,&V4,&V5,&V6,&V7,&V8))
            {
                cout<<"____error in reading the volume mesh____";
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
                cout<<"____error in reading the volume mesh____";
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
        else if(strcmp(eType,"PYRAM")==0)
        {
            anElement.type= PYRAM;
            int V1,V2,V3,V4,V5;
            if(5!=fscanf(filein,"%d%d%d%d%d",&V1,&V2,&V3,&V4,&V5))
            {
                cout<<"____error in reading the volume mesh____";
                return false;
            }
            anElement.theNodeIDs.push_back(V1);
            anElement.theNodeIDs.push_back(V2);
            anElement.theNodeIDs.push_back(V3);
            anElement.theNodeIDs.push_back(V4);
            anElement.theNodeIDs.push_back(V5);
        }
        else if(strcmp(eType,"PRISM")==0)
        {
            anElement.type= PRISM;
            int V1,V2,V3,V4,V5,V6;
            if(6!=fscanf(filein,"%d%d%d%d%d%d",&V1,&V2,&V3,&V4,&V5,&V6))
            {
                cout<<"____error in reading the volume mesh____";
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
    fclose(filein);
    return true;
}

//! ---------------------------
//! function: writeIntoStream
//! details:  use a C++ stream
//! ---------------------------
void Ng_MeshVS_DataSource3D::writeIntoStream(ofstream &os)
{
    if(!os.is_open()) return;

    cout<<"Ng_MeshVS_DataSource3D::writeIntoStream()->____function called____"<<endl;

    //! ----------------
    //! nodes, elements
    //! ----------------
    const TColStd_PackedMapOfInteger &aNodes = this->GetAllNodes();
    const TColStd_PackedMapOfInteger &mapOfElements = this->GetAllElements();

    if(aNodes.Extent()==0)
    {
        cout<<"Ng_MeshVS_DataSource3D::writeIntoStream()->____no node to write____"<<endl;
        return;
    }
    if(mapOfElements.Extent()==0)
    {
        cout<<"Ng_MeshVS_DataSource3D::writeIntoStream()->____no element to write____"<<endl;
        return;
    }

    //! ---------------
    //! mesh dimension
    //! ---------------
    os<<int(3)<<endl;
    os<<aNodes.Extent()<<endl;

    Standard_Real aCoordsBuf[3];
    TColStd_Array1OfReal aCoords(*aCoordsBuf,1,3);
    Standard_Integer nbNodes;
    MeshVS_EntityType aType;
    TColStd_MapIteratorOfPackedMapOfInteger anIter;

    //! ------------------------------------
    //! write the nodID and the coordinates
    //! ------------------------------------
    int k=1;
    for(anIter.Initialize(aNodes);anIter.More();anIter.Next())
    {
        int aKey = anIter.Key();
        if(this->GetGeom(k,Standard_False,aCoords,nbNodes,aType)==true)
        {
            os<<aKey<<"\t"<<aCoords(1)<<"\t"<<aCoords(2)<<"\t"<<aCoords(3)<<endl;
            k++;
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

        MeshVS_EntityType type;
        this->GetGeomType(aKey,true,type);

        char header[16];
        if(type==MeshVS_ET_Face)
        {
            switch(NbNodes)
            {
            case 3: sprintf(header,"TRIG"); break;
            case 4: sprintf(header,"QUAD"); break;
            case 6: sprintf(header,"TRIG6"); break;
            case 8: sprintf(header,"QUAD8"); break;
            default: break;
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
            case 6: sprintf(header,"PRISM"); break;
            case 5: sprintf(header,"PYRAM"); break;
            default: break;
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
}

//! -----------------------------------------
//! function: buildFaceToElementConnectivity
//! details:  optimization => use QHash
//! -----------------------------------------
void Ng_MeshVS_DataSource3D::buildFaceToElementConnectivity()
{
    //cout<<"Ng_MeshVS_DataSource3D::buildFaceToElementConnectivity()->____generating face to element connectivity data____"<<endl;
    myFaceToElements.clear();
    mySurfaceElements.clear();

    for(int localElementID=1; localElementID<=myNumberOfElements; localElementID++)
    {
        ElemType type = static_cast<ElemType>(myElemType->Value(localElementID));
        occHandle(MeshVS_HArray1OfSequenceOfInteger) topology;
        int NbFaces;
        switch(type)
        {
        //! ---------------------
        //! first order elements
        //! ---------------------
        case TET: topology = TET4MeshData; NbFaces = 4; break;
        case PYRAM: topology = PYRAM5MeshData; NbFaces = 5; break;
        case PRISM: topology = PRISM6MeshData; NbFaces = 5; break;
        case HEXA: topology = HEXA8MeshData; NbFaces = 6; break;

        //! ----------------------
        //! second order elements
        //! ----------------------
        case TET10: topology = TET10MeshData; NbFaces = 4; break;
        case HEXA20: topology = HEXA20MeshData; NbFaces = 6; break;
        case PRISM15: topology = PRISM15MeshData; NbFaces = 5; break;
        case PYRAM13: topology = PYRAM13MeshData; NbFaces = 5; break;
        }

        for(int f=1; f<=NbFaces; f++)
        {
            TColStd_SequenceOfInteger indices = topology->Value(f);
            //! -------------------------
            //! create the meshElement2D
            //! -------------------------
            meshElement2D aMeshElement2D;
            int NbIndices = indices.Length();
            for(int n=1; n<=NbIndices; n++)
            {
                int index = indices.Value(n);
                aMeshElement2D.nodeIDs<<myElemNodes->Value(localElementID,index+1);
            }

            //! -----------------------------------------
            //! ::contains() uses "<" in the "weak" form
            //! not aware of nodal ordering
            //! -----------------------------------------
            if(myFaceToElements.contains(aMeshElement2D))
            {
                QList<int> elements = myFaceToElements.value(aMeshElement2D);
                elements<<localElementID;
                myFaceToElements.insert(aMeshElement2D,elements);
            }
            else
            {
                QList<int> elements;
                elements<<localElementID;
                myFaceToElements.insert(aMeshElement2D,elements);
            }
        }
    }

    //! -----------------------------
    //! collect the surface elements
    //! -----------------------------
    int NbInternalElements = 0;
    int NbSurfaceElements = 0;
    for(QMap<meshElement2D,QList<int>>::iterator it = myFaceToElements.begin(); it!=myFaceToElements.end(); ++it)
    {
        QList<int> elements = it.value();
        if(elements.length()==2) NbInternalElements++;
        else
        {
            NbSurfaceElements++;
            const meshElement2D &aMeshElement2D = it.key();
            mySurfaceElements.push_back(aMeshElement2D);
        }
        if(elements.length()==0 || elements.length()>3)
            cerr<<"____something very strange in topological surface mesh generation: "<<elements.length()<<"____"<<endl;
    }
}

//! --------------------------------
//! function: buildElementsTopology
//! details:
//! --------------------------------
void Ng_MeshVS_DataSource3D::buildElementsTopology()
{
    //! ----------------
    //! Linear triangle
    //! ----------------
    TRIGMeshData = new MeshVS_HArray1OfSequenceOfInteger(1,1);
    TRIGMeshData->ChangeValue(1).Append(0);
    TRIGMeshData->ChangeValue(1).Append(1);
    TRIGMeshData->ChangeValue(1).Append(2);

    //! ------------------
    //! Linear quadrangle
    //! ------------------
    QUADMeshData = new MeshVS_HArray1OfSequenceOfInteger(1,1);
    QUADMeshData->ChangeValue(1).Append(0);
    QUADMeshData->ChangeValue(1).Append(1);
    QUADMeshData->ChangeValue(1).Append(2);
    QUADMeshData->ChangeValue(1).Append(3);

    //! -------------------
    //! Quadratic triangle
    //! -------------------
    TRIG6MeshData = new MeshVS_HArray1OfSequenceOfInteger(1,1);
    TRIG6MeshData->ChangeValue(1).Append(0);
    TRIG6MeshData->ChangeValue(1).Append(3);
    TRIG6MeshData->ChangeValue(1).Append(1);
    TRIG6MeshData->ChangeValue(1).Append(4);
    TRIG6MeshData->ChangeValue(1).Append(2);
    TRIG6MeshData->ChangeValue(1).Append(5);

    //! ---------------------
    //! Quadratic quadrangle
    //! ---------------------
    QUAD8MeshData = new MeshVS_HArray1OfSequenceOfInteger(1,1);
    QUAD8MeshData->ChangeValue(1).Append(0);
    QUAD8MeshData->ChangeValue(1).Append(4);
    QUAD8MeshData->ChangeValue(1).Append(1);
    QUAD8MeshData->ChangeValue(1).Append(5);
    QUAD8MeshData->ChangeValue(1).Append(2);
    QUAD8MeshData->ChangeValue(1).Append(6);
    QUAD8MeshData->ChangeValue(1).Append(3);
    QUAD8MeshData->ChangeValue(1).Append(7);

    //! -------------------
    //! Linear tetrahedron
    //! -------------------
    TET4MeshData = new MeshVS_HArray1OfSequenceOfInteger(1,4);
    for(int i=1;i<=4;i++)
    {
        TET4MeshData->ChangeValue(i).Append((i-1)%4);
        TET4MeshData->ChangeValue(i).Append(i%4);
        TET4MeshData->ChangeValue(i).Append((i+1)%4);
    }

    //! ---------------
    //! Linear pyramid
    //! ---------------
    int NbBasePoints = 4;
    PYRAM5MeshData = new MeshVS_HArray1OfSequenceOfInteger(1,NbBasePoints+1);
    for(int i=1; i<=NbBasePoints; i++)
    {
        PYRAM5MeshData->ChangeValue(1).Prepend(i);
        PYRAM5MeshData->ChangeValue(1+i).Append(0);
        PYRAM5MeshData->ChangeValue(1+i).Append(i);
        PYRAM5MeshData->ChangeValue(1+i).Append(i%NbBasePoints+1);
    }

    //! -------------
    //! Linear prism
    //! -------------
    NbBasePoints = 3;
    int Nseq = NbBasePoints + 2;
    PRISM6MeshData = new MeshVS_HArray1OfSequenceOfInteger(1,Nseq);
    int i, next;
    for(i=0; i<NbBasePoints; i++)
    {
        PRISM6MeshData->ChangeValue(1).Prepend(i);
        PRISM6MeshData->ChangeValue(2).Append(i+NbBasePoints);
        PRISM6MeshData->ChangeValue(3+i).Prepend(i);
        PRISM6MeshData->ChangeValue(3+i).Prepend(i+NbBasePoints);
        next = (i+1)%NbBasePoints;
        PRISM6MeshData->ChangeValue(3+i).Prepend(next+NbBasePoints);
        PRISM6MeshData->ChangeValue(3+i).Prepend(next);
    }

    //! ------------------
    //! Linear hexahedron
    //! ------------------
    HEXA8MeshData = new MeshVS_HArray1OfSequenceOfInteger(1,6);
    HEXA8MeshData->ChangeValue(1).Append(3); HEXA8MeshData->ChangeValue(1).Append(2); HEXA8MeshData->ChangeValue(1).Append(1); HEXA8MeshData->ChangeValue(1).Append(0);
    HEXA8MeshData->ChangeValue(2).Append(4); HEXA8MeshData->ChangeValue(2).Append(5); HEXA8MeshData->ChangeValue(2).Append(6); HEXA8MeshData->ChangeValue(2).Append(7);
    HEXA8MeshData->ChangeValue(3).Append(1); HEXA8MeshData->ChangeValue(3).Append(5); HEXA8MeshData->ChangeValue(3).Append(4); HEXA8MeshData->ChangeValue(3).Append(0);
    HEXA8MeshData->ChangeValue(4).Append(2); HEXA8MeshData->ChangeValue(4).Append(6); HEXA8MeshData->ChangeValue(4).Append(5); HEXA8MeshData->ChangeValue(4).Append(1);
    HEXA8MeshData->ChangeValue(5).Append(3); HEXA8MeshData->ChangeValue(5).Append(7); HEXA8MeshData->ChangeValue(5).Append(6); HEXA8MeshData->ChangeValue(5).Append(2);
    HEXA8MeshData->ChangeValue(6).Append(0); HEXA8MeshData->ChangeValue(6).Append(4); HEXA8MeshData->ChangeValue(6).Append(7); HEXA8MeshData->ChangeValue(6).Append(3);

    //! ----------------------
    //! Quadratic tetrahedron
    //! ----------------------
    const int mask[4][6] = {{0,4,1,5,2,6},{1,5,2,7,3,8},{0,6,2,7,3,9},{1,8,3,9,0,4}};
    TET10MeshData = new MeshVS_HArray1OfSequenceOfInteger(1,4);
    TColStd_SequenceOfInteger face[4];
    for(int f=0; f<4; f++)
    {
        for(int j=0; j<6; j++) face[f].Append(mask[f][j]);
        TET10MeshData->SetValue(f+1,face[f]);
    }
}

//! ----------------------------------
//! function: computeNormalAtElements
//! details:
//! ----------------------------------
#define USE_POLYGON
void Ng_MeshVS_DataSource3D::computeNormalAtElements()
{
    cout<<"Ng_MeshVS_DataSource3D::computeNormalAtElements()->____function called____"<<endl;


    TColStd_MapIteratorOfPackedMapOfInteger it;
    for(it.Initialize(myElements); it.More(); it.Next())
    {
        std::vector<polygon::Point> points;

        int globalElementID = it.Key();
        int localElementID = myElementsMap.FindIndex(globalElementID);

        cout<<"Ng_MeshVS_DataSource3D::computeNormalAtElements()->____(global element ID, local element ID) = ("<<globalElementID<<", "<<localElementID<<")____"<<endl;

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
        std::vector<double> n = polygon::getNormal(points);

        myElemNormals->SetValue(localElementID,1,n.at(0));
        myElemNormals->SetValue(localElementID,2,n.at(1));
        myElemNormals->SetValue(localElementID,3,n.at(2));
        cout<<"compute normal at elements____("<<n.at(0)<<", "<<n.at(1)<<", "<<n.at(2)<<")____"<<endl;
    }
}

//! ---------------------------
//! function: changeNodeCoords
//! details:
//! ---------------------------
bool Ng_MeshVS_DataSource3D::changeNodeCoords(int localNodeID, const std::vector<double> &coords)
{
    //! -------------
    //! sanity check
    //! -------------
    if(localNodeID<1 || localNodeID>myNumberOfNodes)
    {
        cerr<<"Ng_MeshVS_DataSource3D::changeNodeCoords()->____local node ID out of range when trying to change coordinates____"<<endl;
        return false;
    }
    if(coords.size()==0) return false;

    myNodeCoords->ChangeValue(localNodeID,1) = coords.at(0);
    myNodeCoords->ChangeValue(localNodeID,2) = coords.at(1);
    myNodeCoords->ChangeValue(localNodeID,3) = coords.at(2);
    return true;
}

//! ------------------------
//! function: removeElement
//! details:  experimental
//! ------------------------
void Ng_MeshVS_DataSource3D::removeElement(int localElementID)
{
    cout<<"Ng_MeshVS_DataSource3D::removeElement()->____function called: removing element ID: "<<
          localElementID<<"____"<<endl;

    //! ---------------------------------------------------------------
    //! copy the map of elements - jump over the element to be removed
    //! ---------------------------------------------------------------
    TColStd_PackedMapOfInteger myElements_new;
    for(TColStd_MapIteratorOfPackedMapOfInteger it(myElements); it.More(); it.Next())
    {
        if(it.Key()==localElementID) continue;
        myElements_new.Add(it.Key());
    }

    //! -----------------------
    //! new number of elements
    //! -----------------------
    int myNumberOfElements_new = myElements_new.Extent();

    occHandle(TColStd_HArray1OfInteger) myElemType_new = new TColStd_HArray1OfInteger(1,myNumberOfElements_new);
    occHandle(TColStd_HArray2OfInteger) myElemNodes_new = new TColStd_HArray2OfInteger(1,myNumberOfElements_new,1,20);
    occHandle(TColStd_HArray2OfReal) myElemNormals_new = new TColStd_HArray2OfReal(1,myNumberOfElements_new,1,3);

    //! -----------------------
    //! clear the map of nodes
    //! -----------------------
    myNodes.Clear();
    myNodesMap.Clear();

    //! -------------------------------------------
    //! redefine the elements and the map of nodes
    //! jump over the element to be removed
    //! -------------------------------------------
    int curElementID_new = 0;
    for(TColStd_MapIteratorOfPackedMapOfInteger it(myElements); it.More(); it.Next())
    {
        int curElementID = it.Key();
        if(curElementID==localElementID) continue;
        curElementID_new++;

        //! --------------------------------------------
        //! retrieve the element type from the old mesh
        //! --------------------------------------------
        ElemType type = static_cast<ElemType>(myElemType->Value(curElementID));

        int NbNodes;
        switch(type)
        {
        case TET:
        {
            myElemType_new->SetValue(curElementID_new,TET);
            NbNodes = 4;
        }
            break;
        case TET10:
        {
            NbNodes = 10;
            myElemType_new->SetValue(curElementID_new,TET10);
        }
            //! other elements ... to do ...
            break;
        }
        for(int k=1; k<=NbNodes; k++)
        {
            //! retrieve the nodes of the element from the old mesh
            int nodeID = myElemNodes->Value(curElementID,k);
            myElemNodes_new->SetValue(curElementID_new,k,nodeID);
            if(!myNodes.Contains(nodeID))
            {
                myNodes.Add(nodeID);
                myNodesMap.Add(nodeID);
            }
        }
    }

    //! ------------------------------
    //! redefine the maps of elements
    //! ------------------------------
    myElements.Clear();
    myElementsMap.Clear();

    for(TColStd_MapIteratorOfPackedMapOfInteger it(myElements_new); it.More(); it.Next())
    {
        myElements.Add(it.Key());
        myElementsMap.Add(it.Key());
    }

    myNumberOfElements = myNumberOfElements_new;

    myElemType->Delete();
    myElemNodes->Delete();
    myElemNormals->Delete();

    myElemType = new TColStd_HArray1OfInteger(1,myNumberOfElements);
    myElemNodes = new TColStd_HArray2OfInteger(1,myNumberOfElements,1,20);
    myElemNormals = new TColStd_HArray2OfReal(1,myNumberOfElements,1,3);

    //! --------------------
    //! new number of nodes
    //! --------------------
    myNumberOfNodes = myNodes.Extent();

    //! ---------------------------------------
    //! redefine the array of node coordinates
    //! ---------------------------------------
    occHandle(TColStd_HArray2OfReal) myNodeCoords_new = new TColStd_HArray2OfReal(1,myNumberOfNodes,1,3);
    int n=0;
    for(TColStd_MapIteratorOfPackedMapOfInteger it(myNodes); it.More(); it.Next())
    {
        n++;
        int nodeID = it.Key();
        double x = myNodeCoords->Value(nodeID,1);
        double y = myNodeCoords->Value(nodeID,2);
        double z = myNodeCoords->Value(nodeID,3);
        myNodeCoords_new->SetValue(myNodesMap.FindKey(n),1,x);
        myNodeCoords_new->SetValue(myNodesMap.FindKey(n) ,2,y);
        myNodeCoords_new->SetValue(myNodesMap.FindKey(n) ,3,z);
    }

    //! ----------------------
    //! redefine the elements
    //! ----------------------
    myElemNodes->Assign(*myElemNodes_new.get());

    //! ------------------------------
    //! redefine the node coordinates
    //! ------------------------------
    myNodeCoords->Delete();
    myNodeCoords = new TColStd_HArray2OfReal(1,myNumberOfNodes,1,3);
    myNodeCoords->Assign(*myNodeCoords_new.get());

    cout<<"Ng_MeshVS_DataSource3D::removeElement()->____exit function___"<<endl;
}


//! -----------------------------
//! function: getNodeCoordinates
//! details:  helper
//! -----------------------------
std::vector<double> Ng_MeshVS_DataSource3D::getNodeCoordinates(int localNodeID)
{
    std::vector<double> coords;
    if(localNodeID<1 || localNodeID>myNumberOfNodes)
    {
        cout<<"Ng_MeshVS_DataSource3D::getNodeCoordinates()->____localNodeID: "<<localNodeID<<" is out of range____"<<endl;
        return coords;
    }
    coords.push_back(myNodeCoords->Value(localNodeID,1));
    coords.push_back(myNodeCoords->Value(localNodeID,2));
    coords.push_back(myNodeCoords->Value(localNodeID,3));
    return coords;
}

//! ------------------------------------------------------------------------------------
//! function: buildCCXFaceToElementConnectivity
//! details:  map
//!           key => meshElement2D (a vector of global node IDs)
//!           value => vector of pairs {(elementID,CCX faceNr),(elementID, CCX faceNr)}
//! ------------------------------------------------------------------------------------
void Ng_MeshVS_DataSource3D::buildCCXFaceToElementConnectivity(std::map<meshElement2D,std::vector<std::pair<int,int>>> &map)
{
    //cout<<"Ng_MeshVS_DataSource3D::buildCCXFaceToElementConnectivity()->____generating face to element connectivity data____"<<endl;

    //! -----------------------------
    //! CCX faces of the TET element
    //! -----------------------------
    std::map<int,int> OCCfaceNrToCCXFaceNr_TET;
    OCCfaceNrToCCXFaceNr_TET[1] = 2;
    OCCfaceNrToCCXFaceNr_TET[2] = 3;
    OCCfaceNrToCCXFaceNr_TET[3] = 4;
    OCCfaceNrToCCXFaceNr_TET[4] = 1;

    //! -----------------------------
    //! CCX faces of the TET element
    //! -----------------------------
    std::map<int,int> OCCfaceNrToCCXFaceNr_TET10;
    OCCfaceNrToCCXFaceNr_TET10[1] = 2;
    OCCfaceNrToCCXFaceNr_TET10[2] = 3;
    OCCfaceNrToCCXFaceNr_TET10[3] = 4;
    OCCfaceNrToCCXFaceNr_TET10[4] = 1;

    //! ------------------------------
    //! CCX faces of the HEXA element
    //! ------------------------------
    std::map<int,int> OCCfaceNrToCCXFaceNr_HEXA;
    OCCfaceNrToCCXFaceNr_HEXA[1] = 1;
    OCCfaceNrToCCXFaceNr_HEXA[2] = 2;
    OCCfaceNrToCCXFaceNr_HEXA[3] = 5;
    OCCfaceNrToCCXFaceNr_HEXA[4] = 4;
    OCCfaceNrToCCXFaceNr_HEXA[5] = 3;
    OCCfaceNrToCCXFaceNr_HEXA[6] = 6;

    //! -------------------------------
    //! CCX faces of the PYRAM element
    //! -------------------------------
    std::map<int,int> OCCfaceNrToCCXFaceNr_PYRAM;
    OCCfaceNrToCCXFaceNr_PYRAM[1] = 1;
    OCCfaceNrToCCXFaceNr_PYRAM[2] = 4;
    OCCfaceNrToCCXFaceNr_PYRAM[3] = 6;
    OCCfaceNrToCCXFaceNr_PYRAM[4] = 5;
    OCCfaceNrToCCXFaceNr_PYRAM[5] = 3;

    //! -------------------------------
    //! CCX faces of the PRISM element
    //! -------------------------------
    std::map<int,int> OCCfaceNrToCCXFaceNr_PRISM;
    OCCfaceNrToCCXFaceNr_PRISM[1] = 1;
    OCCfaceNrToCCXFaceNr_PRISM[2] = 2;
    OCCfaceNrToCCXFaceNr_PRISM[3] = 3;
    OCCfaceNrToCCXFaceNr_PRISM[4] = 4;
    OCCfaceNrToCCXFaceNr_PRISM[5] = 5;

    //! -------------------------
    //! other elements ... to do
    //! -------------------------
    for(int localElementID=1; localElementID<=myNumberOfElements; localElementID++)
    {
        ElemType type = static_cast<ElemType>(myElemType->Value(localElementID));
        occHandle(MeshVS_HArray1OfSequenceOfInteger) topology;
        int NbFaces;
        switch(type)
        {
        case TET: topology = TET4MeshData; NbFaces = 4; break;
        case PYRAM: topology = PYRAM5MeshData; NbFaces = 5; break;
        case PRISM: topology = PRISM6MeshData; NbFaces = 5; break;
        case HEXA: topology = HEXA8MeshData; NbFaces = 6; break;

        case TET10: topology = TET10MeshData; NbFaces = 4; break;
        //case HEXA20: topology = HEXA20MeshData; NbFaces = 8; break;
        //case PRISM15: topology = PRISM15MeshData; NbFaces = 6; break;
        //case PYRAM13: topology = PYRAM13MeshData; NbFaces = 5; break;
        }

        //! --------------------------------------
        //! iterate over the faces of the element
        //! --------------------------------------
        for(int f=1; f<=NbFaces; f++)
        {
            TColStd_SequenceOfInteger indices = topology->Value(f);

            //! -------------------------
            //! create the meshElement2D
            //! -------------------------
            meshElement2D aMeshElement2D;
            int NbIndices = indices.Length();
            for(int n=1; n<=NbIndices; n++)
            {
                int index = indices.Value(n);
                aMeshElement2D.nodeIDs<<myElemNodes->Value(localElementID,index+1);
            }

            //! ----------------------------------------------
            //! check if the current meshElement2D is present
            //! within the map defined at the beginning
            //! ----------------------------------------------
            std::map<meshElement2D,std::vector<std::pair<int,int>>>::iterator it = map.find(aMeshElement2D);
            if(it!=map.end())
            {
                //! ---------------------------------------------------------------
                //! the current mesh element 2D has been already recorded one time
                //! ---------------------------------------------------------------
                std::pair<meshElement2D,std::vector<std::pair<int,int>>> aPair = *it;
                std::vector<std::pair<int,int>> CCXdef = aPair.second;

                //! ------------------------------------------
                //! calculate the CCX element face definition
                //! ------------------------------------------
                std::pair<int,int> apair;
                apair.first = localElementID;
                switch(type)
                {
                case TET: apair.second = OCCfaceNrToCCXFaceNr_TET.at(f); break;
                case TET10: apair.second = OCCfaceNrToCCXFaceNr_TET10.at(f); break;
                case HEXA:  apair.second = OCCfaceNrToCCXFaceNr_HEXA.at(f); break;
                case PYRAM: apair.second = OCCfaceNrToCCXFaceNr_PYRAM.at(f); break;
                case PRISM: apair.second = OCCfaceNrToCCXFaceNr_PRISM.at(f); break;
                }

                CCXdef.push_back(apair);

                //! -----------------
                //! update the value
                //! -----------------
                it->second = CCXdef;
            }
            else
            {
                //! -------------------------------------------------------
                //! the current mesh element has not been already recorded
                //! -------------------------------------------------------
                std::vector<std::pair<int,int>> CCXdef;

                //! ------------------------------------------
                //! calculate the CCX element face definition
                //! ------------------------------------------
                std::pair<int,int> apair;
                apair.first = localElementID;
                switch(type)
                {
                case TET: apair.second = OCCfaceNrToCCXFaceNr_TET.at(f); break;
                case TET10: apair.second = OCCfaceNrToCCXFaceNr_TET10.at(f); break;
                case HEXA: apair.second = OCCfaceNrToCCXFaceNr_HEXA.at(f); break;
                case PYRAM: apair.second = OCCfaceNrToCCXFaceNr_PYRAM.at(f); break;
                case PRISM: apair.second = OCCfaceNrToCCXFaceNr_PRISM.at(f); break;
                }

                CCXdef.push_back(apair);
                std::pair<meshElement2D,std::vector<std::pair<int,int>>> newElement;
                newElement.first = aMeshElement2D;
                newElement.second = CCXdef;

                //! --------------------
                //! insert into the map
                //! --------------------
                map.insert(newElement);
            }
        }
    }
    //cout<<"Ng_MeshVS_DataSource3D::buildCCXFaceToElementConnectivity()->____CCX face to element connectivity data generated____"<<endl;
}


//! -------------------------------------------------------------------
//! function: displaceMySelf
//! details:  displace the nodes using the "displacementField"
//!           the key of the input map should be the global element ID
//!           update normals at elements, at nodes, and 1D boundary
//! -------------------------------------------------------------------
void Ng_MeshVS_DataSource3D::displaceMySelf(const QMap<int, gp_Vec> &displacementField)
{
    for(QMap<int, gp_Vec>::const_iterator it = displacementField.cbegin(); it!= displacementField.cend(); ++it)
    {
        int globalNodeID = it.key();
        int localNodeID = this->myNodesMap.FindIndex(globalNodeID);
        if(localNodeID<0 || localNodeID>myNumberOfNodes) return;
        const gp_Vec &curVec = it.value();
        double dX = curVec.X();
        double dY = curVec.Y();
        double dZ = curVec.Z();
        myNodeCoords->ChangeValue(localNodeID,1) = myNodeCoords->Value(localNodeID,1)+dX;
        myNodeCoords->ChangeValue(localNodeID,2) = myNodeCoords->Value(localNodeID,2)+dY;
        myNodeCoords->ChangeValue(localNodeID,3) = myNodeCoords->Value(localNodeID,3)+dZ;
    }
}

//! -------------------------------------------
//! function: displaceMySelf_asRigidAsPossible
//! details:
//! -------------------------------------------
void Ng_MeshVS_DataSource3D::displaceMySelf_asRigidAsPossible(const QMap<int, gp_Vec> &displacementField,
                                                              const std::vector<int> &nodeGroup,
                                                              int mode)
{
    cout<<"Ng_MeshVS_DataSource3D::displaceMySelf_asRigidAsPossible()->____function called____"<<endl;

    //! -----------------------------------
    //! convert into Eigen format the mesh
    //! -----------------------------------
    Eigen::MatrixXd V(myNumberOfNodes,3);
    for(int row=1; row<=myNumberOfNodes; row++)
        for(int col=1; col<=3; col++)
            V(row-1,col-1) = myNodeCoords->Value(row,col);

    Eigen::MatrixXi T(myNumberOfElements,4);
    for(int i=1; i<=myNumberOfElements; i++)
    {
        int row = i-1;
        T(row,0) = myNodesMap.FindIndex(myElemNodes->Value(i,1))-1;
        T(row,1) = myNodesMap.FindIndex(myElemNodes->Value(i,2))-1;
        T(row,2) = myNodesMap.FindIndex(myElemNodes->Value(i,3))-1;
        T(row,3) = myNodesMap.FindIndex(myElemNodes->Value(i,4))-1;
    }

    int NforAlgo = int(displacementField.size()+nodeGroup.size());
    Eigen::VectorXi b(NforAlgo);
    Eigen::MatrixXd bc(NforAlgo,3);
    Eigen::MatrixXd bc_arap(NforAlgo,3);

    //! ----------------
    //! the fixed nodes
    //! ----------------
    cout<<"Ng_MeshVS_DataSource3D::displaceMySelf_asRigidAsPossible()->____setting up nodes for deformation____"<<endl;
    int i;
    for(i = 0; i<nodeGroup.size(); i++)
    {
        b(i) = myNodesMap.FindIndex(nodeGroup.at(i))-1;
        bc(i,0)= 0;
        bc(i,1)= 0;
        bc(i,2)= 0;
        int localNodeID = myNodesMap.FindIndex(nodeGroup.at(i));
        bc_arap(i,0)= myNodeCoords->Value(localNodeID,1);
        bc_arap(i,1)= myNodeCoords->Value(localNodeID,2);
        bc_arap(i,2)= myNodeCoords->Value(localNodeID,3);
    }

    //! --------------------
    //! the displaced nodes
    //! --------------------
    for(QMap<int,gp_Vec>::const_iterator it = displacementField.cbegin(); it != displacementField.cend(); it++)
    {
        int globalNodeID = it.key();
        int localNodeID = myNodesMap.FindIndex(globalNodeID);
        gp_Vec displ = it.value();

        b(i) = localNodeID-1;

        bc(i,0)= displ.X();
        bc(i,1)= displ.Y();
        bc(i,2)= displ.Z();

        bc_arap(i,0)= myNodeCoords->Value(localNodeID,1)+displ.X();
        bc_arap(i,1)= myNodeCoords->Value(localNodeID,2)+displ.Y();
        bc_arap(i,2)= myNodeCoords->Value(localNodeID,3)+displ.Z();

        i++;
    }
    cout<<"Ng_MeshVS_DataSource3D::displaceMySelf_asRigidAsPossible()->____setting up nodes for deformation: done____"<<endl;

    // -------------------------------------------------------------------
    // Compute k-harmonic weight functions "coordinates".
    //
    // Inputs:
    //   V  #V by dim vertex positions
    //   F  #F by simplex-size list of element indices
    //   b  #b boundary indices into V
    //   bc #b by #W list of boundary values
    //   k  power of harmonic operation (1: harmonic, 2: biharmonic, etc)
    // Outputs:
    //   W  #V by #W list of weights
    // -------------------------------------------------------------------
    cout<<"Ng_MeshVS_DataSource3D::displaceMySelf_asRigidAsPossible()->____harmonic displacement____"<<endl;
    Eigen::MatrixXd U(myNumberOfNodes,3);
    Eigen::MatrixXd D(myNumberOfNodes,3);              //! displacements
    bool isHarmonicDone = igl::harmonic(V,T,b,bc,1,D);
    if(isHarmonicDone) U = V+D;
    else U = V;
    cout<<"Ng_MeshVS_DataSource3D::displaceMySelf_asRigidAsPossible()->____harmonic displacement done____"<<endl;

    if(mode == 1)
    {
        cout<<"Ng_MeshVS_DataSource3D::displaceMySelf_asRigidAsPossible()->____ARAP deformation____"<<endl;
        //! --------------------------------------------------------
        //! https://libigl.github.io/tutorial/#as-rigid-as-possible
        //! "as rigid as possible"
        //! --------------------------------------------------------
        igl::ARAPData arap_data;
        arap_data.energy = igl::ARAP_ENERGY_TYPE_ELEMENTS;
        arap_data.with_dynamics = false;
        arap_data.max_iter = 100;
        arap_data.dim = 3;

        // -----------------------------------------------------------------------------
        // Precompute data needed to efficiently conduct local-global iterations for an
        // arap deformation. This includes the `data` struct employed by
        // `igl::min_quad_with_fixed` to solve the global step" and constructing the
        // bi-linear form `K` that mixes rotation degrees of freedom with unknown
        // positions for preparing the covariance matrices of the local step and the
        // linear term of the global step.
        //
        // Inputs:
        //   V  #V by dim vertex positions
        //   F  #F by simplex-size list of element indices
        //   b  #b indices into V of handle vertices
        // Outputs:
        //   data  pre-factorized system matrix etc. (see `igl::min_quad_with_fixed`)
        //   K  #R*dim by #V
        // -----------------------------------------------------------------------------
        bool isDone = igl::arap_precomputation(V,T,3,b,arap_data);
        if(!isDone)
        {
            cout<<"Ng_MeshVS_DataSource3D::displaceMySelf_asRigidAsPossible()->____pre-computation not done____"<<endl;
            return;
        }
        cout<<"Ng_MeshVS_DataSource3D::displaceMySelf_asRigidAsPossible()->____ARAP pre-computation done____"<<endl;

        // ------------------------------------------------------------------
        // Inputs:
        //   bc  #b by dim list of boundary conditions
        //   data  struct containing necessary precomputation and parameters
        //   U  #V by dim initial guess
        // ------------------------------------------------------------------
        cout<<"Ng_MeshVS_DataSource3D::displaceMySelf_asRigidAsPossible()->____ARAP solving____"<<endl;
        igl::arap_solve(bc_arap,arap_data,U);
        cout<<"Ng_MeshVS_DataSource3D::displaceMySelf_asRigidAsPossible()->____ARAP solving done____"<<endl;
    }

    //! ----------------------------
    //! update the mesh coordinates
    //! ----------------------------
    cout<<"Ng_MeshVS_DataSource3D::displaceMySelf_asRigidAsPossible()->____updatig mesh points____"<<endl;
    int N = U.rows();
    for(int row=0; row<N; row++)
    {
        myNodeCoords->SetValue(row+1,1,U(row,0));
        myNodeCoords->SetValue(row+1,2,U(row,1));
        myNodeCoords->SetValue(row+1,3,U(row,2));
    }
    cout<<"Ng_MeshVS_DataSource3D::displaceMySelf_asRigidAsPossible()->____updatig mesh points done____"<<endl;
}

//! ---------------------------------------------------------------------------------
//! function: computeElementToElementConnectivity
//! details:  element global ID => {vector of attached elements (global element ID)}
//! ---------------------------------------------------------------------------------
void Ng_MeshVS_DataSource3D::buildElementToElementConnectivity(std::map<int,std::vector<int>> &elementToAttachedElements)
{
    //! ------------------------------------------------------------
    //! iterate over the elements - for each face of the element
    //! retrieve the attached elements (this needs that the face to
    //! element connectivity has been generated in advance)
    //! ------------------------------------------------------------
    if(myFaceToElements.isEmpty()) this->buildFaceToElementConnectivity();

    for(int localElementID = 1; localElementID <= myNumberOfElements; localElementID++)
    {
        int globalElementID = myElementsMap.FindKey(localElementID);

        //cout<<"@ -------------------------------"<<endl;
        //cout<<"@ working on element ("<<localElementID<<", "<<globalElementID<<")"<<endl;
        //cout<<"@ -------------------------------"<<endl;

        //! -------------------------------------------
        //! retrieve and scan the faces of the element
        //! -------------------------------------------
        occHandle(MeshVS_HArray1OfSequenceOfInteger) topology;
        switch(myElemType->Value(localElementID))
        {
        case TET: topology = TET4MeshData; break;
        case HEXA: topology = HEXA8MeshData; break;
        case PRISM: topology = PRISM6MeshData; break;
        case PYRAM: topology = PYRAM5MeshData; break;
        }
        int NbFaces = topology->Length();
        for(int n=1; n<=NbFaces; n++)
        {
            //cout<<"@ face number: "<<n<<endl;

            //! ------------------------------------------------------------
            //! define a face - note: the meshElement2D type is unessential
            //! ------------------------------------------------------------
            meshElement2D aMeshElement2D;
            TColStd_SequenceOfInteger indices = topology->Value(n);
            for(int j=1; j<=indices.Length(); j++)
            {
                int index = indices.Value(j);
                int globalNodeID = myElemNodes->Value(localElementID,index+1);
                aMeshElement2D.nodeIDs<<globalNodeID;
            }

            switch(indices.Length())
            {
            case 3: aMeshElement2D.type = TRIG; break;
            case 4: aMeshElement2D.type = QUAD; break;
            }

            //! ----------------------------------------------------------
            //! check the elements connected to the current meshElement2D
            //! (two elements at maximum)
            //! ----------------------------------------------------------
            QList<int> attachedElements = myFaceToElements.value(aMeshElement2D);
            //cout<<"@ number of attached elements: "<<attachedElements.size()<<endl;

            for(int i=0; i<attachedElements.size(); i++)
            {
                //! -------------------------------------------------------------
                //! use find key since the face to element connecticity has been
                //! defined using local element IDs
                //! -------------------------------------------------------------
                int curElementGlobalID = myElementsMap.FindKey(attachedElements.at(i));
                //cout<<"@ attached element: "<<curElementGlobalID;

                if(curElementGlobalID==globalElementID)
                {
                    //cout<<" - jumping over"<<endl;
                    continue;
                }
                //cout<<" - evaluating: ";

                std::map<int,std::vector<int>>::iterator it = elementToAttachedElements.find(globalElementID);
                if(it!=elementToAttachedElements.end())
                {
                    //! --------------------------------------------------------
                    //! an entry for "globalElementID" has already been created
                    //! --------------------------------------------------------
                    std::pair<int,std::vector<int>> aPair = *it;
                    std::vector<int> elementsAttachedToElement_ = aPair.second;
                    //if(std::find(elementsAttachedToElement_.begin(),elementsAttachedToElement_.begin(),curElementGlobalID)==elementsAttachedToElement_.end())
                    //{
                        //exit(8888);
                        //cout<<" not inserted yet into the list of attached"<<endl;
                        elementsAttachedToElement_.push_back(curElementGlobalID);
                        it->second = elementsAttachedToElement_;
                    //}
                }
                else
                {
                    //! --------------------------------------------------------
                    //! an entry for "globalElementID" has not been created yet
                    //! --------------------------------------------------------
                    //cout<<" not inserted yet into the list of attached"<<endl;
                    std::vector<int> elementsAttachedToElement_;
                    elementsAttachedToElement_.push_back(curElementGlobalID);
                    std::pair<int,std::vector<int>> aPair;
                    aPair.first = globalElementID;
                    aPair.second = elementsAttachedToElement_;
                    elementToAttachedElements.insert(aPair);
                }
            }
        }
    }
}

//! ----------------------------------------------------
//! function: getFaceElementsOfElement
//! details:  get the face elements of a volume element
//! ----------------------------------------------------
std::vector<meshElement2D> Ng_MeshVS_DataSource3D::getFaceElementsOfElement(int localElementID)
{
    //! ---------------------------------------------
    //! create the faces of the element remove faces
    //! ---------------------------------------------
    std::vector<meshElement2D> faceElementsOfTheElement;

    occHandle(MeshVS_HArray1OfSequenceOfInteger) topology;
    switch(myElemType->Value(localElementID))
    {
    case TET: topology = TET4MeshData; break;
    case HEXA: topology = HEXA8MeshData; break;
    }

    for(int f=1; f<=topology->Length(); f++)
    {
        TColStd_SequenceOfInteger indices = topology->Value(f);

        //! -------------------------
        //! create the meshElement2D
        //! -------------------------
        meshElement2D aMeshElement2D;
        int NbIndices = indices.Length();
        for(int n=1; n<=NbIndices; n++)
        {
            int index = indices.Value(n);
            aMeshElement2D.nodeIDs<<myElemNodes->Value(localElementID,index+1);

            switch(NbIndices)
            {
            case 3: aMeshElement2D.type = TRIG; break;
            case 4: aMeshElement2D.type = QUAD; break;
            case 6: aMeshElement2D.type = TRIG6; break;
            case 8: aMeshElement2D.type = QUAD8; break;
            }
        }
        faceElementsOfTheElement.push_back(aMeshElement2D);
    }
    return faceElementsOfTheElement;
}

//! ------------------------
//! function: renumberNodes
//! details:
//! ------------------------
void Ng_MeshVS_DataSource3D::renumberNodes(const std::map<std::size_t,int> &mapCoordNodeNumbers,
                                           const std::vector<mesh::tolerantPoint> &vecTolPoints,
                                           const std::vector<int> &GlobalNodeIDs)
{
    cout<<"____renumbering nodes for volume mesh____"<<endl;
    int retrievedByProximity = 0;
    int retrievedByHashNumber = 0;

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
            retrievedByHashNumber++;
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
                cout<<"//--------------------------------------------------------------//"<<endl;
                cout<<"// Ng_MeshVS_DataSourceFace::renumberNodes(): this part of code //"<<endl;
                cout<<"// shoud never be reached. The hash number of the node has not  //"<<endl;
                cout<<"// been found                                                   //"<<endl;
                cout<<"// Cannot found ("<<x<<", "<<y<<", "<<z<<") hash: "<<seed<<endl;
                cout<<"//--------------------------------------------------------------//"<<endl;
                exit(100);
            }
        }
    }

    cout<<"____retrieved by proximity: "<<double(retrievedByProximity)/double(myNumberOfNodes)*100.0<<" %____"<<endl;
    cout<<"____retrieved by hash number: "<<double(retrievedByHashNumber)/double(myNumberOfNodes)*100.0<<" %____"<<endl;
    if((retrievedByHashNumber+retrievedByProximity)==myNumberOfNodes)
        cout<<"____all the nodes have been renumbered____"<<endl;

    //! --------------------------------------
    //! update the definition of the elements
    //! --------------------------------------
    //cout<<"____update the definition of the elements____"<<endl;
    for(int localElementID = 1; localElementID<=myNumberOfElements; localElementID++)
    {
        int NbNodes;
        switch(myElemType->Value(localElementID))
        {
        case TET: NbNodes = 4; break;
        case HEXA: NbNodes = 8; break;
        case PYRAM: NbNodes = 5; break;
        case PRISM: NbNodes = 6; break;
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
            cout<<"//--------------------------------------------------------------//"<<endl;
            cout<<"// Ng_MeshVS_DataSource3D::renumberNodes(): this part of code   //"<<endl;
            cout<<"// shoud never be reached                                       //"<<endl;
            cout<<"//--------------------------------------------------------------//"<<endl;
            exit(102);
        }
    }
}

//! ------------------------------------------------------------------
//! function: buildTolerantPoints
//! details:  return the vector of all the mesh points with tolerance
//! ------------------------------------------------------------------
std::vector<mesh::tolerantPoint> Ng_MeshVS_DataSource3D::buildTolerantPoints(double tolerance)
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
std::vector<mesh::meshSegment> Ng_MeshVS_DataSource3D::getSegmentsOfElement(int elementID, bool isLocal)
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

    //! --------------------------------------------------
    //! iterate over the faces of the element
    //! the following code holds for first order elements
    //! (the "segments" are built using two nodes)
    //! --------------------------------------------------
    const std::vector<meshElement2D> &faceElements = getFaceElementsOfElement(localElementID);
    std::vector<mesh::meshSegment> elementSegments;
    for (int fn = 0; fn<faceElements.size(); fn++)
    {
        const meshElement2D &aFaceElement = faceElements[fn];
        int NbNodes = int(aFaceElement.nodeIDs.size());
        for(int i=0; i<NbNodes; i++)
        {
            int f = i%NbNodes;
            int s = (i+1)%NbNodes;
            mesh::meshSegment aSegment;
            unsigned int n1 = aFaceElement.nodeIDs[f];
            unsigned int n2 = aFaceElement.nodeIDs[s];

            if(n1<n2) aSegment.nodeIDs<<n1<<n2;
            else aSegment.nodeIDs<<n2<<n1;

            if(std::find(elementSegments.begin(),elementSegments.end(),aSegment)==elementSegments.end())
                elementSegments.push_back(aSegment);
        }
    }
    return elementSegments;
}

//! --------------------------------------------------------------
//! function: constructor
//! details:  convert a first order mesh into n-th order elements
//! --------------------------------------------------------------
Ng_MeshVS_DataSource3D::Ng_MeshVS_DataSource3D(const occHandle(Ng_MeshVS_DataSource3D) &aMesh, int order)
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
    myElemNodes = new TColStd_HArray2OfInteger(1,myNumberOfElements,1,20);

    //! -----------------------------
    //! copy the element definitions
    //! -----------------------------
    for(TColStd_MapIteratorOfPackedMapOfInteger eit(aMesh->GetAllElements()); eit.More(); eit.Next())
    {
        int globalElementID = eit.Key();
        int localElementID = aMesh->myElementsMap.FindIndex(globalElementID);
        int NbNodes, buf[20];
        TColStd_Array1OfInteger nodeIDs(*buf,1,20);
        aMesh->GetNodesByElement(globalElementID,nodeIDs,NbNodes);
        for(int i=1; i<=NbNodes; i++) myElemNodes->SetValue(localElementID,i,nodeIDs(i));
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
                std::pair<mesh::meshSegment,mesh::meshPoint> element;
                element.first = curSegment;
                element.second = mp;
                mapSegmentPointBetween.insert(element);

                //! ----------------------
                //! augment the node maps
                //! ----------------------
                myNodes.Add(globalNodeID_newNodes);
                myNodesMap.Add(globalNodeID_newNodes);
            }
        }
    }

    //! -------------------------------------
    //! number of nodes and node coordinates
    //! -------------------------------------
    myNumberOfNodes = aMesh->GetAllNodes().Extent() + additionalNodes;
    myNodeCoords = new TColStd_HArray2OfReal(1,myNumberOfNodes,1,3);

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
        int offset = 0;
        ElemType eType;
        aMesh->GetElementType(eType,globalElementID,false);
        switch(eType)
        {
        case TET:
        {
            myElemType->SetValue(localElementID,TET10);
            offset = 4;
        }
            break;
        case HEXA:
        {
            myElemType->SetValue(localElementID,HEXA20);
            offset = 8;
        }
            break;
        case PYRAM:
        {
            myElemType->SetValue(localElementID,PYRAM13);
            offset = 5;
        }
            break;
        case PRISM:
        {
            myElemType->SetValue(localElementID,PRISM15);
            offset = 6;
        }
            break;
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
            std::pair<mesh::meshSegment,mesh::meshPoint> element = *si;
            mesh::meshPoint aP = element.second;
            int globalNodeID = aP.ID;
            int localNodeID = myNodesMap.FindIndex(globalNodeID);
            myNodeCoords->SetValue(localNodeID,1,aP.x);
            myNodeCoords->SetValue(localNodeID,2,aP.y);
            myNodeCoords->SetValue(localNodeID,3,aP.z);
            myElemNodes->SetValue(localElementID,offset+j+1,globalNodeID);
        }
    }

    //! --------------------
    //! create the topology
    //! --------------------
    this->buildElementsTopology();
}

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
Ng_MeshVS_DataSource3D::Ng_MeshVS_DataSource3D(const occHandle(Ng_MeshVS_DataSource3D) &aMesh, const QMap<int,gp_Vec> &displacements)
{
    myNumberOfElements = aMesh->GetAllElements().Extent();
    myNumberOfNodes = aMesh->GetAllNodes().Extent();

    //! ---------------------------
    //! maps of nodes and elements
    //! ---------------------------
    myElements = aMesh->GetAllElements();
    myElementsMap = aMesh->myNodesMap;
    myNodes = aMesh->GetAllNodes();
    myNodesMap = aMesh->myNodesMap;

    myElemType = new TColStd_HArray1OfInteger(1,myNumberOfElements);
    myNodeCoords = new TColStd_HArray2OfReal(1,myNumberOfNodes,1,3);
    myElemNodes = new TColStd_HArray2OfInteger(1,myNumberOfElements,1,20);

    //! ----------------------------------
    //! elements definition through nodes
    //! ----------------------------------
    int localElementID = 0;
    int NbNodes, b[20];
    TColStd_Array1OfInteger nodeIDs(*b,1,20);
    for(TColStd_MapIteratorOfPackedMapOfInteger it(aMesh->GetAllElements()); it.More(); it.Next())
    {
        localElementID++;
        int globalElementID = it.Key();
        aMesh->GetNodesByElement(globalElementID,nodeIDs,NbNodes);
        for(int c=1; c<=NbNodes; c++) myElemNodes->SetValue(localElementID,c,nodeIDs(c));

        ElemType eType;
        aMesh->GetElementType(eType,globalElementID,false);
        myElemType->SetValue(localElementID,eType);
    }

    //! -----------------
    //! node coordinates
    //! -----------------
    double a[3];
    TColStd_Array1OfReal coords(*a,1,3);
    MeshVS_EntityType aType;
    int localNodeID = 0;
    for(TColStd_MapIteratorOfPackedMapOfInteger it(aMesh->GetAllNodes()); it.More(); it.Next())
    {
        localNodeID ++;
        int globalNodeID = it.Key();
        aMesh->GetGeom(globalNodeID,false,coords,NbNodes,aType);
        const gp_Vec &d = displacements.value(globalNodeID);
        for(int c=0; c<NbNodes; c++)
        {
            int s = 3*c;
            double x = coords(s+1) + d.X();
            double y = coords(s+2) + d.Y();
            double z = coords(s+3) + d.Z();
            myNodeCoords->SetValue(localNodeID,s+1,x);
            myNodeCoords->SetValue(localNodeID,s+2,y);
            myNodeCoords->SetValue(localNodeID,s+3,z);
        }
    }

    //! ------------------
    //! elements topology
    //! ------------------
    this->buildElementsTopology();
}

//! ----------------------
//! function: constructor
//! details:  upcast
//! ----------------------
Ng_MeshVS_DataSource3D::Ng_MeshVS_DataSource3D(const occHandle(MeshVS_DataSource) &aMeshDS)
{
    if(aMeshDS.IsNull()) return;
    if(aMeshDS->GetAllElements().Extent()<1) return;
    if(aMeshDS->GetAllNodes().Extent()<4) return;

    myNodes = aMeshDS->GetAllNodes();
    myNumberOfNodes = myNodes.Extent();
    for(TColStd_MapIteratorOfPackedMapOfInteger it(myNodes); it.More(); it.Next()) myNodesMap.Add(it.Key());

    myElements = aMeshDS->GetAllElements();
    myNumberOfElements = myElements.Extent();
    for(TColStd_MapIteratorOfPackedMapOfInteger it(myElements); it.More(); it.Next()) myElementsMap.Add(it.Key());

    myElemType = new TColStd_HArray1OfInteger(1,myNumberOfElements);
    myNodeCoords = new TColStd_HArray2OfReal(1,myNumberOfNodes,1,3);
    myElemNodes = new TColStd_HArray2OfInteger(1,myNumberOfElements,1,20);

    int NbNodes, buf[20];
    double bufd[3];
    TColStd_Array1OfInteger nodeIDs(*buf,1,20);
    TColStd_Array1OfReal coords(*bufd,1,3);
    MeshVS_EntityType aType;

    int localNodeID = 1;
    for(TColStd_MapIteratorOfPackedMapOfInteger it(myNodes); it.More(); it.Next(), localNodeID++)
    {
        int globalNodeID = it.Key();
        aMeshDS->GetGeom(globalNodeID,false,coords,NbNodes,aType);
        for(int i=1; i<=3; i++) myNodeCoords->SetValue(localNodeID,i,coords(i));
    }

    int localElementID = 1;
    for(TColStd_MapIteratorOfPackedMapOfInteger it(myElements); it.More(); it.Next(),localElementID++)
    {
        int globalElementID = it.Key();
        aMeshDS->GetNodesByElement(globalElementID,nodeIDs,NbNodes);
        for(int i=1; i<=NbNodes; i++) myElemNodes->SetValue(globalElementID,i,nodeIDs(i));
        switch(NbNodes)
        {
        case 4: myElemType->SetValue(localElementID,TET); break;
        case 5: myElemType->SetValue(localElementID,PYRAM); break;
        case 6: myElemType->SetValue(localElementID,PRISM); break;
        case 8: myElemType->SetValue(localElementID,HEXA); break;
        case 10: myElemType->SetValue(localElementID,TET10); break;
        case 13: myElemType->SetValue(localElementID,PYRAM13); break;
        case 15: myElemType->SetValue(localElementID,PRISM15); break;
        case 20: myElemType->SetValue(localElementID,HEXA20); break;
        }
    }

    this->buildElementsTopology();
}
