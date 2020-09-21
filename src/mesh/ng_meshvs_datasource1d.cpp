//! custom includes
#include "ng_meshvs_datasource1d.h"

//! OCC
#include <TColStd_MapIteratorOfPackedMapOfInteger.hxx>
#include <TColStd_HArray2OfReal.hxx>
#include <TColStd_HArray1OfInteger.hxx>
#include <TColStd_IndexedMapOfInteger.hxx>

IMPLEMENT_STANDARD_HANDLE(Ng_MeshVS_DataSource1D,MeshVS_DataSourceExtended)
IMPLEMENT_STANDARD_RTTIEXT(Ng_MeshVS_DataSource1D,MeshVS_DataSourceExtended)

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
Ng_MeshVS_DataSource1D::Ng_MeshVS_DataSource1D(const TColStd_PackedMapOfInteger &theNodes,
                                               const occHandle(TColStd_HArray2OfReal) &theCoords,
                                               Standard_Integer theMeshOrder)
{
#ifdef VERBOSE
    cout<<"Ng_MeshVS_DataSource1D::Ng_MeshVS_DataSource1D()->____CONSTRUCTOR CALLED____"<<endl;
    for(int n=1; n<=theCoords->UpperRow();n++)
        cout<<"n = "<<n<<" x ="<<theCoords->Value(n,1)<<" y= "<<theCoords->Value(n,2)<<" z= "<<theCoords->Value(n,3)<<endl;
#endif

    if(theNodes.Extent()>1)
    {
    //! the mesh order
    myMeshOrder = theMeshOrder;

    //! add the nodes (map of nodes)
    myNodes = theNodes;

    //! the number of nodes
    myNumberOfNodes = theNodes.Extent();

    //! the number of elements
    if(myMeshOrder == 1) myNumberOfEdgeElements = theNodes.Extent()-1;
    else myNumberOfEdgeElements= (theNodes.Extent()-1)/2;

    //! my nodes coordinates
    myNodeCoords = new TColStd_HArray2OfReal(1,myNumberOfNodes,1,3);

    //! fill the node map and add the coordinates
    TColStd_MapIteratorOfPackedMapOfInteger aNodeIter;

    for(aNodeIter.Initialize(myNodes);aNodeIter.More();aNodeIter.Next()) myNodesMap.Add(aNodeIter.Key());
    for(int i=1; i<=myNumberOfNodes; i++)
    {
        Standard_Real x = theCoords->Value(i, 1);
        Standard_Real y = theCoords->Value(i, 2);
        Standard_Real z = theCoords->Value(i, 3);
        myNodeCoords->SetValue(i,1,x);
        myNodeCoords->SetValue(i,2,y);
        myNodeCoords->SetValue(i,3,z);

#ifdef VERBOSE
        int aKey = myNodesMap.FindKey(i);
        cout<<"Ng_MeshVS_DataSource1D::Ng_MeshVS_DataSource1D->____Node index: "<<i<<" ID: "<<aKey<<" x= "<<theCoords->Value(i,1)
           <<" y = "<<theCoords->Value(i,2)<<" z = "<<theCoords->Value(1,3)<<"____"<<endl;
#endif
    }

    std::vector<int> vecID;
    for(aNodeIter.Initialize(myNodes);aNodeIter.More();aNodeIter.Next())
    {
        Standard_Integer aKey = aNodeIter.Key();
        vecID.push_back(aKey);
    }

    if(myMeshOrder==1) myElemNodes = new TColStd_HArray2OfInteger(1,myNumberOfEdgeElements,1,2);
    else myElemNodes = new TColStd_HArray2OfInteger(1,(theNodes.Extent()-1)/2,1,3);

    if(myMeshOrder==1)
    {
        //! first order mesh
        for(int i=0; i<myNumberOfEdgeElements;i++)
        {
#ifdef VERBOSE
            cout<<"("<<vecID.at(i)<<", "<<vecID.at(i+1)<<")"<<endl;
#endif
            myElements.Add(i+1);
            myElementsMap.Add(i+1);
            myElemNodes->SetValue(i+1,1,vecID.at(i));
            myElemNodes->SetValue(i+1,2,vecID.at(i+1));
        }
    }
    else
    {
        //! second order mesh
        for(int i=1; i<=(theNodes.Extent()-1)/2; i++)
        {
            myElements.Add(i);
            myElementsMap.Add(i);
            int l =i-1;
            int m =(myNumberOfNodes-1)/2+i-1;
            int r =i-1+2;

            myElemNodes->SetValue(i,1,vecID.at(l));
            myElemNodes->SetValue(i,2,vecID.at(r));
            myElemNodes->SetValue(i,3,vecID.at(m));
            //cout<<"(l, r, m) = ("<<vecID.at(l)<<", "<<vecID.at(r)<<", "<<vecID.at(m)<<")"<<endl;
        }
    }
#ifdef VERBOSE
    cout<<"Ng_MeshVS_DataSource1D::Ng_MeshVS_DataSource1D->____Number of 1D elements for the edge: "<<myNumberOfEdgeElements<<"____"<<endl;
#endif
    }
}

//! ------------------
//! function: GetGeom
//! details:
//! ------------------
Standard_Boolean Ng_MeshVS_DataSource1D::GetGeom (const Standard_Integer ID,
                          const Standard_Boolean IsElement,
                          TColStd_Array1OfReal& Coords,
                          Standard_Integer& NbNodes,
                          MeshVS_EntityType& Type) const
{
    Standard_Boolean RV = Standard_False;
    if(IsElement==Standard_True)
    {
        if(ID>=1 && ID<=myNumberOfEdgeElements)
        {
            Type = MeshVS_ET_Link;
            myMeshOrder == 1 ? NbNodes = 2: 3;
            for(Standard_Integer i=1,k=1;i<=NbNodes;i++)
            {
                //! i spans the nodes of the element
                //! j spans {x, y, z}
                for(Standard_Integer j=1;j<=3;j++)
                {
                    Coords(k)=myNodeCoords->Value(ID,j);
                    k++;
                }
            }
            RV = Standard_True;
        }
        else RV = Standard_False;
    }
    else if (IsElement==Standard_False)
    {
        if(ID>=myNodes.GetMinimalMapped() && ID<=myNodes.GetMaximalMapped())
        {
            int nodeNr = myNodesMap.FindIndex(ID);
            Type = MeshVS_ET_Node;
            NbNodes = 1;
            Coords(1)=myNodeCoords->Value(nodeNr,1);
            Coords(2)=myNodeCoords->Value(nodeNr,2);
            Coords(3)=myNodeCoords->Value(nodeNr,3);
            RV = Standard_True;
        }
        else RV = Standard_False;
    }
    return RV;
}

//! ---------------------
//! function: GetllNodes
//! details:
//! ---------------------
const TColStd_PackedMapOfInteger& Ng_MeshVS_DataSource1D::GetAllNodes() const
{
    return myNodes;
}

//! -------------------------
//! function: GetAllElements
//! details:
//! -------------------------
const TColStd_PackedMapOfInteger& Ng_MeshVS_DataSource1D::GetAllElements() const
{
    return myElements;
}

//! ----------------------
//! function: GetGeomType
//! details:
//! ----------------------
Standard_Boolean Ng_MeshVS_DataSource1D::GetGeomType(const Standard_Integer ID,
                                             const Standard_Boolean IsElement,
                                             MeshVS_EntityType& Type)
                                             const
{
    //!cout<<"Ng_MeshVS_DataSource1D::GetGeomType()->____called____"<<endl;
    Standard_Boolean RV = Standard_False;
    if(IsElement == Standard_False)
    {
        if(ID>=1 && ID<=myNumberOfNodes)
        {
            RV = Standard_True;
            Type = MeshVS_ET_Node;
        }
    }
    else
    {
        if(ID>=1 && ID<=myNumberOfEdgeElements)
        {
            RV = Standard_True;
            Type = MeshVS_ET_Link;
        }
    }
    return RV;
}

//! ------------------
//! function: GetAddr
//! details:
//! ------------------
Standard_Address Ng_MeshVS_DataSource1D::GetAddr(const Standard_Integer /*ID*/,
                                                 const Standard_Boolean /*IsElement*/) const
{
    return NULL;
}

//! ----------------------------
//! function: GetNodesByElement
//! details:
//! ----------------------------
Standard_Boolean Ng_MeshVS_DataSource1D::GetNodesByElement(const Standard_Integer ID,
                                                          TColStd_Array1OfInteger& theNodeIDs,
                                                          Standard_Integer& theNbNodes) const
{
    if(ID>=1 && ID<=myNumberOfEdgeElements && theNodeIDs.Length()>=2)
    {
        Standard_Integer aLow = theNodeIDs.Lower();
        switch(theNodeIDs.Length())
        {
        case 2:
            theNbNodes=2;
            theNodeIDs(aLow+0)=myElemNodes->Value(ID,1);
            theNodeIDs(aLow+1)=myElemNodes->Value(ID,2);
            break;
        case 3:
            theNbNodes=3;
            theNodeIDs(aLow+0)=myElemNodes->Value(ID,1);
            theNodeIDs(aLow+1)=myElemNodes->Value(ID,2);
            theNodeIDs(aLow+2)=myElemNodes->Value(ID,3);
            break;
        }
        return Standard_True;
    }
    return Standard_False;
}

occHandle(TColStd_HArray2OfReal) Ng_MeshVS_DataSource1D::getNodesCoords() const
{
    return myNodeCoords;
}
