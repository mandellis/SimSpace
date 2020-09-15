//! ----------------
//! custom includes
//! ----------------
#include "ng_meshvs_datasource1d.h"

//! ----
//! OCC
//! ----
#include <TColStd_MapIteratorOfPackedMapOfInteger.hxx>
#include <TColStd_HArray2OfReal.hxx>
#include <TColStd_HArray1OfInteger.hxx>
#include <TColStd_IndexedMapOfInteger.hxx>

//! ----
//! C++
//! ----
#include <vector>

IMPLEMENT_STANDARD_HANDLE(Ng_MeshVS_DataSource1D,MeshVS_DataSourceExtended)
IMPLEMENT_STANDARD_RTTIEXT(Ng_MeshVS_DataSource1D,MeshVS_DataSourceExtended)

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
Ng_MeshVS_DataSource1D::Ng_MeshVS_DataSource1D(const TColStd_PackedMapOfInteger &theNodes,
                                               const occHandle(TColStd_HArray2OfReal) &theCoords,
                                               int theMeshOrder)
{
    if(theNodes.Extent()<2) return;

    myMeshOrder = theMeshOrder;
    myNodes = theNodes;
    myNumberOfNodes = theNodes.Extent();

    if(myMeshOrder == 1) myNumberOfElements = theNodes.Extent()-1;
    else myNumberOfElements= (theNodes.Extent()-1)/2;

    myNodeCoords = new TColStd_HArray2OfReal(1,myNumberOfNodes,1,3);

    for(TColStd_MapIteratorOfPackedMapOfInteger aNodeIter(myNodes);aNodeIter.More();aNodeIter.Next()) myNodesMap.Add(aNodeIter.Key());
    for(int i=1; i<=myNumberOfNodes; i++)
    {
        Standard_Real x = theCoords->Value(i, 1);
        Standard_Real y = theCoords->Value(i, 2);
        Standard_Real z = theCoords->Value(i, 3);
        myNodeCoords->SetValue(i,1,x);
        myNodeCoords->SetValue(i,2,y);
        myNodeCoords->SetValue(i,3,z);
    }

    std::vector<int> vecID;
    for(TColStd_MapIteratorOfPackedMapOfInteger aNodeIter(myNodes);aNodeIter.More();aNodeIter.Next())
    {
        Standard_Integer aKey = aNodeIter.Key();
        vecID.push_back(aKey);
    }

    if(myMeshOrder==1) myElemNodes = new TColStd_HArray2OfInteger(1,myNumberOfElements,1,2);
    else myElemNodes = new TColStd_HArray2OfInteger(1,(theNodes.Extent()-1)/2,1,3);

    if(myMeshOrder==1)
    {
        //! first order mesh
        for(int i=1; i<=myNumberOfElements;i++)
        {
            myElements.Add(i);
            myElementsMap.Add(i);
            myElemNodes->SetValue(i,1,vecID.at(i-1));
            myElemNodes->SetValue(i,2,vecID.at(i));
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
        }
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
    if(IsElement==true)
    {
        if(ID>=1 && ID<=myNumberOfElements)
        {
            Type = MeshVS_ET_Link;
            myMeshOrder == 1 ? NbNodes = 2: 3;
            for(Standard_Integer i=1,k=1;i<=NbNodes;i++)
            {
                for(Standard_Integer j=1;j<=3;j++)
                {
                    Coords(k)=myNodeCoords->Value(ID,j);
                    k++;
                }
            }
            return true;
        }
        else return false;
    }
    else
    {
        if(ID>=myNodes.GetMinimalMapped() && ID<=myNodes.GetMaximalMapped())
        {
            int nodeNr = myNodesMap.FindIndex(ID);
            Type = MeshVS_ET_Node;
            NbNodes = 1;
            Coords(1)=myNodeCoords->Value(nodeNr,1);
            Coords(2)=myNodeCoords->Value(nodeNr,2);
            Coords(3)=myNodeCoords->Value(nodeNr,3);
            return true;
        }
        else return false;
    }
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
    if(IsElement == false)
    {
        if(ID>=1 && ID<=myNumberOfNodes)
        {
            Type = MeshVS_ET_Node;
            return true;
        }
        return false;
    }
    if(ID>=1 && ID<=myNumberOfElements)
    {
        Type = MeshVS_ET_Link;
        return true;
    }
    return false;
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
    if(ID>=1 && ID<=myNumberOfElements && theNodeIDs.Length()>=2)
    {
        int aLow = theNodeIDs.Lower();
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
        return true;
    }
    return false;
}

occHandle(TColStd_HArray2OfReal) Ng_MeshVS_DataSource1D::getNodesCoords() const
{
    return myNodeCoords;
}
