//! custom includes
#include "tg_meshvs_datasource.h"
#include "meshtools.h"

IMPLEMENT_STANDARD_HANDLE(Tg_MeshVS_DataSource,MeshVS_DataSourceExtended)
IMPLEMENT_STANDARD_RTTIEXT(Tg_MeshVS_DataSource,MeshVS_DataSourceExtended)

//! ------------------------------------------------------
//! function: constructor
//! details:
//! ------------------------------------------------------
Tg_MeshVS_DataSource::Tg_MeshVS_DataSource(const tetgenio &aMesh)
{
    cout<<"Tg_MeshVS_DataSource::Tg_MeshVS_DataSource->____function called____"<<endl;

    int NP = aMesh.numberofpoints;
    int NVE = aMesh.numberoftetrahedra;

    myNumberOfVolumeElements = NVE;
    myNumberOfNodes = NP;

    cout<<"Tg_MeshVS_DataSource::Tg_MeshVS_DataSource->____number of nodes: "<<NP<<endl;
    cout<<"Tg_MeshVS_DataSource::Tg_MeshVS_DataSource->____number of elements: "<<NVE<<endl;

    myNodeCoords = new TColStd_HArray2OfReal(1,NP,1,3);

    for(int i=1, S=0; i<=NP; i++)
    {
        myNodes.Add(i);
        myNodesMap.Add(i);
        for(int k=0;k<3;k++)
        {
            double coord = aMesh.pointlist[S];
            myNodeCoords->SetValue(i,k+1,coord);
            //cout<<coord<<" ";
            S++;
        }
        //cout<<endl;
    }
    //cout<<endl;

    myElemType = new TColStd_HArray1OfInteger(1,NVE);
    myElemNodes = new TColStd_HArray2OfInteger(1,NVE,1,10);

    //! ------------------------------------------
    //! map of the volume elements
    //! nodes of the element
    //! ------------------------------------------
    for(int S=0, i=1;i<=NVE;i++)
    {
        myElements.Add(i);
        myElementsMap.Add(i);
        for(int k=0; k<4; k++)
        {
            int nodeID = aMesh.tetrahedronlist[S];
            myElemNodes->SetValue(i,k+1,nodeID);
            //cout<<aMesh.tetrahedronlist[S]<<" ";
            S++;
        }
        //! set the element type
        myElemType->SetValue(i,TRIG);
    }

    //! -------------------------------------------
    //! sequence of nodes for the 3D visualization
    //! 1st order tet: definition of the faces
    //! i= 1 => {0,1,2}
    //! i= 2 => {1,2,3}
    //! i= 3 => {2,3,0}
    //! i= 4 => {3,0,1}
    //! -------------------------------------------
    TET4MeshData = new MeshVS_HArray1OfSequenceOfInteger(1,4);
    for(int i=1;i<=4;i++)
    {
        TET4MeshData->ChangeValue(i).Append((i-1)%4);
        TET4MeshData->ChangeValue(i).Append(i%4);
        TET4MeshData->ChangeValue(i).Append((i+1)%4);
    }
    //cout<<"Tg_MeshVS_DataSource::Tg_MeshVS_DataSource->____constructor OK____"<<endl;
}

//! ---------------------------------------------------------------- //
//! function: GetGeom                                                //
//! details:                                                         //
//! ---------------------------------------------------------------- //
Standard_Boolean Tg_MeshVS_DataSource::GetGeom(const Standard_Integer ID,
                                               const Standard_Boolean IsElement,
                                               TColStd_Array1OfReal& Coords,
                                               Standard_Integer& NbNodes,
                                               MeshVS_EntityType& Type) const
{
    if(IsElement && ID>=1 && ID<=myNumberOfVolumeElements)
    {
        Type = MeshVS_ET_Element;        
        NbNodes = 4;
        for(Standard_Integer i=1,k=1;i<=NbNodes;i++)
        {
            // i spans the nodes of the element
            Standard_Integer IdxNode=myElemNodes->Value(ID,i);
            // j spans {x, y, z}
            for(Standard_Integer j=1;j<=3;j++)
            {
                Coords(k)=myNodeCoords->Value(IdxNode,j);
                //cout<<"\t x"<<k<<"= "<<Coords(k)<<"\t";
                k++;
            }
            //cout<<endl;
        }
        return Standard_True;
    }
    else if(!IsElement && ID>=1 && ID<=myNumberOfNodes)
    {
        Type = MeshVS_ET_Node;
        NbNodes = 1;
        Coords(1)=myNodeCoords->Value(ID,1);
        Coords(2)=myNodeCoords->Value(ID,2);
        Coords(3)=myNodeCoords->Value(ID,3);
        return Standard_True;
    }
    return Standard_False;
}

//! ---------------------------------------------------------------- //
//! function: GetGeomType                                            //
//! details:                                                         //
//! ---------------------------------------------------------------- //
Standard_Boolean Tg_MeshVS_DataSource::GetGeomType(const Standard_Integer ID,
                                                    const Standard_Boolean IsElement,
                                                    MeshVS_EntityType& Type) const
{
    //!cout<<"____Calling GetNGeom for ID "<<ID<<"____"<<endl;
    if(IsElement && ID>=1 && ID<=myNumberOfVolumeElements)
    {
        Type = MeshVS_ET_Volume;
        return Standard_True;
    }
    else if(!IsElement && ID>=1 && ID<=myNumberOfNodes)
    {
        Type = MeshVS_ET_Node;
        return Standard_False;
    }
    return Standard_False;
}

//! ---------------------------------------------------------------- //
//! function: GetAddr                                                //
//! details:                                                         //
//! ---------------------------------------------------------------- //
Standard_Address Tg_MeshVS_DataSource::GetAddr(const Standard_Integer,
                                                const Standard_Boolean) const
{
    return NULL;
}

//! ----------------------------------------------------------------- //
//! function: GetNodesByElement                                       //
//! details: The function returns FALSE if the element does not exist //
//! ----------------------------------------------------------------- //
Standard_Boolean Tg_MeshVS_DataSource::GetNodesByElement(const Standard_Integer ID,
                                                          TColStd_Array1OfInteger& theNodeIDs,
                                                          Standard_Integer& theNbNodes) const
{
    //cout<<"____Calling GetNodesByElement for volume element ID "<<ID<<"____"<<endl;
    if(ID>=1 && ID<=myNumberOfVolumeElements)
    {
        theNbNodes = 4;
        for(int k=1; k<=theNbNodes; k++)
        {
            int nodeID = myElemNodes->Value(ID,k);
            theNodeIDs.SetValue(k,nodeID);
        }
        return Standard_True;
    }
    return Standard_False;
}

//! ---------------------------------------------------------------- //
//! function: GetAllNodes                                            //
//! details:                                                         //
//! ---------------------------------------------------------------- //
const TColStd_PackedMapOfInteger& Tg_MeshVS_DataSource::GetAllNodes() const
{
    return myNodes;
}

//! ---------------------------------------------------------------- //
//! Function: GetAllElements                                         //
//! Details:                                                         //
//! ---------------------------------------------------------------- //
const TColStd_PackedMapOfInteger& Tg_MeshVS_DataSource::GetAllElements() const
{
    return myElements;
}

//! ---------------------------------------------------------------- //
//! function: Get3DGeom                                              //
//! details: Returns information about 3D geometry of an element     //
//! ---------------------------------------------------------------- //
Standard_Boolean Tg_MeshVS_DataSource::Get3DGeom(const Standard_Integer theID,
                                                   Standard_Integer &theNbNodes,
                                                   Handle(MeshVS_HArray1OfSequenceOfInteger) &theData) const
{
    if(theID>=1 && theID<=myNumberOfVolumeElements)
    {
        theNbNodes = 4;
        theData = TET4MeshData;
        return Standard_True;
    }
    return Standard_False;
}
