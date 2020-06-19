#include "ng_meshvs_deformeddatasource2d.h"
#include <gp_Vec.hxx>
#include <MeshVS_Buffer.hxx>
#include <MeshVS_DataSource.hxx>
#include <Standard_Type.hxx>

IMPLEMENT_STANDARD_RTTIEXT(Ng_MeshVS_DeformedDataSource2D,MeshVS_DataSource)
IMPLEMENT_STANDARD_HANDLE(Ng_MeshVS_DataSource2D,MeshVS_DataSource)

//================================================================
// Function : Constructor Ng_MeshVS_DeformedDataSource2D
// Purpose  :
//================================================================
Ng_MeshVS_DeformedDataSource2D::Ng_MeshVS_DeformedDataSource2D( const occHandle(MeshVS_DataSource)& theNonDeformDS,
                                                                const Standard_Real theMagnify )
{
    myNonDeformedDataSource = theNonDeformDS;
    SetMagnify ( theMagnify );
}

//================================================================
// Function : shiftCoord
// Purpose  : auxiliary: shift coordinate of the node on a given vector
//================================================================
static inline void shiftCoord (TColStd_Array1OfReal& Coords, Standard_Integer i, const gp_Vec &aVec)
{
    Coords(3*i-2) = Coords(3*i-2) + aVec.X();
    Coords(3*i-1) = Coords(3*i-1) + aVec.Y();
    Coords(3*i)   = Coords(3*i)   + aVec.Z();
}

//================================================================
// Function : GetGeom
// Purpose  :
//================================================================
Standard_Boolean Ng_MeshVS_DeformedDataSource2D::GetGeom( const Standard_Integer ID,
                                                          const Standard_Boolean IsElement,
                                                          TColStd_Array1OfReal& Coords,
                                                          Standard_Integer& NbNodes,
                                                          MeshVS_EntityType& Type) const
{
    if ( myNonDeformedDataSource.IsNull() ||
         ! myNonDeformedDataSource->GetGeom ( ID, IsElement, Coords, NbNodes, Type ) )
        return Standard_False;

    if ( Type==MeshVS_ET_Node )
    {
        gp_Vec Vect;
        if ( ! GetVector( ID, Vect ) )
            return Standard_False;
        shiftCoord ( Coords, 1, myMagnify * Vect );
    }
    else
    {
        MeshVS_Buffer aNodesBuf (NbNodes * sizeof(Standard_Integer));
        TColStd_Array1OfInteger aNodes (aNodesBuf, 1, NbNodes);
        if ( !myNonDeformedDataSource->GetNodesByElement ( ID, aNodes, NbNodes ) )
            return Standard_False;
        for ( int i=1; i <= NbNodes; i++ )
        {
            gp_Vec Vect;
            if ( ! GetVector( aNodes(i), Vect ) )
                return Standard_False;
            shiftCoord ( Coords, i, myMagnify * Vect );
        }
    }
    return Standard_True;
}

//================================================================
// Function : GetGeomType
// Purpose  :
//================================================================
Standard_Boolean Ng_MeshVS_DeformedDataSource2D::GetGeomType( const Standard_Integer ID,
                                                              const Standard_Boolean IsElement,
                                                              MeshVS_EntityType& Type) const
{
    if (myNonDeformedDataSource.IsNull())
        return Standard_False;
    else
        return myNonDeformedDataSource->GetGeomType(ID,IsElement,Type);
}

//================================================================
// Function : Get3DGeom
// Purpose  :
//================================================================
Standard_Boolean Ng_MeshVS_DeformedDataSource2D::Get3DGeom( const Standard_Integer ID,
                                                            Standard_Integer& NbNodes,
                                                            occHandle( MeshVS_HArray1OfSequenceOfInteger )& Data ) const
{
    if( myNonDeformedDataSource.IsNull() )
        return Standard_False;
    else
        return myNonDeformedDataSource->Get3DGeom( ID, NbNodes, Data );
}

//================================================================
// Function : GetAddr
// Purpose  :
//================================================================
Standard_Address Ng_MeshVS_DeformedDataSource2D::GetAddr( const Standard_Integer ID,
                                                          const Standard_Boolean IsElement ) const
{
    if ( myNonDeformedDataSource.IsNull() )
        return 0;
    else
        return myNonDeformedDataSource->GetAddr( ID, IsElement );
}

//================================================================
// Function : GetNodesByElement
// Purpose  :
//================================================================
Standard_Boolean Ng_MeshVS_DeformedDataSource2D::GetNodesByElement(const Standard_Integer ID,
                                                                   TColStd_Array1OfInteger& NodeIDs,
                                                                   Standard_Integer& NbNodes) const
{
    if ( myNonDeformedDataSource.IsNull() )
        return Standard_False;
    else
        return myNonDeformedDataSource->GetNodesByElement( ID, NodeIDs, NbNodes );
}

//================================================================
// Function : GetAllNodes
// Purpose  :
//================================================================
const TColStd_PackedMapOfInteger& Ng_MeshVS_DeformedDataSource2D::GetAllNodes() const
{
    if ( myNonDeformedDataSource.IsNull() )
        return myEmptyMap;
    else
        return myNonDeformedDataSource->GetAllNodes();
}

//================================================================
// Function : GetAllElements
// Purpose  :
//================================================================
const TColStd_PackedMapOfInteger& Ng_MeshVS_DeformedDataSource2D::GetAllElements() const
{
    if ( myNonDeformedDataSource.IsNull() )
        return myEmptyMap;
    else
        return myNonDeformedDataSource->GetAllElements();
}

//================================================================
// Function : GetVectors
// Purpose  :
//================================================================
const MeshVS_DataMapOfIntegerVector& Ng_MeshVS_DeformedDataSource2D::GetVectors() const
{
    return myVectors;
}

//================================================================
// Function : SetVectors
// Purpose  :
//================================================================
void Ng_MeshVS_DeformedDataSource2D::SetVectors( const MeshVS_DataMapOfIntegerVector& Map )
{
    myVectors = Map;
}

//================================================================
// Function : GetVector
// Purpose  :
//================================================================
Standard_Boolean Ng_MeshVS_DeformedDataSource2D::GetVector( const Standard_Integer ID,
                                                            gp_Vec& Vect ) const
{
    Standard_Boolean aRes = myVectors.IsBound ( ID );
    if ( aRes )
        Vect = myVectors.Find ( ID );
    return aRes;
}

//================================================================
// Function : SetVector
// Purpose  :
//================================================================
void Ng_MeshVS_DeformedDataSource2D::SetVector( const Standard_Integer ID,
                                                const gp_Vec& Vect )
{
    Standard_Boolean aRes = myVectors.IsBound ( ID );
    if ( aRes )
        myVectors.ChangeFind ( ID ) = Vect;
    else
        myVectors.Bind( ID, Vect );
}

//================================================================
// Function : SetNonDeformedDataSource
// Purpose  :
//================================================================
void Ng_MeshVS_DeformedDataSource2D::SetNonDeformedDataSource( const occHandle(MeshVS_DataSource)& theDS )
{
    myNonDeformedDataSource = theDS;
}

//================================================================
// Function : GetNonDeformedDataSource
// Purpose  :
//================================================================
occHandle(MeshVS_DataSource) Ng_MeshVS_DeformedDataSource2D::GetNonDeformedDataSource() const
{
    return myNonDeformedDataSource;
}

//================================================================
// Function : SetMagnify
// Purpose  :
//================================================================
void Ng_MeshVS_DeformedDataSource2D::SetMagnify( const Standard_Real theMagnify )
{
    if ( theMagnify<=0 )
        myMagnify = 1.0;
    else
        myMagnify = theMagnify;
}

//================================================================
// Function : GetMagnify
// Purpose  :
//================================================================
Standard_Real Ng_MeshVS_DeformedDataSource2D::GetMagnify() const
{
    return myMagnify;
}
