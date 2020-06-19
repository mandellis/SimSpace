#ifndef NG_MESHVS_DATASOURCE1D_H
#define NG_MESHVS_DATASOURCE1D_H

//! redefinition of opencascade Handle() macro
#include "occhandle.h"

//! OCC
#include <Standard_DefineHandle.hxx>
#include <Standard_Macro.hxx>
#include <Standard.hxx>
#include <Standard_Type.hxx>
#include <Standard_Address.hxx>
#include <TColStd_PackedMapOfInteger.hxx>
#include <TColStd_IndexedMapOfInteger.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <TColStd_HArray2OfInteger.hxx>
#include <TColStd_HArray2OfReal.hxx>
#include <MeshVS_EntityType.hxx>
#include <MeshVS_HArray1OfSequenceOfInteger.hxx>

//! custom includes
#include <ng_meshvs_datasourceface.h>
#include <ng_meshvs_datasource2d.h>

//#define VERBOSE

class Ng_MeshVS_DataSource1D: public MeshVS_DataSource
{

public:

    DEFINE_STANDARD_RTTIEXT(Ng_MeshVS_DataSource1D,MeshVS_DataSource)

    Standard_EXPORT Ng_MeshVS_DataSource1D(const TColStd_PackedMapOfInteger &theNodes,
                                           const occHandle(TColStd_HArray2OfReal) &theCoords,
                                           Standard_Integer theMeshOrder);


    Standard_EXPORT Standard_Boolean GetGeom (const Standard_Integer ID,
                                              const Standard_Boolean IsElement,
                                              TColStd_Array1OfReal& Coords,
                                              Standard_Integer& NbNodes,
                                              MeshVS_EntityType& Type)
                                              const Standard_OVERRIDE;

    Standard_EXPORT Standard_Boolean GetGeomType(const Standard_Integer ID,
                                                 const Standard_Boolean IsElement,
                                                 MeshVS_EntityType& Type)
                                                 const Standard_OVERRIDE;


    Standard_EXPORT Standard_Address GetAddr(const Standard_Integer ID,
                                             const Standard_Boolean IsElement)
                                             const Standard_OVERRIDE;

    Standard_EXPORT virtual Standard_Boolean GetNodesByElement(const Standard_Integer ID,
                                                               TColStd_Array1OfInteger& NodeIDs,
                                                               Standard_Integer& theNbNodes)
                                                               const Standard_OVERRIDE;

    Standard_EXPORT const TColStd_PackedMapOfInteger& GetAllNodes() const Standard_OVERRIDE;
    Standard_EXPORT const TColStd_PackedMapOfInteger& GetAllElements() const Standard_OVERRIDE;

    occHandle(TColStd_HArray2OfReal) getNodesCoords() const;

protected:

private:

    TColStd_PackedMapOfInteger myNodes;
    TColStd_PackedMapOfInteger myElements;

    occHandle(TColStd_HArray2OfInteger) myElemNodes;
    occHandle(TColStd_HArray2OfReal) myNodeCoords;
    occHandle(TColStd_HArray2OfReal) myElemNormals;

    Standard_Integer myMeshOrder;
    Standard_Integer myNumberOfNodes;
    Standard_Integer myNumberOfEdgeElements;

public:

    TColStd_IndexedMapOfInteger myNodesMap;
    TColStd_IndexedMapOfInteger myElementsMap;
};

DEFINE_STANDARD_HANDLE(Ng_MeshVS_DataSource1D,MeshVS_DataSourceExtended)

#endif // NG_MESHVS_DATASOURCE1D_H
