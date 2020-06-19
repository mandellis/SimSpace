#ifndef NG_MESHVS_DEFORMEDDATASOURCE2D_H
#define NG_MESHVS_DEFORMEDDATASOURCE2D_H

//! redefinition of opencascade Handle() macro
#include "occhandle.h"

#include <Standard.hxx>
#include <Standard_Type.hxx>
#include <Standard_DefineHandle.hxx>

#include <TColStd_PackedMapOfInteger.hxx>
#include <MeshVS_DataMapOfIntegerVector.hxx>
#include <Standard_Real.hxx>
#include <MeshVS_DataSource.hxx>
#include <Standard_Boolean.hxx>
#include <Standard_Integer.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <MeshVS_EntityType.hxx>
#include <MeshVS_HArray1OfSequenceOfInteger.hxx>
#include <Standard_Address.hxx>
#include <TColStd_Array1OfInteger.hxx>
#include <MeshVS_DataSource.hxx>

class MeshVS_DataSource;
class gp_Vec;


class Ng_MeshVS_DeformedDataSource2D;
DEFINE_STANDARD_HANDLE(Ng_MeshVS_DeformedDataSource2D, MeshVS_DataSource)

//! The class provides default class which helps to represent node displacements by deformed mesh
//! This class has an internal handle to canonical non-deformed mesh data source and
//! map of displacement vectors. The displacement can be magnified to useful size.
//! All methods is implemented with calling the corresponding methods of non-deformed data source.
class Ng_MeshVS_DeformedDataSource2D: public MeshVS_DataSource
{

public:

  //! Constructor
  //! theNonDeformDS is canonical non-deformed data source, by which we are able to calculate
  //! deformed mesh geometry
  //! theMagnify is coefficient of displacement magnify
  Standard_EXPORT Ng_MeshVS_DeformedDataSource2D(const occHandle(MeshVS_DataSource)& theNonDeformDS, const Standard_Real theMagnify);

  Standard_EXPORT virtual Standard_Boolean GetGeom (const Standard_Integer ID, const Standard_Boolean IsElement, TColStd_Array1OfReal& Coords, Standard_Integer& NbNodes, MeshVS_EntityType& Type) const Standard_OVERRIDE;

  Standard_EXPORT virtual Standard_Boolean GetGeomType (const Standard_Integer ID, const Standard_Boolean IsElement, MeshVS_EntityType& Type) const Standard_OVERRIDE;

  Standard_EXPORT virtual Standard_Boolean Get3DGeom (const Standard_Integer ID, Standard_Integer& NbNodes, occHandle(MeshVS_HArray1OfSequenceOfInteger)& Data) const Standard_OVERRIDE;

  Standard_EXPORT virtual Standard_Address GetAddr (const Standard_Integer ID, const Standard_Boolean IsElement) const Standard_OVERRIDE;

  Standard_EXPORT virtual Standard_Boolean GetNodesByElement (const Standard_Integer ID, TColStd_Array1OfInteger& NodeIDs, Standard_Integer& NbNodes) const Standard_OVERRIDE;

  Standard_EXPORT virtual const TColStd_PackedMapOfInteger& GetAllNodes() const Standard_OVERRIDE;

  Standard_EXPORT virtual const TColStd_PackedMapOfInteger& GetAllElements() const Standard_OVERRIDE;

  //! This method returns map of nodal displacement vectors
  Standard_EXPORT const MeshVS_DataMapOfIntegerVector& GetVectors() const;

  //! This method sets map of nodal displacement vectors (Map).
  Standard_EXPORT void SetVectors (const MeshVS_DataMapOfIntegerVector& Map);

  //! This method returns vector ( Vect ) assigned to node number ID.
  Standard_EXPORT Standard_Boolean GetVector (const Standard_Integer ID, gp_Vec& Vect) const;

  //! This method sets vector ( Vect ) assigned to node number ID.
  Standard_EXPORT void SetVector (const Standard_Integer ID, const gp_Vec& Vect);

  Standard_EXPORT void SetNonDeformedDataSource (const occHandle(MeshVS_DataSource)& theDS);

  //! With this methods you can read and change internal canonical data source
  Standard_EXPORT occHandle(MeshVS_DataSource) GetNonDeformedDataSource() const;

  Standard_EXPORT void SetMagnify (const Standard_Real theMagnify);

  //! With this methods you can read and change magnify coefficient of nodal displacements
  Standard_EXPORT Standard_Real GetMagnify() const;

  DEFINE_STANDARD_RTTIEXT(Ng_MeshVS_DeformedDataSource2D,MeshVS_DataSource)

protected:

private:

  occHandle(MeshVS_DataSource) myNonDeformedDataSource;
  TColStd_PackedMapOfInteger myEmptyMap;
  MeshVS_DataMapOfIntegerVector myVectors;
  Standard_Real myMagnify;
};


#endif // NG_MESHVS_DEFORMEDDATASOURCE2D_H
