#ifndef Ng_MeshVS_DataSource2D_H
#define Ng_MeshVS_DataSource2D_H

#define USE_POLYGON

//! -------------------------------------------
//! redefinition of opencascade Handle() macro
//! -------------------------------------------
#include "occhandle.h"

//! ----
//! OCC
//! ----
#include <Standard_DefineHandle.hxx>
#include <TColStd_PackedMapOfInteger.hxx>
#include <TColStd_IndexedMapOfInteger.hxx>
#include <TColStd_MapIteratorOfPackedMapOfInteger.hxx>
#include <TColStd_Array1OfInteger.hxx>
#include <TColStd_HArray2OfInteger.hxx>
#include <TColStd_HArray1OfInteger.hxx>
#include <TColStd_HArray2OfReal.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <MeshVS_EntityType.hxx>
#include <MeshVS_DataSource.hxx>
#include <MeshVS_DeformedDataSource.hxx>
#include <MeshVS_HArray1OfSequenceOfInteger.hxx>
#include <StlMesh_Mesh.hxx>

//! ----------------
//! custom includes
//! ----------------
#include <elementtypes.h>
#include <mesh.h>
#include <meshelement2d.h>
#include <ng_meshvs_datasource3d.h>
#include <simulationdata.h>

//! ----
//! C++
//! ----
#include <vector>
#include <memory>

//! ---
//! Qt
//! ---
#include <QMetaType>

namespace nglib {
    #include <nglib.h>
}
using namespace nglib;

#include <tetgen.h>

class Ng_MeshVS_DataSource2D: public MeshVS_DataSource
{
public:

    DEFINE_STANDARD_RTTIEXT(Ng_MeshVS_DataSource2D,MeshVS_DataSource)

    //! destructor
    virtual ~Ng_MeshVS_DataSource2D(){;}

    //! default constructor
    Standard_EXPORT Ng_MeshVS_DataSource2D(){;}

    //! constructor - build a Mesh data source from Netgen
    Standard_EXPORT Ng_MeshVS_DataSource2D(Ng_Mesh *aMesh);

    //! constructor - build a Mesh data source from file
    Standard_EXPORT Ng_MeshVS_DataSource2D(const std::string &fileName);

    //! constructor - build a Mesh data source from Tetgen
    Standard_EXPORT Ng_MeshVS_DataSource2D(const tetgenio &aMesh);

    //! constructor - build a Mesh data source from StlMesh_Mesh
    Standard_EXPORT Ng_MeshVS_DataSource2D(const occHandle(StlMesh_Mesh) &aMesh);

    //! contructor - build a Mesh data source from a .ele Tetgen file
    Standard_EXPORT Ng_MeshVS_DataSource2D(const QString &faceFileName, const QString &nodeFileName);

    //! constructor - constructor from a vector of mesh elements and a vector of mesh points
    Standard_EXPORT Ng_MeshVS_DataSource2D(const std::vector<mesh::meshElement> &elements, std::vector<mesh::meshPoint> &nodes);

    //! constructor from a volume mesh - topological
    Standard_EXPORT Ng_MeshVS_DataSource2D(const occHandle(Ng_MeshVS_DataSource3D) &a3DMeshVS_DataSource);

    //! conversion into second order
    Standard_EXPORT Ng_MeshVS_DataSource2D(const occHandle(Ng_MeshVS_DataSource2D) &aMesh, int order);

    //! constructor: a mesh from another mesh with nodal displacements applied
    Standard_EXPORT Ng_MeshVS_DataSource2D(const occHandle(Ng_MeshVS_DataSource3D) &aMesh, const QMap<int,gp_Vec> &displacements);

    Standard_EXPORT virtual Standard_Boolean GetGeom (const Standard_Integer ID,
                                              const Standard_Boolean IsElement,
                                              TColStd_Array1OfReal& Coords,
                                              Standard_Integer& NbNodes,
                                              MeshVS_EntityType& Type)
                                              const Standard_OVERRIDE;

    Standard_EXPORT virtual Standard_Boolean GetGeomType(const Standard_Integer ID,
                                                 const Standard_Boolean IsElement,
                                                 MeshVS_EntityType& Type)
                                                 const Standard_OVERRIDE;

    Standard_EXPORT virtual Standard_Address GetAddr(const Standard_Integer ID,
                                             const Standard_Boolean IsElement)
                                             const Standard_OVERRIDE;

    Standard_EXPORT virtual Standard_Boolean GetNodesByElement(const Standard_Integer ID,
                                                               TColStd_Array1OfInteger& NodeIDs,
                                                               Standard_Integer& theNbNodes)
                                                               const Standard_OVERRIDE;

    Standard_EXPORT virtual const TColStd_PackedMapOfInteger& GetAllNodes() const Standard_OVERRIDE;

    Standard_EXPORT virtual const TColStd_PackedMapOfInteger& GetAllElements() const Standard_OVERRIDE;

    Standard_EXPORT virtual Standard_Boolean Get3DGeom (const Standard_Integer theID,
                                                        Standard_Integer& NbNodes,
                                                        occHandle(MeshVS_HArray1OfSequenceOfInteger) &Data) const;

    Standard_EXPORT virtual const bool GetElementType(ElemType &eType, int elementID, bool isLocal) const;

    Standard_EXPORT virtual Standard_Boolean GetNormal (const Standard_Integer Id,
                                                        const Standard_Integer Max,
                                                        Standard_Real& nx,
                                                        Standard_Real& ny,
                                                        Standard_Real& nz)
                                                        const Standard_OVERRIDE;

    Standard_EXPORT bool writeMesh(const QString &meshFileName, int dim=2);
    Standard_EXPORT bool readMesh(const QString &meshFileName,
                                  std::vector<mesh::meshPoint> &nodes,
                                  std::vector<mesh::meshElement> &elements);

    Standard_EXPORT void renumberNodes(const std::map<size_t,int> &mapCoordNodeNumbers,
                                       const std::vector<mesh::tolerantPoint> &vecTolPoints,
                                       const std::vector<int> &GlobalNodeIDs);

protected:

private:

    TColStd_PackedMapOfInteger myNodes;
    TColStd_PackedMapOfInteger myElements;
    occHandle(TColStd_HArray1OfInteger) myElemType;
    occHandle(TColStd_HArray2OfInteger) myElemNodes;
    occHandle(TColStd_HArray2OfReal) myNodeCoords;
    occHandle(TColStd_HArray2OfReal) myElemNormals;
    Standard_Integer myNumberOfElements;
    Standard_Integer myNumberOfNodes;

public:

    TColStd_IndexedMapOfInteger myNodesMap;
    TColStd_IndexedMapOfInteger myElementsMap;

    //! -----------------------------------
    //! compute the normal of each element
    //! -----------------------------------
    void computeNormalAtElements();

    //! ------------------------------------
    //! mesh connectivity: node to elements
    //! ------------------------------------
    QMap<int,QList<double>> myNodeNormals;                  //! key: node ID; value: components of the normal vevctor

    //! --------------------------------------------------
    //! deform the mesh using a vector displacement field
    //! --------------------------------------------------
    void displaceMySelf(const QMap<int, gp_Vec> &displacementField);

    //! ------------------
    //! mesh node normals
    //! ------------------
    QMap<int,QList<int>> myNodeToElements;                  //! connectivity. Key: nodeID; value: list of element IDs
    void computeNormalAtNodes();

    //! ------------------------------
    //! node to elements connectivity
    //! ------------------------------
    //QMap<int,QList<int>> myNodeToElements;
    void computeNodeToElementsConnectivity();

    //! -----------------------------------------------
    //! change the coordinates of a mesh node - helper
    //! -----------------------------------------------
    bool changeNodeCoords(int localNodeID, const std::vector<double> &coords);

    //! ---------------------------------------
    //! get the coordinates of a node - helper
    //! ---------------------------------------
    std::vector<double> getNodeCoordinates(int localNodeID);

    //! ----------------------
    //! build tolerant points
    //! ----------------------
    std::vector<mesh::tolerantPoint> buildTolerantPoints(double tolerance = 1e-6);

    //! ------------------------
    //! get segments of element
    //! ------------------------
    std::vector<mesh::meshSegment> getSegmentsOfElement(int elementID, bool isLocal);

    //! ---------------
    //! mid side nodes
    //! ---------------
    std::shared_ptr<std::vector<int>> myMidSideNodes;
    void getMidSideNodes(std::vector<int> *midSideNodes) { midSideNodes=myMidSideNodes.get(); }

    //! ------------------
    //! elements topology
    //! ------------------
    occHandle(MeshVS_HArray1OfSequenceOfInteger) TRIGMeshData;
    occHandle(MeshVS_HArray1OfSequenceOfInteger) QUADMeshData;
    occHandle(MeshVS_HArray1OfSequenceOfInteger) TRIG6MeshData;
    occHandle(MeshVS_HArray1OfSequenceOfInteger) QUAD8MeshData;
    occHandle(MeshVS_HArray1OfSequenceOfInteger) TET4MeshData;
    occHandle(MeshVS_HArray1OfSequenceOfInteger) TET10MeshData;
    occHandle(MeshVS_HArray1OfSequenceOfInteger) HEXA8MeshData;
    occHandle(MeshVS_HArray1OfSequenceOfInteger) HEXA20MeshData;
    occHandle(MeshVS_HArray1OfSequenceOfInteger) PRISM6MeshData;
    occHandle(MeshVS_HArray1OfSequenceOfInteger) PRISM15MeshData;
    occHandle(MeshVS_HArray1OfSequenceOfInteger) PYRAM5MeshData;
    occHandle(MeshVS_HArray1OfSequenceOfInteger) PYRAM13MeshData;
    void buildElementsTopology();
};

DEFINE_STANDARD_HANDLE(Ng_MeshVS_DataSource2D,MeshVS_DataSource)
Q_DECLARE_METATYPE(occHandle(Ng_MeshVS_DataSource2D))

#endif // Ng_MeshVS_DataSource2D_H
