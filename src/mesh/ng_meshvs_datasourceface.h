#ifndef NG_MESHVS_DATASOURCEFACE_H
#define NG_MESHVS_DATASOURCEFACE_H

//! -------------------------------------
//! macro for opencascade::handle<class>
//! -------------------------------------
#include "occhandle.h"

//! ----
//! OCC
//! ----
#include <MeshVS_DataSource.hxx>
#include <TColStd_HArray1OfInteger.hxx>
#include <TColStd_IndexedMapOfInteger.hxx>
#include <TColStd_HArray2OfReal.hxx>
#include <TColStd_HArray2OfInteger.hxx>
#include <Poly_Triangulation.hxx>
#include <Poly_Triangle.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Face.hxx>

//! ---
//! Qt
//! ---
#include <QList>
#include <QMap>
#include <QMetaType>

//! ----------------
//! custom includes
//! ----------------
#include <meshelement2d.h>
#include <StlMesh_Mesh.hxx>
#include "src/utils/hash_c.h"
#include <mesh.h>
#include <ng_meshvs_datasource2d.h>

//! ------
//! Eigen
//! ------
#include <ext/eigen/Eigen/Dense>
#include <ext/eigen/Eigen/Core>
#include <ext/eigen/Eigen/Sparse>

namespace nglib {
    #include <nglib.h>
}
using namespace nglib;

//! -------
//! Tegten
//! -------
#include <tetgen.h>

//! ----
//! C++
//! ----
#include <set>
#include <memory>

class Ng_MeshVS_DataSourceFace: public MeshVS_DataSource
{
public:

    DEFINE_STANDARD_RTTIEXT(Ng_MeshVS_DataSourceFace,MeshVS_DataSource)

    //! ----------------------
    //! constructor - default
    //! ----------------------
    Standard_EXPORT Ng_MeshVS_DataSourceFace();

    //! --------------------------------
    //! constructor from Netgen pointer
    //! --------------------------------
    Standard_EXPORT Ng_MeshVS_DataSourceFace(Ng_Mesh *aMesh, int faceNr);

    //! ---------------------------------------------------------------------------
    //! constructor: build the face mesh datasource from Tetgen pointer and faceNr
    //! ---------------------------------------------------------------------------
    Standard_EXPORT Ng_MeshVS_DataSourceFace(tetgenio *aMesh, int faceNr);

    //! --------------------------------------------------------------
    //! constructor: build the face mesh datasource from StlMesh_Mesh
    //! --------------------------------------------------------------
    Standard_EXPORT Ng_MeshVS_DataSourceFace(const occHandle(StlMesh_Mesh) &anStlMesh_Mesh,
                                             const QList<QMap<int,int>> map_localToGlobal_nodeID,
                                             const QList<QMap<int,int>> map_localToGlobal_elementID,
                                             int faceNr);

    //! ---------------------------------------------------------
    //! constructor - build a Mesh data source from StlMesh_Mesh
    //! ---------------------------------------------------------
    Standard_EXPORT Ng_MeshVS_DataSourceFace(const occHandle(StlMesh_Mesh) &aMesh);

    //! --------------------------------------------------
    //! constructor - from a vector of elements and nodes
    //! --------------------------------------------------
    Standard_EXPORT Ng_MeshVS_DataSourceFace(const QString &faceFileName,
                                             const QMap<int,mesh::meshPoint> indexedMapOfAllMeshPoints,
                                             int faceNr);

    //! -----------------------------------------
    //! constructor from "triangle" library data
    //! -----------------------------------------
    Standard_EXPORT Ng_MeshVS_DataSourceFace(const QString &nodeFileName,
                                             const QString &faceFileName);

    //! ------------------------------
    //! constructor from "Eigen" data
    //! ------------------------------
    Standard_EXPORT Ng_MeshVS_DataSourceFace(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F);

    //! ---------------------------------------------------
    //! constructor: from a list of face mesh data sources
    //! ---------------------------------------------------
    Standard_EXPORT Ng_MeshVS_DataSourceFace(const std::vector<occHandle(Ng_MeshVS_DataSourceFace)> &faceDSList);

    //! --------------------------------------------
    //! constructor: from a volume mesh data source
    //! --------------------------------------------
    Standard_EXPORT Ng_MeshVS_DataSourceFace(const opencascade::handle<Ng_MeshVS_DataSource3D> &aVolumeMeshDS);

    //! -------------------------------------------------------------------
    //! constructor - build the face mesh data source using mesh proximity
    //! -------------------------------------------------------------------
    Standard_EXPORT Ng_MeshVS_DataSourceFace(const TopoDS_Face &aFace,
                                             const occHandle(MeshVS_DataSource) &surfaceMesh,
                                             bool automaticProximity = false,
                                             double proximityValue = 1.0);

    //! --------------------------
    //! constructor - from a file
    //! --------------------------
    Standard_EXPORT Ng_MeshVS_DataSourceFace(const std::string &fileName);

    //! ---------------------------------------------------------------------------------
    //! constructor - from a list of (2D) mesh elements, each one defined by coordinates
    //! ---------------------------------------------------------------------------------
    Standard_EXPORT Ng_MeshVS_DataSourceFace(const std::vector<meshElementByCoords> &meshElements, bool autoRenumberElements=true, bool autoRenumberNodes=true);

    //! ---------------------------------------------------------------
    //! questo è un accessorio per ora indispensabile, ma da rimuovere
    //! ---------------------------------------------------------------
    Standard_EXPORT Ng_MeshVS_DataSourceFace(const occHandle(Ng_MeshVS_DataSource2D) &surfaceMesh);

    //! ---------------------------------------
    //! constructor - clone a mesh data source
    //! ---------------------------------------
    Standard_EXPORT Ng_MeshVS_DataSourceFace(const occHandle(Ng_MeshVS_DataSourceFace) &aFaceMesh);

    //! --------------------------------------------------
    //! constructor - from a vector of elements and nodes
    //! --------------------------------------------------
    Standard_EXPORT Ng_MeshVS_DataSourceFace(const std::vector<mesh::meshElement> &elements, std::vector<mesh::meshPoint> &nodes);

    //! -------------------------------------------------------------------
    //! constructor - retrieve the STL mesh from a TopoDS_Shape and faceNr
    //! -------------------------------------------------------------------
    Standard_EXPORT Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace(const TopoDS_Shape &shape, int faceNr);

    //! --------------------------------
    //! constructor - from a C++ stream
    //! --------------------------------
    Standard_EXPORT Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace(ifstream &stream);

    //! -------------------------------------
    //! convert aMesh into a n-th order mesh
    //! -------------------------------------
    Standard_EXPORT Ng_MeshVS_DataSourceFace(const occHandle(Ng_MeshVS_DataSourceFace) &aMesh, int order);

    //! -----------------------------------------------------------------------
    //! constructor: a mesh from another mesh with nodal displacements applied
    //! -----------------------------------------------------------------------
    Standard_EXPORT Ng_MeshVS_DataSourceFace(const occHandle(Ng_MeshVS_DataSource3D) &aMesh, const QMap<int,gp_Vec> &displacements);

    //! -------------------------------
    //! constructor: meshDSA - meshDSB
    //! -------------------------------
    Standard_EXPORT Ng_MeshVS_DataSourceFace::Ng_MeshVS_DataSourceFace(const occHandle(Ng_MeshVS_DataSourceFace) &meshDSA, occHandle(Ng_MeshVS_DataSourceFace) &meshDSB);

    //! -----------------------
    //! write the mesh on disk
    //! -----------------------
    Standard_EXPORT bool writeMesh(const QString &meshFileName, int dim);
    Standard_EXPORT bool write(ofstream &os, int dim);

    //! ------------------------
    //! read the mesh from disk
    //! ------------------------
    Standard_EXPORT bool readMesh(const QString &meshFileName,
                                  std::vector<mesh::meshPoint> &nodes,
                                  std::vector<mesh::meshElement> &elements);

    //! use for ID the global element/node number
    Standard_EXPORT Standard_Boolean GetGeom(const Standard_Integer ID,
                                             const Standard_Boolean IsElement,
                                             TColStd_Array1OfReal& Coords,
                                             Standard_Integer& NbNodes,
                                             MeshVS_EntityType& Type) const Standard_OVERRIDE;

    //! use for ID the global element/node number
    Standard_EXPORT Standard_Boolean GetGeomType(const Standard_Integer ID,
                                                 const Standard_Boolean IsElement,
                                                 MeshVS_EntityType& Type) const Standard_OVERRIDE;

    Standard_EXPORT Standard_Address GetAddr(const Standard_Integer ID,
                                             const Standard_Boolean IsElement) const Standard_OVERRIDE;

    //! use for ID the global element number
    Standard_EXPORT virtual Standard_Boolean GetNodesByElement(const Standard_Integer ID,
                                                               TColStd_Array1OfInteger& NodeIDs,
                                                               Standard_Integer& theNbNodes) const Standard_OVERRIDE;

    //! get element type
    Standard_EXPORT virtual const bool GetElementType(ElemType &eType, int elementID, bool isLocal) const;

    //! Get3DGeom
    Standard_EXPORT virtual Standard_Boolean Get3DGeom (const Standard_Integer theID,
                                                        Standard_Integer& NbNodes,
                                                        occHandle(MeshVS_HArray1OfSequenceOfInteger) &Data) const;

    //! get all nodes
    Standard_EXPORT virtual const TColStd_PackedMapOfInteger& GetAllNodes() const Standard_OVERRIDE;

    //! get all elements
    Standard_EXPORT virtual const TColStd_PackedMapOfInteger& GetAllElements() const Standard_OVERRIDE;

    //! get normal
    Standard_EXPORT virtual Standard_Boolean GetNormal (const Standard_Integer Id,
                                                        const Standard_Integer Max,
                                                        Standard_Real& nx,
                                                        Standard_Real& ny,
                                                        Standard_Real& nz)
                                                        const Standard_OVERRIDE;

    occHandle(TColStd_HArray2OfReal) getNodesCoordinates() const { return myNodeCoords; }

    std::vector<mesh::meshSegment> getSegmentsOfElement(int elementID, bool isLocal);

protected:

private:

    TColStd_PackedMapOfInteger myNodes;
    TColStd_PackedMapOfInteger myElements;
    Standard_Integer myNumberOfElements;
    Standard_Integer myNumberOfNodes;

    occHandle(TColStd_HArray1OfInteger) myElemType;
    occHandle(TColStd_HArray2OfInteger) myElemNodes;
    occHandle(TColStd_HArray2OfReal) myNodeCoords;
    occHandle(TColStd_HArray2OfReal) myElemNormals;

private:

    void toEigen(Eigen::MatrixXd &V, Eigen::MatrixXi &F);

public:

    TColStd_IndexedMapOfInteger myNodesMap;
    TColStd_IndexedMapOfInteger myElementsMap;

    //! ------------------------------------
    //! mesh connectivity: node to elements
    //! mesh node normals
    //! ------------------------------------
    QMap<int,QList<double>> myNodeNormals;                  //! key: node ID; value: coefficients of the normal
    QMap<int,QList<int>> myNodeToElements;                  //! connectivity. Key: nodeID; value: list of element IDs

    //! -----------------------------
    //! compute the normals at nodes
    //! -----------------------------
    void computeNormalAtNodes();

    //! ---------------------------
    //! compute normal at elements
    //! ---------------------------
    void computeNormalAtElements();

    //! ---------------------------------------
    //! mesh connectivity: segment to elements
    //! free 1D mesh boundaries
    //! ---------------------------------------
    QMap<mesh::meshSegment,QList<int>> mySegmentToElement;
    QList<mesh::meshSegment> myBoundarySegments;
    QList<int> myBoundaryPoints;                                                    //! list of points (ID number) on the 1D boundary
    void computeFreeMeshSegments();

    //! ------------------------------
    //! node to elements connectivity
    //! ------------------------------
    QMap<int,std::set<mesh::meshSegment>> myNodeToSegmentConnectivity;
    void computeNodeToElementsConnectivity();

    //! -----------------------------
    //! node to segment connectivity
    //! -----------------------------
    void computeNodeToSegmentConnectivity();

    //! --------------------------------------------------
    //! deform the mesh using a vector displacement field
    //! --------------------------------------------------
    void displaceMySelf(const QMap<int, QList<double> > &displacementField);

    //! ------------------------------------------------
    //! compute the gradient of the normal at each node
    //! ------------------------------------------------
    QMap<int,QList<double>> myAnglesAtNode;
    void computeAnglesAtNode();

    //! -------------------------------
    //! compute the discrete curvature
    //! -------------------------------
    QMap<int,double> myCurvature;
    QMap<int,double> myCurvatureGradient;
    void computeDiscreteCurvature(int type);

    //! -----------------------------------------------
    //! change the coordinates of a mesh node - helper
    //! -----------------------------------------------
    bool changeNodeCoords(int localNodeID, const std::vector<double> &coords);

    //! --------------------------------------------------------------
    //! get the coordinates of a node: use the local node ID as input
    //! --------------------------------------------------------------
    std::vector<double> getNodeCoordinates(int localNodeID);

    //! ------------------------------------------------------------------
    //! renumber nodes: assign new global node IDs to the face mesh nodes
    //! ------------------------------------------------------------------
    void renumberNodes(const std::map<size_t,int> &mapCoordNodeNumbers,
                       const std::vector<mesh::tolerantPoint> &vecTolPoints = std::vector<mesh::tolerantPoint>(),
                       const std::vector<int> &GlobalNodeIDs = std::vector<int>());

    //! ----------------------------------
    //! setNodeNormal: an access function
    //! ----------------------------------
    void setNodeNormal(int globalNodeID, double *n);

    //! ----------------------
    //! build tolerant points
    //! ----------------------
    std::vector<mesh::tolerantPoint> buildTolerantPoints(double tolerance = 1e-6);

    //! ------------------------
    //! get surrounding point
    //! (uses global numbering)
    //! ------------------------
    bool getSurroundingNodes(int globalNodeID, std::vector<int> &surroundingNodes_global);

    //! -----------------
    //! ARAP deformation
    //! -----------------
    void displaceMySelf_asRigidAsPossible(const QMap<int, gp_Vec> &displacementField);

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

    //! ---------------
    //! add an element
    //! ---------------
    bool addElement(const meshElementByCoords &anElement);
};

DEFINE_STANDARD_HANDLE(Ng_MeshVS_DataSourceFace, MeshVS_DataSource)
Q_DECLARE_METATYPE(Ng_MeshVS_DataSourceFace)

#endif // NG_MESHVS_DATASOURCEFACE_H
