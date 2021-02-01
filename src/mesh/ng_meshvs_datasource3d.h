#ifndef NG_MESHVS_DATASOURCE3D_H
#define NG_MESHVS_DATASOURCE3D_H


//! ----------------
//! custom includes
//! ----------------
#include "occhandle.h"
#include <meshelement2d.h>
#include <mesh.h>
#include <simulationdata.h>
#include <elementtypes.h>

//! ----
//! OCC
//! ----
#include <Standard_DefineHandle.hxx>
#include <TColStd_IndexedMapOfInteger.hxx>
#include <TColStd_HArray1OfInteger.hxx>
#include <TColStd_HArray2OfInteger.hxx>
#include <TColStd_HArray2OfReal.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <TColStd_Array1OfInteger.hxx>
#include <MeshVS_EntityType.hxx>
#include <MeshVS_DataSource.hxx>
#include <MeshVS_HArray1OfSequenceOfInteger.hxx>

//! ---
//! Qt
//! ---
#include <QList>
#include <QTime>
#include <QMap>
#include <QHash>
#include <QMetaType>

//! ----
//! C++
//! ----
#include <algorithm>

//! ------
//! Eigen
//! ------
#include <ext/eigen/Eigen/Dense>

//! ------
//! nglib
//! ------
namespace nglib
{
    #include <nglib.h>
}
using namespace nglib;

//! -------
//! Tetgen
//! -------
#include <tetgen.h>

#include <mesh.h>
#include <meshelementbycoords.h>

class Ng_MeshVS_DataSource3D: public MeshVS_DataSource
{
public:

    DEFINE_STANDARD_RTTIEXT(Ng_MeshVS_DataSource3D,MeshVS_DataSource)

    //! destructor
    virtual ~Ng_MeshVS_DataSource3D(){;}

    //! constructor - default
    Standard_EXPORT Ng_MeshVS_DataSource3D(){;}

    //! constructor I - build a Mesh data source from Netgen
    Standard_EXPORT Ng_MeshVS_DataSource3D(Ng_Mesh *aMesh);

    //! constructor II - build a Mesh data source from file
    Standard_EXPORT Ng_MeshVS_DataSource3D(const std::string &fileName);

    //! constructor III - build a Mesh data source from Tetgen
    Standard_EXPORT Ng_MeshVS_DataSource3D(const tetgenio &aMesh);

    //! contructor IV - build a Mesh data source from a .ele Tetgen file
    Standard_EXPORT Ng_MeshVS_DataSource3D(const QString &tetgenEleFileName, const QString &tetgenNodeFileName);

    //! constructor V - constructor for a list of "meshElementsByCoords"
    Standard_EXPORT Ng_MeshVS_DataSource3D(const std::vector<meshElementByCoords> &vecElements, bool autoNumberElements=true, bool autoNumberNodes=true);

    //! constructor VI - for merging a prismatic mesh and a volume mesh
    Standard_EXPORT Ng_MeshVS_DataSource3D(const occHandle(Ng_MeshVS_DataSource3D) &aMesh1, const occHandle(Ng_MeshVS_DataSource3D) &aMesh2);

    //! constructor VII - from Eigen matrices
    Standard_EXPORT Ng_MeshVS_DataSource3D(const Eigen::MatrixXd &V, const Eigen::MatrixXi &T);

    //! constructor VIII - from vectors of elements and nodes
    Standard_EXPORT Ng_MeshVS_DataSource3D(const std::vector<mesh::meshElement> &elements, const std::vector<mesh::meshPoint> &nodes);

    //! constructor XI - make a n-th order mesh
    Standard_EXPORT Ng_MeshVS_DataSource3D(const occHandle(Ng_MeshVS_DataSource3D) &aMesh, int order);


    //! constructor X - a mesh from another mesh with nodal displacements applied
    Standard_EXPORT Ng_MeshVS_DataSource3D(const occHandle(Ng_MeshVS_DataSource3D) &aMesh, const QMap<int,gp_Vec> &displacements);

    //! upcast
    Standard_EXPORT Ng_MeshVS_DataSource3D(const occHandle(MeshVS_DataSource) &aMeshDS);

    Standard_EXPORT virtual Standard_Boolean GetGeom (const Standard_Integer ID,
                                                      const Standard_Boolean IsElement,
                                                      TColStd_Array1OfReal& Coords,
                                                      Standard_Integer& NbNodes,
                                                      MeshVS_EntityType& Type) const Standard_OVERRIDE;

    Standard_EXPORT virtual Standard_Boolean GetGeomType(const Standard_Integer ID,
                                                         const Standard_Boolean IsElement,
                                                         MeshVS_EntityType& Type) const Standard_OVERRIDE;

    Standard_EXPORT virtual Standard_Address GetAddr(const Standard_Integer ID,
                                                     const Standard_Boolean IsElement) const Standard_OVERRIDE;

    Standard_EXPORT virtual Standard_Boolean GetNodesByElement(const Standard_Integer ID,
                                                               TColStd_Array1OfInteger& NodeIDs,
                                                               Standard_Integer& theNbNodes) const Standard_OVERRIDE;

    Standard_EXPORT virtual const TColStd_PackedMapOfInteger& GetAllNodes() const Standard_OVERRIDE;

    Standard_EXPORT virtual const TColStd_PackedMapOfInteger& GetAllElements() const Standard_OVERRIDE;

    Standard_EXPORT virtual const bool GetElementType(ElemType &eType, int elementID, bool isLocal) const;

    Standard_EXPORT virtual Standard_Boolean Get3DGeom (const Standard_Integer theID,
                                                        Standard_Integer& NbNodes,
                                                        occHandle(MeshVS_HArray1OfSequenceOfInteger) &Data) const;


    Standard_EXPORT virtual Standard_Boolean GetNormalsByElement (const Standard_Integer Id,
                                                                  const Standard_Boolean IsNodal,
                                                                  const Standard_Integer MaxNodes,
                                                                  occHandle(TColStd_HArray1OfReal)& Normals) const Standard_OVERRIDE;

    Standard_EXPORT virtual Standard_Boolean GetNormal (const Standard_Integer Id,
                                                        const Standard_Integer Max,
                                                        Standard_Real& nx,
                                                        Standard_Real& ny,
                                                        Standard_Real& nz) const Standard_OVERRIDE;

    //! -----------------------------
    //! write on disk using C syntax
    //! -----------------------------
    Standard_EXPORT bool writeMesh(const QString &meshFileName, int dim=3);

    //! ---------------------------
    //! write into C++ file stream
    //! ---------------------------
    Standard_EXPORT void writeIntoStream(std::ofstream &os);

    //! ------------------------------
    //! read from file using C syntax
    //! ------------------------------
    Standard_EXPORT bool readMesh(const QString &meshFileName,
                                  std::vector<mesh::meshPoint> &nodes,
                                  std::vector<mesh::meshElement> &elements);

    //! -------------------------------------
    //! read from disk using C++ file stream
    //! -------------------------------------
    Standard_EXPORT bool readFromStream(std::ifstream &stream,
                                        std::vector<mesh::meshPoint> &nodes,
                                        std::vector<mesh::meshElement> &elements);

protected:

private:

    TColStd_PackedMapOfInteger myNodes;
    TColStd_PackedMapOfInteger myElements;
    occHandle(TColStd_HArray1OfInteger) myElemType;
    occHandle(TColStd_HArray2OfInteger) myElemNodes;
    occHandle(TColStd_HArray2OfReal) myElemNormals;

    occHandle(TColStd_HArray2OfReal) myNodeCoords;
    int myNumberOfElements;
    int myNumberOfNodes;

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

    QTime clock();

public:

    TColStd_IndexedMapOfInteger myNodesMap;
    TColStd_IndexedMapOfInteger myElementsMap;

    //! ------------------------------
    //! compute the normal at element
    //! ------------------------------
    void computeNormalAtElements();

    //! -------------------------------------------------
    //! face to element connectivity
    //! it returns also the face elements on the surface
    //! -------------------------------------------------
    QMap<meshElement2D, QList<int>> myFaceToElements;
    std::vector<meshElement2D> mySurfaceElements;
    void buildFaceToElementConnectivity();

    //! -------------
    //! experimental
    //! -------------
    //void buildCCXFaceToElementConnectivity(std::map<meshElement2D, std::vector<std::pair<int, std::string> > > &map);
    void buildCCXFaceToElementConnectivity(std::map<meshElement2D, std::vector<std::pair<int, int> > > &map);

    //! ---------------------------------------
    //! change the coordinates of a mesh point
    //! ---------------------------------------
    bool changeNodeCoords(int localNodeID, const std::vector<double> &coords);

    //! ---------------------------------------
    //! get the coordinates of a node - helper
    //! ---------------------------------------
    std::vector<double> getNodeCoordinates(int localNodeID);

    //! ------------------
    //! remove an element
    //! ------------------
    void removeElement(int localElementID);

    //! ----------------
    //! displace myself
    //! ----------------
    void displaceMySelf(const QMap<int, QList<double> > &displacementField);
    void displaceMySelf_asRigidAsPossible(const QMap<int, QList<double> > &displacementField, const std::vector<int> &nodeGroup, int mode);

    //! ----------------------------------------
    //! compute element to element connectivity
    //! ----------------------------------------
    void buildElementToElementConnectivity(std::map<int,std::vector<int>> &elementToAttachedElements);

    //! ------------------------------------------
    //! get the face elements of a volume element
    //! ------------------------------------------
    std::vector<meshElement2D> getFaceElementsOfElement(int localElementID);

    //! ---------------
    //! renumber nodes
    //! ---------------
    void renumberNodes(const std::map<std::size_t,int> &mapCoordNodeNumbers,
                       const std::vector<mesh::tolerantPoint> &vecTolPoints,
                       const std::vector<int> &GlobalNodeIDs);

    //! ----------------------
    //! build tolerant points
    //! ----------------------
    std::vector<mesh::tolerantPoint> buildTolerantPoints(double tolerance = 1e-6);

    //! ------------------------
    //! get segments of element
    //! ------------------------
    std::vector<mesh::meshSegment> getSegmentsOfElement(int elementID, bool isLocal);
};

DEFINE_STANDARD_HANDLE(Ng_MeshVS_DataSource3D, MeshVS_DataSource)
Q_DECLARE_METATYPE(occHandle(Ng_MeshVS_DataSource3D))

#endif // NG_MESHVS_DATASOURCE3D_H
