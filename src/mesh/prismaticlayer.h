#ifndef PRISMATICLAYER_H
#define PRISMATICLAYER_H

//! ------------------------------------------------------------
//! Definition:
//! "prismatic face" is a face that will undergo mesh inflation
//! ------------------------------------------------------------

//! ----------------
//! custom includes
//! ----------------
#include <meshdatabase.h>
#include <ng_meshvs_datasourceface.h>
#include <ng_meshvs_datasource3d.h>
#include <prismaticlayerparameters.h>
#include "qprogressindicator.h"

//! ----
//! OCC
//! ----
#include <MeshVS_DeformedDataSource.hxx>
#include <MeshVS_Mesh.hxx>

//! ---
//! Qt
//! ---
#include <QMap>
#include <QList>
#include <QString>

//! ----
//! C++
//! ----
#include <vector>
#include <map>
#include <iostream>
using namespace std;

class prismaticLayer
{
public:

    //! ------------
    //! constructor
    //! ------------
    prismaticLayer(meshDataBase *mDB, QProgressIndicator *aProgressIndicator = Q_NULLPTR);

    //! ---------------
    //! set parameters
    //! ---------------
    void setParameters(prismaticLayerParameters parameters);

    //! ------------------------------------------------------------
    //! function merge face meshes: initialize "myOverallSumMeshDS"
    //! and "myPrismaticFacesSumMeshDS"
    //! ------------------------------------------------------------
    bool mergeFaceMeshes();

    //! -------------
    //! set database
    //! -------------
    void setDataBase(meshDataBase *mDB) { myMeshDB = mDB; }

    //! ---------
    //! set body
    //! ---------
    void setBody(int bodyIndex) { myBodyIndex = bodyIndex; }

    //! ---------------------------------------------
    //! set prismatic faces
    //! (bodyIndex, list of prismatic faces on body)
    //! ---------------------------------------------
    //void setPrismaticFaces(const QList<int> &prismaticFaces);
    void setPrismaticFaces(const std::vector<int> &prismaticFaces);


    //! -----------------------------
    //! get all the thickness values
    //! -----------------------------
    inline std::vector<double> getLayerThickness() { return myLayerThickness; }

    //! ------------------------------
    //! build the displaced mesh list
    //! ------------------------------
    bool inflateMesh(QList<occHandle(Ng_MeshVS_DataSourceFace)> &theInflatedMeshes);
    bool inflateMeshAndCompress(QList<occHandle(Ng_MeshVS_DataSourceFace)> &theInflatedMeshes,
                                occHandle(Ng_MeshVS_DataSource3D) &preInflationVolumeMeshDS);

    //! -----------------------------
    //! build the prismatic elements
    //! -----------------------------
    bool buildPrismaticElements(const QList<occHandle(Ng_MeshVS_DataSourceFace)> &theInflatedMeshes,
                                occHandle(Ng_MeshVS_DataSource3D) &prismaticMeshDS3D);

    //! -------------------------
    //! check self intersections
    //! -------------------------
    bool isSelfIntersecting(const occHandle(Ng_MeshVS_DataSourceFace) &aMesh, std::vector<int> &badPairs);

    //! ------------------------------
    //! intersect a mesh with another
    //! ------------------------------
    bool areIntersecting(const occHandle(Ng_MeshVS_DataSourceFace) &aMesh1,
                         const occHandle(Ng_MeshVS_DataSourceFace) &amesh2,
                         std::vector<int> &badNodes1,
                         std::vector<int> &badNodes2);

    //! ---------------------------
    //! set the progress indicator
    //! ---------------------------
    void setProgressIndicator (QProgressIndicator *aProgressIndicator)
    {
        if(aProgressIndicator!=Q_NULLPTR) myProgressIndicator = aProgressIndicator;
    }

    //! --------------------
    //! generateMeshAtWalls
    //! --------------------
    bool generateMeshAtWalls(occHandle(Ng_MeshVS_DataSource3D) &meshAtWalls,
                             occHandle(Ng_MeshVS_DataSourceFace) &lastInflatedMesh,
                             QList<occHandle(Ng_MeshVS_DataSourceFace)> &listOfInflatedMeshes);

    //! ------------------
    //! generateTetLayers
    //! ------------------
    void generateTetLayers(opencascade::handle<Ng_MeshVS_DataSource3D> &meshAtWalls, occHandle(Ng_MeshVS_DataSourceFace) &lastInflatedMesh);

private:

    //! -----
    //! task
    //! -----
    QString myTask;

    //! ------------------
    //! the mesh database
    //! ------------------
    meshDataBase *myMeshDB;

    //! -----------
    //! body index
    //! -----------
    int myBodyIndex;

    //! ---------------------------------------------------
    //! prismatic faces: list of prismatic faces on a body
    //! key: body index; value: list of prismatic faces
    //! ---------------------------------------------------
    //QList<int> myPrismaticFaces;
    std::vector<int> myPrismaticFaces;

    //! -----------
    //! parameters
    //! -----------
    prismaticLayer_sizing myTypeOfSizing;
    int myNbLayers;
    double myTotalThickness;
    double myFirstLayerThickness;
    double myExpRatio;

    //! --------------
    //! to be removed
    //! --------------
    //int myNumberOfModulationDiffusionSteps;             //1
    //double myModulationDiffusionCutoff;                 //2
    //double myModulationCoefficientTransferPercentage;   //3
    //! ---------------------
    //! end of to be removed
    //! ---------------------

    //! ---------------
    //! new parameters
    //! ---------------
    double myCurvatureSensitivityForShrink;
    int myNbGuidingVectorSmoothingSteps;
    int myNbLayerThicknessSmoothingSteps;
    double myCurvatureSensitivityForGuidingVectorsSmoothing;
    double myCurvatureSensitivityForThicknessSmoothing;
    //! ----------------------
    //! end of new parameters
    //! ----------------------

    bool myLockBoundary;
    bool myCheckSelfIntersections;
    bool myCheckMutualIntersections;

    //! --------------
    //! to be removed
    //! --------------
    //int myShrinkFunction;               //4
    //double myMinimumShrink;             //5
    //double myTransition;                //6
    //double myAmplitude;                 //7
    //! ---------------------
    //! end of to be removed
    //! ---------------------

    int myAlgorithm;
    int myBoundaryMeshType;

    //! -----------------------
    //! layer thickness values
    //! -----------------------
    std::vector<double> myLayerThickness;
    std::vector<double> myTotalThicknessAtLayer;

    //! ------------------------------
    //! local layer height modulation
    //! ------------------------------
    QMap<int,double> myLayerHCutOff;

    //! ----------------------
    //! the mesh data sources
    //! ----------------------
    occHandle(Ng_MeshVS_DataSourceFace) myOverallSumMeshDS;
    occHandle(Ng_MeshVS_DataSourceFace) myPrismaticFacesSumMeshDS;
    occHandle(Ng_MeshVS_DataSourceFace) myNonPrismaticFacesSumMeshDS;

    //! ----------------------------------------------
    //! compute the modulation of the vectorial field
    //! ----------------------------------------------
    void computeVecFieldCutOff(bool lockBoundary);

    //! ----------------------------------------
    //! compute the shrink factor for each node
    //! ----------------------------------------
    void computeShrinkFactor(const occHandle(Ng_MeshVS_DataSourceFace) &aMeshVS_DataSource,
                             QMap<int,double> &shrinkFactors);

    //! ------------------------
    //! check self intersection
    //! ------------------------
    bool checkSelfIntersection(const occHandle(Ng_MeshVS_DataSourceFace) &aMeshDS,
                               const occHandle(Ng_MeshVS_DataSourceFace) &aMeshDS_old,
                               bool correct=true);

    //! --------------------------
    //! check mutual intersection
    //! --------------------------
    bool checkMutualIntersection(const occHandle(Ng_MeshVS_DataSourceFace) &aMeshDS,
                                 const occHandle(Ng_MeshVS_DataSourceFace) &aMeshDS_old,
                                 bool correct=true);

    //! -------------------
    //! progress indicator
    //! -------------------
    QProgressIndicator *myProgressIndicator;

    //! -----------
    //! invert tet
    //! -----------
    void invertTet(meshElementByCoords &aTet);

private:

    //! -----------------------------------------------------
    //! map of manifold characteristic and visibility angles
    //! -----------------------------------------------------
    std::map<int,double> betaAverageField;
    std::map<int,double> betaVisibilityField;

    //! --------------------
    //! generateOneTetLayer
    //! --------------------
    void generateOneTetLayer(occHandle(Ng_MeshVS_DataSourceFace) &theMeshToInflate,
                             double displacement,
                             std::vector<meshElementByCoords> &volumeElementsAtWalls);

    //! -------------
    //! compute beta
    //! -------------
    void computeBeta(const occHandle(Ng_MeshVS_DataSourceFace) &aMeshDS);

    //! ---------------------------------------------
    //! check lateral distribution marhing distances
    //! ---------------------------------------------
    void checkLateralDistributionMarchingDistance(const occHandle(Ng_MeshVS_DataSourceFace) &aMeshDS,
                                                  const std::map<int,double> displacementMapOld,
                                                  std::vector<int,double> &displacementMap);

    //! ---------------
    //! local manifold
    //! ---------------
    bool getLocalFanNodes(const occHandle(Ng_MeshVS_DataSourceFace) &aMeshDS, int vertexGlobalNodeID,
                          int t1, int t2,
                          mesh::meshPoint &A,
                          mesh::meshPoint &B,
                          mesh::meshPoint &C,
                          mesh::meshPoint &P);

    //! ---------------
    //! classify nodes
    //! ---------------
    void classifyNodes(const occHandle(Ng_MeshVS_DataSourceFace) &aMeshDS,
                       std::map<int,int> &mapCat1);

    //! --------------------------------------------------
    //! enable/disable the progress indicator stop button
    //! --------------------------------------------------
    void setStopButtonEnabled(bool isEnabled);
};

#endif // PRISMATICLAYER_H
