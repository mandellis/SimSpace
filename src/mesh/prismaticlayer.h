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

    //! --------------------
    //! set prismatic faces
    //! --------------------
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
    std::vector<int> myPrismaticFaces;

    //! -----------
    //! parameters
    //! -----------
    prismaticLayer_sizing myTypeOfSizing;
    int myNbLayers;
    double myTotalThickness;
    double myFirstLayerThickness;
    double myExpRatio;

    bool myLockBoundary;
    bool myCheckSelfIntersections;
    bool myCheckMutualIntersections;

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
    //! surrounding nodes map
    //! ----------------------
    std::map<int,std::vector<int>> mySurroudingNodesMap;

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

    //! -------------------------------------
    //! check tetrahedral front intersection
    //! -------------------------------------
    void checkTetFrontIntersections(const occHandle(Ng_MeshVS_DataSourceFace) &aMeshDS,
                                    const QMap<int,QList<double>> &guidingVectors,
                                    QMap<int,QList<double>> &displacementMap);

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

    //! ------------------------------
    //! map of layer reduction factor
    //! ------------------------------
    std::map<int,double> mapOfReductionFactor;

    //! -------------
    //! compute beta
    //! -------------
    void computeBeta(const occHandle(Ng_MeshVS_DataSourceFace) &aMeshDS);

    //! ---------------------------------------------
    //! check lateral distribution marhing distances
    //! ---------------------------------------------
    void checkLateralDistributionMarchingDistance(const occHandle(Ng_MeshVS_DataSourceFace) &aMeshDS,
                                                  int inflationStep,
                                                  QMap<int,double> &marchingDistanceMap,
                                                  const QMap<int,QList<double>> &guidingDirections);

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

    //! -------------------------
    //! buildSurroundingNodesMap
    //! -------------------------
    void buildSurroundingNodesMap(const occHandle(Ng_MeshVS_DataSourceFace) &aMeshDS, std::map<int,std::vector<int>> &surroundingNodesMap);

    //! --------------
    //! analyze front
    //! --------------
    void analyzeFront(const occHandle(Ng_MeshVS_DataSourceFace) &theMeshToInflate,
                      const QMap<int,QList<double>> &normals,
                      QMap<int,double> &marchingDistancesMap);

    //! ---------------------------
    //! check incomplete manifolds
    //! ---------------------------
    void checkIncompleteManifold(const occHandle(Ng_MeshVS_DataSourceFace) &theMeshToInflate);


    void processMesh(occHandle(Ng_MeshVS_DataSourceFace) &aMeshDS);

    //! --------------------------------------------------
    //! enable/disable the progress indicator stop button
    //! --------------------------------------------------
    void setStopButtonEnabled(bool isEnabled);

    //! -------
    //! helper
    //! -------
    bool getPointCoordinates(const occHandle(MeshVS_DataSource)& aMeshDS, int globalNodeID, double *P)
    {
        int NbNodes;
        MeshVS_EntityType aType;
        double buf[3];
        TColStd_Array1OfReal coords(*buf,1,3);
        bool isDone = aMeshDS->GetGeom(globalNodeID,false,coords,NbNodes,aType);
        P[0] = coords(1);
        P[1] = coords(2);
        P[2] = coords(3);
        if(isDone == false) return false;
        return isDone;
    }

    //! ------------------------------------------------------
    //! rotate vector v about direction k by an angle "angle"
    //! ------------------------------------------------------
    void rotateVector(double *v, double *k, double angle, double *r)
    {
        //! ---------------------
        //! i       j       k
        //! k[0]    k[1]    k[2]
        //! v[0]    v[1]    v[2]
        //! ---------------------
        double sx = k[1]*v[2]-k[2]*v[1];
        double sy = k[2]*v[0]-k[0]*v[2];
        double sz = k[0]*v[1]-k[1]*v[0];

        double dot = v[0]*k[0]+v[1]*k[1]+v[2]*k[2];
        if(dot>1) dot = 1;
        if(dot<-1) dot = -1;
        r[0] = v[0]*cos(angle)+ sx*sin(angle) + k[0]*dot*(1-cos(angle));
        r[1] = v[1]*cos(angle)+ sy*sin(angle) + k[1]*dot*(1-cos(angle));
        r[2] = v[2]*cos(angle)+ sz*sin(angle) + k[2]*dot*(1-cos(angle));
    }
};

#endif // PRISMATICLAYER_H
