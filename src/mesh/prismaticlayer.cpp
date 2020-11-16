//! ----------------
//! custom includes
//! ----------------
#include "prismaticlayer.h"
#include "meshtools.h"
#include "ccout.h"
#include "geomtoolsclass.h"
#include "igtools.h"
#include <ng_meshvs_datasourceface.h>
#include "meshtools.h"
#include "global.h"
#include <tetqualityclass.h>
#include <meshface.h>
#include <OCCface.h>
#include <meshselectnlayers.h>
#include <igtools.h>
#include <smoothingtools.h>
#include <pointtomeshdistance.h>

//! ------
//! Eigen
//! ------
#include <Eigen/Dense>

//! ----------------------
//! used by "displayMesh"
//! ----------------------
#include <occPreGLwidget.h>
#include <AIS_InteractiveContext.hxx>
#include <MeshVS_MeshPrsBuilder.hxx>
#include <MeshVS_Drawer.hxx>
#include <MeshVS_DrawerAttribute.hxx>
#include "arrayofcolors.h"
#include "qprogressevent.h"

//! ---
//! Qt
//! ---
#include <QApplication>

//! ----
//! C++
//! ----
#include <iostream>
#include <memory>
#include <algorithm>

//! ----
//! OCC
//! ----
#include <TopoDS_Face.hxx>
#include <TopoDS.hxx>

//! --------------------------------------------------------
//! function: constructor
//! details:  constructor with default inflation parameters
//! --------------------------------------------------------
prismaticLayer::prismaticLayer(meshDataBase *mDB, QProgressIndicator *aProgressIndicator)
{
    //! ------------------
    //! the mesh database
    //! ------------------
    myMeshDB = mDB;

    //! -----------------------------------------------------
    //! prismatic layers parameters: a set of default values
    //! -----------------------------------------------------
    myTypeOfSizing = prismaticLayer_sizing_FirstLayerThickness;
    myExpRatio = 1.0;
    myNbLayers = 1;
    myTotalThickness = 1;
    myFirstLayerThickness = 1;

    myLockBoundary = true;
    myAlgorithm = 0;            //! generation algorithm

    myProgressIndicator = aProgressIndicator;
    myTask = QString("Building prismatic mesh");
}

//! ----------------------------
//! function: setPrismaticFaces
//! details:
//! ----------------------------
void prismaticLayer::setPrismaticFaces(const std::vector<int> &prismaticFaces)
{
    for(size_t i=0; i<prismaticFaces.size(); i++) myPrismaticFaces.push_back(prismaticFaces.at(i));
}

//! -----------------------------------------------------------------------------------
//! Set the parameters for the prismatic layer generation
//! n parameters are defined:
//! - the thickness of the first layer
//! - the number of layers
//! - the expansion ratio
//! - the total thickness of the prismatic layer mesh
//! One can define the the behavior of the prismatic layer mesh in two ways:
//! 1) At fixed expansion ratio, by entering the first layer thickness, and the number
//!    of thickness. The total thickness of the prismatic mesh is a consequence
//! 2) At fixed expansion ratio, by entering the total thickness, and the number
//!    of layers. The thickness of the first layer is a consequence
//! -----------------------------------------------------------------------------------
void prismaticLayer::setParameters(prismaticLayerParameters parameters)
{
    prismaticLayer_sizing type = parameters.typeOfSizing;
    switch(type)
    {
    //! ------------------------
    //! "First layer thickness"
    //! ------------------------
    case prismaticLayer_sizing_FirstLayerThickness:
    {
        myFirstLayerThickness = parameters.firstLayerThickness;
        myExpRatio = parameters.expRatio;
        myNbLayers = parameters.NbLayers;

        myLayerThickness.push_back(parameters.firstLayerThickness);
        myTotalThicknessAtLayer.push_back(parameters.firstLayerThickness);

        //! -------------------------------------
        //! a consequence: total layer thickness
        //! fill also the vector of layer sizes
        //! -------------------------------------
        double layerThickness_old = myFirstLayerThickness;
        myTotalThickness = myFirstLayerThickness;

        if(myNbLayers>1)
        {
            myTotalThickness = myFirstLayerThickness;

            //! -----------------------------
            //! compute each layer thickness
            //! -----------------------------
            for(int i=1; i<myNbLayers; i++)
            {
                double layerThickness =  layerThickness_old*myExpRatio;
                myLayerThickness.push_back(layerThickness);

                myTotalThickness = myTotalThickness + layerThickness;
                myTotalThicknessAtLayer.push_back(myTotalThickness);
                layerThickness_old = layerThickness;
            }
        }
    }
        break;

    //! ------------------
    //! "Total thickness"
    //! ------------------
    case prismaticLayer_sizing_TotalThickness:
    {
        myNbLayers = parameters.NbLayers;
        myExpRatio = parameters.expRatio;
        myTotalThickness = parameters.totalThickness;

        //! -------------------------------------
        //! a consequence: first layer thickness
        //! and total thickness
        //! -------------------------------------
        myFirstLayerThickness = myTotalThickness*(1-myExpRatio)/(1-pow(myExpRatio,myNbLayers));
        myTotalThicknessAtLayer.push_back(myFirstLayerThickness);

        double layerThickness_old = myFirstLayerThickness;
        myLayerThickness.push_back(myFirstLayerThickness);
        if(myNbLayers>1)
        {
            //! -----------------------------
            //! compute each layer thickness
            //! -----------------------------
            for(int i=1; i<myNbLayers; i++)
            {
                double layerThickness = layerThickness_old*myExpRatio;
                myTotalThicknessAtLayer.push_back(layerThickness+layerThickness_old);
                myLayerThickness.push_back(layerThickness);
                layerThickness_old = layerThickness;
            }
        }
    }
        break;
    }

    myLockBoundary = parameters.lockBoundary;
    myCheckSelfIntersections = parameters.checkSelfIntersections;
    myCheckMutualIntersections = parameters.checkMutualIntersections;

    myAlgorithm = parameters.generationAlgorithm;
    myBoundaryMeshType = parameters.boundaryMeshType;

    //! -----------------
    //! console messages
    //! -----------------
    cout<<"\n@____Prismatic layers parameters____"<<endl;
    cout<<"@____Number of layers: "<<myNbLayers<<"____"<<endl;
    cout<<"@____Expansion ratio: "<<myExpRatio<<"____"<<endl;

    cout<<"@____Prismatic layer: list of heights____"<<endl;
    for(int i=0; i<myLayerThickness.size(); i++)
        cout<<"@____layer nr: "<<i<<"  thickness = "<<myLayerThickness.at(i)<<"; Total thickness @ layer: "<<myTotalThicknessAtLayer.at(i)<<"____"<<endl;

    cout<<"@____Total thickness: "<<myTotalThickness<<"____"<<endl;
    cout<<"@____Boundary behavior "<<(myLockBoundary==true? "Locked":"Free")<<"____"<<endl;
    cout<<"@____Algorithm: "<<myAlgorithm<<"____@"<<endl;
}

//! -----------------------------------------------------------
//! function: computeVecFieldCutOff
//! details:  works with global node IDs
//!
//!           prismatic          prismatic
//!            ______              ____          <- cutoff 1.0
//!           |      |            |    |
//!      _____|      |____________|    |______   <- cutoff 0.0
//!
//! this operation can be performed one sigle time, since
//! the displaced meshes share the same global node IDs
//! -----------------------------------------------------------
void prismaticLayer::computeVecFieldCutOff(bool lockBoundary)
{
    myLayerHCutOff.clear();
    myPrismaticFacesSumMeshDS->computeFreeMeshSegments();

    //! --------------------------------------------
    //! the following list contains global node IDs
    //! --------------------------------------------
    QList<int> nodeIDsOnPrismaticFaces1DBoundary = myPrismaticFacesSumMeshDS->myBoundaryPoints;

    cout<<"____number of points on the prismatic faces boundary: "<<nodeIDsOnPrismaticFaces1DBoundary.length()<<"____"<<endl;

    //! ----------------------------
    //! diagnostic - can be removed
    //! ----------------------------
    double cutoff = 0.0;

    QMap<int,QList<double>> allNodesNormals = myOverallSumMeshDS->myNodeNormals;
    if(allNodesNormals.isEmpty()) cerr<<"____normals at nodes for the overall mesh not defined____"<<endl;

    QMap<int,QList<double>> prismaticNodeNormals = myPrismaticFacesSumMeshDS->myNodeNormals;
    if(prismaticNodeNormals.isEmpty()) cerr<<"____normals at prismatic nodes not defined____"<<endl;

    for(QMap<int,QList<double>>::iterator it = allNodesNormals.begin(); it!=allNodesNormals.end(); ++it)
    {
        int nodeID = it.key();
        if(prismaticNodeNormals.contains(nodeID)) cutoff = 1.0;
        else cutoff = 0.0;

        if(nodeIDsOnPrismaticFaces1DBoundary.contains(nodeID))
        {
            if(lockBoundary==true) cutoff = 0.0;
            else cutoff = 1.0;
        }
        //! -------
        //! cutoff
        //! -------
        myLayerHCutOff.insert(nodeID, cutoff);
    }
}

//! ---------------------------------------------------
//! function: inflateMesh
//! details:  returns the list of the displaced meshes
//! ---------------------------------------------------
bool prismaticLayer::inflateMesh(QList<occHandle(Ng_MeshVS_DataSourceFace)> &inflatedMeshes)
{
    //! -----------------
    //! a progress event
    //! -----------------
    QProgressEvent *progressEvent;

    //! --------------------------------
    //! init the secondary progress bar
    //! --------------------------------
    if(myProgressIndicator!=Q_NULLPTR)
    {
        progressEvent = new QProgressEvent(QProgressEvent_None,0,9999,0,"Start inflating mesh",
                                                QProgressEvent_Init,0,myNbLayers-1,0,"Building prismatic mesh");
        QApplication::postEvent(myProgressIndicator,progressEvent);
        QApplication::processEvents();
    }

    double displacement;

    //! -----------------------------------------------------
    //! set zero displacement on mesh free edge, if required
    //! if "true" the points of the 1D boundary are locked
    //! -----------------------------------------------------
    bool lockBoundary = myLockBoundary;

    //! ------------------------
    //! disable the Stop button
    //! ------------------------
    if(myProgressIndicator!=Q_NULLPTR) this->setStopButtonEnabled(false);

    //! --------------------------------------
    //! decide which node should be displaced
    //! --------------------------------------
    this->computeVecFieldCutOff(lockBoundary);

    //! ---------------------------
    //! make a copy the outer mesh
    //! ---------------------------
    occHandle(Ng_MeshVS_DataSourceFace) theMeshToInflate_old = new Ng_MeshVS_DataSourceFace(myOverallSumMeshDS);
    theMeshToInflate_old->computeNormalAtNodes();

    //! -----------------------------------
    //! generate the surrounding nodes map
    //! -----------------------------------
    mySurroudingNodesMap.clear();
    this->buildSurroundingNodesMap(theMeshToInflate_old, mySurroudingNodesMap);

    //! -----------------------
    //! enable the Stop button
    //! -----------------------
    if(myProgressIndicator!=Q_NULLPTR) myProgressIndicator->enableStop();

    //! -----------------------------------------------------------
    //! insert the outer mesh into the list of the inflated meshes
    //! -----------------------------------------------------------
    inflatedMeshes<<theMeshToInflate_old;

    //! -------------------------------------------
    //! clear the fist layer reduction factors map
    //! -------------------------------------------
    mapOfReductionFactor.clear();

    //! ------------------------------
    //! analyze gaps at the beginning
    //! ------------------------------
    double Sigma = 0;                                           // a parameter for reduction factor calculation
    for(int n=0; n<myNbLayers; n++) Sigma += pow(myExpRatio,n);
    pointToMeshDistance aDistanceMeter;                         // distance meter tool
    aDistanceMeter.init(theMeshToInflate_old);                  // init with mesh
    double firstLayerThickness = myLayerThickness.at(0);        // first layer thickness
    std::map<int,double> mapOfFirstLayerReductionFactor;        // fill a map (nodeID, reduction factor)

    //! ---------------------------------------------
    //! start filling the map of compression factors
    //! ---------------------------------------------
    const QMap<int,QList<double>> &normals = theMeshToInflate_old->myNodeNormals;
    for(QMap<int,QList<double>>::const_iterator itn = normals.cbegin(); itn!=normals.cend(); itn++)
    {
        int globalNodeID = itn.key();
        const QList<double> &normal = itn.value();

        double P[3];
        float distance;
        this->getPointCoordinates(theMeshToInflate_old,globalNodeID,P);

        //! --------------------------------------------------
        //! invert the normal: it points towards the material
        //! for gap calculation
        //! --------------------------------------------------
        double nx = -normal[0];
        double ny = -normal[1];
        double nz = -normal[2];
        double dir[3] {nx,ny,nz};
        aDistanceMeter.distance(P,dir,&distance);   // compute the distance

        bool compress = true;
        if(compress==true)
        {
            //! ---------------------------------------------
            //! compress layers in order to avoid collisions
            //! ---------------------------------------------
            const double compressionFactor = 0.25;
            if(myTotalThickness<distance/3.0)
            {
                //! the node can move with the nominal marching distance
                mapOfFirstLayerReductionFactor.insert(std::make_pair(globalNodeID,1.0));
            }
            else
            {
                //! the node can move with smaller marching distance
                double firstLayerThicknessNew = compressionFactor*distance/Sigma;
                double reductionFactor = firstLayerThicknessNew/firstLayerThickness;
                mapOfFirstLayerReductionFactor.insert(std::make_pair(globalNodeID,reductionFactor));
                mapOfReductionFactor.insert(std::make_pair(globalNodeID,reductionFactor));
            }
        }
        else
        {
            if(myTotalThickness<distance/3.0)
            {
                //! the node can move without compression
                mapOfFirstLayerReductionFactor.insert(std::make_pair(globalNodeID,1.0));
            }
            else
            {
                //! lock point - the node cannot move
                myLayerHCutOff.insert(globalNodeID,0);
                mapOfFirstLayerReductionFactor.insert(std::make_pair(globalNodeID,0.0));
                mapOfReductionFactor.insert(std::make_pair(globalNodeID,0.0));
            }
        }
    }

    for(int n=1; n<=myNbLayers; n++)
    {
        cout<<"prismaticLayer::inflateMesh()->____generating layer: "<<(n-1)<<" - current thickness: "<<myLayerThickness.at(n-1)<<"____"<<endl;

        if(Global::status().code = 0)
        {
            cout<<"\\ ------------------------------------------\\"<<endl;
            cout<<"\\ inflation process interrupted by the user \\"<<endl;
            cout<<"\\ ------------------------------------------\\"<<endl;
            return false;
        }

        //! ----------------------------
        //! the current layer thickness
        //! ----------------------------
        displacement = myLayerThickness.at(n-1);

        //! ---------------------
        //! "best" nodal normals
        //! ---------------------
        occHandle(Ng_MeshVS_DataSourceFace) theMeshToInflate_new = new Ng_MeshVS_DataSourceFace(theMeshToInflate_old);
        theMeshToInflate_new->computeNormalAtNodes();
        this->computeBeta(theMeshToInflate_new);

        //! ---------------
        //! shrink factors
        //! ---------------
        QMap<int,double> shrinkFactors;
        this->computeShrinkFactor(theMeshToInflate_new,shrinkFactors);

        //! ---------------------------------
        //! guiding vectors - "best normals"
        //! ---------------------------------
        QMap<int,QList<double>> normals = theMeshToInflate_new->myNodeNormals;

        //! ---------------
        //! classify nodes
        //! ---------------
        std::map<int,int> mapNodeTypes;
        this->classifyNodes(theMeshToInflate_new,mapNodeTypes);

        //! -----------------------------------------------------------------------------------------------
        //! smooth the guiding vectors directions - use 5/10 smoothing steps
        //! flag normalize == true since the "rotation" of the vector is of interest (change in direction)
        //! -----------------------------------------------------------------------------------------------
        smoothingTools::fieldSmoother(normals,theMeshToInflate_new,betaAverageField,betaVisibilityField,mapNodeTypes,10,true);

        //! ------------------------------------
        //! map of the local marching distances
        //! ------------------------------------
        QMap<int,double> marchingDistanceMap;
        for(QMap<int,QList<double>>::const_iterator it = normals.cbegin(); it!= normals.cend(); ++it)
        {
            int globalNodeID = it.key();
            double shrinkFactor = shrinkFactors.value(globalNodeID);
            double cutOff = myLayerHCutOff.value(globalNodeID);
            double marchingDistance = displacement*(1+shrinkFactor)*cutOff*mapOfFirstLayerReductionFactor.at(globalNodeID);
            marchingDistanceMap.insert(globalNodeID,marchingDistance);
        }

        //! cesere

        //! ------------------------------------------------
        //! check lateral distribution of marhing distances
        //! ------------------------------------------------
        //this->checkLateralDistributionMarchingDistance(theMeshToInflate_new,marchingDistanceMap);

        //! ------------------------------
        //! smooth the marching distances
        //! ------------------------------
        smoothingTools::scalarFieldSmoother(marchingDistanceMap,theMeshToInflate_new,betaAverageField,mapNodeTypes,10);

        //! --------------------------------
        //! generate the displacement field
        //! --------------------------------
        QMap<int,QList<double>> displacementsField;
        for(QMap<int,QList<double>>::const_iterator it = normals.cbegin(); it!= normals.cend(); ++it)
        {
            int globalNodeID = it.key();
            const QList<double> &curNormal = it.value();

            //! --------------------
            //! nodal displacements
            //! --------------------
            double marchingDistance = marchingDistanceMap.value(globalNodeID);
            double vx = -curNormal[0]*marchingDistance;
            double vy = -curNormal[1]*marchingDistance;
            double vz = -curNormal[2]*marchingDistance;

            QList<double> localFieldValue;
            localFieldValue<<vx<<vy<<vz;
            displacementsField.insert(globalNodeID,localFieldValue);
        }

        //! ------------------
        //! displace the mesh
        //! ------------------
        theMeshToInflate_new->displaceMySelf(displacementsField);

        //! ------------------------
        //! check self intersection
        //! ------------------------
        if(myCheckSelfIntersections==true) this->checkSelfIntersection(theMeshToInflate_new,theMeshToInflate_old,true);

        //! ---------------------------------------------------------
        //! start correcting nodes coordinates if the inflated layer
        //! has an intersection with the mesh it originates from
        //! ---------------------------------------------------------
        if(myCheckMutualIntersections==true) this->checkMutualIntersection(theMeshToInflate_new,theMeshToInflate_old,true);

        //! ----------------------------------------------
        //! add the displaced mesh to the list of results
        //! ----------------------------------------------
        inflatedMeshes<<theMeshToInflate_new;

        //! --------------------------------------------
        //! check "Cliff points" - incomplete manifolds
        //! --------------------------------------------
        this->checkIncompleteManifold(theMeshToInflate_old);

        //! --------------------------------------
        //! replace the old mesh with the new one
        //! --------------------------------------
        theMeshToInflate_old = new Ng_MeshVS_DataSourceFace(theMeshToInflate_new);

        //! ---------------------
        //! post an update event
        //! ---------------------
        if(myProgressIndicator!=Q_NULLPTR)
        {
            progressEvent = new QProgressEvent(QProgressEvent_None,-1,-1,-1,QString("Inflation step %1 done").arg(n),
                                    QProgressEvent_Update,-1,-1,n,"Building prismatic mesh");
            QApplication::postEvent(myProgressIndicator,progressEvent);
            QApplication::processEvents();
        }
    }
    //! ---------------------
    //! post an update event
    //! ---------------------
    if(myProgressIndicator!=Q_NULLPTR)
    {
        progressEvent = new QProgressEvent(QProgressEvent_None,-1,-1,-1,"Inflation finished",
                                           QProgressEvent_Reset,0,100,0,"Building prismatic mesh");
        QApplication::postEvent(myProgressIndicator,progressEvent);
        QApplication::processEvents();
    }
    cout<<"prismaticLayer::inflateMesh()->_____inflation finished____"<<endl;
    return true;
}

//! --------------------------
//! function: mergeFaceMeshes
//! details:
//! --------------------------
bool prismaticLayer::mergeFaceMeshes()
{
    cout<<"prismaticLayer::mergeFaceMeshes()->____function called____"<<endl;

    int bodyIndex = myBodyIndex;
    int NbGeometryFaces = myMeshDB->MapOfBodyTopologyMap.value(bodyIndex).faceMap.Extent();

    //! ---------------------------------
    //! faceList0 => any type of face
    //! faceList1 => prismatic faces
    //! faceList2 => non prismatic faces
    //! ---------------------------------
    QList<occHandle(Ng_MeshVS_DataSourceFace)> faceList0, faceList1, faceList2;
    for(int faceNr = 1; faceNr <= NbGeometryFaces; faceNr++)
    {
        const occHandle(Ng_MeshVS_DataSourceFace) &aFaceMeshDS =
                occHandle(Ng_MeshVS_DataSourceFace)::DownCast(myMeshDB->ArrayOfMeshDSOnFaces.getValue(bodyIndex,faceNr));

        faceList0<<aFaceMeshDS;
        if(std::find(myPrismaticFaces.begin(), myPrismaticFaces.end(), faceNr)!=myPrismaticFaces.end()) faceList1<<aFaceMeshDS;
        if(std::find(myPrismaticFaces.begin(), myPrismaticFaces.end(), faceNr)==myPrismaticFaces.end()) faceList2<<aFaceMeshDS;
    }

    cout<<"____faceList0: "<<faceList0.length()<<"____"<<endl;
    cout<<"____faceList1: "<<faceList1.length()<<"____"<<endl;
    cout<<"____faceList2: "<<faceList2.length()<<"____"<<endl;

    myOverallSumMeshDS = new Ng_MeshVS_DataSourceFace(faceList0);
    myPrismaticFacesSumMeshDS = new Ng_MeshVS_DataSourceFace(faceList1);

    myNonPrismaticFacesSumMeshDS = new Ng_MeshVS_DataSourceFace(faceList2);
    //myOverallSumMeshDS->computeNormalAtElements();
    myOverallSumMeshDS->computeNormalAtNodes();
    //myPrismaticFacesSumMeshDS->computeNormalAtElements();
    myPrismaticFacesSumMeshDS->computeNormalAtNodes();
    //myNonPrismaticFacesSumMeshDS->computeNormalAtElements();
    myNonPrismaticFacesSumMeshDS->computeNormalAtNodes();

    return true;
}

//! -----------------------------
//! function: angleBetweenVector
//! details:
//! -----------------------------
double angleBetweenVector(const std::vector<double> &n1, const std::vector<double> &n2, const std::vector<double> &refdir)
{
    double l1 = sqrt(pow(n1[0],2)+pow(n1[1],2)+pow(n1[2],2));
    double l2 = sqrt(pow(n2[0],2)+pow(n2[1],2)+pow(n2[2],2));

    if(l1 * l2 == 0) exit(9999);    // questo non dovrebbe mai succedere

    double dot = (n1[0]*n2[0]+n1[1]*n2[1]+n1[2]*n2[2])/(l1*l2);
    if(dot > 1) dot = 1;        // avoids domain error
    if(dot < -1) dot = -1;      // avoids domain error
    double angle = std::acos(dot);
    double cx = (n1[1]*n2[2]-n1[2]*n2[1]);      //  n1[0]   n1[1]   n1[2]
    double cy = (n1[2]*n2[0]-n1[0]*n2[2]);      //  n2[0]   n2[1]   n2[2]
    double cz = (n1[0]*n2[1]-n1[1]*n2[0]);

    double l_cross = sqrt(cx*cx+cy*cy+cz*cz);
    cx /= l_cross;
    cy /= l_cross;
    cz /= l_cross;

    double v = refdir[0]*cx+refdir[1]*cy+refdir[2]*cz;
    if(v<0.0) return -angle;
    else return angle;
}

//! ------------------------------------
//! check if two elements share an edge
//! ------------------------------------
bool areAdjacent(const std::vector<int> &aFaceElement1, const std::vector<int> &aFaceElement2)
{
    std::set<int> aSet;
    for(std::vector<int>::const_iterator it = aFaceElement1.cbegin(); it!=aFaceElement1.cend(); it++) aSet.insert(*it);
    for(std::vector<int>::const_iterator it = aFaceElement2.cbegin(); it!=aFaceElement2.cend(); it++) aSet.insert(*it);
    if(aFaceElement1.size()+aFaceElement2.size()-aSet.size()==2) return true;
    return false;
}

//! ---------------------------
//! function: getLocalFanNodes
//! details:
//! ---------------------------
bool prismaticLayer::getLocalFanNodes(const occHandle(Ng_MeshVS_DataSourceFace) &aMeshDS, int vertexGlobalNodeID,
                                      int t1, int t2,
                                      mesh::meshPoint &A,
                                      mesh::meshPoint &B,
                                      mesh::meshPoint &C,
                                      mesh::meshPoint &P)
{
    //cout<<"prismaticLayer::getLocalFanNodes()->____function called____"<<endl;

    std::vector<int> supportVector;
    std::vector<int> commonEdge;

    int buf[3], NbNodes;
    TColStd_Array1OfInteger nodeIDs1(*buf,1,3);         // first triangle
    aMeshDS->GetNodesByElement(t1,nodeIDs1,NbNodes);
    for(int i=1; i<=NbNodes; i++)
    {
        int currentNodeID = nodeIDs1(i);
        supportVector.push_back(currentNodeID); // contains A,P,C
        //cout<<"____first wedge triangle: "<<currentNodeID<<"____"<<endl;
    }
    TColStd_Array1OfInteger nodeIDs2(*buf,1,3);         // second triangle
    aMeshDS->GetNodesByElement(t2,nodeIDs2,NbNodes);
    for(int i=1; i<=NbNodes; i++)
    {
        int currentNodeID = nodeIDs2(i);
        //cout<<"____second wedge triangle: "<<currentNodeID<<"____"<<endl;
        std::vector<int>::iterator it = std::find(supportVector.begin(), supportVector.end(), currentNodeID);
        if(it == supportVector.end()) supportVector.push_back(currentNodeID);
        else commonEdge.push_back(currentNodeID);
    }

    //cout<<"____common edge ("<<commonEdge[0]<<", "<<commonEdge[1]<<")____"<<endl;
    std::vector<int>::iterator it_ = std::find(supportVector.begin(), supportVector.end(), commonEdge[0]);
    supportVector.erase(it_);
    it_ = std::find(supportVector.begin(), supportVector.end(), commonEdge[1]);
    supportVector.erase(it_);

    // now the support vector contains A and B
    // here we are not sure if A and B are correct
    // they could be swapped - see at the end
    int globalNodeID_A = supportVector[0];
    int globalNodeID_B = supportVector[1];

    //cout<<"____A: "<<globalNodeID_A<<", B: "<<globalNodeID_B<<"____"<<endl;

    MeshVS_EntityType aType;
    int NbNodes_;
    double bufd[3];
    TColStd_Array1OfReal coordsNodeA(*bufd,1,3);
    TColStd_Array1OfReal coordsNodeB(*bufd,1,3);
    TColStd_Array1OfReal coordsNodeP(*bufd,1,3);

    //! --------------------------
    //! point A, point B, point P
    //! --------------------------
    aMeshDS->GetGeom(globalNodeID_A,false,coordsNodeA,NbNodes_,aType);
    A.ID = globalNodeID_A;
    A.x = coordsNodeA(1);
    A.y = coordsNodeA(2);
    A.z = coordsNodeA(3);
    //cout<<"____A("<<A.x<<", "<<A.y<<", "<<A.z<<")____"<<endl;

    aMeshDS->GetGeom(globalNodeID_B,false,coordsNodeB,NbNodes_,aType);
    B.ID = globalNodeID_B;
    B.x = coordsNodeB(1);
    B.y = coordsNodeB(2);
    B.z = coordsNodeB(3);
    //cout<<"____B("<<B.x<<", "<<B.y<<", "<<B.z<<")____"<<endl;

    aMeshDS->GetGeom(vertexGlobalNodeID,false,coordsNodeP,NbNodes_,aType);
    P.ID = vertexGlobalNodeID;
    P.x = coordsNodeP(1);
    P.y = coordsNodeP(2);
    P.z = coordsNodeP(3);
    //cout<<"____P("<<P.x<<", "<<P.y<<", "<<P.z<<")____"<<endl;

    //! -----------------------
    //! also point C is needed
    //! -----------------------
    int globalNodeID_C;
    if(commonEdge[0] == vertexGlobalNodeID) globalNodeID_C = commonEdge[1];
    else globalNodeID_C = commonEdge[0];
    TColStd_Array1OfReal coordsNodeC(*bufd,1,3);
    aMeshDS->GetGeom(globalNodeID_C,false,coordsNodeC,NbNodes_,aType);
    C.ID = globalNodeID_C;
    C.x = coordsNodeC(1);
    C.y = coordsNodeC(2);
    C.z = coordsNodeC(3);
    //cout<<"____C("<<C.x<<", "<<C.y<<", "<<C.z<<")____"<<endl;

    //! --------------
    //! (A-P) x (C-P)
    //! --------------
    double xAP = A.x-P.x;
    double yAP = A.y-P.y;
    double zAP = A.z-P.z;

    double xCP = C.x-P.x;
    double yCP = C.y-P.y;
    double zCP = C.z-P.z;

    //! --------------
    //! i    j    k
    //! xAP  yAP  zAP
    //! xCP  yCP  zCP
    //! --------------
    double vx = yAP*zCP-zAP*yCP;
    double vy = zAP*xCP-xAP*zCP;
    double vz = xAP*yCP-yAP*xCP;

    //! -----------------------------------------
    //! the normal of the element containing "A"
    //! -----------------------------------------
    double nx, ny, nz;
    aMeshDS->GetNormal(t1,10,nx,ny,nz);

    double scalarP = vx*nx+vy*ny+vz*nz;
    if(scalarP<0)   // swap A and B
    {
        mesh::meshPoint S = B;
        B = mesh::meshPoint(A);
        A = S;
        return true;
    }
    return false;
}

//! ----------------------
//! function: computeBeta
//! details:
//! ----------------------
void prismaticLayer::computeBeta(const occHandle(Ng_MeshVS_DataSourceFace) &aMeshDS)
{
    cout<<"prismaticLayer::computeBeta()->____function called____"<<endl;

    const double PI = 3.1415926538;

    //! ----------------------------------------
    //! check if node normals have been defined
    //! ----------------------------------------
    if(aMeshDS->myNodeNormals.isEmpty()) aMeshDS->computeNormalAtNodes();

    //! --------------------------
    //! clear the private members
    //! --------------------------
    betaAverageField.clear();
    betaVisibilityField.clear();

    //! -----------------------------------
    //! face elements attached to the node
    //! -----------------------------------
    if(aMeshDS->myNodeToElements.isEmpty()) aMeshDS->computeNodeToElementsConnectivity();
    const QMap<int,QList<int>> &nodeToElementsMap = aMeshDS->myNodeToElements;

    //! ------------------------------
    //! do the work on all the points
    //! ------------------------------
    for(TColStd_MapIteratorOfPackedMapOfInteger it(aMeshDS->GetAllNodes()); it.More(); it.Next())
    {
        int globalNodeID = it.Key();
        int localNodeID = aMeshDS->myNodesMap.FindIndex(globalNodeID);

        const QList<int> &attachedElements_localIDs = nodeToElementsMap.value(localNodeID);
        int NbAttachedElements = attachedElements_localIDs.length();

        double bnk_max = -1e10;
        //std::vector<double> bnks;
        std::set<int> surroundingNodes;                             // nodes surrounding the current node
        std::map<int,int> indexedMapOfElements;                     // index, element ID
        std::map<int,std::vector<double>> indexedMapOfNormals;      // index, normal of the element ID
        for(int i=0; i<NbAttachedElements; i++)
        {
            int localElementID = attachedElements_localIDs[i];
            int globalElementID = aMeshDS->myElementsMap.FindKey(localElementID);

            double nx,ny,nz;
            nx = ny = nz = 0;
            aMeshDS->GetNormal(globalElementID,10,nx,ny,nz);

            std::vector<double> aNormal {nx, ny, nz};
            indexedMapOfElements.insert(std::make_pair(i,globalElementID));
            indexedMapOfNormals.insert(std::make_pair(i,aNormal));

            int NbNodes, buf[10];
            TColStd_Array1OfInteger nodeIDs(*buf,1,10);
            aMeshDS->GetNodesByElement(globalElementID,nodeIDs,NbNodes);
            for(int j=1; j<=NbNodes; j++) surroundingNodes.insert(nodeIDs(j));

            //! ---------------------------------------------------------------------------
            //! angle between the current element normal and the best visibility direction
            //! ---------------------------------------------------------------------------
            const QList<double> &nodeNormal = aMeshDS->myNodeNormals.value(globalNodeID);
            //cout<<"____node ID: "<<globalNodeID<<" ("<<nodeNormal[0]<<", "<<nodeNormal[1]<<", "<<nodeNormal[2]<<")____"<<endl;
            double dot = nx*nodeNormal[0]+ny*nodeNormal[1]+nz*nodeNormal[2];
            dot /= sqrt(nx*nx+ny*ny+nz*nz)*sqrt(pow(nodeNormal[0],2)+pow(nodeNormal[1],2)+pow(nodeNormal[2],2));
            if(dot<-1) dot = -1;
            if(dot>1) dot = 1;
            double bnk = std::acos(dot);
            if(bnk>bnk_max) bnk_max = bnk;
        }

        //! -----------------------------------------------------
        //! compute the visibility angle and insert into the map
        //! -----------------------------------------------------
        double betaVisibility = PI/2 - bnk_max;
        betaVisibilityField.insert(std::make_pair(globalNodeID,betaVisibility));

        //cout<<"____betaVisibility (rad,deg) = "<<betaVisibility<<"\t"<<betaVisibility*180/PI<<"____"<<endl;

        std::set<int>::iterator it_ = surroundingNodes.find(globalNodeID);
        surroundingNodes.erase(it_);

        /*
        //! -------------------------------------------------------
        //! given an advancing node P compute the average distance
        //! between the points and its surrounding nodes
        //! -------------------------------------------------------
        double averageDistance = 0;
        MeshVS_EntityType eType;
        int NbNodes;
        double bufd[3];
        TColStd_Array1OfReal coords(*bufd,1,3);
        aMeshDS->GetGeom(globalNodeID,false,coords,NbNodes,eType);
        double xP = coords(1); double yP = coords(2); double zP = coords(3);
        for(std::set<int>::iterator it = surroundingNodes.begin(); it!=surroundingNodes.end(); it++)
        {
            int globalNodeID_surr = *it;
            aMeshDS->GetGeom(globalNodeID_surr,false,coords,NbNodes,eType);
            averageDistance += sqrt(pow(xP-coords(1),2)+pow(yP-coords(2),2)+pow(zP-coords(3),2));
        }
        averageDistance /= surroundingNodes.size();
        */

        //! --------------------------
        //! generate pairs of normals
        //! --------------------------
        std::vector<std::pair<int,int>> vecPairs;
        size_t NbNormals = indexedMapOfNormals.size();

        for(size_t i=0; i<NbNormals-1; i++)
            for(size_t j=i+1; j<NbNormals; j++)
                vecPairs.push_back(std::make_pair(i,j));

        double sumBeta = 0.0;
        double betaAve = 0.0;
        int NbAdjacentPairs = 0;
        for(std::vector<std::pair<int,int>>::iterator it = vecPairs.begin(); it!=vecPairs.end(); it++)
        {
            int fn = it->first;
            int sn = it->second;

            //! ------------------------------
            //! check triangle pair adjacency
            //! ------------------------------
            int buf[3], NbNodes;
            TColStd_Array1OfInteger nodeIDs1(*buf,1,3);
            aMeshDS->GetNodesByElement(indexedMapOfElements.at(fn),nodeIDs1,NbNodes);
            std::vector<int> aFaceElement1;
            for(int i=1; i<=NbNodes; i++) aFaceElement1.push_back(nodeIDs1(i));
            TColStd_Array1OfInteger nodeIDs2(*buf,1,3);
            aMeshDS->GetNodesByElement(indexedMapOfElements.at(sn),nodeIDs2,NbNodes);
            std::vector<int> aFaceElement2;
            for(int i=1; i<=NbNodes; i++) aFaceElement2.push_back(nodeIDs2(i));

            bool shareAnEdge = areAdjacent(aFaceElement1,aFaceElement2);
            if(shareAnEdge==false) continue;

            NbAdjacentPairs++;
            const std::vector<double> &n1 = indexedMapOfNormals.at(fn);
            const std::vector<double> &n2 = indexedMapOfNormals.at(sn);
            //if(sqrt(n1[0]*n1[0]+n1[1]*n1[1]+n1[2]*n1[2])==0) exit(11);   // should never occur
            //if(sqrt(n2[0]*n2[0]+n2[1]*n2[1]+n2[2]*n2[2])==0) exit(12);   // should never occur

            //! ---------------------------------------------------------
            //! get ABCD and define as reference direction for the angle
            //! calculation the common edge (C-P)
            //! ---------------------------------------------------------
            mesh::meshPoint A,B,C,P;
            int element1 = indexedMapOfElements.at(fn);
            int element2 = indexedMapOfElements.at(sn);

            bool swappedAB = this->getLocalFanNodes(aMeshDS,globalNodeID,element1,element2,A,B,C,P);
            double l_CP = sqrt(pow(P.x-C.x,2)+pow(P.y-C.y,2)+pow(P.z-C.z,2));
            std::vector<double> refDir { -(C.x-P.x)/l_CP, -(C.y-P.y)/l_CP, -(C.z-P.z)/l_CP };
            double angle;
            if(swappedAB == false) angle = angleBetweenVector(n1,n2,refDir);
            else angle = angleBetweenVector(n2,n1,refDir);
            sumBeta += angle;
        }
        betaAve = sumBeta/NbAdjacentPairs;          // average angle

        //! ------------------------------------------------
        //! record the manifold characteristic into the map
        //! ------------------------------------------------
        betaAverageField.insert(std::make_pair(globalNodeID,PI-betaAve));
    }
}

//! ------------------------------
//! function: computeShrinkFactor
//! details:
//! ------------------------------
void prismaticLayer::computeShrinkFactor(const occHandle(Ng_MeshVS_DataSourceFace) &aMeshDS,
                                         QMap<int,double> &shrinkFactors)
{
    //FILE *fp = fopen("D:/shrink.txt","w");
    //fprintf(fp,"#NodeID\tx\t\ty\t\tz\t\tbetaAve\t\tbetaVisibility\tshrink\n");

    const double PI = 3.1415926534;
    const double eps = 0.01745329;           // 1 degree

    for(TColStd_MapIteratorOfPackedMapOfInteger it(aMeshDS->GetAllNodes()); it.More(); it.Next())
    {
        int globalNodeID = it.Key();
        double betaAve = betaAverageField.at(globalNodeID);
        double betaVisibility = betaVisibilityField.at(globalNodeID);
        double shrink = 0;
        if(betaAve<PI-eps) shrink = fabs(cos(betaVisibility));      // expansion in concave points
        if(betaAve>PI+eps) shrink =  -fabs(cos(betaVisibility));    // retraction in not concave points
        shrinkFactors.insert(globalNodeID,shrink);

        double P[3];
        this->getPointCoordinates(aMeshDS,globalNodeID,P);
        //fprintf(fp,"%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",globalNodeID,P[0],P[1],P[2],betaAve*180/PI,betaVisibility*180/PI,shrink);
    }
    //fclose(fp);
}

//! ---------------------------------
//! function: buildPrismaticElements
//! details:
//! ---------------------------------
bool prismaticLayer::buildPrismaticElements(const QList<occHandle(Ng_MeshVS_DataSourceFace)> &theInflatedMeshes,
                                            occHandle(Ng_MeshVS_DataSource3D) &prismaticMeshDS3D)
{
    cout<<"prismaticLayer::buildPrismaticElements()->____function called____"<<endl;

    if(theInflatedMeshes.length()==0) exit(44444);

    int NbTet, NbPyram, NbPrism, NbNullElement;
    NbTet = NbPyram = NbPrism = NbNullElement = 0;

    std::vector<meshElementByCoords> listOfMeshElements3D;

    //! -----------------------------------
    //! forward: from external to internal
    //! -----------------------------------
    for(int i=0; i<=theInflatedMeshes.length()-2; i++)
    {
        cout<<"@-----------------------------------@"<<endl;
        cout<<"@  working on mesh pair ("<<i<<","<<i+1<<")"<<endl;
        cout<<"@-----------------------------------@"<<endl;

        const occHandle(Ng_MeshVS_DataSourceFace) &baseMesh = theInflatedMeshes.at(i);
        TColStd_PackedMapOfInteger eMapBase = baseMesh->GetAllElements();

        const occHandle(Ng_MeshVS_DataSourceFace) &topMesh = theInflatedMeshes.at(i+1);
        TColStd_PackedMapOfInteger eMapTop = topMesh->GetAllElements();

        TColStd_MapIteratorOfPackedMapOfInteger eItBase, eItTop;

        for(eItBase.Initialize(eMapBase), eItTop.Initialize(eMapTop);
            eItBase.More(), eItTop.More();
            eItBase.Next(), eItTop.Next())
        {
            //! -------------
            //! base element
            //! -------------
            int globalElementID_base = eItBase.Key();
            int NbNodes;
            MeshVS_EntityType type;
            double buf[24], buf1[24];
            TColStd_Array1OfReal coords(*buf,1,24);
            TColStd_Array1OfReal coords1(*buf1,1,24);

            baseMesh->GetGeom(globalElementID_base,true,coords,NbNodes,type);

            //!cout<<"____working on element nr. "<<globalElementID_base<<"____"<<endl;

            QList<mesh::meshPoint> bottomFace;
            for(int n=0; n<NbNodes; n++)
            {
                int s = 3*n;
                mesh::meshPoint aMeshPoint(coords(s+1),coords(s+2),coords(s+3));
                bottomFace<<aMeshPoint;
            }

            //! ------------
            //! top element
            //! ------------
            int globalElementID_top = eItTop.Key();
            topMesh->GetGeom(globalElementID_top,true,coords1,NbNodes,type);

            QList<mesh::meshPoint> topFace;
            for(int n=0; n<NbNodes; n++)
            {
                int s = 3*n;
                mesh::meshPoint aMeshPoint1(coords1(s+1),coords1(s+2),coords1(s+3));
                topFace<<aMeshPoint1;
            }

            //! -----------------------------
            //! create the prismatic element
            //! -----------------------------
            meshElementByCoords ameshElementByCoords(bottomFace,topFace);
            if(ameshElementByCoords.type!=NULL_ELEMENT) listOfMeshElements3D.push_back(ameshElementByCoords);

            //! ----------------------------
            //! diagnostic - elements count
            //! ----------------------------
            switch(ameshElementByCoords.type)
            {
            case TET: NbTet++; /*cout<<"_TET_"*/; break;
            case PYRAM: NbPyram++; /*cout<<"_PYRAM_"*/; break;
            case PRISM: NbPrism++; /*cout<<"_PRISM_"*/; break;
            case NULL_ELEMENT: NbNullElement++; /*cout<<"_NULL ELEMENT_"*/; break;
            }
        }

        cout<<endl;
        cout<<"@____prismatic mesh summary____@"<<endl;
        cout<<"@____Nb TET: "<<NbTet<<endl;
        cout<<"@____Nb PYRAM: "<<NbPyram<<endl;
        cout<<"@____Nb PRISM: "<<NbPrism<<endl;
        cout<<"@____Nb NULL_: "<<NbNullElement<<endl;
        cout<<"@____end of summary____@"<<endl;
    }

    prismaticMeshDS3D = new Ng_MeshVS_DataSource3D(listOfMeshElements3D);
    if(prismaticMeshDS3D.IsNull()) return false;
    return true;
}

//! --------------------------------------------------------------
//! function: isSelfIntersecting
//! details:  check self intersections (code defined in iglTools)
//! --------------------------------------------------------------
bool prismaticLayer::isSelfIntersecting(const occHandle(Ng_MeshVS_DataSourceFace) &aMesh, std::vector<int> &badTriangles)
{
    bool foundSelfIntersection = iglTools::iglCheckSelfIntersection(aMesh, badTriangles);
    return foundSelfIntersection;
}

//! ---------------------------------------------------------
//! function: areIntersecting
//! details:  check intersection of a mesh with another mesh
//!           (code defined in iglTools)
//! ---------------------------------------------------------
bool prismaticLayer::areIntersecting(const occHandle(Ng_MeshVS_DataSourceFace) &aMesh1,
                                     const occHandle(Ng_MeshVS_DataSourceFace) &aMesh2,
                                     std::vector<int> &badNodes1,
                                     std::vector<int> &badNodes2)
{
    bool foundIntersection = iglTools::iglIntersectMeshes(aMesh1,aMesh2,badNodes1,badNodes2);
    return foundIntersection;
}

//! ---------------------------------------------------------
//! function: postBLBuilder
//! details:  build the displaced mesh list. Here the volume
//!           mesh must have been generated in advance
//! ---------------------------------------------------------
bool prismaticLayer::inflateMeshAndCompress(QList<occHandle(Ng_MeshVS_DataSourceFace)> &theInflatedMeshes,
                                            occHandle(Ng_MeshVS_DataSource3D) &preInflationVolumeMeshDS)
{
    cout<<"****************************************"<<endl;
    cout<<" prismaticLayer::inflateMeshAndCompress "<<endl;
    cout<<"****************************************"<<endl;
    this->mergeFaceMeshes();

    occHandle(Ng_MeshVS_DataSourceFace) prismaticMeshDS = new Ng_MeshVS_DataSourceFace(myPrismaticFacesSumMeshDS);
    occHandle(Ng_MeshVS_DataSourceFace) surfaceMeshDS = new Ng_MeshVS_DataSourceFace(myOverallSumMeshDS);
    occHandle(Ng_MeshVS_DataSourceFace) nonPrismaticMeshDS = new Ng_MeshVS_DataSourceFace(myNonPrismaticFacesSumMeshDS);

    surfaceMeshDS->computeNormalAtNodes();
    prismaticMeshDS->computeNormalAtNodes();
    prismaticMeshDS->computeFreeMeshSegments();

    //! -------------------------------------
    //! make a copy the initial (outer) mesh
    //! -------------------------------------
    occHandle(Ng_MeshVS_DataSourceFace) theMeshToInflate_old = new Ng_MeshVS_DataSourceFace(prismaticMeshDS);
    theMeshToInflate_old->computeNormalAtNodes();

    //! ------------------------------------
    //! store the first mesh (the exterior)
    //! ------------------------------------
    theInflatedMeshes<<theMeshToInflate_old;

    //! -------------------------------------------
    //! reference geometry: use the abstract class
    //! and create it into the heap
    //! -------------------------------------------
    shared_ptr<geometryFace> sharedGeo;

    bool useGeometryForNodeAdjustment = false;
    if(useGeometryForNodeAdjustment)
    {
        //! ----------------------
        //! what is not prismatic
        //! ----------------------
        std::vector<TopoDS_Face> vecFaces;
        for(size_t i=0; i<myPrismaticFaces.size(); i++)
        {
            int faceNr = myPrismaticFaces.at(i);
            TopoDS_Face aFace = TopoDS::Face(myMeshDB->MapOfBodyTopologyMap.value(myBodyIndex).faceMap.FindKey(faceNr));
            if(aFace.IsNull()) continue;
            vecFaces.push_back(aFace);
        }
        sharedGeo = make_shared<OCCFace>();
        static_cast<OCCFace*>(sharedGeo.get())->setGeometry(vecFaces);
    }
    else
    {
        //! --------------------------------------------
        //! what is not defined as prismatic - discrete
        //! this is used as geometry reference
        //! --------------------------------------------
        sharedGeo = make_shared<meshFace>();
        static_cast<meshFace*>(sharedGeo.get())->setGeometry(nonPrismaticMeshDS);
    }

    //! ---------------------------------------------------
    //! the work is done on the prismatic mesh data source
    //! ---------------------------------------------------
    for(int n=1; n<=myNbLayers; n++)
    {
        double displacement = myLayerThickness.at(n-1);
        cout<<"prismaticLayer::inflateMeshAndCompress->____generating layer: "<<n<<" with current thickness: "<<displacement<<"____@"<<endl;

        occHandle(Ng_MeshVS_DataSourceFace) theMeshToInflate_new = new Ng_MeshVS_DataSourceFace(theMeshToInflate_old);
        theMeshToInflate_new->computeNormalAtNodes();
        this->computeBeta(theMeshToInflate_new);

        //! ---------------
        //! shrink factors
        //! ---------------
        QMap<int,double> shrinkFactors;
        this->computeShrinkFactor(theMeshToInflate_new,shrinkFactors);

        //! ---------------------------------
        //! guiding vectors - "best normals"
        //! ---------------------------------
        QMap<int,QList<double>> normals = theMeshToInflate_new->myNodeNormals;

        //! ---------------
        //! classify nodes
        //! ---------------
        std::map<int,int> mapNodeTypes;
        this->classifyNodes(theMeshToInflate_new,mapNodeTypes);

        //! -----------------------------------------------------------------------------------------------
        //! smooth the guiding vectors directions - use 5/10 smoothing steps
        //! flag normalize == true since the "rotation" of the vector is of interest (change in direction)
        //! -----------------------------------------------------------------------------------------------
        smoothingTools::fieldSmoother(normals,theMeshToInflate_new,betaAverageField,betaVisibilityField,mapNodeTypes,10,true);

        //! ------------------------------------
        //! map of the local marching distances
        //! ------------------------------------
        QMap<int,double> marchingDistanceMap;
        for(QMap<int,QList<double>>::const_iterator it = normals.cbegin(); it!= normals.cend(); ++it)
        {
            int globalNodeID = it.key();
            double shrinkFactor = shrinkFactors.value(globalNodeID);
            //double cutOff = myLayerHCutOff.value(globalNodeID);

            //! -----------------------------------------------------
            //! cutoff is 1.0 since all the points will be displaced
            //! component of the displacement field
            //! -----------------------------------------------------
            double cutoff = 1.0;
            double marchingDistance = displacement*(1+shrinkFactor)*cutoff;
            marchingDistanceMap.insert(globalNodeID,marchingDistance);
        }

        //! ------------------------------
        //! smooth the marching distances
        //! ------------------------------
        smoothingTools::scalarFieldSmoother(marchingDistanceMap,theMeshToInflate_new,betaAverageField,mapNodeTypes,10);

        //! -------------------------
        //! check marching distances
        //! 1. Its marching distance is less than a fraction, say 0.60, of the average length of edges sharing p in its manifold
        //! 2.
        //! -------------------------
        //! ...


        //! -------------------
        //! displacement field
        //! -------------------
        QMap<int,QList<double>> displacementsField;
        for(QMap<int,QList<double>>::const_iterator it = normals.cbegin(); it!= normals.cend(); ++it)
        {
            int globalNodeID = it.key();
            const QList<double> &curNormal = it.value();

            //! --------------------
            //! nodal displacements
            //! --------------------
            double marchingDistance = marchingDistanceMap.value(globalNodeID);
            double vx = -curNormal[0]*marchingDistance;
            double vy = -curNormal[1]*marchingDistance;
            double vz = -curNormal[2]*marchingDistance;

            bool correctBoundary = false;
            if(correctBoundary)
            {
                //! -----------------------------------------------------------------
                //! correction of the displacement along the prismatic mesh boundary
                //! -----------------------------------------------------------------
                if(prismaticMeshDS->myBoundaryPoints.contains(globalNodeID))
                {
                    //! -----------------------------
                    //! the current node coordinates
                    //! -----------------------------
                    int localNodeID = prismaticMeshDS->myNodesMap.FindIndex(globalNodeID);
                    std::vector<double> aPoint = prismaticMeshDS->getNodeCoordinates(localNodeID);

                    //! ------------------------------------------------------------------
                    //! the position of the point after translation along the node normal
                    //! ------------------------------------------------------------------
                    double x_displacedAlongNodeNormal = aPoint[0] + vx;
                    double y_displacedAlongNodeNormal = aPoint[1] + vy;
                    double z_displacedAlongNodeNormal = aPoint[2] + vz;

                    double theDisplacedPoint[3];
                    theDisplacedPoint[0] = x_displacedAlongNodeNormal;
                    theDisplacedPoint[1] = y_displacedAlongNodeNormal;
                    theDisplacedPoint[2] = z_displacedAlongNodeNormal;

                    //! -------------------------------------------------------
                    //! the projection of the point onto what is not prismatic
                    //! -------------------------------------------------------
                    double theProjectedPoint[3];
                    bool isDone = sharedGeo->pointProjection(theDisplacedPoint,theProjectedPoint);
                    if(isDone)
                    {
                        //! ---------------------------
                        //! the corrected displacement
                        //! ---------------------------
                        vx = theProjectedPoint[0] - aPoint[0];
                        vy = theProjectedPoint[1] - aPoint[1];
                        vz = theProjectedPoint[2] - aPoint[2];
                    }
                    else
                    {
                        //! --------------------------------------------------------
                        //! here, in case the projection onto the original geometry
                        //! is not done, we chose to apply the displacement. An
                        //! choice could be locking the displacement
                        //! --------------------------------------------------------
                        cerr<<"____projection not done____"<<endl;
                        vx = vy = vz = 0;
                    }
                }
            }
            QList<double> localFieldValue;
            localFieldValue<<vx<<vy<<vz;
            displacementsField.insert(globalNodeID,localFieldValue);
        }

        //! -------------------------------------------
        //! displace the prismatic and the volume mesh
        //! -------------------------------------------
        theMeshToInflate_new->displaceMySelf(displacementsField);
        theInflatedMeshes<<theMeshToInflate_new;

        //! --------------------------------------------------------
        //! translation of nodes without mesh compression: this was
        //! the initial attempt to compress, obviously leading to
        //! mesh elementd folding for large displacements or small
        //! element size - left here for documentation
        //! --------------------------------------------------------
        //preInflationVolumeMeshDS->displaceMySelf(displacementsField);

        //! -------------------------------------------------
        //! "as rigid as possible" volume mesh deformation
        //! this substitutes for the "displaceMySelf" method
        //! -------------------------------------------------
        shared_ptr<meshSelectNLayers> meshSelector(new meshSelectNLayers);
        meshSelector->setVolumeMesh(preInflationVolumeMeshDS);

        occHandle(Ng_MeshVS_DataSource3D) selectedVolumeMeshDS;
        meshSelector->setStartingSurfaceElements(myPrismaticFacesSumMeshDS);

        //! ------------------------------------------------
        //! nodes not affected by the deformation algorithm
        //! nodes affected by the deformation algorithm
        //! ------------------------------------------------
        std::vector<int> constrainedGlobalNodeIDs;
        std::vector<int> groupOfMovingNodes;

        constrainedGlobalNodeIDs.clear();
        groupOfMovingNodes.clear();

        //! -------------------------
        //! 3 steps seems reasonable
        //! -------------------------
        meshSelector->setSteps(3);
        bool isDone = meshSelector->getElements(selectedVolumeMeshDS);
        if(!isDone)
        {
            cout<<"prismaticLayer::inflateMeshAndCompress->____cannot compress mesh____"<<endl;
            return false;
        }
        TColStd_PackedMapOfInteger mapOfGlobalNodeIDs = selectedVolumeMeshDS->GetAllNodes();

        //! --------------------------------------------
        //! nodes affected by the deformation algorithm
        //! --------------------------------------------
        for(TColStd_MapIteratorOfPackedMapOfInteger it(mapOfGlobalNodeIDs);it.More();it.Next())
        {
            int globalNodeID = it.Key();
            groupOfMovingNodes.push_back(globalNodeID);
        }
        //! ------------------------------------------------
        //! nodes not affected by the deformation algorithm
        //! ------------------------------------------------
        for(TColStd_MapIteratorOfPackedMapOfInteger it(preInflationVolumeMeshDS->GetAllNodes()); it.More(); it.Next())
        {
            int globalNodeID = it.Key();

            //! ------------------------------------------------------------------------------
            //! the current node is not contained on the vector of node that are free to move
            //! ------------------------------------------------------------------------------
            if(!mapOfGlobalNodeIDs.Contains(globalNodeID)) constrainedGlobalNodeIDs.push_back(globalNodeID);
        }
        //! ---------------------------------------------
        //! define a minimum number of constrained nodes
        //! ---------------------------------------------
        if(constrainedGlobalNodeIDs.size()<3)
        {
            cout<<"prismaticLayer::inflateMeshAndCompress->____cannot compress mesh____"<<endl;
            return false;
        }

        //! -------------------------------------------------------
        //! mode: "0" harmonic/biharmonic "1" as rigid as possible
        //! -------------------------------------------------------
        int mode = 1;
        preInflationVolumeMeshDS->displaceMySelf_asRigidAsPossible(displacementsField,constrainedGlobalNodeIDs,mode);
        //preInflationVolumeMeshDS->displaceMySelf_asRigidAsPossible(displacementsField,groupOfMovingNodes,mode);

        //! ------------------------------------
        //! update the "old" displacement field
        //! ------------------------------------
        theMeshToInflate_old = theMeshToInflate_new;
    }
    return true;
}

//! ---------------------------------------------------------------
//! function: setStopButtonEnabled
//! details:  enable/disable the progress indicator stop button
//!           here an event is posted to the application main loop
//! ---------------------------------------------------------------
void prismaticLayer::setStopButtonEnabled(bool isEnabled)
{
    QProgressEvent *pe;
    if(isEnabled) pe = new QProgressEvent(QProgressEvent_EnableStop,0,0,0,"",QProgressEvent_None,0,0,0,myTask);
    else pe = new QProgressEvent(QProgressEvent_DisableStop,0,0,0,"",QProgressEvent_None,0,0,0,myTask);
    QApplication::postEvent(myProgressIndicator,pe);
}

//! ------------------------------------------------------
//! function: generateTetLayers
//! details:  return the boundary mesh made of tetrahedra
//!           and the internal deformed surface mesh
//! ------------------------------------------------------
void prismaticLayer::generateTetLayers(occHandle(Ng_MeshVS_DataSource3D) &meshAtWalls,
                                       occHandle(Ng_MeshVS_DataSourceFace) &lastInflatedMesh)
{
    cout<<"prismaticLayers::generateTetLayers()->____function called____"<<endl;

    //! --------------------------------
    //! init the secondary progress bar
    //! --------------------------------
    QProgressEvent *progressEvent;
    if(myProgressIndicator!=Q_NULLPTR)
    {
        progressEvent = new QProgressEvent(QProgressEvent_None,0,9999,0,"Generating boundary mesh",
                                           QProgressEvent_Init,0,myNbLayers,0,"Building prismatic mesh");
        QApplication::postEvent(myProgressIndicator,progressEvent);
        QApplication::processEvents();
    }

    //! ----------------------------
    //! clear reduction factors map
    //! ----------------------------
    mapOfReductionFactor.clear();

    //! -----------------------------
    //! the volume elements at walls
    //! -----------------------------
    std::vector<meshElementByCoords> volumeElementsAtWalls;

    //! ---------------------------------------
    //! decide which nodes should be displaced
    //! ---------------------------------------
    this->computeVecFieldCutOff(myLockBoundary);

    //! ---------------------------
    //! make a copy the outer mesh
    //! ---------------------------
    occHandle(Ng_MeshVS_DataSourceFace) theMeshToInflate = new Ng_MeshVS_DataSourceFace(myOverallSumMeshDS);
    if(myPrismaticFacesSumMeshDS->myBoundarySegments.isEmpty()) myPrismaticFacesSumMeshDS->computeFreeMeshSegments();

    //! ------------------------------------------------------------------
    //! analyze gaps at the beginning: scan all the mesh nodes.
    //! compute the distance from closest point on the mesh => "distance"
    //! apply reduction of first layer thickness if applicable
    //! ------------------------------------------------------------------
    double Sigma = 0;                                           // a parameter for reduction factor calculation
    for(int n=0; n<myNbLayers; n++) Sigma += pow(myExpRatio,n);

    pointToMeshDistance aDistanceMeter;                         // distance meter tool
    aDistanceMeter.init(theMeshToInflate);                      // init with mesh
    double firstLayerThickness = myLayerThickness.at(0);        // first layer thickness
    std::map<int,double> mapOfFirstLayerReductionFactor;        // fill a map (nodeID, reduction factor)

    if(theMeshToInflate->myNodeNormals.isEmpty()) theMeshToInflate->computeNormalAtNodes();
    const QMap<int,QList<double>> &normals = theMeshToInflate->myNodeNormals;

    //! ---------------------------------------------
    //! start filling the map of compression factors
    //! ---------------------------------------------
    for(QMap<int,QList<double>>::const_iterator itn = normals.cbegin(); itn!=normals.cend(); itn++)
    {
        int globalNodeID = itn.key();
        const QList<double> &normal = itn.value();

        double P[3];
        float distance;
        this->getPointCoordinates(theMeshToInflate,globalNodeID,P);

        //! --------------------------------------------------
        //! invert the normal: it points towards the material
        //! for gap calculation
        //! --------------------------------------------------
        double nx = -normal[0];
        double ny = -normal[1];
        double nz = -normal[2];
        double dir[3] {nx,ny,nz};
        aDistanceMeter.distance(P,dir,&distance);   // compute the distance

        bool compress = true;
        if(compress==true)
        {
            //! ---------------------------------------------
            //! compress layers in order to avoid collisions
            //! ---------------------------------------------
            const double compressionFactor = 0.25;
            if(myTotalThickness<distance/3.0)
            {
                //! the node can move with the nominal marching distance
                mapOfFirstLayerReductionFactor.insert(std::make_pair(globalNodeID,1.0));
            }
            else
            {
                //! the node can move with smaller marching distance
                double firstLayerThicknessNew = compressionFactor*distance/Sigma;
                double reductionFactor = firstLayerThicknessNew/firstLayerThickness;
                mapOfFirstLayerReductionFactor.insert(std::make_pair(globalNodeID,reductionFactor));
                mapOfReductionFactor.insert(std::make_pair(globalNodeID,reductionFactor));
            }
        }
        else
        {
            if(myTotalThickness<distance/3.0)
            {
                //! the node can move without compression
                mapOfFirstLayerReductionFactor.insert(std::make_pair(globalNodeID,1.0));
            }
            else
            {
                //! lock point - the node cannot move
                myLayerHCutOff.insert(globalNodeID,0);
                mapOfFirstLayerReductionFactor.insert(std::make_pair(globalNodeID,0.0));
                mapOfReductionFactor.insert(std::make_pair(globalNodeID,0.0));
            }
        }
    }

    //! --------------------------------
    //! build the surrounding nodes map
    //! --------------------------------
    mySurroudingNodesMap.clear();
    this->buildSurroundingNodesMap(theMeshToInflate,mySurroudingNodesMap);

    //! --------------------
    //! generate the layers
    //! --------------------
    for(int n=1; n<=myNbLayers; n++)
    {
        //! --------------
        //! send progress
        //! --------------
        if(myProgressIndicator!=Q_NULLPTR)
        {
            progressEvent = new QProgressEvent(QProgressEvent_None,0,9999,0,QString("Generating boundary mesh layer: %1").arg(n),
                                               QProgressEvent_Update,-1,-1,n,"Building prismatic mesh");
            QApplication::postEvent(myProgressIndicator,progressEvent);
            QApplication::processEvents();
        }

        //! ----------------------------
        //! the current layer thickness
        //! ----------------------------
        double displacement = myLayerThickness.at(n-1);

        //! -----------------------------------------------------------
        //! manifold characteristic for each point of the current mesh
        //! visibility angle for each point of the current mesh
        //! -----------------------------------------------------------
        this->computeBeta(theMeshToInflate);

        //! ------------------------------------------------------------
        //! collapsed into a method the generation of one layer of tets
        //! ------------------------------------------------------------
        this->generateOneTetLayer(theMeshToInflate,displacement,volumeElementsAtWalls,mapOfFirstLayerReductionFactor);
    }

    //! -------------------------------
    //! the last modified surface mesh
    //! -------------------------------
    lastInflatedMesh = theMeshToInflate;
    /*
    //lastInflatedMesh = theMeshToInflate;
    lastInflatedMesh = myOverallSumMeshDS;
    */

    //! ----------------------------------------------
    //! construction of the tetrahedral boundary mesh
    //! ----------------------------------------------
    meshAtWalls = new Ng_MeshVS_DataSource3D(volumeElementsAtWalls,true,true);
    cout<<"____number of generated elements: "<<meshAtWalls->GetAllElements().Extent()<<"____"<<endl;
}

//! ------------------------------
//! function: generateMeshAtWalls
//! details:
//! ------------------------------
bool prismaticLayer::generateMeshAtWalls(occHandle(Ng_MeshVS_DataSource3D) &meshAtWalls,
                                         occHandle(Ng_MeshVS_DataSourceFace) &lastInflatedMesh,
                                         QList<occHandle(Ng_MeshVS_DataSourceFace)> &listOfInflatedMeshes)
{
    //! ---------------------------------------------
    //! merge the prismatic faces into a unique mesh
    //! ---------------------------------------------
    bool isDone = this->mergeFaceMeshes();
    if(!isDone) return false;

    switch(myBoundaryMeshType)
    {
    //! ---------------------
    //! Hybrid mesh at walls
    //! ---------------------
    case 0:
    {
        //! --------------------------
        //! displace the surface mesh
        //! --------------------------
        QList<occHandle(Ng_MeshVS_DataSourceFace)> theInflatedMeshes;
        isDone = this->inflateMesh(theInflatedMeshes);
        if(!isDone) return false;

        //! ---------------------------------------------------------------------
        //! generate the prismatic 3D mesh from the list of the displaced meshes
        //! ---------------------------------------------------------------------
        //cout<<"prismaticLayer::generateMeshAtWalls()->____start generating prismatic 3D elements____"<<endl;
        //cout<<"prismaticLayer::generateMeshAtWalls()->____number of inflated meshes: "<<theInflatedMeshes.length()<<"____"<<endl;
        isDone = this->buildPrismaticElements(theInflatedMeshes,meshAtWalls);

        //cout<<"prismaticLayer::generateMeshAtWalls()->____prismatic mesh: number of elements generated: "<<meshAtWalls->GetAllElements().Extent()<<"____"<<endl;

        if(!isDone) return false;

        //! -----------------------
        //! the last inflated mesh
        //! -----------------------
        lastInflatedMesh = theInflatedMeshes.last();
    }
        break;

    //! -----------------------------------------------------------------------
    //! Tetrahedral mesh at walls. Important note: when using generateTetLayer
    //! the "lastInflatedMesh" is initial, outer surface mesh with advancing
    //! node properly displaced. It can be used for generating a volume mesh
    //! -----------------------------------------------------------------------
    case 1:
    {
        this->generateTetLayers(meshAtWalls,lastInflatedMesh);
    }
        break;
    }
    return true;
}

//! --------------------
//! function: invertTet
//! details:
//! --------------------
void prismaticLayer::invertTet(meshElementByCoords &aTet)
{
    int mask[4] {0,2,1,3};
    QList<mesh::meshPoint> pointListNew;
    for(int i=0; i<4; i++) pointListNew<<aTet.pointList[mask[i]];
    aTet.pointList.clear();
    aTet.pointList<<pointListNew;
}

//! --------------------------------
//! function: checkSelfIntersection
//! details:
//! --------------------------------
bool prismaticLayer::checkSelfIntersection(const occHandle(Ng_MeshVS_DataSourceFace) &aMeshDS,
                                           const occHandle(Ng_MeshVS_DataSourceFace) &aMeshDS_old,
                                           bool correct)
{
    cout<<"prismaticLayer::checkSelfIntersection()->____start checking self intersections____"<<endl;
    std::vector<int> badNodes;
    bool selfIntersection = this->isSelfIntersecting(aMeshDS, badNodes);
    if(correct==false) return selfIntersection;

    if(selfIntersection)
    {
        cout<<"@--------------------------------------"<<endl;
        cout<<"@   a self intersection has been found "<<endl;
        for(int k=0; k<badNodes.size(); k++)
        {
            int localNodeID = badNodes[k]+1;
            std::vector<double> old_coords = aMeshDS_old->getNodeCoordinates(localNodeID);
            aMeshDS->changeNodeCoords(localNodeID,old_coords);
            int globalNodeID = aMeshDS->myNodesMap.FindKey(localNodeID);
            myLayerHCutOff.insert(globalNodeID,0);
        }
        cout<<"@   number of corrected nodes: "<<badNodes.size()<<endl;
        cout<<"@--------------------------------------"<<endl;
    }
    return selfIntersection;
}

//! ----------------------------------
//! function: checkMutualIntersection
//! details:
//! ----------------------------------
bool prismaticLayer::checkMutualIntersection(const occHandle(Ng_MeshVS_DataSourceFace) &aMeshDS,
                                             const occHandle(Ng_MeshVS_DataSourceFace) &aMeshDS_old,
                                             bool correct)
{
    cout<<"prismaticLayer::checkMutualIntersection()->____start checking mutual intersections____"<<endl;
    std::vector<int> badNodes_inflatedMesh, badNodes_sourceMesh;
    bool intersection = this->areIntersecting(aMeshDS,aMeshDS_old,badNodes_inflatedMesh,badNodes_sourceMesh);
    if(correct==false) return intersection;

    if(intersection)
    {
        cout<<"@--------------------------------------------"<<endl;
        cout<<"@  a mutual mesh intersection has been found "<<endl;
        for(int k=0; k<badNodes_inflatedMesh.size(); k++)
        {
            int localNodeID = badNodes_inflatedMesh[k]+1;
            std::vector<double> old_coords = aMeshDS_old->getNodeCoordinates(localNodeID);
            aMeshDS->changeNodeCoords(localNodeID,old_coords);
            int globalNodeID = aMeshDS->myNodesMap.FindKey(localNodeID);
            myLayerHCutOff.insert(globalNodeID,0);
        }
        cout<<"@  number of corrected nodes: "<<badNodes_inflatedMesh.size()<<endl;
        cout<<"@--------------------------------------------"<<endl;
    }
    return intersection;
}

//! ------------------------------
//! function: generateOneTetLayer
//! details:
//! ------------------------------
void prismaticLayer::generateOneTetLayer(occHandle(Ng_MeshVS_DataSourceFace) &theMeshToInflate,
                                         double displacement,
                                         std::vector<meshElementByCoords> &volumeElementsAtWalls,
                                         const std::map<int,double> &mapOfFirstLayerReductionFactor)
{
    cout<<"prismaticLayer::generateOneTetLayer()->____function called____"<<endl;
    const double PI = 3.1415926534;

    //! ----------------
    //! diagnostic info
    //! ----------------
    int NbInvalidTets = 0;
    int blockedByLenght = 0;
    int blockedByAngle = 0;

    //! -----------------------
    //! compute shrink factors
    //! -----------------------
    QMap<int,double> shrinkFactors;
    this->computeShrinkFactor(theMeshToInflate,shrinkFactors);

    //! ---------------------------------------
    //! calculate the nodal displacement field
    //! ---------------------------------------
    theMeshToInflate->computeNormalAtNodes();
    theMeshToInflate->computeFreeMeshSegments();

    //! ---------------------------------
    //! guiding vectors - "best normals"
    //! ---------------------------------
    QMap<int,QList<double>> &normals = theMeshToInflate->myNodeNormals;

    //! --------------------------------------------
    //! check gaps along slightly different normals
    //! --------------------------------------------
    pointToMeshDistance aDistanceMeter;
    aDistanceMeter.init(theMeshToInflate);
    for(TColStd_MapIteratorOfPackedMapOfInteger it(theMeshToInflate->GetAllNodes()); it.More(); it.Next())
    {
        int globalNodeID = it.Key();
        double P[3];
        std::vector<float> gaps;
        float gap;
        QList<double> normal = theMeshToInflate->myNodeNormals.value(globalNodeID);
        double dir[3] { -normal[0], -normal[1], -normal[2] };
        this->getPointCoordinates(theMeshToInflate,globalNodeID,P);

        //! --------------------------
        //! gap along the node normal
        //! --------------------------
        aDistanceMeter.distance(P,dir,&gap);
        gaps.push_back(gap);

        //! ----------------------------------------------------------------------
        //! create a start direction which will be rotated along the node normal:
        //! for this initial deviation the visibility angle (x2) is used
        //! This vector will be rotated by steps around the node normal, creating
        //! a cone of rays
        //! ----------------------------------------------------------------------
        double angleDeviation = 2*betaVisibilityField.at(globalNodeID);
        double v_0[3];
        if(dir[2]!=0)
        {
            v_0[0] = dir[0];
            v_0[1] = dir[1];
            v_0[2] = (cos(angleDeviation)-dir[0]*dir[0]-dir[1]*dir[1])/dir[2];
        }
        else if(dir[1]!=0)
        {
            v_0[0] = (cos(angleDeviation)-dir[0]*dir[0]-dir[2]*dir[2])/dir[1];
            v_0[1] = dir[1];
            v_0[2] = dir[2];
        }
        else if(dir[0]!=0)
        {
            v_0[0] = (cos(angleDeviation)-dir[1]*dir[1]-dir[2]*dir[2])/dir[0];
            v_0[1] = dir[1];
            v_0[2] = dir[2];
        }
        int N = 16;                     // cast 16 rays from P
        double dangle = 2*PI/N;
        for(int i = 0; i<N; i++)
        {
            double dir_rot[3];
            this->rotateVector(v_0,dir,dangle,dir_rot);
            aDistanceMeter.distance(P,dir_rot,&gap);
            gaps.push_back(gap);
            v_0[0] = dir_rot[0];
            v_0[1] = dir_rot[1];
            v_0[2] = dir_rot[2];
        }

        gap = *std::min(gaps.begin(),gaps.end());
        if(displacement>0.25*gap)
        {
            myLayerHCutOff.insert(globalNodeID,0);
            cout<<"prismaticLayer::generateOneTetLayer()->____blocked node ID: "<<globalNodeID<<"____displacement: "<<displacement<<"\t gap: "<<gap<<"____"<<endl;
        }
    }

    //! ---------------
    //! classify nodes
    //! ---------------
    std::map<int,int> mapNodeTypes;
    this->classifyNodes(theMeshToInflate,mapNodeTypes);

    //! -----------------------------------------------------------------------------------------------
    //! smooth the guiding vectors directions - use 5/10 smoothing steps
    //! flag normalize == true since the "rotation" of the vector is of interest (change in direction)
    //! -----------------------------------------------------------------------------------------------
    smoothingTools::fieldSmoother(normals,theMeshToInflate,betaAverageField,betaVisibilityField,mapNodeTypes,10,true);

    //! ------------------------------------
    //! map of the local marching distances
    //! ------------------------------------
    QMap<int,double> marchingDistanceMap;
    for(QMap<int,QList<double>>::const_iterator it = normals.cbegin(); it!= normals.cend(); ++it)
    {
        int globalNodeID = it.key();
        double shrinkFactor = shrinkFactors.value(globalNodeID);
        double cutOff = myLayerHCutOff.value(globalNodeID);
        double marchingDistance = displacement*(1+shrinkFactor)*cutOff*mapOfFirstLayerReductionFactor.at(globalNodeID);
        marchingDistanceMap.insert(globalNodeID,marchingDistance);
    }

    //! ------------------------------
    //! smooth the marching distances
    //! ------------------------------
    smoothingTools::scalarFieldSmoother(marchingDistanceMap,theMeshToInflate,betaAverageField,mapNodeTypes,10);

    //! ---------------------------------------------------------------------------------------------------------------------
    //! front analysis: check the marching distances
    //! 1. Its marching distance is less than a fraction, say 0.60, of the average length of edges sharing p in its manifold
    //! 2. The prismatic element created by marching p has an acceptable quality
    //! ---------------------------------------------------------------------------------------------------------------------
    for(TColStd_MapIteratorOfPackedMapOfInteger it(theMeshToInflate->GetAllNodes()); it.More(); it.Next())
    {
        int globalNodeID = it.Key();
        int localNodeID = theMeshToInflate->myNodesMap.FindIndex(globalNodeID);
        const QList<int> &surroundingElementsIDs_local = theMeshToInflate->myNodeToElements.value(localNodeID);
        int NbSurroundingElements = surroundingElementsIDs_local.length();
        std::set<int> setOfSurroundingNodesIDs_global;
        for(int n=0; n<NbSurroundingElements; n++)
        {
            int curLocalElementID = surroundingElementsIDs_local[n];
            int curGlobalElementID = theMeshToInflate->myElementsMap.FindKey(curLocalElementID);
            int NbNodes, buf[12];
            TColStd_Array1OfInteger nodeIDs(*buf,1,12);
            theMeshToInflate->GetNodesByElement(curGlobalElementID,nodeIDs,NbNodes);
            for(int k=1; k<=NbNodes; k++) setOfSurroundingNodesIDs_global.insert(nodeIDs(k));
        }
        std::set<int>::iterator it_ = setOfSurroundingNodesIDs_global.find(globalNodeID);
        setOfSurroundingNodesIDs_global.erase(it_);

        //! -----------------------------------------------
        //! perform first check - compute average distance
        //! -----------------------------------------------
        double P[3], averageDistance = 0;
        this->getPointCoordinates(theMeshToInflate,globalNodeID,P);
        for(std::set<int>::iterator it__ = setOfSurroundingNodesIDs_global.begin(); it__!=setOfSurroundingNodesIDs_global.end(); it__++)
        {
            int curGlobalNodeID = *it__;
            double Pcurr[3];
            this->getPointCoordinates(theMeshToInflate,curGlobalNodeID,Pcurr);
            double distance = sqrt(pow(P[0]-Pcurr[0],2)+pow(P[1]-Pcurr[1],2)+pow(P[2]-Pcurr[2],2));
            averageDistance += distance;
        }
        averageDistance /= NbSurroundingElements;

        //! -------------------------------------------------------------------
        //! perform first check - stop propagating node if the check is not ok
        //! use 0.5 - 0.6
        //! -------------------------------------------------------------------
        if(marchingDistanceMap.value(globalNodeID)>0.50 * averageDistance)
        {
            myLayerHCutOff.insert(globalNodeID,0);
            blockedByLenght++;
            continue;   // jump over the subsequent check, if any
        }

        //! ------------------------------------------------------------------------------------
        //! second check - compute angles between the guiding vector and the faces at the basis
        //! if the maximum among these angles is greater than 80° block the advancing node
        //! ------------------------------------------------------------------------------------
        double P_star[3];
        this->getPointCoordinates(theMeshToInflate,globalNodeID,P_star);
        const QList<double> &normal = normals.value(globalNodeID);
        double curMarchingDistance = marchingDistanceMap.value(globalNodeID);
        P_star[0] += -normal[0]*curMarchingDistance;
        P_star[1] += -normal[1]*curMarchingDistance;
        P_star[2] += -normal[2]*curMarchingDistance;

        double maxAngle = -1e10;
        double minAngle = 1e10;
        const std::vector<int> &vecSurroundingNodes = mySurroudingNodesMap.at(globalNodeID);
        for(std::vector<int>::const_iterator it = vecSurroundingNodes.cbegin(); it!=vecSurroundingNodes.cend(); it++)
        {
            int globalNodeID_surr = *it;
            double V[3];
            this->getPointCoordinates(theMeshToInflate,globalNodeID_surr,V);
            double lx = -(P_star[0] - V[0]);
            double ly = -(P_star[1] - V[1]);
            double lz = -(P_star[2] - V[2]);
            double dot = (lx*normal[0]+ly*normal[1]*lz*normal[2]);
            if(dot > 1) dot = 1;
            if(dot < -1) dot = -1;
            double angle = std::acos(dot);
            if(angle>maxAngle) maxAngle = angle;
            if(angle<minAngle) minAngle = angle;
        }
        if(maxAngle<70.0*PI/180.0)
        {
            myLayerHCutOff.insert(globalNodeID,0);
            blockedByAngle++;
            continue;   // jump over the subsequent check, if any
        }

        // another check, if any ...
    }

    //! ------------------------------------------------------------------------------------------
    //! 3. check lateral displacements - lateral smoothing can result in additional Cliff points?
    //!    if yes this check should be placed before "checkIncompleteManifold"
    //! ------------------------------------------------------------------------------------------
    //this->checkLateralDistributionMarchingDistance(theMeshToInflate,marchingDistanceMap);

    //! -----------------------------------------------
    //! 4. check "Cliff points" - incomplete manifolds
    //! -----------------------------------------------
    this->checkIncompleteManifold(theMeshToInflate);

    //! -----------------------------
    //! build the displacement field
    //! -----------------------------
    QMap<int,QList<double>> displacementsField;
    for(QMap<int,QList<double>>::const_iterator it = normals.cbegin(); it!= normals.cend(); ++it)
    {
        int globalNodeID = it.key();
        const QList<double> &curNormal = it.value();

        //! -------------------------
        //! jump over boundary nodes
        //! -------------------------
        if(theMeshToInflate->myBoundaryPoints.contains(globalNodeID)) continue;

        //! --------------------
        //! nodal displacements
        //! --------------------
        double marchingDistance = marchingDistanceMap.value(globalNodeID);
        double vx = -curNormal[0]*marchingDistance;
        double vy = -curNormal[1]*marchingDistance;
        double vz = -curNormal[2]*marchingDistance;

        //! -------------------
        //! displacement field
        //! -------------------
        QList<double> localFieldValue; localFieldValue<<vx<<vy<<vz;
        displacementsField.insert(globalNodeID,localFieldValue);
    }

    //! -------------------------
    //! scan the advancing nodes
    //! -------------------------
    for(TColStd_MapIteratorOfPackedMapOfInteger it(theMeshToInflate->GetAllNodes()); it.More(); it.Next())
    {
        int globalNodeID = it.Key();
        if(theMeshToInflate->myBoundaryPoints.contains(globalNodeID)) continue;

        //! ------------------------------
        //! prepare the "displaced point"
        //! ------------------------------
        int localNodeID = theMeshToInflate->myNodesMap.FindIndex(globalNodeID);
        const std::vector<double> &pc = theMeshToInflate->getNodeCoordinates(localNodeID);
        const QList<double> &vecDispl = displacementsField.value(globalNodeID);

        double x_displ = pc[0]+vecDispl[0];
        double y_displ = pc[1]+vecDispl[1];
        double z_displ = pc[2]+vecDispl[2];

        mesh::meshPoint mp_shifted(x_displ,y_displ,z_displ,globalNodeID);

        //! -----------------------------------------------------------
        //! retrieve the surface elements attached to the current node
        //! -----------------------------------------------------------
        const QList<int> &attachedElements = theMeshToInflate->myNodeToElements.value(localNodeID);
        int NbAttachedElements = attachedElements.size();

        //! -----------------------------------------------------------------------
        //! iterate over the surface elements attached to the current node
        //! foreach attached triangle, build a volume element of the boundary mesh
        //! -----------------------------------------------------------------------
        std::vector<meshElementByCoords> vecElementsToAdd;
        bool canMove = true;
        for(int i=0; i<NbAttachedElements; i++)
        {
            int localElementID = attachedElements[i];
            int globalElementID = theMeshToInflate->myElementsMap.FindKey(localElementID);

            //! -------------------------------------------------------------------
            //! a boundary volume element - the volume element ID is not specified
            //! -------------------------------------------------------------------
            meshElementByCoords aVolElement;
            aVolElement.type = TET;
            aVolElement.ID = -1;

            int NbNodes, buf[3];
            TColStd_Array1OfInteger nodeIDs(*buf,1,3);
            theMeshToInflate->GetNodesByElement(globalElementID,nodeIDs,NbNodes);

            for(int k=1; k<=NbNodes; k++)
            {
                int globalNodeID_ = nodeIDs(k);
                int localNodeID_ = theMeshToInflate->myNodesMap.FindIndex(globalNodeID_);
                std::vector<double> &pc_ = theMeshToInflate->getNodeCoordinates(localNodeID_);
                mesh::meshPoint mp_(pc_[0],pc_[1],pc_[2],globalNodeID_);
                aVolElement.pointList<<mp_;
            }

            //! ------------------------------------------
            //! complete the volume element - added after
            //! ------------------------------------------
            aVolElement.pointList<<mp_shifted;

            //! ------------------------
            //! check inverted elements
            //! ------------------------
            TetQualityClass aTetQuality;
            aTetQuality.setPoints(aVolElement.getPoints());

            double V = aTetQuality.Volume();
            if(V<=0.0)
            {
                NbInvalidTets++;
                canMove = false;
                break;
            }

            //! -------------------------------------------------------
            //! check if the volume element is valid, and can be added
            //! - bypass a full element quality check -
            //! -------------------------------------------------------
            //double q0, q1, q2;
            //aTetQuality.setPoints(aVolElement.getPoints());
            //aTetQuality.getQualityMeasure(q0,q1,q2,V);

            //const double QUALITY_LIMIT = 0.01;
            //if(q2<QUALITY_LIMIT)
            //{
            //    cout<<"____CANNOT MOVE____"<<endl;
            //    NbInvalidTets++;
            //    canMove = false;
            //    break;
            //}
            vecElementsToAdd.push_back(aVolElement);
        }
        if(canMove)
        {
            std::vector<double> P {x_displ,y_displ,z_displ};
            theMeshToInflate->changeNodeCoords(localNodeID,P);
            int localNodeID_surfaceMesh = myOverallSumMeshDS->myNodesMap.FindIndex(globalNodeID);

            //! ------------------------------------------------------------------------------------
            //! copy the current advancing mesh: needed by intersection and self intersection algos
            //! ------------------------------------------------------------------------------------
            occHandle(Ng_MeshVS_DataSourceFace) myOverallSumMeshDS_old(myOverallSumMeshDS);

            //! ----------------------------------
            //! move the node of the current mesh
            //! ----------------------------------
            myOverallSumMeshDS->changeNodeCoords(localNodeID_surfaceMesh,P);

            //! -------------------
            //! standard treatment
            //! -------------------
            for(int i=0; i<vecElementsToAdd.size(); i++)
            {
                volumeElementsAtWalls.push_back(vecElementsToAdd[i]);
            }
        }
        else
        {
            //! lock the node for all the subsequent steps
            myLayerHCutOff.insert(globalNodeID,0);
        }
    }
    cout<<"Number of inverted tet: "<<NbInvalidTets<<endl;
    cout<<"Number of shift blocked by lenght criterion: "<<blockedByLenght<<endl;
    cout<<"Number of shift blocked by angle criterion: "<<blockedByAngle<<endl;
}

//! ---------------------------------------------------
//! function: checkLateralDistributionMarchingDistance
//! details:
//! ---------------------------------------------------
void prismaticLayer::checkLateralDistributionMarchingDistance(const occHandle(Ng_MeshVS_DataSourceFace) &aMeshDS,
                                                              QMap<int,double> &marchingDistanceMap)
{
    cout<<"prismaticLayer::checkLateralDistributionMarchingDistance()->____function called____"<<endl;

    QMap<int,double> marchingDistanceMapUpdated;
    for(QMap<int,double>::iterator it = marchingDistanceMap.begin(); it!=marchingDistanceMap.end(); it++)
    {
        int globalNodeID = it.key();
        double distanceToCheck = it.value();
        marchingDistanceMapUpdated.insert(globalNodeID,distanceToCheck);

        //! ----------------------------------------------------------------------
        //! among the surrounding nodes, check the ones which do not satisfy
        //! the inner inequality, and store them into the "referenceNodes" vector
        //! ----------------------------------------------------------------------
        bool foundNodesAround = false;
        const std::vector<int> &surroundingNodes = mySurroudingNodesMap.at(globalNodeID);
        if(surroundingNodes.empty()) exit(9999);
        std::vector<int> referenceNodes;
        for(std::vector<int>::const_iterator it_ = surroundingNodes.cbegin(); it_!=surroundingNodes.cend(); it_++)
        {
            int globalNodeID_surrounding = *it_;
            double curMarchingDistance_surr = marchingDistanceMap.value(globalNodeID_surrounding);
            if(distanceToCheck<0.5*curMarchingDistance_surr || distanceToCheck>2*curMarchingDistance_surr)
            {
                referenceNodes.push_back(globalNodeID_surrounding);
                foundNodesAround = true;
                //distanceToCheck = (5.0/4.0)*curMarchingDistance_surr;
            }
        }
        if(foundNodesAround == true)
        {
            double P[3];
            this->getPointCoordinates(aMeshDS,globalNodeID,P);

            /*
            double dtot_inv = 0;
            double Spartial = 0;
            */

            double dtot = 0;
            double Spartial = 0;

            for(std::vector<int>::iterator it_ = referenceNodes.begin(); it_ != referenceNodes.end(); it_++)
            {
                int nodeID = *it_;
                double B[3];
                this->getPointCoordinates(aMeshDS,nodeID,B);

                /*
                double d_inv = 1/(sqrt(pow(P[0]-B[0],2)+pow(P[1]-B[1],2)+pow(P[2]-B[2],2)));
                Spartial += d_inv * marchingDistanceMap.value(nodeID);
                dtot_inv += 1/d_inv;
                */

                double d = sqrt(pow(P[0]-B[0],2)+pow(P[1]-B[1],2)+pow(P[2]-B[2],2));
                Spartial += d * marchingDistanceMap.value(nodeID);
                dtot += d;
            }
            //double smoothedMarchingDistance = Spartial/dtot_inv;
            double smoothedMarchingDistance = Spartial/dtot;
            marchingDistanceMapUpdated.insert(globalNodeID,smoothedMarchingDistance);
        }
    }

    //! ----------------------------------------------------------------------------------------
    //! update the marching distances - optimize possibile updating only the modified distances
    //! ----------------------------------------------------------------------------------------
    for(QMap<int,double>::iterator it = marchingDistanceMapUpdated.begin(); it != marchingDistanceMapUpdated.end(); it++)
    {
        //cout<<"____"<<it.key()<<"\t"<<marchingDistanceMap.value(it.key())<<"\t"<<it.value()<<"____"<<endl;
        marchingDistanceMap.insert(it.key(), it.value());
    }
}

//! ------------------------------------------------------------------------------
//! function: classifyNodes
//! details:  1 - fixed  2 - direction only 3 - marching distance only - 4 - both
//! ------------------------------------------------------------------------------
void prismaticLayer::classifyNodes(const occHandle(Ng_MeshVS_DataSourceFace) &aMeshDS,
                                   std::map<int,int> &mapCat)
{
    cout<<"prismaticLayer::classifyNodes()->____function called____"<<endl;

    const double PI = 3.145926534;
    const double angleLimit_cat1 = PI/3.0;

    //! --------------
    //! vector of int
    //! --------------
    std::vector<int> temp;

    if(aMeshDS->myNodeToElements.isEmpty()) aMeshDS->computeNodeToElementsConnectivity();

    //cout<<"prismaticLayer::classifyNodes()->____searching for nodes of type 1____"<<endl;
    int NbNodesOfType1 = 0;
    for(TColStd_MapIteratorOfPackedMapOfInteger it(aMeshDS->GetAllNodes()); it.More(); it.Next())
    {
        int globalNodeID = it.Key();
        int localNodeID = aMeshDS->myNodesMap.FindIndex(globalNodeID);
        const QList<int> &surroundingElements_local = aMeshDS->myNodeToElements.value(localNodeID);
        int NbSurroundingElements = surroundingElements_local.length();

        bool isOneAngleGreaterThanLimit = false;
        for(int i=0; i<NbSurroundingElements-1; i++)
        {
            int localElementID_surr_i = surroundingElements_local[i];
            int globalElementID_surr_i = aMeshDS->myElementsMap.FindKey(localElementID_surr_i);
            double nxi,nyi,nzi;
            aMeshDS->GetNormal(globalElementID_surr_i,10,nxi,nyi,nzi);

            for(int j=i+1; j<NbSurroundingElements; j++)
            {
                int localElementID_surr_j = surroundingElements_local[j];
                int globalElementID_surr_j = aMeshDS->myElementsMap.FindKey(localElementID_surr_j);
                double nxj,nyj,nzj;
                aMeshDS->GetNormal(globalElementID_surr_j,10,nxj,nyj,nzj);
                double dot = nxi*nxj+nyi*nyj+nzi*nzj;
                if(dot<-1) dot = -1;
                if(dot>1) dot = 1;
                double angle = std::acos(dot);
                if(fabs(angle)>=angleLimit_cat1)
                {
                    isOneAngleGreaterThanLimit = true;
                    break;
                }
                if(isOneAngleGreaterThanLimit == true) break;
            }
        }
        if(isOneAngleGreaterThanLimit == true)
        {
            //cout<<"____node ID: "<<globalNodeID<<" is of type 1____"<<endl;
            mapCat.insert(std::make_pair(globalNodeID,1));
            mapOfReductionFactor.insert(std::make_pair(globalNodeID,1));
            NbNodesOfType1++;
        }
        else
        {
            mapCat.insert(std::make_pair(globalNodeID,4));
            mapOfReductionFactor.insert(std::make_pair(globalNodeID,4));
            temp.push_back(globalNodeID);
        }
    }
    //cout<<"prismaticLayer::classifyNodes()->____number of nodes of type 1: "<<NbNodesOfType1<<"____"<<endl;

    //! ----------------------------------------------------------------
    //! iterate all over the points =>not<= belonging to category 1
    //! one of this points, if on the prismatic/non prismatic boundary,
    //! belongs to category 2
    //! ----------------------------------------------------------------
    //cout<<"prismaticLayer::classifyNodes()->____residual nodes: "<<temp.size()<<"____"<<endl;
    //cout<<"prismaticLayer::classifyNodes()->____searching for nodes of type 2____"<<endl;
    int NbNodesOfType2 = 0;
    for(std::vector<int>::iterator it = temp.begin(); it!=temp.end();)
    {
        int globalNodeID = *it;
        if(aMeshDS->myBoundaryPoints.contains(globalNodeID))
        {
            std::map<int,int>::iterator itt = mapCat.find(globalNodeID);
            itt->second = 2;
            it = temp.erase(it);
            NbNodesOfType2++;
        }
        else it++;
    }

    //cout<<"prismaticLayer::classifyNodes()->____number of nodes of type 2: "<<NbNodesOfType2<<"____"<<endl;

    //! -----------------------------------------------------------------------
    //! for a given mesh point, check if at least one of the surrounding nodes
    //! is in category 1 (fixed)
    //! -----------------------------------------------------------------------
    //cout<<"prismaticLayer::classifyNodes()->____residual nodes: "<<temp.size()<<"____"<<endl;
    //cout<<"prismaticLayer::classifyNodes()->____searching for nodes of type 3____"<<endl;
    int NbNodesOfType3 = 0;
    for(std::vector<int>::iterator it = temp.begin(); it!=temp.end();)
    {
        int globalNodeID = *it;

        //! ------------------
        //! surrounding nodes
        //! ------------------
        int localNodeID = aMeshDS->myNodesMap.FindIndex(globalNodeID);
        const QList<int> &surroundingElements_local = aMeshDS->myNodeToElements.value(localNodeID);
        int NbSurroundingElements = surroundingElements_local.length();
        std::set<int> setSurroundingNodes_globalIDs;
        for(int i=0; i<NbSurroundingElements; i++)
        {
            int localElementID_surr = surroundingElements_local[i];
            int globalElementID_surr = aMeshDS->myElementsMap.FindKey(localElementID_surr);
            int buf[12], NbNodes_;
            TColStd_Array1OfInteger nodeIDs(*buf,1,12);
            aMeshDS->GetNodesByElement(globalElementID_surr,nodeIDs,NbNodes_);
            for(int n=1; n<=NbNodes_; n++) setSurroundingNodes_globalIDs.insert(nodeIDs(n));
        }
        std::set<int>::iterator it__ = setSurroundingNodes_globalIDs.find(globalNodeID);
        setSurroundingNodes_globalIDs.erase(it__);

        bool aFixedNodeHasBeeFound = false;
        for(std::set<int>::iterator it___ = setSurroundingNodes_globalIDs.begin(); it___ != setSurroundingNodes_globalIDs.end(); it___++)
        {
            int curGlobalNodeID = *it___;
            std::map<int,int>::iterator it____ = mapCat.find(curGlobalNodeID);
            if(it____->second == 1)
            {
                aFixedNodeHasBeeFound = true;
                break;
            }
        }
        if(aFixedNodeHasBeeFound == true)
        {
            mapCat.find(globalNodeID)->second = 3;
            mapOfReductionFactor.find(globalNodeID)->second = 3;
            it = temp.erase(it);
            NbNodesOfType3++;
        }
        else it++;
    }

    cout<<"prismaticLayer::classifyNodes()->____number of nodes of type 3: "<<NbNodesOfType3<<"____"<<endl;
    cout<<"prismaticLayer::classifyNodes()->____number of nodes of type 4: "<<temp.size()<<"____"<<endl;

    //! -----------------------------------------------------------------------------
    //! this section modifies the map of node types: mesh nodes having at least one
    //! neighbor which has "compression" due to wall gap check are put into cat "5":
    //! their marching distance is a weighted average of the marching distance of
    //! the surrounding nodes. Note: a point of category "2" cannot be a point of
    //! category "5"
    //! -----------------------------------------------------------------------------
    for(std::map<int,int>::iterator it = mapCat.begin(); it!= mapCat.end(); it++)
    {
        int globalNodeID = it->first;
        if(mapOfReductionFactor.at(globalNodeID) == 2) continue;    // jump over nodes of category "2"

        //! ----------------------
        //! get surrounding nodes
        //! ----------------------
        int localNodeID = aMeshDS->myNodesMap.FindIndex(globalNodeID);
        const QList<int> &surroundingElements_local = aMeshDS->myNodeToElements.value(localNodeID);

        int NbSurroundingElements = surroundingElements_local.length();
        std::set<int> setSurroundingNodes_globalIDs;
        for(int i=0; i<NbSurroundingElements; i++)
        {
            int localElementID_surr = surroundingElements_local[i];
            int globalElementID_surr = aMeshDS->myElementsMap.FindKey(localElementID_surr);
            int buf[12], NbNodes_;
            TColStd_Array1OfInteger nodeIDs(*buf,1,12);
            aMeshDS->GetNodesByElement(globalElementID_surr,nodeIDs,NbNodes_);
            for(int n=1; n<=NbNodes_; n++) setSurroundingNodes_globalIDs.insert(nodeIDs(n));
        }

        std::set<int>::iterator it__ = setSurroundingNodes_globalIDs.find(globalNodeID);
        setSurroundingNodes_globalIDs.erase(it__);

        //! ------------------------------------------------------------------
        //! check the surrounding nodes: it at least one of the surrounding
        //! nodes has a correction due to wall proximity put the current node
        //! in category 5
        //! ------------------------------------------------------------------
        bool hasCompression = false;
        for(std::set<int>::iterator it = setSurroundingNodes_globalIDs.begin(); it != setSurroundingNodes_globalIDs.end(); it++)
        {
            int globalNodeID_surr = *it;
            if(mapOfReductionFactor.at(globalNodeID_surr)<1.00)
            {
                hasCompression = true;
                break;
            }
        }
        if(hasCompression == true)
        {
            mapOfReductionFactor.find(globalNodeID)->second = 4;
            mapCat.find(globalNodeID)->second = 4;
        }
    }
}

//! ----------------------------------
//! function: checkIncompleteManifold
//! details:
//! ----------------------------------
void prismaticLayer::checkIncompleteManifold(const occHandle(Ng_MeshVS_DataSourceFace) &theMeshToInflate)
{
    cout<<"prismaticLayer::checkIncompleteManifold()->____function called____"<<endl;
    int NbIncompleteManifolds = 0;
    QMap<int,double> layerHCutoff;
    for(TColStd_MapIteratorOfPackedMapOfInteger it = theMeshToInflate->GetAllNodes(); it.More(); it.Next())
    {
        int globalNodeID = it.Key();
        std::vector<int> surroundingNodes = mySurroudingNodesMap.at(globalNodeID);

        bool foundFixedPointAround = false;
        for(std::vector<int>::iterator it_ = surroundingNodes.begin(); it_!=surroundingNodes.end(); it_++)
        {
            int globalNodeID_surrounding = *it_;
            double cutOff = myLayerHCutOff.value(globalNodeID_surrounding);
            if(cutOff==0)
            {
                foundFixedPointAround = true;
                NbIncompleteManifolds++;
                break;
            }
        }
        double newCutOff;
        foundFixedPointAround==true? newCutOff = 0.0: newCutOff = 1.0;
        layerHCutoff.insert(globalNodeID,newCutOff);
    }
    //! --------------------------
    //! update the map of cutoffs
    //! --------------------------
    for(QMap<int,double>::iterator it = layerHCutoff.begin(); it!= layerHCutoff.end(); it++)
    {
        myLayerHCutOff.insert(it.key(),it.value());
    }
    cout<<"prismaticLayer::checkIncompleteManifold()->____number of incomplete manifolds: "<<NbIncompleteManifolds<<"____"<<endl;
}

//! -----------------------------------
//! function: buildSurroundingNodesMap
//! details:
//! -----------------------------------
void prismaticLayer::buildSurroundingNodesMap(const occHandle(Ng_MeshVS_DataSourceFace) &aMeshDS, std::map<int,std::vector<int>> &surroundingNodesMap)
{
    cout<<"prismaticLayer::buildSurroundingNodesMap()->____function called____"<<endl;

    if(aMeshDS->myNodeToElements.isEmpty()) aMeshDS->computeNodeToElementsConnectivity();
    const QMap<int,QList<int>> &nodeToElementMap = aMeshDS->myNodeToElements;
    for(TColStd_MapIteratorOfPackedMapOfInteger it = aMeshDS->GetAllNodes(); it.More(); it.Next())
    {
        int globalNodeID = it.Key();
        int localNodeID = aMeshDS->myNodesMap.FindIndex(globalNodeID);
        const QList<int> &surroundingElements_local = nodeToElementMap.value(localNodeID);
        int NbSurroundingElements = surroundingElements_local.length();
        std::set<int> setOfSurroundingNodes;
        for(int i = 0; i<NbSurroundingElements; i++)
        {
            int localElementID_surrounding = surroundingElements_local[i];
            int globalElemenID_surrounding = aMeshDS->myElementsMap.FindKey(localElementID_surrounding);
            int NbNodes, buf[20];
            TColStd_Array1OfInteger nodeIDs(*buf,1,20);
            aMeshDS->GetNodesByElement(globalElemenID_surrounding,nodeIDs,NbNodes);
            for(int n=1; n<=NbNodes; n++) setOfSurroundingNodes.insert(nodeIDs(n));
        }
        std::set<int>::iterator it__ = setOfSurroundingNodes.find(globalNodeID);
        setOfSurroundingNodes.erase(it__);
        std::vector<int> vecSurroundingNodes;
        for(std::set<int>::iterator it__ = setOfSurroundingNodes.begin(); it__!=setOfSurroundingNodes.end(); it__++)
            vecSurroundingNodes.push_back(*it__);
        surroundingNodesMap.insert(std::make_pair(globalNodeID,vecSurroundingNodes));
    }
    cout<<"prismaticLayer::buildSurroundingNodesMap()->____surrounding nodes map built____"<<endl;
}
