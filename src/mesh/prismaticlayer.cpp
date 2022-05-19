//! ----------------
//! custom includes
//! ----------------
#include "prismaticlayer.h"
#include "src/utils/meshtools.h"
#include "src/utils/ccout.h"
#include "src/utils/geomtoolsclass.h"
#include "igtools.h"
#include <ng_meshvs_datasourceface.h>
#include "src/utils/meshtools.h"
#include "src/utils/global.h"
#include <tetqualityclass.h>
#include <meshface.h>
#include <OCCface.h>
#include <meshselectnlayers.h>
#include <igtools.h>
#include <smoothingtools.h>
#include <pointtomeshdistance.h>
#include <mesh.h>

//! ------
//! Eigen
//! ------
#include <ext/eigen/Eigen/Dense>

//! ----------------------
//! used by "displayMesh"
//! ----------------------
#include <occPreGLwidget.h>
#include <AIS_InteractiveContext.hxx>
#include <MeshVS_MeshPrsBuilder.hxx>
#include <MeshVS_Drawer.hxx>
#include <MeshVS_DrawerAttribute.hxx>
#include "ext/occ_extended/arrayofcolors.h"
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
    cout<<"prismaticLayer::computeVecFieldCutOff()->____function called____"<<endl;
    myLayerHCutOff.clear();

    //! ------------------------------------------------------------------------------
    //! this is very slow on big meshes - done internally within the mesh constructor
    //! from a list of face mesh data sources! Left here for documentation
    //! ------------------------------------------------------------------------------
    //myPrismaticFacesSumMeshDS->computeFreeMeshSegments();

    /*
    //! --------------------------------------------
    //! the following list contains global node IDs
    //! --------------------------------------------
    QList<int> nodeIDsOnPrismaticFaces1DBoundary = myPrismaticFacesSumMeshDS->myBoundaryPoints;

    cout<<"prismaticLayer::computeVecFieldCutOff()->____number of points on the prismatic faces boundary: "<<nodeIDsOnPrismaticFaces1DBoundary.length()<<"____"<<endl;

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
    */

    for(TColStd_MapIteratorOfPackedMapOfInteger it(myOverallSumMeshDS->GetAllNodes()); it.More(); it.Next())
    {
        myLayerHCutOff.insert(it.Key(),0);
    }
    for(TColStd_MapIteratorOfPackedMapOfInteger it(myPrismaticFacesSumMeshDS->GetAllNodes()); it.More(); it.Next())
    {
        myLayerHCutOff.insert(it.Key(),1);
    }
    //cout<<"prismaticLayer::computeVecFieldCutOff()->____number of points on the prismatic faces boundary: "<<myPrismaticFacesSumMeshDS->myBoundaryPoints.length()<<"____"<<endl;
    if(lockBoundary==true)
    {
        const QList<int> &bps = myPrismaticFacesSumMeshDS->myBoundaryPoints;
        for(int i=0; i<bps.length(); i++) myLayerHCutOff.insert(bps[i],0);
    }
    cout<<"prismaticLayer::computeVecFieldCutOff()->____exiting function____"<<endl;
    /* for(QMap<int,double>::iterator it = myLayerHCutOff.begin(); it!= myLayerHCutOff.end(); ++it)
    {
        cout<<"prismaticLayer::myLayerHCutOff()->____"<<it.key()<<"____"<<it.value()<<endl;
    }*/
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
        float gap;
        this->getPointCoordinates(theMeshToInflate_old,globalNodeID,P);

        //! --------------------------------------------------
        //! invert the normal: it points towards the material
        //! for gap calculation
        //! --------------------------------------------------
        double nx = -normal[0];
        double ny = -normal[1];
        double nz = -normal[2];
        double dir[3] {nx,ny,nz};
        aDistanceMeter.distance(P,dir,&gap);   // compute the distance

        bool compress = true;
        if(compress==true)
        {
            //! ---------------------------------------------
            //! compress layers in order to avoid collisions
            //! ---------------------------------------------
            const double compressionFactor = 0.25;
            if(myTotalThickness<gap/4.0)
            {
                //! the node can move with the nominal marching distance
                mapOfFirstLayerReductionFactor.insert(std::make_pair(globalNodeID,1.0));
            }
            else
            {
                //! the node can move with smaller marching distance
                double firstLayerThicknessNew = compressionFactor*gap/Sigma;
                double reductionFactor = firstLayerThicknessNew/firstLayerThickness;
                mapOfFirstLayerReductionFactor.insert(std::make_pair(globalNodeID,reductionFactor));
                mapOfReductionFactor.insert(std::make_pair(globalNodeID,reductionFactor));
            }
        }
        else
        {
            if(myTotalThickness<gap/4.0)
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

        //! -------------------------------------------------
        //! check "Cliff points" - incomplete manifolds
        //! these points could be generated at the beginning
        //! if the option "compress" is true
        //! -------------------------------------------------
        this->checkIncompleteManifold(theMeshToInflate_old);

        //! ----------------------------
        //! the current layer thickness
        //! ----------------------------
        displacement = myLayerThickness.at(n-1);

        //! -------------------
        //! clone the old mesh
        //! -------------------
        occHandle(Ng_MeshVS_DataSourceFace) theMeshToInflate_new = new Ng_MeshVS_DataSourceFace(theMeshToInflate_old);

        //! ---------------------
        //! "best" nodal normals
        //! ---------------------
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
        smoothingTools::fieldSmoother(normals,theMeshToInflate_new,betaAverageField,betaVisibilityField,mapNodeTypes,10,true,n);

        //! ------------------------------------
        //! map of the local marching distances
        //! ------------------------------------
        QMap<int,double> marchingDistanceMap;
        for(QMap<int,QList<double>>::const_iterator it = normals.cbegin(); it!= normals.cend(); ++it)
        {
            int globalNodeID = it.key();
            //cout<<"NodeID "<<globalNodeID<<endl;
            double shrinkFactor = shrinkFactors.value(globalNodeID);
            double cutOff = myLayerHCutOff.value(globalNodeID);
            //cout<<"cutOff "<<cutOff<<" shrink "<<shrinkFactor<<endl;
            double marchingDistance = displacement*(1+shrinkFactor)*cutOff*mapOfFirstLayerReductionFactor.at(globalNodeID);
            marchingDistanceMap.insert(globalNodeID,marchingDistance);
        }

        //! -------------------------------------------------
        //! check lateral distribution of marching distances
        //! -------------------------------------------------
        this->checkLateralDistributionMarchingDistance(theMeshToInflate_new,n,marchingDistanceMap,normals);

        //! ------------------------------
        //! smooth the marching distances
        //! ------------------------------
        smoothingTools::scalarFieldSmoother(marchingDistanceMap,theMeshToInflate_new,betaAverageField,mapNodeTypes,10,n);

        //! --------------
        //! analyze front
        //! --------------
        this->analyzeFront(theMeshToInflate_new,normals,marchingDistanceMap);

        //! --------------------------------------------
        //! check gaps along slightly different normals
        //! --------------------------------------------
        const double PI = 3.1415926538;
        pointToMeshDistance aDistanceMeter;
        aDistanceMeter.init(theMeshToInflate_new);

        for(TColStd_MapIteratorOfPackedMapOfInteger it(theMeshToInflate_new->GetAllNodes()); it.More(); it.Next())
        {
            int globalNodeID = it.Key();
            double P[3];
            this->getPointCoordinates(theMeshToInflate_new,globalNodeID,P);
            float gap;
            QList<double> normal = theMeshToInflate_new->myNodeNormals.value(globalNodeID);
            double dir[3] { -normal[0], -normal[1], -normal[2] };
            double angle = 5.0*PI/180.0;
            aDistanceMeter.distance(P,dir,angle,8,&gap);
            if(marchingDistanceMap.value(globalNodeID)>0.30*gap)
            {
                myLayerHCutOff.insert(globalNodeID,0);
                cout<<"prismaticLayer::inflateMesh()->____blocked node ID: "<<globalNodeID<<"____displacement: "<<displacement<<"\t gap: "<<gap<<"____"<<endl;
            }
        }

        //! ------------------------
        //! check self intersection
        //! ------------------------
        if(myCheckSelfIntersections==true) this->checkSelfIntersection(theMeshToInflate_new,theMeshToInflate_old,true);


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
            double cutOff = myLayerHCutOff.value(globalNodeID);

            //cout<< "marching "<<marchingDistance<<" , cutOFF "<<cutOff<<" , normal "<<curNormal[0]<<endl;

            double vx = -curNormal[0]*marchingDistance*cutOff;
            double vy = -curNormal[1]*marchingDistance*cutOff;
            double vz = -curNormal[2]*marchingDistance*cutOff;
            //cout<< "displacement field "<<vx<<" , "<<vy<<" , "<<vz<<endl;
            QList<double> localFieldValue;
            localFieldValue<<vx<<vy<<vz;
            displacementsField.insert(globalNodeID,localFieldValue);
        }

        //! ------------------
        //! displace the mesh
        //! ------------------
        theMeshToInflate_new->displaceMySelf(displacementsField);

        //! ---------------------------------------------------------
        //! start correcting nodes coordinates if the inflated layer
        //! has an intersection with the mesh it originates from
        //! ---------------------------------------------------------
        if(myCheckMutualIntersections==true) this->checkMutualIntersection(theMeshToInflate_new,theMeshToInflate_old,true);
        //to do -> action for update the displaced mesh

        //! ----------------------------------------------
        //! add the displaced mesh to the list of results
        //! ----------------------------------------------
        inflatedMeshes<<theMeshToInflate_new;

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
        const occHandle(MeshVS_DataSource) &aMeshDS = myMeshDB->ArrayOfMeshDSOnFaces.getValue(bodyIndex,faceNr);
        if(aMeshDS.IsNull())
        {
            //cout<<"____the face mesh data source of face nr: "<<faceNr<<" is NULL____"<<endl;
            continue;
        }
        const occHandle(Ng_MeshVS_DataSourceFace) &aFaceMeshDS = occHandle(Ng_MeshVS_DataSourceFace)::DownCast(aMeshDS);

        faceList0<<aFaceMeshDS;
        if(std::find(myPrismaticFaces.begin(), myPrismaticFaces.end(), faceNr)!=myPrismaticFaces.end()) faceList1<<aFaceMeshDS;
        if(std::find(myPrismaticFaces.begin(), myPrismaticFaces.end(), faceNr)==myPrismaticFaces.end()) faceList2<<aFaceMeshDS;
    }

    cout<<"____number of faces: "<<faceList0.length()<<"____"<<endl;
    cout<<"____number of wall faces: "<<faceList1.length()<<"____"<<endl;
    cout<<"____number of non-wall faces: "<<faceList2.length()<<"____"<<endl;

    myOverallSumMeshDS = new Ng_MeshVS_DataSourceFace(faceList0);               //! surface mesh
    myPrismaticFacesSumMeshDS = new Ng_MeshVS_DataSourceFace(faceList1);        //! prismatic/wall mesh
    myNonPrismaticFacesSumMeshDS = new Ng_MeshVS_DataSourceFace(faceList2);     //! non prismatic/non wall mesh

    myOverallSumMeshDS->computeNormalAtNodes();
    myPrismaticFacesSumMeshDS->computeNormalAtNodes();
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

    if(l1 * l2 == 0)
    {
        cerr<<"angleBetweenVector()->____abnormal termination at line 634____"<<endl;
        exit(9999);    // questo non dovrebbe mai succedere
    }

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
        if(betaVisibility<=0)
        {
            cerr<<"prismaticLayer::computeBeta()->____node ID: "<<globalNodeID<<" negative or zero visibility angle. Value: "
               <<betaVisibility<<"____"<<endl;

            //! --------------
            //! lock the node
            //! --------------
            myLayerHCutOff.insert(globalNodeID,0);
            betaVisibilityField.insert(std::make_pair(globalNodeID,0.0));
        }
        else
        {
            betaVisibilityField.insert(std::make_pair(globalNodeID,betaVisibility));
        }
        //cout<<"____betaVisibility (rad,deg) = "<<betaVisibility<<"\t"<<betaVisibility*180/PI<<"____"<<endl;

        std::set<int>::iterator it_ = surroundingNodes.find(globalNodeID);
        surroundingNodes.erase(it_);

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
            if(sqrt(n1[0]*n1[0]+n1[1]*n1[1]+n1[2]*n1[2])==0)
            {
                //! lock the node
                myLayerHCutOff.insert(globalNodeID,0);
                continue;
            }
            if(sqrt(n2[0]*n2[0]+n2[1]*n2[1]+n2[2]*n2[2])==0)
            {
                //! lock the node
                myLayerHCutOff.insert(globalNodeID,0);
                continue;
            }

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
        if(betaAve>PI+eps) shrink =  -fabs(cos(betaVisibility));    // retraction in convex points
        shrinkFactors.insert(globalNodeID,shrink);

        //double P[3];
        //this->getPointCoordinates(aMeshDS,globalNodeID,P);
        //fprintf(fp,"%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",globalNodeID,P[0],P[1],P[2],betaAve*180/PI,betaVisibility*180/PI,shrink);
    }
    //fclose(fp);
}

//! ------------------------------------------------------------
//! function: buildPrismaticElements
//! details:  build an hydrib volume mesh using inflated meshes
//! ------------------------------------------------------------
bool prismaticLayer::buildPrismaticElements(const QList<occHandle(Ng_MeshVS_DataSourceFace)> &theInflatedMeshes,
                                            occHandle(Ng_MeshVS_DataSource3D) &prismaticMeshDS3D)
{
    cout<<"prismaticLayer::buildPrismaticElements()->____function called____"<<endl;

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
//! function: inflateMeshAndCompress
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

    //! ----------------------
    //! surrounding nodes map
    //! ----------------------
    mySurroudingNodesMap.clear();
    this->buildSurroundingNodesMap(theMeshToInflate_old,mySurroudingNodesMap);

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

    //! ------------------------------------------------------------------
    //! analyze gaps at the beginning: scan all the mesh nodes.
    //! compute the distance from closest point on the mesh => "distance"
    //! apply reduction of first layer thickness if applicable
    //! ------------------------------------------------------------------
    double Sigma = 0;                                           // a parameter for reduction factor calculation
    for(int n=0; n<myNbLayers; n++) Sigma += pow(myExpRatio,n);

    pointToMeshDistance aDistanceMeter;                         // distance meter tool
    aDistanceMeter.init(theMeshToInflate_old);                  // init with mesh
    double firstLayerThickness = myLayerThickness.at(0);        // first layer thickness
    std::map<int,double> mapOfFirstLayerReductionFactor;        // fill a map (nodeID, reduction factor)

    if(theMeshToInflate_old->myNodeNormals.isEmpty()) theMeshToInflate_old->computeNormalAtNodes();
    const QMap<int,QList<double>> &normals = theMeshToInflate_old->myNodeNormals;

    //! ---------------------------------------------
    //! preliminary analysis of gaps
    //! start filling the map of compression factors
    //! ---------------------------------------------
    for(QMap<int,QList<double>>::const_iterator itn = normals.cbegin(); itn!=normals.cend(); itn++)
    {
        int globalNodeID = itn.key();
        const QList<double> &normal = itn.value();

        double P[3];
        float gap;
        this->getPointCoordinates(theMeshToInflate_old,globalNodeID,P);

        //! --------------------------------------------------
        //! invert the normal: it points towards the material
        //! for gap calculation
        //! --------------------------------------------------
        double nx = -normal[0];
        double ny = -normal[1];
        double nz = -normal[2];
        double dir[3] {nx,ny,nz};
        aDistanceMeter.distance(P,dir,&gap);   // compute the distance

        bool compress = true;
        if(compress==true)
        {
            //! ---------------------------------------------
            //! compress layers in order to avoid collisions
            //! ---------------------------------------------
            const double compressionFactor = 0.25;
            if(myTotalThickness<gap/4.0)
            {
                //! the node can move with the nominal marching distance
                mapOfFirstLayerReductionFactor.insert(std::make_pair(globalNodeID,1.0));
            }
            else
            {
                //! the node can move with smaller marching distance
                double firstLayerThicknessNew = compressionFactor*gap/Sigma;
                double reductionFactor = firstLayerThicknessNew/firstLayerThickness;
                mapOfFirstLayerReductionFactor.insert(std::make_pair(globalNodeID,reductionFactor));
                mapOfReductionFactor.insert(std::make_pair(globalNodeID,reductionFactor));
            }
        }
        else
        {
            if(myTotalThickness<gap/4.0)
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
        smoothingTools::fieldSmoother(normals,theMeshToInflate_new,betaAverageField,betaVisibilityField,mapNodeTypes,10,true,n);

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
            //double cutoff = 1.0;
            //double marchingDistance = displacement*(1+shrinkFactor)*cutoff*mapOfReductionFactor.at(globalNodeID);
            double marchingDistance = displacement*(1+shrinkFactor)*mapOfReductionFactor.at(globalNodeID);
            marchingDistanceMap.insert(globalNodeID,marchingDistance);
        }

        //! ---------------------------------------------
        //! check marching distance at P wrt surrounding
        //! ---------------------------------------------
        this->checkLateralDistributionMarchingDistance(theMeshToInflate_new,n,marchingDistanceMap,normals);

        //! ------------------------------
        //! smooth the marching distances
        //! ------------------------------
        smoothingTools::scalarFieldSmoother(marchingDistanceMap,theMeshToInflate_new,betaAverageField,mapNodeTypes,10,n);

        //! --------------
        //! analyze front
        //! --------------
        this->analyzeFront(theMeshToInflate_new,normals,marchingDistanceMap);

        //! --------------
        //! analyze front
        //! --------------
        this->checkIncompleteManifold(theMeshToInflate_new);

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
            double cutOff = myLayerHCutOff.value(globalNodeID);
            double marchingDistance = marchingDistanceMap.value(globalNodeID)*mapOfFirstLayerReductionFactor.at(globalNodeID);
            double vx = -curNormal[0]*marchingDistance*cutOff;
            double vy = -curNormal[1]*marchingDistance*cutOff;
            double vz = -curNormal[2]*marchingDistance*cutOff;

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

        //! -------------------------------------------------------------
        //! translation of nodes without mesh compression: this was
        //! the initial attempt to compress, obviously leading to
        //! mesh element folding in case of large displacements or small
        //! element size - left here for documentation
        //! -------------------------------------------------------------
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
        //! 2 steps seems reasonable
        //! topological distance
        //! -------------------------
        meshSelector->setSteps(2);
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
        theMeshToInflate_old = new Ng_MeshVS_DataSourceFace(theMeshToInflate_new);
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

    const double PI = 3.1415926534;

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
    std::vector<occHandle(Ng_MeshVS_DataSource3D)> vecMeshAtWalls;

    //! ---------------------------------------
    //! decide which nodes should be displaced
    //! ---------------------------------------
    this->computeVecFieldCutOff(myLockBoundary);

    //! ---------------------------
    //! make a copy the outer mesh
    //! ---------------------------
    occHandle(Ng_MeshVS_DataSourceFace) theMeshToInflate = new Ng_MeshVS_DataSourceFace(myOverallSumMeshDS);

    //! --------------------------------------------------------------------
    //! retrieve the points at the boundary (if any) of the prismatic faces
    //! --------------------------------------------------------------------
    if(myPrismaticFacesSumMeshDS->myBoundarySegments.isEmpty())
    {
        //myPrismaticFacesSumMeshDS->computeFreeMeshSegments();
    }
    else
    {
        //exit(1);
    }

    //! ------------------------------------------------------------
    //! precomputation
    //! - preliminary calculation of the shrink factors
    //! - nominal total thickness of the layer at a given point
    //!   including the shrink effect due to curvature. Actually
    //!   map of the shrink factors change at each layer inflation,
    //!   while here it is considered as constant
    //! ------------------------------------------------------------
    cout<<"prismaticLayers::generateTetLayers()->____start precomputations____"<<endl;

    this->computeBeta(theMeshToInflate);
    QMap<int,double> shrinkFactors_;
    this->computeShrinkFactor(theMeshToInflate,shrinkFactors_);
    std::map<int,double> nominalTotalThicknessAtPoint;
    for(TColStd_MapIteratorOfPackedMapOfInteger it(theMeshToInflate->GetAllNodes()); it.More(); it.Next())
    {
        int globalNodeID = it.Key();
        double shrinkFactor = shrinkFactors_.value(globalNodeID);       // the shrink factor is considered as constant
        double cutOff = myLayerHCutOff.value(globalNodeID);

        double totalMarchingDistanceAtPoint = 0;
        for(int n=1; n<=myNbLayers; n++)
            totalMarchingDistanceAtPoint += myLayerThickness.at(n-1)*(1+shrinkFactor)*cutOff;
        nominalTotalThicknessAtPoint.insert(std::make_pair(globalNodeID,totalMarchingDistanceAtPoint));
    }

    //! -----------------------------------------------------------------------
    //! analyze gaps at the beginning: scan all the mesh nodes.
    //! compute the distance from closest point on the mesh => "distance"
    //! apply reduction of first layer thickness if required (flag "compress")
    //! -----------------------------------------------------------------------
    double Sigma = 0;                                                                                  // a parameter for reduction factor calculation
    for(int n=0; n<myNbLayers; n++) Sigma += pow(myExpRatio,n);
    pointToMeshDistance aDistanceMeter_;                                                                // distance meter tool
    aDistanceMeter_.init(theMeshToInflate);                                                             // init with mesh
    double firstLayerThickness = myLayerThickness.at(0);                                               // first layer thickness
    std::map<int,double> mapOfFirstLayerReductionFactor;                                               // fill a map (nodeID, reduction factor)

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
        float gap;
        this->getPointCoordinates(theMeshToInflate,globalNodeID,P);

        //! --------------------------------------------------
        //! invert the normal: it points towards the material
        //! for gap calculation
        //! --------------------------------------------------
        double nx = -normal[0];
        double ny = -normal[1];
        double nz = -normal[2];
        double dir[3] {nx,ny,nz};
        aDistanceMeter_.distance(P,dir,&gap);   // compute the gap

        bool compress = true;
        if(compress==true)
        {
            //! ---------------------------------------------
            //! compress layers in order to avoid collisions
            //! ---------------------------------------------
            const double compressionFactor = 0.25;
            if(nominalTotalThicknessAtPoint.at(globalNodeID)<gap/3.0)
            {
                //! the node can move with the nominal marching distance
                mapOfFirstLayerReductionFactor.insert(std::make_pair(globalNodeID,1.0));
            }
            else
            {
                //! the node can move with smaller marching distance
                double firstLayerThicknessNew = compressionFactor*gap/Sigma;
                double reductionFactor = firstLayerThicknessNew/firstLayerThickness;
                mapOfFirstLayerReductionFactor.insert(std::make_pair(globalNodeID,reductionFactor));
                mapOfReductionFactor.insert(std::make_pair(globalNodeID,reductionFactor));
            }
        }
        else
        {
            if(nominalTotalThicknessAtPoint.at(globalNodeID)<gap/3.0)
            {
                //! the node can move without compression
                mapOfFirstLayerReductionFactor.insert(std::make_pair(globalNodeID,1.0));
            }
            else
            {
                //! lock the point - the node cannot move anymore
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

    cout<<"prismaticLayers::generateTetLayers()->____precomputations done____"<<endl;

    //! --------------------
    //! generate the layers
    //! --------------------
    QMap<int,double> marchingDistanceMapOld;
    int validTetrahedron = 0;
    for(int n=1; n<=myNbLayers; n++)
    {
        cout<<"prismaticLayer::generateTetLayers()->____building layer #"<<n<<"____"<<endl;

        std::map<int,meshElementByCoords> mapOfVolumeElementsAtWalls_;
        std::map<int,std::vector<int>> nodeToTetrahedra;
        //int validTetrahedron = 0;

        int NbBlocked = 0;
        for(QMap<int,double>::iterator it = myLayerHCutOff.begin(); it!=myLayerHCutOff.end(); it++)
            if(it.value()!=0) NbBlocked++;
        if(NbBlocked==0)
        {
            cout<<"******************************************************************"<<endl;
            cout<<" the layer generation has been terminated in advance at step #"<<n<<endl;
            cout<<" since all the nodes has been locked "<<endl;
            cout<<"******************************************************************"<<endl;
            break;
        }

        //! ----------------------------------------------------------------
        //! we store a copy of the "previous" mesh so we can "undo"
        //! some invalid nodal displacements during self intersection check
        //! ----------------------------------------------------------------
        occHandle(Ng_MeshVS_DataSourceFace) theMeshToInflate_old = new Ng_MeshVS_DataSourceFace(theMeshToInflate);

        //! *******************************************
        //! generate the layer of tetrahedra at step n
        //! *******************************************

        //! ---------------------------------------------------
        //! the current layer thickness - nominal displacement
        //! ---------------------------------------------------
        double displacement = myLayerThickness.at(n-1);

        //! -----------------------------------------------------------
        //! manifold characteristic for each point of the current mesh
        //! visibility angle for each point of the current mesh
        //! -----------------------------------------------------------
        this->computeBeta(theMeshToInflate);

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

        //! ----------------------------------------------------------------------------
        //! guiding vectors - "best normals" - we do not use a const reference
        //! since "normals" are passed to the smoothingTools, which will change the map
        //! ----------------------------------------------------------------------------
        theMeshToInflate->computeNormalAtNodes();
        QMap<int,QList<double>> normals = theMeshToInflate->myNodeNormals;

        //! ---------------
        //! classify nodes
        //! ---------------
        std::map<int,int> mapNodeTypes;
        this->classifyNodes(theMeshToInflate,mapNodeTypes);     //! ERR 1 - theMeshToInflate is watertight: it has no 1D boundary

        //! -----------------------------------------------------------------------------------------------
        //! smooth the guiding vectors directions - use 5/10 smoothing steps
        //! flag normalize == true since the "rotation" of the vector is of interest (change in direction)
        //! -----------------------------------------------------------------------------------------------
        smoothingTools::fieldSmoother(normals,theMeshToInflate,betaAverageField,betaVisibilityField,mapNodeTypes,10,true,n);

        //! -------------------------------------------------------------------------------------
        //! create the map of the marching distances
        //! it includes the compression/expansion and the reduction of the first layer thickness
        //! -------------------------------------------------------------------------------------
        QMap<int,double> marchingDistanceMap;
        for(QMap<int,QList<double>>::const_iterator it = normals.cbegin(); it!= normals.cend(); ++it)
        {
            int globalNodeID = it.key();
            double shrinkFactor = shrinkFactors.value(globalNodeID);
            double cutOff = myLayerHCutOff.value(globalNodeID);
            double marchingDistance = displacement*(1+shrinkFactor)*cutOff*mapOfFirstLayerReductionFactor[globalNodeID];
            marchingDistanceMap.insert(globalNodeID,marchingDistance);
        }

        //! -------------------------------------------
        //! it seems able to solve several issues!
        //! note: the marching distance map is changed
        //! -------------------------------------------
        this->checkLateralDistributionMarchingDistance(theMeshToInflate,n,marchingDistanceMap,normals);

        //! ------------------------------
        //! smooth the marching distances
        //! ------------------------------
        smoothingTools::scalarFieldSmoother(marchingDistanceMap,theMeshToInflate,betaAverageField,mapNodeTypes,10,n);

        //! ---------------------------------------------------------------------------------------------------------------------
        //! analyze front - the front analyzer changes the map of the marching distance and also myLayerHCutOff
        //! 1. Its marching distance is less than a fraction, say 0.60, of the average length of edges sharing p in its manifold
        //! 2. The prismatic element created by marching p has an acceptable quality
        //! ---------------------------------------------------------------------------------------------------------------------
        this->analyzeFront(theMeshToInflate,normals,marchingDistanceMap);

        //! --------------------------------------------------------
        //! check the marching distance between two inflation steps
        //! options: lock/adjust
        //! --------------------------------------------------------
        if(n>1)
        {
            cout<<"******************************************************************"<<endl;
            cout<<" comparying marching distance values wrt previous layer"<<endl;

            int NbTooLargeIncrementsWrtOld = 0;
            for(QMap<int,double>::iterator it = marchingDistanceMap.begin(); it != marchingDistanceMap.end(); it++)
            {
                int globalNodeID = it.key();
                if(myLayerHCutOff.value(globalNodeID)==0) continue; // jump over the nodes which should not be displaced
                double curMarchingDistance = it.value();
                double previousMarchingDistance = marchingDistanceMapOld.value(globalNodeID);
                if((curMarchingDistance/previousMarchingDistance)>=myExpRatio*2.0)   // 2.0 is the maximum of (1+shrinkFactor)
                {
                    NbTooLargeIncrementsWrtOld++;
                    cout<<" node ID: "<<globalNodeID<< "has increment of marching distance too large. Ratio: "<<curMarchingDistance/previousMarchingDistance<<endl;

                    //! ------------------------
                    //! this is the option lock
                    //! ------------------------
                    myLayerHCutOff.insert(globalNodeID,0);
                    marchingDistanceMap.insert(globalNodeID,0);
                    marchingDistanceMapOld.insert(globalNodeID,0);
                }
                else
                {
                    marchingDistanceMapOld.insert(globalNodeID,curMarchingDistance);
                }
            }

            cout<<" corrected #"<<NbTooLargeIncrementsWrtOld<<" marching distances"<<endl;
            cout<<"******************************************************************"<<endl;

        }
        else // this occurs @ first step
        {
            //! ---------------------------------------------------------
            //! initialize the marching distance "old" at the first step
            //! ---------------------------------------------------------
            for(QMap<int,double>::iterator it = marchingDistanceMap.begin(); it != marchingDistanceMap.end(); it++)
            {
                int globalNodeID = it.key();
                double curMarchingDistance = it.value();
                marchingDistanceMapOld.insert(globalNodeID,curMarchingDistance);
            }
        }

        //! --------------------------------------------------
        //! final check - check collisions
        //! measure the gaps along slightly different normals
        //! --------------------------------------------------
        cout<<"******************************************************************"<<endl;
        cout<<" layer #"<<n<<" => checking intersections before advancing nodes"<<endl;

        int NbGapNegativeChecks = 0;
        pointToMeshDistance aDistanceMeter;
        aDistanceMeter.init(theMeshToInflate);
        for(TColStd_MapIteratorOfPackedMapOfInteger it(theMeshToInflate->GetAllNodes()); it.More(); it.Next())
        {
            int globalNodeID = it.Key();
            if(myLayerHCutOff.value(globalNodeID)==0) continue;
            float gap;
            const QList<double> &normal = theMeshToInflate->myNodeNormals.value(globalNodeID);
            double dir[3] { -normal[0], -normal[1], -normal[2] };
            double P[3];
            this->getPointCoordinates(theMeshToInflate,globalNodeID,P);

            //! -----------------------------------------
            //! gap measurement along several directions
            //! the returning value is the minimum
            //! -----------------------------------------
            double angleSpan = 2.5*PI/180.0;
            aDistanceMeter.distance(P,dir,angleSpan,5,&gap);

            //! --------------
            //! lock the node
            //! --------------
            if(marchingDistanceMap.value(globalNodeID)>0.3*gap)
            {
                NbGapNegativeChecks++;
                marchingDistanceMap.insert(globalNodeID,0);
                myLayerHCutOff.insert(globalNodeID,0);
            }
        }
        cout<<" intersection check: blocking #"<<NbGapNegativeChecks<<" nodes"<<endl;
        cout<<"******************************************************************"<<endl;

        //! -----------------------------------------------
        //! build the components of the displacement field
        //! -----------------------------------------------
        QMap<int,QList<double>> displacementsField;
        for(QMap<int,QList<double>>::const_iterator it = normals.cbegin(); it!= normals.cend(); ++it)
        {            
            int globalNodeID = it.key();
            if(myLayerHCutOff.contains(globalNodeID)==0) continue;
            const QList<double> &curNormal = it.value();

            //! -------------------------
            //! jump over boundary nodes
            //! -------------------------
            if(myPrismaticFacesSumMeshDS->myBoundaryPoints.contains(globalNodeID)) continue;

            //! --------------------
            //! nodal displacements
            //! --------------------
            double marchingDistance = marchingDistanceMap.value(globalNodeID);
            double cutOff = myLayerHCutOff.value(globalNodeID);
            double vx = -curNormal[0]*marchingDistance*cutOff;
            double vy = -curNormal[1]*marchingDistance*cutOff;
            double vz = -curNormal[2]*marchingDistance*cutOff;

            //! -------------------
            //! displacement field
            //! -------------------
            QList<double> localFieldValue;
            localFieldValue<<vx<<vy<<vz;
            displacementsField.insert(globalNodeID,localFieldValue);
        }

        //! ---------------------------------------------------
        //! scan the advancing nodes and build the thetrahedra
        //! jump over the locked nodes
        //! ---------------------------------------------------
        for(TColStd_MapIteratorOfPackedMapOfInteger it(theMeshToInflate->GetAllNodes()); it.More(); it.Next())
        {
            int globalNodeID = it.Key();
            if(myPrismaticFacesSumMeshDS->myBoundaryPoints.contains(globalNodeID)) continue;

            //! ---------------------------
            //! jump over the locked nodes
            //! ---------------------------
            if(myLayerHCutOff.value(globalNodeID)==0) continue;

            //! ------------------------------
            //! prepare the "displaced point"
            //! ------------------------------
            double pc[3];
            this->getPointCoordinates(theMeshToInflate,globalNodeID,pc);
            int localNodeID = theMeshToInflate->myNodesMap.FindIndex(globalNodeID);
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
                int localElementID_attached = attachedElements[i];
                int globalElementID_attached = theMeshToInflate->myElementsMap.FindKey(localElementID_attached);

                //! -------------------------------------------------------------------
                //! a boundary volume element - the volume element ID is not specified
                //! -------------------------------------------------------------------
                meshElementByCoords aVolElement;
                aVolElement.type = TET;
                aVolElement.ID = -1;

                int NbNodes, buf[3];
                TColStd_Array1OfInteger nodeIDs(*buf,1,3);
                theMeshToInflate->GetNodesByElement(globalElementID_attached,nodeIDs,NbNodes);

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
                /*
                const mesh::meshPoint &P0 = aVolElement.pointList[0];
                const mesh::meshPoint &P1 = aVolElement.pointList[1];
                const mesh::meshPoint &P2 = aVolElement.pointList[2];
                double P10_x = P1.x-P0.x;
                double P10_y = P1.y-P0.y;
                double P10_z = P1.z-P0.z;
                double P20_x = P2.x-P0.x;
                double P20_y = P2.y-P0.y;
                double P20_z = P2.z-P0.z;

                //! ----------------------
                //! i       j       k
                //! P10_x   P10_y   P10_z
                //! P20_x   P20_y   P20_z
                //! ----------------------
                double nx = P10_y*P20_z-P10_z*P20_y;
                double ny = P10_z*P20_x-P10_x*P20_z;
                double nz = P10_x*P20_y-P10_y*P20_x;
                double ln = sqrt(nx*nx+ny*ny+nz*nz);
                nx /= ln; ny /= ln; nz /= ln;

                //! ------
                //! (O-P0)
                //! ------
                double nnx = mp_shifted.x-P0.x;
                double nny = mp_shifted.y-P0.y;
                double nnz = mp_shifted.z-P0.z;

                double dot = nx*nnx+ny*nny+nz*nnz;
                */
                TetQualityClass aTetQuality;
                aTetQuality.setPoints(aVolElement.getPoints());
                double volume = aTetQuality.Volume();
                if(volume<=0.0)
                {
                    //! ---------------------------------------------------------
                    //! remove all the elements generated from this triangle fan
                    //! this avoids abrupt changes in the surface mesh
                    //! ---------------------------------------------------------
                    cerr<<"____found negative volume tet when moving node: "<<globalNodeID<<"\tVolume = "<<volume<<"____"<<endl;
                    //for(int j = 0; j<i; j++) vecElementsToAdd.pop_back();
                    vecElementsToAdd.clear();
                    NbInvalidTets++;
                    canMove = false;
                    break;
                }
                vecElementsToAdd.push_back(aVolElement);
            }
            if(canMove)
            {
                //! -----------------------------------------------------------
                //! move the point forward - update locally the inflation mesh
                //! -----------------------------------------------------------
                std::vector<double> P {x_displ,y_displ,z_displ};
                theMeshToInflate->changeNodeCoords(localNodeID,P);

                //! -------------------
                //! standard treatment
                //! -------------------
                for(int i=0; i<vecElementsToAdd.size(); i++)
                {
                    validTetrahedron++;
                    std::map<int,std::vector<int>>::iterator it = nodeToTetrahedra.find(globalNodeID);
                    if(it == nodeToTetrahedra.end())
                    {
                        //cout<<"____pile up: ("<<globalNodeID<<", "<<validTetrahedron<<")____"<<endl;
                        std::vector<int> v {validTetrahedron};
                        nodeToTetrahedra.insert(std::make_pair(globalNodeID,v));
                    }
                    else
                    {
                        //cout<<"____pile up: ("<<globalNodeID<<", "<<validTetrahedron<<")____"<<endl;
                        it->second.push_back(validTetrahedron);
                    }

                    //! ---------------------------------------------------
                    //! pile up the volume elements building the wall mesh
                    //! ---------------------------------------------------
                    mapOfVolumeElementsAtWalls_.insert(std::make_pair(validTetrahedron,vecElementsToAdd[i]));
                }
            }
            else
            {
                //! do not move the surface mesh point
                //! lock the node for all the subsequent steps
                myLayerHCutOff.insert(globalNodeID,0);
            }
        }

        cout<<"******************************************************************"<<endl;
        cout<<" Number of zero or negative volume tetrahedra: "<<NbInvalidTets<<endl;
        cout<<" Number of shift blocked by lenght criterion: "<<blockedByLenght<<endl;
        cout<<" Number of shift blocked by angle criterion: "<<blockedByAngle<<endl;
        cout<<" Total number of tetrahedra (some could have been removed): "<<validTetrahedron<<endl;
        cout<<"******************************************************************"<<endl;

        //! -------------------------
        //! check self intersections
        //! -------------------------
        if(myCheckSelfIntersections==true)
        {
            std::vector<int> badNodes;
            bool selfIntersection = iglTools::iglCheckSelfIntersection(theMeshToInflate, badNodes);
            if(selfIntersection==true)
            {
                //! ------------------------------------------------------
                //! remove the tetrahedra generated when moving this node
                //! ------------------------------------------------------
                cout<<"@-------------------------------------------------------"<<endl;
                cout<<"@    begin removing elements"<<endl;
                for(size_t k = 0; k<badNodes.size(); k++)
                {
                    int localNodeID_bad = badNodes[k]+1;
                    int globalNodeID_bad = theMeshToInflate->myNodesMap.FindKey(localNodeID_bad);

                    //! -----------------------------------------------
                    //! do not use "nodeToTetrahedra.at(globalNodeID)"
                    //! -----------------------------------------------
                    std::map<int,std::vector<int>>::iterator it = nodeToTetrahedra.find(globalNodeID_bad);
                    if(it==nodeToTetrahedra.end()) continue;
                    const std::vector<int> &elementsToRemove = it->second;
                    if(elementsToRemove.empty()) continue;

                    for(std::vector<int>::const_iterator it__ = elementsToRemove.cbegin(); it__ != elementsToRemove.cend(); it__++)
                    {
                        int elementToRemove = *it__;
                        mapOfVolumeElementsAtWalls_.erase(elementToRemove);
                        cout<<" "<<elementToRemove;
                    }
                }
                //cout<<endl<<"@    end of element removal. Removed: "<<badNodes.size()<<" elements"<<endl;
                cout<<"@    end of element removal. Removed: "<<badNodes.size()<<" elements"<<endl;

                //! ------------------------------------------------------------------------------------
                //! now generate the volume mesh of the layer, using the last modified list of elements
                //! ------------------------------------------------------------------------------------
                cout<<"@    generating the volume mesh of the layer"<<endl;
                std::vector<meshElementByCoords> elementsOfTheLayer;
                for(std::map<int,meshElementByCoords>::iterator it = mapOfVolumeElementsAtWalls_.begin(); it != mapOfVolumeElementsAtWalls_.end(); it++)
                {
                    elementsOfTheLayer.push_back(it->second);
                }
                occHandle(Ng_MeshVS_DataSource3D) volumeMeshOfLayer = new Ng_MeshVS_DataSource3D(elementsOfTheLayer,true,true);
                cout<<"@    mesh of the layer generated"<<endl;

                //! ------------------------------------------------------------------------
                //! generate the hash table of the mesh points of the current boundary mesh
                //! note: if two adjacent elements are removed also the common nodes
                //! disappear from the volume mesh
                //! ------------------------------------------------------------------------
                std::map<size_t,int> pointHashTable;    // point - localNodeID
                for(TColStd_MapIteratorOfPackedMapOfInteger itn = volumeMeshOfLayer->GetAllNodes(); itn.More(); itn.Next())
                {
                    int globalNodeID_ = itn.Key();
                    double P[3];
                    this->getPointCoordinates(volumeMeshOfLayer,globalNodeID_,P);
                    size_t seed = 0;
                    for(int n=0; n<3; n++) hash_c<double>(seed,P[n]);
                    int localNodeID_ = volumeMeshOfLayer->myNodesMap.FindIndex(globalNodeID_);
                    pointHashTable.insert(std::make_pair(seed,localNodeID_));
                }

                //! --------------------------------------------------------------------------------
                //! "reinforce" the point search (a mesh point could be found also using tolerance)
                //! --------------------------------------------------------------------------------
                std::vector<mesh::tolerantPoint> tolerantPoints = volumeMeshOfLayer->buildTolerantPoints(1e-6);

                for(size_t k = 0; k<badNodes.size(); k++)
                {
                    int localNodeID_bad = badNodes[k]+1;
                    int globalNodeID_bad = theMeshToInflate->myNodesMap.FindKey(localNodeID_bad);

                    //! ---------------------------------------------
                    //! the hash code of the current node to correct
                    //! ---------------------------------------------
                    double Pbad[3];
                    this->getPointCoordinates(theMeshToInflate, globalNodeID_bad,Pbad);
                    cout<<"____the node ID: "<<globalNodeID_bad<<" ("<<Pbad[0]<<", "<<Pbad[1]<<", "<<Pbad[2]<<") of the surface mesh is bad____"<<endl;
                    size_t nodeHashCode = 0;
                    for(int n=0; n<3; n++) hash_c<double>(nodeHashCode,Pbad[n]);

                    //! ----------------------------------------------------------------
                    //! get the local node ID of the current wall mesh using hash table
                    //! ----------------------------------------------------------------
                    std::map<size_t,int>::iterator itmap = pointHashTable.find(nodeHashCode);
                    if(itmap==pointHashTable.end())
                    {
                        mesh::tolerantPoint aTolerantPoint(Pbad[0],Pbad[1],Pbad[2]);
                        std::vector<mesh::tolerantPoint>::iterator ittp = std::find(tolerantPoints.begin(),tolerantPoints.end(),aTolerantPoint);
                        if(ittp==tolerantPoints.end())
                        {
                            cout<<"____NOR with tolerant point: this point has desappeared in removing elements____"<<endl;
                            //! this happens when two removed elements share a point
                        }
                        else
                        {
                            cout<<"____OK found point using tolerance____"<<endl;
                            int localNodeID_bad_intoWallMesh = itmap->second;

                            double Pold[3];
                            this->getPointCoordinates(theMeshToInflate_old,globalNodeID_bad,Pold);
                            std::vector<double> newCoords { Pold[0], Pold[1], Pold[2] };
                            volumeMeshOfLayer->changeNodeCoords(localNodeID_bad_intoWallMesh,newCoords);
                        }
                    }
                    else
                    {
                        //cout<<"____node found using hash code____"<<endl;
                        int localNodeID_bad_intoWallMesh = itmap->second;

                        int globalNodeID_bad_intoWallMesh = volumeMeshOfLayer->myNodesMap.FindKey(localNodeID_bad_intoWallMesh);
                        double P[3];
                        this->getPointCoordinates(volumeMeshOfLayer,globalNodeID_bad_intoWallMesh,P);
                        cout<<"____the node ID: "<<globalNodeID_bad_intoWallMesh<<" ("<<P[0]<<", "<<P[1]<<", "<<P[2]<<") of the volume mesh is bad____"<<endl;

                        double Pold[3];
                        this->getPointCoordinates(theMeshToInflate_old,globalNodeID_bad,Pold);
                        std::vector<double> newCoords { Pold[0], Pold[1], Pold[2] };
                        volumeMeshOfLayer->changeNodeCoords(localNodeID_bad_intoWallMesh,newCoords);
                    }
                }
                //! ----------------------
                //! pile up the wall mesh
                //! ----------------------
                vecMeshAtWalls.push_back(volumeMeshOfLayer);

                //! --------------------------------------------
                //! "rewind" the surface mesh in the bad points
                //! --------------------------------------------
                for(int k=0; k<badNodes.size(); k++)
                {
                    int localNodeID_bad = badNodes[k]+1;
                    const std::vector<double> &old_coords = theMeshToInflate_old->getNodeCoordinates(localNodeID_bad);
                    theMeshToInflate->changeNodeCoords(localNodeID_bad,old_coords);

                    //! --------------
                    //! lock the node
                    //! --------------
                    int globalNodeID_bad = theMeshToInflate->myNodesMap.FindKey(localNodeID_bad);
                    myLayerHCutOff.insert(globalNodeID_bad,0);
                }

                cout<<"@    storing a temporary mesh of "<<volumeMeshOfLayer->GetAllNodes().Extent()<<" nodes, "<<volumeMeshOfLayer->GetAllElements().Extent()<<" elements"<<endl;
                cout<<"@-------------------------------------------------------"<<endl;
            }
            else
            {
                //! -----------------------------------------------------------
                //! no self intersection has been found temporary volume mesh
                //! -----------------------------------------------------------
                cout<<"@-------------------------------------------------------"<<endl;
                cout<<"@    generating the volume mesh of the layer"<<endl;
                std::vector<meshElementByCoords> elementsOfTheLayer;
                for(std::map<int,meshElementByCoords>::iterator it = mapOfVolumeElementsAtWalls_.begin(); it != mapOfVolumeElementsAtWalls_.end(); it++)
                {
                    elementsOfTheLayer.push_back(it->second);
                }
                //! ----------------------
                //! pile up the wall mesh
                //! ----------------------
                occHandle(Ng_MeshVS_DataSource3D) volumeMeshOfLayer = new Ng_MeshVS_DataSource3D(elementsOfTheLayer,true,true);
                vecMeshAtWalls.push_back(volumeMeshOfLayer);
                cout<<"@    mesh of the layer generated"<<endl;
                cout<<"@-------------------------------------------------------"<<endl;
            }
        }
        else
        {
            //! ------------------------------------------------------
            //! no self-intersection check required
            //! generate the volume mesh of the for the current layer
            //! ------------------------------------------------------
            cout<<"@-------------------------------------------------------"<<endl;
            cout<<"@    generating the volume mesh of the layer"<<endl;
            std::vector<meshElementByCoords> elementsOfTheLayer;
            for(std::map<int,meshElementByCoords>::iterator it = mapOfVolumeElementsAtWalls_.begin(); it != mapOfVolumeElementsAtWalls_.end(); it++)
            {
                elementsOfTheLayer.push_back(it->second);
            }
            occHandle(Ng_MeshVS_DataSource3D) volumeMeshOfLayer = new Ng_MeshVS_DataSource3D(elementsOfTheLayer,true,true);
            vecMeshAtWalls.push_back(volumeMeshOfLayer);
            cout<<"@    mesh of the layer generated"<<endl;
            cout<<"@-------------------------------------------------------"<<endl;
        }

        //! --------------------------------------------
        //! check "Cliff points" - incomplete manifolds
        //! --------------------------------------------
        this->checkIncompleteManifold(theMeshToInflate);

        //! --------------
        //! send progress
        //! --------------
        if(myProgressIndicator!=Q_NULLPTR)
        {
            progressEvent = new QProgressEvent(QProgressEvent_None,0,9999,0,QString("Generating boundary mesh layer # %1").arg(n),
                                               QProgressEvent_Update,-1,-1,n,"Building prismatic mesh");
            QApplication::postEvent(myProgressIndicator,progressEvent);
            QApplication::processEvents();
        }
    }

    //! --------------------------------------------------------------------------
    //! the last modified surface mesh - should be returned for the volume mesher
    //! --------------------------------------------------------------------------
    lastInflatedMesh = theMeshToInflate;

    //! -----------
    //! diagnostic
    //! -----------
    //MeshTools::saveSTL(theMeshToInflate,"D:/finalMesh.stl");

    //! ------------------------------------------------------
    //! construction of the overall tetrahedral boundary mesh
    //! (sum of the volume mesh of each layer)
    //! ------------------------------------------------------
    std::vector<meshElementByCoords> finalVecElements;
    for(std::vector<occHandle(Ng_MeshVS_DataSource3D)>::iterator itm = vecMeshAtWalls.begin(); itm!=vecMeshAtWalls.end(); itm++)
    {
        const occHandle(Ng_MeshVS_DataSource3D) &curvolumeMeshOfLayer = *itm;
        for(TColStd_MapIteratorOfPackedMapOfInteger ite(curvolumeMeshOfLayer->GetAllElements()); ite.More(); ite.Next())
        {
            meshElementByCoords aMeshElem;
            int globalElementID = ite.Key();
            int NbNodes, buf[20];

            aMeshElem.ID = -1;
            aMeshElem.type = TET;

            TColStd_Array1OfInteger nodeIDs(*buf,1,20);
            curvolumeMeshOfLayer->GetNodesByElement(globalElementID,nodeIDs,NbNodes);
            for(int n=1; n<=NbNodes; n++)
            {
                int globalNodeID = nodeIDs(n);
                double P[3];
                this->getPointCoordinates(curvolumeMeshOfLayer,globalNodeID,P);
                aMeshElem.pointList<<mesh::meshPoint(P[0],P[1],P[2],-1);
            }
            finalVecElements.push_back(aMeshElem);
        }
    }

    meshAtWalls = new Ng_MeshVS_DataSource3D(finalVecElements,true,true);
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

    //! -------------------------------------------------------------------------
    //! Tetrahedral mesh at walls. Important note: when using generateTetLayer
    //! the "lastInflatedMesh" is defirnatuib if the initial, outer surface mesh
    //! which should be used for generating the volume mesh
    //! -------------------------------------------------------------------------
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

//! ---------------------------------------------------
//! function: checkLateralDistributionMarchingDistance
//! details:  does not perform intersection check
//! ---------------------------------------------------
void prismaticLayer::checkLateralDistributionMarchingDistance(const occHandle(Ng_MeshVS_DataSourceFace) &aMeshDS,
                                                              int inflationStep,
                                                              QMap<int,double> &marchingDistanceMap,
                                                              const QMap<int,QList<double>> &guidingDirections)
{
    cout<<"prismaticLayer::checkLateralDistributionMarchingDistance()->____function called____"<<endl;

    const double PI = 3.1415926538;

    //! ---------------------
    //! a gap measuring tool
    //! ---------------------
    pointToMeshDistance aDistanceMeter;
    aDistanceMeter.init(aMeshDS);

    std::map<int,double> marchingDistanceMapUpdated;
    const int NMaxIterations = 100;          //! use this for the option "update"
    //const int NMaxIterations = 1;         //! use this for the option "lock"
    int NbAdjusted_old = 1e10;
    for(int n=0; n<NMaxIterations; n++)
    {
        int NbAdjusted = 0;
        for(QMap<int,double>::iterator it = marchingDistanceMap.begin(); it!=marchingDistanceMap.end(); it++)
        {
            int globalNodeID = it.key();
            if(myLayerHCutOff.value(globalNodeID)==0) continue;
            double distanceToCheck = it.value();

            //! -----------------------------------------------------------------
            //! among the surrounding nodes, check the ones which do not satisfy
            //! the inner inequality. The fixed nodes around must be considered
            //! -----------------------------------------------------------------
            bool foundNodesAround = false;
            const std::vector<int> &surroundingNodes = mySurroudingNodesMap.at(globalNodeID);
            for(std::vector<int>::const_iterator it_ = surroundingNodes.cbegin(); it_!=surroundingNodes.cend(); it_++)
            {
                int globalNodeID_surrounding = *it_;
                double curMarchingDistance_surr = marchingDistanceMap.value(globalNodeID_surrounding);
                if(distanceToCheck<0.5*curMarchingDistance_surr || distanceToCheck>2.0*curMarchingDistance_surr)
                {
                    NbAdjusted++;
                    foundNodesAround = true;
                    break;
                }
            }

            //! --------------------------
            //! this is the option "lock"
            //! --------------------------
            //if(foundNodesAround==true)
            //{
            //    myLayerHCutOff.insert(globalNodeID,0);
            //    marchingDistanceMapUpdated.insert(std::make_pair(globalNodeID,0));
            //}

            //! ----------------------------
            //! this is the option "adjust"
            //! ----------------------------
            if(foundNodesAround == true)
            {
                double betaAve = betaAverageField.at(globalNodeID);
                double smoothedValue;
                smoothingTools::smoothScalarAtPoint(aMeshDS,marchingDistanceMap,globalNodeID,betaAve,inflationStep,smoothedValue);
                std::map<int,double>::iterator it_ = marchingDistanceMapUpdated.find(globalNodeID);
                if(it_!=marchingDistanceMapUpdated.end()) marchingDistanceMapUpdated.insert(std::make_pair(globalNodeID,smoothedValue));
                else it_->second = smoothedValue;
            }
        }
        if(NbAdjusted>=NbAdjusted_old) break;
        NbAdjusted_old = NbAdjusted;
    }

    //! ------------------------------
    //! update the marching distances
    //! ------------------------------
    int NbLocked = 0;
    for(std::map<int,double>::iterator it = marchingDistanceMapUpdated.begin(); it != marchingDistanceMapUpdated.end(); it++)
    {
        int globalNodeID = it->first;
        double updatedMarchingDistance = it->second;
        marchingDistanceMap.insert(globalNodeID,updatedMarchingDistance);
    }
    cout<<"prismaticLayer::checkLateralDistributionMarchingDistance()->____lateral criterion: evaluated "<<NbAdjusted_old<<" nodes____"<<endl;
    cout<<"prismaticLayer::checkLateralDistributionMarchingDistance()->____lateral criterion: locked "<<NbLocked<<" nodes____"<<endl;
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

        //! ---------------------------------------------------------------
        //! a point on the 1D boundary (wall/notwall is always incomplete)
        //! ---------------------------------------------------------------
        if(theMeshToInflate->myBoundaryPoints.contains(globalNodeID)) continue;

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

//! -------------------------------------------------------------------------------------------------------------------------------
//! function: analyzeFront
//! details:  front analysis: check the marching distances
//!           1. Its marching distance is less than a fraction, say 0.50, of the average length of edges sharing p in its manifold
//!           2. The prismatic element created by marching p has an acceptable quality
//! -------------------------------------------------------------------------------------------------------------------------------
void prismaticLayer::analyzeFront(const occHandle(Ng_MeshVS_DataSourceFace) &theMeshToInflate,
                                  const QMap<int,QList<double>> &normals,
                                  QMap<int,double> &marchingDistanceMap)
{
    cout<<"prismaticLayer::analyzeFront()->____analyzing front____"<<endl;

    if(theMeshToInflate.IsNull()) return;

    //! ---------------------------
    //! constant values and a tool
    //! ---------------------------
    const double PI = 3.1415926538;
    pointToMeshDistance aDistanceMeter;
    aDistanceMeter.init(theMeshToInflate);

    int blockedByLenght = 0;
    int blockedByAngle = 0;
    for(TColStd_MapIteratorOfPackedMapOfInteger it(theMeshToInflate->GetAllNodes()); it.More(); it.Next())
    {
        int globalNodeID = it.Key();
        std::vector<int> setOfSurroundingNodesIDs_global = mySurroudingNodesMap.at(globalNodeID);
        size_t NbSurroundingElements = setOfSurroundingNodesIDs_global.size();

        //! ---------------------------------------------------------------------
        //! perform first check - stop propagating node if the marching distance
        //! at point is greater than 0.3 the average length of the segments
        //! connectect to the point
        //! -------------------------------------------------------------------
        double P[3], averageDistance = 0;
        this->getPointCoordinates(theMeshToInflate,globalNodeID,P);
        for(std::vector<int>::iterator it__ = setOfSurroundingNodesIDs_global.begin(); it__!=setOfSurroundingNodesIDs_global.end(); it__++)
        {
            int curGlobalNodeID = *it__;
            double Pcurr[3];
            this->getPointCoordinates(theMeshToInflate,curGlobalNodeID,Pcurr);
            double distance = sqrt(pow(P[0]-Pcurr[0],2)+pow(P[1]-Pcurr[1],2)+pow(P[2]-Pcurr[2],2));
            averageDistance += distance;
        }
        averageDistance /= NbSurroundingElements;

        if(marchingDistanceMap.value(globalNodeID)>0.50*averageDistance)
        {
            cout<<"prismaticLayer::analyzeFront()->____node ID "<<globalNodeID<<": marching distance too large____"<<endl;

            //! --------------
            //! lock the node
            //! --------------
            myLayerHCutOff.insert(globalNodeID,0);
            blockedByLenght++;

            /*
            //! ----------------------------------------
            //! adjust - this option should be improved
            //! ----------------------------------------
            double P[3];
            this->getPointCoordinates(theMeshToInflate,globalNodeID,P);
            float gap = 0.0;
            const QList<double> &curGuidingVector = normals.value(globalNodeID);
            double dir[3] { -curGuidingVector[0], -curGuidingVector[1], -curGuidingVector[2] };
            double span = 2.5*PI/180.0;
            aDistanceMeter.distance(P,dir,span,8,&gap);

            //! ------------------------------------------------------------------
            //! start correcting the marching distance from 0.50*averageDistance
            //! 1. check for proximity in order to avoid collisions
            //! 2. if the candidate value is too large bisect it
            //! 3. if the candidate value is too small break
            //! 4. lock the node if no suitable value has been found
            //! ------------------------------------------------------------------
            bool foundValidMarchingDistance = false;
            double candidateMarchingDistance = 0.50 * averageDistance;
            for(;;)
            {
                if(candidateMarchingDistance<=0.10*0.50*averageDistance)
                {
                    //! the candidate marching distance is too small;
                    break;
                }
                if(candidateMarchingDistance < 0.3*gap)
                {
                    //! can use this distance avoiding collisions
                    marchingDistanceMap.insert(globalNodeID,candidateMarchingDistance);
                    foundValidMarchingDistance = true;
                    break;
                }
                //! "bisection" or something smoother like this
                candidateMarchingDistance *= 0.9;
            }
            if(foundValidMarchingDistance == false)
            {
                //! --------------
                //! lock the node
                //! --------------
                cerr<<"prismaticLayer::analyzeFront()->____node ID "<<globalNodeID<<": cannot update marching distance. Locking node____"<<endl;
                myLayerHCutOff.insert(globalNodeID,0);
                marchingDistanceMap.insert(globalNodeID,0);
                blockedByLenght++;
            }
            */
        }

        //! -----------------------------------------------------------------------
        //! second check - compute angles between the guiding vector
        //! and all the edges generated by P* (P shifted along the guiding vector)
        //! -----------------------------------------------------------------------
        double P_star[3];
        this->getPointCoordinates(theMeshToInflate,globalNodeID,P_star);
        const QList<double> &normal = normals.value(globalNodeID);
        double curMarchingDistance = marchingDistanceMap.value(globalNodeID);
        P_star[0] += -normal[0]*curMarchingDistance;
        P_star[1] += -normal[1]*curMarchingDistance;
        P_star[2] += -normal[2]*curMarchingDistance;

        double maxAngle = -1e10;
        const std::vector<int> &vecSurroundingNodes = mySurroudingNodesMap.at(globalNodeID);
        for(std::vector<int>::const_iterator it = vecSurroundingNodes.cbegin(); it!=vecSurroundingNodes.cend(); it++)
        {
            int globalNodeID_surr = *it;
            double V[3];
            this->getPointCoordinates(theMeshToInflate,globalNodeID_surr,V);
            double lx = V[0] - P_star[0];
            double ly = V[1] - P_star[1];
            double lz = V[2] - P_star[2];

            double dot = (lx*normal[0]+ly*normal[1]+lz*normal[2]);
            double dotNorm = sqrt(lx*lx+ly*ly+lz*lz)*sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
            dot /= dotNorm;
            if(dot > 1) dot = 1;
            if(dot < -1) dot = -1;
            double angle = std::acos(dot);
            if(angle>maxAngle) maxAngle = angle;
        }
        if(maxAngle<=70*PI/180.0)    // this tries to avoid problems at corners, but this is not the solution
        {
            myLayerHCutOff.insert(globalNodeID,0);
            blockedByAngle++;
        }
    }
}

//! --------------------------------------
//! functions: checkTetFrontIntersections
//! details:
//! --------------------------------------
void prismaticLayer::checkTetFrontIntersections(const occHandle(Ng_MeshVS_DataSourceFace) &aMeshDS,
                                                const QMap<int,QList<double>> &guidingVectors,
                                                QMap<int,QList<double>> &displacementMap)
{
    ;
}
