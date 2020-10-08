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

//! ----
//! OCC
//! ----
#include <TopoDS_Face.hxx>
#include <TopoDS.hxx>

//! ------------------------------
//! function: plotAverageNormals
//! details:  diagnostic function
//! ------------------------------
void plotAverageNormals(const occHandle(Ng_MeshVS_DataSourceFace) &aFaceMeshDS)
{
    for(QMap<int,QList<double>>::const_iterator it=aFaceMeshDS->myNodeNormals.cbegin(); it!=aFaceMeshDS->myNodeNormals.cend(); ++it)
    {
        int nodeID = it.key();
        double nx = it.value().at(0);
        double ny = it.value().at(1);
        double nz = it.value().at(2);
        cout<<"____nodeID: "<<nodeID<<"("<<nx<<", "<<ny<<", "<<nz<<")____"<<endl;
    }
}

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

    //! ---------------
    //! new parameters
    //! ---------------
    myCurvatureSensitivityForShrink=100;
    myNbGuidingVectorSmoothingSteps=50;
    myNbLayerThicknessSmoothingSteps=50;
    myCurvatureSensitivityForGuidingVectorsSmoothing =5;
    myCurvatureSensitivityForThicknessSmoothing=50;
    //! ----------------------
    //! end of new parameters
    //! ----------------------

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

    myCurvatureSensitivityForShrink = parameters.curvatureSensitivityForShrink;
    myCurvatureSensitivityForGuidingVectorsSmoothing = parameters.curvatureSensitivityForGuidingVectorsSmoothing;
    myNbGuidingVectorSmoothingSteps = parameters.NbGuidingVectorSmoothingSteps;
    myCurvatureSensitivityForThicknessSmoothing = parameters.curvatureSensitivityForThicknessSmoothing;
    myNbLayerThicknessSmoothingSteps = parameters.NbLayerThicknessSmoothingSteps;
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
    if(myProgressIndicator!=Q_NULLPTR)
    {
        this->setStopButtonEnabled(false);
    }

    //! --------------------------------------
    //! decide which node should be displaced
    //! --------------------------------------
    this->computeVecFieldCutOff(lockBoundary);

    //! ---------------------------
    //! make a copy the outer mesh
    //! ---------------------------
    occHandle(Ng_MeshVS_DataSourceFace) theMeshToInflate_old = new Ng_MeshVS_DataSourceFace(myOverallSumMeshDS);
    theMeshToInflate_old->computeNormalAtElements();
    theMeshToInflate_old->computeNormalAtNodes();
    theMeshToInflate_old->computeAnglesAtNode();

    //! -----------------------
    //! enable the Stop button
    //! -----------------------
    if(myProgressIndicator!=Q_NULLPTR) myProgressIndicator->enableStop();

    //! -----------------------------------------------
    //! check self intersection at the very first step
    //! -----------------------------------------------
    //cout<<"____checking self intersections for the surface mesh____"<<endl;
    //this->isSelfIntersecting(theMeshToInflate_old);

    //! -----------------------------------------------------------
    //! insert the outer mesh into the list of the inflated meshes
    //! -----------------------------------------------------------
    inflatedMeshes<<theMeshToInflate_old;
    for(int n=1; n<=myNbLayers; n++)
    {
        if(Global::status().code = 0)
        {
            cout<<"\\ ------------------------------------------\\"<<endl;
            cout<<"\\ inflation process interrupted by the user \\"<<endl;
            cout<<"\\ ------------------------------------------\\"<<endl;
            return false;
        }

        displacement = myLayerThickness.at(n-1);

        //! ----------------
        //! console message
        //! ----------------
        cout<<"@____generating layer: "<<(n-1)<<" - current thickness: "<<displacement<<"____@"<<endl;

        occHandle(Ng_MeshVS_DataSourceFace) theMeshToInflate_new = new Ng_MeshVS_DataSourceFace(theMeshToInflate_old);

        //! -------------------------
        //! inflation needs
        //! 1) normal at nodes
        //! 2) curvature information
        //! -------------------------
        theMeshToInflate_new->computeNormalAtNodes();
        theMeshToInflate_new->computeAnglesAtNode();

        //! ---------------------------------------------------
        //! curvature before displacing: "0" mean "1" gaussian
        //! ---------------------------------------------------
        int mode = 1; //! Gaussian curvature
        theMeshToInflate_new->computeDiscreteCurvature(mode);
        QMap<int,double> shrinkFactors;
        this->computeShrinkFactor(theMeshToInflate_new,shrinkFactors);

        //! --------------------------------------------------------------------
        //! create the displacements vectorial field
        //! Important note:
        //! - "myNodeNormals" is internally built using the global node IDs
        //! - "displacementsField" is also built using the global node IDs
        //! - Ng_MeshVS_DataSourceFace::displaceMySelf uses the global node IDs
        //! --------------------------------------------------------------------
        QMap<int,QList<double>> normals = theMeshToInflate_new->myNodeNormals;

        //! ---------------
        //! a sanity check
        //! ---------------
        if(normals.isEmpty()) return false;

        //! ----------------------------------
        //! smooth guiding vectors directions
        //! ----------------------------------
        int NbSmoothingSteps = myNbGuidingVectorSmoothingSteps;
        double cs = myCurvatureSensitivityForGuidingVectorsSmoothing;
        this->fieldSmoother(normals,theMeshToInflate_new,cs,NbSmoothingSteps);

        QMap<int,gp_Vec> displacementsField;
        for(QMap<int,QList<double>>::const_iterator it = normals.cbegin(); it!= normals.cend(); ++it)
        {
            int globalNodeID = it.key();
            const QList<double> &curNormal = it.value();

            double cutoff = myLayerHCutOff.value(globalNodeID);
            double shrinkFactor = shrinkFactors.value(globalNodeID);
            double C = cutoff*shrinkFactor;
            double vx = -curNormal[0]*C;
            double vy = -curNormal[1]*C;
            double vz = -curNormal[2]*C;

            gp_Vec normVec(vx,vy,vz);
            displacementsField.insert(globalNodeID,normVec);
        }

        //! -------------------------------
        //! smooth the displacements field
        //! -------------------------------
        this->smoothDisplacementField(displacementsField,normals,theMeshToInflate_new);

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

        //if(myPrismaticFaces.contains(faceNr)) faceList1<<aFaceMeshDS;
        //if(!myPrismaticFaces.contains(faceNr)) faceList2<<aFaceMeshDS;
    }

    cout<<"____faceList0: "<<faceList0.length()<<"____"<<endl;
    cout<<"____faceList1: "<<faceList1.length()<<"____"<<endl;
    cout<<"____faceList2: "<<faceList2.length()<<"____"<<endl;

    myOverallSumMeshDS = new Ng_MeshVS_DataSourceFace(faceList0);
    myPrismaticFacesSumMeshDS = new Ng_MeshVS_DataSourceFace(faceList1);

    myNonPrismaticFacesSumMeshDS = new Ng_MeshVS_DataSourceFace(faceList2);
    myOverallSumMeshDS->computeNormalAtElements();
    myOverallSumMeshDS->computeNormalAtNodes();
    myPrismaticFacesSumMeshDS->computeNormalAtElements();
    myPrismaticFacesSumMeshDS->computeNormalAtNodes();
    myNonPrismaticFacesSumMeshDS->computeNormalAtElements();
    myNonPrismaticFacesSumMeshDS->computeNormalAtNodes();

    return true;
}

//! ------------------------------
//! function: computeShrinkFactor
//! details:
//! ------------------------------
void prismaticLayer::computeShrinkFactor(const occHandle(Ng_MeshVS_DataSourceFace) &aMeshDS,
                                         QMap<int,double> &shrinkFactors)
{
    //! ----------------
    //! diagnostic file
    //! ----------------
    //FILE *f = fopen("D:/shrink field.txt","w");

    int mode = 1; // Gaussian curvature
    if(aMeshDS->myCurvatureGradient.isEmpty()) aMeshDS->computeDiscreteCurvature(mode);
    for(TColStd_MapIteratorOfPackedMapOfInteger it(aMeshDS->GetAllNodes()); it.More(); it.Next())
    {
        int globalNodeID = it.Key();
        double curvature = aMeshDS->myCurvature.value(globalNodeID);

        //! ---------------------------
        //! use the gaussian curvature
        //! ---------------------------
        const double maxShrink = 1.0;
        double k = myCurvatureSensitivityForShrink;
        double shrink = maxShrink*(1./(1.+exp(-k*sqrt(fabs(curvature)))));
        if(curvature<0) shrink = maxShrink - shrink;
        shrinkFactors.insert(globalNodeID,shrink);

        //fprintf(f,"%d\t%lf\t%lf\n",globalNodeID,curvature,shrink);
    }
    //fclose(f);
}

//! ---------------------------------
//! function: buildPrismaticElements
//! details:
//! ---------------------------------
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
            //! ---------------
            //! end diagnostic
            //! ---------------
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
    this->mergeFaceMeshes();

    occHandle(Ng_MeshVS_DataSourceFace) prismaticMeshDS = new Ng_MeshVS_DataSourceFace(myPrismaticFacesSumMeshDS);
    occHandle(Ng_MeshVS_DataSourceFace) surfaceMeshDS = new Ng_MeshVS_DataSourceFace(myOverallSumMeshDS);
    occHandle(Ng_MeshVS_DataSourceFace) nonPrismaticMeshDS = new Ng_MeshVS_DataSourceFace(myNonPrismaticFacesSumMeshDS);

    surfaceMeshDS->computeNodeToElementsConnectivity();
    surfaceMeshDS->computeNormalAtNodes();

    prismaticMeshDS->computeNormalAtNodes();
    prismaticMeshDS->computeFreeMeshSegments();

    //! -------------------------------------
    //! make a copy the initial (outer) mesh
    //! -------------------------------------
    occHandle(Ng_MeshVS_DataSourceFace) theMeshToInflate_old = new Ng_MeshVS_DataSourceFace(prismaticMeshDS);
    theMeshToInflate_old->computeNormalAtNodes();
    theMeshToInflate_old->computeAnglesAtNode();

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
        cout<<"@____generating layer: "<<n<<" with current thickness: "<<displacement<<"____@"<<endl;

        occHandle(Ng_MeshVS_DataSourceFace) theMeshToInflate_new = new Ng_MeshVS_DataSourceFace(theMeshToInflate_old);

        //! -------------------------------------------
        //! inflation needs: normal at nodes, gradient
        //! of the normal at nodes, shrink factor
        //! -------------------------------------------
        theMeshToInflate_new->computeNormalAtNodes();
        theMeshToInflate_new->computeAnglesAtNode();

        //! ------------------
        //! compute curvature
        //! ------------------
        int curvatureType = 1;      //! gaussian
        theMeshToInflate_new->computeDiscreteCurvature(curvatureType);

        //! -----------------------
        //! compute shrink factors
        //! -----------------------
        QMap<int,double> shrinkFactors;
        this->computeShrinkFactor(theMeshToInflate_new,shrinkFactors);

        //! -----------------------------------------
        //! create the displacements vectorial field
        //! -----------------------------------------
        QMap<int,QList<double>> normals = theMeshToInflate_new->myNodeNormals;
        if(normals.isEmpty()) return false;

        //! -------------------
        //! smooth the normals
        //! -------------------
        const int NbSmoothingSteps = 50;
        this->fieldSmoother(normals,theMeshToInflate_new,5,NbSmoothingSteps);

        QMap<int,gp_Vec> displacementsField;
        for(QMap<int,QList<double>>::const_iterator it = normals.cbegin(); it!= normals.cend(); ++it)
        {
            int globalNodeID = it.key();
            const QList<double> &curNormal = it.value();

            //! -----------------------------------------------------
            //! cutoff is 1.0 since all the points will be displaced
            //! component of the displacement field
            //! -----------------------------------------------------
            double cutoff = 1.0;

            double shrinkFactor = shrinkFactors.value(globalNodeID);
            double vx = -curNormal.at(0)*displacement*cutoff*shrinkFactor;
            double vy = -curNormal.at(1)*displacement*cutoff*shrinkFactor;
            double vz = -curNormal.at(2)*displacement*cutoff*shrinkFactor;            

            bool correctBoundary = true;
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
            gp_Vec displacementVector(vx,vy,vz);
            displacementsField.insert(globalNodeID,displacementVector);
        }

        //! ------------------------------
        //! smooth the displacement field
        //! ------------------------------
        this->smoothDisplacementField(displacementsField,normals,theMeshToInflate_new);

        //! -------------------------------------------
        //! displace the prismatic and the volume mesh
        //! -------------------------------------------
        theMeshToInflate_new->displaceMySelf(displacementsField);
        theInflatedMeshes<<theMeshToInflate_new;

        //! ----------------------------------------------
        //! translation of nodes without mesh compression
        //! ----------------------------------------------
        //preInflationVolumeMeshDS->displaceMySelf(displacementsField);

        //! ------------------------------------------------------------
        //! experimental - for testing as rigid as possible deformation
        //! this substitutes for the "displaceMySelf" method
        //! ------------------------------------------------------------
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
        bool canApplyDeformation = false;   // dummy initialization;

        for(int P = myNbLayers; P>0; P--)
        {
            canApplyDeformation = false;
            constrainedGlobalNodeIDs.clear();
            groupOfMovingNodes.clear();

            meshSelector->setSteps(P);
            bool isDone = meshSelector->getElements(selectedVolumeMeshDS);

            //! ---------------------------------------------------------------
            //! for some reason cannot retrieve a volume submesh using P steps
            //! ---------------------------------------------------------------
            if(!isDone) continue;

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
            if(constrainedGlobalNodeIDs.size()>3)
            {
                canApplyDeformation = true;
                break;
            }
        }

        if(!canApplyDeformation) return false;

        //! -------------------------------------------------------
        //! mode: "0" harmonic/biharmonic "1" as rigig as possible
        //! -------------------------------------------------------
        int mode = 1;
        preInflationVolumeMeshDS->displaceMySelf_asRigidAsPossible(displacementsField,constrainedGlobalNodeIDs,mode);
        //preInflationVolumeMeshDS->displaceMySelf_asRigidAsPossible(displacementsField,groupOfMovingNodes,mode);
        //! ------------
        //! end testing
        //! ------------

        //! ------------------------------------
        //! update the "old" displacement field
        //! ------------------------------------
        theMeshToInflate_old = theMeshToInflate_new;
    }

    //! ---------------------------------------------------------------
    //! insert the deformed (compressed) volume mesh into the database
    //! ---------------------------------------------------------------
    myMeshDB->ArrayOfMeshDS.insert(myBodyIndex,preInflationVolumeMeshDS);                   // to be commented
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

//! -------------------------------------------------------
//! function: generateTetLayers
//! details:  return the boundary mesh made of tetrahedral
//!           and the deformed surface mesh
//! -------------------------------------------------------
void prismaticLayer::generateTetLayers(occHandle(Ng_MeshVS_DataSource3D) &meshAtWalls,
                                       occHandle(Ng_MeshVS_DataSourceFace) &lastInflatedMesh)
{
    cout<<"prismaticLayers::generateTetLayers()->____function called____"<<endl;

    //! -----------------------
    //! number of inverted tet
    //! -----------------------
    int NumberOfInvertedTet = 0;

    //! -----------------------------
    //! the volume elements at walls
    //! -----------------------------
    std::vector<meshElementByCoords> volumeElementsAtWalls;

    //! --------------------------------------
    //! decide which node should be displaced
    //! --------------------------------------
    this->computeVecFieldCutOff(myLockBoundary);

    //! ---------------------------
    //! make a copy the outer mesh
    //! ---------------------------
    occHandle(Ng_MeshVS_DataSourceFace) theMeshToInflate = new Ng_MeshVS_DataSourceFace(myPrismaticFacesSumMeshDS);
    //occHandle(Ng_MeshVS_DataSourceFace) theMeshToInflate = new Ng_MeshVS_DataSourceFace(myOverallSumMeshDS);

    if(myPrismaticFacesSumMeshDS->myBoundarySegments.isEmpty()) myPrismaticFacesSumMeshDS->computeFreeMeshSegments();

    //! --------------------
    //! generate the layers
    //! --------------------
    for(int n=1; n<=myNbLayers; n++)
    {
        //! ----------------------------
        //! the current layer thickness
        //! ----------------------------
        double displacement = myLayerThickness.at(n-1);

        //! ------------------------------------------------------------
        //! collapsed into a method the generation of one layer of tets
        //! ------------------------------------------------------------
        this->generateOneTetLayer(theMeshToInflate,displacement,volumeElementsAtWalls);
    }

    cout<<"Number of inverted tet: "<<NumberOfInvertedTet<<endl;

    //! -------------------------------
    //! the last modified surface mesh
    //! -------------------------------
    //lastInflatedMesh = theMeshToInflate;
    lastInflatedMesh = myOverallSumMeshDS;

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
        cout<<"prismaticLayer::generateMeshAtWalls()->____start generating prismatic 3D elements____"<<endl;
        isDone = this->buildPrismaticElements(listOfInflatedMeshes,meshAtWalls);
        if(!isDone) return false;

        //! -----------------------
        //! the last inflated mesh
        //! -----------------------
        lastInflatedMesh = listOfInflatedMeshes.last();
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
            int localNodeID = badNodes.at(k)+1;
            std::vector<double> old_coords = aMeshDS_old->getNodeCoordinates(localNodeID);
            aMeshDS->changeNodeCoords(localNodeID,old_coords);
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
            int localNodeID = badNodes_inflatedMesh.at(k)+1;
            std::vector<double> old_coords = aMeshDS_old->getNodeCoordinates(localNodeID);
            aMeshDS->changeNodeCoords(localNodeID,old_coords);
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
                                         std::vector<meshElementByCoords> &volumeElementsAtWalls)
{
    //! ----------------------------------------------------------------------
    //! calculate the discrete curvature and the curvature gradient amplitude
    //! compute the shrink factor
    //! ----------------------------------------------------------------------
    int curvatureType = 1;  //! gaussian curvature
    theMeshToInflate->computeDiscreteCurvature(curvatureType);
    QMap<int,double> shrinkFactors;
    this->computeShrinkFactor(theMeshToInflate,shrinkFactors);

    //! ---------------------------------------
    //! calculate nodal the displacement field
    //! ---------------------------------------
    theMeshToInflate->computeNormalAtElements();
    theMeshToInflate->computeNormalAtNodes();
    theMeshToInflate->computeFreeMeshSegments();

    QMap<int,QList<double>> normals = theMeshToInflate->myNodeNormals;

    //! --------------------------------------
    //! smooth the guiding vectors directions
    //! --------------------------------------
    int NbSmoothingSteps = myNbGuidingVectorSmoothingSteps;
    double cs = myCurvatureSensitivityForGuidingVectorsSmoothing;
    this->fieldSmoother(normals,theMeshToInflate,cs,NbSmoothingSteps);

    QMap<int,gp_Vec> displacementsField;
    for(QMap<int,QList<double>>::const_iterator it = normals.cbegin(); it!= normals.cend(); ++it)
    {
        int globalNodeID = it.key();
        const QList<double> &curNormal = it.value();

        //! -------------------------
        //! jump over boundary nodes
        //! -------------------------
        if(theMeshToInflate->myBoundaryPoints.contains(globalNodeID)) continue;

        //! --------------
        //! shrink factor
        //! --------------
        double shrinkFactor = shrinkFactors.value(globalNodeID);
        double cutOff = myLayerHCutOff.value(globalNodeID);
        double C = displacement*shrinkFactor*cutOff;
        double vx = -curNormal.at(0)*C;
        double vy = -curNormal.at(1)*C;
        double vz = -curNormal.at(2)*C;

        //! -------------------
        //! displacement field
        //! -------------------
        gp_Vec normVec(vx,vy,vz);
        displacementsField.insert(globalNodeID,normVec);
    }

    //! ------------------------------
    //! smooth the displacement field
    //! ------------------------------
    this->smoothDisplacementField(displacementsField,normals,theMeshToInflate);

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
        const gp_Vec &vecDispl = displacementsField.value(globalNodeID);

        double x_displ = pc[0]+vecDispl.X();
        double y_displ = pc[1]+vecDispl.Y();
        double z_displ = pc[2]+vecDispl.Z();

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
            if(V<0.0)
            {
                //NumberOfInvertedTet++;
                cout<<"____INVERTED ELEMENT: CANNOT MOVE____"<<endl;
                canMove = false;
                break;

                //this->invertTet(aVolElement);
            }

            //! -------------------------------------------------------
            //! check if the volume element is valid, and can be added
            //! -------------------------------------------------------
            double q0, q1, q2;
            aTetQuality.setPoints(aVolElement.getPoints());
            aTetQuality.getQualityMeasure(q0,q1,q2,V);

            const double QUALITY_LIMIT = 0.01;
            if(q2<QUALITY_LIMIT)
            {
                cout<<"____CANNOT MOVE____"<<endl;
                canMove = false;
                break;
            }
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

            //! ------------------------
            //! check self intersection
            //! ------------------------
            //bool selfIntersection = this->checkSelfIntersection(myOverallSumMeshDS,myOverallSumMeshDS_old,true);
            //if(selfIntersection == false)
            //{
            //    for(int i=0; i<vecElementsToAdd.size(); i++)
            //    {
            //        volumeElementsAtWalls.push_back(vecElementsToAdd[i]);
            //    }
            //}
            //! ---------------------------
            //! check mutual intersections
            //! ---------------------------
            //bool mutualIntersection = this->checkMutualIntersection(myOverallSumMeshDS,myOverallSumMeshDS_old,true);
            //if(mutualIntersection==false)
            //{
            //    for(int i=0; i<vecElementsToAdd.size(); i++)
            //    {
            //        volumeElementsAtWalls.push_back(vecElementsToAdd[i]);
            //    }
            //}

            //! -------------------
            //! standard treatment
            //! -------------------
            for(int i=0; i<vecElementsToAdd.size(); i++)
            {
                volumeElementsAtWalls.push_back(vecElementsToAdd[i]);
            }
        }
    }
}

//! ------------------------
//! function: fieldSmoother
//! details:
//! ------------------------
void prismaticLayer::fieldSmoother(QMap<int,QList<double>> &field,
                                   const occHandle(Ng_MeshVS_DataSourceFace) &aMeshDS,
                                   double k,
                                   int NbSteps)
{
    cout<<"prismaticLayer::fieldSmoother()->____function called____"<<endl;

    if(aMeshDS->myCurvature.isEmpty())
    {
        int mode = 1;       //! use gaussian curvature
        aMeshDS->computeDiscreteCurvature(mode);
    }

    //! ----------------
    //! smoothing steps
    //! ----------------
    for(int q = 1; q<=NbSteps; q++)
    {
        QMap<int,QList<double>> smoothedField;
        for(QMap<int,QList<double>>::iterator it = field.begin(); it!= field.end(); ++it)
        {
            int globalNodeID = it.key();
            int localNodeID = aMeshDS->myNodesMap.FindIndex(globalNodeID);

            const std::vector<double> &P = aMeshDS->getNodeCoordinates(localNodeID);

            //! -------------------------------------------
            //! value of the field at the current position
            //! -------------------------------------------
            const QList<double> &localValue = it.value();

            //! --------------------------
            //! get the surrounding nodes
            //! --------------------------
            std::vector<int> surroundingNodes;
            QList<int> elements = aMeshDS->myNodeToElements.value(localNodeID);
            for(int i=0; i<elements.length(); i++)
            {
                int localElementID = elements[i];
                int globalElementID = aMeshDS->myElementsMap.FindKey(localElementID);
                int NbNodes, buf[20];
                TColStd_Array1OfInteger nodeIDs(*buf,1,20);
                aMeshDS->GetNodesByElement(globalElementID,nodeIDs,NbNodes);
                for(int n=1; n<=NbNodes; n++)
                {
                    int globalNodeID_ = nodeIDs(n);

                    if(globalNodeID_==globalNodeID) continue;
                    if(std::find(surroundingNodes.begin(),surroundingNodes.end(),globalNodeID_)==surroundingNodes.end())
                        surroundingNodes.push_back(globalNodeID_);
                }
            }
            //bool isDone = aMeshDS->getSurroundingNodes(localNodeID,surroundingNodes);
            //if(!isDone) exit(22);

            int NbSurrounding = int(surroundingNodes.size());

            //cout<<"@____begin smoothing at node____@"<<endl;
            //cout<<"@____nodeID: "<<globalNodeID<<" number of surrounding: "<<NbSurrounding<<"____@"<<endl;

            double denom = 0;
            double num_x, num_y, num_z;
            num_x = num_y = num_z = 0;

            for(int n=0; n<NbSurrounding; n++)
            {
                //cout<<"____n-th surrounding: "<<n+1<<"____"<<endl;
                int globalNodeID_surrounding = surroundingNodes[n];
                int localNodeID_surrounding = aMeshDS->myNodesMap.FindIndex(globalNodeID_surrounding);
                const std::vector<double> &S = aMeshDS->getNodeCoordinates(localNodeID_surrounding);
                double d = sqrt(pow(S[0]-P[0],2)+pow(S[1]-P[1],2)+pow(S[2]-P[2],2));
                //cout<<"____distance ("<<globalNodeID<<", "<<globalNodeID_surrounding<<") ="<<d<<"____"<<endl;
                denom += 1.0/d;
                const QList<double> &localValue_surrounding = field.value(globalNodeID_surrounding);
                num_x += (1.0/d)*localValue_surrounding[0];
                num_y += (1.0/d)*localValue_surrounding[1];
                num_z += (1.0/d)*localValue_surrounding[2];
                //cout<<"____"<<num_x<<", "<<num_y<<", "<<num_z<<"____"<<endl;
            }
            if(denom == 0)
            {
                cerr<<"____zero denominator during laplacian smoothing____"<<endl;
                exit(20);
            }

            double curvature = aMeshDS->myCurvature.value(globalNodeID);
            const double maxOmega = 1.0;
            double omega = maxOmega*(1-1/(1+exp(-k*sqrt(fabs(curvature)))));
            if(curvature<0) omega = maxOmega-omega;

            double snx = omega*localValue[0]+(1-omega)*(num_x/denom);
            double sny = omega*localValue[1]+(1-omega)*(num_y/denom);
            double snz = omega*localValue[2]+(1-omega)*(num_z/denom);

            double l = sqrt(snx*snx+sny*sny+snz*snz);
            snx /= l;
            sny /= l;
            snz /= l;

            //cout<<"____old field value ("<<field.value(globalNodeID).at(0)<<", "<<field.value(globalNodeID).at(1)<<", "<<field.value(globalNodeID).at(2)<<")____"<<endl;
            //cout<<"____new field value ("<<snx<<", "<<sny<<", "<<snz<<")____"<<endl;
            //cout<<"@____end smoothing at node____@"<<endl<<endl;

            QList<double> localFieldValueSmoothed;
            localFieldValueSmoothed<<snx<<sny<<snz;
            smoothedField.insert(globalNodeID,localFieldValueSmoothed);
        }

        //! -------------------------
        //! replace the field values
        //! -------------------------
        for(QMap<int,QList<double>>::iterator it = field.begin(); it!= field.end(); ++it)
        {
            int globalNodeID = it.key();
            it.value() = smoothedField.value(globalNodeID);
        }
    }
}

//! ----------------------------------
//! function: smoothDisplacementField
//! details:  helper
//! ----------------------------------
void prismaticLayer::smoothDisplacementField(QMap<int,gp_Vec> &displacementsField,
                                             const QMap<int,QList<double>> &normals,
                                             const occHandle(Ng_MeshVS_DataSourceFace) &theMeshToInflate)
{
    cout<<"smoothDisplacementField()->____function called____"<<endl;

    //! ---------------------
    //! smoothing parameters
    //! ---------------------
    int NbSmoothingSteps = myNbLayerThicknessSmoothingSteps;
    double cs = myCurvatureSensitivityForThicknessSmoothing;

    //! ----------------------------------
    //! copy into QMap<int,QList<double>>
    //! ----------------------------------
    QMap<int,QList<double>> field;
    for(QMap<int,gp_Vec>::iterator it = displacementsField.begin(); it!=displacementsField.end(); it++)
    {
        int globalNodeID = it.key();
        gp_Vec v = it.value();

        QList<double> fieldValue;
        double d = sqrt(v.X()*v.X()+v.Y()*v.Y()+v.Z()*v.Z());
        fieldValue<<d<<0.01<<0.01;
        field.insert(globalNodeID,fieldValue);
    }

    //! -------
    //! smooth
    //! -------
    this->fieldSmoother(field,theMeshToInflate,cs,NbSmoothingSteps);

    //! -------------------------------
    //! converto into QMap<int,gp_Vec>
    //! -------------------------------
    for(QMap<int,QList<double>>::iterator it = field.begin(); it!=field.end(); it++)
    {
        int globalNodeID = it.key();
        double smoothedValue = it.value().at(0);
        QList<double> curNormal = normals.value(globalNodeID);
        double vx = -smoothedValue*curNormal.at(0);
        double vy = -smoothedValue*curNormal.at(1);
        double vz = -smoothedValue*curNormal.at(2);
        gp_Vec vec(vx,vy,vz);
        displacementsField.insert(globalNodeID,vec);
    }
}
