//! ----------------
//! custom includes
//! ----------------
#include "particlesinfieldssolver.h"
#include "particlesemitter.h"
#include "poissonsolver.h"
#include "ng_meshvs_datasource3d.h"
#include "simulationnodeclass.h"
#include "simulationdatabase.h"
#include "indexedmapofmeshdatasources.h"
#include "geometrytag.h"

//! ---
//! Qt
//! ---
#include <QVector>

//! -------
//! libigl
//! -------
#include <libigl/include/igl/in_element.h>
#include <libigl/include/igl/grad.h>
#include <libigl/include/igl/massmatrix.h>
#include <libigl/include/igl/invert_diag.h>

//! -------------------------------------------------
//! function: tetVol
//! details:  return the determinant of a 4x4 matrix
//! -------------------------------------------------
double tetVol(double *P0,double* P1, double *P2, double *P3)
{
    Eigen::MatrixXd J(4,4);
    J<<P0[0],P0[1],P0[2],1.0,
            P1[0],P1[1],P1[2],1.0,
            P2[0],P2[1],P2[2],1.0,
            P3[0],P3[1],P3[2],1.0;
    return (1.0/6.0)*J.determinant();
}

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
particlesInFieldsSolver::particlesInFieldsSolver(simulationDataBase *sDB, QStandardItem *simulationRoot, QObject *parent): QObject(parent)
{
    mySimulationRoot = simulationRoot;
    mySDB =  sDB;
    this->init(mySimulationRoot);

    //! ---------------------
    //! computational domain
    //! ---------------------
    occHandle(Ng_MeshVS_DataSource3D) volumeMesh = occHandle(Ng_MeshVS_DataSource3D)::DownCast(sDB->ArrayOfMeshDS.first());

    //! ------------------------
    //! setup mesh: mesh points
    //! ------------------------
    int NbPoints = volumeMesh->GetAllNodes().Extent();
    V.resize(NbPoints,3);
    for(int row=1; row<=NbPoints; row++)
    {
        const std::vector<double> &P = volumeMesh->getNodeCoordinates(row);
        for(int col=1; col<=3; col++) V(row-1,col-1) = P[col];
    }

    //! --------------
    //! mesh elements
    //! --------------
    int NbElements = volumeMesh->GetAllElements().Extent();
    T.resize(NbElements,4);
    int NbNodes = -1, buf[4];
    TColStd_Array1OfInteger nodeIDs(*buf,1,4);
    TColStd_IndexedMapOfInteger nodesMap = volumeMesh->myNodesMap;
    for(int i=1; i<=NbElements; i++)
    {
        int globalElementID = volumeMesh->myElementsMap.FindIndex(i);
        volumeMesh->GetNodesByElement(globalElementID,nodeIDs,NbNodes);

        int row = i-1;
        T(row,0) = nodesMap.FindIndex(nodeIDs(1))-1;
        T(row,1) = nodesMap.FindIndex(nodeIDs(2))-1;
        T(row,2) = nodesMap.FindIndex(nodeIDs(3))-1;
        T(row,3) = nodesMap.FindIndex(nodeIDs(4))-1;

        std::vector<double> P0 = volumeMesh->getNodeCoordinates(nodeIDs(1));
        std::vector<double> P1 = volumeMesh->getNodeCoordinates(nodeIDs(2));
        std::vector<double> P2 = volumeMesh->getNodeCoordinates(nodeIDs(3));
        std::vector<double> P3 = volumeMesh->getNodeCoordinates(nodeIDs(4));

        myVolumes(row) = tetVol(P0.data(),P1.data(),P2.data(),P3.data());
    }

    //! ----------------------------------------
    //! initialize AABB tree for point location
    //! ----------------------------------------
    mySearchTree.init(V,T);

    //! ----------------------------------------
    //! allocate space for density distribution
    //! and reset distribution
    //! ----------------------------------------
    myDensity.resize(V.rows());
    for(int i=0; i<V.rows(); i++) myDensity(i) = 0.0;

    //! --------------------------------
    //! the same for the electric field
    //! --------------------------------
    myElectricField.resize(V.rows(),3);
    for(int i=0; i<V.rows(); i++) myElectricField(i,0) = myElectricField(i,1) = myElectricField(i,2) = 0.0;

    //! --------------------------
    //! create the Poisson solver
    //! --------------------------
    myPoissonSolver = new PoissonSolver(volumeMesh,this);

    //! -----------------------------
    //! init the number of particles
    //! -----------------------------
    myNbParticles = 0;

    //! ----------------------------
    //! mass matrix and its inverse
    //! ----------------------------
    igl::massmatrix(V,T,igl::MASSMATRIX_TYPE_BARYCENTRIC, M);
    igl::invert_diag(M,Minv);

    //! ------------------
    //! gradient operator
    //! ------------------
    igl::grad(V,T,G);
}

//! ---------------
//! function: init
//! details:
//! ---------------
bool particlesInFieldsSolver::init(QStandardItem* simulationRoot)
{
    //! ----------------
    //! simulation root
    //! ----------------
    mySimulationRoot = simulationRoot;

    //! -----------------------------
    //! init the number of particles
    //! -----------------------------
    myNbParticles = 0;

    //! ---------------------------------
    //! retrieve the "Analysis settings"
    //! ---------------------------------
    SimulationNodeClass *node = mySimulationRoot->data(Qt::UserRole).value<SimulationNodeClass*>();
    if(node->getType()!= SimulationNodeClass::nodeType_particlesInFieldsAnalysis) return false;

    //! ----------------------------
    //! simulation time - time step
    //! ----------------------------
    myFinalTime = node->getPropertyValue<double>("Step end time");
    myTimeStep = node->getPropertyValue<double>("Time step size");

    //! --------------------------------
    //! for each face a potential patch
    //! --------------------------------
    std::map<int,occHandle(Ng_MeshVS_DataSourceFace)> patches;

    //! --------------
    //! scan the tree
    //! --------------
    for(int n=1; n<simulationRoot->rowCount()-1; n++)
    {
        QStandardItem *item = simulationRoot->child(n,0);
        SimulationNodeClass *curNode = item->data(Qt::UserRole).value<SimulationNodeClass*>();
        switch(curNode->getType())
        {
        case SimulationNodeClass::nodeType_electrostaticPotential:
        {
            const std::vector<GeometryTag> &vecLoc = curNode->getPropertyValue<std::vector<GeometryTag>>("Tags");
            int NbPatches = vecLoc.size();
            QList<occHandle(Ng_MeshVS_DataSourceFace)> listOfFaceMeshDS;
            for(int n=0; n<NbPatches; n++)
            {
                std::pair<int,occHandle(Ng_MeshVS_DataSourceFace)> apair;
                int faceNr = vecLoc[0].subTopNr;
                const occHandle(Ng_MeshVS_DataSourceFace) &aFaceMeshDS =
                        occHandle(Ng_MeshVS_DataSourceFace)::DownCast(mySDB->ArrayOfMeshDSOnFaces.getValue(1,faceNr));
                patches.insert(std::make_pair(faceNr,aFaceMeshDS));

                //! -------------------------------
                //! sum the face mesh data sources
                //! -------------------------------
                listOfFaceMeshDS.append(aFaceMeshDS);
            }

            double val = curNode->getPropertyValue<double>("Potential");
            bool isEmitter = curNode->getPropertyValue<bool>("Is emitter");

            //! -----------------------
            //! the wall is an emitter
            //! -----------------------
            if(isEmitter)
            {
                double intensity = curNode->getPropertyValue<double>("Intensity");
                double mass = curNode->getPropertyValue<double>("Particle mass");
                double charge = curNode->getPropertyValue<double>("Electric charge");

                occHandle(Ng_MeshVS_DataSourceFace) emitterFaceMeshDS = new Ng_MeshVS_DataSourceFace(listOfFaceMeshDS);
                std::shared_ptr<particlesEmitter> anEmitter(new particlesEmitter(emitterFaceMeshDS,intensity));
                myEmitters.push_back(anEmitter);
            }
        }
            break;
        }

        //! ------------------------------------------------
        //! define the patches and initialize the potential
        //! on all the surface mesh nodes
        //! ------------------------------------------------
        myPoissonSolver->definePatches(patches);
        myPoissonSolver->initPotentialOnBoundary(0.0);
    }
    return true;
}

//! --------------
//! function: run
//! details:
//! --------------
void particlesInFieldsSolver::run()
{
    cout<<"@ -------------------------------------"<<endl;
    cout<<"@ Simulation started @ time "<<this->getTimeStamp()<<endl;
    cout<<"@ -------------------------------------"<<endl;

    //! -----------------------------------------------------
    //! init the simulation time and the number of particles
    //! -----------------------------------------------------
    myNbParticles = 0;
    double time = -myTimeStep;
    for(;;)
    {
        if(time>myFinalTime)
        {
            cout<<"@ -------------------------------"<<endl;
            cout<<"@ Simulation ended @ time "<<this->getTimeStamp()<<endl;
            cout<<"@ -------------------------------"<<endl;
        }
        break;
        time += myTimeStep;

        //! ---------------------------------------
        //! load the particles at the current time
        //! ---------------------------------------
        for(int i=0; i<myEmitters.size(); i++)
        {
            std::shared_ptr<particlesEmitter> anEmitter = myEmitters.at(i);
            std::vector<particle> newParticles;
            anEmitter->generateParticles(newParticles,time);
            for(int n=0; n<newParticles.size(); n++)
            {
                myNbParticles++;
                const particle &aP = newParticles[n];
                std::shared_ptr<particle> aParticle(new particle(aP));
                std::pair<int,std::shared_ptr<particle>> apair;
                apair.first= myNbParticles;
                apair.second = aParticle;
                myParticles.insert(apair);
            }
        }

        //! -------------------------
        //! compute density at nodes
        //! -------------------------
        this->computeChargeDensityAtNodes();

        //! --------------------------------
        //! compute electrostatic potential
        //! --------------------------------
        //... to do

        //! ---------------------------------------------
        //! compute electric field at particles position
        //! ---------------------------------------------


        //! ---------------
        //! move particles
        //! ---------------
        this->moveParticles();

    }
}

//! --------------------------------
//! function: computeDensityAtNodes
//! details:
//! --------------------------------
void particlesInFieldsSolver::computeChargeDensityAtNodes()
{
    Eigen::VectorXd charge(V.rows());

    const int sequences[4][3] { {1,2,3}, {3,2,0}, {0,1,3}, {0,1,2} };

    //! ------------------
    //! reset the density
    //! ------------------
    //this->initChargeDensityAtNodes(0.0);

    Eigen::MatrixXd Q(1,3);
    for(std::map<int,std::shared_ptr<particle>>::iterator it = myParticles.begin(); it!= myParticles.end(); it++)
    {
        std::pair<int,std::shared_ptr<particle>> apair = *it;
        std::shared_ptr<particle> aParticle = apair.second;
        double xP = aParticle->x[0];
        double yP = aParticle->x[1];
        double zP = aParticle->x[2];

        //! ---------------------------------------------------
        //! search for the tetrahedron the particle belongs to
        //! ---------------------------------------------------
        Q(0,0) = xP; Q(0,1) = yP; Q(0,2) = zP;
        Eigen::VectorXi I;
        igl::in_element(V,T,Q,mySearchTree,I);
        if(I.rows()==0)
        {
            cout<<"particlesInFieldsSolver::computeChargeDensityAtNodes()->____particle outside domain____"<<endl;
            continue;
        }

        int tetNumber = I(0);
        double cellVolume = myVolumes(tetNumber);

        //! ------------------------------
        //! particle point and tet points
        //! ------------------------------
        double P[4][3];
        P[3][0] = xP; P[3][1] = yP; P[3][2] = zP;   // position of the particle

        //! -------------------------------------------------------------------
        //! iterate over the subdivision tetrahedra: use volume coordinates
        //! for distributing the electric charge onto the nodes of the current
        //! mesh element
        //! -------------------------------------------------------------------
        for(int n=0; n<=3; n++)
        {
            for(int i=0; i<=2; i++)
            {
                int index = sequences[n][i];
                int pnum =  T(tetNumber,index);
                for(int c=0; c<3; c++) P[i][c] = V(pnum,c);
            }
            double r = tetVol(&P[n][0],&P[n][1],&P[n][2],&P[n][3]);
            r /= cellVolume;

            //! ---------------
            //! charge at node
            //! ---------------
            charge(T(tetNumber,n)) += aParticle->q*r;
        }

        //! --------------------------------------
        //! charge density: divide by cell volume
        //! --------------------------------------
        myDensity = (Minv*myDensity).eval();
    }
}

//! ------------------------
//! function: moveParticles
//! details:
//! ------------------------
void particlesInFieldsSolver::moveParticles()
{
    for(std::map<int,std::shared_ptr<particle>>::iterator it = myParticles.begin(); it!= myParticles.end(); it++)
    {
        std::pair<int,shared_ptr<particle>> apair = *it;
        double x = apair.second->x[0];
        double y = apair.second->x[1];
        double z = apair.second->x[2];

        double vx = apair.second->v[0];
        double vy = apair.second->v[1];
        double vz = apair.second->v[2];

        //! ------------
        //! Runge-Kutta
        //! ------------

    }
}

//! -------------------------------
//! function: computeElectricField
//! details:
//! -------------------------------
void particlesInFieldsSolver::computeElectricField()
{
    myElectricField = -Eigen::Map<const Eigen::MatrixXd>((G*myElectrostaticPotential).eval().data(),T.rows(),3);
}

//! -----------------------------
//! function: getGradientAtPoint
//! details:
//! -----------------------------
void particlesInFieldsSolver::getElectricFieldAtPoint(double x, double y, double z, double &Ex, double &Ez, double &Ey)
{
    int tetNumber = -1;
    bool isFound = this->getCellOfAPoint(x,y,z,tetNumber);
    if(!isFound) Ex = Ey = Ez = 0.0;

    Ex = myElectricField(tetNumber,0);
    Ey = myElectricField(tetNumber,1);
    Ez = myElectricField(tetNumber,2);
}

//! --------------------------
//! function: getCellOfAPoint
//! details:  helper
//! --------------------------
bool particlesInFieldsSolver::getCellOfAPoint(double x, double y, double z, int& tetNumber)
{
    Eigen::MatrixXd Q;
    Q<<x,y,z;
    Eigen::VectorXi I;
    igl::in_element(V,T,Q,mySearchTree,I);
    if(I.rows()==0) return false;
    {
        cout<<"particlesInFieldsSolver::computeChargeDensityAtNodes()->____particle outside domain____"<<endl;
        return false;
    }
    tetNumber = I(0);
    return true;
}

//! -----------------------
//! function: getTimeStamp
//! details:
//! -----------------------
std::string particlesInFieldsSolver::getTimeStamp()
{
    auto theNow = std::chrono::system_clock::now();
    std::time_t theNow_time = std::chrono::system_clock::to_time_t(theNow);
    std::string timeStamp(std::ctime(&theNow_time));
    return timeStamp;
}
