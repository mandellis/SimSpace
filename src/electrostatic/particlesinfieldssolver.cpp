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
#include <meshelementbycoords.h>
#include "customtablemodel.h"
#include "tabulardatacolumns.h"

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

//! ----
//! C++
//! ----
#include <utility>

//! -------------------------
//! function: writeInputFile
//! details:
//! -------------------------
bool particlesInFieldsSolver::writeInputFile(simulationDataBase *sDB, QStandardItem *simulationRoot, const string &inputFilePath)
{
    //! --------------
    //! open a stream
    //! --------------
    ofstream os(inputFilePath);
    if(os.is_open()==false) return false;

    //! ---------------------
    //! write the body index
    //! ---------------------
    int bodyIndex = 1;
    os<<"BODYINDEX "<<bodyIndex<<endl;

    //! ----------------------
    //! write the volume mesh
    //! ----------------------
    os<<"VOLUMEMESH"<<endl;
    occHandle(Ng_MeshVS_DataSource3D) volumeMesh = occHandle(Ng_MeshVS_DataSource3D)::DownCast(sDB->ArrayOfMeshDS.value(bodyIndex));
    volumeMesh->writeIntoStream(os);

    //! --------------------------
    //! write the number of faces
    //! --------------------------
    int NbFaces = sDB->MapOfBodyTopologyMap.value(bodyIndex).faceMap.Extent();
    os<<"NUMBEROFFACES "<<NbFaces<<endl;

    for(int faceNr = 1; faceNr<=NbFaces; faceNr++)
    {
        //! ----------------------
        //! write the face number
        //! ----------------------
        os<<"FACENUMBER "<<faceNr<<endl;

        //! -----------------------------------------------------------------------------
        //! write the definition of the face mesh through the ID of the surface elements
        //! -----------------------------------------------------------------------------
        occHandle(Ng_MeshVS_DataSourceFace) aFaceMeshDS = occHandle(Ng_MeshVS_DataSourceFace)::DownCast(sDB->ArrayOfMeshDSOnFaces.getValue(bodyIndex,faceNr));
        int NbFaceElements = aFaceMeshDS->GetAllElements().Extent();
        os<<"NUMBEROFELEMENTS "<<NbFaceElements<<endl;

        for(TColStd_MapIteratorOfPackedMapOfInteger it(aFaceMeshDS->GetAllElements()); it.More(); it.Next())
        {
            int NbNodes, buf[8];
            TColStd_Array1OfInteger nodeIDs(*buf,1,8);
            aFaceMeshDS->GetNodesByElement(it.Key(),nodeIDs,NbNodes);

            switch(NbNodes)
            {
            case 3:
            {
                os<<"TRIG"<<endl;
                os<<nodeIDs(1)<<"\t"<<nodeIDs(2)<<"\t"<<nodeIDs(3)<<endl;
            }
                break;
            case 4:
            {
                os<<"QUAD"<<endl;
                os<<nodeIDs(1)<<"\t"<<nodeIDs(2)<<"\t"<<nodeIDs(3)<<"\t"<<nodeIDs(4)<<endl;
            }
                break;
            }
        }
    }

    //! ----------------------------------------
    //! write the number of boundary conditions
    //! ----------------------------------------
    int NbBC = simulationRoot->rowCount()-2;
    os<<"BOUNDARYCONDITIONSNUMBER "<<NbBC<<endl;

    //! -------------------------------------------------------------
    //! now scan the simulation root (start after "Analysis setting"
    //! and finish an item before "Solution"). Count also the number
    //! of emitter associated with the boundary condition
    //! -------------------------------------------------------------
    int NbEmitters = 0;
    for(int n=1; n<simulationRoot->rowCount()-1; n++)
    {
        QStandardItem *item = simulationRoot->child(n,0);
        SimulationNodeClass *curNode = item->data(Qt::UserRole).value<SimulationNodeClass*>();
        if(curNode->getPropertyItem("Emitter")!=Q_NULLPTR)
        {
            bool isEmitter = curNode->getPropertyValue<bool>("Emitter");
            if(isEmitter) NbEmitters++;
        }
        switch(curNode->getType())
        {
        case SimulationNodeClass::nodeType_electrostaticPotential:
        {
            //! ---------------
            //! write the type
            //! ---------------
            os<<"ELECTROSTATICPOTENTIAL"<<endl;

            //! --------------------------------------
            //! write the time behavior and the value
            //! --------------------------------------
            os<<"CONSTANT"<<endl;
            double val = curNode->getPropertyValue<double>("Potential");
            os<<val<<endl;

            //! -----------------------------------------
            //! build a composite face mesh on whith the
            //! current boundary condition is applied
            //! -----------------------------------------
            const QVector<GeometryTag> &vecLoc = curNode->getPropertyValue<QVector<GeometryTag>>("Tags");
            QList<occHandle(Ng_MeshVS_DataSourceFace)> listOfFaceMeshDS;
            for(int n=0; n<vecLoc.length(); n++)
            {
                int faceNr = vecLoc[n].subTopNr;
                occHandle(Ng_MeshVS_DataSourceFace) aFaceMeshDS = occHandle(Ng_MeshVS_DataSourceFace)::DownCast(sDB->ArrayOfMeshDSOnFaces.getValue(bodyIndex,faceNr));
                listOfFaceMeshDS.append(aFaceMeshDS);
            }
            occHandle(Ng_MeshVS_DataSourceFace) compositeFaceDS = new Ng_MeshVS_DataSourceFace(listOfFaceMeshDS);

            //! ---------------------------
            //! write the surface elements
            //! ---------------------------
            os<<"BOUNDARYPATCHDEFINITION"<<endl;

            int NbFaceElements = compositeFaceDS->GetAllElements().Extent();
            os<<"NUMBEROFELEMENTS "<<NbFaceElements<<endl;
            for(TColStd_MapIteratorOfPackedMapOfInteger it(compositeFaceDS->GetAllElements()); it.More(); it.Next())
            {
                int NbNodes, buf[8];
                TColStd_Array1OfInteger nodeIDs(*buf,1,8);
                compositeFaceDS->GetNodesByElement(it.Key(),nodeIDs,NbNodes);
                switch(NbNodes)
                {
                case 3: os<<"TRIG"<<endl; os<<nodeIDs(1)<<"\t"<<nodeIDs(2)<<"\t"<<nodeIDs(3)<<endl; break;
                case 4: os<<"QUAD"<<endl; os<<nodeIDs(1)<<"\t"<<nodeIDs(2)<<"\t"<<nodeIDs(3)<<"\t"<<nodeIDs(4)<<endl; break;
                }
            }
        }
            break;
        }
    }

    //! -----------------------------
    //! write the number of emitters
    //! -----------------------------
    os<<"NUMBEROFEMITTERS "<<NbEmitters<<endl;

    //! --------------------------------------------
    //! if emitter faces have been found write them
    //! --------------------------------------------
    if(NbEmitters!=0)
    {
        //! --------------------------------------------------------------
        //! rescan the simulation root (start after "Analysis setting"
        //! and finish an item before "Solution") searching for particles
        //! emitters
        //! --------------------------------------------------------------
        for(int n=1; n<simulationRoot->rowCount()-1; n++)
        {
            QStandardItem *item = simulationRoot->child(n,0);
            SimulationNodeClass *curNode = item->data(Qt::UserRole).value<SimulationNodeClass*>();
            QStandardItem* itemEmitter = curNode->getPropertyItem("Emitter");
            if(itemEmitter==Q_NULLPTR) continue;
            bool emitterOn = curNode->getPropertyValue<bool>("Emitter");
            if(emitterOn==false) continue;

            os<<"EMITTER"<<endl;
            os<<"LOCATION"<<endl;

            //! ----------------------------
            //! build a composite face mesh
            //! ----------------------------
            const QVector<GeometryTag> &vecLoc = curNode->getPropertyValue<QVector<GeometryTag>>("Tags");
            QList<occHandle(Ng_MeshVS_DataSourceFace)> listOfFaceMeshDS;
            for(int n=0; n<vecLoc.length(); n++)
            {
                int faceNr = vecLoc[n].subTopNr;
                occHandle(Ng_MeshVS_DataSourceFace) aFaceMeshDS = occHandle(Ng_MeshVS_DataSourceFace)::DownCast(sDB->ArrayOfMeshDSOnFaces.getValue(bodyIndex,faceNr));
                listOfFaceMeshDS.append(aFaceMeshDS);
            }
            occHandle(Ng_MeshVS_DataSourceFace) compositeFaceDS = new Ng_MeshVS_DataSourceFace(listOfFaceMeshDS);

            //! ------------------------------
            //! write the composite face mesh
            //! ------------------------------
            int NbFaceElements = compositeFaceDS->GetAllElements().Extent();
            os<<"NUMBEROFELEMENTS "<<NbFaceElements<<endl;
            for(TColStd_MapIteratorOfPackedMapOfInteger it(compositeFaceDS->GetAllElements()); it.More(); it.Next())
            {
                int NbNodes, buf[8];
                TColStd_Array1OfInteger nodeIDs(*buf,1,8);
                compositeFaceDS->GetNodesByElement(it.Key(),nodeIDs,NbNodes);
                switch(NbNodes)
                {
                case 3: os<<"TRIG"<<endl; os<<nodeIDs(1)<<"\t"<<nodeIDs(2)<<"\t"<<nodeIDs(3)<<endl; break;
                case 4: os<<"QUAD"<<endl; os<<nodeIDs(1)<<"\t"<<nodeIDs(2)<<"\t"<<nodeIDs(3)<<"\t"<<nodeIDs(4)<<endl; break;
                }
            }

            double intensity = curNode->getPropertyValue<double>("Intensity");
            os<<"INTENSITY"<<endl;
            os<<"CONSTANT"<<endl;
            os<<intensity<<endl;

            double mass = curNode->getPropertyValue<double>("Particle mass");
            os<<"MASS"<<endl;
            os<<"CONSTANT"<<endl;
            os<<mass<<endl;

            double charge = curNode->getPropertyValue<double>("Electric charge");
            os<<"CHARGE"<<endl;
            os<<"CONSTANT"<<endl;
            os<<charge<<endl;

            double radius = curNode->getPropertyValue<double>("Radius");
            os<<"RADIUS"<<endl;
            os<<"CONSTANT"<<endl;
            os<<radius<<endl;
        }
    }

/*
    //! ----------------
    //! particles packs
    //! ----------------
    int NbParticlesPacks = 0;
    for(int n=1; n<simulationRoot->rowCount()-1; n++)
    {
        QStandardItem *item = simulationRoot->child(n,0);
        SimulationNodeClass *curNode = item->data(Qt::UserRole).value<SimulationNodeClass*>();
        if(curNode==SimulationNodeClass::nodeType_particlesPack) NbParticlesPacks++;
    }

    os<<"NUMBEROFPARTICLESPACKS "<<NbParticlesPacks<<endl;

  */

    os.close();
    return true;
}

//! ------------------------
//! function: readInputFile
//! details:
//! ------------------------
bool particlesInFieldsSolver::readInputFile(const std::string &inputFilePath,
                                            occHandle(Ng_MeshVS_DataSource3D) &volumeMesh,
                                            std::map<int,occHandle(Ng_MeshVS_DataSourceFace)> &allFacesMeshDSMap,
                                            std::vector<particlesEmitter> &theEmitters)
{
    //! --------------
    //! open a stream
    //! --------------
    ifstream is(inputFilePath);
    if(is.is_open()==false) return false;

    std::string val;
    char tmp[512];

    //! --------------------
    //! read the body index
    //! --------------------
    int bodyIndex;
    std::getline(is,val);
    sscanf(val.c_str(),"%s%d",tmp,&bodyIndex);
    cout<<val<<endl;

    //! ---------------------
    //! read the volume mesh
    //! ---------------------
    std::getline(is,val);
    cout<<val<<endl;

    volumeMesh = new Ng_MeshVS_DataSource3D(is);

    //! -------------
    //! blank line ?
    //! -------------
    std::getline(is,val);
    cout<<val<<endl;

    //! -------------------------
    //! read the number of faces
    //! -------------------------
    int NbFaces;
    std::getline(is,val);
    sscanf(val.c_str(),"%s%d",tmp,&NbFaces);
    cout<<val<<endl;

    for(int n = 1; n<=NbFaces; n++)
    {
        //! ---------------------
        //! read the face number
        //! ---------------------
        int faceNr;
        std::getline(is,val);
        sscanf(val.c_str(),"%s%d",tmp,&faceNr);
        cout<<val<<endl;

        //! ------------------------------------
        //! read the number of surface elements
        //! ------------------------------------
        int NbElements;
        std::getline(is,val);
        sscanf(val.c_str(),"%s%d",tmp,&NbElements);
        cout<<val<<endl;

        std::vector<meshElementByCoords> vecElem;
        for(int n=1; n<=NbElements; n++)
        {
            std::getline(is,val);
            cout<<val<<endl;
            if(strcmp(val.c_str(),"TRIG")==0)
            {
                int V[3];
                std::getline(is,val);
                sscanf(val.c_str(),"%d%d%d",&V[0],&V[1],&V[2]);
                cout<<val<<endl;

                //! -----------------------------
                //! build a surface mesh element
                //! -----------------------------
                meshElementByCoords anElement;
                anElement.ID = n;

                for(int n=0; n<3; n++)
                {
                    double c[3];
                    int NbNodes;
                    TColStd_Array1OfReal coords(*c,1,3);
                    MeshVS_EntityType aType;
                    volumeMesh->GetGeom(V[n],false,coords,NbNodes,aType);
                    anElement.pointList<<mesh::meshPoint(coords(1),coords(2),coords(3),V[n]);
                }
                vecElem.push_back(anElement);
            }
        }

        //! -----------------------------------
        //! build the surface mesh data source
        //! and insert it into the map
        //! -----------------------------------
        occHandle(Ng_MeshVS_DataSourceFace) aFaceMeshDS = new Ng_MeshVS_DataSourceFace(vecElem,false,false);
        std::pair<int,occHandle(Ng_MeshVS_DataSourceFace)> apair;
        apair.first = faceNr;
        apair.second = aFaceMeshDS;
        allFacesMeshDSMap.insert(apair);
    }

    //! ---------------------------------------
    //! read the number of boundary conditions
    //! ---------------------------------------
    int NbBC;
    std::getline(is,val);
    sscanf(val.c_str(),"%s%d",tmp,&NbBC);
    cout<<val<<endl;

    //! ----------------------------------------
    //! read the boundary conditions one by one
    //! ----------------------------------------
    for(int i=0; i<NbBC; i++)
    {
        //! --------------
        //! read the type
        //! --------------
        std::string BCType;
        std::getline(is,BCType);
        cout<<BCType<<endl;

        //! -----------------------
        //! read the time behavior
        //! -----------------------
        std::string value;
        if(BCType=="CONSTANT")
        {
            std::getline(is,value);
            cout<<value<<endl;
        }
        if(BCType=="FUNCTIONOFTIME")
        {
            std::getline(is,value);
            cout<<value<<endl;
        }

        //! ----------------
        //! a blank line ??
        //! ----------------
        std::getline(is,val);
        cout<<val<<endl;

        //! -------------------------------
        //! read "BOUNDARYPATCHDEFINITION"
        //! -------------------------------
        std::getline(is,val);
        cout<<val<<endl;

        //! ----------------
        //! a blank line ??
        //! ----------------
        std::getline(is,val);
        cout<<val<<endl;

        //! ---------------------------------
        //! read the number of face elements
        //! ---------------------------------
        int NbElements;
        std::getline(is,val);
        sscanf(val.c_str(),"%s%d",tmp,&NbElements);
        cout<<val<<endl;

        std::vector<meshElementByCoords> vecElem;
        for(int n=1; n<=NbElements; n++)
        {
            std::getline(is,val);
            cout<<val<<endl;
            if(strcmp(val.c_str(),"TRIG")==0)
            {
                int V[3];
                std::getline(is,val);
                sscanf(val.c_str(),"%d%d%d",&V[0],&V[1],&V[2]);
                cout<<val<<endl;

                //! -----------------------------
                //! build a surface mesh element
                //! -----------------------------
                meshElementByCoords anElement;
                anElement.ID = n;

                for(int n=0; n<3; n++)
                {
                    double c[3];
                    int NbNodes;
                    TColStd_Array1OfReal coords(*c,1,3);
                    MeshVS_EntityType aType;
                    volumeMesh->GetGeom(V[n],false,coords,NbNodes,aType);
                    anElement.pointList<<mesh::meshPoint(coords(1),coords(2),coords(3),V[n]);
                }
                vecElem.push_back(anElement);
            }
        }

        //! -----------------------------------
        //! build the surface mesh data source
        //! and insert it into the map
        //! -----------------------------------
        occHandle(Ng_MeshVS_DataSourceFace) aCompositeFaceDS = new Ng_MeshVS_DataSourceFace(vecElem,false,false);

        //! ---------------------
        //! for testing purposes
        //! ---------------------
        aCompositeFaceDS->writeMesh("D:/Work/WBtest/BCMesh.txt",2);
    }

    //! ----------------------------
    //! read the number of emitters
    //! ----------------------------
    int NbEmitters;
    std::getline(is,val);
    sscanf(val.c_str(),"%s%d",tmp,&NbEmitters);
    cout<<val<<endl;

    //! ------------------
    //! read the emitters
    //! ------------------
    for(int i=0; i<NbEmitters; i++)
    {
        //! --------------------------
        //! read the header "EMITTER"
        //! --------------------------
        std::getline(is,val);
        cout<<val<<endl;

        //! ---------------------------
        //! read the header "LOCATION"
        //! ---------------------------
        std::getline(is,val);
        cout<<val<<endl;

        //! ------------------------------------------------------
        //! read the composite face mesh
        //! read the number of surface elements, the elements one
        //! by one, finally build the composite mesh data source
        //! ------------------------------------------------------
        int NbElements;
        std::getline(is,val);
        sscanf(val.c_str(),"%s%d",tmp,&NbElements);
        cout<<val<<endl;

        std::vector<meshElementByCoords> vecElem;
        for(int n=1; n<=NbElements; n++)
        {
            cout<<n<<endl;
            std::getline(is,val);
            cout<<val<<endl;
            if(strcmp(val.c_str(),"TRIG")==0)
            {
                int V[3];
                std::getline(is,val);
                sscanf(val.c_str(),"%d%d%d",&V[0],&V[1],&V[2]);
                cout<<val<<endl;

                //! -----------------------------
                //! build a surface mesh element
                //! -----------------------------
                meshElementByCoords anElement;
                anElement.ID = n;

                for(int n=0; n<3; n++)
                {
                    double c[3];
                    int NbNodes;
                    TColStd_Array1OfReal coords(*c,1,3);
                    MeshVS_EntityType aType;
                    volumeMesh->GetGeom(V[n],false,coords,NbNodes,aType);
                    anElement.pointList<<mesh::meshPoint(coords(1),coords(2),coords(3),V[n]);
                }
                vecElem.push_back(anElement);
            }
        }

        //! -----------------------------------
        //! build the surface mesh data source
        //! and insert it into the map
        //! -----------------------------------
        occHandle(Ng_MeshVS_DataSourceFace) aCompositeFaceDS = new Ng_MeshVS_DataSourceFace(vecElem,false,false);

        //! ---------------------
        //! for testing purposes
        //! ---------------------
        aCompositeFaceDS->writeMesh(QString("D:/Work/WBtest/EmitterMesh_%1.txt").arg(i+1),2);

        //! read the header "INTENSITY"
        std::getline(is,val);
        cout<<val<<endl;

        //! read the header "CONSTANT" or "FUNCTION"
        std::getline(is,val);
        cout<<val<<endl;

        double intensity = 0;
        std::string timeDependentIntensity;
        if(val=="CONSTANT")
        {
            //! -------------------
            //! read the intensity
            //! -------------------
            std::getline(is,val);
            sscanf(val.c_str(),"%lf",&intensity);
            cout<<val<<endl;
        }
        if(val=="FUNCTION")
        {
            cerr<<"particlesInFieldsSolver::readInputFile()->____intensity defined by function not supported____"<<endl;
            exit(999999);
        }

        //! read the header "MASS"
        std::getline(is,val);
        cout<<val<<endl;

        //! read the header "CONSTANT" or "FUNCTION"
        std::getline(is,val);
        cout<<val<<endl;

        double mass;
        std::string massDistribution;
        if(val=="CONSTANT")
        {
            //! -------------------
            //! read the intensity
            //! -------------------
            std::getline(is,val);
            sscanf(val.c_str(),"%lf",&mass);
            cout<<val<<endl;
        }
        if(val=="FUNCTION")
        {
            cerr<<"particlesInFieldsSolver::readInputFile()->____mass defined by function not supported____"<<endl;
            exit(999999);
        }

        //! read the header "CHARGE"
        std::getline(is,val);
        cout<<val<<endl;

        //! read the header "CONSTANT" or "FUNCTION"
        std::getline(is,val);
        cout<<val<<endl;

        double charge;
        std::string chargeValueDistribution;
        if(val=="CONSTANT")
        {
            //! -------------------
            //! read the intensity
            //! -------------------
            std::getline(is,val);
            sscanf(val.c_str(),"%lf",&charge);
            cout<<val<<endl;
        }
        if(val=="FUNCTION")
        {
            cerr<<"particlesInFieldsSolver::readInputFile()->____charge distribution defined by function not supported____"<<endl;
            exit(999999);
        }

        //! read the header "RADIUS"
        std::getline(is,val);
        cout<<val<<endl;

        //! read the header "CONSTANT" or "FUNCTION"
        std::getline(is,val);
        cout<<val<<endl;

        double radius;
        std::string radiusDistribution;
        if(val=="CONSTANT")
        {
            //! -------------------
            //! read the intensity
            //! -------------------
            std::getline(is,val);
            sscanf(val.c_str(),"%lf",&radius);
            cout<<val<<endl;
        }
        if(val=="FUNCTION")
        {
            cerr<<"particlesInFieldsSolver::readInputFile()->____radius distribution defined by function not supported____"<<endl;
            exit(999999);
        }

        //! -----------------
        //! build an emitter
        //! -----------------
        particlesEmitter aParticleEmitter;
        aParticleEmitter.setLocation(aCompositeFaceDS);
        aParticleEmitter.setIntensity(intensity);
        aParticleEmitter.setChange(charge);
        aParticleEmitter.setMass(mass);
        aParticleEmitter.setRadius(radius);

        //! ------------------------------
        //! add to the vector of emitters
        //! ------------------------------
        theEmitters.push_back(aParticleEmitter);
    }

    return true;
}


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
//! details:  from file
//! ----------------------
particlesInFieldsSolver::particlesInFieldsSolver(const string &inputFilePath)
{
    ifstream is(inputFilePath);
    if(is.is_open()==false) return;

    //! ---------------------------------------------
    //! read the mesh, the face meshes, the emitters
    //! ---------------------------------------------
    occHandle(Ng_MeshVS_DataSource3D) volumeMesh;
    std::map<int,occHandle(Ng_MeshVS_DataSourceFace)> allFacesMeshDS;
    std::vector<particlesEmitter> vecEmitters;
    this->readInputFile(inputFilePath,volumeMesh,allFacesMeshDS,vecEmitters);

    //! -----------------------------------------------------------------
    //! - configure the domain mesh in Eigen format
    //! - initialize AABB tree for point location
    //! - allocate space for density distribution and reset distribution
    //! - allocate the space for the electric field
    //! - init the number of particles to zero
    //! - mass matrix and its inverse
    //! - gradient operator
    //! -----------------------------------------------------------------
    this->setUpDomainMesh(volumeMesh);

    //! --------------------------
    //! create the Poisson solver
    //! --------------------------
    myPoissonSolver = std::make_shared<PoissonSolver>(volumeMesh);
    myPoissonSolver->definePatches(allFacesMeshDS);
    myPoissonSolver->initPotentialOnBoundary(0.0);

    //! --------------------
    //! set up the emitters
    //! --------------------
    for(int i=0; i<vecEmitters.size(); i++) myEmitters.push_back(vecEmitters[i]);


    //! -----------------------------
    //! init the number of particles
    //! -----------------------------
    myNbParticles = 0;
}

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
particlesInFieldsSolver::particlesInFieldsSolver(simulationDataBase *sDB, QStandardItem *simulationRoot, QObject *parent): QObject(parent)
{
    mySimulationRoot = simulationRoot;
    mySDB =  sDB;

    //! -----------------------
    //! set up the domain mesh
    //! -----------------------
    occHandle(Ng_MeshVS_DataSource3D) volumeMesh = occHandle(Ng_MeshVS_DataSource3D)::DownCast(sDB->ArrayOfMeshDS.first());

    //! -----------------------------------------------------------------
    //! - configure the domain mesh in Eigen format
    //! - initialize AABB tree for point location
    //! - allocate space for density distribution and reset distribution
    //! - allocate the space for the electric field
    //! - init the number of particles
    //! - mass matrix and its inverse
    //! - gradient operator
    //! -----------------------------------------------------------------
    this->setUpDomainMesh(volumeMesh);

    //! --------------------------
    //! create the Poisson solver
    //! --------------------------
    myPoissonSolver = std::make_shared<PoissonSolver>(volumeMesh);
    std::map<int,occHandle(Ng_MeshVS_DataSourceFace)> allFacesMap;
    int bodyIndex = -1;
    int NbFaces = sDB->MapOfBodyTopologyMap.value(bodyIndex).faceMap.Extent();
    for(int faceNr = 1; faceNr<NbFaces; faceNr++)
    {
        const occHandle(Ng_MeshVS_DataSourceFace) &aFaceMeshDS = occHandle(Ng_MeshVS_DataSourceFace)::DownCast(sDB->ArrayOfMeshDSOnFaces.getValue(bodyIndex,faceNr));
        std::pair<int,occHandle(Ng_MeshVS_DataSourceFace)> apair;
        apair.first = faceNr;
        apair.second = aFaceMeshDS;
        allFacesMap.insert(apair);
    }
    myPoissonSolver->definePatches(allFacesMap);
    myPoissonSolver->initPotentialOnBoundary(0.0);

    //! -----------------------------
    //! init the number of particles
    //! -----------------------------
    myNbParticles = 0;

    //! --------------------------------------------------------------
    //! simulation time - time step from the "Analysis settings" item
    //! --------------------------------------------------------------
    SimulationNodeClass *nodeAnalysisSetting = mySimulationRoot->child(0,0)->data(Qt::UserRole).value<SimulationNodeClass*>();

    CustomTableModel *tabData = nodeAnalysisSetting->getTabularDataModel();
    int lastTableRow = tabData->rowCount();
    int col = TABULAR_DATA_STEP_END_TIME_COLUMN;

    myFinalTime = tabData->dataRC(lastTableRow,col).toDouble();
    myTimeStep = nodeAnalysisSetting->getPropertyValue<double>("Time step size");
}

//! ---------------
//! function: init
//! details:
//! ---------------
bool particlesInFieldsSolver::init(QStandardItem* simulationRoot)
{
    //! ---------------------------------
    //! model faces on which a BC is put
    //! ---------------------------------
    std::map<int,occHandle(Ng_MeshVS_DataSourceFace)> BoundaryConditionsFaceDSMap;

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
            const QVector<GeometryTag> &vecLoc = curNode->getPropertyValue<QVector<GeometryTag>>("Tags");
            int NbPatches = vecLoc.length();
            QList<occHandle(Ng_MeshVS_DataSourceFace)> listOfFaceMeshDS;
            for(int n=0; n<NbPatches; n++)
            {
                std::pair<int,occHandle(Ng_MeshVS_DataSourceFace)> apair;
                int faceNr = vecLoc[0].subTopNr;
                const occHandle(Ng_MeshVS_DataSourceFace) &aFaceMeshDS =
                        occHandle(Ng_MeshVS_DataSourceFace)::DownCast(mySDB->ArrayOfMeshDSOnFaces.getValue(1,faceNr));
                apair.first = faceNr;
                apair.second = aFaceMeshDS;
                BoundaryConditionsFaceDSMap.insert(apair);

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
                double radius = curNode->getPropertyValue<double>("Radius");

                occHandle(Ng_MeshVS_DataSourceFace) emitterFaceMeshDS = new Ng_MeshVS_DataSourceFace(listOfFaceMeshDS);
                particlesEmitter anEmitter(emitterFaceMeshDS,intensity);
                myEmitters.push_back(anEmitter);
            }
        }
            break;
        }

        //! ------------------------------------------------
        //! define the patches and initialize the potential
        //! on all the surface mesh nodes
        //! ------------------------------------------------
        myPoissonSolver->definePatches(BoundaryConditionsFaceDSMap);

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
            particlesEmitter anEmitter = myEmitters.at(i);
            std::vector<particle> newParticles;
            anEmitter.generateParticles(newParticles,time);
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

//! -----------------------
//! function: loadTimeInfo
//! details:
//! -----------------------
bool particlesInFieldsSolver::loadTimeInfo(ifstream &inputFileStream)
{
    std::string val;

    //! header
    std::getline(inputFileStream,val);
    if(strcmp(val.c_str(),"Final time")!=0) return false;

    //! final simulation time
    inputFileStream>>myFinalTime;

    //! header
    std::getline(inputFileStream,val);
    if(strcmp(val.c_str(),"Time step")!=0) return false;

    //! time step size
    inputFileStream>>myTimeStep;

    return true;
}

//! --------------------------
//! function: setUpDomainMesh
//! details:
//! --------------------------
bool particlesInFieldsSolver::setUpDomainMesh(occHandle(Ng_MeshVS_DataSource3D) volumeMesh)
{
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

    return true;
}

//! -----------------------
//! function: loadParticle
//! details:
//! -----------------------
int particlesInFieldsSolver::loadParticles(ifstream &inputFileStream)
{
    int N;
    inputFileStream>>N;
    for(int n=1; n<=N; n++)
    {
        double x,y,z,vx,vy,vz,m,q,r;
        std::string val;
        std::getline(inputFileStream,val);
        sscanf(val.c_str(),"%lf%lf%lf%lf%lf%lf%lf%lf%lf",&x,&y,&z,&vx,&vy,&vz,&m,&q,&r);
        particle aParticle(x,y,z,vx,vy,m,q,r);
        std::pair<int,particle> apair;
        apair.first = 1;
        apair.second = aParticle;
        myParticles1->insert(apair);
    }
    myNbParticles++;
    return N;
}
