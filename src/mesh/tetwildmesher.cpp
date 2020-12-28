#define TETWILD_PATH "D:/Work/Qt/build_simSpace/release/FloatTetwild_bin.exe"

//! -------
//! pyMesh
//! -------
#include <MshSaver.h>

//! ---------------
//! tetgen library
//! ---------------
#include <tetgen/inc/tetgen.h>

//! ----------------
//! custom includes
//! ----------------
#include <mshconvert.h>
#include <tetwildmesher.h>
#include <meshdatabase.h>

//! ---
//! Qt
//! ---
#include <QProcess>
#include <QString>
#include <QDir>

//! ----
//! C++
//! ----
#include <direct.h>
#include <random>
#include <iostream>
using namespace std;

//! ----
//! OCC
//! ----
#include <Bnd_Box.hxx>
#include <BRepBndLib.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepTools.hxx>
#include <BRep_Tool.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <CPnts_AbscissaPoint.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Vertex.hxx>

//! --------
//! Windows
//! --------
#include <sys/stat.h>
#include <windows.h>

#include <stdio.h>
#include <dirent.h>
#define SUCCESS_STAT 0

//! --------------------------------------
//! function: dirExists
//! details:  check if a directory exists
//! --------------------------------------
bool dirExists(std::string dir_path)
{
    struct stat sb;

    if (stat(dir_path.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode))
        return true;
    else
        return false;
}

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
tetWildMesher::tetWildMesher(QObject *parent):QObject(parent),
    myWillRunOnDisk(true),
    myAbsoluteOutputFilePath("")
{
    myTetWildProcess = new QProcess(this);
    this->createSupportFilesDirectory();
}

//! ---------------------
//! function: destructor
//! details:
//! ---------------------
tetWildMesher::~tetWildMesher()
{
    Sleep(2000);
    cout<<"tetWildMesher::~tetWildMesher()->____DESTRUCTOR CALLED____"<<endl;
    this->removeSupportFilesDirectory();
}

//! ---------------------------------------------------------------------
//! function: setMeshDataBase
//! details:  mesh sizing controls are contained into the mesh data base
//! ---------------------------------------------------------------------
void tetWildMesher::setDataBase(meshDataBase *mDB)
{
    myMeshDB = mDB;
}

//! -------------------
//! function: setInput
//! details:
//! -------------------
void tetWildMesher::setInput(const Eigen::MatrixXd &VI, const Eigen::MatrixXi &FI)
{
    myVertices = VI;
    myFaces = FI;
}

//! -------------------------------------------------
//! function: setParameters
//! details:  meshParam contains default parameters.
//!           Using this function they are by-passed
//! -------------------------------------------------
void tetWildMesher::setGlobalParameters(float l_rel, float eps_rel)
{
    myMeshParam.eps_rel = eps_rel;
    myMeshParam.initial_edge_len_rel = l_rel;
}

//! ------------------------------------------------------
//! function: performVolume
//! details:  VO - output vertices, TO - ouput tetrahedra
//! ------------------------------------------------------
bool tetWildMesher::perform_inMemory(Eigen::MatrixXd &VO, Eigen::MatrixXi &TO)
{
    Q_UNUSED(VO)
    Q_UNUSED(TO)
    //Eigen::VectorXd AO;
    //tetwild::tetrahedralization(myVertices,myFaces,VO,TO,AO,myMeshParam);
    return false;
}

//! ------------------------------------------------------
//! function: performVolume_onDisk
//! details:  VO - output vertices, TO - ouput tetrahedra
//! ------------------------------------------------------
bool tetWildMesher::perform_onDisk(const QString &absoluteInputFilePath)
{
    cout<<"\n@____Experimental mesher will run on disk____@"<<endl;

    //! ------------------
    //! preliminary check
    //! ------------------
    QFile file;
    file.setFileName(absoluteInputFilePath);
    if(!file.exists()) return false;
    if(absoluteInputFilePath.isEmpty()) return false;

    QString tmp(absoluteInputFilePath);
    QString relativeInputFileName = tmp.split("/").last();
    tmp.chop(relativeInputFileName.length());
    QString workingDir(tmp);
    QString extensionFreeRelativeInputFileName = relativeInputFileName;
    extensionFreeRelativeInputFileName.chop(4);
    QString absoluteOuputFilePath = workingDir+extensionFreeRelativeInputFileName+"_tetWild.msh";

    myAbsoluteOutputFilePath = absoluteOuputFilePath;

    //! -----------------------------------------------------------------------
    //! search for fTetWild_bin.exe
    //! as a general rule the file should placed into the executable directory
    //! -----------------------------------------------------------------------
    /*
    std::string exePath = this->getPathOfExecutable()+"\\FloatTetwild_bin.exe";
    cout<<"@____fTetWildPath: "<<exePath<<"____"<<endl;
    struct stat buf;
    int fileExist = stat(exePath.c_str(),&buf);
    if(fileExist!=0) return false;
    myTetWildProcess->setProgram(QString::fromStdString(exePath));
    */
    myTetWildProcess->setProgram(TETWILD_PATH);

    //! --------------
    //! set arguments
    //! --------------
    QList<QString> arguments;
    arguments<<"-l"<<QString("%1").arg(myMeshParam.initial_edge_len_rel);
    arguments<<"-e"<<QString("%1").arg(myMeshParam.eps_rel);
    arguments<<"--input"<<absoluteInputFilePath;
    arguments<<"--output"<<absoluteOuputFilePath;
    myTetWildProcess->setArguments(arguments);

    //! ------------------------------
    //! diagnostic - console messages
    //! ------------------------------
    cout<<"@____running fTetWild mesher with arguments: ";
    for(int i=0; i<arguments.length(); i++) cout<<arguments.at(i).toStdString()<<" ";
    cout<<endl;

    disconnect(myTetWildProcess,SIGNAL(readyReadStandardOutput()),this,SLOT(readTetWildProcess()));
    connect(myTetWildProcess,SIGNAL(readyReadStandardOutput()),this,SLOT(readTetWildProcess()));

    //! ---------------
    //! start fTetWild
    //! ---------------
    myTetWildProcess->start();
    int exitCode = myTetWildProcess->exitCode();
    myTetWildProcess->waitForFinished(-1);

    cout<<"@____\"Experimental mesher\" process exit code: "<<exitCode<<endl<<endl;

    if(exitCode==0) return true;
    return false;
}

//! -----------------------------------------
//! function: retrieveMeshDataSourceFromDisk
//! details:
//! -----------------------------------------
bool tetWildMesher::retrieveMeshDataSourceFromDisk(opencascade::handle<Ng_MeshVS_DataSource3D> &mesh3D_DataSource)
{
    cout<<"tetWildMesher::retrieveMeshDataSourceFromDisk()->____function called____"<<endl;
    if(myAbsoluteOutputFilePath.isEmpty()) return false;

    cout<<"tetWildMesher::retrieveMeshDataSourceFromDisk()->____calling .msh converter____"<<endl;
    mshConvert mshConverter(myAbsoluteOutputFilePath);
    Eigen::MatrixXd V;
    Eigen::MatrixXi T;

    bool isDone = mshConverter.toEigen(V,T);
    if(!isDone) return false;

    cout<<"tetWildMesher::retrieveMeshDataSourceFromDisk()->____build the 3D mesh data source____"<<endl;
    mesh3D_DataSource = new Ng_MeshVS_DataSource3D(V,T);
    if(mesh3D_DataSource.IsNull()) return false;
    return true;
}

//! -------------------------
//! function: extractSurface
//! details:
//! -------------------------
bool tetWildMesher::extractSurface(const Eigen::MatrixXd &VI, const Eigen::MatrixXi &TI,
                                   Eigen::MatrixXd &VS, Eigen::MatrixXi &FS)
{
    Q_UNUSED(VI)
    Q_UNUSED(TI)
    Q_UNUSED(VS)
    Q_UNUSED(FS)
    return false;
}

//! -------------------------------------
//! function: readTetWildProcess
//! details:  send message to std output
//! -------------------------------------
void tetWildMesher::readTetWildProcess()
{
    QByteArray msg = myTetWildProcess->readAllStandardOutput();
    cout<<msg.toStdString()<<endl;
}

//! ------------------------------
//! function: getPathOfExecutable
//! details:  helper
//! ------------------------------
std::string tetWildMesher::getPathOfExecutable()
{
    LPWSTR buffer;
    GetModuleFileName(NULL, buffer, MAX_PATH);
    std::string aString;
    aString.reserve(wcslen(buffer));
    for (;*buffer; buffer++) aString += (char)*buffer;
    std::string::size_type pos = aString.find_last_of("\\/");
    if (pos == std::string::npos) return "";
    return aString.substr(0, pos);
}

//! --------------------------------------------------------------------
//! function: sampleGeometry
//! details:  a set of points is sampled on a sub-geometry and returned
//! --------------------------------------------------------------------
void tetWildMesher::sampleGeometry(const TopoDS_Shape &aShape,
                                   void *parametersForSampling,
                                   std::vector<tetWildMesher::point> &sampledPoints)
{
    cout<<"sampleGeometry::sampleGeometry()->____function called____"<<endl;
    if(aShape.IsNull())
    {
        cout<<"sampleGeometry::sampleGeometry()->____the input shape is null____"<<endl;
        return;
    }

    //! ------------------------------
    //! the bounding box of the shape
    //! ------------------------------
    Bnd_Box BB;
    BRepBndLib::Add(aShape, BB);
    double xmin, ymin, zmin, xmax, ymax, zmax;
    BB.Get(xmin,ymin,zmin,xmax,ymax,zmax);

    switch(aShape.ShapeType())
    {
    case TopAbs_FACE:
    {
        //! ------------------------------------------------------------
        //! one parameter in "samplingParameters": face sizing at point
        //! ------------------------------------------------------------
        std::vector<double>* samplingParameters =(std::vector<double>*)(parametersForSampling);
        double elementSizeOnFace = samplingParameters->at(0);

        //! -----------------------------------------------
        //! a new bounding box for the face - a bit larger
        //! -----------------------------------------------
        double dx = (xmax-xmin)*0.025;
        double dy = (ymax-ymin)*0.025;
        double dz = (zmax-zmin)*0.025;

        double xmin_bb = xmin - dx; double xmax_bb = xmax + dx;
        double ymin_bb = ymin - dy; double ymax_bb = ymax + dy;
        double zmin_bb = zmin - dz; double zmax_bb = zmax + dz;

        //! ---------
        //! the face
        //! ---------
        const TopoDS_Face &aFace = TopoDS::Face(aShape);

        //! --------------------
        //! surface of the face
        //! --------------------
        BRepAdaptor_Surface adaptor(aFace,true);
        const GeomAdaptor_Surface &s_adaptor = adaptor.Surface();
        occHandle(Geom_Surface) CurGeomSurface = s_adaptor.Surface();

        //! --------------------------------------------------
        //! retrieve the u,v bounds
        //! enclose the into a rectangle and perform sampling
        //! for the moment a 10x10 grid is used ...
        //! --------------------------------------------------
        double umin,umax,vmin,vmax;
        BRepTools::UVBounds(aFace,umin,umax,vmin,vmax);
        double Lu = umax-umin;
        double Lv = vmax-vmin;
        double du = Lu/10;
        double dv = Lv/10;

        for(double ucur = umin; ucur<=umax; ucur += du)
        {
            for(double vcur = vmin; vcur<=vmax; vcur +=dv)
            {
                const gp_Pnt &curPoint = CurGeomSurface->Value(ucur,vcur);
                double x = curPoint.X();
                double y = curPoint.Y();
                double z = curPoint.Z();

                //! --------------------------------------------------
                //! check if the point is whithin the BB of the shape
                //! --------------------------------------------------
                if((x>=xmin_bb && x<=xmax_bb) && (y>=ymin_bb && y<=ymax_bb) && (z>=zmin_bb && z<=zmax_bb))
                {
                    point aPoint(x,y,z,elementSizeOnFace);
                    sampledPoints.push_back(aPoint);
                }
            }
        }
    }
        break;

    case TopAbs_EDGE:
    {
        //! ------------------------------------------------------------------------------------------------------
        //! parameter for the edge sizing: two parameters in the vector
        //! first parameter: type of sizing:
        //!                 "0" => through size of the element => the second parameter is the element size
        //!                 "1" => through number of divisions => the second parameter is the number of divisions
        //! ------------------------------------------------------------------------------------------------------
        std::vector<double> *samplingParameters =  (std::vector<double>*)(parametersForSampling);

        double typeOfSizing = samplingParameters->at(0);
        double secondParameter = samplingParameters->at(1);

        //! ---------
        //! the edge
        //! ---------
        const TopoDS_Edge &anEdge = TopoDS::Edge(aShape);
        if(BRep_Tool::Degenerated(anEdge)) break;
        BRepAdaptor_Curve BRepAdaptor(anEdge);
        GeomAdaptor_Curve curve = BRepAdaptor.Curve();

        //! -------------------------------
        //! compute the length of the edge
        //! -------------------------------
        CPnts_AbscissaPoint CP;
        CP.Init(curve);
        double s_old = BRepAdaptor.FirstParameter();
        const gp_Pnt &P_onEdge = BRepAdaptor.Value(s_old);
        double x = P_onEdge.X();
        double y = P_onEdge.Y();
        double z = P_onEdge.Z();

        double Lcurrent = 0;
        double L = CP.Length(curve,1e-2);
        double size = 0;
        if(typeOfSizing == 0) size = secondParameter;
        else if(typeOfSizing == 1) size = L/secondParameter;

        sampledPoints.push_back(point(x,y,z,size)); // first point of the current edge
        for(double s=s_old-size;Lcurrent<=L+size;)
        {
            CP.Perform(size,s,1e-2);
            const gp_Pnt &P1_onEdge = BRepAdaptor.Value(s);
            x = P1_onEdge.X();
            y = P1_onEdge.Y();
            z = P1_onEdge.Z();
            sampledPoints.push_back(point(x,y,z,size)); // other points of the edge
            s_old=s;
            s=CP.Parameter();
            Lcurrent += fabs(s-s_old);
        }
    }
        break;

    case TopAbs_VERTEX:
    {
        //! --------------------------------------------------
        //! parameter for the mesh sizing at the vertex
        //! only one parameters in the vector: pinball radius
        //! --------------------------------------------------
        std::vector<double> *samplingParameters = (std::vector<double>*)(parametersForSampling);

        //! radial, theta, phi divisions
        const double PI = 3.141592654;
        const int radialDiv = 10;
        const int thetaDiv = 20;
        const int phiDiv = 10;

        //! pinball
        double R = samplingParameters->at(0);
        double dr = R/radialDiv;

        //! mesh sizing
        double meshSize = samplingParameters->at(1);

        //! coordinates of the selected vertex
        const gp_Pnt &V = BRep_Tool::Pnt(TopoDS::Vertex(aShape));

        double deltaTheta = 2*PI/thetaDiv;
        double deltaPhi = PI/phiDiv;
        for(double theta =0; theta<2*PI; theta +=deltaTheta)
        {
            for(double phi = -PI/2; phi<=PI/2; phi +=deltaPhi)
            {
                for(double r = dr; r<=R; r+=dr)
                {
                    double x = V.X()+r*cos(phi)*cos(theta);
                    double y = V.Y()+r*cos(phi)*sin(theta);
                    double z = V.Z()+r*sin(phi);
                    sampledPoints.push_back(point(x,y,z,meshSize));
                }
            }
        }
    }
        break;
    }
}

//! -------------------------------------------------------------------------
//! function: writeMeshSizingField
//! details:  the mesh sizing fiels is defined through a vector of locations
//!           having a value (the mesh size at location)
//! -------------------------------------------------------------------------
void tetWildMesher::computeMeshSizingField(int bodyIndex, std::vector<tetWildMesher::point> &meshSizingField)
{
    //! -----------------------------------------------------
    //! retrieve the mesh size parameters from the data base
    //! -----------------------------------------------------
    if(myMeshDB==nullptr)
    {
        cout<<"********************************************************************************"<<endl;
        cout<<" cannot retrieve the mesh sizing parameters: the mesh data base is null/not set "<<endl;
        cout<<"********************************************************************************"<<endl;
        return;
    }

    //! -----------------------
    //! mesh controls on faces
    //! -----------------------
    int NbFaces = myMeshDB->MapOfBodyTopologyMap.value(bodyIndex).faceMap.Extent();
    for(int faceNr = 1; faceNr<=NbFaces; faceNr++)
    {
        bool isMeshControlOnFace = myMeshDB->MapOfIsFaceModified.getValue(bodyIndex,faceNr);
        if(isMeshControlOnFace == false) continue;

        //! ------------------------------------------------------------
        //! retrieve the mesh face sizing for the current body and face
        //! from the mesh data base
        //! ------------------------------------------------------------
        std::vector<double> parametersFaces;
        double faceMeshSizing = myMeshDB->MapOfElementSizeOnFace.getValue(bodyIndex,faceNr);
        parametersFaces.push_back(faceMeshSizing);
        void *parametersForSamplingFaces = (void*)(&parametersFaces);

        TopoDS_Shape aFace = myMeshDB->MapOfBodyTopologyMap.value(bodyIndex).faceMap.FindKey(faceNr);
        if(aFace.IsNull()) continue;
        std::vector<tetWildMesher::point> meshSizingFieldOnFace;
        this->sampleGeometry(aFace,parametersForSamplingFaces,meshSizingFieldOnFace);

        //! ---------------------------
        //! fill the mesh sizing field
        //! ---------------------------
        for(std::vector<tetWildMesher::point>::iterator it = meshSizingFieldOnFace.begin(); it!= meshSizingFieldOnFace.end(); it++)
        {
            const tetWildMesher::point &aSampledPoint = *it;
            meshSizingField.push_back(aSampledPoint);
        }
    }

    //! -----------------------
    //! mesh controls on edges
    //! -----------------------
    int NbEdges = myMeshDB->MapOfBodyTopologyMap.value(bodyIndex).edgeMap.Extent();
    for(int edgeNr = 1; edgeNr<=NbEdges; edgeNr++)
    {
        bool isMeshControlOnEdge = myMeshDB->MapOfIsEdgeModified.getValue(bodyIndex,edgeNr);
        if(isMeshControlOnEdge == false) continue;
        TopoDS_Shape anEdge = myMeshDB->MapOfBodyTopologyMap.value(bodyIndex).edgeMap.FindKey(edgeNr);
        if(anEdge.IsNull()) continue;

        //! ------------------------------------------------------------
        //! retrieve the mesh face sizing for the current body and face
        //! from the mesh data base
        //! ------------------------------------------------------------
        double typeOfSizing = myMeshDB->MapOfSizingTypeOnEdge.getValue(bodyIndex,edgeNr);
        double secondParameter = myMeshDB->MapOfElementSizeOnEdge.getValue(bodyIndex,edgeNr);
        std::vector<double> parametersEdge { typeOfSizing, secondParameter };
        void *parametersForSamplingEdge = (void*)(&parametersEdge);
        std::vector<tetWildMesher::point> meshSizingFieldOnEdge;
        this->sampleGeometry(anEdge,parametersForSamplingEdge,meshSizingFieldOnEdge);

        //! ---------------------------
        //! fill the mesh sizing field
        //! ---------------------------
        for(std::vector<tetWildMesher::point>::iterator it = meshSizingFieldOnEdge.begin(); it!= meshSizingFieldOnEdge.end(); it++)
        {
            const tetWildMesher::point &aSampledPoint = *it;
            meshSizingField.push_back(aSampledPoint);
        }
    }

    //! -------------------------------------------
    //! mesh controls on points
    //! first parameter: pinball radius
    //! second parameter: mesh size within pinball
    //! -------------------------------------------
    int NbVertex = myMeshDB->MapOfBodyTopologyMap.value(bodyIndex).vertexMap.Extent();
    for(int vertexNr = 1; vertexNr<=NbVertex; vertexNr++)
    {
        bool isMeshControlOnVertex = myMeshDB->MapOfIsVertexModified.getValue(bodyIndex,vertexNr);
        if(isMeshControlOnVertex == false) continue;
        TopoDS_Shape aVertex = myMeshDB->MapOfBodyTopologyMap.value(bodyIndex).vertexMap.FindKey(vertexNr);
        if(aVertex.IsNull()) continue;

        double pinBallRadius = myMeshDB->MapOfVertexPinball.getValue(bodyIndex,vertexNr);
        double meshPointSize = myMeshDB->MapOfElementSizeOnVertex.getValue(bodyIndex,vertexNr);
        std::vector<double> parameterPoint { pinBallRadius, meshPointSize };
        void *parametersForSamplingPoint = (void*)(&parameterPoint);
        std::vector<tetWildMesher::point> meshSizingFieldOnVertex;
        this->sampleGeometry(aVertex,parametersForSamplingPoint,meshSizingFieldOnVertex);

        //! ---------------------------
        //! fill the mesh sizing field
        //! ---------------------------
        for(std::vector<tetWildMesher::point>::iterator it = meshSizingFieldOnVertex.begin(); it!= meshSizingFieldOnVertex.end(); it++)
        {
            const tetWildMesher::point &aSampledPoint = *it;
            meshSizingField.push_back(aSampledPoint);
        }
    }
}

void ReplaceStringInPlace(std::string& subject, const std::string& search,
                          const std::string& replace) {
    size_t pos = 0;
    while((pos = subject.find(search, pos)) != std::string::npos) {
         subject.replace(pos, search.length(), replace);
         pos += replace.length();
    }
}

//! -------------------------------
//! function: writeMeshSizingField
//! details:
//! -------------------------------
bool tetWildMesher::writeMeshSizingField(const std::string &backgroundMesh, const std::vector<tetWildMesher::point> &inputMeshSizingField)
{
    if(backgroundMesh.empty()) return false;
    if(inputMeshSizingField.empty()) return false;

    //! -----------------------------------------------------------------------------------------------
    //! filter the mesh sizing field: the value of the field at a duplicated location will be averaged
    //! -----------------------------------------------------------------------------------------------
    std::map<tetWildMesher::point,std::vector<double>> mapOfPointsWithValue;
    for(std::vector<tetWildMesher::point>::const_iterator it = inputMeshSizingField.cbegin(); it != inputMeshSizingField.cend(); it++)
    {
        const tetWildMesher::point &aP = *it;
        //cout<<aP.x<<"\t"<<aP.y<<"\t"<<aP.z<<"\t"<<aP.value<<endl;
        std::map<tetWildMesher::point,std::vector<double>>::iterator it1 = mapOfPointsWithValue.find(aP);
        if(it1==mapOfPointsWithValue.end())
        {
            std::vector<double> v { aP.value };
            mapOfPointsWithValue.insert(std::make_pair(aP,v));
        }
        else
        {
            it1->second.push_back(aP.value);
        }
    }

    //! --------------------------------------------------
    //! a new sizing field with slighlty perturbed points
    //! --------------------------------------------------
    std::random_device rd;      //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd());     //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(-0.1, 0.1);

    //! ----------------------------------------
    //! average the values at duplicated points
    //! ----------------------------------------
    std::vector<tetWildMesher::point> meshSizingField;
    for(std::map<tetWildMesher::point,std::vector<double>>::iterator it = mapOfPointsWithValue.begin(); it!=mapOfPointsWithValue.end(); it++)
    {
        tetWildMesher::point aPoint = it->first;
        std::vector<double> values = it->second;
        double averageValue = 0;
        int N = 0;
        for(std::vector<double>::iterator it1 = values.begin(); it1!=values.end(); it1++)
        {
            averageValue += *it1;
            N++;
        }
        //cout<<"____"<<values.size()<<"____"<<endl;
        averageValue /= N;
        aPoint.x += dis(gen);
        aPoint.y += dis(gen);
        aPoint.z += dis(gen);
        aPoint.value = averageValue;
        meshSizingField.push_back(aPoint);
    }

    int NbPoints = int(meshSizingField.size());

    //! ---------
    //! use disk
    //! ---------
    char tetgenNodeFilePath[512];
    sprintf(tetgenNodeFilePath,"%s\\%s",mySupportFilesDirectory.c_str(),"meshSizingField.node");
    cout<<"____"<<tetgenNodeFilePath<<"____"<<endl;

    FILE *fp = fopen(tetgenNodeFilePath,"w");
    fprintf(fp,"# Node count, 3 dim, no attribute, no boundary marker\n");
    fprintf(fp,"%d\t%d\t%d\t%d\n",NbPoints,3,0,0);

    for(int n=0; n<NbPoints; n++)
    {
        const tetWildMesher::point &aP = meshSizingField[n];
        fprintf(fp,"%d\t%.6lf\t%.6lf\t%.6lf\n",n+1,aP.x,aP.y,aP.z);
    }
    fclose(fp);


    //! --------------------
    //! tetgen.exe location
    //! --------------------
    char tetgenExeLocation[128];
    sprintf(tetgenExeLocation,"D:/tetgen1.5.0/build/Release/tetgen.exe");

    //! ------------------------------------
    //! execute tetgen throgh shell command
    //! ------------------------------------
    char command[256];
    sprintf(command, "%s %s", tetgenExeLocation, tetgenNodeFilePath);
    system((char *)command);

    char eleFileLocation[512];
    sprintf(eleFileLocation,"%s\\%s",mySupportFilesDirectory.c_str(),"meshSizingField1.ele");
    cout<<"____"<<eleFileLocation<<"____"<<endl;
    char nodeFileLocation[512];
    sprintf(nodeFileLocation,"%s\\%s",mySupportFilesDirectory.c_str(),"meshSizingField1.node");
    cout<<"____"<<nodeFileLocation<<"____"<<endl;

    std::string a(nodeFileLocation);
    ReplaceStringInPlace(a,std::string("\\"),std::string("/"));
    std::string b(eleFileLocation);
    ReplaceStringInPlace(b,std::string("\\"),std::string("/"));

    cout<<"____"<<a<<"____"<<endl;
    cout<<"____"<<b<<"____"<<endl;

    occHandle(Ng_MeshVS_DataSource3D) OCCBkgMeshDS = new Ng_MeshVS_DataSource3D(QString::fromStdString(b),
                                                                                QString::fromStdString(a));

    if(OCCBkgMeshDS.IsNull()) return false;
    if(OCCBkgMeshDS->GetAllElements().Extent()==0) return false;
    if(OCCBkgMeshDS->GetAllNodes().Extent()<3) return false;

    cout<<"____"<<OCCBkgMeshDS->GetAllElements().Extent()<<"____"<<endl;
    cout<<"____"<<OCCBkgMeshDS->GetAllNodes().Extent()<<"____"<<endl;
    exit(1);

    //! ------------------------------------------------------
    //! convert the OCC volume mesh into Eigen representation
    //! ------------------------------------------------------
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    V.resize(NbPoints,4);
    TColStd_MapIteratorOfPackedMapOfInteger it;
    int NbNodes;
    MeshVS_EntityType type;
    double buf[3];
    TColStd_Array1OfReal coords(*buf,1,3);
    for(TColStd_MapIteratorOfPackedMapOfInteger it(OCCBkgMeshDS->GetAllNodes());it.More();it.Next())
    {
        int globalID = it.Key();
        int localNodeID = OCCBkgMeshDS->myNodesMap.FindIndex(globalID);
        OCCBkgMeshDS->GetGeom(globalID,false,coords,NbNodes,type);
        int pos = localNodeID-1;
        V(pos,0) = coords(1);
        V(pos,1) = coords(2);
        V(pos,2) = coords(3);
    }

    int NE = OCCBkgMeshDS->GetAllElements().Extent();
    F.resize(NE,4);

    int bufn[3];
    TColStd_Array1OfInteger nodeIDs(*bufn,1,3);
    for(TColStd_MapIteratorOfPackedMapOfInteger it(OCCBkgMeshDS->GetAllElements());it.More();it.Next())
    {
        int globalElementID = it.Key();
        int localElementID = OCCBkgMeshDS->myElementsMap.FindIndex(globalElementID);
        OCCBkgMeshDS->GetNodesByElement(globalElementID,nodeIDs,NbNodes);

        int pos = localElementID-1;
        for(int i=1; i<=NbNodes; i++)
        {
            int curGlobalNodeID = nodeIDs.Value(i);
            int curLocalNodeID = OCCBkgMeshDS->myNodesMap.FindIndex(curGlobalNodeID);
            F(pos,i-1) = curLocalNodeID-1;
        }
    }

    PyMesh::MshSaver aMshSaver("D:/sizingFunctionMsh.msh",false);
    aMshSaver.save_header();
    aMshSaver.save_nodes(V);
    aMshSaver.save_elements(F,PyMesh::MshSaver::TET);
    //aMshSaver.save_scalar_field("Mesh sizing function",);
    return true;
}

//! -------------------------------------
//! function: createSupportFileDirectory
//! details:
//! -------------------------------------
void tetWildMesher::createSupportFilesDirectory()
{
    //! -------------------------------
    //! get the path of the executable
    //! -------------------------------
    TCHAR buffer[MAX_PATH] = { 0 };
    GetModuleFileName( NULL, buffer, MAX_PATH );
    std::wstring::size_type pos = std::wstring(buffer).find_last_of(L"\\/");
    std::wstring currentPathWS = std::wstring(buffer).substr(0, pos);
    string currentPath = string(currentPathWS.begin(), currentPathWS.end());

    //! -------------------------------------------
    //! create the directory for the support files
    //! -------------------------------------------
    mySupportFilesDirectory = currentPath+"\\supportFiles";

    //! ---------------------------
    //! remove if already existing
    //! ---------------------------
    experimental::filesystem::path path(mySupportFilesDirectory.c_str());
    if(dirExists(mySupportFilesDirectory)==true) experimental::filesystem::remove_all(path);
    experimental::filesystem::create_directory(path);
}

//! --------------------------------------
//! function: removeSupportFilesDirectory
//! detail:
//! --------------------------------------
void tetWildMesher::removeSupportFilesDirectory()
{
    experimental::filesystem::path path(mySupportFilesDirectory.c_str());
    experimental::filesystem::remove_all(path);
}
