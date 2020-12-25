#define TETWILD_PATH "D:/Work/Qt/build_simSpace/release/FloatTetwild_bin.exe"

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

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
tetWildMesher::tetWildMesher(QObject *parent):QObject(parent),
    myWillRunOnDisk(true),
    myAbsoluteOutputFilePath("")
{
    myTetWildProcess = new QProcess(this);
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
void tetWildMesher::setParameters(float l_rel, float eps_rel)
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

    //cout<<"@____input file absolute  path: "<<absoluteInputFilePath.toStdString()<<endl;
    //cout<<"@____output file absolute path: "<<absoluteOuputFilePath.toStdString()<<endl;
    //cout<<"@____envelope relative size: "<<myMeshParam.eps_rel<<endl;
    //cout<<"@____initial relative edge length: "<<myMeshParam.initial_edge_len_rel<<endl;

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
                                   std::vector<std::vector<double>> &sampledPoints)
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

        std::vector<double>* samplingParameters =(std::vector<double>*)(parametersForSampling);
        double elementSizeOnFace = samplingParameters->at(0);
        Q_UNUSED(elementSizeOnFace);

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
        //! for the moment a 100x100 grid is used ...
        //! --------------------------------------------------
        double umin,umax,vmin,vmax;
        BRepTools::UVBounds(aFace,umin,umax,vmin,vmax);
        double Lu = umax-umin;
        double Lv = vmax-vmin;
        double du = Lu/100;
        double dv = Lv/100;

        //FILE *fp = fopen("D:\\sampled_face.txt","w");
        for(double ucur = umin; ucur<=umax; ucur += du)
        {
            for(double vcur = vmin; vcur<=vmax; vcur +=dv)
            {
                const gp_Pnt &curPoint = CurGeomSurface->Value(ucur,vcur);
                double x = curPoint.X();
                double y = curPoint.Y();
                double z = curPoint.Z();
                //fprintf(fp,"%lf\t%lf\t%lf\n",x,y,z);

                //! --------------------------------------------------
                //! check if the point is whithin the BB of the shape
                //! --------------------------------------------------
                if((x>=xmin_bb && x<=xmax_bb) && (y>=ymin_bb && y<=ymax_bb) && (z>=zmin_bb && z<=zmax_bb))
                {
                    std::vector<double> P{x,y,z};
                    sampledPoints.push_back(P);
                }
            }
        }
        //fclose(fp);
    }
        break;

    case TopAbs_EDGE:
    {
        //! ------------------------------------------------
        //! parameter for the edge sizing
        //! at the moment only one parameter in the vector:
        //! the element size on the edge (length)
        //! ------------------------------------------------
        std::vector<double> *samplingParameters =  (std::vector<double>*)(parametersForSampling);
        double size = samplingParameters->at(0);

        //FILE *fp = fopen("D:\\sampled_edge.txt","w");

        //! ---------
        //! the edge
        //! ---------
        const TopoDS_Edge &anEdge = TopoDS::Edge(aShape);
        if(BRep_Tool::Degenerated(anEdge)) break;
        BRepAdaptor_Curve BRepAdaptor(anEdge);
        GeomAdaptor_Curve curve = BRepAdaptor.Curve();

        //! -------------------
        //! length of the edge
        //! -------------------
        CPnts_AbscissaPoint CP;
        CP.Init(curve);
        double s_old = BRepAdaptor.FirstParameter();
        const gp_Pnt &P_onEdge = BRepAdaptor.Value(s_old);
        double x = P_onEdge.X();
        double y = P_onEdge.Y();
        double z = P_onEdge.Z();
        std::vector<double> P{x,y,z};
        sampledPoints.push_back(P);

        double Lcurrent = 0;
        double L = CP.Length(curve,1e-2);

        for(double s=s_old-size;Lcurrent<=L+size;)
        {
            CP.Perform(size,s,1e-2);

            const gp_Pnt &P1_onEdge = BRepAdaptor.Value(s);
            x = P1_onEdge.X();
            y = P1_onEdge.Y();
            z = P1_onEdge.Z();
            std::vector<double> P1{x,y,z};
            //fprintf(fp,"%lf\t%lf\t%lf\n",x,y,z);
            sampledPoints.push_back(P1);

            s_old=s;
            s=CP.Parameter();
            Lcurrent += fabs(s-s_old);
        }
        //fclose(fp);
    }
        break;

    case TopAbs_VERTEX:
    {
        //! ---------------------------------------------
        //! parameter for the mesh sizing at the vertex.
        //! one parameters in the vector: pinball radius
        //! ---------------------------------------------
        std::vector<double> *samplingParameters = (std::vector<double>*)(parametersForSampling);

        //! radial, theta, phi divisions
        const double PI = 3.141592654;
        const int radialDiv = 10;
        const int thetaDiv = 20;
        const int phiDiv = 10;

        //! pinball
        double R = samplingParameters->at(0);
        double dr = R/radialDiv;

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
                    std::vector<double> P1{x,y,z};
                    sampledPoints.push_back(P1);
                }
            }
        }
    }
        break;
    }
}


//! -------------------------------
//! function: writeMeshSizingField
//! details:
//! -------------------------------
//! -------
//! pymesh
//! -------
#include <MshSaver.h>
void tetWildMesher::writeMeshSizingField(int bodyIndex)
{
    struct point
    {
        double x,y,z;
        double value;
    };

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
    std::vector<double> parametersFaces;
    void *parametersForSamplingFaces = (void*)(&parametersFaces);
    //! to be filled ...
    //!
    int NbFaces = myMeshDB->MapOfBodyTopologyMap.value(bodyIndex).faceMap.Extent();
    for(int faceNr = 1; faceNr<=NbFaces; faceNr++)
    {
        bool isMeshControlOnFace = myMeshDB->MapOfIsFaceModified.getValue(bodyIndex,faceNr);
        if(isMeshControlOnFace == false) continue;
        std::vector<std::vector<double>> sampledPoints;
        TopoDS_Shape aFace = myMeshDB->MapOfBodyTopologyMap.value(bodyIndex).faceMap.FindKey(faceNr);
        if(aFace.IsNull()) continue;
        this->sampleGeometry(aFace,parametersForSamplingFaces,sampledPoints);

    }

    //! -----------------------
    //! mesh controls on edges
    //! -----------------------
    std::vector<double> parametersEdges;
    void *parametersForSamplingEdges = (void*)(&parametersEdges);
    //! to be filled ....

    int NbEdges = myMeshDB->MapOfBodyTopologyMap.value(bodyIndex).edgeMap.Extent();
    for(int edgeNr = 1; edgeNr<=NbEdges; edgeNr++)
    {
        bool isMeshControlOnEdge = myMeshDB->MapOfIsEdgeModified.getValue(bodyIndex,edgeNr);
        if(isMeshControlOnEdge == false) continue;
        std::vector<std::vector<double>> sampledPoints;
        TopoDS_Shape anEdge = myMeshDB->MapOfBodyTopologyMap.value(bodyIndex).edgeMap.FindKey(edgeNr);
        if(anEdge.IsNull()) continue;
        this->sampleGeometry(anEdge,parametersForSamplingEdges,sampledPoints);

    }

    //! ------------------------
    //! mesh controls on points
    //! ------------------------
    std::vector<double> parametersPoints;
    void *parametersForSamplingPoints = (void*)(&parametersPoints);
    //! to be filled ...

    int NbVertex = myMeshDB->MapOfBodyTopologyMap.value(bodyIndex).vertexMap.Extent();
    for(int vertexNr = 1; vertexNr<=NbEdges; vertexNr++)
    {
        bool isMeshControlOnVertex = myMeshDB->MapOfIsVertexModified.getValue(bodyIndex,vertexNr);
        if(isMeshControlOnVertex == false) continue;
        std::vector<std::vector<double>> sampledPoints;
        TopoDS_Shape aVertex = myMeshDB->MapOfBodyTopologyMap.value(bodyIndex).vertexMap.FindKey(vertexNr);
        if(aVertex.IsNull()) continue;
        this->sampleGeometry(aVertex,parametersForSamplingPoints,sampledPoints);

    }

}
