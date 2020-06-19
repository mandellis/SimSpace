//#define TETWILD_PATH "D:/Work/Qt/pro_25.0_OCC7.3.0/ftetwild/bin/FloatTetwild_bin.exe"
//#define TETWILD_PATH "D:/Work/Qt/pro_25.0_OCC7.3.0/tetwild/bin/TetWild.exe"
#define TETWILD_PATH "C:/ExperimentalMesher/ExpMesher.exe"
#include "tetwildmesher.h"

//! ----------------
//! custom includes
//! ----------------
#include "ccout.h"
#include <mshconvert.h>

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

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
tetWildMesher::tetWildMesher(QObject *parent):QObject(parent)
{
    //! -------------------------------------
    //! QProcess for running TetWild on disk
    //! -------------------------------------
    myWillRunOnDisk = true;
    myTetWildProcess = new QProcess(this);
    myTetWildProcess->setProgram(TETWILD_PATH);
    myAbsoluteOutputFilePath = "";
}

//! --------------------------
//! function: setTriangleSoup
//! details:
//! --------------------------
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
    QFile file; file.setFileName(absoluteInputFilePath);
    if(!file.exists() || absoluteInputFilePath.isEmpty()) return false;

    QString tmp(absoluteInputFilePath);
    QString relativeInputFileName = tmp.split("/").last();
    tmp.chop(relativeInputFileName.length());
    QString workingDir(tmp);
    QString extensionFreeRelativeInputFileName = relativeInputFileName;
    extensionFreeRelativeInputFileName.chop(4);
    QString absoluteOuputFilePath = workingDir+extensionFreeRelativeInputFileName+"_tetWild.msh";

    myAbsoluteOutputFilePath = absoluteOuputFilePath;

    QList<QString> arguments;
    arguments<<"-l"<<QString("%1").arg(myMeshParam.initial_edge_len_rel);
    arguments<<"-e"<<QString("%1").arg(myMeshParam.eps_rel);
    arguments<<"--input"<<absoluteInputFilePath;
    arguments<<"--output"<<absoluteOuputFilePath;

    //! ------------------------------
    //! diagnostic - console messages
    //! ------------------------------
    cout<<"@____running Experimental mesher with arguments: ";
    for(int i=0; i<arguments.length(); i++) cout<<arguments.at(i).toStdString()<<" ";
    cout<<endl;

    cout<<"@____input file absolute  path: "<<absoluteInputFilePath.toStdString()<<endl;
    cout<<"@____output file absolute path: "<<absoluteOuputFilePath.toStdString()<<endl;
    cout<<"@____envelope relative size: "<<myMeshParam.eps_rel<<endl;
    cout<<"@____initial relative edge length: "<<myMeshParam.initial_edge_len_rel<<endl;

    disconnect(myTetWildProcess,SIGNAL(readyReadStandardOutput()),this,SLOT(readTetWildProcess()));
    connect(myTetWildProcess,SIGNAL(readyReadStandardOutput()),this,SLOT(readTetWildProcess()));

    myTetWildProcess->setArguments(arguments);
    myTetWildProcess->start();
    int exitCode = myTetWildProcess->exitCode();
    myTetWildProcess->waitForFinished(-1);

    cout<<"@____\"Experimental mesher\" process exit code: "<<exitCode<<endl<<endl;

    if(exitCode==0) return true;
    else return false;
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

//! -----------------------------
//! function: readTetWildProcess
//! details:
//! -----------------------------
void tetWildMesher::readTetWildProcess()
{
    QByteArray msg = myTetWildProcess->readAllStandardOutput();

    //! send to std output
    cout<<msg.toStdString()<<endl;
}
