//! ----------------
//! custom includes
//! ----------------
#include "simplifymesh.h"
#include "stldoctor.h"
#include "qprogressindicator.h"
#include "qprogressevent.h"

//! ---
//! Qt
//! ---
#include <QByteArray>
#include <QProcess>
#include <QFile>
#include <QApplication>

//! ----
//! C++
//! ----
#include <iostream>
using namespace std;

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
STLdoctor::STLdoctor(QObject *parent): QObject(parent)
{
    //! ---------------------------
    //! create the ADMesh QProcess
    //! ---------------------------
    ADMesh = new QProcess(this);
    connect(ADMesh,SIGNAL(readyReadStandardOutput()),this,SLOT(readADMeshStdOutput()));

    //! ----------------------------
    //! create the MeshFix QProcess
    //! ----------------------------
    MeshFix = new QProcess(this);
    connect(MeshFix,SIGNAL(readyReadStandardOutput()),this,SLOT(readMeshFixStdOutput()));

    //! ----------------------------
    //! init the progress indicator
    //! ----------------------------
    myProgressIndicator = Q_NULLPTR;
}

//! -------------------------------
//! function: setProgressIndicator
//! details:
//! -------------------------------
void STLdoctor::setProgressIndicator(QProgressIndicator *aProgressIndicator)
{
    myProgressIndicator = aProgressIndicator;
}

//! ------------------
//! function: perform
//! details:
//! ------------------
bool STLdoctor::perform(const QString &inputFile, const QString &outputFile)
{
    //cout<<"/-----------------------------------/"<<endl;
    //cout<<"/ STLdoctor: performing exact check /"<<endl;
    //cout<<"/-----------------------------------/"<<endl;

    //! ---------------------------------
    //! set the program for the QProcess
    //! ---------------------------------
    QString program = QString(std::getenv("ADMESH_PATH")).append("\\admesh.exe");
    ADMesh->setProgram(program);

    //! ----------------------------
    //! arguments - all checks done
    //! ----------------------------
    QStringList arguments;
    arguments<<inputFile<<QString("-a")<<outputFile;

    //! -----------------------------------------------
    //! if you want to perform a custom set of checks
    //! using ADMesh
    //!
    //! --exact --nearby --tolerance=<tolerance value>
    //! --iterations=<number of iterations> --exact
    //! --normal-directions --normal-values
    //! -----------------------------------------------

    //! ------------------
    //! set the arguments
    //! ------------------
    ADMesh->setArguments(arguments);

    cout<<"STLdoctor::perform()->____ADMesh starting____"<<endl;
    ADMesh->start();
    bool isDone = ADMesh->waitForFinished(-1);

    return isDone;
}

//! ------------------------------
//! function: performADMeshNearby
//! details:
//! ------------------------------
bool STLdoctor::performADMeshNearby(const QString &inputFile, const QString &outputFile)
{
    cout<<"/------------------------------/"<<endl;
    cout<<"/ STLdoctor: performing nearby /"<<endl;
    cout<<"/------------------------------/"<<endl;

    //! define the program: for the moment use command line through ADMesh
    QString program = QString(ADMESH_PATH);

    //! set the program for the QProcess
    ADMesh->setProgram(program);

    //! ------------------
    //! ADMesh parameters
    //! ------------------
    QStringList arguments;
    arguments<<inputFile<<
               QString("--nearby")<<
               QString("--tolerance=1e-10")<<
               QString("--increment=1e-10")<<
               QString("--iteration=10000")<<
               QString("--normal-directions")<<
               QString("-a")<<outputFile;

    //! ------------------
    //! set the arguments
    //! ------------------
    ADMesh->setArguments(arguments);
    ADMesh->start();
    bool isDone = ADMesh->waitForFinished(-1);
    return isDone;
}

//! ------------------------------
//! function: readADMeshStdOutput
//! details:
//! ------------------------------
void STLdoctor::readADMeshStdOutput()
{
    cout<<"STLdoctor::readStdOutput()->____function called____"<<endl;
    QByteArray msg = ADMesh->readAllStandardOutput();
    cout<<msg.toStdString()<<endl;
}

//! ------------------------------
//! function: readADMeshStdOutput
//! details:
//! ------------------------------
void STLdoctor::readMeshFixStdOutput()
{
    //cout<<"STLdoctor::readMeshFixOutput()->____function called____"<<endl;
    QByteArray msg = MeshFix->readAllStandardOutput();
    cout<<msg.toStdString()<<endl;
    if(myProgressIndicator!=Q_NULLPTR)
    {
        QProgressEvent *e = new QProgressEvent(QProgressEvent_None,-1,-1,0,QString::fromStdString(msg.toStdString()),QProgressEvent_Init,0,100,0);
        QApplication::postEvent(myProgressIndicator,e);
        QApplication::processEvents();
    }
}

//! --------------------------------------
//! function: simplifyMesh
//! details:  interface to "SimplifyMesh"
//! --------------------------------------
void STLdoctor::simplifyMesh(const QString &meshFileIn,
                             const QString &meshFileOut,
                             int method,
                             double surfaceMeshSizeFactor,
                             double pairDistance,
                             const std::vector<mesh::tolerantPoint> &fixedPoints)
{
    SimplifyMesh sm;
    int nodes, elements;
    //! ---------------------------------------------
    //! set the fixed points before loading the .stl
    //! ---------------------------------------------
    sm.setFixedPoints(fixedPoints);
    sm.load_Stl(meshFileIn.toStdString().c_str(), nodes, elements);

    switch(method)
    {
    case 0:
    {
        //! use surface mesh reduction factor
        int count = int(surfaceMeshSizeFactor*nodes);
        double aggressiveness = 5.0;
        sm.simplify_mesh(count,aggressiveness,true);
        sm.write_Stl(meshFileOut.toStdString().c_str());
    }
        break;

    case 1:
    {
        //! use pair distance
        sm.simplify_mesh_lossless(pairDistance,true);
        sm.write_Stl(meshFileOut.toStdString().c_str());
    }
        break;
    }
}

//! -------------------------
//! function: performMeshFix
//! details:
//! -------------------------
bool STLdoctor::performMeshFix(const QString &inputFile, const QString &outputFile)
{
    cout<<"/-------------------------------/"<<endl;
    cout<<"/ STLdoctor: perform MeshFix_VT /"<<endl;
    cout<<"/ with geometry tag propagation /"<<endl;
    cout<<"/-------------------------------/"<<endl;

    //! ---------------------------------
    //! set the program for the QProcess
    //! ---------------------------------
    QString program = QString::fromLatin1(std::getenv("MESHFIX_VT_PATH")).append("\\MeshFix_WT.exe");
    MeshFix->setProgram(program);

    //! ---------------------------------------------------
    //! arguments - use "-j" for generating an .stl output
    //! ---------------------------------------------------
    QStringList arguments;
    arguments<<inputFile<<outputFile<<QString("-j");

    //! ------------------
    //! set the arguments
    //! ------------------
    MeshFix->setArguments(arguments);

    cout<<"STLdoctor::perform()->____MeshFix starting____"<<endl;
    MeshFix->start();
    bool isDone = MeshFix->waitForFinished(-1);

    return isDone;
}
