#ifndef STLDOCTOR_H
#define STLDOCTOR_H

//! -------------------------------------------
//! redefinition of opencascade Handle() macro
//! -------------------------------------------
#include "occhandle.h"

//! must contain the .exe file
//#define ADMESH_PATH "C:/MinGW/msys/1.0/home/User/admesh-0.98.2/admesh.exe"
//#define MESHFIX_PATH "D:/Work/MeshFIX_WT/bin64/MeshFix_WT.exe"
#define ADMESH_PATH "C:/ADMesh_22.0/admesh.exe"
#define MESHFIX_PATH "C:/MeshFix/MeshFix_WT.exe"

//! ---
//! Qt
//! ---
#include <QObject>
#include <QString>

//! ----
//! OCC
//! ----
#include <StlMesh_Mesh.hxx>
#include <StlMesh_SequenceOfMeshTriangle.hxx>
#include <gp_Pnt.hxx>

//! ----------------
//! custom includes
//! ----------------
#include <mesh.h>

class QProcess;

class STLdoctor: public QObject
{
    Q_OBJECT

public:

    STLdoctor(QObject *parent=0);

private:

    //! ADMesh QProcess
    QProcess *ADMesh;

    //! MeshFix QProcess
    QProcess *MeshFix;

public:

    //! perform ADMesh exact check
    bool perform(const QString &inputFile, const QString &outputFile);

    //! perform ADMesh nearby
    bool performADMeshNearby(const QString &inputFile, const QString &outputFile);

    //! perform MeshFix
    bool performMeshFix(const QString &inputFile, const QString &outputFile);

    //! perform mesh simplification through edge collapse
    void simplifyMesh(const QString &meshFileIn,
                      const QString &meshFileOut,
                      int method=0,
                      double surfaceMeshSizeFactor=0.95,
                      double pairDistance=0.01,
                      const std::vector<mesh::tolerantPoint> &fixedPoints = std::vector<mesh::tolerantPoint>());

    //! list small triangles
    //static std::list<occHandle(StlMesh_MeshTriangle)> getSmallTriangles(const occHandle(StlMesh_Mesh) &aSTLMesh_Mesh, double AreaThreshold);

private slots:

    void readADMeshStdOutput();
    void readMeshFixStdOutput();
};

#endif // STLDOCTOR_H
