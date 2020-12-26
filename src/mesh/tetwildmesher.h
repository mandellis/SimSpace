#ifndef TETWILDMESHER_H
#define TETWILDMESHER_H

//! -----------------------------------------
//! interface for the TetWild mesher library
//! -----------------------------------------
#include <tetwild.h>
#include <Args.h>

//! Eigen
#include <Eigen/Dense>

//! ---
//! Qt
//! ---
#include <QObject>
#include <QString>

//! ----------------
//! custom includes
//! ----------------
#include <ng_meshvs_datasource3d.h>

//! ----
//! OCC
//! ----
#include <TopoDS_Shape.hxx>

class QProcess;
class meshDataBase;

class tetWildMesher: public QObject
{
    Q_OBJECT

private:

    struct point
    {
        point (double aX, double aY, double aZ, double aValue): x(aX),y(aY),z(aZ),value(aValue){;}
        point (const point &aPoint) { x = aPoint.x; y = aPoint.y; z = aPoint.z; value = aPoint.value; }
        double x,y,z,value;
    };

    meshDataBase *myMeshDB;

    Eigen::MatrixXd myVertices;
    Eigen::MatrixXi myFaces;
    tetwild::Args myMeshParam;
    bool myWillRunOnDisk;

    QProcess *myTetWildProcess;
    QString myAbsoluteOutputFilePath;

public:

    //! constructor
    tetWildMesher(QObject *parent=Q_NULLPTR);

    //! constructor
    tetWildMesher(const Eigen::MatrixXd &VI, const Eigen::MatrixXi &FI): myVertices(VI), myFaces(FI){;}

    //! set input triangle soup
    void setInput(const Eigen::MatrixXd &VI, const Eigen::MatrixXi &FI);

    //! set arguments
    void setParameters(float l_rel, float eps_rel);

    //! set mesh data base
    void setDataBase(meshDataBase *mDB);

    //! setRunOnDisk
    void setRunOnDisk(bool flag) { myWillRunOnDisk = flag; }

    //! perform volume
    bool perform_inMemory(Eigen::MatrixXd &VO, Eigen::MatrixXi &TO);

    //! perform volume on disk
    bool perform_onDisk(const QString &absoluteInputFilePath);

    //! retrieve the 3D mesh data source from disk
    bool retrieveMeshDataSourceFromDisk(opencascade::handle<Ng_MeshVS_DataSource3D> &mesh3D_DataSource);

    //! extract surface
    bool extractSurface(const Eigen::MatrixXd &VI, const Eigen::MatrixXi &TI,
                        Eigen::MatrixXd &VS, Eigen::MatrixXi &FS);

private:

    //! write mesh sizing field
    void writeMeshSizingField(int bodyIndex);

    //! sample geometry
    void sampleGeometry(const TopoDS_Shape &aShape,
                        void *parametersForSampling,
                        std::vector<tetWildMesher::point> &sampledPoints);


    //! get path of the executable - helper
    std::string getPathOfExecutable();

private slots:

    //! capture process output
    void readTetWildProcess();
};

#endif // TETWILDMESHER_H
