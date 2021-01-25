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
#include <hash_c.h>

//! ----
//! OCC
//! ----
#include <TopoDS_Shape.hxx>

//! ----
//! C++
//! ----
#include <string>
#include <experimental/filesystem>

class QProcess;
class meshDataBase;

class tetWildMesher: public QObject
{
    Q_OBJECT

private:

    meshDataBase *myMeshDB;

    Eigen::MatrixXd myVertices;
    Eigen::MatrixXi myFaces;
    tetwild::Args myMeshParam;
    bool myWillRunOnDisk;

    QProcess *myTetWildProcess;
    QString myAbsoluteOutputFilePath;

    //! directory for the support files
    std::string mySupportFilesDirectory;

public:

    struct point
    {
        point (double aX, double aY, double aZ, double aValue): x(aX),y(aY),z(aZ),value(aValue){;}
        point (const point &aPoint) { x = aPoint.x; y = aPoint.y; z = aPoint.z; value = aPoint.value; }
        point operator = (const point &aPoint) { x = aPoint.x; y = aPoint.y; z = aPoint.z; value = aPoint.value; return *this; }
        bool operator < (const point &rhs) const
        {
            size_t seed0, seed1;
            seed0 = 0; seed1 = 0;
            hash_c<double>(seed0,x); hash_c<double>(seed1,rhs.x);
            hash_c<double>(seed0,y); hash_c<double>(seed1,rhs.y);
            hash_c<double>(seed0,z); hash_c<double>(seed1,rhs.z);
            if(seed0<seed1) return true;
            return false;
        }
        bool operator == (const point &rhs) const
        {
            if( x == rhs.x && y == rhs.y && z == rhs.z) return true;
            return false;
        }
        double x,y,z,value;
    };

    //! constructor
    tetWildMesher(QObject *parent=Q_NULLPTR);

    //! constructor
    tetWildMesher(const Eigen::MatrixXd &VI, const Eigen::MatrixXi &FI): myVertices(VI), myFaces(FI){;}

    //! destructor
    ~tetWildMesher();

    //! set input triangle soup
    void setInput(const Eigen::MatrixXd &VI, const Eigen::MatrixXi &FI);

    //! set arguments
    void setGlobalParameters(float l_rel, float eps_rel);

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

    //! write mesh sizing field
    void computeMeshSizingField(int bodyIndex, std::vector<tetWildMesher::point> &meshSizingField);

    //! write the mesh sizing field
    bool writeMeshSizingField(const std::string &backgroundMesh, const std::vector<tetWildMesher::point> &inputMeshSizingField);


private:

    //! sample geometry
    void sampleGeometry(const TopoDS_Shape &aShape,
                        void *parametersForSampling,
                        std::vector<tetWildMesher::point> &sampledPoints);

    //! get path of the executable - helper
    std::string getPathOfExecutable();

    //! create support file directory
    void createSupportFilesDirectory();
    void removeSupportFilesDirectory();

private slots:

    //! capture process output
    void readTetWildProcess();
};

#endif // TETWILDMESHER_H
