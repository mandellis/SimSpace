#ifndef NETTGENTOOL2_H
#define NETTGENTOOL2_H
#define HANDLE_NGLIB_SYSTEM_EXCEPTIONS

//! ------
//! nglib
//! ------
#include <vector>
namespace nglib
{
    #include <nglib.h>
}
using namespace nglib;

//! --------------------------------------
//! replacement for opencascade::handle<>
//! --------------------------------------
#include <occhandle.h>

//! ----------------
//! custom includes
//! ----------------
#include <userMessage.h>
#include "meshdatabase.h"

//! ---
//! Qt
//! ---
#include <QVariant>
#include <QObject>
#include <QTime>

//! ----
//! C++
//! ----
#include <vector>

#ifdef HANDLE_NGLIB_SYSTEM_EXCEPTIONS
//! ----------------
//! system specific
//! ----------------
#include <stdio.h>
#include <windows.h>
#include <eh.h>
#include "se_exception.h"

extern void trans_func_N(unsigned int, EXCEPTION_POINTERS*);
#endif

class MeshVS_DataSource;
class Ng_MeshVS_DataSource2D;
class Ng_MeshVS_DataSource3D;
class QProgressIndicator;

class NettgenTool2: public QObject
{
    Q_OBJECT

private:

#ifdef HANDLE_NGLIB_SYSTEM_EXCEPTIONS
    struct STLMeshing_InOutParameters
    {
        bool isVolume;
        void *netgenTessellation;
        void *externallyDefinedEdges;
        Ng_Meshing_Parameters *meshingParameters;
        void *localMeshSizingControls;
        void *surfaceMeshPointsP;
        void *surfaceMeshElementsP;
        void *volumeMeshPointsP;
        void *volumeMeshElementsP;
        void *surfaceMeshDataP;
        void *unassociatedSurfaceElements;
    };

    struct BRepMeshing_InOutParameters
    {
        bool isVolume;
        void *shape;
        void *validFaceTags;    //! unused
        Ng_Meshing_Parameters *mp;
        void *localMeshSizingControls;
        void *surfaceMeshPointsP;
        void *surfaceMeshElementsP;
        void *volumeMeshPointsP;
        void *volumeMeshElementsP;
        void *faceMeshesP;
    };

    struct meshing_InOutParameters
    {
        void *inputSurfaceMeshPoints;
        void *inputSurfaceMeshElements;
        Ng_Meshing_Parameters *mp;
        void *localMeshSizingControls;
        void *surfaceMeshPointsP;
        void *surfaceMeshElementsP;
        void *volumeMeshPointsP;
        void *volumeMeshElementsP;
    };
#endif

private:

    meshDataBase *myMDB;
    QProgressIndicator *myProgressIndicator;

public:

    //! ------------
    //! constructor
    //! ------------
    NettgenTool2(QObject* parent=0);

    //! -----------
    //! destructor
    //! -----------
    ~NettgenTool2()
    {
        cout<<"NettgenTool2::~NettgenTool2()->____DESTRUCTOR CALLED____"<<endl;
    }

    //! ------------------
    //! set mesh database
    //! ------------------
    void setMeshDataBase(meshDataBase *mDB)
    {
        myMDB = mDB;
        //checkDataBase();
    }

    //! --------
    //! perform
    //! --------
    userMessage perform(int bodyIndex,
                          bool isVolume,
                          opencascade::handle<MeshVS_DataSource> &outputSurfaceMeshDS,
                          opencascade::handle<MeshVS_DataSource> &outputVolumeMeshDS,
                            std::vector<vector<double>> &vecUnassociatedElements);

    //! --------------------------------------
    //! performBrep - works on a CAD geometry
    //! --------------------------------------
    userMessage performBrep(int bodyIndex,
                           bool isVolume,
                           occHandle(MeshVS_DataSource) &outputSurfaceMeshDS,
                           occHandle(MeshVS_DataSource) &outputVolumeMeshDS);

    //! ------------------------------
    //! mesh volume from surface mesh
    //! ------------------------------
    userMessage meshVolumeFromSurfaceMesh(int bodyIndex,
                                            const opencascade::handle<MeshVS_DataSource> &aSurfaceMeshDS,
                                            occHandle(MeshVS_DataSource) &outputVolumeMeshDS);


    //! -----------------------
    //! set progress indicator
    //! -----------------------
    void setProgressIndicator(QProgressIndicator *aProgressIndicator);

private:

    //! ----------------
    //! sample geometry
    //! ----------------
    void sampleGeometry(const TopoDS_Shape &aShape,
                        const QVariant &parametersForSampling,
                        std::vector<std::vector<double>> &sampledPoints);

    //! -------------------------------
    //! init netgen meshing parameters
    //! -------------------------------
    Ng_Meshing_Parameters initNetgenMeshingParameters();

    //! ----------------
    //! data base check
    //! ----------------
    bool checkDataBase();

#ifdef HANDLE_NGLIB_SYSTEM_EXCEPTIONS
    //! ---------------------------------------------
    //! system specific - meshing from an .stl input
    //! ---------------------------------------------
    bool SystemErrorFunction(STLMeshing_InOutParameters inOutParameters);

    //! ---------------------------------------------
    //! system specific - meshing from an Brep input
    //! ---------------------------------------------
    bool SystemErrorFunction_1(BRepMeshing_InOutParameters inOutParameters);

    //! -----------------------------------------------------
    //! system specific - meshing from a netgen surface mesh
    //! -----------------------------------------------------
    bool SystemErrorFunction_2(meshing_InOutParameters inOutParameters);
#endif

signals:

    //! ----------------------------------------
    //! stop and start the netgen enquire timer
    //! signals collected by the mesherClass
    //! ----------------------------------------
    void requestStoppingNetgenEnquireTimer();
    void requestStartingNetgenEnquireTimer();

    //! ---------------------------------------------------------
    //! signals for driving the button of the progress indicator
    //! ---------------------------------------------------------
    void requestEnableStop();
    void requestDisableStop();

private:

    //! --------------------------
    //! function: clock
    //! details:  time diagnostic
    //! --------------------------
    QTime clock()
    {
        QTime time = QTime::currentTime();
        QString text = time.toString("hh::mm:ss:zzz");
        cout<<text.toStdString()<<endl;
        return time;
    }
};

#endif // NETTGENTOOL2_H
