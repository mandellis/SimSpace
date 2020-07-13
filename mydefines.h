#ifndef MYDEFINES
#define MYDEFINES

//! C++
#include <iostream>
#include <string>

//! Qt
#include <QMap>
#include <Quantity_NameOfColor.hxx>
#include <QMetaType>

//! OCC
#include <TopAbs_ShapeEnum.hxx>

//! custom includes
#include "ais_colorscaleextended.h"
#include "listofmesh.h"
#include "hash_c.h"

#define SYSTEM_PROGRAM_DATA "C:/ProgramData"
#define GEOMETRY_DIR "D:/Work/Qt/geometries"
#define USERDATA_DIR "D:/Work/Qt/user data"
#define IGES_FILES "Iges files (*.iges *.igs )"
#define STEP_FILES "Step files (*.step *.stp)"
#define STL_FILES "Stl files (*.stl)"
#define INP_FILES "Solver input (*.inp)"
#define CCX_FILES "CCX files (*.frd)"

#define PNG_FILES ".png"
#define JPG_FILES ".jpg"
#define BMP_FILES ".bmp"
#define GIF_FILES ".gif"

#define DEF_R1  255
#define DEF_G1  255
#define DEF_B1  255

#define DEF_R2  255
#define DEF_G2  255
#define DEF_B2  255

#define GIL_FILES "Simulation data (*.gil)"
#define NODE_FILES "Node file (*.node)"
#define APPNAME "SimSpace - Simulation Environment"

#define TRANSPARENCY_IN_WORKING_MODE_MODEL 0.0
#define TRANSPARENCY_IN_WORKING_MODE_MESH 0.8
#define TRANSPARENCY_IN_WORKING_MODE_CONTACT 0.9
#define TRANSPARENCY_IN_WORKING_MODE_CONTACT_PARTS_IN_CONTACT 0.1
#define TRANSIENT_MESSAGE_TIMEOUT 1000

#define MESH_CONTROL_COLOR Quantity_NOC_VIOLET

#define FRICTIONLESS_SUPPORT_COLOR Quantity_NOC_BLUE1
#define CYLINDRICAL_SUPPORT_COLOR Quantity_NOC_BLUE1
#define FIXED_SUPPORT_COLOR Quantity_NOC_BLUE1
#define FORCE_COLOR Quantity_NOC_RED
#define REMOTE_FORCE_COLOR Quantity_NOC_ORANGE
#define MOMENT_COLOR Quantity_NOC_RED
#define PRESSURE_COLOR Quantity_NOC_RED
#define DISPLACEMENT_COLOR Quantity_NOC_YELLOW
#define REMOTE_DISPLACEMENT_COLOR Quantity_NOC_BLUE4
#define MASTER_COLOR Quantity_NOC_BLUE1
#define SLAVE_COLOR Quantity_NOC_RED
#define NAMED_SELECTION_COLOR Quantity_NOC_GREEN1
#define THERMALCONDITION_COLOR Quantity_NOC_RED
#define CONNECTION_GROUP_COLOR Quantity_NOC_BLUE1
#define THERMAL_BOUNDARY_COLOR Quantity_NOC_RED
#define COMPRESSION_ONLY_SUPPORT_COLOR Quantity_NOC_BLUE2

#define SURFACE_MESH_ONLY_EDGECOLOR Quantity_NOC_BLUE1
#define VOLUME_MESH_EDGECOLOR Quantity_NOC_BLACK
#define PRISMATIC_MESH_EDGE_COLOR Quantity_NOC_RED
#define NODE_MESH_COLOR Quantity_NOC_BLACK

#define INITIAL_NUMBER_OF_COLORBOX_LEVELS 10
#define EXPFORMAT_PRECISION 9

#define OMF_TEMPLIC "OMF-7.1          PUZZLE_DIE  06.07.2017 31.07.2017      0    0 185a648e4b32827fb33d594c28222bb5dc64d96953a3db1d5f2a8c2f41a7dc3e Co.Stamp s.r.l. / Temporary license"
#define EMESH_TEMPLIC "EMESH-7.1        PUZZLE_DIE  06.07.2017 31.07.2017      0    0 cad0081eeb6c8f1d14270117deab6d2b0e09c48898d6599833471730d9d1d4b2 Co.Stamp s.r.l. / Temporary license"
//#define OMF_DEFLIC "          PUZZLE_DIE         0        0      0    0 8d19e677b16fb8e2e0c2301df07656dee5d2ba44cc03eebeeeca55d64ed3ed2c Co.Stamp s.r.l. / Master key"
//#define EMESH_DEFLIC "        PUZZLE_DIE         0        0      0    0 3f85a38fd58cab86c9f3831e6025957ba6f63394b841eaf92f88da8e13ddd1f7 Co.Stamp s.r.l. / Master key"
#define OMF_DEFLIC =  "          PUZZLE_DIE         0        0      0    0 613699cc791fa76a08d0c5f8534c2db8eb89565b94cfcdcbc600ba4639221061 Co.Stamp s.r.l. / Master key";
#define EMESH_DEFLIC  "        PUZZLE_DIE         0        0      0    0 5ba43f346b9dfbd4ece27b57443b309961676c9bb98db05e041047c22e93a907 Co.Stamp s.r.l. / Master key";

//! ----------
//! meshParam
//! ----------
struct meshParam
{
    //! ----------------------------------
    //! global controls (at "body" level)
    //! ----------------------------------
    double grading = 0.3;
    double minElementSize = 1;
    double maxElementSize = 100;
    int isSurfaceQuad = 0;
    int isVolumeHexa = 0;
    bool isSecondOrder = false;
    bool isElementStraight = false;
    bool isVolume = true;
    bool isHealingOn;
    bool isDefeaturingOn;
    bool defeaturingBy;                        //! 0 => by "Surface mesh size" 1 => by "Pair distance"
    bool defeaturingParameterValue;            //! "Surface mesh size" or "Pair distance"

    //! ---------------
    //! local controls
    //! ---------------
    QMap<int, double> faceSize;
    QMap<int, double> edgeSize;
    QMap<int, int> edgeNbDivisions;
    QMap<int, int> edgeTypeOfSizing;
    QMap<int, double> vertexSize;
    QMap<int, double> pinballSize;

    //! -----------------
    //! prismatic layers
    //! -----------------
    QMap<int,QList<int>> prismaticFaces;

    //! -----------
    //! operator =
    //! -----------
    meshParam operator=(const meshParam &other)
    {
        grading=other.grading;
        minElementSize=other.minElementSize;
        maxElementSize=other.maxElementSize;
        isSurfaceQuad=other.isSurfaceQuad;
        isVolumeHexa=other.isVolumeHexa;
        isSecondOrder=other.isSecondOrder;
        isElementStraight=other.isElementStraight;
        isVolume=other.isVolume;

        isHealingOn = other.isHealingOn;
        isDefeaturingOn = other.isDefeaturingOn;
        defeaturingBy = other.defeaturingBy;
        defeaturingParameterValue = other.defeaturingParameterValue;

        faceSize=other.faceSize;
        edgeSize=other.edgeSize;
        edgeNbDivisions=other.edgeNbDivisions;
        edgeTypeOfSizing=other.edgeTypeOfSizing;

        vertexSize=other.vertexSize;
        pinballSize=other.pinballSize;

        for(QMap<int,QList<int>>::const_iterator it=other.prismaticFaces.begin(); it!=other.prismaticFaces.cend(); ++it)
        {
            int bodyIndex = it.key();
            const QList<int> &faceNbs = it.value();
            for(int i=0; i<faceNbs.length(); i++) prismaticFaces.insert(bodyIndex,faceNbs);
        }

        return *this;
    }
};

#endif // MYDEFINES
