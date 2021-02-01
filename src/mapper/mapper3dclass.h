#ifndef MAPPER3DCLASS_H
#define MAPPER3DCLASS_H

//! -------------------------------------------
//! redefinition of opencascade Handle() macro
//! -------------------------------------------
#include "occhandle.h"

//! ---
//! Qt
//! ---
#include <QObject>
#include <QMap>

//! ----------------
//! custom includes
//! ----------------
#include <mesh.h>

//! ----
//! C++
//! ----
#include <vector>
#include <set>
#include <iostream>
using namespace std;

//! ----
//! OCC
//! ----
#include <MeshVS_DataSource.hxx>

class QProgressIndicator;

class Mapper3DClass: public QObject
{
    Q_OBJECT

public:

    //! --------------------------------------------------------
    //! a source point (x,y,z) with an attached list of scalars
    //! ---------------------------------------------------------
    struct nodeSource
    {
        double x[3];
        std::vector<double> vecVal;

        //! --------------------
        //! default constructor
        //! --------------------
        nodeSource()
        {
            x[0]=0; x[1]=0; x[2]=0;
            vecVal = std::vector<double>();
        }

        //! -----------
        //! operator =
        //! -----------
        nodeSource operator = (const nodeSource &rhs)
        {
            x[0] = rhs.x[0];
            x[1] = rhs.x[1];
            x[2] = rhs.x[2];
            vecVal.clear();
            vecVal = rhs.vecVal;
            return *this;
        }

        //! -----------------
        //! copy constructor
        //! -----------------
        inline nodeSource(const nodeSource &rhs)
        {
            x[0] = rhs.x[0];
            x[1] = rhs.x[1];
            x[2] = rhs.x[2];
            vecVal.clear();
            vecVal = rhs.vecVal;
        }

        //! ------------
        //! operator ==
        //! ------------
        bool operator == (const nodeSource &other)
        {
            //! check the coordinates
            if(x[0] != other.x[0]) return false;
            if(x[1] != other.x[1]) return false;
            if(x[2] != other.x[2]) return false;

            //! check the values
            for(int i=0; i<other.vecVal.size(); i++)
            {
                if (vecVal[i]!=other.vecVal[i]) return false;
            }
            return true;
        }
    };

    //! -----------------
    //! target mesh node
    //! -----------------
    struct targetMeshNode
    {
        double x,y,z;
        int nodeID;
        std::vector<double> values;
        double pinball;
        bool hasValue;

        //! -----------------
        //! copy constructor
        //! -----------------
        targetMeshNode(const targetMeshNode &other)
        {
            x = other.x; y = other.y; z = other.z;
            nodeID = other.nodeID;
            values.clear();
            values = other.values;
            pinball = other.pinball;
            hasValue = other.hasValue;
        }

        //! -----------
        //! contructor
        //! -----------
        targetMeshNode(double xp, double yp, double zp, int aNodeID)
        {
            x = xp; y = yp; z = zp;
            nodeID = aNodeID;
            values = std::vector<double>();
            pinball = 1e80;
            hasValue = false;
        }

        //! ---------------------
        //! function: setPinball
        //! ---------------------
        void setPinball(double aPinball) { pinball = aPinball; }

        //! ------------------------------------------------------------------
        //! function: put
        //! details:  set "values" usign source node data, if the source node
        //!           is within the region of influence (pinball)
        //! ------------------------------------------------------------------
        void put(const nodeSource &aNodeSource)
        {
            double d = sqrt(pow(x-aNodeSource.x[0],2)+pow(y-aNodeSource.x[1],2)+pow(z-aNodeSource.x[2],2));
            if(d<=pinball)
            {
                pinball = d;
                values.clear();
                values = aNodeSource.vecVal;
                hasValue = true;
            }
        }
    };


private:

    //! --------------------------------------------
    //! inline definition of a 3D array containing
    //! vectors of source elements
    //! --------------------------------------------
    struct bucket
    {
        std::vector<nodeSource> *vecNodes;
        int myNx, myNy, myNz;

        //! ------------------------------
        //! function: default constructor
        //! details:
        //! ------------------------------
        bucket()
        {
            myNx = myNy = myNz = 0;
            vecNodes = NULL;
        }

        //! -------------------
        //! function: allocate
        //! details:
        //! -------------------
        void allocate(int Nx, int Ny, int Nz)
        {
            if(Nx>=1 && Ny>=1 && Nz>=1)
            {
                myNx = Nx;
                myNy = Ny;
                myNz = Nz;
                vecNodes = new std::vector<Mapper3DClass::nodeSource>[Nx*Ny*Nz];
            }
            else vecNodes = NULL;
        }

        //! ---------------------------------------------
        //! function: insert
        //! details:  insert into the (i, j, k) position
        //! ---------------------------------------------
        void setValue(int i, int j, int k, Mapper3DClass::nodeSource aSourceNode)
        {
            //! x + y * WIDTH + z * WIDTH * DEPTH
            //! arr[x + width * (y + depth * z)]
            (vecNodes+i*myNy*myNz+j*myNz+k)->push_back(aSourceNode);
        }

        //! ---------------------------------------------
        //! function: value
        //! details:  return the vector of source points
        //! ---------------------------------------------
        std::vector<Mapper3DClass::nodeSource> value(int i, int j, int k)
        {
            return *(vecNodes+i*myNy*myNz+j*myNz+k);
        }

        //! ----------------------
        //! function: deallocate
        //! details:  free memory
        //! ----------------------
        void deallocate()
        {
            delete [] vecNodes;
        }

        //! -----------
        //! operator <
        //! -----------
        bool operator <(const bucket &other) const
        {
            size_t seed1, seed2;
            seed1 = seed2 = 0;
            hash_c<int>(seed1,myNx);
            hash_c<int>(seed1,myNy);
            hash_c<int>(seed1,myNz);

            hash_c<int>(seed2,other.myNx);
            hash_c<int>(seed2,other.myNy);
            hash_c<int>(seed2,other.myNz);

            return seed1<seed2;
        }
    };

public:

    //! --------------
    //! constructor I
    //! --------------
    Mapper3DClass(const occHandle(MeshVS_DataSource) &theTargetMeshDS, QObject *parent=0);

    //! --------------
    //! constructor II
    //! --------------
    Mapper3DClass(QObject *parent=0);

    //! -----------
    //! destructor
    //! -----------
    ~Mapper3DClass();

    //! --------------------
    //! set the target mesh
    //! --------------------
    void setTargetMesh(const occHandle(MeshVS_DataSource) &theTargetMeshDS);

    //! -------------
    //! setNbBuckets
    //! -------------
    void setNbBuckets(int NbBucketsX, int NbBucketsY, int NbBucketsZ);

    //! ---------------------
    //! set the source nodes
    //! ---------------------
    void setSource(const std::vector<Mapper3DClass::nodeSource> &vecSourceNodes);

    //! remap by target elements ... to do
    int remapByTargetElements();

    //! -------------------------------------
    //!  put the list of values on each node
    //! -------------------------------------
    void putScalarOnTargetNode(int interpolationTypeFunction=0);

    //! --------------------------------------
    //! retrieve resMap from multiInterpolate
    //! --------------------------------------
    void retrieveResMap(int pos);

    //! -----------
    //! get result
    //! -----------
    std::map<int,double> getResults();

    //! --------------------------
    //! return min and max values
    //! --------------------------
    std::pair<double, double> getMinMax() const;

    //! --------------
    //! divide domain
    //! --------------
    void splitSourceIntoBuckets();

    //! set remap flag
    void setRemap(bool remapFlag) {myRemapFlag = remapFlag;}

    //! set the number of remapping steps
    void setRemappingSteps(int rs){ myRemappingSteps = rs; }
/*
    //! --------------------------------------
    //! getMultiresult for multiInterpolation
    //! --------------------------------------
    std::map<int, std::vector<double> > getMultiResults();
*/
    //! ---------------------------
    //! set the progress indicator
    //! ---------------------------
    void setProgressIndicator(QProgressIndicator *aProgressIndicator);

private:

    //! -------------------
    //! progress indicator
    //! -------------------
    QProgressIndicator *myProgressIndicator;

    //! ----------------
    //! the target mesh
    //! ----------------
    occHandle(MeshVS_DataSource) myTargetMesh;

    //! -------------------------------------------
    //! check if a point is inside an element
    //! if Y return the interpolation coefficients
    //! -------------------------------------------
    bool isInner(const mesh::meshElement &anElement, const nodeSource &aSourceNode, std::vector<double> &IC);

    //! -------------
    //! source nodes
    //! -------------
    std::vector<nodeSource> mySourceNodes;
    
    //! ---------------------
    //! target mesh elements
    //! ---------------------
    std::vector<mesh::meshElement> myVecElements;
    std::vector<mesh::meshElement> myVecEmptyElements;

    std::map<int,std::vector<double>> myMapOfNodes;         //! target mesh nodes - key->nodeID, value[3]->coords(x, y, z)
    std::map<int,double> myRes;                             //! target mesh nodes - key->nodeID, averaged interpolated scalar value

    std::multimap<int,std::vector<double>> myMultiResNodes; //! target mesh nodes - intermediate result of the multiInterpolation
    std::map<int,std::vector<double>> myMultiRes;           //! target mesh nodes - key->nodeID, averaged interpolated scalar value for multiinterpolate

    //! -------------------------
    //! grouping of source nodes
    //! -------------------------
    bucket myBuckets;
    int myNbucketsX;
    int myNbucketsY;
    int myNbucketsZ;

    //! ---------------------
    //! remapping parameters
    //! ---------------------
    bool myRemapFlag;
    int myRemappingSteps;

    //! ------------------
    //! mesh bounding box
    //! ------------------
    double myXmin, myYmin, myZmin, myXmax, myYmax, myZmax;

    //! -------------------
    //! max and min values
    //! -------------------
    double myMinValue;
    double myMaxValue;

    //! --------------------
    //! the biggest element
    //! --------------------
    double theBiggestElement();

    //! clock
    static void clock();

public:

    //! ------------------------------
    //! perform using shape functions
    //! ------------------------------
    void performShapeFunctions();

    //! ----------------------------------------
    //! perform using "in pinball" source point
    //! ----------------------------------------
    void performNearest(double pinball = 100);
    void performNearestNeighboring(double pinball = 100);

public slots:

    //! -----
    //! slot
    //! -----
    void perform(int theAlgo, const occHandle(MeshVS_DataSource) &theTargetMesh = occHandle(MeshVS_DataSource)());

signals:

    //! -----------------------------
    //! emit interpolation finisched
    //! -----------------------------
    void interpolationFinished();
};

#endif // MAPPER3DCLASS_H
