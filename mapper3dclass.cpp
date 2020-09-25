#define OVERLAP_BUCKETS
#define CORRECT_BUCKETS

//! ---
//! Qt
//! ---
#include <QLabel>
#include <QTime>
#include <QDir>
#include <QMessageBox>
#include <QApplication>

//! ----------------
//! custom includes
//! ----------------
#include "mapper3dclass.h"
#include "mathtools.h"
#include "mydefines.h"
#include "qprogressevent.h"
#include "tools.h"
#include "qprogressindicator.h"
#include "tools.h"

//! ----
//! OCC
//! ----
#include <MeshVS_DataSource.hxx>
#include <TColStd_PackedMapOfInteger.hxx>
#include <TColStd_MapIteratorOfPackedMapOfInteger.hxx>
#include <TColStd_Array1OfInteger.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <Bnd_Box.hxx>

//! ----
//! C++
//! ----
#include <algorithm>

//! -------
//! system
//! -------
#include <Windows.h>

//! ---------------------
//! function: destructor
//! details:
//! ---------------------
Mapper3DClass::~Mapper3DClass()
{
    cout<<"Mapper3DClass::~Mapper3DClass()->____DESTRUCTOR CALLED____"<<endl;
    myBuckets.deallocate();
}

//! -------------------------
//! function: constructor I
//! details:
//! -------------------------
Mapper3DClass::Mapper3DClass(const opencascade::handle<MeshVS_DataSource> &theTargetMeshDS,QObject *parent):
    myTargetMesh(theTargetMeshDS),
    QObject(parent)
{
    cout<<"Mapped3DClass::Mapper3DClass->____constructor I called____"<<endl;
    myNbucketsX = 1;
    myNbucketsY = 1;
    myNbucketsZ = 1;
    myRemapFlag = false;
    myRemappingSteps = 1;
    myMinValue = 1e80;
    myMaxValue = -1e80;
    myProgressIndicator = Q_NULLPTR;

    this->setTargetMesh(theTargetMeshDS);

    //cout<<"Mapper3DClass->____target mesh bounding box size: ("<<myXmax-myXmin<<", "<<myYmax-myYmin<<", "<<myZmax-myZmin<<")____"<<endl;
}

//! -------------------------
//! function: constructor II
//! details:
//! -------------------------
Mapper3DClass::Mapper3DClass(QObject *parent):QObject(parent)
{
    cout<<"Mapped3DClass::Mapper3DClass->____constructor I called____"<<endl;
    myNbucketsX = 1;
    myNbucketsY = 1;
    myNbucketsZ = 1;
    myRemapFlag = false;
    myRemappingSteps = 1;
    myMinValue = 1e80;
    myMaxValue = -1e80;
    myProgressIndicator = Q_NULLPTR;
}

//! ------------------------
//! function: setTargetMesh
//! details:
//! ------------------------
void Mapper3DClass::setTargetMesh(const opencascade::handle<MeshVS_DataSource> &theTargetMeshDS)
{
    //! ---------------------------------------------
    //! initialize with OCC mesh
    //! prepare the elements as std::vector<element>
    //! element {int ID, vecNodes[], enum Type}
    //! ---------------------------------------------
    TColStd_PackedMapOfInteger eMap = theTargetMeshDS->GetAllElements();
    TColStd_MapIteratorOfPackedMapOfInteger eIter;
    int buf[20];
    for(eIter.Initialize(eMap);eIter.More();eIter.Next())
    {
        int elementID = eIter.Key();

        int NbNodes;
        TColStd_Array1OfInteger nodeIDs(*buf,1,20);
        theTargetMeshDS->GetNodesByElement(elementID,nodeIDs,NbNodes);

        mesh::meshElement anElem;
        anElem.ID=elementID;
        for(int k=0;k<NbNodes;k++)
        {
            int nodeID = nodeIDs.Value(k+1);
            anElem.theNodeIDs.push_back(nodeID);
        }
        switch(NbNodes)
        {
        case 4: anElem.type = TET; break;
        case 10: anElem.type = TET10; break;
        case 8: anElem.type = HEXA; break;
        case 20: anElem.type = HEXA20; break;
        }
        myVecElements.push_back(anElem);
        myVecEmptyElements.push_back(anElem);
    }
    //cout<<"Mapper3DClass::setTargetMesh()->____"<<myVecElements.size()<<" elements prepared"<<endl;

    //! -------------------------------------------------------
    //! prepare the nodes as std::map<int,std::vector<double>>
    //! <int> => nodeID, double[3] => spatial coords
    //! -------------------------------------------------------
    //cout<<"Mapper3DClass::setTargetMesh()->____preparing nodes____"<<endl;

    double bufn[3];
    for(TColStd_MapIteratorOfPackedMapOfInteger nIter(theTargetMeshDS->GetAllNodes());nIter.More();nIter.Next())
    {
        int globalNodeID = nIter.Key();
        int NbNodes;
        MeshVS_EntityType type;
        TColStd_Array1OfReal Coords(*bufn,1,3);
        theTargetMeshDS->GetGeom(globalNodeID,false,Coords,NbNodes,type);
        std::vector<double> coor {Coords.Value(1), Coords.Value(2), Coords.Value(3)};
        //myMapOfNodes.insert(globalNodeID,coor);
        myMapOfNodes.insert(std::make_pair(globalNodeID,coor));
    }

    //! -----------------------------
    //! the target mesh bounding box
    //! increase the size a bit
    //! -----------------------------
    theTargetMeshDS->GetBoundingBox().Get(myXmin,myYmin,myZmin,myXmax,myYmax,myZmax);
    double Lx = myXmax-myXmin;
    double Ly = myYmax-myYmin;
    double Lz = myZmax-myZmin;

    myXmin = myXmin - Lx*0.10;
    myYmin = myYmin - Ly*0.10;
    myZmin = myZmin - Lz*0.10;
    myXmax = myXmax + Lx*0.10;
    myYmax = myYmax + Ly*0.10;
    myZmax = myZmax + Lz*0.10;

    //cout<<"Mapper3DClass::setTargetMesh()->____target mesh bounding box size: ("<<myXmax-myXmin<<", "<<myYmax-myYmin<<", "<<myZmax-myZmin<<")____"<<endl;
}

//! -----------------------
//! function: setNbBuckets
//! details:
//! -----------------------
void Mapper3DClass::setNbBuckets(int NbBucketsX, int NbBucketsY, int NbBucketsZ)
{
    myNbucketsX = NbBucketsX;
    myNbucketsY = NbBucketsY;
    myNbucketsZ = NbBucketsZ;
}

//! ----------------------------
//! function: theBiggestElement
//! details:
//! ----------------------------
double Mapper3DClass::theBiggestElement()
{
    double Vmax = -1e80;
    for(std::vector<mesh::meshElement>::iterator it = myVecElements.begin() ;it!=myVecElements.end(); ++it)
    {
        const mesh::meshElement &curElem = *it;
        switch(curElem.type)
        {
        case TET:
        {
            const int n = 4;
            MatDoub a(n,n);
            for(int row=0; row<n; row++)
            {
                int nodeID = curElem.theNodeIDs.at(row);
                //const std::vector<double> &x = myMapOfNodes.value(nodeID);
                const std::vector<double> &x = myMapOfNodes.at(nodeID);
                for(int col=0; col<n-1; col++) a[row][col]=x[col];
            }
            for(int row=0; row<n; row++) a[row][n-1]=1.0;
            LUdcmp alu(a);
            double V = alu.det();
            if(V>Vmax) Vmax = V;
        }
            break;
            //! --------------------
            //! other element types
            //! --------------------
        }
    }
    return Vmax;
}

//! -------------------------------------------------------
//! function: isInner
//! details: if a source node is inside the target element
//! -------------------------------------------------------
#include <polyhedron.h>
bool Mapper3DClass::isInner(const mesh::meshElement &anElement, const nodeSource &aSourceNode, std::vector<double> &IC)
{
    bool containsNode = false;
    switch(anElement.type)
    {
    case TET:
    {
        //cout<<"Mapper3DClass::isInner()->____TET found____"<<endl;
        const Int n = 4;

        //! -----------------------------------------------------------------
        //! define and fill the matrix for the element volume calculation
        //! init the matrix
        //!
        //!     | x1 y1 z1 1 |                  | x1 y1 z1 1 |
        //!     | x2 y2 z2 1 |                  | x  y  z  1 |
        //!  V= | x3 y3 z3 1 |      V(voln=1) = | x3 y3 z3 1 |  [...]
        //!     | x4 y4 z4 1 |                  | x4 y4 z4 1 |
        //!
        //! voln => sub volume number
        //!
        //! If voln = 0 the matrix is initialized to return the volume of V0
        //! If voln = 1 the matrix is initialized to return the volume of V1
        //! -----------------------------------------------------------------

        MatDoub a(n,n);
        for(int row=0; row<n; row++)
        {
            int nodeID = anElement.theNodeIDs.at(row);
            //const std::vector<double> &x = myMapOfNodes.value(nodeID);
            const std::vector<double> &x = myMapOfNodes.at(nodeID);
            for(int col=0; col<n-1; col++)
            {
                double val = x[col];
                a[row][col]=val;
            }
        }
        //! fill the last column with "1"
        for(int row=0; row<n; row++) a[row][n-1]=1.0;

        //! calculate the volume V of the element
        LUdcmp alu(a);
        double V = alu.det();
        if(V<=0) return false;

        //! define the matrix for the sub-volumes calculation
        MatDoub a_(n,n);
        double D[4];

        //! calculate the sub-volumes
        for(int voln=0; voln<n; voln++)
        {
            //! fill one matrix
            for(int row=0; row<n; row++)
            {
                int nodeID = anElement.theNodeIDs.at(row);
                //const std::vector<double> &x = myMapOfNodes.value(nodeID);
                const std::vector<double> &x = myMapOfNodes.at(nodeID);
                for(int col=0; col<n-1; col++)
                {
                    double val;
                    if(row==voln) val = aSourceNode.x[col];
                    else val = x[col];
                    a_[row][col]=val;
                }
            }
            //! fill the last column with "1"
            for(int row=0; row<n; row++) a_[row][n-1]=1.0;

            //! calculate the natural coordinates (normalized sub-volume in case of a TET)
            LUdcmp alu_(a_);
            D[voln] = alu_.det()/V;
            //cout<<"Mapper3DClass::isInner->____element sub volume ("<<voln<<"): "<<V<<"____"<<endl;
            //cout<<"Mapper3DClass::isInner->____volume ratio for subvol "<<voln<<" = "<<D[voln]<<"____"<<endl;

            if(D[voln]<-gp::Resolution() || D[voln]>1.0+gp::Resolution())
            {
                containsNode = false;
                IC.push_back(V);
                break;
            }
            else
            {
                IC.push_back(D[voln]);
                containsNode=true;
            }
        }
    }
        break;

    case PYRAM:
    {
        /*
        //cout<<"Mapper3DClass::isInner()->____PYRAM found____"<<endl;

        //! ---------------------
        //! volume of th element
        //! ---------------------
        std::vector<int> nodeIDs = anElement.theNodeIDs;

        int apexIndex = nodeIDs.at(0);
        double xV = myMapOfNodes.value(apexIndex).at(0);
        double yV = myMapOfNodes.value(apexIndex).at(1);
        double zV = myMapOfNodes.value(apexIndex).at(2);
        polygon::Point P0(xV,yV,zV);

        int index1 = nodeIDs.at(1);
        double x1 = myMapOfNodes.value(index1).at(0);
        double y1 = myMapOfNodes.value(index1).at(1);
        double z1 = myMapOfNodes.value(index1).at(2);
        polygon::Point P1(x1,y1,z1);

        int index2 = nodeIDs.at(2);
        double x2 = myMapOfNodes.value(index2).at(0);
        double y2 = myMapOfNodes.value(index2).at(1);
        double z2 = myMapOfNodes.value(index2).at(2);
        polygon::Point P2(x2,y2,z2);

        int index3 = nodeIDs.at(3);
        double x3 = myMapOfNodes.value(index3).at(0);
        double y3 = myMapOfNodes.value(index3).at(1);
        double z3 = myMapOfNodes.value(index3).at(2);
        polygon::Point P3(x3,y3,z3);

        int index4 = nodeIDs.at(4);
        double x4 = myMapOfNodes.value(index4).at(0);
        double y4 = myMapOfNodes.value(index4).at(1);
        double z4 = myMapOfNodes.value(index4).at(2);
        polygon::Point P4(x4,y4,z4);

        //! ------------------------------
        //! volume of the pyramid element
        //! ------------------------------
        std::vector<polygon::Point> faceQuad,faceTrig1,faceTrig2,faceTrig3,faceTrig4;
        faceQuad.push_back(P4); faceQuad.push_back(P3); faceQuad.push_back(P2); faceQuad.push_back(P1);
        faceTrig1.push_back(P0); faceTrig1.push_back(P1); faceTrig1.push_back(P2);
        faceTrig2.push_back(P0); faceTrig2.push_back(P2); faceTrig2.push_back(P3);
        faceTrig3.push_back(P0); faceTrig3.push_back(P3); faceTrig3.push_back(P4);
        faceTrig3.push_back(P0); faceTrig3.push_back(P4); faceTrig3.push_back(P1);

        std::vector<std::vector<polygon::Point>> faces;
        faces.push_back(faceQuad);
        faces.push_back(faceTrig1);
        faces.push_back(faceTrig2);
        faces.push_back(faceTrig3);
        faces.push_back(faceTrig4);
        polyhedron elementPoly(faces);

        //! ----------------------------------
        //! the volume of the pyramid element
        //! ----------------------------------
        double elementVol = elementPoly.volume();

        //! ---------------------------------
        //! the source point to be evaluated
        //! ---------------------------------
        double xP = aSourceNode.x[0];
        double yP = aSourceNode.x[1];
        double zP = aSourceNode.x[2];
        std::vector<polygon::Point> P(xP,yP,zP);

        //! ------------------------------------------------------------------
        //! volume of a pyramid having apex in "P"
        //! use P,1,2,3,4
        //!    [0,1,2,3,4]
        //! remember the definition {4,3,2,1} {0,2,1} {0,2,3} {0,3,4} {0,4,1}
        //! ------------------------------------------------------------------
        std::vector<polygon::Point> faceQuad_sub;
        std::vector<polygon::Point> faceTrig1_sub, faceTrig2_sub,faceTrig3_sub, faceTrig4_sub;
        faceQuad_sub.push_back(P4); faceQuad_sub.push_back(P3); faceQuad_sub.push_back(P2); faceQuad_sub.push_back(P1);
        faceTrig1_sub.push_back(P); faceTrig1_sub.push_back(2); faceTrig1_sub.push_back(1);
        faceTrig2_sub.push_back(P); faceTrig2_sub.push_back(2); faceTrig2_sub.push_back(3);
        faceTrig3_sub.push_back(P); faceTrig3_sub.push_back(3); faceTrig2_sub.push_back(4);
        faceTrig4_sub.push_back(P); faceTrig4_sub.push_back(4); faceTrig4_sub.push_back(1);

        std::vector<std::vector<polygon::Point>> pyrFaces;
        pyrFaces<<faceQuad_sub;
        pyrFaces<<faceTrig1_sub;
        pyrFaces.push_back(faceTrig2_sub);
        pyrFaces.push_back(faceTrig3_sub);
        pyrFaces.push_back(faceTrig4_sub);
        polyhedron aPyr(pyrFaces);
        double vol0 = elementVol-aPyr.volume();

        //! ------------
        //! use P,0,2,3
        //! ------------
        double vol1 = elementVol;
        std::vector<polygon::Point> tet1;
        //polyhedron::tetrahedronVolume();

        //! ------------
        //! use P,1,2,3
        //! ------------
        std::vector<polygon::Point> tet2;
        double vol2 = elementVol;

        //! ------------
        //! use P,0,1,4
        //! ------------
        std::vector<polygon::Point> tet3;
        double vol3 = elementVol;

        */
    }
        break;

    case PRISM:
    {
        //cout<<"Mapper3DClass::isInner()->____PRISM found____"<<endl;
    }
        break;
    }
    return containsNode;
}

//! -----------------------------------------------
//! function: splitSourceIntoBuckets
//! details:  split the source points into buckets
//! -----------------------------------------------
void Mapper3DClass::splitSourceIntoBuckets()
{
    //cout<<"Mapper3DClass::splitSourceIntoBuckets()->____function called____"<<endl;

#ifdef CORRECT_BUCKETS
    //! --------------------------------
    //! first check the biggest element
    //! --------------------------------
    double maxVolume = this->theBiggestElement();
    if(maxVolume!=0)
    {
        int recalc_NbBucketsX, recalc_NbBucketsY, recalc_NbBucketsZ;
        double typicalSize = 1.0*pow(maxVolume,1.0/3.0);
        //cout<<"Mapper3DClass::splitSourceIntoBuckets()->____maximum linear element size: "<<typicalSize<<"____"<<endl;

        //! --------------
        //! check along x
        //! --------------
        recalc_NbBucketsX = floor((myXmax-myXmin)/typicalSize);
        if(recalc_NbBucketsX == 0)  recalc_NbBucketsX++;
        if(myNbucketsX>recalc_NbBucketsX)
        {
            myNbucketsX = recalc_NbBucketsX;
        }
        //! --------------
        //! check along y
        //! --------------
        recalc_NbBucketsY = floor((myYmax-myYmin)/typicalSize);
        if(recalc_NbBucketsY == 0)  recalc_NbBucketsY++;
        if(myNbucketsY>recalc_NbBucketsY)
        {
            myNbucketsY = recalc_NbBucketsY;
        }
        //! --------------
        //! check along z
        //! --------------
        recalc_NbBucketsZ = floor((myZmax-myZmin)/typicalSize);
        if(recalc_NbBucketsZ == 0)  recalc_NbBucketsZ++;
        if(myNbucketsZ>recalc_NbBucketsZ)
        {
            myNbucketsZ = recalc_NbBucketsZ;
        }
        //cout<<"Mapper3DClass::splitSourceIntoBuckets()->____("<<myNbucketsX<<", "<<myNbucketsY<<", "<<myNbucketsZ<<")____"<<endl;
    }
#endif

    //! -----------
    //! allocation
    //! -----------
    myBuckets.allocate(myNbucketsX,myNbucketsY,myNbucketsZ);

    //! ---------------------
    //! width of the buckets
    //! ---------------------
    double deltaX = (myXmax-myXmin)/myNbucketsX;
    double deltaY = (myYmax-myYmin)/myNbucketsY;
    double deltaZ = (myZmax-myZmin)/myNbucketsZ;

    //! --------------
    //! bucket number
    //! --------------
    int bucketX, bucketY, bucketZ;

    int npointsoutsidebuckets = 0;
    int npointinsidebuckets = 0;
    for(std::vector<Mapper3DClass::nodeSource>::iterator it=mySourceNodes.begin(); it!=mySourceNodes.end(); ++it)
    {
        const Mapper3DClass::nodeSource &curSN = *it;
        double xsource = curSN.x[0];
        double ysource = curSN.x[1];
        double zsource = curSN.x[2];

        //! ----------------
        //! find the bucket
        //! ----------------
        bucketX = int((xsource-myXmin)/deltaX);
        bucketY = int((ysource-myYmin)/deltaY);
        bucketZ = int((zsource-myZmin)/deltaZ);

        if(bucketX==myNbucketsX) bucketX--;     //! point on the right x boundary of the bucket
        if(bucketY==myNbucketsY) bucketY--;     //! point on the right y boundary of the bucket
        if(bucketZ==myNbucketsZ) bucketZ--;     //! point on the right z boundary of the bucket

        if(bucketX<myNbucketsX && bucketY<myNbucketsY && bucketZ<myNbucketsZ && bucketX>=0 && bucketY>=0 && bucketZ>=0)
        {
            myBuckets.setValue(bucketX,bucketY,bucketZ,curSN);
            npointinsidebuckets++;
        }
#ifndef OVERLAP_BUCKETS
        else
        {
            //cout<<"____the point is outside the bucket____"<<endl;
            npointsoutsidebuckets++;
        }
#endif
#ifdef OVERLAP_BUCKETS
        //! ---------------------------------------
        //! add some points to the current buckets
        //! using left and right data
        //! ---------------------------------------
        double k = 0.25;
        double halfWidthX = deltaX*k;
        double halfWidthY = deltaY*k;
        double halfWidthZ = deltaZ*k;
        double x_min = bucketX*deltaX;
        double y_min = bucketY*deltaY;
        double z_min = bucketZ*deltaZ;
        double x_max = (bucketX+1)*deltaX+halfWidthX;
        double y_max = (bucketY+1)*deltaY+halfWidthY;
        double z_max = (bucketZ+1)*deltaZ+halfWidthZ;
        if(xsource>=x_min && xsource<=x_max && ysource>=y_min && ysource<=y_max && zsource>=z_min && zsource<=z_max)
        {
            myBuckets.setValue(bucketX,bucketY,bucketZ,curSN);
            npointinsidebuckets++;
        }
        else
        {
            //cout<<"____the point is outside the bucket____"<<endl;
            npointsoutsidebuckets++;
        }
#endif
    }
    //cout<<"Mapper3DClass::splitSourceIntoBuckets()->____total number of source nodes: "<<mySourceNodes.size()<<"____"<<endl;
    //cout<<"Mapper3DClass::splitSourceIntoBuckets()->____total number of source nodes into buckects: "<<npointinsidebuckets<<"____"<<endl;
    //cout<<"Mapper3DClass::splitSourceIntoBuckets()->____total number of source nodes outside buckets: "<<npointsoutsidebuckets<<"____"<<endl;
}

//! ---------------------
//! function: set source
//! details:
//! ---------------------
void Mapper3DClass::setSource(const std::vector<nodeSource> &vecSourceNodes)
{
    mySourceNodes = vecSourceNodes;
}

//!--------------------------------
//! function: setProgressIndicator
//! details:
//!--------------------------------
void Mapper3DClass::setProgressIndicator(QProgressIndicator *aProgressIndicator)
{
    myProgressIndicator = aProgressIndicator;
}

//! --------------------------------
//! function: performShapeFunctions
//! details:
//! --------------------------------
void Mapper3DClass::performShapeFunctions()
{
    cout<<"Mapper3DClass::interpolate()->____function called____"<<endl;

    //! -----------------
    //! a progress event
    //! -----------------
    QProgressEvent *aPrgEvvent;

    if(myProgressIndicator!=Q_NULLPTR)
    {
        //! ----------------------------
        //! init the progress indicator
        //! ----------------------------
        QString msg("Start interpolating");
        int NbEvents = int(myVecElements.size());
        aPrgEvvent = new QProgressEvent(QProgressEvent_None,0,9999,0,msg,QProgressEvent_Init,0,NbEvents-1,0);
        QApplication::postEvent(myProgressIndicator,aPrgEvvent);
        QApplication::processEvents();
        Sleep(500);
    }

    //! -----------------------------------
    //! calculate the width of the buckets
    //! -----------------------------------
    //cout<<"Mapper3DClass::interpolate()->____Number of buckets: ("<<myNbucketsX<<", "<<myNbucketsY<<", "<<myNbucketsZ<<")____"<<endl;
    double dx = (myXmax-myXmin)/myNbucketsX;
    double dy = (myYmax-myYmin)/myNbucketsY;
    double dz = (myZmax-myZmin)/myNbucketsZ;
    //cout<<"Mapper3DClass::interpolate()->____Buckets width ("<<dx<<", "<<dy<<", "<<dz<<")____"<<endl;

    //! --------------------------------------
    //! iterate over the target mesh elements
    //! --------------------------------------
    int en=1;
    int nn=0;
    for(std::vector<mesh::meshElement>::iterator it = myVecElements.begin() ;it!=myVecElements.end(); ++it, en++)
    {        
        const mesh::meshElement &anElement = *it;

        //! ------------------------------------------------
        //! center of the element: once the center is known
        //! it is possible to get its buckets
        //! ------------------------------------------------
        double xc, yc, zc;
        if(anElement.type==TET)
        {
            //! ------------------------------------
            //! calculate the center of the element
            //! ------------------------------------
            xc= yc= zc = 0.0;
            const std::vector<int> &nodeIDs = anElement.theNodeIDs;
            for(std::vector<int>::const_iterator aNodeIDIter = nodeIDs.cbegin(); aNodeIDIter!=nodeIDs.cend(); ++aNodeIDIter)
            {
                int nodeID = *aNodeIDIter;
                //xc = xc + myMapOfNodes.value(nodeID).at(0);
                //yc = yc + myMapOfNodes.value(nodeID).at(1);
                //zc = zc + myMapOfNodes.value(nodeID).at(2);
                xc += myMapOfNodes.at(nodeID)[0];
                yc += myMapOfNodes.at(nodeID)[1];
                zc += myMapOfNodes.at(nodeID)[2];
            }
            xc = 0.25*xc;
            yc = 0.25*yc;
            zc = 0.25*zc;
        }

        //! ------------------------------------------------
        //! extract the source nodes of the selected bucket
        //! ------------------------------------------------
        int bucketX_targetElement=floor((xc-myXmin)/dx);
        int bucketY_targetElement=floor((yc-myYmin)/dy);
        int bucketZ_targetElement=floor((zc-myZmin)/dz);

        if(bucketX_targetElement==myNbucketsX) bucketX_targetElement--;
        if(bucketY_targetElement==myNbucketsY) bucketY_targetElement--;
        if(bucketZ_targetElement==myNbucketsZ) bucketZ_targetElement--;

        std::vector<nodeSource> vecSourceNodesFiltered = myBuckets.value(bucketX_targetElement,bucketY_targetElement,bucketZ_targetElement);

        //! ------------------------------
        //! iterate over the source nodes
        //! ------------------------------
        bool isEmpty = true;

        for(std::vector<nodeSource>::iterator itSourceNodes = vecSourceNodesFiltered.begin();itSourceNodes!=vecSourceNodesFiltered.end();)
        {
            const nodeSource &ns = *itSourceNodes;
            std::vector<double> IC;
            bool isInner = this->isInner(anElement,ns,IC);
            if(isInner==true)
            {
                isEmpty = false;
                nn++;
                std::vector<double> val = ns.vecVal;
                for(int i=0;i<IC.size();i++)
                {
                    int nodeID = anElement.theNodeIDs[i];
                    //myMultiResNodes.insert(nodeID,val);
                    myMultiResNodes.insert(std::make_pair(nodeID,val));
                }
                //! --------------------------------
                //! remove the node while iterating
                //! --------------------------------
                itSourceNodes = vecSourceNodesFiltered.erase(itSourceNodes);

                //! -------------------------------------------------------------------
                //! erase the current element from the vector of empty target elements
                //! -------------------------------------------------------------------
                myVecEmptyElements.erase(std::remove(myVecEmptyElements.begin(), myVecEmptyElements.end(), anElement), myVecEmptyElements.end());
            }
            else
            {
                ++itSourceNodes;
            }
        }
        if(en%250==0)
        {
            if(myProgressIndicator!=Q_NULLPTR)
            {
                //! --------------------
                //! update the progress
                //! --------------------
                QString msg = QString("Working on element: %1").arg(en);
                aPrgEvvent = new QProgressEvent(QProgressEvent_Update,0,100,1,msg,QProgressEvent_Update,0,9999,en);
                QApplication::postEvent(myProgressIndicator,aPrgEvvent);
                QApplication::processEvents();
            }
        }
    }

    //cout<<"Mapper3DClass::interpolate()->____number of empty elements: "<<myVecEmptyElements.size()<<"____"<<endl;

    //! --------------------------
    //! try to remap if requested
    //! --------------------------
    int nr = 0;
    if(myRemapFlag==true) nr=this->remapByTargetElements();

    //! --------------------------------------------
    //! put scalar on target node: "0" uses average
    //! --------------------------------------------
    this->putScalarOnTargetNode(0);

    //! -----------------------------------------
    //! notify end of interpolation
    //! the progress bar is closed by the caller
    //! -----------------------------------------
    emit interpolationFinished();
}

//! --------------------------------------------
//! function: putScalarOnTargetNode
//! details:  interpolation type: 0 => average
//!                               1 => one node
//! --------------------------------------------
void Mapper3DClass::putScalarOnTargetNode(int interpolationTypeFunction)
{
    cout<<"Mapper3DClass::putScalarOnTargetNode()->____function called____"<<endl;
    //cout<<"Mapper3DClass::putScalarOnTargetNode()->____myMultiresNode size "<<myMultiResNodes.size()<<"____"<<endl;

    for(std::multimap<int,std::vector<double>>::iterator it = myMultiResNodes.begin(); it!=myMultiResNodes.end(); ++it)
    //for(QMultiMap<int,std::vector<double>>::iterator it = myMultiResNodes.begin(); it!=myMultiResNodes.end(); ++it)
    {
        std::vector<double> vecVal;               //! list of values at node at time
        //const QList<std::vector<double>> &values = myMultiResNodes.values(it.key());
        QList<std::vector<double>> values;
        int curKey = it->first;
        auto itr1 = myMultiResNodes.lower_bound(curKey);
        auto itr2 = myMultiResNodes.upper_bound(curKey);
        while (itr1 != itr2)
        {
            if (itr1 -> first == curKey) values.push_back(itr1->second);
            itr1++;
        }

        int sizeList = values.length();     //! number of times
        int m = values.at(0).size();        //!

        switch(interpolationTypeFunction)
        {
        case 0:
        {
            //! -------------------
            //! arithmetic average
            //! -------------------
            for(int i=0;i<m;i++)
            {
                double S=0;
                for(QList<std::vector<double>>::const_iterator li=values.cbegin();li!=values.cend();++li)
                {
                    std::vector<double> value=*li;
                    S+=value.at(i);
                }
                vecVal.push_back(S/sizeList);
            }

            /*
            for(int timeIndex = 0; timeIndex<sizeList; timeIndex++)
            {
                const QList<double> &valuesAtNodeAtTime = values.at(timeIndex);
                double S = 0;
                for(int j = 0; j<valuesAtNodeAtTime.length(); j++)
                {
                    S = S+valuesAtNodeAtTime.at(j);
                }
                S = S/valuesAtNodeAtTime.length();
                vecVal<<S;
            }
            */
        }
            break;

        case 1:
        {
            //! ---------------
            //! a single value
            //! ---------------
            vecVal = values.at(0);
        }
            break;
        }
        //myMultiRes.insert(it.key(),vecVal);
        myMultiRes.insert(std::make_pair(it->first,vecVal));
    }
}

//! --------------------
//! function: getMinMax
//! details:
//! --------------------
std::pair<double,double> Mapper3DClass::getMinMax() const
{
    std::pair<double,double> pair;
    pair.first=myMinValue;
    pair.second=myMaxValue;
    return pair;
}

//! --------------------------
//! function: clock
//! details:  time diagnostic
//! --------------------------
void Mapper3DClass::clock()
{
    QTime time = QTime::currentTime();
    QString text = time.toString("mm:ss:zzz");
    cout<<text.toStdString()<<endl;
}

//! -----------------------------------------------------------
//! function: remapByTargetElements
//! details:  scan the list of the empty elements and map them
//! -----------------------------------------------------------
int Mapper3DClass::remapByTargetElements()
{
    cout<<"Mapper3DClass::remapByTargetElements->____function called. Number of empty elements: "<<myVecEmptyElements.size()<<"____"<<endl;

    int nremapped = 0;
    int tryCount = 1;
    double pinball = 0;
    for(std::vector<mesh::meshElement>::iterator itElem = myVecEmptyElements.begin(); itElem!=myVecEmptyElements.end(); ++itElem)
    {
        const mesh::meshElement &curElem = *itElem;        
        const std::vector<int> &nodeIDs = curElem.theNodeIDs;

        //! ----------------------------------------------
        //! the typical size of the unmapped element
        //! (element which does not contain source nodes)
        //! ----------------------------------------------
        double L=0;

        //! ----------------------
        //! center of the element
        //! ----------------------
        double xc = 0, yc = 0, zc = 0;

        if(curElem.type==TET)
        {
            for(std::vector<int>::const_iterator itIDs = nodeIDs.cbegin(); itIDs!= nodeIDs.cend(); ++itIDs)
            {
                int nodeID = *itIDs;
                xc += xc + myMapOfNodes.at(nodeID)[0];
                yc += yc + myMapOfNodes.at(nodeID)[1];
                zc += zc + myMapOfNodes.at(nodeID)[2];
            }
            xc = xc*0.25;
            yc = yc*0.25;
            zc = zc*0.25;

            //! --------------------------------------------------------------------
            //! calculate the volume of the element, and from that, its typcal size
            //! --------------------------------------------------------------------
            const int n = 4;
            MatDoub a(n,n);
            for(int row=0; row<n; row++)
            {
                int nodeID = curElem.theNodeIDs.at(row);
                const std::vector<double> &x = myMapOfNodes.at(nodeID);
                for(int col=0; col<n-1; col++) a[row][col]=x[col];
            }
            for(int row=0; row<n; row++) a[row][n-1]=1.0;
            LUdcmp alu(a);
            double V = alu.det();
            L = pow(V,1.0/3.0);
        }

        //! --------------------------------------------------------
        //! do not use all the source nodes for remapping, but only
        //! the ones contained in the element bucket
        //! --------------------------------------------------------
        int bucketX_targetElement=floor((xc-myXmin)/(myXmax-myXmin));
        int bucketY_targetElement=floor((yc-myYmin)/(myYmax-myYmin));
        int bucketZ_targetElement=floor((zc-myZmin)/(myZmax-myZmin));

        if(bucketX_targetElement==myNbucketsX) bucketX_targetElement--;
        if(bucketY_targetElement==myNbucketsY) bucketY_targetElement--;
        if(bucketZ_targetElement==myNbucketsZ) bucketZ_targetElement--;
        std::vector<Mapper3DClass::nodeSource> vecSourceNodesFiltered = myBuckets.value(bucketX_targetElement,bucketY_targetElement,bucketZ_targetElement);

        //! ---------------------------------------------------
        //! once the center is known scan all the source nodes
        //! init the pinball value
        //! ---------------------------------------------------
        if(tryCount==1) pinball = 1.0*L;
        else pinball = pinball*1.1;

        bool hasBeenRemappedDone=false;
        for(std::vector<Mapper3DClass::nodeSource>::iterator itSourceNodes = vecSourceNodesFiltered.begin(); itSourceNodes!=vecSourceNodesFiltered.end();)
        {
            const Mapper3DClass::nodeSource &aSN = *itSourceNodes;

            //! -----------------------------------------------------------------------
            //! distance between the current source node and the center of the element
            //! -----------------------------------------------------------------------
            double d = sqrt(pow(xc-aSN.x[0],2)+pow(yc-aSN.x[1],2)+pow(zc-aSN.x[2],2));
            if(d<=pinball)
            {
                //! --------------------------------------
                //! the source node is inside the pinball
                //! --------------------------------------
                for(int i=0;i<curElem.theNodeIDs.size();i++)
                {
                    int nodeID = curElem.theNodeIDs.at(i);
                    myMultiResNodes.insert(std::make_pair(nodeID,aSN.vecVal));
                }
                itSourceNodes = vecSourceNodesFiltered.erase(itSourceNodes);
                nremapped++;
                hasBeenRemappedDone = true;
                break;
            }
            else
            {
                hasBeenRemappedDone = false;
                ++itSourceNodes;
            }
        }

        if(hasBeenRemappedDone == false)
        {
            if(tryCount<=myRemappingSteps)
            {
                tryCount++;
                itElem--;
            }
            else
            {
                //cout<<"Mapper3DClass::remapByTargetElements->____element: "<<curElem.ID<<" cannot be renapped in "<<myRemappingSteps<<" steps____"<<endl;
                tryCount = 1;
                pinball = 0;
                continue;
            }
        }
        else
        {
            //cout<<"Mapper3DClass::remapByTargetElements->____element: "<<curElem.ID<<" remapped using pinball: "<<pinball<<"____"<<endl;
            tryCount = 1;
            pinball = 0;
        }
    }

    //cout<<"Mapper3DClass::remapByTargetElements()->____succesfully remapped % "<<double(nremapped)/double(myVecEmptyElements.size())*100<<" elemenst____"<<endl;
    return nremapped;
}

/*
//! ---------------------------------------------------------
//! function: getMultiResults
//! details:  retrieve the multiRes map for multiInterpolate
//! added 13/03/19
//! ---------------------------------------------------------
std::map<int,std::vector<double>> Mapper3DClass::getMultiResults()
{
    return myMultiRes;
}*/

//! ---------------------
//! function: getResults
//! details:
//! ---------------------
std::map<int,double> Mapper3DClass::getResults()
{
    return myRes;
}

//! ------------------------
//! function: retriveResMap
//! details:
//! ------------------------
void Mapper3DClass::retrieveResMap(int pos)
{
    //! -------------------------------
    //! key => nodeID
    //! value => list of values @ node
    //! -------------------------------
    myRes.clear();
    for(std::map<int,std::vector<double>>::iterator it = myMultiRes.begin(); it!=myMultiRes.end(); ++it)
    {
        int nodeID = it->first;
        const std::vector<double> &vec = myMultiRes.at(nodeID);
        double scalarVal = vec[pos];
        if(scalarVal>=myMaxValue) myMaxValue=scalarVal;
        if(scalarVal<=myMinValue) myMinValue=scalarVal;
        myRes.insert(std::make_pair(nodeID,scalarVal));
    }
}

//! -------------------------
//! function: performNearest
//! details:
//! -------------------------
void Mapper3DClass::performNearest(double pinball)
{
    cout<<"Mapper3DClass::performNearest()->____function called____"<<endl;

    //! -------------------
    //! read some settings
    //! -------------------
    int NbMappingSteps = 1;
    if(myRemapFlag == true) NbMappingSteps = myRemappingSteps+1;
    else NbMappingSteps = 1;

    //! -------------------------
    //! initialize with OCC mesh
    //! -------------------------
    TColStd_PackedMapOfInteger nMap = myTargetMesh->GetAllNodes();

    //! --------------
    //! send progress
    //! --------------
    if(myProgressIndicator!=Q_NULLPTR)
    {
        myProgressIndicator->setSecondaryBarVisible(true);
        QProgressEvent *aPrgEvent = new QProgressEvent(QProgressEvent_None,0,9999,0,"Start intepolating",QProgressEvent_Init,0,nMap.Extent()-1,0);
        QApplication::postEvent(myProgressIndicator,aPrgEvent);
        QApplication::processEvents();
        Sleep(500);
    }

    //! --------------------------
    //! progress bar increment 5%
    //! --------------------------
    int aProgressStep = int(nMap.Extent()*0.05);

    //! -----------------
    //! the target nodes
    //! -----------------
    int NbMappedNodes = 0;
    MeshVS_EntityType type;
    int NbNodes;
    double buf[3];
    int n=0;
    for(TColStd_MapIteratorOfPackedMapOfInteger it(nMap);it.More();it.Next(), n++)
    {
        int nodeID = it.Key();
        TColStd_Array1OfReal coords(*buf,1,3);
        myTargetMesh->GetGeom(nodeID,false,coords,NbNodes,type);
        targetMeshNode aTargetMeshNode(coords(1),coords(2),coords(3),nodeID);

        bool success = false;
        double localPinball = 0;
        for(int NbSteps = 1; NbSteps<=NbMappingSteps; NbSteps++)
        {
            if(NbSteps==1) localPinball = pinball;
            else
            {
                //cout<<"Mapper3DClass::performNearest()->____cannot map nodeID: "<<nodeID<<" increasing pinball from: "<<localPinball;
                localPinball = localPinball*1.1;
                //cout<<" to: "<<localPinball<<"____"<<endl;
            }

            //! ------------------------------
            //! iterate over the source nodes
            //! ------------------------------
            for(int n=0; n<mySourceNodes.size(); n++)
            {
                const nodeSource &aNodeSource = mySourceNodes[n];
                aTargetMeshNode.put(aNodeSource);
                if(aTargetMeshNode.pinball<localPinball)
                {
                    int nodeID = aTargetMeshNode.nodeID;
                    const std::vector<double> &vecVal = aTargetMeshNode.values;
                    //myMultiRes.insert(nodeID,vecVal);
                    myMultiRes.insert(std::make_pair(nodeID,vecVal));
                    success = true;
                    break;
                }
            }
            if(success)
            {
                NbMappedNodes++;
                break;
            }
        }

        //! --------------
        //! send progress
        //! --------------
        if(n%aProgressStep==0)
        {
            if(myProgressIndicator!=Q_NULLPTR)
            {
                QProgressEvent *aPrgEvent = new QProgressEvent(QProgressEvent_None,0,9999,0,"",QProgressEvent_Update,0,9999,n);
                QApplication::postEvent(myProgressIndicator,aPrgEvent);
                QApplication::processEvents();
            }
        }
    }

    //! -----------------------------------------------------------------
    //! check the number of mapped target nodes with the current pinball
    //! -----------------------------------------------------------------
    double p = 100.0*NbMappedNodes/nMap.Extent();
    cout<<"Mapper3DClass::performNearest()->_____mapping finished. Mapped nodes: "<<p<<" %____"<<endl;
}

//! -----------------------------------------------------
//! function: interpolate
//! details:  switch among several interpolation methods
//! -----------------------------------------------------
void Mapper3DClass::perform(int theAlgo, const opencascade::handle<MeshVS_DataSource> &theTargetMesh)
{
    if(theTargetMesh.IsNull()) this->setTargetMesh(theTargetMesh);

    switch(theAlgo)
    {
    case 0:
        //! "in pinball" algo
        this->performNearest();
        break;

    case 1:
        //! nearest
        this->performNearestNeighboring();
        break;

    case 2: break;
        //! shape functions
        this->performShapeFunctions();
        break;
    }
}

//! -------------------------------------
//! function: performNearest Neighboring
//! details:
//! -------------------------------------
void Mapper3DClass::performNearestNeighboring(double pinball)
{
    cout<<"Mapper3DClass::performNearestNeighboring()->____function called____"<<endl;

    //! -------------------
    //! read some settings
    //! -------------------
    int NbMappingSteps = 1;
    if(myRemapFlag == true) NbMappingSteps = myRemappingSteps+1;
    else NbMappingSteps = 1;

    double dx = (myXmax-myXmin)/myNbucketsX;
    double dy = (myYmax-myYmin)/myNbucketsY;
    double dz = (myZmax-myZmin)/myNbucketsZ;

    //! -------------------------
    //! initialize with OCC mesh
    //! -------------------------
    TColStd_PackedMapOfInteger nMap = myTargetMesh->GetAllNodes();

    //! --------------
    //! send progress
    //! --------------
    if(myProgressIndicator!=Q_NULLPTR)
    {
        myProgressIndicator->setSecondaryBarVisible(true);
        QProgressEvent *aPrgEvent = new QProgressEvent(QProgressEvent_None,0,9999,0,"Start intepolating",QProgressEvent_Init,0,nMap.Extent()-1,0);
        QApplication::postEvent(myProgressIndicator,aPrgEvent);
        QApplication::processEvents();
        Sleep(500);
    }

    //! ----------------------------
    //! progress bar increment 2.5%
    //! ----------------------------
    int aProgressStep = int(nMap.Extent()*0.025);

    //! -----------------
    //! the target nodes
    //! -----------------
    int NbMappedNodes = 0;
    MeshVS_EntityType type;
    int NbNodes;
    double buf[3];
    int n=0;

    for(TColStd_MapIteratorOfPackedMapOfInteger it(nMap);it.More();it.Next(), n++)
    {
        int nodeID = it.Key();
        TColStd_Array1OfReal coords(*buf,1,3);
        myTargetMesh->GetGeom(nodeID,false,coords,NbNodes,type);
        targetMeshNode aTargetMeshNode(coords(1),coords(2),coords(3),nodeID);

        bool success = false;
        double localPinball = 0;
        for(int NbSteps = 1; NbSteps<=NbMappingSteps; NbSteps++)
        {
            if(NbSteps==1) localPinball = pinball;
            else localPinball = localPinball*1.5;
            aTargetMeshNode.setPinball(localPinball);
            double minDistance=1e80;

            //! ------------------------------------------------
            //! center of the element: once the center is known
            //! it is possible to get its buckets
            //! ------------------------------------------------
            double xc, yc, zc;
            xc = aTargetMeshNode.x;
            yc = aTargetMeshNode.y;
            zc = aTargetMeshNode.z;

            //! ------------------------------------------------
            //! extract the source nodes of the selected bucket
            //! ------------------------------------------------
            int bucketX_targetElement=floor((xc-myXmin)/dx);
            int bucketY_targetElement=floor((yc-myYmin)/dy);
            int bucketZ_targetElement=floor((zc-myZmin)/dz);

            if(bucketX_targetElement==myNbucketsX) bucketX_targetElement--;
            if(bucketY_targetElement==myNbucketsY) bucketY_targetElement--;
            if(bucketZ_targetElement==myNbucketsZ) bucketZ_targetElement--;

            std::vector<nodeSource> vecSourceNodesFiltered = myBuckets.value(bucketX_targetElement,bucketY_targetElement,bucketZ_targetElement);

            //! ------------------------------
            //! iterate over the source nodes
            //! ------------------------------
            for(int n=0; n<vecSourceNodesFiltered.size(); n++)
            {
                double curDistance;
                const nodeSource &aNodeSource = vecSourceNodesFiltered[n];
                aTargetMeshNode.put(aNodeSource);
                if(aTargetMeshNode.hasValue)
                {
                    curDistance = aTargetMeshNode.pinball;
                    if(curDistance<minDistance)
                    {
                        int nodeID = aTargetMeshNode.nodeID;
                        const std::vector<double> &vecVal = aTargetMeshNode.values;
                        //myMultiRes.insert(nodeID,vecVal);
                        myMultiRes.insert(std::make_pair(nodeID,vecVal));
                        minDistance = curDistance;
                        success = true;
                    }
                }
            }
            if(success)
            {
                NbMappedNodes++;
                break;              //! exit remapping cycle
            }
        }

        //! --------------
        //! send progress
        //! --------------
        if(n%aProgressStep==0)
        {
            if(myProgressIndicator!=Q_NULLPTR)
            {
                QProgressEvent *aPrgEvent = new QProgressEvent(QProgressEvent_None,0,9999,0,"",QProgressEvent_Update,0,9999,n);
                QApplication::postEvent(myProgressIndicator,aPrgEvent);
                QApplication::processEvents();
            }
        }
    }

    //! -----------------------------------------------------------------
    //! check the number of mapped target nodes with the current pinball
    //! -----------------------------------------------------------------
    double p = 100.0*NbMappedNodes/nMap.Extent();
    cout<<"Mapper3DClass::performNearest()->_____mapping finished. Mapped nodes: "<<p<<" %____"<<endl;
}
