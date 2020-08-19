//! ----------------
//! custom includes
//! ----------------
#include "isostripbuilder.h"
#include "isostrip.h"
#include <meshelement2d.h>  // an attempt to optimize: used in ::perform1()

//! ----
//! OCC
//! ----
#include <MeshVS_DataSource.hxx>
#include <MeshVS_HArray1OfSequenceOfInteger.hxx>
#include <MeshVS_EntityType.hxx>
#include <TColStd_SequenceOfInteger.hxx>
#include <TColStd_Array1OfInteger.hxx>
#include <TColStd_Array1OfReal.hxx>

//! ----
//! C++
//! ----
#include <chrono>
#include <memory>
#include <iostream>
using namespace std;

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
isoStripBuilder::isoStripBuilder(const occHandle(MeshVS_DataSource) &aMeshDS, QObject *parent):
    myMeshDS(aMeshDS),
    QObject(parent)
{    
    ;
}

//! ----------------------------
//! function: setMeshDataSource
//! details:
//! ----------------------------
void isoStripBuilder::setMeshDataSource(const occHandle(MeshVS_DataSource) &aMeshDS)
{
    myMeshDS = aMeshDS;
}

//! -----------------------
//! function: setIsoStrips
//! details:
//! -----------------------
void isoStripBuilder::setIsoStrips(const std::vector<isoStrip> &theIsoStrips)
{
    myIsoStrips = theIsoStrips;
}

//! --------------------
//! function: setValues
//! details:
//! --------------------
void isoStripBuilder::setValues(const std::map<int, double> &values)
{
    myValues = values;
}

//! ------------------
//! function: perform
//! details:
//! ------------------
bool isoStripBuilder::perform(std::vector<meshElementByCoords> &vecMeshElements)
{
    if(myIsoStrips.size()==0) return false;
    if(myMeshDS.IsNull()) return false;
    if(myValues.size()==0)
    {
        cerr<<"isoStripBuilder::perform()->____empty results____"<<endl;
        return false;
    }

    //! ------------------------
    //! performance measurement
    //! ------------------------
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    //! ---------------------------------------------------------
    //! classify the nodes according to the isostrip definitions
    //! ---------------------------------------------------------
    myMultiMap<int,int> pointToIsoStrip;
    this->classifyNodes(pointToIsoStrip);

    //! ----------------------
    //! vector of face tables
    //! ----------------------
    std::vector<faceTable> vecFaces;

    //! -------------------------------------------------------------
    //! iterate over the elements and over the faces of each element
    //! -------------------------------------------------------------
    TColStd_MapIteratorOfPackedMapOfInteger it(myMeshDS->GetAllElements());
    int NbElements = myMeshDS->GetAllElements().Extent();
    for(int localElementID = 1; localElementID<=NbElements; localElementID++,it.Next())
    {
        int globalElementID = it.Key();
        //cout<<"____local element ID: "<<localElementID<<" globalNodeID: "<<globalElementID<<"____"<<endl;

        int NbNodes;
        occHandle(MeshVS_HArray1OfSequenceOfInteger) topology;
        myMeshDS->Get3DGeom(globalElementID,NbNodes,topology);

        int NbFaces = topology->Length();
        //cout<<"____NbFaces: "<<NbFaces<<"____"<<endl;

        int NbNodesOfElement, b[20];
        TColStd_Array1OfInteger nodeIDs(*b,1,20);
        myMeshDS->GetNodesByElement(globalElementID,nodeIDs,NbNodesOfElement);

        for(int f = 1; f<=NbFaces; f++)
        {
            faceTable aFaceTable;

            //! --------------------------------------------------
            //! put the existing nodes of the face into the table
            //! --------------------------------------------------
            TColStd_SequenceOfInteger aSeq = topology->Value(f);
            int N = aSeq.Length();
            for(int i = 0; i<N; i++)
            {
                //! -------------------------------------------------------
                //! put the face points into the columns of the face table
                //! -------------------------------------------------------
                int index = aSeq.Value(i+1)+1;
                int globalNodeID = nodeIDs(index);
                //cout<<"____(index, globalNodeID)____("<<index<<", "<<globalNodeID<<")____"<<endl;

                const isoStripList &isoStripOf = pointToIsoStrip.values(globalNodeID);

                int maxIsoStrip1 = -1e10;
                int minIsoStrip1 = 1e10;
                for(int i=0; i<isoStripOf.size(); i++)
                {
                    int curIsoStripNb = isoStripOf[i];
                    double c[3];
                    this->pointCoord(c,globalNodeID);
                    isoStripPoint aP(c[0],c[1],c[2],myValues.at(globalNodeID));
                    aFaceTable.appendAtCol(curIsoStripNb,aP);
                    if(curIsoStripNb>=maxIsoStrip1) maxIsoStrip1 = curIsoStripNb;
                    if(curIsoStripNb<=minIsoStrip1) minIsoStrip1 = curIsoStripNb;
                }
                //cout<<"____(sx,dx) = ("<<minIsoStrip1<<", "<<maxIsoStrip1<<")____"<<endl;
            }

            //! ------------------------------
            //! scan the segments of the face
            //! ------------------------------
            for(int i = 0; i<N; i++)
            {
                //cout<<"____scanning segment n: "<<i+1<<"____"<<endl;

                //! --------------
                //! recompute ...
                //! --------------
                int index1 = aSeq.Value(i%N+1)+1;
                int globalNodeID1 = nodeIDs(index1);

                const isoStripList &isoStripOf1 = pointToIsoStrip.values(globalNodeID1);

                int maxIsoStrip1 = -1e10;
                int minIsoStrip1 = 1e10;
                for(int i=0; i<isoStripOf1.size(); i++)
                {
                    int curIsoStripNb = isoStripOf1[i];
                    if(curIsoStripNb<=minIsoStrip1) minIsoStrip1 = curIsoStripNb;
                    if(curIsoStripNb>=maxIsoStrip1) maxIsoStrip1 = curIsoStripNb;
                }

                //! --------------
                //! recompute ...
                //! --------------
                int index2 = aSeq.Value((i+1)%N+1)+1;
                int globalNodeID2 = nodeIDs(index2);

                const isoStripList &isoStripOf2 = pointToIsoStrip.values(globalNodeID2);

                int maxIsoStrip2 = -1e10;
                int minIsoStrip2 = 1e10;

                for(int i=0; i<isoStripOf2.size(); i++)
                {
                    int curIsoStripNb = isoStripOf2[i];
                    if(curIsoStripNb<=minIsoStrip2) minIsoStrip2 = curIsoStripNb;
                    if(curIsoStripNb>=maxIsoStrip2) maxIsoStrip2 = curIsoStripNb;
                }

                double valueOfPoint1 = myValues.at(globalNodeID1);
                double valueOfPoint2 = myValues.at(globalNodeID2);

                //cout<<"____"<<globalNodeID1<<", "<<valueOfPoint1<<"____"<<endl;
                //cout<<"____"<<globalNodeID2<<", "<<valueOfPoint2<<"____"<<endl;

                //! ---------------------------------
                //! zero gradient between two points
                //! ---------------------------------
                if(valueOfPoint1 == valueOfPoint2) continue;

                //! ------------------
                //! increasing values
                //! ------------------
                if(valueOfPoint2-valueOfPoint1>0)
                {
                    int deltaStrip = minIsoStrip2 - maxIsoStrip1;

                    if(deltaStrip==0) continue;

                    double csp[3];
                    this->pointCoord(csp,globalNodeID1);
                    isoStripPoint sP(csp[0],csp[1],csp[2],myValues.at(maxIsoStrip1));

                    double cep[3];
                    this->pointCoord(cep,globalNodeID2);
                    isoStripPoint eP(cep[0],cep[1],cep[2],myValues.at(globalNodeID2));

                    int startIndex = maxIsoStrip1;
                    for(int k=0; k<abs(deltaStrip); k++,startIndex++)
                    {
                        //cout<<"____startIndex case >: "<<startIndex<<"____"<<endl;
                        double newPointValue = myIsoStrips.at(startIndex).vmax;
                        //cout<<"____"<<newPointValue<<"____"<<endl;
                        double t = (newPointValue-valueOfPoint1)/(valueOfPoint2-valueOfPoint1);
                        //cout<<"____t = "<<t<<"____"<<endl;
                        double P1[3], P2[3];
                        this->pointCoord(P1,globalNodeID1);
                        this->pointCoord(P2,globalNodeID2);
                        double xP = P1[0] + t * (P2[0]-P1[0]);
                        double yP = P1[1] + t * (P2[1]-P1[1]);
                        double zP = P1[2] + t * (P2[2]-P1[2]);

                        isoStripPoint aP(xP,yP,zP,newPointValue);

                        //cout<<"_>___trying inserting after____"<<endl;
                        aFaceTable.insertAfter(startIndex,sP,aP);
                        //cout<<"_>___insert after done____"<<endl;

                        //cout<<"_>___trying inserting before____"<<endl;
                        bool isDone = aFaceTable.insertBefore(startIndex+1,eP,aP);
                        if(isDone==false) aFaceTable.appendAtCol(startIndex+1,aP);
                        //cout<<"_>___insert before done____"<<endl;
                        sP = aP;
                    }
                }
                //! ------------------
                //! decreasing values
                //! ------------------
                else
                {
                    int deltaStrip = minIsoStrip1 - maxIsoStrip2;
                    if(deltaStrip==0) continue;

                    double csp[3];
                    this->pointCoord(csp,globalNodeID2);
                    isoStripPoint sP(csp[0],csp[1],csp[2],myValues.at(globalNodeID2));

                    double cep[3];
                    this->pointCoord(cep,globalNodeID1);
                    isoStripPoint eP(cep[0],cep[1],cep[2],myValues.at(globalNodeID1));

                    int startIndex = maxIsoStrip2;
                    for(int k=0; k<abs(deltaStrip); k++, startIndex++)
                    {                        
                        //cout<<"____startIndex  case <: "<<startIndex<<"____"<<endl;
                        double newPointValue = myIsoStrips.at(startIndex).vmax;
                        double t = (newPointValue-valueOfPoint2)/(valueOfPoint1-valueOfPoint2);
                        double P1[3], P2[3];
                        this->pointCoord(P1,globalNodeID2);
                        this->pointCoord(P2,globalNodeID1);
                        double xP = P1[0] + t * (P2[0]-P1[0]);
                        double yP = P1[1] + t * (P2[1]-P1[1]);
                        double zP = P1[2] + t * (P2[2]-P1[2]);

                        //cout<<"____case deltastrip < 0: t = "<<t<<" newPoint("<<xP<<", "<<yP<<", "<<zP<<")____"<<endl;

                        isoStripPoint aP(xP,yP,zP,newPointValue);

                        //cout<<"_<___trying inserting before____"<<endl;
                        aFaceTable.insertBefore(startIndex,sP,aP);
                        //cout<<"_<___insert before done____"<<endl;

                        //cout<<"_<___trying inserting after____"<<endl;
                        bool isDone = aFaceTable.insertAfter(startIndex+1,eP,aP);
                        if(isDone==false) aFaceTable.appendAtCol(startIndex+1,aP);
                        //cout<<"_<___insert after done____"<<endl;
                        sP = aP;
                    }
                }
            }
            //! ------------------------
            //! pile up the face tables
            //! ------------------------
            vecFaces.push_back(aFaceTable);
        }
    }

    //! ------------------------
    //! performance measurement
    //! ------------------------
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[microseconds]" << std::endl;

    //! --------
    //! testing
    //! --------
    this->getAllElements(vecFaces,vecMeshElements);

    return true;
}

//! -------------------------
//! function: getAllElements
//! details:
//! -------------------------
void isoStripBuilder::getAllElements(const std::vector<faceTable> &vecFaceTables, std::vector<meshElementByCoords> &vecMeshElements)
{
    cout<<"isoStripBuilder::getAllElements()->____function called. Number of face tables: "<<vecFaceTables.size()<<"____"<<endl;
    for(int i=0; i<vecFaceTables.size(); i++)
    {
        std::vector<meshElementByCoords> vecMeshElementsOfFace;
        vecFaceTables[i].getElements(vecMeshElementsOfFace);
        for(int i=0; i<vecMeshElementsOfFace.size(); i++) vecMeshElements.push_back(vecMeshElementsOfFace[i]);
    }
    cout<<"isoStripBuilder::getAllElements()->____number of elements: "<<vecMeshElements.size()<<"____"<<endl;
    /*
    int N = vecFaceTables.size();
    const int NbThreads = 17;
    int NbFacesForThread = N/NbThreads;
    int NbFaceTables[NbThreads];
    for(int n=0; n<NbThreads-1; n++) NbFaceTables[n] = NbFacesForThread;
    NbFaceTables[NbThreads-1] = N-(NbThreads-1)*NbFacesForThread;
    */
}

//! --------------------------------------------------------------------------------------
//! function: getIsoStripElements
//! details:  input all the face tables and output the mesh elements grouped by isostrips
//! --------------------------------------------------------------------------------------
void isoStripBuilder::getIsoStripElements(const std::vector<faceTable> &allFaceTables,
                                          std::multimap<int,meshElementByCoords> &meshElementsByIsoStripNb)
{
    int NbFaces = (int)allFaceTables.size();
    int NbIsostrip = (int)myIsoStrips.size();

    //! ---------------------------
    //! iterate over all the faces
    //! ---------------------------
    for(int i=0; i<NbFaces; i++)
    {
        for(int c = 0; c<NbIsostrip; c++)
        {
            meshElementByCoords me;
            allFaceTables[i].getElementOfColumn(c,me);
            std::pair<int,meshElementByCoords> apair;
            apair.first = c;
            apair.second = me;
            meshElementsByIsoStripNb.insert(apair);
        }
    }
}

//! ----------------------
//! function: pointCoord
//! details:  helper
//! ----------------------
void isoStripBuilder::pointCoord(double *c, int globalNodeID)
{
    MeshVS_EntityType aType;
    int NbNodes;
    double b[3];
    TColStd_Array1OfReal cn(*b,1,3);
    myMeshDS->GetGeom(globalNodeID, false, cn, NbNodes, aType);
    for(int i=0; i<3; i++) c[i] = cn(i+1);
}

//! ------------------------
//! function: classifyNodes
//! details:
//! ------------------------
void isoStripBuilder::classifyNodes(myMultiMap<int,int> &pointToIsoStrip)
{
    //cout<<"isoStripBuilder::classifyNodes()->____function called____"<<endl;
    TColStd_MapIteratorOfPackedMapOfInteger it(myMeshDS->GetAllNodes());
    int NbNodes = myMeshDS->GetAllNodes().Extent();
    for(int localNodeID = 1; localNodeID<=NbNodes; localNodeID++,it.Next())
    {
        int globalNodeID = it.Key();
        double val = myValues.at(globalNodeID);
        for(int isoStripNb = 0; isoStripNb<myIsoStrips.size(); isoStripNb++)
        {
            isoStrip curIsoStrip = myIsoStrips.at(isoStripNb);
            if(curIsoStrip.contains(val)) pointToIsoStrip.insertMulti(globalNodeID,isoStripNb);
        }
    }
}

//! ------------------
//! function: perform
//! details:
//! ------------------
bool isoStripBuilder::perform1(std::vector<meshElementByCoords> &vecMeshElements)
{
    if(myIsoStrips.size()==0) return false;
    if(myMeshDS.IsNull()) return false;
    if(myValues.size()==0)
    {
        cerr<<"isoStripBuilder::perform()->____empty results____"<<endl;
        return false;
    }

    //! ------------------------
    //! performance measurement
    //! ------------------------
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    //! ---------------------------------------------------------
    //! classify the nodes according to the isostrip definitions
    //! ---------------------------------------------------------
    myMultiMap<int,int> pointToIsoStrip;
    this->classifyNodes(pointToIsoStrip);

    //! ----------------------
    //! vector of face tables
    //! ----------------------
    std::vector<faceTable> vecFaces;
    std::map<meshElement2D,int> mapProcessedFaces;  // "1" processed

    //! -------------------------------------------------------------
    //! iterate over the elements and over the faces of each element
    //! -------------------------------------------------------------
    TColStd_MapIteratorOfPackedMapOfInteger it(myMeshDS->GetAllElements());
    int NbElements = myMeshDS->GetAllElements().Extent();
    for(int localElementID = 1; localElementID<=NbElements; localElementID++,it.Next())
    {
        int globalElementID = it.Key();

        int NbNodes;
        occHandle(MeshVS_HArray1OfSequenceOfInteger) topology;
        myMeshDS->Get3DGeom(globalElementID,NbNodes,topology);

        int NbFaces = topology->Length();

        int NbNodesOfElement, b[20];
        TColStd_Array1OfInteger nodeIDs(*b,1,20);
        myMeshDS->GetNodesByElement(globalElementID,nodeIDs,NbNodesOfElement);

        for(int f = 1; f<=NbFaces; f++)
        {
            faceTable aFaceTable;
            QList<int> nodeIDs_;

            //! --------------------------------------------------
            //! put the existing nodes of the face into the table
            //! --------------------------------------------------
            TColStd_SequenceOfInteger aSeq = topology->Value(f);
            int N = aSeq.Length();
            for(int i = 0; i<N; i++)
            {
                //! -------------------------------------------------------
                //! put the face points into the columns of the face table
                //! -------------------------------------------------------
                int index = aSeq.Value(i+1)+1;
                int globalNodeID = nodeIDs(index);
                nodeIDs_<<globalNodeID;

                const isoStripList &isoStripOf = pointToIsoStrip.values(globalNodeID);

                for(int i=0; i<isoStripOf.size(); i++)
                {
                    int curIsoStripNb = isoStripOf[i];
                    double c[3];
                    this->pointCoord(c,globalNodeID);
                    isoStripPoint aP(c[0],c[1],c[2],myValues.at(globalNodeID));
                    aFaceTable.appendAtCol(curIsoStripNb,aP);
                }
            }

            meshElement2D aFaceMesh_(nodeIDs_);
            if(mapProcessedFaces.count(aFaceMesh_)!=0) continue;
            mapProcessedFaces.insert(std::make_pair(aFaceMesh_,1));

            //! ------------------------------
            //! scan the segments of the face
            //! ------------------------------
            for(int i = 0; i<N; i++)
            {
                //cout<<"____scanning segment n: "<<i+1<<"____"<<endl;

                //! --------------
                //! recompute ...
                //! --------------
                int index1 = aSeq.Value(i%N+1)+1;
                int globalNodeID1 = nodeIDs(index1);

                const isoStripList &isoStripOf1 = pointToIsoStrip.values(globalNodeID1);
                int maxIsoStrip1 = -1e10;
                int minIsoStrip1 = 1e10;
                for(int i=0; i<isoStripOf1.size(); i++)
                {
                    int curIsoStripNb = isoStripOf1[i];
                    if(curIsoStripNb<=minIsoStrip1) minIsoStrip1 = curIsoStripNb;
                    if(curIsoStripNb>=maxIsoStrip1) maxIsoStrip1 = curIsoStripNb;
                }

                //! --------------
                //! recompute ...
                //! --------------
                int index2 = aSeq.Value((i+1)%N+1)+1;
                int globalNodeID2 = nodeIDs(index2);

                const isoStripList &isoStripOf2 = pointToIsoStrip.values(globalNodeID2);

                int maxIsoStrip2 = -1e10;
                int minIsoStrip2 = 1e10;

                for(int i=0; i<isoStripOf2.size(); i++)
                {
                    int curIsoStripNb = isoStripOf2[i];
                    if(curIsoStripNb<=minIsoStrip2) minIsoStrip2 = curIsoStripNb;
                    if(curIsoStripNb>=maxIsoStrip2) maxIsoStrip2 = curIsoStripNb;
                }

                double valueOfPoint1 = myValues.at(globalNodeID1);
                double valueOfPoint2 = myValues.at(globalNodeID2);

                //cout<<"____"<<globalNodeID1<<", "<<valueOfPoint1<<"____"<<endl;
                //cout<<"____"<<globalNodeID2<<", "<<valueOfPoint2<<"____"<<endl;

                //! ---------------------------------
                //! zero gradient between two points
                //! ---------------------------------
                if(valueOfPoint1 == valueOfPoint2) continue;

                //! ------------------
                //! increasing values
                //! ------------------
                if(valueOfPoint2-valueOfPoint1>0)
                {
                    int deltaStrip = minIsoStrip2 - maxIsoStrip1;

                    if(deltaStrip==0) continue;

                    double csp[3];
                    this->pointCoord(csp,globalNodeID1);
                    isoStripPoint sP(csp[0],csp[1],csp[2],myValues.at(maxIsoStrip1));

                    double cep[3];
                    this->pointCoord(cep,globalNodeID2);
                    isoStripPoint eP(cep[0],cep[1],cep[2],myValues.at(globalNodeID2));

                    int startIndex = maxIsoStrip1;
                    for(int k=0; k<abs(deltaStrip); k++,startIndex++)
                    {
                        //cout<<"____startIndex case >: "<<startIndex<<"____"<<endl;
                        double newPointValue = myIsoStrips.at(startIndex).vmax;
                        //cout<<"____"<<newPointValue<<"____"<<endl;
                        double t = (newPointValue-valueOfPoint1)/(valueOfPoint2-valueOfPoint1);
                        //cout<<"____t = "<<t<<"____"<<endl;
                        double P1[3], P2[3];
                        this->pointCoord(P1,globalNodeID1);
                        this->pointCoord(P2,globalNodeID2);
                        double xP = P1[0] + t * (P2[0]-P1[0]);
                        double yP = P1[1] + t * (P2[1]-P1[1]);
                        double zP = P1[2] + t * (P2[2]-P1[2]);

                        isoStripPoint aP(xP,yP,zP,newPointValue);

                        //cout<<"_>___trying inserting after____"<<endl;
                        aFaceTable.insertAfter(startIndex,sP,aP);
                        //cout<<"_>___insert after done____"<<endl;

                        //cout<<"_>___trying inserting before____"<<endl;
                        bool isDone = aFaceTable.insertBefore(startIndex+1,eP,aP);
                        if(isDone==false) aFaceTable.appendAtCol(startIndex+1,aP);
                        //cout<<"_>___insert before done____"<<endl;
                        sP = aP;
                    }
                }
                //! ------------------
                //! decreasing values
                //! ------------------
                else
                {
                    int deltaStrip = minIsoStrip1 - maxIsoStrip2;
                    if(deltaStrip==0) continue;

                    double csp[3];
                    this->pointCoord(csp,globalNodeID2);
                    isoStripPoint sP(csp[0],csp[1],csp[2],myValues.at(globalNodeID2));

                    double cep[3];
                    this->pointCoord(cep,globalNodeID1);
                    isoStripPoint eP(cep[0],cep[1],cep[2],myValues.at(globalNodeID1));

                    int startIndex = maxIsoStrip2;
                    for(int k=0; k<abs(deltaStrip); k++, startIndex++)
                    {
                        //cout<<"____startIndex  case <: "<<startIndex<<"____"<<endl;
                        double newPointValue = myIsoStrips.at(startIndex).vmax;
                        double t = (newPointValue-valueOfPoint2)/(valueOfPoint1-valueOfPoint2);
                        double P1[3], P2[3];
                        this->pointCoord(P1,globalNodeID2);
                        this->pointCoord(P2,globalNodeID1);
                        double xP = P1[0] + t * (P2[0]-P1[0]);
                        double yP = P1[1] + t * (P2[1]-P1[1]);
                        double zP = P1[2] + t * (P2[2]-P1[2]);

                        //cout<<"____case deltastrip < 0: t = "<<t<<" newPoint("<<xP<<", "<<yP<<", "<<zP<<")____"<<endl;

                        isoStripPoint aP(xP,yP,zP,newPointValue);

                        //cout<<"_<___trying inserting before____"<<endl;
                        aFaceTable.insertBefore(startIndex,sP,aP);
                        //cout<<"_<___insert before done____"<<endl;

                        //cout<<"_<___trying inserting after____"<<endl;
                        bool isDone = aFaceTable.insertAfter(startIndex+1,eP,aP);
                        if(isDone==false) aFaceTable.appendAtCol(startIndex+1,aP);
                        //cout<<"_<___insert after done____"<<endl;
                        sP = aP;
                    }
                }
            }
            //! ------------------------
            //! pile up the face tables
            //! ------------------------
            vecFaces.push_back(aFaceTable);
        }
    }

    //! ------------------------
    //! performance measurement
    //! ------------------------
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[microseconds]" << std::endl;

    //! --------
    //! testing
    //! --------
    this->getAllElements(vecFaces,vecMeshElements);

    //! ----------
    //! testing 1
    //! ----------
    //myMultiMap<int,meshElementByCoords> mm;
    //this->getIsoStripElements(vecFaces,mm);

    return true;
}

//! ----------------------------
//! function: performIsoSurface
//! details:
//! ----------------------------
bool isoStripBuilder::computeFaceTables()
{
    if(myIsoStrips.size()==0) return false;
    if(myMeshDS.IsNull()) return false;
    if(myValues.size()==0)
    {
        cerr<<"isoStripBuilder::perform()->____empty results____"<<endl;
        return false;
    }

    myMapElementFaceTables.clear();

    //! ------------------------
    //! performance measurement
    //! ------------------------
    //std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    //! ---------------------------------------------------------
    //! classify the nodes according to the isostrip definitions
    //! ---------------------------------------------------------
    myMultiMap<int,int> pointToIsoStrip;
    this->classifyNodes(pointToIsoStrip);

    //! ----------------------
    //! vector of face tables
    //! ----------------------
    std::vector<faceTable> vecFaces;

    //! -------------------------------------------------------------
    //! iterate over the elements and over the faces of each element
    //! -------------------------------------------------------------
    TColStd_MapIteratorOfPackedMapOfInteger it(myMeshDS->GetAllElements());
    int NbElements = myMeshDS->GetAllElements().Extent();
    for(int localElementID = 1; localElementID<=NbElements; localElementID++,it.Next())
    {
        int globalElementID = it.Key();

        int NbNodes;
        occHandle(MeshVS_HArray1OfSequenceOfInteger) topology;
        myMeshDS->Get3DGeom(globalElementID,NbNodes,topology);

        int NbFaces = topology->Length();
        int NbNodesOfElement, b[20];
        TColStd_Array1OfInteger nodeIDs(*b,1,20);
        myMeshDS->GetNodesByElement(globalElementID,nodeIDs,NbNodesOfElement);

        std::vector<faceTable> elementFaceTables;
        for(int f = 1; f<=NbFaces; f++)
        {
            faceTable aFaceTable;

            //! --------------------------------------------------
            //! put the existing nodes of the face into the table
            //! --------------------------------------------------
            TColStd_SequenceOfInteger aSeq = topology->Value(f);
            int N = aSeq.Length();
            for(int i = 0; i<N; i++)
            {
                //! -------------------------------------------------------
                //! put the face points into the columns of the face table
                //! -------------------------------------------------------
                int index = aSeq.Value(i+1)+1;
                int globalNodeID = nodeIDs(index);
                const isoStripList &isoStripOf = pointToIsoStrip.values(globalNodeID);

                int maxIsoStrip1 = -1e10;
                int minIsoStrip1 = 1e10;
                for(int i=0; i<isoStripOf.size(); i++)
                {
                    int curIsoStripNb = isoStripOf[i];

                    double c[3];
                    this->pointCoord(c,globalNodeID);
                    isoStripPoint aP(c[0],c[1],c[2],myValues[globalNodeID],0);

                    //! ---------------------------------------------------------------------
                    //! determine the type of point
                    //! The code below has been commented since an existing mesh point with
                    //! value could belong to two isostrips if the value is suitable
                    //! ---------------------------------------------------------------------
                    if(myValues[globalNodeID]==myIsoStrips[curIsoStripNb].vmin)
                    {
                        isoStripPoint aP(c[0],c[1],c[2],myValues[globalNodeID],-1);
                        aFaceTable.appendAtCol(curIsoStripNb,aP);
                        if(curIsoStripNb>=maxIsoStrip1) maxIsoStrip1 = curIsoStripNb;
                        if(curIsoStripNb<=minIsoStrip1) minIsoStrip1 = curIsoStripNb;
                        continue;
                    }
                    if(myValues[globalNodeID]==myIsoStrips[curIsoStripNb].vmax)
                    {
                        isoStripPoint aP(c[0],c[1],c[2],myValues[globalNodeID],1);
                        aFaceTable.appendAtCol(curIsoStripNb,aP);
                        if(curIsoStripNb>=maxIsoStrip1) maxIsoStrip1 = curIsoStripNb;
                        if(curIsoStripNb<=minIsoStrip1) minIsoStrip1 = curIsoStripNb;
                        continue;
                    }
                    aFaceTable.appendAtCol(curIsoStripNb,aP);
                    if(curIsoStripNb>=maxIsoStrip1) maxIsoStrip1 = curIsoStripNb;
                    if(curIsoStripNb<=minIsoStrip1) minIsoStrip1 = curIsoStripNb;
                }
            }

            //! ------------------------------
            //! scan the segments of the face
            //! ------------------------------
            for(int i = 0; i<N; i++)
            {
                //! --------------
                //! recompute ...
                //! --------------
                int index1 = aSeq.Value(i%N+1)+1;
                int globalNodeID1 = nodeIDs(index1);
                const isoStripList &isoStripOf1 = pointToIsoStrip.values(globalNodeID1);
                int maxIsoStrip1 = -1e10;
                int minIsoStrip1 = 1e10;
                for(int i=0; i<isoStripOf1.size(); i++)
                {
                    int curIsoStripNb = isoStripOf1[i];
                    if(curIsoStripNb<=minIsoStrip1) minIsoStrip1 = curIsoStripNb;
                    if(curIsoStripNb>=maxIsoStrip1) maxIsoStrip1 = curIsoStripNb;
                }

                //! --------------
                //! recompute ...
                //! --------------
                int index2 = aSeq.Value((i+1)%N+1)+1;
                int globalNodeID2 = nodeIDs(index2);
                const isoStripList &isoStripOf2 = pointToIsoStrip.values(globalNodeID2);
                int maxIsoStrip2 = -1e10;
                int minIsoStrip2 = 1e10;
                for(int i=0; i<isoStripOf2.size(); i++)
                {
                    int curIsoStripNb = isoStripOf2[i];
                    if(curIsoStripNb<=minIsoStrip2) minIsoStrip2 = curIsoStripNb;
                    if(curIsoStripNb>=maxIsoStrip2) maxIsoStrip2 = curIsoStripNb;
                }
                double valueOfPoint1 = myValues[globalNodeID1];
                double valueOfPoint2 = myValues[globalNodeID2];

                //! ---------------------------------
                //! zero gradient between two points
                //! ---------------------------------
                if(valueOfPoint1 == valueOfPoint2) continue;

                //! ------------------
                //! increasing values
                //! ------------------
                if(valueOfPoint2-valueOfPoint1>0)
                {
                    int deltaStrip = minIsoStrip2 - maxIsoStrip1;
                    if(deltaStrip==0) continue;

                    double csp[3];
                    this->pointCoord(csp,globalNodeID1);
                    isoStripPoint sP(csp[0],csp[1],csp[2],myValues[maxIsoStrip1]);

                    double cep[3];
                    this->pointCoord(cep,globalNodeID2);
                    isoStripPoint eP(cep[0],cep[1],cep[2],myValues[globalNodeID2]);

                    int startIndex = maxIsoStrip1;
                    for(int k=0; k<abs(deltaStrip); k++,startIndex++)
                    {
                        double newPointValue = myIsoStrips[startIndex].vmax;
                        double t = (newPointValue-valueOfPoint1)/(valueOfPoint2-valueOfPoint1);
                        double P1[3], P2[3];
                        this->pointCoord(P1,globalNodeID1);
                        this->pointCoord(P2,globalNodeID2);
                        double xP = P1[0] + t * (P2[0]-P1[0]);
                        double yP = P1[1] + t * (P2[1]-P1[1]);
                        double zP = P1[2] + t * (P2[2]-P1[2]);

                        isoStripPoint aP(xP,yP,zP,newPointValue,1);
                        aFaceTable.insertAfter(startIndex,sP,aP);
                        isoStripPoint aP_(xP,yP,zP,newPointValue,-1);   // clone of the previous at bottom
                        bool isDone = aFaceTable.insertBefore(startIndex+1,eP,aP_);
                        if(isDone==false) aFaceTable.appendAtCol(startIndex+1,aP_);
                        sP = aP;
                    }
                }
                //! ------------------
                //! decreasing values
                //! ------------------
                else
                {
                    int deltaStrip = minIsoStrip1 - maxIsoStrip2;
                    if(deltaStrip==0) continue;

                    double csp[3];
                    this->pointCoord(csp,globalNodeID2);
                    isoStripPoint sP(csp[0],csp[1],csp[2],myValues[globalNodeID2]);

                    double cep[3];
                    this->pointCoord(cep,globalNodeID1);
                    isoStripPoint eP(cep[0],cep[1],cep[2],myValues[globalNodeID1]);

                    int startIndex = maxIsoStrip2;
                    for(int k=0; k<abs(deltaStrip); k++, startIndex++)
                    {
                        double newPointValue = myIsoStrips[startIndex].vmax;
                        double t = (newPointValue-valueOfPoint2)/(valueOfPoint1-valueOfPoint2);
                        double P1[3], P2[3];
                        this->pointCoord(P1,globalNodeID2);
                        this->pointCoord(P2,globalNodeID1);
                        double xP = P1[0] + t * (P2[0]-P1[0]);
                        double yP = P1[1] + t * (P2[1]-P1[1]);
                        double zP = P1[2] + t * (P2[2]-P1[2]);

                        isoStripPoint aP(xP,yP,zP,newPointValue,1);
                        aFaceTable.insertBefore(startIndex,sP,aP);
                        isoStripPoint aP_(xP,yP,zP,newPointValue,-1);    // clone of the previous at bottom
                        bool isDone = aFaceTable.insertAfter(startIndex+1,eP,aP_);
                        if(isDone==false) aFaceTable.appendAtCol(startIndex+1,aP_);
                        sP = aP;
                    }
                }
            }
            //! ------------------------
            //! pile up the face tables
            //! ------------------------
            vecFaces.push_back(aFaceTable);
            elementFaceTables.push_back(aFaceTable);
        }
        myMapElementFaceTables.insert(std::make_pair(globalElementID,elementFaceTables));
    }
    return true;
}

//! --------------------------------------------------------------------
//! function: performIsoSurface
//! details:  the returning map binds an element global ID with a level
//! --------------------------------------------------------------------
bool isoStripBuilder::performIsoSurface(int NbLevels, std::vector<meshElementByCoords> &vecMeshElements,
                                        std::map<int,int> &mapElementLevel)
{
    //! ----------------------------------
    //! strip table for each element face
    //! ----------------------------------
    bool isDone = this->computeFaceTables();
    if(isDone==false) return false;

    //! ---------------------------
    //! isostrip selection - dummy
    //! ---------------------------
    int aType = -1;  // bottom
    int localElementID  = 0;
    for(int level = 2; level<NbLevels+2; level++)
    {
        //! -------------------------------
        //! the elements of an iso-surface
        //! -------------------------------
        for(std::map<int,std::vector<faceTable>>::const_iterator it = myMapElementFaceTables.cbegin(); it!=myMapElementFaceTables.cend(); it++)
        {
            //int globalElementID = it->first;
            const std::vector<faceTable> &elementTables = it->second;
            if(elementTables.size()<4) continue;
            QList<mesh::meshPoint> pointList;
            for(std::vector<faceTable>::const_iterator it_ = elementTables.cbegin(); it_!=elementTables.cend(); it_++)
            {
                const faceTable &aFaceTable = *it_;
                const std::vector<isoStripPoint> &points = aFaceTable.getColumn(level);
                for(std::vector<isoStripPoint>::const_iterator it__ = points.cbegin(); it__!=points.cend(); it__++)
                {
                    const isoStripPoint &aPoint = *it__;
                    mesh::meshPoint aMeshPoint(aPoint.x,aPoint.y,aPoint.z);
                    if(aPoint.type()!=aType) continue;
                    //if(pointList.contains(aMeshPoint)) continue;      // this question should be clarified
                    pointList.push_back(aMeshPoint);
                }
            }
            int NbPoints = (int)pointList.size();
            if(NbPoints>=3)
            {
                localElementID++;
                meshElementByCoords ame;
                for(int n=0; n<NbPoints; n++) ame.pointList<<pointList[n];
                switch(NbPoints)
                {
                case 3: ame.type = TRIG; break;
                case 4: ame.type = QUAD; break;
                case 5: ame.type = PENTA; break;
                case 6: ame.type = TRIG6; break;
                case 7: ame.type = EPTA; break;
                case 8: ame.type = QUAD8; break;
                }
                vecMeshElements.push_back(ame);
                mapElementLevel.insert(std::make_pair(localElementID,level));
            }
        }
    }
    return true;
}
