//! ----------------
//! custom includes
//! ----------------
#include "isostripbuilder.h"
#include "isostrip.h"

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
#include <memory>
#include <iostream>
using namespace std;

//! ---
//! Qt
//! ---
#include <QMap>

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
void isoStripBuilder::setValues(const QMap<int, double> &values)
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
    if(myValues.isEmpty())
    {
        cerr<<"isoStripBuilder::perform()->____empty results____"<<endl;
        return false;
    }

    std::vector<faceTable> vecFaces;

    //! ---------------------------------------------------------
    //! classify the nodes according to the isostrip definitions
    //! ---------------------------------------------------------
    QMap<int,int> pointToIsoStrip;
    this->classifyNodes(pointToIsoStrip);

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
            //cout<<"____iterating____"<<endl;
            faceTable aFaceTable;
            aFaceTable.resize(myIsoStrips.size());

            //! --------------------------------------------------
            //! put the existing nodes of the face into the table
            //! --------------------------------------------------
            TColStd_SequenceOfInteger aSeq = topology->Value(f);
            //cout<<"____aSeq length: "<<aSeq.Length()<<"____"<<endl;
            int N = aSeq.Length();
            for(int i = 0; i<aSeq.Length(); i++)
            {
                //! -------------------------------------
                //! put the first point into the columns
                //! -------------------------------------
                int index1 = aSeq.Value(i%N+1)+1;
                //cout<<"____index1: "<<index1<<"____"<<endl;
                int globalNodeID1 = nodeIDs(index1);
                //cout<<"____(index1, globalNodeID1)____("<<index1<<", "<<globalNodeID1<<")____"<<endl;

                const QList<int> &isoStripOf1 = pointToIsoStrip.values(globalNodeID1);
                int maxIsoStrip1 = -1e10;
                int minIsoStrip1 = 1e10;
                for(int i=0; i<isoStripOf1.length(); i++)
                {
                    int curIsoStripNb = isoStripOf1.at(i);
                    double c[3];
                    this->pointCoord(c,globalNodeID1);
                    isoStripPoint aP(c[0],c[1],c[2],myValues.value(globalNodeID1));
                    aFaceTable[curIsoStripNb].push_back(aP);
                    if(curIsoStripNb>=maxIsoStrip1) maxIsoStrip1 = curIsoStripNb;
                    if(curIsoStripNb<=minIsoStrip1) minIsoStrip1 = curIsoStripNb;
                }
                //cout<<"____(sx,dx) = ("<<minIsoStrip1<<", "<<maxIsoStrip1<<")____"<<endl;

                //! --------------------------------------
                //! put the second point into the columns
                //! --------------------------------------
                int index2 = aSeq.Value((i+1)%N+1)+1;
                //cout<<"____index2: "<<index2<<"____"<<endl;
                int globalNodeID2 = nodeIDs(index2);
                //cout<<"____(index2, globalNodeID2)____("<<index2<<", "<<globalNodeID2<<")____"<<endl;

                const QList<int> &isoStripOf2 = pointToIsoStrip.values(globalNodeID2);
                int maxIsoStrip2 = -1e10;
                int minIsoStrip2 = 1e10;
                for(int i=0; i<isoStripOf2.length(); i++)
                {
                    int curIsoStripNb = isoStripOf2.at(i);
                    double c[3];
                    this->pointCoord(c,globalNodeID2);
                    isoStripPoint aP(c[0],c[1],c[2],myValues.value(globalNodeID2));
                    aFaceTable[curIsoStripNb].push_back(aP);
                    if(curIsoStripNb<=minIsoStrip2) minIsoStrip2 = curIsoStripNb;
                    if(curIsoStripNb>=maxIsoStrip2) maxIsoStrip2 = curIsoStripNb;
                }
                //cout<<"____(sx,dx) = ("<<minIsoStrip2<<", "<<maxIsoStrip2<<")____"<<endl;
            }

            //! ------------------------------
            //! scan the segments of the face
            //! ------------------------------
            for(int i = 0; i<aSeq.Length(); i++)
            {
                //! recompute ...
                int index1 = aSeq.Value(i%N+1)+1;
                int globalNodeID1 = nodeIDs(index1);

                const QList<int> &isoStripOf1 = pointToIsoStrip.values(globalNodeID1);
                int maxIsoStrip1 = -1e10;
                int minIsoStrip1 = 1e10;
                for(int i=0; i<isoStripOf1.length(); i++)
                {
                    int curIsoStripNb = isoStripOf1.at(i);
                    if(curIsoStripNb<=minIsoStrip1) minIsoStrip1 = curIsoStripNb;
                    if(curIsoStripNb>=maxIsoStrip1) maxIsoStrip1 = curIsoStripNb;
                }

                //! recompute ...
                int index2 = aSeq.Value((i+1)%N+1)+1;
                int globalNodeID2 = nodeIDs(index2);

                const QList<int> &isoStripOf2 = pointToIsoStrip.values(globalNodeID2);
                int maxIsoStrip2 = -1e10;
                int minIsoStrip2 = 1e10;

                for(int i=0; i<isoStripOf2.length(); i++)
                {
                    int curIsoStripNb = isoStripOf2.at(i);
                    if(curIsoStripNb<=minIsoStrip2) minIsoStrip2 = curIsoStripNb;
                    if(curIsoStripNb>=maxIsoStrip2) maxIsoStrip2 = curIsoStripNb;
                }

                double csp[3];
                this->pointCoord(csp,globalNodeID1);
                isoStripPoint sP(csp[0],csp[1],csp[2]);

                double valueOfPoint1 = myValues.value(globalNodeID1);
                double valueOfPoint2 = myValues.value(globalNodeID2);
                if(valueOfPoint2-valueOfPoint1>=0)
                {
                    int deltaStrip = minIsoStrip2 - maxIsoStrip1;
                    if(deltaStrip==0) continue;
                    for(int k=1; k<=abs(deltaStrip); k++)
                    {
                        int startIndex = maxIsoStrip1;
                        int newStripIndex = startIndex+k;
                        double newPointValue = myIsoStrips.at(newStripIndex).vmin;
                        double t = (newPointValue-valueOfPoint1)/(valueOfPoint2-valueOfPoint1);
                        double P1[3], P2[3];
                        this->pointCoord(P1,globalNodeID1);
                        this->pointCoord(P2,globalNodeID2);
                        double xP = P1[0] + t * (P2[0]-P1[0]);
                        double yP = P1[1] + t * (P2[1]-P1[1]);
                        double zP = P1[2] + t * (P2[2]-P1[2]);

                        cout<<"____case deltastrip>0: t = "<<t<<" newPoint("<<xP<<", "<<yP<<", "<<zP<<")____"<<endl;

                        isoStripPoint aP(xP,yP,zP,newPointValue);
                        aFaceTable.insertAfter(startIndex,sP,aP);
                        aFaceTable.insertBefore(startIndex+1,sP,aP);
                        sP = aP;
                    }
                }
                else
                {
                    int deltaStrip = minIsoStrip1 - maxIsoStrip2;
                    if(deltaStrip==0) continue;
                    for(int k=1; k<=abs(deltaStrip); k++)
                    {
                        int startIndex = maxIsoStrip2;
                        int newStripIndex = startIndex+k;
                        double newPointValue = myIsoStrips.at(newStripIndex).vmin;
                        double t = (newPointValue-valueOfPoint2)/(valueOfPoint1-valueOfPoint2);
                        double P1[3], P2[3];
                        this->pointCoord(P1,globalNodeID2);
                        this->pointCoord(P2,globalNodeID1);
                        double xP = P1[0] + t * (P2[0]-P1[0]);
                        double yP = P1[1] + t * (P2[1]-P1[1]);
                        double zP = P1[2] + t * (P2[2]-P1[2]);

                        cout<<"____case deltastrip<0: t = "<<t<<" newPoint("<<xP<<", "<<yP<<", "<<zP<<")____"<<endl;

                        isoStripPoint aP(xP,yP,zP,newPointValue);
                        //aFaceTable.insertAfter(startIndex,sP,aP);
                        //aFaceTable.insertBefore(startIndex+1,sP,aP);
                        aFaceTable.insertBefore(startIndex,sP,aP);
                        aFaceTable.insertAfter(startIndex+1,sP,aP);
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
void isoStripBuilder::getAllElements(const std::vector<faceTable> &vecFaces, std::vector<meshElementByCoords> &vecMeshElements)
{
    cout<<"isoStripBuilder::getAllElements()->____function called. Number of face tables: "<<vecFaces.size()<<"____"<<endl;
    for(std::vector<faceTable>::const_iterator it = vecFaces.cbegin(); it!=vecFaces.cend(); it++)
    {
        std::vector<meshElementByCoords> vecMeshElementsOfFace;
        it->getElements(vecMeshElementsOfFace);
        for(int i=0; i<vecMeshElementsOfFace.size(); i++)
            vecMeshElements.push_back(vecMeshElementsOfFace[i]);
    }
    cout<<"isoStripBuilder::getAllElements()->____number of elements: "<<vecMeshElements.size()<<"____"<<endl;
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
void isoStripBuilder::classifyNodes(QMap<int,int> &pointToIsoStrip)
{
    cout<<"isoStripBuilder::classifyNodes()->____function called____"<<endl;
    TColStd_MapIteratorOfPackedMapOfInteger it(myMeshDS->GetAllNodes());
    int NbNodes = myMeshDS->GetAllNodes().Extent();
    for(int localNodeID = 1; localNodeID<=NbNodes; localNodeID++,it.Next())
    {
        int globalNodeID = it.Key();
        //cout<<"____inserting "<<globalNodeID<<" into strip: ";
        double val = myValues.value(globalNodeID);
        for(int isoStripNb = 0; isoStripNb<myIsoStrips.size(); isoStripNb++)
        {
            isoStrip curIsoStrip = myIsoStrips.at(isoStripNb);
            if(curIsoStrip.contains(val))
            {
                pointToIsoStrip.insertMulti(globalNodeID,isoStripNb);
                //cout<<"("<<curIsoStrip.vmin<<","<<curIsoStrip.vmax<<")";
            }
        }
        //cout<<"\n"<<endl;
    }
    cout<<"isoStripBuilder::classifyNodes()->____exiting function____"<<endl;
}
