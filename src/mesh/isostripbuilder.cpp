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

//! ------
//! Eigen
//! ------
#include <Eigen/Dense>

//! ----
//! C++
//! ----
#include <map>
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
isoStripBuilder::isoStripBuilder(const occHandle(MeshVS_DataSource) &aMeshDS, int NbStrips, QObject *parent):
    myMeshDS(aMeshDS),
    myNbStrips(NbStrips),
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

//! ----------------------
//! function: setNbStrips
//! details:
//! ----------------------
void isoStripBuilder::setNbStrip(int NbStrips)
{
    myNbStrips = NbStrips;
}

//! -----------------------
//! function: setIsoStrips
//! details:
//! -----------------------
void isoStripBuilder::setIsoStrips(std::vector<isoStrip> *theIsoStrips)
{
    myIsoStrips = theIsoStrips;
}


//! ------------------
//! function: perform
//! details:
//! ------------------
bool isoStripBuilder::perform()
{
    if(myNbStrips==0) return false;
    if(myMeshDS.IsNull()) return false;
    if(myValues==Q_NULLPTR) return false;
    if(myValues->isEmpty())
    {
        cerr<<"soStripBuilder::perform()->____empty results____"<<endl;
        return false;
    }

    //! ---------------------------------------------------------
    //! classify the nodes according to the isostrip definitions
    //! ---------------------------------------------------------
    QMap<int,int> nodeToIsoStrip;
    this->classifyNodes(nodeToIsoStrip);

    //! ---------------------------
    //! process faces and segments
    //! ---------------------------
    std::vector<std::vector<int>> faceTable;
    faceTable.resize(myNbStrips);

    MIPMI it(myMeshDS->GetAllElements());
    int NbElements = myMeshDS->GetAllElements().Extent();
    for(int localElementID = 1; localElementID<=NbElements; localElementID++,it.Next())
    {
        int globalElementID = it.Key();
        int NbNodes;
        occHandle(MeshVS_HArray1OfSequenceOfInteger) topology;
        myMeshDS->Get3DGeom(globalElementID,NbNodes,topology);

        int NbNodesOfElement, b[20];
        TColStd_Array1OfInteger nodeIDs(*b,1,20);
        myMeshDS->GetNodesByElement(globalElementID,nodeIDs,NbNodesOfElement);
        int NbFaces = topology->Length();
        for(int f = 0; f<NbFaces; f++)
        {
            //! ----------------------------------------
            //! the sequence of nodes defining the face
            //! ----------------------------------------
            TColStd_SequenceOfInteger aSeq = topology->Value(f);
            for(int i = 0; i<NbNodes; i++)
            {
                //! -------------------------------------------
                //! put the first point into the right columns
                //! -------------------------------------------
                int index1 = aSeq.Value(i%NbNodes)+1;
                int globalNodeID1 = nodeIDs(index1);
                const QList<int> &isoStripOf1 = nodeToIsoStrip.values(globalNodeID1);
                for(int i=0; i<isoStripOf1.length(); i++)
                {
                    int curIsoStripNb = isoStripOf1.at(i);
                    faceTable[curIsoStripNb].push_back(globalNodeID1);
                }



                //! --------------------------------------------
                //! put the second point into the right columns
                //! --------------------------------------------
                int index2 = aSeq.Value((i+1)%NbNodes)+1;
                int globalNodeID2 = nodeIDs(index2);
                const QList<int> &isoStripOf2 = nodeToIsoStrip.values(globalNodeID2);
                for(int i=0; i<isoStripOf1.length(); i++)
                {
                    int curIsoStripNb = isoStripOf2.at(i);
                    faceTable[curIsoStripNb].push_back(globalNodeID2);
                }

                //! ----------------------------
                //! point coordinates if needed
                //! ----------------------------
                MeshVS_EntityType aType;
                int NbNodes1;
                double b1[3];
                TColStd_Array1OfReal cn1(*b1,1,3);
                myMeshDS->GetGeom(globalNodeID1, false, cn1, NbNodes1, aType);
                TColStd_Array1OfReal cn2(*b1,1,3);
                myMeshDS->GetGeom(globalNodeID1, false, cn2, NbNodes1, aType);



            }
        }

    }
    return true;
}


//! ------------------------
//! function: classifyNodes
//! details:
//! ------------------------
void isoStripBuilder::classifyNodes(QMap<int,int> &nodeToIsoStrip)
{
    MIPMI it(myMeshDS->GetAllNodes());
    int NbNodes = myMeshDS->GetAllNodes().Extent();
    for(int localNodeID = 1; localNodeID<=NbNodes; localNodeID++,it.Next())
    {
        int globalNodeID = it.Key();
        double val = myValues->value(globalNodeID);
        for(int isoStripNb = 0; isoStripNb<myIsoStrips->size(); isoStripNb++)
        {
            isoStrip curIsoStrip = myIsoStrips->at(isoStripNb);
            if(curIsoStrip.contains(val)) nodeToIsoStrip.insertMulti(globalNodeID,isoStripNb);
        }
    }
}
