#include "occtoof.h"
#include "elementtypes.h"
#include "meshelement2d.h"
#include "ng_meshvs_datasource3d.h"

//! ---
//! Std
//! ---
#include <set>
#include <map>
#include <fstream>
#include <string>

//! ---
//! Qt
//! ---
#include <QList>
#include <QMap>
#include <QString>

//! ----
//! OCC
//! ----
#include <TColStd_PackedMapOfInteger.hxx>
#include <TColStd_MapIteratorOfPackedMapOfInteger.hxx>


namespace of {

    bool occToOF(const occHandle(MeshVS_DataSource) &sSMesh,
                 const std::map<std::string,std::map<int,occHandle(Ng_MeshVS_DataSourceFace)>> &allBound,
                 const std::string &workingPath)
    {
        if(!sSMesh.IsNull())
        {
            //!initialize variable for function GetGeom

            const TColStd_PackedMapOfInteger &aNodes = sSMesh->GetAllNodes();
            double buf1[3];
            TColStd_Array1OfReal coords(*buf1,1,3);
            int nbNodes1;
            MeshVS_EntityType aType;

            //!open points ofstream

            const std::string pointsPath = workingPath + "/constant/polyMesh/points";
            std::ofstream points;
            points.open(pointsPath.c_str());
            if(!points.is_open()) return false;

            //!print points file heading

            of::printMark(points);
            of::printHeading(points, std::string("vectorField"), std::string("points"), std::vector<int>());
            points<<aNodes.Extent()<<endl;
            points<<"("<<endl;

            //!print points file

            for (TColStd_MapIteratorOfPackedMapOfInteger anIter(aNodes); anIter.More();anIter.Next())
            {
                int globalNodeID = anIter.Key();
                if (!sSMesh->GetGeom(globalNodeID,false,coords,nbNodes1,aType)) return false;
                points<<"("<<coords(1)<<" "<<coords(2)<<" "<<coords(3)<<")"<<endl;
            }

            of::printEnd(points);

            //!open faces, owner, neighbour ofstream

            const std::string facesPath = workingPath + "/constant/polyMesh/faces";
            const std::string ownerPath = workingPath + "/constant/polyMesh/owner";
            const std::string neighbourPath = workingPath + "/constant/polyMesh/neighbour";
            const std::string boundaryPath = workingPath + "/constant/polyMesh/boundary";
            std::ofstream faces, owner,neighbour, boundary;
            faces.open(facesPath.c_str());
            owner.open(ownerPath.c_str());
            neighbour.open(neighbourPath.c_str());
            boundary.open(boundaryPath.c_str());
            if(!faces.is_open() || !owner.is_open() || !neighbour.is_open() || !boundary.is_open()) return false;

            //!create a set of surface elements

            const occHandle(Ng_MeshVS_DataSource3D) &volumeSSMesh =
                    occHandle(Ng_MeshVS_DataSource3D)::DownCast(sSMesh);
            if (volumeSSMesh->myFaceToElements.isEmpty())
                volumeSSMesh->buildFaceToElementConnectivity();
            std::set<meshElement2D> surFaces;
            for (int k=0; k<volumeSSMesh->mySurfaceElements.size(); k++) {
                surFaces.insert(volumeSSMesh->mySurfaceElements.at(k));
            }

            //!print faces, owner and neighbour heading

            std::vector<int> notes;
            const TColStd_PackedMapOfInteger &aElem = sSMesh->GetAllElements();
            notes.push_back(aNodes.Extent());
            notes.push_back(aElem.Extent());
            notes.push_back(volumeSSMesh->myFaceToElements.size());
            notes.push_back(volumeSSMesh->myFaceToElements.size()- int(volumeSSMesh->mySurfaceElements.size()));

            of::printMark(faces);
            of::printHeading(faces, std::string("faceList"), std::string("faces"), std::vector<int>());
            faces<<notes.at(2)<<endl;
            faces<<"("<<endl;
            of::printMark(owner);
            of::printHeading(owner, std::string("labelList"), std::string("owner"), notes);
            owner<<notes.at(2)<<endl;
            owner<<"("<<endl;
            of::printMark(neighbour);
            of::printHeading(neighbour, std::string("labelList"), std::string("owner"), notes);
            neighbour<<notes.at(3)<<endl;
            neighbour<<"("<<endl;
            of::printMark(boundary);
            of::printHeading(boundary, std::string("polyBoundaryMesh"), std::string("boundary"), std::vector<int>());
            boundary<<allBound.size()<<endl;
            boundary<<"("<<endl;

            //!print internal surfaces in faces, owner and neighbour files

            QMap<meshElement2D, QList<int>>::const_iterator it1;
            for (it1 = volumeSSMesh->myFaceToElements.begin(); it1!=volumeSSMesh->myFaceToElements.end(); ++it1) {
                if(surFaces.count(it1.key())==0) {
                    const meshElement2D &tempMeshElement1 = it1.key();
                    int nbNodes2 = tempMeshElement1.nodeIDs.size();
                    faces<<nbNodes2<<"(";
                    QList<int>::const_iterator qIt1;
                    for (qIt1 = tempMeshElement1.nodeIDs.begin(); qIt1!=tempMeshElement1.nodeIDs.end(); ++qIt1)
                        faces<<*qIt1-1<<" ";
                    faces<<")"<<endl;
                    const QList<int> &ownerNeightbours = it1.value();
                    int st = ownerNeightbours.at(0)-1;
                    int nd = ownerNeightbours.at(1)-1;
                    if (st<nd) {owner<<st<<endl; neighbour<<nd<<endl;}
                    else if (nd<st) {owner<<nd<<endl; neighbour<<st<<endl;}
                    else return false;
                }
            }

            //!create a set of boundary surface elements

            std::map<std::string,std::set<meshElement2D>> mapOfBoundarySets;
            for(std::map<std::string,std::map<int,occHandle(Ng_MeshVS_DataSourceFace)>>::const_iterator
                it2 = allBound.cbegin(); it2!=allBound.cend(); it2++)
            {
                const std::string &boundaryName = it2->first;
                const std::map<int,occHandle(Ng_MeshVS_DataSourceFace)> &mapOfBoundaryMeshDS = it2->second;
                std::set<meshElement2D> setOfMeshElement2D;
                for(std::map<int,occHandle(Ng_MeshVS_DataSourceFace)>::const_iterator
                    it3 = mapOfBoundaryMeshDS.cbegin(); it3!=mapOfBoundaryMeshDS.cend(); it3++)
                {
                    const occHandle(Ng_MeshVS_DataSourceFace) &curMeshDSFace = it3->second;
                    for(TColStd_MapIteratorOfPackedMapOfInteger it4(curMeshDSFace->GetAllElements());
                        it4.More(); it4.Next())
                    {
                        int globalElementID = it4.Key();
                        int nbNodes3, buf[20];
                        TColStd_Array1OfInteger nodeIDs(*buf,1,20);
                        if (!curMeshDSFace->GetNodesByElement(globalElementID,nodeIDs,nbNodes3)) return false;
                        meshElement2D aMeshElement2D;
                        for(int p=1; p<=nbNodes3; p++) aMeshElement2D.nodeIDs.push_back(nodeIDs.Value(p));
                        setOfMeshElement2D.insert(aMeshElement2D);
                    }
                }
                mapOfBoundarySets.insert(std::make_pair(boundaryName, setOfMeshElement2D));

            }

            //!print boundary surfaces in faces, owner and boundary files

            int gap = notes.at(3);

            for(std::map<std::string,std::set<meshElement2D>>::const_iterator
                it5 = mapOfBoundarySets.cbegin(); it5!=mapOfBoundarySets.cend(); it5++)
            {
                const std::string &bName = it5->first;
                const std::set<meshElement2D> &boundarySet = it5->second;
                for (std::set<meshElement2D>::const_iterator it6 = boundarySet.begin(); it6!=boundarySet.end(); ++it6)
                {
                    meshElement2D tempMeshElement2 = *it6;
                    QList<int>::const_iterator qIt2;
                    faces<<tempMeshElement2.nodeIDs.size()<<"(";
                    for (qIt2 = tempMeshElement2.nodeIDs.begin(); qIt2 != tempMeshElement2.nodeIDs.end(); ++qIt2)
                        faces<<*qIt2-1<<" ";
                    faces<<")"<<endl;
                    const QList<int> &surfaceOwner = volumeSSMesh->myFaceToElements.value(tempMeshElement2);
                    owner<<surfaceOwner.at(0)-1<<endl;
                }
                printMshBoundary(boundary, bName, int(boundarySet.size()), gap);
                gap += int(boundarySet.size());
            }


            /*for (std::set<meshElement2D>::iterator it2=surFaces.begin(); it2!=surFaces.end(); ++it2) {
                meshElement2D tempMeshElement2 = *it2;
                int nbNodes4 = tempMeshElement2.nodeIDs.size();
                faces<<nbNodes4<<"(";
                QList<int>::iterator qIt2;
                for (qIt2 = tempMeshElement2.nodeIDs.begin(); qIt2 != tempMeshElement2.nodeIDs.end(); ++qIt2)
                    faces<<*qIt2-1<<" ";
                faces<<")"<<endl;
                QList<int> surfaceOwner = volumeSSMesh->myFaceToElements.value(tempMeshElement2);
                owner<<surfaceOwner.at(0)-1<<endl;
            }*/


            //!print the end lines in faces, owner and neighbour files

            of::printEnd(faces);
            of::printEnd(owner);
            of::printEnd(neighbour);
            of::printEnd(boundary);

        }

        return true;
    }



    void printMark(std::ofstream &textFile)
    {
        if (!textFile.is_open()) return;
        textFile<<"/*--------------------------------*-  C++ -*----------------------------------*\\"<<endl;
        textFile<<"| =========                   |                                                 |"<<endl;
        textFile<<"| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |"<<endl;
        textFile<<"|  \\\\    /   O peration     | Version:  v2012                                 |"<<endl;
        textFile<<"|   \\\\  /    A nd           | Website:  www.openfoam.com                      |"<<endl;
        textFile<<"|    \\\\/     M anipulation  |                                                 |"<<endl;
        textFile<<"\\*----------------------------------------------------------------------------*/"<<endl;
        return;
    }

    void printHeading(std::ofstream &textFile, std::string clas, std::string object, std::vector<int> &notes)
    {
        if (!textFile.is_open()) return;
        textFile<<"FoamFile"<<endl;
        textFile<<"{"<<endl;
        textFile<<"    version     2.0;"<<endl;
        textFile<<"    format      ascii;"<<endl;
        textFile<<"    class       "<<clas<<";"<<endl;
        if (notes.size()>0) textFile<<"    note        "<<"\"nPoints:"<<notes.at(0)<<"  nCells:"<<notes.at(1)
               <<"  nFaces:"<<notes.at(2)<<"  nInternalFaces:"<<notes.at(3)<<"\";"<<endl;
        textFile<<"    location    \"constant/polyMesh\";"<<endl;
        textFile<<"    object      "<<object<<";"<<endl;
        textFile<<"}"<<endl;
        textFile<<"// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //"<<endl;
        textFile<<endl;
        textFile<<endl;
        return;
    }

    void printEnd(std::ofstream &textFile)
     {
         if (!textFile.is_open()) return;
         textFile<<")"<<endl;
         textFile<<endl;
         textFile<<endl;
         textFile<<"// ************************************************************************* //"<<endl;
         return;
     }

    void printMshBoundary(std::ofstream &textFile, std::string name, int nFaces, int startFace)
    {
        if (!textFile.is_open()) return;
        textFile<<"    "<<name<<endl;
        textFile<<"    {"<<endl;
        QString qName = QString::fromStdString(name);
        if (qName.startsWith("inlet") || qName.startsWith("outlet")) {
            textFile<<"        type            patch;"<<endl;
        } else if (qName.startsWith("wall")) {
            textFile<<"        type            wall;"<<endl;
            textFile<<"        inGroups        1(wall);"<<endl;
        } else if (qName.startsWith("simmetry")) {
            textFile<<"        type            symmetryPlane;"<<endl;
            textFile<<"        inGroups        1(symmetryPlane);"<<endl;
        } else {
            textFile<<"        type            empty;"<<endl;
            textFile<<"        inGroups        1(empty);"<<endl;
        }
        textFile<<"        nFaces          "<<nFaces<<";"<<endl;
        textFile<<"        startFace       "<<startFace<<";"<<endl;
        textFile<<"    }"<<endl;
        return;
    }
    void printBoundary(std::ofstream &textFile, std::string stype, std::vector<double> value)
    {
        if (!textFile.is_open()) return;
        textFile<<"    "<<stype<<endl;
        textFile<<"    {"<<endl;
        textFile<<" type    "<<stype<<endl;
        if(stype=="fixedValue")
        textFile<<" value       uniform (";
        if(value.size()>1)
                  textFile<<value.at(0)<<" "<<value.at(1)<<" "<<value.at(2)<<endl;
        else textFile<<value.at(0)<<");"<<endl;
        textFile<<"    }"<<endl;
        return;
    }
    void printBoundary(std::ofstream &textFile, std::string stype, std::string value)
    {
        if (!textFile.is_open()) return;
        textFile<<"    "<<stype<<endl;
        textFile<<"    {"<<endl;
        textFile<<" type    "<<stype<<endl;
        if(stype=="fixedValue")
        textFile<<" value       "<<value<<endl;
        textFile<<"    }"<<endl;
        return;
    }

    void closeBoundary(std::ofstream &textFile)
    {
        if (!textFile.is_open()) return;
        textFile<<"}"<<endl;
        textFile<<endl;
        textFile<<"// ************************************************************************* //"<<endl;
        return;
    }
}
