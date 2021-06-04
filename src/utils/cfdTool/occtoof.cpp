//! ----
//! Custom
//! ----
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
#include <filesystem>

//! ---
//! Qt
//! ---
#include <QList>
#include <QMap>
#include <QString>
#include <QDir>

//! ----
//! OCC
//! ----
#include <TColStd_PackedMapOfInteger.hxx>
#include <TColStd_MapIteratorOfPackedMapOfInteger.hxx>

namespace fs = std::experimental::filesystem;

bool of::occToOF(const occHandle(MeshVS_DataSource) &sSMesh,
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

        const std::string polyMeshPath = workingPath + "/constant/polyMesh";
        fs::path myPolyMeshPath = fs::path(polyMeshPath);

        if (fs::exists(myPolyMeshPath)) {
            myPolyMeshPath.clear();
        } else {
            fs::create_directories(myPolyMeshPath);
        }

        const std::string pointsPath = workingPath + "/constant/polyMesh/points";
        std::ofstream points;
        points.open(pointsPath.c_str());
        if(!points.is_open()) return false;

        //!print points file heading

        of::printMark(points);
        of::printHeading(points, std::string("vectorField"), std::string("points"), std::vector<int>());
        points<<aNodes.Extent()<<endl;
        points<<"("<<endl;

        //!vector of points

        std::vector<polygon::Point> vecPoints;

        //!print points file

        for (TColStd_MapIteratorOfPackedMapOfInteger anIter(aNodes); anIter.More();anIter.Next())
        {
            int globalNodeID = anIter.Key();
            if (!sSMesh->GetGeom(globalNodeID,false,coords,nbNodes1,aType)) return false;
            points<<"("<<coords(1)<<" "<<coords(2)<<" "<<coords(3)<<")"<<endl;
            vecPoints.push_back(polygon::Point(coords(1), coords(2), coords(3)));
        }

        of::printEnd(points);

        //!open faces, owner, neighbour, boundary ofstream

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

        //!build multimap of faces and "owner&neighbour"

        std::multimap<int, std::pair<int, vector<int>>> ownNeighFacesMap;
        QMap<meshElement2D, QList<int>>::const_iterator it1;
        for (it1 = volumeSSMesh->myFaceToElements.begin(); it1!=volumeSSMesh->myFaceToElements.end(); ++it1) {
            if(surFaces.count(it1.key())==0) {
                const meshElement2D &tempMeshElement1 = it1.key();
                QList<int>::const_iterator qIt1;
                vector<int> facePoints;
                for (qIt1 = tempMeshElement1.nodeIDs.begin(); qIt1!=tempMeshElement1.nodeIDs.end(); ++qIt1) {
                    facePoints.push_back(*qIt1-1);
                }
                const QList<int> &ownerNeightbours = it1.value();
                vector<int> ownNeigh;
                int st = ownerNeightbours.at(0)-1;
                int nd = ownerNeightbours.at(1)-1;
                if (st<nd) {
                    ownNeighFacesMap.insert(std::make_pair(st, std::make_pair(nd, facePoints)));
                } else if (nd<st) {
                    ownNeighFacesMap.insert(std::make_pair(nd, std::make_pair(st, facePoints)));
                } else return false;
            }
        }

        //!create cell vector and pointsElements map
        //!(this one can have duplicates points with elements connection)

        std::vector<polyhedron::cell> vecCell;
        for (int w=0; w<notes.at(1); w++) vecCell.push_back(polyhedron::cell());
        std::multimap<int, int> pointsElements;

        //!print internal surfaces in faces, owner and neighbour files
        //!create cell vector.

        for(int j=0; j<notes.at(1); j++)
        {
            std::pair<std::multimap<int, std::pair<int, vector<int>>>::iterator,
                      std::multimap<int, std::pair<int, vector<int>>>::iterator> range;
            range = ownNeighFacesMap.equal_range(j);
            for(std::multimap<int, std::pair<int, vector<int>>>::iterator it5 = range.first; it5 != range.second; ++it5)
            {
                owner<<j<<endl;
                neighbour<<it5->second.first<<endl;
                std::vector<int> facce = it5->second.second;
                polyhedron::face polyFace1;
                faces<<facce.size()<<"(";
                for(int k=0; k<facce.size(); k++)
                {
                    faces<<facce.at(k)<<" ";
                    polyFace1.push_back(vecPoints.at(facce.at(k)));
                    pointsElements.insert(std::make_pair(facce.at(k), j));
                    pointsElements.insert(std::make_pair(facce.at(k), it5->second.first));
                }
                faces<<")"<<endl;
                vecCell.at(j).addFace(polyFace1);
                vecCell.at(it5->second.first).addFace(polyFace1);
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
                const QList<int> &surfaceOwner = volumeSSMesh->myFaceToElements.value(tempMeshElement2);
                polyhedron::face polyFace2;
                faces<<tempMeshElement2.nodeIDs.size()<<"(";
                for (qIt2 = tempMeshElement2.nodeIDs.begin(); qIt2 != tempMeshElement2.nodeIDs.end(); ++qIt2)
                {
                    faces<<*qIt2-1<<" ";
                    polyFace2.push_back(vecPoints.at(*qIt2-1));
                    pointsElements.insert(std::make_pair(*qIt2-1, surfaceOwner.at(0)-1));
                }
                faces<<")"<<endl;
                vecCell.at(surfaceOwner.at(0)-1).addFace(polyFace2);
                owner<<surfaceOwner.at(0)-1<<endl;
            }
            of::printMshBoundary(boundary, bName, int(boundarySet.size()), gap);
            gap += int(boundarySet.size());
        }

        //!print the end lines in faces, owner and neighbour files

        of::printEnd(faces);
        of::printEnd(owner);
        of::printEnd(neighbour);
        of::printEnd(boundary);

        std::vector<std::map<int, double>> pointsElementsWeights;
        of::evaluatePointsWeights(pointsElementsWeights, pointsElements, vecCell, vecPoints);
        volumeSSMesh->pointsElementsWeights = pointsElementsWeights;
    }

    return true;
}

bool of::oFToOcc(const opencascade::handle<MeshVS_DataSource> &sSMesh,
                 const std::string &targetPath, const std::string &sourcePath)
{
    //!initialize directories

    fs::path myTargetPath = fs::path(targetPath);
    fs::path mySourcePath = fs::path(sourcePath);

    if (fs::exists(myTargetPath)) {
        myTargetPath.clear();
    } else {
        fs::create_directories(myTargetPath);
    }

    std::vector<std::string> folders;
    for (auto &it5: fs::directory_iterator(mySourcePath))
    {
        std::string fold = it5.path().stem().string();
        bool isDir = fs::is_directory(it5.path());
        bool isConst = fold == "constant";
        bool isSys = fold == "system";
        if(isDir && !isConst && !isSys) folders.push_back(fold);
    }

    //last time, all time, one time

    const occHandle(Ng_MeshVS_DataSource3D) &volumeSSMesh =
            occHandle(Ng_MeshVS_DataSource3D)::DownCast(sSMesh);

    std::vector<std::map<int, double>> pointsElementsWeights = volumeSSMesh->pointsElementsWeights;

    //come gestisco la scelta dei tempi? DOMANDA
    std::string mode = "last time";
    std::string time = "20";

    if (mode == "last time")
    {
        std::string step = folders.back();
        evaluateAndPrintResults(pointsElementsWeights, mySourcePath, myTargetPath, step);
    }
    else if (mode == "all time")
    {
        for(int i=1; i<folders.size(); i++)
        {
            evaluateAndPrintResults(pointsElementsWeights, mySourcePath, myTargetPath, folders.at(i));
        }
    }
    else if (mode == "time point")
    {
        evaluateAndPrintResults(pointsElementsWeights, mySourcePath, myTargetPath, time);
    }

    return true;
}

void of::printMark(std::ofstream &textFile)
{
    if (!textFile.is_open()) return;
    textFile<<"/*--------------------------------*- C++ -*----------------------------------*\\"<<endl;
    textFile<<"| =========                 |                                                 |"<<endl;
    textFile<<"| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |"<<endl;
    textFile<<"|  \\\\    /   O peration     | Version:  v2012                                 |"<<endl;
    textFile<<"|   \\\\  /    A nd           | Website:  www.openfoam.com                      |"<<endl;
    textFile<<"|    \\\\/     M anipulation  |                                                 |"<<endl;
    textFile<<"\\*---------------------------------------------------------------------------*/"<<endl;
    return;
}

void of::printHeading(std::ofstream &textFile, std::string clas, std::string object, std::vector<int> &notes)
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

void of::printEnd(std::ofstream &textFile)
 {
     if (!textFile.is_open()) return;
     textFile<<")"<<endl;
     textFile<<endl;
     textFile<<endl;
     textFile<<"// ************************************************************************* //"<<endl;
     return;
 }

void of::printMshBoundary(std::ofstream &textFile, std::string name, int nFaces, int startFace)
{
    if (!textFile.is_open()) return;
    textFile<<"    "<<name<<endl;
    textFile<<"    {"<<endl;
    QString qName = QString::fromStdString(name);
    if (qName.startsWith("velocity") || qName.startsWith("pressure")) {
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

void of::printBoundaryIntro(std::ofstream &textFile, int meter, int second, std::vector<double> value)
{
    if (!textFile.is_open()) return;
    textFile<<endl;
    textFile<<"dimensions      [0 "<<meter<<" "<<second<<" 0 0 0 0];"<<endl;
    textFile<<endl;
    textFile<<"internalField   uniform ("<<value.at(0)<<" "<<value.at(1)<<" "<<value.at(2)<<");"<<endl;
    textFile<<endl;
    textFile<<"boundaryField"<<endl;
    textFile<<"{"<<endl;
    return;
}

void of::printBoundaryIntro(std::ofstream &textFile, int meter, int second, double value)
{
    if (!textFile.is_open()) return;
    textFile<<endl;
    textFile<<"dimensions      [0 "<<meter<<" "<<second<<" 0 0 0 0];"<<endl;
    textFile<<endl;
    textFile<<"internalField   uniform "<<value<<";"<<endl;
    textFile<<endl;
    textFile<<"boundaryField"<<endl;
    textFile<<"{"<<endl;
    return;
}

void of::printBoundary(std::ofstream &textFile, std::string name, std::string stype, std::vector<double> value)
{
    if (!textFile.is_open()) return;
    textFile<<"    "<<name<<endl;
    textFile<<"    {"<<endl;
    textFile<<"        type            "<<stype<<";"<<endl;
    if(stype=="fixedValue")
    {
        textFile<<"        value           uniform ";
        if (value.size()>1)
        {
            textFile<<"("<<value.at(0)<<" "<<value.at(1)<<" "<<value.at(2)<<");"<<endl;
        }
        else
        {
            textFile<<value.at(0)<<";"<<endl;
        }
    }
    textFile<<"    }"<<endl;
    textFile<<endl;
    return;
}

void of::printBoundary(std::ofstream &textFile, std::string name, std::string stype, std::string value)
{
    if (!textFile.is_open()) return;
    textFile<<"    "<<name<<endl;
    textFile<<"    {"<<endl;
    textFile<<"        type            "<<stype<<";"<<endl;
    if(stype=="fixedValue") textFile<<"        value           "<<value<<";"<<endl;
    textFile<<"    }"<<endl;
    textFile<<endl;
    return;
}

void of::closeBoundary(std::ofstream &textFile)
{
    if (!textFile.is_open()) return;
    textFile<<"}"<<endl;
    textFile<<endl;
    textFile<<"// ************************************************************************* //"<<endl;
    return;
}

void of::evaluatePointsWeights(std::vector<std::map<int, double>> &pointsElementsWeights,
                               std::multimap<int,int> &pointsElements, std::vector<polyhedron::cell> &vecCell,
                               std::vector<polygon::Point> vecPoints)
{
    std::vector<polygon::Point> centroids;
    for(int j=0; j<vecCell.size(); j++)
    {
        polygon::Point p;
        vecCell.at(j).getCentroid(p.x, p.y, p.z);
        centroids.push_back(p);
    }
    for(int i=0; i<vecPoints.size(); i++)
    {
        std::pair<std::multimap<int, int>::iterator, std::multimap<int, int>::iterator> range;
        range = pointsElements.equal_range(i);
        std::set<int> elementi;
        for(std::multimap<int, int>::iterator it1 = range.first; it1 != range.second; ++it1)
        {
            int elID = it1->second;
            elementi.insert(elID);
        }
        std::map<int, double> elemCentrDist;
        double total = 0;
        polygon::Point myP = vecPoints.at(i);
        for(std::set<int>::iterator it4 = elementi.begin(); it4 != elementi.end(); ++it4)
        {
            int elID = *it4;
            polygon::Point myC = centroids.at(elID);
            double xdist = myP.x - myC.x;
            double ydist = myP.y - myC.y;
            double zdist = myP.z - myC.z;
            double distsum = std::pow(xdist,2) + std::pow(ydist,2) + std::pow(zdist,2);
            double dist = std::sqrt(distsum);
            elemCentrDist.insert(std::make_pair(elID, dist));
            total += dist;
        }
        std::map<int,double> elemAbsoluteWeights;
        double absTotal = 0;
        std::map<int,double>::iterator it2;
        for(it2 = elemCentrDist.begin(); it2 != elemCentrDist.end(); ++it2)
        {
            int elID = it2->first;
            double distance = it2->second;
            double absWeight = total/distance;
            elemAbsoluteWeights.insert(std::make_pair(elID, absWeight));
            absTotal += absWeight;
        }

        std::map<int,double> elemWeights;
        std::map<int,double>::iterator it3;
        for(it3 = elemAbsoluteWeights.begin(); it3 != elemAbsoluteWeights.end(); ++it3)
        {
            int elID = it3->first;
            double absW = it3->second;
            double weight = absW/absTotal;
            elemWeights.insert(std::make_pair(elID, weight));
        }
        pointsElementsWeights.push_back(elemWeights);
    }
}

void of::evaluateAndPrintResults(std::vector<std::map<int, double>> &pointsElementsWeights,
                                 fs::path sourcePath,fs::path targetPath, std::string folder)
{
    fs::path myPath = sourcePath / folder;

    if (!fs::is_directory(myPath))
    {
        cout<<"occtoof::oFToOcc()->____\"0\" directory not found____"<<endl;
        return;
    }

    std::ifstream iS;
    iS.open(myPath.string());
    if(!iS.is_open())
    {
        cout<<"occtoof::evaluateAndPrintResults()->____input stream isn't open____"<<endl;
        return;
    }

    std::vector<double> cellPressure;
    std::vector<double> pointPressure;
    of::readValue(cellPressure, iS);
    iS.close();

    of::evaluatePointsValue(pointsElementsWeights, cellPressure, pointPressure);

    fs::path finalPath = targetPath / folder / "p";
    std::ofstream myfile;
    myfile.open(finalPath.c_str());
    if(!myfile.is_open())
    {
        cout<<"occtoof::evaluateAndPrintResults()->____output stream isn't open____"<<endl;
        return;
    }
    myfile<<pointPressure.size()<<endl;
    myfile<<"("<<endl;
    for(int i=0; i<pointPressure.size(); i++)
    {
        myfile<<pointPressure.at(i)<<endl;
    }
    myfile<<")"<<endl;

}

void of::readValue(std::vector<double> &cellPressure, ifstream &iS)
{
    //! ----------------
    //! skip the header
    //! ----------------

    std::string val;
    std::getline(iS,val);

    while(val.compare("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //")!=0)
        std::getline(iS,val);

    //! ----------------
    //! get pressure type
    //! ----------------

    std::getline(iS,val);
    while(val.empty()) std::getline(iS,val);
    std::getline(iS,val);
    while(val.empty()) std::getline(iS,val);
    char str1[128], str2[128];
    sscanf(val.c_str(),"%s%s",str1, str2);

    //! ----------------
    //! fill vector of pressure
    //! ----------------

    if(std::strcmp("nonuniform", str2)==0)
    {
        std::getline(iS,val);
        while(val.empty()) std::getline(iS,val);
        int eNumber;
        sscanf(val.c_str(),"%d",&eNumber);
        std::getline(iS,val);

        double pressure;
        for(int k=0; k<eNumber; k++)
        {
            std::getline(iS,val);
            sscanf(val.c_str(),"%lf", &pressure);
            cellPressure.push_back(pressure);
        }
    }
}

void of::evaluatePointsValue(std::vector<std::map<int, double>> &pointsElementsWeights,
                             std::vector<double> &cellPressure, std::vector<double> &pointPressure)
{
    for(int i=0; i<pointsElementsWeights.size(); i++)
    {
        std::map<int, double> nodeElemMap = pointsElementsWeights.at(i);
        double pointP = 0;
        std::map<int, double>::iterator it1;
        for (it1 = nodeElemMap.begin(); it1 != nodeElemMap.end(); ++it1)
        {
            int elID = it1->first;
            double weight = it1->second;
            pointP += weight * cellPressure.at(elID);
        }
        pointPressure.push_back(pointP);
    }
}
