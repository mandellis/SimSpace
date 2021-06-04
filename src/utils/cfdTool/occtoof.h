#ifndef OCCTOOF_H
#define OCCTOOF_H

#include <stdio.h>
#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <filesystem>

#include <MeshVS_DataSource.hxx>
#include <ng_meshvs_datasourceface.h>
#include <occhandle.h>
#include "polygon.h"
#include "polyhedron.h"

namespace fs = std::experimental::filesystem;

namespace of {

    struct BB
    {
        char boundaryName[256];
        int startFace;
        int nFaces;
    };

    bool occToOF(const occHandle(MeshVS_DataSource) &sSMesh,
                  const std::map<std::string,std::map<int,occHandle(Ng_MeshVS_DataSourceFace)>> &allBound,
                  const std::string &workingPath);

    bool oFToOcc(const occHandle(MeshVS_DataSource) &sSMesh, const std::string &targetPath, const std::string &sourcePath);

    void printMark(std::ofstream &textFile);

    void printHeading(std::ofstream &textFile, std::string clas, std::string object, std::vector<int> &notes);

    void printEnd(std::ofstream &textFile);

    void printBoundaryIntro(std::ofstream &textFile, int meter, int second, std::vector<double> value);

    void printBoundaryIntro(std::ofstream &textFile, int meter, int second, double value);

    void printMshBoundary(std::ofstream &textFile, std::string name, int nFaces, int startFace);

    void printBoundary(std::ofstream &textFile, std::string name, std::string stype, std::vector<double> value);

    void printBoundary(std::ofstream &textFile, string name, std::string stype, std::string value);

    void closeBoundary(std::ofstream &textFile);

    /*void readBoundary(std::vector<of::BB> &vecBoundaries, ifstream &iS);

    void readPoints(std::vector<polygon::Point> &vecPoints, ifstream &iS);

    void readFaces(std::vector<polyhedron::face> &vecFaces,
                   std::vector<polygon::Point> &vecPoints,
                   std::vector<std::vector<int> > &vecFacesPointsIDs, ifstream &iS);

    void readCell(std::vector<polyhedron::cell> &vecCell,
                   std::vector<polyhedron::face> &vecFaces, std::vector<std::vector<int> > &vecCellFacesIDs,
                   ifstream &iS1, ifstream &iS2);

    void buildPointsToElementsConnectivity(std::multimap<int,int> &pointsElements,
                                           std::vector<std::vector<int>> &vecCellFacesIDs,
                                           std::vector<std::vector<int>> vecFacesPointsIDs);*/

    void readValue(std::vector<double> &cellPressure, ifstream &iS);

    void evaluatePointsValue(std::vector<std::map<int, double>> &pointsElementsWeights,
                                 std::vector<double> &cellPressure, std::vector<double> &pointPressure);

    void evaluatePointsWeights(std::vector<std::map<int, double>> &pointsElementsWeights,
                                   std::multimap<int,int> &pointsElements, std::vector<polyhedron::cell> &vecCell,
                                   std::vector<polygon::Point> vecPoints);

    void evaluateAndPrintResults(std::vector<std::map<int, double>> &pointsElementsWeights,
                                     fs::path sourcePath,fs::path targetPath, std::string folder);


}
#endif // OCCTOOF_H
