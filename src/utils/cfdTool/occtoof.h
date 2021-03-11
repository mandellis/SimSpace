#include <stdio.h>
#include <vector>
#include <map>
#include <string>

#include <MeshVS_DataSource.hxx>
#include <ng_meshvs_datasourceface.h>
#include <occhandle.h>

namespace of {

    bool occToOF (const occHandle(MeshVS_DataSource) &sSMesh,
                  const std::map<std::string,std::map<int,occHandle(Ng_MeshVS_DataSourceFace)>> &allBound,
                  const std::string &workingPath); 

    void printMark(std::ofstream &textFile);

    void printHeading(std::ofstream &textFile, std::string clas, std::string object, std::vector<int> &notes);

    void printEnd(std::ofstream &textFile);

    void printMshBoundary(std::ofstream &textFile, std::string name, int nFaces, int startFace);

    void printBoundary(std::ofstream &textFile, std::string stype, std::vector<double> value);
    void printBoundary(std::ofstream &textFile, std::string stype, std::string value);

    void closeBoundary(std::ofstream &textFile);

}

/* simulationmanager(generateboundaryconditiondatasource)
 * datasourcebuilder*/
