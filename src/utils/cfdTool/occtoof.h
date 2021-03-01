#include <stdio.h>
#include <vector>
#include <map>
#include <string>

#include <MeshVS_DataSource.hxx>
#include <ng_meshvs_datasourceface.h>
#include <occhandle.h>

namespace opFo {

    bool occToOF (const occHandle(MeshVS_DataSource) &sSMesh,
                  const std::map<std::string,std::map<int,occHandle(Ng_MeshVS_DataSourceFace)>> &allBound,
                  const std::string &workingPath);

    void printMark(std::ofstream &textFile);

    void printHeading(std::ofstream &textFile, std::string classe, std::string object, std::vector<int> &notes);

    void printEnd(std::ofstream &textFile);

    void printBoundary(std::ofstream &textFile, std::string name, int nFaces, int startFace);

}

/* simulationmanager(generateboundaryconditiondatasource)
 * datasourcebuilder*/
