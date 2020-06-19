#ifndef MESHSELECTNLAYERS_H
#define MESHSELECTNLAYERS_H

#include "occhandle.h"
#include <ng_meshvs_datasource3d.h>
#include <ng_meshvs_datasourceface.h>

class meshSelectNLayers
{
public:

    meshSelectNLayers(const occHandle(Ng_MeshVS_DataSource3D) &aVolumeMeshDS = occHandle(Ng_MeshVS_DataSource3D)(), int N=1);
    void setVolumeMesh(const occHandle(Ng_MeshVS_DataSource3D) &aVolumeMesh);
    void setStartingSurfaceElements(const occHandle(Ng_MeshVS_DataSourceFace) &facesDS);
    void setSteps(int NbSteps);
    bool getElements(occHandle(Ng_MeshVS_DataSource3D) &theElements);

private:

    occHandle(Ng_MeshVS_DataSource3D) myVolumeMeshDS;
    occHandle(Ng_MeshVS_DataSourceFace) myStartFace;
    int myN;
    int myNbSteps;
};

#endif // MESHSELECTNLAYERS_H
