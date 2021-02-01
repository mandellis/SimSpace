#ifndef CONTACTPARAMETERS_H
#define CONTACTPARAMETERS_H

#include <ng_meshvs_datasourceface.h>
#include <geometrydatabase.h>

class contactParameters
{
public:

    contactParameters();

    static double polygon_area (const QList<QList<double>> &points);
    static double calc_discretizedArea(const occHandle(Ng_MeshVS_DataSourceFace) &faceMesh);
    static double calc_K(QList<occHandle(Ng_MeshVS_DataSourceFace)> master,QList<occHandle(Ng_MeshVS_DataSourceFace)> slave);
    static double calc_aveDiscretizedArea(const occHandle(Ng_MeshVS_DataSourceFace) &faceMesh);

};

#endif // CONTACTPARAMETERS_H
