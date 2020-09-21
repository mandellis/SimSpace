#ifndef MESHSCRIPT_H
#define MESHSCRIPT_H

#include <QObject>
#include <QString>

#include "mainwindow.h"
#include "qextendedstandarditem.h"

class SimulationNodeClass;

class MeshScript : public QObject
{

    Q_OBJECT

public:

    explicit MeshScript(QObject *parent);

public slots:

    //Item & node creation

    void generateMesh();
    void generateSurfaceMesh();
    void insertMethod(QString name);
    void insertMeshType(QString name);
    void insertBodySizing(QString name);
    void insertFaceSizing(QString name);
    void insertEdgeSizing(QString name);
    void insertVertexSizing(QString name);
    void insertPrismaticLayer(QString name);
    void insertMeshMetric(QString name);

    //Mesh

    void setShowMeshNodes(QString value);
    void setRelevance(int value);
    void setInitialSizeSeed(QString value);
    void setSmoothing(QString value);
    void setElementMidsideNodes(QString value);
    void setSubmeshes(QString value);
    void setStraightSidedElements(QString value);

    //General

    void setScopingMethod(QString nodeName, QString value);
    void setVolumeScope(QString nodeName, int pSNr);
    void setFaceScope(QString nodeName, int pSNr, int sTNr);
    void setEdgeScope(QString nodeName, int pSNr, int sTNr);
    void setNamedSelection(QString nodeName, QString value);
    void setSuppressed(QString nodeName, QString value);

    //Method

    void setPatchConforming(QString nodeName, QString value);
    void setTessellator(QString nodeName, QString value);
    void setAngularDeflection(QString nodeName, double value);
    void setLinearDeflection(QString nodeName, double value);
    void setMinFaceSize(QString nodeName, double value);
    void setMaxFaceSize(QString nodeName, double value);
    void setMeshOrder(QString nodeName, QString value);
    void setSurfaceMesher(QString nodeName, QString value);
    void setVolumeMesher(QString nodeName, QString value);
    void setRunInMemory(QString nodeName, QString value);

    //Mesh type

    void setMeshType(QString nodeName, QString value);

    //Body sizing

    void setMinElementSize(QString nodeName, double value);
    void setMaxElementSize(QString nodeName, double value);
    void setGrading(QString nodeName, double value);

    //Face sizing

    void setFaceSizing(QString nodeName, double value);

    //Edge sizing

    void setSizingType(QString nodeName, QString value);
    void setEdgeElementSize(QString nodeName, double value);
    void setNumberOfDivisions(QString nodeName, int value);

    //Vertex sizing

    void setPinball(QString nodeName, double value);
    void setVertexElementSize(QString nodeName, double value);

    //Prismatic layer

    void setBoundaryScopingMethod(QString nodeName, QString value);
    void setBoundary(QString nodeName, int pSNr, int sTNr);
    void setBoundaryNamedSelection(QString nodeName, QString value);
    void setOptions(QString nodeName, QString value);
    void setNumberOfLayers(QString nodeName, int value);
    void setFirstLayerHeight(QString nodeName, double value);
    void setTotalThickness(QString nodeName, double value);
    void setExpansionRatio(QString nodeName, double value);
    void setAlgorithm(QString nodeName, QString value);
    void setBoundaryMeshType(QString nodeName, QString value);
    void setGuidingVectorsSmoothingSteps(QString nodeName, int value);
    void setThicknessSmoothingSteps(QString nodeName, int value);
    void setGVSCurvatureSensivity(QString nodeName, double value);
    void setTSCurvatureSensivity(QString nodeName, double value);
    void setCurvatureSensivity(QString nodeName, double value);
    void setLockBoundary(QString nodeName, QString value);
    void setCheckSelfIntersections(QString nodeName, QString value);
    void setCheckMutualIntersections(QString nodeName, QString value);

    //Mesh metric

    void setMetricType(QString nodeName, QString value);

private:

    SimulationNodeClass* getNode(QString nodeName);

    QExtendedStandardItem* getNSItem(QString nodeName);

    QExtendedStandardItem* getMeshItem(QString nodeName);

    MainWindow *mw;
};

#endif // MESHSCRIPT_H
