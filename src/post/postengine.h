#ifndef POSTENGINE_H
#define POSTENGINE_H

//! ---
//! Qt
//! ---
#include <QObject>

//! ----------------
//! custom includes
//! ----------------
#include "frdreader.h"
#include <meshdatabase.h>
#include "postobject.h"
#include "global.h"
#include "qprogressindicator.h"
#include "qprogressevent.h"

//! ----
//! C++
//! ----
#include <iostream>
#include <fstream>
#include <map>
#include <memory>

class OCCMeshToCCXmesh;

enum fatigueModelType
{
    fatigueModel_BCM, // Basquin Coffin Manson
    fatigueModel_ESR  // Effective Strain Range (ASME VIII Div2 Assestment)
};
Q_DECLARE_METATYPE(fatigueModelType)

struct fatigueModel
{
    fatigueModelType type;
    QList<double> coeffs;
};
Q_DECLARE_METATYPE(fatigueModel)

class postEngine : public QObject
{
    Q_OBJECT

public:

    enum TypeOfResult
    {
        TypeOfResult_U,
        TypeOfResult_S,
        TypeOfResult_TOSTRAIN,
        TypeOfResult_MESTRAIN,
        TypeOfResult_NT,
        TypeOfResult_UDR,
        TypeOfResult_F,
        TypeOfResult_EPS,
        TypeOfResult_CONT,
        TypeOfResult_HFL
    };

public:

    explicit postEngine(QObject *parent = 0);

private:

    meshDataBase *myMeshDataBase;
    QString myResultsFilePath;
    QMap<QString,TypeOfResult> m;

    void buildMap();

signals:

public slots:

    //! set fatigue model
    void setFatigueModel (int fatigueAlgo);

    //! set mesh data base
    void setDataBase(meshDataBase *mDB);

    //! set results file
    void setResultsFile(QString resultsFilePath);

    //! it calls internally an FrdReader object
    bool perform();

    //! evaluate results
    std::map<GeometryTag,std::vector<std::map<int,double>>> evaluateResult(const QString &resultKeyName,
                                                                            int requiredSubStepNb,
                                                                            int requiredStepNb,
                                                                            int requiredMode,
                                                                            const std::vector<GeometryTag> &vecLoc,
                                                                            double &requiredTime);

    //! evaluate fatigue results
    bool buildFatiguePostObject(int type, const std::vector<GeometryTag> &locs, std::vector<double> times, QMap<int,int> materialBodyMap, int nCycle, sharedPostObject &aPostObject);

    //! build a post object
    bool buildPostObject(const QString &keyName,
                         int component,
                         int requiredSubStepNb,
                         int requiredStepNb,
                         int requiredMode,
                         const std::vector<GeometryTag> &vecLoc,
                         sharedPostObject &aPostObject);

private:

    //! create a string title for the colorbox, carrying the type of data
    QString resultName(const QString &keyName, int component, int step, int subStep, double time);

    //! time stamp
    QString timeStamp();

    //! read fatigue results
    std::map<GeometryTag, std::map<int, std::vector<double>>> readFatigueResults(int type,
                                                                                 const std::vector<GeometryTag> &vecLoc,
                                                                                 std::vector<double> times);

    //! fatigue model
    fatigueModel myFatigueModel;

    //! discrete time map
    QMap<double,QVector<int>> myDTM;

public:

    void setDiscreteTimeMap(const QMap<double,QVector<int>> &dtm);
    void updateResultsPresentation(QList<sharedPostObject> &postObjectList);
    void updateIsostrips(sharedPostObject &aPostObject, int scaleType, double minValue, double maxValue, int NbIntervals);

private:

    void groupDeformationFieldByBodies(const std::map<GeometryTag,std::map<int,gp_Vec>> &mapDisplMap,
                                       std::map<GeometryTag,std::map<int,gp_Vec>> &mapDisplMap_byBodies);

    void groupTagsByBodies(const std::vector<GeometryTag> &vecLoc, std::vector<GeometryTag> &vecLoc_byBodies);

    void groupResultsByBodies(const std::map<GeometryTag,std::vector<std::map<int,double>>> &resMap,
                              std::map<GeometryTag,std::vector<std::map<int,double>>> &resMap_byBody);

    void groupAndMergeMeshDataSourcesByBodies(const std::vector<GeometryTag> &vecLoc,
                                              std::map<GeometryTag, opencascade::handle<MeshVS_DataSource> > &meshDSforResults);

};

#endif // POSTENGINE_H
