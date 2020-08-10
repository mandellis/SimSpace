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

//! ----
//! C++
//! ----
#include <iostream>
#include <fstream>
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

    //! ---------------------------------------------------
    //! from CCX node to OCC mesh node
    //! location (m,n) <=> QMap<CCXnodeID, OCCnodeID>
    //! ---------------------------------------------------
    OCCMeshToCCXmesh *OCCtoCCXinterface;


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
    QMap<GeometryTag,QList<QMap<int,double>>> evaluateResult(const QString &resultKeyName,
                                                             int requiredSubStepNb,
                                                             int requiredStepNb,
                                                             int requiredMode,
                                                             const QVector<GeometryTag> &vecLoc,
                                                             double &requiredTime);

    //! evaluate fatigue results
    bool evaluateFatigueResults(int type, QVector<GeometryTag> locs, const QList<double> &times, QMap<int,int> materialBodyMap, int nCycle, sharedPostObject &aPostObject);

    //! build a post object
    bool buildPostObject(const QString &keyName,
                         int component,
                         int requiredSubStepNb,
                         int requiredStepNb,
                         int requiredMode,
                         const QVector<GeometryTag> &vecLoc,
                         sharedPostObject &aPostObject);

private:

    //! create a string title for the colorbox, carrying the type of data
    QString resultName(const QString &keyName, int component, int step, int subStep, double time);

    //! time stamp
    QString timeStamp();

    //! read fatigue results
    QMap<GeometryTag,QMap<int,QList<double>>> readFatigueResults(int type,
                                                                 const QVector<GeometryTag> &vecLoc,
                                                                 const QList<double> &times);

    //! fatigue model
    fatigueModel myFatigueModel;

    //! discrete time map
    QMap<double,QVector<int>> myDTM;

public:

    void setDiscreteTimeMap(const QMap<double,QVector<int>> &dtm);
    void updateResultsPresentation(QList<sharedPostObject> &postObjectList);
    void updateIsostrips(sharedPostObject &aPostObject, int scaleType, double minValue, double maxValue, int NbIntervals);
};

#endif // POSTENGINE_H
