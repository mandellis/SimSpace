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

//! ----
//! C++
//! ----
#include <iostream>
#include <fstream>

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

    void buildMap()
    {
        m.insert("DISP",TypeOfResult_U);
        m.insert("STRESS",TypeOfResult_S);
        m.insert("TOSTRAIN",TypeOfResult_TOSTRAIN);
        m.insert("MESTRAIN",TypeOfResult_MESTRAIN);
        m.insert("NDTEMP",TypeOfResult_NT);
        m.insert("UDR",TypeOfResult_UDR);
        m.insert("FORC",TypeOfResult_F);
        m.insert("PE",TypeOfResult_EPS);
        m.insert("FLUX",TypeOfResult_HFL);
        m.insert("CONTACT",TypeOfResult_CONT);
    }

signals:

public slots:

    //! set fatigue model
    void setFatigueModel (int fatigueAlgo);

    void setDataBase(meshDataBase *mDB);

    void setResultsFile(QString resultsFilePath);

    //! it calls internally an FrdReader which reads
    //!
    bool perform();

    QMap<GeometryTag,QList<QMap<int,double>>> evaluateResult(const QString &resultKeyName,
                                                             int requiredSubStepNb,
                                                             int requiredStepNb,
                                                             int requiredMode,
                                                             const QVector<GeometryTag> &vecLoc,
                                                             double &requiredTime);

    postObject evaluateFatigueResults(int type, QVector<GeometryTag> locs, const QList<double> &times, QMap<int,int> materialBodyMap, int nCycle);
    //postObject evaluateFatigueResults1(int type, QVector<GeometryTag> locs, const QList<double> &times, QMap<int,int> materialBodyMap, int nCycle);

    postObject buildPostObject(const QString &keyName,
                               int component,
                               int requiredSubStepNb,
                               int requiredStepNb,
                               int requiredMode,
                               const QVector<GeometryTag> &vecLoc);

private:

    //! create a string title for the colorbox, carrying the type of data
    QString resultName(const QString &keyName, int component, int step, int subStep, double time);

    //! time stamp
    QString timeStamp();

    //! read fatigue results
    QMap<GeometryTag,QMap<int,QList<double>>> readFatigueResults(int type,
                                                                    const QVector<GeometryTag> &vecLoc,
                                                                    const QList<double> &times);

    fatigueModel myFatigueModel;
    QMap<double,QVector<int>> myDTM;

public:

    //! retrieve the discrete time map
    //QMap<double,QVector<int>> getDTM();

    void setDiscreteTimeMap(const QMap<double,QVector<int>> &dtm);

    //! experimental
    void updateResultScale(postObject &aPostObject, int scaleType, double minValue, double maxValue, int NbIntervals);
};

#endif // POSTENGINE_H
