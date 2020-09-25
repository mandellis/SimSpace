#ifndef RAINFLOW_H
#define RAINFLOW_H

#include <geometrytag.h>

#include <cstdlib>
#include <iostream>
#include <fstream>

#include <vector>

//#include <math.h>

#include <stdio.h>

#include <string.h>
#include <stdlib.h>
#include <conio.h>

#include <QVector>
#include <QVariant>
#include <QList>

#include <time.h>
#include <sys/types.h>
#include <sys/timeb.h>

#define MAX 51000000

#include <QObject>
#include <QMetaType>
#include "mydefines.h"
#include "postengine.h"

using namespace std;
/*
enum fatigueModelType
{
    fatigueModel_BCM // Basquin Coffin Manson
    fatigueModel_ESR // Effective Strain Range (ASME VIII Div2 Assestment)
};
Q_DECLARE_METATYPE(fatigueModelType)

struct fatigueModel
{
    fatigueModelType type;
    QList<double> coeffs;
};
Q_DECLARE_METATYPE(fatigueModel)
*/

class rainflow: public QObject
{
    Q_OBJECT

public:

    //! constructor
    rainflow(QObject *parent = 0);

    //! constructor I
    rainflow(GeometryTag loc, QObject *parent = 0);

    //! constructor II
    rainflow(GeometryTag loc, fatigueModel fm, QObject *parent = 0);

    //! set location
    inline void setLocation(GeometryTag aLoc) { myLoc = aLoc; }

    //! set fatigue model
    inline void setFatigueModel (fatigueModel fm) { myFatigueModel = fm; }

    //! perform
    //bool perform(QMap<int, QList<double> > strainDistTimeHistory, QMap<int, double> &damageDist);
    bool perform(std::map<int, std::vector<double> > strainDistTimeHistory, std::map<int, double> &damageDist);

    //void read_data();

signals:


public slots:


private:

    //!tolerance on solveEquation
    double maxErr = 1.0e+0;
    //! tolerance on deltaEps value
    double tol2=3.5e-3;

    //! loc
    GeometryTag myLoc;

    //! fatigue model
    fatigueModel myFatigueModel;

    //! damage index
    double damage_index(std::vector<double> y);

    //! called within damage_index()
    int solve(double eps);

    //! ---
    std::vector<double> rainflow_engine(std::vector<double> y);
};

#endif // RAINFLOW_H
