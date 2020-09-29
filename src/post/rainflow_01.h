#ifndef RAINFLOW_01_H
#define RAINFLOW_01_H

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

class rainflow_01: public QObject
{
    Q_OBJECT

public:

    //! constructor
    rainflow_01(QObject *parent = 0);

    //! constructor I
    rainflow_01(GeometryTag loc, QObject *parent = 0);

    //! constructor II
    rainflow_01(GeometryTag loc, fatigueModel fm, QObject *parent = 0);

    //! set location
    inline void setLocation(GeometryTag aLoc) { myLoc = aLoc; }

    //! set fatigue model
    inline void setFatigueModel (fatigueModel fm) { myFatigueModel = fm; }

    //! perform
    bool perform(std::map<int, std::vector<double> > strainDistTimeHistory, std::map<int, double> &damageDist);

    //void read_data();





signals:


public slots:


private:

    int rf3(double *array_ext, int nr, double *array_out);

    int rf5(double *array_ext, int nr, double *array_t, double *array_out);

    int sig2ext(double *sig, double *time_sig, long n, int clsn,
                double *ext, double *exttime);

    double arr_min(double *sig, int n, int *pos);
    double arr_max(double *sig, int n, int *pos);
    #define NNEW(a,b) (a *)calloc((b),sizeof(a))
    #define RENEW(a,b,c) a=(b *) realloc((b *)(a),(c)*sizeof(b))
    double *diff(double *vec, int n);
    int repl(double *x, int *filt, int n, double *x_repl);

    double damage_index(double *y,long timeSize);


    //!tolerance on solveEquation
    const double maxErr = 1.0e+0;

    //! tolerance on deltaEps value
    const double tol2 = 3.5e-3;

    //! loc
    GeometryTag myLoc;

    //! fatigue model
    fatigueModel myFatigueModel;

    //! damage index
    double damage_index(const std::vector<double> &y);

    //! called within damage_index()
    double solve(double eps, double epsF, double c, double sigmaF, double E, double b);
    double solve_exact(double eps, double epsF, double c, double sigmaF, double E, double b);

    //! ---
    std::vector<double> rainflow_engine_01(std::vector<double> y);
};

#endif // RAINFLOW_H
