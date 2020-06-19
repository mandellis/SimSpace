#ifndef VECTORTOOL_H
#define VECTORTOOL_H

#include <QVector>

class vectorTool
{
public:

    vectorTool();
    static QVector<double> getComponents(const double &magnitude, const QVector<double> &cosines);
    static QVector<double> getComponents(const QVector<double> components, const QVector<QVector<double>> &directionalData);
};

#endif // VECTORTOOL_H
