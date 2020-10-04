#ifndef TOOLS_H
#define TOOLS_H

//! C++
#include <iostream>
#include <fstream>

//! Qt
#include <QString>
#include <QVariant>
#include <QColor>
#include <QMetaType>
#include <QWidget>
#include <QTime>
#include <QApplication>
#include <QEventLoop>

//! custom includes
#include "simulationnodeclass.h"
#include "geometrytag.h"

class tools
{
public:

    tools();

    //! new
    template<class T>
    static void writeQVector2(const QVector<T> &vec, std::ofstream &os)
    {
        int Nelements = vec.size();
        os<<Nelements<<std::endl;
        for(int i=0; i<Nelements; i++)
        {
            os<<vec.at(i)<<std::endl;
        }
    }

    template<class T>
    static void writeQVector(const T &vec, std::ofstream &os)
    {
        int Nelements = vec.size();
        os<<Nelements<<std::endl;
        for(int i=0; i<Nelements; i++)
        {
            os<<vec.at(i)<<std::endl;
        }
    }

    template<class T>
    static void writeTensor2(const T &tensor, std::ofstream &os)
    {
        int Nrow = tensor.size();
        int Ncol = tensor.at(0).size();
        os<<Nrow<<endl;
        os<<Ncol<<endl;

        for(int i=0;i<Nrow;i++)
            for(int k=0; k<Ncol; k++)
                os<<tensor.at(i).at(k)<<std::endl;
    }

    template<class T>
    static void writeTensor(const QVector<QVector<T>> &tensor, std::ofstream &os)
    {
        int Nrow = tensor.size();
        int Ncol = tensor.at(0).size();
        os<<Nrow<<endl;
        os<<Ncol<<endl;

        for(int i=0;i<Nrow;i++)
            for(int k=0; k<Ncol; k++)
                os<<tensor.at(i).at(k)<<std::endl;
    }

    template<class T>
    static QVector<QVector<T>> readTensor2(std::ifstream &is)
    {
        int Nrow, Ncol;
        is>>Nrow;
        is>>Ncol;

        QVector<QVector<T>> tensor;
        for(int i=0;i<Nrow;i++)
        {
            QVector<T> vec;
            for(int j=0;j<Ncol;j++)
            {
                T val;
                is>>val;
                vec.push_back(val);
            }
            tensor.push_back(vec);
        }
        return tensor;
    }

    template<class T>
    static QVector<T> readQVector(std::ifstream &is)
    {
        int Nelements;
        is>>Nelements;
        T val;
        QVector<T> vec;
        for(int i=0; i<Nelements; i++)
        {
            is>>val;
            vec.push_back(val);
        }
        return vec;
    }

    static void clearDir(const QString &path);

    static void writeQVariant(const QVariant &var, std::ofstream &os);
    static QVariant readQVariant(std::ifstream &is);

    static std::vector<GeometryTag> getScope(SimulationNodeClass *aNode);
    static std::vector<GeometryTag> getScope(QExtendedStandardItem *item);

    static void writeVectorOfLocations(const std::vector<GeometryTag> &vecLocs, std::ofstream &os);
    static std::vector<GeometryTag> readVectorOfLocations(std::ifstream &is);

    static void writeColor(QColor color, std::ofstream &os);
    static QColor readColor(std::ifstream &is);

    static QWidget* getWidgetByName(const QString &widgetName);

    //! get the working directoty
    static QString getWorkingDir();

    //! eliminate duplicated values from a vector
    static std::vector<int> clearFromDuplicates(const std::vector<int> &aVec);

    //! change the opacity of an icon
    static void changeIconOpacity(QExtendedStandardItem *item, bool isOpaque=false);

    //! time stamp
    static QString timeStamp();

    //! get path of executable
    static std::string getPathOfExecutable();
};

#endif // TOOLS_H
