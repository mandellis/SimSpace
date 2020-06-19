#include "vectortool.h"

vectorTool::vectorTool()
{
    ;
}

//! --------------------------------------------------------------------------//
//! function: getComponents(...)                                              //
//! details:  for a given vector defined by magnitude and direction corines   //
//!           returns the components with respect to the global coordinates   //
//! --------------------------------------------------------------------------//
QVector<double> vectorTool::getComponents(const double &magnitude, const QVector<double> &cosines)
{
    double Xglobal_component = magnitude*cosines.at(0);
    double Yglobal_component = magnitude*cosines.at(1);
    double Zglobal_component = magnitude*cosines.at(2);
    QVector<double> vec;
    vec.push_back(Xglobal_component);
    vec.push_back(Yglobal_component);
    vec.push_back(Zglobal_component);
    return vec;
}

//! --------------------------------------------------------------------------//
//! function: getComponents()                                                 //
//! details:  for a given vector V=QVector<double> components, in a given     //
//!           system of reference directionalData=QVector<QVector<double>>    //
//!           returns the components of V in the global coordinate system     //
//! --------------------------------------------------------------------------//
QVector<double> vectorTool::getComponents(const QVector<double> components, const QVector<QVector<double> > &directionalData)
{
    double a11 = directionalData.at(0).at(0);
    double a12 = directionalData.at(0).at(1);
    double a13 = directionalData.at(0).at(2);
    double a21 = directionalData.at(1).at(0);
    double a22 = directionalData.at(1).at(1);
    double a23 = directionalData.at(1).at(2);
    double a31 = directionalData.at(2).at(0);
    double a32 = directionalData.at(2).at(1);
    double a33 = directionalData.at(2).at(2);

    double Xglobal_component, Yglobal_component, Zglobal_component;
    Xglobal_component = components.at(0)*a11+components.at(1)*a21+components.at(2)*a31;
    Yglobal_component = components.at(0)*a12+components.at(1)*a22+components.at(2)*a32;
    Zglobal_component = components.at(0)*a13+components.at(1)*a23+components.at(2)*a33;

    QVector<double> vec;
    vec.push_back(Xglobal_component);
    vec.push_back(Yglobal_component);
    vec.push_back(Zglobal_component);
    return vec;
}
