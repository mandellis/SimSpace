//! ----------------
//! custom includes
//! ----------------
#include "load.h"
#include "tools.h"

//! -----------------------------------------
//! function: constructor
//! details:  creates a load vector with one
//!           value ("0") and load type none
//! -----------------------------------------
load::load()
{
    QVector<QVariant> valueVec;
    QVariant data;
    data.setValue(0);
    valueVec.push_back(data);

    myLoadType = Property::loadType_none;
    for(int i=0; i<valueVec.size(); i++) myValues.push_back(valueVec.at(i));
}

//! ------------------------
//! function: constructor I
//! details:
//! ------------------------
load::load(QVector<QVariant> values,Property::loadType type)
{
    //cout<<"load::load()->____constructor called____"<<endl;
    myLoadType = type;
    for(int k=0; k<values.size();k++) myValues.push_back(values.at(k));
}

//! ------------------
//! function: setData
//! details:
//! ------------------
void load::setData(const QVector<QVariant> &values)
{
    myValues.clear();
    for(int k=0; k<values.size();k++) myValues.push_back(values.at(k));
}

//! ------------------
//! function: NbTimes
//! details:
//! ------------------
int load::NbTimes() const
{
    return myValues.size();
}

//! -----------------
//! function: values
//! details:
//! -----------------
QVector<QVariant> load::values() const
{
    return myValues;
}

//! ---------------
//! function: type
//! details:
//! ---------------
Property::loadType load::type() const
{
    return myLoadType;
}

//! ------------------
//! function: setType
//! details:
//! ------------------
void load::setType(Property::loadType theType)
{
    myLoadType = theType;
}

//! ----------------
//! function: write
//! details:
//! ----------------
void load::write(std::ofstream &out) const
{
    //! -----------------------
    //! write the type of load
    //! -----------------------
    QString loadType = this->getLoadType();
    QVariant var;
    var.setValue(loadType);
    tools::writeQVariant(var,out);

    switch(myLoadType)
    {
    //! --------------------------------
    //! write a vector of vector of int
    //! --------------------------------
    case Property::loadType_autoTimeStepping:
    case Property::loadType_timeIncrementationParameters:
    case Property::loadType_outputSettings:
    case Property::loadType_storeResultsAt:
    {
        //! write the type of content
        QVariant data;
        data.setValue(QString("vec_int"));
        tools::writeQVariant(data,out);

        //! write the number of components
        tools::writeQVariant(QVariant(myValues.size()),out);

        //! write the components
        QVector<QVector<int>> vec2Int;
        for(int i=0;i<myValues.size();i++)
        {
            vec2Int.push_back(myValues.at(i).value<QVector<int>>());
            tools::writeQVector2<int> (myValues.at(i).value<QVector<int>>(),out);
        }
    }
        break;
    //! -------------------------
    //! write a vector of double
    //! -------------------------
    case Property::loadType_fieldParameters:
    case Property::loadType_cutBackParameters:
    case Property::loadType_lineSearchParameters:
    {
        //! write the type of content
        QVariant data;
        data.setValue(QString("vec_double"));
        tools::writeQVariant(data,out);

        //! write the number of components
        tools::writeQVariant(QVariant(myValues.size()),out);

        //! write the components
        QVector<QVector<double>> vec2Double;
        for(int i=0;i<myValues.size();i++)
        {
            vec2Double.push_back(myValues.at(i).value<QVector<double>>());
            tools::writeQVector2<double> (myValues.at(i).value<QVector<double>>(),out);
        }
    }
        break;
    //! ----------------------
    //! write a vector of int
    //! ----------------------
    case Property::loadType_stepNumber:
    case Property::loadType_solverType:
    case Property::loadType_boltStatusDefinedBy:
    case Property::loadType_fluxConvergence:
    case Property::loadType_solutionConvergence:
    case Property::loadType_timeIncrementation:
    case Property::loadType_cutBack:
    case Property::loadType_lineSearch:
    case Property::loadType_modelChange:
    case Property::loadType_analysisType:
    case Property::loadType_timeIntegration:
    {
        //! write the type of content
        QVariant data;
        data.setValue(QString("int"));
        tools::writeQVariant(data,out);

        //! build the vector and write it
        QVector<int> vecInt;
        for(int i=0;i<myValues.size();i++)
        {
            if(myLoadType == Property::loadType_boltStatusDefinedBy)
            {
                cout<<"____bolt status: "<<myValues.at(i).toInt()<<endl;
            }
            vecInt.push_back(myValues.at(i).toInt());
        }
        tools::writeQVector2<int>(vecInt,out);
    }
        break;
    //! -------------------------
    //! write a vector of double
    //! -------------------------
    case Property::loadType_temperatureMagnitude:
    case Property::loadType_thermalFluxMagnitude:
    case Property::loadType_thermalFlowMagnitude:
    case Property::loadType_thermalPowerMagnitude:
    case Property::loadType_thermalConvectionFilmCoefficientMagnitude:
    case Property::loadType_thermalConvectionReferenceTemperatureMagnitude:
    case Property::loadType_displacementMagnitude:
    case Property::loadType_displacementX:
    case Property::loadType_displacementY:
    case Property::loadType_displacementZ:
    case Property::loadType_remoteRotationX:
    case Property::loadType_remoteRotationY:
    case Property::loadType_remoteRotationZ:
    case Property::loadType_remoteRotationMagnitude:
    case Property::loadType_forceMagnitude:
    case Property::loadType_forceX:
    case Property::loadType_forceY:
    case Property::loadType_forceZ:
    case Property::loadType_momentMagnitude:
    case Property::loadType_momentX:
    case Property::loadType_momentY:
    case Property::loadType_momentZ:
    case Property::loadType_rotationalVelocityX:
    case Property::loadType_rotationalVelocityY:
    case Property::loadType_rotationalVelocityZ:
    case Property::loadType_rotationalVelocityMagnitude:
    case Property::loadType_accelerationMagnitude:
    case Property::loadType_accelerationX:
    case Property::loadType_accelerationY:
    case Property::loadType_accelerationZ:
    case Property::loadType_pressureMagnitude:
    case Property::loadType_thermalConditionTemperature:
    case Property::loadType_stepEndTime:
    case Property::loadType_time:
    case Property::loadType_boltForce:
    case Property::loadType_boltAdjustment:
    {
        //! write the type of content
        QVariant data;
        data.setValue(QString("double"));
        tools::writeQVariant(data,out);

        //! build the vector and write it
        QVector<double> vecDouble;
        for(int i=0;i<myValues.size();i++) vecDouble.push_back(myValues.at(i).toDouble());
        tools::writeQVector2<double>(vecDouble,out);
    }
        break;
    }
}

//! ------------------------------------------
//! function: readLoad
//! details:  return a load - static function
//! ------------------------------------------
load load::readLoad(std::ifstream &in)
{
    //! ----------------------
    //! read the type of load
    //! ----------------------
    QVariant var = tools::readQVariant(in);
    QString loadType = var.toString();
    cout<<"load::readLoad()->____reading load of type: "<<loadType.toStdString()<<"____"<<endl;

    QVector<QVariant> vecVariant;
    QVariant data;

    //! -------------------------
    //! read the type of content
    //! -------------------------
    var = tools::readQVariant(in);
    QString type = var.toString();
    cout<<"load::readLoad->reading type of value: "<<type.toStdString()<<"<----"<<endl;

    if(loadType == "loadType_stepEndTime" ||
            loadType == "loadType_forceX" ||
            loadType == "loadType_forceY" ||
            loadType == "loadType_forceZ" ||
            loadType == "loadType_forceMagnitude" ||
            loadType == "loadType_rotationalVelocityMagnitude" ||
            loadType == "loadType_rotationalVelocityX" ||
            loadType == "loadType_rotationalVelocityY" ||
            loadType == "loadType_rotationalVelocityZ" ||
            loadType == "loadType_accelerationMagnitude" ||
            loadType == "loadType_accelerationX" ||
            loadType == "loadType_accelerationY" ||
            loadType == "loadType_accelerationZ" ||
            loadType == "loadType_momentX" ||
            loadType == "loadType_momentY" ||
            loadType == "loadType_momentZ" ||
            loadType == "loadType_momentMagnitude" ||
            loadType == "loadType_displacementX" ||
            loadType == "loadType_displacementY" ||
            loadType == "loadType_displacementZ" ||
            loadType == "loadType_displacementMagnitude" ||
            loadType == "loadType_pressureMagnitude" ||
            loadType == "loadType_thermalConditionTemperature" ||
            loadType == "loadType_magnitude" ||
            loadType == "loadType_temperatureMagnitude" ||
            loadType == "loadType_stepEndTime" ||
            loadType == "loadType_time" ||
            loadType == "loadType_boltForce" ||
            loadType == "loadType_boltAdjustment" ||
            loadType == "loadType_thermalFlowMagnitude" ||
            loadType == "loadType_thermalFluxMagnitude" ||
            loadType == "loadType_thermalConvectionFilmCoefficientMagnitude" ||
            loadType == "loadType_thermalConvectionReferenceTemperatureMagnitude")
    {
        QVector<double> vecDouble = tools::readQVector<double>(in);
        for(int i=0; i<vecDouble.size();i++)
        {
            cout<<"____"<<type.toStdString()<<"_____value: "<<vecDouble[i]<<"____"<<endl;
            data.setValue(vecDouble.at(i));
            vecVariant.push_back(data);
        }
    }
    else if (loadType == "loadType_solverType" ||
             loadType==  "loadType_stepNumber" ||
             loadType == "loadType_boltStatusDefinedBy" ||
             loadType == "loadType_solutionConvergence" ||
             loadType == "loadType_fluxConvergence" ||
             loadType == "loadType_timeIncrementation" ||
             loadType == "loadType_cutBack" ||
             loadType == "loadType_lineSearch" ||
             loadType == "loadType_modelChange" ||
             loadType == "loadType_analysisType" ||
             loadType == "loadType_timeIntegration")
    {
        QVector<int> vecInt = tools::readQVector<int>(in);
        for(int i=0; i<vecInt.size();i++)
        {
            int val = vecInt.at(i);
            if(loadType == "loadType_boltStatusDefinedBy")
            {
                cout<<"_____reading bolt type: "<<val<<"____"<<endl;
                Property::boltStatusDefinedBy boltStatus = static_cast<Property::boltStatusDefinedBy>(val);
                data.setValue(boltStatus);
                vecVariant.push_back(data);
                continue;
            }
            Property::loadType lt = static_cast<Property::loadType>(val);
            data.setValue(lt);
            vecVariant.push_back(data);

        }
    }
    else if(loadType == "loadType_autoTimeStepping" ||
            loadType == "loadType_timeIncrementationParameters" ||
            loadType == "loadType_outputSettings" ||
            loadType == "loadType_storeResultsAt")
    {
        int N = tools::readQVariant(in).toInt();
        QVector<int> vecInt;
        for(int i=0; i<N; i++)
        {
            vecInt = tools::readQVector<int>(in);
            data.setValue(vecInt);
            vecVariant.push_back(data);
        }
    }
    else if(loadType == "loadType_fieldParameters" ||
            loadType == "loadType_cutBackParameters" ||
            loadType == "loadType_lineSearchParameters")
    {
        int N = tools::readQVariant(in).toInt();
        QVector<double> vecDouble;
        for(int i=0; i<N; i++)
        {
            vecDouble = tools::readQVector<double>(in);
            data.setValue(vecDouble);
            vecVariant.push_back(data);
        }
    }

    load aLoad(vecVariant);
    aLoad.setLoadType(loadType);
    return aLoad;
}

//! ----------------------
//! function: getLoadType
//! details:
//! ----------------------
const QString load::getLoadType() const
{
    const QMetaObject &mo = Property::staticMetaObject;
    int index = mo.indexOfEnumerator("loadType");
    //cout<<"load::getLoadType()->____index: "<<index<<endl;
    QMetaEnum metaEnum = mo.enumerator(index);
    return metaEnum.valueToKey(myLoadType);
}

//! ----------------------
//! function: setLoadType
//! details:
//! ----------------------
void load::setLoadType(const QString &theLoadType)
{
    const QMetaObject &mo = Property::staticMetaObject;
    int index = mo.indexOfEnumerator("loadType");
    //cout<<"load::setLoadType->____index: "<<index<<endl;
    QMetaEnum metaEnum = mo.enumerator(index);
    int value = metaEnum.keyToValue(qPrintable(theLoadType));
    myLoadType = static_cast<Property::loadType>(value);
}
