#include "customtablemodel.h"
#include "load.h"

#include <QVector>
#include <QHeaderView>
#include <QRect>
#include <QColor>

#include <iostream>
using namespace std;

//! ----------------------------------
//! function: constructor
//! details:  for "Analysis settings"
//! ----------------------------------
CustomTableModel::CustomTableModel(const QVector<load> &vecLoads, bool addFirstRow, QObject *parent):QAbstractTableModel(parent)
{
    //! the list of load type
    for(int k=0; k<vecLoads.size(); k++)
        m_loadTypes.append(vecLoads.at(k).type());

    m_columnCount = vecLoads.size();
    m_rowCount = vecLoads.at(0).NbTimes();  //! number of values

    for(int i=0; i<m_rowCount; i++)
    {
        QVector<QVariant> *dataVec = new QVector<QVariant>;
        for(int col=0; col<m_columnCount; col++)
        {
            //!cout<<"type: "<<vecLoads.at(col).values().at(i).typeName()<<endl;
            dataVec->push_back(vecLoads.at(col).values().at(i));
        }
        m_data.append(dataVec);
    }

    //! create the additional first row (0-th row)
    //! put inside the if statement... to do...
    QVector<QVariant> *dataVec_first = new QVector<QVariant>;

    if(addFirstRow==true)
    {
        //! step number
        dataVec_first->push_back(1);

        //! time
        dataVec_first->push_back(0.0);

        for(int i=2; i<m_data.at(0)->size();i++)
            dataVec_first->push_back(m_data.at(0)->at(i));
        m_data.insert(0,dataVec_first);
    }
}

//! -------------------------------------
//! function: autoUpdateTable
//! details:  used for "Bolt pretension"
//! -------------------------------------
void CustomTableModel::autoUpdateTable(QModelIndex topLeftIndex, QModelIndex /*bottomRightIndex*/)
{
    //cout<<"CustomTableModel::autoUpdateTable()->____function called____"<<endl;
    cout<<"CustomTableModel::autoUpdateTable()->____function called. Row: "<<topLeftIndex.row()<<" column: "<<topLeftIndex.column()<<"____"<<endl;

    if(strcmp(topLeftIndex.data(Qt::EditRole).typeName(),"Property::boltStatusDefinedBy")==0)
    //if(strcmp(topLeftIndex.data(Qt::EditRole).typeName(),"Property::defineBy")==0)
    {

        Property::boltStatusDefinedBy boltStatusDefineBy = topLeftIndex.data(Qt::EditRole).value<Property::boltStatusDefinedBy>();
        //Property::defineBy boltStatusDefineBy = topLeftIndex.data(Qt::EditRole).value<Property::defineBy>();
        QVariant value;
        int currentRow = topLeftIndex.row();
        int loadColumn = topLeftIndex.column()+1;
        int adjustmentColumn = topLeftIndex.column()+2;
        cerr<<"CustomTableModel::autoUpdateTable()->____{"<<topLeftIndex.column()<<", "<<loadColumn<<", "<<adjustmentColumn<<"}____"<<endl;
        QModelIndex indexAdjustment = this->createIndex(currentRow,adjustmentColumn);

        switch(boltStatusDefineBy)
        {
        //case Property::defineBy_load:
        case Property::boltStatusDefinedBy_load:
        {
            cerr<<"CustomTableModel::autoUpdateTable()->____disabling adjustment, enablig load____"<<endl;
            //! ------------------------
            //! disable adjustment cell
            //! ------------------------
            value.setValue(QString("N/A"));
            m_data[currentRow]->replace(adjustmentColumn, value);
            //! -----------------------
            //! enable load cell
            //! -----------------------
            value.setValue(0.0);
            m_data[currentRow]->replace(loadColumn, value);

            emit dataChanged(topLeftIndex,indexAdjustment);
            emit boltStatusDefineByChangedThroughTabularView();
        }
            break;

        case Property::boltStatusDefinedBy_adjustment:
        //case Property::defineBy_adjustment:
        {
            cerr<<"CustomTableModel::autoUpdateTable()->____disabling load, enabling adjustment____"<<endl;
            //! -----------------------
            //! disable load cell
            //! -----------------------
            value.setValue(QString("N/A"));
            m_data[currentRow]->replace(loadColumn, value);
            //! -----------------------
            //! enable adjustment cell
            //! -----------------------
            value.setValue(0.0);
            m_data[currentRow]->replace(adjustmentColumn, value);

            emit dataChanged(topLeftIndex,indexAdjustment);
            emit boltStatusDefineByChangedThroughTabularView();
        }
            break;

        case Property::boltStatusDefinedBy_open:
        case Property::boltStatusDefinedBy_lock:
        //case Property::defineBy_open:
        //case Property::defineBy_lock:
        {
            cerr<<"CustomTableModel::autoUpdateTable()->____disabling load and adjustment____"<<endl;
            //! ----------------------------------
            //! disable load and adjustment cells
            //! ----------------------------------
            value.setValue(QString("N/A"));
            m_data[currentRow]->replace(loadColumn, value);
            m_data[currentRow]->replace(adjustmentColumn, value);
            emit dataChanged(topLeftIndex,indexAdjustment);
            emit boltStatusDefineByChangedThroughTabularView();
        }
            break;
        }
    }
}

//! ---------------------
//! function: headerData
//! details:
//! ---------------------
QVariant CustomTableModel::headerData(int section, Qt::Orientation orientation, int role) const
{
    if (role!= Qt::DisplayRole) return QVariant();

    if (orientation == Qt::Horizontal)
    {
        return this->getHeaderString(section);
    }
    else
    {
        switch (section)
        {
        case 0: return 1; break;
        case 1: return 1; break;
        default: return section; break;
        }
    }
}

//! -----------------------------------------
//! function: dataRC
//! details:  get data through (row, column)
//! -----------------------------------------
QVariant CustomTableModel::dataRC(int row, int column, int role) const
{
    if (role == Qt::EditRole)
    {
        return m_data[row]->at(column);
    }
    return QVariant();
}

//! -----------------------------------------
//! function: data
//! details:  setup also the Qt::DisplayRole
//! -----------------------------------------
QVariant CustomTableModel::data(const QModelIndex &index, int role) const
{
    if (role == Qt::DisplayRole)
    {
        QVariant data = m_data[index.row()]->at(index.column());
        if(m_loadTypes.at(index.column())==Property::loadType_stepNumber)
        {
            QString cellString = QString("%1").arg(data.toInt());
            QVariant cellStringVariant;
            cellStringVariant.setValue(cellString);
            return cellStringVariant;
        }
        else if(m_loadTypes.at(index.column())==Property::loadType_stepEndTime ||
                m_loadTypes.at(index.column())==Property::loadType_time)
        {
            QString cellString = QString("%1").arg(data.toDouble());
            QVariant cellStringVariant;
            cellStringVariant.setValue(cellString);
            return cellStringVariant;
        }
        else if(m_loadTypes.at(index.column())==Property::loadType_solverType)
        {
            QString cellString;
            switch(data.value<Property::solverType>())
            {
            case Property::solverType_direct: cellString="Direct"; break;
            case Property::solverType_iterative: cellString="Iterative"; break;
            case Property::solverType_programControlled: cellString="Auto"; break;
            }
            QVariant cellStringVariant;
            cellStringVariant.setValue(cellString);
            return cellStringVariant;
        }
        else if(m_loadTypes.at(index.column())==Property::loadType_fieldParameters)
        {
            const QVector<double> &fieldParameters = data.value<QVector<double>>();
            QString cellString = QString("(%1, %2, %3, %4, %5, %6, %7, %8)")
                    .arg(fieldParameters.at(0))
                    .arg(fieldParameters.at(1))
                    .arg(fieldParameters.at(2))
                    .arg(fieldParameters.at(3))
                    .arg(fieldParameters.at(4))
                    .arg(fieldParameters.at(5))
                    .arg(fieldParameters.at(6))
                    .arg(fieldParameters.at(7));
            QVariant cellStringVariant;
            cellStringVariant.setValue(cellString);
            return cellStringVariant;
        }
        else if(m_loadTypes.at(index.column())==Property::loadType_lineSearch)
        {
            QString cellString;
            int val = data.toInt();
            switch(val)
            {
            case 0: cellString ="Program controlled"; break;
            case 1: cellString ="Custom"; break;
            }
            QVariant cellStringVariant;
            cellStringVariant.setValue(cellString);
            return cellStringVariant;
        }
        else if(m_loadTypes.at(index.column())==Property::loadType_lineSearchParameters)
        {
            const QVector<double> &fieldParameters = data.value<QVector<double>>();
            QString cellString = QString("(%1, %2)")
                    .arg(fieldParameters.at(0))
                    .arg(fieldParameters.at(1));
            QVariant cellStringVariant;
            cellStringVariant.setValue(cellString);
            return cellStringVariant;
        }
        else if(m_loadTypes.at(index.column())==Property::loadType_outputSettings)
        {
            const QVector<int> &outputControls = data.value<QVector<int>>();
            QString cellString = QString("(%1, %2, %3, %4)")
                    .arg(outputControls.at(0))
                    .arg(outputControls.at(1))
                    .arg(outputControls.at(2))
                    .arg(outputControls.at(3));
            QVariant cellStringVariant;
            cellStringVariant.setValue(cellString);
            return cellStringVariant;
        }
        else if(m_loadTypes.at(index.column())==Property::loadType_timeIncrementationParameters)
        {
            const QVector<int> &timeIncrementationParameters = data.value<QVector<int>>();
            QString cellString = QString("(%1, %2, %3, %4, %5, %6, %7, %8, %9, %10)")
                    .arg(timeIncrementationParameters.at(0))
                    .arg(timeIncrementationParameters.at(1))
                    .arg(timeIncrementationParameters.at(2))
                    .arg(timeIncrementationParameters.at(3))
                    .arg(timeIncrementationParameters.at(4))
                    .arg(timeIncrementationParameters.at(5))
                    .arg(timeIncrementationParameters.at(6))
                    .arg(timeIncrementationParameters.at(7))
                    .arg(timeIncrementationParameters.at(8))
                    .arg(timeIncrementationParameters.at(9));
            QVariant cellStringVariant;
            cellStringVariant.setValue(cellString);
            return cellStringVariant;
        }
        else if(m_loadTypes.at(index.column())==Property::loadType_storeResultsAt)
        {
            const QVector<int> &storeResultsAt = data.value<QVector<int>>();
            QString cellString = QString("(%1, %2)")
                    .arg(storeResultsAt.at(0))
                    .arg(storeResultsAt.at(1));
            QVariant cellStringVariant;
            cellStringVariant.setValue(cellString);
            return cellStringVariant;
        }
        else if(m_loadTypes.at(index.column())==Property::loadType_cutBackParameters)
        {
            const QVector<double> &cutBackParameters = data.value<QVector<double>>();
            QString cellString = QString("(%1, %2, %3, %4, %5, %6, %7, %8)")
                    .arg(cutBackParameters.at(0))
                    .arg(cutBackParameters.at(1))
                    .arg(cutBackParameters.at(2))
                    .arg(cutBackParameters.at(3))
                    .arg(cutBackParameters.at(4))
                    .arg(cutBackParameters.at(5))
                    .arg(cutBackParameters.at(6))
                    .arg(cutBackParameters.at(7));
            QVariant cellStringVariant;
            cellStringVariant.setValue(cellString);
            return cellStringVariant;
        }
        else if(m_loadTypes.at(index.column())==Property::loadType_fluxConvergence || m_loadTypes.at(index.column())==Property::loadType_solutionConvergence)
        {
            QString cellString;
            int val = data.toInt();
            switch(val)
            {
            case 0: cellString = QString("Off"); break;
            case 1: cellString = QString("On"); break;
            case 2: cellString = QString("Program controlled"); break;
            }
            QVariant cellStringVariant;
            cellStringVariant.setValue(cellString);
            return cellStringVariant;
        }
        else if(m_loadTypes.at(index.column())==Property::loadType_timeIncrementation ||
                m_loadTypes.at(index.column())==Property::loadType_cutBack)
        {
            QString cellString;
            int val = data.toInt();
            switch(val)
            {
            case 0: cellString = QString("Program controlled"); break;
            case 1: cellString = QString("Custom"); break;
            }
            QVariant cellStringVariant;
            cellStringVariant.setValue(cellString);
            return cellStringVariant;
        }
        else if(m_loadTypes.at(index.column())==Property::loadType_boltStatusDefinedBy)
        {
            QString cellString;
            switch(data.value<Property::boltStatusDefinedBy>())
            //switch(data.value<Property::defineBy>())
            {
            case Property::boltStatusDefinedBy_load: cellString="Load"; break;
            case Property::boltStatusDefinedBy_adjustment: cellString="Adjustment"; break;
            case Property::boltStatusDefinedBy_open: cellString="Open"; break;
            case Property::boltStatusDefinedBy_lock: cellString="Lock"; break;

            /*
            case Property::defineBy_load: cellString="Load"; break;
            case Property::defineBy_adjustment: cellString="Adjustment"; break;
            case Property::defineBy_open: cellString="Open"; break;
            case Property::defineBy_lock: cellString="Lock"; break;
            */
            }
            QVariant cellStringVariant;
            cellStringVariant.setValue(cellString);
            return cellStringVariant;
        }
        else if(m_loadTypes.at(index.column())==Property::loadType_forceMagnitude ||
                m_loadTypes.at(index.column())==Property::loadType_forceX ||
                m_loadTypes.at(index.column())==Property::loadType_forceY ||
                m_loadTypes.at(index.column())==Property::loadType_forceZ ||
                m_loadTypes.at(index.column())==Property::loadType_rotationalVelocityMagnitude ||
                m_loadTypes.at(index.column())==Property::loadType_rotationalVelocityX ||
                m_loadTypes.at(index.column())==Property::loadType_rotationalVelocityY ||
                m_loadTypes.at(index.column())==Property::loadType_rotationalVelocityZ ||
                m_loadTypes.at(index.column())==Property::loadType_accelerationMagnitude ||
                m_loadTypes.at(index.column())==Property::loadType_accelerationX ||
                m_loadTypes.at(index.column())==Property::loadType_accelerationY ||
                m_loadTypes.at(index.column())==Property::loadType_accelerationZ ||
                m_loadTypes.at(index.column())==Property::loadType_momentMagnitude ||
                m_loadTypes.at(index.column())==Property::loadType_momentX ||
                m_loadTypes.at(index.column())==Property::loadType_momentY ||
                m_loadTypes.at(index.column())==Property::loadType_momentZ ||
                m_loadTypes.at(index.column())==Property::loadType_pressureMagnitude ||
                m_loadTypes.at(index.column())==Property::loadType_displacementX ||
                m_loadTypes.at(index.column())==Property::loadType_displacementY ||
                m_loadTypes.at(index.column())==Property::loadType_displacementZ ||
                m_loadTypes.at(index.column())==Property::loadType_displacementMagnitude ||
                m_loadTypes.at(index.column())==Property::loadType_remoteRotationX ||
                m_loadTypes.at(index.column())==Property::loadType_remoteRotationY ||
                m_loadTypes.at(index.column())==Property::loadType_remoteRotationZ ||
                m_loadTypes.at(index.column())==Property::loadType_remoteRotationMagnitude ||
                m_loadTypes.at(index.column())==Property::loadType_temperatureMagnitude ||
                m_loadTypes.at(index.column())==Property::loadType_thermalConvectionFilmCoefficientMagnitude ||
                m_loadTypes.at(index.column())==Property::loadType_thermalConvectionReferenceTemperatureMagnitude ||
                m_loadTypes.at(index.column())==Property::loadType_thermalFluxMagnitude ||
                m_loadTypes.at(index.column())==Property::loadType_thermalFlowMagnitude ||
                m_loadTypes.at(index.column())==Property::loadType_thermalPowerMagnitude)
        {
            QString cellString = QString("%1").arg(data.toDouble());
            QVariant cellStringVariant;
            cellStringVariant.setValue(cellString);
            return cellStringVariant;
        }
        else if(m_loadTypes.at(index.column())==Property::loadType_boltForce ||
                m_loadTypes.at(index.column())==Property::loadType_boltAdjustment)
        {
            QString cellString = data.toString();
            QVariant cellStringVariant;
            cellStringVariant.setValue(cellString);
            return cellStringVariant;
        }
        else if(m_loadTypes.at(index.column())==Property::loadType_autoTimeStepping)
        {
            QVector<int> N = data.value<QVector<int>>();
            QString cellString = QString("(%1, %2, %3, %4)").arg(N[0]).arg(N[1]).arg(N[2]).arg(N[3]);
            QVariant cellStringVariant;
            cellStringVariant.setValue(cellString);
            return cellStringVariant;
        }
        else if(m_loadTypes.at(index.column())==Property::loadType_modelChange)
        {
            QString cellString;
            int val = data.toInt();
            switch(val)
            {
            case 0: cellString = QString("Inactive N.A."); break;
            case 1: cellString = QString("Add"); break;
            case -1: cellString = QString("Remove"); break;
            }
            QVariant cellStringVariant;
            cellStringVariant.setValue(cellString);
            return cellStringVariant;
        }
        else if(m_loadTypes.at(index.column())==Property::loadType_analysisType)
        {
            QString cellString;
            Property::analysisType analysisType = data.value<Property::analysisType>();
            switch(analysisType)
            {
            case Property::analysisType_structural: cellString = QString("Structural"); break;
            case Property::analysisType_thermal: cellString = QString("Thermal"); break;
            case Property::analysisType_modal: cellString = QString("Modal"); break;
            case Property::analysisType_frequencyResponse: cellString = QString("Frequency response"); break;
            case Property::analysisType_uncoupledTemperatureDisplacement: cellString = QString("Uncoupled temperature displacement"); break;
            case Property::analysisType_coupledTemperatureDisplacement: cellString = QString("Temperature displacement"); break;
            }
            QVariant cellStringVariant;
            cellStringVariant.setValue(cellString);
            return cellStringVariant;
        }
        else if(m_loadTypes.at(index.column())==Property::loadType_timeIntegration)
        {
            QString cellString;
            Property::timeIntegration timeIntegration = data.value<Property::timeIntegration>();
            switch(timeIntegration)
            {
            case Property::timeIntegration_steadyState: cellString = QString("Static"); break;
            case Property::timeIntegration_transient: cellString = QString("Transient"); break;
            }
            QVariant cellStringVariant;
            cellStringVariant.setValue(cellString);
            return cellStringVariant;
        }
    }
    else if (role == Qt::EditRole)
    {
        return m_data[index.row()]->at(index.column());
    }
    else if (role == Qt::TextAlignmentRole)
    {
        return Qt::AlignCenter;
    }
    return QVariant();
}

//! ------------------
//! function: setData
//! details:
//! ------------------
bool CustomTableModel::setData(const QModelIndex &index, const QVariant &value, int role)
{
    if(index.isValid() && role == Qt::EditRole)
    {
        m_data[index.row()]->replace(index.column(), value);
        cout<<"CustomTableModel::setData()->____emitting dataChanged()____"<<endl;
        emit dataChanged(index, index);
        //this->autoUpdateTable(index,index);
        return true;
    }
    return false;
}

//! --------------------
//! function: setDataRC
//! details:
//! --------------------
bool CustomTableModel::setDataRC(const QVariant &value, int row, int column, int role)
{
    if(row>=0 && column >=0)
    {
        if(role == Qt::EditRole)
        {
            m_data[row]->replace(column, value);
            QModelIndex index = this->createIndex(row,column);
            emit dataChanged(index, index);
            return true;
        }
    }
    return false;
}

//! --------------------------------------------
//! function: strings for the horizontal header
//! details:
//! --------------------------------------------
QString CustomTableModel::getHeaderString(int section) const
{
    QString horizontalHeaderString = QString("");
    switch(m_loadTypes.at(section))
    {
    case Property::loadType_temperatureMagnitude: horizontalHeaderString = "T"; break;
    case Property::loadType_thermalConvectionFilmCoefficientMagnitude: horizontalHeaderString = "Film coeff"; break;
    case Property::loadType_thermalConvectionReferenceTemperatureMagnitude: horizontalHeaderString = "Ref temperature"; break;
    case Property::loadType_thermalFluxMagnitude: horizontalHeaderString = "Flux"; break;
    case Property::loadType_thermalFlowMagnitude: horizontalHeaderString = "Flow"; break;
    case Property::loadType_thermalPowerMagnitude: horizontalHeaderString = "Heat generation"; break;
    case Property::loadType_stepNumber: horizontalHeaderString = "Step #"; break;
    case Property::loadType_displacementX: horizontalHeaderString = "dX"; break;
    case Property::loadType_displacementY: horizontalHeaderString = "dY"; break;
    case Property::loadType_displacementZ: horizontalHeaderString = "dZ"; break;
    case Property::loadType_displacementMagnitude: horizontalHeaderString = "displ"; break;
    case Property::loadType_forceX: horizontalHeaderString = "Fx"; break;
    case Property::loadType_forceY: horizontalHeaderString = "Fy"; break;
    case Property::loadType_forceZ: horizontalHeaderString = "Fz"; break;
    case Property::loadType_forceMagnitude: horizontalHeaderString = "F"; break;
    case Property::loadType_rotationalVelocityX: horizontalHeaderString = "Rx"; break;
    case Property::loadType_rotationalVelocityY: horizontalHeaderString = "Ry"; break;
    case Property::loadType_rotationalVelocityZ: horizontalHeaderString = "Rz"; break;
    case Property::loadType_remoteRotationMagnitude: horizontalHeaderString = "Rot [°]"; break;
    case Property::loadType_remoteRotationX: horizontalHeaderString = "Rot X [°]"; break;
    case Property::loadType_remoteRotationY: horizontalHeaderString = "Rot Y [°]"; break;
    case Property::loadType_remoteRotationZ: horizontalHeaderString = "Rot Z [°]"; break;
    case Property::loadType_rotationalVelocityMagnitude: horizontalHeaderString = "R"; break;
    case Property::loadType_accelerationX: horizontalHeaderString = "Ax"; break;
    case Property::loadType_accelerationY: horizontalHeaderString = "Ay"; break;
    case Property::loadType_accelerationZ: horizontalHeaderString = "Az"; break;
    case Property::loadType_accelerationMagnitude: horizontalHeaderString = "A"; break;
    case Property::loadType_momentX: horizontalHeaderString = "Mx"; break;
    case Property::loadType_momentY: horizontalHeaderString = "My"; break;
    case Property::loadType_momentZ: horizontalHeaderString = "Mz"; break;
    case Property::loadType_momentMagnitude: horizontalHeaderString = "M"; break;
    case Property::loadType_pressureMagnitude: horizontalHeaderString = "P"; break;
    case Property::loadType_stepEndTime: horizontalHeaderString ="End time[s]"; break;
    case Property::loadType_time: horizontalHeaderString ="Time[s]"; break;
    case Property::loadType_solverType: horizontalHeaderString ="Solver type"; break;
    case Property::loadType_autoTimeStepping: horizontalHeaderString ="Time stepping policy"; break;
    case Property::loadType_thermalConditionTemperature: horizontalHeaderString ="T [K]"; break;
    case Property::loadType_boltStatusDefinedBy: horizontalHeaderString = "Define by"; break;
    case Property::loadType_boltForce: horizontalHeaderString = "Force"; break;
    case Property::loadType_boltAdjustment: horizontalHeaderString = "Adjustment"; break;
    case Property::loadType_fieldParameters: horizontalHeaderString = "Field parameters"; break;
    case Property::loadType_fluxConvergence: horizontalHeaderString = "Flux convergence"; break;
    case Property::loadType_solutionConvergence: horizontalHeaderString = "Solution convergence"; break;
    case Property::loadType_timeIncrementationParameters: horizontalHeaderString = "Time incrementation parameters"; break;
    case Property::loadType_timeIncrementation: horizontalHeaderString = "Time incrementation type"; break;
    case Property::loadType_cutBack: horizontalHeaderString = "Cutback type"; break;
    case Property::loadType_cutBackParameters: horizontalHeaderString = "Cutback parameters"; break;
    case Property::loadType_lineSearch: horizontalHeaderString = "Line search type"; break;
    case Property::loadType_lineSearchParameters: horizontalHeaderString = "Line search parameters"; break;
    case Property::loadType_outputSettings: horizontalHeaderString = "Output controls"; break;
    case Property::loadType_storeResultsAt: horizontalHeaderString = "Store results at"; break;
    case Property::loadType_modelChange: horizontalHeaderString = "Model change"; break;
    case Property::loadType_analysisType: horizontalHeaderString = "Analysis type"; break;
    case Property::loadType_timeIntegration: horizontalHeaderString = "Static/Transient"; break;
        //! empty string for "loadType_none"
    case Property::loadType_none: horizontalHeaderString = ""; break;
    }
    horizontalHeaderString.prepend(QString("[%1]").arg(section));
    return horizontalHeaderString;
}

//! ------------------------
//! function: append rows
//! details:  to be removed
//! ------------------------
bool CustomTableModel::appendRow(int count, const QModelIndex &parent)
{
    Q_UNUSED(parent);

    //cout<<"CustomTableModel::appendRow()->____function called____"<<endl;
    this->beginInsertRows(QModelIndex(),1,2);
    for(int i=0;i<count;i++)
    {
        QVector<QVariant> *dataVec = new QVector<QVariant>(this->columnCount());
        m_data.append(dataVec);
    }
    this->endInsertRows();
    return true;
}

//! ---------------------
//! function: insertRows
//! details:
//! ---------------------
bool CustomTableModel::insertRows(int row, int count, const QModelIndex &parent)
{
    cout<<"\\--------------------------------------------------------------\\"<<endl;
    cout<<"\\ CustomTableModel::insertRows->____function called____________\\"<<endl;

    Q_UNUSED(parent)
    this->beginInsertRows(QModelIndex(),row,row+count-1);

    //! retrieve the content of the last row
    QVector<QVariant> *theLastDataVec = m_data.at(row-1);

    //! a time increment
    double theTimeIncrement = 1.0;

    for(int c=row; c<row+count; c++)
    {
        //! the copied row
        QVector<QVariant> *dataVec = new QVector<QVariant>;
        QVariant data;
        //! increment the time step number
        //cout<<"Pushing back step no: "<<theLastDataVec->at(0).toInt()+c-row+1<<endl;
        data.setValue(theLastDataVec->at(0).toInt()+c-row+1);
        dataVec->push_back(data);
        //! increment the step end time
        data.setValue(theLastDataVec->at(1).toDouble()+(c-row+1)*theTimeIncrement);
        //cout<<"Pushing back time: "<<theLastDataVec->at(1).toDouble()+(c-row+1)*theTimeIncrement<<endl;
        dataVec->push_back(data);
        for(int i=2;i<theLastDataVec->length();i++)
        {
            //cout<<"Pushing back: "<<theLastDataVec->at(i).typeName()<<endl;
            dataVec->push_back(theLastDataVec->at(i));
        }
        //! adding the row to the data
        m_data.append(dataVec);
    }
    this->endInsertRows();

    //! ----------------------------------------------------------------
    //! new code - inform the viewer about the position of the new data
    //! ----------------------------------------------------------------
    QModelIndex topLeftIndex = this->createIndex(row,theLastDataVec->length()-1);
    QModelIndex bottomRightIndex = this->createIndex(row+count-1,theLastDataVec->length()-1);

    cout<<"\\ top left(row = "<<row<<", col = "<<theLastDataVec->length()-1<<") bottom right(row = "<<row+count-1<<", col = "<<theLastDataVec->length()-1<<")"<<endl;
    cout<<"\\--------------------------------------------------------------\\"<<endl;

    emit dataChanged(topLeftIndex,bottomRightIndex);
    return true;
}

//! ---------------------
//! function: removeRows
//! details:
//! ---------------------
bool CustomTableModel::removeRows(int first, int last, const QModelIndex &parent)
{
    Q_UNUSED(parent);
    //cout<<"CustomTableModel::removeRows->____function called____"<<endl;
    this->beginRemoveRows(QModelIndex(),first,last);
    for(int k=last; k>=first; k--)m_data.removeAt(k);
    this->endRemoveRows();
    return true;
}

//! ---------------------------------------
//! function: minXVal
//! details:  this is for the graph viewer
//! ---------------------------------------
double CustomTableModel::minXVal() const
{
    double Xmin = 1e20;
    for(int i=0; i<m_data.length(); i++)
    {
        //! at(1) for times
        double value = m_data.at(i)->at(1).toDouble();
        if(value<=Xmin)Xmin=value;
    }
    return Xmin;
}

//! --------------------------------
//! function: maxXVal
//! details:  this is for the graph
//! --------------------------------
double CustomTableModel::maxXVal() const
{
    double Xmax = -1e20;
    for(int i=0; i<m_data.length(); i++)
    {
        //! at(1) for times
        double value = m_data.at(i)->at(1).toDouble();
        if(value>=Xmax)Xmax=value;
    }
    return Xmax;
}

//! --------------------------------
//! function: maxYVal
//! details:  this is for the graph
//! --------------------------------
double CustomTableModel::maxVal() const
{
    //! this is for the graph
    double max = -1e20;
    for(int i=0; i<m_data.length(); i++)
    {
        int ncol = m_data.at(i)->length();
        for(int j = 3; j<ncol; j++)
        {
            double value = m_data.at(i)->at(j).toDouble();
            if(value>=max)max = value;
        }
    }
    cout<<"____=>The max value is: "<<max<<"<=____"<<endl;
    return max;
}

//! -----------------------------
//! function: minVal
//! details:  this for the graph
//! -----------------------------
double CustomTableModel::minVal() const
{
    //! this is for the graph
    double min = 1e20;
    for(int i=0; i<m_data.length(); i++)
    {
        int ncol = m_data.at(i)->length();
        for(int j = 3; j<ncol; j++)
        {
            double value = m_data.at(i)->at(j).toDouble();
            if(value<=min) min = value;
        }
    }
    return min;
}

//! -------------------
//! function: rowCount
//! details:
//! -------------------
int CustomTableModel::rowCount(const QModelIndex &parent) const
{
    Q_UNUSED(parent)
    return m_data.count();
}

//! ----------------------
//! function: columnCount
//! details:
//! ----------------------
int CustomTableModel::columnCount(const QModelIndex &parent) const
{
    Q_UNUSED(parent)
    return m_data.at(0)->size();
    //return m_columnCount;
}

//! ----------------------------------------------------------------
//! function: flags
//! details:  a not editable cell is created by setting the QString
//!           value "N\A"
//! ----------------------------------------------------------------
Qt::ItemFlags CustomTableModel::flags(const QModelIndex &index) const
{
    if(index.column()==0)
    {
        //! the first column (containing step numbers) is not editable
        return QAbstractItemModel::flags(index) | Qt::ItemIsSelectable;
    }
    else
    {
        //! if the cell contains the QString "N/A" set the cell as not editable
        if(index.data(Qt::EditRole).canConvert<QString>())
        {
            if(index.data(Qt::EditRole).toString()=="N/A")
                return QAbstractItemModel::flags(index) | Qt::ItemIsSelectable;
        }
        return QAbstractItemModel::flags(index) | Qt::ItemIsEditable;
    }
}

//! ------------------------
//! function: insertColumns
//! details:
//! ------------------------
bool CustomTableModel::insertColumns(int column, int count, const QModelIndex &parent)
{
    cout<<"\\--------------------------------------------------------------\\"<<endl;
    cout<<"\\ CustomTableModel::insertColumns()->____function called"<<endl;
    cout<<"\\ inserting column starting from col position: "<<column<<"____"<<endl;

    Q_UNUSED(parent)

    //! append the new load type
    m_loadTypes.insert(column,loadToAppend.type());

    //! calculate the position
    int first = column;
    int last = column+count-1;
    this->beginInsertColumns(QModelIndex(),first,last);

    //! -----------------------------------------------
    //! handle data - some definitions:
    //! Nrow = m_data.size(); N = aLoad.values().size;
    //! 1) N = 0;
    //! A load with an empty data vector is provided
    //! -----------------------------------------------
    int Nrows = this->rowCount();
    if(loadToAppend.values().isEmpty())
    {
        double value = 0.0;
        QVariant data;
        data.setValue(value);
        for(int k=0; k<Nrows; k++)
        {
            QVector<QVariant> *dataVec = m_data.at(k);
            dataVec->insert(column,data);
            m_data.replace(k,dataVec);
        }
    }
    //! -------------------------------------------
    //! 2) N >= Nrow;
    //! A load with a number of rows equal or
    //! greater than the m_data.size() is provided
    //! -------------------------------------------
    else if(loadToAppend.values().size()>=m_data.length())
    {
        for(int k=0; k<Nrows; k++)
        {
            QVector<QVariant> *dataVec = m_data.at(k);
            dataVec->insert(column,loadToAppend.values().at(k));
            m_data.replace(k,dataVec);
        }
    }
    //! ---------------------------------------
    //! 3) N>=0; N<Nrow
    //! A load with a nonzero row number but
    //! smaller than m_data.size() is provided
    //! ---------------------------------------
    else if(loadToAppend.values().size()!=0 && loadToAppend.values().size()<m_data.length())
    {
        int k;
        for(k=0; k<loadToAppend.values().size(); k++)
        {
            QVector<QVariant> *dataVec = m_data.at(k);
            dataVec->insert(column,loadToAppend.values().at(k));
            m_data.replace(k,dataVec);
        }
        QVariant lastAvailableValue = loadToAppend.values().at(k-1);
        for(int k=loadToAppend.values().size(); k<Nrows; k++)
        {
            QVector<QVariant> *dataVec = m_data.at(k);
            dataVec->insert(column,lastAvailableValue);
            m_data.replace(k,dataVec);
        }
    }

    //! -----------------------------
    //! update the number of columns
    //! -----------------------------
    //m_columnCount++;
    m_columnCount = m_data.at(0)->size();
    this->endInsertColumns();

    //! -----------------------------------------------------
    //! inform the viewer about the position of the new data
    //! -----------------------------------------------------
    QModelIndex topLeftIndex = this->createIndex(0,first);
    QModelIndex bottomRightIndex = this->createIndex(this->rowCount()-1,last);

    cout<<"\\ top left(row = "<<0<<", col = "<<first<<") bottom right(row = "<<this->rowCount()-1<<", col = "<<last<<")"<<endl;
    cout<<"\\--------------------------------------------------------------\\"<<endl;

    emit dataChanged(topLeftIndex,bottomRightIndex);
    return true;
}

//! --------------------------
//! function: setLoadToInsert
//! details:
//! --------------------------
void CustomTableModel::setLoadToInsert(load aload)
{
    //cout<<"CustomTableModel::setLoadToInsert->____function called____"<<endl;
    loadToAppend = aload;
}

//! ------------------------
//! function: removeColumns
//! details:
//! ------------------------
bool CustomTableModel::removeColumns(int column, int count, const QModelIndex &parent)
{
    cout<<"CustomTableModel::removeColumns()->____function called: removing____start: "<<column<<" count: "<<count<<"____"<<endl;
    Q_UNUSED(parent)
    this->beginRemoveColumns(QModelIndex(),column,column+count-1);

    //! -----------------------
    //! modify the data
    //! [1] erase the loadType
    //! -----------------------
    //cout<<"CustomTableModel::removeColumns()->____initial number of load types: "<<m_loadTypes.length()<<"____"<<endl;
    for(int i=column;i<column+count;i++)
    {
        cout<<"CustomTableModel::removeColumns()->____removing load type____"<<endl;
        //! -------------------------------------------------------------------
        //! warning: remove the element of the list keeping the position fixed
        //! -------------------------------------------------------------------
        m_loadTypes.removeAt(column);
    }

    //! --------------------
    //! [2] change the data
    //! --------------------
    for(int i=0; i<m_data.size(); i++)
    {
        QVector<QVariant> *dataVec = m_data.at(i);
        dataVec->remove(column,count);
        m_data.replace(i,dataVec);
    }
    m_columnCount = m_data.at(0)->size();
    this->endRemoveColumns();
    cout<<"CustomTableModel::removeColumns()->____final number of columns: "<<m_data.at(0)->size()<<"____"<<endl;
    return true;
}

//! ----------------------------------------------
//! function: appendColumn
//! details:  append a column to the tabular data
//! ----------------------------------------------
bool CustomTableModel::appendColumn(const load &aLoad, const QModelIndex &parent)
{
    cout<<"\\--------------------------------------------------------------\\"<<endl;
    cout<<"\\ CustomTableModel::appendColumn()->____function called________\\"<<endl;

    Q_UNUSED(parent)

    //! append the new load type
    m_loadTypes.append(aLoad.type());

    //! calculate the position
    int position = this->columnCount();
    int Nrows = this->rowCount();

    //cout<<"CustomTableModel::appendColumn()->____number of row in analysis settings: "<<Nrows<<"____"<<endl;
    //cout<<"CustomTableModel::appendColumn()->____new column number: "<<position<<"____"<<endl;

    this->beginInsertColumns(QModelIndex(),position,position);

    //! -----------------------------------------------
    //! handle data - some definitions:
    //! Nrow = m_data.size(); N = aLoad.values().size;
    //! 1) N = 0;
    //! A load with an empty data vector is provided
    //! -----------------------------------------------
    if(aLoad.values().isEmpty())
    {
        double value = 0.0;
        QVariant data;
        data.setValue(value);
        for(int k=0; k<Nrows; k++)
        {
            QVector<QVariant> *dataVec = m_data.at(k);
            dataVec->push_back(data);
            m_data.replace(k,dataVec);
            //cout<<"____dataVec size after insertion: "<<m_data.at(k)->size()<<endl;
        }
    }
    //! -------------------------------------------
    //! 2) N >= Nrow;
    //! A load with a number of rows equal or
    //! greater than the m_data.size() is provided
    //! -------------------------------------------
    else if(aLoad.values().size()>=m_data.length())
    {
        for(int k=0; k<Nrows; k++)
        {
            QVector<QVariant> *dataVec = m_data.at(k);
            dataVec->push_back(aLoad.values().at(k));
            m_data.replace(k,dataVec);
            //cout<<"____dataVec size after insertion: "<<m_data.at(k)->size()<<endl;
        }
    }
    //! ---------------------------------------
    //! 3) N>=0; N<Nrow
    //! A load with a nonzero row number but
    //! smaller than m_data.size() is provided
    //! ---------------------------------------
    else if(aLoad.values().size()!=0 && aLoad.values().size()<m_data.length())
    {
        int k;
        for(k=0; k<aLoad.values().size(); k++)
        {
            QVector<QVariant> *dataVec = m_data.at(k);
            dataVec->push_back(aLoad.values().at(k));
            m_data.replace(k,dataVec);
        }
        QVariant lastAvailableValue = aLoad.values().at(k-1);
        for(int k=aLoad.values().size(); k<Nrows; k++)
        {
            QVector<QVariant> *dataVec = m_data.at(k);
            dataVec->push_back(lastAvailableValue);
            m_data.replace(k,dataVec);
        }
    }

    //! -----------------------------
    //! update the number of columns
    //! -----------------------------
    m_columnCount = m_data.at(0)->size();
    this->endInsertColumns();

    //! -----------------------------------------------------
    //! inform the viewer about the position of the new data
    //! -----------------------------------------------------
    QModelIndex topLeftIndex = this->createIndex(0,position);
    QModelIndex bottomRightIndex = this->createIndex(Nrows-1,position);

    cout<<"\\ top left(row = "<<0<<", col = "<<position<<") bottom right(row = "<<Nrows-1<<", col = "<<position<<")"<<endl;
    cout<<"\\--------------------------------------------------------------\\"<<endl;

    emit dataChanged(topLeftIndex,bottomRightIndex);
    return true;
}

//! ----------------------------------------------------------------
//! function: -----
//! details:
//! ----------------------------------------------------------------
void CustomTableModel::addMapping(QString color, QRect area)
{
    m_mapping.insertMulti(color, area);
}

//! --------------------
//! function: getColumn
//! details:
//! --------------------
load CustomTableModel::getColumn(int col) const
{
    QVector<QVariant> vecVariant;
    int NValues = this->rowCount();
    bool isInvalid = true;
    for(int i=0;i<NValues;i++)
    {
        QVector<QVariant> *vecData = m_data.at(i);
        QVariant data = vecData->at(col);

        if(!data.isValid())
        {
            //! this could occur in case of a corrupted database
            //! the "null" or "none" loadType, is defined as
            //! type => Property::loadType_none
            //! values => QVector<QVariant> in which the value of QVariant is int "0"
            vecVariant.push_back(0);
            isInvalid = true;
        }
        else
        {
            vecVariant.push_back(data);
            isInvalid = false;
        }
    }

    //! ----------------
    //! returning value
    //! ----------------
    if(isInvalid)
    {
        load theLoad(vecVariant, Property::loadType_none);
        return theLoad;
    }
    else
    {
        load theLoad(vecVariant, m_loadTypes.at(col));
        return theLoad;
    }
}
