//! ----------------
//! custom includes
//! ----------------
#include "tableviewclassitemdelegate.h"
#include "tabulardatacolumns.h"
#include "simulationnodeclass.h"
#include "customtablemodel.h"
#include "property.h"

//! ---
//! Qt
//! ---
#include <QLineEdit>
#include <QComboBox>
#include <QTableView>

//! -----------------------
//! function: constructor
//! details:
//! -----------------------
tableViewClassItemDelegate::tableViewClassItemDelegate(QWidget *parent): QStyledItemDelegate(parent)
{
    //!cout<<"tableViewClassItemDelegate::tableViewClassItemDelegate()->____constructor called____"<<endl;
}

//! -----------------------
//! function: createEditor
//! details:
//! -----------------------
QWidget* tableViewClassItemDelegate::createEditor(QWidget *parent, const QStyleOptionViewItem &/*option*/, const QModelIndex &index) const
{
    if(!index.isValid()) return 0;

    QVariant data = index.model()->data(index,Qt::EditRole);
    cout<<"tableViewClassItemDelegate::createEditor()->____creating editor for type: "<<data.typeName()<<"____"<<endl;

    if(strcmp(data.typeName(),"Property::timeIntegration")==0)
    {
        QComboBox *editor = new QComboBox(parent);
        editor->addItem("Static",Property::timeIntegration_steadyState);
        editor->addItem("Transient",Property::timeIntegration_transient);
        connect(editor,SIGNAL(currentIndexChanged(int)),this, SLOT(commitAndCloseComboBox()));
        return editor;
    }
    if(strcmp(data.typeName(),"Property::analysisType")==0)
    {
        QComboBox *editor = new QComboBox(parent);
        editor->addItem("Structural",Property::analysisType_structural);
        editor->addItem("Thermal",Property::analysisType_thermal);
        editor->addItem("Modal",Property::analysisType_modal);
        editor->addItem("Frequency response",Property::analysisType_frequencyResponse);
        editor->addItem("Uncoupled temperature displacement",Property::analysisType_uncoupledTemperatureDisplacement);
        editor->addItem("Coupled temperature displacement",Property::analysisType_coupledTemperatureDisplacement);
        connect(editor,SIGNAL(currentIndexChanged(int)),this, SLOT(commitAndCloseComboBox()));
        return editor;
    }
    if(strcmp(data.typeName(),"Property::modelChangeActivationStatus")==0)
    {
        QComboBox *editor = new QComboBox(parent);
        QVariant data;
        data.setValue(Property::modelChangeActivationStatus_Remove);        //! index "0"
        editor->addItem("Remove",data);
        data.setValue(Property::modelChangeActivationStatus_Inactive);      //! index "1"
        editor->addItem("Inactive",data);
        data.setValue(Property::modelChangeActivationStatus_Add);           //! index "2"
        editor->addItem("Add",data);
        connect(editor,SIGNAL(currentIndexChanged(int)),this, SLOT(commitAndCloseComboBox()));
        return editor;
    }
    if(strcmp(data.typeName(),"Property::solverType")==0)
    {
        QComboBox *editor = new QComboBox(parent);
        //editor->clear();
        QVariant data;
        data.setValue(Property::solverType_programControlled);
        editor->addItem("Auto",data);
        data.setValue(Property::solverType_direct);
        editor->addItem("Direct",data);
        data.setValue(Property::solverType_iterative);
        editor->addItem("Iterative",data);
        connect(editor,SIGNAL(currentIndexChanged(int)),this, SLOT(commitAndCloseComboBox()));
        return editor;
    }
    if(strcmp(data.typeName(),"Property::boltStatusDefinedBy")==0)
    {
        QComboBox *editor = new QComboBox(parent);
        QVariant data;
        data.setValue(Property::boltStatusDefinedBy_load);
        editor->addItem("Load",data);
        data.setValue(Property::boltStatusDefinedBy_adjustment);
        editor->addItem("Adjustment",data);
        data.setValue(Property::boltStatusDefinedBy_open);
        editor->addItem("Open",data);
        data.setValue(Property::boltStatusDefinedBy_lock);
        editor->addItem("Lock",data);
        return editor;
    }
    if(strcmp(data.typeName(),"double")==0)
    {
        QLineEdit *editor = new QLineEdit(parent);
        QDoubleValidator *doubleValidator = new QDoubleValidator();
        editor->setValidator(doubleValidator);
        return editor;
    }
    else return Q_NULLPTR;
}

//! ------------------------
//! function: setEditorData
//! details:
//! ------------------------
void tableViewClassItemDelegate::setEditorData(QWidget *editor, const QModelIndex &index) const
{
    QVariant data = index.data(Qt::EditRole);

    if(strcmp(data.typeName(),"Property::timeIntegration")==0)
    {
        Property::timeIntegration value = data.value<Property::timeIntegration>();
        QComboBox *comboBox = static_cast<QComboBox*>(editor);
        switch(value)
        {
        case Property::timeIntegration_steadyState: comboBox->setCurrentIndex(0); break;
        case Property::timeIntegration_transient: comboBox->setCurrentIndex(1); break;
        }
    }
    if(strcmp(data.typeName(),"Property::analysisType")==0)
    {
        Property::analysisType value = data.value<Property::analysisType>();
        QComboBox *comboBox = static_cast<QComboBox*>(editor);
        switch(value)
        {
        case Property::analysisType_structural: comboBox->setCurrentIndex(0); break;
        case Property::analysisType_thermal: comboBox->setCurrentIndex(1); break;
        case Property::analysisType_modal: comboBox->setCurrentIndex(2); break;
        case Property::analysisType_frequencyResponse: comboBox->setCurrentIndex(3); break;
        case Property::analysisType_uncoupledTemperatureDisplacement: comboBox->setCurrentIndex(4); break;
        case Property::analysisType_coupledTemperatureDisplacement: comboBox->setCurrentIndex(5); break;
        }
    }
    if(strcmp(data.typeName(),"Property::modelChangeActivationStatus")==0)
    {
        Property::modelChangeActivationStatus value = data.value<Property::modelChangeActivationStatus>();
        QComboBox *comboBox = static_cast<QComboBox*>(editor);
        switch(value)
        {
        case Property::modelChangeActivationStatus_Remove: comboBox->setCurrentIndex(0); break;     //! "0" <=> "Remove"
        case Property::modelChangeActivationStatus_Inactive: comboBox->setCurrentIndex(1); break;   //! "1" <=> "Inactive"
        case Property::modelChangeActivationStatus_Add: comboBox->setCurrentIndex(2); break;        //! "2" <=> "Active"
        }
    }
    if(strcmp(data.typeName(),"Property::solverType")==0)
    {
        Property::solverType value = data.value<Property::solverType>();
        QComboBox *comboBox = static_cast<QComboBox*>(editor);
        switch(value)
        {
        case Property::solverType_programControlled: comboBox->setCurrentIndex(0); break;
        case Property::solverType_direct: comboBox->setCurrentIndex(1); break;
        case Property::solverType_iterative: comboBox->setCurrentIndex(2); break;
        }
    }
    if(strcmp(data.typeName(),"Property::boltStatusDefinedBy")==0)
    {
        Property::boltStatusDefinedBy value = data.value<Property::boltStatusDefinedBy>();

        QComboBox *comboBox = static_cast<QComboBox*>(editor);
        switch(value)
        {
        case Property::boltStatusDefinedBy_load: comboBox->setCurrentIndex(0); break;
        case Property::boltStatusDefinedBy_adjustment: comboBox->setCurrentIndex(1); break;
        case Property::boltStatusDefinedBy_open: comboBox->setCurrentIndex(2); break;
        case Property::boltStatusDefinedBy_lock: comboBox->setCurrentIndex(3); break;
        }
        connect(editor,SIGNAL(currentIndexChanged(int)),this, SLOT(commitAndCloseBoltStatusDefineBy()));
    }
    if(strcmp(data.typeName(),"double")==0)
    {
        double value = index.data(Qt::EditRole).toDouble();
        QLineEdit *lineEdit = static_cast<QLineEdit*>(editor);
        lineEdit->setText(QString("%1").arg(value));
        connect(editor,SIGNAL(editingFinished()),this, SLOT(commitAndCloseLineEdit()));
    }
}

//! -----------------------
//! function: setModelData
//! details:
//! -----------------------
void tableViewClassItemDelegate::setModelData(QWidget *editor, QAbstractItemModel *model, const QModelIndex &index) const
{
    if(!index.isValid()) return;

    QVariant dataContent = index.model()->data(index, Qt::EditRole);

    QVariant data;
    if(strcmp(dataContent.typeName(),"Property::timeIntegration")==0)
    {
        QComboBox *comboBox = static_cast<QComboBox*>(editor);
        switch(comboBox->currentIndex())
        {
        case 0: data.setValue(Property::timeIntegration_steadyState); break;
        case 1: data.setValue(Property::timeIntegration_transient); break;
        }
    }
    if(strcmp(dataContent.typeName(),"Property::analysisType")==0)
    {
        QComboBox *comboBox = static_cast<QComboBox*>(editor);
        switch(comboBox->currentIndex())
        {
        case 0: data.setValue(Property::analysisType_structural); break;
        case 1: data.setValue(Property::analysisType_thermal); break;
        case 2: data.setValue(Property::analysisType_modal); break;
        case 3: data.setValue(Property::analysisType_frequencyResponse); break;
        case 4: data.setValue(Property::analysisType_uncoupledTemperatureDisplacement); break;
        case 5: data.setValue(Property::analysisType_coupledTemperatureDisplacement); break;
        }
    }
    if(strcmp(dataContent.typeName(),"Property::solverType")==0)
    {
        QComboBox *comboBox = static_cast<QComboBox*>(editor);
        switch(comboBox->currentIndex())
        {
        case 0: data.setValue(Property::solverType_programControlled); break;
        case 1: data.setValue(Property::solverType_direct); break;
        case 2: data.setValue(Property::solverType_iterative); break;
        }
    }
    if(strcmp(dataContent.typeName(),"Property::modelChangeActivationStatus")==0)
    {
        QComboBox *comboBox = static_cast<QComboBox*>(editor);
        switch(comboBox->currentIndex())
        {
        case 0: data.setValue(Property::modelChangeActivationStatus_Remove); break;
        case 1: data.setValue(Property::modelChangeActivationStatus_Inactive); break;
        case 2: data.setValue(Property::modelChangeActivationStatus_Add); break;
        }
    }
    if(strcmp(dataContent.typeName(),"Property::boltStatusDefinedBy")==0)
    {
        QComboBox *comboBox = static_cast<QComboBox*>(editor);
        switch(comboBox->currentIndex())
        {
        case 0: data.setValue(Property::boltStatusDefinedBy_load); break;
        case 1: data.setValue(Property::boltStatusDefinedBy_adjustment); break;
        case 2: data.setValue(Property::boltStatusDefinedBy_open); break;
        case 3: data.setValue(Property::boltStatusDefinedBy_lock); break;
        }
    }
    if(strcmp(dataContent.typeName(),"double")==0)
    {
        QLineEdit *lineEdit = static_cast<QLineEdit*>(editor);
        double value = lineEdit->text().toDouble();
        data.setValue(value);
    }
    model->setData(index,data,Qt::EditRole);
}

//! ---------------------------------
//! function: commitAndCloseLineEdit
//! details:
//! ---------------------------------
void tableViewClassItemDelegate::commitAndCloseLineEdit()
{
    QLineEdit *editor = qobject_cast<QLineEdit *>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
}

//! ---------------------------------
//! function: commitAndCloseComboBox
//! details:
//! ---------------------------------
void tableViewClassItemDelegate::commitAndCloseComboBox()
{
    QComboBox *editor = qobject_cast<QComboBox *>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
}

//! -------------------------------------------
//! function: commitAndCloseBoltStatusDefineBy
//! details:
//! -------------------------------------------
void tableViewClassItemDelegate::commitAndCloseBoltStatusDefineBy()
{
    QComboBox *editor = qobject_cast<QComboBox *>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit boltStatusDefineByChanged();
}
