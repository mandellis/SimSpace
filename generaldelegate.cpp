//! ----------------
//! custom includes
//! ----------------
#include "generaldelegate.h"
#include "property.h"
#include "mydefines.h"
#include "myenumvariables.h"
#include "qextendedstandarditem.h"
#include "shapeselector.h"
#include "meshselector.h"
#include "directionselector.h"
#include "lineedit.h"
#include "qbackgroundevent.h"
#include "tools.h"
#include "optionsWidget/colorselector.h"
#include <itemselector.h>
#include <qfileselect.h>

//! ---
//! Qt
//! ---
#include <QStandardItemModel>
#include <QComboBox>
#include <QLineEdit>
#include <QPainter>
#include <QTreeView>
#include <QFileDialog>
#include <QSpinBox>
#include <QPushButton>
#include <QSlider>


#define MIN_MESH_REDUCTION  0.40
#define MAX_MESH_REDUCTION  0.99

using namespace std;

const QString sliderStyleSheet = QString("QSlider::groove:horizontal { "
                                 "border: 1px solid #999999; "
                                 "background: white; "
                                 "} "
                                 "QSlider::handle:horizontal { "
                                 "background: qlineargradient(x1:0, y1:0, x2:1, y2:1, stop:0 #b4b4b4, stop:1 #8f8f8f); "
                                 "border: 1px solid #5c5c5c; "
                                 "width: 18px; "
                                 "margin: 0px 0px; "
                                 "} ");

//! ------------------------
//! function: getWorkingDir
//! details:
//! ------------------------
QString GeneralDelegate::getWorkingDir() const
{
    ifstream is;
    QString settingsFileName = QString(SYSTEM_PROGRAM_DATA).append("\\WB\\settings.txt");
    cout<<"Settings file name: "<<settingsFileName.toStdString()<<endl;
    is.open(settingsFileName.toStdString());
    if(is.is_open())
    {
        std::string val;
        char tmp[512],wd[512];
        std::getline(is,val);
        is.close();
        sscanf(val.c_str(),"%s%s",tmp,wd);
        return QString::fromLatin1(wd);
    }
    return "";
}

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
GeneralDelegate::GeneralDelegate(QWidget *parent): QStyledItemDelegate(parent)
{
    //cout<<"GeneralDelegate::GeneralDelegate()->____constructor called____"<<endl;
}

//! ------------------------
//! function: constructor I
//! details:
//! ------------------------
GeneralDelegate::GeneralDelegate(const occHandle(AIS_InteractiveContext) &aCTX, QWidget *parent):
    QStyledItemDelegate(parent),myCTX(aCTX)
{
    //cout<<"GeneralDelegate::GeneralDelegate()->____constructor called I____"<<endl;
}

//! -----------------------
//! function: createEditor
//! details:
//! -----------------------
QWidget* GeneralDelegate::createEditor(QWidget *parent, const QStyleOptionViewItem &option, const QModelIndex &index) const
{
    Q_UNUSED(option)

    if(!index.isValid()) return 0;
    QVariant data = index.model()->data(index, Qt::UserRole);

    //! -----------------------------------------------------------
    //! type of data contained in the item
    //! this is used to check if the second colum has been clicked
    //! -----------------------------------------------------------
    QString typeName = QString("%1").arg(data.typeName());

    //! ---------------------------------------------
    //! create the editor according to the data type
    //! ---------------------------------------------
    if(typeName=="Property")
    {
        QString propertyName = index.model()->data(index, Qt::UserRole).value<Property>().getName();

#ifdef COSTAMP_VERSION
        //! --------------------------------------
        //! closure force && third phase pressure
        //! --------------------------------------
        if(propertyName=="Intensification pressure" || propertyName=="Force value")
        {
            QLineEdit *editor = new QLineEdit(parent);
            QDoubleValidator *validator = new QDoubleValidator();
            validator->setBottom(0.0);
            editor->setValidator(validator);
            return editor;
        }

        //! -------------------------------
        //! Time step builder - "Activate"
        //! -------------------------------
        if(propertyName =="Activate")
        {
            QPushButton *editor = new QPushButton("Go",parent);
            editor->setContentsMargins(0,0,0,0);
            connect(editor,SIGNAL(pressed()),this,SLOT(commitAndCloseTimeStepBuilderButton()));
            return editor;
        }

        //! --------------------
        //! "Time history file"
        //! --------------------
        if(propertyName =="Time history file")
        {
            QString timeHistoryFileLoc = this->getCurrentNode()->getPropertyValue<QString>("Time history file");
            QString relativeFileName = timeHistoryFileLoc.split("/").last();
            timeHistoryFileLoc.chop(relativeFileName.length()+1);

            QFileSelect *editor = new QFileSelect("",parent);
            connect(editor,SIGNAL(editingFinished()),this,SLOT(commitAndCloseTimeStepBuilderFileSelector()));

            return editor;
        }
#endif

        //! ----------------------
        //! "Jx" "Jy" "Jz" "Mass"
        //! ----------------------
        if(propertyName =="Jx" || propertyName =="Jy" || propertyName =="Jz" || propertyName =="Mass")
        {
            QLineEdit *le = new QLineEdit(parent);
            QDoubleValidator *dval = new QDoubleValidator();
            if(propertyName == "Mass") dval->setBottom(0.0);
            le->setValidator(dval);
            return le;
        }
        //! ------------------
        //! "Static/Transient"
        //! ------------------
        if(propertyName =="Static/Transient")
        {
            QComboBox *cb = new QComboBox(parent);
            data.setValue(Property::timeIntegration_steadyState);
            cb->addItem("Static",data);
            data.setValue(Property::timeIntegration_transient);
            cb->addItem("Transient",data);
            return cb;
        }
        //! ----------------
        //! "Analysis type"
        //! ----------------
        if(propertyName =="Analysis type")
        {
            QComboBox *cb = new QComboBox(parent);
            QStandardItem *curMainTreeItem = this->getCurrenItem();
            QStandardItem *curAnalysisRoot = curMainTreeItem->parent();
            SimulationNodeClass *curAnalysisNode = curAnalysisRoot->data(Qt::UserRole).value<SimulationNodeClass*>();
            switch(curAnalysisNode->getType())
            {
            case SimulationNodeClass::nodeType_structuralAnalysis:
            {
                return Q_NULLPTR;
            }
                break;
            case SimulationNodeClass::nodeType_thermalAnalysis:
            {
                return Q_NULLPTR;
            }
                break;
            case SimulationNodeClass::nodeType_combinedAnalysis:
            {
                data.setValue(Property::analysisType_structural);
                cb->addItem("Structural",data);
                data.setValue(Property::analysisType_thermal);
                cb->addItem("Thermal",data);
                data.setValue(Property::analysisType_modal);
                cb->addItem("Modal",data);
                data.setValue(Property::analysisType_frequencyResponse);
                cb->addItem("Frequency response",data);
                data.setValue(Property::analysisType_uncoupledTemperatureDisplacement);
                cb->addItem("Uncoupled temperature displacement",data);
                data.setValue(Property::analysisType_coupledTemperatureDisplacement);
                cb->addItem("Coupled temperature displacement",data);
            }
                break;
            }
            return cb;
        }
        //! -----------------------
        //! "Structural time step"
        //! -----------------------
        if(propertyName == "Structural time step")
        {
            QLineEdit *editor = new QLineEdit(parent);
            QIntValidator *validator = new QIntValidator();
            validator->setBottom(0);
            return editor;
        }
        //! ---------------------
        //! "Coupling time step"
        //! ---------------------
        if(propertyName =="Coupling time step")
        {
            QComboBox *cb = new QComboBox(parent);
            SimulationNodeClass *curNode = this->getCurrentNode();
            QString parentTimeTag = curNode->getPropertyValue<QString>("Time tag");

            //! -------------------------------
            //! retrieve the "Model" root item
            //! -------------------------------
            SimulationManager* sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
            QStandardItem *itemModelRoot = sm->getTreeItem(SimulationNodeClass::nodeType_root);

            QStandardItem *curAnalysisRoot = Q_NULLPTR;
            for(int n=0; n<itemModelRoot->rowCount(); n++)
            {
                SimulationNodeClass *node = itemModelRoot->child(n,0)->data(Qt::UserRole).value<SimulationNodeClass*>();
                if(node->isAnalysisRoot()==false) continue;
                QString timeTag = node->getPropertyValue<QString>("Time tag");
                if(timeTag!=parentTimeTag) continue;
                curAnalysisRoot = itemModelRoot->child(n,0);
            }
            if(curAnalysisRoot==Q_NULLPTR) return Q_NULLPTR;
            QStandardItem *itemAnalysisSettings = curAnalysisRoot->child(0,0);
            SimulationNodeClass *analysisSettingsNode = itemAnalysisSettings->data(Qt::UserRole).value<SimulationNodeClass*>();
            int NbSteps = analysisSettingsNode->getPropertyValue<int>("Number of steps");

            for(int n=1; n<=NbSteps; n++)
            {
                data.setValue(n);
                cb->addItem(QString("%1").arg(n),data);
            }
            return cb;
        }
        //! ----------------------------
        //! "Imported body temperature"
        //! ----------------------------
        if(propertyName =="Imported body temperature")
        {
            QComboBox *cb = new QComboBox(parent);

            //! -----------------------------
            //! retrieve the "Analysis" item
            //! -----------------------------
            SimulationNodeClass *curNode = this->getCurrentNode();

            void *p = curNode->getPropertyValue<void*>("Analysis");
            QStandardItem *item = (QStandardItem*)(p);
            SimulationNodeClass *node = item->data(Qt::UserRole).value<SimulationNodeClass*>();
            if(node->getType()==SimulationNodeClass::nodeType_NULL)
            {
                cout<<"____Thermal analysis root selected yet____"<<endl;
                return Q_NULLPTR;
            }
            cout<<"____Analysis root recorded into this item: \""<<item->data(Qt::DisplayRole).toString().toStdString()<<"\"____"<<endl;

            //! -----------------------------------------------------------------------
            //! the void pointer stored into this item =>after reload<= contains
            //! an alias of the (void*)(QStandardItem*) stored into thermal analysis
            //! so it has no child, and at line 284 the number of children is required
            //! -----------------------------------------------------------------------
            QString timeTag = node->getPropertyValue<QString>("Time tag");
            SimulationManager* sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
            QStandardItem *itemModelRoot = sm->getTreeItem(SimulationNodeClass::nodeType_root);
            //cout<<"____"<<itemModelRoot->data(Qt::DisplayRole).value<QString>().toStdString()<<"____"<<endl;
            bool found = false;
            QStandardItem *itemAnalysis;
            for(int i=0; itemModelRoot->rowCount()-1; i++)
            {
                itemAnalysis = itemModelRoot->child(i,0);
                SimulationNodeClass *nodeAnalysisRoot = itemAnalysis->data(Qt::UserRole).value<SimulationNodeClass*>();
                if(nodeAnalysisRoot->isAnalysisRoot()==false) continue;
                QString nodeAnalysisRootTimeTag = nodeAnalysisRoot->getPropertyValue<QString>("Time tag");
                if(nodeAnalysisRootTimeTag==timeTag)
                {
                    found = true;
                    break;
                }
            }
            if(!found) return Q_NULLPTR;
            item = itemAnalysis;

            //! -----------------------------------------------------------
            //! retrieve the "Solution" item of the selected analysis root
            //! -----------------------------------------------------------
            cout<<"____number of thermal results: "<<item->rowCount()<<"____"<<endl;
            QStandardItem *itemOtherAnalysisSolution = item->child(item->rowCount()-1);
            bool foundTemperatureDistribution = false;
            for(int n=1; n<itemOtherAnalysisSolution->rowCount(); n++)
            {
                //cout<<"____n: "<<n<<"____"<<endl;
                QStandardItem *curResultItem = itemOtherAnalysisSolution->child(n,0);
                SimulationNodeClass *curResultNode = curResultItem->data(Qt::UserRole).value<SimulationNodeClass*>();

                //! -------------
                //! sanity check
                //! -------------
                if(curResultNode->isAnalysisResult()==false) continue;
                if(curResultNode->getType()!=SimulationNodeClass::nodeType_solutionThermalTemperature) continue;

                void *pp = (void*)(curResultItem);
                data.setValue(pp);
                QString name = curResultNode->getName();
                cb->addItem(name,data);
                foundTemperatureDistribution = true;
            }
            if(foundTemperatureDistribution==true)
            {
                cout<<"____found result temperature distribution____"<<endl;
                return cb;
            }
            else
            {
                cout<<"____Thermal analysis selected but no imported body item selected yet____"<<endl;
                return Q_NULLPTR;
            }
        }
        //! -------------------------------------------
        //! "Analysis" - for imported body temperature
        //! -------------------------------------------
        if(propertyName =="Analysis")
        {
            QComboBox *cb = new QComboBox(parent);

            //! -------------------------------
            //! retrieve the "Model" root item
            //! -------------------------------
            SimulationManager* sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
            QStandardItem *itemModelRoot = sm->getTreeItem(SimulationNodeClass::nodeType_root);
            cout<<"____"<<itemModelRoot->data(Qt::DisplayRole).value<QString>().toStdString()<<"____"<<endl;

            //! -----------------------------------
            //! retrieve the current analysis root
            //! -----------------------------------
            QStandardItem *activeAnalysisRoot = sm->getActiveAnalysisBranch();
            QString timeTag = activeAnalysisRoot->data(Qt::UserRole).value<SimulationNodeClass*>()->getPropertyValue<QString>("Time tag");
            cout<<"____Active analysis branch: "<<activeAnalysisRoot->data(Qt::DisplayRole).value<QString>().toStdString()<<"____"<<endl;
            cout<<"____Time tag: "<<timeTag.toStdString()<<"____"<<endl;

            //! --------------------
            //! scan the model root
            //! --------------------
            bool found = false;
            for(int n=0; n<itemModelRoot->rowCount(); n++)
            {
                QStandardItem *curAnalysisRootItem = itemModelRoot->child(n,0);
                SimulationNodeClass *curAnalysisRootNode = curAnalysisRootItem->data(Qt::UserRole).value<SimulationNodeClass*>();
                if(curAnalysisRootNode->isAnalysisRoot()==false)
                {
                    //cout<<"____not a simulation root____"<<endl;
                    continue;
                }

                //! ------------------------------------------------
                //! jump over what is not a "Thermal analysis" root
                //! ------------------------------------------------
                if(curAnalysisRootNode->getType()!=SimulationNodeClass::nodeType_thermalAnalysis)
                {
                    cout<<"____\""<<curAnalysisRootNode->getName().toStdString()<<"\" not a thermal analysis branch____"<<endl;
                    continue;
                }

                cout<<"____found a thermal analysis root____"<<endl;

                QString curAnalysisRootNodeTimeTag = curAnalysisRootNode->getPropertyValue<QString>("Time tag");
                if(timeTag == curAnalysisRootNodeTimeTag)
                {
                    cout<<"____jumping over the same root____"<<endl;
                    continue;
                }

                //! ----------------------------
                //! add an item to the combobox
                //! ----------------------------
                cout<<"____add item to the combobox____"<<endl;
                QString name = curAnalysisRootItem->data(Qt::DisplayRole).toString();
                void *p = (void*)(curAnalysisRootItem);
                data.setValue(p);
                cb->addItem(name,data);
                found = true;
            }
            if(found==true)
            {
                cout<<"____returning an editor____"<<endl;
                return cb;
            }
            return Q_NULLPTR;
        }
        //! ---------------
        //! "Fatigue algo"
        //! ---------------
        if(propertyName =="Fatigue algo")
        {
            QComboBox *cb = new QComboBox(parent);
            cb->addItem("Basquin Coffin Manson",0);
            cb->addItem("Equivalent strain range",1);
            return cb;
        }
        //! --------------------
        //! Material assignment
        //! --------------------
        if(propertyName =="Assignment")
        {
            //! ---------------------------
            //! material names definitions
            //! ---------------------------
            QList<QString> matNames;
            matNames<<"Structural steel"<<
                      "Bilinear steel"<<
                      "H11 fatigue"<<
                      "F22 fatigue"<<
                      "B16_fatigue"<<
                      "F6NM_fatigue"<<
                      "F92_fatigue"<<
                      "A479_fatigue"<<
                      "SA479_XM19_fatigue"<<
                      "SA182-B8M_CL2"<<
                      "SA182-F316"<<
                      "SA352-LCB";

            QComboBox *cb = new QComboBox(parent);
            for(int index = 0; index<matNames.length(); index++)
            cb->addItem(matNames[index],index);
            return cb;
        }
        //! ----------------------------------------------
        //! "Metric type" - currently implemented for tet
        //! 0 => "Modified Mean Ratio (MMR)"
        //! 1 => "Modified Condition Number (MCN)"
        //! 2 => "Modified volume-length (MIVL)"
        //! ----------------------------------------------
        if(propertyName =="Metric type")
        {
            QComboBox *cb = new QComboBox(parent);
            cb->addItem("Modified Mean Ratio",0);
            cb->addItem("Modified Condition Number",1);
            cb->addItem("Modified volume-length",2);
            cb->addItem("None",3);
            return cb;
        }
        //! -----------
        //! "Grouping"
        //! -----------
        if(propertyName =="Grouping")
        {
            QComboBox *cb = new QComboBox(parent);
            cb->addItem("By bodies",0);
            cb->addItem("By master faces",1);
            cb->addItem("By slave faces",2);
            cb->addItem("By master body",3);
            cb->addItem("By slave body",4);
            cb->addItem("Ungrouped",5);
            return cb;
        }
        //! ---------------
        //! "Element list"
        //! ---------------
        if(propertyName =="Element list")
        {
            QLineEdit *editor = new QLineEdit(parent);
            //QRegExpValidator *validator = new QRegExpValidator();
            //validator->setRegExp(QRegExp("\\d,?+"));
            //editor->setValidator(validator);
            return editor;
        }
        //! ------------------
        //! "Selection method"
        //! ------------------
        if(propertyName =="Selection method")
        {
            QComboBox *editor = new QComboBox(parent);
            editor->addItem("Elements numbers list",0);
            editor->addItem("Picking",1);
            return editor;
        }
        //! -----------------------------------
        //! "Activation status" (model change)
        //! -----------------------------------
        if(propertyName =="Activation status")
        {
            QComboBox *editor = new QComboBox(parent);
            editor->addItem("Remove",0);
            editor->addItem("Inactive",1);
            editor->addItem("Add",2);
            return editor;
        }
        //! ---------------------------
        //! "Item type" (model change)
        //! ---------------------------
        if(propertyName =="Item type")
        {
             QComboBox *editor = new QComboBox(parent);
             editor->addItem("Bodies",0);
             editor->addItem("Contact pairs",1);
             return editor;
        }
        //! -------------------
        //! "Source file path"
        //! -------------------
        if(propertyName =="Source file path")
        {
            QFileSelect *editor = new QFileSelect(QDir::currentPath(),parent);
            return editor;
        }
        //! ----------------------------------------------
        //! "Solid bodies" "Surface bodies" "Line bodies"
        //! ----------------------------------------------
        if(propertyName =="Solid bodies" || propertyName == "Surface bodies" || propertyName == "Line bodies")
        {
            QComboBox *editor = new QComboBox(parent);
            data.setValue(false);
            editor->addItem("No",data);
            data.setValue(true);
            editor->addItem("Yes",data);
            return editor;
        }
        //! ------------------
        //! "Tolerance value"
        //! ------------------
        if(propertyName=="Tolerance value")
        {
            QLineEdit *editor = new QLineEdit(parent);
            QDoubleValidator *validator = new QDoubleValidator();
            validator->setBottom(1e-6);
            validator->setTop(100);
            return editor;
        }
        //! --------------------
        //! "Angular criterion"
        //! --------------------
        if(propertyName =="Angular criterion")
        {
            QLineEdit *editor = new QLineEdit(parent);
            QDoubleValidator *validator = new QDoubleValidator();
            validator->setBottom(1.0);
            validator->setTop(89.0);
            return editor;
        }
        //! -------------------
        //! "Geometry healing"
        //! -------------------
        if(propertyName =="Geometry healing")
        {
            QComboBox *editor = new QComboBox(parent);
            data.setValue(false);
            editor->addItem("No",data);
            data.setValue(true);
            editor->addItem("Yes",data);
            return editor;
        }
        //! -----------------------------------------
        //! "Defeaturing" "Healing" "Simplification"
        //! -----------------------------------------
        if(propertyName=="Defeaturing" || propertyName == "Healing" || propertyName == "Simplification" ||
                propertyName =="Preserve boundary conditions edges" || propertyName == "Project points on geometry")
        {
            QComboBox *editor = new QComboBox(parent);
            data.setValue(false);
            editor->addItem("Off",data);
            data.setValue(true);
            editor->addItem("On",data);
            return editor;
        }
        //! --------------------------------
        //! "Min face size" "Max face size"
        //! --------------------------------
        if(propertyName =="Min face size" || propertyName=="Max face size")
        {
            QLineEdit *editor = new QLineEdit(parent);
            QDoubleValidator *validator = new QDoubleValidator();
            validator->setBottom(1e-6);
            return editor;
        }
        //! -----------------------------------------
        //! "Angular deflection" "Linear deflection"
        //! -----------------------------------------
        if(propertyName == "Angular deflection" || propertyName == "Linear deflection")
        {
            QLineEdit *editor = new QLineEdit(parent);
            QDoubleValidator *validator = new QDoubleValidator();
            if(propertyName=="Angular deflection")
            {
                const double minAngDeflection = 0.001;      //! 0.1 deg
                const double maxAngDeflection = 0.7854;     //! 45 deg
                validator->setBottom(minAngDeflection);
                validator->setTop(maxAngDeflection);
            }
            if(propertyName=="Linear deflection")
            {
                const double minLinearDeflection = 0.001;    //! check value
                const double maxLinearDeflection = 100;      //! check value
                validator->setBottom(minLinearDeflection);
                validator->setTop(maxLinearDeflection);
            }
            editor->setValidator(validator);
            connect(editor,SIGNAL(editingFinished()),this,SLOT(commitAndCloseLineEdit()));
            return editor;
        }
        //! --------------
        //! "Tessellator"
        //! --------------
        if(propertyName =="Tessellator")
        {
            QVariant data;
            QComboBox *editor = new QComboBox(parent);
            data.setValue(Property::meshEngine2D_OCC_STL); editor->addItem("Standard STL",data);
            data.setValue(Property::meshEngine2D_OCC_ExpressMesh); editor->addItem("Express mesh",data);
            return editor;
        }
        //! -----------
        //! "Patch conforming"
        //! -----------
        if(propertyName =="Patch conforming")
        {
            QComboBox *editor = new QComboBox(parent);
            data.setValue(false);
            editor->addItem("Off",data);
            data.setValue(true);
            editor->addItem("On",data);
            return editor;
        }
        //! -----------------
        //! "Surface mesher"
        //! -----------------
        if(propertyName == "Surface mesher")
        {
            QComboBox *editor = new QComboBox(parent);
            if(this->getCurrentNode()->getPropertyValue<bool>("Patch conforming")==true)
            {
                data.setValue(Property::meshEngine2D_Netgen);
                editor->addItem("Netgen",data);
                data.setValue(Property::meshEngine2D_OCC_ExpressMesh);
                editor->addItem("Express mesh",data);
            }
            else
            {
                data.setValue(Property::meshEngine2D_Netgen_STL);
                editor->addItem("Netgen STL",data);
                //data.setValue(Property::meshEngine2D_TetgenBR);
                //editor->addItem("Tetgen BR",data);
            }
            connect(editor,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseComboBox()));
            return editor;
        }
        //! --------------
        //! "Volume mesh"
        //! --------------
        if(propertyName == "Volume mesher")
        {
            QComboBox *editor = new QComboBox(parent);
            if(this->getCurrentNode()->getPropertyValue<bool>("Patch conforming")==true)
            {
                data.setValue(Property::meshEngine3D_Netgen);   // index "0"
                editor->addItem("Netgen",data);
                data.setValue(Property::meshEngine3D_Tetgen);   // index "1"
                editor->addItem("Tetgen",data);
            }
            else    // "Patch conforming" = false
            {
                data.setValue(Property::meshEngine3D_Netgen_STL);   // index "0"
                editor->addItem("Netgen STL",data);
                data.setValue(Property::meshEngine3D_Tetgen);       // index "1"
                editor->addItem("Tetgen",data);
                data.setValue(Property::meshEngine3D_Tetgen_BR);    // index "2"
                editor->addItem("Tetgen BR",data);
                //! ---------------
                //! blindo TetWild
                //! ---------------
                data.setValue(Property::meshEngine3D_TetWild);      // index "3"
                editor->addItem("Experimental mesher",data);
            }
            return editor;
        }
        //! ----------------------------------------
        //! "Envelope sizing" "Ideal length sizing"
        //! "0" => "Relative"
        //! "1" => "Absolute"
        //! ----------------------------------------
        if(propertyName == "Envelope sizing" || propertyName =="Ideal length sizing")
        {
            QComboBox *editor = new QComboBox(parent);
            editor->addItem("Relative",int(0));
            editor->addItem("Absolute",int(1));
            return editor;
        }
        //! --------------------------------------------------------------------
        //! "Relative size" "Absolute Size" "Relative length" "Absolute length"
        //! --------------------------------------------------------------------
        if(propertyName =="Relative size" || propertyName =="Relative length"
                || propertyName =="Absolute size" || propertyName =="Absolute length")
        {
            QLineEdit *editor = new QLineEdit(parent);
            QDoubleValidator *validator = new QDoubleValidator();
            validator->setBottom(0.0001);
            editor->setValidator(validator);
            return editor;
        }
        //! ----------------------------------------
        //! stress/strain source - for fatigue tool
        //! ----------------------------------------
        if(propertyName =="Stress/strain source")
        {
            QComboBox *editor = new QComboBox(parent);
            editor->addItem("Equivalent strain",0);
            editor->addItem("Equivalent plastic strain",1);
            return editor;
        }
        //! --------
        //! Mapping
        //! --------
        if(propertyName=="Mapping")
        {
            QComboBox *editor = new QComboBox(parent);
            editor->addItem("Linear",0);
            editor->addItem("Logarithmic",1);
            return editor;
        }
        //! -------------------------------
        //! "Component" - for fatigue tool
        //! -------------------------------
        if(propertyName =="Component")
        {
            if(this->getCurrentNode()->getPropertyValue<int>("Fatigue algo")==1) return Q_NULLPTR;

            QComboBox *editor = new QComboBox(parent);
            editor->addItem("Total",0);
            editor->addItem("Mechanical",1);
            editor->addItem("Thermal",2);
            return editor;
        }
        //! --------------------------------------
        //! "Number of cycles" - for fatigue tool
        //! --------------------------------------
        if(propertyName =="Number of cycles")
        {
            QLineEdit *editor = new QLineEdit(parent);
            QIntValidator *validator = new QIntValidator();
            validator->setBottom(10);
            editor->setValidator(validator);
            return editor;
        }
        //! ------------------------------------------
        //! triaxiality correction - for fatigue tool
        //! ------------------------------------------
        if(propertyName =="Triaxiality correction")
        {
            QComboBox *editor = new QComboBox(parent);
            editor->addItem("Off",0);
            editor->addItem("On",1);
            return editor;
        }
        //! -----------------------------------------------
        //! Basquin Coffin Mason fatigue model coefficients
        //! -----------------------------------------------
        if(propertyName =="Rupture stress" || propertyName =="Rupture strain" ||
                propertyName=="n'" || propertyName =="b")
        {
            QLineEdit *editor = new QLineEdit(parent);
            QDoubleValidator *validator = new QDoubleValidator();
            validator->setBottom(0.0);
            editor->setValidator(validator);
            return editor;
        }
        //! ---------------------
        //! "Volume mesh engine"
        //! ---------------------
        if(propertyName =="Volume mesh engine")
        {
            QComboBox *editor = new QComboBox(parent);
            editor->addItem("Netgen");
            editor->addItem("Tetgen");
            return editor;
        }

        //! ----------------
        //! "Lock boundary"
        //! ----------------
        if(propertyName =="Lock boundary")
        {
            QComboBox *editor = new QComboBox(parent);
            editor->addItem("Locked",true);
            editor->addItem("Free",false);
            return editor;
        }
        //! --------------------------------------------------------------
        //! "Guiding vectors smoothing steps" "Thickness smoothing steps"
        //! --------------------------------------------------------------
        if(propertyName=="Guiding vectors smoothing steps" || propertyName=="Thickness smoothing steps")
        {
            QLineEdit *editor = new QLineEdit(parent);
            QIntValidator *validator = new QIntValidator();
            validator->setBottom(1);
            editor->setValidator(validator);
            return editor;
        }
        //! ------------------------
        //! "Curvature sensitivity"
        //!  -----------------------
        if(propertyName =="Curvature sensitivity" ||
                propertyName =="Guiding vector smoothing - curvature sensitivity" ||
                propertyName =="Thickness smoothing - curvature sensitivity")
        {
            QLineEdit *editor = new QLineEdit(parent);
            QDoubleValidator *validator = new QDoubleValidator();
            validator->setBottom(1e-5);
            editor->setValidator(validator);
            return editor;
        }
        //! -------------
        //! "Transition"
        //! -------------
        if(propertyName =="Transition")
        {
            QLineEdit *editor = new QLineEdit(parent);
            QDoubleValidator *validator = new QDoubleValidator();
            validator->setBottom(0.1);
            validator->setTop(100.0);
            editor->setValidator(validator);
            return editor;
        }
        //! --------------------------------------------
        //! "Options" (for prismatic layers generation)
        //! 0 => "First layer thickness"
        //! 1 => "Total thickness"
        //! --------------------------------------------
        if(propertyName =="Options")
        {
            QVariant data;
            QComboBox *editor = new QComboBox(parent);
            data.setValue(0); editor->addItem("First layer thickness",data);    //! data=0, index "0"
            data.setValue(1); editor->addItem("Total thickness",data);          //! data=1, index "1"
            return editor;
        }
        //! -------------------
        //! "Number of layers"
        //! -------------------
        if(propertyName=="Number of layers")
        {
            QLineEdit *editor = new QLineEdit(parent);
            QIntValidator *validator = new QIntValidator();
            validator->setBottom(1);
            editor->setValidator(validator);
            return editor;
        }
        //! -----------------------------------------------------
        //! "Expansion ratio" "Total thickness" "First layer height"
        //! -----------------------------------------------------
        if(propertyName =="Expansion ratio" || propertyName =="Total thickness" || propertyName =="First layer height")
        {
            QLineEdit *editor = new QLineEdit(parent);
            QDoubleValidator *validator = new QDoubleValidator();
            if(propertyName =="Expansion ratio") validator->setBottom(1.0);
            if(propertyName =="Total thickness") validator->setBottom(0.0);
            if(propertyName =="First layer height") validator->setBottom(0.0);
            editor->setValidator(validator);
            return editor;
        }
        //! --------------------------------------------------------
        //! "Check self intersections" "Check mutual intersections"
        //! --------------------------------------------------------
        if(propertyName=="Check mutual intersections" || propertyName=="Check self intersections")
        {
            QVariant data;
            QComboBox *editor = new QComboBox(parent);
            data.setValue(false); editor->addItem("Off",data);          //! data=0, index "0"
            data.setValue(true); editor->addItem("On",data);            //! data=1, index "1"
            return editor;
        }
        //! ----------------
        //! "Run in memory"
        //! ----------------
        if(propertyName =="Run in memory")
        {
            QVariant data;
            QComboBox *editor = new QComboBox(parent);
            data.setValue(true); editor->addItem("Active",data);    //! index "0"
            data.setValue(false); editor->addItem("Off",data);      //! index "1"
            return editor;
        }
        //! ---------------------------
        //! "Surface mesh size" "Pair"
        //! ---------------------------
        if(propertyName =="Pair distance")
        {
            QLineEdit *editor = new QLineEdit(parent);
            QDoubleValidator *doubleValidator = new QDoubleValidator();
            if(propertyName=="Surface mesh size")
            {
                doubleValidator->setBottom(0.01);
                doubleValidator->setTop(100.0);
            }
            editor->setValidator(doubleValidator);
            return editor;
        }
        //! --------
        //! "Level"
        //! --------
        if(propertyName =="Level")
        {
            QSlider *editor = new QSlider(parent);
            editor->setOrientation(Qt::Horizontal);
            editor->setStyleSheet(sliderStyleSheet);
            editor->setRange(0,100);

            //connect(editor,SIGNAL(sliderPressed()),this,SLOT(displayLevelSliderValue(int)));
            connect(editor,SIGNAL(sliderMoved(int)),this,SLOT(displayLevelSliderValue(int)));
            return editor;
        }
        //! ----------------
        //! New mesh method
        //! ----------------
        if(propertyName =="Method")
        {
            QVariant data;
            QComboBox *editor = new QComboBox(parent);
            data.setValue(0); editor->addItem("Automatic",data);
            data.setValue(1); editor->addItem("Netgen Tetgen",data);
            data.setValue(2); editor->addItem("Full Netgen",data);
            data.setValue(3); editor->addItem("Emesh Tetgen",data);
            data.setValue(4); editor->addItem("Emesh Netgen",data);
            return editor;
        }
        //! ------------
        //! "Min" "Max"
        //! ------------
        if(propertyName =="Min" || propertyName =="Max")
        {
            QLineEdit *editor = new QLineEdit(parent);
            QDoubleValidator *doubleValidator = new QDoubleValidator();
            editor->setValidator(doubleValidator);
            return editor;
        }
        //! --------------
        //! "# intervals"
        //! --------------
        if(propertyName =="# intervals")
        {
            QSpinBox *editor = new QSpinBox(parent);
            editor->setMinimum(3);
            editor->setMaximum(32);
            return editor;
        }
        //! -------------
        //! "Scale type"
        //! -------------
        else if(propertyName =="Scale type")
        {
            QComboBox *editor = new QComboBox(parent);
            QVariant data;
            data.setValue(int(0));
            editor->addItem("Autoscale",data);
            data.setValue(int(1));
            editor->addItem("Custom",data);
            return editor;
        }
        //! ----------------
        //! "Remote points"
        //! ----------------
        else if(propertyName =="Remote points")
        {
            //! -----------------------------------------------------------------
            //! the editor (Combo box) for the remote point lists can be created
            //! only if the remote point root exixts
            //! -----------------------------------------------------------------
            QComboBox *editor = new QComboBox(parent);
            SimulationManager* sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
            QStandardItem *itemRemotePointRoot = sm->getTreeItem(SimulationNodeClass::nodeType_remotePointRoot);
            if(itemRemotePointRoot!=NULL)
            {
                for(int k=0; k<itemRemotePointRoot->rowCount(); k++)
                {
                    QStandardItem *itemCurrentRemotePoint = itemRemotePointRoot->child(k,0);
                    QString name = itemCurrentRemotePoint->data(Qt::DisplayRole).toString();
                    void *currRemotePointPointer = (void*)itemCurrentRemotePoint;
                    data.setValue(currRemotePointPointer);
                    editor->addItem(name,data);
                }
                return editor;
            }
            return 0;
        }
        //! ----------------------------------------------------------
        //! update frequency for writing results: "--Recurrence rate"
        //! ----------------------------------------------------------
        else if(propertyName =="--Recurrence rate")
        {
            QLineEdit *editor = new QLineEdit(parent);
            QIntValidator *validator = new QIntValidator();
            validator->setBottom(1);
            editor->setValidator(validator);
            return editor;
        }
        //! -----------------------------------
        //! Output settings "Store results at"
        //! -----------------------------------
        else if(propertyName =="Store results at")
        {
            QVariant data;
            QComboBox *editor = new QComboBox(parent);
            data.setValue(int(0)); editor->addItem("All time points",data);
            data.setValue(int(1)); editor->addItem("Last time point",data);
            data.setValue(int(2)); editor->addItem("Specified recurrence rate",data);
            return editor;
        }
        //! ------------------------------------------------------------------
        //! Ouput settings "Stress" "Strain" "Reaction forces" "Contact data"
        //! ------------------------------------------------------------------
        else if(propertyName =="Stress" || propertyName=="Strain" || propertyName =="Reaction forces" || propertyName=="Contact data")
        {
            QVariant data;
            QComboBox *editor = new QComboBox(parent);
            data.setValue(int(0)); editor->addItem("No",data);
            data.setValue(int(1)); editor->addItem("Yes",data);
            return editor;
        }
        //! --------------------------------------
        //! Advanced properties for contact group
        //! --------------------------------------
        else if(propertyName =="K" || propertyName=="Sigma infty" /*|| propertyName =="C0"*/ ||
                propertyName =="Lambda" || propertyName =="P0")
        {
            SimulationNodeClass *curNode = this->getCurrentNode();
            if(curNode->getType()==SimulationNodeClass::nodeType_connectionPair)
            {
                Property::contactBehavior behavior = this->getCurrentNode()->getPropertyItem("Behavior")->data(Qt::UserRole).value<Property>().getData().value<Property::contactBehavior>();
                Property::overpressureFunction overpressure = this->getCurrentNode()->getPropertyItem("Overpressure")->data(Qt::UserRole).value<Property>().getData().value<Property::overpressureFunction>();

                switch(behavior)
                {
                case Property::contactBehavior_asymmetric:
                {
                    switch(overpressure)
                    {
                    case Property::overpressureFunction_linear:
                    {
                        //! ------------------------------------------------
                        //! "C0" "K" "Lambda" "Sigma infty" can be modified
                        //! ------------------------------------------------
                        if(propertyName=="C0" || propertyName =="Lambda" || propertyName =="Sigma infty" || propertyName =="K")
                        {
                            QLineEdit *editor = new QLineEdit(parent);
                            QDoubleValidator *doubleValidator = new QDoubleValidator();
                            doubleValidator->setBottom(0);
                            editor->setValidator(doubleValidator);
                            return editor;
                        }
                        else return 0;
                    }
                        break;

                    case Property::overpressureFunction_exponential:
                    {
                        //! ---------------------------------
                        //! "P0" "C0" "Lambda" can be edited
                        //! ---------------------------------
                        if(propertyName =="C0" || propertyName =="Lambda" || propertyName =="P0")
                        {
                            QLineEdit *editor = new QLineEdit(parent);
                            QDoubleValidator *doubleValidator = new QDoubleValidator();
                            doubleValidator->setBottom(0);
                            editor->setValidator(doubleValidator);
                            return editor;
                        }
                        else return 0;
                    }
                        break;
                    }
                }
                    break;

                case Property::contactBehavior_symmetric:
                {
                    switch(overpressure)
                    {
                    case Property::overpressureFunction_linear:
                    {
                        //! -----------------------------------------------------------------
                        //! "K" "Lambda" "can be modified. For a symmetric face to face
                        //! contact "C0" is set to zero by CCX, and "Sigma infty" is not
                        //! defined (the overpressure function is truly bilinear)
                        //! From CCX routine "springfc_f2f.f" appears that the overpressure
                        //! equation is the same for the two cases. For a "Tied" overpressure
                        //! do not return the editor for "Lambda" (since "irrelevant"
                        //! according to the guide)
                        //! ------------------------------------------------------------------
                        if(propertyName =="Lambda" || propertyName =="K")
                        {
                            QLineEdit *editor = new QLineEdit(parent);
                            QDoubleValidator *doubleValidator = new QDoubleValidator();
                            doubleValidator->setBottom(0);
                            editor->setValidator(doubleValidator);
                            return editor;
                        }
                        else return 0;
                    }
                        break;

                    case Property::overpressureFunction_tied:
                    {
                        //! -------------------------
                        //! see the previous comment
                        //! -------------------------
                        if(propertyName =="Lambda" || propertyName =="K")
                        {
                            QLineEdit *editor = new QLineEdit(parent);
                            QDoubleValidator *doubleValidator = new QDoubleValidator();
                            doubleValidator->setBottom(0);
                            editor->setValidator(doubleValidator);
                            return editor;
                        }
                        else return 0;
                    }
                        break;

                    case Property::overpressureFunction_exponential:
                    {
                        //! ---------------------------------
                        //! "P0 "C0" "Lambda" can be edited
                        //! ---------------------------------
                        if(propertyName =="P0" || propertyName =="C0" || propertyName =="Lambda")
                        {
                            QLineEdit *editor = new QLineEdit(parent);
                            QDoubleValidator *doubleValidator = new QDoubleValidator();
                            doubleValidator->setBottom(0);
                            editor->setValidator(doubleValidator);
                            return editor;
                        }
                        else return 0;
                    }
                        break;
                    }
                }
                    break;
                }
            }

            //! ---------------------------------------------
            //! "K" and "Sigma infty" are also properties of
            //! the "Compression only support"
            //! ---------------------------------------------
            if(curNode->getType()==SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_CompressionOnlySupport)
            {
                QLineEdit *editor = new QLineEdit(parent);
                QDoubleValidator *doubleValidator = new QDoubleValidator();
                doubleValidator->setBottom(0);
                editor->setValidator(doubleValidator);
                return editor;
            }
        }
        //! -----
        //! "C0"
        //! -----
        else if(propertyName == "C0")
        {
            QLineEdit *editor = new QLineEdit(parent);
            QDoubleValidator *doubleValidator = new QDoubleValidator();
            doubleValidator->setBottom(0);
            editor->setValidator(doubleValidator);
            return editor;
        }
        //! ----------------
        //! "Small sliding"
        //! ----------------
        else if(propertyName=="Small sliding")
        {
            SimulationNodeClass *curNode = this->getCurrentNode();
            Property::contactBehavior val = curNode->getPropertyItem("Behavior")->data(Qt::UserRole).value<Property>().getData().value<Property::contactBehavior>();
            if(val == Property::contactBehavior_asymmetric)
            {
                QComboBox *editor = new QComboBox(parent);
                QVariant data;
                data.setValue(0);
                editor->addItem("Off",data);
                data.setValue(1);
                editor->addItem("On",data);
                return editor;
            }
            else return 0;
        }
        //! --------------------------------------------------------
        //! Line search "0" => "Program controlled" "1" => "Custom"
        //! --------------------------------------------------------
        else if(propertyName =="Line search")
        {
            QComboBox *editor = new QComboBox(parent);
            editor->addItem("Program controlled"); data.setValue(0); editor->setItemData(0,data);
            editor->addItem("Custom"); data.setValue(1); editor->setItemData(1,data);
            return editor;
        }
        //! -------------------------
        //! "Line search parameters"
        //! -------------------------
        else if(propertyName == "Min value" || propertyName == "Max value")
        {
            int lineSearch = this->getCurrentNode()->getPropertyItem("Line search")->data(Qt::UserRole).value<Property>().getData().toInt();
            double curMin = this->getCurrentNode()->getPropertyItem("Min value")->data(Qt::UserRole).value<Property>().getData().toDouble();
            double curMax = this->getCurrentNode()->getPropertyItem("Max value")->data(Qt::UserRole).value<Property>().getData().toDouble();
            if(lineSearch!=0)
            {
                QLineEdit *editor = new QLineEdit(parent);
                QDoubleValidator *doubleValidator = new QDoubleValidator();
                if(propertyName == "Min value")
                {
                    doubleValidator->setBottom(0.25);
                    doubleValidator->setTop(curMax);
                }
                if(propertyName == "Max value")
                {
                    doubleValidator->setTop(1.01);
                    doubleValidator->setBottom(curMin);
                }
                editor->setValidator(doubleValidator);
                return editor;
            }
            else return 0;
        }
        //! ------------------
        //! "Cutback factors"
        //! ------------------
        else if(propertyName =="D_f" || propertyName =="D_C" || propertyName =="D_B" || propertyName =="D_A"
                || propertyName =="D_S" || propertyName =="D_H" || propertyName =="D_D" || propertyName =="W_G")
        {
            int cutBackType = this->getCurrentNode()->getPropertyItem("Cutback factors")->data(Qt::UserRole).value<Property>().getData().toInt();
            if(cutBackType==1)
            {
                if(propertyName !="D_S" && propertyName!="D_H" && propertyName!="W_G")
                {
                    QLineEdit *editor = new QLineEdit(parent);
                    QDoubleValidator *doubleValidator = new QDoubleValidator();
                    doubleValidator->setBottom(0.0);
                    editor->setValidator(doubleValidator);
                    return editor;
                }
                else return 0;
            }
            else return 0;
        }
        //! ----------------------------------------------
        //! "Flux convergence" and "Solution convergence"
        //! ----------------------------------------------
        else if(propertyName == "Flux convergence" || propertyName =="Solution convergence")
        {
            QComboBox *editor = new QComboBox(parent);
            data.setValue(0); editor->addItem("Remove"); editor->setItemData(0,data);
            data.setValue(1); editor->addItem("On"); editor->setItemData(1,data);
            data.setValue(2); editor->addItem("Program controlled"); editor->setItemData(2,data);
            return editor;
        }
        //! ---------------------------------------
        //! "Time incrementation", "Cutback" types
        //! ---------------------------------------
        else if(propertyName == "Time incrementation" || propertyName=="Cutback factors")
        {
            QComboBox *editor = new QComboBox(parent);
            editor->addItem("Program controlled"); data.setValue(0); editor->setItemData(0,data);
            editor->addItem("Custom"); data.setValue(1); editor->setItemData(1,data);
            return editor;
        }
        //! ---------------------------------
        //! "Time incrementation" parameters
        //! ---------------------------------
        else if(propertyName =="I_0" || propertyName =="I_R" || propertyName =="I_P" || propertyName =="I_C" || propertyName =="I_L"
              || propertyName =="I_G" || propertyName =="I_S" || propertyName =="I_A" || propertyName =="I_J" || propertyName =="I_T")
        {
            int timeIncrementationType = this->getCurrentNode()->getPropertyItem("Time incrementation")->data(Qt::UserRole).value<Property>().getData().toInt();
            if(timeIncrementationType==1)
            {
                QLineEdit *editor = new QLineEdit(parent);
                QIntValidator *intValidator = new QIntValidator();
                intValidator->setBottom(1);
                editor->setValidator(intValidator);
                if(propertyName =="I_S" || propertyName =="I_J" || propertyName=="I_T") editor = 0;
                return editor;
            }
            else return 0;
        }
        //! --------------------------------------------------
        //! --R_alpha_n --R_alpha_P R_alpha_l --epsilon_alpha
        //! --------------------------------------------------
        else if(propertyName =="--R_alpha_n" || propertyName =="--R_alpha_P"  || propertyName == "--R_alpha_l" || propertyName=="--epsilon_alpha")
        {
            int fluxConvergenceFlag = this->getCurrentNode()->getPropertyItem("Flux convergence")->data(Qt::UserRole).value<Property>().getData().toInt();
            if(fluxConvergenceFlag==1)
            {
                cout<<"____returning editor____"<<endl;
                QLineEdit *editor = new QLineEdit(parent);
                QDoubleValidator *doubleValidator = new QDoubleValidator();
                editor->setValidator(doubleValidator);
                return editor;
            }
            else return 0;
        }
        //! --------
        //! --Value
        //! --------
        else if(propertyName=="--Value")
        {
            int fluxConvergenceFlag = this->getCurrentNode()->getPropertyItem("Flux convergence")->data(Qt::UserRole).value<Property>().getData().toInt();
            if(fluxConvergenceFlag==1)
            {
                cout<<"____returning editor____"<<endl;
                QLineEdit *editor = new QLineEdit(parent);
                QDoubleValidator *doubleValidator = new QDoubleValidator();
                editor->setValidator(doubleValidator);
                return editor;
            }
            else return 0;
        }
        //! ------------
        //! --q_alpha_0
        //! ------------
        else if(propertyName=="--q_alpha_0")
        {
            int fluxConvergenceFlag = this->getCurrentNode()->getPropertyItem("Flux convergence")->data(Qt::UserRole).value<Property>().getData().toInt();
            if(fluxConvergenceFlag==1)
            {
                cout<<"____returning editor____"<<endl;
                QLineEdit *editor = new QLineEdit(parent);
                QDoubleValidator *doubleValidator = new QDoubleValidator();
                editor->setValidator(doubleValidator);
                return editor;
            }
            else return 0;
        }
        //! ------------------------------
        //! --C_alpha_n --C_alpha_epsilon
        //! ------------------------------
        else if(propertyName == "--C_alpha_n" || propertyName == "--C_alpha_epsilon")
        {
            int solutionConvergenceFlag = this->getCurrentNode()->getPropertyItem("Solution convergence")->data(Qt::UserRole).value<Property>().getData().toInt();
            if(solutionConvergenceFlag==1)
            {
                QLineEdit *editor = new QLineEdit(parent);
                QDoubleValidator *doubleValidator = new QDoubleValidator();
                editor->setValidator(doubleValidator);
                return editor;
            }
            else return 0;
        }
        //! ---------------
        //! solution items
        //! ---------------
        else if (propertyName =="Large deflection")
        {
            QComboBox *editor = new QComboBox(parent);
            data.setValue(0); editor->addItem("Off",data);
            data.setValue(1); editor->addItem("On",data);
            return editor;
        }
        //! ----------------
        //! update interval
        //! ----------------
        else if(propertyName =="Update interval")
        {
            QLineEdit *editor = new QLineEdit(parent);
            QDoubleValidator *doubleValidator = new QDoubleValidator();
            doubleValidator->setBottom(0.0);
            editor->setValidator(doubleValidator);
            return editor;
        }
        //! ---------------
        //! "Display time"
        //! ---------------
        else if(propertyName =="Display time")
        {
            QLineEdit *editor = new QLineEdit(parent);
            QDoubleValidator *doubleValidator = new QDoubleValidator();
            doubleValidator->setBottom(0.0);
            editor->setValidator(doubleValidator);
            return editor;
        }
        //! ---------------
        //! "Mode number"
        //! ---------------
        else if(propertyName =="Mode number")
        {
            QLineEdit *editor = new QLineEdit(parent);
            QIntValidator *intValidator = new QIntValidator();
            intValidator->setBottom(0);
            editor->setValidator(intValidator);
            return editor;
        }
        //! ------------
        //! "Tolerance"
        //! ------------
        else if(propertyName =="Tolerance")
        {
            SimulationNodeClass::nodeType nodeType = this->getCurrentNode()->getType();
            if(nodeType==SimulationNodeClass::nodeType_import)
            {
                QComboBox *editor = new QComboBox(parent);
                editor->addItem("User defined",0);
                editor->addItem("Normal",1);
                editor->addItem("Loose",2);
                return editor;
            }
            else
            {
                QLineEdit *editor = new QLineEdit(parent);
                QDoubleValidator *doubleValidator = new QDoubleValidator();
                doubleValidator->setBottom(0.0);
                editor->setValidator(doubleValidator);
                return editor;
            }
        }
        //! -------------
        //! "Set number"
        //! -------------
        else if(propertyName=="Set number")
        {
            QLineEdit *editor = new QLineEdit(parent);
            QIntValidator *intValidator = new QIntValidator();
            intValidator->setBottom(1);
            editor->setValidator(intValidator);
            return editor;
        }
        //! -----
        //! "By"
        //! -----
        else if(propertyName=="By")
        {
            QComboBox *editor = new QComboBox(parent);
            QVariant data;
            switch(this->getCurrentNode()->getType())
            {
            case SimulationNodeClass::nodeType_solutionThermalFlux:
            case SimulationNodeClass::nodeType_solutionThermalTemperature:
            case SimulationNodeClass::nodeType_solutionStructuralMechanicalStrain:
            case SimulationNodeClass::nodeType_solutionStructuralNodalDisplacement:
            case SimulationNodeClass::nodeType_solutionStructuralStress:
            case SimulationNodeClass::nodeType_solutionStructuralTemperature:
            case SimulationNodeClass::nodeType_solutionStructuralThermalStrain:
            case SimulationNodeClass::nodeType_solutionStructuralTotalStrain:
            case SimulationNodeClass::nodeType_solutionStructuralEquivalentPlasticStrain:
            case SimulationNodeClass::nodeType_solutionStructuralNodalForces:
            {
                data.setValue(0); editor->addItem("Time",data);
                data.setValue(1); editor->addItem("Set",data);
            }
                break;

            case SimulationNodeClass::nodeType_meshMethod:
            {
                data.setValue(0); editor->addItem("Level",data);
                data.setValue(1); editor->addItem("Pair distance",data);
            }
                break;
            }
            return editor;
        }
        //! --------
        //! "Type "
        //! --------
        else if(propertyName == "Type ")
        {
            DetailViewer *theDetailViewer = static_cast<DetailViewer*>(parent->parent());
            SimulationNodeClass *node = theDetailViewer->getNode();
            SimulationNodeClass::nodeType type = node->getType();

            QComboBox *editor;
            QVariant data;

            switch(type)
            {
            case SimulationNodeClass::nodeType_solutionThermalTemperature:
            case SimulationNodeClass::nodeType_solutionThermalFlux:
                return Q_NULLPTR;
                break;

            case SimulationNodeClass::nodeType_solutionStructuralEquivalentPlasticStrain:
            {
                editor = new QComboBox(parent);
                data.setValue(0); editor->addItem("Equivalent plastic strain",data);
            }
                break;

            case SimulationNodeClass::nodeType_solutionStructuralTemperature:
            {
                editor = new QComboBox(parent);
                data.setValue(0); editor->addItem("Temperature",data);
            }
                break;
            case SimulationNodeClass::nodeType_solutionStructuralNodalDisplacement:
            {
                editor = new QComboBox(parent);
                data.setValue(0); editor->addItem("Total displacement",data);
                data.setValue(1); editor->addItem("Directional displacement X",data);
                data.setValue(2); editor->addItem("Directional displacement Y",data);
                data.setValue(3); editor->addItem("Directional displacement Z",data);
            }
                break;

            case SimulationNodeClass::nodeType_solutionStructuralNodalForces:
            {
                editor = new QComboBox(parent);
                data.setValue(0); editor->addItem("Total",data);
                data.setValue(1); editor->addItem("X direction",data);
                data.setValue(2); editor->addItem("Y direction",data);
                data.setValue(3); editor->addItem("Z direction",data);
            }
                break;
            case SimulationNodeClass::nodeType_solutionStructuralStress:
            {
                editor = new QComboBox(parent);
                data.setValue(0); editor->addItem("Equivalent stress",data);
                data.setValue(1); editor->addItem("Stress intensity",data);
                data.setValue(2); editor->addItem("Maximum principal stress",data);
                data.setValue(3); editor->addItem("Middle principal stress",data);
                data.setValue(4); editor->addItem("Minimum principal stress",data);
                data.setValue(5); editor->addItem("Normal stress X",data);
                data.setValue(6); editor->addItem("Normal stress Y",data);
                data.setValue(7); editor->addItem("Normal stress Z",data);
                data.setValue(8); editor->addItem("Shear stress XY",data);
                data.setValue(9); editor->addItem("Shear stress YZ",data);
                data.setValue(10); editor->addItem("Shear stress ZX",data);
            }
                break;

            case SimulationNodeClass::nodeType_solutionStructuralTotalStrain:
            {
                editor = new QComboBox(parent);
                data.setValue(0); editor->addItem("Equivalent total strain",data);
                data.setValue(1); editor->addItem("Total strain intensity",data);
                data.setValue(2); editor->addItem("Total maximum principal strain",data);
                data.setValue(3); editor->addItem("Total middle principal strain",data);
                data.setValue(4); editor->addItem("Total minimum principal strain",data);
            }
                break;

            case SimulationNodeClass::nodeType_solutionStructuralMechanicalStrain:
            {
                editor = new QComboBox(parent);
                data.setValue(0); editor->addItem("Equivalent mechanical strain",data);
                data.setValue(1); editor->addItem("Mechanical strain intensity",data);
                data.setValue(2); editor->addItem("Mechanical maximum principal strain",data);
                data.setValue(3); editor->addItem("Mechanical middle principal strain",data);
                data.setValue(4); editor->addItem("Mechanical minimum principal strain",data);
            }
                break;

            case SimulationNodeClass::nodeType_solutionStructuralThermalStrain:
            {
                editor = new QComboBox(parent);
                data.setValue(0); editor->addItem("Equivalent thermal strain",data);
                data.setValue(1); editor->addItem("Thermal strain intensity",data);
                data.setValue(2); editor->addItem("Thermal maximum principal strain",data);
                data.setValue(3); editor->addItem("Thermal middle principal strain",data);
                data.setValue(4); editor->addItem("Thermal minimum principal strain",data);
            }
                break;

            case SimulationNodeClass::nodeType_solutionStructuralContact:
            {
                editor = new QComboBox(parent);
                data.setValue(0); editor->addItem("Contact pressure",data);
                data.setValue(1); editor->addItem("Frictional stress",data);
                data.setValue(2); editor->addItem("Penetration",data);
                data.setValue(3); editor->addItem("Sliding",data);
            }
                break;
            }
            return editor;
        }
        //! -----------------------
        //! "Solution information"
        //! -----------------------
        else if(propertyName =="Solution information")
        {
            //! ------------------------------
            //! 0 => solver txt messages
            //! 1 => force convergence
            //! 2 => displacement convergence
            //! 3 => line search
            //! 4 => time step size
            //! ------------------------------
            QComboBox *editor = new QComboBox(parent);
            QVariant data;
            data.setValue(0); editor->addItem("Solver output",data);
            data.setValue(1); editor->addItem("Force convergence",data);
            data.setValue(2); editor->addItem("Displacement convergence",data);
            data.setValue(3); editor->addItem("Line search",data);
            data.setValue(4); editor->addItem("Time step size",data);
            return editor;
        }
        //! ------------------------------------------------------------------------------
        //! "Coupling" type
        //! if a "Remote force" or a "Remote displacement" has a "Scoping method" defined
        //! through "Remote point", do not return the editor
        //! ------------------------------------------------------------------------------
        else if(propertyName == "Coupling")
        {
            SimulationNodeClass *curNode = this->getCurrentNode();
            SimulationNodeClass::nodeType nodeType = curNode->getType();
            if(nodeType == SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteForce ||
                    nodeType == SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteDisplacement ||
                    nodeType == SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteRotation)
            {
                Property::ScopingMethod scopingMethod = curNode->getPropertyValue<Property::ScopingMethod>("Scoping method");
                if(scopingMethod == Property::ScopingMethod_RemotePoint) return 0;
            }

            QComboBox *editor = new QComboBox(parent);
            QVariant data;
            data.setValue(0); editor->addItem("Kinematic",data);
            data.setValue(1); editor->addItem("Distributed",data);
            return editor;
        }
        //! ------------------------------------------------------------------------------
        //! "DOFs selection"
        //! if a "Remote force" or a "Remote displacement" has a "Scoping method" defined
        //! through "Remote point", do not return the editor
        //! ------------------------------------------------------------------------------
        else if(propertyName == "DOFs selection")
        {
            SimulationNodeClass *curNode = this->getCurrentNode();
            SimulationNodeClass::nodeType nodeType = curNode->getType();
            if(nodeType == SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteForce ||
                    nodeType == SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteDisplacement ||
                    nodeType == SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteRotation)
            {
                Property::ScopingMethod scopingMethod = curNode->getPropertyValue<Property::ScopingMethod>("Scoping method");
                if(scopingMethod == Property::ScopingMethod_RemotePoint) return 0;
            }

            QComboBox *editor = new QComboBox(parent);
            QVariant data;
            data.setValue(0); editor->addItem("Program controlled",data);
            data.setValue(1); editor->addItem("Manual",data);
            return editor;
        }
        //! ------------------------------------------------------------------------------
        //! coupling DOFs
        //! if a "Remote force" or a "Remote displacement" has a "Scoping method" defined
        //! through "Remote point", do not return the editor
        //! ------------------------------------------------------------------------------
        else if(propertyName == "X component " || propertyName == "Y component " || propertyName == "Z component " ||
                propertyName == "X rotation" || propertyName == "Y rotation" || propertyName == "Z rotation")
        {
            SimulationNodeClass *curNode = this->getCurrentNode();
            SimulationNodeClass::nodeType nodeType = curNode->getType();
            if(nodeType == SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteForce ||
                    nodeType == SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteDisplacement ||
                    nodeType == SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteRotation)
            {
                Property::ScopingMethod scopingMethod = curNode->getPropertyValue<Property::ScopingMethod>("Scoping method");
                if(scopingMethod == Property::ScopingMethod_RemotePoint) return 0;
            }

            QComboBox *editor = new QComboBox(parent);
            QVariant data;
            data.setValue(0); editor->addItem("Inactive",data);
            data.setValue(1); editor->addItem("Active",data);
            return editor;
        }
        //! --------------------
        //! "Number of threads"
        //! --------------------
        else if(propertyName == "Number of threads")
        {
            QSpinBox *editor = new QSpinBox(parent);
            editor->setMinimum(1);
            editor->setMaximum(256);
            return editor;
        }
        //! -----------------------------
        //! "First color" "Second color"
        //! -----------------------------
        else if(propertyName == "First color" || propertyName == "Second color")
        {
            ColorSelector *editor = new ColorSelector(parent);
            return editor;
        }
        //! -----------
        //! "Gradient"
        //! -----------
        else if(propertyName == "Gradient")
        {
            QComboBox *editor = new QComboBox(parent);
            QVariant data;
            data.setValue(0); editor->addItem("Uniform");
            data.setValue(1); editor->addItem("Horizontal");
            data.setValue(2); editor->addItem("Vertical");
            return editor;
        }
        //! ----------------
        //! "Analysis time"
        //! ----------------
        else if(propertyName == "Analysis time")
        {
            QLineEdit *editor = new QLineEdit(parent);
            QDoubleValidator *doubleValidator = new QDoubleValidator();
            doubleValidator->setBottom(0.0);
            editor->setValidator(doubleValidator);
            return editor;
        }
        //! ----------------------------
        //! "Buckets" "Remapping steps"
        //! ----------------------------
        else if(propertyName =="X buckets" ||
                propertyName =="Y buckets" ||
                propertyName =="Z buckets" ||
                propertyName == "Remapping steps")
        {
            QLineEdit *editor = new QLineEdit(parent);
            QIntValidator *intValidator = new QIntValidator();
            intValidator->setBottom(1);
            editor->setValidator(intValidator);
            return editor;
        }
        //! ------------
        //! "Algorithm"
        //! ------------
        else if(propertyName =="Algorithm")
        {
            QComboBox *editor = new QComboBox(parent);
            QVariant data;
            switch(this->getCurrentNode()->getType())
            {
            case SimulationNodeClass::nodeType_importedBodyScalar:
            {
                data.setValue(0); editor->addItem("Nearest point 1",data);
                data.setValue(1); editor->addItem("Nearest point 2",data);
                data.setValue(2); editor->addItem("Shape functions",data);
            }
                break;
            case SimulationNodeClass::nodeType_meshPrismaticLayer:
            {
                data.setValue(0); editor->addItem("Pre",data);
                data.setValue(1); editor->addItem("Post",data);
            }
                break;
            }
            return editor;
        }
        //! ---------------------
        //! "Boundary mesh type"
        //! ---------------------
        if(propertyName =="Boundary mesh type")
        {
            QComboBox *editor = new QComboBox(parent);
            QVariant data;
            data.setValue(0); editor->addItem("Hybrid",data);
            data.setValue(1); editor->addItem("Tetrahedral",data);
            return editor;
        }
        //! -------------
        //! "Split data"
        //! -------------
        else if(propertyName == "Split data")
        {
            QComboBox *editor = new QComboBox(parent);
            QVariant data;
            data.setValue(0); editor->addItem("Single file",data);
            data.setValue(1); editor->addItem("Multiple files",data);
            return editor;
        }
        //! ------------
        //! "Relevance"
        //! ------------
        else if(propertyName == "Relevance")
        {
            QSlider *editor = new QSlider(parent);
            editor->setStyleSheet(sliderStyleSheet);
            editor->setOrientation(Qt::Horizontal);
            editor->setRange(-100,100);
            return editor;
        }
        //! --------------------
        //! "Initial size seed"
        //! --------------------
        else if(propertyName =="Initial size seed")
        {
            QComboBox *editor = new QComboBox(parent);
            QVariant data;
            data.setValue(0); editor->addItem("Program controlled",data);
            data.setValue(1); editor->addItem("Assembly",data);
            data.setValue(2); editor->addItem("Part");
            return editor;
        }
        else if(propertyName =="Element midside nodes")
        {
            QComboBox *editor = new QComboBox(parent);
            QVariant data;
            data.setValue(0); editor->addItem("Program controlled",data);
            data.setValue(1); editor->addItem("Dropped",data);
            data.setValue(2); editor->addItem("Kept");
            return editor;
        }
        else if(propertyName =="Straight sided elements")
        {
            QComboBox *editor = new QComboBox(parent);
            data.setValue(0); editor->addItem("No");
            data.setValue(1); editor->addItem("Yes");
            return editor;
        }
        else if(propertyName =="Remap")
        {
            QComboBox *editor = new QComboBox(parent);
            QVariant data;
            data.setValue(true); editor->addItem("Yes",data);
            data.setValue(false); editor->addItem("No",data);
            //connect(editor,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseRemapFlagControl()));
            return editor;
        }
        //! ------------
        //! "Submeshes"
        //! ------------
        else if(propertyName =="Submeshes")
        {
            QComboBox *editor = new QComboBox(parent);
            QVariant data;
            data.setValue(true);
            editor->addItem("At meshing time",data);
            data.setValue(false);
            editor->addItem("Deferred",data);
            connect(editor,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseSubmeshesControl()));
            return editor;
        }
        else if(propertyName =="Generate")
        {
            QPushButton *editor = new QPushButton(parent);
            editor->setText("Go!");
            connect(editor,SIGNAL(pressed()),this,SLOT(commitAndCloseGenerate()));
            return editor;
        }
        else if(propertyName =="Step number")
        {
            QSpinBox *editor = new QSpinBox(parent);
            editor->setMinimum(1);

            if(getCurrentNode()->getType()==SimulationNodeClass::nodeType_structuralAnalysisSettings)
            {
                //! --------------------
                //! set a maximum value
                //! --------------------
                DetailViewer *theDetailViewer = static_cast<DetailViewer*>(parent->parent());
                SimulationManager *sm = theDetailViewer->parent()->parent()->findChild<SimulationManager*>();
                int max = sm->getAnalysisSettingsNodeFromCurrentItem()->getPropertyItem("Number of steps")->
                        data(Qt::UserRole).value<Property>().getData().toInt();
                editor->setMaximum(max);
                return editor;
            }
            return 0;
        }
        else if(propertyName =="Step selection mode")
        {
            QComboBox *editor = new QComboBox(parent);
            editor->addItem("All steps",0);
            editor->addItem("First",1);
            editor->addItem("Last",2);
            editor->addItem("By number",3);
            editor->addItem("Automatic time history", 4);
            return editor;
        }
        else if(propertyName =="Source file")
        {
            QFileDialog *editor = new QFileDialog(parent);
            editor->setWindowTitle("Select file");
            editor->setAcceptMode(QFileDialog::AcceptOpen);
            editor->setViewMode(QFileDialog::Detail);
            editor->setFileMode(QFileDialog::ExistingFiles);
            editor->setNameFilter("Text files (*.txt *.dat)");
            //cout<<"GeneralDelegate::getWorkingDir()->____function called____"<<this->getWorkingDir().toStdString()<<"____"<<endl;
            editor->setDirectory(this->getWorkingDir());
            return editor;
        }
        else if(propertyName =="Source directory" || propertyName == "Target directory")
        {
            QFileDialog *editor = new QFileDialog(parent);
            editor->setWindowTitle("Select file");
            editor->setAcceptMode(QFileDialog::AcceptOpen);
            editor->setViewMode(QFileDialog::Detail);
            editor->setFileMode(QFileDialog::Directory);
            editor->setNameFilter("Text files (*.txt *.dat)");
            //cout<<"GeneralDelegate::getWorkingDir()->____function called____"<<this->getWorkingDir().toStdString()<<"____"<<endl;
            editor->setDirectory(this->getWorkingDir());
            return editor;
        }
        else if(propertyName =="Smoothing")
        {
            QComboBox *editor = new QComboBox(parent);
            QVariant data;
            data.setValue(0);
            editor->addItem("Off",data);
            data.setValue(1);
            editor->addItem("Low",data);
            data.setValue(2);
            editor->addItem("Medium",data);
            data.setValue(3);
            editor->addItem("High",data);
            return editor;
        }
        else if(propertyName =="Show mesh nodes")
        {
            QComboBox *editor = new QComboBox(parent);
            QVariant data;
            data.setValue(false);
            editor->addItem("Off",data);
            data.setValue(true);
            editor->addItem("Active",data);
            return editor;
        }
        else if(propertyName =="Environment temperature" || propertyName =="X coordinate"
                || propertyName =="Y coordinate" || propertyName == "Z coordinate")
        {
            QLineEdit *editor = new QLineEdit(parent);
            QDoubleValidator *doubleValidator = new QDoubleValidator();
            editor->setValidator(doubleValidator);
            return editor;
        }
        else if(propertyName =="Initial substeps" || propertyName =="Minimum substeps" || propertyName =="Maximum substeps"
                || propertyName =="Number of substeps")
        {
            QLineEdit *editor = new QLineEdit(parent);
            QIntValidator *intValidator = new QIntValidator();
            if(propertyName == "Number of substeps") intValidator->setBottom(1);
            if(propertyName == "Minimum time step")
            {
                 intValidator->setBottom(1);
                 int maxSubSteps = this->getCurrentNode()->getPropertyValue<int>("Maximum substeps");
                 intValidator->setTop(maxSubSteps);
            }
            if(propertyName == "Maximum substeps")
            {
                int minSubSteps = this->getCurrentNode()->getPropertyValue<int>("Minimum substeps");
                intValidator->setBottom(minSubSteps);
            }
            if(propertyName =="Initial substeps")
            {
                int minSubSteps = this->getCurrentNode()->getPropertyValue<int>("Minimum substeps");
                intValidator->setBottom(minSubSteps);
            }
            editor->setValidator(intValidator);
            connect(editor,SIGNAL(editingFinished()),this,SLOT(commitAndCloseSubstepEditors()));
            return editor;
        }
        else if(propertyName =="Offset X" || propertyName =="Offset Y" || propertyName =="Offset Z" ||
                propertyName =="Rotation X" || propertyName =="Rotation Y" || propertyName =="Rotation Z")
        {
            QLineEdit *editor = new QLineEdit(parent);
            QDoubleValidator *doubleValidator = new QDoubleValidator();
            editor->setValidator(doubleValidator);
            return editor;
        }
        else if((propertyName =="Origin X" || propertyName =="Origin Y" || propertyName =="Origin Z") &&
                (index.model()->parent(index).child(0,1)).data(Qt::UserRole).value<Property>().getData().value<Property::defineBy>()==Property::defineBy_globalCoordinates)
        {
            //! return the editor for the coordinates of the origin only if "Define by" is "Global coordinates"
            QLineEdit *editor = new QLineEdit(parent);
            QDoubleValidator *doubleValidator = new QDoubleValidator();
            editor->setValidator(doubleValidator);
            return editor;
        }
        //! -------------------------------------------
        //! "Film coefficient" "Reference temperature"
        //! -------------------------------------------
        else if(propertyName =="Film coefficient" || propertyName =="Reference temperature")
        {
            bool editorHasFreeOption = false;
            LineEdit *editor = new LineEdit(editorHasFreeOption, parent);
            return editor;
        }
        //! ------------
        //! "Magnitude"
        //! ------------
        else if (propertyName =="Magnitude")
        {
            SimulationNodeClass *theNode = this->getCurrentNode();
            LineEdit *editor;
            if(theNode->getType()==SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Displacement ||
                    theNode->getType()==SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteDisplacement ||
                    theNode->getType()==SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteRotation)
            {
                //! -----------------------------------------------------------
                //! the editing tool has the option "free" for "Displacement",
                //! but only when "defineBy" is "defineBy_components"
                //! -----------------------------------------------------------
                Property::defineBy theDefineBy = theNode->getPropertyValue<Property::defineBy>("Define by");
                if(theDefineBy == Property::defineBy_components)
                {
                    //! activate the option "free"
                    editor = new LineEdit(true, parent);
                }
                else
                {
                    //! deactivate the option "free"
                    editor = new LineEdit(false, parent);
                }
            }
            else
            {
                //! the editing tool has not the option "free"
                editor = new LineEdit(false, parent);
            }
            return editor;
        }
        //! --------------
        //! "X component"
        //! --------------
        else if (propertyName =="X component")
        {
            SimulationNodeClass *theNode = this->getCurrentNode();
            LineEdit *editor;
            if(theNode->getType()==SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Displacement ||
                    theNode->getType()==SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteDisplacement ||
                    theNode->getType()==SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteRotation)
            {
                //! the editing tool has the option "free"
                editor = new LineEdit(true, parent);
            }
            else
            {
                //! the editing tool has not the option "free"
                editor = new LineEdit(false, parent);
            }
            return editor;
        }
        //! --------------
        //! "Y component"
        //! --------------
        else if (propertyName =="Y component")
        {
            SimulationNodeClass *theNode = this->getCurrentNode();
            LineEdit *editor;
            if(theNode->getType()==SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Displacement ||
                    theNode->getType()==SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteDisplacement ||
                    theNode->getType()==SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteRotation)
            {
                //! the editing tool has the option "free"
                editor = new LineEdit(true, parent);
            }
            else
            {
                //! the editing tool has not the option "free"
                editor = new LineEdit(false, parent);
            }
            return editor;
        }
        //! --------------
        //! "Z component"
        //! --------------
        else if (propertyName =="Z component")
        {
            SimulationNodeClass *theNode = this->getCurrentNode();
            LineEdit *editor;
            if(theNode->getType()==SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Displacement ||
                    theNode->getType()==SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteDisplacement ||
                    theNode->getType()==SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteRotation)
            {
                //! the editing tool has the option "free"
                editor = new LineEdit(true, parent);
            }
            else
            {
                //! the editing tool has not the option "free"
                editor = new LineEdit(false, parent);
            }
            return editor;
        }
        //! ------------------------------------------------------------------------------
        //! Overpressure
        //! "Overpressure" can be edited for a frictional or frictionless contact
        //! If the contact is "Bonded" no overpressure function is defined
        //! If the contact is "Tied" only the "tied" overpressure function can be defined
        //! ------------------------------------------------------------------------------
        else if(propertyName =="Overpressure")
        {
            SimulationNodeClass *curNode = this->getCurrentNode();
            Property::contactType type = curNode->getPropertyItem("Type")->data(Qt::UserRole).value<Property>().getData().value<Property::contactType>();
            if(type!= Property::contactType_tied && type!= Property::contactType_bonded)
            {
                QComboBox *editor = new QComboBox(parent);
                QVariant data;
                data.setValue(Property::overpressureFunction_linear);
                editor->addItem("Linear",data);
                data.setValue(Property::overpressureFunction_exponential);
                editor->addItem("Exponential",data);
                return editor;
            }
            else return 0;
        }
        //! ---------------------------------------------------------------------
        //! "Behavior"
        //! Do not return an editor for "Tied", since the only behavior defined
        //! for that case is "Symmetric"
        //! ---------------------------------------------------------------------
        else if(propertyName =="Behavior")
        {
            SimulationNodeClass *curNode = this->getCurrentNode();
            Property::contactType type = curNode->getPropertyItem("Type")->data(Qt::UserRole).value<Property>().getData().value<Property::contactType>();
            if(type!=Property::contactType_tied)
            {
                QComboBox *editor = new QComboBox(parent);
                editor->clear();
                QVariant data;
                data.setValue(Property::contactBehavior_asymmetric);
                editor->addItem("Asymmetric",data);
                data.setValue(Property::contactBehavior_symmetric);
                editor->addItem("Symmetric",data);
                return editor;
            }
            else return 0;
        }
        else if(propertyName =="Sizing type")
        {
            cout<<"Creating editor for edge mesh sizing"<<endl;
            QComboBox *editor = new QComboBox(parent);
            editor->clear();
            data.setValue(0); editor->addItem("Element size",data);
            data.setValue(1); editor->addItem("Number of divisions",data);
            return editor;
        }
        //! -------
        //! "Type"
        //! -------
        else if(propertyName =="Type")
        {
            //! -----------------------------------------------------------------
            //! non standard treatment of the dummy node for multiple selections
            //! -----------------------------------------------------------------
            DetailViewer *theDetailViewer = static_cast<DetailViewer*>(parent->parent());
            SimulationNodeClass *aNode = theDetailViewer->getCurrentMultipleSelectionNode();
            if(aNode!=Q_NULLPTR)
            {
                if(aNode->getPropertyItem("Type")->data(Qt::UserRole).value<Property>().getData().isValid()==false)
                {
                    //gildotta
                }
            }
            //! -----------------
            //! end experimental
            //! -----------------

            QVariant data;
            QComboBox *editor = new QComboBox(parent);
            SimulationNodeClass *node = this->getCurrentNode();
            switch(node->getType())
            {
            case SimulationNodeClass::nodeType_connectionPair:
            {
                cout<<"Creating editor for \"Contact type\""<<endl;
                editor->clear();
                data.setValue(Property::contactType_bonded);
                editor->addItem("Bonded",data);
                data.setValue(Property::contactType_frictional);
                editor->addItem("Frictional",data);
                data.setValue(Property::contactType_frictionless);
                editor->addItem("Frictionless",data);
                data.setValue(Property::contactType_tied);
                editor->addItem("Tied",data);
            }
                break;
            }
            return editor;
        }
        //! --------------------
        //! "Coordinate system"
        //! --------------------
        else if(propertyName =="Coordinate system")
        {
            cout<<"Creating editor for Coordinate system"<<endl;
            SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
            QExtendedStandardItem *itemCSRoot = sm->getTreeItem(SimulationNodeClass::nodeType_coordinateSystems);

            QComboBox *editor = new QComboBox(parent);
            for(int k=0; k<itemCSRoot->rowCount(); k++)
            {
                QStandardItem *itemCScurrent = (QStandardItem*)(itemCSRoot->child(k,0));
                QString name = itemCScurrent->data(Qt::DisplayRole).toString();
                void *itemCScurrentVoid = (void*)itemCScurrent;
                data.setValue(itemCScurrentVoid);
                editor->addItem(name,data);
            }
            cout<<"Editor created"<<endl;
            return editor;
        }
        //! ---------------------------------------------
        //! "Named selection" "Boundary named selection"
        //! ---------------------------------------------
        else if(propertyName =="Named selection" || propertyName =="Boundary named selection")
        {
            SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
            QExtendedStandardItem *itemNSRoot = sm->getTreeItem(SimulationNodeClass::nodeType_namedSelection);

            QComboBox *editor = new QComboBox(parent);
            for(int k=0; k<itemNSRoot->rowCount(); k++)
            {
                itemNSRoot->child(k,0);
                QExtendedStandardItem *itemNScurrent = dynamic_cast<QExtendedStandardItem*>(itemNSRoot->child(k,0));
                QString name = itemNScurrent->data(Qt::DisplayRole).toString();
                void *itemNScurrentVoid = (void*)itemNScurrent;
                data.setValue(itemNScurrentVoid);
                editor->addItem(name,data);
            }
            return editor;
        }
        //! ------------------------------
        //! "Contact pair" (model change)
        //! ------------------------------
        else if(propertyName =="Contact pair")
        {
            SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
            QExtendedStandardItem *itemConnectionRoot = sm->getTreeItem(SimulationNodeClass::nodeType_connection);

            QComboBox *editor = new QComboBox(parent);

            //! ----------------------------------------------------------------
            //! add the dummy item "Select from list", which is the first child
            //! ----------------------------------------------------------------
            QStandardItem *itemDummyContactPair = itemConnectionRoot->child(0,0);
            QString name = itemDummyContactPair->data(Qt::DisplayRole).toString();
            void *p = (void*)itemDummyContactPair;
            data.setValue(p);
            editor->addItem(name,data);

            //! ---------------------------------------------
            //! we know that by definition the other contact
            //! items are placed under a connection group
            //! ---------------------------------------------
            for(int k=1; k<itemConnectionRoot->rowCount(); k++)
            {
                QStandardItem *childItem = itemConnectionRoot->child(k,0);
                for(int j=0; j<childItem->rowCount(); j++)
                {
                    QStandardItem *itemContactPair = childItem->child(j,0);
                    QString name = itemContactPair->data(Qt::DisplayRole).toString();
                    void *p = (void*)itemContactPair;
                    data.setValue(p);
                    editor->addItem(name,data);
                }
            }
            return editor;
        }
        //! -----------------
        //! "Master" "Slave"
        //! -----------------
        else if (propertyName =="Master" || propertyName =="Slave")
        {
            cout<<"____creating editor for master or slave____"<<endl;
            //! -----------------------------------------------------------------
            //! non standard treatment of the dummy node for multiple selections
            //! -----------------------------------------------------------------
            DetailViewer *theDetailViewer = static_cast<DetailViewer*>(parent->parent());
            SimulationNodeClass *aNode = theDetailViewer->getCurrentMultipleSelectionNode();
            if(aNode!=Q_NULLPTR)
            {
                if(aNode->getPropertyItem("Master")->data(Qt::UserRole).value<Property>().getData().isValid()==false) return Q_NULLPTR;
                if(aNode->getPropertyItem("Slave")->data(Qt::UserRole).value<Property>().getData().isValid()==false) return Q_NULLPTR;
            }
            //! -----------------
            //! end experimental
            //! -----------------
            SimulationNodeClass *node = this->getCurrentNode();
            Property::ScopingMethod theScopingMethod = node->getPropertyValue<Property::ScopingMethod>("Scoping method");
            switch(theScopingMethod)
            {
            case Property::ScopingMethod_GeometrySelection:
            case Property::ScopingMethod_Automatic:
            {
                //! -------------------------------------------------------------------------
                //! the scope of master/slave is defined through a direct geometry selection
                //! -------------------------------------------------------------------------
                ShapeSelector *editor = new ShapeSelector(myCTX,parent);
                return editor;
            }
                break;

            case Property::ScopingMethod_NamedSelection:
            {
                //! -----------------------------------------------------------------
                //! the scope of master/slave is defined through a "Named selection"
                //! -----------------------------------------------------------------
                SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
                QExtendedStandardItem *itemNSRoot = sm->getTreeItem(SimulationNodeClass::nodeType_namedSelection);

                QComboBox *editor = new QComboBox(parent);
                for(int k=0; k<itemNSRoot->rowCount(); k++)
                {
                    itemNSRoot->child(k,0);
                    QExtendedStandardItem *itemNScurrent = dynamic_cast<QExtendedStandardItem*>(itemNSRoot->child(k,0));
                    QString name = itemNScurrent->data(Qt::DisplayRole).toString();
                    void *itemNScurrentVoid = (void*)itemNScurrent;
                    data.setValue(itemNScurrentVoid);
                    editor->addItem(name,data);
                }
                return editor;
            }
                break;
            }
        }
        //! ------------
        //! "Direction"
        //! ------------
        else if(propertyName =="Direction")
        {
            DirectionSelector *editor = new DirectionSelector(myCTX,parent);
            //! warning: parent->parent() -  10/12/2017
            DetailViewer *theDetailViewer = static_cast<DetailViewer*>(parent->parent());
            //!cout<<theDetailViewer->metaObject()->className()<<endl; // 10/12/2017 - identity of the parent
            connect(theDetailViewer,SIGNAL(requestHandleSelectionChanged()),editor,SLOT(showMarker()));
            connect(editor,SIGNAL(editingFinished()),this,SLOT(commitAndCloseDirectionSelector()));
            return editor;
        }
        //! --------------
        //! "Bolt status"
        //! --------------
        else if(propertyName =="Bolt status")
        {
            QComboBox *editor = new QComboBox(parent);
            editor->clear();
            QVariant data;

            data.setValue(Property::boltStatusDefinedBy_load);
            editor->addItem("Load",data);
            data.setValue(Property::boltStatusDefinedBy_adjustment);
            editor->addItem("Adjustment",data);
            data.setValue(Property::boltStatusDefinedBy_open);
            editor->addItem("Open");
            data.setValue(Property::boltStatusDefinedBy_lock);
            editor->addItem("Lock");

            return editor;
        }
        //! ------------
        //! "Define by"
        //! ------------
        else if(propertyName =="Define by")
        {
            SimulationNodeClass *node = this->getCurrentNode();
            SimulationNodeClass::nodeType type = node->getType();

            if(type!=SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Pressure &&
                    type!=SimulationNodeClass::nodeType_structuralAnalysisBoltPretension)
            {
                QComboBox *editor = new QComboBox(parent);
                QVariant data;
                data.setValue(Property::defineBy_components);
                editor->addItem("Components",data);
                data.setValue(Property::defineBy_vector);
                editor->addItem("Vector",data);

                //! ---------------------------------------------------------------------------
                //! do not add the option "Normal to" in case of a "Force" (or "Remote force")
                //! since physically it would be a "Pressure", which has it own nummerical
                //! representation (CCX => surface defined by elements"
                //! ---------------------------------------------------------------------------
                if(type!=SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Force &&
                        type!=SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteForce)
                {
                    data.setValue(Property::defineBy_normal);
                    editor->addItem("Normal to");
                }
                return editor;
            }
            /*
            if(type==SimulationNodeClass::nodeType_structuralAnalysisBoltPretension)
            {
                QComboBox *editor = new QComboBox(parent);
                editor->clear();
                QVariant data;

                data.setValue(Property::boltStatusDefinedBy_load);
                editor->addItem("Load",data);
                data.setValue(Property::boltStatusDefinedBy_adjustment);
                editor->addItem("Adjustment",data);
                data.setValue(Property::boltStatusDefinedBy_open);
                editor->addItem("Open");
                data.setValue(Property::boltStatusDefinedBy_lock);
                editor->addItem("Lock");

                return editor;
            }
            */
            return 0;
        }
        //! -------------
        //! "Define by "
        //! -------------
        else if(propertyName =="Define by ") //! warning: observe the space at the end
        {
            QComboBox *editor = new QComboBox(parent);
            editor->clear();
            QVariant data;
            data.setValue(Property::defineBy_geometrySelection);
            editor->addItem("Geometry selection",data);
            data.setValue(Property::defineBy_globalCoordinates);
            editor->addItem("Global coordinates",data);
            return editor;
        }
        //! --------------------
        //! "Load" "Adjustment"
        //! --------------------
        else if (propertyName=="Load" || propertyName == "Adjustment")
        {
            QLineEdit *editor = new QLineEdit(parent);
            QDoubleValidator *doubleValidator = new QDoubleValidator();
            doubleValidator->setBottom(0.0);
            editor->setValidator(doubleValidator);
            return editor;
        }
        //! ----------------------------------------------
        //! "Ambient" "Diffuse" "Specular" "Transparency"
        //! ----------------------------------------------
        else if (propertyName=="Ambient" || propertyName=="Diffuse" || propertyName=="Specular" || propertyName=="Transparency")
        {
            QLineEdit *editor = new QLineEdit(parent);
            QDoubleValidator *doubleValidator = new QDoubleValidator();
            doubleValidator->setDecimals(1);
            editor->setValidator(doubleValidator);
            return editor;
        }
        //! ---------------
        //! "Transparency"
        //! ---------------
        else if (propertyName=="Transparency")
        {
            QLineEdit *editor = new QLineEdit(parent);
            QDoubleValidator *doubleValidator = new QDoubleValidator();
            doubleValidator->setRange(0,1,1);
            editor->setValidator(doubleValidator);
            //connect(editor,SIGNAL(editingFinished()),this,SLOT(commitAndCloseLineEditTransparency()));
            return editor;
        }
        //! ------------------
        //! "Element control"
        //! ------------------
        else if(propertyName =="Element control")
        {
            QComboBox *editor = new QComboBox(parent);
            editor->clear();
            QVariant data;
            data.setValue(Property::elementControl_programControlled);
            editor->addItem("Program controlled",data);
            data.setValue(Property::elementControl_manual);
            editor->addItem("Manual",data);
            return editor;
        }
        else if(propertyName =="Integration scheme")
        {
            QComboBox *editor = new QComboBox(parent);
            editor->clear();
            QVariant data;
            data.setValue(Property::integrationScheme_full);
            editor->addItem("Full",data);
            data.setValue(Property::integrationScheme_reduced);
            editor->addItem("Reduced",data);
            return editor;
        }
        //! ------------------------------
        //! "Radial" "Axial" "Tangential"
        //! ------------------------------
        else if(propertyName =="Radial" || propertyName =="Axial" || propertyName =="Tangential")
        {
            QComboBox *editor = new QComboBox(parent);
            editor->clear();
            QVariant data;
            data.setValue(Property::DOFfreedom_fixed);
            editor->addItem("Fixed",data);
            data.setValue(Property::DOFfreedom_free);
            editor->addItem("Free",data);
            return editor;
        }
        //! ----------------------------------------------------
        //! "Min element size" "Max element size" "Grading" ...
        //! ----------------------------------------------------
        else if(propertyName=="Min element size" || propertyName=="Max element size" || propertyName == "Grading"
                || propertyName=="Face sizing" || propertyName == "Element size" || propertyName =="Pinball")
        {
            QLineEdit *editor = new QLineEdit(parent);
            QDoubleValidator *doubleValidator = new QDoubleValidator();
            doubleValidator->setBottom(0.0);
            editor->setValidator(doubleValidator);
            return editor;
        }
        //! -----------------------
        //! "Number of divisions"
        //! -----------------------
        else if(propertyName=="Number of divisions")
        {
            QSpinBox *editor = new QSpinBox(parent);
            editor->setMinimum(1);
            return editor;
        }
        //! ------------
        //! "Mesh type"
        //! ------------
        else if(propertyName =="Mesh type")
        {
            QComboBox *editor = new QComboBox(parent);
            editor->addItem("Full tri-tet",0);
            editor->addItem("Full quad-hexa",1);
            editor->addItem("Quad-hexa dominant",2);
            return editor;
        }
        //! -----------------
        //! "Mesh engine 2D"
        //! Program controlled
        //! Netgen
        //! Express mesh
        //! OCC STL
        //! -----------------
        else if (propertyName=="Mesh engine 2D")
        {
            QComboBox *editor = new QComboBox(parent);
            editor->clear();
            QVariant data;
            data.setValue(Property::meshEngine2D_ProgramControlled);
            editor->addItem("Program controlled",data);
            data.setValue(Property::meshEngine2D_Netgen);
            editor->addItem("Netgen 2D",data);
            data.setValue(Property::meshEngine2D_OCC_ExpressMesh);
            editor->addItem("Express mesh",data);
            data.setValue(Property::meshEngine2D_OCC_STL);
            editor->addItem("OCC STL mesh",data);
            return editor;
        }
        else if (propertyName=="Mesh engine 3D")
        {
            QComboBox *editor = new QComboBox(parent);
            editor->addItem("Netgen",0);
            editor->addItem("Tetgen",1);
            return editor;
        }
        else if (propertyName=="Surface mesh type")
        {
            QComboBox *editor = new QComboBox(parent);
            editor->addItem("All triangles",0);
            editor->addItem("Quad dominant",1);
            return editor;
        }
        else if (propertyName=="Volume mesh type")
        {
            QComboBox *editor = new QComboBox(parent);
            editor->addItem("All tet",0);
            return editor;
        }
        else if (propertyName=="Mesh order")
        {
            QComboBox *editor = new QComboBox(parent);
            editor->addItem("First",0);
            editor->addItem("Second",1);
            return editor;
        }
        else if (propertyName=="Surface mesh type")
        {
            QComboBox *editor = new QComboBox(parent);
            editor->addItem("All triangles",0);
            editor->addItem("Quad dominant",1);
            return editor;
        }
        //! ---------------
        //! Scoping method
        //! ---------------
        else if (propertyName=="Scoping method")
        {
            //! -----------------------------------------------------------------
            //! non standard treatment of the dummy node for multiple selections
            //! -----------------------------------------------------------------
            DetailViewer *theDetailViewer = static_cast<DetailViewer*>(parent->parent());
            SimulationNodeClass *aNode = theDetailViewer->getCurrentMultipleSelectionNode();
            if(aNode!=Q_NULLPTR)
            {
                if(aNode->getPropertyItem("Scoping method")->data(Qt::UserRole).value<Property>().getData().isValid()==false) return Q_NULLPTR;
            }
            //! -----------------
            //! end experimental
            //! -----------------

            SimulationNodeClass *curNode = this->getCurrentNode();
            if(curNode->getPropertyItem("Scoping method")->data(Qt::UserRole).isValid()==false) return Q_NULLPTR;

            SimulationNodeClass::nodeType nodeType = curNode->getType();

            QComboBox *editor = new QComboBox(parent);
            editor->addItem("Geometry Selection",0);
            editor->addItem("Named Selection",1);

            //! -----------------------------------------------------------------------------------
            //! add the item "Remote point" only if the remote point root has been already created
            //! (it automatically contains the "dummy" remote point "Select from list"
            //! -----------------------------------------------------------------------------------
            SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
            if(sm->getTreeItem(SimulationNodeClass::nodeType_remotePointRoot)!=NULL)
            {
                if(nodeType == SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteForce ||
                        nodeType == SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteDisplacement ||
                        nodeType == SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteRotation)
                {
                    editor->addItem("Remote point",2);
                }
            }

            //! --------------------------------------------------------------
            //! add the "Automatic" scoping method if the item is a "Contact"
            //! --------------------------------------------------------------
            if(nodeType==SimulationNodeClass::nodeType_connectionPair)
            {
                editor->addItem("Automatic",2);
            }
            return editor;
        }
        //! --------------------------
        //! "Boundary scoping method"
        //! --------------------------
        else if(propertyName =="Boundary scoping method")
        {
            QComboBox *editor = new QComboBox(parent);
            editor->addItem("Geometry Selection",0);
            editor->addItem("Named Selection",1);
            return editor;
        }
        //! -------------
        //! "Suppressed"
        //! -------------
        else if (propertyName=="Suppressed")
        {
            QComboBox *editor = new QComboBox(parent);
            editor->addItem("No, active",0);
            editor->addItem("Suppressed",1);
            return editor;
        }
        else if (propertyName=="Visible")
        {
            QComboBox *editor = new QComboBox(parent);
            editor->addItem("Yes",0);
            editor->addItem("No",1);
            return editor;
        }
        //! ----------------------
        //! "Geometry" "Location"
        //! ----------------------
        else if (propertyName=="Geometry" || propertyName=="Location")
        {
            //! ----------------------------------------------------------------------
            //! the editor is created only when "Geometry" can be mofied by the user.
            //! Geometry root, a body, and a repair tool have not editable geometries
            //! by consequence no editor is returned whem the cell is clicked
            //! ----------------------------------------------------------------------
            DetailViewer *theDetailViewer = static_cast<DetailViewer*>(parent->parent());
            SimulationNodeClass *node = theDetailViewer->getNode();
            SimulationNodeClass::nodeType type = node->getType();
            if(type!=SimulationNodeClass::nodeType_geometry &&
                    type!=SimulationNodeClass::nodeType_geometryBody &&
                    type!=SimulationNodeClass::nodeType_repairTool)
            {
                ShapeSelector *editor = new ShapeSelector(myCTX,parent);
                return editor;
            }
            else return 0;
        }
        //! ----------------
        //! "Mesh entities"
        //! ----------------
        else if (propertyName=="Mesh entities")
        {
            MeshSelector *aMeshSelector = new MeshSelector(parent);
            return aMeshSelector;
        }
        //! ----------------------------------
        //! "Boundary" - for prismatic layers
        //! ----------------------------------
        else if (propertyName =="Boundary")
        {
            ShapeSelector *editor = new ShapeSelector(myCTX,parent);
            return editor;
        }
        //! ----------------------------------------
        //! "Number of steps" "Current step number"
        //! ----------------------------------------
        else if (propertyName=="Number of steps" || propertyName=="Current step number")
        {
            QSpinBox *editor = new QSpinBox(parent);
            editor->setMinimum(1);
            editor->setMaximum(1e6);
            if(propertyName =="Current step number")
            {
                cout<<"GeneralDelegate::createEditor()->____set maximum for Current step number____"<<endl;
                int maximum = (index.model()->parent(index).child(0,1)).data(Qt::UserRole).value<Property>().getData().toInt();
                editor->setMaximum(maximum);
            }
            return editor;
        }
        else if (propertyName=="Step end time")
        {
            QLineEdit *editor = new QLineEdit(parent);
            //connect(editor,SIGNAL(editingFinished()),this,SLOT(StepEndTimeLineEditClosed()));
            return editor;
        }
        else if (propertyName =="Friction coefficient" || propertyName =="C0" ||
                propertyName =="Normal stiffness" || propertyName =="Tau")
        {
            QLineEdit *editor = new QLineEdit(parent);
            QDoubleValidator *doubleValidator = new QDoubleValidator();
            doubleValidator->setDecimals(3);
            editor->setValidator(doubleValidator);
            return editor;
        }
        else if (propertyName=="Auto time stepping")
        {
            QComboBox *editor = new QComboBox(parent);
            editor->clear();
            QVariant data;
            data.setValue(Property::autoTimeStepping_ProgramControlled);
            editor->addItem("Program Controlled",data);
            data.setValue(Property::autoTimeStepping_ON);
            editor->addItem("On",data);
            data.setValue(Property::autoTimeStepping_OFF);
            editor->addItem("Off",data);
            return editor;
        }
        //! --------------
        //! "Solver type"
        //! --------------
        else if (propertyName=="Solver type")
        {
            QComboBox *editor = new QComboBox(parent);
            editor->clear();
            QVariant data;
            data.setValue(Property::solverType_programControlled);
            editor->addItem("Program Controlled",data);
            data.setValue(Property::solverType_direct);
            editor->addItem("Direct",data);
            data.setValue(Property::solverType_iterative);
            editor->addItem("Iterative",data);
            //connect(editor,SIGNAL(currentIndexChanged(int)),this, SLOT(commitAndCloseComboBox()));
            return editor;
        }
        //! -------
        //! "Name"
        //! -------
        else if(propertyName=="Name")
        {
            cout<<"GeneralDelegate::createEditor()->____creating the name editor____"<<endl;
            QLineEdit *editor = new QLineEdit(parent);
            return editor;
        }
        //! -----------------
        //! "Time step size"
        //! -----------------
        else if(propertyName=="Time step size")
        {
            QLineEdit *editor = new QLineEdit(parent);
            QDoubleValidator *validator = new QDoubleValidator();
            validator->setBottom(1e-20);
            editor->setValidator(validator);
            return editor;
        }
        //! ------------
        //! "Potential"
        //! ------------
        else if(propertyName =="Potential")
        {
            QLineEdit *editor = new QLineEdit(parent);
            QDoubleValidator *validator = new QDoubleValidator();
            editor->setValidator(validator);
            return editor;
        }
        //! ----------
        //! "Emitter"
        //! ----------
        else if(propertyName =="Emitter")
        {
            QComboBox *cb = new QComboBox(parent);
            cb->addItem("Off",0);
            cb->addItem("On",1);
            return cb;
        }
        //! ----------------------------------------------
        //! "Particle mass" "Electric charge" "Intensity"
        //! ----------------------------------------------
        else if(propertyName =="Electric charge" || propertyName =="Particle mass" || propertyName =="Intensity")
        {
            QLineEdit *editor = new QLineEdit(parent);
            QDoubleValidator *validator = new QDoubleValidator();
            editor->setValidator(validator);
            return editor;
        }
        else
        {
            return 0;
        }
    }
    else
    {
        return 0;
    }
}

//! ------------------------
//! function: setEditorData
//! details:
//! ------------------------
void GeneralDelegate::setEditorData(QWidget *editor, const QModelIndex &index) const
{
    //cout<<"GeneralDelegate::setEditorData()->____function called____"<<endl;
    QVariant data = index.model()->data(index, Qt::UserRole);
    QString propertyName = index.model()->data(index, Qt::UserRole).value<Property>().getName();

#ifdef COSTAMP_VERSION
    if(propertyName =="Activate")
    {
        // ...
    }

    //! ----------------------------------------------
    //! "Intensification pressure" && "Closure force"
    //! ----------------------------------------------
    if(propertyName =="Intensification pressure" || propertyName=="Force value")
    {
        QLineEdit *le = static_cast<QLineEdit*>(editor);
        double val = data.value<Property>().getData().toDouble();
        le->setText(QString("%1").arg(val));
    }
    if(propertyName =="Time history file")
    {
        QFileSelect *fs = static_cast<QFileSelect*>(editor);
        QString text = data.value<Property>().getData().toString();
        fs->setText(text);
    }
#endif

    //! ----------------------
    //! "Jx" "Jy" "Jz" "Mass"
    //! ----------------------
    if(propertyName =="Jx" || propertyName =="Jy" || propertyName =="Jz" || propertyName =="Mass")
    {
        QLineEdit *le =static_cast<QLineEdit*>(editor);
        double val = data.value<Property>().getData().toDouble();
        le->setText(QString("%1").arg(val));
        connect(le,SIGNAL(returnPressed()),this,SLOT(commitAndCloseLineEdit()));
    }
    //! -------------------
    //! "Static/Transient"
    //! -------------------
    if(propertyName =="Static/Transient")
    {
         QComboBox *cb = static_cast<QComboBox*>(editor);
         Property::timeIntegration ti = data.value<Property>().getData().value<Property::timeIntegration>();
         switch(ti)
         {
         case Property::timeIntegration_steadyState: cb->setCurrentIndex(0); break;
         case Property::timeIntegration_transient: cb->setCurrentIndex(1); break;
         }
         connect(cb,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseTimeIntegrationChanged()));
    }
    //! ----------------
    //! "Analysis type"
    //! ----------------
    if(propertyName =="Analysis type")
    {
         QComboBox *cb = static_cast<QComboBox*>(editor);
         Property::analysisType at = data.value<Property>().getData().value<Property::analysisType>();

         QStandardItem *curMainTreeItem = this->getCurrenItem();
         QStandardItem *curAnalysisRoot = curMainTreeItem->parent();
         SimulationNodeClass *curAnalysisNode = curAnalysisRoot->data(Qt::UserRole).value<SimulationNodeClass*>();

         switch(curAnalysisNode->getType())
         {
         case SimulationNodeClass::nodeType_structuralAnalysis:
         {
             switch(at)
             {
             case Property::analysisType_structural: cb->setCurrentIndex(0); break;
             }
         }
             break;
         case SimulationNodeClass::nodeType_thermalAnalysis:
         {
             switch(at)
             {
             case Property::analysisType_thermal: cb->setCurrentIndex(0); break;
             }
         }
             break;
         case SimulationNodeClass::nodeType_combinedAnalysis:
         {
             switch(at)
             {
             case Property::analysisType_structural: cb->setCurrentIndex(0); break;
             case Property::analysisType_thermal: cb->setCurrentIndex(1); break;
             case Property::analysisType_modal: cb->setCurrentIndex(2); break;
             case Property::analysisType_frequencyResponse: cb->setCurrentIndex(3); break;
             case Property::analysisType_uncoupledTemperatureDisplacement: cb->setCurrentIndex(4); break;
             case Property::analysisType_coupledTemperatureDisplacement: cb->setCurrentIndex(5); break;
             }
         }
             break;
         }
         connect(cb,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseAnalysisType()));
    }

    //! -----------------------
    //! "Structural time step"
    //! -----------------------
    if(propertyName == "Structural time step")
    {
        QLineEdit *le = static_cast<QLineEdit*>(editor);
        int timeStep = data.value<Property>().getData().toInt();
        le->setText(QString("%1").arg(timeStep));
        connect(le,SIGNAL(returnPressed()),this,SLOT(commitAndCloseLineEdit()));
    }
    //! ---------------------
    //! "Coupling time step"
    //! ---------------------
    if(propertyName =="Coupling time step")
    {
        ;
    }
    //! ----------------------------
    //! "Imported body temperature"
    //! ----------------------------
    if(propertyName =="Imported body temperature")
    {
        void *p = data.value<Property>().getData().value<void*>();
        QComboBox *comboBox = static_cast<QComboBox*>(editor);
        QStandardItem *item = (QStandardItem*)(p);
        if(item->data(Qt::UserRole).value<SimulationNodeClass*>()->getPropertyItem("Time tag")==Q_NULLPTR) return;

        QString timeTag = item->data(Qt::UserRole).value<SimulationNodeClass*>()->getPropertyValue<QString>("Time tag");

        QStandardItemModel *model = (QStandardItemModel*)(comboBox->model());
        //cout<<"____number of named selections: "<<model->rowCount()<<"____"<<endl;
        bool found = false;
        int n=0;
        for(; n<model->rowCount(); n++)
        {
            QModelIndex mi = model->index(n,0,QModelIndex());
            void *pp = model->data(mi,Qt::UserRole).value<void*>();
            QStandardItem *curItem = (QStandardItem*)(pp);
            SimulationNodeClass *curNode = curItem->data(Qt::UserRole).value<SimulationNodeClass*>();
            cout<<"____"<<curNode->getPropertyValue<QString>("Time tag").toStdString()<<"____"<<endl;
            if(timeTag == curNode->getPropertyValue<QString>("Time tag"))
            {
                found = true;
                break;
            }
        }
        if(found == true) comboBox->setCurrentIndex(n);
        else
        {
            return;
            //comboBox->setCurrentIndex(0);
        }
        connect(comboBox,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseComboBox()));
    }
    //! -------------------------------------------
    //! "Analysis" - for imported body temperature
    //! -------------------------------------------
    if(propertyName =="Analysis")
    {
        void *p = data.value<Property>().getData().value<void*>();
        QComboBox *comboBox = static_cast<QComboBox*>(editor);
        QStandardItem *item = (QStandardItem*)(p);

        if(item->data(Qt::UserRole).value<SimulationNodeClass*>()->getPropertyItem("Time tag")==Q_NULLPTR) return;
        QString timeTag = item->data(Qt::UserRole).value<SimulationNodeClass*>()->getPropertyValue<QString>("Time tag");

        QStandardItemModel *model = (QStandardItemModel*)(comboBox->model());

        cout<<"____number of analysis branches: "<<model->rowCount()<<"____"<<endl;

        bool found = false;
        int n=0;
        for(; n<model->rowCount(); n++)
        {
            QModelIndex mi = model->index(n,0,QModelIndex());
            void *pp = model->data(mi,Qt::UserRole).value<void*>();
            QStandardItem *curItem = (QStandardItem*)(pp);
            SimulationNodeClass *curNode = curItem->data(Qt::UserRole).value<SimulationNodeClass*>();
            //if(curNode->getPropertyItem("Time tag")==Q_NULLPTR) continue;
            cout<<"____"<<curNode->getPropertyValue<QString>("Time tag").toStdString()<<"____"<<endl;
            if(timeTag == curNode->getPropertyValue<QString>("Time tag"))
            {
                found = true;
                break;
            }
        }
        if(found == true) comboBox->setCurrentIndex(n);
        else
        {
            comboBox->setCurrentIndex(0);
        }
        connect(comboBox,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseComboBox()));
    }
    //! ---------------
    //! "Fatigue algo"
    //! ---------------
    if(propertyName =="Fatigue algo")
    {
        QComboBox *cb = static_cast<QComboBox*>(editor);
        int val = data.value<Property>().getData().toInt();
        cb->setCurrentIndex(val);
        connect(cb,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseFatigueAlgo()));
    }
    //! --------------------
    //! Material assignment
    //! --------------------
    if(propertyName =="Assignment")
    {
        QComboBox *cb = static_cast<QComboBox*>(editor);
        int val = data.value<Property>().getData().toInt();
        cb->setCurrentIndex(val);
        connect(cb,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseComboBox()));
    }
    //! ----------------------------------------------
    //! "Metric type" - currently implemented for tet
    //! 0 => "Modified Mean Ratio (MMR)"
    //! 1 => "Modified Condition Number (MCN)"
    //! 2 => "Modified volume-length (MIVL)"
    //! ----------------------------------------------
    if(propertyName == "Metric type")
    {
        QComboBox *cb = static_cast<QComboBox*>(editor);
        cb->setCurrentIndex(data.value<Property>().getData().toInt());
        connect(cb,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseMeshMetricCombobox()));
    }
    //! -----------
    //! "Grouping"
    //! -----------
    if(propertyName =="Grouping")
    {
        QComboBox *cb = static_cast<QComboBox*>(editor);
        cb->setCurrentIndex(data.value<Property>().getData().toInt());
        connect(cb,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseComboBox()));
    }
    //! ---------------
    //! "Element list"
    //! ---------------
    if(propertyName =="Element list")
    {
        QLineEdit *le = static_cast<QLineEdit*>(editor);
        QString val = data.value<Property>().getData().toString();
        le->setText(val);
        //connect(le,SIGNAL(editingFinished()),this,SLOT(commitAndCloseLineEditElementList()));
        connect(le,SIGNAL(returnPressed()),this,SLOT(commitAndCloseLineEditElementList()));
    }
    //! -------------------
    //! "Selection method"
    //! -------------------
    if(propertyName =="Selection method")
    {
        QComboBox *cb = static_cast<QComboBox*>(editor);
        int selectionMethod = data.value<Property>().getData().toInt();
        cb->setCurrentIndex(selectionMethod);
        connect(cb,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseSelectionMethod()));
    }
    //! -----------------------------------
    //! "Activation status" (model change)
    //! -----------------------------------
    if(propertyName == "Activation status")
    {
        //! -------------------------------------------------------------------------------------
        //! the compact form is simpler, but as a rule we are preferring readability when coding
        //! cb->setCurrentIndex(static_cast<int>(val)+1);
        //! -------------------------------------------------------------------------------------
        QComboBox *cb = static_cast<QComboBox*>(editor);
        Property::modelChangeActivationStatus val = data.value<Property>().getData().value<Property::modelChangeActivationStatus>();
        switch(val)
        {
        case Property::modelChangeActivationStatus_Remove: cb->setCurrentIndex(0); break;
        case Property::modelChangeActivationStatus_Inactive: cb->setCurrentIndex(1); break;
        case Property::modelChangeActivationStatus_Add: cb->setCurrentIndex(2); break;
        }
        connect(cb,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseModelChangeActivationStatus()));
    }
    //! ------------
    //! "Item type"
    //! ------------
    if(propertyName =="Item type")
    {
        QComboBox *cb = static_cast<QComboBox*>(editor);
        int val = data.value<Property>().getData().toInt();
        cb->setCurrentIndex(val);
        connect(cb,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseItemType()));
    }
    //! -------------------
    //! "Source file path"
    //! -------------------
    if(propertyName == "Source file path")
    {
        QFileSelect *fileSelect = static_cast<QFileSelect*>(editor);
        QString val = data.value<Property>().getData().toString();
        fileSelect->setText(val);
        connect(fileSelect,SIGNAL(editingFinished()),this,SLOT(commitAndCloseFileSelect()));
    }
    //! ----------------------------------------------
    //! "Solid bodies" "Surface bodies" "Line bodies"
    //! ----------------------------------------------
    if(propertyName == "Solid bodies" || propertyName == "Surface bodies" || propertyName =="Line bodies")
    {
        QComboBox *cb = static_cast<QComboBox*>(editor);
        bool val = data.value<Property>().getData().toBool();
        if(val==false) cb->setCurrentIndex(0);
        else cb->setCurrentIndex(1);
        connect(editor,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseComboBox()));
    }
    //! ------------------
    //! "Tolerance value"
    //! ------------------
    if(propertyName =="Tolerance value")
    {
        QLineEdit *le = static_cast<QLineEdit*>(editor);
        double val = data.value<Property>().getData().toDouble();
        le->setText(QString("%1").arg(val));
        connect(le,SIGNAL(editingFinished()),this,SLOT(commitAndCloseLineEdit()));
    }
    //! --------------------
    //! "Angular criterion"
    //! --------------------
    if(propertyName =="Angular criterion")
    {
        QLineEdit *le = static_cast<QLineEdit*>(editor);
        double val = data.value<Property>().getData().toDouble();
        le->setText(QString("%1").arg(val));
        connect(le,SIGNAL(editingFinished()),this,SLOT(commitAndCloseLineEdit()));
    }
    //! -------------------
    //! "Geometry healing"
    //! -------------------
    if(propertyName == "Geometry healing")
    {
        QComboBox *cb = static_cast<QComboBox*>(editor);
        bool val = data.value<Property>().getData().toBool();
        if(val==false) cb->setCurrentIndex(0);
        else cb->setCurrentIndex(1);
        connect(editor,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseGeometryHealing()));
    }
    //! -------------------------------------
    //! "Preserve boundary conditions edges"
    //! -------------------------------------
    if(propertyName =="Preserve boundary conditions edges")
    {
        QComboBox *cb = static_cast<QComboBox*>(editor);
        bool val = data.value<Property>().getData().toBool();
        if(val==false) cb->setCurrentIndex(0);
        else cb->setCurrentIndex(1);
        connect(editor,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseComboBox()));
    }
    //! -----------------------------
    //! "Project points on geometry"
    //! -----------------------------
    if(propertyName =="Project points on geometry")
    {
        QComboBox *cb = static_cast<QComboBox*>(editor);
        bool val = data.value<Property>().getData().toBool();
        if(val==false) cb->setCurrentIndex(0);
        else cb->setCurrentIndex(1);
        connect(editor,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseComboBox()));
    }
    //! -----------------
    //! "Simplification"
    //! -----------------
    if(propertyName=="Simplification")
    {
        QComboBox *cb = static_cast<QComboBox*>(editor);
        bool val = data.value<Property>().getData().toBool();
        if(val==false) cb->setCurrentIndex(0);
        else cb->setCurrentIndex(1);
        connect(editor,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseSimplificationControl()));
    }
    //! ----------
    //! "Healing"
    //! ----------
    if(propertyName=="Healing")
    {
        QComboBox *cb = static_cast<QComboBox*>(editor);
        bool val = data.value<Property>().getData().toBool();
        if(val==false) cb->setCurrentIndex(0);
        else cb->setCurrentIndex(1);
        connect(editor,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseComboBox()));
    }
    //! --------------
    //! "Defeaturing"
    //! --------------
    if(propertyName=="Defeaturing")
    {
        QComboBox *cb = static_cast<QComboBox*>(editor);
        bool val = data.value<Property>().getData().toBool();
        if(val==false) cb->setCurrentIndex(0);
        else cb->setCurrentIndex(1);
        connect(editor,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseDefeaturingControl()));
    }
    //! ----------------------
    //! "Tessellator"
    //! ----------------------
    if(propertyName=="Tessellator")
    {
        QComboBox *cb = static_cast<QComboBox*>(editor);
        Property::meshEngine2D val = data.value<Property>().getData().value<Property::meshEngine2D>();
        switch(val)
        {
        case Property::meshEngine2D_OCC_STL: cb->setCurrentIndex(0); break;
        case Property::meshEngine2D_OCC_ExpressMesh: cb->setCurrentIndex(1); break;
        }
        //connect(cb,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseComboBox()));
        connect(cb,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseTessellator()));
    }
    //! -----------------------------------------
    //! "Angular deflection" "Linear deflection"
    //! -----------------------------------------
    if(propertyName =="Angular deflection" || propertyName =="Linear deflection")
    {
        QLineEdit *le = static_cast<QLineEdit*>(editor);
        double val = data.value<Property>().getData().toDouble();
        le->setText(QString("%1").arg(val));
    }
    //! --------------------------------------------------------------------
    //! "Min face size" "Max face size" - force values to follow some rules
    //! --------------------------------------------------------------------
    if(propertyName =="Min face size" || propertyName =="Max face size")
    {
        QLineEdit *le = static_cast<QLineEdit*>(editor);
        double val = data.value<Property>().getData().toDouble();
        if(propertyName =="Min face size")
        {
            double curMaxFaceSize = this->getCurrentNode()->getPropertyValue<double>("Max face size");
            if(val>=curMaxFaceSize) val = curMaxFaceSize;
        }
        if(propertyName =="Max face size")
        {
            double curMinFaceSize = this->getCurrentNode()->getPropertyValue<double>("Min face size");
            if(val<=curMinFaceSize) val = curMinFaceSize;
        }
        le->setText(QString("%1").arg(val));
    }
    //! -------------------
    //! "Patch conforming"
    //! -------------------
    if(propertyName =="Patch conforming")
    {
        QComboBox *cb = static_cast<QComboBox*>(editor);
        bool val = data.value<Property>().getData().toBool();
        if(val==false) cb->setCurrentIndex(0);
        else cb->setCurrentIndex(1);
        connect(editor,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseUseBRep()));
    }
    //! -----------------
    //! "Surface mesher"
    //! -----------------
    if(propertyName =="Surface mesher")
    {
        QComboBox *cb = static_cast<QComboBox*>(editor);
        Property::meshEngine2D val = data.value<Property>().getData().value<Property::meshEngine2D>();
        bool useBRep = this->getCurrentNode()->getPropertyValue<bool>("Patch conforming");
        if(useBRep)
        {
            switch(val)
            {
            case Property::meshEngine2D_Netgen: cb->setCurrentIndex(0); break;
            case Property::meshEngine2D_OCC_ExpressMesh: cb->setCurrentIndex(1); break;
            }
        }
        else
        {
            switch(val)
            {
            case Property::meshEngine2D_Netgen_STL: cb->setCurrentIndex(0); break;
            //case Property::meshEngine2D_TetgenBR: cb->setCurrentIndex(1); break;
            }
        }
    }
    //! ----------------
    //! "Volume mesher"
    //! ----------------
    if(propertyName =="Volume mesher")
    {
        QComboBox *cb = static_cast<QComboBox*>(editor);
        Property::meshEngine3D val = data.value<Property>().getData().value<Property::meshEngine3D>();
        bool useBRep = this->getCurrentNode()->getPropertyValue<bool>("Patch conforming");

        if(useBRep)
        {
            switch(val)
            {
            case Property::meshEngine3D_Netgen: cb->setCurrentIndex(0); break;
            case Property::meshEngine3D_Tetgen: cb->setCurrentIndex(1); break;
            }
        }
        else    // "Patch conforming" = false
        {
            if(this->getCurrentNode()->getPropertyItem("Surface mesher")!=NULL)
            {
                switch(val)
                {
                case Property::meshEngine3D_Netgen_STL: cb->setCurrentIndex(0); break;
                case Property::meshEngine3D_Tetgen: cb->setCurrentIndex(1); break;
                //case Property::meshEngine3D_Tetgen_BR: cb->setCurrentIndex(2); break;
                //case Property::meshEngine3D_TetWild: cb->setCurrentIndex(3); break;
                }
            }
            else
            {
                //if(val==Property::meshEngine3D_Tetgen_BR) cb->setCurrentIndex(2);
                switch(val)
                {
                case Property::meshEngine3D_Tetgen_BR: cb->setCurrentIndex(2); break;
                case Property::meshEngine3D_TetWild: cb->setCurrentIndex(3); break;
                }
            }
        }
        connect(editor,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseVolumeMesher()));
    }
    //! ------------------
    //! "Envelope sizing"
    //! ------------------
    if(propertyName =="Envelope sizing")
    {
        QComboBox *cb = static_cast<QComboBox*>(editor);
        int val = data.value<Property>().getData().toInt();
        if(val ==0) cb->setCurrentIndex(0);
        else cb->setCurrentIndex(1);
        connect(cb,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseEnvelopeSizing()));
    }
    //! ----------------------
    //! "Ideal length sizing"
    //! ----------------------
    if(propertyName == "Ideal length sizing")
    {
        QComboBox *cb = static_cast<QComboBox*>(editor);
        int val = data.value<Property>().getData().toInt();
        if(val ==0) cb->setCurrentIndex(0);
        else cb->setCurrentIndex(1);
        connect(cb,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseIdealLength()));
    }
    //! --------------------------------------------------------------------
    //! "Relative size" "Absolute size" "Relative length" "Absolute length"
    //! --------------------------------------------------------------------
    if(propertyName =="Relative size" || propertyName =="Absolute size" ||
            propertyName =="Relative length" || propertyName =="Absolute length")
    {
        QLineEdit *le = static_cast<QLineEdit*>(editor);
        double val = data.value<Property>().getData().toDouble();
        le->setText(QString("%1").arg(val));

        //! the slot does not emit signals (it only closes the editor)
        connect(le,SIGNAL(editingFinished()),this,SLOT(commitAndCloseLineEdit()));
    }
    //! --------
    //! Mapping
    //! --------
    if(propertyName=="Mapping")
    {
        QComboBox *cb = static_cast<QComboBox*>(editor);
        int val = data.value<Property>().getData().toInt();
        cb->setCurrentIndex(val);
        connect(cb,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseComboBox()));
    }
    //! ---------------------------------------------
    //! stress strain/strain source for fatigue tool
    //! ---------------------------------------------
    if(propertyName =="Stress/strain source")
    {
        QComboBox *cb = static_cast<QComboBox*>(editor);
        int val = data.value<Property>().getData().toInt();
        cb->setCurrentIndex(val);
        connect(cb,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseComboBox()));
    }
    //! ------------------------------------------------
    //! stress strain/strain component for fatigue tool
    //! ------------------------------------------------
    if(propertyName =="Component")
    {
        QComboBox *cb = static_cast<QComboBox*>(editor);
        int val = data.value<Property>().getData().toInt();
        cb->setCurrentIndex(val);
        connect(cb,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseComboBox()));
    }
    //! ----------------------------------------
    //! triaxiality correction for fatigue tool
    //! ----------------------------------------
    if(propertyName =="Triaxiality correction")
    {
        QComboBox *cb = static_cast<QComboBox*>(editor);
        int val = data.value<Property>().getData().toInt();
        cb->setCurrentIndex(val);
        connect(cb,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseComboBox()));
    }
    //! --------------------------------------
    //! "Number of cycles" - for fatigue tool
    //! --------------------------------------
    if(propertyName =="Number of cycles")
    {
        QLineEdit *le = static_cast<QLineEdit*>(editor);
        int val = data.value<Property>().getData().toInt();
        le->setText(QString("%1").arg(val));
        connect(le,SIGNAL(returnPressed()),this,SLOT(commitAndCloseLineEdit()));
    }
    //! ----------------
    //! "Lock boundary"
    //! ----------------
    if(propertyName =="Lock boundary")
    {
        QComboBox *cb = static_cast<QComboBox*>(editor);
        bool val = data.value<Property>().getData().toBool();
        if(val==true) cb->setCurrentIndex(0);
        else cb->setCurrentIndex(1);
    }
    //! --------------------------------------------------------------
    //! "Guiding vectors smoothing steps" "Thickness smoothing steps"
    //! --------------------------------------------------------------
    if(propertyName=="Guiding vectors smoothing steps" || propertyName=="Thickness smoothing steps")
    {
        QLineEdit *le = static_cast<QLineEdit*>(editor);
        int val = data.value<Property>().getData().toInt();
        le->setText(QString("%1").arg(val));
        connect(le,SIGNAL(returnPressed()),this,SLOT(commitAndCloseLineEdit()));
    }
    //!  --------------------------------------------------------------------------
    //! "Curvature sensitivity" "Guiding vector smoothing - curvature sensitivity"
    //! "Thickness smoothing - curvature sensitivity"
    //!  --------------------------------------------------------------------------
    if(propertyName =="Curvature sensitivity" ||
            propertyName =="Guiding vector smoothing - curvature sensitivity" ||
            propertyName =="Thickness smoothing - curvature sensitivity")
    {
        QLineEdit *le = static_cast<QLineEdit*>(editor);
        double val = data.value<Property>().getData().toDouble();
        le->setText(QString("%1").arg(val));
        connect(le,SIGNAL(returnPressed()),this,SLOT(commitAndCloseLineEdit()));
    }
    //! ---------------------
    //! "Volume mesh engine"
    //! ---------------------
    if(propertyName =="Volume mesh engine")
    {
        QComboBox *cb = static_cast<QComboBox*>(editor);
        Property::meshEngine3D volumeMeshEngine = data.value<Property>().getData().value<Property::meshEngine3D>();
        switch(volumeMeshEngine)
        {
        case Property::meshEngine3D_Netgen: cb->setCurrentIndex(0); break;
        case Property::meshEngine3D_Tetgen: cb->setCurrentIndex(1); break;
        }
    }
    //! ---------------------------------------------
    //! "Modulation coefficient transfer percentage"
    //! ---------------------------------------------
    if(propertyName =="Modulation coefficient transfer percentage")
    {
        QLineEdit *le = static_cast<QLineEdit*>(editor);
        double val = data.value<Property>().getData().toDouble();
        le->setText(QString("%1").arg(val));
    }
    //! -----------------------------
    //! "Modulation diffusion steps"
    //! -----------------------------
    if(propertyName =="Modulation diffusion cutoff")
    {
        QLineEdit *le = static_cast<QLineEdit*>(editor);
        double val = data.value<Property>().getData().toDouble();
        le->setText(QString("%1").arg(val));
    }
    //! -----------------------------
    //! "Modulation diffusion steps"
    //! -----------------------------
    if(propertyName =="Modulation diffusion steps")
    {
        QLineEdit *le = static_cast<QLineEdit*>(editor);
        int val = data.value<Property>().getData().toInt();
        le->setText(QString("%1").arg(val));
    }
    //! ------------------------------------------
    //! "Transition" "Minimum shrink" "Amplitude"
    //! ------------------------------------------
    if(propertyName =="Transition" || propertyName =="Minimum shrink" || propertyName =="Amplitude")
    {
        QLineEdit *le = static_cast<QLineEdit*>(editor);
        double val = data.value<Property>().getData().toDouble();
        le->setText(QString("%1").arg(val));
    }
    //! -------------------------------------------
    //! "Options" (for prismatic layer generation)
    //! 0 => "First layer thickness"
    //! 1 => "Total thickness"
    //! -------------------------------------------
    if(propertyName =="Options")
    {
        QComboBox *cb = static_cast<QComboBox*>(editor);
        int val = data.value<Property>().getData().toInt();
        switch(val)
        {
        case 0: cb->setCurrentIndex(0); break;
        case 1: cb->setCurrentIndex(1); break;
        }
        connect(cb,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndClosePrismaticLayerOptions()));
    }
    //! -------------------
    //! "Number of layers"
    //! -------------------
    if(propertyName =="Number of layers")
    {
        QLineEdit *le = static_cast<QLineEdit*>(editor);
        int val = data.value<Property>().getData().toInt();
        le->setText(QString("%1").arg(val));
    }
    //! -----------------------------------------------------
    //! "Expansion ratio" "Total thickness" "First layer height"
    //! -----------------------------------------------------
    if(propertyName =="Expansion ratio" || propertyName =="Total thickness" || propertyName =="First layer height")
    {
        QLineEdit *le = static_cast<QLineEdit*>(editor);
        double val = data.value<Property>().getData().toDouble();
        le->setText(QString("%1").arg(val));
    }
    //! -------------------------------------------------------
    //! "Check self intersections" "Check mutual intersections
    //! -------------------------------------------------------
    if(propertyName =="Check self intersections" || propertyName =="Check mutual intersections")
    {
        QComboBox *cb = static_cast<QComboBox*>(editor);
        bool val = data.value<Property>().getData().toBool();
        if(val==true) cb->setCurrentIndex(1);
        else cb->setCurrentIndex(0);
        connect(cb,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseComboBox()));
    }
    //! ----------------
    //! "Run in memory"
    //! ----------------
    if(propertyName =="Run in memory")
    {
        QComboBox *cb = static_cast<QComboBox*>(editor);
        bool val = data.value<Property>().getData().toBool();
        if(val==true) cb->setCurrentIndex(0);
        else cb->setCurrentIndex(1);
        connect(cb,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseComboBox()));
    }
    //! ----------------
    //! "Pair distance"
    //! ----------------
    if(propertyName =="Pair distance")
    {
        QLineEdit *le = static_cast<QLineEdit*>(editor);
        double val = data.value<Property>().getData().toDouble();
        le->setText(QString("%1").arg(val));
        connect(le,SIGNAL(editingFinished()),this,SLOT(commitAndCloseMeshDefeaturingParameterValueControl()));
    }
    //! --------
    //! "Level"
    //! --------
    if(propertyName =="Level")
    {
        double level = data.value<Property>().getData().toDouble();
        int value = (level-MIN_MESH_REDUCTION)/(MAX_MESH_REDUCTION-MIN_MESH_REDUCTION)*100.0;
        QSlider *slider = static_cast<QSlider*>(editor);
        slider->setValue(value);
        //! the signal "sliderReleased()" is not OK. Please, change it ... to do
        connect(slider,SIGNAL(sliderReleased()),this,SLOT(commitAndCloseMeshDefeaturingParameterValueControl()));
    }
    //! ----------------
    //! New mesh method
    //! ----------------
    if(propertyName =="Method")
    {
        QComboBox *cb = static_cast<QComboBox*>(editor);
        int val = data.value<Property>().getData().toInt();
        cb->setCurrentIndex(val);
        connect(cb,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseMeshMethodSelector()));
    }
    //! ------------
    //! "Min" "Max"
    //! ------------
    if(propertyName =="Min" || propertyName =="Max")
    {
        QLineEdit *le = static_cast<QLineEdit*>(editor);
        double val = data.value<Property>().getData().toDouble();
        le->setText(QString("%1").arg(val));
        connect(le,SIGNAL(editingFinished()),this,SLOT(commitAndCloseMinMaxControls()));
    }
    //! -------------
    //! "# intervals
    //! -------------
    if(propertyName =="# intervals")
    {
        QSpinBox *spinBox = static_cast<QSpinBox*>(editor);
        int val = data.value<Property>().getData().toInt();
        spinBox->setValue(val);
        connect(editor,SIGNAL(editingFinished()),this,SLOT(commitAndCloseNumberOfInterval()));
    }
    //! -------------
    //! "Scale type"
    //! -------------
    else if(propertyName =="Scale type")
    {
        QComboBox *cb = static_cast<QComboBox*>(editor);
        int val = data.value<Property>().getData().toInt();
        cb->setCurrentIndex(val);
        connect(cb,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseScaleTypeSelector()));
    }
    //! --------------------
    //! "--Recurrence rate"
    //! --------------------
    else if(propertyName =="--Recurrence rate")
    {
        QLineEdit *le = static_cast<QLineEdit*>(editor);
        int value = data.value<Property>().getData().toInt();
        le->setText(QString("%1").arg(value));
    }
    //! ---------------------------------------
    //! Output settings for "Store results at"
    //! ---------------------------------------
    else if(propertyName =="Store results at")
    {
        QComboBox *cb = static_cast<QComboBox*>(editor);
        int value = data.value<Property>().getData().toInt();
        switch(value)
        {
        case 0: cb->setCurrentIndex(0); break;
        case 1: cb->setCurrentIndex(1); break;
        case 2: cb->setCurrentIndex(2); break;
        }
        connect(cb,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseStoreResultsAt()));
    }
    //! -----------------------------------------------------------------------
    //! Output settings for "Stress" "Strain" "Reaction forces" "Contact data"
    //! -----------------------------------------------------------------------
    else if(propertyName =="Stress" || propertyName=="Strain" || propertyName =="Reaction forces" || propertyName=="Contact data")
    {
        QComboBox *cb = static_cast<QComboBox*>(editor);
        int value = data.value<Property>().getData().toInt();
        switch(value)
        {
        case 0: cb->setCurrentIndex(0); break;
        case 1: cb->setCurrentIndex(1); break;
        }
        connect(cb,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseComboBoxForOutputSettings()));
    }
    //! ----------------------------------------
    //! Advanced controls for the contact group
    //! ----------------------------------------
    else if(propertyName =="K" || propertyName=="Sigma infty" || propertyName =="C0" ||
            propertyName =="Lambda" || propertyName =="P0")
    {
        //SimulationNodeClass *curNode = this->getCurrentNode();
        //Property::contactBehavior val = curNode->getPropertyItem("Behavior")->data(Qt::UserRole).value<Property>().getData().value<Property::contactBehavior>();
        //if(val == Property::contactBehavior_symmetric && propertyName != "C0" && propertyName == "Sigma infty")
        //{
            double value = data.value<Property>().getData().toDouble();
            QLineEdit *le = static_cast<QLineEdit*>(editor);
            le->setText(QString("%1").arg(value));
        //}
    }
    //! ----------------
    //! "Small sliding"
    //! ----------------
    else if(propertyName =="Small sliding")
    {
        SimulationNodeClass *curNode = this->getCurrentNode();
        Property::contactBehavior val = curNode->getPropertyItem("Behavior")->data(Qt::UserRole).value<Property>().getData().value<Property::contactBehavior>();
        if(val == Property::contactBehavior_asymmetric)
        {
            int value = data.value<Property>().getData().toInt();
            QComboBox *cb = static_cast<QComboBox*>(editor);
            cb->setCurrentIndex(value);
            connect(cb,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseSmallSlidingControl()));
        }
    }

    //! ----------------
    //! "Adjust to touch"
    //! ----------------
    else if(propertyName =="Adjust to touch")
    {
        SimulationNodeClass *curNode = this->getCurrentNode();
        Property::contactBehavior val = curNode->getPropertyItem("Behavior")->data(Qt::UserRole).value<Property>().getData().value<Property::contactBehavior>();
        Property::contactType type = curNode->getPropertyItem("Type")->data(Qt::UserRole).value<Property>().getData().value<Property::contactType>();
        if(val == Property::contactBehavior_asymmetric || type == Property::contactType_bonded)
        {
            int value = data.value<Property>().getData().toInt();
            QComboBox *cb = static_cast<QComboBox*>(editor);
            cb->setCurrentIndex(value);
            connect(cb,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseAdjustControl()));
        }
    }

    //! --------------
    //! "Line search"
    //! --------------
    if(propertyName == "Line search")
    {
        int value = data.value<Property>().getData().toInt();
        QComboBox *cb = static_cast<QComboBox*>(editor);
        cb->setCurrentIndex(value);
        connect(cb,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseLineSearchEditor()));
    }
    //! -------------------------
    //! "Line search parameters"
    else //! -------------------------
    if(propertyName == "Min value" || propertyName == "Max value")
    {
        int lineSearch = this->getCurrentNode()->getPropertyItem("Line search")->data(Qt::UserRole).value<Property>().getData().toInt();
        if(lineSearch!=0)
        {
            double val = data.value<Property>().getData().toDouble();
            QLineEdit *le = static_cast<QLineEdit*>(editor);
            le->setText(QString("%1").arg(val));
            connect(le,SIGNAL(editingFinished()),this,SLOT(commitAndCloseLineParametersChanged()));
        }
    }
    //! ------------------
    //! "Cutback factors"
    //! ------------------
    else if(propertyName =="Cutback factors")
    {
        int value = data.value<Property>().getData().toInt();
        QComboBox *cb = static_cast<QComboBox*>(editor);
        cb->setCurrentIndex(value);
        connect(cb,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseCutBackEditor()));
    }
    //! ---------------------
    //! "Cutback parameters"
    //! ---------------------
    else if(propertyName =="D_f" || propertyName =="D_C" || propertyName =="D_B" || propertyName =="D_A"
            || propertyName =="D_S" || propertyName =="D_H" || propertyName =="D_D" || propertyName =="W_G")
    {
        int cutBackType = this->getCurrentNode()->getPropertyItem("Cutback factors")->data(Qt::UserRole).value<Property>().getData().toInt();
        if(cutBackType==1)
        {
            if(propertyName!="D_S" && propertyName!="D_H" && propertyName!="W_G")
            {
                double value = data.value<Property>().getData().toDouble();
                QLineEdit *le = static_cast<QLineEdit*>(editor);
                le->setText(QString("%1").arg(value));
                connect(le,SIGNAL(editingFinished()),this,SLOT(commitAndCloseCutBackParametersLineEdit()));
            }
        }
    }
    //! ----------------------
    //! "Time incrementation"
    //! ----------------------
    else if(propertyName =="Time incrementation")
    {
        int value = data.value<Property>().getData().toInt();
        QComboBox *cb = static_cast<QComboBox*>(editor);
        cb->setCurrentIndex(value);
        connect(cb,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseTimeIncrementationControl()));
    }
    //! ---------------------------------
    //! "Time incrementation" parameters
    //! ---------------------------------
    else if(propertyName =="I_0" || propertyName =="I_R" || propertyName =="I_P" || propertyName =="I_C" || propertyName =="I_L"
          || propertyName =="I_G" || propertyName =="I_S" || propertyName =="I_A" || propertyName =="I_J" || propertyName =="I_T")
    {
        int timeIncrementationType = this->getCurrentNode()->getPropertyItem("Time incrementation")->data(Qt::UserRole).value<Property>().getData().toInt();
        if(timeIncrementationType==1)
        {
            if(propertyName !="I_S" && propertyName !="I_J" && propertyName!="I_T")
            {
                int value = data.value<Property>().getData().toInt();
                QLineEdit *le = static_cast<QLineEdit*>(editor);
                le->setText(QString("%1").arg(value));
                connect(le,SIGNAL(editingFinished()),this,SLOT(commitAndCloseTimeIncrementationEditor()));
            }
        }
    }
    //! -------------------
    //! "Flux convergence"
    //! -------------------
    else if(propertyName =="Flux convergence")
    {
        int value = data.value<Property>().getData().toInt();
        QComboBox *cb = static_cast<QComboBox*>(editor);
        cb->setCurrentIndex(value);
        connect(cb,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseComboBoxFluxConvergence()));
    }
    //! ----------------------------------------------------
    //! --R_alpha_n --R_alpha_P __R_alpha_l --epsilon_alpha
    //! ----------------------------------------------------
    else if(propertyName =="--R_alpha_n" || propertyName =="--R_alpha_P" || propertyName == "--R_alpha_l" || propertyName=="--epsilon_alpha")
    {
        int fluxConvergenceFlag = this->getCurrentNode()->getPropertyItem("Flux convergence")->data(Qt::UserRole).value<Property>().getData().toInt();
        if(fluxConvergenceFlag==1)
        {
            double val = data.value<Property>().getData().toDouble();
            QLineEdit *le = static_cast<QLineEdit*>(editor);
            le->setText(QString("%1").arg(val));
            connect(le,SIGNAL(editingFinished()),this,SLOT(commitAndCloseFieldParametersLineEditEditor()));
        }
    }
    //! ------------------------------
    //! --C_alpha_n --C_alpha_epsilon
    //! ------------------------------
    else if(propertyName == "--C_alpha_n" || propertyName == "--C_alpha_epsilon")
    {
        int solutionConvergenceFlag = this->getCurrentNode()->getPropertyItem("Solution convergence")->data(Qt::UserRole).value<Property>().getData().toInt();
        if(solutionConvergenceFlag==1)
        {
            double val = data.value<Property>().getData().toDouble();
            QLineEdit *le = static_cast<QLineEdit*>(editor);
            le->setText(QString("%1").arg(val));
            connect(le,SIGNAL(editingFinished()),this,SLOT(commitAndCloseFieldParametersLineEditEditor()));
        }
    }
    //! -------------------------
    //! "--Value" or --q_alpha_u
    //! -------------------------
    else if(propertyName=="--Value")
    {
        int fluxConvergenceFlag = this->getCurrentNode()->getPropertyItem("Flux convergence")->data(Qt::UserRole).value<Property>().getData().toInt();
        if(fluxConvergenceFlag==1)
        {
            double val = data.value<Property>().getData().toDouble();
            QLineEdit *le = static_cast<QLineEdit*>(editor);
            if(val=0.0) le->setText(QString("Calculated by solver"));
            else le->setText(QString("%1").arg(val));
            connect(le,SIGNAL(editingFinished()),this,SLOT(commitAndCloseFieldParametersLineEditEditor()));
        }
    }
    //! ------------
    //! --q_alpha_0
    //! ------------
    else if(propertyName=="--q_alpha_0")
    {
        int fluxConvergenceFlag = this->getCurrentNode()->getPropertyItem("Flux convergence")->data(Qt::UserRole).value<Property>().getData().toInt();
        if(fluxConvergenceFlag==1)
        {
            double val = data.value<Property>().getData().toDouble();
            QLineEdit *le = static_cast<QLineEdit*>(editor);
            if(val=0.0) le->setText(QString("Calculated by solver"));
            else le->setText(QString("%1").arg(val));
            connect(le,SIGNAL(editingFinished()),this,SLOT(commitAndCloseFieldParametersLineEditEditor()));
        }
    }
    //! ------------------------
    //! --R_alpha_n --R_alpha_P
    //! ------------------------
    else if(propertyName =="--R_alpha_n" || propertyName =="--R_alpha_P")
    {
        int fluxConvergenceFlag = this->getCurrentNode()->getPropertyItem("Flux convergence")->data(Qt::UserRole).value<Property>().getData().toInt();
        if(fluxConvergenceFlag==1)
        {
            double val = data.value<Property>().getData().toDouble();
            QLineEdit *le = static_cast<QLineEdit*>(editor);
            le->setText(QString("%1").arg(val));
            connect(le,SIGNAL(editingFinished()),this,SLOT(commitAndCloseFieldParametersLineEditEditor()));
        }
    }
    //! -----------------------
    //! "Solution convergence"
    //! -----------------------
    else if(propertyName=="Solution convergence")
    {
        int value = data.value<Property>().getData().toInt();
        QComboBox *cb = static_cast<QComboBox*>(editor);
        cb->setCurrentIndex(value);
        connect(cb,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseComboBoxSolutionConvergence()));
    }
    //! -------------------
    //! "Large deflection"
    //! -------------------
    else if(propertyName=="Large deflection")
    {
        int value = data.value<Property>().getData().toInt();
        QComboBox *cb = static_cast<QComboBox*>(editor);
        cb->setCurrentIndex(value);
        connect(cb,SIGNAL(currentIndexChanged(int)),this,SLOT(closeComboBox()));
    }
    //! ---------------
    //! "Display time"
    //! ---------------
    else if(propertyName =="Display time")
    {
        double value = data.value<Property>().getData().toDouble();
        QLineEdit *le = static_cast<QLineEdit*>(editor);
        le->setText(QString("%1").arg(value));
    }
    //! ---------------
    //! "Mode number"
    //! ---------------
    else if(propertyName =="Mode number")
    {
        int value = data.value<Property>().getData().toInt();
        QLineEdit *le = static_cast<QLineEdit*>(editor);
        le->setText(QString("%1").arg(value));
    }
    //! ------------------
    //! "Update interval"
    //! ------------------
    else if(propertyName =="Update interval")
    {
        double value = data.value<Property>().getData().toDouble();
        QLineEdit *le = static_cast<QLineEdit*>(editor);
        le->setText(QString("%1").arg(value));
        connect(le,SIGNAL(editingFinished()),this,SLOT(commitAndCloseUpdateInterval()));
    }
    //! ------------
    //! "Tolerance"
    //! ------------
    else if(propertyName =="Tolerance")
    {
        SimulationNodeClass::nodeType nodeType = this->getCurrentNode()->getType();
        if(nodeType==SimulationNodeClass::nodeType_import)
        {
            QComboBox *cb = static_cast<QComboBox*>(editor);
            int toleranceType = data.value<Property>().getData().toInt();
            cb->setCurrentIndex(toleranceType);
            connect(cb,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseToleranceType()));
        }
        else
        {
            double value = data.value<Property>().getData().toDouble();
            QLineEdit *le = static_cast<QLineEdit*>(editor);
            le->setText(QString("%1").arg(value));
            connect(editor,SIGNAL(editingFinished()),this, SLOT(commitAndCloseContactToleranceEditor()));
        }
    }
    //! -------------
    //! "Set number"
    //! -------------
    else if(propertyName =="Set number")
    {
        int value = data.value<Property>().getData().toInt();
        QLineEdit *le = static_cast<QLineEdit*>(editor);
        le->setText(QString("%1").arg(value));
    }
    //! -----
    //! "By"
    //! -----
    else if(propertyName == "By")
    {
        switch(this->getCurrentNode()->getType())
        {
        case SimulationNodeClass::nodeType_solutionThermalFlux:
        case SimulationNodeClass::nodeType_solutionThermalTemperature:
        case SimulationNodeClass::nodeType_solutionStructuralMechanicalStrain:
        case SimulationNodeClass::nodeType_solutionStructuralNodalDisplacement:
        case SimulationNodeClass::nodeType_solutionStructuralStress:
        case SimulationNodeClass::nodeType_solutionStructuralTemperature:
        case SimulationNodeClass::nodeType_solutionStructuralThermalStrain:
        case SimulationNodeClass::nodeType_solutionStructuralTotalStrain:
        case SimulationNodeClass::nodeType_solutionStructuralEquivalentPlasticStrain:
        case SimulationNodeClass::nodeType_solutionStructuralNodalForces:
        {
            int val = data.value<Property>().getData().toInt();
            QComboBox *cb = static_cast<QComboBox*>(editor);
            cb->setCurrentIndex(val);
            connect(cb,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseBySelector()));
        }
            break;

        case SimulationNodeClass::nodeType_meshMethod:
        {
            int val = data.value<Property>().getData().toInt();
            QComboBox *cb = static_cast<QComboBox*>(editor);
            cb->setCurrentIndex(val);
            connect(cb,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseBySelector()));
        }
            break;
        }
    }
    //! --------
    //! "Type "
    //! --------
    else if(propertyName == "Type ")
    {
        int val = data.value<Property>().getData().toInt();
        QComboBox *cb = static_cast<QComboBox*>(editor);
        cb->setCurrentIndex(val);
        connect(cb,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseSolutionComponentSelector()));
    }

    //! -----------------------
    //! "Solution information"
    //! -----------------------
    else if(propertyName =="Solution information")
    {
        //! -----------------------------------------------
        //! 0 => solver txt messages
        //! 1 => force convergence
        //! 2 => displacement convergence
        //! 3 => line search
        //! 4 => time step size
        //! -----------------------------------------------
        int val = data.value<Property>().getData().toInt();
        QComboBox *cb = static_cast<QComboBox*>(editor);
        cb->setCurrentIndex(val);
        connect(cb,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseSolutionInformation()));
    }
    //! ---------
    //! Coupling
    //! ---------
    else if(propertyName=="Coupling")
    {
        int val = data.value<Property>().getData().toInt();
        QComboBox *cb = static_cast<QComboBox*>(editor);
        cb->setCurrentIndex(val);
        connect(cb,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseCouplingSelector()));
    }
    //! ---------------
    //! DOFs selection
    //! ---------------
    else if(propertyName=="DOFs selection")
    {
        int val = data.value<Property>().getData().toInt();
        QComboBox *cb = static_cast<QComboBox*>(editor);
        cb->setCurrentIndex(val);
        connect(cb,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseDOFswitchSelector()));
    }
    //! ------------------------------------------------
    //! "X component" "Y component " "Z component " ...
    //! ------------------------------------------------
    else if(propertyName == "X component " || propertyName == "Y component " || propertyName == "Z component " ||
            propertyName == "X rotation" || propertyName == "Y rotation" || propertyName == "Z rotation")
    {
        int val = data.value<Property>().getData().toInt();
        QComboBox *cb = static_cast<QComboBox*>(editor);
        cb->setCurrentIndex(val);
        connect(cb,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseRPDOFselector()));
    }
    //! ------------------
    //! number of threads
    //! ------------------
    else if(propertyName =="Number of threads")
    {
        //cerr<<"____set editor data for \"Number of threads\"____"<<endl;
        //if(data.canConvert<Property>()) cerr<<"____can convert data to \"Property\"____"<<endl;
        //if(data.value<Property>().getData().canConvert<int>()) cerr<<"____can convert property value to \"int\"____"<<endl;
        int val = data.value<Property>().getData().toInt();
        QSpinBox *spinBox = static_cast<QSpinBox*>(editor);
        spinBox->setValue(val);
        disconnect(spinBox,SIGNAL(valueChanged(int)),this,SLOT(commitAndCloseNbThreads()));
        connect(spinBox,SIGNAL(valueChanged(int)),this,SLOT(commitAndCloseNbThreads()));
    }
    //! -----------------------------
    //! "First color" "Second color"
    //! -----------------------------
    else if(propertyName == "First color" || propertyName =="Second color")
    {
        QVector<int> rgb = data.value<Property>().getData().value<QVector<int>>();
        QColor color(rgb.at(0),rgb.at(1),rgb.at(2));
        color.setRed(rgb.at(0));
        color.setGreen(rgb.at(1));
        color.setBlue(rgb.at(2));
        ColorSelector *cs = static_cast<ColorSelector*>(editor);
        cs->setColor(color);
        connect(editor,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseBackgroundControl()));
    }
    //! -----------
    //! "Gradient"
    //! -----------
    else if(propertyName =="Gradient")
    {
        int value = data.value<Property>().getData().toInt();
        QComboBox *comboBox = static_cast<QComboBox*>(editor);
        comboBox->setCurrentIndex(value);
        connect(editor,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseBackgroundControl()));
    }
    //! ----------------
    //! "Analysis time"
    //! ----------------
    else if(propertyName =="Analysis time")
    {
        double value = data.value<Property>().getData().toDouble();
        QLineEdit *le = static_cast<QLineEdit*>(editor);
        le->setText(QString("%1").arg(value));
    }
    //! ----------------------------
    //! "Buckets" "Remapping steps"
    //! ----------------------------
    else if(propertyName =="X buckets" ||
            propertyName =="Y buckets" ||
            propertyName =="Z buckets" ||
            propertyName =="Remapping steps")
    {
        int value = data.value<Property>().getData().toInt();
        QLineEdit *le = static_cast<QLineEdit*>(editor);
        le->setText(QString("%1").arg(value));
    }
    //! -------------
    //! "Split data"
    //! -------------
    else if(propertyName =="Split data")
    {
        int value = data.value<Property>().getData().toInt();
        QComboBox *comboBox = static_cast<QComboBox*>(editor);
        comboBox->setCurrentIndex(value);
        connect(editor,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseSplitFileModeSelector()));
    }
    //! ------------
    //! "Relevance"
    //! ------------
    else if(propertyName =="Relevance")
    {
        int value = data.value<Property>().getData().toInt();
        QSlider *slider = static_cast<QSlider*>(editor);
        slider->setValue(value);
        //! the signal "sliderReleased()" is not OK. Please, change it ... to do
        connect(slider,SIGNAL(sliderReleased()),this,SLOT(commitAndCloseRelevanceControl()));
    }
    //! --------------------
    //! "Initial size seed"
    //! --------------------
    else if(propertyName=="Initial size seed")
    {
        int value = data.value<Property>().getData().toInt();
        QComboBox *comboBox = static_cast<QComboBox*>(editor);
        comboBox->setCurrentIndex(value);
        connect(editor,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseInitialSizeSeedSelector()));
    }
    //! ------------------------
    //! "Element midside nodes"
    //! ------------------------
    else if(propertyName=="Element midside nodes")
    {
        int value = data.value<Property>().getData().toInt();
        QComboBox *comboBox = static_cast<QComboBox*>(editor);
        comboBox->setCurrentIndex(value);
        connect(editor,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseMidsideNodesSelector()));
    }
    else if(propertyName =="Straight sided elements")
    {
        bool val = data.value<Property>().getData().toBool();
        QComboBox *comboBox = static_cast<QComboBox*>(editor);
        comboBox->setCurrentIndex(val);
        connect(editor,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseStraightSidedElementsControl()));
    }
    else if(propertyName =="Remap")
    {
        bool value = data.value<Property>().getData().toBool();
        QComboBox *comboBox = static_cast<QComboBox*>(editor);
        switch(value)
        {
        case true: comboBox->setCurrentIndex(0); break;
        case false: comboBox->setCurrentIndex(1); break;
        }
        connect(editor,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseRemapFlagControl()));
    }
    //! ------------
    //! "Algorithm"
    //! ------------
    else if(propertyName =="Algorithm")
    {
        switch(this->getCurrentNode()->getType())
        {
        case SimulationNodeClass::nodeType_importedBodyScalar:
        {
            int value = data.value<Property>().getData().toInt();
            QComboBox *comboBox = static_cast<QComboBox*>(editor);
            switch(value)
            {
            case 0: comboBox->setCurrentIndex(0); break;
            case 1: comboBox->setCurrentIndex(1); break;
            case 2: comboBox->setCurrentIndex(2); break;
            case 3: comboBox->setCurrentIndex(3); break;
            }
            connect(editor,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseInterpolationAlgorithm()));
        }
            break;

        case SimulationNodeClass::nodeType_meshPrismaticLayer:
        {
            int value = data.value<Property>().getData().toInt();
            QComboBox *comboBox = static_cast<QComboBox*>(editor);
            switch(value)
            {
            case 0: comboBox->setCurrentIndex(0); break;
            case 1: comboBox->setCurrentIndex(1); break;
            }
            connect(editor,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseComboBox()));
        }
            break;
        }
    }
    //! ---------------------
    //! "Boundary mesh type"
    //! ---------------------
    else if(propertyName =="Boundary mesh type")
    {
        int value = data.value<Property>().getData().toInt();
        QComboBox *comboBox = static_cast<QComboBox*>(editor);
        comboBox->setCurrentIndex(value);
        connect(editor,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseComboBox()));
    }
    //! ------------
    //! "Submeshes"
    //! ------------
    else if(propertyName =="Submeshes")
    {
        bool value = data.value<Property>().getData().toBool();
        QComboBox *comboBox = static_cast<QComboBox*>(editor);
        switch(value)
        {
        case true: comboBox->setCurrentIndex(0); break;
        case false: comboBox->setCurrentIndex(1); break;
        }
    }
    else if(propertyName=="Generate")
    {
        ;
    }
    else if(propertyName=="Step number")
    {
        int value = data.value<Property>().getData().toInt();
        QSpinBox *spinBox = static_cast<QSpinBox*>(editor);
        spinBox->setValue(value);
    }
    else if(propertyName =="Step selection mode")
    {
        int value = data.value<Property>().getData().toInt();
        QComboBox *comboBox = static_cast<QComboBox*>(editor);
        switch(value)
        {
            case 0: comboBox->setCurrentIndex(0); break;
            case 1: comboBox->setCurrentIndex(1); break;
            case 2: comboBox->setCurrentIndex(2); break;
            case 3: comboBox->setCurrentIndex(3); break;
            case 4: comboBox->setCurrentIndex(4); break;
        }
        connect(editor,SIGNAL(currentIndexChanged(int)),this, SLOT(commitAndCloseStepSelectionMode()));
    }
    else if(propertyName =="Source file")
    {
        ;   //! do nothing
    }
    else if(propertyName =="Source directory" || propertyName =="Target directory")
    {
        ;   //! do nothing
    }
    else if(propertyName =="Smoothing")
    {
        int value = data.value<Property>().getData().toInt();
        QComboBox *comboBox = static_cast<QComboBox*>(editor);
        switch(value)
        {
            case 0: comboBox->setCurrentIndex(0); break;
            case 1: comboBox->setCurrentIndex(1); break;
            case 2: comboBox->setCurrentIndex(2); break;
            case 3: comboBox->setCurrentIndex(3); break;
        }
        connect(editor,SIGNAL(currentIndexChanged(int)),this, SLOT(commitAndCloseSmoothingControl()));
    }
    else if(propertyName =="Show mesh nodes")
    {
        bool value = data.value<Property>().getData().toBool();
        QComboBox *comboBox = static_cast<QComboBox*>(editor);
        if(value==false) comboBox->setCurrentIndex(0);
        else comboBox->setCurrentIndex(1);
        connect(editor,SIGNAL(currentIndexChanged(int)),this, SLOT(commitAndCloseShowMeshNodes()));
    }
    else if(propertyName =="Environment temperature" || propertyName =="X coordinate"
            || propertyName =="Y coordinate" || propertyName =="Z coordinate")
    {
        double value = data.value<Property>().getData().toDouble();
        QLineEdit *lineEdit = static_cast<QLineEdit*>(editor);
        lineEdit->setText(QString("%1").arg(value));
    }
    else if(propertyName =="Initial substeps" || propertyName =="Minimum substeps" || propertyName =="Maximum substeps"
            || propertyName =="Number of substeps")
    {
        int N = data.value<Property>().getData().toInt();
        QLineEdit *lineEdit = static_cast<QLineEdit*>(editor);
        lineEdit->setText(QString("%1").arg(N));
    }
    else if(propertyName=="Origin X" || propertyName=="Origin Y" || propertyName=="Origin Z" ||
            propertyName=="Offset X" || propertyName=="Offset Y" || propertyName=="Offset Z" ||
            propertyName=="Rotation X" || propertyName=="Rotation Y" || propertyName=="Rotation Z")
    {
        double value = data.value<Property>().getData().toDouble();
        QLineEdit *lineEdit = static_cast<QLineEdit*>(editor);
        lineEdit->setText(QString("%3").arg(value));
        connect(editor,SIGNAL(editingFinished()),this,SLOT(commitAndCloseLineEdiCSTransformation()));
    }
    else if(propertyName =="Define by ") //! warning: observe the space at the end
    {
        Property::defineBy value = data.value<Property>().getData().value<Property::defineBy>();
        QComboBox *comboBox = static_cast<QComboBox*>(editor);
        switch(value)
        {
        case Property::defineBy_geometrySelection: comboBox->setCurrentIndex(0); break;
        case Property::defineBy_globalCoordinates: comboBox->setCurrentIndex(1); break;
        }
        connect(editor,SIGNAL(currentIndexChanged(int)),this, SLOT(commitAndCloseComboBoxDefineBy_()));
    }
    //! ------------
    //! "Direction"
    //! ------------
    else if(propertyName == "Direction")
    {
        QVector<double> vec = data.value<Property>().getData().value<QVector<double>>();
        DirectionSelector *dirSelector = static_cast<DirectionSelector*>(editor);
        dirSelector->setDirection(vec);
    }
    //! --------------
    //! "X component"
    //! --------------
    else if(propertyName =="X component")
    {
        Property::loadDefinition value = data.value<Property>().getData().value<Property::loadDefinition>();
        LineEdit *le = static_cast<LineEdit*>(editor);
        le->setData(value);
        connect(le,SIGNAL(editingFinished(QString)),this,SLOT(commitAndCloseLe_Xcomponent()));
    }
    //! --------------
    //! "Y component"
    //! --------------
    else if(propertyName =="Y component")
    {
        Property::loadDefinition value = data.value<Property>().getData().value<Property::loadDefinition>();
        LineEdit *le = static_cast<LineEdit*>(editor);
        le->setData(value);
        connect(le,SIGNAL(editingFinished(QString)),this,SLOT(commitAndCloseLe_Ycomponent()));
    }
    //! --------------
    //! "Z component"
    //! --------------
    else if(propertyName =="Z component")
    {
        Property::loadDefinition value = data.value<Property>().getData().value<Property::loadDefinition>();
        LineEdit *le = static_cast<LineEdit*>(editor);
        le->setData(value);
        connect(le,SIGNAL(editingFinished(QString)),this,SLOT(commitAndCloseLe_Zcomponent()));
    }
    //! -------------------------------------------
    //! "Film coefficient" "Reference temperature"
    //! -------------------------------------------
    else if(propertyName == "Film coefficient" || propertyName =="Reference temperature")
    {
        Property::loadDefinition value = data.value<Property>().getData().value<Property::loadDefinition>();
        LineEdit *le = static_cast<LineEdit*>(editor);
        le->setData(value);
        if(propertyName == "Film coefficient") connect(le,SIGNAL(editingFinished(QString)),this,SLOT(commitAndCloseLe_filmCoefficient()));
        if(propertyName == "Reference temperature") connect(le,SIGNAL(editingFinished(QString)),this,SLOT(commitAndCloseLe_referenceTemperature()));
    }
    //! ------------
    //! "Magnitude"
    //! ------------
    else if(propertyName =="Magnitude")
    {
        Property::loadDefinition value = data.value<Property>().getData().value<Property::loadDefinition>();
        LineEdit *le = static_cast<LineEdit*>(editor);
        le->setData(value);
        connect(le,SIGNAL(editingFinished(QString)),this,SLOT(commitAndCloseLe_magnitude()));
    }
    //! ---------------
    //! "Overpressure"
    //! ---------------
    else if(propertyName=="Overpressure")
    {
        SimulationNodeClass *curNode = this->getCurrentNode();
        Property::contactType type = curNode->getPropertyItem("Type")->data(Qt::UserRole).value<Property>().getData().value<Property::contactType>();
        if(type!= Property::contactType_tied && type!= Property::contactType_bonded)
        {
            Property::overpressureFunction value = data.value<Property>().getData().value<Property::overpressureFunction>();
            QComboBox *comboBox = static_cast<QComboBox*>(editor);
            switch(value)
            {
            case Property::overpressureFunction_linear: comboBox->setCurrentIndex(0); break;
            case Property::overpressureFunction_exponential: comboBox->setCurrentIndex(1); break;
            }
            connect(editor,SIGNAL(currentIndexChanged(int)),this, SLOT(commitAndCloseComboBox()));
        }
    }
    //! --------------
    //! to be removed
    //! --------------
    else if(propertyName=="Formulation")
    {
        Property::contactFormulation value = data.value<Property>().getData().value<Property::contactFormulation>();
        QComboBox *comboBox = static_cast<QComboBox*>(editor);
        switch(value)
        {
        case Property::contactFormulation_penalty: comboBox->setCurrentIndex(0); break;
        case Property::contactFormulation_MPC: comboBox->setCurrentIndex(1); break;
        }
        connect(editor,SIGNAL(currentIndexChanged(int)),this, SLOT(commitAndCloseComboBox()));
    }
    //! ---------
    //! Behavior
    //! ---------
    else if(propertyName=="Behavior")
    {
        SimulationNodeClass *curNode = this->getCurrentNode();
        Property::contactType type = curNode->getPropertyItem("Type")->data(Qt::UserRole).value<Property>().getData().value<Property::contactType>();
        if(type!=Property::contactType_tied)
        {
            Property::contactBehavior value = data.value<Property>().getData().value<Property::contactBehavior>();
            QComboBox *comboBox = static_cast<QComboBox*>(editor);
            switch(value)
            {
            case Property::contactBehavior_asymmetric: comboBox->setCurrentIndex(0); break;
            case Property::contactBehavior_symmetric: comboBox->setCurrentIndex(1); break;
            }
            connect(editor,SIGNAL(currentIndexChanged(int)),this, SLOT(commitAndCloseComboBox()));
        }
    }
    //! -------
    //! "Type"
    //! -------
    else if(propertyName=="Type")
    {
        //SimulationNodeClass *node = this->getCurrentNode();
        //SimulationNodeClass::nodeType nodeType = node->getType();
        //switch(nodeType)
        //{
        //case SimulationNodeClass::nodeType_connectionPair:
        //{
            Property::contactType value = data.value<Property>().getData().value<Property::contactType>();
            QComboBox *comboBox = static_cast<QComboBox*>(editor);
            switch(value)
            {
            case Property::contactType_bonded: comboBox->setCurrentIndex(0); break;
            case Property::contactType_frictional: comboBox->setCurrentIndex(1); break;
            case Property::contactType_frictionless: comboBox->setCurrentIndex(2); break;
            case Property::contactType_tied: comboBox->setCurrentIndex(3); break;
            }
            connect(editor,SIGNAL(currentIndexChanged(int)),this, SLOT(commitAndCloseComboBox()));
        //}
        //    break;
        //}
    }
    //! --------------
    //! "Sizing type"
    //! --------------
    else if(propertyName =="Sizing type")
    {
        int typeOfDivision = data.value<Property>().getData().toInt();
        QComboBox *comboBox = static_cast<QComboBox*>(editor);
        comboBox->setCurrentIndex(typeOfDivision);
        connect(comboBox,SIGNAL(currentIndexChanged(int)),this, SLOT(commitAndCloseComboBoxTypeOfSizing()));
    }
    //! --------------------
    //! "Coordinate system"
    //! --------------------
    else if(propertyName == "Coordinate system")
    {
        cout<<"Set editor data for \"Coordinate system\""<<endl;

        //void *p = data.value<Property>().getData().value<void*>();
        //QComboBox *comboBox = static_cast<QComboBox*>(editor);
        //QVariant t;
        //t.setValue(p);
        //int index = comboBox->findData(t);
        //if(index!=-1) comboBox->setCurrentIndex(index);
        //else comboBox->setCurrentIndex(0);

        void *p = data.value<Property>().getData().value<void*>();
        QComboBox *comboBox = static_cast<QComboBox*>(editor);
        QStandardItem *item = (QStandardItem*)(p);
        //if(item==Q_NULLPTR) exit(1);
        //cout<<item->data(Qt::DisplayRole).toString().toStdString()<<"____"<<endl;
        //exit(2);
        if(item->data(Qt::UserRole).value<SimulationNodeClass*>()->getPropertyItem("Time tag")==Q_NULLPTR) exit(5555);
        QString timeTag = item->data(Qt::UserRole).value<SimulationNodeClass*>()->getPropertyValue<QString>("Time tag");

        //cout<<"____time tag: "<<timeTag.toStdString()<<"____"<<endl;

        QStandardItemModel *model = (QStandardItemModel*)(comboBox->model());
        //cout<<"____number of coordinate systems defined: "<<model->rowCount()<<"____"<<endl;
        bool found = false;
        int n=0;
        for(; n<model->rowCount(); n++)
        {
            QModelIndex mi = model->index(n,0,QModelIndex());
            void *pp = model->data(mi,Qt::UserRole).value<void*>();
            QStandardItem *curItem = (QStandardItem*)(pp);
            SimulationNodeClass *curNode = curItem->data(Qt::UserRole).value<SimulationNodeClass*>();
            //cout<<"____"<<curNode->getPropertyValue<QString>("Time tag").toStdString()<<"____"<<endl;
            if(timeTag == curNode->getPropertyValue<QString>("Time tag"))
            {
                found = true;
                break;
            }
        }
        if(found == true) comboBox->setCurrentIndex(n);
        else comboBox->setCurrentIndex(0);

        connect(editor,SIGNAL(currentIndexChanged(int)),SLOT(commitAndCloseCSSelector()));
    }
    //! ----------------
    //! "Remote points"
    //! ----------------
    else if(propertyName =="Remote points")
    {
        void *p = data.value<Property>().getData().value<void*>();
        QComboBox *comboBox = static_cast<QComboBox*>(editor);
        QStandardItem *item = (QStandardItem*)(p);
        QString timeTag = item->data(Qt::UserRole).value<SimulationNodeClass*>()->getPropertyValue<QString>("Time tag");

        QStandardItemModel *model = (QStandardItemModel*)(comboBox->model());
        //cout<<"____number of named selections: "<<model->rowCount()<<"____"<<endl;
        bool found = false;
        int n=0;
        for(; n<model->rowCount(); n++)
        {
            QModelIndex mi = model->index(n,0,QModelIndex());
            void *pp = model->data(mi,Qt::UserRole).value<void*>();
            QStandardItem *curItem = (QStandardItem*)(pp);
            SimulationNodeClass *curNode = curItem->data(Qt::UserRole).value<SimulationNodeClass*>();
            //cout<<"____"<<curNode->getPropertyValue<QString>("Time tag").toStdString()<<"____"<<endl;
            if(timeTag == curNode->getPropertyValue<QString>("Time tag"))
            {
                found = true;
                break;
            }
        }
        if(found == true) comboBox->setCurrentIndex(n);
        else comboBox->setCurrentIndex(0);

        connect(editor,SIGNAL(currentIndexChanged(int)),SLOT(commitAndCloseRemotePointsSelector()));
    }
    //! ------------------
    //! "Named selection"
    //! ------------------
    else if(propertyName == "Named selection")
    {
        void *p = data.value<Property>().getData().value<void*>();
        QComboBox *comboBox = static_cast<QComboBox*>(editor);
        QStandardItem *item = (QStandardItem*)(p);
        QString timeTag = item->data(Qt::UserRole).value<SimulationNodeClass*>()->getPropertyValue<QString>("Time tag");

        QStandardItemModel *model = (QStandardItemModel*)(comboBox->model());
        //cout<<"____number of named selections: "<<model->rowCount()<<"____"<<endl;
        bool found = false;
        int n=0;
        for(; n<model->rowCount(); n++)
        {
            QModelIndex mi = model->index(n,0,QModelIndex());
            void *pp = model->data(mi,Qt::UserRole).value<void*>();
            QStandardItem *curItem = (QStandardItem*)(pp);
            SimulationNodeClass *curNode = curItem->data(Qt::UserRole).value<SimulationNodeClass*>();
            //cout<<"____"<<curNode->getPropertyValue<QString>("Time tag").toStdString()<<"____"<<endl;
            if(timeTag == curNode->getPropertyValue<QString>("Time tag"))
            {
                found = true;
                break;
            }
        }
        if(found == true) comboBox->setCurrentIndex(n);
        else comboBox->setCurrentIndex(0);
        connect(editor,SIGNAL(currentIndexChanged(int)),SLOT(commitAndCloseNSSelector()));
    }
    //! -------------------------------
    //! "Contact pairs" (model change)
    //! -------------------------------
    else if(propertyName == "Contact pair")
    {
        void *p = data.value<Property>().getData().value<void*>();
        QComboBox *comboBox = static_cast<QComboBox*>(editor);
        QStandardItem *item = (QStandardItem*)(p);
        QString timeTag = item->data(Qt::UserRole).value<SimulationNodeClass*>()->getPropertyValue<QString>("Time tag");

        QStandardItemModel *model = (QStandardItemModel*)(comboBox->model());
        //cout<<"____number of named selections: "<<model->rowCount()<<"____"<<endl;
        bool found = false;
        int n=0;
        for(; n<model->rowCount(); n++)
        {
            QModelIndex mi = model->index(n,0,QModelIndex());
            void *pp = model->data(mi,Qt::UserRole).value<void*>();
            QStandardItem *curItem = (QStandardItem*)(pp);
            SimulationNodeClass *curNode = curItem->data(Qt::UserRole).value<SimulationNodeClass*>();
            //cout<<"____"<<curNode->getPropertyValue<QString>("Time tag").toStdString()<<"____"<<endl;
            if(timeTag == curNode->getPropertyValue<QString>("Time tag"))
            {
                found = true;
                break;
            }
        }
        if(found == true) comboBox->setCurrentIndex(n);
        else comboBox->setCurrentIndex(0);
        connect(editor,SIGNAL(currentIndexChanged(int)),SLOT(commitAndCloseContactPairSelector()));
    }
    //! ---------------------------
    //! "Boundary named selection"
    //! ---------------------------
    else if(propertyName == "Boundary named selection")
    {
        void *p = data.value<Property>().getData().value<void*>();
        QComboBox *comboBox = static_cast<QComboBox*>(editor);
        QStandardItem *item = (QStandardItem*)(p);
        QString timeTag = item->data(Qt::UserRole).value<SimulationNodeClass*>()->getPropertyValue<QString>("Time tag");

        QStandardItemModel *model = (QStandardItemModel*)(comboBox->model());
        //cout<<"____number of named selections: "<<model->rowCount()<<"____"<<endl;
        bool found = false;
        int n=0;
        for(; n<model->rowCount(); n++)
        {
            QModelIndex mi = model->index(n,0,QModelIndex());
            void *pp = model->data(mi,Qt::UserRole).value<void*>();
            QStandardItem *curItem = (QStandardItem*)(pp);
            SimulationNodeClass *curNode = curItem->data(Qt::UserRole).value<SimulationNodeClass*>();
            //cout<<"____"<<curNode->getPropertyValue<QString>("Time tag").toStdString()<<"____"<<endl;
            if(timeTag == curNode->getPropertyValue<QString>("Time tag"))
            {
                found = true;
                break;
            }
        }
        if(found == true) comboBox->setCurrentIndex(n);
        else comboBox->setCurrentIndex(0);
        connect(editor,SIGNAL(currentIndexChanged(int)),SLOT(commitAndCloseShapeSelectorEditor_boundaryPrismaticLayer()));
    }
    //! -----------------
    //! "Slave" "Master"
    //! -----------------
    else if(propertyName =="Slave" || propertyName =="Master")
    {
        //SimulationNodeClass *node = this->getCurrentNode();
        //Property::ScopingMethod theScopingMethod = node->getPropertyValue<Property::ScopingMethod>("Scoping method");
        //if(theScopingMethod==Property::ScopingMethod_Automatic)
        //{
        //    cout<<"____automatic method____"<<endl;
        //    return;
        //}

        if(data.value<Property>().getData().canConvert<std::vector<GeometryTag>>())
        {
            ShapeSelector *shapeSelector = static_cast<ShapeSelector*>(editor);
            std::vector<GeometryTag> vecLoc = data.value<Property>().getData().value<std::vector<GeometryTag>>();
            shapeSelector->setShape(vecLoc);
            //connect(editor, SIGNAL(editingFinished()),this, SLOT(commitAndCloseShapeSelectorEditor()));
            connect(shapeSelector,SIGNAL(editingSelectionFinished()),this, SLOT(commitAndCloseShapeSelectorEditor()));
        }
        else
        {

            void *p = data.value<Property>().getData().value<void*>();
            QComboBox *comboBox = static_cast<QComboBox*>(editor);
            QStandardItem *item = (QStandardItem*)(p);
            QString timeTag = item->data(Qt::UserRole).value<SimulationNodeClass*>()->getPropertyValue<QString>("Time tag");

            QStandardItemModel *model = (QStandardItemModel*)(comboBox->model());
            //cout<<"____number of named selections: "<<model->rowCount()<<"____"<<endl;
            bool found = false;
            int n=0;
            for(; n<model->rowCount(); n++)
            {
                QModelIndex mi = model->index(n,0,QModelIndex());
                void *pp = model->data(mi,Qt::UserRole).value<void*>();
                QStandardItem *curItem = (QStandardItem*)(pp);
                SimulationNodeClass *curNode = curItem->data(Qt::UserRole).value<SimulationNodeClass*>();
                //cout<<"____"<<curNode->getPropertyValue<QString>("Time tag").toStdString()<<"____"<<endl;
                if(timeTag == curNode->getPropertyValue<QString>("Time tag"))
                {
                    found = true;
                    break;
                }
            }
            if(found == true) comboBox->setCurrentIndex(n);
            else comboBox->setCurrentIndex(0);
            connect(editor,SIGNAL(currentIndexChanged(int)),SLOT(commitAndCloseNSSelector()));
        }
    }
    //! --------------
    //! "Bolt status"
    //! --------------
    else if(propertyName =="Bolt status")
    {
        cout<<"GeneralDelegate::setEditorData()->____set editor data for \"Bolt pretension\"____"<<endl;
        Property::defineBy value = data.value<Property>().getData().value<Property::defineBy>();
        QComboBox *comboBox = static_cast<QComboBox*>(editor);
        switch(value)
        {
        case Property::boltStatusDefinedBy_load: comboBox->setCurrentIndex(0); break;
        case Property::boltStatusDefinedBy_adjustment: comboBox->setCurrentIndex(1); break;
        case Property::boltStatusDefinedBy_open: comboBox->setCurrentIndex(2); break;
        case Property::boltStatusDefinedBy_lock: comboBox->setCurrentIndex(3); break;
        }
        connect(editor,SIGNAL(currentIndexChanged(int)),this, SLOT(commitAndCloseBoltStatusDefinedBy()));
    }
    //! ------------
    //! "Define by"
    //! ------------
    else if(propertyName == "Define by")
    {
        SimulationNodeClass *node = this->getCurrentNode();
        SimulationNodeClass::nodeType nodeType = node->getType();

        if(nodeType !=SimulationNodeClass::nodeType_structuralAnalysisBoltPretension)
        {
            Property::defineBy value = data.value<Property>().getData().value<Property::defineBy>();
            QComboBox *comboBox = static_cast<QComboBox*>(editor);
            switch(value)
            {
            case Property::defineBy_components: comboBox->setCurrentIndex(0); break;
            case Property::defineBy_vector: comboBox->setCurrentIndex(1); break;
            case Property::defineBy_normal: comboBox->setCurrentIndex(2); break;
            }
            connect(editor,SIGNAL(currentIndexChanged(int)),this, SLOT(commitAndCloseDefineByControlComboBox()));
        }
        /*
        else
        {
            cout<<"GeneralDelegate::setEditorData()->____set editor data for \"Bolt pretension\"____"<<endl;
            Property::defineBy value = data.value<Property>().getData().value<Property::defineBy>();
            QComboBox *comboBox = static_cast<QComboBox*>(editor);
            switch(value)
            {
            case Property::boltStatusDefinedBy_load: comboBox->setCurrentIndex(0); break;
            case Property::boltStatusDefinedBy_adjustment: comboBox->setCurrentIndex(1); break;
            case Property::boltStatusDefinedBy_open: comboBox->setCurrentIndex(2); break;
            case Property::boltStatusDefinedBy_lock: comboBox->setCurrentIndex(3); break;
            }
            connect(editor,SIGNAL(currentIndexChanged(int)),this, SLOT(commitAndCloseBoltStatusDefinedBy()));
        }
        */
    }
    //! --------------------
    //! "Load" "Adjustment"
    //! --------------------
    else if(propertyName == "Load" || propertyName == "Adjustment")
    {
        double value = data.value<Property>().getData().toDouble();
        QLineEdit *lineEdit = static_cast<QLineEdit*>(editor);
        lineEdit->setText(QString("%1").arg(value));
        if(propertyName == "Load") connect(lineEdit,SIGNAL(editingFinished()),this,SLOT(commitAndCloseBoltLoadControl()));
        if(propertyName == "Adjustment") connect(lineEdit,SIGNAL(editingFinished()),this,SLOT(commitAndCloseBoltAdjustmentControl()));
    }
    //! -------------------------------
    //! "Ambient" "Diffuse" "Specular"
    //! -------------------------------
    else if(propertyName == "Ambient" || propertyName == "Diffuse" || propertyName == "Specular")
    {
        double value = data.value<Property>().getData().toDouble();
        QLineEdit *lineEdit = static_cast<QLineEdit*>(editor);
        lineEdit->setText(QString("%1").arg(value));
        connect(editor,SIGNAL(editingFinished()),this, SLOT(commitAndCloseLineEdit()));
    }
    //! ---------------
    //! "Transparency"
    //! ---------------
    else if(propertyName=="Transparency")
    {
        double value = data.value<Property>().getData().toDouble();
        QLineEdit *lineEdit = static_cast<QLineEdit*>(editor);
        lineEdit->setText(QString("%1").arg(value));
        connect(editor,SIGNAL(editingFinished()),this,SLOT(commitAndCloseLineEditTransparency()));
    }
    //! ------------------
    //! "Element control"
    //! ------------------
    else if(propertyName == "Element control")
    {
        Property::elementControl value = data.value<Property>().getData().value<Property::elementControl>();
        QComboBox *comboBox = static_cast<QComboBox*>(editor);
        switch(value)
        {
        case Property::elementControl_programControlled: comboBox->setCurrentIndex(0); break;
        case Property::elementControl_manual: comboBox->setCurrentIndex(1); break;
        }
        connect(editor,SIGNAL(currentIndexChanged(int)),this, SLOT(commitAndCloseElementControlComboBox()));
    }
    //! ---------------------
    //! "Integration scheme"
    //! ---------------------
    else if(propertyName == "Integration scheme")
    {
        Property::integrationScheme value = data.value<Property>().getData().value<Property::integrationScheme>();
        QComboBox *comboBox = static_cast<QComboBox*>(editor);
        switch(value)
        {
        case Property::integrationScheme_full: comboBox->setCurrentIndex(0); break;
        case Property::integrationScheme_reduced: comboBox->setCurrentIndex(1); break;
        }
        connect(comboBox,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseComboBox()));
    }
    //! ------------------------------
    //! "Radial" "Axial" tangential"
    //! ------------------------------
    else if(propertyName == "Radial" || propertyName == "Axial" || propertyName == "Tangential")
    {
        Property::DOFfreedom value = data.value<Property>().getData().value<Property::DOFfreedom>();
        QComboBox *comboBox = static_cast<QComboBox*>(editor);
        switch(value)
        {
        case Property::DOFfreedom_fixed: comboBox->setCurrentIndex(0); break;
        case Property::DOFfreedom_free: comboBox->setCurrentIndex(1); break;
        }
        connect(editor,SIGNAL(currentIndexChanged(int)),this, SLOT(commitAndCloseComboBox()));
    }
    //! ----------------------------------------------------------------------------------------
    //! "Min element size" "Max elemenet size" "Grading" "Face sizing" "Element size" "Pinball"
    //! ----------------------------------------------------------------------------------------
    else if(propertyName == "Min element size" || propertyName == "Max element size" || propertyName == "Grading"
            || propertyName =="Face sizing" || propertyName =="Element size" || propertyName =="Pinball")
    {
        double value = data.value<Property>().getData().toDouble();
        QLineEdit *lineEdit = static_cast<QLineEdit*>(editor);
        lineEdit->setText(QString("%1").arg(value));
        connect(editor,SIGNAL(editingFinished()),this, SLOT(commitAndCloseLineEdit()));
    }
    //! ----------------------
    //! "Number of divisions"
    //! ----------------------
    else if(propertyName =="Number of divisions")
    {
        int Ndiv = data.value<Property>().getData().toInt();
        QSpinBox *spinBox = static_cast<QSpinBox*>(editor);
        spinBox->setValue(Ndiv);
    }
    //! ------------
    //! "Mesh type"
    //! ------------
    else if(propertyName =="Mesh type")
    {
        int val = data.value<Property>().getData().toInt();
        QComboBox *comboBox = static_cast<QComboBox*>(editor);
        comboBox->setCurrentIndex(val);
        connect(editor,SIGNAL(currentIndexChanged(int)),this, SLOT(commitAndCloseComboBox()));
    }
    //! -----------------
    //! "Mesh engine 2D"
    //! Program controlled
    //! Netgen
    //! Express mesh
    //! OCC STL
    //! -----------------
    else if(propertyName == "Mesh engine 2D")
    {
        Property::meshEngine2D value = data.value<Property>().getData().value<Property::meshEngine2D>();
        QComboBox *comboBox = static_cast<QComboBox*>(editor);
        switch(value)
        {
        case Property::meshEngine2D_ProgramControlled: comboBox->setCurrentIndex(0); break;
        case Property::meshEngine2D_Netgen: comboBox->setCurrentIndex(1); break;
        case Property::meshEngine2D_OCC_ExpressMesh: comboBox->setCurrentIndex(2); break;
        case Property::meshEngine2D_OCC_STL: comboBox->setCurrentIndex(3); break;
        }
        connect(editor,SIGNAL(currentIndexChanged(int)),this, SLOT(commitAndCloseComboBox()));
    }
    //! -----------------
    //! "Mesh engine 3D"
    //! -----------------
    else if(propertyName == "Mesh engine 3D")
    {
        Property::meshEngine3D value = data.value<Property>().getData().value<Property::meshEngine3D>();
        QComboBox *comboBox = static_cast<QComboBox*>(editor);
        switch(value)
        {
        case Property::meshEngine3D_Netgen: comboBox->setCurrentIndex(0); break;
        case Property::meshEngine3D_Tetgen: comboBox->setCurrentIndex(1); break;
        }
        connect(editor,SIGNAL(currentIndexChanged(int)),this, SLOT(commitAndCloseComboBox()));
    }
    //! --------------------
    //! "Surface mesh type"
    //! --------------------
    else if(propertyName == "Surface mesh type")
    {
        Property::meshType_Surface value = data.value<Property>().getData().value<Property::meshType_Surface>();
        QComboBox *comboBox = static_cast<QComboBox*>(editor);
        if(value==Property::meshType_Surface_AllTrig)
        {
            comboBox->setCurrentIndex(0);
        }
        else
        {
            comboBox->setCurrentIndex(1);
        }
        connect(editor,SIGNAL(currentIndexChanged(int)),this, SLOT(commitAndCloseComboBox()));
    }
    //! -------------------
    //! "Volume mesh type"
    //! -------------------
    else if(propertyName == "Volume mesh type")
    {
        Property::meshType_Volume value = data.value<Property>().getData().value<Property::meshType_Volume>();
        QComboBox *comboBox = static_cast<QComboBox*>(editor);
        if(value==Property::meshType_Volume_AllTet)
        {
            comboBox->setCurrentIndex(0);
        }
        //else comboBox->setCurrentIndex(1);
        connect(editor,SIGNAL(currentIndexChanged(int)),this, SLOT(commitAndCloseComboBox()));
    }
    //! -------------
    //! "Mesh order"
    //! -------------
    else if(propertyName == "Mesh order")
    {
        Property::meshOrder value = data.value<Property>().getData().value<Property::meshOrder>();
        QComboBox *comboBox = static_cast<QComboBox*>(editor);
        switch(value)
        {
        case Property::meshOrder_First: comboBox->setCurrentIndex(0); break;
        case Property::meshOrder_Second: comboBox->setCurrentIndex(1); break;
        }
        connect(editor,SIGNAL(currentIndexChanged(int)),this, SLOT(commitAndCloseComboBox()));
    }
    //! -----------------
    //! "Scoping method"
    //! -----------------
    else if(propertyName == "Scoping method")
    {
        Property::ScopingMethod value = data.value<Property>().getData().value<Property::ScopingMethod>();
        QComboBox *comboBox = static_cast<QComboBox*>(editor);
        switch(value)
        {
        case Property::ScopingMethod_GeometrySelection: comboBox->setCurrentIndex(0); break;
        case Property::ScopingMethod_NamedSelection: comboBox->setCurrentIndex(1); break;

        //! -------------------------------------------------------------------------
        //! the scoping method "Remote point" and "Automatic" are mutually exclusive
        //! -------------------------------------------------------------------------
        case Property::ScopingMethod_RemotePoint: comboBox->setCurrentIndex(2); break;
        case Property::ScopingMethod_Automatic: comboBox->setCurrentIndex(2); break;
        }
        connect(editor,SIGNAL(currentIndexChanged(int)),this, SLOT(commitAndCloseScopingMethodComboBox()));
    }
    //! --------------------------
    //! "Boundary scoping method"
    //! --------------------------
    else if(propertyName =="Boundary scoping method")
    {
        Property::ScopingMethod value = data.value<Property>().getData().value<Property::ScopingMethod>();
        QComboBox *comboBox = static_cast<QComboBox*>(editor);
        switch(value)
        {
        case Property::ScopingMethod_GeometrySelection: comboBox->setCurrentIndex(0); break;
        case Property::ScopingMethod_NamedSelection: comboBox->setCurrentIndex(1); break;
        }
        connect(comboBox,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseBoundaryScopingMethod()));
    }
    //! -------------
    //! "Suppressed"
    //! -------------
    else if(propertyName == "Suppressed")
    {
        Property::SuppressionStatus value = data.value<Property>().getData().value<Property::SuppressionStatus>();
        QComboBox *comboBox = static_cast<QComboBox*>(editor);
        switch(value)
        {
        case Property::SuppressionStatus_Active: comboBox->setCurrentIndex(0); break;
        case Property::SuppressionStatus_Suppressed: comboBox->setCurrentIndex(1); break;
        }
        connect(editor,SIGNAL(currentIndexChanged(int)),this, SLOT(commitAndCloseComboBox()));
    }
    //! ----------
    //! "Visible"
    //! ----------
    else if(propertyName == "Visible")
    {
        bool isVisible = data.value<Property>().getData().toBool();
        QComboBox *comboBox = static_cast<QComboBox*>(editor);
        if(isVisible) comboBox->setCurrentIndex(0);
        else comboBox->setCurrentIndex(1);
        connect(editor,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseVisibilityComboBox()));
    }
    //! ----------------------
    //! "Geometry" "Location"
    //! ----------------------
    else if(propertyName == "Geometry" || propertyName== "Location")
    {
        ShapeSelector *shapeSelector = static_cast<ShapeSelector*>(editor);
        std::vector<GeometryTag> vecLoc = data.value<Property>().getData().value<std::vector<GeometryTag>>();
        shapeSelector->setShape(vecLoc);

        connect(shapeSelector, SIGNAL(editingSelectionFinished()),this, SLOT(commitAndCloseShapeSelectorEditor()));

        //! -------------------------------------------------------------------
        //! when modifying a coordinate system, connect to the originChanged()
        //! when modifying a remote point, connect to
        //! -------------------------------------------------------------------
        SimulationNodeClass::nodeType type = this->getCurrentNode()->getType();
        switch(type)
        {
        case SimulationNodeClass::nodeType_coordinateSystem:
        {
            if(propertyName=="Geometry")
            {
                //! cesere
                //cout<<"GeneralDelegate::createEditor()->____creating editor for property ->Geometry<-____"<<endl;
                //connect(editor,SIGNAL(editingFinished()),this,SLOT(emitOriginAndDirectionChanged()));
            }
            else    //! "Location"
            {
                //cout<<"GeneralDelegate::createEditor()->____creating editor for property ->Location<-____"<<endl;
                //connect(editor,SIGNAL(editingFinished()),this,SLOT(emitOriginChanged()));
            }
        }
            break;

        case SimulationNodeClass::nodeType_remotePoint:
        {
            if(propertyName == "Location")
            {
                connect(editor,SIGNAL(editingFinished()),this,SLOT(emitRemotePointChangedByLocation()));
            }
        }
            break;
        }
    }
    //! ----------------------------------
    //! "Boundary" - for prismatic layers
    //! ----------------------------------
    else if(propertyName=="Boundary")
    {
        cerr<<"____setting editor data for \"Boundary\": case geometry selection____"<<endl;
        ShapeSelector *shapeSelector = static_cast<ShapeSelector*>(editor);
        std::vector<GeometryTag> vecLoc = data.value<Property>().getData().value<std::vector<GeometryTag>>();
        shapeSelector->setShape(vecLoc);
        //connect(editor, SIGNAL(editingFinished()), this, SLOT(commitAndCloseShapeSelectorEditor_boundaryPrismaticLayer()));
        connect(shapeSelector, SIGNAL(editingSelectionFinished()), this, SLOT(commitAndCloseShapeSelectorEditor_boundaryPrismaticLayer()));
    }
    //! ------------------
    //! "Number of steps"
    //! ------------------
    else if(propertyName =="Number of steps")
    {
        int value = data.value<Property>().getData().toInt();
        QSpinBox *spinBox = static_cast<QSpinBox*>(editor);
        spinBox->setValue(value);
        //connect(editor,SIGNAL(editingFinished()),this,SLOT(NumberOfStepSpinBoxClosed()));
    }
    //! ----------------------
    //! "Current step number"
    //! ----------------------
    else if(propertyName =="Current step number")
    {
        int value = data.value<Property>().getData().toInt();
        QSpinBox *spinBox = static_cast<QSpinBox*>(editor);
        spinBox->setValue(value);
        //connect(editor,SIGNAL(editingFinished()),this,SLOT(CurrentStepSpinBoxClosed()));
    }
    //! ----------------
    //! "Step end time"
    //! ----------------
    else if(propertyName == "Step end time")
    {
        double value = data.value<Property>().getData().toDouble();
        QLineEdit *lineEdit = static_cast<QLineEdit*>(editor);
        lineEdit->setText(QString("%1").arg(value));
        //connect(lineEdit,SIGNAL(editingFinished()),this,SLOT(StepEndTimeLineEditClosed()));
    }
    //! ----------------------------
    //! "Friction coefficient" "C0"
    //! ----------------------------
    else if(propertyName == "Friction coefficient"  || propertyName =="C0")
    {
        double value = data.value<Property>().getData().toDouble();
        QLineEdit *lineEdit = static_cast<QLineEdit*>(editor);
        lineEdit->setText(QString("%1").arg(value));
        connect(editor,SIGNAL(editingFinished()),this, SLOT(commitAndCloseLineEdit()));
    }
    //! ---------------------
    //! "Auto time stepping"
    //! ---------------------
    else if(propertyName == "Auto time stepping")
    {
        Property::autoTimeStepping theTimeStepping = data.value<Property>().getData().value<Property::autoTimeStepping>();
        QComboBox *comboBox = static_cast<QComboBox*>(editor);
        switch(theTimeStepping)
        {
        case Property::autoTimeStepping_ProgramControlled: comboBox->setCurrentIndex(0); break;
        case Property::autoTimeStepping_ON: comboBox->setCurrentIndex(1); break;
        case Property::autoTimeStepping_OFF: comboBox->setCurrentIndex(2); break;
        }
        connect(editor,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseAutoTimeStepping()));
    }
    //! --------------
    //! "Solver type"
    //! --------------
    else if(propertyName == "Solver type")
    {
        Property::solverType theSolverType = data.value<Property>().getData().value<Property::solverType>();
        QComboBox *comboBox = static_cast<QComboBox*>(editor);
        switch(theSolverType)
        {
        case Property::solverType_programControlled: comboBox->setCurrentIndex(0); break;
        case Property::solverType_direct: comboBox->setCurrentIndex(1); break;
        case Property::solverType_iterative: comboBox->setCurrentIndex(2); break;
        }
        connect(editor,SIGNAL(currentIndexChanged(int)),this, SLOT(commitAndCloseSolverTypeComboBox()));
    }
    //! -------
    //! "Name"
    //! -------
    else if(propertyName == "Name")
    {
        QLineEdit *lineEdit = static_cast<QLineEdit*>(editor);
        QString name = data.value<Property>().getData().toString();
        lineEdit->setText(name);
        connect(editor,SIGNAL(editingFinished()),this, SLOT(commitAndCloseLineEdit()));
    }
    //! -----------------
    //! "Time step size"
    //! -----------------
    else if(propertyName =="Time step size")
    {
        QLineEdit *lineEdit = static_cast<QLineEdit*>(editor);
        QString name = data.value<Property>().getData().toString();
        lineEdit->setText(name);
        connect(editor,SIGNAL(returnPressed()),this, SLOT(commitAndCloseLineEdit()));
    }
    //! ------------
    //! "Potential"
    //! ------------
    else if(propertyName =="Potential")
    {
        QLineEdit *lineEdit = static_cast<QLineEdit*>(editor);
        QString name = data.value<Property>().getData().toString();
        lineEdit->setText(name);
        connect(editor,SIGNAL(returnPressed()),this, SLOT(commitAndCloseLineEdit()));
    }
    //! ----------
    //! "Emitter"
    //! ----------
    else if(propertyName=="Emitter")
    {
        int value = data.value<Property>().getData().toInt();
        QComboBox *comboBox = static_cast<QComboBox*>(editor);
        comboBox->setCurrentIndex(value);
        connect(comboBox,SIGNAL(currentIndexChanged(int)),this,SLOT(commitAndCloseComboBoxEmitter()));
    }
    //! ----------------------------------------------
    //! "Particle mass" "Electric charge" "Intensity"
    //! ----------------------------------------------
    else if(propertyName =="Particle mass" || propertyName =="Electric charge" || propertyName =="Intensity")
    {
        QLineEdit *lineEdit = static_cast<QLineEdit*>(editor);
        QString name = data.value<Property>().getData().toString();
        lineEdit->setText(name);
        connect(editor,SIGNAL(returnPressed()),this, SLOT(commitAndCloseLineEdit()));
    }
}

//! -----------------------
//! function: setModelData
//! detail:
//! -----------------------
void GeneralDelegate::setModelData(QWidget *editor, QAbstractItemModel *model, const QModelIndex &index) const
{
    if(!index.isValid()) return;

    //! "content" is the QVariant cast of "Property"
    QVariant content = index.model()->data(index, Qt::UserRole);

    //! ----------------------------------------------
    //! modifications are allowed only for "Property"
    //! ----------------------------------------------
    if(content.canConvert<Property>())
    {
        Property theProp = content.value<Property>();
        QString propertyName = theProp.getName();
        QVariant data;

#ifdef COSTAMP_VERSION
        //! -------------------------------
        //! Time step builder - "Activate"
        //! -------------------------------
        if(propertyName =="Activate")
        {
            // ...
        }
        //! --------------------------------------------
        //! "Intensification pressure" && "Force value"
        //! --------------------------------------------
        if(propertyName=="Intensification pressure" || propertyName=="Force value")
        {
            QLineEdit *le = static_cast<QLineEdit*>(editor);
            double val = le->text().toDouble();
            data.setValue(val);
        }
        //! --------------------
        //! "Time history file"
        //! --------------------
        if(propertyName =="Time history file")
        {
            QFileSelect *fd = static_cast<QFileSelect*>(editor);
            QString timeHistoryFileLoc = fd->getText();
            data.setValue(timeHistoryFileLoc);
        }
#endif

        //! ----------------------
        //! "Jx" "Jy" "Jz" "Mass"
        //! ----------------------
        if(propertyName =="Jx" || propertyName =="Jy" || propertyName =="Jz" || propertyName =="Mass")
        {
            QLineEdit *le =static_cast<QLineEdit*>(editor);
            double val = le->text().toDouble();
            data.setValue(val);
        }
        //! -------------------
        //! "Static/Transient"
        //! -------------------
        if(propertyName =="Static/Transient")
        {
            QComboBox *cb = static_cast<QComboBox*>(editor);
            data.setValue(cb->currentData());
        }
        //! ----------------
        //! "Analysis type"
        //! ----------------
        if(propertyName =="Analysis type")
        {
            QComboBox *cb = static_cast<QComboBox*>(editor);
            data.setValue(cb->currentData());
        }
        //! -----------------------
        //! "Structural time step"
        //! -----------------------
        if(propertyName=="Structural time step")
        {
            QLineEdit *le = static_cast<QLineEdit*>(editor);
            int timeStep = le->text().toInt();
            data.setValue(timeStep);
        }
        //! ---------------------
        //! "Coupling time step"
        //! ---------------------
        if(propertyName =="Coupling time step")
        {
            ;
        }
        //! ----------------------------
        //! "Imported body temperature"
        //! ----------------------------
        if(propertyName =="Imported body temperature")
        {
            QComboBox *cb = static_cast<QComboBox*>(editor);
            data = cb->currentData();
        }
        //! -----------
        //! "Analysis"
        //! -----------
        if(propertyName =="Analysis")
        {
            QComboBox *cb = static_cast<QComboBox*>(editor);
            data = cb->currentData();
        }
        //! ---------------
        //! "Fatigue algo"
        //! ---------------
        if(propertyName =="Fatigue algo")
        {
            QComboBox *cb = static_cast<QComboBox*>(editor);
            int val = cb->currentIndex();
            data.setValue(val);
        }
        //! --------------------
        //! Material assignment
        //! --------------------
        if(propertyName =="Assignment")
        {
            QComboBox *cb = static_cast<QComboBox*>(editor);
            int val = cb->currentIndex();
            data.setValue(val);
        }
        //! ----------------------------------------------
        //! "Metric type" - currently implemented for tet
        //! 0 => "Modified Mean Ratio (MMR)"
        //! 1 => "Modified Condition Number (MCN)"
        //! 2 => "Modified volume-length (MIVL)"
        //! ----------------------------------------------
        if(propertyName =="Metric type")
        {
            QComboBox *cb = static_cast<QComboBox*>(editor);
            data.setValue(cb->currentData().toInt());
        }
        //! -----------
        //! "Grouping"
        //! -----------
        if(propertyName =="Grouping")
        {
            QComboBox *cb = static_cast<QComboBox*>(editor);
            data.setValue(cb->currentData().toInt());
        }
        //! ---------------
        //! "Element list"
        //! ---------------
        if(propertyName =="Element list")
        {
            QLineEdit *le = static_cast<QLineEdit*>(editor);
            data.setValue(le->text());
        }
        //! -------------------
        //! "Selection method"
        //! -------------------
        if(propertyName =="Selection method")
        {
            QComboBox *cb = static_cast<QComboBox*>(editor);
            data.setValue(cb->currentIndex());
        }
        //! --------------------
        //! "Activation status"
        //! --------------------
        if(propertyName =="Activation status")
        {
            QComboBox *cb = static_cast<QComboBox*>(editor);
            switch(cb->currentIndex())
            {
            case 0: data.setValue(Property::modelChangeActivationStatus_Remove); break;
            case 1: data.setValue(Property::modelChangeActivationStatus_Inactive); break;
            case 2: data.setValue(Property::modelChangeActivationStatus_Add); break;
            }
        }
        //! ------------
        //! "Item type"
        //! ------------
        if(propertyName =="Item type")
        {
             QComboBox *cb = static_cast<QComboBox*>(editor);
             switch(cb->currentIndex())
             {
             case 0: data.setValue(0); break;
             case 1: data.setValue(1); break;
             }
        }
        //! --------------------
        //!  "Source file path"
        //! --------------------
        if(propertyName =="Source file path")
        {
            QFileSelect *fileSelect = static_cast<QFileSelect*>(editor);
            QString val = fileSelect->getText();
            data.setValue(val);
        }
        //! -------------------
        //! "Geometry healing"
        //! -------------------
        if(propertyName =="Geometry healing")
        {
            QComboBox *cb = static_cast<QComboBox*>(editor);
            switch(cb->currentIndex())
            {
            case 0: data.setValue(false); break;
            case 1: data.setValue(true); break;
            }
        }
        //! ------------------
        //! "Tolerance value"
        //! ------------------
        if(propertyName =="Tolerance value")
        {
            QLineEdit *le = static_cast<QLineEdit*>(editor);
            double val = le->text().toDouble();
            data.setValue(val);
        }
        //! --------------------
        //! "Angular criterion"
        //! --------------------
        if(propertyName =="Angular criterion")
        {
            QLineEdit *le = static_cast<QLineEdit*>(editor);
            double val = le->text().toDouble();
            data.setValue(val);
        }
        //! ----------------------------------------------
        //! "Solid bodies" "Surface bodies" "Line bodies"
        //! ----------------------------------------------
        if(propertyName =="Solid bodies" || propertyName == "Surface bodies" || propertyName == "Line bodies")
        {
            QComboBox *cb = static_cast<QComboBox*>(editor);
            switch(cb->currentIndex())
            {
            case 0: data.setValue(false); break;
            case 1: data.setValue(true); break;
            }
        }
        //! -------------------
        //! "Geometry healing"
        //! -------------------
        if(propertyName =="Geometry healing")
        {
            QComboBox *cb = static_cast<QComboBox*>(editor);
            switch(cb->currentIndex())
            {
            case 0: data.setValue(false); break;
            case 1: data.setValue(true); break;
            }
        }
        //! --------------
        //! "Defeaturing"
        //! --------------
        if(propertyName =="Defeaturing" || propertyName =="Healing" || propertyName =="Simplification" ||
                propertyName == "Preserve boundary conditions edges" || propertyName == "Project points on geometry")
        {
            QComboBox *cb = static_cast<QComboBox*>(editor);
            switch(cb->currentIndex())
            {
            case 0: data.setValue(false); break;
            case 1: data.setValue(true); break;
            }
        }
        //! ----------------------
        //! "Tessellator"
        //! ----------------------
        if(propertyName=="Tessellator")
        {
            QComboBox *cb = static_cast<QComboBox*>(editor);
            switch(cb->currentIndex())
            {
            case 0: data.setValue(Property::meshEngine2D_OCC_STL); break;
            case 1: data.setValue(Property::meshEngine2D_OCC_ExpressMesh); break;
            }
        }
        //! -----------------------------------------
        //! "Angular deflection" "Linear deflection"
        //! -----------------------------------------
        if(propertyName =="Angular deflection" || propertyName=="Linear deflection")
        {
            QLineEdit *le = static_cast<QLineEdit*>(editor);
            double val = le->text().toDouble();
            data.setValue(val);
        }
        //! --------------------------------
        //! "Min face size" "Max face size"
        //! --------------------------------
        if(propertyName =="Min face size" || propertyName =="Max face size")
        {
            QLineEdit *le = static_cast<QLineEdit*>(editor);
            double val = le->text().toDouble();
            data.setValue(val);
        }
        //! -----------
        //! "Patch conforming"
        //! -----------
        if(propertyName =="Patch conforming")
        {
            QComboBox *cb = static_cast<QComboBox*>(editor);
            bool val = cb->currentData().toBool();
            data.setValue(val);
        }
        //! -----------------
        //! "Surface mesher"
        //! -----------------
        if(propertyName =="Surface mesher")
        {
            QComboBox *cb = static_cast<QComboBox*>(editor);
            bool useBRep = this->getCurrentNode()->getPropertyValue<bool>("Patch conforming");
            if(useBRep)
            {
                Property::meshEngine2D val = cb->currentData().value<Property::meshEngine2D>();
                data.setValue(val);
            }
            else
            {
                switch(cb->currentIndex())
                {
                case 0: data.setValue(Property::meshEngine2D_Netgen_STL); break;
                //case 1: data.setValue(Property::meshEngine2D_TetgenBR); break;
                }
            }
        }
        //! ----------------
        //! "Volume mesher"
        //! ----------------
        if(propertyName =="Volume mesher")
        {
            QComboBox *cb = static_cast<QComboBox*>(editor);
            Property::meshEngine3D val = cb->currentData().value<Property::meshEngine3D>();
            data.setValue(val);
        }
        //! ----------------------------------------
        //! "Envelope sizing" "Ideal length sizing"
        //! ----------------------------------------
        if(propertyName =="Envelope sizing" || propertyName =="Ideal length sizing")
        {
            QComboBox *cb = static_cast<QComboBox*>(editor);
            int val = cb->currentData().toInt();
            data.setValue(val);
        }
        //! --------------------------------------------------------------------
        //! "Relative size" "Absolute size" "Relative length" "Absolute length"
        //! --------------------------------------------------------------------
        if(propertyName=="Relative size" || propertyName=="Absolute size" ||
                propertyName=="Relative length" || propertyName=="Absolute length")
        {
            QLineEdit *le = static_cast<QLineEdit*>(editor);
            double val = le->text().toDouble();
            data.setValue(val);
        }
        //! --------
        //! Mapping
        //! --------
        if(propertyName=="Mapping")
        {
            QComboBox *cb = static_cast<QComboBox*>(editor);
            int val = cb->currentData().toInt();
            data.setValue(val);
        }
        //! -----------------------
        //! "Stress/strain source"
        //! -----------------------
        if(propertyName == "Stress/strain source")
        {
            QComboBox *cb = static_cast<QComboBox*>(editor);
            int val = cb->currentData().toInt();
            data.setValue(val);
        }
        //! ------------------------------------------------
        //! stress strain/strain component for fatigue tool
        //! ------------------------------------------------
        if(propertyName =="Component")
        {
            QComboBox *cb = static_cast<QComboBox*>(editor);
            int val = cb->currentData().toInt();
            data.setValue(val);
        }
        //! ---------------------------------------
        //! triaxality correction for fatigue tool
        //! ---------------------------------------
        if(propertyName =="Triaxiality correction")
        {
            QComboBox *cb = static_cast<QComboBox*>(editor);
            int val = cb->currentData().toInt();
            data.setValue(val);
        }
        //! ----------------------------------
        //! "Number of cycles" - fatigue tool
        //! ----------------------------------
        if(propertyName=="Number of cycles")
        {
            QLineEdit *le = static_cast<QLineEdit*>(editor);
            int val = le->text().toInt();
            data.setValue(val);
        }
        //! ---------------------
        //! "Volume mesh engine"
        //! ---------------------
        if(propertyName =="Volume mesh engine")
        {
            QComboBox *cb = static_cast<QComboBox*>(editor);
            Property::meshEngine3D val = cb->currentData().value<Property::meshEngine3D>();
            data.setValue(val);
        }
        //! ----------------
        //! "Lock boundary"
        //! ----------------
        if(propertyName =="Lock boundary")
        {
            QComboBox *cb = static_cast<QComboBox*>(editor);
            bool val = cb->currentData().toBool();
            data.setValue(val);
        }
        //! --------------------------------------------------------------
        //! "Guiding vectors smoothing steps" "Thickness smoothing steps"
        //! --------------------------------------------------------------
        if(propertyName=="Guiding vectors smoothing steps" || propertyName=="Thickness smoothing steps")
        {
            QLineEdit *le = static_cast<QLineEdit*>(editor);
            int val = le->text().toInt();
            data.setValue(val);
        }
        //!  --------------------------------------------------------------------------
        //! "Curvature sensitivity" "Guiding vector smoothing - curvature sensitivity"
        //! "Thickness smoothing - curvature sensitivity"
        //!  --------------------------------------------------------------------------
        if(propertyName =="Curvature sensitivity" ||
                propertyName =="Guiding vector smoothing - curvature sensitivity" ||
                propertyName =="Thickness smoothing - curvature sensitivity")
        {
            QLineEdit *le = static_cast<QLineEdit*>(editor);
            double val = le->text().toDouble();
            data.setValue(val);
        }
        //! -------------------------------------------
        //! "Options" (for prismatic layer generation)
        //! -------------------------------------------
        if(propertyName =="Options")
        {
            QComboBox *cb = static_cast<QComboBox*>(editor);
            data.setValue(cb->currentData().toInt());
        }
        //! -------------------
        //! "Number of layers"
        //! -------------------
        if(propertyName=="Number of layers")
        {
            QLineEdit *le = static_cast<QLineEdit*>(editor);
            int val = le->text().toInt();
            data.setValue(val);
        }
        //! -----------------------------------------------------
        //! "Expansion ratio" "Total thickness" "First layer height"
        //! -----------------------------------------------------
        if(propertyName=="Expansion ratio" || propertyName =="Total thickness" || propertyName =="First layer height")
        {
            QLineEdit *le = static_cast<QLineEdit*>(editor);
            double val = le->text().toDouble();
            data.setValue(val);
        }
        //! --------------------------------------------------------
        //! "Check self intersections" "Check mutual intersections"
        //! --------------------------------------------------------
        if(propertyName =="Check self intersections" || propertyName =="Check mutual intersections")
        {
            QComboBox *cb = static_cast<QComboBox*>(editor);
            data.setValue(cb->currentData());
        }
        /*
        //! -----------------
        //! "Split pyramids"
        //! -----------------
        if(propertyName =="Split pyramids")
        {
            QComboBox *cb = static_cast<QComboBox*>(editor);
            data.setValue(cb->currentData());
        }
        */
        //! ----------------
        //! "Run in memory"
        //! ----------------
        if(propertyName =="Run in memory")
        {
            QComboBox *cb = static_cast<QComboBox*>(editor);
            if(cb->currentIndex()==0) data.setValue(true);
            else data.setValue(false);
        }
        //! ----------------
        //! "Pair distance"
        //! ----------------
        if(propertyName =="Pair distance")
        {
            QLineEdit *le = static_cast<QLineEdit*>(editor);
            double val = le->text().toDouble();
            data.setValue(val);
        }
        //! --------
        //! "Level"
        //! --------
        if(propertyName =="Level")
        {
            QSlider *slider = static_cast<QSlider*>(editor);
            int val = slider->value();
            double level = MIN_MESH_REDUCTION+double(val)*(MAX_MESH_REDUCTION-MIN_MESH_REDUCTION)/100.0;
            data.setValue(level);
        }
        //! ---------
        //! "Method"
        //! ---------
        if(propertyName=="Method")
        {
            QComboBox *cb = static_cast<QComboBox*>(editor);
            data.setValue(cb->currentData());
        }
        //! ------------
        //! "Min" "Max"
        //! ------------
        if(propertyName =="Min" || propertyName =="Max")
        {
            QLineEdit *le = static_cast<QLineEdit*>(editor);
            double val = le->text().toDouble();
            data.setValue(val);
        }
        //! --------------
        //! "# intervals"
        //! --------------
        if(propertyName =="# intervals")
        {
            QSpinBox *spinBox = static_cast<QSpinBox*>(editor);
            int val = spinBox->value();
            data.setValue(val);
        }
        //! -------------
        //! "Scale type"
        //! -------------
        else if(propertyName =="Scale type")
        {
            QComboBox *cb = static_cast<QComboBox*>(editor);
            int val = cb->currentIndex();
            data.setValue(val);
        }
        //! --------------------
        //! "--Recurrence rate"
        //! --------------------
        else if(propertyName =="--Recurrence rate")
        {
            QLineEdit *le = static_cast<QLineEdit*>(editor);
            int val = le->text().toInt();
            data.setValue(val);
        }
        //! -------------------
        //! "Store results at"
        //! -------------------
        else if(propertyName=="Store results at")
        {
            QComboBox *cb = static_cast<QComboBox*>(editor);
            int val = cb->currentIndex();
            data.setValue(val);
        }
        //! -----------------------------------------------------------------------
        //! Output settings for "Stress" "Strain" "Reaction forces" "Contact data"
        //! -----------------------------------------------------------------------
        else if(propertyName =="Stress" || propertyName =="Strain" || propertyName =="Reaction forces" || propertyName=="Contact data")
        {
            QComboBox *cb = static_cast<QComboBox*>(editor);
            int val = cb->currentIndex();
            data.setValue(val);
        }
        //! ----------------------------------------
        //! Advanced controls for the contact group
        //! ----------------------------------------
        else if(propertyName =="K" || propertyName=="Sigma infty" || /*propertyName =="C0" ||*/
                propertyName =="Lambda" || propertyName =="P0")
        {
            SimulationNodeClass *curNode = this->getCurrentNode();

            if(curNode->getType()==SimulationNodeClass::nodeType_connectionPair)
            {
                Property::contactBehavior val = curNode->getPropertyValue<Property::contactBehavior>("Behavior");
                Property::overpressureFunction val1 = curNode->getPropertyValue<Property::overpressureFunction>("Overpressure");

                switch(val)
                {
                case Property::contactBehavior_asymmetric:
                {
                    //! -------------------------------------
                    //! asymmetric node to face contact pair
                    //! -------------------------------------
                    switch (val1)
                    {
                    case Property::overpressureFunction_linear:
                    {
                        //! linear: "K", "C0", "Sigma infty", "Lambda" can be modified
                        QLineEdit *le = static_cast<QLineEdit*>(editor);
                        data.setValue(le->text().toDouble());
                    }
                        break;

                    case Property::overpressureFunction_exponential:
                    {
                        //! exponential: "C0", "P0", "Lambda" can be modified
                        if(/*propertyName =="C0" || */propertyName =="P0" || propertyName =="Lambda")
                        {
                            QLineEdit *le = static_cast<QLineEdit*>(editor);
                            data.setValue(le->text().toDouble());
                        }
                    }
                        break;
                    }
                }
                    break;

                case Property::contactBehavior_symmetric:
                {
                    //! ------------------------------------
                    //! symmetric face to face contact pair
                    //! ------------------------------------
                    switch(val1)
                    {
                    case Property::overpressureFunction_linear:
                    {
                        //! linear: C0 cannot be set: CCX sets C0 = 0.0
                        if(propertyName !="C0")
                        {
                            QLineEdit *le = static_cast<QLineEdit*>(editor);
                            data.setValue(le->text().toDouble());
                        }
                    }
                        break;

                    case Property::overpressureFunction_tied:
                    {
                        if(propertyName =="K" || propertyName =="Lambda")
                        {
                            QLineEdit *le = static_cast<QLineEdit*>(editor);
                            data.setValue(le->text().toDouble());
                        }
                    }
                        break;

                    case Property::overpressureFunction_exponential:
                    {
                        if(propertyName=="P0" /*|| propertyName=="C0" */|| propertyName=="Lambda")
                        {
                            QLineEdit *le = static_cast<QLineEdit*>(editor);
                            data.setValue(le->text().toDouble());
                        }
                    }
                        break;
                    }
                }
                    break;
                }
            }
            else if(propertyName =="C0")
            {
                QLineEdit *le = static_cast<QLineEdit*>(editor);
                data.setValue(le->text().toDouble());
            }
            //! ------------------------------------------
            //! "K" and "Sigma infty" are also properties
            //! of "Compression only support"
            //! ------------------------------------------
            if(curNode->getType()==SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_CompressionOnlySupport)
            {
                QLineEdit *le = static_cast<QLineEdit*>(editor);
                data.setValue(le->text().toDouble());
            }
        }
        //! ----------------
        //! "Small sliding"
        //! ----------------
        else if(propertyName =="Small sliding")
        {
            SimulationNodeClass *curNode = this->getCurrentNode();
            Property::contactBehavior val = curNode->getPropertyValue<Property::contactBehavior>("Behavior");
            if(val == Property::contactBehavior_asymmetric)
            {
                QComboBox *cb = static_cast<QComboBox*>(editor);
                data.setValue(cb->currentData().toInt());
            }
        }
        //! ------------
        //! Line search
        //! ------------
        else if(propertyName =="Line search")
        {
            QComboBox *cb = static_cast<QComboBox*>(editor);
            data.setValue(cb->currentData().toInt());
        }
        //! -----------------------
        //! Line search parameters
        //! -----------------------
        else if(propertyName =="Min value" || propertyName =="Max value")
        {
            int lineSearch = this->getCurrentNode()->getPropertyValue<int>("Line search");
            if(lineSearch!=0)
            {
                QLineEdit *le = static_cast<QLineEdit*>(editor);
                data.setValue(le->text().toDouble());
            }
        }
        //! ------------------
        //! "Cutback" factors
        //! ------------------
        else if(propertyName =="D_f" || propertyName =="D_C" || propertyName =="D_B" || propertyName =="D_A"
                || propertyName =="D_S" || propertyName =="D_H" || propertyName =="D_D" || propertyName =="W_G")
        {
            int cutBackType = this->getCurrentNode()->getPropertyValue<int>("Cutback factors");
            if(cutBackType==1)
            {
                if(propertyName !="D_S" && propertyName!="D_H" && propertyName!="W_G")
                {
                    QLineEdit *le = static_cast<QLineEdit*>(editor);
                    data.setValue(le->text().toDouble());
                }
            }
        }
        //! ----------------------
        //! "Time incrementation"
        //! ----------------------
        else if(propertyName == "Time incrementation" || propertyName =="Cutback factors")
        {
            QComboBox *cb = static_cast<QComboBox*>(editor);
            data.setValue(cb->currentIndex());
        }
        //! -------------------------------
        //! Time incrementation parameters
        //! -------------------------------
        else if(propertyName =="I_0" || propertyName =="I_R" || propertyName =="I_P" || propertyName =="I_C" || propertyName =="I_L"
              || propertyName =="I_G" || propertyName =="I_S" || propertyName =="I_A" || propertyName =="I_J" || propertyName =="I_T")
        {
            int timeIncrementationType = this->getCurrentNode()->getPropertyValue<int>("Time incrementation");
            if(timeIncrementationType==1)
            {
                if(propertyName !="I_S" && propertyName !="I_J" && propertyName!="I_T")
                {
                    QLineEdit *le = static_cast<QLineEdit*>(editor);
                    data.setValue(le->text().toInt());
                }
            }
        }
        //! ------------------------------------------
        //! "Flux convergence" "Solution convergence"
        //! ------------------------------------------
        else if(propertyName =="Flux convergence" || propertyName =="Solution convergence")
        {
            QComboBox *cb = static_cast<QComboBox*>(editor);
            data.setValue(cb->currentIndex());
        }
        //! ---------------------------------------------------
        //! --R_alpha_n --R_alpha_P -R_alpha_l --epsilon_alpha
        //! ---------------------------------------------------
        else if(propertyName =="--R_alpha_n" || propertyName =="--R_alpha_P" || propertyName =="--R_alpha_l" || propertyName=="--epsilon_alpha")
        {
            int fluxConvergenceFlag = this->getCurrentNode()->getPropertyValue<int>("Flux convergence");
            if(fluxConvergenceFlag==1)
            {
                QLineEdit *le = static_cast<QLineEdit*>(editor);
                data.setValue(le->text().toDouble());
            }
        }
        //! -----------------------
        //! "--Value" or q_alpha_u
        //! -----------------------
        else if(propertyName=="--Value")
        {
            int fluxConvergenceFlag = this->getCurrentNode()->getPropertyValue<int>("Flux convergence");
            if(fluxConvergenceFlag==1)
            {
                QLineEdit *le = static_cast<QLineEdit*>(editor);
                data.setValue(le->text().toDouble());
            }
            else data.setValue(QString("Calculated by solver"));
        }
        //! ------------
        //! --q_alpha_0
        //! ------------
        else if(propertyName=="--q_alpha_0")
        {
            int fluxConvergenceFlag = this->getCurrentNode()->getPropertyValue<int>("Flux convergence");
            if(fluxConvergenceFlag==1)
            {
                QLineEdit *le = static_cast<QLineEdit*>(editor);
                data.setValue(le->text().toDouble());
            }
            else data.setValue(QString("Calculated by solver"));
        }
        //! ------------------------------------------
        //! --C_alpha_n --C_alpha_epsilon --R_alpha_l
        //! ------------------------------------------
        else if(propertyName == "--C_alpha_n" || propertyName == "--C_alpha_epsilon")
        {
            int solutionConvergenceFlag = this->getCurrentNode()->getPropertyValue<int>("Solution convergence");
            if(solutionConvergenceFlag==1)
            {
                QLineEdit *le = static_cast<QLineEdit*>(editor);
                data.setValue(le->text().toDouble());
            }
        }
        //! ------------------
        //! "Large deflection
        //! ------------------
        else if(propertyName=="Large deflection")
        {
            QComboBox *cb = static_cast<QComboBox*>(editor);
            data.setValue(cb->currentIndex());
        }
        //! ---------------
        //! "Display time"
        //! ---------------
        else if(propertyName=="Display time")
        {
            QLineEdit *le = static_cast<QLineEdit*>(editor);
            data.setValue(le->text().toDouble());
        }
        //! ---------------
        //! "Mode number"
        //! ---------------
        else if(propertyName=="Mode number")
        {
            QLineEdit *le = static_cast<QLineEdit*>(editor);
            data.setValue(le->text().toInt());
        }
        //! ------------------
        //! "Update interval"
        //! ------------------
        else if(propertyName=="Update interval")
        {
            QLineEdit *le = static_cast<QLineEdit*>(editor);
            double val = le->text().toDouble();
            if(val<1.5) val = 1.5;
            data.setValue(val);
        }
        //! ------------
        //! "Tolerance"
        //! ------------
        else if(propertyName=="Tolerance")
        {
            SimulationNodeClass::nodeType nodeType = this->getCurrentNode()->getType();
            if(nodeType==SimulationNodeClass::nodeType_import)
            {
                QComboBox *cb = static_cast<QComboBox*>(editor);
                data.setValue(cb->currentIndex());
            }
            else
            {
                QLineEdit *le = static_cast<QLineEdit*>(editor);
                data.setValue(le->text().toDouble());
            }
        }
        //! -------------
        //! "Set number"
        //! -------------
        else if(propertyName=="Set number")
        {
            QLineEdit *le = static_cast<QLineEdit*>(editor);
            data.setValue(le->text().toInt());
        }
        //! -----
        //! "By"
        //! -----
        else if(propertyName=="By")
        {
            QComboBox *cb = static_cast<QComboBox*>(editor);
            int val = cb->currentIndex();
            data.setValue(val);
        }
        //! -------
        //! "Type"
        //! -------
        else if(propertyName =="Type ")
        {
            QComboBox *cb = static_cast<QComboBox*>(editor);
            int val = cb->currentIndex();
            data.setValue(val);
        }
        else if(propertyName =="Solution information")
        {
            QComboBox *cb = static_cast<QComboBox*>(editor);
            int val = cb->currentIndex();
            data.setValue(val);
        }
        //! ----------------
        //! "Coupling" type
        //! ----------------
        else if(propertyName =="Coupling")
        {
            QComboBox *cb = static_cast<QComboBox*>(editor);
            int val = cb->currentIndex();
            data.setValue(val);
        }
        //! -----------------
        //! "DOFs selection"
        //! -----------------
        else if(propertyName =="DOFs selection")
        {
            QComboBox *cb = static_cast<QComboBox*>(editor);
            int val = cb->currentIndex();
            data.setValue(val);
        }
        //! --------------
        //! coupling DOFs
        //! --------------
        else if(propertyName =="X component " || propertyName =="Y component " || propertyName =="Z component " ||
                propertyName =="X rotation" || propertyName =="Y rotation" || propertyName =="Z rotation")
        {
            QComboBox *cb = static_cast<QComboBox*>(editor);
            int val = cb->currentIndex();
            data.setValue(val);
        }
        //! -----------------------------
        //! "First color" "Second color"
        //! -----------------------------
        else if(propertyName =="First color" || propertyName =="Second color")
        {
            ColorSelector *cs = static_cast<ColorSelector*>(editor);
            QVector<int> rgb;
            rgb.push_back(cs->color().red());
            rgb.push_back(cs->color().green());
            rgb.push_back(cs->color().blue());
            data.setValue(rgb);
        }
        //! -----------
        //! "Gradient"
        //! -----------
        else if(propertyName=="Gradient")
        {
            QComboBox *cb = static_cast<QComboBox*>(editor);
            data.setValue(cb->currentIndex());
        }
        //! --------------------
        //! "Number of threads"
        //! --------------------
        else if(propertyName =="Number of threads")
        {
            cerr<<"____set model data for \"Number of threads\"____"<<endl;
            QSpinBox *sb = static_cast<QSpinBox*>(editor);
            data.setValue(sb->value());
        }
        //! ----------------
        //! "Analysis time"
        //! ----------------
        else if(propertyName =="Analysis time")
        {
            QLineEdit *le = static_cast<QLineEdit*>(editor);
            data.setValue(le->text().toDouble());
        }
        //! ----------------------------
        //! "Buckets" "Remapping steps"
        //! ----------------------------
        else if(propertyName =="X buckets" ||
                propertyName =="Y buckets" ||
                propertyName =="Z buckets" ||
                propertyName =="Remapping steps")
        {
            QLineEdit *le = static_cast<QLineEdit*>(editor);
            data.setValue(le->text().toInt());
        }
        else if(propertyName == "Split data")
        {
            QComboBox *comboBox = static_cast<QComboBox*>(editor);
            data.setValue(comboBox->currentIndex());
        }
        //! ------------
        //! "Relevance"
        //! ------------
        else if(propertyName =="Relevance")
        {
            QSlider *slider = static_cast<QSlider*>(editor);
            data.setValue(slider->value());
        }
        else if(propertyName =="Initial size seed")
        {
            QComboBox *comboBox = static_cast<QComboBox*>(editor);
            data.setValue(comboBox->currentIndex());
        }
        else if(propertyName =="Element midside nodes")
        {
            QComboBox *comboBox = static_cast<QComboBox*>(editor);
            data.setValue(comboBox->currentIndex());
        }
        /*
        else if(propertyName=="Triangle surface mesher")
        {
            QComboBox *comboBox = static_cast<QComboBox*>(editor);
            data.setValue(comboBox->currentIndex());
        }
        */
        else if(propertyName=="Straight sided elements")
        {
            QComboBox *comboBox = static_cast<QComboBox*>(editor);
            data.setValue(comboBox->currentIndex());
        }
        else if(propertyName =="Remap")
        {
            QComboBox *comboBox = static_cast<QComboBox*>(editor);
            switch(comboBox->currentIndex())
            {
            case 0: data.setValue(true); break;
            case 1: data.setValue(false); break;
            }
        }
        //! ---------------------
        //! "Boundary mesh type"
        //! ---------------------
        else if(propertyName =="Boundary mesh type")
        {
            QComboBox *comboBox = static_cast<QComboBox*>(editor);
            data.setValue(comboBox->currentIndex());
        }
        //! ------------
        //! "Algorithm"
        //! ------------
        else if(propertyName =="Algorithm")
        {
            switch(this->getCurrentNode()->getType())
            {
            case SimulationNodeClass::nodeType_importedBodyScalar:
            {
                QComboBox *comboBox = static_cast<QComboBox*>(editor);
                switch(comboBox->currentIndex())
                {
                case 0: data.setValue(0); break;
                case 1: data.setValue(1); break;
                case 2: data.setValue(2); break;
                case 3: data.setValue(3); break;
                }
            }
                break;

            case SimulationNodeClass::nodeType_meshPrismaticLayer:
            {
                QComboBox *comboBox = static_cast<QComboBox*>(editor);
                data.setValue(comboBox->currentIndex());
            }
                break;
            }
        }
        //! ------------
        //! "Submeshes"
        //! ------------
        else if(propertyName =="Submeshes")
        {
            QComboBox *comboBox = static_cast<QComboBox*>(editor);
            switch(comboBox->currentIndex())
            {
            case 0: data.setValue(true); break;
            case 1: data.setValue(false); break;
            }
        }
        else if(propertyName=="Generate")
        {
            //! to be implemented [?]
            ;
        }
        //else if(propertyName=="Translate")
        //{
            //! to be implemented [?]
        //    ;
        //}
        else if(propertyName =="Step number")
        {
            QSpinBox *spinBox = static_cast<QSpinBox*>(editor);
            data.setValue(spinBox->value());
        }
        else if(propertyName =="Step selection mode")
        {
            QComboBox *comboBox = static_cast<QComboBox*>(editor);
            switch(comboBox->currentIndex())
            {
            case 0: data.setValue(0); break;
            case 1: data.setValue(1); break;
            case 2: data.setValue(2); break;
            case 3: data.setValue(3); break;
            case 4: data.setValue(4); break;
            }
        }
        else if(propertyName == "Source file")
        {
            QFileDialog *aDialog = static_cast<QFileDialog*>(editor);
            if(!aDialog->selectedFiles().isEmpty())
            {
                QString fileName = aDialog->selectedFiles().at(0);
                data.setValue(fileName);
            }
        }
        else if(propertyName =="Source directory" || propertyName =="Target directory")
        {
            QFileDialog *aDialog = static_cast<QFileDialog*>(editor);
            if(!aDialog->selectedFiles().isEmpty())
            {
                QString fileName = aDialog->selectedFiles().at(0);
                data.setValue(fileName);
            }
        }
        else if(propertyName== "Smoothing")
        {
            QComboBox *comboBox = static_cast<QComboBox*>(editor);
            int nsteps = comboBox->currentIndex();
            switch(nsteps)
            {
            case 0: data.setValue(0); break;
            case 1: data.setValue(1); break;
            case 2: data.setValue(2); break;
            case 3: data.setValue(3); break;
            }
        }
        else if(propertyName =="Show mesh nodes")
        {
            QComboBox *comboBox = static_cast<QComboBox*>(editor);
            bool meshNodeVisible = comboBox->currentIndex()==0? false:true;
            data.setValue(meshNodeVisible);
        }
        else if(propertyName =="Environment temperature" || propertyName =="X coordinate"
                || propertyName =="Y coordinate" || propertyName =="Z coordinate")
        {
            QLineEdit *lineEdit = static_cast<QLineEdit*>(editor);
            double value = lineEdit->text().toDouble();
            data.setValue(value);
        }
        else if(propertyName =="Initial substeps" || propertyName =="Minimum substeps" || propertyName =="Maximum substeps"
                || propertyName =="Number of substeps")
        {
            QLineEdit *lineEdit = static_cast<QLineEdit*>(editor);
            double value = lineEdit->text().toInt();
            data.setValue(value);
        }
        else if(propertyName=="Origin X" || propertyName=="Origin Y" || propertyName=="Origin Z" ||
           propertyName=="Offset X" || propertyName=="Offset Y" || propertyName=="Offset Z" ||
           propertyName=="Rotation X" || propertyName=="Rotation Y" || propertyName=="Rotation Z")
        {
            QLineEdit *lineEdit = static_cast<QLineEdit*>(editor);
            double value = lineEdit->text().toDouble();
            data.setValue(value);
        }
        else if(propertyName=="Direction")
        {
            DirectionSelector *dirSelector =static_cast<DirectionSelector*>(editor);
            QVector<double> values = dirSelector->getDirection();
            data.setValue(values);
        }
        //! -------------------------------------------
        //! "Film coefficient" "Reference temperature"
        //! -------------------------------------------
        else if(propertyName =="Film coefficient" || propertyName =="Reference temperature")
        {
            LineEdit *le = static_cast<LineEdit*>(editor);
            data.setValue(le->getLoadDefinition());
        }
        //! ------------------------------------------------------
        //! "Magnitude" "X component" "Y component" "Z component"
        //! ------------------------------------------------------
        else if(propertyName=="Magnitude" || propertyName=="X component" ||
                propertyName =="Y component" || propertyName =="Z component")
        {
            LineEdit *le = static_cast<LineEdit*>(editor);
            data.setValue(le->getLoadDefinition());

            //! -----------
            //! diagnostic
            //! -----------
            switch(le->getLoadDefinition())
            {
            case Property::loadDefinition_constant: cout<<"____constant____"<<endl; break;
            case Property::loadDefinition_tabularData : cout<<"____tabular data____"<<endl; break;
            case Property::loadDefinition_free: cout<<"____free____"<<endl; break;
            }
        }
        //! -------------
        //! Overpressure
        //! -------------
        else if(propertyName=="Overpressure")
        {
            SimulationNodeClass *curNode = this->getCurrentNode();
            Property::contactType type = curNode->getPropertyItem("Type")->data(Qt::UserRole).value<Property>().getData().value<Property::contactType>();
            if(type!= Property::contactType_tied &&  type!= Property::contactType_bonded)
            {
                QComboBox *comboBox = static_cast<QComboBox*>(editor);
                data.setValue(comboBox->currentData(Qt::UserRole).value<Property::overpressureFunction>());
            }
        }
        //! --------------
        //! to be removed
        //! --------------
        //else if(propertyName=="Formulation")
        //{
        //    QComboBox *comboBox = static_cast<QComboBox*>(editor);
        //    data.setValue(comboBox->currentData(Qt::UserRole).value<Property::contactFormulation>());
        //}
        //! ---------
        //! Behavior
        //! ---------
        else if(propertyName=="Behavior")
        {
            SimulationNodeClass *curNode = this->getCurrentNode();
            Property::contactType type = curNode->getPropertyItem("Type")->data(Qt::UserRole).value<Property>().getData().value<Property::contactType>();
            if(type!=Property::contactType_tied)
            {
                QComboBox *comboBox = static_cast<QComboBox*>(editor);
                data.setValue(comboBox->currentData(Qt::UserRole).value<Property::contactBehavior>());
            }
        }
        //! --------------
        //! "Sizing type"
        //! --------------
        else if(propertyName=="Sizing type")
        {
            QComboBox *comboBox = static_cast<QComboBox*>(editor);
            int sizingType = comboBox->currentIndex();
            data.setValue(sizingType);
        }
        //! -------
        //! "Type"
        //! -------
        else if(propertyName=="Type")
        {
            QComboBox *comboBox = static_cast<QComboBox*>(editor);
            SimulationNodeClass *node = this->getCurrentNode();
            switch(node->getType())
            {
            case SimulationNodeClass::nodeType_connectionPair:
            {
                Property::contactType contactType = comboBox->currentData(Qt::UserRole).value<Property::contactType>();
                data.setValue(contactType);
            }
                break;
            }
        }
        //! --------------------
        //! "Coordinate system"
        //! --------------------
        else if(propertyName=="Coordinate system")
        {
            QComboBox *comboBox = static_cast<QComboBox*>(editor);
            QVariant data1 = comboBox->currentData(Qt::UserRole);
            data.setValue(data1);
        }
        //! ---------------------------------------------
        //! "Named selection" "Boundary named selection"
        //! ---------------------------------------------
        else if(propertyName=="Named selection" || propertyName =="Boundary named selection")
        {
            QComboBox *comboBox = static_cast<QComboBox*>(editor);
            QVariant data1 = comboBox->currentData(Qt::UserRole);
            data.setValue(data1);
        }
        //! ------------------------------
        //! "Contact pair" (model change)
        //! ------------------------------
        else if(propertyName =="Contact pair")
        {
            QComboBox *comboBox = static_cast<QComboBox*>(editor);
            QVariant data1 = comboBox->currentData(Qt::UserRole);
            data.setValue(data1);
        }
        //! ---------------------------
        //! "Boundary named selection"
        //! ---------------------------
        else if(propertyName=="Boundary named selection")
        {
            QComboBox *comboBox = static_cast<QComboBox*>(editor);
            QVariant data1 = comboBox->currentData(Qt::UserRole);
            data.setValue(data1);
        }
        //! ----------------
        //! "Remote points"
        //! ----------------
        else if(propertyName=="Remote points")
        {
            QComboBox *comboBox = static_cast<QComboBox*>(editor);
            QVariant data1 = comboBox->currentData(Qt::UserRole);
            data.setValue(data1);
        }
        //! -----------------
        //! "Master" "Slave"
        //! -----------------
        else if(propertyName=="Master" || propertyName =="Slave")
        {
            if(content.value<Property>().getData().canConvert<std::vector<GeometryTag>>())
            {
                ShapeSelector *shapeSelector = static_cast<ShapeSelector*>(editor);
                data.setValue(shapeSelector->getVecLoc());
            }
            else
            {
                QComboBox *comboBox = static_cast<QComboBox*>(editor);
                QVariant data1 = comboBox->currentData(Qt::UserRole);
                data.setValue(data1);
            }
        }
        //! -------------
        //! "Define by "
        //! -------------
        else if(propertyName=="Define by ") //! this is for coordinate systems
        {
            QComboBox *comboBox = static_cast<QComboBox*>(editor);
            switch(comboBox->currentIndex())
            {
            case 0: data.setValue(Property::defineBy_geometrySelection); break;
            case 1: data.setValue(Property::defineBy_globalCoordinates); break;
            }
        }
        //! --------------
        //! "Bolt status"
        //! --------------
        else if(propertyName =="Bolt status")
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
        //! ------------
        //! "Define by"
        //! ------------
        else if(propertyName=="Define by")
        {
            SimulationNodeClass *theNode = this->getCurrentNode();
            SimulationNodeClass::nodeType theNodeType = theNode->getType();

            if(theNodeType!=SimulationNodeClass::nodeType_structuralAnalysisBoltPretension)
            {
                QComboBox *comboBox = static_cast<QComboBox*>(editor);
                switch(comboBox->currentIndex())
                {
                case 0: data.setValue(Property::defineBy_components); break;
                case 1: data.setValue(Property::defineBy_vector); break;
                case 2: data.setValue(Property::defineBy_normal); break;
                }
            }
            /*
            else
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
            */
        }
        else if(propertyName =="Load" || propertyName =="Adjustment" || propertyName=="Ambient" || propertyName=="Diffuse" || propertyName =="Specular")
        {
            QLineEdit *lineEdit = static_cast<QLineEdit*>(editor);
            double value = lineEdit->text().toDouble();
            data.setValue(value);
        }
        else if(propertyName=="Transparency")
        {
            QLineEdit *lineEdit = static_cast<QLineEdit*>(editor);
            double value = lineEdit->text().toDouble();
            data.setValue(value);
        }
        else if(propertyName =="Element control")
        {
            QComboBox *comboBox = static_cast<QComboBox*>(editor);
            data.setValue(comboBox->currentData(Qt::UserRole).value<Property::elementControl>());
        }
        else if(propertyName =="Integration scheme")
        {
            QComboBox *comboBox = static_cast<QComboBox*>(editor);
            data.setValue(comboBox->currentData(Qt::UserRole).value<Property::integrationScheme>());
        }
        else if(propertyName =="Radial" || propertyName =="Axial" || propertyName =="Tangential")
        {
            QComboBox *comboBox = static_cast<QComboBox*>(editor);
            data.setValue(comboBox->currentData(Qt::UserRole).value<Property::DOFfreedom>());
        }
        else if(propertyName=="Max element size" || propertyName=="Min element size" || propertyName =="Grading"
                || propertyName=="Face sizing" || propertyName =="Element size" || propertyName =="Pinball")
        {
            QLineEdit *lineEdit = static_cast<QLineEdit*>(editor);
            double value = lineEdit->text().toDouble();
            data.setValue(value);
        }
        else if(propertyName=="Number of divisions")
        {
            QSpinBox *spinBox = static_cast<QSpinBox*>(editor);
            int Ndiv = spinBox->value();
            data.setValue(Ndiv);
        }
        //! ------------
        //! "Mesh type"
        //! ------------
        else if(propertyName=="Mesh type")
        {
            QComboBox *comboBox = static_cast<QComboBox*>(editor);
            data.setValue(comboBox->currentData(Qt::UserRole).toInt());
        }
        else if(propertyName=="Mesh engine 2D")
        {
            QComboBox *comboBox = static_cast<QComboBox*>(editor);
            data.setValue(comboBox->currentData(Qt::UserRole).value<Property::meshEngine2D>());
        }
        else if(propertyName=="Mesh engine 3D")
        {
            QComboBox *comboBox = static_cast<QComboBox*>(editor);
            int currentIndex = comboBox->currentIndex();
            switch(currentIndex)
            {
            case 0: data.setValue(Property::meshEngine3D_Netgen); break;
            case 1: data.setValue(Property::meshEngine3D_Tetgen); break;
            }
        }
        else if(propertyName=="Surface mesh type")
        {
            QComboBox *comboBox = static_cast<QComboBox*>(editor);
            int currentIndex = comboBox->currentIndex();
            switch(currentIndex)
            {
            case 0: data.setValue(Property::meshType_Surface_AllTrig); break;
            case 1: data.setValue(Property::meshType_Surface_QuadDominant); break;
            }
        }
        else if(propertyName=="Volume mesh type")
        {
            QComboBox *comboBox = static_cast<QComboBox*>(editor);
            int currentIndex = comboBox->currentIndex();
            switch(currentIndex)
            {
            case 0: data.setValue(Property::meshType_Volume_AllTet); break;
            default: break; //! hexa dominant mesh is not available
            }
        }
        else if(propertyName=="Mesh order")
        {
            QComboBox *comboBox = static_cast<QComboBox*>(editor);
            int currentIndex = comboBox->currentIndex();
            switch(currentIndex)
            {
            case 0: data.setValue(Property::meshOrder_First); break;
            case 1: data.setValue(Property::meshOrder_Second); break;
            }
        }
        else if(propertyName=="Scoping method")
        {
            QComboBox *comboBox = static_cast<QComboBox*>(editor);
            int currentIndex = comboBox->currentIndex();
            switch(currentIndex)
            {
            case 0: data.setValue(Property::ScopingMethod_GeometrySelection); break;
            case 1: data.setValue(Property::ScopingMethod_NamedSelection); break;
            case 2:
            {
                SimulationNodeClass::nodeType type = this->getCurrentNode()->getType();
                if(type==SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteDisplacement ||
                        type ==SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteRotation ||
                        type==SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteForce)
                {
                    data.setValue(Property::ScopingMethod_RemotePoint);
                }
                if(type==SimulationNodeClass::nodeType_connectionPair)
                {
                    data.setValue(Property::ScopingMethod_Automatic);
                }
            }
                break;
            }
        }
        else if(propertyName =="Boundary scoping method")
        {
            QComboBox *comboBox = static_cast<QComboBox*>(editor);
            int currentIndex = comboBox->currentIndex();
            switch(currentIndex)
            {
            case 0: data.setValue(Property::ScopingMethod_GeometrySelection); break;
            case 1: data.setValue(Property::ScopingMethod_NamedSelection); break;
            }
        }
        else if(propertyName=="Suppressed")
        {
            QComboBox *comboBox = static_cast<QComboBox*>(editor);
            int currentIndex = comboBox->currentIndex();
            if(currentIndex==0)
            {
                data.setValue(Property::SuppressionStatus_Active);
                emit suppressionChanged(Property::SuppressionStatus_Active);
            }
            else
            {
                data.setValue(Property::SuppressionStatus_Suppressed);
                emit suppressionChanged(Property::SuppressionStatus_Suppressed);
            }
        }
        else if(propertyName=="Visible")
        {
            QComboBox *comboBox = static_cast<QComboBox*>(editor);
            int currentIndex = comboBox->currentIndex();
            if(currentIndex==0)
            {
                data.setValue(true);
                //static int i;
                //cout<<"GeneralDelegate::setModelData->____visibility changed in VISIBLE: "<<i++<<"____"<<endl;
                //emit visibilityChanged(true);
            }
            else
            {
                //static int i;
                //cout<<"GeneralDelegate::setModelData->____visibility changed in HIDDEN: "<<i++<<"____"<<endl;
                data.setValue(false);
                //emit visibilityChanged(false);
            }
        }
        //! ----------------------
        //! "Geometry" "Location"
        //! ----------------------
        else if(propertyName=="Geometry" || propertyName=="Location")
        {
            ShapeSelector *shapeSelector = static_cast<ShapeSelector*>(editor);
            data.setValue(shapeSelector->getVecLoc());
        }
        //! ----------------------------------
        //! "Boundary" - for prismatic layers
        //! ----------------------------------
        else if(propertyName=="Boundary")
        {
            ShapeSelector *shapeSelector = static_cast<ShapeSelector*>(editor);
            data.setValue(shapeSelector->getVecLoc());
        }
        //! ------------------
        //! "Number of steps"
        //! ------------------
        else if(propertyName=="Number of steps")
        {
            QSpinBox *spinBox = static_cast<QSpinBox*>(editor);
            int value = spinBox->value();
            data.setValue(value);
            connect(editor,SIGNAL(editingFinished()),this,SLOT(NumberOfStepSpinBoxClosed()));
        }
        else if(propertyName=="Current step number")
        {
            QSpinBox *spinBox = static_cast<QSpinBox*>(editor);
            int value = spinBox->value();
            data.setValue(value);
            connect(editor,SIGNAL(editingFinished()),this,SLOT(CurrentStepSpinBoxClosed()));
        }
        else if(propertyName=="Step end time")
        {
            QLineEdit *lineEdit = static_cast<QLineEdit*>(editor);
            double value = lineEdit->text().toDouble();
            data.setValue(value);
            connect(editor,SIGNAL(editingFinished()),this,SLOT(StepEndTimeLineEditClosed()));
        }
        else if(propertyName == "Friction coefficient" || propertyName=="C0"
                || propertyName =="Normal stiffness" || propertyName =="Tau")
        {
            QLineEdit *lineEdit = static_cast<QLineEdit*>(editor);
            double value = lineEdit->text().toDouble();
            data.setValue(value);
        }
        //! --------------------
        //! Auto time stepping"
        //! --------------------
        else if(propertyName=="Auto time stepping")
        {
            QComboBox *comboBox = static_cast<QComboBox*>(editor);
            int currentIndex = comboBox->currentIndex();
            switch(currentIndex)
            {
            case 0: data.setValue(Property::autoTimeStepping_ProgramControlled); break;
            case 1: data.setValue(Property::autoTimeStepping_ON); break;
            case 2: data.setValue(Property::autoTimeStepping_OFF); break;
            }
        }
        //! --------------
        //! "Solver type"
        //! --------------
        else if(propertyName=="Solver type")
        {
            QComboBox *comboBox = static_cast<QComboBox*>(editor);
            int currentIndex = comboBox->currentIndex();
            switch(currentIndex)
            {
            case 0: data.setValue(Property::solverType_programControlled); break;
            case 1: data.setValue(Property::solverType_direct); break;
            case 2: data.setValue(Property::solverType_iterative); break;
            }
        }
        //! -------
        //! "Name"
        //! -------
        else if(propertyName=="Name")
        {
            QLineEdit *lineEdit = static_cast<QLineEdit*>(editor);
            data.setValue(lineEdit->text());
        }
        //! -----------------
        //! "Time step size"
        //! -----------------
        else if(propertyName =="Time step size")
        {
            QLineEdit *lineEdit = static_cast<QLineEdit*>(editor);
            data.setValue(lineEdit->text());
        }
        //! ------------
        //! "Potential"
        //! ------------
        else if(propertyName=="Potential")
        {
            QLineEdit *lineEdit = static_cast<QLineEdit*>(editor);
            data.setValue(lineEdit->text());
        }
        //! ----------
        //! "Emitter"
        //! ----------
        else if(propertyName =="Emitter")
        {
            QComboBox *cb = static_cast<QComboBox*>(editor);
            data.setValue(cb->currentIndex());
        }
        //! ----------------------------------------------
        //! "Particle mass" "Electric charge" "Intensity"
        //! ----------------------------------------------
        else if(propertyName =="Particle mass" || propertyName =="Electric charge" || propertyName =="Intensity")
        {
            QLineEdit *lineEdit = static_cast<QLineEdit*>(editor);
            data.setValue(lineEdit->text());
        }

        theProp.setData(data);
        QVariant thePropVariant;
        thePropVariant.setValue(theProp);
        model->setData(index, thePropVariant, Qt::UserRole);
    }
}

//! --------------------------------------------------------
//! function: updateEditorGeometry
//! details:  give the editor the info on size and location
//! --------------------------------------------------------
void GeneralDelegate::updateEditorGeometry(QWidget *editor,const QStyleOptionViewItem &option,const QModelIndex &index) const
{
    Q_UNUSED(index)
    editor->setGeometry(option.rect);
}

//! --------------------
//! function: size hint
//! details:
//! --------------------
QSize GeneralDelegate::sizeHint(const QStyleOptionViewItem & option, const QModelIndex & index) const
{
    QSize size = QStyledItemDelegate::sizeHint(option, index);
    size.setHeight(20);
    return size;
}

//! ----------------
//! function: paint
//! details:
//! ----------------
void GeneralDelegate::paint(QPainter *painter, const QStyleOptionViewItem &opt, const QModelIndex &index) const
{
    QStyleOptionViewItem option = opt;

    //! retrieve the item from the index using the index.internalPointer()
    QExtendedStandardItem *item = static_cast<QExtendedStandardItem*>(index.internalPointer());

    //! font color
    option.palette.setColor(QPalette::Text, option.font.Black);

    //! check if the item is a separator
    bool isSeparator = index.model()->data(index, Qt::UserRole+1).toString() == "separator"? true:false;

    //! calculate indentation
    auto view = qobject_cast<const QTreeView*>(option.widget);
    int indent = view->indentation();

    if(!isSeparator)
    {
        //! if the item is not a separator remove the indentation
        option.rect.adjust(-indent,0,0,0);
        option.state &= QStyle::State_Selected;
        if(index.column()==1)
        {
            option.rect.adjust(indent,0,0,0);
            option.state &= ~QStyle::State_Selected;
        }
    }
    else
    {
        //! if the item is a separator make the first item larger
        //! and bring the length of the second to zero
        int W = qobject_cast<QTreeView *>(this->parent())->size().width();
        if(index.column()==0)
        {
            option.rect.setWidth(W);
        }
        else
        {
            option.rect.adjust(W,0,W,0);
            option.rect.setWidth(0);
        }
    }

    //! get the style to draw with
    QStyle* style = option.widget->style();
    painter->setPen(Qt::lightGray);
    style->drawControl(QStyle::CE_ItemViewItem, &option, painter, option.widget);

    int L = option.rect.x();
    int B = option.rect.y();
    int H = option.rect.height();
    if(item && index.column()!=1)painter->drawLine(L,B,L,B+H);       //! | left vertical line
    if(item && isSeparator)
    {
        int W = qobject_cast<QTreeView *>(this->parent())->size().width();

        //QBrush bg = option.palette.dark();
        //option.rect.setWidth(W);
        //painter->fillRect(option.rect,bg);

        painter->drawLine(L,B,L+W,B);       //! ________    long bottom line
        painter->drawLine(L,B+H,L+W,B+H);   //! --------    long top line
        painter->drawLine(L+W,B,L+W,B+H);   //! | left vertical line
    }
    else if(item && !isSeparator)
    {
        int W = option.rect.width();
        painter->drawLine(L,B,L+W,B);       //! ________    linea bassa corta
        painter->drawLine(L,B+H,L+W,B+H);   //! --------    linea alta  corta
        painter->drawLine(L+W,B,L+W,B+H);   //! | a sinistra
    }

    //! comment: drawing a rectangle would be easier but
    //! the separator row has the central segment which
    //! splits the single column into two columns

    if(item && isSeparator)
    {
        //option.palette.setColor(QPalette::Text, option.palette.color(QPalette::BrightText))
        option.font.setBold(true);
        painter->setPen(Qt::black);
        //option.state &= ~QStyle::State_Selected;
        option.state &= QStyle::State_Selected;
    }
    option.state &= ~QStyle::State_HasFocus;
    QStyledItemDelegate::paint(painter, option, index);
}

//! ---------------------
//! function: setContext
//! details:
//! ---------------------
void GeneralDelegate::setContext(const occHandle(AIS_InteractiveContext) &aCTX)
{
    cout<<"GeneralDelegate::setContext()->____function called____"<<endl;
    myCTX = aCTX;
    if(myCTX.IsNull()==false) cout<<"GeneralDelegate::setContext()->____context OK____"<<endl;
}

//! -------------------------
//! function: setMeshContext
//! details:
//! -------------------------
void GeneralDelegate::setMeshContext(const occHandle(AIS_InteractiveContext) &aMeshCTX)
{
    cout<<"GeneralDelegate::setMeshContext()->____function called____"<<endl;
    myMeshCTX = aMeshCTX;
    if(myCTX.IsNull()==false) cout<<"GeneralDelegate::setMeshContext()->____mesh context OK____"<<endl;
}

//! --------------------------------------------
//! function: commitAndCloseShapeSelectorEditor
//! details:
//! --------------------------------------------
void GeneralDelegate::commitAndCloseShapeSelectorEditor()
{
    ShapeSelector *editor = qobject_cast<ShapeSelector *>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit scopeChanged();

    if(this->getCurrentNode()->getType()==SimulationNodeClass::nodeType_coordinateSystem)
    {
        //emit originChanged();
        emit originAndDirectionChanged();
    }
}

//! -------------------------------------------------------------------
//! function: commitAndCloseShapeSelectorEditor_boundaryPrismaticLayer
//! details:
//! -------------------------------------------------------------------
void GeneralDelegate::commitAndCloseShapeSelectorEditor_boundaryPrismaticLayer()
{
    static int i;
    cout<<"commitAndCloseShapeSelectorEditor_boundaryPrismaticLayer()____function called: "<<i++<<"____"<<endl;
    ShapeSelector *editor = qobject_cast<ShapeSelector *>(sender());
    if(editor!=Q_NULLPTR)
    {
        ShapeSelector *editor = qobject_cast<ShapeSelector *>(sender());
        emit commitData(editor);
        emit closeEditor(editor);
        emit boundaryScopeChanged();
    }

    else
    {
        QComboBox *editor = qobject_cast<QComboBox *>(sender());
        emit commitData(editor);
        emit closeEditor(editor);
        emit boundaryScopeChanged();
    }
}

//! ---------------------------------
//! function: commitAndCloseLineEdit
//! details:  no signal emitted
//! ---------------------------------
void GeneralDelegate::commitAndCloseLineEdit()
{
    QLineEdit *editor = qobject_cast<QLineEdit *>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
}

//! ---------------------------------------------
//! function: commitAndCloseLineEditTransparency
//! details:
//! ---------------------------------------------
static int s=0;
void GeneralDelegate::commitAndCloseLineEditTransparency()
{
    cout<<"commitAndCloseLineEditTransparency()____function called: "<<s++<<"____"<<endl;

    QLineEdit *editor = qobject_cast<QLineEdit*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);

    //! -------------------------------------
    //! signal emitted 3 times! Please check
    //! -------------------------------------
    if(s==2) emit transparencyChanged();
}

//! --------------------------------------------------------------------
//! function: commitAndCloseComboBox
//! details:  this closes a combobox editor without emitting any signal
//! --------------------------------------------------------------------
void GeneralDelegate::commitAndCloseComboBox()
{
    QComboBox *editor = qobject_cast<QComboBox *>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
}

//! ----------------------------------------
//! function: commitAndCloseComboBoxEmitter
//! details:
//! ----------------------------------------
void GeneralDelegate::commitAndCloseComboBoxEmitter()
{
    QComboBox *editor = qobject_cast<QComboBox *>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit EmitterStatusChanged();
}

//! -------------------------------------------
//! function: commitAndCloseMeshMetricCombobox
//! details:
//! -------------------------------------------
void GeneralDelegate::commitAndCloseMeshMetricCombobox()
{
    QComboBox *editor = qobject_cast<QComboBox *>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit meshMetricChanged();
}

//! -------------------------------------------
//! function: commitAndCloseVisibilityComboBox
//! details:
//! -------------------------------------------
void GeneralDelegate::commitAndCloseVisibilityComboBox()
{
    static int i;
    cout<<"GeneralDelegate::commitAndVisibilityCloseComboBox()->____function called: "<<i++<<"____"<<endl;
    QComboBox *editor = qobject_cast<QComboBox *>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit visibilityChanged_();
}

//! ------------------------------------------------
//! function: commitAndCloseComboBox of "Define by"
//! details:
//! ------------------------------------------------
void GeneralDelegate::commitAndCloseDefineByControlComboBox()
{
    QComboBox *editor = qobject_cast<QComboBox *>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit defineByChanged();
    switch(this->getCurrentNode()->getType())
    {
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Acceleration:
        emit accelerationChanged();
        break;
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Moment:
        emit momentChanged();
        break;
    }
}

//! ------------------------------------------------------
//! function: commitAndCloseComboBox of "Element control"
//! details:
//! ------------------------------------------------------
void GeneralDelegate::commitAndCloseElementControlComboBox()
{
    //static int i;
    //cout<<"GeneralDelegate::commitAndCloseElementControlComboBox()->____Element control changed: "<<i++<<"____"<<endl;
    QComboBox *editor = qobject_cast<QComboBox*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit elementControlChanged();
}

//! -----------------------------------------------------
//! function: commitAndCloseScopingMethodComboBox
//! details:  commitAndCloseComboBox of "Scoping method"
//! -----------------------------------------------------
void GeneralDelegate::commitAndCloseScopingMethodComboBox()
{
    //static int i;
    //cout<<"GeneralDelegate::commitAndCloseScopingMethodComboBox()->____function called: "<<i++<<"____"<<endl;
    QComboBox *editor = qobject_cast<QComboBox*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit scopingMethodChanged();
}

//! ----------------------------------------------
//! function: commitAndCloseBoundaryScopingMethod
//! details:
//! ----------------------------------------------
void GeneralDelegate::commitAndCloseBoundaryScopingMethod()
{
    QComboBox *editor = qobject_cast<QComboBox*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit boundaryScopingMethodChanged();
}

//! ---------------------------------------------
//! function: commitAndCloseNSSelector
//! details:  close the names selection selector
//! ---------------------------------------------
void GeneralDelegate::commitAndCloseNSSelector()
{
    QComboBox *editor = qobject_cast<QComboBox *>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit scopeChanged();
}

//! -----------------------------------------------
//! function: commitAndCloseCSSelector
//! details:  close the coordinate system selector
//! -----------------------------------------------
void GeneralDelegate::commitAndCloseCSSelector()
{
    QComboBox *editor = qobject_cast<QComboBox *>(sender());
    emit commitData(editor);
    //emit scopeChanged();
    emit closeEditor(editor);
    switch(this->getCurrentNode()->getType())
    {
    case SimulationNodeClass::nodeType_remotePoint: emit remotePointSystemOfReferenceChanged(); break;
    case SimulationNodeClass::nodeType_structuralAnalysisBoltPretension: emit boltCSChanged(); break;
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Acceleration: emit accelerationChanged(); break;
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Moment: emit momentChanged(); break;
    }
}

//! ------------------------------------
//! function: NumberOfStepSpinBoxClosed
//! details:
//! ------------------------------------
void GeneralDelegate::NumberOfStepSpinBoxClosed()
{
    //static int i;
    //cout<<"GeneralDelegate::NumberOfStepSpinBox()->____Number of steps changed____"<<i++<<endl;
    emit numberOfStepChanged();
}

//! -----------------------------------
//! function: CurrentStepSpinBoxClosed
//! details:
//! -----------------------------------
void GeneralDelegate::CurrentStepSpinBoxClosed()
{
    emit currentStepNumberChanged();
}

//! -------------------------------------------
//! function: commitAndCloseSolverTypeComboBox
//! details:
//! -------------------------------------------
void GeneralDelegate::commitAndCloseSolverTypeComboBox()
{
    QComboBox *editor = qobject_cast<QComboBox *>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit solverTypeChanged();
}

//! ------------------------------------
//! function: StepEndTimeLineEditClosed
//! details:
//! ------------------------------------
void GeneralDelegate::StepEndTimeLineEditClosed()
{
    QLineEdit *editor = qobject_cast<QLineEdit *>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    //static int i;
    //cout<<"GeneralDelegate::StepEndTimeLineEditClosed: "<<i++<<endl;
    emit StepEndTimeChanged();
}

/*
//! -----------------------------
//! function: commitAndCloseLe()
//! details:
//! -----------------------------
void GeneralDelegate::commitAndCloseLe()
{
    //static int i;
    //cout<<"GeneralDelegate::commitAndCloseLe(): "<<i++<<endl;
    LineEdit *editor = qobject_cast<LineEdit *>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit loadDefinitionChanged();
}
*/

//! ----------------------------------------------------
//! function: commitAndCloseModelChangeActivationStatus
//! details:
//! ----------------------------------------------------
void GeneralDelegate::commitAndCloseModelChangeActivationStatus()
{
    cout<<"GeneralDelegate::commitAndCloseModelChangeActivationStatus()->____function called____"<<endl;
    QComboBox *editor = qobject_cast<QComboBox *>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit modelChangeActivationStatusChanged();
}

//! -------------------------------------------
//! function: commitAndCloseLe_filmCoefficient
//! details:
//! -------------------------------------------
void GeneralDelegate::commitAndCloseLe_filmCoefficient()
{
    LineEdit *editor = qobject_cast<LineEdit *>(sender());
    QString text = editor->getText();
    emit commitData(editor);
    emit closeEditor(editor);
    emit loadDefinitionFilmCoefficientChanged(text);
}

//! ------------------------------------------------
//! function: commitAndCloseLe_referenceTemperature
//! details:
//! ------------------------------------------------
void GeneralDelegate::commitAndCloseLe_referenceTemperature()
{
    LineEdit *editor = qobject_cast<LineEdit *>(sender());
    QString text = editor->getText();
    emit commitData(editor);
    emit closeEditor(editor);
    emit loadDefinitionReferenceTemperatureChanged(text);
}

//! -------------------------------------
//! function: commitAndCloseLe_magnitude
//! details:
//! -------------------------------------
void GeneralDelegate::commitAndCloseLe_magnitude()
{
    LineEdit *editor = qobject_cast<LineEdit *>(sender());
    QString text = editor->getText();
    //emit loadDefinitionMagnitudeChanged(text);
    emit commitData(editor);
    emit closeEditor(editor);
    emit loadDefinitionMagnitudeChanged(text);  //! 5/11/2019

    switch(this->getCurrentNode()->getType())
    {
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Acceleration:
        emit accelerationChanged();
        break;
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Moment:
        emit momentChanged();
        break;
    }
}

//! --------------------------------------
//! function: commitAndCloseLe_Xcomponent
//! details:
//! --------------------------------------
void GeneralDelegate::commitAndCloseLe_Xcomponent()
{
    LineEdit *editor = qobject_cast<LineEdit *>(sender());
    QString text = editor->getText();
    emit commitData(editor);
    emit closeEditor(editor);
    emit loadDefinitionXChanged(text);
    switch(this->getCurrentNode()->getType())
    {
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Acceleration:
        emit accelerationChanged();
        break;
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Moment:
        emit momentChanged();
        break;
    }
}

//! --------------------------------------
//! function: commitAndCloseLe_Ycomponent
//! details:
//! --------------------------------------
void GeneralDelegate::commitAndCloseLe_Ycomponent()
{
    LineEdit *editor = qobject_cast<LineEdit *>(sender());
    QString text = editor->getText();
    emit commitData(editor);
    emit closeEditor(editor);
    emit loadDefinitionYChanged(text);
    switch(this->getCurrentNode()->getType())
    {
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Acceleration:
        emit accelerationChanged();
        break;
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Moment:
        emit momentChanged();
        break;
    }
}

//! --------------------------------------
//! function: commitAndCloseLe_Zcomponent
//! details:
//! --------------------------------------
void GeneralDelegate::commitAndCloseLe_Zcomponent()
{
    LineEdit *editor = qobject_cast<LineEdit *>(sender());
    QString text = editor->getText();
    emit commitData(editor);
    emit closeEditor(editor);
    emit loadDefinitionZChanged(text);
    switch(this->getCurrentNode()->getType())
    {
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Acceleration:
        emit accelerationChanged();
        break;
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Moment:
        emit momentChanged();
        break;
    }
}

//! -------------------------------------------
//! function: commitAndCloseDirectionSelector
//! details:
//! -------------------------------------------
void GeneralDelegate::commitAndCloseDirectionSelector()
{
    DirectionSelector *editor = qobject_cast<DirectionSelector *>(sender());
    emit commitData(editor);
    emit closeEditor(editor);

    switch(this->getCurrentNode()->getType())
    {
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Acceleration:
        emit accelerationChanged();
        break;
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Moment:
        emit momentChanged();
        break;
    }
}

void GeneralDelegate::commitAndCloseComboBoxDefineBy_()
{
    static int i;
    cout<<"GeneralDelegate::commitAndCloseComboBoxDefineBy_: "<<i++<<endl;
    QComboBox *editor = qobject_cast<QComboBox *>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit defineBy_Changed();
}

void GeneralDelegate::emitOriginChanged()
{
    emit originChanged();
}

void GeneralDelegate::emitOriginAndDirectionChanged()
{
    emit originAndDirectionChanged();
}

//! -------------------------------------------
//! function: emitRemotePointChangedByLocation
//! details:
//! -------------------------------------------
void GeneralDelegate::emitRemotePointChangedByLocation()
{
    emit remotePointChangedByLocation(); 
}

//! -----------------------------------------
//! function: commitAndCloseLineEditCSOrigin
//! details:
//! -----------------------------------------
void GeneralDelegate::commitAndCloseLineEditCSOrigin()
{
    static int i;
    cout<<"GeneralDelegate::commitAndCloseLineEditCSOrigin()->____function called: "<<i++<<endl;
    QLineEdit *editor = qobject_cast<QLineEdit*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
}

//! ------------------------------------------------
//! function: commitAndCloseLineEdiCSTransformation
//! details:
//! ------------------------------------------------
void GeneralDelegate::commitAndCloseLineEdiCSTransformation()
{
    static int i;
    cout<<"GeneralDelegate::commitAndCloseLineEdiCSTransformation()->____function called: "<<i++<<"____"<<endl;
    QLineEdit *editor = qobject_cast<QLineEdit *>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    // this is a workaround, since this function is called three (3!) times [???]
    if(i%3==2)
        emit transformationsChanged();
}

void GeneralDelegate::commitAndCloseAutoTimeStepping()
{
    static int i;
    cout<<"GeneralDelegate::commitAndCloseAutoTimeStepping()->____function called: "<<i++<<"____"<<endl;
    QComboBox *editor = qobject_cast<QComboBox *>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit autoTimeSteppingChanged();
}

void GeneralDelegate::commitAndCloseSubstepEditors()
{
    static int i;
    cout<<"GeneralDelegate::commitAndCloseSubstepEditors()->____function called: "<<i++<<"____"<<endl;
    QComboBox *editor = qobject_cast<QComboBox *>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit timeDivisionChanged();
}

void GeneralDelegate::commitAndCloseShowMeshNodes()
{
    QComboBox *editor = qobject_cast<QComboBox *>(sender());
    emit commitData(editor);
    bool meshNodesVisible = editor->currentIndex()==0? false:true;
    emit closeEditor(editor);
    emit meshNodesVisibilityChanged(meshNodesVisible);
}

void GeneralDelegate::commitAndCloseSmoothingControl()
{
    QComboBox *editor = qobject_cast<QComboBox *>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit meshSmoothingChanged();
}

void GeneralDelegate::commitAndCloseStepSelectionMode()
{
    QComboBox *editor = qobject_cast<QComboBox *>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit stepSelectionModeChanged();
}

void GeneralDelegate::commitAndCloseGenerate()
{
    cout<<"GeneralDelegate::commitAndCloseGenerate()->____function called____"<<endl;
    QPushButton *editor = qobject_cast<QPushButton*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit requestStartInterpolator();
}

void GeneralDelegate::commitAndCloseSubmeshesControl()
{
    QComboBox *editor = qobject_cast<QComboBox*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
}

void GeneralDelegate::commitAndCloseRemapFlagControl()
{
    QComboBox *editor = qobject_cast<QComboBox*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit remapFlagChanged();
}

/*
void GeneralDelegate::commitAndCloseSurfaceMesherSelector()
{
    QComboBox *editor = qobject_cast<QComboBox*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit globalMeshControlChanged();
}
*/

void GeneralDelegate::commitAndCloseMidsideNodesSelector()
{
    QComboBox *editor = qobject_cast<QComboBox*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit globalMeshControlChanged();
}

void GeneralDelegate::commitAndCloseInitialSizeSeedSelector()
{
    QComboBox *editor = qobject_cast<QComboBox*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit globalMeshControlChanged();
}

void GeneralDelegate::commitAndCloseRelevanceControl()
{
    QSlider *editor = qobject_cast<QSlider*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit globalMeshControlChanged();
}

void GeneralDelegate::commitAndCloseSplitFileModeSelector()
{
    QComboBox *editor = qobject_cast<QComboBox*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    cout<<editor->currentIndex()<<endl;
    emit splitModeChanged();
}
/*
void GeneralDelegate::commitAndCloseTranslate()
{
    cout<<"GeneralDelegate::commitAndCloseTranslate()->____function called____"<<endl;
    QPushButton *editor = qobject_cast<QPushButton*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit requestStartOpenFoamScalarDataTranslator();
}
*/
//! ----------------------------------------------------
//! function: commitAndCloseBackgroundChanged()
//! details:  first or second color or gradient changed
//! ----------------------------------------------------
void GeneralDelegate::commitAndCloseBackgroundControl()
{
    cout<<"GeneralDelegate::commitAndCloseBackgroundControl()->____function called____"<<endl;
    QComboBox *editor = qobject_cast<QComboBox*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit backgroundChanged();
}

//! -------------------------------------------
//! function: commitAndCloseCouplingSelector()
//! details:  handle "Kinematic" "Distributed"
//! -------------------------------------------
void GeneralDelegate::commitAndCloseCouplingSelector()
{
    QComboBox *editor = qobject_cast<QComboBox*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit couplingChanged();
}

//! -----------------------------------------------
//! function: commitAndCloseRPDOFselector
//! details:  handle "Manual" "Program controlled"
//! -----------------------------------------------
void GeneralDelegate::commitAndCloseRPDOFselector()
{
    QComboBox *editor = qobject_cast<QComboBox*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
}

//! ------------------------------------------
//! function: commitAndCloseDOFswitchSelector
//! details:  handle "Inactive" "Active"
//! ------------------------------------------
void GeneralDelegate::commitAndCloseDOFswitchSelector()
{
    cout<<"GeneralDelegate::commitAndCloseDOFswitchSelector()->____function called____"<<endl;
    QComboBox *editor = qobject_cast<QComboBox*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit DOFselectorChanged();
}

//! ----------------------------------
//! function: commitAndCloseNbThreads
//! details:
//! ----------------------------------
void GeneralDelegate::commitAndCloseNbThreads()
{
    QSpinBox *editor = qobject_cast<QSpinBox*>(sender());
    if(editor!=NULL) cerr<<"____editor OK____"<<endl;
    emit commitData(editor);

    //! do not close the editor. Commit only
    emit NbThreadsChanged();
}

//! ------------------------------------------------------
//! function: commitAndCloseSolutionInformation
//! details:  handle the solution information type switch
//! ------------------------------------------------------
void GeneralDelegate::commitAndCloseSolutionInformation()
{
    cout<<"GeneralDelegate::commitAndCloseSolutionInformation()->____function called____"<<endl;
    QComboBox *editor = qobject_cast<QComboBox*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit solutionInformationChanged();
}

//! ----------------------------------------------------
//! function: commitAndCloseSolutionComponentSelector()
//! details:  postprocessing: component changed
//! ----------------------------------------------------
void GeneralDelegate::commitAndCloseSolutionComponentSelector()
{
    cout<<"GeneralDelegate::commitAndCloseSolutionComponentSelector()->____function called____"<<endl;
    QComboBox *editor = qobject_cast<QComboBox*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit solutionComponentChanged();
}

//! --------------------------------------------------
//! function: commitAndCloseBySelector
//! details:  postprocessing: solution by time or set
//! --------------------------------------------------
void GeneralDelegate::commitAndCloseBySelector()
{
    cout<<"GeneralDelegate::commitAndCloseBySelector()->____function called____"<<endl;
    QComboBox *editor = qobject_cast<QComboBox*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);

    switch(this->getCurrentNode()->getType())
    {
    case SimulationNodeClass::nodeType_solutionThermalFlux:
    case SimulationNodeClass::nodeType_solutionThermalTemperature:
    {
        emit byChanged();
    }
        break;

    case SimulationNodeClass::nodeType_solutionStructuralEquivalentPlasticStrain:
    case SimulationNodeClass::nodeType_solutionStructuralMechanicalStrain:
    case SimulationNodeClass::nodeType_solutionStructuralNodalDisplacement:
    case SimulationNodeClass::nodeType_solutionStructuralStress:
    case SimulationNodeClass::nodeType_solutionStructuralTemperature:
    case SimulationNodeClass::nodeType_solutionStructuralThermalStrain:
    case SimulationNodeClass::nodeType_solutionStructuralTotalStrain:
    {
        emit byChanged();
    }
        break;

    case SimulationNodeClass::nodeType_meshMethod:
    {
        emit meshSimplificationByChanged();
    }
        break;
    }
}

//! ---------------------------------------------
//! function: closeComboBox
//! details:  this is used by "Large deflection"
//! ---------------------------------------------
void GeneralDelegate::closeComboBox()
{
    //cout<<"GeneralDelegate::commitAndCloseComboBox()->____function called____"<<endl;
    QComboBox *editor = qobject_cast<QComboBox*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
}

//! ---------------------------------------------------------
//! function: commitAndCloseContactToleranceEditor
//! details:  closes the editor for the contact tolerance
//!           parameter (automatic contact search algorithm)
//! ---------------------------------------------------------
void GeneralDelegate::commitAndCloseContactToleranceEditor()
{
    //cout<<"commitAndCloseContactToleranceEditor()->____function called____"<<endl;
    QLineEdit *editor = qobject_cast<QLineEdit*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit automaticContactSearchToleranceChanged();
}

//! ---------------------------------------------------
//! function: commitAndCloseComboBoxTypeOfSizing
//! details:  closes the editor for the type of sizing
//!           (mesh controls: division/element size...
//! ---------------------------------------------------
void GeneralDelegate::commitAndCloseComboBoxTypeOfSizing()
{
    //cout<<"GeneralDelegate::commitAndCloseComboBoxTypeOfSizing()->____function called____"<<endl;
    QComboBox *editor = qobject_cast<QComboBox*>(sender());
    emit commitData(editor);

    //! ----------------------------------------------------------
    //! important note: emit typeOfSizingChanged() before closing
    //! the editor, for synch reasons. Otherwise the program will
    //! crash
    //! ----------------------------------------------------------
    emit typeOfSizingChanged();

    emit closeEditor(editor);
}

//! -------------------------------------------------------------
//! function: commitAndCloseStraightSidedElementsControl
//! details:  control for the generation of super/iso-parametric
//!           elements
//! -------------------------------------------------------------
void GeneralDelegate::commitAndCloseStraightSidedElementsControl()
{
    QComboBox *editor = qobject_cast<QComboBox*>(sender());
     emit commitData(editor);
    emit closeEditor(editor);
    emit globalMeshControlChanged();
}

//! -------------------------------------------
//! function: commitAndCloseUpdateInterval
//! details:  solution information update rate
//! -------------------------------------------
void GeneralDelegate::commitAndCloseUpdateInterval()
{
    QLineEdit *editor = qobject_cast<QLineEdit*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit solutionInformationUpdateIntervalChanged();
}

//! --------------------------------------------
//! function: commitAndCloseBoltStatusDefinedBy
//! details:  bolt status definition changed
//! --------------------------------------------
void GeneralDelegate::commitAndCloseBoltStatusDefinedBy()
{
    QComboBox *editor = qobject_cast<QComboBox*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    //! do not move signal "BoltStatusDefinedByChanged()" from this position
    emit BoltStatusDefinedByChanged();
}

//! ----------------------------------------
//! function: commitAndCloseBoltLoadControl
//! details:
//! ----------------------------------------
void GeneralDelegate::commitAndCloseBoltLoadControl()
{
    QLineEdit *editor = qobject_cast<QLineEdit*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit BoltLoadChanged();
}

//! ---------------------------------------
//! function: commitAndCloseBoltAdjustment
//! details:
//! ---------------------------------------
void GeneralDelegate::commitAndCloseBoltAdjustmentControl()
{
    QLineEdit *editor = qobject_cast<QLineEdit*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit BoltAdjustmentChanged();
}

//! -------------------------------------------
//! function: commitAndCloseBoxFluxConvergence
//! details:
//! -------------------------------------------
void GeneralDelegate::commitAndCloseComboBoxFluxConvergence()
{
    QComboBox *editor = qobject_cast<QComboBox*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit fluxConvegernceChanged();
}

void GeneralDelegate::commitAndCloseComboBoxSolutionConvergence()
{
    QComboBox *editor = qobject_cast<QComboBox*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit solutionConvegernceChanged();
}

//! ----------------------------------------------------------------------
//! function: commitAndCloseFieldParametersLineEditEditor()
//! details:  this closes the QLineEdit for changing the field parameters
//! ----------------------------------------------------------------------
void GeneralDelegate::commitAndCloseFieldParametersLineEditEditor()
{
    QLineEdit *editor = qobject_cast<QLineEdit*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit fieldParametersChanged();
}

//! ------------------------------------------------------------------------------------
//! function: commitAndCloseTimeIncrementationEditor()
//! details:  this closes the QLineEdit for changing the time incrementation parameters
//! ------------------------------------------------------------------------------------
void GeneralDelegate::commitAndCloseTimeIncrementationEditor()
{
    QLineEdit *editor = qobject_cast<QLineEdit*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit timeIncrementationChanged();
}

//! --------------------------------------------------------------
//! function: commitAndCloseTimeIncrementationControl()
//! details:  this closes the "Time incrementation" type combobox
//! --------------------------------------------------------------
void GeneralDelegate::commitAndCloseTimeIncrementationControl()
{
    QComboBox *editor = qobject_cast<QComboBox*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit timeIncrementationChanged();
}

//! -------------------------------------------------------------------
//! function: commitAndCloseCutBackParametersLineEdit()
//! details:  this closes one the line edit editor for cutback factors
//! -------------------------------------------------------------------
void GeneralDelegate::commitAndCloseCutBackParametersLineEdit()
{
    QLineEdit *editor = qobject_cast<QLineEdit*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit cutbackFactorsChanged();
}

//! ---------------------------------------------------
//! function: commitAndCloseCutBackEditor
//! details:  this closes the "Cutback factors" editor
//! ---------------------------------------------------
void GeneralDelegate::commitAndCloseCutBackEditor()
{
    QComboBox *editor = qobject_cast<QComboBox*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit cutbackFactorsChanged();
}

//! -----------------------------------------------
//! function: commitAndCloseLineSearchEditor()
//! details:  this closes the "Line search" editor
//! -----------------------------------------------
void GeneralDelegate::commitAndCloseLineSearchEditor()
{
    QComboBox *editor = qobject_cast<QComboBox*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit lineSearchChanged();
}

//! ---------------------------------------------
//! function: commitAndCloseLineSearchEditor()
//! details:  this closes the line search editor
//! ---------------------------------------------
void GeneralDelegate::commitAndCloseLineParametersChanged()
{
    QLineEdit *editor = qobject_cast<QLineEdit*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit lineSearchParametersChanged();
}

//! --------------------------------------------------
//! function: commitAndCloseSmallSlidingControl()
//! details:  this closes the "Small sliding" control
//! --------------------------------------------------
void GeneralDelegate::commitAndCloseSmallSlidingControl()
{
    QComboBox *editor = qobject_cast<QComboBox*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
}

//! --------------------------------------------------
//! function: commitAndCloseAdjustControl()
//! details:  this closes the "Small sliding" control
//! --------------------------------------------------
void GeneralDelegate::commitAndCloseAdjustControl()
{
    QComboBox *editor = qobject_cast<QComboBox*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
}

//! ---------------------------------------------------------------------
//! function: commitAndCloseComboBoxForOutputSettings()
//! details:  this closes one of the combobox among of "Output settings"
//! ---------------------------------------------------------------------
void GeneralDelegate::commitAndCloseComboBoxForOutputSettings()
{
    QComboBox *editor = qobject_cast<QComboBox*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit outputControlsChanged();
}

//! --------------------------------------------------
//! function: commitAndCloseComboBoxStoreResultsAt()
//! details:  this closes the "Small sliding" control
//! --------------------------------------------------
void GeneralDelegate::commitAndCloseStoreResultsAt()
{
    QComboBox *editor = qobject_cast<QComboBox*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit storeResultsAtChanged();
}

//! ----------------------------------------------------
//! function: commitAndCloseRemotePointsSelector()
//! details:  this closes the selector of remote points
//! ----------------------------------------------------
void GeneralDelegate::commitAndCloseRemotePointsSelector()
{
    QComboBox *editor = qobject_cast<QComboBox*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit scopeChanged();
}

//! -------------------------------------------------------
//! function: commitAndCloseScaleTypeSelector()
//! details:  this closes the scale selector (result item)
//! -------------------------------------------------------
void GeneralDelegate::commitAndCloseScaleTypeSelector()
{
    //cout<<"GeneralDelegate::commitAndCloseScaleTypeSelector()->____function called____"<<endl;
    QComboBox *editor = qobject_cast<QComboBox*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit colorBoxScaleChanged();
}

//! --------------------------------------------------------------
//! function: commitAndCloseNumberOfInterval()
//! details:  this closes "# intervals" control (color scale in a
//!           result item)
//! --------------------------------------------------------------
void GeneralDelegate::commitAndCloseNumberOfInterval()
{
    QSpinBox *editor = qobject_cast<QSpinBox*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit colorBoxScaleChanged();
}

//! --------------------------------------------------
//! function: commitAndCloseMinMaxControls
//! details:  this closes the "Min" or "Max" controls
//!           (color scale in a result item)
//! --------------------------------------------------
void GeneralDelegate::commitAndCloseMinMaxControls()
{
    QLineEdit *editor = qobject_cast<QLineEdit*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit colorBoxScaleChanged();
}

void GeneralDelegate::commitAndCloseMeshMethodSelector()
{
    QComboBox *editor = qobject_cast<QComboBox*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit scopeChanged();
    emit meshMethodChanged();
}

//! -------------------------------------------------
//! function: commitAndCloseMeshDefeaturingControl()
//! details:  to be removed
//! -------------------------------------------------
void GeneralDelegate::commitAndCloseMeshDefeaturingControl()
{
    QComboBox *editor = qobject_cast<QComboBox*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit scopeChanged();        // why? check ... to do...
    emit meshDefeaturingChanged();
}

//! to be removed
void GeneralDelegate::commitAndCloseMeshSimplificationControl()
{
    QComboBox *editor = qobject_cast<QComboBox*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit scopeChanged();        // why? check ... to do...
    emit meshSimplificationChanged();
}

void GeneralDelegate::commitAndCloseMeshHealingControl()
{
    QComboBox *editor = qobject_cast<QComboBox*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit scopeChanged();        // why? check ... to do...
    emit meshHealingChanged();
}

//! -------------------------------------------------------------
//! function: commitAndCloseMeshDefeaturingParameterValueControl
//! details:  mesh defeaturing parameters changed
//! -------------------------------------------------------------
void GeneralDelegate::commitAndCloseMeshDefeaturingParameterValueControl()
{
    QWidget *widget;
    widget = qobject_cast<QLineEdit*>(sender());
    if(widget!=Q_NULLPTR)
    {
        QLineEdit *editor = qobject_cast<QLineEdit*>(sender());
        emit commitData(editor);
        emit closeEditor(editor);
    }
    widget = qobject_cast<QSlider*>(sender());
    if(widget!=Q_NULLPTR)
    {
        QSlider *editor = qobject_cast<QSlider*>(sender());
        emit commitData(editor);
        emit closeEditor(editor);
    }
    emit meshSimplificationParameterValueChanged();
}

//! ----------------------------------------------------------------
//! function: commitAndClosePrismaticLayerOptions()
//! details:  close the editor for "Options" in boundary layer item
//! ----------------------------------------------------------------
void GeneralDelegate::commitAndClosePrismaticLayerOptions()
{
    QComboBox *editor = qobject_cast<QComboBox*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit prismaticLayerOptionsChanged();
}

//! --------------------------------
//! function: commitAndCloseUseBRep
//! details:
//! --------------------------------
void GeneralDelegate::commitAndCloseUseBRep()
{
    cout<<"GeneralDelegate::commitAndCloseUseBRep()->____function called____"<<endl;
    QComboBox *editor = qobject_cast<QComboBox*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit BRepFlagChanged();
}

//! -------------------------------------------
//! function: commitAndCloseDefeaturingControl
//! details:
//! -------------------------------------------
void GeneralDelegate::commitAndCloseDefeaturingControl()
{
    QComboBox *editor = qobject_cast<QComboBox*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit defeaturingFlagChanged();
}

//! ----------------------------------------------
//! function: commitAndCloseSimplificationControl
//! details:
//! ----------------------------------------------
void GeneralDelegate::commitAndCloseSimplificationControl()
{
    QComboBox *editor = qobject_cast<QComboBox*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit simplificationFlagChanged();
}

//! ------------------------------------
//! function: commitAndCloseTessellator
//! details:
//! ------------------------------------
void GeneralDelegate::commitAndCloseTessellator()
{
    QComboBox *editor = qobject_cast<QComboBox*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit TessellatorChanged();
}

//! -------------------------------------
//! function: commitAndCloseVolumeMesher
//! details:
//! -------------------------------------
void GeneralDelegate::commitAndCloseVolumeMesher()
{
    QComboBox *editor = qobject_cast<QComboBox*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit volumeMesherChanged();
}

//! ----------------------------------------
//! function: commitAndCloseGeometryHealing
//! details:
//! ----------------------------------------
void GeneralDelegate::commitAndCloseGeometryHealing()
{
    QComboBox *editor = qobject_cast<QComboBox*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit geometryHealingChanged();
}

#include <QToolTip>
void GeneralDelegate::displayLevelSliderValue(int pos)
{
    double level = MIN_MESH_REDUCTION+pos*(MAX_MESH_REDUCTION-MIN_MESH_REDUCTION)/100.0;
    cout<<"____slider value: "<<level<<"____"<<endl;
}

//! -----------------------------------
//! function: commitAndCloseFileSelect
//! details:
//! -----------------------------------
void GeneralDelegate::commitAndCloseFileSelect()
{
    QFileSelect *editor = qobject_cast<QFileSelect*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
}

//! --------------------------------------
//! function: commitAndCloseToleranceType
//! details:  geometry healing section
//! --------------------------------------
void GeneralDelegate::commitAndCloseToleranceType()
{
    QComboBox *editor = qobject_cast<QComboBox*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit geometryHealingChanged();
}

#ifdef COSTAMP_VERSION
//! ------------------------------------------------
//! function: commitAndCloseTimeStepBuilderButton()
//! details:
//! ------------------------------------------------
void GeneralDelegate::commitAndCloseTimeStepBuilderButton()
{
    cout<<"GeneralDelegate::commitAndCloseTimeStepBuilder()->____function called____"<<endl;
    QPushButton *editor = qobject_cast<QPushButton*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit requestStartTimeStepBuilder();
}

void GeneralDelegate::commitAndCloseTimeStepBuilderFileSelector()
{
    QFileSelect *editor = qobject_cast<QFileSelect*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
}
#endif

//! -------------------------------------
//! TetWild parameters - envelope sizing
//! -------------------------------------
void GeneralDelegate::commitAndCloseEnvelopeSizing()
{
    QComboBox *editor = qobject_cast<QComboBox*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit envelopeSizingChanged();
}

//! -------------------------------------------
//! function: commitAndCloseIdealLength
//! details: TetWild parameters - ideal length
//! -------------------------------------------
void GeneralDelegate::commitAndCloseIdealLength()
{
    QComboBox *editor = qobject_cast<QComboBox*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit idealLengthChanged();
}

//! -----------------------------------------------
//! function: commitAndCloseInterpolationAlgorithm
//! details:
//! -----------------------------------------------
void GeneralDelegate::commitAndCloseInterpolationAlgorithm()
{
    QComboBox *editor = qobject_cast<QComboBox*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit interpolationAlgorithmChanged();
}

//! ---------------------------------
//! function: commitAndCloseItemType
//! details:
//! ---------------------------------
void GeneralDelegate::commitAndCloseItemType()
{
    QComboBox *editor = qobject_cast<QComboBox*>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit modelChangeScopingMethodChanged();
}

//! --------------------------------------------
//! function: commitAndCloseContactPairSelector
//! details:
//! --------------------------------------------
void GeneralDelegate::commitAndCloseContactPairSelector()
{
    QComboBox *editor = qobject_cast<QComboBox *>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
}

//! ----------------------------------------
//! function: commitAndCloseSelectionMethod
//! details:
//! ----------------------------------------
void GeneralDelegate::commitAndCloseSelectionMethod()
{
    QComboBox *editor = qobject_cast<QComboBox *>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit selectionMethodChanged();
}

//! --------------------------------------------
//! function: commitAndCloseLineEditElementList
//! details:
//! --------------------------------------------
void GeneralDelegate::commitAndCloseLineEditElementList()
{
    QLineEdit *editor = qobject_cast<QLineEdit *>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit elementListChanged();
}

//! ------------------------------------
//! function: commitAndCloseFatigueAlgo
//! details:
//! ------------------------------------
void GeneralDelegate::commitAndCloseFatigueAlgo()
{
    QComboBox *editor = qobject_cast<QComboBox *>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit fatigueAlgoChanged();
}

//! -------------------------------------
//! function: commitAndCloseAnalysisType
//! details:
//! -------------------------------------
void GeneralDelegate::commitAndCloseAnalysisType()
{
    QComboBox *editor = qobject_cast<QComboBox *>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit analysisTypeChangedChanged();
}

//! -----------------------------------------------
//! function: commitAndCloseTimeIntegrationChanged
//! details:
//! -----------------------------------------------
void GeneralDelegate::commitAndCloseTimeIntegrationChanged()
{
    QComboBox *editor = qobject_cast<QComboBox *>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
    emit timeIntegrationChanged();
}
