//! ----------------
//! custom includes
//! ----------------
#include "qextendedstandarditem.h"
#include "mydefines.h"
#include "postobject.h"
#include "ng_meshvs_datasource2d.h"
#include "ng_meshvs_datasource3d.h"
#include "tools.h"
#include "listofmesh.h"
#include "detailviewer.h"
#include "simulationmanager.h"
#include "customtablemodel.h"
#include "maintreetools.h"
#include "ccxsolvermessage.h"
#include "solutioninfo.h"

//! ----
//! OCC
//! ----
#include "ais_colorscaleextended.h"

//! ---
//! Qt
//! ---
#include <QDebug>
#include <QPalette>
#include <QTableWidget>

//! ----
//! C++
//! ----
#include <iostream>
using namespace std;

//! -----------------------------------------------------------------------//
//! function: setData                                                      //
//! details:  this method has been redefined in order to introcuce         //
//!           the automatic handling of the item icons (Qt::DisplayRole).  //
//!           when the item contains a SimulationNodeClass* as value       //
//!           (Qt::UserRole). Algorithm:                                   //
//!           1) check if the passed value can be converted into a         //
//!              SimulationNodeClass*                                      //
//!           2) check if the role is Qt::UserRole                         //
//!           3) if 1) && 2) and && the SuppressionStatus is "Suppressed"  //
//!              use the icon representing the suppression                 //
//!           4) if 1) && 2) && the SuppressionStatus is "No, Active"      //
//!              use a specific item icon                                  //
//! -----------------------------------------------------------------------//
void QExtendedStandardItem::setData(const QVariant &value, int role)
{
    if(value.canConvert<SimulationNodeClass*>() && role==Qt::UserRole)
    {
        //!myCurNode = value.value<SimulationNodeClass*>();
        //!QMainWindow *theMainWindow = qobject_cast<QMainWindow*>(myCurNode->parent()->parent()->parent());
        //!cout<<"QExtendedStandardItem::setData()->____function called on node____"<<endl;

        QStandardItem::setData(value,role);

        SimulationNodeClass::nodeType theNodeType = value.value<SimulationNodeClass*>()->getType();
        QVariant data;
        data.setValue(this->getIcon(theNodeType));
        QStandardItem::setData(data,Qt::DecorationRole);
    }
    else
    {
        //! other types of items which do not contain a SimulationNodeClass*
        QStandardItem::setData(value,role);
    }

    //! ----------------------------------------------------------
    //! post processing items - creation with the checkable role
    //! ----------------------------------------------------------
    if(value.canConvert<SimulationNodeClass*>() && role == Qt::UserRole)
    {
        SimulationNodeClass *node = value.value<SimulationNodeClass*>();
        SimulationNodeClass::nodeType theFamily = node->getFamily();
        SimulationNodeClass::nodeType theType = node->getType();
        //cout<<"family "<<theFamily<<" type: "<<theType<<endl;
        if(theFamily == SimulationNodeClass::nodeType_StructuralAnalysisSolution &&
                theType != SimulationNodeClass::nodeType_StructuralAnalysisSolutionInformation &&
                theType != SimulationNodeClass::nodeType_StructuralAnalysisSolution)
        {
            QExtendedStandardItem *item = node->getPropertyItem("Post object");
            if(item == NULL)
            {
                //cout<<"no post object has been inserted"<<endl;
                QStandardItem::setData(true,Qt::CheckStateRole);
                QStandardItem::setCheckState(Qt::Unchecked);
            }
            else
            {
                //cout<<"the item has a postObject: check if data are present or not"<<endl;
                bool isEmpty = item->data(Qt::UserRole).value<Property>().getData().value<postObject>().isEmpty();
                if(isEmpty) QStandardItem::setCheckState(Qt::Unchecked);
                else QStandardItem::setCheckState(Qt::Checked);
            }
        }
    }
}

//! ---------------
//! function: data
//! details:
//! ---------------
#include "handle_ais_doublearrowmarker_reg.h"
#include "meshvs_mesh_handle_reg.h"
#include "occhandle.h"
#include <indexedmapofmeshdatasources.h>

QVariant QExtendedStandardItem::data(int role) const
{
    if(role==Qt::DisplayRole)
    {
        QString name = data(Qt::UserRole).value<Property>().getName();
        QString name1 = data(Qt::UserRole+1).toString();

        QVariant data;

        //! ----------------------
        //! "Mass" "Jx" "Jy" "Jz"
        //! ----------------------
        if(name =="Mass" || name =="Jx" || name =="Jy" || name == "Jz")
        {
            double val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toDouble();
            data.setValue(val);
            return data;
        }
        //! -----------------
        //! Convergence data
        //! -----------------
        if(name=="Convergence data")
        {
            const QList<solutionInfo> solInfoList = QStandardItem::data(Qt::UserRole).value<Property>().getData().value<QList<solutionInfo>>();
            if(solInfoList.isEmpty()) data.setValue(QString("Empty"));
            else data.setValue(QString("%1 substeps").arg(solInfoList.length()));
            return data;
        }
        //! ------------------
        //! Time intergration
        //! ------------------
        if(name=="Static/Transient")
        {
            Property::timeIntegration val = QStandardItem::data(Qt::UserRole).value<Property>().getData().value<Property::timeIntegration>();
            switch(val)
            {
            case 0: data.setValue(QString("Static")); break;
            case 1: data.setValue(QString("Transient")); break;
            }
            return data;
        }
        //! --------------
        //! Analysis type
        //! --------------
        if(name=="Analysis type")
        {
            Property::analysisType val = QStandardItem::data(Qt::UserRole).value<Property>().getData().value<Property::analysisType>();
            switch(val)
            {
            case Property::analysisType_structural: data.setValue(QString("Structural")); break;
            case Property::analysisType_thermal: data.setValue(QString("Thermal")); break;
            case Property::analysisType_modal: data.setValue(QString("Modal")); break;
            case Property::analysisType_frequencyResponse: data.setValue(QString("Frequency response")); break;
            case Property::analysisType_uncoupledTemperatureDisplacement: data.setValue(QString("Uncoupled temperature displacement")); break;
            case Property::analysisType_coupledTemperatureDisplacement: data.setValue(QString("Coupled temperature displacement")); break;
            }
            return data;
        }
        //! --------------------
        //! "Discrete time map"
        //! --------------------
        if(name=="Discrete time map")
        {
            QMap<double,QVector<int>> dtm = QStandardItem::data(Qt::UserRole).value<Property>().getData().value<QMap<double,QVector<int>>>();
            if(dtm.isEmpty()) data.setValue(QString("Empty"));
            else data.setValue(QString("%1").arg(dtm.size()));
            return data;
        }
        //! ---------------
        //! "Solver output
        //! ---------------
        if(name =="Solver output")
        {
            CCXSolverMessage val = QStandardItem::data(Qt::UserRole).value<Property>().getData().value<CCXSolverMessage>();
            if(val.myText.isEmpty()) data.setValue(QString("No solver output"));
            else data.setValue(QString("Shown in worksheet"));
            return data;
        }
        //! -----------------------
        //! "Structural time step"
        //! -----------------------
        if(name =="Structural time step")
        {
            int val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toInt();
            data.setValue(QString("%1").arg(val));
            return data;
        }
        //! ---------------------
        //! "Coupling time step"
        //! ---------------------
        if(name =="Coupling time step")
        {
            int val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toInt();
            data.setValue(QString("%1").arg(val));
            return data;
        }
        //! ----------------------------
        //! "Imported body temperature"
        //! ----------------------------
        if(name =="Imported body temperature")
        {
            void *p = QStandardItem::data(Qt::UserRole).value<Property>().getData().value<void*>();
            QStandardItem *item = (QStandardItem*)(p);
            if(item!=Q_NULLPTR) return item->data(Qt::DisplayRole);
            data.setValue(QString("No body temperature item selected"));
            return data;
        }
        //! -----------
        //! "Analysis"
        //! -----------
        if(name =="Analysis")
        {
            void *p = QStandardItem::data(Qt::UserRole).value<Property>().getData().value<void*>();
            QStandardItem *item = (QStandardItem*)(p);
            if(item!=Q_NULLPTR) return item->data(Qt::DisplayRole);
            data.setValue(QString("No analysis root selected"));
            return data;
        }
        //! --------------
        //! "Assignement"
        //! --------------
        if(name == "Assignment")
        {
            int val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toInt();
            QString matName;
            switch(val)
            {
            case 0: matName = "Structural steel"; break;
            case 1: matName = "Bilinear steel"; break;
            case 2: matName = "H11 fatigue"; break;
            case 3: matName = "F22 fatigue"; break;
            case 4: matName = "B16_fatigue"; break;
            case 5: matName = "F6NM_fatigue"; break;
            case 6: matName = "F92_fatigue"; break;
            case 7: matName = "A479_fatigue"; break;
            case 8: matName = "SA479_XM19_fatigue"; break;
            case 9: matName = "SA182-B8M_CL2"; break;
            case 10: matName = "SA182-F316"; break;
            case 11: matName = "SA352-LCB"; break;
            }
            data.setValue(matName);
            return data;
        }
        //! -----------
        //! "Grouping"
        //! -----------
        if(name == "Grouping")
        {
            int val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toInt();
            switch(val)
            {
            case 0: data.setValue(QString("By bodies")); break;
            case 1: data.setValue(QString("By master faces")); break;
            case 2: data.setValue(QString("By slave faces")); break;
            case 3: data.setValue(QString("By master body")); break;
            case 4: data.setValue(QString("By slave body")); break;
            case 5: data.setValue(QString("Ungrouped")); break;
            }
            return data;
        }
        //! --------------------
        //! "Selected elements"
        //! --------------------
        if(name == "Selected elements")
        {
            const occHandle(Ng_MeshVS_DataSource3D) &meshDS = QStandardItem::data(Qt::UserRole).value<Property>().getData().value<occHandle(Ng_MeshVS_DataSource3D)>();
            int NbNodes = meshDS->GetAllNodes().Extent();
            int NbElements = meshDS->GetAllElements().Extent();
            QString val = QString("Nodes: %1 Elements: %2").arg(NbNodes).arg(NbElements);
            data.setValue(val);
            return data;
        }
        //! ---------------
        //! "Element list"
        //! ---------------
        if(name =="Element list")
        {
            return QStandardItem::data(Qt::UserRole).value<Property>().getData();
        }
        //! -------------------
        //! "Selection method"
        //! -------------------
        if(name =="Selection method")
        {
            int val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toInt();
            switch(val)
            {
            case 0: data.setValue(QString("Elements numbers list")); break;
            case 1: data.setValue(QString("Picking")); break;
            }
            return data;
        }
        //! ------------
        //! "Mesh type"
        //! ------------
        if(name =="Mesh type")
        {
            int val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toInt();
            switch(val)
            {
            case 0: data.setValue(QString("Full tri-tet")); break;
            case 1: data.setValue(QString("Full quad-hexa")); break;
            case 2: data.setValue(QString("Quad-hexa dominant")); break;
            }
            return data;
        }
        //! -----------------------------------
        //! "Activation status" (model change)
        //! -----------------------------------
        if(name =="Activation status")
        {
            Property::modelChangeActivationStatus val = QStandardItem::data(Qt::UserRole).value<Property>().getData().value<Property::modelChangeActivationStatus>();
            switch(val)
            {
            case Property::modelChangeActivationStatus_Remove: data.setValue(QString("Remove")); break;
            case Property::modelChangeActivationStatus_Inactive: data.setValue(QString("Inactive")); break;
            case Property::modelChangeActivationStatus_Add: data.setValue(QString("Add")); break;
            }
            return data;
        }
        //! ---------------------------
        //! "Item type" (model change)
        //! ---------------------------
        if(name =="Item type")
        {
            int val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toInt();
            if(val==0) data.setValue(QString("Bodies"));
            else data.setValue(QString("Contact pairs"));
            return data;
        }
        //! -----------------------------------------------------------
        //! mesh data sources resulting from a geometry user selection
        //! write in the cell the number of mesh data sources
        //! -----------------------------------------------------------
        if(name =="Mesh data sources" || name == "Master mesh data source" || name == "Slave mesh data source")
        {
            const IndexedMapOfMeshDataSources &ds =
                    QStandardItem::data(Qt::UserRole).value<Property>().getData().value<IndexedMapOfMeshDataSources>();            
            data.setValue(QString("Computed on %1 bodies").arg(ds.size()));
            return data;
        }
        //! -------------------------------------------------------------
        //! TetWild parameters - "Envelope sizing" "Ideal length sizing"
        //! -------------------------------------------------------------
        if(name =="Envelope sizing" || name =="Ideal length sizing")
        {
            int val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toInt();
            if(val==0) data.setValue(QString("Relative"));
            else data.setValue(QString("Absolute"));
            return data;
        }
        //! -----------------------------------------------------------------------------------------
        //! TetWild parameters - "Relative size" "Absolute size" "Relative length" "Absolute length"
        //! -----------------------------------------------------------------------------------------
        if(name =="Relative size" || name =="Absolute size" || name =="Relative length" || name =="Absolute length")
        {
            double val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toDouble();
            data.setValue(QString("%1").arg(val));
            return data;
        }

#ifdef COSTAMP_VERSION
        if(name =="Activate")
        {
            bool hasRun = QStandardItem::data(Qt::UserRole).value<Property>().getData().toBool();
            switch(hasRun)
            {
            case false: data.setValue(QString("Press to run")); break;
            case true: data.setValue(QString("Running")); break;
            }
            return data;
        }

        if(name == "Time history file")
        {
            QString timeHistoryFileLoc = QStandardItem::data(Qt::UserRole).value<Property>().getData().toString();
            data.setValue(timeHistoryFileLoc);
            if(timeHistoryFileLoc.isEmpty()) data.setValue(QString("Select file"));
            return data;
        }

        if(name =="Intensification pressure" || name=="Force value")
        {
            double val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toDouble();
            data.setValue(QString("%1").arg(val));
            return data;
        }
#endif

        if(name =="Source file path")
        {
            QString val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toString();
            if(val.isEmpty() || val.isNull()) data.setValue(QString("No file loaded"));
            else data.setValue(val);
            return data;
        }
        if(name =="Tolerance value")
        {
            double val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toDouble();
            data.setValue(val);
            return data;
        }
        if(name =="Geometry healing")
        {
            bool val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toBool();
            if(val==false) data.setValue(QString("No"));
            else data.setValue(QString("Yes"));
            return data;
        }
        if(name =="Solid bodies" || name =="Surface bodies" || name =="Line bodies")
        {
            bool val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toBool();
            if(val==false) data.setValue(QString("No"));
            else data.setValue(QString("Yes"));
            return data;
        }
        if(name =="Tessellator")
        {
            Property::meshEngine2D sd = QStandardItem::data(Qt::UserRole).value<Property>().getData().value<Property::meshEngine2D>();
            switch(sd)
            {
            case Property::meshEngine2D_OCC_STL: data.setValue(QString("Standard STL")); break;
            case Property::meshEngine2D_OCC_ExpressMesh: data.setValue(QString("Express mesh")); break;
            }
            return data;
        }
        if(name =="Patch conforming" || name =="Defeaturing" || name =="Healing" || name =="Simplification"
                || name == "Preserve boundary conditions edges" || name == "Project points on geometry")
        {
            bool val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toBool();
            if(val==false) data.setValue(QString("Off"));
            else data.setValue(QString("On"));
            return data;
        }
        if(name =="Surface mesher")
        {
            Property::meshEngine2D val = QStandardItem::data(Qt::UserRole).value<Property>().getData().value<Property::meshEngine2D>();
            switch(val)
            {
            case Property::meshEngine2D_Netgen: data.setValue(QString("Netgen")); break;
            case Property::meshEngine2D_Netgen_STL: data.setValue(QString("Netgen STL")); break;
            case Property::meshEngine2D_OCC_ExpressMesh : data.setValue(QString("Express mesh")); break;
            }
            return data;
        }
        if(name =="Volume mesher")
        {
            Property::meshEngine3D val = QStandardItem::data(Qt::UserRole).value<Property>().getData().value<Property::meshEngine3D>();
            switch(val)
            {
            case Property::meshEngine3D_Netgen_STL: data.setValue(QString("Netgen STL")); break;
            case Property::meshEngine3D_Netgen: data.setValue(QString("Netgen")); break;
            case Property::meshEngine3D_Tetgen: data.setValue(QString("Tetgen")); break;
            case Property::meshEngine3D_Tetgen_BR: data.setValue(QString("Tetgen BR")); break;
            case Property::meshEngine3D_TetWild: data.setValue(QString("Experimental mesher")); break;
            }
            return data;
        }
        if(name =="Angular deflection" || name =="Linear deflection" || name =="Min face size" || name =="Max face size")
        {
            double val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toDouble();
            data.setValue(QString("%1").arg(val));
            return data;
        }
        if(name =="Fatigue algo")
        {
            int val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toInt();
            switch(val)
            {
            case 0: data.setValue(QString("Basquin Coffin Manson")); break;
            case 1: data.setValue(QString("Equivalent strain range")); break;
            }
            return data;
        }
        if(name=="Mapping")
        {
            int val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toInt();
            switch(val)
            {
            case 0: data.setValue(QString("Linear")); break;
            case 1: data.setValue(QString("Logarithmic")); break;
            }
            return data;
        }
        if(name=="Stress/strain source")
        {
            int val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toInt();
            switch(val)
            {
            case 0: data.setValue(QString("Equivalent strain")); break;
            case 1: data.setValue(QString("Equivalent plastic strain")); break;
            }
            return data;
        }
        if(name=="Component")
        {
            int val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toInt();
            switch(val)
            {
            case 0:
            {
                int fatigueAlgo = this->getCurrentNode()->getPropertyValue<int>("Fatigue algo");
                data.setValue(QString("Total"));
                if(fatigueAlgo==1) data.setValue(QString("Not available"));
            }
                break;
            case 1: data.setValue(QString("Mechanical")); break;
            case 2: data.setValue(QString("Thermal")); break;
            }
            return data;
        }
        if(name=="Triaxiality correction")
        {
            int val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toInt();
            switch(val)
            {
            case 0: data.setValue(QString("Off")); break;
            case 1: data.setValue(QString("On")); break;
            }
            return data;
        }
        if(name=="Number of cycles")
        {
            int val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toInt();
            data.setValue(QString("%1").arg(val));
            return data;
        }
        if(name =="Rupture stress" || name =="Rupture strain" || name =="n'" || name=="b")
        {
            double val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toDouble();
            data.setValue(QString("%1").arg(val));
            return data;
        }
        if(name =="Import solid bodies" || name =="Import surface bodies" || name =="Import line bodies")
        {
            bool val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toBool();
            switch(val)
            {
            case true: data.setValue(QString("Yes")); break;
            case false: data.setValue(QString("No")); break;
            }
            return data;
        }
        if(name =="Shape" || name=="Source geometry")
        {
            //const TopoDS_Shape_Reg &shape = QStandardItem::data(Qt::UserRole).value<Property>().getData().value<TopoDS_Shape_Reg>();
            const TopoDS_Shape &shape = QStandardItem::data(Qt::UserRole).value<Property>().getData().value<TopoDS_Shape>();
            if(shape.IsNull())
            {
                data.setValue(QString("no shape loaded"));
                return data;
            }
            switch(shape.ShapeType())
            {
            case TopAbs_COMPOUND: data.setValue(QString("Compound")); break;
            case TopAbs_COMPSOLID: data.setValue(QString("Composite solid")); break;
            case TopAbs_SOLID: data.setValue(QString("Solid")); break;
            case TopAbs_SHELL: data.setValue(QString("Shell")); break;
            case TopAbs_FACE: data.setValue(QString("Face")); break;
            case TopAbs_WIRE: data.setValue(QString("Wire")); break;
            case TopAbs_EDGE: data.setValue(QString("Edge")); break;
            }
            return data;
        }
        if(name == "3D mesh" || name == "2D mesh")
        {
            const occHandle(MeshVS_Mesh) &val = QStandardItem::data(Qt::UserRole).value<Property>().getData().value<MeshVS_Mesh_handle_reg>();
            if(val!=occHandle(MeshVS_Mesh)())
                data.setValue(QString("Nodes: %1 Elements: %2").
                              arg(val->GetDataSource()->GetAllNodes().Extent()).
                              arg(val->GetDataSource()->GetAllElements().Extent()));
            else data.setValue(QString("not meshed"));
            return data;
        }
        if(name == "Graphic object")
        {
            const occHandle(AIS_Shape) &val = QStandardItem::data(Qt::UserRole).value<Property>().getData().value<AIS_DoubleArrowMarker_handle_reg>();
            data.setValue(QString::fromLatin1(val->get_type_name()));
            return data;
        }
        if(name == "Volume mesh engine")
        {
            Property::meshEngine3D volumeMeshEngine= QStandardItem::data(Qt::UserRole).value<Property>().getData().value<Property::meshEngine3D>();
            switch(volumeMeshEngine)
            {
            case Property::meshEngine3D_Netgen: data.setValue(QString("Netgen")); break;
            case Property::meshEngine3D_Tetgen: data.setValue(QString("Tetgen")); break;
            }
            return data;
        }
        if(name =="Transition")
        {
            double val = QStandardItem::data(Qt::UserRole).value<Property>().getData().value<double>();
            data.setValue(QString("%1").arg(val));
            return data;
        }
        if(name=="Guiding vectors smoothing steps" || name =="Thickness smoothing steps")
        {
            int val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toInt();
            data.setValue(QString("%1").arg(val));
            return data;
        }
        if(name=="Curvature sensitivity" || name=="Guiding vector smoothing - curvature sensitivity" || name=="Thickness smoothing - curvature sensitivity")
        {
            double val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toDouble();
            data.setValue(QString("%1").arg(val));
            return data;
        }
        if(name =="Lock boundary")
        {
            bool val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toBool();
            if(val==true) data.setValue(QString("Locked"));
            else data.setValue(QString("Free"));
            return data;
        }
        if(name=="Options")
        {
            bool val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toInt();
            switch(val)
            {
            case 0: data.setValue(QString("First layer thickness")); break;
            case 1: data.setValue(QString("Total thickness")); break;
            }
            return data;
        }
        if(name== "Number of layers")
        {
            int val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toInt();
            data.setValue(QString("%1").arg(val));
            return data;
        }
        if(name== "Expansion ratio")
        {
            double val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toDouble();
            data.setValue(QString("%1").arg(val));
            return data;
        }
        if(name =="Check self intersections" || name =="Check mutual intersections")
        {
            bool val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toBool();
            if(val == true) data.setValue(QString("On")); else data.setValue(QString("Off"));
            return data;
        }
        if(name== "Total thickness")
        {
            double val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toDouble();
            data.setValue(QString("%1").arg(val));
            return data;
        }
        if(name== "First layer height")
        {
            double val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toDouble();
            data.setValue(QString("%1").arg(val));
            return data;
        }
        if(name=="Run in memory")
        {
            bool val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toInt();
            if(val==true) data.setValue(QString("Active"));
            else data.setValue(QString("Off"));
            return data;
        }
        //! -----------------
        //! "Surface mesher"
        //! -----------------
        if(name=="Surface mesher")
        {
            int val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toInt();
            switch(val)
            {
            case Property::meshEngine2D_OCC_STL: data.setValue(QString("OCC STL")); break;
            case Property::meshEngine2D_OCC_ExpressMesh: data.setValue(QString("Express mesh")); break;
            case Property::meshEngine2D_Netgen_STL: data.setValue(QString("Netgen STL")); break;
            }
            return data;
        }
        //! ----------------
        //! "Pair distance"
        //! ----------------
        if(name =="Pair distance")
        {
            double val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toDouble();
            data.setValue(QString("%1").arg(val));
            return data;
        }
        //! ---------------------------------
        //! "Level" (of mesh simplification)
        //! ---------------------------------
        if(name =="Level")
        {
            double val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toDouble();
            data.setValue(QString("%1").arg(val));
            return data;
        }
        //! ----------------------------------------
        //! color box property: type of scale
        //! ----------------------------------------
        if(name =="Mesh healing" || name == "Mesh simplification" || name =="Mesh defeaturing")
        {
            int val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toInt();
            if(val==0) data.setValue(QString("Off"));
            else data.setValue(QString("On"));
            return data;
        }
        if(name =="Method")
        {            
            Property::meshMethod mm = QStandardItem::data(Qt::UserRole).value<Property>().getData().value<Property::meshMethod>();
            switch(mm)
            {
            case Property::meshMethod_automatic:
                data.setValue(QString("Automatic"));
                break;
            case Property::meshMethod_NetgenTetgen:
                data.setValue(QString("Netgen Tetgen"));
                break;
            case Property::meshMethod_NetgenNetgen:
                data.setValue(QString("Full Netgen"));
                break;
            case Property::meshMethod_EMeshTetgen:
                data.setValue(QString("Emesh Tetgen"));
                break;
            case Property::meshMethod_EMeshNetgen:
                data.setValue(QString("Emesh Netgen"));
                break;
            }
            return data;
        }
        //! ----------------------------------------------
        //! "Metric type" - currently implemented for tet
        //! 0 => "Modified Mean Ratio (MMR)"
        //! 1 => "Modified Condition Number (MCN)"
        //! 2 => "Modified volume-length (MIVL)"
        //! ----------------------------------------------
        if(name =="Metric type")
        {
            int val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toInt();
            switch(val)
            {
            case 0: data.setValue(QString("Modified Mean Ratio")); break;
            case 1: data.setValue(QString("Modified Condition Number")); break;
            case 2: data.setValue(QString("Modified volume-length")); break;
            case 3: data.setValue(QString("None")); break;
            }
            return data;
        }
        //! --------------
        //! "Metric data"
        //! --------------
        if(name =="Metric data")
        {
            histogramData hisData = QStandardItem::data(Qt::UserRole).value<Property>().getData().value<histogramData>();
            data.setValue(QString("%1 bins").arg(hisData.size()));
            return data;
        }
        //! -------------
        //! "Scale type"
        //! -------------
        if(name =="Scale type")
        {
            int n = QStandardItem::data(Qt::UserRole).value<Property>().getData().toInt();
            switch(n)
            {
            case 0: data.setValue(QString("Autoscale")); break;
            case 1: data.setValue(QString("Custom")); break;
            }
            return data;
        }
        //! ----------------------------------------
        //! color box property: number of intervals
        //! ----------------------------------------
        if(name =="# intervals")
        {
            int n = QStandardItem::data(Qt::UserRole).value<Property>().getData().toInt();
            data.setValue(n);
            if(n==0) data.setValue(QString("Program controlled"));
            return data;
        }
        //! ----------------------------------------
        //! color box property: min and max values
        //! ----------------------------------------
        if(name =="Min" || name =="Max")
        {
            double val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toDouble();
            data.setValue(QString("%1").arg(val));
            return data;
        }
        //! ----------------------------------------------------------
        //! reference point for remote force and remote displacements
        //! ----------------------------------------------------------
        if(name =="Reference point")
        {
            QVector<double> RP = QStandardItem::data(Qt::UserRole).value<Property>().getData().value<QVector<double>>();
            if(!RP.isEmpty()) data.setValue(QString("(%1; %2; %3)").arg(RP.at(0)).arg(RP.at(1)).arg(RP.at(2)));
            else data.setValue(QString("not defined yet"));
            return data;
        }
        if(name =="--Recurrence rate")
        {
            int val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toInt();
            data.setValue(QString("%1").arg(val));
            return data;
        }
        if(name =="Store results at")
        {
            int val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toInt();
            switch(val)
            {
            case 0: data.setValue(QString("All time points")); break;
            case 1: data.setValue(QString("Last time point")); break;
            case 2: data.setValue(QString("Specified recurrence rate")); break;
            }
            return data;
        }
        if(name =="Stress" || name =="Strain" || name =="Reaction forces" || name=="Contact data")
        {
            int val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toInt();
            switch(val)
            {
            case 0: data.setValue(QString("No")); break;
            case 1: data.setValue(QString("Yes")); break;
            }
            return data;
        }
        if(name =="Small sliding")
        {
            int val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toInt();
            switch(val)
            {
            case 0: data.setValue(QString("Off")); break;
            case 1: data.setValue(QString("On")); break;
            }
            return data;
        }
        if(name == "K" || name =="Sigma infinity" || name =="C0" || name =="Lambda"
                || name =="P0" || name =="Beta" || name == "Thermal conductance")
        {
            double val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toDouble();
            if(val==0.0) data.setValue(QString("Program controlled"));
            else data.setValue(QString("%1").arg(val));
            return data;
        }
        if(name == "Line search")
        {
            int val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toInt();
            switch(val)
            {
            case 0: data.setValue(QString("Program controlled")); break;
            case 1: data.setValue(QString("Custom")); break;
            }
            return data;
        }
        if(name == "Max value" || name == "Min value")
        {
            double val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toDouble();
            data.setValue(QString("%1").arg(val));
            return data;
        }
        if(name == "Time incrementation" || name == "Cutback factors")
        {
            int val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toInt();
            switch(val)
            {
            case 0: data.setValue(QString("Program controlled")); break;
            case 1: data.setValue(QString("Custom")); break;
            }
            return data;
        }
        if(name == "D_f" || name == "D_C" || name == "D_B"|| name == "D_A" || name == "D_S" || name == "D_H"
                || name == "D_D" || name == "W_G")
        {
            double val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toDouble();
            data.setValue(QString("%1").arg(val));
            if(name =="D_S" || name =="D_H" || name =="W_G") data.setValue(QString("Unused"));
            return data;
        }
        if(name =="I_0" || name =="I_R" || name =="I_P" || name =="I_C" || name =="I_L"
                || name =="I_G" || name =="I_S" || name =="I_A" || name =="I_J" || name =="I_T")
        {
            int val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toInt();
            data.setValue(QString("%1").arg(val));
            if(name =="I_S" || name =="I_J" || name =="I_T") data.setValue(QString("Unused"));
            return data;
        }
        if(name =="Flux convergence" || name=="Solution convergence")
        {
            int val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toInt();
            switch(val)
            {
            case 0: data.setValue(QString("Remove")); break;
            case 1: data.setValue(QString("On")); break;
            case 2: data.setValue(QString("Program controlled")); break;
            }
            return data;
        }
        if(name =="--Value")
        {
            double val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toDouble();
            if(val==0.0) data.setValue(QString("Calculated by solver"));
            else data.setValue(QString("%1").arg(val));
            return data;
        }
        if(name =="--q_alpha_0")
        {
            double val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toDouble();
            if(val==0.0) data.setValue(QString("Calculated by solver"));
            else data.setValue(QString("%1").arg(val));
            return data;
        }
        if(name == "--R_alpha_n" || name == "--R_alpha_P" || name == "--R_alpha_l" || name == "--epsilon_alpha" ||
                name == "--C_alpha_n" || name == "--C_alpha_epsilon")
        {
            double val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toDouble();
            data.setValue(QString("%1").arg(val));
            return data;
        }
        if(name =="Pinball")
        {
            double radius = QStandardItem::data(Qt::UserRole).value<Property>().getData().toDouble();
            data.setValue(QString("%1").arg(radius));
            return data;
        }
        if(name =="Surface mesh")
        {
            const occHandle(Ng_MeshVS_DataSource2D) &aMesh = QStandardItem::data(Qt::UserRole).value<Property>().getData().value<occHandle(Ng_MeshVS_DataSource2D)>();
            if(!aMesh.IsNull())
            {
                int NbNodes = aMesh->GetAllNodes().Extent();
                int NbElements = aMesh->GetAllElements().Extent();
                data.setValue(QString("Nodes: %1 Elements: %2").arg(NbNodes).arg(NbElements));
            }
            else data.setValue(QString("Not meshed"));
            return data;
        }
        else if(name =="Volume mesh")
        {
            const occHandle(Ng_MeshVS_DataSource3D) &aMesh = QStandardItem::data(Qt::UserRole).value<Property>().getData().value<occHandle(Ng_MeshVS_DataSource3D)>();
            if(!aMesh.IsNull())
            {
                int NbNodes = aMesh->GetAllNodes().Extent();
                int NbElements = aMesh->GetAllElements().Extent();
                data.setValue(QString("Nodes: %1 Elements: %2").arg(NbNodes).arg(NbElements));
            }
            else data.setValue(QString("Not meshed"));
            return data;
        }
        else if(name =="Project files dir")
        {
            QString val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toString();
            data.setValue(val);
            return data;
        }
        else if (name == "Simulation time")
        {
            double val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toDouble();
            data.setValue(val);
            return data;
        }
        else if (name == "Step" || name == "Substep")
        {
            int val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toInt();
            data.setValue(val);
            return data;
        }
        else if (name == "Simulation status")
        {
            int val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toInt();
            switch(val)
            {
            case 0: data.setValue(QString("No convergence")); break;
            case 1: data.setValue(QString("Convergence")); break;
            }
            return data;
        }
        //! post processing
        else if(name =="Large deflection")
        {
            int val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toInt();
            if(val == 0) data.setValue(QString("Off"));
            else data.setValue(QString("On"));
            return data;
        }
        else if(name=="Display time")
        {
            double val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toDouble();
            char v[16];
            sprintf(v,"%g",val);
            data.setValue(QString::fromLatin1(v));
            return data;
        }
        else if(name=="Mode number")
        {
            int val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toInt();
            data.setValue(val);
            return data;
        }
        else if(name=="Set number")
        {
            int val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toInt();
            char v[16];
            sprintf(v,"%d",val);
            data.setValue(QString::fromLatin1(v));
            return data;
        }
        else if(name =="By")
        {
            SimulationNodeClass::nodeType nodeType = this->getCurrentNode()->getType();
            switch(nodeType)
            {
            case SimulationNodeClass::nodeType_solutionThermalFlux:
            case SimulationNodeClass::nodeType_solutionThermalTemperature:
            case SimulationNodeClass::nodeType_solutionStructuralContact:
            case SimulationNodeClass::nodeType_solutionStructuralNodalForces:
            case SimulationNodeClass::nodeType_solutionStructuralEquivalentPlasticStrain:
            case SimulationNodeClass::nodeType_solutionStructuralMechanicalStrain:
            case SimulationNodeClass::nodeType_solutionStructuralNodalDisplacement:
            case SimulationNodeClass::nodeType_solutionStructuralStress:
            case SimulationNodeClass::nodeType_solutionStructuralTemperature:
            case SimulationNodeClass::nodeType_solutionStructuralThermalStrain:
            case SimulationNodeClass::nodeType_solutionStructuralTotalStrain:
            case SimulationNodeClass::nodeType_solutionStructuralGamma:
            case SimulationNodeClass::nodeType_solutionStructuralReactionForce:
            {
                int val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toInt();
                switch(val)
                {
                case 0: data.setValue(QString("Time")); break;
                case 1: data.setValue(QString("Set")); break;
                }
            }
                break;

            case SimulationNodeClass::nodeType_meshMethod:
            {
                int val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toInt();
                switch(val)
                {
                case 0: data.setValue(QString("Level")); break;
                case 1: data.setValue(QString("Pair distance")); break;
                }
            }
                break;
            }
            return data;
        }
        else if(name =="Type ")
        {
            int val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toInt();
            SimulationNodeClass *curNode = this->getCurrentNode();
            SimulationNodeClass::nodeType type = curNode->getType();
            switch (type)
            {
            case SimulationNodeClass::nodeType_solutionThermalTemperature:
                data.setValue(QString("Temperature"));
                break;

            case SimulationNodeClass::nodeType_solutionThermalFlux:
                data.setValue(QString("Thermal flux"));
                break;

            case SimulationNodeClass::nodeType_solutionStructuralEquivalentPlasticStrain:
                data.setValue(QString("Equivalent plastic strain"));
                break;

            case SimulationNodeClass::nodeType_solutionStructuralTemperature:
                data.setValue(QString("Temperature"));
                break;

            case SimulationNodeClass::nodeType_solutionStructuralNodalDisplacement:
                switch (val)
                {
                case 0: data.setValue(QString("Total displacement")); break;
                case 1: data.setValue(QString("Directional displacement X")); break;
                case 2: data.setValue(QString("Directional displacement Y")); break;
                case 3: data.setValue(QString("Directional displacement Z")); break;
                }
                break;

            case SimulationNodeClass::nodeType_solutionStructuralStress:
                switch(val)
                {
                case 0: data.setValue(QString("Equivalent stress")); break;
                case 1: data.setValue(QString("Stress intensity")); break;
                case 2: data.setValue(QString("Maximum principal")); break;
                case 3: data.setValue(QString("Middle principal stress")); break;
                case 4: data.setValue(QString("Minimum principal stress")); break;
                case 5: data.setValue(QString("Normal stress X")); break;
                case 6: data.setValue(QString("Normal stress Y")); break;
                case 7: data.setValue(QString("Normal stress Z")); break;
                case 8: data.setValue(QString("Shear stress XY")); break;
                case 9: data.setValue(QString("Shear stress YZ")); break;
                case 10: data.setValue(QString("Shear stress ZX")); break;
                }

                break;

            case SimulationNodeClass::nodeType_solutionStructuralTotalStrain:
                switch(val)
                {
                case 0: data.setValue(QString("Equivalent strain")); break;
                case 1: data.setValue(QString("Strain intensity")); break;
                case 2: data.setValue(QString("Maximum principal")); break;
                case 3: data.setValue(QString("Middle principal")); break;
                case 4: data.setValue(QString("Minimum principal")); break;
                }
                break;

            case SimulationNodeClass::nodeType_solutionStructuralNodalForces:
                switch(val)
                {
                case 0: data.setValue(QString("Total")); break;
                case 1: data.setValue(QString("X direction")); break;
                case 2: data.setValue(QString("Y direction")); break;
                case 3: data.setValue(QString("Z direction")); break;
                }
                break;

            case SimulationNodeClass::nodeType_solutionStructuralContact:
                switch(val)
                {
                case 0: data.setValue(QString("Contact pressure")); break;
                case 1: data.setValue(QString("Frictional stress")); break;
                case 2: data.setValue(QString("Penetration")); break;
                case 3: data.setValue(QString("Sliding")); break;
                }
                break;
            }
            return data;
        }

        else if(name =="Solution information")
        {
            int val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toInt();
            switch(val)
            {
            case 0: data.setValue(QString("Solver output")); break;
            case 1: data.setValue(QString("Force convergence")); break;
            case 2: data.setValue(QString("Displacement convergence")); break;
            case 3: data.setValue(QString("Line search")); break;
            case 4: data.setValue(QString("Time step size")); break;
            }
            return data;
        }
        else if(name =="X rotation" || name =="Y rotation" || name =="Z rotation" ||
                name =="X component " || name =="Y component " || name =="Z component ")
        {
            //! --------------------------------------------------------------------------------
            //! the keys "X component " "Y component " "Z component " have one space at the end
            //! in order to distinguish them from "X component" "Y component" "Z component"
            //! --------------------------------------------------------------------------------
            int val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toInt();
            switch(val)
            {
            case 0: data.setValue(QString("Inactive")); break;
            case 1: data.setValue(QString("Active")); break;
            }
            return data;
        }
        else if(name =="DOFs selection")
        {
            int val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toInt();
            switch(val)
            {
            case 0: data.setValue(QString("Program controlled")); break;
            case 1: data.setValue(QString("Manual")); break;
            }
            return data;
        }
        else if(name =="Coupling")
        {
            int val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toInt();
            switch(val)
            {
            case 0: data.setValue(QString("Kinematic")); break;
            case 1: data.setValue(QString("Distributed")); break;
            }
            return data;
        }
        else if(name =="First color" || name =="Second color")
        {
            QVector<int> rgb = QStandardItem::data(Qt::UserRole).value<Property>().getData().value<QVector<int>>();
            QColor color(rgb.at(0),rgb.at(1),rgb.at(2),255);
            QString s = color.name();
            data.setValue(s);
            return data;
        }
        else if(name =="Gradient")
        {
            int val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toInt();
            switch(val)
            {
            case 0: data.setValue(QString("Uniform")); break;
            case 1: data.setValue(QString("Horizontal")); break;
            case 2: data.setValue(QString("Vertical")); break;
            }
            return data;
        }
        else if(name =="Number of threads")
        {
            int val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toInt();
            data.setValue(QString("%1").arg(val));
            return data;
        }
        if (name =="Analysis time" || name == "Source time")
        {
            double value = QStandardItem::data(Qt::UserRole).value<Property>().getData().toDouble();
            data.setValue(QString("%1").arg(value));
            return data;
        }
        else if (name =="X buckets" ||
                 name =="Y buckets" ||
                 name =="Z buckets" ||
                 name == "Remapping steps")
        {
            int value = QStandardItem::data(Qt::UserRole).value<Property>().getData().toInt();
            data.setValue(QString("%1").arg(value));
            return data;
        }
        //! ------------------------------------
        //! "Algorithm" interpolation algorithm
        //! ------------------------------------
        else if(name == "Algorithm")
        {
            switch(this->getCurrentNode()->getType())
            {
                case SimulationNodeClass::nodeType_importedBodyScalar:
                {
                    int value = QStandardItem::data(Qt::UserRole).value<Property>().getData().toInt();
                    switch(value)
                    {
                    case 0: data.setValue(QString("Nearest point 1")); break;
                    case 1: data.setValue(QString("Nearest point 2")); break;
                    case 2: data.setValue(QString("Shape functions")); break;
                    }
                }
                break;

            case SimulationNodeClass::nodeType_meshPrismaticLayer:
            {
                int value = QStandardItem::data(Qt::UserRole).value<Property>().getData().toInt();
                switch(value)
                {
                case 0: data.setValue(QString("Pre")); break;
                case 1: data.setValue(QString("Post")); break;
                }
            }
                break;
            }
            return data;
        }
        //! ---------------------
        //! "Boundary mesh type"
        //! ---------------------
        else if(name =="Boundary mesh type")
        {
            int value = QStandardItem::data(Qt::UserRole).value<Property>().getData().toInt();
            switch(value)
            {
            case 0: data.setValue(QString("Hybrid")); break;
            case 1: data.setValue(QString("Tetrahedral")); break;
            }
            return data;
        }
        //! -------------------
        //! "Target directory"
        //! -------------------
        else if(name =="Target directory")
        {
            if(QStandardItem::data(Qt::UserRole).value<Property>().getData().toString().isEmpty())
                data.setValue(QString("Select folder"));
            else data.setValue(QStandardItem::data(Qt::UserRole).value<Property>().getData().toString());
            return data;
        }
        else if(name =="Split data")
        {
            int value = QStandardItem::data(Qt::UserRole).value<Property>().getData().toInt();
            switch(value)
            {
            case 0: data.setValue(QString("Single file")); break;
            case 1: data.setValue(QString("Multiple files")); break;
            }
            return data;
        }
        else if(name =="Relevance")
        {
            int value = QStandardItem::data(Qt::UserRole).value<Property>().getData().toDouble();
            data.setValue(QString("%1").arg(value));
            return data;
        }
        else if(name =="Minimum edge length")
        {
            double value = QStandardItem::data(Qt::UserRole).value<Property>().getData().toDouble();
            char v[32];
            sprintf(v,"%.6g",value);
            data.setValue(QString::fromLatin1(v));
            return data;
        }
        if(name =="Initial size seed")
        {
            int value = QStandardItem::data(Qt::UserRole).value<Property>().getData().toInt();
            switch(value)
            {
            case 0: data.setValue(QString("Program controlled")); break;
            case 1: data.setValue(QString("Assembly")); break;
            case 2: data.setValue(QString("Part")); break;
            }
            return data;
        }
        else if(name =="Element midside nodes")
        {
            int value = QStandardItem::data(Qt::UserRole).value<Property>().getData().toInt();
            switch(value)
            {
            case 0: data.setValue(QString("Program controlled")); break;
            case 1: data.setValue(QString("Dropped")); break;
            case 2: data.setValue(QString("Kept")); break;
            }
            return data;
        }
        else if(name == "Post object")
        {
            //const postObject &aPostObject = QStandardItem::data(Qt::UserRole).value<Property>().getData().value<postObject>();
            const sharedPostObject &aPostObject = QStandardItem::data(Qt::UserRole).value<Property>().getData().value<sharedPostObject>();
            data.setValue(QString("%1 mesh objects").arg(aPostObject->NbMeshes()));
            return data;
        }
        else if(name=="Remap")
        {
            bool val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toBool();
            if(val==true) data.setValue(QString("Yes"));
            else if(val==false) data.setValue(QString("No"));
            return data;
        }
        //! ------------
        //! "Submeshes"
        //! ------------
        else if(name=="Submeshes")
        {
            bool val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toBool();
            if(val==true) data.setValue(QString("At meshing time"));
            else data.setValue(QString("Deferred"));
            return data;
        }
        else if(name=="Generate")
        {
            data.setValue(QString("Press to start"));
            return data;
        }
        //else if(name =="Translate")
        //{
        //    data.setValue(QString("Press to start"));
        //    return data;
        //}
        else if(name=="Step number")
        {
            int step = QStandardItem::data(Qt::UserRole).value<Property>().getData().toInt();
            data.setValue(QString("%1").arg(step));
            return data;
        }
        else if(name=="Step selection mode")
        {
            switch(QStandardItem::data(Qt::UserRole).value<Property>().getData().toInt())
            {
            case 0: data.setValue(QString("All steps")); break;
            case 1: data.setValue(QString("First")); break;
            case 2: data.setValue(QString("Last")); break;
            case 3: data.setValue(QString("By number")); break;
            case 4: data.setValue(QString("Automatic time history")); break;
            }
            return data;
        }
        else if(name=="Source file")
        {
            if(QStandardItem::data(Qt::UserRole).value<Property>().getData().toString().isEmpty())
                data.setValue(QString("Select file"));
            else data.setValue(QStandardItem::data(Qt::UserRole).value<Property>().getData().toString());
            return data;
        }
        else if(name=="Source directory")
        {
            if(QStandardItem::data(Qt::UserRole).value<Property>().getData().toString().isEmpty())
                data.setValue(QString("Select folder"));
            else data.setValue(QStandardItem::data(Qt::UserRole).value<Property>().getData().toString());
            return data;
        }
        else if(name=="NS index")
        {
            int val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toInt();
            data.setValue(QString("%1").arg(val));
            return data;
        }
        else if(name =="Smoothing")
        {
            int val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toInt();
            switch(val)
            {
            case 0: data.setValue(QString("Off")); break;
            case 1: data.setValue(QString("Low")); break;
            case 2: data.setValue(QString("Medium")); break;
            case 3: data.setValue(QString("High")); break;
            }
            return data;
        }
        else if(name == "Time tag" || name == "Parent time tag")
        {
            data.setValue(QStandardItem::data(Qt::UserRole).value<Property>().getData().toString());
            return data;
        }
        else if(name =="Tolerance")
        {
            if(this->getCurrentNode()->getType()==SimulationNodeClass::nodeType_import)
            {
                int toleranceType = QStandardItem::data(Qt::UserRole).value<Property>().getData().toInt();
                switch(toleranceType)
                {
                case 0: data.setValue(QString("User defined")); break;
                case 1: data.setValue(QString("Normal")); break;
                case 2: data.setValue(QString("Loose")); break;
                }
            }
            else
            {
                double val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toDouble();
                data.setValue(QString("%1").arg(val));
            }
            return data;
        }
        else if(name =="Angular criterion")
        {
            double val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toDouble();
            data.setValue(QString("%1").arg(val));
            return data;
        }
        else if(name == "Tags" || name == "Tags master" || name =="Tags slave" || name =="Boundary tags")
        {
            if(!QStandardItem::data(Qt::UserRole).value<Property>().getData().isValid())
            {
                data.setValue(QString("Multiple selection"));
                return data;
            }
            std::vector<GeometryTag> vecLocs = QStandardItem::data(Qt::UserRole).value<Property>().getData().value<std::vector<GeometryTag>>();
            QString val;
            std::vector<GeometryTag>::iterator it;
            if(vecLocs.size()!=0)
                for(it=vecLocs.begin();it!=vecLocs.end();++it)
                {
                    const GeometryTag &loc = *it;
                    int parentShapeNumber = loc.parentShapeNr;
                    int childShapeNumber = loc.subTopNr;
                    bool isParent = loc.isParent;
                    TopAbs_ShapeEnum type = loc.subShapeType;
                    QString shapeType;
                    switch(type)
                    {
                    case TopAbs_VERTEX: shapeType = "Points"; break;
                    case TopAbs_EDGE: shapeType = "Edged"; break;
                    case TopAbs_FACE: shapeType = "Faces"; break;
                    case TopAbs_WIRE: shapeType = "Wires"; break;
                    case TopAbs_SOLID: shapeType = "Solids"; break;
                    case TopAbs_SHELL: shapeType = "Shells"; break;
                    case TopAbs_COMPSOLID: shapeType = "CSolids"; break;
                    case TopAbs_COMPOUND: shapeType = "Compounds"; break;
                    }
                    val.append(QString("(%1, %2 - %3)%4 ").arg(parentShapeNumber).arg(childShapeNumber).arg(isParent).arg(shapeType));
                }
            else val ="Empty tag list";
            data.setValue(val);
            return val;
        }
        else if(name == "Map index")
        {
            QString val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toString();
            data.setValue(val);
            return data;
        }
        else if(name=="Show mesh nodes")
        {
            bool showNodes = QStandardItem::data(Qt::UserRole).value<Property>().getData().toBool();
            QString val;
            showNodes==true? val="Active":val="Off";
            data.setValue(val);
            return data;
        }
        else if(name =="X coordinate" || name =="Y coordinate" || name =="Z coordinate" ||
                name =="X abs coordinate" || name =="Y abs coordinate" || name =="Z abs coordinate")
        {
            double val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toDouble();
            data.setValue(QString("%1").arg(val));
            return data;
        }
        else if(name =="Environment temperature")
        {
            QString name = QStandardItem::data(Qt::UserRole).value<Property>().getData().toString();
            data.setValue(name.append(" K"));
            return data;
        }
        else if(name =="Name" || name =="Assignment")
        {
            QString name = QStandardItem::data(Qt::UserRole).value<Property>().getData().toString();
            data.setValue(name);
            return data;
        }
        else if(name =="Is valid")
        {
            bool isValid = QStandardItem::data(Qt::UserRole).value<Property>().getData().toBool();
            if(isValid) data.setValue(QString("Valid"));
            else data.setValue(QString("Invalid"));
            return data;
        }
        else if(name =="Initial substeps" || name =="Minimum substeps" || name =="Maximum substeps"
                || name =="Number of substeps")
        {
            int N = QStandardItem::data(Qt::UserRole).value<Property>().getData().toInt();
            data.setValue(QString("%1").arg(N));
            return data;
        }
        else if(name == "Status")
        {
            Property::solutionInformation theSolutionInformation = QStandardItem::data(Qt::UserRole).value<Property>().getData().value<Property::solutionInformation>();
            switch(theSolutionInformation)
            {
            case Property::solutionInformation_solveRequired: data.setValue(QString("Solve required")); break;
            case Property::solutionInformation_running: data.setValue(QString("Running")); break;
            case Property::solutionInformation_failed: data.setValue(QString("Solve failed")); break;
            case Property::solutionInformation_interrupted: data.setValue(QString("Solve interrupted")); break;
            }
            return data;
        }
        else if(name == "Update interval")
        {
            double updateInterval = QStandardItem::data(Qt::UserRole).value<Property>().getData().toDouble();
            data.setValue(QString("%1 s").arg(updateInterval));
            return data;
        }
        else if(name=="Base directional data")
        {
            QVector<QVector<double>> dvec = QStandardItem::data(Qt::UserRole).value<Property>().getData().value<QVector<QVector<double>>>();
            data.setValue(QString("[%1; %2; %3][%4; %5; %6][%7; %8; %9]").
                          arg(dvec.at(0).at(0)).arg(dvec.at(0).at(1)).arg(dvec.at(0).at(2)).
                          arg(dvec.at(1).at(0)).arg(dvec.at(1).at(1)).arg(dvec.at(1).at(2)).
                          arg(dvec.at(2).at(0)).arg(dvec.at(2).at(1)).arg(dvec.at(2).at(2)));
            return data;
        }
        else if(name=="Base origin")
        {
            QVector<double> vec = QStandardItem::data(Qt::UserRole).value<Property>().getData().value<QVector<double>>();
            data.setValue(QString("[%1; %2; %3]").arg(vec.at(0)).arg(vec.at(1)).arg(vec.at(2)));
            return data;
        }
        else if (name =="X axis data" || name =="Y axis data" || name =="Z axis data")
        {
            QVector<double> vec = QStandardItem::data(Qt::UserRole).value<Property>().getData().value<QVector<double>>();
            data.setValue(QString("[%1; %2; %3]").arg(vec.at(0)).arg(vec.at(1)).arg(vec.at(2)));
            return data;
        }
        else if(name =="Geometry selection")
        {
            data.setValue(QString("Geometry selection"));
            return data;
        }
        else if (name =="Global coordinates")
        {
            data.setValue(QString("Global coordinates"));
            return data;
        }
        else if (name=="Direction")
        {
            QVector<double> vec = QStandardItem::data(Qt::UserRole).value<Property>().getData().value<QVector<double>>();
            data.setValue(QString("[%1; %2; %3]").arg(vec.at(3)).arg(vec.at(4)).arg(vec.at(5)));
            return data;
        }
        //! -------------------
        //! "Film coefficient"
        //! -------------------
        else if(name =="Film coefficient")
        {
            //! --------------------------
            //! retrieve the tabular data
            //! --------------------------
            SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
            SimulationNodeClass *nodeAnalysisSettings = sm->getActiveAnalysisBranch()->child(0,0)->data(Qt::UserRole).value<SimulationNodeClass*>();
            if(nodeAnalysisSettings->isAnalysisSettings()==false)
            {
                cout<<"@-------------------------------------------------@"<<endl;
                cout<<"@- Analysis settings not found - check interface -@"<<endl;
                cout<<"@-------------------------------------------------@"<<endl;
                data.setValue(QString(""));
                return data;
            }
            int row = nodeAnalysisSettings->getPropertyValue<int>("Current step number");
            int SC = mainTreeTools::calculateStartColumn(sm->myTreeView);
            CustomTableModel *tabularDataModel = nodeAnalysisSettings->getTabularDataModel();
            double val = tabularDataModel->dataRC(row,SC,Qt::EditRole).toDouble();
            Property::loadDefinition theLoadDefinition = QStandardItem::data(Qt::UserRole).value<Property>().getData().value<Property::loadDefinition>();
            switch(theLoadDefinition)
            {
            case Property::loadDefinition_constant: data.setValue(QString("%1 (ramped)").arg(val)); break;
            case Property::loadDefinition_tabularData: data.setValue(QString("Tabular data")); break;
            }
            return data;
        }
        //! ------------------------
        //! "Reference temperature"
        //! ------------------------
        else if(name =="Reference temperature")
        {
            //! --------------------------
            //! retrieve the tabular data
            //! --------------------------
            SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
            SimulationNodeClass *nodeAnalysisSettings = sm->getActiveAnalysisBranch()->child(0,0)->data(Qt::UserRole).value<SimulationNodeClass*>();
            if(nodeAnalysisSettings->isAnalysisSettings()==false)
            {
                cout<<"@-------------------------------------------------@"<<endl;
                cout<<"@- Analysis settings not found - check interface -@"<<endl;
                cout<<"@-------------------------------------------------@"<<endl;
                data.setValue(QString(""));
                return data;
            }
            int row = nodeAnalysisSettings->getPropertyValue<int>("Current step number");
            int SC = mainTreeTools::calculateStartColumn(sm->myTreeView);
            CustomTableModel *tabularDataModel = nodeAnalysisSettings->getTabularDataModel();
            double val = tabularDataModel->dataRC(row,SC+1,Qt::EditRole).toDouble();
            Property::loadDefinition theLoadDefinition = QStandardItem::data(Qt::UserRole).value<Property>().getData().value<Property::loadDefinition>();
            switch(theLoadDefinition)
            {
            case Property::loadDefinition_constant: data.setValue(QString("%1 (ramped)").arg(val)); break;
            case Property::loadDefinition_tabularData: data.setValue(QString("Tabular data")); break;
            }
            return data;
        }
        //! ------------------------------------------------------
        //! "Magniture" "X component" "Y component" "Z component"
        //! ------------------------------------------------------
        else if(name =="Magnitude" || name =="X component" || name =="Y component" || name =="Z component")
        {
            //! --------------------------
            //! retrieve the tabular data
            //! --------------------------
            SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
            SimulationNodeClass *nodeAnalysisSettings = sm->getActiveAnalysisBranch()->child(0,0)->data(Qt::UserRole).value<SimulationNodeClass*>();
            if(nodeAnalysisSettings->isAnalysisSettings()==false)
            {
                cout<<"@-------------------------------------------------@"<<endl;
                cout<<"@- Analysis settings not found - check interface -@"<<endl;
                cout<<"@-------------------------------------------------@"<<endl;
                data.setValue(QString(""));
                return data;
            }
            /*
            CustomTableModel *tabularDataModel = nodeAnalysisSettings->getTabularDataModel();

            //! -----------------------------------------------------------------
            //! retrieve the numerical value of "Magnitude" of "X/Y/Z component"
            //! at the current step number
            //! -----------------------------------------------------------------
            int row = nodeAnalysisSettings->getPropertyValue<int>("Current step number");
            int SC = mainTreeTools::calculateStartColumn(sm->myTreeView);
            int col;

            double val;
            if(name =="Magnitude")
            {
                col = SC;
                val = tabularDataModel->dataRC(row,col,Qt::EditRole).toDouble();
            }
            else
            {
                SimulationNodeClass *curNode = this->getCurrentNode();
                bool b0 = curNode->getPropertyValue<Property::loadDefinition>("X component")!=Property::loadDefinition_free? true:false;
                bool b1 = curNode->getPropertyValue<Property::loadDefinition>("Y component")!=Property::loadDefinition_free? true:false;
                bool b2 = curNode->getPropertyValue<Property::loadDefinition>("Z component")!=Property::loadDefinition_free? true:false;

                //! --------------------------------------------
                //! one column for the load component is active
                //! --------------------------------------------
                if((b0==true && b1==false && b2==false) || (b0==false && b1==true && b2==false) || (b0==false && b1==false &&b2 ==true))
                {
                    col = SC;
                    val = tabularDataModel->dataRC(row,col,Qt::EditRole).toDouble();
                }
                //! -----------------------------------------------
                //! two columns for the load components are active
                //! (x,y) =>    (SC, SC+1)
                //! -----------------------------------------------
                if((b0==true && b1==true && b2 ==false) || (b0==true && b1==false && b2==true) || (b0==false && b1==true && b2==true))
                {
                    if(name =="X component")
                    {
                        col = SC;
                        val = tabularDataModel->dataRC(row,col,Qt::EditRole).toDouble();
                    }
                    if(name =="Y component")
                    {
                        if(b0==false) col = SC;
                        else col = SC+1;
                        val = tabularDataModel->dataRC(row,col,Qt::EditRole).toDouble();
                    }
                    if(name =="Z component")
                    {
                        col = SC+1;
                        val = tabularDataModel->dataRC(row,col,Qt::EditRole).toDouble();
                    }
                }
                //! -------------------------------------------------
                //! three columns for the load components are active
                //! -------------------------------------------------
                if(b0==true && b1==true && b2==true)
                {
                    if(name=="X component") col = SC;
                    if(name=="Y component") col = SC+1;
                    if(name=="Z component") col = SC+2;
                    val = tabularDataModel->dataRC(row,col,Qt::EditRole).toDouble();
                };
            }
            */
            Property::loadDefinition theLoadDefinition = QStandardItem::data(Qt::UserRole).value<Property>().getData().value<Property::loadDefinition>();
            switch(theLoadDefinition)
            {
            case Property::loadDefinition_constant: data.setValue(QString("Ramped")); break;
            case Property::loadDefinition_tabularData: data.setValue(QString("Tabular data")); break;
            case Property::loadDefinition_free: data.setValue(QString("Free")); break;
            }
            return data;
        }
        //! ----------------------
        //! Overpressure function
        //! ----------------------
        else if(name=="Overpressure")
        {
            Property::overpressureFunction theOverpressure = QStandardItem::data(Qt::UserRole).value<Property>().getData().value<Property::overpressureFunction>();
            switch(theOverpressure)
            {
            case Property::overpressureFunction_exponential: data.setValue(QString("Exponential")); break;
            case Property::overpressureFunction_linear: data.setValue(QString("Linear")); break;
            case Property::overpressureFunction_tied: data.setValue(QString("Tied")); break;
            }
            return data;
        }
        //! -----------
        //! "Formulation"
        //! -----------
        else if(name=="Formulation")
        {
            Property::contactFormulation theFormulation = QStandardItem::data(Qt::UserRole).value<Property>().getData().value<Property::contactFormulation>();
            switch(theFormulation)
            {
            case Property::contactFormulation_lagrange: data.setValue(QString("Lagrange")); break;
            case Property::contactFormulation_penalty: data.setValue(QString("Pure penalty")); break;
            case Property::contactFormulation_MPC: data.setValue(QString("MPC")); break;
            }
            return data;
        }
        //! -----------
        //! "Behavior"
        //! -----------
        else if(name=="Behavior")
        {
            Property::contactBehavior theBehavior = QStandardItem::data(Qt::UserRole).value<Property>().getData().value<Property::contactBehavior>();
            switch(theBehavior)
            {
            case Property::contactBehavior_asymmetric: data.setValue(QString("Asymmetric")); break;
            case Property::contactBehavior_symmetric: data.setValue(QString("Symmetric")); break;
            }
            return data;
        }
        //! --------------
        //! "Sizing type"
        //! --------------
        else if(name =="Sizing type")
        {
            int sizingType = QStandardItem::data(Qt::UserRole).value<Property>().getData().toInt();
            switch(sizingType)
            {
            case 0: data.setValue(QString("Element size")); break;
            case 1: data.setValue(QString("Number of divisions")); break;
            }
            return data;
        }
        //! -----
        //! Type
        //! -----
        else if(name=="Type")
        {
            switch(this->getCurrentNode()->getType())
            {
            case SimulationNodeClass::nodeType_connectionPair:
            {
                if(QStandardItem::data(Qt::UserRole).value<Property>().getData().isValid()==false)
                {
                    data.setValue(QString("Multiple selection"));
                    return data;
                }

                Property::contactType theType = QStandardItem::data(Qt::UserRole).value<Property>().getData().value<Property::contactType>();
                switch(theType)
                {
                case Property::contactType_bonded: data.setValue(QString("Bonded")); break;
                case Property::contactType_frictional: data.setValue(QString("Frictional")); break;
                case Property::contactType_frictionless: data.setValue(QString("Frictionless")); break;
                case Property::contactType_tied: data.setValue(QString("Tied")); break;
                }
            }
                break;
            }
            return data;
        }
        else if(name =="Number of divisions" && this->getCurrentNode()->getType() == SimulationNodeClass::nodeType_meshEdgeSize)
        {

            int Ndiv = QStandardItem::data(Qt::UserRole).value<Property>().getData().toInt();
            data.setValue(QString("%1").arg(Ndiv));
            return data;
        }
        else if(name =="Element size" && (this->getCurrentNode()->getType() == SimulationNodeClass::nodeType_meshEdgeSize ||
                                          this->getCurrentNode()->getType() == SimulationNodeClass::nodeType_meshFaceSize))
        {
            double eSize = QStandardItem::data(Qt::UserRole).value<Property>().getData().toDouble();
            data.setValue(QString("%1").arg(eSize));
            return data;
        }
        else if(name =="Named selection" || name =="Boundary named selection" || name== "Contact pair")
        {
            void *p = QStandardItem::data(Qt::UserRole).value<Property>().getData().value<void*>();
            QExtendedStandardItem *item = static_cast<QExtendedStandardItem*>(p);
            QString name = item->data(Qt::DisplayRole).toString();
            data.setValue(name);
            return data;
        }
        else if(name == "Coordinate system")
        {
            void *p = QStandardItem::data(Qt::UserRole).value<Property>().getData().value<void*>();
            QExtendedStandardItem *item = static_cast<QExtendedStandardItem*>(p);
            QString name = item->data(Qt::DisplayRole).toString();
            data.setValue(name);
            return data;
        }
        else if(name =="Remote points")
        {
            void *p = QStandardItem::data(Qt::UserRole).value<Property>().getData().value<void*>();
            QExtendedStandardItem *item = static_cast<QExtendedStandardItem*>(p);
            QString name = item->data(Qt::DisplayRole).toString();
            data.setValue(name);
            return data;
        }
        else if (name =="Define by ")   //! this is for coordinate systems
        {
            if(QStandardItem::data(Qt::UserRole).value<Property>().getData().value<Property::defineBy>() == Property::defineBy_geometrySelection)
            {
                data.setValue(QString("Geometry selection"));
            }
            else if(QStandardItem::data(Qt::UserRole).value<Property>().getData().value<Property::defineBy>() == Property::defineBy_globalCoordinates)
            {
                data.setValue(QString("Global coordinates"));
            }
            return data;
        }
        else if (name =="Bolt status")
        {
            switch(QStandardItem::data(Qt::UserRole).value<Property>().getData().value<Property::boltStatusDefinedBy>())
            {
            case Property::boltStatusDefinedBy_load: data.setValue(QString("Load")); break;
            case Property::boltStatusDefinedBy_adjustment: data.setValue(QString("Adjustment")); break;
            case Property::boltStatusDefinedBy_open: data.setValue(QString("Open")); break;
            case Property::boltStatusDefinedBy_lock: data.setValue(QString("Lock")); break;
            }
            return data;
        }
        else if (name =="Define by")
        {
            /*
            switch(this->getCurrentNode()->getType())
            {
            case SimulationNodeClass::nodeType_structuralAnalysisBoltPretension:
                //switch(QStandardItem::data(Qt::UserRole).value<Property>().getData().value<Property::boltStatusDefinedBy>())
                switch(QStandardItem::data(Qt::UserRole).value<Property>().getData().value<Property::defineBy>())
                {
                case Property::boltStatusDefinedBy_load: data.setValue(QString("Load")); break;
                case Property::boltStatusDefinedBy_adjustment: data.setValue(QString("Adjustment")); break;
                case Property::boltStatusDefinedBy_open: data.setValue(QString("Open")); break;
                case Property::boltStatusDefinedBy_lock: data.setValue(QString("Lock")); break;
                }
                break;
            default:
                switch(QStandardItem::data(Qt::UserRole).value<Property>().getData().value<Property::defineBy>())
                {
                case Property::defineBy_components: data.setValue(QString("Components")); break;
                case Property::defineBy_vector: data.setValue(QString("Vector")); break;
                case Property::defineBy_normal: data.setValue(QString("Normal to")); break;
                }
                break;
            }
            */
            switch(QStandardItem::data(Qt::UserRole).value<Property>().getData().value<Property::defineBy>())
            {
            case Property::defineBy_components: data.setValue(QString("Components")); break;
            case Property::defineBy_vector: data.setValue(QString("Vector")); break;
            case Property::defineBy_normal: data.setValue(QString("Normal to")); break;
            }
            return data;
        }
        else if(name =="Load" || name =="Adjustment")
        {
            QString value;
            value.sprintf("%g",QStandardItem::data(Qt::UserRole).value<Property>().getData().toDouble());   
            data.setValue(value);
            return data;
        }
        else if(name =="Ambient" || name=="Diffuse" || name=="Specular" || name =="Transparency")
        {
            QString value;
            value.sprintf("%.1f",QStandardItem::data(Qt::UserRole).value<Property>().getData().toDouble());
            data.setValue(value);
            return data;
        }
        else if(name == "Element control")
        {
            if(QStandardItem::data(Qt::UserRole).value<Property>().getData().value<Property::elementControl>() == Property::elementControl_programControlled)
            {
                data.setValue(QString("Program controlled"));
            }
            else
            {
                data.setValue(QString("Manual"));
            }
            return data;
        }
        else if(name == "Integration scheme")
        {
            if(QStandardItem::data(Qt::UserRole).value<Property>().getData().value<Property::integrationScheme>() == Property::integrationScheme_full)
            {
                data.setValue(QString("Full"));
            }
            else
            {
                data.setValue(QString("Reduced"));
            }
            return data;
        }
        else if(name == "Radial" || name == "Axial" || name == "Tangential")
        {
            if(QStandardItem::data(Qt::UserRole).value<Property>().getData().value<Property::DOFfreedom>() == Property::DOFfreedom_fixed)
            {
                data.setValue(QString("Fixed"));
            }
            else
            {
                data.setValue(QString("Free"));
            }
            return data;
        }
        else if(name =="Center of mass x" || name=="Center of mass y" || name=="Center of mass z")
        {
            QString value;
            value.sprintf("%.3f",QStandardItem::data(Qt::UserRole).value<Property>().getData().toDouble());
            data.setValue(value);
            return data;
        }        
        else if(name=="Moment of inertia Ixx" || name=="Moment of inertia Iyy" || name=="Moment of inertia Izz")
        {
            double d = QStandardItem::data(Qt::UserRole).value<Property>().getData().toDouble();
            QString value = QString("%1").arg(d,0,'E',3);
            data.setValue(value);
            return data;
        }
        else if(name=="Origin X" || name == "Origin Y" || name == "Origin Z" || name =="Offset X" || name =="Offset Y" ||
                name =="Offset Z" || name =="Rotation X" || name =="Rotation Y" || name =="Rotation Z")
        {
            double v = QStandardItem::data(Qt::UserRole).value<Property>().getData().toDouble();
            QString value;
            if(name=="Origin X" || name == "Origin Y" || name == "Origin Z" || name =="Offset X" || name =="Offset Y" || name =="Offset Z")
            {
                value = QString("%1").arg(v);
            }
            else
            {
                value = QString("%1").arg(v).append(" °");
            }
            data.setValue(value);
            return data;
        }
        else if(name=="Length X" || name=="Length Y" || name=="Length Z")
        {
            QString value;
            value.sprintf("%.3f",QStandardItem::data(Qt::UserRole).value<Property>().getData().toDouble());
            data.setValue(value);
            return data;
        }
        else if (name== "Number of nodes")
        {
            //cout<<"Number of nodes: "<<QStandardItem::data(Qt::UserRole).value<Property>().getData().toString().toStdString()<<endl;
            data.setValue(QStandardItem::data(Qt::UserRole).value<Property>().getData().toString());
            if(QStandardItem::data(Qt::UserRole).value<Property>().getData().toInt()==0)
            {
                data.setValue(QString("not meshed"));
            }
            return data;
        }
        else if (name== "Number of elements")
        {
            //cout<<"Number of elements: "<<QStandardItem::data(Qt::UserRole).value<Property>().getData().toString().toStdString()<<endl;
            data.setValue(QStandardItem::data(Qt::UserRole).value<Property>().getData().toString());
            if(QStandardItem::data(Qt::UserRole).value<Property>().getData().toInt()==0)
            {
                data.setValue(QString("not meshed"));
            }
            return data;
        }
        else if(name == "Min element size" || name =="Max element size" || name =="Grading" || name == "Face sizing" || name == "Element size")
        {
            data.setValue(QStandardItem::data(Qt::UserRole).value<Property>().getData().toString());
            return data;
        }        
        else if(name == "Mesh engine 2D")
        {
            switch(QStandardItem::data(Qt::UserRole).value<Property>().getData().value<Property::meshEngine2D>())
            {
            case Property::meshEngine2D_ProgramControlled: data.setValue(QString("Program controlled")); break;
            case Property::meshEngine2D_Netgen: data.setValue(QString("Netgen 2D")); break;
            case Property::meshEngine2D_OCC_ExpressMesh: data.setValue(QString("Express mesh")); break;
            case Property::meshEngine2D_OCC_STL: data.setValue(QString("OCC STL mesh")); break;
            case Property::meshEngine2D_Netgen_STL: data.setValue(QString("Netgen STL")); break;
            }
            return data;
        }
        else if(name == "Mesh engine 3D")
        {
            Property::meshEngine3D meshEng3D = QStandardItem::data(Qt::UserRole).value<Property>().getData().value<Property::meshEngine3D>();
            switch(meshEng3D)
            {
            case Property::meshEngine3D_Netgen: data.setValue(QString("Netgen")); break;
            case Property::meshEngine3D_Tetgen: data.setValue(QString("Tetgen")); break;
            }
            return data;
        }
        else if(name == "Surface mesh type")
        {
            switch(QStandardItem::data(Qt::UserRole).value<Property>().getData().value<Property::meshType_Surface>())
            {
            case Property::meshType_Surface_AllTrig: data.setValue(QString("All triangles")); break;
            case Property::meshType_Surface_QuadDominant: data.setValue(QString("Quad dominant")); break;
            }
            return data;
        }
        else if(name == "Volume mesh type")
        {
            switch(QStandardItem::data(Qt::UserRole).value<Property>().getData().value<Property::meshType_Volume>())
            {
            case Property::meshType_Volume_AllTet: data.setValue(QString("All tet")); break;
            default: data.setValue(QString("not implemented yet")); break;
            }
            return data;
        }
        else if(name == "Mesh order")
        {
            if(QStandardItem::data(Qt::UserRole).value<Property>().getData().value<Property::meshOrder>() == Property::meshOrder_First)
                data.setValue(QString("First"));
            else data.setValue(QString("Second"));
            return data;
        }
        else if(name == "Scoping method")
        {
            if(!QStandardItem::data(Qt::UserRole).value<Property>().getData().isValid())
            {
                data.setValue(QString("Multiple selection"));
                return data;
            }
            switch(QStandardItem::data(Qt::UserRole).value<Property>().getData().value<Property::ScopingMethod>())
            {
            case Property::ScopingMethod_GeometrySelection: data.setValue(QString("Geometry selection")); break;
            case Property::ScopingMethod_NamedSelection: data.setValue(QString("Named selection")); break;
            case Property::ScopingMethod_RemotePoint: data.setValue(QString("Remote point")); break;
            case Property::ScopingMethod_Automatic: data.setValue(QString("Automatic")); break;
            }
            return data;
        }
        else if(name == "Boundary scoping method")
        {
            switch(QStandardItem::data(Qt::UserRole).value<Property>().getData().value<Property::ScopingMethod>())
            {
            case Property::ScopingMethod_GeometrySelection: data.setValue(QString("Geometry selection")); break;
            case Property::ScopingMethod_NamedSelection: data.setValue(QString("Named selection")); break;
            }
            return data;
        }
        else if(name == "Straight sided elements")
        {
            if(QStandardItem::data(Qt::UserRole).value<Property>().getData().toBool()) return QString("Yes");
            else return QString("No");
        }
        else if(name == "Suppressed")
        {
            if(QStandardItem::data(Qt::UserRole).value<Property>().getData().value<Property::SuppressionStatus>() == Property::SuppressionStatus_Active)
                data.setValue(QString("Active"));
            else data.setValue(QString("Suppressed"));
            return data;
        }
        //! -----------------
        //! "Master" "Slave"
        //! -----------------
        else if(name =="Master" || name =="Slave")
        {
            if(!QStandardItem::data(Qt::UserRole).value<Property>().getData().isValid())
            {
                data.setValue(QString("Multiple selection"));
                return data;
            }
            //! -------------------------------------------------------------
            //! a map of mesh data sources is stored into the item: show the
            //! number of mesh
            //! -------------------------------------------------------------
            if(QStandardItem::data(Qt::UserRole).value<Property>().getData().canConvert<IndexedMapOfMeshDataSources>())
            {
                QVariant data;
                const IndexedMapOfMeshDataSources &ds =
                        QStandardItem::data(Qt::UserRole).value<Property>().getData().value<IndexedMapOfMeshDataSources>();
                data.setValue(QString("%1 mesh data sources").arg(ds.values().length()));
                return data;
            }
            //! ------------------------------------------------
            //! a vector of GeometryTag is stored into the item
            //! ------------------------------------------------
            if(QStandardItem::data(Qt::UserRole).value<Property>().getData().canConvert<std::vector<GeometryTag>>())
            {
                QString msg;
                std::vector<GeometryTag> vecLoc = QStandardItem::data(Qt::UserRole).value<Property>().getData().value<std::vector<GeometryTag>>();
                if(vecLoc.size()!=0)
                {
                    msg = QString("%1").arg(vecLoc.size());
                    switch(vecLoc.begin()->subShapeType)
                    {
                    case TopAbs_VERTEX:
                        if(vecLoc.size()>1)msg.append(" vertexes");
                        else msg.append(" vertex");
                        break;
                    case TopAbs_EDGE:
                        if(vecLoc.size()>1)msg.append(" edges");
                        else msg.append(" edge");
                        break;
                    case TopAbs_FACE:
                        if(vecLoc.size()>1)msg.append(" faces");
                        else msg.append(" face");
                        break;
                    case TopAbs_SOLID:
                        if(vecLoc.size()>1)msg.append(" solids");
                        else msg.append(" solid");
                        break;
                    default:
                        break;
                    }
                    data.setValue(msg);
                }
                else
                {
                    msg ="No selection";
                    data.setValue(msg);
                }
                return data;
            }
            else //if(QStandardItem::data(Qt::UserRole).value<Property>().getData().canConvert<void*>())
            {
                //! -----------------------------------------------------------
                //! indirect reference to the list of items in Named Selection
                //! -----------------------------------------------------------
                void *p = QStandardItem::data(Qt::UserRole).value<Property>().getData().value<void*>();
                QExtendedStandardItem *item = static_cast<QExtendedStandardItem*>(p);
                QString name = item->data(Qt::DisplayRole).toString();
                data.setValue(name);
                return data;
            }
        }
        else if(name == "Location")
        {
            QString msg ="Click to change";
            data.setValue(msg);
            return data;
        }
        else if(name == "Geometry")
        {
            std::vector<GeometryTag> vecLoc = QStandardItem::data(Qt::UserRole).value<Property>().getData().value<std::vector<GeometryTag>>();
            if(vecLoc.size()!=0)
            {
                QString shapeType;
                for(int i=0; i<vecLoc.size();i++)
                {
                    switch(vecLoc[0].subShapeType)
                    {
                    case TopAbs_COMPOUND: shapeType = "CC"; break;
                    case TopAbs_COMPSOLID: shapeType = "CS"; break;
                    case TopAbs_SOLID: shapeType = "SO"; break;
                    case TopAbs_FACE: shapeType = "FA"; break;
                    case TopAbs_EDGE: shapeType = "ED"; break;
                    case TopAbs_VERTEX: shapeType = "VE"; break;
                    }
                    data.setValue(QString("%1 %2").arg(shapeType).arg(vecLoc.size()));
                }
            }
            else
            {
                data.setValue(QString("No selection"));
            }
            //data.setValue(val);
            return data;
        }
        else if(name =="Boundary")
        {
            Property::ScopingMethod scopingMethod = this->getCurrentNode()->getPropertyValue<Property::ScopingMethod>("Boundary scoping method");
            switch(scopingMethod)
            {
            case Property::ScopingMethod_GeometrySelection:
            {
                std::vector<GeometryTag> vecLoc = QStandardItem::data(Qt::UserRole).value<Property>().getData().value<std::vector<GeometryTag>>();
                if(vecLoc.size()!=0)
                {
                    QString shapeType;
                    for(int i=0; i<vecLoc.size();i++)
                    {
                        switch(vecLoc[0].subShapeType)
                        {
                        case TopAbs_FACE: shapeType = "FA"; break;
                        }
                        data.setValue(QString("%1 %2").arg(shapeType).arg(vecLoc.size()));
                    }
                }
                else
                {
                    data.setValue(QString("No selection"));
                }
            }
                break;

            case Property::ScopingMethod_NamedSelection:
            {
                void *p = QStandardItem::data(Qt::UserRole).value<Property>().getData().value<void*>();
                QExtendedStandardItem *item = static_cast<QExtendedStandardItem*>(p);
                QString name = item->data(Qt::DisplayRole).toString();
                data.setValue(name);
            }
                break;
            }
            return data;
        }
        else if(name=="Source")
        {
            return QStandardItem::data(Qt::UserRole).value<Property>().getData().toString();
        }
        else if(name=="Visible")
        {
            if(QStandardItem::data(Qt::UserRole).value<Property>().getData().toBool()== true)
                return QString("Yes");
            else return QString("No");
        }
        else if(name=="Number of steps" || name == "Current step number")
        {
            QString value = QString("%1").arg(QStandardItem::data(Qt::UserRole).value<Property>().getData().toInt());
            data.setValue(value);
            return value;
        }
        else if(name=="Step end time")
        {
            QString value = QString("%1").arg(QStandardItem::data(Qt::UserRole).value<Property>().getData().toDouble());
            data.setValue(value);
            return value;
        }
        else if(name=="Friction coefficient" || name =="C0" || name =="Normal stiffness" || name =="Tau")
        {
            QString value;
            value.sprintf("%.3f",QStandardItem::data(Qt::UserRole).value<Property>().getData().toDouble());
            data.setValue(value);
            return data;
        }
        else if(name == "Auto time stepping")
        {
            if(QStandardItem::data(Qt::UserRole).value<Property>().getData().value<Property::autoTimeStepping>() == Property::autoTimeStepping_ProgramControlled)
            {
                data.setValue(QString("Program Controlled"));
            }
            else if (QStandardItem::data(Qt::UserRole).value<Property>().getData().value<Property::autoTimeStepping>() == Property::autoTimeStepping_ON)
            {
                data.setValue(QString("On"));
            }
            else if (QStandardItem::data(Qt::UserRole).value<Property>().getData().value<Property::autoTimeStepping>() == Property::autoTimeStepping_OFF)
            {
                data.setValue(QString("Off"));
            }
            return data;
        }
        else if(name == "Solver type")
        {
            switch(QStandardItem::data(Qt::UserRole).value<Property>().getData().value<Property::solverType>())
            {
            case Property::solverType_programControlled: data.setValue(QString("Program Controlled")); break;
            case Property::solverType_direct: data.setValue(QString("Direct")); break;
            case Property::solverType_iterative: data.setValue(QString("Iterative")); break;
            }
            return data;
        }
        else if(name == "Time step size")
        {
            double val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toDouble();
            data.setValue(QString("%1").arg(val, 0, 'E', 3));
            return data;
        }
        else if(name == "Particle mass" || name =="Intensity" || name == "Electric charge")
        {
            double val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toDouble();
            data.setValue(QString("%1").arg(val, 0, 'E', 3));
            return data;
        }
        else if(name =="Potential")
        {
            double val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toDouble();
            data.setValue(QString("%1").arg(val));
            return data;
        }
        else if(name =="Emitter")
        {
            int val = QStandardItem::data(Qt::UserRole).value<Property>().getData().toInt();
            data.setValue(val==0? QString("Off"):QString("On"));
            return data;
        }
        else if (name1 == "separator")
        {
            //!cout<<"name1: "<<name1.toStdString()<<endl;
            return QStandardItem::data(Qt::DisplayRole).toString();
        }
        else return QStandardItem::data(Qt::DisplayRole).toString();
    }
    else return QStandardItem::data(role);
}


//! ------------------
//! function: getIcon
//! details:
//! ------------------
QIcon QExtendedStandardItem::getIcon(SimulationNodeClass::nodeType theNodeType) const
{
    switch(theNodeType)
    {
    case SimulationNodeClass::nodeType_combinedAnalysis: return QIcon(":/icons/icon_combined analysis.png"); break;
    case SimulationNodeClass::nodeType_particlesInFieldsAnalysis: return QIcon(":/icons/icon_a field.png"); break;
    case SimulationNodeClass::nodeType_meshMeshMetric: return QIcon(":/icons/icon_metric.png"); break;
    case SimulationNodeClass::nodeType_meshMeshType: return QIcon(":/icons/icon_type of mesh.png"); break;
    case SimulationNodeClass::nodeType_modelChange: return QIcon(":/icons/icon_balordo.png"); break;
    case SimulationNodeClass::nodeType_postObject: return QIcon(":/icons/icon_locator.png"); break;
    case SimulationNodeClass::nodeType_solutionStructuralTemperature: return QIcon(":/icons/icon_insert body temperature dist.png"); break;
    case SimulationNodeClass::nodeType_geometry: return QIcon(":/icons/icon_geometry root item.png"); break;
    case SimulationNodeClass::nodeType_import: return QIcon(":/icons/icon_mesh method.png"); break;
    case SimulationNodeClass::nodeType_geometryBody:
    case SimulationNodeClass::nodeType_geometryPart:
        return QIcon(":/icons/icon_body3D.png"); break;
    case SimulationNodeClass::nodeType_meshControl: return QIcon(":/icons/icon_volume mesh.png"); break;
    case SimulationNodeClass::nodeType_meshBodyMeshControl: return QIcon(":/icons/icon_volume mesh.png"); break;
    case SimulationNodeClass::nodeType_meshBodyMeshMethod: return QIcon(":/icons/icon_mesh method.png"); break;
    case SimulationNodeClass::nodeType_meshMethod: return QIcon(":/icons/icon_mesh method.png"); break;
    case SimulationNodeClass::nodeType_meshFaceSize: return QIcon(":/icons/icon_mesh face sizing.png"); break;
    case SimulationNodeClass::nodeType_meshEdgeSize: return QIcon(":/icons/icon_edge sizing.png"); break;
    case SimulationNodeClass::nodeType_meshVertexSize: return QIcon(":/icons/icon_point.png"); break;
    case SimulationNodeClass::nodeType_thermalAnalysis: return QIcon(":/icons/icon_thermal analysis.png"); break;

        //! --------------------
        //! "Analysis settings"
        //! --------------------
    case SimulationNodeClass::nodeType_thermalAnalysisSettings:
    case SimulationNodeClass::nodeType_structuralAnalysisSettings:
    case SimulationNodeClass::nodeType_combinedAnalysisSettings:
    case SimulationNodeClass::nodeType_particlesInFieldsAnalysisSettings:
        return QIcon(":/icons/icon_analysis settings.png"); break;

    case SimulationNodeClass::nodeType_root: return QIcon(":/icons/icon_model root item.png"); break;
    case SimulationNodeClass::nodeType_connection:
    case SimulationNodeClass::nodeType_connectionPair: return QIcon(":/icons/icon_insert contact.png"); break;
    case SimulationNodeClass::nodeType_connectionGroup: return QIcon(":/icons/icon_folder.png"); break;
    case SimulationNodeClass::nodeType_structuralAnalysis: return QIcon(":/icons/icon_static structural.png"); break;
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_FrictionlessSupport: return QIcon(":/icons/icon_BC frictionless support.png"); break;
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_CylindricalSupport: return QIcon(":/icons/icon_BC cylindrical support.png"); break;
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryContidion_FixedSupport: return QIcon(":/icons/icon_BC fixed support.png"); break;
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_CompressionOnlySupport: return QIcon(":/icons/icon_compression only support.png"); break;
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Displacement: return QIcon(":/icons/icon_BC displacement.png"); break;
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Moment: return QIcon(":/icons/icon_BC moment.png"); break;
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Force:
    case SimulationNodeClass::nodeType_solutionStructuralNodalForces: return QIcon(":/icons/icon_BC force.png"); break;
    case SimulationNodeClass::nodeType_solutionStructuralReactionForce: return QIcon(":/icons/icon_BC force.png"); break;
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Pressure: return QIcon(":/icons/icon_BC pressure.png"); break;
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_ImportedTemperatureDistribution: return QIcon(":/icons/icon_insert body temperature dist.png");
    case SimulationNodeClass::nodeType_namedSelection: return QIcon(":/icons/icon_named selection root item.png");break;
    case SimulationNodeClass::nodeType_namedSelectionElement: return QIcon(":/icons/icon_volume mesh.png"); break;
    case SimulationNodeClass::nodeType_namedSelectionGeometry: return QIcon(":/icons/icon_named selection geometry.png"); break;
    case SimulationNodeClass::nodeType_coordinateSystems: case SimulationNodeClass::nodeType_coordinateSystem:
    case SimulationNodeClass::nodeType_coordinateSystem_global: return QIcon(":/icons/icon_system of reference.png"); break;
    case SimulationNodeClass::nodeType_solutionStructuralGamma: return QIcon(":/icons/icon_BC force.png"); break;
        //! -----------
        //! "Solution"
        //! -----------
    case SimulationNodeClass::nodeType_StructuralAnalysisSolution:
    case SimulationNodeClass::nodeType_thermalAnalysisSolution:
    case SimulationNodeClass::nodeType_combinedAnalysisSolution:
    case SimulationNodeClass::nodeType_particlesInFieldsSolution:
        return QIcon(":/icons/icon_solution.png"); break;

    case SimulationNodeClass::nodeType_structuralAnalysisThermalCondition:return QIcon(":/icons/icon_thermal analysis.png"); break;
    case SimulationNodeClass::nodeType_remotePoint:
    case SimulationNodeClass::nodeType_remotePointRoot: return QIcon(":/icons/icon_remote point.png"); break;
    case SimulationNodeClass::nodeType_mapper: return QIcon(":/icons/icon_mapping.png"); break;
    case SimulationNodeClass::nodeType_OpenFoamScalarData: return QIcon(":/icons/icon_openfoam.png"); break;
    case SimulationNodeClass::nodeType_importedBodyScalar: return QIcon(":/icons/icon_insert body temperature dist.png"); break;
    case SimulationNodeClass::nodeType_solutionStructuralNodalDisplacement: return QIcon(":/icons/icon_deformation.png"); break;
    case SimulationNodeClass::nodeType_solutionStructuralStress: return QIcon(":/icons/icon_spring.png"); break;
    case SimulationNodeClass::nodeType_solutionStructuralTotalStrain: return QIcon(":/icons/icon_spring.png"); break;

        //! -----------------------
        //! "Solution information"
        //! -----------------------
    case SimulationNodeClass::nodeType_StructuralAnalysisSolutionInformation:
    case SimulationNodeClass::nodeType_thermalAnalysisSolutionInformation:
    case SimulationNodeClass::nodeType_combinedAnalysisSolutionInformation:
    case SimulationNodeClass::nodeType_particlesInFieldsSolutionInformation:
        return QIcon(":/icons/icon_solution information.png"); break;

    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Acceleration: return QIcon(":/icons/icon_acceleration.png"); break;
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RotationalVelocity: return QIcon(":/icons/icon_rotational velocity.png"); break;
    case SimulationNodeClass::nodeType_structuralAnalysisBoltPretension: return QIcon(":/icons/icon_bolt.png"); break;
    case SimulationNodeClass::nodeType_repairTool: return QIcon(":/icons/icon_repair tool.png"); break;
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteForce: return QIcon(":/icons/icon_remote force.png"); break;
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteDisplacement: return QIcon(":/icons/icon_remote displacement.png"); break;
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteRotation: return QIcon(":/icons/icon_rotation.png"); break;
    case SimulationNodeClass::nodeType_meshPrismaticLayer: return QIcon(":/icons/icon_prismatic layer.png"); break;
    case SimulationNodeClass::nodeType_solutionStructuralFatigueTool:return QIcon(":/icons/icon_crack.png"); break;
    case SimulationNodeClass::nodeType_thermalAnalysisTemperature: return QIcon(":/icons/icon_temperature.png"); break;
    case SimulationNodeClass::nodeType_thermalAnalysisConvection: return QIcon(":/icons/icon_convection.png"); break;
    case SimulationNodeClass::nodeType_thermalAnalysisThermalFlux: return QIcon(":/icons/icon_thermal flux.png"); break;
    case SimulationNodeClass::nodeType_thermalAnalysisThermalFlow: return QIcon(":/icons/icon_thermal flow.png"); break;
    case SimulationNodeClass::nodeType_thermalAnalysisThermalPower: return QIcon(":/icons/icon_thermal power.png"); break;
    case SimulationNodeClass::nodeType_thermalAnalysisAdiabaticWall: return QIcon(":/icons/icon_adiabatic wall.png"); break;
    case SimulationNodeClass::nodeType_solutionThermalTemperature: return QIcon(":/icons/icon_temperature.png"); break;
    case SimulationNodeClass::nodeType_solutionThermalFlux: return QIcon(":/icons/icon_thermal flux.png"); break;
    case SimulationNodeClass::nodeType_electrostaticPotential: return QIcon(":/icons/icon_electrostatic potential.png"); break;
    case SimulationNodeClass::nodeType_pointMass: return QIcon(":/icons/icon_point mass.png"); break;

#ifdef COSTAMP_VERSION
    case SimulationNodeClass::nodeType_timeStepBuilder: return QIcon(":/icons/icon_clock.png"); break;
    case SimulationNodeClass::nodeType_processParametersClosureForce: return QIcon(":/icons/icon_closure force.png"); break;
    case SimulationNodeClass::nodeType_processParametersPressure: return QIcon(":/icons/icon_BC pressure.png"); break;
#endif
    default:
        return QIcon();
        break;
    }    
}

Q_DECLARE_OPAQUE_POINTER(QExtendedStandardItem*)
