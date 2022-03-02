//! ----------------
//! custom includes
//! ----------------
#include "markerbuilder.h"
#include "markers.h"
#include "src/utils/geomtoolsclass.h"
#include "qextendedstandarditem.h"
#include "ext/occ_extended/handle_ais_doublearrowmarker_reg.h"
#include "ext/occ_extended/handle_ais_customtrihedron_reg.h"
#include "src/main/simulationmanager.h"
#include "src/gui/tabularData/customtablemodel.h"
#include "src/main/maintreetools.h"
#include "occPreGLwidget.h"

//! ---
//! Qt
//! ---
#include <Qvector>
#include <QVariant>
#include <QStandardItemModel>

//! ----
//! OCC
//! ----
#include <Bnd_Box.hxx>
#include <BRepBndLib.hxx>
#include <TopoDS_Builder.hxx>
#include <gp_Pnt.hxx>
#include <GProp_GProps.hxx>
#include <BRepGProp.hxx>

//! ------------------------------------------------------------
//! function: addMarker
//! details:  add a marker to the node as "Graphic object"
//!           returns true if the marker has been added/updated
//! ------------------------------------------------------------
bool markerBuilder::addMarker(SimulationNodeClass *node, geometryDataBase *gDB)
{
    cout<<"markerBuilder::addMarker()->____function called____"<<endl;
    SimulationNodeClass::nodeType nodeType = node->getType();
    switch(nodeType)
    {
    case SimulationNodeClass::nodeType_structuralAnalysisBoltPretension:
    {
        //! -------------------------------------------
        //! retrieve the coordinate system of the bolt
        //! -------------------------------------------
        void *CSpointer = node->getPropertyValue<void*>("Coordinate system");
        QStandardItem *itemCS = static_cast<QStandardItem*>(CSpointer);
        SimulationNodeClass *nodeCS = itemCS->data(Qt::UserRole).value<SimulationNodeClass*>();

        double x0,y0,z0;
        double nx,ny,nz;
        double D;

        std::vector<GeometryTag> locs = node->getPropertyValue<std::vector<GeometryTag>>("Tags");
        if(locs.size()==0) return false;

        //! -----------------------------------------------
        //! only solids are acceptable for bolt pretension
        //! -----------------------------------------------
        if(locs.size()==1 && locs[0].subShapeType==TopAbs_SOLID)
        {
            const TopoDS_Shape &theBoltsShape = gDB->bodyMap.value(locs.at(0).subTopNr);
            Bnd_Box BB;
            BRepBndLib::Add(theBoltsShape, BB);
            double diag = sqrt(BB.SquareExtent());
            D = diag/40.0;

            if(nodeCS->getType()==SimulationNodeClass::nodeType_coordinateSystem_global)
            {
                x0 = y0 = z0 = 0.0;
                nx = ny = 0.0; nz = 1.0;
            }
            else
            {
                //! --------------------------------
                //! middle point =><= of the marker
                //! --------------------------------
                const QVector<double> &baseOrigin = nodeCS->getPropertyValue<QVector<double>>("Base origin");
                x0 = baseOrigin.at(0);
                y0 = baseOrigin.at(1);
                z0 = baseOrigin.at(2);

                //! -----------------
                //! axis of the bolt
                //! -----------------
                const QVector<QVector<double>> &baseDirectionalData =
                        nodeCS->getPropertyValue<QVector<QVector<double>>>("Base directional data");
                nx = baseDirectionalData.at(2).at(0);
                ny = baseDirectionalData.at(2).at(1);
                nz = baseDirectionalData.at(2).at(2);
            }

            //cout<<"____Diameter: "<<D<<"____";
            //cout<<"____Center: ("<<x0<<", "<<y0<<", "<<z0<<")____"<<endl;
            //cout<<"____Direction: ("<<nx<<", "<<ny<<", "<<nz<<")____"<<endl;

            //! ---------------------------
            //! actually create the marker
            //! ---------------------------
            gp_Pnt boltOrigin(x0,y0,z0);
            gp_Dir boltDirection(nx,ny,nz);
            const AIS_DoubleArrowMarker_handle_reg &boltMarker = markers::buildDoubleArrowMarker(boltOrigin,boltDirection,D);
            QVariant data;
            data.setValue(boltMarker);
            Property prop_marker("Graphic object",data,Property::PropertyGroup_GraphicObjects);

            //! ---------------------------------------
            //! add the graphic object if not existing
            //! otherwise replace it
            //! ---------------------------------------
            if(node->getPropertyItem("Graphic object")==Q_NULLPTR) node->addProperty(prop_marker);
            else node->replaceProperty("Graphic object",prop_marker);
            return true;
        }
    }
        break;

    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Acceleration:
    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Force:
    {
        const std::vector<GeometryTag> &locs = node->getPropertyValue<std::vector<GeometryTag>>("Geometry");
        //! ---------------------------------------------------------
        //! create the marker only if something is present in "Tags"
        //! ---------------------------------------------------------
        if(locs.size()==0)
        {
            cout<<"markerBuilder::addMarker()->____locs empty: returning false____"<<endl;
            return false;
        }
        //! -----------------------------------------------
        //! retrieve the center of mass - build a compound
        //! using the shapes recorded in tags
        //! -----------------------------------------------
        TopoDS_Compound compound;
        TopoDS_Builder theBuilder;
        theBuilder.MakeCompound(compound);
        for(int i=0; i<locs.size(); i++)
        {
            const GeometryTag &loc = locs.at(i);
            const TopoDS_Shape &curShape = gDB->bodyMap.value(loc.parentShapeNr);
            theBuilder.Add(compound,curShape);
        }

        gp_Pnt CM = GeomToolsClass::getCenterOfMass(compound);

        //! -------------------------------------------
        //! retrieve the direction of the acceleration
        //! -------------------------------------------
        gp_Dir dir;

        //! ------------------------------------------
        //! read the components from the tabular data
        //! ------------------------------------------
        SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
        QTreeView *theTree = sm->myTreeView;
        SimulationNodeClass *nodeAnalysisSetting = mainTreeTools::getAnalysisSettingsNodeFromCurrentItem(theTree);
        CustomTableModel *tabModel = nodeAnalysisSetting->getTabularDataModel();
        QList<int> tableColumns = mainTreeTools::getColumnsToRead(theTree,tabModel->getColumnBeforeBC());

        Property::defineBy defineBy = node->getPropertyValue<Property::defineBy>("Define by");
        int curStepNumber = nodeAnalysisSetting->getPropertyValue<int>("Current step number");

        switch(defineBy)
        {
        case Property::defineBy_components:
        {
            void *CS = node->getPropertyValue<void*>("Coordinate system");

            QStandardItem *itemCS = static_cast<QStandardItem*>(CS);
            SimulationNodeClass *nodeCS = itemCS->data(Qt::UserRole).value<SimulationNodeClass*>();
            cout<<"____"<<nodeCS->getName().toStdString()<<"____"<<endl;
            cout<<"____"<<nodeCS->type().toStdString()<<"____"<<endl;
            cout<<nodeCS->getType()<<endl;
            QVector<QVector<double>> directionalData;
            switch(nodeCS->getType())
            {
            case SimulationNodeClass::nodeType_coordinateSystem_global:
            {
                QVector<double> XAxisData = nodeCS->getPropertyValue<QVector<double>>("X axis data");
                QVector<double> YAxisData = nodeCS->getPropertyValue<QVector<double>>("Y axis data");
                QVector<double> ZAxisData = nodeCS->getPropertyValue<QVector<double>>("Z axis data");

                directionalData.push_back(XAxisData);
                directionalData.push_back(YAxisData);
                directionalData.push_back(ZAxisData);
                cout<<XAxisData.at(0)<<" "<<XAxisData.at(1)<<" "<<XAxisData.at(2)<<endl;
                cout<<YAxisData.at(0)<<" "<<YAxisData.at(1)<<" "<<YAxisData.at(2)<<endl;
                cout<<ZAxisData.at(0)<<" "<<ZAxisData.at(1)<<" "<<ZAxisData.at(2)<<endl;
            }
                break;

            case SimulationNodeClass::nodeType_coordinateSystem:
            {
                if(nodeCS->getPropertyItem("Base directional data")==NULL) exit(1);
                directionalData = nodeCS->getPropertyValue<QVector<QVector<double>>>("Base directional data");
            }
                break;
            }

            double Xcomp_local = tabModel->dataRC(curStepNumber,tableColumns.at(0),Qt::EditRole).toDouble();
            double Ycomp_local = tabModel->dataRC(curStepNumber,tableColumns.at(1),Qt::EditRole).toDouble();
            double Zcomp_local = tabModel->dataRC(curStepNumber,tableColumns.at(2),Qt::EditRole).toDouble();

            double Xcomp = Xcomp_local*directionalData.at(0).at(0)+Ycomp_local*directionalData.at(1).at(0)+Zcomp_local*directionalData.at(2).at(0);
            double Ycomp = Xcomp_local*directionalData.at(0).at(1)+Ycomp_local*directionalData.at(1).at(1)+Zcomp_local*directionalData.at(2).at(1);
            double Zcomp = Xcomp_local*directionalData.at(0).at(2)+Ycomp_local*directionalData.at(1).at(2)+Zcomp_local*directionalData.at(2).at(2);

            //! ---------------------------------------------
            //! check if the acceleration module is not zero
            //! ---------------------------------------------
            if(sqrt(pow(Xcomp,2)+pow(Ycomp,2)+pow(Zcomp,2))<1e-6)
            {
                if(node->getPropertyItem("Graphic object")!=NULL) node->removeProperty("Graphic object");
                return false;
            }
            gp_Vec v(Xcomp,Ycomp,Zcomp);
            dir = gp_Dir(v);
        }
            break;

        case Property::defineBy_vector:
        {
            //! ---------------------------------------------
            //! check if the acceleration module is not NULL
            //! ---------------------------------------------
            double magnitude = tabModel->dataRC(curStepNumber,tableColumns.at(0),Qt::EditRole).toDouble();
            if(magnitude<1e-6 && magnitude>-1e-6)
            {
                cerr<<"____zero magnitude: returning false____"<<endl;
                if(node->getPropertyItem("Graphic object")!=NULL) node->removeProperty("Graphic object");
                return false;
            }

            QVector<double> direction = node->getPropertyValue<QVector<double>>("Direction");
            double nx = direction.at(3);
            double ny = direction.at(4);
            double nz = direction.at(5);
            gp_Vec v(nx,ny,nz);
            dir = gp_Dir(v);
        }
            break;
        }

        //! ----------------------
        //! diameter of the arrow
        //! ----------------------
        Bnd_Box BB;
        BRepBndLib::Add(compound, BB);
        double diag = sqrt(BB.SquareExtent());
        double D = diag/30.0;

        //! --------------------------
        //! actually create the arrow
        //! --------------------------
        AIS_ArrowMarker_handle_reg theArrow =  markers::buildArrowMarker(CM,dir,D);
        QVariant data;
        data.setValue(theArrow);
        Property prop_marker("Graphic object",data,Property::PropertyGroup_GraphicObjects);

        //! ---------------------------------------
        //! add the graphic object if not existing
        //! otherwise replace it
        //! ---------------------------------------
        if(node->getPropertyItem("Graphic object")==NULL) node->addProperty(prop_marker);
        else node->replaceProperty("Graphic object",prop_marker);
    }
        break;

    case SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Moment:
    {
        const std::vector<GeometryTag> &locs = node->getPropertyValue<std::vector<GeometryTag>>("Geometry");
        //! ---------------------------------------------------------
        //! create the marker only if something is present in "Tags"
        //! ---------------------------------------------------------
        if(locs.size()==0)
        {
            cout<<"markerBuilder::addMarker()->____locs empty: returning false____"<<endl;
            return false;
        }

        //! ---------------------------
        //! area of the selected faces
        //! ---------------------------
        TopoDS_Compound compound;
        TopoDS_Builder theBuilder;
        theBuilder.MakeCompound(compound);
        cout<<"markerBuilder::addMarker()->____tag00____"<<endl;
        //double Area;
        for(int i=0; i<locs.size(); i++)
        {
            int bodyIndex = locs.at(i).parentShapeNr;
            int faceIndex = locs.at(i).subTopNr;
            cout<<"markerBuilder::addMarker()->____BI____"<<bodyIndex<<endl;
            cout<<"markerBuilder::addMarker()->____FI____"<<faceIndex<<endl;

            const TopoDS_Shape &curFace = gDB->MapOfBodyTopologyMap.value(bodyIndex).faceMap.FindKey(faceIndex);
            //GProp_GProps prop;
            //BRepGProp::SurfaceProperties(compound/*curFace*/,prop);
            cout<<"markerBuilder::addMarker()->____tag01"<<endl;

            //Area += prop.Mass();
            theBuilder.Add(compound,curFace);
            //gp_Pnt p = prop.CentreOfMass();

        }
        cout<<"markerBuilder::addMarker()->____tag01"<<endl;
        GProp_GProps prop;
        BRepGProp::SurfaceProperties(compound,prop);
        cout<<"markerBuilder::addMarker()->____tag02"<<endl;

        double Area = prop.Mass();
        cout<<"markerBuilder::addMarker()->____mass____"<<Area<<endl;

        //! ----------------------------------------------
        //! this sets the size of the curved arrow marker
        //! ----------------------------------------------
        double Rin = sqrt(Area)/10.0;

        //! -----------------
        //! marker placement
        //! -----------------
        gp_Ax2 axes;
        gp_Pnt P = GeomToolsClass::getCenterOfMass(compound)/* prop.CentreOfMass()*/;
        axes.SetLocation(P);

        //! ------------------------------------------------------
        //! this for reading the components from the tabular data
        //! ------------------------------------------------------
        SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
        QTreeView *theTree = sm->myTreeView;
        SimulationNodeClass *nodeAnalysisSetting = mainTreeTools::getAnalysisSettingsNodeFromCurrentItem(theTree);
        CustomTableModel *tabModel = nodeAnalysisSetting->getTabularDataModel();
        QList<int> tableColumns = mainTreeTools::getColumnsToRead(theTree,tabModel->getColumnBeforeBC());
        Property::defineBy defineBy = node->getPropertyValue<Property::defineBy>("Define by");

        int curStepNumber = nodeAnalysisSetting->getPropertyValue<int>("Current step number");
        cout<<"markerBuilder::addMarker()->____ts____"<<curStepNumber<<endl;

        switch(defineBy)
        {
        case Property::defineBy_vector:
        {
            //! -------------------------------------
            //! check that the magnitude is not zero
            //! -------------------------------------
            int tableColumn = tableColumns.at(0);
            double magnitude = tabModel->dataRC(curStepNumber,tableColumn,Qt::EditRole).toDouble();
            if(fabs(magnitude)<1e-6)
            {
                cerr<<"____magnitude zero: returning false____"<<endl;
                if(node->getPropertyItem("Graphic object")!=NULL) node->removeProperty("Graphic object");
                return false;
            }
            QVector<double> direction = node->getPropertyValue<QVector<double>>("Direction");
            double nx = direction.at(3);
            double ny = direction.at(4);
            double nz = direction.at(5);

            gp_Vec v(nx,ny,nz);
            gp_Dir dir = gp_Dir(v);
            axes.SetDirection(dir);
        }
            break;

        case Property::defineBy_components:
        {
            double Xcomp_local = tabModel->dataRC(curStepNumber,tableColumns.at(0),Qt::EditRole).toDouble();
            double Ycomp_local = tabModel->dataRC(curStepNumber,tableColumns.at(1),Qt::EditRole).toDouble();
            double Zcomp_local = tabModel->dataRC(curStepNumber,tableColumns.at(2),Qt::EditRole).toDouble();

            //! -----------------------------------
            //! check if the magnitude is not zero
            //! -----------------------------------
            if(sqrt(pow(Xcomp_local,2)+pow(Ycomp_local,2)+pow(Zcomp_local,2))<1e-6)
            {
                cerr<<"____the magnitude is zero: returning false____"<<endl;
                if(node->getPropertyItem("Graphic object")!=NULL) node->removeProperty("Graphic object");
                return false;
            }

            QVector<QVector<double>> directionalData;

            void *CS = node->getPropertyValue<void*>("Coordinate system");
            QStandardItem *itemCS = static_cast<QStandardItem*>(CS);
            SimulationNodeClass *nodeCS = itemCS->data(Qt::UserRole).value<SimulationNodeClass*>();
            switch(nodeCS->getType())
            {
            case SimulationNodeClass::nodeType_coordinateSystem_global:
            {
                QVector<double> XAxisData = nodeCS->getPropertyValue<QVector<double>>("X axis data");
                QVector<double> YAxisData = nodeCS->getPropertyValue<QVector<double>>("Y axis data");
                QVector<double> ZAxisData = nodeCS->getPropertyValue<QVector<double>>("Z axis data");

                directionalData.push_back(XAxisData);
                directionalData.push_back(YAxisData);
                directionalData.push_back(ZAxisData);
            }
                break;

            case SimulationNodeClass::nodeType_coordinateSystem:
            {
                const QVector<QVector<double>> &baseDirectionalData = nodeCS->getPropertyValue<QVector<QVector<double>>>("Base directional data");
                directionalData.push_back(baseDirectionalData.at(0));
                directionalData.push_back(baseDirectionalData.at(1));
                directionalData.push_back(baseDirectionalData.at(2));
            }
                break;
            }

            double Xcomp = Xcomp_local*directionalData.at(0).at(0)+Ycomp_local*directionalData.at(1).at(0)+Zcomp_local*directionalData.at(2).at(0);
            double Ycomp = Xcomp_local*directionalData.at(0).at(1)+Ycomp_local*directionalData.at(1).at(1)+Zcomp_local*directionalData.at(2).at(1);
            double Zcomp = Xcomp_local*directionalData.at(0).at(2)+Ycomp_local*directionalData.at(1).at(2)+Zcomp_local*directionalData.at(2).at(2);

            gp_Vec v(Xcomp,Ycomp,Zcomp);
            gp_Dir dir = gp_Dir(v);

            axes.SetDirection(dir);
        }
            break;
        }
        cout<<"markerBuilder::addMarker()->____tag01____"<<endl;

        //! ---------------------------
        //! actually create the marker
        //! ---------------------------
        AIS_CurvedArrowMarker_handle_reg marker = markers::buildCurvedArrow(axes,Rin,false);

        QVariant data;
        data.setValue(marker);
        Property prop_marker("Graphic object",data,Property::PropertyGroup_GraphicObjects);

        //! ---------------------------------------
        //! add the graphic object if not existing
        //! otherwise replace it
        //! ---------------------------------------
        if(node->getPropertyItem("Graphic object")==NULL) node->addProperty(prop_marker);
        else node->replaceProperty("Graphic object",prop_marker);
        cout<<"markerBuilder::addMarker()->____tag02____"<<endl;

    }
        break;

        //! -----------
        //! Point mass
        //! -----------
    case SimulationNodeClass::nodeType_pointMass:
    {
        //! -------------------------------
        //! coordinate of the remote point
        //! -------------------------------
        double x = node->getPropertyValue<double>("X coordinate");
        double y = node->getPropertyValue<double>("Y coordinate");
        double z = node->getPropertyValue<double>("Z coordinate");

        //! ---------------------
        //! calculate the radius
        //! ---------------------
        double X,Y;
        occPreGLWidget *viewPort = static_cast<occPreGLWidget*>(tools::getWidgetByName("maingwindow"));
        viewPort->getView()->Size(X,Y);
        double L = sqrt(X*Y);
        double radius = L/100.0;

        AIS_SphereMarker_handle_reg marker = markers::buildSphereMarker(gp_Pnt(x,y,z),radius);
        QVariant data;
        data.setValue(marker);
        Property prop_marker("Graphic object",data,Property::PropertyGroup_GraphicObjects);

        //! ---------------------------------------
        //! add the graphic object if not existing
        //! otherwise replace it
        //! ---------------------------------------
        if(node->getPropertyItem("Graphic object")==Q_NULLPTR) node->addProperty(prop_marker);
        else node->replaceProperty("Graphic object",prop_marker);
    }
        break;

        //! -------------
        //! remote point
        //! -------------
    case SimulationNodeClass::nodeType_remotePoint:
    {
        //! -------------------------------
        //! coordinate of the remote point
        //! -------------------------------
        double x = node->getPropertyValue<double>("X abs coordinate");
        double y = node->getPropertyValue<double>("Y abs coordinate");
        double z = node->getPropertyValue<double>("Z abs coordinate");

        //cout<<"markerBuilder::addMarker()->____adding remote point marker at ("<<x<<", "<<y<<", "<<z<<")____"<<endl;

        //! ---------------------
        //! calculate the radius
        //! ---------------------
        double X,Y;
        occPreGLWidget *viewPort = static_cast<occPreGLWidget*>(tools::getWidgetByName("maingwindow"));
        viewPort->getView()->Size(X,Y);
        double L = sqrt(X*Y);
        double radius = L/100.0;

        AIS_SphereMarker_handle_reg marker = markers::buildSphereMarker(gp_Pnt(x,y,z),radius);
        QVariant data;
        data.setValue(marker);
        Property prop_marker("Graphic object",data,Property::PropertyGroup_GraphicObjects);

        //! ---------------------------------------
        //! add the graphic object if not existing
        //! otherwise replace it
        //! ---------------------------------------
        if(node->getPropertyItem("Graphic object")==NULL) node->addProperty(prop_marker);
        else node->replaceProperty("Graphic object",prop_marker);
    }
        break;
    }
    cout<<"markerBuilder::addMarker()->____exiting function____"<<endl;
    return true;
}
