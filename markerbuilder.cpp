//! ----------------
//! custom includes
//! ----------------
#include "markerbuilder.h"
#include "markers.h"
#include "geomtoolsclass.h"
#include "qextendedstandarditem.h"
#include "handle_ais_doublearrowmarker_reg.h"
#include "handle_ais_customtrihedron_reg.h"
#include "simulationmanager.h"
#include "customtablemodel.h"
#include "maintreetools.h"
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
/*
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

        QVector<GeometryTag> locs = node->getPropertyValue<QVector<GeometryTag>>("Tags");
        if(locs.isEmpty()) return false;

        //! -----------------------------------------------
        //! only solids are acceptable for bolt pretension
        //! -----------------------------------------------
        if(locs.length()==1 && locs.at(0).subShapeType==TopAbs_SOLID)
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
    {
        const QVector<GeometryTag> &locs = node->getPropertyValue<QVector<GeometryTag>>("Geometry");
        //! ---------------------------------------------------------
        //! create the marker only if something is present in "Tags"
        //! ---------------------------------------------------------
        if(locs.isEmpty())
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
        for(int i=0; i<locs.length(); i++)
        {
            const GeometryTag &loc = locs.at(i);
            const TopoDS_Shape &curShape = gDB->bodyMap.value(loc.parentShapeNr);
            theBuilder.Add(compound,curShape);
        }
        gp_Pnt CM = GeomToolsClass::getCenterOfMass(compound);
        //cout<<"____CM ("<<CM.X()<<", "<<CM.Y()<<", "<<CM.Z()<<")____"<<endl;

        //! -------------------------------------------
        //! retrieve the direction of the acceleration
        //! -------------------------------------------
        gp_Dir dir;

        //! ------------------------------------------
        //! read the components from the tabular data
        //! ------------------------------------------
        SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
        QTreeView *theTree = sm->myTreeView;
        QList<int> tableColumns = mainTreeTools::getColumnsToRead(theTree);
        SimulationNodeClass *nodeAnalysisSetting = sm->getAnalysisSettingsNodeFromCurrentItem();
        CustomTableModel *tabModel = nodeAnalysisSetting->getTabularDataModel();
        Property::defineBy defineBy = node->getPropertyValue<Property::defineBy>("Define by");
        int curStepNumber = nodeAnalysisSetting->getPropertyValue<int>("Current step number");

        switch(defineBy)
        {
        case Property::defineBy_components:
        {
            void *CS = node->getPropertyValue<void*>("Coordinate system");

            QStandardItem *itemCS = static_cast<QStandardItem*>(CS);
            SimulationNodeClass *nodeCS = itemCS->data(Qt::UserRole).value<SimulationNodeClass*>();
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
             }
                break;

            case SimulationNodeClass::nodeType_coordinateSystem:
                directionalData = nodeCS->getPropertyValue<QVector<QVector<double>>>("Base directional data");
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
        const QVector<GeometryTag> &locs = node->getPropertyValue<QVector<GeometryTag>>("Geometry");
        //! ---------------------------------------------------------
        //! create the marker only if something is present in "Tags"
        //! ---------------------------------------------------------
        if(locs.isEmpty())
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
        for(int i=0; i<locs.length(); i++)
        {
            int bodyIndex = locs.at(i).parentShapeNr;
            int faceIndex = locs.at(i).subTopNr;
            const TopoDS_Shape &curFace = gDB->MapOfBodyTopologyMap.value(bodyIndex).faceMap.FindKey(faceIndex);
            theBuilder.Add(compound,curFace);
        }
        GProp_GProps prop;
        BRepGProp::SurfaceProperties(compound,prop);
        double Area = prop.Mass();

        //! ----------------------------------------------
        //! this sets the size of the curved arrow marker
        //! ----------------------------------------------
        double Rin = sqrt(Area)/10.0;

        //! -----------------
        //! marker placement
        //! -----------------
        gp_Ax2 axes;
        gp_Pnt P = prop.CentreOfMass();
        axes.SetLocation(P);

        //! ------------------------------------------------------
        //! this for reading the components from the tabular data
        //! ------------------------------------------------------
        SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
        QTreeView *theTree = sm->myTreeView;
        QList<int> tableColumns = mainTreeTools::getColumnsToRead(theTree);
        SimulationNodeClass *nodeAnalysisSetting = sm->getAnalysisSettingsNodeFromCurrentItem();
        CustomTableModel *tabModel = nodeAnalysisSetting->getTabularDataModel();
        Property::defineBy defineBy = node->getPropertyValue<Property::defineBy>("Define by");

        cout<<"____tag02____"<<endl;
        int curStepNumber = nodeAnalysisSetting->getPropertyValue<int>("Current step number");
        cout<<"____tag03____"<<endl;

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
    }
        break;
*/
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

    default:
    {
        return false;
    }
        break;
    }
    cout<<"markerBuilder::addMarker()->____exiting function____"<<endl;
    return true;
}
