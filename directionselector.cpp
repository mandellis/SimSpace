//! custom includes
#include "directionselector.h"
#include "geomtoolsclass.h"
#include "markers.h"

//! Qt
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QGridLayout>
#include <QPushButton>

//! OCC
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Edge.hxx>
#include <gp_Ax1.hxx>
#include <gp_Pnt.hxx>
#include <Graphic3d_ZLayerId.hxx>
#include <Bnd_Box.hxx>
#include <BRepBndLib.hxx>

//! C++
#include <iostream>
using namespace std;

//! ------------------------
//! function: constructor I
//! details:
//! ------------------------
DirectionSelector::DirectionSelector(const occHandle(AIS_InteractiveContext) &aCTX, QWidget *parent):
    ShapeSelectorBox(parent),myCTX(aCTX)
{
    //!cout<<"DirectionSelector::DirectionSelector()->____constructor I called____"<<endl;
    this->createContent();
}

//! -------------------------
//! function: constructor II
//! details:
//! -------------------------
DirectionSelector::DirectionSelector(QWidget *parent):ShapeSelectorBox(parent)
{
    //!cout<<"DirectionSelector::DirectionSelector()->____constructor II called____"<<endl;
    this->createContent();
    this->hidePushButtons();
}

//! ---------------------
//! function: destructor
//! details:
//! ---------------------
DirectionSelector::~DirectionSelector()
{
    if(myCTX->NbSelected())myCTX->ClearSelected(false);
    if(!myMarker.IsNull())myCTX->Remove(myMarker,false);
    myCTX->UpdateCurrentViewer();

}

//! -------------------------
//! function: showPushButton
//! details:
//! -------------------------
void DirectionSelector::showPushButtons()
{
    //!cout<<"DirectionSelector::hidePushButtons()->____function called____"<<endl;
    if(!(bap->isVisible() && bcn->isVisible()))
    {
        bap->show();
        bcn->show();
    }
}

//! -----------------------------------------------------------------//
//! function: hidePushButton                                         //
//! details:                                                         //
//! -----------------------------------------------------------------//
void DirectionSelector::hidePushButtons()
{
    //!cout<<"DirectionSelector::hidePushButtons()->____function called____"<<endl;
    if((bap->isVisible() && bcn->isVisible()))
    {
        bap->hide();
        bcn->hide();
    }
}

//! ----------------------
//! function: setAccepted
//! details:
//! ----------------------
void DirectionSelector::setAccepted()
{
    //!cout<<"DirectionSelector::setAccepted()->____set accepted called____"<<endl;
    hidePushButtons();
    myCTX->ClearSelected(true);

    //! -----------------------------------
    //! remove the marker from the context
    //! -----------------------------------
    if(!myMarker.IsNull()) myCTX->Remove(myMarker,true);
    emit editingFinished();
}

//! ----------------------
//! function: setRejected
//! details:
//! ----------------------
void DirectionSelector::setRejected()
{
    //!cout<<"DirectionSelector::setAccepted()->____set rejected called____"<<endl;
    hidePushButtons();
    myCTX->ClearSelected(true);

    //! -----------------------------------
    //! remove the marker from the context
    //! -----------------------------------
    if(!myMarker.IsNull()) myCTX->Remove(myMarker,true);
    emit editingFinished();
}

//! -----------------------
//! function: getDirection
//! details:
//! -----------------------
QVector<double> DirectionSelector::getDirection()
{
    //!cout<<"DirectionSelector::getDirection()->____function called____"<<endl;
    return myDirection;
}

//! ------------------------------------------------
//! function: setDirection
//! details:  set myShape and activate the selector
//! ------------------------------------------------
void DirectionSelector::setDirection(QVector<double> vec)
{
    //cout<<"DirectionSelector::setDirection()->____function called____"<<endl;

    myDirection.clear();
    for(int i=0;i<6;i++) myDirection.push_back(vec.at(i));

    this->showPushButtons();
    gp_Ax1 theAxis;
    gp_Vec V1;
    gp_Pnt P;
    P.SetX(vec[0]);
    P.SetY(vec[1]);
    P.SetZ(vec[2]);
    theAxis.SetLocation(P);
    for(int i=3; i<6; i++) V1.SetCoord(i-2,vec[i]);
    //! warning: first create a vector for the direction,
    //! then normalize it
    gp_Dir V(V1);
    theAxis.SetDirection(V);

    myMarker = markers::buildArrowMarker(P, theAxis.Direction(), 5);
    myMarker->SetColor(DIRECTION_SELECTOR_COLOR);
    myCTX->SetZLayer(myMarker,Graphic3d_ZLayerId_Top);
    myCTX->Display(myMarker,1,-1,false);
    myCTX->UpdateCurrentViewer();
}

//! ------------------------
//! function: createContent
//! details:
//! ------------------------
void DirectionSelector::createContent()
{
    //! hide the blinking cursor
    this->setReadOnly(true);

    h = new QHBoxLayout(this);
    h->setMargin(0);
    h->setContentsMargins(0,0,0,0);
    h->setSpacing(0);

    g = new QGridLayout();
    g->setHorizontalSpacing(0);
    g->setMargin(0);
    g->setContentsMargins(0,0,0,0);

    h->addLayout(g);

    bap = new QPushButton("apply");
    bcn = new QPushButton("cancel");
    g->addWidget(bap,0,0);
    g->addWidget(bcn,0,1);

    /*
    QString msg;
    if(myShape.Extent()==0)
        msg = "Click to change";
    this->setText(msg);
    */

    //! connections
    connect(bap,SIGNAL(clicked()),this,SLOT(setAccepted()));
    connect(bcn,SIGNAL(clicked()),this,SLOT(setRejected()));
}

//! ---------------------
//! function: setContext
//! details:
//! ---------------------
void DirectionSelector::setContext(const opencascade::handle<AIS_InteractiveContext> &aCTX)
{
    myCTX = aCTX;
}

//! ---------------------
//! function: showMarker
//! details:
//! ---------------------
void DirectionSelector::showMarker()
{
    //cout<<"DirectionSelector::showMarker()->____function called____"<<endl;
    if(!myMarker.IsNull()) myCTX->Remove(myMarker,true);

    TopoDS_Shape theShape;

    //! ------------------------------------------
    //! scan the context and go to the last shape
    //! ------------------------------------------
    for(myCTX->InitSelected(); myCTX->MoreSelected(); myCTX->NextSelected())
    {
        theShape = myCTX->SelectedShape();
    }
    if(!theShape.IsNull() && theShape.ShapeType()!=TopAbs_SOLID && theShape.ShapeType()!=TopAbs_VERTEX)
    {
        //! ----------------------------------
        //! get the bounding box of the shape
        //! ----------------------------------
        Bnd_Box boundingBox;
        BRepBndLib::Add(theShape, boundingBox);
        double D = sqrt(boundingBox.SquareExtent());
        double size = D/30;

        gp_Ax1 theAxis;
        gp_Pnt P;
        //! ------------------------------------------------------------
        //! generate the direction when the selected topology is a face
        //! ------------------------------------------------------------
        if(theShape.ShapeType()==TopAbs_FACE)
        {
            GeomAbs_SurfaceType type;
            TopoDS_Face theCurFace = TopoDS::Face(theShape);
            GeomToolsClass::getFaceType(theCurFace,type);
            if(type == GeomAbs_Plane)
            {
                GeomToolsClass::getPlanarFaceInfo(theCurFace, theAxis, P);
                myMarker = markers::buildArrowMarker(P, theAxis.Direction(), size);
            }
            else if(type == GeomAbs_Cylinder)
            {
                GeomToolsClass::getCylindricalFaceInfo(theCurFace, theAxis, P);
                myMarker = markers::buildArrowMarker(P, theAxis.Direction(), size);
            }
            else if(type == GeomAbs_Sphere)
            {
                GeomToolsClass::getSphericalFaceInfo(theCurFace, theAxis, P);
                myMarker = markers::buildArrowMarker(P, theAxis.Direction(), size);
            }
            else if(type== GeomAbs_Cone)
            {
                GeomToolsClass::getConicalFaceInfo(theCurFace, theAxis, P);
                myMarker = markers::buildArrowMarker(P, theAxis.Direction(), size);
            }
        }
        //! -------------------------------------------------------------
        //! generate the direction when the selected topology is an edge
        //! -------------------------------------------------------------
        else if(theShape.ShapeType()==TopAbs_EDGE)
        {
            TopoDS_Edge theEdge = TopoDS::Edge(theShape);
            GeomAbs_CurveType type;
            GeomToolsClass::getCurveType(theEdge,type);
            if(type == GeomAbs_Line)
            {
                //!cout<<"straight edge selected"<<endl;
                GeomToolsClass::getEdgeInfo(theEdge, theAxis, P);
                myMarker = markers::buildArrowMarker(P, theAxis.Direction(), 5);
            }
        }
        if(!myMarker.IsNull())
        {
            myMarker->SetColor(DIRECTION_SELECTOR_COLOR);
            myCTX->SetZLayer(myMarker,Graphic3d_ZLayerId_Top);
            myCTX->Display(myMarker,1,-1,true);
        }

        //myDirection.clear();
        myDirection.replace(0,P.X());
        myDirection.replace(1,P.Y());
        myDirection.replace(2,P.Z());
        myDirection.replace(3,theAxis.Direction().X());
        myDirection.replace(4,theAxis.Direction().Y());
        myDirection.replace(5,theAxis.Direction().Z());
    }
}
