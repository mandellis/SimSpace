//! ----------------
//! custom includes
//! ----------------
#include "shapeselector.h"
#include "src/utils/topologytools.h"
#include "src/utils/tools.h"
#include <geometrydatabase.h>
#include "src/main/simulationmanager.h"
#include "detailviewer.h"
#include "src/utils/ccout.h"

//! ----
//! OCC
//! ----
#include <SelectMgr_IndexedMapOfOwner.hxx>
#include <TopTools_ListOfShape.hxx>
#include <SelectMgr_EntityOwner.hxx>
#include <StdSelect_BRepOwner.hxx>

//! ---
//! Qt
//! ---
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QGridLayout>
#include <QPushButton>

//! ----
//! C++
//! ----
#include <iostream>
using namespace std;

//! ------------------------
//! function: constructor I
//! details:
//! ------------------------
ShapeSelector::ShapeSelector(const occHandle(AIS_InteractiveContext) &aCTX, QWidget *parent):
    ShapeSelectorBox(parent), myCTX(aCTX)
{
    cout<<"ShapeSelector::ShapeSelector()->____constructor I called____"<<endl;
    if(myCTX.IsNull()) { cout<<"ShapeSelector::ShapeSelector()->____THE CONTEXT IS NOT INIZIALIZED____"<<endl; }
    this->createContent();
}

//! -------------------------
//! function: constructor II
//! details:
//! -------------------------
ShapeSelector::ShapeSelector(QWidget *parent):ShapeSelectorBox(parent)
{
    cout<<"ShapeSelector::ShapeSelector()->____constructor II called____"<<endl;

    //! shape content
    myShapes.Clear();
    myVecLoc.clear();

    this->createContent();
    this->hidePushButtons();
}

//! ---------------------
//! function: destructor
//! details:
//! ---------------------
ShapeSelector::~ShapeSelector()
{
    cout<<"ShapeSelector::~ShapeSelector()->____destructor called____"<<endl;
    this->clearContext();
    cout<<"ShapeSelector::~ShapeSelector()->____Number of selected shapes (must be zero): "<<myCTX->NbSelected()<<"____"<<endl;
}

//! -------------------------
//! function: showPushButton
//! details:
//! -------------------------
void ShapeSelector::showPushButtons()
{
    //cout<<"ShapeSelector::showPushButtons()->____function called____"<<endl;
    if(!(bap->isVisible() && bcn->isVisible()))
    {
        bap->show();
        bcn->show();
    }
}

//! -------------------------
//! function: hidePushButton
//! details:
//! -------------------------
void ShapeSelector::hidePushButtons()
{
    //cout<<"ShapeSelector::hidePushButtons()->____function called____"<<endl;
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
void ShapeSelector::setAccepted()
{
    cout<<"ShapeSelector::setAccepted()->____function called____"<<endl;
    this->hidePushButtons();

    myShapes.Clear();
    myVecLoc.clear();
    for(myCTX->InitSelected();myCTX->MoreSelected();myCTX->NextSelected())
    {
        // DEPRECATED METHOD
        //const TopoDS_Shape &curSelectedShape = myCTX->SelectedShape();

        const occHandle(StdSelect_BRepOwner) &sb  = occHandle(StdSelect_BRepOwner)::DownCast(myCTX->SelectedOwner());
        if(sb.IsNull())
        {
            cout<<"ShapeSelector::setAccepted()->____NULL BRepOwner - case handled____"<<endl;
            continue;
        }
        const TopoDS_Shape &ownerShape = sb->Shape();
        TopoDS_Shape selectedShape = ownerShape.Located(sb->Location() * ownerShape.Location());
        myShapes.Append(selectedShape);
    }

    myVecLoc = TopologyTools::generateLocationPairs(gdb,myShapes);

    //for(int i=0; i<myVecLoc.size(); i++) cout<<"____("<<myVecLoc.at(i).parentShapeNr<<", "<<myVecLoc.at(i).subTopNr<<")____"<<endl;
    //cout<<"ShapeSelector::setAccepted()->____number of current selected shapes: "<<myShapes.Extent()<<"____"<<endl;
    //cout<<"ShapeSelector::setAccepted()->____number of shapes in vecLoc: "<<myVecLoc.size()<<"____"<<endl;

    //emit editingSelectionFinished();    //new position
    this->clearContext();
    emit editingSelectionFinished(); //old position

    //cout<<"ShapeSelector::setAccepted()->____exiting____"<<endl;
}

//! ----------------------
//! function: setRejected
//! details:
//! ----------------------
void ShapeSelector::setRejected()
{
    cout<<"ShapeSelector::setRejected()->____set rejected called____"<<endl;
    myVecLoc = TopologyTools::generateLocationPairs(gdb,myShapes);
    cout<<"ShapeSelector::setAccepted()->____number of shapes in vecLoc: "<<myVecLoc.size()<<"____"<<endl;
    this->hidePushButtons();
    this->clearContext();
    emit editingSelectionFinished();
}

//! -------------------
//! function: getShape
//! details:
//! -------------------
ListOfShape ShapeSelector::getShape() const
{
    return myShapes;
}

//! --------------------
//! function: getVecLoc
//! details:
//! --------------------
std::vector<GeometryTag> ShapeSelector::getVecLoc() const
{
    return myVecLoc;
}

//! -------------------
//! function: setShape
//! details:
//! -------------------
void ShapeSelector::setShape(const std::vector<GeometryTag> &vecLoc)
{
    cout<<"ShapeSelector::setShape()->____set shape called____"<<endl;

    myVecLoc.clear();
    myVecLoc = vecLoc;

    myShapes.Clear();

    //! -----------------------------------------------
    //! rebuild the shape list from the vector of tags
    //! -----------------------------------------------
    for(int i=0; i<vecLoc.size(); i++)
    {
        int parentShapeNr = vecLoc[i].parentShapeNr;
        int subTopNr = vecLoc[i].subTopNr;
        bool isParent = vecLoc[i].isParent;
        TopAbs_ShapeEnum type = vecLoc[i].subShapeType;
        if(isParent) myShapes.Append(gdb->bodyMap.value(parentShapeNr));
        else
        {
            switch(type)
            {
            case TopAbs_FACE: myShapes.Append(gdb->MapOfBodyTopologyMap.value(parentShapeNr).faceMap.FindKey(subTopNr)); break;
            case TopAbs_EDGE: myShapes.Append(gdb->MapOfBodyTopologyMap.value(parentShapeNr).edgeMap.FindKey(subTopNr)); break;
            case TopAbs_VERTEX: myShapes.Append(gdb->MapOfBodyTopologyMap.value(parentShapeNr).vertexMap.FindKey(subTopNr)); break;
            }
        }
    }

    //! -----------------------------------------------------------------------------------
    //! trying to remove the deprecated method AIS_InteractiveContext::AddOrRemoveSelected
    //! -----------------------------------------------------------------------------------
    AIS_ListOfInteractive listIO;
    myCTX->ObjectsInside(listIO,AIS_KOI_Shape,0);
    for(AIS_ListIteratorOfListOfInteractive it(listIO); it.More(); it.Next())
    {
        const occHandle(AIS_Shape) &curAISShape = occHandle(AIS_Shape)::DownCast(it.Value());

        //! ----------------------------------------
        //! all the owners of the current AIS_Shape
        //! ----------------------------------------
        occHandle(SelectMgr_IndexedMapOfOwner) m;
        myCTX->EntityOwners(m,curAISShape,-1);

        //! ----------------------------
        //! iterate over all the owners
        //! ----------------------------
        for(int n=1; n<=m->Extent(); n++)
        {
            occHandle(SelectMgr_EntityOwner) s = m->FindKey(n);
            if(s.IsNull())
            {
                cout<<"ShapeSelector::setShape()->____NULL owner - handled case____"<<endl;
                continue;
            }
            occHandle(StdSelect_BRepOwner) sb  = occHandle(StdSelect_BRepOwner)::DownCast(s);
            if(sb.IsNull())
            {
                cout<<"ShapeSelector::setShape()->____NULL BRep owner - handled case____"<<endl;
                continue;
            }
            const TopoDS_Shape &ownerShape = sb->Shape();
            TopoDS_Shape locatedShape = ownerShape.Located(sb->Location() * ownerShape.Location());
            if(myShapes.Contains(locatedShape)==false) continue;

            myCTX->AddOrRemoveSelected(sb,false);
        }
    }
    myCTX->UpdateCurrentViewer();
    this->showPushButtons();

    cout<<"ShapeSelector::setShape()->____myShapes size: "<<myShapes.Extent()<<"____"<<endl;
    cout<<"ShapeSelector::setShape()->____myVecLoc size: "<<myVecLoc.size()<<"____"<<endl;
    cout<<"ShapeSelector::setShape()->____exiting function____"<<endl;
}

//! ------------------------
//! function: createContent
//! details:
//! ------------------------
void ShapeSelector::createContent()
{
    //! ------------------------------------------------------------
    //! link to the simulationManager and get the geometry database
    //! ------------------------------------------------------------
    SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
    gdb = sm->getDataBase();

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

    QString msg;
    if(myShapes.Extent()!=0)
    {
        msg = QString("%1").arg(myShapes.Extent());
        switch(myShapes.First().ShapeType())
        {
        case TopAbs_VERTEX: msg.append(" vertexes"); break;
        case TopAbs_EDGE: msg.append(" edges"); break;
        case TopAbs_FACE: msg.append(" faces"); break;
        case TopAbs_SOLID: msg.append(" solids"); break;
        default: break;
        }
    }
    else
    {
        msg = "No selection";
    }
    this->setText(msg);

    //! connections logic
    //disconnect(bap,SIGNAL(clicked()),this,SLOT(setAccepted()));
    //disconnect(bcn,SIGNAL(clicked()),this,SLOT(setRejected()));
    connect(bap,SIGNAL(clicked()),this,SLOT(setAccepted()));
    connect(bcn,SIGNAL(clicked()),this,SLOT(setRejected()));
}

//! ---------------------
//! function: setContext
//! details:
//! ---------------------
void ShapeSelector::setContext(const occHandle(AIS_InteractiveContext) &aCTX)
{
    myCTX = aCTX;
}

//! -----------------------
//! function: clearContext
//! details:
//! -----------------------
void ShapeSelector::clearContext(bool updateViewer)
{
    //if(myCTX->NbSelected()>0) myCTX->ClearSelected(updateViewer);
    myCTX->ClearSelected(updateViewer);
}
