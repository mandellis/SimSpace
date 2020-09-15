//! ----------------
//! custom includes
//! ----------------
#include "meshselector.h"

//! ----
//! OCC
//! ----
#include <MeshVS_MeshEntityOwner.hxx>
#include <MeshVS_Mesh.hxx>
#include <MeshVS_DataSource.hxx>

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
MeshSelector::MeshSelector(const occHandle(AIS_InteractiveContext) &aCTX, QWidget *parent):
    ShapeSelectorBox(parent), myCTX(aCTX)
{
    cout<<"ShapeSelector::ShapeSelector()->____constructor I called____"<<endl;
    this->createContent();
}

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
MeshSelector::MeshSelector(QWidget *parent):ShapeSelectorBox(parent)
{
    myElements.clear();
    myIDs.clear();

    this->createContent();
    this->hidePushButtons();
}

//! ---------------------
//! function: destructor
//! details:
//! ---------------------
MeshSelector::~MeshSelector()
{
    cout<<"MeshSelector::~MeshSelector()->____destructor called____"<<endl;
    this->clearContext();

    // for checking what remains in selection after destroying object
    cout<<"MeshSelector::~MeshSelector()->____Number of selected elements: "<<myCTX->NbSelected()<<"____"<<endl;
}

//! -------------------------
//! function: showPushButton
//! details:
//! -------------------------
void MeshSelector::showPushButtons()
{
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
void MeshSelector::hidePushButtons()
{
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
void MeshSelector::setAccepted()
{
    cout<<"MeshSelector::setAccepted()->____function called____"<<endl;
    hidePushButtons();

    myElements.clear();
    myIDs.clear();

    //! -------------------------------
    //! retrieve the selected elements
    //! -------------------------------
    for(myCTX->InitSelected();myCTX->MoreSelected();myCTX->NextSelected())
    {
        const occHandle(MeshVS_MeshEntityOwner) &anEntityOwner = occHandle(MeshVS_MeshEntityOwner)::DownCast(myCTX->DetectedOwner());
        const occHandle(MeshVS_Mesh) &aMeshObject = occHandle(MeshVS_Mesh)::DownCast(anEntityOwner->Selectable());
        const occHandle(MeshVS_DataSource) &aMeshDS = aMeshObject->GetDataSource();
        int ID = anEntityOwner->ID();

        //! ---------------
        //! a mesh element
        //! ---------------
        meshElementByCoords aMeshElement;
        aMeshElement.ID = ID;

        int NbNodes, buf[20];
        TColStd_Array1OfInteger nodeIDs(*buf,1,20);
        aMeshDS->GetNodesByElement(ID,nodeIDs,NbNodes);
        for(int n=1; n<=NbNodes; n++)
        {
            int globalNodeID = nodeIDs(n);
            MeshVS_EntityType aType;
            int NbNodes1;
            double bufd[3];
            TColStd_Array1OfReal coords(*bufd,1,3);
            aMeshDS->GetGeom(globalNodeID,false,coords,NbNodes1,aType);
            aMeshElement.pointList<<mesh::meshPoint(coords(1),coords(2),coords(3),globalNodeID);
        }
        MeshVS_EntityType Type;
        bool isElement = true; // <= this constraint the type of elements that can be selected
        aMeshDS->GetGeomType(ID,isElement,Type);
        switch(Type)
        {
        case MeshVS_ET_Face:
        {
            switch(NbNodes)
            {
            case 3: aMeshElement.type = TRIG; break;
            case 4: aMeshElement.type = QUAD; break;
            case 5: aMeshElement.type = PENTA; break;
            case 6: aMeshElement.type = TRIG6; break;
            case 7: aMeshElement.type = EPTA; break;
            case 8: aMeshElement.type = QUAD8; break;
            }
        }
            break;
        case MeshVS_ET_Volume:
        {
            switch(NbNodes)
            {
            case 4: aMeshElement.type = TET; break;
            case 5: aMeshElement.type = PYRAM; break;
            case 6: aMeshElement.type = PRISM; break;
            case 8: aMeshElement.type = HEXA; break;
            case 10: aMeshElement.type = TET10; break;
            case 13: aMeshElement.type = PYRAM13; break;
            case 15: aMeshElement.type = PRISM15; break;
            case 20: aMeshElement.type = HEXA20; break;
            }
        }
            break;
        case MeshVS_ET_Link:
        {
            // not supported for the moment
        }
            break;
        }
        myElements.push_back(aMeshElement);
        myIDs.push_back(ID);
    }

    cout<<"MeshSelector::setAccepted()->____number of current selected elements: "<<myElements.size()<<"____"<<endl;

    this->clearContext();
    emit editingSelectionFinished(); //old position
}

//! ----------------------
//! function: setRejected
//! details:
//! ----------------------
void MeshSelector::setRejected()
{
    cout<<"MeshSelector::setRejected()->____set rejected called____"<<endl;
    this->hidePushButtons();
    this->clearContext();
    emit editingSelectionFinished();
}

//! ----------------------
//! function: getElements
//! details:
//! ----------------------
std::vector<meshElementByCoords> MeshSelector::getElements() const
{
    return myElements;
}

//! ----------------------
//! function: setElements
//! details:
//! ----------------------
void MeshSelector::setElements(const std::vector<meshElementByCoords> &vecElements)
{

}

//! ------------------------
//! function: createContent
//! details:
//! ------------------------
void MeshSelector::createContent()
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

    this->setText(QString("%1 entity selected").arg(myElements.size()));

    //! connections logic
    //disconnect(bap,SIGNAL(clicked()),this,SLOT(setAccepted()));
    //disconnect(bcn,SIGNAL(clicked()),this,SLOT(setRejected()));
    connect(bap,SIGNAL(clicked()),this,SLOT(setAccepted()));
    connect(bcn,SIGNAL(clicked()),this,SLOT(setRejected()));
}

//! -----------------------------------------------
//! function: clearContext
//! details:  clear the context from the selection
//! -----------------------------------------------
void MeshSelector::clearContext(bool updateViewer)
{
    myCTX->ClearSelected(updateViewer);
}
