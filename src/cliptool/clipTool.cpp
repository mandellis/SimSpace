//! ----------------
//! custom includes
//! ----------------
#include "clipTool.h"
#include "cliptooldelegate.h"
#include "simulationnodeclass.h"
#include "qextendedstandarditem.h"
#include <meshslicer.h>

//! ---
//! Qt
//! ---
#include <QAction>
#include <QMenu>
#include <QStyledItemDelegate>
#include <QMessageBox>
#include <QHeaderView>
#include <QScrollBar>

//! ----
//! C++
//! ----
#include <time.h>
#include <random>

//! ----
//! OCC
//! ----
#include <MeshVS_MeshPrsBuilder.hxx>
#include <MeshVS_Drawer.hxx>
#include <MeshVS_DrawerAttribute.hxx>
#include <AIS_Plane.hxx>
#include <TColStd_HPackedMapOfInteger.hxx>
#include <TColStd_MapIteratorOfPackedMapOfInteger.hxx>
#include <TColStd_PackedMapOfInteger.hxx>

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
clipTool::clipTool(QWidget *parent):QTableView(parent),
    myCurPoint(0,0),
    myMaxClipPlanes(0),
    myCurNumberOfClipPlanes(0),
    myWorkingMode(2)
{
    //! ---------------------------------
    //! init the random number generator
    //! ---------------------------------
    srand(time(NULL));

    //! ----------------
    //! set object name
    //! ----------------
    this->setObjectName("clipTool");

    //! ----------------------
    //! enable mouse tracking
    //! ----------------------
    this->setMouseTracking(true);

    //! --------------------------
    //! create the internal model
    //! --------------------------
    internalModel = new QStandardItemModel(this);
    this->setModel(internalModel);

    //! -----------------
    //! set the delegate
    //! -----------------
    clipToolDelegate *theDelegate = new clipToolDelegate(this);
    this->setItemDelegate(theDelegate);

    //! -------------------
    //! enable custom menu
    //! -------------------
    this->setContextMenuPolicy(Qt::CustomContextMenu);
    connect(this,SIGNAL(customContextMenuRequested(QPoint)),this,SLOT(showContextMenu(QPoint)));

    //! ----------------------------
    //! hide the horizontal headers
    //! ----------------------------
    this->horizontalHeader()->hide();

    //! -----------------------------------------------------
    //! automatic update of the plane data in the table cell
    //! -----------------------------------------------------
    connect(theDelegate,SIGNAL(currentCSChanged()),this,SLOT(updateClipPlane()));

    //! --------------------------------
    //! update of the clip plane status
    //! --------------------------------
    //connect(theDelegate,SIGNAL(currentCSStatusChanged()),this,SLOT(updateClipPlane()));

    //! --------------------------------------------
    //! dynamic update after clip plane translation
    //! --------------------------------------------
    connect(theDelegate,SIGNAL(currentCSTranslationApplied(int)),this,SLOT(updateCSTranslation(int)));

    //! ------------------------
    //! the table fits the view
    //! ------------------------
    this->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
}

//! --------------------------
//! function: setMeshDataBase
//! details:
//! --------------------------
void clipTool::setMeshDataBase(meshDataBase *aMeshDataBase)
{
    myMDB = aMeshDataBase;
}

//! -----------------------------------------
//! function: isCurrentRowActive
//! details:  get the flag on/off of a row
//!           (status of the clipping plane)
//! -----------------------------------------
bool clipTool::isCurrentRowActive()
{
    QModelIndex index = this->currentIndex();
    if(index.isValid())
    {
        QVariant data = internalModel->item(index.row(), CLIPPLANE_STATUS_COLUMN)->data(Qt::UserRole);
        return data.toBool();
    }
    return false;
}

//! ----------------------------------------------------------
//! function: setWorkingMode
//! details:  0 - mesh; 1 - contacts; 2 - model; 3 - solution
//! ----------------------------------------------------------
void clipTool::setWorkingMode(int workingMode)
{
    cout<<"clipTool::setWorkingMode()->____current working mode: "<<workingMode<<"____"<<endl;

    myWorkingMode = workingMode;

    //! --------------------------------
    //! prevent null pointer exceptions
    //! --------------------------------
    if(myOCCViewer==Q_NULLPTR) return;
    if(myMDB==Q_NULLPTR) return;
    if(internalModel==Q_NULLPTR) return;

    //! ----------------------------------------
    //! retrieve the the active clipping planes
    //! ----------------------------------------
    std::vector<int> activeClipPlanes;
    for(int row=0; row<internalModel->rowCount(); row++)
    {
        bool isActive = internalModel->index(row,CLIPPLANE_STATUS_COLUMN).data(Qt::UserRole).toBool();
        int planeID = internalModel->index(row,CLIPPLANE_ID_COLUMN).data(Qt::UserRole).toInt();
        if(!isActive) continue;
        activeClipPlanes.push_back(planeID);
    }

    switch(myWorkingMode)
    {
    case 0: // on mesh
    {
        //! ----------------------------------------------------------
        //! recompute the hidden elements and send them to the viewer
        //! ----------------------------------------------------------
        this->computeHiddenElements();
        myOCCViewer->setHiddenElements(myHiddenElements);
    }
        break;

    default:    // other working modes
    {
        if(!myOCCViewer->getMeshContext().IsNull()) myOCCViewer->getMeshContext()->RemoveAll(true);
    }
        break;
    }
}

//! --------------------------
//! function: showContextMenu
//! details:
//! --------------------------
void clipTool::showContextMenu(QPoint aPoint)
{
    //! ---------------------------------------------------------------
    //! create the context menu only if coordinate systems are present
    //! ---------------------------------------------------------------
    if(myCoordinateSystemRoot!=NULL)
    {
        QPoint pos = this->mapToGlobal(aPoint);
        QMenu *ctxMenu= new QMenu(this);

        QAction *actionAdd = ctxMenu->addAction("Add");
        actionAdd->setIcon(QIcon(":/icons/icon_add.png"));
        actionAdd->setData(0);

        //! ----------------------------------------------------
        //! add the option "Remove" only when items are present
        //! ----------------------------------------------------
        if(internalModel->rowCount()>0)
        {
            QAction *actionRemove = ctxMenu->addAction("Remove");
            actionRemove->setIcon(QIcon(":/icons/icon_delete.png"));
            actionRemove->setData(1);
        }

        QAction* selectedItem = ctxMenu->exec(pos);

        if(selectedItem)
        {
            switch(selectedItem->data().toInt())
            {
            case 0:
            {
                this->addItemToTable();
                emit clipPlaneAdded();
            }
                break;

            case 1:
            {
                this->removeItemFromTable();
                emit clipPlaneRemoved();
            }
                break;
            }
        }
    }
}

//! -------------------------
//! function: addItemToTable
//! details:
//! -------------------------
void clipTool::addItemToTable()
{
    if(myCoordinateSystemRoot!=Q_NULLPTR && myOCCViewer!=Q_NULLPTR)
    {
        if(myCurNumberOfClipPlanes>myMaxClipPlanes)
        {
            QMessageBox::critical(this,APPNAME,"Maximum number of clippig planes reached",QMessageBox::Ok);
            return;
        }

        //! ----------------------------------
        //! generate the ID of the clip plane
        //! ----------------------------------
        int clipPlaneID = rand();

        //! ------------------------------
        //! set the current clip plane ID
        //! ------------------------------
        myOCCViewer->setCurrentClipPlane(clipPlaneID);

        cout<<"clipTool::addItemToTable()->____set current plane ID: "<<clipPlaneID<<"____"<<endl;
        cout<<"clipTool::addItemToTable()->____function called: clip plane ID: "<<clipPlaneID<<"____"<<endl;

        //! -------
        //! create
        //! -------
        QList<QStandardItem*> itemList;
        QStandardItem *item;
        QVariant data;

        //! --------------
        //! name: col = 0
        //! --------------
        item = new QStandardItem();
        QString clipPlaneName = QString("Clip plane %1").arg(clipPlaneID);
        data.setValue(clipPlaneName);
        item->setData(data,Qt::UserRole);
        item->setData(data,Qt::DisplayRole);
        itemList<<item;

        //! ------------
        //! ID: col = 1
        //! ------------
        item = new QStandardItem();
        data.setValue(clipPlaneID);
        item->setData(clipPlaneID,Qt::UserRole);
        data.setValue(QString("%1").arg(clipPlaneID));
        item->setData(data,Qt::DisplayRole);
        itemList<<item;

        //! ------------------------------------
        //! status: "0" => Off => "1" On: col 2
        //! ------------------------------------
        item = new QStandardItem();
        data.setValue(true);
        item->setData(data,Qt::UserRole);
        data.setValue(QString("On"));
        item->setData(data,Qt::DisplayRole);
        itemList<<item;

        //! --------------------------------
        //! "base" coordinate system: col 3
        //! --------------------------------
        item = new QStandardItem();
        QExtendedStandardItem *itemCSGlobal = static_cast<QExtendedStandardItem*>(myCoordinateSystemRoot->child(0,0));
        void *p = (void*)itemCSGlobal;
        data.setValue(p);
        item->setData(data,Qt::UserRole);
        data.setValue(itemCSGlobal->data(Qt::UserRole).value<SimulationNodeClass*>()->getName());
        item->setData(data,Qt::DisplayRole);
        itemList<<item;

        //! -------------------------------------------------------------------------------
        //! create an adjacent cell with the data of the (global) coordinate system: col 4
        //! -------------------------------------------------------------------------------
        QVector<double> coeff = this->getPlaneCoefficients(static_cast<QExtendedStandardItem*>(myCoordinateSystemRoot->child(0,0)));
        data.setValue(coeff);

        item = new QStandardItem();
        item->setData(data,Qt::UserRole);
        double A,B,C,D;
        A = coeff.at(0); B = coeff.at(1); C = coeff.at(2); D = coeff.at(3);
        data.setValue(QString("(%1, %2, %3, %4)").arg(A).arg(B).arg(C).arg(D));
        item->setData(data,Qt::DisplayRole);
        itemList<<item;

        //! ---------------------
        //! "Translation": col 5
        //! ---------------------
        item = new QStandardItem();
        double deltaZ = 0.0;
        data.setValue(deltaZ);
        item->setData(data,Qt::UserRole);
        data.setValue(QString("%1").arg(deltaZ));
        item->setData(data,Qt::DisplayRole);
        itemList<<item;

        //! -----------------------------------
        //! append to to the table the new row
        //! -----------------------------------
        internalModel->appendRow(itemList);

        //! ---------------------------------
        //! add the clip plane to the viewer
        //! ---------------------------------
        myOCCViewer->addClipPlane(A,B,C,D,clipPlaneID,true);

        //! --------------------------------------------------
        //! recompute the hidden elements
        //! send to the viewer the map of element IDs to hide
        //! --------------------------------------------------
        this->computeHiddenElements();
        myOCCViewer->setHiddenElements(myHiddenElements);

        //! ------------------------------
        //! labels for horizontal headers
        //! ------------------------------
        QList<QString> horizontalHeaderLabels;
        horizontalHeaderLabels<<"Clip plane"<<"ID"<<"Status"<<"Coordinate system"<<"Plane coeff"<<"Z position";
        internalModel->setHorizontalHeaderLabels(horizontalHeaderLabels);

        this->verticalHeader()->hide();
        this->horizontalHeader()->setVisible(true);
        this->verticalScrollBar()->setVisible(true);
        this->setAlternatingRowColors(true);

        //! ------------------
        //! hide some columns
        //! ------------------
        //this->setColumnHidden(CLIPPLANE_ID_COLUMN,true);
        //this->setColumnHidden(CLIPPLANE_BASE_PLANE_DATA_COLUMN,true);

        myCurNumberOfClipPlanes++;
    }
    //! --------------------------------
    //! disable all the selection modes
    //! --------------------------------
    this->setSelectionMode(QAbstractItemView::SingleSelection);
    //this->setSelectionMode(QAbstractItemView::NoSelection);

    this->resizeColumnsToContents();
    this->resizeRowsToContents();

    //this->horizontalHeader()->setStretchLastSection(true);
}

//! -----------------------------
//! function: retrieveClipPlanes
//! details:
//! -----------------------------
int clipTool::retrieveActiveClipPlanes(std::map<int,std::vector<double>> &mapOfClipPlanes)
{
    double lx,ly,lz;
    myOCCViewer->getSceneBoundingBox(lx,ly,lz);
    double BBXdiagonal = sqrt(lx*lx+ly*ly+lz*lz);

    int NbClipPlanes = internalModel->rowCount();
    if(NbClipPlanes==0) return 0;

    for(int row = 0; row<NbClipPlanes; row++)
    {
        QModelIndex index1 = internalModel->index(row,CLIPPLANE_STATUS_COLUMN);
        bool isActive = index1.data(Qt::UserRole).toBool();
        if(isActive==false) continue;

        QModelIndex index2 = internalModel->index(row,CLIPPLANE_BASE_PLANE_DATA_COLUMN);
        const QVector<double> &coeffs = index2.data(Qt::UserRole).value<QVector<double>>();

        QModelIndex index3 = internalModel->index(row,CLIPPLANE_TRANSLATION_COLUMN);
        int sliderValue = index3.data(Qt::UserRole).toInt();

        QModelIndex index4 = internalModel->index(row,CLIPPLANE_ID_COLUMN);
        int planeID = index4.data(Qt::UserRole).toInt();

        double k = double(sliderValue/100.0)*(BBXdiagonal/1.0);
        std::vector<double> aPlaneCoeffs {coeffs[0],coeffs[1],coeffs[2],coeffs[3]+k};
        mapOfClipPlanes.insert(std::make_pair(planeID,aPlaneCoeffs));
    }
    return (int)mapOfClipPlanes.size();
}

//! ------------------------------
//! function: removeItemFromTable
//! details:  remove a clip plane
//! ------------------------------
void clipTool::removeItemFromTable()
{
    cout<<"TableView::removeItemFromTable()->____function called____"<<endl;
    QModelIndex curClickedIndex = this->indexAt(myCurPoint);
    if(curClickedIndex.isValid())
    {
        int row = curClickedIndex.row();
        int clipPlaneID = internalModel->index(row,CLIPPLANE_ID_COLUMN).data(Qt::UserRole).toInt();
        cout<<"____removing: "<<clipPlaneID<<"____"<<endl;
        internalModel->removeRow(curClickedIndex.row());
        myOCCViewer->removeClipPlane(clipPlaneID);
        myCurNumberOfClipPlanes--;
    }

    //! --------------------------
    //! recompute hidden elements
    //! --------------------------
    this->computeHiddenElements();
    myOCCViewer->setHiddenElements(myHiddenElements);

    //! ------------------------------
    //! set the current clip plane ID
    //! ------------------------------
    //myOCCViewer->setCurrentClipPlane(-1);
}

//! -----------------------------
//! function: updateCSDefinition
//! details:
//! -----------------------------
void clipTool::updateCSDefinition()
{
    cout<<"clipTool::updateCSDefinition()->____function called____"<<endl;
    this->computeHiddenElements();
    myOCCViewer->setHiddenElements(myHiddenElements);
    /*
    QModelIndex index = internalModel->index(this->currentIndex().row(),this->currentIndex().column());
    void *p = index.data(Qt::UserRole).value<void*>();
    QStandardItem *curItemCS = static_cast<QStandardItem*>(p);

    const QVector<double> &coeffs = this->getPlaneCoefficients(curItemCS);
    QVariant data;
    data.setValue(coeffs);
    QModelIndex indexPlaneCoeffs = internalModel->index(this->currentIndex().row(),CLIPPLANE_BASE_PLANE_DATA_COLUMN);
    internalModel->setData(indexPlaneCoeffs,data,Qt::UserRole);
    double A = coeffs[0];
    double B = coeffs[1];
    double C = coeffs[2];
    double D = coeffs[3];
    data.setValue(QString("(%1, %2, %3, %4").arg(A).arg(B).arg(C).arg(D));
    internalModel->setData(indexPlaneCoeffs,data,Qt::DisplayRole);

    cout<<"____(A,B,C,D) = ("<<A<<", "<<B<<", "<<C<<", "<<D<<")____"<<endl;

    QModelIndex indexID = internalModel->index(this->currentIndex().row(),CLIPPLANE_ID_COLUMN);
    int ID = indexID.data(Qt::UserRole).toInt();
    //cout<<"clipTool::updateCSDefinition()->____ID: "<<ID<<"____"<<endl;

    QModelIndex indexZTranslation = internalModel->index(this->currentIndex().row(),CLIPPLANE_TRANSLATION_COLUMN);
    int zVal = indexZTranslation.data(Qt::UserRole).toInt();
    //cout<<"clipTool::updateCSDefinition()->____z val: "<<zVal<<"____"<<endl;

    myOCCViewer->addClipPlane(A,B,C,D,ID,true);
    myOCCViewer->updateClipPlaneTranslation(ID,zVal,coeffs);
    */
}

//! ------------------------------------------------------------
//! function: updateCSDataByExternalCSChange
//! details:  the definition of a coordinate system is changed:
//!           update all the clipPlanes containing that CS
//! ------------------------------------------------------------
void clipTool::updateCSDataByExternalCSChange(QStandardItem *theCurrentModifiedCS)
{
    SimulationNodeClass *nodeCS = theCurrentModifiedCS->data(Qt::UserRole).value<SimulationNodeClass*>();
    cerr<<"clipTool::updateCSDataByExternalCSChange()->____current changed CS: "<<nodeCS->getName().toStdString()<<"____"<<endl;

    //! --------------
    //! scan the rows
    //! --------------
    for(int row = 0; row< this->internalModel->rowCount(); row++)
    {
        QModelIndex index = this->internalModel->index(row,CLIPPLANE_BASE_COORDINATE_SYSTEM_COLUMN);
        void *p = index.data(Qt::UserRole).value<void*>();
        QStandardItem *curItem = static_cast<QStandardItem*>(p);
        if(curItem == theCurrentModifiedCS)
        {
            cerr<<"+++++++++++++++"<<endl;
            this->updateCSDefinition();
        }
    }
}

//! -------------------------------
//! function: getPlaneCoefficients
//! details:
//! -------------------------------
QVector<double> clipTool::getPlaneCoefficients(QStandardItem *aCSItem)
{
    //cout<<"clipTool::getPlaneCoefficients()->____function called____"<<endl;

    //! -----------------------------------------
    //! item from the coordinate system selector
    //! -----------------------------------------
    SimulationNodeClass *curCSNode = aCSItem->data(Qt::UserRole).value<SimulationNodeClass*>();
    double X_origin, Y_origin, Z_origin;
    QVector<double> XaxisData, YaxisData, ZaxisData;

    switch(curCSNode->getType())
    {
    case SimulationNodeClass::nodeType_coordinateSystem_global:
    {
        cout<<"clipTool::getPlaneCoefficients()->____global coordinate system selected____"<<endl;
        X_origin = curCSNode->getPropertyValue<double>("Origin X");
        Y_origin = curCSNode->getPropertyValue<double>("Origin Y");
        Z_origin = curCSNode->getPropertyValue<double>("Origin Z");
        XaxisData = curCSNode->getPropertyValue<QVector<double>>("X axis data");
        YaxisData = curCSNode->getPropertyValue<QVector<double>>("Y axis data");
        ZaxisData = curCSNode->getPropertyValue<QVector<double>>("Z axis data");
    }
        break;

    case SimulationNodeClass::nodeType_coordinateSystem:
    {
        cout<<"clipTool::getPlaneCoefficients()->____coordinate system selected____"<<endl;
        QVector<double> baseOrigin = curCSNode->getPropertyValue<QVector<double>>("Base origin");
        X_origin = baseOrigin[0];
        Y_origin = baseOrigin[1];
        Z_origin = baseOrigin[2];
        QVector<QVector<double>> axisData = curCSNode->getPropertyValue<QVector<QVector<double>>>("Base directional data");
        XaxisData = axisData[0];
        YaxisData = axisData[1];
        ZaxisData = axisData[2];
    }
        break;
    }

    double axx = XaxisData[0]; double axy = XaxisData[1]; double axz = XaxisData[2];
    //double ayx = YaxisData[0]; double ayy = YaxisData[1]; double ayz = YaxisData[2];
    double azx = ZaxisData[0]; double azy = ZaxisData[1]; double azz = ZaxisData[2];

    //! -------------------------------------------------------
    //! definition of N - "Main direction" normal to the plane
    //! -------------------------------------------------------
    double xN = X_origin + azx;
    double yN = Y_origin + azy;
    double zN = Z_origin + azz;
    gp_Vec vec_N (gp_Pnt(X_origin,Y_origin,Z_origin),gp_Pnt(xN,yN,zN));

    //! --------------------------------------
    //! definition of X - X axis of the plane
    //! --------------------------------------
    double x1 = X_origin + axx;
    double y1 = X_origin + axy;
    double z1 = X_origin + axz;
    gp_Vec vec_X (gp_Pnt(X_origin,X_origin,X_origin),gp_Pnt(x1,y1,z1));

    gp_Ax3 localCSSys(gp_Pnt(X_origin,X_origin,X_origin), vec_N, vec_X);
    gp_Pln plane(localCSSys);

    //! ---------------------------
    //! get the plane coefficients
    //! ---------------------------
    double A,B,C,D;
    plane.Coefficients(A,B,C,D);
    QVector<double> coeff{A,B,C,D};
    return coeff;
}

//! --------------------------
//! function: updateClipPlane
//! details:
//! --------------------------
void clipTool::updateClipPlane()
{
    cout<<"clipTool::updateClipPlane()->____function called____"<<endl;
    bool isOn = internalModel->index(this->currentIndex().row(),CLIPPLANE_STATUS_COLUMN).data(Qt::UserRole).toBool();
    if(isOn==false) return;

    void* p = this->currentIndex().data(Qt::UserRole).value<void*>();
    QStandardItem *itemCS = (QStandardItem*)p;
    SimulationNodeClass *nodeCS = itemCS->data(Qt::UserRole).value<SimulationNodeClass*>();
    QVector<double> ZAxisData, origin;
    switch(nodeCS->getType())
    {
    case SimulationNodeClass::nodeType_coordinateSystem_global:
        ZAxisData = nodeCS->getPropertyValue<QVector<double>>("Z axis data");
        origin.push_back(0); origin.push_back(0); origin.push_back(0);
        break;
    case SimulationNodeClass::nodeType_coordinateSystem:
        ZAxisData = nodeCS->getPropertyValue<QVector<QVector<double>>>("Base directional data")[2];
        origin = nodeCS->getPropertyValue<QVector<double>>("Base origin");
        break;
    }

    //! ------------------------------
    //! normalized plane coefficients
    //! ------------------------------
    double a = ZAxisData[0]; double b = ZAxisData[1]; double c = ZAxisData[2];
    double d = -(a*origin[0]+b*origin[1]+c*origin[2]);
    double l = sqrt(a*a+b*b+c*c);
    a /= l; b /=l; c /=l; d/=l;

    QModelIndex index = this->currentIndex();   // cell of clip plane selector

    //! -------------------------------------------------------------------
    //! sei andato a scuola, per favore scrivi le equazioni esplicitamente
    //! -------------------------------------------------------------------
    gp_Pln aPlane(a,b,c,d);
    gp_Dir aDir(a,b,c);
    gp_Vec aVec(aDir);
    QModelIndex indexTranslation = internalModel->index(index.row(),CLIPPLANE_STATUS_COLUMN);
    int sliderPosition = indexTranslation.data(Qt::UserRole).toInt();

    double lx,ly,lz;
    myOCCViewer->getSceneBoundingBox(lx,ly,lz);
    double translation = sqrt(lx*lx+ly*ly+lz*lz)*double(sliderPosition/100.0);
    aVec.Scale(translation);

    aPlane.Translate(aVec);
    aPlane.Coefficients(a,b,c,d);
    l = sqrt(a*a+b*b+c*c);
    a /= l; b /=l; c /=l; d/=l;
    QVector<double> coeffs {a,b,c,d};

    //! ----------------
    //! update the cell
    //! ----------------
    QModelIndex index1 = index.sibling(index.row(),CLIPPLANE_BASE_PLANE_DATA_COLUMN);

    QVariant value;
    value.setValue(coeffs);
    internalModel->setData(index1,value,Qt::UserRole);
    value.setValue(QString("(%1 ,%2 ,%3 , %4").arg(a).arg(b).arg(c).arg(d));
    internalModel->setData(index1,value,Qt::DisplayRole);

    int clipPlaneID = internalModel->index(this->currentIndex().row(),CLIPPLANE_ID_COLUMN).data(Qt::UserRole).toInt();
    myOCCViewer->addClipPlane(a,b,c,d,clipPlaneID,true);

    //! --------------------------
    //! recompute hidden elements
    //! --------------------------
    this->computeHiddenElements();
    myOCCViewer->setHiddenElements(myHiddenElements);
}

//! ------------------------------
//! function: updateCSTranslation
//! details:
//! ------------------------------
void clipTool::updateCSTranslation(int sliderPosition)
{
    //! ------------------------------------------
    //! get the current CS undergoing translation
    //! and the plane coefficients
    //! ------------------------------------------
    int ID = internalModel->index(this->currentIndex().row(),CLIPPLANE_ID_COLUMN).data(Qt::UserRole).toInt();
    const QVector<double> &coeffs = internalModel->index(this->currentIndex().row(),CLIPPLANE_BASE_PLANE_DATA_COLUMN).data(Qt::UserRole).
            value<QVector<double>>();

    //! -------------------------------------------------
    //! this moves the plane and redraw the clipped view
    //! -------------------------------------------------
    myOCCViewer->updateClipPlaneTranslation(ID,sliderPosition,coeffs);

    if(myMDB==NULL)
    {
        cout<<"clipTool::setWorkingMode()->____cannot change the cliptool working mode: the mesh data base is null____"<<endl;
        return;
    }

    if(myWorkingMode==0)
    {
        const QMap<int,occHandle(Graphic3d_ClipPlane)> &clipPlanes = myOCCViewer->getClipPlanes();
        for(int row =0; row<internalModel->rowCount(); row++)
        {
            int clipPlaneID = internalModel->index(row,CLIPPLANE_ID_COLUMN).data(Qt::UserRole).toInt();
            bool clipPlaneIsOn = internalModel->index(row,CLIPPLANE_STATUS_COLUMN).data(Qt::UserRole).toBool();
            if(!clipPlaneIsOn) continue;
            const occHandle(Graphic3d_ClipPlane) &curClipPlane = clipPlanes.value(clipPlaneID);
            Graphic3d_ClipPlane::Equation planeEquation = curClipPlane->GetEquation();
            double a = *planeEquation.GetData();
            double b = *(planeEquation.GetData()+1);
            double c = *(planeEquation.GetData()+2);
            double d = *(planeEquation.GetData()+3);

            //! --------------------------------
            //! compute the element IDs to hide
            //! --------------------------------
            const QMap<int,occHandle(AIS_InteractiveObject)> &meshObjects = myOCCViewer->getMeshObjects();
            for(QMap<int,occHandle(AIS_InteractiveObject)>::const_iterator it = meshObjects.cbegin(); it!=meshObjects.cend(); it++)
            {
                int bodyIndex = it.key();
                const occHandle(MeshVS_Mesh) &aMeshObject = occHandle(MeshVS_Mesh)::DownCast(it.value());
                const occHandle(MeshVS_DataSource) &aMeshDS = aMeshObject->GetDataSource();
                if(aMeshDS.IsNull()) continue;

                meshSlicer aMeshSlicer(aMeshDS);
                occHandle(TColStd_HPackedMapOfInteger) hiddenElementIDs;
                aMeshSlicer.perform(a,b,c,d,hiddenElementIDs);
                std::map<int,occHandle(TColStd_HPackedMapOfInteger)>::iterator itt = myHiddenElements.find(bodyIndex);
                if(itt==myHiddenElements.end()) myHiddenElements.insert(std::make_pair(bodyIndex,hiddenElementIDs));
                else itt->second = hiddenElementIDs;
            }
        }
        //! --------------------------------------------------
        //! send to the viewer the map of element IDs to hide
        //! --------------------------------------------------
        myOCCViewer->setHiddenElements(myHiddenElements);
    }
}

//! --------------------------------
//! function: updateClippedMeshView
//! details:
//! --------------------------------
void clipTool::updateClippedMeshView(bool onlyExterior)
{
    cout<<"clipTool::updateClippedMeshView()->____"<<onlyExterior<<"____"<<endl;
}

/*
//! -----------------------
//! function: displayMesh
//! details:  display tool
//! -----------------------
void clipTool::displayMesh(const occHandle(MeshVS_DataSource) &aMeshDS,
                           Quantity_NameOfColor aColorName,
                           bool showMeshEdges)
{
    cout<<"occPreGLWidget::displayMesh()->____function called____"<<endl;

    static int bodyIndex;
    bodyIndex++;

    //! ----------------------------
    //! the mesh interactive object
    //! ----------------------------
    occHandle(MeshVS_Mesh) aMeshIO = new MeshVS_Mesh();
    aMeshIO->SetDataSource(aMeshDS);

    //! ----------------------------------------
    //! create and add the presentation builder
    //! ----------------------------------------
    occHandle(MeshVS_MeshPrsBuilder) aBuilder = new MeshVS_MeshPrsBuilder(aMeshIO);
    aMeshIO->AddBuilder(aBuilder,Standard_False);

    //! --------------------------------------
    //! cosmetic for displaying the face mesh
    //! --------------------------------------
    Graphic3d_MaterialAspect myAspect(Graphic3d_NOM_GOLD);
    myAspect.SetColor(static_cast<Quantity_Color>(aColorName));

    aMeshIO->GetDrawer()->SetMaterial(MeshVS_DA_FrontMaterial, myAspect);
    aMeshIO->GetDrawer()->SetBoolean(MeshVS_DA_DisplayNodes, Standard_False);
    aMeshIO->GetDrawer()->SetBoolean(MeshVS_DA_ShowEdges, showMeshEdges);
    aMeshIO->GetDrawer()->SetColor(MeshVS_DA_EdgeColor,Quantity_NOC_BLACK);
    aMeshIO->SetDisplayMode(MeshVS_DMF_Shading);

    myOCCViewer->getContext()->Display(aMeshIO,Standard_True);
    mySlicedMeshIO.insert(bodyIndex,aMeshIO);
}
*/

//! --------------------------------
//! function: computeHiddenElements
//! details:
//! --------------------------------
void clipTool::computeHiddenElements()
{
    cout<<"clipTool::computeHiddenElements()->____function called____"<<endl;

    //! ------------------------------------------------------
    //! in case no clip plane is active or there are not clip
    //! planes defined, show all elements
    //! ------------------------------------------------------
    std::map<int,std::vector<double>> mapOfClipPlanes;
    int NbActiveClipPlanes = this->retrieveActiveClipPlanes(mapOfClipPlanes);
    if(NbActiveClipPlanes==0)
    {
        cout<<"____no clip plane left: showing all elements and returning____"<<endl;

        TColStd_PackedMapOfInteger emptyMap;
        occHandle(TColStd_HPackedMapOfInteger) emptyHMap = new TColStd_HPackedMapOfInteger;
        emptyHMap->ChangeMap()=emptyMap;
        const QMap<int,occHandle(AIS_InteractiveObject)> &mapOfMeshObjects = myOCCViewer->getMeshObjects();
        for(QMap<int,occHandle(AIS_InteractiveObject)>::const_iterator it = mapOfMeshObjects.cbegin(); it != mapOfMeshObjects.cend(); it++)
        {
            int bodyIndex = it.key();
            std::map<int,occHandle(TColStd_HPackedMapOfInteger)>::iterator itt = myHiddenElements.find(bodyIndex);
            if(itt==myHiddenElements.end()) myHiddenElements.insert(std::make_pair(bodyIndex,emptyHMap));
            else itt->second = emptyHMap;
        }
        return;
    }

    //! ---------------------
    //! mesh slicer instance
    //! ---------------------
    meshSlicer aMeshSlicer;

    //! ------------------------
    //! iterate over the meshes
    //! ------------------------
    const QMap<int,occHandle(AIS_InteractiveObject)> &mapOfMeshObjects = myOCCViewer->getMeshObjects();
    for(QMap<int,occHandle(AIS_InteractiveObject)>::const_iterator it = mapOfMeshObjects.cbegin(); it != mapOfMeshObjects.cend(); it++)
    {
        const occHandle(MeshVS_Mesh) &aMeshObject = occHandle(MeshVS_Mesh)::DownCast(it.value());
        if(aMeshObject.IsNull()) continue;
        const occHandle(MeshVS_DataSource) &aMeshDS = aMeshObject->GetDataSource();
        if(aMeshDS.IsNull()) continue;

        //! -------------------
        //! current body index
        //! -------------------
        int bodyIndex = it.key();

        //! ---------------------------
        //! initialize the mesh slicer
        //! ---------------------------
        aMeshSlicer.setMeshDataSource(aMeshDS);

        //! -----------------------------------------
        //! the hidden elements for the current mesh
        //! -----------------------------------------
        TColStd_PackedMapOfInteger hiddenElementsForCurrentMesh;

        //! ------------------------
        //! iterate over the planes
        //! ------------------------
        for(std::map<int,std::vector<double>>::const_iterator itplanes = mapOfClipPlanes.cbegin(); itplanes !=mapOfClipPlanes.cend() ; itplanes++)
        {
            //int planeID = itPlanes->first;
            const std::vector<double> &c = itplanes->second;
            occHandle(TColStd_HPackedMapOfInteger) hiddenElementsIDs;
            bool isDone = aMeshSlicer.perform(c[0],c[1],c[2],c[3],hiddenElementsIDs);
            if(isDone==false) continue;
            hiddenElementsForCurrentMesh.Unite(hiddenElementsIDs->Map());
        }
        occHandle(TColStd_HPackedMapOfInteger) hiddenHElementsForCurrentMesh = new TColStd_HPackedMapOfInteger;
        hiddenHElementsForCurrentMesh->ChangeMap() = hiddenElementsForCurrentMesh;

        std::map<int,occHandle(TColStd_HPackedMapOfInteger)>::iterator itt = myHiddenElements.find(bodyIndex);
        if(itt==myHiddenElements.end()) myHiddenElements.insert(std::make_pair(bodyIndex,hiddenHElementsForCurrentMesh));
        else itt->second = hiddenHElementsForCurrentMesh;
    }
}
