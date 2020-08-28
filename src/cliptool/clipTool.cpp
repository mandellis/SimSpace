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
    connect(theDelegate,SIGNAL(currentCSChanged()),this,SLOT(updateClipPlaneOfRow()));

    //! --------------------------------
    //! update of the clip plane status
    //! --------------------------------
    connect(theDelegate,SIGNAL(currentCSStatusChanged()),this,SLOT(handleClipPlaneOfRowStatus()));

    //! --------------------------------------------
    //! dynamic update after clip plane translation
    //! --------------------------------------------
    connect(theDelegate,SIGNAL(currentCSTranslationApplied(int)),this,SLOT(updateCSTranslation(int)));
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

    switch(myWorkingMode)
    {
    case 0: // on mesh
    {
        //! ----------------------------------------------------------
        //! recompute the hidden elements and send them to the viewer
        //! ----------------------------------------------------------
        //this->computeHiddenElements();
        //myOCCViewer->setHiddenElements(myHiddenElements);

        myOCCViewer->clipMesh();
    }
        break;

    case 3: // on solution
    {
        myOCCViewer->clipResult();
    }
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
    if(myCoordinateSystemRoot==Q_NULLPTR) return;
    if(myOCCViewer==Q_NULLPTR) return;

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
    QString clipPlaneName = QString("Clip plane %1").arg(myCurNumberOfClipPlanes+1);
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

    //! ----------------------------------------
    //! "Status": "0" => Off => "1" On: col = 2
    //! ----------------------------------------
    item = new QStandardItem();
    data.setValue(true);
    item->setData(data,Qt::UserRole);
    data.setValue(QString("On"));
    item->setData(data,Qt::DisplayRole);
    itemList<<item;

    //! ----------------------------------
    //! "Base" coordinate system: col = 3
    //! ----------------------------------
    item = new QStandardItem();
    QExtendedStandardItem *itemCSGlobal = static_cast<QExtendedStandardItem*>(myCoordinateSystemRoot->child(0,0));
    void *p = (void*)itemCSGlobal;
    data.setValue(p);
    item->setData(data,Qt::UserRole);
    data.setValue(itemCSGlobal->data(Qt::UserRole).value<SimulationNodeClass*>()->getName());
    item->setData(data,Qt::DisplayRole);
    itemList<<item;

    //! ---------------------------------------------------------------------------------
    //! create an adjacent cell with the data of the (global) coordinate system: col = 4
    //! ---------------------------------------------------------------------------------
    const QVector<double> &coeff = this->getPlaneCoefficients(myCoordinateSystemRoot->child(0,0));
    data.setValue(coeff);

    item = new QStandardItem();
    item->setData(data,Qt::UserRole);
    double A,B,C,D;
    A = coeff[0]; B = coeff[1]; C = coeff[2]; D = coeff[3];
    data.setValue(QString("(%1, %2, %3, %4)").arg(A).arg(B).arg(C).arg(D));
    item->setData(data,Qt::DisplayRole);
    itemList<<item;

    //! -----------------------
    //! "Translation": col = 5
    //! -----------------------
    item = new QStandardItem();
    double deltaZ = 0.0;
    data.setValue(deltaZ);
    item->setData(data,Qt::UserRole);
    data.setValue(QString("%1").arg(deltaZ));
    item->setData(data,Qt::DisplayRole);
    itemList<<item;

    //! --------------------------------------
    //! "Shifted plane coefficients": col = 6
    //! At the beginning the same of 4-th col
    //! since translation is zero
    //! --------------------------------------
    item = new QStandardItem();
    data.setValue(coeff);

    item = new QStandardItem();
    item->setData(data,Qt::UserRole);
    data.setValue(QString("(%1, %2, %3, %4)").arg(A).arg(B).arg(C).arg(D));
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
    horizontalHeaderLabels<<"Plane"<<"ID"<<"Status"<<"Csys"<<"Base plane"<<"Shift"<<"Shifted plane";
    internalModel->setHorizontalHeaderLabels(horizontalHeaderLabels);

    this->verticalHeader()->hide();
    this->horizontalHeader()->setVisible(true);
    this->verticalScrollBar()->setVisible(true);
    this->setAlternatingRowColors(true);

    //! ------------------
    //! hide some columns
    //! ------------------
    this->setColumnHidden(CLIPPLANE_ID_COLUMN,true);
    this->setColumnHidden(CLIPPLANE_BASE_PLANE_DATA_COLUMN,true);
    this->setColumnHidden(CLIPPLANE_SHIFTED_PLANE_COEFFICIENTS_COLUMN,true);

    myCurNumberOfClipPlanes++;

    //! -----------------
    //! single selection
    //! -----------------
    this->setSelectionMode(QAbstractItemView::SingleSelection);

    //! ---------------------------
    //! hide the horizontal header
    //! ---------------------------
    this->horizontalHeader()->setVisible(false);
    this->resizeColumnsToContents();
    this->setColumnWidth(CLIPPLANE_TRANSLATION_COLUMN,80);
}

//! -----------------------------
//! function: retrieveClipPlanes
//! details:
//! -----------------------------
int clipTool::retrieveActiveClipPlanes(std::map<int,std::vector<double>> &mapOfClipPlanes)
{
    int NbClipPlanes = internalModel->rowCount();
    if(NbClipPlanes==0) return 0;

    for(int row = 0; row<NbClipPlanes; row++)
    {
        QModelIndex index1 = internalModel->index(row,CLIPPLANE_STATUS_COLUMN);
        bool isActive = index1.data(Qt::UserRole).toBool();
        if(isActive==false) continue;

        QModelIndex index4 = internalModel->index(row,CLIPPLANE_ID_COLUMN);
        int planeID = index4.data(Qt::UserRole).toInt();

        QModelIndex index = internalModel->index(row,CLIPPLANE_SHIFTED_PLANE_COEFFICIENTS_COLUMN);
        const QVector<double> &coeffs = index.data(Qt::UserRole).value<QVector<double>>();

        std::vector<double> aPlaneCoeffs {coeffs[0],coeffs[1],coeffs[2],coeffs[3]};
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
        internalModel->removeRow(curClickedIndex.row());
        myOCCViewer->removeClipPlane(clipPlaneID);
        myCurNumberOfClipPlanes--;
    }

    //! --------------------------
    //! recompute hidden elements
    //! --------------------------
    myOCCViewer->clipMesh();
    myOCCViewer->clipResult();
}

//! -----------------------------
//! function: updateCSDefinition
//! details:
//! -----------------------------
void clipTool::updateCSDefinition()
{
    cout<<"clipTool::updateCSDefinition()->____function called____"<<endl;
    myOCCViewer->clipMesh();
    myOCCViewer->clipResult();
}

//! -------------------------------------------
//! function: updateClipPlanesByExternalChange
//! details:
//! -------------------------------------------
void clipTool::updateCSDataByExternalCSChange(QStandardItem *theCurrentModifiedCS)
{    
    cout<<"clipTool::updateCSDataByExternalCSChange()->____function called____"<<endl;
    for(int row = 0; row< this->internalModel->rowCount(); row++)
    {
        QModelIndex index = internalModel->index(row,CLIPPLANE_BASE_COORDINATE_SYSTEM_COLUMN);
        if(index.data(Qt::UserRole).canConvert<void*>() == false) cout<<"____ERROR IN CONVERSION____"<<endl;
        else
        {
            cout<<"____CONVERSION OK____"<<endl;
            void *p = index.data(Qt::UserRole).value<void*>();
            QStandardItem *item = static_cast<QStandardItem*>(p);
            cout<<"____"<<item->data(Qt::DisplayRole).toString().toStdString()<<"____"<<endl;
        }
    }

    SimulationNodeClass *nodeCS = theCurrentModifiedCS->data(Qt::UserRole).value<SimulationNodeClass*>();
    const QString &extTimeTag = nodeCS->getPropertyValue<QString>("Time tag");

    cout<<"____modified CS tag: "<<extTimeTag.toStdString()<<"____"<<endl;

    //! -------------------------------------------------
    //! coefficients of CLIPPLANE_BASE_PLANE_DATA_COLUMN
    //! -------------------------------------------------
    const QVector<double> &coeffs = this->getPlaneCoefficients(theCurrentModifiedCS);
    cout<<"____new coefficients("<<coeffs[0]<<", "<<coeffs[1]<<", "<<coeffs[2]<<", "<<coeffs[3]<<")____"<<endl;

    //! ---------------------------
    //! scan the rows of the table
    //! ---------------------------
    for(int row = 0; row< this->internalModel->rowCount(); row++)
    {
        QModelIndex index = internalModel->index(row,CLIPPLANE_BASE_COORDINATE_SYSTEM_COLUMN);
        void *p = index.data(Qt::UserRole).value<void*>();
        QStandardItem *curItem = (QStandardItem*)p;
        if(p==NULL) cout<<"____nullptr____"<<endl;
        SimulationNodeClass *curNode = curItem->data(Qt::UserRole).value<SimulationNodeClass*>();
        const QString &curTimeTag = curNode->getPropertyValue<QString>("Time tag");
        cout<<"____scanning CS with tag: "<<curTimeTag.toStdString()<<"____"<<endl;
        if(curTimeTag!=extTimeTag) continue;
        cout<<"____FOUND____"<<endl;

        //! ----------------------------
        //! update the row of the table
        //! ----------------------------
        QModelIndex indexBasePlaneCoeffs = internalModel->index(row,CLIPPLANE_BASE_COORDINATE_SYSTEM_COLUMN);
        QVariant value;
        value.setValue(coeffs);
        internalModel->setData(indexBasePlaneCoeffs,value,Qt::UserRole);
        value.setValue(QString("(%1, %2, %3, %4)").arg(coeffs[0]).arg(coeffs[1]).arg(coeffs[2]).arg(coeffs[3]));
        cout<<"____old coefficients: "<<value.toString().toStdString()<<"____"<<endl;
        internalModel->setData(indexBasePlaneCoeffs,value,Qt::DisplayRole);

        //! --------------------
        //! translate the plane
        //! --------------------
        QModelIndex indexPlaneTranslation = internalModel->index(row,CLIPPLANE_TRANSLATION_COLUMN);
        int sliderValue = indexPlaneTranslation.data(Qt::UserRole).toInt();
        double lx,ly,lz;
        myOCCViewer->getSceneBoundingBox(lx,ly,lz);
        double D = sqrt(lx*lx+ly*ly+lz*lz);
        double translation = double(sliderValue/100.0)*(D/1.0);
        cout<<"____slider bar position: "<<sliderValue<<"____"<<endl;
        double a = coeffs[0];
        double b = coeffs[1];
        double c = coeffs[2];
        double d = coeffs[3];
        this->translatePlane(a,b,c,d,translation);

        QVector<double> coeffs_trans {a,b,c,d};
        QModelIndex indexShiftedPlane = internalModel->index(row,CLIPPLANE_SHIFTED_PLANE_COEFFICIENTS_COLUMN);
        value.setValue(coeffs_trans);
        internalModel->setData(indexShiftedPlane,value,Qt::UserRole);
        value.setValue(QString("(%1, %2, %3, %4").arg(a).arg(b).arg(c).arg(d));
        cout<<"____updated coefficients: "<<QString("(%1, %2, %3, %4").arg(a).arg(b).arg(c).arg(d).toStdString()<<"____"<<endl;
        internalModel->setData(indexShiftedPlane,value,Qt::DisplayRole);
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

//! -------------------------------
//! function: updateClipPlaneOfRow
//! details:
//! -------------------------------
void clipTool::updateClipPlaneOfRow()
{
    //cout<<"clipTool::updateClipPlaneOfRow()->____function called____"<<endl;
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
    a /= l; b /=l; c /=l; d /=l;

    QModelIndex index = this->currentIndex();   // cell of clip plane selector

    //! ----------------------------------
    //! update the cell of the base plane //cesere
    //! ----------------------------------
    QModelIndex index1 = index.sibling(index.row(),CLIPPLANE_BASE_PLANE_DATA_COLUMN);
    QVariant value;
    QVector<double> coeffs_base_plane{a,b,c,d};
    value.setValue(coeffs_base_plane);
    internalModel->setData(index1,value,Qt::UserRole);
    value.setValue(QString("(%1 ,%2 ,%3 , %4)").arg(a).arg(b).arg(c).arg(d));
    internalModel->setData(index1,value,Qt::DisplayRole);

    //! -------------------------------------
    //! update the cell of the shifted plane
    //! -------------------------------------
    QModelIndex indexTranslation = internalModel->index(index.row(),CLIPPLANE_TRANSLATION_COLUMN);
    int sliderPosition = indexTranslation.data(Qt::UserRole).toInt();

    double lx,ly,lz;
    myOCCViewer->getSceneBoundingBox(lx,ly,lz);
    double translation = sqrt(lx*lx+ly*ly+lz*lz)*double(sliderPosition/100.0);
    this->translatePlane(a,b,c,d,translation);
    QVector<double> coeffs_shifted_plane {a,b,c,d};

    QModelIndex index2 = index.sibling(index.row(),CLIPPLANE_SHIFTED_PLANE_COEFFICIENTS_COLUMN);
    value.setValue(coeffs_shifted_plane);
    internalModel->setData(index2,value,Qt::UserRole);
    value.setValue(QString("(%1 ,%2 ,%3 , %4").arg(a).arg(b).arg(c).arg(d));
    internalModel->setData(index2,value,Qt::DisplayRole);

    int clipPlaneID = internalModel->index(this->currentIndex().row(),CLIPPLANE_ID_COLUMN).data(Qt::UserRole).toInt();
    myOCCViewer->addClipPlane(a,b,c,d,clipPlaneID,true);

    //! --------------------------
    //! recompute hidden elements
    //! --------------------------
    myOCCViewer->clipMesh();
    myOCCViewer->clipResult();
}

//! ------------------------------
//! function: updateCSTranslation
//! details:
//! ------------------------------
void clipTool::updateCSTranslation(int sliderPosition)
{
    //cout<<"clipTool::updateCSTranslation()->____slider position: "<<sliderPosition<<"____"<<endl;
    double lx,ly,lz;
    myOCCViewer->getSceneBoundingBox(lx,ly,lz);
    double BBXDiagonal = sqrt(lx*lx+ly*ly+lz*lz);

    //! ------------------------------------------
    //! get the current CS undergoing translation
    //! and the plane coefficients
    //! ------------------------------------------
    int ID = internalModel->index(this->currentIndex().row(),CLIPPLANE_ID_COLUMN).data(Qt::UserRole).toInt();
    const QVector<double> &coeffs = internalModel->index(this->currentIndex().row(),CLIPPLANE_BASE_PLANE_DATA_COLUMN).data(Qt::UserRole).
            value<QVector<double>>();

    //! --------------------
    //! compute translation
    //! --------------------
    double translation = double(sliderPosition/100.0)*(BBXDiagonal/1.0);

    double a = coeffs[0];
    double b = coeffs[1];
    double c = coeffs[2];
    double d = coeffs[3];
    this->translatePlane(a,b,c,d,translation);

    //! --------------------------------------------
    //! update the table shifted plane coefficients
    //! --------------------------------------------
    QVariant value;
    QVector<double> shifted_plane_coeffs {a,b,c,d};
    QModelIndex index = this->currentIndex().sibling(this->currentIndex().row(),CLIPPLANE_SHIFTED_PLANE_COEFFICIENTS_COLUMN);
    value.setValue(shifted_plane_coeffs);
    internalModel->setData(index,value,Qt::UserRole);
    value.setValue(QString("%1, %2, %3, %4").arg(a).arg(b).arg(c).arg(d));
    internalModel->setData(index,value,Qt::DisplayRole);

    //! -------------------------------------------------
    //! this moves the plane and redraw the clipped view
    //! -------------------------------------------------
    myOCCViewer->updateClipPlaneCoefficients(ID,shifted_plane_coeffs);

    switch(myWorkingMode)
    {
    case 0:
    {
        if(myMDB==Q_NULLPTR) return;
        myOCCViewer->clipMesh();
    }
        break;
    case 3:
    {
        myOCCViewer->clipMesh();
        myOCCViewer->clipResult();
    }
        break;
    }
}

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

//! -------------------------------------
//! function: handleClipPlaneOfRowStatus
//! details:
//! -------------------------------------
void clipTool::handleClipPlaneOfRowStatus()
{
    //! --------------------------
    //! handle the view of shapes
    //! --------------------------
    int clipPlaneID = internalModel->index(this->currentIndex().row(),CLIPPLANE_ID_COLUMN).data(Qt::UserRole).toInt();
    bool isOn = internalModel->index(this->currentIndex().row(),this->currentIndex().column()).data(Qt::UserRole).toBool();
    myOCCViewer->setClipPlaneOn(clipPlaneID,isOn);

    //! ----------------------
    //! handle the mesh slice
    //! ----------------------
    myOCCViewer->clipMesh();
    myOCCViewer->clipResult();
}

//! -------------------------
//! function: translatePlane
//! details:
//! -------------------------
void clipTool::translatePlane(double &a, double &b,double &c,double &d, const double &t)
{
    double l = sqrt(a*a+b*b+c*c);
    a /= l; b/=l; c/=l; d/=l;
    d += t;
}
