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
    connect(theDelegate,SIGNAL(currentCSChanged()),this,SLOT(updateCSDefinition()));

    //! --------------------------------
    //! update of the clip plane status
    //! --------------------------------
    connect(theDelegate,SIGNAL(currentCSStatusChanged()),this,SLOT(setClipPlaneActive()));

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

    //! ----------------------------------
    //! avoid crash at the very beginning
    //! (click on a mesh item)
    //! ----------------------------------
    if(myMDB==Q_NULLPTR)
    {
        cout<<"clipTool::setWorkingMode()->____cannot change the cliptool working mode: the mesh data base is null____"<<endl;
        return;
    }
    if(internalModel==Q_NULLPTR)
    {
        cerr<<"____the internal model pointer of the clip tool is NULL____"<<endl;
        return;
    }

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
    case 0:
    {
        if(myOCCViewer==Q_NULLPTR)
        {
            cerr<<"____the viewer pointer for the clip tool has not been set____"<<endl;
            return;
        }
        QMap<int,occHandle(Graphic3d_ClipPlane)> mapOfClipPlanes = myOCCViewer->getClipPlanes();
        for(QMap<int,occHandle(Graphic3d_ClipPlane)>::iterator it = mapOfClipPlanes.begin(); it!= mapOfClipPlanes.end(); ++it)
        {
            int clipPlaneNr = it.key();
            const occHandle(Graphic3d_ClipPlane) &curClipPlane = mapOfClipPlanes.value(clipPlaneNr);
            const double *coeff = curClipPlane->GetEquation().GetData();
            double a = *coeff;
            double b = *(coeff+1);
            double c = *(coeff+2);
            double d = *(coeff+3);

            const QMap<int,occHandle(AIS_InteractiveObject)> &meshObjects = myOCCViewer->getMeshObjects();
            for(QMap<int,occHandle(AIS_InteractiveObject)>::const_iterator it = meshObjects.cbegin(); it!=meshObjects.cend(); it++)
            {
                const occHandle(MeshVS_Mesh) &aMeshObject = occHandle(MeshVS_Mesh)::DownCast(it.value());
                const occHandle(MeshVS_DataSource) &aMeshDS = aMeshObject->GetDataSource();
                if(aMeshDS.IsNull()) continue;
                int bodyIndex = it.key();
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
        break;

    default:
    {
        if(myOCCViewer==Q_NULLPTR)
        {
            cerr<<"____the viewer pointer for the clip tool has not been set____"<<endl;
            return;
        }
        if(!myOCCViewer->getMeshContext().IsNull())
        {
            cout<<"____removing sliced meshes____"<<endl;
            myOCCViewer->getMeshContext()->EraseAll(true);
        }
    }
        break;
    }
    //! ----------------------------------------------------
    //! change the status on/off of the defined clip planes
    //! ----------------------------------------------------
    myOCCViewer->updateClipPlanes(activeClipPlanes);
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

        //! ---------------------------------
        //! add the clip plane to the viewer
        //! ---------------------------------
        myOCCViewer->addClipPlane(A,B,C,D,clipPlaneID,true);

        internalModel->appendRow(itemList);

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
void clipTool::retrieveClipPlanes()
{
    int Nb = internalModel->rowCount();
    for(int row = 0; row<Nb; row++)
    {
        QModelIndex index2 = internalModel->index(row,CLIPPLANE_STATUS_COLUMN);
        bool isActive = index2.data(Qt::UserRole).toBool();
        if(isActive==false) continue;

        QModelIndex index = internalModel->index(row,CLIPPLANE_BASE_COORDINATE_SYSTEM_COLUMN);
        void *CSp = internalModel->data(index,Qt::UserRole).value<void*>();
        QStandardItem* CS = (QStandardItem*)CSp;
        SimulationNodeClass *aNode = CS->data(Qt::UserRole).value<SimulationNodeClass*>();

        QModelIndex index1 = internalModel->index(row,CLIPPLANE_TRANSLATION_COLUMN);
    }
}

//! ------------------------------
//! function: removeItemFromTable
//! details:
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

    //! ------------------------------
    //! set the current clip plane ID
    //! ------------------------------
    myOCCViewer->setCurrentClipPlane(-1);
}

//! -----------------------------
//! function: updateCSDefinition
//! details:
//! -----------------------------
void clipTool::updateCSDefinition()
{
    cout<<"clipTool::updateCSDefinition()->____function called____"<<endl;

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

//! -----------------------------
//! function: switchOffClipPlane
//! details:
//! -----------------------------
void clipTool::setClipPlaneActive()
{
    int ID = internalModel->index(this->currentIndex().row(),CLIPPLANE_ID_COLUMN).data(Qt::UserRole).toInt();
    bool isOn = internalModel->index(this->currentIndex().row(),CLIPPLANE_STATUS_COLUMN).data(Qt::UserRole).toBool();
    myOCCViewer->setClipPlaneOn(ID, isOn);
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
