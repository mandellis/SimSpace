//! custom includes
#include "qcustomplotextended.h"

//! Qt
#include <QRubberBand>

#define TRACER_SIZE 10

//! ----------------
//! file extensions
//! ----------------
#define PNG_FILES ".png"
#define JPG_FILES ".jpg"
#define BMP_FILES ".bmp"
#define PDF_FILES ".pdf"

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
QCustomPlotExtended::QCustomPlotExtended(QWidget *parent): QCustomPlot(parent),
    myX(0),myY(0),myCurAction(action::action_none),isTracingEnabled(true)
{
    this->setContextMenuPolicy(Qt::CustomContextMenu);
    this->setInteraction(QCP::iSelectItems);

    //! --------
    //! actions
    //! --------
    actionSaveImage = new QAction("Save image",this);
    actionSaveImage->setIcon(QIcon(":/icons/icon_save image.png"));

    connect(actionSaveImage,SIGNAL(triggered(bool)),this,SLOT(saveImage()));
    cout<<"____tag022____"<<endl;
    //! -------------------
    //! add a phase tracer
    //! -------------------
    myTracer = new QCPItemTracer(this);
    cout<<"____tag023____"<<endl;

    myTracer->setVisible(false);
    myTracer->setInterpolating(false);
    myTracer->setPen(QPen(Qt::red));
    myTracer->setBrush(Qt::red);
    myTracer->setSize(TRACER_SIZE); // the size is in pixels
    myTracer->setStyle(QCPItemTracer::tsCircle);
    myTracer->setClipToAxisRect(false);

    cout<<"____tag03____"<<endl;

    //! ------------
    //! rubber band
    //! ------------
    myRectBand = new QRubberBand(QRubberBand::Rectangle, this);
    myRectBand->setStyle(QStyleFactory::create("windows"));

    cout<<"____tag04____"<<endl;

    //! ------------------------
    //! create the context menu
    //! ------------------------
    myContextMenu = new QMenu(this);
    connect(this,SIGNAL(customContextMenuRequested(QPoint)),this,SLOT(ShowContextMenu(QPoint)));

    //! ---------------------
    //! the very custom menu
    //! ---------------------
    myVeryCustomMenu = new QMenu(this);
}

//! -------------------------
//! function: showPointValue
//! details:
//! -------------------------
void QCustomPlotExtended::showPointValue(QMouseEvent* event)
{
    QCPGraph *graph = this->graph();
    if(graph!=NULL)
    {
        //! ----------------------
        //! Setup the item tracer
        //! ----------------------
        if(myTracer->visible()==false) myTracer->setVisible(true);
        this->myTracer->setGraph(graph);
        myTracer->setGraphKey(this->xAxis->pixelToCoord(event->pos().x()));

        this->replot();

        QPointF temp = myTracer->position->coords();
        emit tracerInterceptsTime(temp.x());

        /*
        //! --------------------------------------
        //! Show a tooltip which tells the values
        //! --------------------------------------
        QPointF temp = myTracer->position->coords();
        QToolTip::showText(event->globalPos(),
                           tr("<h4>%L1</h4>"
                              "<table>"
                              "<tr>"
                              "<td><h5>X: %L2</h5></td>" "<td>  ,  </td>" "<td><h5>Y: %L3</h5></td>"
                              "</tr>"
                              "</table>").arg(graph->name()).arg(QString::number(temp.x(),'f',0))
                           .arg(QString::number(temp.y(),'f',1)),this,this->rect());
        */
    }
}

//! --------------------------
//! function: ShowContextMenu
//! details:
//! --------------------------
void QCustomPlotExtended::ShowContextMenu(QPoint pos)
{
    //cout<<"____context menu activated____"<<endl;
    QPoint globalPos = this->mapToGlobal(pos);

    myContextMenu->addAction(actionSaveImage);
    myContextMenu->exec(globalPos);
}

//! --------------------------
//! function: mousePressEvent
//! details:
//! --------------------------
void QCustomPlotExtended::mousePressEvent(QMouseEvent *event)
{
    xinitmin = this->xAxis->range().lower;
    yinitmin = this->yAxis->range().lower;
    xinitmax = this->xAxis->range().upper;
    yinitmax = this->yAxis->range().upper;

    myX_start = event->pos().x();
    myY_start = event->pos().y();

    myX = event->pos().x();
    myY = event->pos().y();

    if(event->button()==Qt::LeftButton)
    {
        //cout<<"____lb pressed"<<endl;

        QPoint P(myX,myY);
        QCPAbstractItem *theItem = this->itemAt(P,true);
        if(theItem==NULL)
        {
            //cerr<<"___NULL selection: item not updated____"<<endl;
            myCurSelectedItem = NULL; //this resets the previous selection
        }
        else
        {
            //cout<<"____lb pressed____"<<endl;
            myCurSelectedItem = theItem;

            //! ----------------------------------
            //! the item tracer has been clicked
            //! activate the very custom menu
            //! ----------------------------------
            QCPItemTracer *it = qobject_cast<QCPItemTracer*>(theItem);
            if(it)
            {
                //! ---------------------------
                //! create the interactive menu
                //! ---------------------------
                //cout<<"____build menu____"<<endl;
                this->fillVeryCustomMenu(it);

                //! -----------------------------
                //! handling the selected action
                //! check if the action is valid
                //! -----------------------------
                QAction *sa = myVeryCustomMenu->exec(event->globalPos());
                if(sa)
                {
                    TimeStepType ts = sa->data().value<TimeStepType>();
                    switch(ts)
                    {
                    case TimeStepType::ToBeRemoved:
                    {
                        //! ------------------
                        //! removing a marker
                        //! ------------------
                        double time = myTracer->position->key();
                        for(int i=0; i<TimeStepMarkerList.length(); i++)
                        {
                            const TimeStepMarker &tsm = TimeStepMarkerList.at(i);
                            double curTime = tsm.time;
                            if(curTime==time)
                            {
                                //cout<<"____REMOVING A MARKER____"<<endl;

                                //! remove the vertical line from the screen
                                this->removeItem(vLineMap.value(tsm));

                                this->replot();

                                //! remove the vertical line from the list
                                vLineMap.remove(tsm);

                                //! remove the time step marker from the list
                                TimeStepMarkerList.removeAt(i);

                                //! update boxes
                                this->updateBoxes();
                            }
                        }
                        //cout<<"____NUMBER OF MARKERS: "<<TimeStepMarkerList.length()<<"____"<<endl;
                    }
                        break;

                    case TimeStepType::ClosedAssemblyWithoutPressure:
                    case TimeStepType::ClosedAssemblyWithPressure:
                    case TimeStepType::OpenAssembly:
                    {
                        //! --------------------------
                        //! create a time step marker
                        //! --------------------------
                        TimeStepMarker aMarker;
                        aMarker.confirmed = true;
                        aMarker.time = myTracer->position->key();
                        aMarker.value = myTracer->position->value();
                        aMarker.tsType = ts;

                        //cout<<"____ADDING TIME STEP MARKER: ____"<<aMarker.time<<"____"<<endl;
                        TimeStepMarkerList<<aMarker;

                        //! -----------------------
                        //! create a vertical line
                        //! -----------------------
                        QCPItemLine *vLine = new QCPItemLine(this);
                        vLine->setPen(QPen(Qt::blue, 0.5, Qt::DashLine));

                        double y_start = -1e80; //this->yAxis->pixelToCoord(0);
                        double y_end = 1e80; //this->yAxis->pixelToCoord(this->size().height());

                        vLine->start->setCoords(aMarker.time,y_start);
                        vLine->end->setCoords(aMarker.time,y_end);

                        vLineMap.insert(aMarker,vLine);

                        //! -------------
                        //! update boxes
                        //! -------------
                        this->updateBoxes();

                        //cout<<"____NUMBER OF MARKERS: "<<TimeStepMarkerList.length()<<"____"<<endl;
                        //cout<<"____NUMBER OF RECT: "<<rectMap.size()<<"____"<<endl;
                    }
                        break;
                    }
                }
            }
        }
    }
    else if(event->button()==Qt::RightButton)
    {
        //cout<<"____rb pressed: unset action drag____"<<endl;
        myCurAction = action::action_none;
    }
}

//! ----------------------------
//! function: mouseReleaseEvent
//! details:
//! ----------------------------
void QCustomPlotExtended::mouseReleaseEvent(QMouseEvent *event)
{
    switch(event->button())
    {
    case Qt::LeftButton:
    {
        //cout<<"____lb released____"<<endl;

        myX = event->pos().x();
        myY = event->pos().y();
        if(myCurAction==action::action_zoom)
        {
            //! ---------------------
            //! hide the rubber band
            //! ---------------------
            if(myRectBand) myRectBand->hide();

            //! ------------------------------------------
            //! avoid the zoom at a single press=>release
            //! ------------------------------------------
            if(myX_start==myX && myY_start==myY) return;

            //cout<<"____NOW ZOOMING____"<<endl;

            double x_min, x_max, y_min, y_max;
            x_min = this->xAxis->pixelToCoord(myX_start);
            x_max = this->xAxis->pixelToCoord(event->pos().x());
            y_min = this->yAxis->pixelToCoord(myY_start);
            y_max = this->yAxis->pixelToCoord(event->pos().y());

            this->xAxis->setRangeLower(x_min);
            this->xAxis->setRangeUpper(x_max);
            this->yAxis->setRangeLower(y_min);
            this->yAxis->setRangeUpper(y_max);

            this->replot();
        }
        isTracingEnabled = true;
        //myCurAction = action::action_none;
    }
        break;
    }
}

//! -------------------------
//! function: mouseMoveEvent
//! details:
//! -------------------------
void QCustomPlotExtended::mouseMoveEvent(QMouseEvent *event)
{
    myX = event->pos().x();
    myY = event->pos().y();

    if(event->buttons() & Qt::LeftButton)
    {
        isTracingEnabled = false;
        if(myCurAction == action::action_drag)
        {
            //! avoid crashes when calling objectName on a NULL object
            if(myCurSelectedItem!=NULL)
            {
                ;
            }
        }
        if(myCurAction == action::action_zoom)
        {
            double deltaX = myX-myX_start;
            double deltaY = myY-myY_start;
            if(deltaX> 10 && deltaY > 10)
            {
                cout<<"___deltaX: "<<myX-myX_start<<" deltaY: "<<myY-myY_start<<"____"<<endl;
                this->drawRubberBand(myX_start,myY_start,myX,myY);
            }
        }
        if(myCurAction == action::action_pan)
        {
            double deltaX = myX_start-myX;
            double deltaY = myY_start-myY;

            double kx = (this->xAxis->pixelToCoord(myX_start)-this->xAxis->pixelToCoord(myX))/deltaX;
            double ky = (this->yAxis->pixelToCoord(myY_start)-this->yAxis->pixelToCoord(myY))/deltaY;

            deltaX =deltaX*kx;
            deltaY =deltaY*ky;

            if(fabs(deltaX)>=1 && fabs(deltaY)>=1)
            {
                this->xAxis->setRangeLower(xinitmin+deltaX);
                this->xAxis->setRangeUpper(xinitmax+deltaX);
                this->yAxis->setRangeLower(yinitmin+deltaY);
                this->yAxis->setRangeUpper(yinitmax+deltaY);
                this->replot();
            }
        }
    }

    if(isTracingEnabled)
    {
        showPointValue(event);
    }
}

//! --------------
//! function: fit
//! details:
//! --------------
void QCustomPlotExtended::fit()
{
    //cout<<"____fit____"<<endl;
    myCurAction = action_fit;
    isTracingEnabled = true;
    if(this->graph())
    {
        bool foundRange;
        QCPRange rangeX = this->graph(0)->getKeyRange(foundRange);
        QCPRange rangeY = this->graph(0)->getValueRange(foundRange);
        if(foundRange)
        {
            double xmin = rangeX.lower;
            double xmax = rangeX.upper;
            double ymin = rangeY.lower;
            double ymax = rangeY.upper;
            double deltax = xmax-xmin;
            double deltay = ymax-ymin;

            this->xAxis->setRangeLower(xmin-deltax*0.025);
            this->yAxis->setRangeLower(ymin-deltay*0.025);
            this->xAxis->setRangeUpper(xmax+deltax*0.025);
            this->yAxis->setRangeUpper(ymax+deltay*0.025);

            this->yAxis->scaleRange(1.1, this->yAxis->range().center());
            this->xAxis->scaleRange(1.1, this->xAxis->range().center());

            //! ----------------------
            //! update vertical lines
            //! ----------------------
            for(QMap<TimeStepMarker,QCPItemLine*>::iterator it = vLineMap.begin(); it!=vLineMap.end(); ++it)
            {
                QCPItemLine* curVLine =it.value();

                double y_start = this->yAxis->pixelToCoord(0);
                double y_end = this->yAxis->pixelToCoord(this->size().height());

                double x_start = curVLine->start->coords().x();
                curVLine->start->setCoords(x_start,y_start);
                curVLine->end->setCoords(x_start,y_end);
            }

            this->replot();
        }
    }
}

//! ---------------
//! function: zoom
//! details:
//! ---------------
void QCustomPlotExtended::zoom()
{
    //cout<<"____zoom____"<<endl;
    myCurAction = action_zoom;
    isTracingEnabled = false;
}

//! --------------
//! function: pan
//! details:
//! --------------
void QCustomPlotExtended::pan()
{
    //cout<<"____pan____"<<endl;
    myCurAction = action_pan;
    isTracingEnabled = false;
}

//! -----------------------------
//! function: fillVeryCustomMenu
//! details:
//! -----------------------------
void QCustomPlotExtended::fillVeryCustomMenu(QCPItemTracer *aTimeTracer)
{
    myVeryCustomMenu->clear();

    //! -----------------------------------------------
    //! add actions only if a time marker has not been
    //! previously set at this time point
    //! -----------------------------------------------
    bool found = false;
    double time =aTimeTracer->position->key();
    for(int i=0; i<TimeStepMarkerList.length(); i++)
    {
        if(TimeStepMarkerList.at(i).time==time)
        {
            found = true;
            break;
        }
    }
    if(!found)
    {
        QList<QAction*> listOfActions;

        //! --------------------
        //! very custom actions
        //! --------------------
        QAction *actionInsertTimeStepType_AssemblyClosedNoPressure = new QAction("Closed assembly no pressure",this);
        actionInsertTimeStepType_AssemblyClosedNoPressure->setIcon(QIcon(":/icons/icon_closed.png"));
        actionInsertTimeStepType_AssemblyClosedNoPressure->setData(TimeStepType::ClosedAssemblyWithoutPressure);

        QAction *actionInsertTimeStepType_AssemblyClosedAndPressure = new QAction("Closed assembly plus pressure",this);
        actionInsertTimeStepType_AssemblyClosedAndPressure->setIcon(QIcon(":/icons/icon_closed.png"));
        actionInsertTimeStepType_AssemblyClosedAndPressure->setData(TimeStepType::ClosedAssemblyWithPressure);

        QAction *actionInsertTimeStepType_AssemblyOpen = new QAction("Open assembly",this);
        actionInsertTimeStepType_AssemblyOpen->setIcon(QIcon(":/icons/icon_opened.png"));
        actionInsertTimeStepType_AssemblyOpen->setData(TimeStepType::OpenAssembly);

        listOfActions<<actionInsertTimeStepType_AssemblyClosedNoPressure<<
                       actionInsertTimeStepType_AssemblyClosedAndPressure<<
                       actionInsertTimeStepType_AssemblyOpen;

        myVeryCustomMenu->addActions(listOfActions);
    }

    //! ------------------------------------------
    //! add the clear action only if at least one
    //! time marker has been added
    //! ------------------------------------------
    if(!TimeStepMarkerList.isEmpty())
    {
        //! add a separator
        myVeryCustomMenu->addSeparator();

        QAction *actionClear = new QAction("Clear marker",this);
        actionClear->setIcon(QIcon(":/icons/icon_clear.png"));
        actionClear->setData(TimeStepType::ToBeRemoved);

        myVeryCustomMenu->addAction(actionClear);
    }
}

//! ---------------------
//! function: clearPanel
//! details:
//! ---------------------
void QCustomPlotExtended::clearPanel()
{
    this->clearGraphs();

    //! ----------------
    //! hide the tracer
    //! ----------------
    myTracer->setVisible(false);

    //! ----------------------------------------------
    //! remove all the vertical lines from the screen
    //! ----------------------------------------------
    for(QMap<TimeStepMarker,QCPItemLine*>::iterator it = vLineMap.begin(); it!=vLineMap.end(); ++it)
    {
        QCPItemLine *curLine = it.value();
        this->removeItem(curLine);
    }
    //! --------------
    //! clear the map
    //! --------------
    vLineMap.clear();

    //! ---------------------
    //! remove all the boxes
    //! ---------------------
    TimeStepMarkerList.clear();

    //! -------------
    //! clear labels
    //! -------------
    this->xAxis->setLabel("");
    this->yAxis->setLabel("");

    //! -------------
    //! update boxes
    //! -------------
    updateBoxes();

    this->replot();
}

//! -------------------------
//! function: drawRubberBand
//! details:
//! -------------------------
void QCustomPlotExtended::drawRubberBand(int tlx, int tly, int brx, int bry)
{
    QRect aRect;

    //! ------------------
    //! set the rectangle
    //! ------------------
    aRect.setCoords(tlx,tly,brx,bry);

    myRectBand->setGeometry(aRect);
    myRectBand->show();
}

//! -------------------
//! function: getBrush
//! details:
//! -------------------
QBrush QCustomPlotExtended::getBrush(const TimeStepMarker &aTimeStepMarker)
{
    QBrush aBrush;
    aBrush.setStyle(Qt::SolidPattern);
    switch(aTimeStepMarker.tsType)
    {
    case ClosedAssemblyWithoutPressure:

        aBrush.setColor(QColor(0,0,255,20));
        break;

    case ClosedAssemblyWithPressure:
        aBrush.setColor(QColor(0,255,0,20));
        break;

    case OpenAssembly:
        aBrush.setColor(QColor(255,0,0,20));
        break;
    }
    return aBrush;
}

//! ----------------------
//! function: updateBoxes
//! details:
//! ----------------------
void QCustomPlotExtended::updateBoxes()
{
    //cout<<"____updateBoxes____"<<endl;

    double ymax = 1e80; //this->yAxis->range().upper;
    double ymin = -1e80; //this->yAxis->range().lower;

    //! ---------------------------------------------
    //! remove all the colored boxes from the screen
    //! ---------------------------------------------
    for(QMap<TimeStepMarker,QCPItemRect*>::iterator it = rectMap.begin(); it!=rectMap.end(); ++it)
    {
        //! remove the box from the graph
        this->removeItem(it.value());
    }
    rectMap.clear();
    //cout<<"____number of rectangles: "<<rectMap.size()<<"____"<<endl;
    this->replot();
    Sleep(125);

    double previousTime = 0; // minimum at left
    double nextTime = 0;

    std::sort(TimeStepMarkerList.begin(),TimeStepMarkerList.end());

    for(int i=0; i<TimeStepMarkerList.length(); i++)
    {
        const TimeStepMarker &curTimeStepMarker = TimeStepMarkerList.at(i);
        nextTime = curTimeStepMarker.time;

        QCPItemRect *aRect = new QCPItemRect(this);
        const QBrush &curBrush = this->getBrush(curTimeStepMarker);
        aRect->setPen(curBrush.color());
        aRect->setBrush(curBrush);
        aRect->topLeft->setCoords(previousTime,ymax);
        aRect->bottomRight->setCoords(nextTime,ymin);

        rectMap.insert(curTimeStepMarker,aRect);

#ifdef DEBUG
        QString name;
        switch(curTimeStepMarker.tsType)
        {
        case ClosedAssemblyWithoutPressure: name = "Closed assembly without pressure"; break;
        case ClosedAssemblyWithPressure: name ="Closed assembly with pressure"; break;
        case OpenAssembly: name = "Open assembly"; break;
        }
        cout<<"____"<<i+1<<"____("<<previousTime<<", "<<nextTime<<") "<<name.toStdString()<<"____"<<endl;
#endif

        previousTime = nextTime;
    }
    this->replot();
}

//! --------------------
//! function: saveImage
//! details:
//! --------------------
void QCustomPlotExtended::saveImage()
{
    QString selectedFilter;
    QString fileName("D:/file.png");
    //QString fileName = QFileDialog::getSaveFileName(0,"Save image as ",QDir::currentPath(),PNG_FILES";;"JPG_FILES";;"BMP_FILES,&selectedFilter,0);

    if(!fileName.isEmpty() && this->isVisible())
    {
        fileName.append(selectedFilter);
        if(selectedFilter==".jpg") this->saveJpg(fileName);
        if(selectedFilter==".bmp") this->saveBmp(fileName);
        if(selectedFilter==".png") this->savePng(fileName);
    }
}

//! -------------------
//! function: fitRange
//! details:
//! -------------------
void QCustomPlotExtended::fitRange()
{
    this->rescaleAxes();
    this->replot();
}
