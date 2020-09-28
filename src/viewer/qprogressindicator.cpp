//! ----------------
//! custom includes
//! ----------------
#include "qprogressindicator.h"
#include "qprogressevent.h"
#include <tetgenmesher.h>

//! ---
//! Qt
//! ---
#include <QVBoxLayout>
#include <QGridLayout>
#include <QProgressBar>
#include <QPushButton>
#include <QLabel>
#include <QProcess>
#include <QThread>

//! ----------------
//! system specific
//! ----------------
#include <windows.h>

//! ----
//! C++
//! ----
#include <iostream>
using namespace std;

//! -------
//! global
//! -------
#include "global.h"

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
QProgressIndicator::QProgressIndicator(bool additionalBar, QWidget *parent):QWidget(parent)
{
    cout<<"QProgressIndicator::QProgressIndicator()->____constructor called____"<<endl;
    this->createContent(additionalBar);
    this->installEventFilter(this);
    this->hide();
    this->setWindowFlags(Qt::WindowStaysOnTopHint);
}

//! -------------------------------------------
//! function: destructor
//! details:  defined for checking destruction
//! -------------------------------------------
QProgressIndicator::~QProgressIndicator()
{
    cout<<"QProgressIndicator::~QProgressIndicator()->____DESTRUCTOR CALLED____"<<endl;
}

//! ------------------------
//! function: createContent
//! details:
//! ------------------------
void QProgressIndicator::createContent(bool additionalBar)
{
    this->setWindowTitle("Progress indicator");
    this->setObjectName("progressIndicator");
    this->setFixedWidth(400);

    //! --------------------------
    //! style of the progress bar
    //! --------------------------
    QString progressBarStyleSheet = "QProgressBar"
                                    "{"
                                    "border: 2px solid grey;"
                                    "border-radius: 5px;"
                                    "text-align: center;"
                                    "}";

    myProgressBar = new QProgressBar(this);
    myProgressBar->setFixedWidth(300);
    myProgressBar->setStyleSheet(progressBarStyleSheet);

    //! -------------------
    //! the additional bar
    //! -------------------
    myProgressBar1 = new QProgressBar(this);
    myProgressBar1->setFixedWidth(300);
    myProgressBar1->setStyleSheet(progressBarStyleSheet);
    if(!additionalBar) myProgressBar1->hide();

    //! ------------------------------
    //! label for displaying messages
    //! ------------------------------
    myLabel = new QLabel(this);
    myLabel->setAlignment(Qt::AlignCenter);

    //! --------------
    //! "Stop" button
    //! --------------
    myStopButton = new QPushButton(this);
    myStopButton->setFixedWidth(100);
    myStopButton->setText("Stop");

    //! ----------------
    //! vertical layout
    //! ----------------
    QVBoxLayout *vl = new QVBoxLayout();
    vl->addWidget(myLabel);
    vl->addWidget(myProgressBar);
    if(myProgressBar1!=Q_NULLPTR) vl->addWidget(myProgressBar1);

    //! ------------------------------
    //! horizontal layout for buttons
    //! ------------------------------
    QHBoxLayout *hl = new QHBoxLayout();
    hl->addWidget(myStopButton);

    vl->addLayout(hl);
    this->setLayout(vl);

    if(!myStopButton->isEnabled()) myStopButton->setEnabled(true);
    connect(myStopButton,SIGNAL(pressed()),this,SLOT(handleStopPressed()));

    //! fix the window size // ???
    this->layout()->setSizeConstraint( QLayout::SetFixedSize );

    //! -------------------------------------------------------
    //! disable the "X" (close) and "?" (context help) buttons
    //! -------------------------------------------------------
    setWindowFlags(((this->windowFlags() | Qt::CustomizeWindowHint) & ~Qt::WindowCloseButtonHint) & ~Qt::WindowContextHelpButtonHint);    
}

//! --------------------------------------
//! function: setSecondatybarVisible
//! details:  show/hide the secondary bar
//! --------------------------------------
void QProgressIndicator::setSecondaryBarVisible(bool isVisible)
{
    cout<<"QProgressIndicator::setSecondaryBarVisible()->____function called____"<<endl;
    if(myProgressBar1!=Q_NULLPTR)
    {
        myProgressBar1->setVisible(isVisible);

        //! ------------------------------------------------------------------------
        //! from Qt documentation
        //!
        //! Updates the widget unless updates are disabled or the widget is hidden.
        //! This function does not cause an immediate repaint; instead it schedules
        //! a paint event for processing when Qt returns to the main event loop.
        //! This permits Qt to optimize for more speed and less flicker than a call
        //! to repaint() does. Calling update() several times normally results in
        //! just one paintEvent() call. Qt normally erases the widget's area before
        //! the paintEvent() call. If the Qt::WA_OpaquePaintEvent widget attribute
        //! is set, the widget is responsible for painting all its pixels with an
        //! opaque color.
        //! ------------------------------------------------------------------------
        this->update();
    }
}

//! -------------------------------
//! function: setPrimarybarVisible
//! details:
//! -------------------------------
void QProgressIndicator::setPrimaryBarVisible(bool isVisible)
{
    if(myProgressBar!=Q_NULLPTR)
    {
        myProgressBar->setVisible(isVisible);
        this->update();
    }
}

//! -----------------------
//! function: event filter
//! details:
//! -----------------------
bool QProgressIndicator::eventFilter(QObject *object, QEvent *event)
{
    if(event->type()==QProgressEvent::type())
    {
        QProgressEvent *e = static_cast<QProgressEvent*>(event);

        //! -----------------------------------------------------
        //! the current caller object
        //! at the moment used for stopping Tetgen sengin kill()
        //! -----------------------------------------------------
        if(e->getObject()!=NULL)
        {
            myCurCallerObject = e->getObject();
            //const QString &callerName =  e->getObject()->objectName();
            //cout<<"QProgressIndicator::eventFilter()->____caller: \""<<callerName.toStdString()<<"\"____"<<endl;
            //exit(9999);
        }

        //! -----------------------------------------
        //! update the current running task variable
        //! -----------------------------------------
        myCurrentRunningTask = e->getTask();

        //! -------------------------------
        //! handling the main progress bar
        //! -------------------------------
        switch(e->action())
        {
        case QProgressEvent_None:
        {
            ;   //! do nothing
        }
            break;

        case QProgressEvent_Init:
        {
            //! -------------------------------
            //! show the widget if not visible
            //! -------------------------------
            if(!this->isVisible()) this->setVisible(true);
            if(this->isHidden()) this->show();

            //! --------------------------------------------
            //! enable the "Stop" button if it was disabled
            //! --------------------------------------------
            if(!myStopButton->isEnabled()) myStopButton->setEnabled(true);

            myLabel->setText(e->getMessage());

            //! ----------------------------------------------
            //! set the range and set the progress to the min
            //! ----------------------------------------------
            myProgressBar->setRange(e->min(),e->max());
            myProgressBar->setValue(e->min());
        }
            break;

        case QProgressEvent_Update:
        {
            //! -------------------------------
            //! show the widget if not visible
            //! -------------------------------
            if(!this->isVisible()) this->setVisible(true);

            //! --------------------
            //! update the progress
            //! --------------------
            myProgressBar->setValue(e->val());
            myLabel->setText(e->getMessage());
        }
            break;

        case QProgressEvent_Reset:
        {
            //! -----------------
            //! reset to default
            //! -----------------
            myProgressBar->setRange(0,100);
            myProgressBar->setValue(0);
            myLabel->setText("");

            //! ----------------
            //! hide the widget
            //! ----------------
            if(!this->isHidden()) this->hide();
        }
            break;

        case QProgressEvent_DisableStop:
        {
            if(myStopButton->isEnabled()) myStopButton->setEnabled(false);
        }
            break;

        case QProgressEvent_EnableStop:
        {
            if(!myStopButton->isEnabled()) myStopButton->setEnabled(true);
        }
            break;
        }

        //! ------------------------------------
        //! handling the secondary progress bar
        //! ------------------------------------
        if(myProgressBar1!=Q_NULLPTR)
        {
            switch(e->action1())
            {
            case QProgressEvent_None:
            {
                ; //! do nothing
            }
                break;

            case QProgressEvent_Init:
            {
                //! -------------------------------
                //! show the widget if not visible
                //! -------------------------------
                if(!this->isVisible()) this->setVisible(true);

                //! ----------------------------------------------
                //! set the range and put the progress to the min
                //! ----------------------------------------------
                myProgressBar1->setRange(e->min1(),e->max1());
                myProgressBar1->setValue(e->min1());
            }
                break;

            case QProgressEvent_Update:
            {
                //! -------------------------------
                //! show the widget if not visible
                //! -------------------------------
                if(!this->isVisible()) this->setVisible(true);

                //! --------------------
                //! update the progress
                //! --------------------
                myProgressBar1->setValue(e->val1());
                myLabel->setText(e->getMessage());
            }
                break;

            case QProgressEvent_Reset:
            {
                //! -------------------
                //! reset to a default
                //! -------------------
                myProgressBar1->setRange(0,100);
                myProgressBar1->setValue(0);
            }
                break;

            case QProgressEvent_DisableStop:
            {
                if(myStopButton->isEnabled()) myStopButton->setEnabled(false);
            }
                break;

            case QProgressEvent_EnableStop:
            {
                if(!myStopButton->isEnabled()) myStopButton->setEnabled(true);
            }
                break;
            }
        }
        //! ---------------------
        //! the event stops here
        //! ---------------------
        return true;
    }
    else
    {
        return QObject::eventFilter(object, event);
    }
}

//! ----------------------------
//! function: handleStopPressed
//! details:
//! ----------------------------
void QProgressIndicator::handleStopPressed()
{
    cout<<"@-----------------------------------------@"<<endl;
    cout<<"@-         Stop has been pressed         -@"<<endl;
    cout<<"@-----------------------------------------@"<<endl;

    if(myCurrentRunningTask == "Netgen meshing")
    {
        Global::status().code = 0;

        //! -----------------------------
        //! this will stop the nglib.dll
        //! -----------------------------
        emit abortNetgenSTLPressed();
    }
    if(myCurrentRunningTask == "Tetgen meshing running on disk")
    {
        Global::status().code = 0;

        //! ------------------------------------------------------------------------------------
        //! interrupt tetgen
        //! From Qt documentation: "Kills the current process, causing it to exit immediately.
        //! On Windows, kill() uses TerminateProcess, and on Unix and macOS, the SIGKILL signal
        //! is sent to the process."
        //! ------------------------------------------------------------------------------------
        static_cast<TetgenMesher*>(myCurCallerObject)->getProcess()->kill();
    }
    if(myCurrentRunningTask == "Netgen rebuilding associativity" ||
            myCurrentRunningTask == "Tetgen building face mesh datasources from disk" ||
            myCurrentRunningTask == "MeshTools building PLC on disk" ||
            myCurrentRunningTask == "Creating automatic connections" ||
            myCurrentRunningTask == "Building prismatic mesh" ||
            myCurrentRunningTask == "TetHex conversion" ||
            myCurrentRunningTask == "Projecting mesh points on geometry" ||
            myCurrentRunningTask == "Writing CCX input file")
    {
        Global::status().code = 0;
    }
    if(myCurrentRunningTask == "Running CCX")
    {
        Global::status().code = 0;
        emit requestStopCCX();
    }

    QThread::msleep(125);
    this->hide();
}

//! ------------------
//! function: getText
//! details:
//! ------------------
QString QProgressIndicator::getText()
{
    return myLabel->text();
}

//! ------------------
//! function: setText
//! details:
//! ------------------
void QProgressIndicator::setText(const QString &text)
{
    myLabel->setText(text);
}

//! --------------------------------
//! function: setAbortButtonVisible
//! details:
//! --------------------------------
void QProgressIndicator::setAbortButtonVisible(bool isVisible)
{
    if(isVisible) myStopButton->setVisible(false);
    else myStopButton->setVisible(true);
}

//! -----------------------------
//! function: getCurrentProgress
//! details:
//! -----------------------------
int QProgressIndicator::getCurrentProgress()
{
    return myProgressBar1->value();
}

//! ------------------------------
//! function: getCurrentProgress1
//! details:
//! ------------------------------
int QProgressIndicator::getCurrentProgress1()
{
    return myProgressBar1->value();
}

//! ----------------------
//! function: enableAbort
//! details:
//! ----------------------
void QProgressIndicator::enableStop()
{
    if(!myStopButton->isEnabled()) myStopButton->setEnabled(true);
}

//! -----------------------
//! function: disableAbort
//! details:
//! -----------------------
void QProgressIndicator::disableStop()
{
    if(myStopButton->isEnabled()) myStopButton->setEnabled(false);
}

//! ------------------
//! function: setTask
//! details:
//! ------------------
void QProgressIndicator::setTask(const QString &task)
{
    myCurrentRunningTask = task;
}
