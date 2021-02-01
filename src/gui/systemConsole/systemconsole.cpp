//! custom includes
#include "systemconsole.h"
#include "qconsoleevent.h"
#include "src/utils/tools.h"
#include "src/utils/ccout.h"
//#include "ccxconsoletofile.h"
#include "qprogressevent.h"

//! Qt
#include <QVBoxLayout>
#include <QMenuBar>
#include <QMenu>
#include <QFileDialog>
#include <QAction>
#include <QEvent>
#include <QApplication>

//! C++
#include <fstream>
#include <iostream>

using namespace std;

//! ---------------------------------------
//! function: eventFilter
//! details:
//! ---------------------------------------
bool systemConsole::eventFilter(QObject* object, QEvent* event)
{
    if(object == this && event->type() == QConsoleEvent::type())
    {
        QString text = static_cast<QConsoleEvent*>(event)->getMessage();
        if(this->isRunning)
        {
            //! first append the text
            myTextEdit->appendPlainText(text);
        }
        if(this->myContinuoslyLogging==true)
        {
            *myOutFileStream<<text.toStdString()<<endl;
        }
        //! the event is stopped here
        return true;
    }
    else
    {
        // standard event processing
        return QObject::eventFilter(object, event);
    }
}

//! ---------------------------------------
//! function: constructor
//! details:
//! ---------------------------------------
systemConsole::systemConsole(bool aIsCreatedAsRunning,
                             bool isMenuBarVisible,
                             bool isContinuosLog,
                             const QString &logFile, QWidget *parent):QWidget(parent)
{
    //! init private members
    myIsMenuBarVisible = isMenuBarVisible;
    myLogFile = logFile;
    myContinuoslyLogging = isContinuosLog;

    //! open a log file
    if(myContinuoslyLogging==true)
    {
        myOutFileStream = new ofstream;
        myOutFileStream->open(myLogFile.toStdString(),ios::ate);
    }
    isCreatedAsRunning = aIsCreatedAsRunning;

    this->createContent();
    this->installEventFilter(this);

    if(aIsCreatedAsRunning==true) isRunning=true;
    else isRunning = false;
}

//! --------------------------------
//! function: setContinuoslyLogging
//! details:
//! ---------------------------------
void systemConsole::setContinuoslyLogging(bool isOn, const QString &logFile)
{
    myContinuoslyLogging = isOn;
    if(isOn)
    {
        myLogFile = logFile;
        myOutFileStream = new ofstream;
        myOutFileStream->open(myLogFile.toStdString(),ios::ate);
    }
    else
    {
        myLogFile = QString();
        if(myOutFileStream->is_open()) myOutFileStream->close();
        if(myOutFileStream!=NULL) delete myOutFileStream;
    }
}

//! ---------------------
//! function: destructor
//! details:
//! ---------------------
systemConsole::~systemConsole()
{
    cout<<"systemConsole::~systemConsole()->____DESTRUCTOR CALLED____"<<endl;
    if(myOutFileStream->is_open()) myOutFileStream->close();
    emit systemConsoleClosed();
}

//! ------------------------
//! function: createContent
//! details:
//! ------------------------
void systemConsole::createContent()
{
    //! -------------------
    //! create the content
    //! -------------------
    vl = new QVBoxLayout(this);
    vl->setMargin(0);
    vl->setContentsMargins(0,0,0,0);
    vl->setSpacing(0);

    myTextEdit = new QPlainTextEdit;

    //! --------------------
    //! context menu policy
    //! --------------------
    myTextEdit->setContextMenuPolicy(Qt::CustomContextMenu);
    connect(myTextEdit,SIGNAL(customContextMenuRequested(QPoint)),this,SLOT(showContextMenu(QPoint)));

    menuBar = new QMenuBar(this);

    //! ------------
    //! "File" menu
    //! ------------
    QMenu *menuFile = new QMenu("File",this);
    actionStart = new QAction("Start",this);
    actionStart->setIcon(QIcon(":/icons/icon_solve.png"));
    connect(actionStart,SIGNAL(triggered(bool)),this,SLOT(start()));

    actionSave = new QAction("Save",this);
    actionSave->setIcon(QIcon(":/icons/icon_save.png"));
    connect(actionSave,SIGNAL(triggered(bool)),this,SLOT(save()));

    actionSaveAs = new QAction("Save as",this);
    actionSaveAs->setIcon(QIcon(":/icons/icon_save as.png"));
    connect(actionSaveAs,SIGNAL(triggered(bool)),this,SLOT(saveAs()));

    actionStop = new QAction("Stop",this);
    actionStop->setIcon(QIcon(":/icons/icon_pause.png"));
    connect(actionStop,SIGNAL(triggered(bool)),this,SLOT(stop()));

    //actionExit = new QAction("Exit",this);
    //actionExit->setIcon(QIcon(":/icons/icon_exit.png"));
    //connect(actionExit,SIGNAL(triggered(bool)),this,SLOT(close()));

    //! -------------------
    //! fill the main menu
    //! -------------------
    menuFile->addAction(actionStart);
    menuFile->addAction(actionStop);

    if(isCreatedAsRunning==false)
    {
        actionStart->setEnabled(true);
        actionStop->setEnabled(false);
    }
    else
    {
        actionStart->setEnabled(false);
        actionStop->setEnabled(true);
    }

    //! ------------
    //! "Edit" menu
    //! ------------
    QMenu *menuEdit = new QMenu("Edit",this);
    actionClear = new QAction("Clear",this);
    actionClear->setIcon(QIcon(":/icons/icon_clear.png"));
    connect(actionClear,SIGNAL(triggered(bool)),this,SLOT(clear()));
    menuEdit->addAction(actionSave);
    menuEdit->addAction(actionSaveAs);
    menuEdit->addAction(actionClear);

    menuBar->addMenu(menuFile);
    menuBar->addMenu(menuEdit);

    vl->addWidget(menuBar);
    vl->addWidget(myTextEdit);

    if(myIsMenuBarVisible==true) menuBar->show();
    else menuBar->hide();

    //! -------------------
    //! set a default font
    //! -------------------
    myFont.setStyleHint(QFont::Helvetica);
    myFont.setPointSize(10);
    myTextEdit->setFont(myFont);
}

//! ---------------
//! function: stop
//! details:
//! ---------------
void systemConsole::stop()
{
    isRunning = false;

    actionStart->setEnabled(true);
    actionStop->setEnabled(false);

    emit systemConsoleStopped();
}

//! ----------------
//! function: start
//! details:
//! ----------------
void systemConsole::start()
{
    isRunning = true;

    actionStart->setEnabled(false);
    actionStop->setEnabled(true);
    emit systemConsoleStarted();
}

//! ---------------
//! function: save
//! details:
//! ---------------
void systemConsole::save()
{
    myTextEdit->appendPlainText("systemConsole::save->____function called____");
    QString workDir = tools::getWorkingDir();
    cout<<workDir.toStdString()<<endl;
    if(!workDir.isNull() && !workDir.isEmpty())
    {
        QString fileName = workDir.append("/default_log.txt");
        myTextEdit->appendPlainText(QString("systemConsole::save->____").append(fileName).append(QString("____")));
        ccout(QString("systemConsole::save()->____").append(fileName).append("____"));
        ofstream file;
        file.open(fileName.toStdString());
        if(file.is_open())
        {
            QString text = myTextEdit->toPlainText();
            QList<QString> textLine = text.split("/n");
            for(QList<QString>::iterator it = textLine.begin(); it!=textLine.end(); ++it)
            {
                QString aLine = *it;
                file<<aLine.toStdString()<<endl;
            }
            file.close();
        }
    }
}

//! -----------------
//! function: saveAs
//! details:
//! -----------------
void systemConsole::saveAs()
{
    QString filter;
    QString startDir;
    QString fileName = QFileDialog::getSaveFileName(this,"Save as",startDir,".txt",&filter,0);
    if(!fileName.isEmpty() && !fileName.isNull())
    {
        ofstream file;
        file.open(fileName.append(filter).toStdString().c_str());
        if(file.is_open())
        {
            QString text = myTextEdit->toPlainText();
            QList<QString> textLine = text.split("/n");
            for(QList<QString>::iterator it = textLine.begin(); it!=textLine.end(); ++it)
            {
                QString aLine = *it;
                file<<aLine.toStdString()<<endl;
            }
        }
    }
}

//! ----------------
//! function: clear
//! details:
//! ----------------
void systemConsole::clear()
{
    myTextEdit->clear();
}

//! --------------------------
//! function: showContextMenu
//! details:
//! --------------------------
void systemConsole::showContextMenu(QPoint p)
{
    QPoint pos = this->mapToGlobal(p);

    //! ----------
    //! main menu
    //! ----------
    QMenu* ctxMenu = new QMenu(this);

    //! -------
    //! "File"
    //! -------
    //QMenu *subMenu = new QMenu("File",this);

    ctxMenu->addAction(actionStart);
    ctxMenu->addAction(actionStop);
    //ctxMenu->addAction(actionExit);
    ctxMenu->addSeparator();

    //! -------
    //! "Edit"
    //! -------
    //QMenu *subMenu1 = new QMenu("Edit",this);

    ctxMenu->addSeparator();
    ctxMenu->addAction(actionSave);
    ctxMenu->addAction(actionSaveAs);
    ctxMenu->addAction(actionClear);

    //ctxMenu->addMenu(subMenu);
    //ctxMenu->addMenu(subMenu1);

    ctxMenu->addSeparator();

    //! ---------------------------------------
    //! "show/hide" menubar
    //! ---------------------------------------
    QAction *actionShowMenuBar = new QAction("Show menu bar",this);
    actionShowMenuBar->setIcon(QIcon(":/icons/icon_show menu.png"));
    connect(actionShowMenuBar,SIGNAL(triggered(bool)),this,SLOT(setMenuVisible()));
    if(myIsMenuBarVisible==false) ctxMenu->addAction(actionShowMenuBar);

    QAction *actionHideMenuBar = new QAction("Hide menu bar",this);
    actionHideMenuBar->setIcon(QIcon(":/icons/icon_hide menu.png"));
    connect(actionHideMenuBar,SIGNAL(triggered(bool)),this,SLOT(setMenuHidden()));
    if(myIsMenuBarVisible==true) ctxMenu->addAction(actionHideMenuBar);

    QAction* selectedItem = ctxMenu->exec(pos);
    if(selectedItem)
    {
        ;
    }
}

//! ------------------
//! function: setText
//! details:
//! ------------------
void systemConsole::setText(const QString &text)
{
    myTextEdit->setPlainText(text);
}
