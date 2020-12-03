//! -------------------------------------------------
//! to do list
//! 1) progress bar
//! 2) hide abort button from the progress indicator
//! -------------------------------------------------

#define A_SMALL_PROGRESS_TIMEOUT 100
//#define USE_FACE_CENTERS

//! ----------------
//! custom includes
//! ----------------
#include "openfoamreader.h"
#include <elementtypes.h>
#include <polygon_triangulate.h>
#include "tools.h"
#include "qprogressevent.h"
#include <qprogressindicator.h>
#include "geomtoolsclass.h"

//! ----
//! OCC
//! ----
#include <TColStd_PackedMapOfInteger.hxx>
#include <TColStd_MapIteratorOfPackedMapOfInteger.hxx>

//! ---
//! Qt
//! ---
#include <QDirIterator>
#include <QMessageBox>
#include <QLabel>
#include <QPushButton>
#include <QApplication>
#include <QTimer>
#include <QTreeView>

//! ----
//! C++
//! ----
#include <iostream>

//! ----------------
//! system specific
//! ----------------
#include <Windows.h>

using namespace std;
const double TINY = 1e-20;
const double TINY1 = 1e-12;
const double SCALE = 1000.0;

//! ---------------------
//! function: destructor
//! details:
//! ---------------------
OpenFoamReader::~OpenFoamReader()
{
    cout<<"OpenFoamReader::~OpenFoamReader()->____DESTRUCTOR CALLED____"<<endl;
    this->freeMemory();
}

//! ----------------------
//! function: deallocate
//! details:  free memory
//! ----------------------
void OpenFoamReader::freeMemory()
{
    cout<<"OpenFoamReader::freeMemory()->____function called____"<<endl;

    //! ----------------------
    //! delete the point list
    //! ----------------------
    if(myPointList.size()>0)
    {
        //cout<<"OpenFoamReader::freeMemory()->____deleting points____"<<endl;
        for(std::vector<double*>::iterator it = myPointList.begin(); it != myPointList.end(); ++it)
        {
            double *point = *it;
            delete [] point;
        }
        myPointList.clear();
    }

    //! ------------------------------------
    //! delete the map of rotation matrices
    //! ------------------------------------
    if(!myMapOfRotationMatrices.isEmpty())
    {
        //cout<<"OpenFoamReader::freeMemory()->____deleting the map of rotation matrices____"<<endl;
        for(QMap<int,double*>::iterator it = myMapOfRotationMatrices.begin(); it != myMapOfRotationMatrices.end(); ++it)
        {
            double *rotationMatrix = it.value();
            delete [] rotationMatrix;
        }
        myMapOfRotationMatrices.clear();
    }

    //! -------------------------------
    //! delete the map of face centers
    //! -------------------------------
    if(!myMapOfFaceCenters.isEmpty())
    {
        //cout<<"OpenFoamReader::freeMemory()->____deleting the map of face centers____"<<endl;
        for(QMap<int,double*>::iterator it = myMapOfFaceCenters.begin(); it != myMapOfFaceCenters.end(); ++it)
        {
            double *center = it.value();
            delete [] center;
        }
        myMapOfFaceCenters.clear();
    }

    /*
    //! -----------------------------------
    //! delete the map of the cell centers
    //! -----------------------------------
    if(!myMapOfCellCenters.isEmpty())
    {
        cout<<"OpenFoamReader::freeMemory()->____deleting the map of cell centers____"<<endl;
        for(QMap<int, double*>::iterator it = myMapOfCellCenters.begin(); it != myMapOfCellCenters.end(); ++it)
        {
            double *cellCenter = it.value();
            delete [] cellCenter;
        }
        myMapOfCellCenters.clear();
    }
    */
    cout<<"OpenFoamReader::freeMemory()->____free memory done____"<<endl;
}

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
OpenFoamReader::OpenFoamReader(/*const QString &sourceDirPath, const QString &targetDirPath,*/ int fileMode, QObject *parent):QObject(parent),myFileMode(fileMode)
{
    //cout<<"OpenFoamReader::OpenFoamReader()->____constructor called. Source file: "<<sourceDirPath.toStdString()<<"____"<<endl;
    //mySourceDir = sourceDirPath;
    //myTargetDir = targetDirPath;
    this->setObjectName("openFoamReader");

    //! ---------------------
    //! a timer - diagnostic
    //! ---------------------
    QTimer *timer = new QTimer(this);
    timer->setInterval(500);
    connect(timer,SIGNAL(timeout()),this,SLOT(updateTimer()));
    timer->start();
}

//! -----------------------------------------------------------------
//! function: updateTimer
//! details:  a message for checking if the thread is running or not
//! -----------------------------------------------------------------
void OpenFoamReader::updateTimer()
{
    static int ct = 0;
    ct +=500;
    cout<<"openFoamReader::updateTimer()->____Elapsed time: "<<ct<<"____"<<endl;
}

//! -------------------------------
//! function: setProgressIndicator
//! details:
//! -------------------------------
void OpenFoamReader::setProgressIndicator(QProgressIndicator *aProgressIndicator)
{
    myProgressIndicator = aProgressIndicator;
}
/*
//! ---------------------------------------------
//! function: setTargetDir
//! details:  set the location for storing files
//! ---------------------------------------------
void OpenFoamReader::setTargetDir(const QString &theTargetDir)
{
    cout<<"OpenFoamReader::setTargetDir()->____setting the target dir: "<<theTargetDir.toStdString()<<"____"<<endl;
    myTargetDir = theTargetDir;
    QDir dir(theTargetDir);
    if(dir.exists())
    {
        cout<<"OpenFoamReader::setTargetDir()->____the target directory exists: deleting content____"<<endl;
        tools::clearDir(myTargetDir);
    }
    else
    {
        cout<<"OpenFoamReader::setTargetDir()->____the target directory does not exist: creating a new one____"<<endl;
        dir.mkdir(myTargetDir);
    }
}
//! ---------------------------------------------------
//! function: setSourceDir
//! details:  set the location where reading the files
//! ---------------------------------------------------
void OpenFoamReader::setSourceDir(const QString &theSourceDir)
{
    cout<<"OpenFoamReader::setSourceDir()->____setting the source dir: "<<theSourceDir.toStdString()<<"____"<<endl;
    mySourceDir = theSourceDir;
}
*/
//! ------------------
//! function: perform
//! details:
//! ------------------
bool OpenFoamReader::perform(SimulationNodeClass *OFnode)
{
    cout<<"OpenFoamReader::perform()->____function called____"<<endl;

    mySourceDir = OFnode->getPropertyValue<QString>("Source directory");
    myTargetDir = OFnode->getPropertyValue<QString>("Target directory");
    QDir dir(myTargetDir);
    if(dir.exists())
    {
        cout<<"OpenFoamReader::setTargetDir()->____the target directory exists: deleting content____"<<endl;
        tools::clearDir(myTargetDir);
    }
    else
    {
        cout<<"OpenFoamReader::setTargetDir()->____the target directory does not exist: creating a new one____"<<endl;
        dir.mkdir(myTargetDir);
    }
    cout<<"OpenFoamReader::perform()->____"<<mySourceDir.toStdString()<<"____"<<endl;

#ifdef COSTAMP_VERSION
    //! retrieve the time list from the node
    const QVector<double> timeList = OFnode->getPropertyValue<QVector<double>>("Time list");
    myTimeFolders = timeList.toStdVector();
#endif
    QDir curDir(mySourceDir);

    cout<<(QString("OpenFoamReader::OpenFoamReader->____starting dir: ").append(curDir.absolutePath()).append("____")).toStdString()<<endl;
    cout<<(QString("OpenFoamReader::OpenFoamReader->____Results split mode %1____").arg(myFileMode)).toStdString()<<endl;
    cout<<(QString("OpenFoamReader::OpenFoamReader->____Starting directory ").append(curDir.absolutePath()).append("____")).toStdString()<<endl;
    //! ----------------------------------------------------
    //! [1] searching for directories containing the data
    //! entering the directory "0" for retrieving the names
    //! the directory "0" always exists
    //! ----------------------------------------------------
    if(!curDir.cd("0"))
    {
        cout<<"OpenFoamReader::OpenFoamReader()->____\"0\" directory not found____"<<endl;
        QMessageBox::warning(0,APPNAME,"OpenFoamReader: directory \"0\" is missing."
                                       "\nOpenfoam data archive is not valid", QMessageBox::Ok);
        return false;
    }
    //cout<<"OpenFoamReader::OpenFoamReader()->____Entering directory \"0\": retrieving the body names____"<<endl;
    QList<QFileInfo> entriesInfo = curDir.entryInfoList();
    QList<QString> directoryList = curDir.entryList();
    QList<QString> directoryListFiltered;

    for(int i=0; i<entriesInfo.length(); i++)
    {
        if(entriesInfo.at(i).isDir() && directoryList.at(i)!="."
                && directoryList.at(i)!=".." && directoryList.at(i)!="polyMesh"
                && directoryList.at(i)!="Casting0")
        {
            directoryListFiltered.append(directoryList.at(i));
            cout<<"OpenFoamReader::OpenFoamReader()->____found body: "<<directoryList.at(i).toStdString()<<"____"<<endl;
        }
    }

    //! ---------------------------------
    //! calculation steps for each body
    //! [1] reading boundary
    //! [2] reading points
    //! [3] reading faces
    //! [4] rebuilding mesh
    //! [5] building maps
    //! [6] translate data for each time
    //! ---------------------------------

    //! ---------------------
    //! a progress indicator
    //! ---------------------
    QProgressEvent *e;
    int NbBodies = directoryListFiltered.length();
    //cout<<"____tag01____"<<NbBodies<<endl;

    if(myProgressIndicator!=Q_NULLPTR)
    {
        e = new QProgressEvent(QProgressEvent_Init,0,NbBodies-1,0,"Start converting",QProgressEvent_Init,1,6,1);
        QApplication::postEvent(myProgressIndicator,e);
        QApplication::processEvents();
    }

    curDir.cdUp();
    curDir.cd("constant");
    //cout<<"____tag02____"<<endl;

    //! ------------------------------
    //! entering the data directories
    //! ------------------------------
    for(int i=0; i<NbBodies; i++)
    {
        //! ------------------------------------------------
        //! the directory name is also the name of the body
        //! ------------------------------------------------
        QString dirName = directoryListFiltered.at(i);

        curDir.cd(dirName);
        curDir.cd("polyMesh");
        //! ----------------------------------
        //! [1] read boundary fields settings
        //! ----------------------------------
        char fileName[512];
        sprintf(fileName,"%s",curDir.filePath("boundary").toStdString().c_str());

        //cout<<"OpenFoamReader::perform()->____Opening file: "<<fileName<<"____"<<endl;
        //cout<<"OpenFoamReader::perform()->____Entering directory: "<<dirName.toStdString()<<"____"<<endl;

        if(myProgressIndicator!=Q_NULLPTR)
        {
            QString msg = QString("Reading boundary of: ").append(QString::fromStdString(dirName.toStdString()));
            e = new QProgressEvent(QProgressEvent_Update,0,NbBodies-1,i,msg,QProgressEvent_Update,1,9999,1);
            QApplication::postEvent(myProgressIndicator,e);
            QApplication::processEvents();
        }

        Sleep(A_SMALL_PROGRESS_TIMEOUT);

        //cout<<"OpenFoamReader::perform()->____Reading boundary definitions____"<<endl;

        fstream is;
        openWithLock(is,fileName, std::ios_base::in);
        std::vector<BB> vecBoundaries;
        this->readBoundarySection(vecBoundaries, is);
        closeAndRemoveLock(is,fileName);

        //! ---------------------
        //! [2] read mesh points
        //! ---------------------
        if(myProgressIndicator!=Q_NULLPTR)
        {
            QString msg = QString("Reading mesh points of body: ").append(QString::fromStdString(dirName.toStdString()));
            e = new QProgressEvent(QProgressEvent_Update,0,NbBodies-1,i,msg,QProgressEvent_Update,1,9999,2);
            QApplication::postEvent(myProgressIndicator,e);
            QApplication::processEvents();
        }

        Sleep(A_SMALL_PROGRESS_TIMEOUT);

        sprintf(fileName,"%s",curDir.filePath("points").toStdString().c_str());
        //cout<<"OpenFoamReader::perform()->____Opening: "<<fileName<<"____"<<endl;

        //! ------------------------------
        //! uses the optimized overload
        //! "myPointList" is created here
        //! ------------------------------
        openWithLock(is,fileName,std::ios_base::in);
        this->readPoints(is);
        closeAndRemoveLock(is,fileName);

        //! ---------------
        //! [3] read faces
        //! ---------------
        if(myProgressIndicator!=Q_NULLPTR)
        {
            QString msg = QString("Reading faces of: ").append(QString::fromStdString(dirName.toStdString()));
            e = new QProgressEvent(QProgressEvent_Update,0,NbBodies-1,i,msg,QProgressEvent_Update,0,9999,3);
            QApplication::postEvent(myProgressIndicator,e);
            QApplication::processEvents();
        }

        Sleep(A_SMALL_PROGRESS_TIMEOUT);

        sprintf(fileName,"%s",curDir.filePath("faces").toStdString().c_str());
        //cout<<"OpenFoamReader::perform()->____Opening file: "<<fileName<<"____"<<endl;
        //cout<<"OpenFoamReader::perform()->____Reading faces____"<<endl;

        openWithLock(is,fileName, std::ios_base::in);

        //! -----------------------------------------------------------
        //! the definition of one face is through the list of node IDs
        //! -----------------------------------------------------------
        std::vector<std::vector<int>> listOfFaceDefinitions;

        //! -----------------------------------------
        //! uses the optimized overload
        //! here "myMapOfRotationMatrices" is filled
        //! and "myMapOfFaceCenters" is filled also
        //! -----------------------------------------
        this->readFaces(is, listOfFaceDefinitions);

        //! ---------------------------------------
        //! create an indexed map of all the faces
        //! ---------------------------------------
        TColStd_PackedMapOfInteger allFacesMap;

        for(int faceNr = 0; faceNr<listOfFaceDefinitions.size(); faceNr++) allFacesMap.Add(faceNr); //cese
        //for(int faceNr = 0; faceNr<=listOfFaceDefinitions.size(); faceNr++) allFacesMap.Add(faceNr);  //! versione originale

        closeAndRemoveLock(is,fileName);

        //! ----------------------------------------------
        //! read cells: calculate the center of the cells
        //! ----------------------------------------------
        if(myProgressIndicator!=Q_NULLPTR)
        {
            QString msg = QString("Calculating the centers of the cells for body: ").append(QString::fromStdString(dirName.toStdString()));
            e = new QProgressEvent(QProgressEvent_Update,0,NbBodies-1,i,msg,QProgressEvent_Update,1,9999,4);
            QApplication::postEvent(myProgressIndicator,e);
            QApplication::processEvents();
        }

        Sleep(A_SMALL_PROGRESS_TIMEOUT);

        //cout<<"OpenFoamReader::perform()->____Calculating the centers of the cells____"<<endl;

        sprintf(fileName,"%s",curDir.filePath("owner").toStdString().c_str());
        //cout<<"OpenFoamReader::perform()->____Opening: "<<fileName<<"____"<<endl;

        openWithLock(is,fileName, std::ios_base::in);

        //cout<<"OpenFoamReader::perform()->____Opening file: "<<fileName<<"____"<<endl;

        char fileName1[128];
        sprintf(fileName1,"%s",curDir.filePath("neighbour").toStdString().c_str());

        //cout<<"OpenFoamReader::perform()->____Opening: "<<fileName<<"____"<<endl;

        fstream is1;
        openWithLock(is1,fileName1, std::ios_base::in);

        //cout<<"OpenFoamReader::perform()->____Opening file: "<<fileName<<"____"<<endl;

        //! ----------------------------------------
        //! calculate the center of each polyhedron
        //! here "myCellCenters" is filled
        //! ----------------------------------------
        this->getCellCenters(is,is1,listOfFaceDefinitions);

        closeAndRemoveLock(is,fileName);
        closeAndRemoveLock(is1,fileName1);

        //! ---------------------------
        //! [5] building the face maps
        //! ---------------------------
        if(myProgressIndicator!=Q_NULLPTR)
        {
            QString msg = QString("Building the face maps for body: ").append(QString::fromStdString(dirName.toStdString()));
            e = new QProgressEvent(QProgressEvent_Update,0,NbBodies-1,i,msg,QProgressEvent_Update,1,9999,5);
            QApplication::postEvent(myProgressIndicator,e);
            QApplication::processEvents();
        }

        Sleep(A_SMALL_PROGRESS_TIMEOUT);

        //cout<<"OpenFoamReader::perform()->____Grouping the element faces____"<<endl;

        QMap<std::string, QList<int>> localMaps;
        for(std::vector<BB>::iterator it = vecBoundaries.begin(); it!= vecBoundaries.end(); ++it)
        {
            const BB &curBoundary = *it;

            char bName[512];
            sprintf(bName,"%s",curBoundary.boundaryName);
            int startF = curBoundary.startFace;
            int endF = curBoundary.nFaces+startF-1;

            TColStd_PackedMapOfInteger tempFacesMap;

            QList<int> localToGlobalFaceNr;
            for(int globalFaceNr = startF; globalFaceNr<=endF; globalFaceNr++)
            {
                localToGlobalFaceNr.append(globalFaceNr);
                tempFacesMap.Add(globalFaceNr);
            }
            QString tmp = QString::fromLatin1(bName);
            std::string boundaryName = tmp.toStdString();

            //! -----------------------------------------------------------------
            //! (name of the boundary, link between local and global face index)
            //! -----------------------------------------------------------------
            localMaps.insert(boundaryName, localToGlobalFaceNr);
            allFacesMap.Subtract(tempFacesMap);
        }

        //cout<<"OpenFoamReader::perform()->____number of internal faces: "<<allFacesMap.Extent()<<"____"<<endl;

        //! ---------------------------------------------------
        //! inserting into the map the group of internal faces
        //! ---------------------------------------------------
        QList<int> localToGlobalInternalFaceNr;
        for(TColStd_MapIteratorOfPackedMapOfInteger anIt(allFacesMap); anIt.More(); anIt.Next())
        {
            localToGlobalInternalFaceNr.append(anIt.Key());
        }
        localMaps.insert("internalField",localToGlobalInternalFaceNr);
        //cout<<"OpenFoamReader::perform()->____number of internal faces: "<<localToGlobalInternalFaceNr.size()<<"____"<<endl;

        //! ----------------
        //! position @ root
        //! ----------------
        curDir.cdUp();
        curDir.cdUp();
        curDir.cdUp();

        //! -----------------------------------------
        //! search for directories with a time label
        //! -----------------------------------------
        QList<QFileInfo> entriesInfo = curDir.entryInfoList();
        QList<QString> directoryList = curDir.entryList();
        QList<QString> directoryListFiltered;

        for(int i=0; i<entriesInfo.length(); i++)
        {
#ifdef COSTAMP_VERSION
        int nBstep = int(myTimeFolders.size());
        for(int j=0;j<nBstep;j++)
        {
            double endTime = myTimeFolders.at(j);
            QString dirTime = QString("%1").arg(endTime,0,'g',-1);
            //cout<<"time "<<endTime<<" dirTime "<<directoryList.at(i).toStdString()<<endl;
            if(entriesInfo.at(i).isDir() && directoryList.at(i) == dirTime)
            {
                directoryListFiltered.append(directoryList.at(i));
            }
        }
#endif
#ifndef COSTAMP_VERSION
            if(entriesInfo.at(i).isDir() && directoryList.at(i)!="." && directoryList.at(i)!=".."
                    && directoryList.at(i)!="constant" && directoryList.at(i)!="system")
            {
                directoryListFiltered.append(directoryList.at(i));
                //cout<<"____"<<directoryList.at(i).toStdString()<<"____"<<endl;
            }
#endif
        }

        //! ------------------------------
        //! [6] working on "time" folders
        //! ------------------------------
        for(int t=0; t<directoryListFiltered.length(); t++)
        {
            //! the "Time" directory (curDirName is ->relative<-)
            const QString &curDirName = directoryListFiltered.at(t);

            //cout<<"OpenFoamReader::perform()->____entering time directory: \""<<curDirName.toStdString()<<"\"____"<<endl;

            if(myProgressIndicator!=Q_NULLPTR)
            {
                QString msg = QString("Body: %1 Translating data at time %2").arg(dirName).arg(atof(curDirName.toLatin1()));
                e = new QProgressEvent(QProgressEvent_Update,0,NbBodies-1,i,msg,QProgressEvent_Update,1,9999,6);
                QApplication::postEvent(myProgressIndicator,e);
                QApplication::processEvents();
            }

            Sleep(A_SMALL_PROGRESS_TIMEOUT);

            //cout<<"OpenFoamReader::perform()->____Body: \""<<dirName.toStdString()<<"\" Translating data at time: "<<curDirName.toStdString()<<"____"<<endl;

            curDir.cd(curDirName);  // entering the time directory
            curDir.cd(dirName);     // entering the body directory

            //! ------------------------------------------------------
            //! open a file for writing results
            //! output data for the current body, at the current time
            //! ------------------------------------------------------
            fstream fout;
            fstream foutOverall;

            QString foutName = myTargetDir+QString(("/Body_%1_Time_%2.txt")).arg(dirName).arg(curDirName);
            QString foutNameOverall = myTargetDir+QString(("/Overall_%1.txt")).arg(curDirName);

            //cout<<"OpenFoamReader::perform()->____target directory: "<<myTargetDir.toStdString()<<"____"<<endl;

            if(myFileMode==1)
            {
                //cout<<"OpenFoamReader::perform()->____Creating file: \""<<foutName.toStdString()<<"\"____"<<endl;
                openWithLock(fout,foutName.toStdString().c_str(), std::ios_base::out);
            }
            else
            {
                //cout<<"OpenFoamReader::perform()->____Creating file: \""<<foutNameOverall.toStdString()<<"\"____"<<endl;
                openWithLock(foutOverall,foutNameOverall.toStdString().c_str(), std::ios_base::app);
            }

            //! -----------------------------------------------
            //! reading the scalar dist from the openfoam file
            //! "T" can be changed into another scalar
            //! -----------------------------------------------
            sprintf(fileName,"%s",curDir.absoluteFilePath("TCel").toStdString().c_str());
            //cout<<"OpenFoamReader::perform()->____Opening file: \""<<fileName<<"\"____"<<endl;
            is.open(fileName);

            //! ----------------
            //! skip the header
            //! ----------------
            std::string val;
            for(int n=0;n<=100;n++)
            {
                std::getline(is,val);
                if(strcmp("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //",val.c_str())==0) break;
            }

            std::getline(is,val);
            std::getline(is,val);
            std::getline(is,val);

            //! -------------------------------------------------------------------------
            //! handling "internalField"
            //! the value of the scalar is given at the center of the cell (by default?)
            //! or at the center of each internal face
            //! -------------------------------------------------------------------------
            std::getline(is,val);
            char type1[128], type2[128];
            sscanf(val.c_str(),"%s%s",type1,type2);

            if(strcmp("uniform",type2)==0)      //! case "uniform"
            {
                double scalarValue = 0;
                sscanf(val.c_str(),"%s%s%lf",type1,type2,&scalarValue);
                //cout<<"OpenFoamReader::perform()->____uniform scalar value: "<<scalarValue<<endl;
                //scalarValue+=273.15;

                //! ------------------------------
                //! this works on element centers
                //! ------------------------------
                //for(QMap<int,double*>::iterator it = myMapOfCellCenters.begin(); it!=myMapOfCellCenters.end(); ++it)
                for(QMap<int,std::vector<double>>::iterator it = myMapOfCellCenters.begin(); it!=myMapOfCellCenters.end(); ++it)
                {
                    //double *aCenter = it.value();
                    std::vector<double> aCenter = it.value();

                    char s[512];
                    sprintf(s,"%.9e\t%.9e\t%.9e\t%.9e\n",aCenter[0]*SCALE,aCenter[1]*SCALE,aCenter[2]*SCALE,scalarValue+273.15);

                    //! --------------
                    //! write on disk
                    //! --------------
                    if(myFileMode==1)
                    {
                        fout<<s;
                    }
                    else
                    {
                        foutOverall<<s;
                    }
                    //! or store in memory ...
                }

                /* this piece of code works on internal faces
                QList<int> internalFieldIndexes = localMaps.value("internalField");
                for(int i=0; i<internalFieldIndexes.length(); i++)
                {
                    int faceGlobalNumber = internalFieldIndexes.at(i);
                    QList<double> sourcePoint = map.value(faceGlobalNumber);
                    //! write on disk
                    if(myFileMode==1) fprintf(foutData,"%.9e\t%.9e\t%.9e\t%.9e\n",sourcePoint.at(0)*SCALE,sourcePoint.at(1)*SCALE,sourcePoint.at(2)*SCALE,value);
                    //! or store in memory ...
                }
                */
            }
            else    //! case "non uniform"
            {
                std::getline(is,val);
                int Ndata;
                sscanf(val.c_str(),"%d",&Ndata);
                //cout<<"OpenFoamReader::perform()->____Number of scalar values: "<<Ndata<<endl;

                std::getline(is,val);   //! capture "("
                char value[32];
                sscanf(val.c_str(),"%s",value);

                //! ------------------------------
                //! this works on element centers
                //! ------------------------------
                for(int row=0;;row++)
                {
                    double scalarValue = 0;

                    std::getline(is,val);
                    sscanf(val.c_str(),"%s",value);

                    //! check for termination of data: capture ")"
                    if(strcmp(value,")")==0)
                    {
                        std::getline(is,val);   //! then capture ";"
                        break;
                    }
                    sscanf(val.c_str(),"%lf",&scalarValue);
                    //cout<<"non uniform scalar value: "<<scalarValue<<endl;  //! diagnostic

                    //scalarValue+=273.15;

                    //double* aCenter = myMapOfCellCenters.value(row);
                    std::vector<double> aCenter = myMapOfCellCenters.value(row);

                    //! --------------
                    //! write on disk
                    //! --------------
                    char s[512];
                    sprintf(s,"%.9e\t%.9e\t%.9e\t%.9e\n",aCenter[0]*SCALE,aCenter[1]*SCALE,aCenter[2]*SCALE,scalarValue+273.15);
                    if(myFileMode==1)
                    {
                        fout<<s;
                    }
                    else
                    {
                        foutOverall<<s;
                    }
                    //! or store in memory ...
                }

                /* this works on internal faces
                //QList<int> internalFieldIndexes = localMaps.value("internalField");
                for(int row=0;;row++)
                {
                    double scalarValue=0;
                    std::getline(is,val);
                    sscanf(val.c_str(),"%s",value);
                    if(strcmp(value,")")==0)
                    {
                        std::getline(is,val);   //! capture ";"
                        //cout<<"------------------>"<<val.c_str()<<endl;
                        break;
                    }
                    sscanf(val.c_str(),"%lf",&scalarValue);
                    //cout<<"non uniform scalar value: "<<scalarValue<<endl;  //! diagnostic
                    //!int faceGlobalNumber = internalFieldIndexes.at(row);
                    //!QList<double> sourcePoint = map.value(faceGlobalNumber);
                    //! write on disk
                    if(myFileMode==1) fprintf(foutData,"%.9e\t%.9e\t%.9e\t%.9e\n",sourcePoint.at(0)*SCALE,sourcePoint.at(1)*SCALE,sourcePoint.at(2)*SCALE,scalarValue+273.15);
                    //! or store in memory ...
                }
                */
            }

#ifdef USE_FACE_CENTERS
            //! --------------------------
            //! handling "Boundary field"
            //! --------------------------
            std::getline(is,val);   //! capture an empty line
            std::getline(is,val);
            std::getline(is,val);   //! capture "{"
            char s[512];

            for(int i=0; i<vecBoundaries.size(); i++)
            {
                //! name
                char bName[512];
                const BB &bn = vecBoundaries.at(i);
                sprintf(bName,"%s",bn.boundaryName);
                //cout<<"OpenFoamReader::perform()->____boundary name: "<<bn.boundaryName<<"____"<<endl;
                std::string aKey(bName);
                //! end name

                std::getline(is,val);
                sscanf(val.c_str(),"%s",s);
                while(strcmp(s,"}")!=0)
                {
                    if(strcmp(s,"value")==0)
                    {
                        sscanf(val.c_str(),"%s%s",type1,type2);
                        if(strcmp(type2,"uniform")==0)
                        {
                            //! handling uniform as before
                            double value = 0;
                            sscanf(val.c_str(),"%s%s%lf",type1,type2,&value);
                            cout<<"OpenFoamReader::perform()->____uniform: value = "<<value<<____"<<endl;
                            //value=value+273.15;
                            QList<int> indexes = localMaps.value(aKey);
                            for(int i=0; i<indexes.length(); i++)
                            {
                                int faceGlobalNumber = indexes.at(i);
                                double *sourcePoint = myMapOfFaceCenters.value(faceGlobalNumber);
                                //! write on disk
                                char s[512];
                                sprintf(s,"%.9e\t%.9e\t%.9e\t%.9e\n",sourcePoint[0]*SCALE,sourcePoint[1]*SCALE,sourcePoint[2]*SCALE,value+273.15);
                                if(myFileMode==1)
                                {
                                    fout<<s;
                                }
                                else
                                {
                                    foutOverall<<s;
                                }
                                //! or store in memory ...
                            }
                        }
                        else
                        {
                            //!----------------------------------------------------------------
                            //! Handling nonuniform value per column
                            //! details:
                            //! value           nonuniform List<scalar> 2(90.293947 90.283879);
                            //!-----------------------------------------------------------------
                            char type3[128],type4[128];
                            cout<<"OpenFoamReader::perform()->____reading boundary block nonuniform scalarValue value = "<<val<<endl;
                            sscanf(val.c_str(),"%s%s%s%s",type1,type2,type3,type4);
                            cout<<"OpenFoamReader::perform()->____reading boundary block type4  = "<<type4<<endl;
                            //! -----------------------------------------
                            //! small record of data collected "in line"
                            //! jump over it
                            //! -----------------------------------------
                            if(type4==NULL) continue;
                            //! handling non-uniform as before
                            cout<<"OpenFoamReader::perform()->____type: non uniform____"<<endl;
                            std::getline(is,val);
                            std::getline(is,val);   //! caputure N data
                            cout<<"OpenFoamReader::perform()->____number of data: "<<val<<endl;
                            char s[512];
                            int n=0;
                            for(int row=0;;row++)
                            {
                                std::getline(is,val);   //! capture "("
                                sscanf(val.c_str(),"%s",s);
                                if(strcmp(s,")")==0)
                                {
                                    //cout<<"@row: "<<row<<endl;
                                    //cout<<val<<endl;
                                    break;
                                }
                                else
                                {
                                    double value=0;
                                    n++;
                                    sscanf(val.c_str(),"%lf",&value);
                                    //cout<<"-------->"<<value<<endl;
                                    //value+=273.15;
                                    QList<int> indexes = localMaps.value(aKey);
                                    if(n==1) cout<<indexes.length()<<endl;
                                    //cout<<"@row: "<<row<<endl;
                                    int faceGlobalNumber = indexes[row];   //! check carefully
                                    double *sourcePoint = myMapOfFaceCenters.value(faceGlobalNumber);
                                    //! --------------
                                    //! write on disk
                                    //! --------------
                                    char s[512];
                                    sprintf(s,"%.9e\t%.9e\t%.9e\t%.9e\n",sourcePoint[0]*SCALE,sourcePoint[1]*SCALE,sourcePoint[2]*SCALE,value+273.15);
                                    if(myFileMode==1)
                                    {
                                        fout<<s;
                                    }
                                    else
                                    {
                                        foutOverall<<s;
                                    }
                                    //! or store in memory
                                }
                            }
                            cout<<"OpenFoamReader::perform()->____"<<n<<" non uniform values written___"<<endl;
                        }
                    }
                    std::getline(is,val);
                    sscanf(val.c_str(),"%s",s);
                }
            }
#endif
            is.close();
            curDir.cdUp();

            //! --------------------------
            //! exiting the "time" folder
            //! --------------------------
            curDir.cdUp();
            if(myFileMode==1)
            {
                //cout<<"OpenFoamReader::perform()->____Closing file: "<<foutName.toStdString()<<"____"<<endl;
                closeAndRemoveLock(fout,foutName.toStdString().c_str());
            }
            else
            {
                //cout<<"OpenFoamReader::perform()->____Closing file: "<<foutNameOverall.toStdString()<<"____"<<endl;
                closeAndRemoveLock(foutOverall,foutNameOverall.toStdString().c_str());
            }
            //cout<<"OpenFoamReader::perform()->____cur directory after exiting time folder: "<<curDir.absolutePath().toStdString()<<"\"____"<<endl;
        }

        curDir.cd("constant");

        //! --------------------------------
        //! delete the maps and the vectors
        //! --------------------------------
        this->freeMemory();
    }

    //! ---------------------------
    //! resetting the progress bar
    //! ---------------------------
    if(myProgressIndicator!=Q_NULLPTR)
    {
        e = new QProgressEvent(QProgressEvent_Reset,0,9999,0,"",QProgressEvent_Reset,0,9999,0);
        QApplication::postEvent(myProgressIndicator,e);
        QApplication::processEvents();
    }

    //! -------------------------
    //! emit conversion finished
    //! -------------------------
    emit conversionFinished();
    return true;
}

//! -------------------------------
//! function: readBoundarySection
//! details:  helper
//! -------------------------------
void OpenFoamReader::readBoundarySection(std::vector<BB> &vecBoundaries, fstream &is)
{
    //! ----------------
    //! skip the header
    //! ----------------
    std::string val;
    for(int n=0;n<=100;n++)
    {
        std::getline(is,val);
        if(strcmp("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //",val.c_str())==0) break;
    }
    std::getline(is,val);

    //! ---------------------------
    //! reading "boundary" section
    //! ---------------------------
    std::getline(is,val);
    int Nb;
    sscanf(val.c_str(),"%d",&Nb);
    //cout<<"OpenFoamReader::readBoundarySection()->____number of boundaries: "<<Nb<<"____"<<endl;

    std::getline(is,val);   //! (

    for(int boundaryNumber=1; boundaryNumber<=Nb; boundaryNumber++)
    {
        BB curBoundary;
        int nFaces,startFace;

        std::getline(is,val);   //! read the <name> of the boundary
        char curName[512];
        sscanf(val.c_str(),"%s",curName);
        //cout<<"OpenFoamReader::readBoundarySection()->____name of the boundary: "<<curName<<"____"<<endl;

        QString domainCheckString = QString::fromLatin1(curName);
        sprintf(curBoundary.boundaryName,"%s",curName);
        while(strcmp(curName,"}")!=0)
        {
            std::getline(is,val);      //! capture {
            sscanf(val.c_str(),"%s",curName);
            if(strcmp(curName,"startFace")==0)
            {
                sscanf(val.c_str(),"%s%d",curName,&startFace);
                curBoundary.startFace = startFace;
                //cout<<"startFace field found: "<<startFace<<endl;
            }
            if(strcmp(curName,"nFaces")==0)
            {
                sscanf(val.c_str(),"%s%d",curName,&nFaces);
                curBoundary.nFaces = nFaces;
                //cout<<"nFaces field found: "<<nFaces<<endl;
            }
            //std::getline(is,val);      //! capture }
            //cout<<curName<<endl;
        }
        if(!domainCheckString.contains("_to_"))
        {
            vecBoundaries.push_back(curBoundary);
        }
        else
        {
            //cout<<"____skipping since defined as \"interface\"____"<<endl;
        }
    }
    std::getline(is,val);   //! ")"
}

//! -----------------------------
//! function: readPoints
//! details:  helper - optimized
//! -----------------------------
void OpenFoamReader::readPoints(fstream &is)
{
    //! ----------------
    //! skip the header
    //! ----------------
    std::string val;
    for(int n=0;n<=100;n++)
    {
        std::getline(is,val);
        if(strcmp("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //",val.c_str())==0)
            break;
    }

    //! ----------------
    //! capture Npoints
    //! ----------------
    int Nent;
    for(;;)
    {
        std::getline(is,val);
        if(sscanf(val.c_str(),"%d",&Nent)==1) break;
    }

    //! ------------
    //! capture "("
    //! ------------
    std::getline(is,val);

    double x,y,z;
    for(int i=1; i<=Nent; i++)
    {
        std::getline(is,val);
        sscanf(val.c_str(),"(%lf %lf %lf)",&x,&y,&z);
        double *point = new double[3];
        point[0] = x;
        point[1] = y;
        point[2] = z;
        myPointList.push_back(point);
    }
}

//! --------------------
//! function: readFaces
//! details:  optimized
//! --------------------
void OpenFoamReader::readFaces(fstream &is, std::vector<std::vector<int>> &listOfFaceDefinitions)
{
    //cout<<"____start reading faces: "<<this->clock()<<"____"<<endl;

    //! ip: inserted centers
    //! bad: failures
    int bad=0;
    std::string val;
    for(int n=0;n<=30;n++)
    {
        std::getline(is,val);
        if(strcmp("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //",val.c_str())==0)
            break;
    }

    //! -------------------------------
    //! capture Nent (number of faces)
    //! -------------------------------
    int Nent;
    for(;;)
    {
        std::getline(is,val);
        if(sscanf(val.c_str(),"%d",&Nent)==1) break;
    }

    //! ------------
    //! capture "("
    //! ------------
    std::getline(is,val);

    //! -----------------------------------------
    //! int => position within the list
    //! std::vector<double> => point coordinates
    //! -----------------------------------------
    char c;
    int Ncomp;

    //! -----------------------
    //! iterate over the faces
    //! -----------------------
    for(int faceNr=0; faceNr<Nent; faceNr++)
    {
        //if(faceNr%1000==0)//cout<<"____reading face nr: "<<faceNr<<" from file____"<<endl;

        //! -----------------------------------
        //! a list of node IDs defining a face
        //! -----------------------------------
        std::vector<int> nodeIDs;

        is>>Ncomp;
        is>>c;
        for(int k=1; k<=Ncomp; k++)
        {
            int nodeID;
            is>>nodeID;
            nodeIDs.push_back(nodeID);
        }
        listOfFaceDefinitions.push_back(nodeIDs);

        is>>c;

        //! ------------------
        //! working on a face
        //! ------------------
        double *P1 = myPointList[nodeIDs[0]];
        double *P2 = myPointList[nodeIDs[1]];
        double *P3;

        //cout<<"____P1("<<P1[0]<<", "<<P1[1]<<", "<<P1[2]<<")____"<<endl;
        //cout<<"____P2("<<P2[0]<<", "<<P2[1]<<", "<<P2[2]<<")____"<<endl;

        //! ---------------------------
        //! three points
        //! {1,0,2} {2,1,3} {3,2,4} ...
        //! ---------------------------
        bool flag = false;
        for(int h=2; h<nodeIDs.size(); h++)
        {
            P3 = myPointList.at(nodeIDs[h]);
            //! -----------------------------
            //! uses the "optimized" version
            //! -----------------------------
            bool isCollinear = GeomToolsClass::isCollinear(P1,P2,P3);

            if(isCollinear==false)
            {
                flag = true;
                break;
            }
            else
            {
                flag = false;
            }
        }
        if(flag==false)
        {
            //cout<<"____bad face: all points are aligned____"<<endl;
            bad++;
            continue;
        }

        //! ------------
        //! three nodes
        //! ------------
        double *rotationMatrix = new double[9];
        bool isDone = GeomToolsClass::getRotationMatrix(P1,P2,P3,rotationMatrix);

        if(isDone)
        {
            //! ------------------------------------------------------
            //! this is built but not used here: it will be used when
            //! calculating the cell centers
            //! ------------------------------------------------------
            myMapOfRotationMatrices.insert(faceNr,rotationMatrix);

#ifdef USE_FACE_CENTERS
            //! ---------------------------
            //! get the center of the face
            //! ---------------------------
            double *C = new double[3];
            double faceArea;
            bool isFaceOK;
            isFaceOK = GeomToolsClass::getFaceCenter(myPointList, rotationMatrix, C, faceArea);
            if(isFaceOK)
            {
                myMapOfFaceCenters.insert(faceNr,C);
            }
            else
            {
                bad++;
                continue;
            }
#endif
        }
        else
        {
            bad++;
            continue;
        }
    }

    //cout<<"____Number bad faces: "<<bad<<"____"<<endl;
    //double failures = double(bad)/double(Nent);
    //cout<<"____Failures ratio: "<<failures<<endl;
    //cout<<"____end reading faces: "<<this->clock()<<"____"<<endl;
}

//! -------------------------
//! function: getCellCenters
//! details:
//! -------------------------
void OpenFoamReader::getCellCenters(fstream &is,
                                    fstream &is1,
                                    const std::vector<std::vector<int>> &listOfFaceDefinitions)
{
    //cout<<"OpenFoamReader::getCellCenters()->____function called: rebuilding elements____"<<endl;
    //! ------------------
    //! read "owner" file
    //! ------------------

    //! ----------------
    //! skip the header
    //! ----------------
    std::string val;
    for(int n=0;n<=100;n++)
    {
        std::getline(is,val);
        if(strcmp("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //",val.c_str())==0) break;
    }
    int NE;
    for(;;)
    {
        std::getline(is,val);
        if(sscanf(val.c_str(),"%d",&NE)==1) break;
    }
    //cout<<"OpenFoamReader::getCellCenters()->____Number of owner entries: "<<NE<<"____"<<endl;

    //! capture "("
    std::getline(is,val);

    QMultiMap<int,int> elementsByFaces;

    //! capture first number
    int n;
    int pos;

    for(pos=0;;pos++)
    {
        std::getline(is,val);

        //! check termination
        if(strcmp(val.c_str(),")")==0) break;
        sscanf(val.c_str(),"%d", &n);
        elementsByFaces.insert(n,pos);
    }

    //cout<<"OpenFoamReader::getCellCenters()->____\"owner\" file read: number of entries: "<<NE<<"____"<<endl;

    //! ----------------------
    //! read "neighbour" file
    //! skip the header
    //! ----------------------
    for(int n=0;n<=30;n++)
    {
        std::getline(is1,val);
        if(strcmp("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //",val.c_str())==0) break;
    }

    int NN;
    for(;;)
    {
        std::getline(is1,val);
        if(sscanf(val.c_str(),"%d",&NN)==1) break;
    }

    //cout<<"OpenFoamReader::getCellCenters()->____Number of \"neighbour\" entries: "<<NN<<"____"<<endl;

    //! capture "("
    std::getline(is1,val);
    for(pos=0;;pos++)
    {
        std::getline(is1,val);

        //! check termination
        if(strcmp(val.c_str(),")")==0) break;
        sscanf(val.c_str(),"%d", &n);
        elementsByFaces.insert(n,pos);
    }

    //! ------------------------------------------------------
    //! aim: decompose the polyhedron boundary into triangles
    //! iterate over the cells
    //! ------------------------------------------------------
    //cout<<"OpenFoamReader::getCellCenters()->____decomposing cell faces into triangles____"<<endl;

    FILE *fc = fopen("D:/centers.txt","w");
    for(NE=elementsByFaces.firstKey(); NE<=elementsByFaces.lastKey(); NE++)
    {
        //! -------------------------------------------------------------------------
        //! this is the ->final<- decomposition of a polyhedron into triangles
        //! std::vector<double> => array of 9 triangle coordinates
        //! std::vector<std::vector<double*>> => the triangulation of the polyhedron
        //! -------------------------------------------------------------------------
        std::vector<std::vector<double>> polyhedronTriangulation;

        //! -----------------------------------------
        //! the current element defined by its faces
        //! -----------------------------------------
        QList<int> theElementByFacePoss = elementsByFaces.values(NE);

        //! ----------------------------------------------
        //! bad cell faces - jump over non valid elements
        //! ----------------------------------------------
        if(theElementByFacePoss.size()<4)
        {
            cerr<<"____bad face found____"<<endl;
            continue;
        }

        //! ---------------
        //! build the cell
        //! ---------------
        for(QList<int>::iterator itFaces = theElementByFacePoss.begin(); itFaces!=theElementByFacePoss.end(); ++itFaces)
        {
            int facePos = *itFaces;

            //! -----------------------------
            //! get the node IDs of the face
            //! -----------------------------
            std::vector<int> nodeIDs = listOfFaceDefinitions[facePos];

            //cout<<"____working on face: (";
            //for(int i=0; i<nodeIDs.size()-1; i++) cout<<nodeIDs[i]<<", ";
            //cout<<nodeIDs[nodeIDs.size()-1]<<")____"<<endl;

            //! -------------------------------------------------
            //! transform the coordinates of the nodes of a face
            //! -------------------------------------------------
            double* rotationMatrix = myMapOfRotationMatrices.value(facePos);
            double a11 = rotationMatrix[0], a12 = rotationMatrix[1], a13 = rotationMatrix[2];
            double a21 = rotationMatrix[3], a22 = rotationMatrix[4], a23 = rotationMatrix[5];
            double a31 = rotationMatrix[6], a32 = rotationMatrix[7], a33 = rotationMatrix[8];

            double *P1 = myPointList[nodeIDs[0]];
            double xP1 = P1[0], yP1 = P1[1], zP1 = P1[2];

            //! ----------------------
            //! the transformed nodes
            //! ----------------------
            //std::vector<double*> transfPointList;
            std::vector<std::vector<double>> transfPointList;

            //cout<<"____transforming the nodes of the face____"<<endl;

            //! ------------------
            //! working on a face
            //! ------------------
            for(std::vector<int>::iterator it1 = nodeIDs.begin(); it1!=nodeIDs.end(); ++it1)
            {
                int curNodeID = *it1;
                double *curPoint = myPointList[curNodeID];

                //! ---------------------------------------------------
                //! the current point into the local coordinate system
                //! ---------------------------------------------------
                std::vector<double> transPoint;

                //! ---------------------------------
                //! from global to local coordinates
                //! ---------------------------------
                double dxP = curPoint[0]-xP1;
                double dyP = curPoint[1]-yP1;
                double dzP = curPoint[2]-zP1;

                double XP = dxP*a11+dyP*a12+dzP*a13;
                double YP = dxP*a21+dyP*a22+dzP*a23;
                double ZP = dxP*a31+dyP*a32+dzP*a33;

                transPoint.push_back(XP);
                transPoint.push_back(YP);
                transPoint.push_back(ZP);   //! must be = 0 by definition
                transfPointList.push_back(transPoint);
            }

            //cout<<"____triangulate the face in 2D____"<<endl;

            //! -----------------------------------------------------------------------
            //! [2] triangulate the polygon in 2D
            //! definition: P[j] = {X[j], Y[j]};  Z[j] = 0.0 foreach j, by definition
            //! -----------------------------------------------------------------------
            int numberOfVertexes = int(nodeIDs.size());   //! number of polygon vertexes
            double *X = new double[numberOfVertexes];
            double *Y = new double[numberOfVertexes];

            for(int j=0; j<numberOfVertexes; j++)
            {
                std::vector<double> transformedPoint = transfPointList[j];
                X[j] = transformedPoint[0];
                Y[j] = transformedPoint[1];
            }

            //! --------------------------------------------------------------------
            //! triangulate the current face of the polyhedron
            //!
            //! Input, int N, the number of vertices.
            //! Input, double X[N], Y[N], the coordinates of each vertex.
            //! Output, int TRIANGLES[3*(N-2)], the triangles of the triangulation
            //!
            //! --------------------------------------------------------------------
            int *trs = new int[3*(numberOfVertexes-2)];
            bool isDone = polygon_triangulate(numberOfVertexes, X, Y, trs);
            if(!isDone) continue;

            //cout<<"____face triangulation done____"<<endl;

            std::vector<std::vector<int>> listOfTriads;
            int NbTriads = numberOfVertexes-2;
            //cout<<"*****************************"<<endl;
            for(int i=0; i<NbTriads; i++)
            {
                //! -------------------------------------------
                //! definition of a triangle through local IDs
                //! -------------------------------------------
                std::vector<int> aTriad;
                int s = 3*i;
                aTriad.push_back(trs[s  ]);
                aTriad.push_back(trs[s+1]);
                aTriad.push_back(trs[s+2]);
                listOfTriads.push_back(aTriad);
                //cout<<"____triad("<<aTriad[0]<<", "<<aTriad[1]<<", "<<aTriad[2]<<")____"<<endl;
            }
            //cout<<"*****************************"<<endl;

            //cout<<"____transforming the face triangulation backward____"<<endl;

            //! -------------------------------------------------------
            //! transform back each triangle of the face triangulation
            //! -------------------------------------------------------
            for(int k=0; k<NbTriads; k++)
            {
                std::vector<int> aTriad = listOfTriads[k];

                //cout<<"*****************************"<<endl;
                std::vector<double> aFaceTriangle;
                for(int j=0; j<3; j++)
                {
                    //! -------------------
                    //! transform backward
                    //! -------------------
                    int index = aTriad[j];

                    //! -----------------------------------------------
                    //! a point of the triad (defined into a 2D plane)
                    //! -----------------------------------------------
                    double X_ = transfPointList[index][0];
                    double Y_ = transfPointList[index][1];

                    //! ---------------------------------
                    //! the same point into the 3D space
                    //! ---------------------------------
                    double x = xP1+a11*X_+a21*Y_;
                    double y = yP1+a12*X_+a22*Y_;
                    double z = zP1+a13*X_+a23*Y_;

                    aFaceTriangle.push_back(x);
                    aFaceTriangle.push_back(y);
                    aFaceTriangle.push_back(z);

                    //cout<<"___P("<<j+1<<")"<<"\t("<<x<<", "<<y<<", "<<z<<")____"<<endl;
                }
                //cout<<"*****************************"<<endl;

                polyhedronTriangulation.push_back(aFaceTriangle);
            }

            //! ----------------------------------------------
            //! delete inputs & output of polygon_triangulate
            //! ----------------------------------------------
            //cout<<"____deleting the polygon_triangulate inputs and output____"<<endl;
            delete [] trs;
            delete [] X;
            delete [] Y;
        }

        //cout<<"____calculating the center of the cell____"<<endl;

        //! ------------------------------------
        //! calculate the center of the element
        //! ------------------------------------
        double Area_i = 0.0, AreaTot = 0.0;
        double AiRi_x = 0.0, AiRi_y = 0.0, AiRi_z = 0.0;

        //! -------------------------------------
        //! working on the decomposed polyhedron
        //! -------------------------------------
        int NbTriangles = int(polyhedronTriangulation.size());
        for(int j=0; j<NbTriangles; j++)
        {
            std::vector<double> curTrig = polyhedronTriangulation[j];

            //! ---------------------------------------
            //! the three points defining the triangle
            //! ---------------------------------------
            double xAi = curTrig[0]; double yAi = curTrig[1]; double zAi = curTrig[2];
            double xBi = curTrig[3]; double yBi = curTrig[4]; double zBi = curTrig[5];
            double xCi = curTrig[6]; double yCi = curTrig[7]; double zCi = curTrig[8];

            Area_i = fabs(xAi*(yBi-yCi)+xBi*(yCi-yAi)+xCi*(yAi-yBi));
            AreaTot = AreaTot + Area_i;

            double Ri_x = (xAi+xBi+xCi)/3.0;
            double Ri_y = (yAi+yBi+yCi)/3.0;
            double Ri_z = (zAi+zBi+zCi)/3.0;

            AiRi_x = AiRi_x + Ri_x * Area_i;
            AiRi_y = AiRi_y + Ri_y * Area_i;
            AiRi_z = AiRi_z + Ri_z * Area_i;
        }
        if(AreaTot!=0.0)
        {
            double Cx = AiRi_x/AreaTot;
            double Cy = AiRi_y/AreaTot;
            double Cz = AiRi_z/AreaTot;
            //double *aCenter = new double[3];
            //aCenter[0] = Cx; aCenter[1] = Cy; aCenter[2] = Cz;
            std::vector<double> aCenter{Cx,Cy,Cz};

            fprintf(fc,"%lf\t%lf\t%lf\n",Cx,Cy,Cz);
            //cout<<"____CENTER("<<Cx<<", "<<Cy<<", "<<Cz<<")____"<<endl;
            myMapOfCellCenters.insert(NE,aCenter);
        }
        else
        {
            //! ----------------
            //! console messges
            //! ----------------
            cout<<"OpenFoamReader::getCellCenters()->____face area is zero or polyhedron triangle(s) is(are) missing: jumping over the element____"<<endl;
            continue;
        }
    }
    fclose(fc);
}

//! -----------------------
//! function: openWithLock
//! details:
//! -----------------------
void OpenFoamReader::openWithLock(fstream &stream, const char *fileName, std::ios_base::openmode theOpenMode)
{
    cout<<"OpenFoamReader::openWithLock()->____function called on file"<<fileName<<"____"<<endl;

    //! --------------
    //! open the file
    //! --------------
    stream.open(fileName,theOpenMode);

    //! -----------------------------------
    //! create the corresponding lock file
    //! the lock file is created only if
    //! the open mode is append
    //! -----------------------------------
    if(theOpenMode==std::ios_base::app)
    {
        std::string lockFileName = std::string(fileName)+".lck";
        ofstream lockFile;
        lockFile.open(lockFileName);
    }
}
//! -----------------------------
//! function: closeAndRemoveLock
//! details:
//! -----------------------------
void OpenFoamReader::closeAndRemoveLock(fstream &stream, const char *fileName)
{
    stream.close();

    //! ----------------------
    //! name of the lock file
    //! ----------------------
    std::string lockFileName = std::string(fileName)+".lck";
    remove(lockFileName.c_str());
}
