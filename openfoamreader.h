#ifndef OPENFOAMREADER_H
#define OPENFOAMREADER_H

//! Qt
#include <QObject>
#include <QDir>
#include <QString>
#include <QMap>
#include <QTime>

//! OCC
#include <TColStd_PackedMapOfInteger.hxx>

//! C++
#include <vector>
#include <iostream>
#include <fstream>
//! Custom
#include <customtablemodel.h>
#include <simulationnodeclass.h>

using namespace std;

class QProgressIndicator;

class OpenFoamReader: public QObject
{
    Q_OBJECT

public:

    //! ---------------------------
    //! boundary block in OpenFoam
    //! Ex:
    //! ShotSleeve
    //! {
    //!    type            wall;
    //!    inGroups        1(wall);
    //!    nFaces          11300;
    //!    startFace       132636;
    //! }
    //! ---------------------------
    struct BB
    {
        char boundaryName[256];
        int startFace;
        int nFaces;
    };

    //! --------------------------------------------------------
    //! fileMode= 0 - a single output file, for each time
    //! fileMode= 1 - multiple files, one for each time
    //! fileMode= 2 - performed in memory - not implemented yet
    //! --------------------------------------------------------
    OpenFoamReader(/*const QString &sourceDirPath, const QString &targetDirPath,*/ int fileMode=0, QObject *parent=0);

    //! -----------------------------------------
    //! set the target and source directory path
    //! -----------------------------------------
    //void setTargetDir(const QString &theTargetDir);
    //void setSourceDir(const QString &theSourceDir);

    //! ---------------------------
    //! set the progress indicator
    //! ---------------------------
    void setProgressIndicator(QProgressIndicator *aProgressIndicator);

    //! -----------
    //! destructor
    //! -----------
    ~OpenFoamReader();

private:

    //! -----------------------------------
    //! source and target file directories
    //! -----------------------------------
    QString mySourceDir;
    QString myTargetDir;

    //! -------------------------
    //! map of rotation matrices
    //! -------------------------
    QMap<int, double*> myMapOfRotationMatrices;

    //! -----------
    //! point list
    //! -----------
    std::vector<double*> myPointList;

    //! -------------
    //! face centers
    //! -------------
    QMap<int, double*> myMapOfFaceCenters;

    //! -------------
    //! cell centers
    //! -------------
    QMap<int,std::vector<double>> myMapOfCellCenters;

    //! --------------------
    //! file splitting mode
    //! --------------------
    int myFileMode;

    //! -----------------------
    //! the progress indicator
    //! -----------------------
    QProgressIndicator *myProgressIndicator;

#ifdef COSTAMP_VERSION
    std::vector<double> myTimeFolders;
#endif

private:

    //! --------------------------
    //! read the boundary section
    //! --------------------------
    void readBoundarySection(std::vector<BB> &vecBoundaries, fstream &is);

    //! ------------
    //! read points
    //! ------------
    void readPoints(fstream &is);

    //! -----------
    //! read faces
    //! -----------
    void readFaces(fstream &is, std::vector<std::vector<int>> &listOfFaceDefinitions);

    //! -------------------------
    //! get the center of a cell
    //! -------------------------
    void getCellCenters(fstream &is,
                        fstream &is1,
                        const std::vector<std::vector<int>> &listOfFaceDefinitions
                        /*, const QMap<int, double*> &mapRotationMatrices, */
                        /*const std::vector<double*> &pointList, */
                        /*QMap<int, double*> &centers*/);

    //! -----------
    //! freeMemory
    //! -----------
    void freeMemory();

    //! ------
    //! clock
    //! ------
    std::string clock()
    {
        QTime time = QTime::currentTime();
        QString text = time.toString("mm:ss:zzz");
        return text.toStdString();
    }

    //! -------------------------
    //! open and close with lock
    //! -------------------------
    void openWithLock(fstream &stream, const char *fileName, std::ios_base::openmode theOpenMode);
    void closeAndRemoveLock(fstream &stream, const char *fileName);

signals:

    //! signal for results handling
    void conversionFinished();

public slots:

    //! perform
    bool perform(SimulationNodeClass *OFnode);

    //! timer
    void updateTimer();

#ifdef COSTAMP_VERSION
    void setTimeFolders(const std::vector<double> timeFolders)
    {myTimeFolders = timeFolders;}
#endif
};

#endif // OPENFOAMREADER_H
