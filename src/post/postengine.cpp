//! ----------------
//! custom includes
//! ----------------
#include "postengine.h"
#include "posttools.h"
#include "occmeshtoccxmesh.h"
#include <rainflow.h>

//! ---
//! Qt
//! ---
#include <QDirIterator>
#include <QDateTime>

//! ----
//! C++
//! ----
#include <iostream>
#include <fstream>

//! -------
//! global
//! -------
#include "global.h"

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
postEngine::postEngine(QObject *parent) : QObject(parent),myDTM(QMap<double,QVector<int>>())
{
    this->buildMap();
}

//! -----------------------------
//! function: setDiscreteTimeMap
//! details:
//! -----------------------------
void postEngine::setDiscreteTimeMap(const QMap<double,QVector<int>> &dtm)
{
    myDTM = dtm;
}

//! ------------------
//! function: perform
//! details:
//! ------------------
bool postEngine::perform()
{
    if(myMeshDataBase==NULL) return false;
    if(myResultsFilePath.isNull()) return false;
    if(myResultsFilePath.isEmpty()) return false;

    FrdReader anFrdReader(myResultsFilePath);
    bool isDone = anFrdReader.perform();
    return isDone;
}

//! ----------------------
//! function: setDataBase
//! details:
//! ----------------------
void postEngine::setDataBase(meshDataBase *mDB)
{
    myMeshDataBase = mDB;
}

//! -------------------------
//! function: setResultsFile
//! details:
//! -------------------------
void postEngine::setResultsFile(QString resultsFilePath)
{
    myResultsFilePath  = resultsFilePath;
}

//! -------------------------
//! function: evaluateResult
//! details:
//! --------------------------
QMap<GeometryTag,QList<QMap<int,double>>> postEngine::evaluateResult(const QString &resultKeyName,
                                                                     int requiredSubStepNb,
                                                                     int requiredStepNb,
                                                                     int requiredMode,
                                                                     const QVector<GeometryTag> &vecLoc)
{
    cout<<"postEngine::evaluateResult()->____function called for variable: "<<resultKeyName.toStdString()<<"____"<<endl;
    cout<<"postEngine::evaluateResult()->____on Nr: "<<vecLoc.length()<<" locations____"<<endl;

    //! ----------------------------------------------------
    //! generate the results on all the requested locations
    //! ----------------------------------------------------
    QMap<GeometryTag,QList<QMap<int,double>>> resMap;

    QVector<GeometryTag>::const_iterator it;
    for(it = vecLoc.cbegin(); it!= vecLoc.cend(); ++it)
    {
        GeometryTag loc = *it;

        //! -------------------------------------------------------------------------
        //! node conversion map: (Calculix mesh nodeID,nodeID for MeshVS_DataSource)
        //! -------------------------------------------------------------------------
        QMap<int,int> indexedMapOfNodes = OCCMeshToCCXmesh::perform(loc,myMeshDataBase);

        //! ------------------------------------------------
        //! enter <...>/SolutionData/ResultsData
        //! ------------------------------------------------
        QString tmp = myResultsFilePath.split("/").last();
        QString path = myResultsFilePath;
        path.chop(tmp.length());

        QDir curDir(path);
        //cout<<"____"<<curDir.absolutePath().toStdString()<<"____"<<endl;
        curDir.cd("ResultsData");
        //cout<<"____"<<curDir.absolutePath().toStdString()<<"____"<<endl;
        QFileInfoList entriesInfo = curDir.entryInfoList();
        //cout<<"____"<<entriesInfo.length()<<"____"<<endl;

        QList<QString> entryList = curDir.entryList();
        QList<QString> fileList;

        //! ------------------------------------------------------
        //! retrieve the files (discard directories)
        //! it could be used to setup the range of a progress bar
        //! ------------------------------------------------------
        for(int k=0; k<entryList.length(); k++)
        {
            if(entriesInfo.at(k).isFile())
            {
                QString fileName = curDir.absolutePath()+"/"+entryList.at(k);
                fileList.append(fileName);
                //cout<<"adding file name: "<<entryList.at(k).toStdString()<<endl;
            }
        }

        //! --------------------------
        //! the results on a location
        //! --------------------------
        QList<QMap<int,double>> res;

        //! ---------------
        //! scan the files
        //! ---------------
        for(int i=0; i<fileList.length(); i++)
        {
            //cout<<"scanning file name: "<<fileList.at(i).toStdString()<<endl;
            QString filePath = fileList.at(i);
            fstream curFile(filePath.toStdString(),ios::in);
            if(!curFile.is_open())
            {
                cout<<"postEngine::evaluateResult()->____cannot open the results file: "<<resultKeyName.toStdString()<<"____"<<endl;
                break;
            }
            std::string val;
            std::getline(curFile,val);
            char analysisType[24];
            double mode;
            sscanf(val.c_str(),"%s%lf",&analysisType,&mode);

            std::getline(curFile,val);
            double time;
            sscanf(val.c_str(),"Time= %lf",&time);
            std::getline(curFile,val);
            int subStepNb, stepNb;
            sscanf(val.c_str(),"Substep n=%d Step n=%d",&subStepNb,&stepNb);
            std::getline(curFile,val);
            char tdata[32];
            sscanf(val.c_str(),"%s",tdata);

            //printf("File %s Time= %lf Substep n=%d Step n=%d\n",tdata,time,subStepNb,stepNb);
            //printf("compare to File %s Time= %lf Substep n=%d Step n=%d\n",resultKeyName.toStdString().c_str(),time,requiredSubStepNb,requiredStepNb);

            if(strcmp(tdata,resultKeyName.toStdString().c_str())==0 && subStepNb==requiredSubStepNb && stepNb == requiredStepNb && mode == requiredMode)
            {
                //printf("file @ required time found\n");
                TypeOfResult tor = m.value(tdata);
                switch(tor)
                {
                case TypeOfResult_HFL:
                {
                    QMap<int,double> resComp_normal;
                    occHandle(Ng_MeshVS_DataSourceFace) curFaceDS = occHandle(Ng_MeshVS_DataSourceFace)::
                            DownCast(myMeshDataBase->ArrayOfMeshDSOnFaces.getValue(loc.parentShapeNr,loc.subTopNr));

                    if(curFaceDS->myNodeNormals.isEmpty()) curFaceDS->computeNormalAtNodes();

                    //! <>::eof(): call getline before while, then inside {}, @ as last instruction
                    std::getline(curFile,val);
                    while(curFile.eof()!=true)
                    {
                        int ni;
                        double cxx,cyy,czz;
                        sscanf(val.c_str(),"%d%lf%lf%lf",&ni,&cxx,&cyy,&czz);

                        //! nodeIDs defining the MeshVS_dataSource
                        int OCCnodeID = indexedMapOfNodes.value(ni,-1);
                        if(OCCnodeID!=-1)
                        {
                            const QList<double> &normal = curFaceDS->myNodeNormals.value(OCCnodeID);
                            double normalFlux = cxx*normal[0]+cyy*normal[1]+czz*normal[2];
                            resComp_normal.insert(OCCnodeID,normalFlux);
                        }
                        std::getline(curFile,val);
                    }

                    //! result
                    res<<resComp_normal;
                    //cout<<"postEngine::evaluateResult()->____Number of components: "<<res.length()<<"____"<<endl;
                }
                    break;
                case TypeOfResult_U:
                case TypeOfResult_F:
                //case TypeOfResult_HFL:
                {
                    QMap<int,double> resComp_X,resComp_Y,resComp_Z,resComp_Total;

                    //! <>::eof(): call getline before while, then inside {}, @ as last instruction
                    std::getline(curFile,val);
                    while(curFile.eof()!=true)
                    {
                        int ni;
                        double cxx,cyy,czz,total;
                        sscanf(val.c_str(),"%d%lf%lf%lf",&ni,&cxx,&cyy,&czz);

                        //! nodeIDs defining the MeshVS_dataSource
                        int OCCnodeID = indexedMapOfNodes.value(ni,-1);
                        if(OCCnodeID!=-1)
                        {
                            total = sqrt(pow(cxx,2)+pow(cyy,2)+pow(czz,2));
                            resComp_Total.insert(OCCnodeID,total);
                            resComp_X.insert(OCCnodeID,cxx);
                            resComp_Y.insert(OCCnodeID,cyy);
                            resComp_Z.insert(OCCnodeID,czz);
                        }
                        std::getline(curFile,val);
                    }

                    //! result
                    res<<resComp_Total<<resComp_X<<resComp_Y<<resComp_Z;
                    //cout<<"postEngine::evaluateResult()->____Number of components: "<<res.length()<<"____"<<endl;
                }
                    break;

                case TypeOfResult_S:
                case TypeOfResult_TOSTRAIN:
                case TypeOfResult_MESTRAIN:
                {
                    //!                   0       1        2      3       4        5       6       7       8      9      10
                    QMap<int,double> resMISES, resSINT, resSI, resSII, resSIII, resSXX, resSYY, resSZZ, resSXY,resSYZ,resSXZ;

                    //! <>::eof(): call getline before while, then inside {}, @ as last instruction
                    std::getline(curFile,val);
                    while(curFile.eof()!=true)
                    {
                        //! read the components of the 3x3 data
                        int ni;
                        double cxx,cyy,czz,cxy,cyz,cxz;
                        sscanf(val.c_str(),"%d%lf%lf%lf%lf%lf%lf",&ni,&cxx,&cyy,&czz,&cxy,&cyz,&cxz);

                        //! nodeIDs defining the MeshVS_dataSource
                        int OCCnodeID = indexedMapOfNodes.value(ni,-1);
                        if(OCCnodeID!=-1)
                        {
                            //! --------------------------------------------
                            //! compute the equivalent stress/strain
                            //! --------------------------------------------
                            double vonMises;
                            if(tor==TypeOfResult_TOSTRAIN || tor == TypeOfResult_MESTRAIN)   // equivalent strain formula
                            {
                                vonMises = (2.0/3.0)*sqrt((3.0/2.0)*(cxx*cxx+cyy*cyy+czz*czz)+(3.0/4.0)*(cxy*cxy+cyz*cyz+cxz*cxz));
                            }
                            else    // equivalent stress formula
                            {
                                //vonMises = sqrt(cxx*cxx+cyy*cyy+czz*czz-cxx*cyy-cyy*czz-czz*cxx+3*(cxz*cxz+cyz*cyz+cxz*cxz));
                                vonMises = sqrt(0.5*(pow(cxx-cyy,2)+pow(cxx-czz,2)+pow(cyy-czz,2))+3*(cxz*cxz+cyz*cyz+cxz*cxz));
                            }
                            resMISES.insert(OCCnodeID,vonMises);

                            //! --------------------------------------------
                            //! compute the principal components
                            //! --------------------------------------------
                            QList<double> sik;
                            sik<<cxx<<cyy<<czz<<cxy<<cyz<<cxz;
                            QList<double> s = postTools::principalComponents(sik);

                            resSI.insert(OCCnodeID,s.at(2));        //! maximum
                            resSII.insert(OCCnodeID,s.at(1));       //! middle
                            resSIII.insert(OCCnodeID,s.at(0));      //! minimum

                            //! -----------------------------------------------
                            //! compute the stress/strain intensity (2*Tresca)
                            //! maximum shear stress/strain
                            //! -----------------------------------------------
                            double sint = fabs(s.at(2)-s.at(0));
                            resSINT.insert(OCCnodeID,sint);

                            resSXX.insert(OCCnodeID,cxx);
                            resSYY.insert(OCCnodeID,cyy);
                            resSZZ.insert(OCCnodeID,czz);
                            resSXY.insert(OCCnodeID,cxy);
                            resSYZ.insert(OCCnodeID,cyz);
                            resSXZ.insert(OCCnodeID,cxz);
                        }
                        std::getline(curFile,val);
                    }
                    //! result
                    res<<resMISES<<resSINT<<resSI<<resSII<<resSIII<<resSXX<<resSYY<<resSZZ<<resSXY<<resSYZ<<resSXZ;
                    //cout<<"postEngine::evaluateResult()->____Number of components: "<<res.length()<<"____"<<endl;;
                }
                    break;

                case TypeOfResult_NT:
                case TypeOfResult_EPS:
                {
                    QMap<int,double> resT;

                    //! <>::eof(): call getline before while, then inside {}, @ as last instruction
                    std::getline(curFile,val);
                    while(curFile.eof()!=true)
                    {
                        int ni;
                        double v;
                        sscanf(val.c_str(),"%d%lf",&ni,&v);

                        //! nodeIDs defining the MeshVS_dataSource
                        int OCCnodeID = indexedMapOfNodes.value(ni,-1);
                        if(OCCnodeID!=-1)
                        {
                            resT.insert(OCCnodeID,v);
                        }
                        std::getline(curFile,val);
                    }
                    //! result
                    res<<resT;
                    //cout<<"Number of components: "<<res.length();
                }
                    break;

                case TypeOfResult_CONT:
                {
                    //!                         col 1           col 2+3     col 4           col 5+6
                    QMap<int,double>  resContPenetration, resContSliding,resContPress, resContFrictStress;

                    //! <>::eof(): call getline before while, then inside {}, @ as last instruction
                    std::getline(curFile,val);
                    while(curFile.eof()!=true)
                    {
                        //! read the components of the 3x3 data
                        int ni;
                        double cxx,cyy,czz,cxy,cyz,cxz;
                        sscanf(val.c_str(),"%d%lf%lf%lf%lf%lf%lf",&ni,&cxx,&cyy,&czz,&cxy,&cyz,&cxz);

                        //! nodeIDs defining the MeshVS_dataSource
                        int OCCnodeID = indexedMapOfNodes.value(ni,-1);
                        if(OCCnodeID!=-1)
                        {
                            //! --------------------------------------------------
                            //! compute the frictional stress and contact sliding
                            //! -------------------------------------------------
                            double frictStress = sqrt(cyz*cyz+cxz*cxz);
                            double contSliding = sqrt(cyy*cyy+czz*czz);

                            resContFrictStress.insert(OCCnodeID,frictStress);
                            resContSliding.insert(OCCnodeID,contSliding);
                            resContPress.insert(OCCnodeID,cxy);
                            resContPenetration.insert(OCCnodeID,cxx);
                        }
                        std::getline(curFile,val);
                    }
                    //! result
                    res<<resContPress<<resContFrictStress<<resContPenetration<<resContSliding;
                }
                    break;
                }
                resMap.insert(loc,res);
                curFile.close();
                break;
            }
            else
            {
                curFile.close();
            }
            curFile.close();
        }
        //! --------------------------------------------------------
        //! resMap => QMap<GeometryTag,QList<QMap<int,double>>
        //! res => QList<QMap<int,double>>
        //! --------------------------------------------------------
    }

    //! diagnostic function
    return resMap;
}

//! ----------------------------------------------------------
//! function: colorBoxTitle
//! details:  build a title for the colorbox according to the
//!           type of result
//! ----------------------------------------------------------
QString postEngine::colorBoxTitle(const QString &keyName, int component, int step, int subStep)
{
    QString timeInfo = QString("\nStep %1\nSubstep %2").arg(step).arg(subStep);

    QString colorBoxTitle;

    if(keyName =="Damage")
    {
        colorBoxTitle = keyName;
        return colorBoxTitle.append("\n");
    }

    TypeOfResult tor = m.value(keyName);
    switch(tor)
    {
    case TypeOfResult_HFL:
        switch(component)
        {
        case 0: colorBoxTitle="Thermal flux"; break;
        }
        break;
    case TypeOfResult_CONT:
        switch(component)
        {
        case 0: colorBoxTitle="Contact pressure"; break;
        case 1: colorBoxTitle="Frictional stress"; break;
        case 2: colorBoxTitle="Contact penetration"; break;
        case 3: colorBoxTitle="Contact sliding"; break;
        }
        break;
    case TypeOfResult_U:
        switch(component)
        {
        //case 0: colorBoxTitle="Total displacement"; break;
        case 0: colorBoxTitle="Total displacement"; break;
        case 1: colorBoxTitle="Directional displacement X"; break;
        case 2: colorBoxTitle="Directional displacement Y"; break;
        case 3: colorBoxTitle="Directional displacement Z"; break;
        }
        break;
    case TypeOfResult_S:
        switch(component)
        {
        case 0: colorBoxTitle="Equivalent stress"; break;
        case 1: colorBoxTitle = "Stress intensity"; break;
        case 2: colorBoxTitle = "Maximum principal stress"; break;
        case 3: colorBoxTitle = "Middle principal stress"; break;
        case 4: colorBoxTitle = "Minimum principal stress"; break;
        case 5: colorBoxTitle = "Normal stress X"; break;
        case 6: colorBoxTitle = "Normal stress Y"; break;
        case 7: colorBoxTitle = "Normal stres Z"; break;
        case 9: colorBoxTitle = "Shear stress XY"; break;
        case 10: colorBoxTitle = "Shear stress YZ"; break;
        case 11: colorBoxTitle = "Shear stress ZX"; break;
        }
        break;
    case TypeOfResult_TOSTRAIN:
        switch(component)
        {
        case 0: colorBoxTitle="Equivalent strain"; break;
        case 1: colorBoxTitle = "Strain intensity"; break;
        case 2: colorBoxTitle = "Maximum principal strain"; break;
        case 3: colorBoxTitle = "Middle principal strain"; break;
        case 4: colorBoxTitle = "Minimum principal strain"; break;
        }
        break;
    case TypeOfResult_MESTRAIN:
        switch(component)
        {
        case 0: colorBoxTitle="Equivalent strain"; break;
        case 1: colorBoxTitle = "Strain intensity"; break;
        case 2: colorBoxTitle = "Maximum principal strain"; break;
        case 3: colorBoxTitle = "Middle principal strain"; break;
        case 4: colorBoxTitle = "Minimum principal strain"; break;
        }
        break;
    case TypeOfResult_NT:
        switch(component)
        {
        case 0: colorBoxTitle = "Temperature"; break;
        }
        break;
    case TypeOfResult_F:
        switch(component)
        {
        case 0: colorBoxTitle="Total force"; break;
        case 1: colorBoxTitle="Directional Force X"; break;
        case 2: colorBoxTitle="Directional Force Y"; break;
        case 3: colorBoxTitle="Directional Force Z"; break;
        }
        break;
    case TypeOfResult_EPS:
        switch(component)
        {
        case 0: colorBoxTitle="Equivalent Plastic Strain"; break;
        }
        break;
    }
    return this->timeStamp().append("\n").append(colorBoxTitle).append(timeInfo).append("\n");
}

//! --------------------------
//! function: buildPostObject
//! details:
//! --------------------------
postObject postEngine::buildPostObject(const QString &keyName,
                                       int component,
                                       int requiredSubStepNb,
                                       int requiredStepNb,
                                       int requiredMode,
                                       const QVector<GeometryTag> &vecLoc)
{
    //! -------------------------
    //! build the colorBox title
    //! -------------------------
    QString aColorBoxTitle = this->colorBoxTitle(keyName, component, requiredStepNb, requiredSubStepNb);

    //! --------------------
    //! call the postEngine
    //! --------------------
    QMap<GeometryTag,QList<QMap<int,double>>> resMap = this->evaluateResult(keyName, requiredSubStepNb, requiredStepNb, requiredMode, vecLoc);

    //! ------------------------------------------------------------------------------------------------------------
    //! create the map of nodal vectorial displacements for the deformed mesh presentation. Here:
    //! QMap<int,gp_Vec> displMap                    => map of nodal vectorial displacements
    //! QMap<GeometryTag,QList<QMap<int,gp_Vec>>> => each location has its own map of nodal vectorial displacements
    //! ------------------------------------------------------------------------------------------------------------
    QMap<int,gp_Vec> displMap;
    QMap<GeometryTag,QMap<int,gp_Vec>> mapDisplMap;
    QMap<GeometryTag,QList<QMap<int,double>>> nodalDisplacements = this->evaluateResult("DISP", requiredSubStepNb, requiredStepNb,requiredMode, vecLoc);

    for(QMap<GeometryTag,QList<QMap<int,double>>>::iterator it = nodalDisplacements.begin(); it!=nodalDisplacements.end(); ++it)
    {
        const GeometryTag &aLoc= it.key();

        QList<QMap<int,double>> nodalDisplacementsComponents = it.value();
        QMap<int,double> displX = nodalDisplacementsComponents[1];
        QMap<int,double> displY = nodalDisplacementsComponents[2];
        QMap<int,double> displZ = nodalDisplacementsComponents[3];

        QMap<int,double>::iterator itX = displX.begin();
        QMap<int,double>::iterator itY = displY.begin();
        QMap<int,double>::iterator itZ = displZ.begin();

        for(;itX!=displX.end() && itY!=displY.end() && itZ!=displZ.end(); ++itX, ++itY, ++itZ)
        {
            int nodeID = itX.key();
            double x = itX.value();
            double y = itY.value();
            double z = itZ.value();
            gp_Vec aVec(x,y,z);
            displMap.insert(nodeID,aVec);
        }
        mapDisplMap.insert(aLoc,displMap);
    }

    //! ----------------------------------------------------------------
    //! create the postObject
    //! the last options create the post object using the surface mesh,
    //! if the underlying mesh is a volume mesh
    //! ----------------------------------------------------------------
    bool showSolidMeshAsSurface = Global::status().isVolumeMeshShownAsSurface;
    postObject aPostObject(resMap,vecLoc,mapDisplMap,aColorBoxTitle,showSolidMeshAsSurface);
    return aPostObject;
}

//! ------------------------------
//! function: plotDataSummary
//! details:  diagnostic function
//! ------------------------------
void postEngine::plotDataSummary(QMap<GeometryTag, QList<QMap<int,double>>> data)
{
    if(!data.isEmpty())
    {
        QMap<GeometryTag,QList<QMap<int,double>>>::const_iterator anIt;
        for(anIt = data.cbegin(); anIt!= data.cend(); ++anIt)
        {
            GeometryTag loc = anIt.key();
            QList<QMap<int,double>> l = anIt.value();
        }
    }
    else
    {
        cout<<"---------->no result has been found<-----------------"<<endl;
    }
}

//! --------------------
//! function: timeStamp
//! details:
//! --------------------
QString postEngine::timeStamp()
{
    //! builds the timestamp
    QDateTime dateTime;
    QString dateFormat = "dd/MM/yyyy";
    QString dateString = dateTime.currentDateTime().toString(dateFormat);
    //QString timeFormat = "hh:mm";
    //QString timeString = dateTime.currentDateTime().toString(timeFormat);
    QString timeStamp;
    //timeStamp.append("Date: ").append(dateString).append("\n").append("Time: ").append(timeString);
    timeStamp.append("Date: ").append(dateString);
    return timeStamp;
}

//! ----------------------------
//! function: updateResultScale
//! details:
//! ----------------------------
void postEngine::updateResultScale(postObject &aPostObject, int scaleType, double minValue, double maxValue, int NbIntervals)
{
    cout<<"postEngine::updateResultScale()->____function called____"<<endl;

    //! ------------------------------------------------------------
    //! retrieve the mesh data sources from the current MeshVS_Mesh
    //! ------------------------------------------------------------
    const QMap<GeometryTag,occHandle(MeshVS_Mesh)> &theMeshes = aPostObject.getColoredMeshes();
    QMap<GeometryTag,opencascade::handle<MeshVS_DataSource>> myMapOfMeshDataSources;
    for(QMap<GeometryTag,occHandle(MeshVS_Mesh)>::const_iterator it = theMeshes.cbegin(); it!=theMeshes.cend(); ++it)
    {
        const GeometryTag &loc = it.key();
        const opencascade::handle<MeshVS_Mesh> &curMeshVS_mesh = it.value();
        myMapOfMeshDataSources.insert(loc,curMeshVS_mesh->GetDataSource());
    }

    //! -------------------
    //! retrieve the scale
    //! -------------------
    //double scale = aPostObject.getScale();

    //! -----------------------
    //! retrieve the component
    //! -----------------------
    int solutionDataComponent = aPostObject.getSolutionDataComponent();

    //! -----------------------------------------
    //! update the colored mesh and the colorbox
    //! the scale is left unchanged
    //! -----------------------------------------
    switch (scaleType)
    {
    case 0:
        //! autoscale min max, custom number of levels
        aPostObject.buildMeshIO(myMapOfMeshDataSources,0,0,NbIntervals,true,solutionDataComponent);
        break;
    case 1:
        //! custom scale (custom min, max, numbre of levels
        aPostObject.buildMeshIO(myMapOfMeshDataSources,minValue,maxValue,NbIntervals,false,solutionDataComponent);
        break;
    }
}

//! ---------------------------------
//! function: evaulateFatigueResults
//! details:
//! ---------------------------------
postObject postEngine::evaluateFatigueResults(int type, QVector<GeometryTag> locs, const QList<double> &times, QMap<int,int> materialBodyMap, int nCycle)
{
    QMap<GeometryTag,QList<QMap<int,double>>> fatigueResults;

    switch(myFatigueModel.type)
    {
    case fatigueModel_BCM:
    {
        cout<<"postEngine::evaluateResult()->____fatigue model BCM called___"<<endl;

        QMap<GeometryTag,QMap<int,QList<double>>> r = readFatigueResults(type,locs,times);
        rainflow rf;
        for(QMap<GeometryTag,QMap<int,QList<double>>>::iterator it = r.begin(); it!=r.end(); ++it)
        {
            GeometryTag curLoc = it.key();
            rf.setLocation(curLoc);

            QMap<int,QList<double>> strainDistTimeHistory = it.value();
            rf.setFatigueModel(myFatigueModel);
            QMap<int,double> damageDist;

            bool isDone = rf.perform(strainDistTimeHistory,damageDist);
            if(isDone)
            {
                QList<QMap<int,double>>damageDistList;
                damageDistList<<damageDist;
                fatigueResults.insert(curLoc,damageDistList);
            }
        }
    }
        break;

    case fatigueModel_ESR:
    {
        cout<<"postEngine::evaluateResult()->____fatigue model ESR called___"<<endl;
        int step,substep;
        //QMap<double,QVector<int>> dtm= getDTM();
        //QMap<double,QVector<int>> dtm = myDTM;

        postTools::getStepSubStepByTimeDTM(myDTM,times.last(),step,substep);
        QString tor_eps = m.key(TypeOfResult_EPS);
        QString tor_mises = m.key(TypeOfResult_S);
        int mode =0;
        QMap<GeometryTag,QList<QMap<int,double>>> pe = this->evaluateResult(tor_eps,substep,step,mode,locs);
        QMap<GeometryTag,QList<QMap<int,double>>> stress = this->evaluateResult(tor_mises,substep,step,mode,locs);

        for(QVector<GeometryTag>::iterator it=locs.begin();it!=locs.end();it++)
        {
            GeometryTag curLoc = *it;
            int bodyIndex = curLoc.parentShapeNr;
            const QList<QMap<int, double>> &listOfResPe = pe.value(curLoc);
            const QMap<int, double> &curPe = listOfResPe.first();
            const QList<QMap<int, double>> &listOfResMises = stress.value(curLoc);
            const QMap<int, double> &curMises = listOfResMises.first();
            QList<QMap<int,double>> damageIndex;
            QMap<int,double> damageIndexData;

            double elasticModulusMedium, elasticModulusMin,r,a,b,c,d,e,f,g,h;
            int material = materialBodyMap.value(bodyIndex);

            for(QMap<int,double>::const_iterator itt=curPe.cbegin(); itt!=curPe.cend(); itt++)
            {
                double altStress,X,Y;
                int curPos = itt.key();
                double eps = itt.value();
                double mises = curMises.value(curPos);

                switch(material)
                {
                case 5:
                case 6:
                case 0:
                case 1:
                case 2:
                case 3:
                case 4:
                case 7:
                case 8:
                case 9:
                {
                    elasticModulusMedium = 1.76000000e+005;
                    altStress = 0.5*(mises+eps*elasticModulusMedium);
                    Y = log10(28.3*pow(10,3)*altStress/elasticModulusMedium);

                    r=35.9;
                    a=9.030556;
                    b=8.1906623;
                    c=0.36077181;
                    d=0.4706984;
                    e=42.08579;
                    f=12.514054;
                    g=4.3290016;
                    h=0.60540862;

                    if(pow(10,Y)>r) X = (a-b*Y)/(1-c*Y-d*Y*Y);
                    else  X = (-e+f*Y)/(1-g*Y+h*Y*Y);
                }
                    break;
                }
                double damage;
                double min,max;
                min=1;
                max=8;
                if(X>min && X<max) damage = nCycle/(pow(10,X));
                else damage = 0;
                damageIndexData.insert(curPos,damage);
            }
            damageIndex<<damageIndexData;
            fatigueResults.insert(curLoc,damageIndex);
        }
    }
        break;
    }

    //! -----------------------
    //! create the post object
    //! -----------------------
    QString label = this->colorBoxTitle("Damage",0,1,1);
    postObject aPostObject(fatigueResults,locs,label);
    //postObject aPostObject(fatigueResults,locs);
    return aPostObject;
}

//! ---------------------------------------------------
//! function: readFatigueResults
//! details:  type = 1 => equivalent mechanical strain
//!           type = 0 => equivalent total strain
//!           ... other types
//! ---------------------------------------------------
QMap<GeometryTag,QMap<int,QList<double>>> postEngine::readFatigueResults(int type,
                                                                         const QVector<GeometryTag> &vecLoc,
                                                                         const QList<double> &times)
{
    QString resultKeyName;
    switch(type)
    {
    case 0: resultKeyName = "TOSTRAIN"; break;
    case 1: resultKeyName = "MESTRAIN"; break;
    //case 2: resultKeyName = "THSTRAIN"; break;
    }

    //! ----------------------------------------------------
    //! generate the results on all the requested locations
    //! ----------------------------------------------------
    QMap<GeometryTag,QMap<int,QList<double>>> resMap;
    QVector<GeometryTag>::const_iterator it;
    for(it = vecLoc.cbegin(); it!= vecLoc.cend(); ++it)
    {
        GeometryTag loc = *it;
        //! -------------------------------------------------------------------------
        //! node conversion map: (Calculix mesh nodeID,nodeID for MeshVS_DataSource)
        //! -------------------------------------------------------------------------
        QMap<int,int> indexedMapOfNodes = OCCMeshToCCXmesh::perform(loc,myMeshDataBase);

        //! -------------------------------------
        //! enter <...>/SolutionData/ResultsData
        //! -------------------------------------
        QString tmp = myResultsFilePath.split("/").last();
        QString path = myResultsFilePath;
        path.chop(tmp.length());

        QDir curDir(path);
        curDir.cd("ResultsData");
        QFileInfoList entriesInfo = curDir.entryInfoList();

        QList<QString> entryList = curDir.entryList();
        QList<QString> fileList;

        //! ------------------------------------------------------
        //! retrieve the files (discard directories)
        //! it could be used to setup the range of a progress bar
        //! ------------------------------------------------------
        for(int k=0; k<entryList.length(); k++)
        {
            if(entriesInfo.at(k).isFile())
            {
                QString fileName = curDir.absolutePath()+"/"+entryList.at(k);
                fileList.append(fileName);
            }
        }

        //! ---------------
        //! scan the files
        //! ---------------
        int n=0;    //! time index
        QMap<int,QList<double>> resMISES;
        for(int i=0; i<fileList.length(); i++)
        {
            QString filePath = fileList.at(i);
            ifstream curFile(filePath.toStdString());

            std::string val;
            std::getline(curFile,val);
            double time;
            sscanf(val.c_str(),"Time= %lf",&time);
            std::getline(curFile,val);
            int subStepNb, stepNb;
            sscanf(val.c_str(),"Substep n=%d Step n=%d",&subStepNb,&stepNb);
            std::getline(curFile,val);
            char tdata[32];
            sscanf(val.c_str(),"%s",tdata);

            int step,substep;
            postTools::getStepSubStepByTimeDTM(myDTM,times.at(n),step,substep);

            QList<double> timeHistory;
            if(strcmp(tdata,resultKeyName.toStdString().c_str())==0 && subStepNb==substep && stepNb == step)
            {
                n++;
                //! ----------------------------------------------------------------------------
                //! <>::eof(): call getline before while, then inside {}, @ as last instruction
                //! ----------------------------------------------------------------------------
                std::getline(curFile,val);
                while(curFile.eof()!=true)
                {
                    //! read the components of the 3x3 data
                    int ni;
                    double cxx,cyy,czz,cxy,cyz,cxz;
                    sscanf(val.c_str(),"%d%lf%lf%lf%lf%lf%lf",&ni,&cxx,&cyy,&czz,&cxy,&cyz,&cxz);

                    //! nodeIDs defining the MeshVS_dataSource
                    int OCCnodeID = indexedMapOfNodes.value(ni,-1);
                    if(OCCnodeID!=-1 && n==0)
                    {
                        //! -------------------------------------
                        //! compute the equivalent stress/strain
                        //! -------------------------------------
                        double vonMises;
                        vonMises = (2.0/3.0)*sqrt((3.0/2.0)*(cxx*cxx+cyy*cyy+czz*czz)+(3.0/4.0)*(cxy*cxy+cyz*cyz+cxz*cxz));
                        timeHistory<<vonMises;
                        resMISES.insert(OCCnodeID,timeHistory);

                        //! -----------------------------------------------------
                        //! compute the principal components: index 0 is minimum
                        //! -----------------------------------------------------
                        QList<double> sik;
                        sik<<cxx<<cyy<<czz<<cxy<<cyz<<cxz;
                        QList<double> s = postTools::principalComponents(sik);

                    }
                    else if(OCCnodeID!=-1)
                    {
                        //! -------------------------------------
                        //! compute the equivalent stress/strain
                        //! -------------------------------------
                        double vonMises;
                        vonMises = (2.0/3.0)*sqrt((3.0/2.0)*(cxx*cxx+cyy*cyy+czz*czz)+(3.0/4.0)*(cxy*cxy+cyz*cyz+cxz*cxz));
                        timeHistory = resMISES.value(OCCnodeID);
                        timeHistory<<vonMises;
                        resMISES.insert(OCCnodeID,timeHistory);

                        //! ---------------------------------
                        //! compute the principal components
                        //! ---------------------------------
                        QList<double> sik;
                        sik<<cxx<<cyy<<czz<<cxy<<cyz<<cxz;
                        QList<double> s = postTools::principalComponents(sik);

                    }
                    std::getline(curFile,val);
                }
                curFile.close();
            }
            else
            {
                curFile.close();
            }
        }
        resMap.insert(loc,resMISES);
    }
    return resMap;
}

//! --------------------------
//! function: setFatigueModel
//! details:
//! --------------------------
void postEngine::setFatigueModel(int fatigueAlgo)
{
    switch(fatigueAlgo)
    {
    case 0:
    {
        myFatigueModel.type = fatigueModel_BCM;
        myFatigueModel.coeffs<<0.6<<-0.69<<1407.64<<198609.5<<-0.065;
    }
        break;

    case 1:
    {
        myFatigueModel.type = fatigueModel_ESR;
        myFatigueModel.coeffs<<0<<1; // to DO
    }
        break;
    }
}

/*
double Rect(double x, double l, double r)
{
    if(x<l) return 0;
    if(x>=l && x<=r) return 1.0;
    return 0;
}

void coefficients(double x, double *c)
{
    c[0] = 195000;
    c[1] = Rect(x,103,248)*5.562508+Rect(x,248,4881)*1.55481e1;
    c[2] = Rect(x,103,248)*-1.014634e1+Rect(x,248,4881)*6.229821e-2;
    c[3] = Rect(x,103,248)*-5.738073e1+Rect(x,248,4881)*-8.425030e-2;
    c[4] = Rect(x,103,248)*7.152267e-1+Rect(x,248,4881)*-8.596020e-4;
    c[5] = Rect(x,103,248)*4.578432+Rect(x,248,4881)*1.029439e-4;
    c[6] = Rect(x,103,248)*3.584816e-3+Rect(x,248,4881)*8.030748e-6;
    c[7] = Rect(x,103,248)*0.0+Rect(x,248,4881)*1.603119e-5;
    c[8] = Rect(x,103,248)*0.0+Rect(x,248,4881)*5.051589e-9;
    c[9] = Rect(x,103,248)*0.0+Rect(x,248,4881)*-7.849028e-9;
    c[10]= Rect(x,103,248)*0.0+Rect(x,248,4881)*0.0;
    c[11]= Rect(x,103,248)*0.0+Rect(x,248,4881)*0.0;
}

//! ---------------------------------
//! function: evaulateFatigueResults
//! details:
//! ---------------------------------
postObject postEngine::evaluateFatigueResults(int type, QVector<GeometryTag> locs, const QList<double> &times, QMap<int,int> materialBodyMap, int nCycle)
{
    QMap<GeometryTag,QList<QMap<int,double>>> fatigueResults;

    switch(myFatigueModel.type)
    {
    case fatigueModel_BCM:
    {
        cout<<"postEngine::evaluateResult()->____fatigue model BCM called___"<<endl;

        QMap<GeometryTag,QMap<int,QList<double>>> r = readFatigueResults(type,locs,times);
        rainflow rf;
        for(QMap<GeometryTag,QMap<int,QList<double>>>::iterator it = r.begin(); it!=r.end(); ++it)
        {
            GeometryTag curLoc = it.key();
            rf.setLocation(curLoc);

            QMap<int,QList<double>> strainDistTimeHistory = it.value();
            rf.setFatigueModel(myFatigueModel);
            QMap<int,double> damageDist;

            bool isDone = rf.perform(strainDistTimeHistory,damageDist);
            if(isDone)
            {
                QList<QMap<int,double>>damageDistList;
                damageDistList<<damageDist;
                fatigueResults.insert(curLoc,damageDistList);
            }
        }
    }
        break;

    case fatigueModel_ESR:
    {
        cout<<"postEngine::evaluateResult()->____fatigue model ESR called___"<<endl;
        int step,substep;
        //QMap<double,QVector<int>> dtm= getDTM();
        //QMap<double,QVector<int>> dtm = myDTM;

        postTools::getStepSubStepByTimeDTM(myDTM,times.last(),step,substep);
        QString tor_eps = m.key(TypeOfResult_EPS);
        QString tor_mises = m.key(TypeOfResult_S);

        QMap<GeometryTag,QList<QMap<int,double>>> pe = this->evaluateResult(tor_eps,substep,step,locs);
        QMap<GeometryTag,QList<QMap<int,double>>> stress = this->evaluateResult(tor_mises,substep,step,locs);

        for(QVector<GeometryTag>::iterator it=locs.begin();it!=locs.end();it++)
        {
            GeometryTag curLoc = *it;
            int bodyIndex = curLoc.parentShapeNr;
            const QList<QMap<int, double>> &listOfResPe = pe.value(curLoc);
            const QMap<int, double> &curPe = listOfResPe.first();
            const QList<QMap<int, double>> &listOfResMises = stress.value(curLoc);
            const QMap<int, double> &curMises = listOfResMises.first();
            QList<QMap<int,double>> damageIndex;
            QMap<int,double> damageIndexData;

            int material = materialBodyMap.value(bodyIndex);

            for(QMap<int,double>::const_iterator itt=curPe.cbegin(); itt!=curPe.cend(); itt++)
            {
                int curPos = itt.key();
                double eps = itt.value();
                double mises = curMises.value(curPos);
                double N = 0;

                switch(material)
                {
                case 5:
                case 6:
                case 0:
                case 1:
                case 2:
                case 3:
                case 4:
                case 7:
                case 8:
                case 9:
                {
                    double elasticModulusAve = 176000;
                    double altStress = 0.5*(mises+eps*elasticModulusAve);
                    double c[12];
                    coefficients(altStress,c);
                    const double Cu = 6.894757;
                    double Y = (altStress/Cu)*(c[0]/elasticModulusAve);
                    double Xnum = c[1]+c[3]*pow(Y,1)+c[5]*pow(Y,2)+c[7]*pow(Y,3)+c[9]*pow(Y,4)+c[11]*pow(Y,5);
                    double Xden = 1+c[2]*pow(Y,1)+c[4]*pow(Y,2)+c[6]*pow(Y,3)+c[8]*pow(Y,4)+c[10]*pow(Y,5);
                    double X = Xnum/Xden;

                    N = pow(10.0,X);
                    if(altStress<103) N = 1e20;
                    if(altStress>4881) N = 1e20;
                }
                    break;
                }
                double damage = double(nCycle)/N;
                damageIndexData.insert(curPos,damage);
                //damageIndexData.insert(curPos,log10l(damage+1e-4));
            }
            damageIndex<<damageIndexData;
            fatigueResults.insert(curLoc,damageIndex);
        }
    }
        break;
    }

    //! -----------------------
    //! create the post object
    //! -----------------------
    QString label = this->colorBoxTitle("Damage",0,1,1);
    postObject aPostObject(fatigueResults,locs,label);
    //postObject aPostObject(fatigueResults,locs);
    return aPostObject;
}
*/
