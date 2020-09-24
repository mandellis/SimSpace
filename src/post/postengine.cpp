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

//! -------------------
//! function: buildMap
//! details:
//! -------------------
void postEngine::buildMap()
{
    m.insert("DISP",TypeOfResult_U);
    m.insert("STRESS",TypeOfResult_S);
    m.insert("TOSTRAIN",TypeOfResult_TOSTRAIN);
    m.insert("MESTRAIN",TypeOfResult_MESTRAIN);
    m.insert("NDTEMP",TypeOfResult_NT);
    m.insert("UDR",TypeOfResult_UDR);
    m.insert("FORC",TypeOfResult_F);
    m.insert("PE",TypeOfResult_EPS);
    m.insert("FLUX",TypeOfResult_HFL);
    m.insert("CONTACT",TypeOfResult_CONT);
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
    if(myMeshDataBase==Q_NULLPTR) return false;
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
std::map<GeometryTag,std::vector<std::map<int,double>>> postEngine::evaluateResult(const QString &resultKeyName,
                                                                                   int requiredSubStepNb,
                                                                                   int requiredStepNb,
                                                                                   int requiredMode,
                                                                                   const std::vector<GeometryTag> &vecLoc,
                                                                                   double &requiredTime)
{
    cout<<"@ -------------------------------------------------"<<endl;
    cout<<"@ - postEngine::evaluateResult "<<endl;
    cout<<"@ - variable "<<resultKeyName.toStdString()<<" on "<<vecLoc.size()<<" locations"<<endl;
    cout<<"@ -------------------------------------------------"<<endl;

    //! ----------------------------------------------------
    //! generate the results on all the requested locations
    //! ----------------------------------------------------
    std::map<GeometryTag,std::vector<std::map<int,double>>> resMap;
    for(std::vector<GeometryTag>::const_iterator it = vecLoc.cbegin(); it!= vecLoc.cend(); ++it)
    {
        GeometryTag loc = *it;

        //! -------------------------------------------------------------------------
        //! node conversion map: (Calculix mesh nodeID,nodeID for MeshVS_DataSource)
        //! -------------------------------------------------------------------------
        std::map<int,int> indexedMapOfNodes = OCCMeshToCCXmesh::performCCXtoOCC(loc,myMeshDataBase);

        //! ------------------------------------------------
        //! enter <...>/SolutionData/ResultsData
        //! ------------------------------------------------
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

        std::vector<std::map<int,double>> res;

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
            int mode;
            sscanf(val.c_str(),"%s%d",&analysisType,&mode);

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

            bool eval=false;
            switch (requiredMode)
            {
            case 0:
            {
                if(strcmp(tdata,resultKeyName.toStdString().c_str())==0 && subStepNb==requiredSubStepNb && stepNb == requiredStepNb && mode == requiredMode)
                {
                    requiredTime = time;
                    eval = true;
                }
            }
                break;
            default:
            {
                if(strcmp(tdata,resultKeyName.toStdString().c_str())==0 && mode == requiredMode)
                {
                    requiredTime = time;
                    eval = true;
                }
            }
                break;
            }
            if(eval)
            {
                //printf("file @ required time found\n");
                TypeOfResult tor = m.value(tdata);
                switch(tor)
                {
                case TypeOfResult_HFL:
                {
                    std::map<int,double> resComp_normal;
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
                        std::map<int,int>::iterator it = indexedMapOfNodes.find(ni);

                        if(it!=indexedMapOfNodes.end())
                        {
                            int OCCnodeID = it->second;
                            const QList<double> &normal = curFaceDS->myNodeNormals.value(OCCnodeID);
                            double normalFlux = cxx*normal[0]+cyy*normal[1]+czz*normal[2];
                            resComp_normal.insert(std::make_pair(OCCnodeID,normalFlux));
                        }
                        std::getline(curFile,val);
                    }

                    //! result
                    res.push_back(resComp_normal);
                    //cout<<"postEngine::evaluateResult()->____Number of components: "<<res.length()<<"____"<<endl;
                }
                    break;
                case TypeOfResult_U:
                case TypeOfResult_F:
                //case TypeOfResult_HFL:
                {
                    std::map<int,double> resComp_X,resComp_Y,resComp_Z,resComp_Total;

                    //! <>::eof(): call getline before while, then inside {}, @ as last instruction
                    std::getline(curFile,val);
                    while(curFile.eof()!=true)
                    {
                        int ni;
                        double cxx,cyy,czz,total;
                        sscanf(val.c_str(),"%d%lf%lf%lf",&ni,&cxx,&cyy,&czz);

                        //! nodeIDs defining the MeshVS_dataSource
                        std::map<int,int>::iterator it = indexedMapOfNodes.find(ni);
                        if(it!=indexedMapOfNodes.end())
                        {
                            int OCCnodeID = it->second;
                            total = sqrt(pow(cxx,2)+pow(cyy,2)+pow(czz,2));
                            resComp_Total.insert(std::make_pair(OCCnodeID,total));
                            resComp_X.insert(std::make_pair(OCCnodeID,cxx));
                            resComp_Y.insert(std::make_pair(OCCnodeID,cyy));
                            resComp_Z.insert(std::make_pair(OCCnodeID,czz));

                        }
                        std::getline(curFile,val);
                    }

                    //! result
                    res.push_back(resComp_Total);
                    res.push_back(resComp_X);
                    res.push_back(resComp_Y);
                    res.push_back(resComp_Z);

                    //cout<<"postEngine::evaluateResult()->____Number of components: "<<res.length()<<"____"<<endl;
                }
                    break;

                case TypeOfResult_S:
                case TypeOfResult_TOSTRAIN:
                case TypeOfResult_MESTRAIN:
                {
                    //!                   0       1        2      3       4        5       6       7       8      9      10
                    std::map<int,double> resMISES, resSINT, resSI, resSII, resSIII, resSXX, resSYY, resSZZ, resSXY,resSYZ,resSXZ;

                    //! <>::eof(): call getline before while, then inside {}, @ as last instruction
                    std::getline(curFile,val);
                    while(curFile.eof()!=true)
                    {
                        //! read the components of the 3x3 data
                        int ni;
                        double cxx,cyy,czz,cxy,cyz,cxz;
                        sscanf(val.c_str(),"%d%lf%lf%lf%lf%lf%lf",&ni,&cxx,&cyy,&czz,&cxy,&cyz,&cxz);

                        //! nodeIDs defining the MeshVS_dataSource
                        std::map<int,int>::iterator it = indexedMapOfNodes.find(ni);
                        if(it!=indexedMapOfNodes.end())
                        {
                            int OCCnodeID = it->second;

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
                            //resMISES.insert(OCCnodeID,vonMises);
                            resMISES.insert(std::make_pair(OCCnodeID,vonMises));

                            //! ---------------------------------
                            //! compute the principal components
                            //! ---------------------------------
                            double sik[6] {cxx,cyy,czz,cxy,cyz,cxz};
                            double s[3];
                            postTools::principalComponents(sik,s);

                            resSI.insert(std::make_pair(OCCnodeID,s[2]));        //! maximum
                            resSII.insert(std::make_pair(OCCnodeID,s[1]));       //! middle
                            resSIII.insert(std::make_pair(OCCnodeID,s[0]));      //! minimum

                            //! -----------------------------------------------
                            //! compute the stress/strain intensity (2*Tresca)
                            //! maximum shear stress/strain
                            //! -----------------------------------------------
                            double sint = fabs(s[2]-s[0]);
                            resSINT.insert(std::make_pair(OCCnodeID,sint));
                            resSXX.insert(std::make_pair(OCCnodeID,cxx));
                            resSYY.insert(std::make_pair(OCCnodeID,cyy));
                            resSZZ.insert(std::make_pair(OCCnodeID,czz));
                            resSXY.insert(std::make_pair(OCCnodeID,cxy));
                            resSYZ.insert(std::make_pair(OCCnodeID,cyz));
                            resSXZ.insert(std::make_pair(OCCnodeID,cxz));

                        }
                        std::getline(curFile,val);
                    }
                    //! result
                    res.push_back(resMISES);
                    res.push_back(resSINT);
                    res.push_back(resSI);
                    res.push_back(resSII);
                    res.push_back(resSIII);
                    res.push_back(resSXX);
                    res.push_back(resSYY);
                    res.push_back(resSZZ);
                    res.push_back(resSXY);
                    res.push_back(resSYZ);
                    res.push_back(resSXZ);
                    //cout<<"postEngine::evaluateResult()->____Number of components: "<<res.length()<<"____"<<endl;;
                }
                    break;

                case TypeOfResult_NT:
                case TypeOfResult_EPS:
                {
                    std::map<int,double> resT;

                    //! <>::eof(): call getline before while, then inside {}, @ as last instruction
                    std::getline(curFile,val);
                    while(curFile.eof()!=true)
                    {
                        int ni;
                        double v;
                        sscanf(val.c_str(),"%d%lf",&ni,&v);

                        //! nodeIDs defining the MeshVS_dataSource
                        std::map<int,int>::iterator it = indexedMapOfNodes.find(ni);
                        if(it!=indexedMapOfNodes.end())
                        {
                            int OCCnodeID = it->second;
                            //resT.insert(OCCnodeID,v);
                            resT.insert(std::make_pair(OCCnodeID,v));
                        }
                        std::getline(curFile,val);
                    }
                    //! result
                    res.push_back(resT);
                    //cout<<"Number of components: "<<res.length();
                }
                    break;

                case TypeOfResult_CONT:
                {
                    //!                         col 1           col 2+3     col 4           col 5+6
                    std::map<int,double> resContPenetration, resContSliding, resContPress, resContFrictStress;

                    //! <>::eof(): call getline before while, then inside {}, @ as last instruction
                    std::getline(curFile,val);
                    while(curFile.eof()!=true)
                    {
                        //! read the components of the 3x3 data
                        int ni;
                        double cxx,cyy,czz,cxy,cyz,cxz;
                        sscanf(val.c_str(),"%d%lf%lf%lf%lf%lf%lf",&ni,&cxx,&cyy,&czz,&cxy,&cyz,&cxz);

                        //! nodeIDs defining the MeshVS_dataSource
                        std::map<int,int>::iterator it = indexedMapOfNodes.find(ni);
                        if(it!=indexedMapOfNodes.end())
                        {
                            int OCCnodeID = it->second;

                            //! --------------------------------------------------
                            //! compute the frictional stress and contact sliding
                            //! -------------------------------------------------
                            double frictStress = sqrt(cyz*cyz+cxz*cxz);
                            double contSliding = sqrt(cyy*cyy+czz*czz);

                            resContFrictStress.insert(std::make_pair(OCCnodeID,frictStress));
                            resContSliding.insert(std::make_pair(OCCnodeID,contSliding));
                            resContPress.insert(std::make_pair(OCCnodeID,cxy));
                            resContPenetration.insert(std::make_pair(OCCnodeID,cxx));
                        }
                        std::getline(curFile,val);
                    }
                    //! result
                    res.push_back(resContPress);
                    res.push_back(resContFrictStress);
                    res.push_back(resContPenetration);
                    res.push_back(resContSliding);
                }
                    break;
                }
                resMap.insert(std::make_pair(loc,res));
                curFile.close();
                break;
            }
            else
            {
                curFile.close();
            }
            //curFile.close();
        }
        //! --------------------------------------------------------------------
        //! resMap => std::vector<GeometryTag,std::vector<std::map<int,double>>
        //! res => std::vector<std::map<int,double>>
        //! --------------------------------------------------------------------
    }
    return resMap;
}

//! ---------------------------------
//! function: resultName
//! details:  title for the colorbox
//! ---------------------------------
QString postEngine::resultName(const QString &keyName, int component, int step, int subStep, double time)
{
    //QString timeInfo = QString("\nTime %1\nStep %2").arg(time).arg(step);
    //QString timeInfo = QString("\nStep %1\nSub Step %2").arg(step).arg(subStep);
    QString timeInfo = QString("\nTime %1\nStep %2\nSub Step %3").arg(time).arg(step).arg(subStep);

    QString resultName;

    if(keyName =="Damage")
    {
        resultName = keyName;
        return resultName.append("\n");
    }

    TypeOfResult tor = m.value(keyName);
    switch(tor)
    {
    case TypeOfResult_HFL:
        switch(component)
        {
        case 0: resultName="Thermal flux"; break;
        }
        break;
    case TypeOfResult_CONT:
        switch(component)
        {
        case 0: resultName="Contact pressure"; break;
        case 1: resultName="Frictional stress"; break;
        case 2: resultName="Contact penetration"; break;
        case 3: resultName="Contact sliding"; break;
        }
        break;
    case TypeOfResult_U:
        switch(component)
        {
        case 0: resultName="Total displacement"; break;
        case 1: resultName="Directional displacement X"; break;
        case 2: resultName="Directional displacement Y"; break;
        case 3: resultName="Directional displacement Z"; break;
        }
        break;
    case TypeOfResult_S:
        switch(component)
        {
        case 0: resultName="Equivalent stress"; break;
        case 1: resultName = "Stress intensity"; break;
        case 2: resultName = "Maximum principal stress"; break;
        case 3: resultName = "Middle principal stress"; break;
        case 4: resultName = "Minimum principal stress"; break;
        case 5: resultName = "Normal stress X"; break;
        case 6: resultName = "Normal stress Y"; break;
        case 7: resultName = "Normal stres Z"; break;
        case 9: resultName = "Shear stress XY"; break;
        case 10: resultName = "Shear stress YZ"; break;
        case 11: resultName = "Shear stress ZX"; break;
        }
        break;
    case TypeOfResult_TOSTRAIN:
        switch(component)
        {
        case 0: resultName="Equivalent strain"; break;
        case 1: resultName = "Strain intensity"; break;
        case 2: resultName = "Maximum principal strain"; break;
        case 3: resultName = "Middle principal strain"; break;
        case 4: resultName = "Minimum principal strain"; break;
        }
        break;
    case TypeOfResult_MESTRAIN:
        switch(component)
        {
        case 0: resultName="Equivalent strain"; break;
        case 1: resultName = "Strain intensity"; break;
        case 2: resultName = "Maximum principal strain"; break;
        case 3: resultName = "Middle principal strain"; break;
        case 4: resultName = "Minimum principal strain"; break;
        }
        break;
    case TypeOfResult_NT:
        switch(component)
        {
        case 0: resultName = "Temperature"; break;
        }
        break;
    case TypeOfResult_F:
        switch(component)
        {
        case 0: resultName="Total force"; break;
        case 1: resultName="Directional Force X"; break;
        case 2: resultName="Directional Force Y"; break;
        case 3: resultName="Directional Force Z"; break;
        }
        break;
    case TypeOfResult_EPS:
        switch(component)
        {
        case 0: resultName="Equivalent Plastic Strain"; break;
        }
        break;
    default:
        resultName = "Unnamed result"; break;
    }
    return this->timeStamp().append("\n").append(resultName).append(timeInfo).append("\n");
}

//! --------------------------------------------------------------------
//! function: buildPostObject
//! details:  this method reads the data from the .frd file (from disk)
//! --------------------------------------------------------------------
bool postEngine::buildPostObject(const QString &keyName,
                                 int component,
                                 int requiredSubStepNb,
                                 int requiredStepNb,
                                 int requiredMode,
                                 const std::vector<GeometryTag> &vecLoc,
                                 sharedPostObject &aPostObject)
{
    //! --------------------
    //! call the postEngine
    //! --------------------
    double time;
    const std::map<GeometryTag,std::vector<std::map<int,double>>> &resMap = this->evaluateResult(keyName, requiredSubStepNb, requiredStepNb, requiredMode, vecLoc, time);

    //! -------------------------
    //! build the colorBox title
    //! -------------------------
    QString aResultName = this->resultName(keyName, component, requiredStepNb, requiredSubStepNb, time);

    //! --------------------------------------------------------------------------------------------------------------------------
    //! create the map of nodal vectorial displacements for the deformed mesh presentation. Here:
    //! std::map<int,gp_Vec> displMap                    => map of nodal vectorial displacements
    //! std::map<GeometryTag,std:.vector<std::map<int,gp_Vec>>> => each location has its own map of nodal vectorial displacements
    //! --------------------------------------------------------------------------------------------------------------------------
    std::map<int,gp_Vec> displMap;
    std::map<GeometryTag,std::map<int,gp_Vec>> mapDisplMap;
    const std::map<GeometryTag,std::vector<std::map<int,double>>> &nodalDisplacements = this->evaluateResult("DISP", requiredSubStepNb, requiredStepNb,requiredMode, vecLoc, time);
    for(std::map<GeometryTag,std::vector<std::map<int,double>>>::const_iterator it = nodalDisplacements.cbegin(); it!=nodalDisplacements.cend(); ++it)
    {
        const GeometryTag &aLoc= it->first;

        const std::vector<std::map<int,double>> &nodalDisplacementsComponents = it->second;
        const std::map<int,double> &displX = nodalDisplacementsComponents[1];
        const std::map<int,double> &displY = nodalDisplacementsComponents[2];
        const std::map<int,double> &displZ = nodalDisplacementsComponents[3];

        std::map<int,double>::const_iterator itX = displX.cbegin();
        std::map<int,double>::const_iterator itY = displY.cbegin();
        std::map<int,double>::const_iterator itZ = displZ.cbegin();
        for(;itX!=displX.cend() && itY!=displY.cend() && itZ!=displZ.cend(); ++itX, ++itY, ++itZ)
        {
            int nodeID = itX->first;
            gp_Vec aVec(itX->second,itY->second,itZ->second);
            displMap[nodeID] = aVec;        // do not use "insert"
        }
        mapDisplMap[aLoc]=displMap;         // do not use "insert"
    }

    //! ----------------------
    //! create the postObject
    //! ----------------------
    bool useSurfaceMeshForVolumeResults = Global::status().myResultPresentation.useExteriorMeshForVolumeResults;
    aPostObject = std::make_shared<postObject>(resMap,vecLoc,mapDisplMap,aResultName,useSurfaceMeshForVolumeResults);
    aPostObject->init(myMeshDataBase);
    double magnifyFactor = Global::status().myResultPresentation.theScale;
    bool isDone = aPostObject->buildMeshIO(-1,-1,10,true,component,magnifyFactor);
    return isDone;
}

//! --------------------
//! function: timeStamp
//! details:
//! --------------------
QString postEngine::timeStamp()
{
    QDateTime dateTime;
    QString dateFormat = "dd/MM/yyyy";
    QString dateString = dateTime.currentDateTime().toString(dateFormat);
    QString timeStamp;
    timeStamp.append("Date: ").append(dateString);
    return timeStamp;
}

//! --------------------------
//! function: updateIsostrips
//! details:
//! --------------------------
void postEngine::updateIsostrips(sharedPostObject &aPostObject, int scaleType, double minValue, double maxValue, int NbIntervals)
{
    cout<<"postEngine::updateIsostrips()->____function called____"<<endl;

    //! ---------------------------------------------------------------
    //! retrieve the mesh data sources and the solution data component
    //! ---------------------------------------------------------------
    int solutionDataComponent = aPostObject->getSolutionDataComponent();

    //! -----------------------------------------
    //! update the colored mesh and the colorbox
    //! the scale is left unchanged
    //! -----------------------------------------
    switch (scaleType)
    {
    case 0:
        //! autoscale min max, custom number of levels
        aPostObject->buildMeshIO(-1,-1,NbIntervals,true,solutionDataComponent);
        break;
    case 1:
        //! custom scale (custom min, max, number of levels)
        aPostObject->buildMeshIO(minValue,maxValue,NbIntervals,false,solutionDataComponent);
        break;
    }
}

//! ---------------------------------
//! function: evaulateFatigueResults
//! details:  overload
//! ---------------------------------
bool postEngine::evaluateFatigueResults(int type, std::vector<GeometryTag> locs, const QList<double> &times, QMap<int,int> materialBodyMap, int nCycle, sharedPostObject &aPostObject)
{
    std::map<GeometryTag,std::vector<std::map<int,double>>> fatigueResults;

    switch(myFatigueModel.type)
    {
    case fatigueModel_BCM:
    {
        cout<<"postEngine::evaluateFatigueResult()->____fatigue model BCM called___"<<endl;
        std::map<GeometryTag,std::map<int,QList<double>>> r = readFatigueResults(type,locs,times);

        rainflow rf;
        for(std::map<GeometryTag,std::map<int,QList<double>>>::iterator it = r.begin(); it!=r.end(); ++it)
        {
            GeometryTag curLoc = it->first;

            rf.setLocation(curLoc);

            std::map<int,QList<double>> strainDistTimeHistory = it->second;

            rf.setFatigueModel(myFatigueModel);
            std::map<int,double> damageDist;

            bool isDone = rf.perform(strainDistTimeHistory,damageDist);
            if(isDone)
            {
                std::vector<std::map<int,double>> damageDistList;

                damageDistList.push_back(damageDist);
                fatigueResults.insert(std::make_pair(curLoc,damageDistList));
            }
        }
    }
        break;

    case fatigueModel_ESR:
    {
        cout<<"postEngine::evaluateFatigueResult()->____fatigue model ESR called___"<<endl;
        int step,substep;
        double requiredTime;
        postTools::getStepSubStepByTimeDTM(myDTM,times.last(),step,substep);
        QString tor_eps = m.key(TypeOfResult_EPS);
        QString tor_mises = m.key(TypeOfResult_S);
        int mode =0;

        std::map<GeometryTag,std::vector<std::map<int,double>>> pe = this->evaluateResult(tor_eps,substep,step,mode,locs,requiredTime);
        std::map<GeometryTag,std::vector<std::map<int,double>>> stress = this->evaluateResult(tor_mises,substep,step,mode,locs,requiredTime);

        for(std::vector<GeometryTag>::iterator it=locs.begin();it!=locs.end();it++)
        {
            GeometryTag curLoc = *it;
            int bodyIndex = curLoc.parentShapeNr;
            const std::vector<std::map<int, double>> &listOfResPe = pe.at(curLoc);
            const std::map<int, double> &curPe = listOfResPe[0];
            const std::vector<std::map<int, double>> &listOfResMises = stress.at(curLoc);
            const std::map<int, double> &curMises = listOfResMises[0];
            std::vector<std::map<int,double>> damageIndex;
            std::map<int,double> damageIndexData;

            double elasticModulusAve, elasticModulusMin,r,a,b,c,d,e,f,g,h;
            Q_UNUSED (elasticModulusMin)
            int material = materialBodyMap.value(bodyIndex);

            for(std::map<int,double>::const_iterator itt=curPe.cbegin(); itt!=curPe.cend(); itt++)
            {
                double altStress,X,Y;
                int curPos = itt->first;
                double eps = itt->second;
                double mises = curMises.at(curPos);

                switch(material)
                {
                case 5:case 6:case 0:case 1:case 2:case 3:case 4:case 7:case 8:case 9:
                {
                    elasticModulusAve = 1.76e5;
                    altStress = 0.5*(mises+eps*elasticModulusAve);
                    Y = log10(28300*altStress/elasticModulusAve);

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
                damageIndexData.insert(std::make_pair(curPos,damage));
            }
            damageIndex.push_back(damageIndexData);
            fatigueResults.insert(std::make_pair(curLoc,damageIndex));
        }
    }
        break;
    }

    //! -----------------------
    //! create the post object
    //! -----------------------
    bool useSurfaceMeshForVolumeResults = Global::status().myResultPresentation.useExteriorMeshForVolumeResults;
    QString label = this->resultName("Damage",0,1,1,0);
    std::map<GeometryTag,std::map<int,gp_Vec>> mapDisplMap;
    aPostObject = std::make_shared<postObject>(fatigueResults,locs,mapDisplMap,label,useSurfaceMeshForVolumeResults);
    aPostObject->init(myMeshDataBase);
    double magnifyFactor = Global::status().myResultPresentation.theScale;
    int component = 0;
    bool isDone = aPostObject->buildMeshIO(-1,-1,10,true,component,magnifyFactor);
    return isDone;
}

//! ---------------------------------------------------
//! function: readFatigueResults
//! details:  type = 1 => equivalent mechanical strain
//!           type = 0 => equivalent total strain
//!           ... other types
//! ---------------------------------------------------
std::map<GeometryTag,std::map<int,QList<double>>> postEngine::readFatigueResults(int type,
                                                                                 const std::vector<GeometryTag> &vecLoc,
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
    std::map<GeometryTag,std::map<int,QList<double>>> resMap;

    std::vector<GeometryTag>::const_iterator it;
    for(it = vecLoc.cbegin(); it!= vecLoc.cend(); ++it)
    {
        //! -------------------------------------------------------------------------
        //! node conversion map: (Calculix mesh nodeID,nodeID for MeshVS_DataSource)
        //! -------------------------------------------------------------------------
        GeometryTag loc = *it;
        std::map<int,int> indexedMapOfNodes = OCCMeshToCCXmesh::performCCXtoOCC(loc,myMeshDataBase);

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
        //QMap<int,QList<double>> resMISES;
        std::map<int,QList<double>> resMISES;

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
                    std::map<int,int>::iterator it = indexedMapOfNodes.find(ni);
                    if(it!=indexedMapOfNodes.end() && n==0)
                    //if(OCCnodeID!=-1 && n==0)
                    {
                        int OCCnodeID = it->second;

                        //! -------------------------------------
                        //! compute the equivalent stress/strain
                        //! -------------------------------------
                        double vonMises = (2.0/3.0)*sqrt((3.0/2.0)*(cxx*cxx+cyy*cyy+czz*czz)+(3.0/4.0)*(cxy*cxy+cyz*cyz+cxz*cxz));
                        timeHistory<<vonMises;
                        //resMISES.insert(OCCnodeID,timeHistory);
                        resMISES.insert(std::make_pair(OCCnodeID,timeHistory));

                        //! -----------------------------------------------------
                        //! compute the principal components: index 0 is minimum
                        //! -----------------------------------------------------
                        double sik[6] {cxx,cyy,czz,cxy,cyz,cxz};
                        double s[3];
                        postTools::principalComponents(sik,s);

                    }
                    else if(it!=indexedMapOfNodes.end())
                    //else if(OCCnodeID!=-1)
                    {
                        int OCCnodeID = indexedMapOfNodes.find(ni)->second;

                        //! -------------------------------------
                        //! compute the equivalent stress/strain
                        //! -------------------------------------
                        double vonMises;
                        vonMises = (2.0/3.0)*sqrt((3.0/2.0)*(cxx*cxx+cyy*cyy+czz*czz)+(3.0/4.0)*(cxy*cxy+cyz*cyz+cxz*cxz));
                        //timeHistory = resMISES.value(OCCnodeID);
                        timeHistory = resMISES.at(OCCnodeID);

                        timeHistory<<vonMises;
                        //resMISES.insert(OCCnodeID,timeHistory);
                        resMISES.insert(std::make_pair(OCCnodeID,timeHistory));

                        //! ---------------------------------
                        //! compute the principal components
                        //! ---------------------------------
                        double sik[6] {cxx,cyy,czz,cxy,cyz,cxz};
                        double s[3];
                        postTools::principalComponents(sik,s);
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
        resMap.insert(std::make_pair(loc,resMISES));
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

//! ------------------------------------
//! function: updateResultsPresentation
//! details:
//! ------------------------------------
void postEngine::updateResultsPresentation(QList<sharedPostObject> &postObjectList)
{
    static resultPresentation previousResultPresentation;
    resultPresentation newResultsPresentation = Global::status().myResultPresentation;
    for(QList<sharedPostObject>::iterator it = postObjectList.begin(); it!=postObjectList.end(); ++it)
    {
        sharedPostObject aPostObject = *it;

        //! ----------------------------------------------------------------------------
        //! the surface mesh=>volume mesh/volume mesh=>surface mesh for viewing results
        //! the presentation should be fully rebuilt
        //! ----------------------------------------------------------------------------
        /* commented: heavy computations (isostrips, isosurfaces, isolines should not be
         * re-done when only some features of the results presentation has been changed
         * left here for documentation
        if(previousResultPresentation != newResultsPresentation)
        {
            double min = aPostObject->getMin();
            double max = aPostObject->getMax();
            int NbLevels = aPostObject->getNbLevels();
            int component = aPostObject->getSolutionDataComponent();
            bool isAutoScale = aPostObject->IsAutoscale();
            int magnifyFactor = newResultsPresentation.theScale;

            aPostObject->setMode(newResultsPresentation.useExteriorMeshForVolumeResults);
            aPostObject->buildMeshIO(min,max,NbLevels,isAutoScale,component,magnifyFactor);
        }
        */
        bool toBeUpdated = false;
        if(previousResultPresentation.theTypeOfPresentation != newResultsPresentation.theTypeOfPresentation) toBeUpdated = true;
        if(previousResultPresentation.theScale != newResultsPresentation.theScale) toBeUpdated = true;
        if(toBeUpdated)
        {
            double min = aPostObject->getMin();
            double max = aPostObject->getMax();
            int NbLevels = aPostObject->getNbLevels();
            int component = aPostObject->getSolutionDataComponent();
            bool isAutoScale = aPostObject->IsAutoscale();
            int magnifyFactor = newResultsPresentation.theScale;

            aPostObject->setMode(newResultsPresentation.useExteriorMeshForVolumeResults);
            aPostObject->buildMeshIO(min,max,NbLevels,isAutoScale,component,magnifyFactor);
        }
    }
    previousResultPresentation = newResultsPresentation;
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
postObject postEngine::evaluateFatigueResults(int type, std::vector<GeometryTag> locs, const QList<double> &times, QMap<int,int> materialBodyMap, int nCycle)
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

        for(std::vector<GeometryTag>::iterator it=locs.begin();it!=locs.end();it++)
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
    QString label = this->resultName("Damage",0,1,1);
    postObject aPostObject(fatigueResults,locs,label);
    //postObject aPostObject(fatigueResults,locs);
    return aPostObject;
}
*/
