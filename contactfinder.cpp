//! -------
//! global
//! -------
#include "global.h"

//! ----------------
//! custom includes
//! ----------------
#include "contactFinder.h"
#include "occhandle.h"
#include <mesh.h>
#include <ng_meshvs_datasource3d.h>
#include <ng_meshvs_datasource2d.h>
#include <ng_meshvs_datasourceface.h>
#include <polygon.h>
#include "qprogressindicator.h"
#include "qprogressevent.h"
#include "extendedrwstl.h"
#include "stlapiwriter.h"
#include "meshtools.h"
#include "triangletotriangleintersection.h"
#include "hash_c.h"

//! --------------
//! compatibility
//! --------------
#include <StlTransfer.hxx>
#include <StlMesh_Mesh.hxx>

//! ----
//! OCC
//! ----
#include <MeshVS_DataSource.hxx>
#include <TColStd_MapIteratorOfPackedMapOfInteger.hxx>
#include <TColStd_Array1OfInteger.hxx>
#include <BRepTools.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <NCollection_Array1.hxx>
#include <Bnd_Box.hxx>
#include <BRepBndLib.hxx>
#include <TopoDS_Face.hxx>
#include <TopAbs.hxx>
//! ---
//! Qt
//! ---
#include <QList>
#include <QWidget>
#include <QApplication>
#include <QThread>

//! ----
//! C++
//! ----
#include <set>
#include <map>
#include <vector>
#include <algorithm>
#include <chrono>

//! ----
//! PQP
//! ----
#include "TriDist.h"

//! ----------------------
//! function: constructor
//! details:  default
//! ----------------------
contactFinder::contactFinder():myAngularCriterion(30.0)
{
    ;
}

//! ------------------------------
//! function: setAngularCriterion
//! details:
//! ------------------------------
void contactFinder::setAngularCriterion(double angle)
{
    if(angle>0 && angle<90) myAngularCriterion = angle;
}

//! -------------------------------
//! function: setProgressIndicator
//! details:
//! -------------------------------
void contactFinder::setProgressIndicator(QProgressIndicator *aProgressIndicator)
{
    if(aProgressIndicator!=Q_NULLPTR) myProgressIndicator = aProgressIndicator;
}

//! ----------------------
//! function: setDataBase
//! details:
//! ----------------------
void contactFinder::setDataBase(meshDataBase *mDB)
{
    myMDB = mDB;
}

//! -----------------------
//! function: setBodyPairs
//! details:
//! -----------------------
void contactFinder::setBodyPairs(const std::vector<std::pair<GeometryTag,GeometryTag>> &vectorOfBodyPairs)
{
    myVectorOfBodyPairs = vectorOfBodyPairs;
}

//! ----------------------
//! function: findFaceTag
//! details:
//! ----------------------
GeometryTag contactFinder::findFaceTag(const TopoDS_Face &aFace)
{
    if(myMDB==NULL)
    {
        cout<<"contactFinder::findFaceTag()->____cannot retrieve the face tag: the mesh data base has not been set____"<<endl;
        return GeometryTag();
    }
    GeometryTag aGeometryTag;
    for(QMap<int,TopoDS_Shape>::iterator it = myMDB->bodyMap.begin(); it!=myMDB->bodyMap.end(); ++it)
    {
        int bodyIndex = it.key();
        const TopTools_IndexedMapOfShape &faceMap = myMDB->MapOfBodyTopologyMap.value(bodyIndex).faceMap;
        if(faceMap.Contains(aFace))
        {
            int faceNr = faceMap.FindIndex(aFace);
            aGeometryTag.parentShapeNr = bodyIndex;
            aGeometryTag.subTopNr = faceNr;
            aGeometryTag.subShapeType = TopAbs_FACE;
            aGeometryTag.isParent = false;
            break;
        }
    }
    return aGeometryTag;
}

//! ----------------------------------------------------------------
//! function: AngularCriterion
//! details:  input a vector in 3D space with components (xA,yA,zA)
//!           input a vector in 3D space with components (xB,yB,zB)
//! ----------------------------------------------------------------
bool contactFinder::AngularCriterion(double xA, double yA, double zA, double xB, double yB, double zB)
{
    const double pi = 3.141592654;
    double dot = xA*xB+yA*yB+zA*zB;
    if(dot>1.0)
    {
        dot =1.0;
    }
    if(dot<-1.0)
    {
        dot =-1;
    }
    double angle = 180.0*acos(dot)/pi;
    //if(myAngularCriterion+angle<=180 && angle>=180-myAngularCriterion)
    if(angle<=180 && angle>=180-myAngularCriterion)
    {
        return true;
    }
    else
    {
        return false;
    }
}

//! -----------------------
//! function: triangleHash
//! details:
//! -----------------------
std::size_t contactFinder::triangleHash(const triangle &aTriangle)
{
    std::size_t seed = 0;
    for(int i=0; i<3; i++)
    {
        hash_c<double>(seed,aTriangle[i].x);
        hash_c<double>(seed,aTriangle[i].y);
        hash_c<double>(seed,aTriangle[i].z);
    }
    return seed;
}

//! ------------------
//! function: perform
//! details:
//! ------------------
bool contactFinder::perform(const std::vector<std::pair<GeometryTag,GeometryTag>> &vectorOfBodyPairs,
                            std::vector<std::pair<std::vector<GeometryTag>,std::vector<GeometryTag>>> &allContactsPairs,
                            double tolerance,
                            int grouping)
{
    //! -----------------
    //! has been stopped
    //! -----------------
    bool stopped = false;

    //! ----------------
    //! disable "Abort"
    //! ----------------
    myProgressIndicator->disableStop();

    //! ------------------------
    //! set the status: running
    //! ------------------------
    Global::status().code = 1;

    //! ----------------------------------------------------------
    //! create an init event: let the user read what is happening
    //! ----------------------------------------------------------
    int NbPairs = int(vectorOfBodyPairs.size());
    QString task = "Creating automatic connections";
    QProgressEvent *aProgressEvent = new QProgressEvent(QProgressEvent_Init,0,NbPairs-1,0,"Searching for contact pairs",
                                                   QProgressEvent_Init,0,100,0,task);
    QApplication::postEvent(myProgressIndicator,aProgressEvent);
    QApplication::processEvents();
    QThread::msleep(500);

    //! -----------------------------
    //! prepare the STL tessellation
    //! -----------------------------
    std::map<int,occHandle(MeshVS_DataSource)> mapOfSurfaceMeshDataSources;
    std::map<int,std::map<int,occHandle(MeshVS_DataSource)>> mapOfMapOfFaceMeshDataSources;

    bool rebuildTessellations = true;
    this->prepareSTLTessellation_useDisk(vectorOfBodyPairs,rebuildTessellations,mapOfSurfaceMeshDataSources,mapOfMapOfFaceMeshDataSources);

    //! ---------------------------------------------------
    //! build the maps <bodyIndex, <triangleHash, faceNr>>
    //! ---------------------------------------------------
    std::map<int,std::map<size_t,int>> trianglesMap;
    for(std::map<int,std::map<int,occHandle(MeshVS_DataSource)>>::iterator it = mapOfMapOfFaceMeshDataSources.begin(); it!=mapOfMapOfFaceMeshDataSources.end(); ++it)
    {
        std::pair<int,std::map<int,occHandle(MeshVS_DataSource)>> aPair = *it;
        int curBodyIndex = aPair.first;
        std::map<int,occHandle(MeshVS_DataSource)> faceMeshDSOfBody = aPair.second;

        //! -----------------------------------------------------------
        //! iterate over the face mesh data source of the current body
        //! -----------------------------------------------------------
        std::map<size_t,int> faceTriangleMap;
        for(std::map<int,occHandle(MeshVS_DataSource)>::iterator it_ = faceMeshDSOfBody.begin(); it_!=faceMeshDSOfBody.end(); ++it_)
        {
            std::pair<int,occHandle(MeshVS_DataSource)> apair_ =*it_;
            int curFaceNr = apair_.first;
            occHandle(MeshVS_DataSource) curFaceMeshDS = apair_.second;

            //! ----------------------------------------------------------------
            //! iterate over the triangles of the current face mesh data source
            //! ----------------------------------------------------------------
            for(TColStd_MapIteratorOfPackedMapOfInteger eit(curFaceMeshDS->GetAllElements()); eit.More(); eit.Next())
            {
                int globalElementID = eit.Key();
                MeshVS_EntityType aType;
                int NbNodes;
                double buf[9];
                TColStd_Array1OfReal coords(*buf,1,9);
                curFaceMeshDS->GetGeom(globalElementID,true,coords,NbNodes,aType);

                //! ------------------------
                //! hash code of a triangle
                //! ------------------------
                size_t seed = 0;
                for(int n=0; n<NbNodes; n++)
                {
                    int s = 3*n;
                    hash_c<double>(seed,coords(s+1));
                    hash_c<double>(seed,coords(s+2));
                    hash_c<double>(seed,coords(s+3));
                }
                std::pair<size_t,int> mapElement;
                mapElement.first = seed;
                mapElement.second = curFaceNr;
                faceTriangleMap.insert(mapElement);
            }
        }
        std::pair<int,std::map<size_t,int>> outerElementMap;
        outerElementMap.first = curBodyIndex;
        outerElementMap.second = faceTriangleMap;
        trianglesMap.insert(outerElementMap);
    }

    //! --------------------------------------
    //! iterate over the vector of body pairs
    //! --------------------------------------
    int done = 0;
    double TotalTimeForSearch = 0;
    for(std::vector<std::pair<GeometryTag,GeometryTag>>::const_iterator it = vectorOfBodyPairs.cbegin(); it!=vectorOfBodyPairs.cend(); ++it)
    {
        //! --------------
        //! enable "Stop"
        //! --------------
        myProgressIndicator->enableStop();

        //! -----------------
        //! check the status
        //! -----------------
        if(Global::status().code ==0)
        {
            cout<<"contactFinder::perform()->____process interrupted by the used____"<<endl;
            Global::status().code = 1;

            //! -----------------------------
            //! reset the progress indicator
            //! -----------------------------
            aProgressEvent = new QProgressEvent(QProgressEvent_Reset,0,100,0,"Process interrupted by the user",
                                                QProgressEvent_Reset,0,100,0,"");
            QApplication::postEvent(myProgressIndicator,aProgressEvent);
            QApplication::processEvents();

            stopped = true;
            return stopped;
        }

        const std::pair<GeometryTag,GeometryTag> &curPair = *it;
        int masterBodyIndex = curPair.first.parentShapeNr;
        int slaveBodyIndex = curPair.second.parentShapeNr;

        cout<<"contactFinder::perform()->____working on pair (masterIndex, slaveIndex) = ("<<masterBodyIndex<<", "<<slaveBodyIndex<<")____"<<endl;

        //! --------------------------------------------
        //! preliminary check of bounding box distances
        //! --------------------------------------------
        Bnd_Box boundingBoxMaster, boundingBoxSlave;
        const TopoDS_Shape &masterBody = myMDB->bodyMap.value(masterBodyIndex);
        const TopoDS_Shape &slaveBody = myMDB->bodyMap.value(slaveBodyIndex);
        BRepBndLib::Add(masterBody, boundingBoxMaster);
        BRepBndLib::Add(slaveBody, boundingBoxSlave);
        double bbdistance = boundingBoxSlave.Distance(boundingBoxMaster);

        if(bbdistance>tolerance)
        {
            //! ---------------------
            //! post an update event
            //! ---------------------
            aProgressEvent = new QProgressEvent(QProgressEvent_Update,0,9999,done,"Shape pair BBXs too far",
                                                QProgressEvent_None,0,9999,0,task);
            QApplication::postEvent(myProgressIndicator,aProgressEvent);
            QApplication::processEvents();

            done++;
            continue;
        }

        //! ------------------------------------
        //! source and target mesh for the algo
        //! ------------------------------------
        std::pair<int,occHandle(MeshVS_DataSource)> p1 = *mapOfSurfaceMeshDataSources.find(masterBodyIndex);
        if(p1.second.IsNull())
        {
            cerr<<"contactFinder::perform()->____error: cannot retrieve the tessellation of master body"<<endl;
            exit(1);
        }
        std::pair<int,occHandle(MeshVS_DataSource)> p2 = *mapOfSurfaceMeshDataSources.find(slaveBodyIndex);
        if(p2.second.IsNull())
        {
            cerr<<"contactFinder::perform()->____error: cannot retrieve the tessellation of slave body"<<endl;
            exit(2);
        }

        //! ---------------------------------------
        //! "master" and "slave" mesh data sources
        //! ---------------------------------------
        occHandle(Ng_MeshVS_DataSource2D) aMasterMesh2D = occHandle(Ng_MeshVS_DataSource2D)::DownCast(p1.second);
        occHandle(Ng_MeshVS_DataSource2D) aSlaveMesh2D = occHandle(Ng_MeshVS_DataSource2D)::DownCast(p2.second);

        //! -----------------------------------------------------------
        //! source (centers) to target (elements)
        //! here the master is used as source and, the slave as target
        //! -----------------------------------------------------------
        int start = getTimePoint();

        std::map<size_t,int> trianglesMapForBodySource = trianglesMap.at(masterBodyIndex);
        std::map<size_t,int> trianglesMapForBodyTarget = trianglesMap.at(slaveBodyIndex);
        std::map<int,occHandle(MeshVS_DataSource)> mapOfFaceMeshDSForBodySource = mapOfMapOfFaceMeshDataSources.at(masterBodyIndex);
        std::map<int,occHandle(MeshVS_DataSource)> mapOfFaceMeshDSForBodyTarget = mapOfMapOfFaceMeshDataSources.at(slaveBodyIndex);

        std::vector<std::pair<int,int>>  vecSourceTargetFaceNr;
        stopped = this->sourceNodesToTargetElement1(aMasterMesh2D,
                                                    aSlaveMesh2D,
                                                    mapOfFaceMeshDSForBodyTarget,
                                                    mapOfFaceMeshDSForBodySource,
                                                    tolerance,
                                                    trianglesMapForBodyTarget,
                                                    trianglesMapForBodySource,
                                                    vecSourceTargetFaceNr);

        int end = getTimePoint();
        for(std::vector<std::pair<int,int>>::iterator it = vecSourceTargetFaceNr.begin(); it!= vecSourceTargetFaceNr.end(); it++)
        {
            std::vector<GeometryTag> tagsMaster,tagsSlave;

            std::pair<int,int> apair = *it;

            GeometryTag aMasterTag;
            aMasterTag.isParent = false;
            aMasterTag.parentShapeNr = masterBodyIndex;
            aMasterTag.subTopNr = apair.first;
            aMasterTag.subShapeType = TopAbs_FACE;
            tagsMaster.push_back(aMasterTag);

            GeometryTag aSlaveTag;
            aSlaveTag.isParent = false;
            aSlaveTag.parentShapeNr = slaveBodyIndex;
            aSlaveTag.subTopNr = apair.second;
            aSlaveTag.subShapeType = TopAbs_FACE;

            tagsSlave.push_back(aSlaveTag);

            std::pair<std::vector<GeometryTag>,std::vector<GeometryTag>> apair_;
            apair_.first = tagsMaster;
            apair_.second = tagsSlave;

            allContactsPairs.push_back(apair_);
        }

        TotalTimeForSearch += double((end-start)/1000.0);

        //cout<<"@------------------------------------@"<<endl;
        //cout<<"@- time for search [s] "<<double((end-start)/1000.0)<<endl;
        //cout<<"@------------------------------------@"<<endl;

        if(stopped)
        {
            cout<<"contactFinder::perform()->____process interrupted by the used during mesh proximity analysis____"<<endl;
            Global::status().code = 1;

            //! -----------------------------
            //! reset the progress indicator
            //! -----------------------------
            aProgressEvent = new QProgressEvent(QProgressEvent_Reset,0,100,0,"Process interrupted by the user",
                                                QProgressEvent_Reset,0,100,0,"");
            QApplication::postEvent(myProgressIndicator,aProgressEvent);
            QApplication::processEvents();
            return stopped;
        }

        bool swap = true;
        if(swap)
        {
            //! -----------------------------------------------------------
            //! target (centers) to source (elements)
            //! here the slave is used as source, and the master as target
            //! -----------------------------------------------------------
            start = getTimePoint();

            std::map<size_t,int> trianglesMapForBodySource = trianglesMap.at(slaveBodyIndex);
            std::map<size_t,int> trianglesMapForBodyTarget = trianglesMap.at(masterBodyIndex);
            std::map<int,occHandle(MeshVS_DataSource)> mapOfFaceMeshDSForBodySource = mapOfMapOfFaceMeshDataSources.at(slaveBodyIndex);
            std::map<int,occHandle(MeshVS_DataSource)> mapOfFaceMeshDSForBodyTarget = mapOfMapOfFaceMeshDataSources.at(masterBodyIndex);

            stopped = this->sourceNodesToTargetElement1(aSlaveMesh2D,
                                                        aMasterMesh2D,
                                                        mapOfFaceMeshDSForBodyTarget,
                                                        mapOfFaceMeshDSForBodySource,
                                                        tolerance,
                                                        trianglesMapForBodyTarget,
                                                        trianglesMapForBodySource,
                                                        vecSourceTargetFaceNr);

            end = getTimePoint();
            TotalTimeForSearch += double((end-start)/1000.0);

            //cout<<"@------------------------------------@"<<endl;
            //cout<<"@- time for search [s] "<<double((end-start)/1000.0)<<endl;
            //cout<<"@------------------------------------@"<<endl;

            if(stopped)
            {
                cout<<"contactFinder::perform()->____process interrupted by the used during mesh proximity analysis____"<<endl;
                Global::status().code = 1;

                //! -----------------------------
                //! reset the progress indicator
                //! -----------------------------
                aProgressEvent = new QProgressEvent(QProgressEvent_Reset,0,100,0,"Process interrupted by the user",
                                                    QProgressEvent_Reset,0,100,0,"");
                QApplication::postEvent(myProgressIndicator,aProgressEvent);
                QApplication::processEvents();
                return stopped;
            }
        }

        //! ---------------------
        //! post an update event
        //! ---------------------
        aProgressEvent = new QProgressEvent(QProgressEvent_Update,0,9999,done,"Contact pairs found",
                                            QProgressEvent_None,0,9999,0,task);
        QApplication::postEvent(myProgressIndicator,aProgressEvent);
        QApplication::processEvents();

        done++;
    }

    cout<<"@------------------------------------@"<<endl;
    cout<<"@- total time for search [s] "<<TotalTimeForSearch<<endl;
    cout<<"@------------------------------------@"<<endl;

    //! ---------------------------------------------------
    //! post an event: let the user see what has been done
    //! ---------------------------------------------------
    aProgressEvent = new QProgressEvent(QProgressEvent_Reset,0,100,0,"Contact pairs found",
                                        QProgressEvent_Reset,0,100,0,"");
    QApplication::postEvent(myProgressIndicator,aProgressEvent);
    QApplication::processEvents();
    QThread::msleep(500);

    std::vector<std::pair<std::vector<GeometryTag>,std::vector<GeometryTag>>> reordered;
    this->groupBy(allContactsPairs,reordered,grouping);
    allContactsPairs.clear();       //! unessential: left here for showing the intention
    allContactsPairs = reordered;

    //! ----------------------------------------------
    //! return true: the process has not been stopped
    //! ----------------------------------------------
    return stopped;
}

//! ------------------------------------------------------------------------------
//! function: isPointProjectionInsideTriangle
//! details:  input a triangle T, a point P
//!           output the coordinates of the point PP projected onto the triangle,
//!           the (abs value) of the distance d(P,PP), true is PP lays inside T
//! ------------------------------------------------------------------------------
bool contactFinder::isPointProjectionInsideTriangle(const std::vector<polygon::Point> &aTriangle,
                                                     const polygon::Point &P,
                                                     polygon::Point &PP,
                                                     double &distance)
{
    //! ---------------------------
    //! the points of the triangle
    //! ---------------------------
    const polygon::Point &V1 = aTriangle[0];
    const polygon::Point &V2 = aTriangle[1];
    const polygon::Point &V3 = aTriangle[2];

    //! ------------------------------------------------------------------
    //! check degenerancy: the STL triangulation could not be well formed
    //! ------------------------------------------------------------------
    if(polygon::triangleArea(V1,V2,V3)<1e-6) return false;

    //! ---
    //! eX
    //! ---
    double A = V2.x - V1.x;
    double B = V2.y - V1.y;
    double C = V2.z - V1.z;
    double l21 = sqrt(A*A + B*B + C*C);

    if(l21==0) exit(21);

    double eXx = A / l21; double eXy = B / l21; double eXz = C / l21;

    //! ----------
    //! eZ
    //!	i	j	k
    //!	A	B	C
    //!	D	E	F
    //! ----------
    double D = V3.x - V1.x;
    double E = V3.y - V1.y;
    double F = V3.z - V1.z;
    double eZx = B*F-C*E;
    double eZy = C*D-A*F;
    double eZz = A*E-B*D;
    double l31 = sqrt(eZx*eZx + eZy*eZy + eZz*eZz);

    if(l31==0) exit(22);

    eZx /= l31; eZy /= l31; eZz /= l31;

    //! ------------
    //! eY
    //!	i	j	k
    //!	eZx	eZy	eZz
    //!	eXx	eXy	eXz
    //! ------------
    double eYx = eZy*eXz-eZz*eXy;
    double eYy = eZz*eXx-eZx*eXz;
    double eYz = eZx*eXy-eZy*eXx;
    double ly = sqrt(eYx*eYx + eYy*eYy + eYz*eYz);

    if(ly==0) exit(23);

    eYx /= ly; eYy /= ly; eYz /= ly;

    //! -------------------
    //! the local triangle
    //! -------------------
    double X1 = 0.0; double Y1 = 0.0; double Z1 = 0.0;
    double X2 = l21; double Y2 = 0.0; //double Z2 = 0.0;
    double Dx = V3.x - V1.x;
    double Dy = V3.y - V1.y;
    double Dz = V3.z - V1.z;
    double X3 = eXx*Dx + eXy*Dy + eXz*Dz;
    double Y3 = eYx*Dx + eYy*Dy + eYz*Dz;
    //double Z3 = eZx*Dx + eZy*Dy + eZz*Dz;		// should be exactly zero by definition

    //! -------------------------------------------
    //! P point projection onto the triangle plane
    //! -------------------------------------------
    Dx = P.x - V1.x;
    Dy = P.y - V1.y;
    Dz = P.z - V1.z;
    double XP = eXx*Dx + eXy*Dy + eXz*Dz;
    double YP = eYx*Dx + eYy*Dy + eYz*Dz;
    double ZP = eZx*Dx + eZy*Dy + eZz*Dz;		// "elevation" of P over the triangle plane => distance from the triangle plane
    //cout << "\\____PP projection - local CS - (" << XP << ", " << YP << ", " << ZP << ")____" << endl;

    //! -------------------------------------------------------------------
    //! return value - projected point into the original coordinate system
    //! -------------------------------------------------------------------
    PP.x = XP*eXx + YP*eYx + ZP*eZx + V1.x;
    PP.y = XP*eXy + YP*eYy + ZP*eZy + V1.y;
    PP.z = XP*eXz + YP*eYz + ZP*eZz + V1.z;

    polygon::Point PPp(XP, YP, 0.0);	//! projected point in the local coordinate system
    polygon::Point P1(X1, Y1, Z1);		//! (0,0,0) local origin
    polygon::Point P2(X2, Y2, 0.0);
    polygon::Point P3(X3, Y3, 0.0);

    double Area = triangleArea(P1,P2,P3);
    if(Area<1e-6) exit(20);

    double lambda1 = triangleArea(PPp, P2, P3)/Area;	//!
    double lambda2 = triangleArea(PPp, P3, P1)/Area;	//!
    const double lamtol = 1e-6;
    distance = fabs(ZP);

    //if (lambda1 > -lamtol && lambda2 > -lamtol && lambda1 < (1+lamtol) && lambda2 < (1+lamtol) && (lambda1 + lambda2) < (1 + lamtol))
    if (lambda1 > -lamtol && lambda2 > -lamtol && (1 - lambda1 - lambda2) > lamtol)
    {
        return true;
    }
    else return false;
}

//! -----------------------
//! function: getTimePoint
//! details:
//! -----------------------
int contactFinder::getTimePoint()
{
    std::chrono::time_point<std::chrono::system_clock> now = std::chrono::system_clock::now();
    auto duration = now.time_since_epoch();
    auto millis = std::chrono::duration_cast<std::chrono::milliseconds>(duration).count();
    return millis;
}

//! --------------------------------------------------------------------------------------------
//! function: prepareSTLTessellation_useDisk
//! details:  given a vector of body pairs (actually of GeometryTag) returns:
//!           1. a map whose index is a body index, and whose value is the STL mesh data source
//!           2. a map whose index is a body index, and whose value is a map, whose index is
//!              an STL triangle hash code and whose value the corresponding face number
//! --------------------------------------------------------------------------------------------
int contactFinder::prepareSTLTessellation_useDisk(const std::vector<std::pair<GeometryTag,GeometryTag>> &vectorOfBodyPairs,
                                                   bool rebuildTessellation,
                                                   std::map<int,occHandle(MeshVS_DataSource)> &mapOfSurfaceMeshDataSources,
                                                   std::map<int,std::map<int,occHandle(MeshVS_DataSource)>> &mapOfMapOfFaceMeshDataSources)
{
    cout<<"contactFinder::prepareSTLTessellation_useDisk()->____function called____"<<endl;

    std::set<int> parentShapeNumbers;
    for(std::vector<std::pair<GeometryTag,GeometryTag>>::const_iterator it = vectorOfBodyPairs.cbegin();
        it != vectorOfBodyPairs.cend(); it++)
    {
        std::pair<GeometryTag,GeometryTag> aPair = *it;
        parentShapeNumbers.insert(aPair.first.parentShapeNr);
        parentShapeNumbers.insert(aPair.second.parentShapeNr);
    }

    //! -------------------------
    //! rebuild the tessellation
    //! -------------------------
    if(rebuildTessellation)
    {
        for(std::set<int>::iterator it = parentShapeNumbers.begin(); it!= parentShapeNumbers.end(); it++)
        {
            int bodyIndex = *it;
            const TopoDS_Shape &curShape = myMDB->bodyMap.value(bodyIndex,TopoDS_Shape());
            if(curShape.IsNull()) continue;
            double ad = 0.4; double ld = 0.1;
            BRepTools::Clean(curShape);
            BRepMesh_IncrementalMesh aMesher(curShape,ld,true,ad,true,false);
            Q_UNUSED(aMesher);
        }
    }

    //! -------------------------------------------------
    //! convert the tessellations into mesh data sources
    //! -------------------------------------------------
    for(std::set<int>::iterator it = parentShapeNumbers.begin(); it!= parentShapeNumbers.end(); it++)
    {
        int bodyIndex = *it;
        const TopoDS_Shape &curShape = myMDB->bodyMap.value(bodyIndex,TopoDS_Shape());
        if(curShape.IsNull()) continue;

        QString fileName("D:/test.stl");

        //! ----------------------------
        //! write the extended stl file
        //! ----------------------------
        STLAPIWriter aWriter;
        aWriter.WriteExtended(curShape,fileName.toStdString().c_str());

        //! ---------------------------
        //! read the extended stl file
        //! ---------------------------
        int NbFaces = myMDB->MapOfBodyTopologyMap.value(bodyIndex).faceMap.Extent();
        occHandle(Ng_MeshVS_DataSource2D) anSTLMesh;

        //NCollection_Array1<occHandle(Ng_MeshVS_DataSourceFace)> faceSubMeshes = NCollection_Array1<occHandle(Ng_MeshVS_DataSourceFace)>(0,NbFaces);
        //MeshTools::toSTLMesh(curShape,fileName,anSTLMesh,faceSubMeshes,Q_NULLPTR,0);

        std::vector<occHandle(Ng_MeshVS_DataSourceFace)> faceSubMeshes;
        for(int i=0; i<=NbFaces; i++) faceSubMeshes.push_back(occHandle(Ng_MeshVS_DataSourceFace)());
        MeshTools::toSTLMesh1(curShape,fileName,anSTLMesh,faceSubMeshes);

        //! ----------------------------------------------
        //! fill the map of the surface mesh data sources
        //! ----------------------------------------------
        std::pair<int,occHandle(MeshVS_DataSource)> aPair;
        aPair.first = bodyIndex;
        aPair.second = anSTLMesh;
        mapOfSurfaceMeshDataSources.insert(aPair);

        //! -----------------------------------------------
        //! fill the map of the face mesh stl data sources
        //! -----------------------------------------------
        std::map<int,occHandle(MeshVS_DataSource)> mapOfFaceMeshDataSources;
        for(int faceNr = 0; faceNr<=NbFaces; faceNr++)
        //for (int faceNr = faceSubMeshes.Lower(); faceNr<=faceSubMeshes.Upper(); faceNr++)
        {
            //occHandle(Ng_MeshVS_DataSourceFace) aFaceMeshDS = faceSubMeshes.Value(faceNr);
            occHandle(Ng_MeshVS_DataSourceFace) aFaceMeshDS = faceSubMeshes[faceNr];
            if(aFaceMeshDS.IsNull())
            {
                cerr<<"____faceNr: "<<faceNr<<" NULL face mesh data source____"<<endl;
                continue;
            }
            std::pair<int,occHandle(MeshVS_DataSource)> aPair_;
            aPair_.first = faceNr;
            aPair_.second = aFaceMeshDS;
            mapOfFaceMeshDataSources.insert(aPair_);
        }
        std::pair<int,std::map<int,occHandle(MeshVS_DataSource)>> aPair__;
        aPair__.first = bodyIndex;
        aPair__.second = mapOfFaceMeshDataSources;
        mapOfMapOfFaceMeshDataSources.insert(aPair__);
    }
    return int(mapOfSurfaceMeshDataSources.size());
}

//! ----------------------------------------------
//! function: TriangleToTriangleTest
//! details:  test triangle to triangle proximity
//! ----------------------------------------------
void contactFinder::TriangleToTriangleTest(const triangle &trig1, const triangle &trig2, bool &intersecting, double &distance)
{
    //! ----------------------------------------------------------------------------------------------------
    //! first check if triangles are intersecting. It uses:
    //! tri_tri_intersection_test_3d(real p1[3], real q1[3], real r1[3], real p2[3], real q2[3], real r2[3],
    //! int * coplanar, real source[3], real target[3] )
    //! -----------------------------------------------------------------------------------------------------

    //! ---------------------------------------------
    //! define the input coordinates: first triangle
    //! ---------------------------------------------
    TriangleToTriangleIntersection::real p1[3], q1[3], r1[3];
    const polygon::Point &V1 = trig1[0];
    const polygon::Point &V2 = trig1[1];
    const polygon::Point &V3 = trig1[2];
/*
    p1[0] = V1.x; p1[1] = V2.x; p1[2] = V3.x;
    q1[0] = V1.y; q1[1] = V2.y; q1[2] = V3.y;
    r1[0] = V1.z; r1[1] = V2.z; r1[2] = V3.z;
*/
    p1[0] = V1.x; p1[1] = V1.y; p1[2] = V1.z;
    q1[0] = V2.x; q1[1] = V2.y; q1[2] = V2.z;
    r1[0] = V3.x; r1[1] = V3.y; r1[2] = V3.z;

    //! ----------------------------------------------
    //! define the input coordinates: second triangle
    //! ----------------------------------------------
    TriangleToTriangleIntersection::real p2[3], q2[3], r2[3];
    const polygon::Point &Z1 = trig2[0];
    const polygon::Point &Z2 = trig2[1];
    const polygon::Point &Z3 = trig2[2];
/*
    p2[0] = Z1.x; p2[1] = Z2.x; p2[2] = Z3.x;
    q2[0] = Z1.y; q2[1] = Z2.y; q2[2] = Z3.y;
    r2[0] = Z1.z; r2[1] = Z2.z; r2[2] = Z3.z;
*/
    p2[0] = Z1.x; p2[1] = Z1.y; p2[2] = Z1.z;
    q2[0] = Z2.x; q2[1] = Z2.y; q2[2] = Z2.z;
    r2[0] = Z3.x; r2[1] = Z3.y; r2[2] = Z3.z;

    //! ------------------------------------------------------------
    //! define the results: if res == 1 the two triangles intersect
    //! ------------------------------------------------------------
    int coplanar;
    TriangleToTriangleIntersection::real source[3], target[3];

    int res = TriangleToTriangleIntersection::tri_tri_intersection_test_3d(p1, q1, r1, p2, q2, r2, &coplanar, source, target);
    if(res==1)
    {
        distance = 0.0;
        intersecting = true;
        return;
    }

    //Q_UNUSED(coplanar)
    Q_UNUSED(source)
    Q_UNUSED(target)

    intersecting = false;

    //cout<<"@ triangles do not intersect "<<endl;
    //cout<<"@ triangles are"<<(coplanar==1? " ":" not")<<" coplanar"<<endl;

    //! -----------------------------------------------------------------------------------------------
    //! if the triangles are not intersecting perform a distance query
    //!
    //! PQP_REAL TriDist(PQP_REAL p[3], PQP_REAL q[3], const PQP_REAL s[3][3], const PQP_REAL t[3][3]);
    //!
    //! computes the closest points on two triangles, and returns the distance between them.
    //! s and t are the triangles, stored tri[point][dimension].
    //! If the triangles are disjoint, p and q give the closest points of s and t respectively.
    //! However, if the triangles overlap, p and q are basically a random pair of points from
    //! the triangles, not coincident points on the intersection of the triangles, as might be expected
    //!
    //! struct Tri
    //! {
    //!   PQP_REAL p1[3];
    //!   PQP_REAL p2[3];
    //!   PQP_REAL p3[3];
    //!   int id;
    //! };
    //! ------------------------------------------------------------------------------------------------

    //! ---------------------------
    //! define the input triangles
    //! ---------------------------
    PQP_REAL s[3][3], t[3][3];
    s[0][0] = V1.x; s[0][1] = V1.y; s[0][2] = V1.z;
    s[1][0] = V2.x; s[1][1] = V2.y; s[1][2] = V2.z;
    s[2][0] = V3.x; s[2][1] = V3.y; s[2][2] = V3.z;

    t[0][0] = Z1.x; t[0][1] = Z1.y; t[0][2] = Z1.z;
    t[1][0] = Z2.x; t[1][1] = Z2.y; t[1][2] = Z2.z;
    t[2][0] = Z3.x; t[2][1] = Z3.y; t[2][2] = Z3.z;

    //! -------------------
    //! define the results
    //! -------------------
    /*
    PQP_REAL R[3][3];
    R[0][0] = 1.0; R[0][1] = 0.0; R[0][2] = 0.0;
    R[1][0] = 0.0; R[1][1] = 1.0; R[1][2] = 0.0;
    R[2][0] = 0.0; R[2][1] = 0.0; R[2][2] = 1.0;
    PQP_REAL T[3]{0,0,0};
    Tri t1,t2;
    t1.p1[0] = V1.x; t1.p1[1] = V1.y; t1.p1[2] = V1.z;
    t1.p2[0] = V2.x; t1.p2[1] = V2.y; t1.p2[2] = V2.z;
    t1.p3[0] = V3.x; t1.p3[1] = V3.y; t1.p3[2] = V3.z;

    t2.p1[0] = Z1.x; t2.p1[1] = Z1.y; t2.p1[2] = Z1.z;
    t2.p2[0] = Z2.x; t2.p2[1] = Z2.y; t2.p2[2] = Z2.z;
    t2.p3[0] = Z3.x; t2.p3[1] = Z3.y; t2.p3[2] = Z3.z;


    cout<<"____("<<V1.x<<", "<<V1.y<<", "<<V1.z<<")____"<<endl;
    cout<<"____("<<V2.x<<", "<<V2.y<<", "<<V2.z<<")____"<<endl;
    cout<<"____("<<V3.x<<", "<<V3.y<<", "<<V3.z<<")____"<<endl;

    cout<<"____("<<Z1.x<<", "<<Z1.y<<", "<<Z1.z<<")____"<<endl;
    cout<<"____("<<Z2.x<<", "<<Z2.y<<", "<<Z2.z<<")____"<<endl;
    cout<<"____("<<Z3.x<<", "<<Z3.y<<", "<<Z3.z<<")____"<<endl;

    PQP_REAL a[3], b[3];
    distance = TriDistance(R,T,&t1,&t2,a,b);
    */

    PQP_REAL a[3], b[3];
    distance = TriDist(a,b,s,t);
}

//! --------------------------------------
//! function: sourceNodesToTargetElement1
//! details:
//! --------------------------------------
bool contactFinder::sourceNodesToTargetElement1(const occHandle(Ng_MeshVS_DataSource2D) &aSourceMesh,
                                                 const occHandle(Ng_MeshVS_DataSource2D) &aTargetMesh,
                                                 const std::map<int,occHandle(MeshVS_DataSource)> &mapOfFaceMeshDSForBodyTarget,
                                                 const std::map<int,occHandle(MeshVS_DataSource)> &mapOfFaceMeshDSForBodySource,
                                                 double tolerance,
                                                 const std::map<size_t,int> &trianglesMapForBodyTarget,
                                                 const std::map<size_t,int> &trianglesMapForBodySource,
                                                 std::vector<std::pair<int,int>> &vecSourceTargetFaceNr)
{
    cout<<"contactFinder::sourceNodesToTargetElement1()->____function called____"<<endl;

    //! ---------------
    //! execution flag
    //! ---------------
    bool stopped = false;

    //! ---------------------------------------------
    //! iterate over the elements of the source mesh
    //! ---------------------------------------------
    int NbElements = aSourceMesh->GetAllElements().Extent();
    int NbScannedElements = 0;

    //! -------------------
    //! send an init event
    //! -------------------
    QProgressEvent *aProgressEvent;
    QString msg;
    if(myProgressIndicator!=Q_NULLPTR)
    {
        msg = QString("Analyzing %1 surface mesh elements").arg(aSourceMesh->GetAllElements().Extent());
        QProgressEvent *aProgressEvent = new QProgressEvent(QProgressEvent_None,0,9999,0,msg,
                                                            QProgressEvent_Init,0,NbElements,0,"Creating automatic connections");
        QApplication::postEvent(myProgressIndicator,aProgressEvent);
        QApplication::processEvents();
    }

    //! ----------------------------------------
    //! iterate over the "source" mesh elements
    //! ----------------------------------------
    TColStd_PackedMapOfInteger mapOfTargetElements = aTargetMesh->GetAllElements();
    TColStd_MapIteratorOfPackedMapOfInteger elIt;

    bool pairFound = false;
    int sourceFaceNr, targetFaceNr;

    TColStd_MapIteratorOfPackedMapOfInteger it;
    TColStd_PackedMapOfInteger mapOfSourceElements = aSourceMesh->GetAllElements();
    it.Initialize(mapOfSourceElements);

    //for(TColStd_MapIteratorOfPackedMapOfInteger it(aSourceMesh->GetAllElements()); it.More(); it.Next())
    for(; it.More(); it.Next())
    {
        NbScannedElements++;

        //! ----------------------
        //! check the global flag
        //! ----------------------
        if(Global::status().code==0)
        {
            //! --------------------------
            //! send an interrupt message
            //! --------------------------
            aProgressEvent = new QProgressEvent(QProgressEvent_None,0,100,0,"Operation interrupted by the user",
                                                QProgressEvent_Update,0,NbElements,NbScannedElements,"");
            QApplication::postEvent(myProgressIndicator,aProgressEvent);
            QApplication::processEvents();
            stopped = true;
            return stopped;
        }

        //! -----------------------------------------------------------------
        //! retrieve the coordinates of the center "C" of the source element
        //! -----------------------------------------------------------------
        int globalElementID_source = it.Key();

        MeshVS_EntityType aType_s;
        int NbNodes_s;
        double buf_s[9];
        TColStd_Array1OfReal coords_s(*buf_s,1,9);
        aSourceMesh->GetGeom(globalElementID_source,true,coords_s,NbNodes_s,aType_s);

        //! -------------------------
        //! create a source triangle
        //! -------------------------
        triangle aSourceTriangle;
        for(int k=0; k<NbNodes_s; k++)
        {
            int s = 3*k;
            polygon::Point aPoint(coords_s(s+1),coords_s(s+2),coords_s(s+3));
            aSourceTriangle.push_back(aPoint);
        }
        size_t t_source_hash = this->triangleHash(aSourceTriangle);

        //! ---------------------------------------
        //! the face number of the source triangle
        //! ---------------------------------------
        sourceFaceNr = trianglesMapForBodySource.at(t_source_hash);

        double area_s = polygon::triangleArea(aSourceTriangle.at(0),aSourceTriangle.at(1),aSourceTriangle.at(2));
        if(area_s <1e-6) continue;

        //! -------------------------------------------------
        //! iterate over the elements of the target mesh:
        //! check the proximity of the center of the current
        //! source element to the current target element:
        //! if the proximity is OK exit cycle
        //! -------------------------------------------------
        int j = 1;
        for(elIt.Initialize(mapOfTargetElements); elIt.More(); elIt.Next(), j++)
        {
            int globalElementID_target = elIt.Key();

            int NbNodes_t;
            double buf_t[9];
            MeshVS_EntityType aType_t;
            TColStd_Array1OfReal coords_t(*buf_t,1,9);
            bool isDone = aTargetMesh->GetGeom(globalElementID_target,true,coords_t,NbNodes_t,aType_t);
            if(!isDone) continue;

            //! -------------------------
            //! create a target triangle
            //! -------------------------
            triangle aTargetTriangle;
            for(int k=0; k<NbNodes_t; k++)
            {
                int s = 3*k;
                polygon::Point aPoint(coords_t(s+1),coords_t(s+2),coords_t(s+3));
                aTargetTriangle.push_back(aPoint);
            }

            double area_t = polygon::triangleArea(aTargetTriangle.at(0),aTargetTriangle.at(1),aTargetTriangle.at(2));
            if(area_t <1e-6) continue;

            double distance = 1e20;     //! dummy initialization
            bool intersecting = false;  //! dummy initialization

            TriangleToTriangleTest(aSourceTriangle,aTargetTriangle,intersecting,distance);
            if(intersecting == false && distance>tolerance) continue;

            //! -------------------------------------------------
            //! angular criterion: discard the "source" element
            //! if it is "almost" normal to the "target" element
            //! -------------------------------------------------
            double sx,sy,sz,tx,ty,tz;
            bool normalOK = aSourceMesh->GetNormal(globalElementID_source,8,sx,sy,sz);
            if(!normalOK) continue;
            normalOK = aTargetMesh->GetNormal(globalElementID_target,8,tx,ty,tz);
            if(!normalOK) continue;
            if(AngularCriterion(sx,sy,sz,tx,ty,tz)==false) continue;

            //! ------------------------------------------------------------
            //! the (source element, target element) pair has been accepted
            //! ------------------------------------------------------------

            //! ----------------------------------------------------------
            //! - retrieve the face number of the detected target element
            //! - retrieve all the elements of the face
            //! - remove from the target surface mesh
            //! ----------------------------------------------------------
            size_t t_hash = triangleHash(aTargetTriangle);
            if(trianglesMapForBodyTarget.count(t_hash)==1) cout<<"____found in target____"<<endl;

            //! -----------------------------------
            //! the face number of the target face
            //! -----------------------------------
            targetFaceNr = trianglesMapForBodyTarget.at(t_hash);
            pairFound = true;

            //! ------------------------------------------------------------------
            //! speed-up the search: remove the target elements of the whole face
            //! ------------------------------------------------------------------
            if(pairFound)
            {
                cout<<"____("<<sourceFaceNr<<", "<<targetFaceNr<<")____"<<endl;
                std::map<int,occHandle(MeshVS_DataSource)>::const_iterator iit = mapOfFaceMeshDSForBodyTarget.find(targetFaceNr);

                if(iit==mapOfFaceMeshDSForBodyTarget.cend())
                {
                    cerr<<"___the target face mesh is NULL____"<<endl;
                    continue;
                }
                const occHandle(MeshVS_DataSource) &faceMeshDS = (*iit).second;

                //cout<<"_____________________________________"<<endl;
                //cout<<"____initial target map key dimension: "<<mapOfTargetElements.Extent()<<"____"<<endl;
                for(TColStd_MapIteratorOfPackedMapOfInteger anit(faceMeshDS->GetAllElements()); anit.More(); anit.Next())
                {
                    int aKey = anit.Key();
                    bool keyRemoved = mapOfTargetElements.Remove(aKey);
                    if(keyRemoved==false) cout<<"____cannot remove____"<<endl;
                    //cout<<"____target face Nr: "<<targetFaceNr<<"____removing: "<<aKey<<"____"<<endl;
                }
                //cout<<"____final target map key dimension: "<<mapOfTargetElements.Extent()<<"____"<<endl;
                //cout<<"_____________________________________"<<endl;

                std::pair<int,int> sourceTargetPair;
                sourceTargetPair.first = sourceFaceNr;
                sourceTargetPair.second = targetFaceNr;
                vecSourceTargetFaceNr.push_back(sourceTargetPair);

                pairFound = false;

                //! ----------------------------------------------------------------
                //! experimental - this is a source of error for complex geometries
                //! has been removed with this - it remains to understand why
                //! ----------------------------------------------------------------
                elIt.Initialize(mapOfTargetElements);
            }

            //break;
        }
        /*
        //! ------------------------------------------------------------------
        //! speed-up the search: remove the target elements of the whole face
        //! ------------------------------------------------------------------
        if(pairFound)
        {
            cout<<"____("<<sourceFaceNr<<", "<<targetFaceNr<<")____"<<endl;
            std::map<int,occHandle(MeshVS_DataSource)>::const_iterator iit = mapOfFaceMeshDSForBodyTarget.find(targetFaceNr);

            if(iit==mapOfFaceMeshDSForBodyTarget.cend())
            {
                cerr<<"___the target face mesh is NULL____"<<endl;
                continue;
            }
            const occHandle(MeshVS_DataSource) &faceMeshDS = (*iit).second;

            //cout<<"_____________________________________"<<endl;
            //cout<<"____initial target map key dimension: "<<mapOfTargetElements.Extent()<<"____"<<endl;
            for(TColStd_MapIteratorOfPackedMapOfInteger anit(faceMeshDS->GetAllElements()); anit.More(); anit.Next())
            {
                int aKey = anit.Key();
                bool keyRemoved = mapOfTargetElements.Remove(aKey);
                if(keyRemoved==false) cout<<"____cannot remove____"<<endl;
                //cout<<"____target face Nr: "<<targetFaceNr<<"____removing: "<<aKey<<"____"<<endl;
            }
            //cout<<"____final target map key dimension: "<<mapOfTargetElements.Extent()<<"____"<<endl;
            //cout<<"_____________________________________"<<endl;

            std::pair<int,int> sourceTargetPair;
            sourceTargetPair.first = sourceFaceNr;
            sourceTargetPair.second = targetFaceNr;
            vecSourceTargetFaceNr.push_back(sourceTargetPair);

            pairFound = false;

            //! ----------------------------------------------------------------
            //! experimental - this is a source of error for complex geometries
            //! has been removed with this - it remains to understand why
            //! ----------------------------------------------------------------
            elIt.Initialize(mapOfTargetElements);
        }
        */

        //! --------------------------
        //! remove the source element
        //! --------------------------
        //mapOfSourceElements.Remove(globalElementID_source);

        //! ----------------------------------
        //! update the secondary progress bar
        //! ----------------------------------
        if(myProgressIndicator!=Q_NULLPTR)
        {
            if(NbScannedElements%100==0)
            {
                aProgressEvent = new QProgressEvent(QProgressEvent_None,0,100,0,msg,
                                                    QProgressEvent_Update,0,NbElements-1,NbScannedElements,
                                                    "Creating automatic connections");
                QApplication::postEvent(myProgressIndicator,aProgressEvent);
                QApplication::processEvents();
            }
        }
    }

    //! --------------------------------------------------------
    //! the process has not been manually stopped: return false
    //! --------------------------------------------------------
    return stopped;
}

//! -------------------
//! function: groupBy
//! details:
//! -------------------
void contactFinder::groupBy(const std::vector<std::pair<std::vector<GeometryTag>,std::vector<GeometryTag>>> &allContactsPairs,
                             std::vector<std::pair<std::vector<GeometryTag>,std::vector<GeometryTag>>> &reordered,
                             int grouping)
{
    switch(grouping)
    {
    //! ---------------
    //! by bodies pair
    //! ---------------
    case 0:
    {
        cout<<"____FILTERING BY BODY PAIR____"<<endl;
        for(std::vector<std::pair<GeometryTag,GeometryTag>>::iterator it = myVectorOfBodyPairs.begin(); it!=myVectorOfBodyPairs.end(); ++it)
        {
            std::pair<std::vector<GeometryTag>,std::vector<GeometryTag>> reorderedPair;
            std::vector<GeometryTag> tagsMasterReordered;
            std::vector<GeometryTag> tagsSlaveReordered;

            std::pair<GeometryTag,GeometryTag> &curPair = *it;
            int masterBodyIndex = curPair.first.parentShapeNr;
            int slaveBodyIndex = curPair.second.parentShapeNr;

            int Nb = int(allContactsPairs.size());
            for(int i=0; i<Nb; i++)
            {
                std::pair<std::vector<GeometryTag>,std::vector<GeometryTag>> aContactPair = allContactsPairs[i];
                std::vector<GeometryTag> tagsMaster = aContactPair.first;
                std::vector<GeometryTag> tagsSlave = aContactPair.second;
                GeometryTag aMtag = tagsMaster.at(0);
                GeometryTag aStag = tagsSlave.at(0);
                if(aMtag.parentShapeNr == masterBodyIndex && aStag.parentShapeNr == slaveBodyIndex)
                {
                    tagsMasterReordered.push_back(aMtag);
                    tagsSlaveReordered.push_back(aStag);
                }
            }
            if(tagsMasterReordered.size()>0 && tagsSlaveReordered.size()>0)
            {
                reorderedPair.first = tagsMasterReordered;
                reorderedPair.second = tagsSlaveReordered;
                reordered.push_back(reorderedPair);
            }
        }
    }
        break;

    //! ----------------
    //! by master faces
    //! ----------------
    case 1:
    {
        std::map<GeometryTag,std::vector<GeometryTag>> supportMap;
        int Nb = int(allContactsPairs.size());
        for(int i=0; i<Nb; i++)
        {
            std::pair<std::vector<GeometryTag>,std::vector<GeometryTag>> aContact = allContactsPairs[i];
            GeometryTag masterTag = aContact.first[0];
            std::map<GeometryTag,std::vector<GeometryTag>>::iterator it = supportMap.find(masterTag);
            if(it==supportMap.end())
            {
                std::pair<GeometryTag,std::vector<GeometryTag>> element;
                element.first = masterTag;
                element.second = aContact.second;
                supportMap.insert(element);
            }
            else
            {
                //(*it).second.append(aContact.second);
                (*it).second.insert((*it).second.end(),aContact.second.begin(),aContact.second.end());
            }
        }

        for(std::map<GeometryTag,std::vector<GeometryTag>>::iterator it = supportMap.begin(); it!=supportMap.end(); it++)
        {
            std::pair<GeometryTag,std::vector<GeometryTag>> aPair = *it;
            std::vector<GeometryTag> masterTags;
            masterTags.push_back(aPair.first);
            std::vector<GeometryTag> slaveTags;
            slaveTags = aPair.second;

            std::pair<std::vector<GeometryTag>,std::vector<GeometryTag>> aPair_;
            aPair_.first = masterTags;
            aPair_.second = slaveTags;
            reordered.push_back(aPair_);
        }

        //! --------
        //! testing
        //! --------
        cout<<"@____AGGREGATION____@"<<endl;
        for(std::map<GeometryTag,std::vector<GeometryTag>>::iterator it = supportMap.begin(); it!=supportMap.end(); ++it)
        {
            GeometryTag mt = (*it).first;
            std::vector<GeometryTag> sts = (*it).second;
            cout<<"____master: ("<<mt.parentShapeNr<<", "<<mt.subTopNr<<")____"<<endl;
            for(int i=0; i<sts.size(); i++)
            {
                cout<<"_______slave: ("<<sts[i].parentShapeNr<<", "<<sts[i].subTopNr<<")_______"<<endl;
            }
        }
        cout<<"@____END AGGREGATION____@"<<endl;
    }
        break;

    //! ---------------
    //! by slave faces
    //! ---------------
    case 2:
    {
        //! -----------------------------------------
        //! key => slave face; value => master faces
        //! -----------------------------------------
        std::map<GeometryTag,std::vector<GeometryTag>> supportMap;
        int Nb = int(allContactsPairs.size());
        for(int i=0; i<Nb; i++)
        {
            std::pair<std::vector<GeometryTag>,std::vector<GeometryTag>> aContact = allContactsPairs[i];
            GeometryTag slaveTag = aContact.second[0];
            std::map<GeometryTag,std::vector<GeometryTag>>::iterator it = supportMap.find(slaveTag);
            if(it==supportMap.end())
            {
                std::pair<GeometryTag,std::vector<GeometryTag>> element;
                element.first = slaveTag;
                element.second = aContact.first;
                supportMap.insert(element);
            }
            else
            {
                //(*it).second.append(aContact.first);
                it->second.insert(it->second.end(),aContact.first.begin(),aContact.first.end());
            }
        }

        for(std::map<GeometryTag,std::vector<GeometryTag>>::iterator it = supportMap.begin(); it!=supportMap.end(); it++)
        {
            std::pair<GeometryTag,std::vector<GeometryTag>> aPair = *it;
            std::vector<GeometryTag> slaveTags;
            slaveTags.push_back(aPair.first);
            std::vector<GeometryTag> masterTags;
            masterTags = aPair.second;

            std::pair<std::vector<GeometryTag>,std::vector<GeometryTag>> aPair_;
            aPair_.first = masterTags;
            aPair_.second = slaveTags;
            reordered.push_back(aPair_);
        }
    }
        break;

    //! ---------------
    //! by master body
    //! ---------------
    case 3:
    {
        cout<<"____FILTERING BY MASTER BODY____"<<endl;
        //! ----------------------------------------
        //! key => index of master body
        //! value => index of the array of contacts
        //! ----------------------------------------
        std::map<int,std::vector<int>> supportMap;
        int Nb = int(allContactsPairs.size());
        for(int i=0; i<Nb; i++)
        {
            std::pair<std::vector<GeometryTag>,std::vector<GeometryTag>> aContact = allContactsPairs[i];

            //! ----------------------------------------------------------------
            //! by definition the "masterTags" vector contains only one element
            //! ----------------------------------------------------------------
            std::vector<GeometryTag> masterTags = aContact.first;

            int masterBodyNumber = masterTags[0].parentShapeNr;
            std::map<int,std::vector<int>>::iterator it = supportMap.find(masterBodyNumber);
            if(it == supportMap.end())
            {
                std::pair<int,std::vector<int>> element;
                element.first = masterBodyNumber;
                std::vector<int> indexList;
                indexList.push_back(i);
                element.second = indexList;
                supportMap.insert(element);
            }
            else { (*it).second.push_back(i); }
        }

        for(std::map<int,std::vector<int>>::iterator it = supportMap.begin(); it!=supportMap.end(); it++)
        {
            std::vector<int> listOfIndexOfContactPairs = (*it).second;
            std::vector<GeometryTag> masterTags, slaveTags;
            for(int k=0; k<listOfIndexOfContactPairs.size(); k++)
            {
                int indexOfContactPair = listOfIndexOfContactPairs[k];
                std::pair<std::vector<GeometryTag>,std::vector<GeometryTag>> apair = allContactsPairs[indexOfContactPair];
                masterTags.insert(masterTags.end(),apair.first.begin(),apair.first.end());
                slaveTags.insert(slaveTags.end(),apair.second.begin(),apair.second.end());
            }

            if(masterTags.size()==0) { cerr<<"____EMPTY MASTER TAGS____"<<endl; exit(1); }
            if(slaveTags.size()==0) { cerr<<"____EMPTY SLAVE TAGS____"<<endl; exit(2); }

            //! --------------------------
            //! define a new contact pair
            //! --------------------------
            std::pair<std::vector<GeometryTag>,std::vector<GeometryTag>> aContact;
            aContact.first = masterTags;
            aContact.second = slaveTags;
            reordered.push_back(aContact);
        }
    }
        break;

    //! ----------------
    //! By slave bodies
    //! ----------------
    case 4:
    {
        cout<<"____FILTERING BY SLAVE BODY____"<<endl;
        //! ----------------------------------------
        //! key => index of slave body
        //! value => index of the array of contacts
        //! ----------------------------------------
        std::map<int,std::vector<int>> supportMap;
        int Nb = int(allContactsPairs.size());
        for(int i=0; i<Nb; i++)
        {
            std::pair<std::vector<GeometryTag>,std::vector<GeometryTag>> aContact = allContactsPairs[i];

            //! ---------------------------------------------------------------
            //! by definition the "slaveTags" vector contains only one element
            //! ---------------------------------------------------------------
            std::vector<GeometryTag> slaveTags = aContact.second;

            int slaveBodyNumber = slaveTags[0].parentShapeNr;
            std::map<int,std::vector<int>>::iterator it = supportMap.find(slaveBodyNumber);
            if(it == supportMap.end())
            {
                std::pair<int,std::vector<int>> element;
                element.first = slaveBodyNumber;
                std::vector<int> indexList;
                indexList.push_back(i);
                element.second = indexList;
                supportMap.insert(element);
            }
            else { (*it).second.push_back(i); }
        }

        for(std::map<int,std::vector<int>>::iterator it = supportMap.begin(); it!=supportMap.end(); it++)
        {
            std::vector<int> listOfIndexOfContactPairs = (*it).second;

            std::vector<GeometryTag> masterTags, slaveTags;
            for(int k=0; k<listOfIndexOfContactPairs.size(); k++)
            {
                int indexOfContactPair = listOfIndexOfContactPairs[k];
                std::pair<std::vector<GeometryTag>,std::vector<GeometryTag>> apair = allContactsPairs[indexOfContactPair];
                masterTags.insert(masterTags.end(),apair.first.begin(),apair.first.end());
                slaveTags.insert(slaveTags.end(),apair.second.begin(),apair.second.end());
            }

            if(masterTags.size()==0) { cerr<<"____EMPTY MASTER TAGS____"<<endl; exit(1); }
            if(slaveTags.size()==0) { cerr<<"____EMPTY SLAVE TAGS____"<<endl; exit(2); }

            //! --------------------------
            //! define a new contact pair
            //! --------------------------
            std::pair<std::vector<GeometryTag>,std::vector<GeometryTag>> aContact;
            aContact.first = masterTags;
            aContact.second = slaveTags;
            reordered.push_back(aContact);
        }
    }
        break;

    //! -----------------------
    //! no grouping ("case 5")
    //! -----------------------
    default:
    {
        reordered = allContactsPairs;
    }
        break;
    }
}
