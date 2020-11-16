#include "smoothingtools.h"

//! ----
//! C++
//! ----
#include <vector>

//! -------------------------------------------------------------------------------------------------------
//! function: scalarFieldSmoother
//! details:  laplacian smoother of a scalar field as found in:
//!           Hybrid Grid Generation for Viscous Flow Simulations in Complex Geometries
//!           Hongfei Ye, Yang Liuc, Bo Chenc, Zhiwei Liua, Jianjing Zhenga, Yufei Pangc, Jianjun Chena
//!           Here the definition of omega relies on local dicrete gaussian curvature and on a sentivitity
//!           real parameter k
//! -------------------------------------------------------------------------------------------------------

//! ----------------------------------
//! definition of auxiliary functions
//! ----------------------------------
double fomega(double betaAve)
{
    double val = 0.5*(1+cos(betaAve/2));
    //const double PI = 3.1415926538;
    //double val = (PI-betaAve/2)*(1/(PI));
    return val;
}

double kij(double n, const std::vector<double> &P, const std::vector<double> &S, double Sdij)
{
    double dij = sqrt(pow(P[0]-S[0],2)+pow(P[1]-S[1],2)+pow(P[2]-S[2],2));
    double val = n*dij/Sdij;
    return val;
}

double wij(double betaAve, const std::vector<double> &P, const std::vector<double> &S, double Sdij, double n)
{
    const double PI = 3.1415926538;
    //const double eps = 5.0*PI/180.0;
    //double val = 1.0;
    double val = 0.0;
    if(betaAve<PI) val = 1/pow(kij(n,P,S,Sdij),2);         // "concave" point
    if(betaAve>=PI) val = pow(kij(n,P,S,Sdij),2);          // "convex" point
    return val;
}

//! ------------------------------
//! function: scalarFieldSmoother
//! details:
//! ------------------------------
void smoothingTools::scalarFieldSmoother(QMap<int,double> &field,
                                         const occHandle(Ng_MeshVS_DataSourceFace) &aMeshDS,
                                         const std::map<int,double> &betaAveField,
                                         const std::map<int,int> &mapOfNodeTypes,
                                         int NbSteps)
{
    cout<<"prismaticLayer::scalarFieldSmoother()->____function called____"<<endl;

    QMap<int,double> smoothedField;
    for(int n = 1; n<=NbSteps; n++)
    {
        for(QMap<int,double>::iterator it = field.begin(); it!= field.end(); ++it)
        {
            double snx = 0;
            int globalNodeID = it.key();

            //! -----------------------
            //! check the type of node
            //! -----------------------
            int nodeType = mapOfNodeTypes.find(globalNodeID)->second;
            switch(nodeType)
            {
                //! ---------------------------------------------------------------------------
                //! node fixed - nor the guiding vector nor the marching distance are smoothed
                //! use the current value of the field as smoothed value
                //! "p is in this category if p 2 F with its manifold  and if (ti,tj)>=60°
                //! The marching vector and distance of p are fixed"
                //! ---------------------------------------------------------------------------
            case 1:
            {
                double localValue = it.value();
                smoothedField.insert(globalNodeID,localValue);
            }
                break;

                //! ---------------------------------------------------------------------------------
                //! category 2: node on the prismatic/non prismatic boundary
                //! "Its marching vector and distance are weighted averages of
                //! the marching vector and step of the neighbouring points belonging to Category 2"
                //! ---------------------------------------------------------------------------------
            case 2:
            {
                //! ------------------
                //! surrounding nodes
                //! ------------------
                std::set<int> setOfSurroundingNodes;
                smoothingTools::surroundingNodes(aMeshDS,globalNodeID,false,setOfSurroundingNodes);
                /*
                int localNodeID = aMeshDS->myNodesMap.FindIndex(globalNodeID);
                const QList<int> &surroundingElements_local = aMeshDS->myNodeToElements.value(localNodeID);
                int NbSurroundingElements = surroundingElements_local.length();
                std::set<int> setOfSurroundingNodes;
                for(int n=0; n<NbSurroundingElements; n++)
                {
                    int localElementID_surr = surroundingElements_local[n];
                    int globalElementID_surr = aMeshDS->myElementsMap.FindKey(localElementID_surr);
                    int NbNodes, buf[12];
                    TColStd_Array1OfInteger nodeIDs(*buf,1,10);
                    aMeshDS->GetNodesByElement(globalElementID_surr,nodeIDs,NbNodes);
                    for(int n=1; n<=NbNodes; n++) setOfSurroundingNodes.insert(nodeIDs(n));
                }
                std::set<int>::iterator it_ = setOfSurroundingNodes.find(globalNodeID);
                setOfSurroundingNodes.erase(it_);
                */
                std::vector<int> surroundingNodesOfType2;
                for(std::set<int>::iterator it__ = setOfSurroundingNodes.begin(); it__!=setOfSurroundingNodes.end(); it__++)
                {
                    int curGlobalNodeID_surr = *it__;
                    int typeOfNode = mapOfNodeTypes.at(curGlobalNodeID_surr);
                    if(typeOfNode!=2) continue;
                    surroundingNodesOfType2.push_back(curGlobalNodeID_surr);
                }

                //! in case no surrounding node of type 2 exists, use the current
                //! value of the field as smoothed value
                if(surroundingNodesOfType2.size()==0)
                {
                    double localValue = it.value();
                    smoothedField.insert(globalNodeID,localValue);
                    break;
                }

                //! perform weighted average of displacements - use distances as weight
                //! this could be changed according to the answer of emailed question
                double bufd[3];
                int NbNodes_;
                TColStd_Array1OfReal coords(*bufd,1,3);
                MeshVS_EntityType aType;
                aMeshDS->GetGeom(globalNodeID,false,coords,NbNodes_,aType);
                double xP = coords(1);
                double yP = coords(2);
                double zP = coords(3);
                double Scurr = 0, S = 0;
                for(std::vector<int>::iterator it = surroundingNodesOfType2.begin(); it!=surroundingNodesOfType2.end(); it++)
                {
                    int curGlobalNodeID_surr_type2 = *it;
                    aMeshDS->GetGeom(curGlobalNodeID_surr_type2,false,coords,NbNodes_,aType);
                    double xn = coords(1);
                    double yn = coords(2);
                    double zn = coords(3);
                    double curDist = sqrt(pow(xP-xn,2)+pow(yP-yn,2)+pow(zP-zn,2));
                    S += sqrt(pow(xP-xn,2)+pow(yP-yn,2)+pow(zP-zn,2));
                    Scurr += curDist*field.value(curGlobalNodeID_surr_type2);
                }
                Scurr /= S;
                smoothedField.insert(globalNodeID,Scurr);

            }
                break;

                //! --------------------------------------------------------------------------------
                //! The marching distance is not smoothed, and only its marching vector is smoothed
                //! Use the current value of the field as smoothed value
                //! --------------------------------------------------------------------------------
            case 3:
            {
                double localValue = it.value();
                smoothedField.insert(globalNodeID,localValue);
            }
                break;

                //! -------------------------------------------------------------------------------
                //! free node - both the guiding vector and both the marching distance are updated
                //! "All the neighbouring points of p in its manifold  contribute to the smoothing
                //! of its marching vector and distance."
                //! -------------------------------------------------------------------------------
            case 4:
            {
                int localNodeID = aMeshDS->myNodesMap.FindIndex(globalNodeID);
                const std::vector<double> &P = aMeshDS->getNodeCoordinates(localNodeID);

                //! -------------------------------------------
                //! value of the field at the current position
                //! -------------------------------------------
                double localValue = it.value();

                //! ------------------------
                //! manifold characteristic
                //! ------------------------
                double betaAve = betaAveField.at(globalNodeID);

                //! --------------------------
                //! get the surrounding nodes
                //! --------------------------
                std::set<int> surroundingNodes;
                smoothingTools::surroundingNodes(aMeshDS,localNodeID,true,surroundingNodes);
                /*
                std::set<int> surroundingNodes;
                const QList<int> &elements = aMeshDS->myNodeToElements.value(localNodeID);
                for(int i=0; i<elements.length(); i++)
                {
                    int localElementID = elements[i];
                    int globalElementID = aMeshDS->myElementsMap.FindKey(localElementID);
                    int NbNodes, buf[20];
                    TColStd_Array1OfInteger nodeIDs(*buf,1,20);
                    aMeshDS->GetNodesByElement(globalElementID,nodeIDs,NbNodes);
                    for(int j=1; j<=NbNodes; j++) surroundingNodes.insert(nodeIDs(j));
                }
                std::set<int>::iterator itt = surroundingNodes.find(globalNodeID);
                if(itt!=surroundingNodes.end()) surroundingNodes.erase(itt);
                */

                //! ---------------------------------------------------------
                //! sum of all the distances from the current advancing node
                //! ---------------------------------------------------------
                double Sdij = 0;
                for(std::set<int>::iterator it = surroundingNodes.begin(); it!=surroundingNodes.end(); it++)
                {
                    int globalNodeID_surrounding = *it;
                    int localNodeID_surrounding = aMeshDS->myNodesMap.FindIndex(globalNodeID_surrounding);
                    const std::vector<double> &S = aMeshDS->getNodeCoordinates(localNodeID_surrounding);
                    double d = sqrt(pow(S[0]-P[0],2)+pow(S[1]-P[1],2)+pow(S[2]-P[2],2));
                    Sdij += d;
                }

                double Swij = 0;
                size_t n = surroundingNodes.size();
                for(std::set<int>::iterator it = surroundingNodes.begin(); it!=surroundingNodes.end(); it++)
                {
                    int globalNodeID_surrounding = *it;
                    int localNodeID_surrounding = aMeshDS->myNodesMap.FindIndex(globalNodeID_surrounding);
                    const std::vector<double> &S = aMeshDS->getNodeCoordinates(localNodeID_surrounding);
                    Swij += wij(betaAve,P,S,Sdij,n);
                }

                double a0 = 0.0;
                for(std::set<int>::iterator it = surroundingNodes.begin(); it!=surroundingNodes.end(); it++)
                {
                    int globalNodeID_surrounding = *it;
                    int localNodeID_surrounding = aMeshDS->myNodesMap.FindIndex(globalNodeID_surrounding);
                    const std::vector<double> &S = aMeshDS->getNodeCoordinates(localNodeID_surrounding);
                    double localValue_surrounding = field.value(globalNodeID_surrounding);
                    a0 += localValue_surrounding*wij(betaAve,P,S,Sdij,n);
                }

                //! ------------------
                //! relaxation factor
                //! ------------------
                double omega = fomega(betaAve);
                snx = (1-omega)*localValue+(omega/Swij)*a0;
                smoothedField.insert(globalNodeID,snx);
            }
                break;

                //! ---------------------------------------------------------------
                //! nodes having a neighbor which has compression due to proximity
                //! rule: use the an averaged value of the marching distance
                //! (average by distances? Inverse of distances?)
                //! ---------------------------------------------------------------
            case 5:
            {
                //! ------------------------------
                //! get current point coordinates
                //! ------------------------------
                double buff[3];
                int NbNodes;
                TColStd_Array1OfReal coordsP(*buff,1,3);
                MeshVS_EntityType aType;
                aMeshDS->GetGeom(globalNodeID,false,coordsP,NbNodes,aType);

                //! -----------------------------
                //! search for surrounding nodes
                //! -----------------------------
                std::set<int> setOfSurroundingNodes;
                smoothingTools::surroundingNodes(aMeshDS,globalNodeID,false,setOfSurroundingNodes);
                /*
                int localNodeID = aMeshDS->myNodesMap.FindIndex(globalNodeID);
                const QList<int> &surroundingElements_local = aMeshDS->myNodeToElements.value(localNodeID);
                int NbSurroundingElements = surroundingElements_local.length();
                std::set<int> setOfSurroundingNodes;
                for(int n=0; n<NbSurroundingElements; n++)
                {
                    int localElementID_surr = surroundingElements_local[n];
                    int globalElementID_surr = aMeshDS->myElementsMap.FindKey(localElementID_surr);
                    int NbNodes, buf[10];
                    TColStd_Array1OfInteger nodeIDs(*buf,1,10);
                    aMeshDS->GetNodesByElement(globalElementID_surr,nodeIDs,NbNodes);
                    for(int n=1; n<=NbNodes; n++) setOfSurroundingNodes.insert(nodeIDs(n));
                }
                setOfSurroundingNodes.erase(globalNodeID);
                */

                //! -----------------------------------
                //! iterate over the surrounding nodes
                //! -----------------------------------
                double buff_[3];
                TColStd_Array1OfReal coordsS(*buff_,1,3);
                double sumOfWeigths = 0;
                double partialSum = 0;
                for(std::set<int>::iterator it = setOfSurroundingNodes.begin(); it!=setOfSurroundingNodes.end(); it++)
                {
                    int globalNodeID_surrounding = *it;
                    aMeshDS->GetGeom(globalNodeID_surrounding,false,coordsS,NbNodes,aType);
                    double d_PS = sqrt(pow(coordsP(1)-coordsS(1),2)+pow(coordsP(2)-coordsS(2),2)+pow(coordsP(3)-coordsS(3),2));
                    //cout<<"____case5: P("<<coordsP(1)<<", "<<coordsP(2)<<", "<<coordsP(3)<<")____"<<endl;
                    //cout<<"____case5: S("<<coordsS(1)<<", "<<coordsS(2)<<", "<<coordsS(3)<<")____"<<endl;
                    //cout<<"____case 5: "<<d_PS<<"____"<<endl;
                    sumOfWeigths += d_PS;
                    double localValue = field.value(globalNodeID_surrounding);
                    partialSum += d_PS*localValue;
                }
                double localField_updated = partialSum/sumOfWeigths;
                smoothedField.insert(globalNodeID,localField_updated);
            }
                break;
            }
        }

        //! --------------------------------
        //! replace the scalar field values
        //! --------------------------------
        for(QMap<int,double>::iterator it = field.begin(); it!= field.end(); ++it)
        {
            int globalNodeID = it.key();
            it.value() = smoothedField.value(globalNodeID);
        }
    }
}

//! ------------------------------------------------------------------------------
//! function: fieldSmoother
//! details:  same of previous but working on a vectorial field
//!           provide also normalization in case only vector "rotation" is needed
//! ------------------------------------------------------------------------------
void smoothingTools::fieldSmoother(QMap<int,QList<double>> &field,
                                   const occHandle(Ng_MeshVS_DataSourceFace) &aMeshDS,
                                   const std::map<int,double> &betaAveField,
                                   const std::map<int,double> &betaVisibilityField,
                                   const std::map<int,int> &mapOfNodeTypes,
                                   int NbSteps,
                                   bool normalize)
{
    cout<<"prismaticLayer::fieldSmoother()->____function called____"<<endl;

    const double PI = 3.14159236534;
    const double dmax = 0.8;

    //! ----------------
    //! smoothing steps
    //! ----------------
    QMap<int,QList<double>> smoothedField;
    for(int n = 1; n<=NbSteps; n++)
    {
        for(QMap<int,QList<double>>::iterator it = field.begin(); it!= field.end(); ++it)
        {
            double snx, sny, snz;
            snx = sny = snz = 0;

            int globalNodeID = it.key();

            //! -----------------------
            //! check the type of node
            //! -----------------------
            int nodeType = mapOfNodeTypes.find(globalNodeID)->second;
            //cout<<"prismaticLayer::fieldSmoother()->____node category: "<<nodeType<<"____"<<endl;

            switch(nodeType)
            {
                //! ------------------------
                //! category 1: fixed nodes
                //! ------------------------
            case 1:
            {
                const QList<double> &localFieldValue = smoothedField.value(globalNodeID);
                smoothedField.insert(globalNodeID,localFieldValue);
            }
                break;

                //! --------------------------------------------------
                //! category 2: node on prismatic/non prismatic walls
                //! --------------------------------------------------
            case 2:
            {
                std::vector<int> surroundingNodesOfType2;

                //! ------------------
                //! surrounding nodes
                //! ------------------
                int localNodeID = aMeshDS->myNodesMap.FindIndex(globalNodeID);
                std::set<int> setOfSurroundingNodes;
                smoothingTools::surroundingNodes(aMeshDS,localNodeID,true,setOfSurroundingNodes);
                /*
                const QList<int> &surroundingElements_local = aMeshDS->myNodeToElements.value(localNodeID);
                int NbSurroundingElements = surroundingElements_local.length();
                std::set<int> setOfSurroundingNodes;
                for(int n=0; n<NbSurroundingElements; n++)
                {
                    int localElementID_surr = surroundingElements_local[n];
                    int globalElementID_surr = aMeshDS->myElementsMap.FindKey(localElementID_surr);
                    int NbNodes, buf[10];
                    TColStd_Array1OfInteger nodeIDs(*buf,1,10);
                    aMeshDS->GetNodesByElement(globalElementID_surr,nodeIDs,NbNodes);
                    for(int n=1; n<=NbNodes; n++) setOfSurroundingNodes.insert(nodeIDs(n));
                }
                std::set<int>::iterator it_ = setOfSurroundingNodes.find(globalNodeID);
                setOfSurroundingNodes.erase(it_);
                */

                for(std::set<int>::iterator it_ = setOfSurroundingNodes.begin(); it_!=setOfSurroundingNodes.end(); it_++)
                {
                    int curGlobalNodeID_surr = *it_;
                    int typeOfNode = mapOfNodeTypes.at(curGlobalNodeID_surr);
                    if(typeOfNode!=2) continue;
                    surroundingNodesOfType2.push_back(curGlobalNodeID_surr);
                }

                //! --------------------------------------------------------------
                //! in case no surrounding node of type 2 exists, use the current
                //! value of the field as smoothed value and break
                //! --------------------------------------------------------------
                if(surroundingNodesOfType2.size()==0)
                {
                    const QList<double> &localFieldValue = smoothedField.value(globalNodeID);
                    smoothedField.insert(globalNodeID,localFieldValue);
                    break;  // immediately exit from this case
                }

                //! ------------------------------------------------------------------------
                //! perform the weighted average of displacements - use distances as weight
                //! this could be changed according to the answer of emailed question
                //! ------------------------------------------------------------------------
                double bufd[3];
                int NbNodes_;
                TColStd_Array1OfReal coords(*bufd,1,3);
                MeshVS_EntityType aType;
                aMeshDS->GetGeom(globalNodeID,false,coords,NbNodes_,aType);
                double xP = coords(1);
                double yP = coords(2);
                double zP = coords(3);
                double Scurrx = 0, Scurry = 0, Scurrz = 0, S = 0;
                for(std::vector<int>::iterator it = surroundingNodesOfType2.begin(); it!=surroundingNodesOfType2.end(); it++)
                {
                    int curGlobalNodeID_surr_type2 = *it;
                    aMeshDS->GetGeom(curGlobalNodeID_surr_type2,false,coords,NbNodes_,aType);
                    double xn = coords(1);
                    double yn = coords(2);
                    double zn = coords(3);
                    double curDist = sqrt(pow(xP-xn,2)+pow(yP-yn,2)+pow(zP-zn,2));
                    S += sqrt(pow(xP-xn,2)+pow(yP-yn,2)+pow(zP-zn,2));

                    Scurrx += curDist*field.value(curGlobalNodeID_surr_type2)[0];
                    Scurry += curDist*field.value(curGlobalNodeID_surr_type2)[1];
                    Scurrz += curDist*field.value(curGlobalNodeID_surr_type2)[2];
                }
                Scurrx /= S;
                Scurry /= S;
                Scurrz /= S;

                QList<double> localFieldValue;
                if(normalize==true)
                {
                    double l = sqrt(pow(Scurrx,2)+pow(Scurry,2)+pow(Scurrz,2));
                    Scurrx /= l;
                    Scurry /= l;
                    Scurrz /= l;
                }
                localFieldValue<<Scurrx<<Scurry<<Scurrz;
                smoothedField.insert(globalNodeID,localFieldValue);
            }
                break;

                //! -------------------------------------------------
                //! category 3: only the marching vector is smoothed
                //! category 4: free nodes
                //! -------------------------------------------------
            case 3:
            case 4:
            {
                int localNodeID = aMeshDS->myNodesMap.FindIndex(globalNodeID);
                const std::vector<double> &P = aMeshDS->getNodeCoordinates(localNodeID);

                //! ---------------------------------------------
                //! manifold characteristic at the current point
                //! ---------------------------------------------
                double betaAve = betaAveField.at(globalNodeID);

                //! -------------------------------------------
                //! value of the field at the current position
                //! -------------------------------------------
                const QList<double> &localValue = it.value();

                //! --------------------------
                //! get the surrounding nodes
                //! --------------------------
                std::set<int> surroundingNodes;
                smoothingTools::surroundingNodes(aMeshDS,localNodeID,true,surroundingNodes);
                /*
                QList<int> elements = aMeshDS->myNodeToElements.value(localNodeID);
                for(int i=0; i<elements.length(); i++)
                {
                    int localElementID = elements[i];
                    int globalElementID = aMeshDS->myElementsMap.FindKey(localElementID);
                    int NbNodes, buf[20];
                    TColStd_Array1OfInteger nodeIDs(*buf,1,20);
                    aMeshDS->GetNodesByElement(globalElementID,nodeIDs,NbNodes);
                    for(int j=1; j<=NbNodes; j++) surroundingNodes.insert(nodeIDs(j));

                }
                std::set<int>::iterator itt = surroundingNodes.find(globalNodeID);
                if(itt!=surroundingNodes.end()) surroundingNodes.erase(itt);
                */

                //! ---------------------------------------------------------
                //! sum of all the distances from the current advancing node
                //! ---------------------------------------------------------
                double Sdij = 0;
                for(std::set<int>::iterator it = surroundingNodes.begin(); it!=surroundingNodes.end(); it++)
                {
                    int globalNodeID_surrounding = *it;
                    int localNodeID_surrounding = aMeshDS->myNodesMap.FindIndex(globalNodeID_surrounding);
                    const std::vector<double> &S = aMeshDS->getNodeCoordinates(localNodeID_surrounding);
                    double d = sqrt(pow(S[0]-P[0],2)+pow(S[1]-P[1],2)+pow(S[2]-P[2],2));
                    Sdij += d;
                }

                double Swij = 0;
                size_t n = surroundingNodes.size();
                for(std::set<int>::iterator it = surroundingNodes.begin(); it!=surroundingNodes.end(); it++)
                {
                    int globalNodeID_surrounding = *it;
                    int localNodeID_surrounding = aMeshDS->myNodesMap.FindIndex(globalNodeID_surrounding);
                    const std::vector<double> &S = aMeshDS->getNodeCoordinates(localNodeID_surrounding);
                    Swij += wij(betaAve,P,S,Sdij,n);
                }

                double a0,a1,a2;
                a0=a1=a2=0.0;
                for(std::set<int>::iterator it = surroundingNodes.begin(); it!=surroundingNodes.end(); it++)
                {
                    int globalNodeID_surrounding = *it;
                    int localNodeID_surrounding = aMeshDS->myNodesMap.FindIndex(globalNodeID_surrounding);
                    const std::vector<double> &S = aMeshDS->getNodeCoordinates(localNodeID_surrounding);
                    const QList<double> &localValue_surrounding = field.value(globalNodeID_surrounding);
                    a0 += localValue_surrounding[0]*wij(betaAve,P,S,Sdij,n);
                    a1 += localValue_surrounding[1]*wij(betaAve,P,S,Sdij,n);
                    a2 += localValue_surrounding[2]*wij(betaAve,P,S,Sdij,n);
                }

                //! ----------------------
                //! overrelaxation factor
                //! ----------------------
                double omega = fomega(betaAve);
                double oneminusomega = 1-omega;

                snx = oneminusomega*localValue[0] + (omega/Swij)*a0;
                sny = oneminusomega*localValue[1] + (omega/Swij)*a1;
                snz = oneminusomega*localValue[2] + (omega/Swij)*a2;

                if(normalize == true)
                {
                    double l = sqrt(snx*snx+sny*sny+snz*snz);
                    snx /= l;
                    sny /= l;
                    snz /= l;
                }
                QList<double> localFieldValueSmoothed;
                localFieldValueSmoothed<<snx<<sny<<snz;
                smoothedField.insert(globalNodeID,localFieldValueSmoothed);
            }
                break;
            }
        }
    }

    //! -------------------------
    //! replace the field values
    //! -------------------------
    const double maxDeviationFromPreviousGuiding = PI/6.0;
    int NbNotSmoothed = 0;
    for(QMap<int,QList<double>>::iterator it = field.begin(); it!= field.end(); ++it)
    {
        int globalNodeID = it.key();
        int localNodeID = aMeshDS->myNodesMap.FindIndex(globalNodeID);

        const QList<double> &localFieldValue = smoothedField.value(globalNodeID);
        double snx = localFieldValue[0];
        double sny = localFieldValue[1];
        double snz = localFieldValue[2];

        //! ------------------------------------------------------
        //! check deviation from the not smoothed marching vector
        //! ------------------------------------------------------
        //const QList<double> &localFieldValueNotSmoothed = field.value(globalNodeID);
        //double nsnx = localFieldValueNotSmoothed[0];
        //double nsny = localFieldValueNotSmoothed[1];
        //double nsnz = localFieldValueNotSmoothed[2];

        //double dot = snx*nsnx+sny*nsny+snz*nsnz;
        //if(dot<-1) dot = -1;
        //if(dot>1) dot = 1;
        //double deviationFromPrevious = std::acos(dot);

        //! ------------------------------------------------------------------
        //! if the smoothed vector deviates from the non smoothed vector more
        //! than the limit (typically 30°) an average value is considered
        //! ------------------------------------------------------------------
        //if(deviationFromPrevious>maxDeviationFromPreviousGuiding)
        //{
        //    snx = (snx+nsnx)/2.0; sny = (sny+nsny)/2.0; snz = (snz+nsnz)/2.0;
        //    double l = sqrt(snx*snx+sny*sny+snz*snz);
        //    snx /= l; sny /= l; snz /= l;
        //    localFieldValue.clear();
        //    localFieldValue<<snx<<sny<<snz;
        //}

        //! ------------------------------------------------------
        //! check if the smoothed vectors stand within visibility
        //! ------------------------------------------------------
        const QList<double> &aNodeNormal = aMeshDS->myNodeNormals.value(globalNodeID);          // best normal
        const QList<int> &localElementIDs_surr = aMeshDS->myNodeToElements.value(localNodeID);
        bool isSmoothedDirectionOK = true;
        for(int n=0; n<localElementIDs_surr.length(); n++)
        {
            int localElementID_surr = localElementIDs_surr[n];
            int globalElementID_surr = aMeshDS->myElementsMap.FindKey(localElementID_surr);
            double nx,ny,nz;
            aMeshDS->GetNodeNormal(globalElementID_surr,10,nx,ny,nz);
            double dot = aNodeNormal[0]*snx+aNodeNormal[1]*sny+aNodeNormal[2]*snz;      // no need to normalize here
            if(dot<-1) dot = -1;
            if(dot>1) dot = 1;
            double dpmin = cos(PI/2-(1-dmax)*betaVisibilityField.at(globalNodeID));
            if(dot<dpmin)
            {
                isSmoothedDirectionOK = false;
                break;
            }
        }
        if(isSmoothedDirectionOK==true) it.value() = smoothedField.value(globalNodeID);
    }

    cout<<"**************************************************"<<endl;
    cout<<" fraction of guiding vectors not smoothed: "<<double(NbNotSmoothed/aMeshDS->GetAllNodes().Extent())<<endl;
    cout<<"**************************************************"<<endl;
}

//! ----------------------------
//! function: surroundingPoints
//! details:  helper
//! ----------------------------
void smoothingTools::surroundingNodes(const occHandle(Ng_MeshVS_DataSourceFace) &aMeshDS, int nodeID, bool isLocal, std::set<int> &setOfSurroundingNodes)
{
    int localNodeID, globalNodeID;
    if(isLocal==true)
    {
        localNodeID = nodeID;
        globalNodeID = aMeshDS->myNodesMap.FindKey(nodeID);
    }
    else
    {
        localNodeID = aMeshDS->myNodesMap.FindIndex(nodeID);
        globalNodeID = nodeID;
    }

    const QList<int> &surroundingElements_local = aMeshDS->myNodeToElements.value(localNodeID);
    int NbSurroundingElements = surroundingElements_local.length();
    for(int n=0; n<NbSurroundingElements; n++)
    {
        int localElementID_surr = surroundingElements_local[n];
        int globalElementID_surr = aMeshDS->myElementsMap.FindKey(localElementID_surr);
        int NbNodes, buf[10];
        TColStd_Array1OfInteger nodeIDs(*buf,1,10);
        aMeshDS->GetNodesByElement(globalElementID_surr,nodeIDs,NbNodes);
        for(int n=1; n<=NbNodes; n++) setOfSurroundingNodes.insert(nodeIDs(n));
    }
    std::set<int>::iterator it_ = setOfSurroundingNodes.find(globalNodeID);
    setOfSurroundingNodes.erase(it_);
}
