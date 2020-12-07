#include "smoothingtools.h"

//! ----
//! C++
//! ----
#include <vector>

//! ------------------------------
//! function: getPointCoordinates
//! details:
//! ------------------------------
bool smoothingTools::getPointCoordinate(const occHandle(Ng_MeshVS_DataSourceFace) &aMeshDS, int globalNodeID, double *P)
{
    if(aMeshDS.IsNull()) return false;
    int NbNodes;
    MeshVS_EntityType aType;
    double buf[3];
    TColStd_Array1OfReal coords(*buf,1,3);
    bool isDone = aMeshDS->GetGeom(globalNodeID,false,coords,NbNodes,aType);
    P[0] = coords(1); P[1] = coords(2); P[2] = coords(3);
    return isDone;
}

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
    double val = 0.5*(1+cos(betaAve/2));    // (PI-betaAve/2)*(1/(PI))
    return val;
}

double kij(int n, const std::vector<double> &P, const std::vector<double> &S, double Sdij)
{
    double dij = sqrt(pow(P[0]-S[0],2)+pow(P[1]-S[1],2)+pow(P[2]-S[2],2));
    double val = n*dij/Sdij;
    return val;
}

double wij(double betaAve, const std::vector<double> &P, const std::vector<double> &S, double Sdij, int n)
{
    const double PI = 3.1415926538;
    const double power = 2.0;
    double val = 0.0;
    if(betaAve<PI) val = 1/pow(kij(n,P,S,Sdij),power);         // "concave" point
    if(betaAve>=PI) val = pow(kij(n,P,S,Sdij),-power);          // "convex" point
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
                                         int NbSteps,
                                         int inflationStep)
{
    cout<<"smoothingTools::scalarFieldSmoother()->____function called____"<<endl;

    //! --------------------------------------------------------
    //! record the initial field: these "best" normals are used
    //! for checking the visibility constraint after smoothing
    //! --------------------------------------------------------
    QMap<int,double> initialField(field);

    QMap<int,double> smoothedField;
    for(int n = 1; n<=NbSteps; n++)
    {
        for(QMap<int,double>::iterator it = field.begin(); it!= field.end(); ++it)
        {
            int globalNodeID = it.key();

            //! -----------------------
            //! check the type of node
            //! -----------------------
            int nodeCategory = mapOfNodeTypes.find(globalNodeID)->second;
            switch(nodeCategory)
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
                std::vector<int> surroundingNodesOfType2;
                for(std::set<int>::iterator it__ = setOfSurroundingNodes.begin(); it__!=setOfSurroundingNodes.end(); it__++)
                {
                    int curGlobalNodeID_surr = *it__;
                    int typeOfNode = mapOfNodeTypes.at(curGlobalNodeID_surr);
                    if(typeOfNode!=2) continue;
                    surroundingNodesOfType2.push_back(curGlobalNodeID_surr);
                }

                //! --------------------------------------------------------------
                //! in case no surrounding node of type 2 exists, use the current
                //! value of the field as smoothed value
                //! --------------------------------------------------------------
                if(surroundingNodesOfType2.size()==0)
                {
                    double localValue = it.value();
                    smoothedField.insert(globalNodeID,localValue);
                    break;
                }

                //! --------------------------------------------------------------------
                //! perform weighted average of displacements - use distances as weight
                //! this could be changed according to the answer of emailed question
                //! --------------------------------------------------------------------
                double P[3];
                smoothingTools::getPointCoordinate(aMeshDS,globalNodeID,P);
                double xP = P[0];
                double yP = P[1];
                double zP = P[2];
                double Scurr = 0, S = 0;
                for(std::vector<int>::iterator it = surroundingNodesOfType2.begin(); it!=surroundingNodesOfType2.end(); it++)
                {
                    int curGlobalNodeID_surr_type2 = *it;
                    double N[3];
                    smoothingTools::getPointCoordinate(aMeshDS,curGlobalNodeID_surr_type2,N);
                    double xn = N[0];
                    double yn = N[1];
                    double zn = N[2];
                    double curDist = sqrt(pow(xP-xn,2)+pow(yP-yn,2)+pow(zP-zn,2));
                    S+= 1/curDist;
                    Scurr += (1/curDist)*field.value(curGlobalNodeID_surr_type2);
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

                for(std::set<int>::iterator it = surroundingNodes.begin(); it!=surroundingNodes.end(); it++)
                {
                    int globalNodeID_surrounding = *it;
                    int localNodeID_surrounding = aMeshDS->myNodesMap.FindIndex(globalNodeID_surrounding);
                    const std::vector<double> &S = aMeshDS->getNodeCoordinates(localNodeID_surrounding);
                    Swij += wij(betaAve,P,S,Sdij,inflationStep);
                }

                double a0 = 0.0;
                for(std::set<int>::iterator it = surroundingNodes.begin(); it!=surroundingNodes.end(); it++)
                {
                    int globalNodeID_surrounding = *it;
                    int localNodeID_surrounding = aMeshDS->myNodesMap.FindIndex(globalNodeID_surrounding);
                    const std::vector<double> &S = aMeshDS->getNodeCoordinates(localNodeID_surrounding);
                    double localValue_surrounding = field.value(globalNodeID_surrounding);
                    a0 += localValue_surrounding*wij(betaAve,P,S,Sdij,inflationStep);
                }

                //! -----------
                //! relaxation
                //! -----------
                double omega = fomega(betaAve);
                double smoothedValue = (1-omega)*localValue+(omega/Swij)*a0;
                smoothedField.insert(globalNodeID,smoothedValue);
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
                //double P[3];
                //smoothingTools::getPointCoordinate(aMeshDS,globalNodeID,P);

                double smoothedValue;
                double betaAve = betaAveField.at(globalNodeID);
                smoothingTools::smoothScalarAtPoint(aMeshDS,field,globalNodeID,betaAve,inflationStep,smoothedValue);
                field.insert(globalNodeID,smoothedValue);

                /*
                //! -----------------------------
                //! search for surrounding nodes
                //! -----------------------------
                std::set<int> setOfSurroundingNodes;
                smoothingTools::surroundingNodes(aMeshDS,globalNodeID,false,setOfSurroundingNodes);

                //! -----------------------------------
                //! iterate over the surrounding nodes
                //! -----------------------------------
                double sumOfWeigths = 0;
                double weight = 0;
                double partialSum = 0;
                for(std::set<int>::iterator it = setOfSurroundingNodes.begin(); it!=setOfSurroundingNodes.end(); it++)
                {
                    int globalNodeID_surrounding = *it;
                    double S[3];
                    smoothingTools::getPointCoordinate(aMeshDS,globalNodeID_surrounding,S);
                    double d_PS = sqrt(pow(P[0]-S[0],2)+pow(P[1]-S[1],2)+pow(P[2]-S[2],2));
                    weight = pow(d_PS,-2);
                    sumOfWeigths += weight;
                    double localValue = field.value(globalNodeID_surrounding);
                    partialSum += weight*localValue;
                }                

                double localField_updated = partialSum/sumOfWeigths;
                smoothedField.insert(globalNodeID,localField_updated);
                */
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
    cout<<"smoothingTools::scalarFieldSmoother()->____smoothing done____"<<endl;
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
                                   bool normalize,
                                   int inflationStep)
{
    cout<<"smoothingTools::fieldSmoother()->____function called____"<<endl;

    const double PI = 3.14159236534;
    const double dmax = 0.8;

    //! -------------------------
    //! record the initial field
    //! -------------------------
    QMap<int,QList<double>> initialField(field);

    //! ---------------
    //! smoothed field
    //! ---------------
    QMap<int,QList<double>> smoothedField;

    //! ----------------
    //! smoothing steps
    //! ----------------
    for(int n = 1; n<=NbSteps; n++)
    {
        for(QMap<int,QList<double>>::iterator it = field.begin(); it!= field.end(); ++it)
        {
            int globalNodeID = it.key();

            //! -----------------------
            //! check the type of node
            //! -----------------------
            int nodeCategory = mapOfNodeTypes.find(globalNodeID)->second;
            switch(nodeCategory)
            {
                //! ------------------------
                //! category 1: fixed nodes
                //! ------------------------
            case 1:
            {
                //const QList<double> &localFieldValue = smoothedField.value(globalNodeID);     // ERR - left here for documentation
                const QList<double> &localFieldValue = initialField.value(globalNodeID);
                smoothedField.insert(globalNodeID,localFieldValue);
            }
                break;

                //! --------------------------------------------------
                //! category 2: node on prismatic/non prismatic walls
                //! --------------------------------------------------
            case 2:
            {
                //! ------------------
                //! surrounding nodes
                //! ------------------
                int localNodeID = aMeshDS->myNodesMap.FindIndex(globalNodeID);
                std::set<int> setOfSurroundingNodes;
                smoothingTools::surroundingNodes(aMeshDS,localNodeID,true,setOfSurroundingNodes);

                //! -------------------------------------------------------
                //! nodes surrounding the current, belonging to category 2
                //! -------------------------------------------------------
                std::vector<int> surroundingNodesOfType2;
                for(std::set<int>::iterator it_ = setOfSurroundingNodes.begin(); it_!=setOfSurroundingNodes.end(); it_++)
                {
                    int curGlobalNodeID_surr = *it_;
                    int typeOfNode = mapOfNodeTypes.at(curGlobalNodeID_surr);
                    if(typeOfNode!=2) continue;
                    surroundingNodesOfType2.push_back(curGlobalNodeID_surr);
                }

                //! --------------------------------------------------------------
                //! in case no surrounding node of type 2 exists, use the current
                //! value of the field as smoothed value and immediately break
                //! --------------------------------------------------------------
                if(surroundingNodesOfType2.size()==0)
                {
                    const QList<double> &localFieldValue = smoothedField.value(globalNodeID);
                    smoothedField.insert(globalNodeID,localFieldValue);
                    break;
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
                    S += 1/sqrt(pow(xP-xn,2)+pow(yP-yn,2)+pow(zP-zn,2));

                    Scurrx += (1/curDist)*field.value(curGlobalNodeID_surr_type2)[0];
                    Scurry += (1/curDist)*field.value(curGlobalNodeID_surr_type2)[1];
                    Scurrz += (1/curDist)*field.value(curGlobalNodeID_surr_type2)[2];
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
                //double P[3];
                //smoothingTools::getPointCoordinate(aMeshDS,globalNodeID,P);
                std::vector<double> &P = aMeshDS->getNodeCoordinates(localNodeID);

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

                //! ---------------------------------------------------------
                //! sum of all the distances from the current advancing node
                //! ---------------------------------------------------------
                double Sdij = 0;
                for(std::set<int>::iterator it = surroundingNodes.begin(); it!=surroundingNodes.end(); it++)
                {
                    int globalNodeID_surrounding = *it;
                    double S[3];
                    smoothingTools::getPointCoordinate(aMeshDS,globalNodeID_surrounding,S);
                    double d = sqrt(pow(S[0]-P[0],2)+pow(S[1]-P[1],2)+pow(S[2]-P[2],2));
                    Sdij += d;
                }

                double Swij = 0;
                for(std::set<int>::iterator it = surroundingNodes.begin(); it!=surroundingNodes.end(); it++)
                {
                    int globalNodeID_surrounding = *it;
                    int localNodeID_surrounding = aMeshDS->myNodesMap.FindIndex(globalNodeID_surrounding);
                    const std::vector<double> &S = aMeshDS->getNodeCoordinates(localNodeID_surrounding);
                    Swij += wij(betaAve,P,S,Sdij,inflationStep);
                }

                double a0,a1,a2;
                a0=a1=a2=0.0;
                for(std::set<int>::iterator it = surroundingNodes.begin(); it!=surroundingNodes.end(); it++)
                {
                    int globalNodeID_surrounding = *it;
                    int localNodeID_surrounding = aMeshDS->myNodesMap.FindIndex(globalNodeID_surrounding);
                    const std::vector<double> &S = aMeshDS->getNodeCoordinates(localNodeID_surrounding);
                    const QList<double> &localValue_surrounding = field.value(globalNodeID_surrounding);
                    a0 += localValue_surrounding[0]*wij(betaAve,P,S,Sdij,inflationStep);
                    a1 += localValue_surrounding[1]*wij(betaAve,P,S,Sdij,inflationStep);
                    a2 += localValue_surrounding[2]*wij(betaAve,P,S,Sdij,inflationStep);
                }

                //! ----------------------
                //! overrelaxation factor
                //! ----------------------
                double omega = fomega(betaAve);
                double oneminusomega = 1-omega;

                double snx = oneminusomega*localValue[0] + (omega/Swij)*a0;
                double sny = oneminusomega*localValue[1] + (omega/Swij)*a1;
                double snz = oneminusomega*localValue[2] + (omega/Swij)*a2;

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

        //! -------------------------
        //! replace the field values
        //! -------------------------
        for(QMap<int,QList<double>>::iterator it = field.begin(); it!= field.end(); ++it)
        {
            int globalNodeID = it.key();
            field.insert(globalNodeID,smoothedField.value(globalNodeID));
        }
        /*
        const double maxDeviationFromPreviousGuiding = PI/6.0;
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
            //it.value() = smoothedField.value(globalNodeID);
        }
        */
    }

    //! -------------------------------------------------------------------------------------
    //! final checks
    //! 1. the smoothed guiding direction cannot deviate from the previous more than a limit
    //! -------------------------------------------------------------------------------------
    int NbTooDeviated = 0;
    for(QMap<int,QList<double>>::iterator it = field.begin(); it!= field.end(); ++it)
    {
        int globalNodeID = it.key();
        const QList<double> &previousDir = initialField.value(globalNodeID);
        const QList<double> &smoothedDir = it.value();

        double nx_old = previousDir[0];
        double ny_old = previousDir[1];
        double nz_old = previousDir[2];
        double nx = smoothedDir[0];
        double ny = smoothedDir[1];
        double nz = smoothedDir[2];

        double l = sqrt(nx_old*nx_old+ny_old*ny_old+nz_old*nz_old)*sqrt(nx*nx+ny*ny+nz*nz);
        double dot = nx*nx_old+ny*ny_old+nz*nz_old;
        dot /= l;
        if(dot>1) dot = 1;
        if(dot<-1) dot = -1;
        double curDeviation = std::acos(dot);
        if(curDeviation>30*PI/180)
        {
            //cout<<"____deviation: "<<curDeviation*180/PI<<"____"<<endl;
            NbTooDeviated++;
            nx = (nx+nx_old)*0.5;
            ny = (ny+ny_old)*0.5;
            nz = (nz+nz_old)*0.5;
            double l = sqrt(nx*nx+ny*ny+nz*nz);
            nx /= l;
            ny /= l;
            nz /= l;
            QList<double> averaged; averaged<<nx<<ny<<nz;
            field.insert(globalNodeID,averaged);
        }
    }
    cout<<"smoothingTools::fieldSmoother()->____guiding vectors too deviated: "<<double(NbTooDeviated)/aMeshDS->GetAllNodes().Extent()<<" %____"<<endl;

    //! --------------------
    //! final checks
    //! 2. visibility check
    //! --------------------
    int NbNotSmoothed = 0;
    for(QMap<int,QList<double>>::iterator it = field.begin(); it!= field.end(); ++it)
    {
        int globalNodeID = it.key();
        const QList<double> &previousDir = initialField.value(globalNodeID);
        const QList<double> &smoothedDir = it.value();
        double nx_old = previousDir[0];
        double ny_old = previousDir[1];
        double nz_old = previousDir[2];
        double nx = smoothedDir[0];
        double ny = smoothedDir[1];
        double nz = smoothedDir[2];

        double l = sqrt(nx_old*nx_old+ny_old*ny_old+nz_old*nz_old)*sqrt(nx*nx+ny*ny+nz*nz);
        double dot = nx*nx_old+ny*ny_old+nz*nz_old;
        dot /= l;
        if(dot>1) dot = 1;
        if(dot<-1) dot = -1;
        double curDeviation = std::acos(dot);
        if(curDeviation>betaVisibilityField.at(globalNodeID))
        {
            NbNotSmoothed++;
            field.insert(globalNodeID,initialField.value(globalNodeID));
        }
    }
    cout<<"smoothingTools::fieldSmoother()->____guiding vectors not smoothed: "<<double(NbNotSmoothed)/aMeshDS->GetAllNodes().Extent()<<" %____"<<endl;
    cout<<"smoothingTools::fieldSmoother()->____smoothing done____"<<endl;
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

//! ------------------------------
//! function: smoothScalarAtPoint
//! details:
//! ------------------------------
void smoothingTools::smoothScalarAtPoint(const occHandle(Ng_MeshVS_DataSourceFace) &aMeshDS,
                                         const QMap<int,double> &field,
                                         int globalNodeID,
                                         double betaAve,
                                         int inflationStep,
                                         double &fieldValue)
{
    //cout<<"smoothingTools::smoothScalarAtPoint()->_____smoothing the marching distance at node ID: "<<globalNodeID<<"____"<<endl;

    //! ---------------------------------
    //! coordinates of the current point
    //! ---------------------------------
    int localNodeID = aMeshDS->myNodesMap.FindIndex(globalNodeID);
    const std::vector<double> &P = aMeshDS->getNodeCoordinates(localNodeID);

    std::set<int> surroundingNodes;
    smoothingTools::surroundingNodes(aMeshDS,globalNodeID,false,surroundingNodes);

    //! ---------------------------------------------------------
    //! sum of all the distances from the current advancing node
    //! ---------------------------------------------------------
    double Sdij = 0;
    for(std::set<int>::iterator it = surroundingNodes.begin(); it!=surroundingNodes.end(); it++)
    {
        int globalNodeID_surrounding = *it;
        double S[3];
        smoothingTools::getPointCoordinate(aMeshDS,globalNodeID_surrounding,S);
        double d = sqrt(pow(S[0]-P[0],2)+pow(S[1]-P[1],2)+pow(S[2]-P[2],2));
        Sdij += d;
    }

    double Swij = 0;
    for(std::set<int>::iterator it = surroundingNodes.begin(); it!=surroundingNodes.end(); it++)
    {
        int globalNodeID_surrounding = *it;
        int localNodeID_surrounding = aMeshDS->myNodesMap.FindIndex(globalNodeID_surrounding);
        const std::vector<double> &S = aMeshDS->getNodeCoordinates(localNodeID_surrounding);
        Swij += wij(betaAve,P,S,Sdij,inflationStep);
    }

    //! -----------
    //! relaxation
    //! -----------
    double omega = fomega(betaAve);
    double oneminusomega = 1-omega;
    for(int n=0; n<10; n++)
    {
        double a0=0.0;
        for(std::set<int>::iterator it = surroundingNodes.begin(); it!=surroundingNodes.end(); it++)
        {
            int globalNodeID_surrounding = *it;
            int localNodeID_surrounding = aMeshDS->myNodesMap.FindIndex(globalNodeID_surrounding);
            const std::vector<double> &S = aMeshDS->getNodeCoordinates(localNodeID_surrounding);
            double localValue_surrounding = field.value(globalNodeID_surrounding);
            a0 += localValue_surrounding*wij(betaAve,P,S,Sdij,inflationStep);
        }
        fieldValue = oneminusomega*fieldValue+(omega/Swij)*a0;
    }
}
