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
    //double val = 0.5*(1-cos(betaAve));
    const double PI = 3.1415926538;
    double val = (PI-betaAve/2)*(1/(PI));
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
    const double eps = 5.0*PI/180.0;
    double val = 1.0;
    if(betaAve<PI-eps) val = 1/pow(kij(n,P,S,Sdij),2);        // "concave" point
    if(betaAve>PI+eps) val = pow(kij(n,P,S,Sdij),2);          // "convex" point
    return val;
}

//! ------------------------------
//! function: scalarFieldSmoother
//! details:
//! ------------------------------
void smoothingTools::scalarFieldSmoother(QMap<int,double> &field,
                                         const occHandle(Ng_MeshVS_DataSourceFace) &aMeshDS,
                                         const std::map<int,double> &betaAveField,
                                         int NbSteps)
{
    cout<<"prismaticLayer::scalarFieldSmoother()->____function called____"<<endl;

    for(int n = 1; n<=NbSteps; n++)
    {
        QMap<int,double> smoothedField;
        for(QMap<int,double>::iterator it = field.begin(); it!= field.end(); ++it)
        {
            double snx = 0;
            int globalNodeID = it.key();
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

            //! -----------------------
            //! over relaxation factor
            //! -----------------------
            double omega = fomega(betaAve);
            snx = (1-omega)*localValue+(omega/Swij)*a0;
            smoothedField.insert(globalNodeID,snx);
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

//! ----------------------------------------------------------------------------
//! function: fieldSmoother
//! details:  same of previous but working on a vectorial field
//!           provide also normalization in case only vector rotation is needed
//! ----------------------------------------------------------------------------
void smoothingTools::fieldSmoother(QMap<int,QList<double>> &field,
                                   const occHandle(Ng_MeshVS_DataSourceFace) &aMeshDS,
                                   const std::map<int,double> &betaAveField,
                                   int NbSteps,
                                   bool normalize)
{
    cout<<"prismaticLayer::fieldSmoother()->____function called____"<<endl;

    //! ---------------------------------------------------------------------
    //! check if the mesh has curvature data inside - use gaussian curvature
    //! ---------------------------------------------------------------------
    if(aMeshDS->myCurvature.isEmpty()) aMeshDS->computeDiscreteCurvature(1);

    //! ----------------
    //! smoothing steps
    //! ----------------
    for(int n = 1; n<=NbSteps; n++)
    {
        QMap<int,QList<double>> smoothedField;
        for(QMap<int,QList<double>>::iterator it = field.begin(); it!= field.end(); ++it)
        {
            double snx, sny, snz;
            snx = sny = snz = 0;

            int globalNodeID = it.key();
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

            snx = (1-omega)*localValue[0] + (omega/Swij)*a0;
            sny = (1-omega)*localValue[1] + (omega/Swij)*a1;
            snz = (1-omega)*localValue[2] + (omega/Swij)*a2;

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

        //! -------------------------
        //! replace the field values
        //! -------------------------
        for(QMap<int,QList<double>>::iterator it = field.begin(); it!= field.end(); ++it)
        {
            int globalNodeID = it.key();
            it.value() = smoothedField.value(globalNodeID);
        }
    }
}
