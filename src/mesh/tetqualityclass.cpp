//! ----------------
//! custom includes
//! ----------------
#include "tetqualityclass.h"

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
TetQualityClass::TetQualityClass()
{
    //! default initialization
    x0 = y0 = z0 = x1 = y1 = z1 = x2 = y2 = z2 = x3 = y3 = z3 = 0;
}

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
TetQualityClass::TetQualityClass(double ax0, double ay0, double az0,
                                 double ax1, double ay1, double az1,
                                 double ax2, double ay2, double az2,
                                 double ax3, double ay3, double az3):
    x0(ax0),y0(ay0),z0(az0),x1(ax1),y1(ay1),z1(az1),x2(ax2),y2(ay2),z2(az2),x3(ax3),y3(ay3),z3(az3)
{
    ;
}

//! --------------------
//! function: setPoints
//! details:
//! --------------------
void TetQualityClass::setPoints(const std::vector<mesh::meshPoint> &aTet)
{
    x0 = aTet[0].x; y0 = aTet[0].y; z0 = aTet[0].z;
    x1 = aTet[1].x; y1 = aTet[1].y; z1 = aTet[1].z;
    x2 = aTet[2].x; y2 = aTet[2].y; z2 = aTet[2].z;
    x3 = aTet[3].x; y3 = aTet[3].y; z3 = aTet[3].z;
}

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
TetQualityClass::TetQualityClass(const std::vector<mesh::meshPoint> &aTet)
{
    //! ----------------------------------
    //! Jacobian matrix of the affine map
    //! ----------------------------------
    x0 = aTet[0].x; y0 = aTet[0].y; z0 = aTet[0].z;
    x1 = aTet[1].x; y1 = aTet[1].y; z1 = aTet[1].z;
    x2 = aTet[2].x; y2 = aTet[2].y; z2 = aTet[2].z;
    x3 = aTet[3].x; y3 = aTet[3].y; z3 = aTet[3].z;
}

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
TetQualityClass::TetQualityClass(const std::vector<polygon::Point> &aTet)
{
    //! ----------------------------------
    //! Jacobian matrix of the affine map
    //! ----------------------------------
    x0 = aTet[0].x; y0 = aTet[0].y; z0 = aTet[0].z;
    x1 = aTet[1].x; y1 = aTet[1].y; z1 = aTet[1].z;
    x2 = aTet[2].x; y2 = aTet[2].y; z2 = aTet[2].z;
    x3 = aTet[3].x; y3 = aTet[3].y; z3 = aTet[3].z;
}

//! ----------------------------------
//! function: invertTetForCalculation
//! details:
//! ----------------------------------
void TetQualityClass::invertTetForCalculation()
{
    ;
}


//! ----------------------------
//! function: getQualityMeasure
//! details:
//! ----------------------------
void TetQualityClass::getQualityMeasure(double &q0, double &q1, double &q2, double &V)
{
    //! -----------------------------------------
    //! Jacobiam matrix of the reference element
    //! -----------------------------------------
    Eigen::MatrixXd W(3,3);
    W<<1.0,0.5,0.5,0.0,sqrt(3.0)/2.0,sqrt(3.0)/6.0,0.0,0.0,sqrt(2.0/3.0);

    Eigen::Matrix<double,3,3> A;
    A<<x1-x0,x2-x0,x3-x0,y1-y0,y2-y0,y3-y0,z1-z0,z2-z0,z3-z0;

    //! ----------------
    //! Jacobian matrix
    //! ----------------
    Eigen::MatrixXd S(3,3);
    S = A*(W.inverse());

    double F = S.norm();
    double F1 = S.adjoint().norm();

    //! -----------------------------------
    //! determinant of the Jacobian matrix
    //! -----------------------------------
    double sigma = S.determinant();

    //! --------------------------
    //! Modified mean-ratio (MMR)
    //! --------------------------
    q0 = (3.0*pow(h(sigma),2.0/3.0))/pow(F,2);

    //! --------------------------------
    //! Modified condition number (MCN)
    //! --------------------------------
    q1 = (3.0*h(sigma))/(F*F1);    

    //! ---------------------------------
    //! root mean square of edge lengths
    //! ---------------------------------
    double l01 = pow(x1-x0,2)+pow(y1-y0,2)+pow(z1-z0,2);
    double l02 = pow(x2-x0,2)+pow(y2-y0,2)+pow(z2-z0,2);
    double l03 = pow(x3-x0,2)+pow(y3-y0,2)+pow(z3-z0,2);
    double l12 = pow(x2-x1,2)+pow(y2-y1,2)+pow(z2-z1,2);
    double l13 = pow(x3-x1,2)+pow(y3-y1,2)+pow(z3-z1,2);
    double l23 = pow(x3-x2,2)+pow(y3-y2,2)+pow(z3-z2,2);

    double lrms = sqrt((1.0/6.0)*(l01+l02+l03+l12+l13+l23));

    //! ------------------------------
    //! Modified volume-length (MIVL)
    //! ------------------------------
    q2 = (12.0*h(Volume()))/(sqrt(2.0)*(lrms*lrms*lrms));

    //! -------
    //! Volume
    //! -------
    V = Volume();

    //cout<<"TetQualityClass::getQualityMeasure()->____("<<q0<<", "<<q1<<", "<<q2<<")____"<<endl;
}

//! --------------------------------------------------------------------------------------
//! function: h
//! details:  "Simultaneous untangling and smoothing of tetrahedral meshes"
//!           J.M. Escobar, E. Rodriguez, R. Montenegro∗, G. Montero, J.M. Gonzalez-Yuste
//! --------------------------------------------------------------------------------------
double TetQualityClass::h(double sigma)
{
    //! ----------------
    //! machine epsilon
    //! ----------------
    const double gamma = std::numeric_limits<double>::epsilon();
    double delta_min = 0;
    if(sigma<gamma)
    {
        delta_min = sqrt(gamma*(gamma-sigma));
    }
    double val = 0.5*(sigma+sqrt(sigma*sigma+4*delta_min*delta_min));
    return val;
}

//! -----------------
//! function: Volume
//! details:
//! -----------------
double TetQualityClass::Volume()
{
    //! --------------
    //! x0  y0  z0  1
    //! x1  y1  z1  1
    //! x2  y2  z2  1
    //! x3  y3  z3  1
    //! --------------
    Eigen::MatrixXd M(4,4);
    M<<x0,y0,z0,1.0,x1,y1,z1,1.0,x2,y2,z2,1.0,x3,y3,z3,1.0;
    return M.determinant()/6.0;
}
