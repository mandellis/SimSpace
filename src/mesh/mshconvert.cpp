//! ----------------
//! custom includes
//! ----------------
#include "mshConvert.h"

//! -------
//! libigl
//! -------
#include <ext/libigl/include/igl/readMSH.h>
//#include <libigl/include/igl/writeMSH.h>

//! -------
//! pymesh
//! -------
#include <ext/pymesh/MshLoader.h>

//! ----
//! C++
//! ----
#include <iostream>
#include <fstream>
using namespace std;

//! ---------------------------------------------------------------------------------
//! function: constructor
//! details:  simply a container for libraries able to read .msh files (text or bin)
//!           here libigl is used
//! ---------------------------------------------------------------------------------
mshConvert::mshConvert(const QString &absoluteFilePath):myMshInputFile(absoluteFilePath)
{
    ;
}

//! ------------------------------------
//! function: toEigen
//! details:  convert to Eigen matrices
//! ------------------------------------
bool mshConvert::toEigen(Eigen::MatrixXd &V, Eigen::MatrixXi &T)
{
    try
    {
        //! -------------------
        //! nodes and elements
        //! -------------------
        igl::readMSH(myMshInputFile.toStdString(),V,T);
    }
    catch(...)
    {
        cerr<<"____error in converting the .msh file____"<<endl;
        return false;
    }
    return true;
}
