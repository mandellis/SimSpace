//! ----------------
//! custom includes
//! ----------------
#include "frdreader.h"
#include "mydefines.h"
#include "meshdatabase.h"
#include "tools.h"

//! ---
//! Qt
//! ---
#include <QProgressBar>
#include <QMessageBox>
#include <QMap>

//! ----
//! C++
//! ----
#include <iostream>
#include <fstream>
#include <string>
using namespace std;

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
FrdReader::FrdReader(QObject *parent):QObject(parent),myFRDfileName(""),myPath("")
{
    ;
}

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
FrdReader::FrdReader(const QString &FRDfileName, QObject *parent):
    QObject(parent),
    myFRDfileName(FRDfileName),
    myPath("")
{
}

//! ------------------------------------------------------------
//! function: splitFrdFile
//! details:  splt the .frd file up to the final available time
//! ------------------------------------------------------------
bool FrdReader::perform()
{
    cout<<"FrdReader::splitFrdFile()->____function called____"<<endl;

    QString FRDname = myFRDfileName.split("/").last();
    QString resPath = myFRDfileName;
    resPath.chop(FRDname.length());

    QDir dir(resPath);
    if(dir.cd("ResultsData")==true)
    {
        resPath = dir.absolutePath();
        //! empty the previous
        tools::clearDir(resPath);
    }
    else
    {
        dir.mkdir("ResultsData");
        dir.cd("ResultsData");
        resPath = dir.absolutePath();
    }

    //! open the .frd file
    cout<<"FrdReader::splitFrdFile()->____opening file: "<<myFRDfileName.toStdString()<<"____"<<endl;
    ifstream is;
    is.open(myFRDfileName.toStdString());
    if(is.is_open()==false)
    {
        cerr<<"FrdReader::splitFrdFile()->____cannot open the .frd results file: \""<<myFRDfileName.toStdString()<<"\"____"<<endl;
    }

    //! ----------------------------------------------------------------
    //! 2C = begin node list
    //! -1 = in node list is the node number and its coordinates
    //! 3C = begin element list
    //! -1 = information about element number and type of element
    //! -2 = element composition (node number)
    //! -3 = end of section
    //! >    1Ctest             !''1C'' defines a new calc of name ''test''
    //!   1UDATE  26.01.2000 !''1U'' stores user job-informations
    //!    2C                 !''2C'' starts a block of node coordinates
    //! -1    1 0.00000E+00 0.00000E+00 0.00000E+00 ! 1. node
    //! -1    2 0.10000E+01 0.00000E+00 0.00000E+00 ! 2. node ....
    //! -3                    !end of the current block
    //!    3C                 !''3C'' starts a block of elem definitions
    //! -1    1    4    0    0!first elem, type of that elem is 4 (he20)
    //! -2    1    2    3    4   13   14   15   16    5    6    7    8 ..
    //! -2   12   17   18   19   20 !twenty nodes defining that element
    //! -1 ...
    //! -2 ...
    //! -2 ...
    //! -3                    !end of the current block
    //!   1PHID      10      !defines a parameter HID value 10
    //! 100CL101            !''100C'' starts a user defined result block
    //! -4  DISP        3    1         !Attribute Header Record (Dataset)
    //! -5  D1          1    2    1    0 !Component Def. Record  (Entity)
    //! -5  D2          1    2    2    0
    //! -5  D3          1    2    3    0
    //! -1    1 0.00000E+00 1.00000E+00 1.00000E+00  !Nodal Values
    //! -1    2 1.00000E+00 0.00000E+00 0.00000E+00
    //! -3                      !end of the current block
    //! 9999                    !end of data
    //! -----------------------------------------------------------------

    //! skip the header
    std::string val;

    int Nnodes;
    for(int n=0;;n++)
    {
        std::getline(is,val);
        if(is.eof()) return false;

        int m2;
        char c[24];
        sscanf(val.c_str(),"%s%d",c,&m2);

        //! capture Nnodes
        if(strcmp("2C",c)==0)
        {
            Nnodes=m2;
            break;
        }
    }

    //! reading nodes
    for (int n=0;n<(Nnodes+1);n++)
    {
        std::getline(is,val);
        int m1,i;
        double x,y,z;
        sscanf(val.c_str(),"%d%d%lf%lf%lf",&m1,&i,&x,&y,&z);
        std::string a=std::to_string(m1);
        if(strcmp("-3",a.c_str())==0)
        {
            break;
        }
    }

    //! capture the number of elements
    int b;
    char s[24];
    std::getline(is,val);
    sscanf(val.c_str(),"%s%d",s,&b);
    int Nel;
    if (strcmp("3C",s)==0)
    {
        Nel=b;
    }
    else
    {
        return false;
    }

    //! -----------------------------------------------
    //! reading element section (skip element section)
    //! -----------------------------------------------
    int NN=Nel*2+1;
    for (int n=0;n<NN;n++)
    {
        std::getline(is,val);
        int m1;
        sscanf(val.c_str(),"%d",&m1);
        if(m1==-3)
        {
            break;
        }
    }

    //! -----------------------------------------------------------------------------------------
    //! create the discreteTimeMap for results evaluation by time or by set (used in "postools")
    //! details: to remove
    //! -----------------------------------------------------------------------------------------
    //int subStep,step;
    while(readResults(is,resPath/*,subStep,step,time*/) == true)
    {
        std::getline(is,val);
        int tl; //termination line, =-3, begin new result section
        sscanf(val.c_str(),"%d",&tl);
        if (tl!=-3)
        {
            cout<<"FrdReader::splitFRD()->____error in reading result block____"<<endl;
            break;
        }
    }

    //! --------------------
    //! close the .frd file
    //! --------------------
    cout<<"FrdReader::splitFrdFile()->____closing the .frd file____"<<endl;
    is.close();
    return true;
}

//! -----------------------
//! function: read results
//! details:
//! -----------------------
bool FrdReader::readResults(ifstream &is, QString path/*, int &sb,int &st,double &t*/) //! modified Mandelli 05/11/2018
{
    cout<<"FrdReader::readResults()->____function called____"<<endl;

    //! --------------------------------------------------------------------------------
    //! read 1st line  "1PSTEP    1(dataset n°)           1(substep n°)     1(step n°)"
    //! --------------------------------------------------------------------------------
    std::string val;
    std::getline(is,val);
    int c,sb,st;
    double t;
    char s[24],dtN[4];
    std::string s_end1,s_end2;
    s_end1="9999";
    s_end2="99999999";
    sscanf(val.c_str(),"%s%s%d%d",s,dtN,&sb,&st);

    if(strcmp(s,s_end1.c_str())==0 || strcmp(s,s_end2.c_str())==0)
    {
        cout<<"end of results data"<<endl;
        return false;
    }

    //! ----------------------------------------------
    //! read Frequency/modal (block before results)
    //!     1PGM                1.000000E+00
    //!     1PGK                1.322220E+09
    //!     1PHID                         -1
    //!     1PSUBC                         0
    //!     1PMODE                         1
    //! -------------------------------------------------
    std::getline(is,val);
    std::vector<double> dval;
    double d;
    char cd[24];
    if(sscanf(val.c_str(),"%s%lf",cd,&d)==2)
    {
        if(QString(cd) =="1PGM")
        {
            dval.push_back(d);
            for(int i=0;i<4;i++)
            {
                double e;
                std::getline(is,val);
                sscanf(val.c_str(),"%s%lf",cd,&e);
                dval.push_back(e);
            }
            std::getline(is,val);
        }
    }
    //! ------------------------------------------------------------------------------------------------------------------
    //! read 2nd line  "100CL  101(dataset n°) 2.00000E-01 (current time)        9 (n° of data == Nnodes)    0    1    1"
    //! ------------------------------------------------------------------------------------------------------------------
    //std::getline(is,val);
    int nd,b1,b2,b3;
    char ct[24];
    if(sscanf(val.c_str(),"%s%d%lf%d%d%d%d",ct,&b1,&t,&nd,&b1,&b2,&b3)==7)
    {
        //! ------------------------------------------------------------------------
        //! read 3rd line  "-4  DISP (type of data)       4 (n° of variables)    1 "
        //! ------------------------------------------------------------------------
        std::getline(is,val);
        int c1,nv,c2;
        char tdata[24];
        if(sscanf(val.c_str(),"%d%s%d%d",&c1,tdata,&nv,&c2)==4)
        {
            //!-------------------------------------------
            //! Possibile type of data:
            //! 1) Displacement, 3 variables
            //! 2) Stress, 6 variables
            //! 3) Temperature, 1 variables
            //! 4) Equivalent Plastic strain, 1 variables
            //! 5) Strain, 6 variables
            //! 6) Force, 3 variables
            //! -------------------------------------------

            //! nodal results
            int ni;
            double cxx,cyy,czz,cxy,cxz,cyz;
            double v;

            QString resultFileName=path.append("/").append(dtN).append(".rst");
            cout<<"FrdReader::readResults()->____creating file: "<<resultFileName.toStdString()<<"____"<<endl;

            fstream aResultFile(resultFileName.toStdString(),ios::out);
            aResultFile.setf(ios::scientific);
            aResultFile.precision(EXPFORMAT_PRECISION);

            aResultFile<<"Time="<<t<<endl;
            aResultFile<<"Substep n="<<sb<<" Step n="<<st<<endl;
            aResultFile<<tdata<<endl;

            //! -----------------------------------
            //! read displacement, force reactions
            //! -----------------------------------
            if(strcmp(tdata,"DISP")==0 || strcmp(tdata,"FORC")==0 || strcmp(tdata,"FLUX")==0)
            {
                for(int i=0;i<nv;i++)
                {
                    std::getline(is,val);
                }
                bool isOK = false;
                if(nd!=0)
                {
                for(int i=0;i<nd;i++)
                {
                    std::getline(is,val);
                    if(sscanf(val.c_str(),"%d%d%lf%lf%lf",&c,&ni,&cxx,&cyy,&czz)==5)
                    {
                        isOK = true;
                        aResultFile<<ni<<" "<<cxx<<" "<<cyy<<" "<<czz<<endl;
                    }
                    else
                    {
                        //! incomplete data
                        //! consider using throw()
                        isOK = false;
                        break;
                    }
                }
                }
                else
                {
                    isOK=true;
                }
                aResultFile.close();
                return isOK;
            }
            //! ---------------------------------------------------------
            //! read nodal temperature, energy, equivalent plastic strain
            //! ---------------------------------------------------------
            if(strcmp(tdata,"NDTEMP")==0 ||
                    strcmp(tdata,"ENER")==0 ||
                    strcmp(tdata,"PE")==0 ||
                    strcmp(tdata,"CELS")==0 ||
                    strcmp(tdata,"ERROR")==0 ||
                    strcmp(tdata,"HERROR")==0)
            {
                cout<<"FrdReader::readResults()->____reading a scalar from .frd____"<<endl;
                for(int i=0;i<nv;i++)
                {
                    std::getline(is,val);
                }
                bool isOK = false;
                if(nd!=0)
                {
                    for(int i=0;i<nd;i++)
                    {
                        std::getline(is,val);
                        if(sscanf(val.c_str(),"%d%d%lf",&c,&ni,&v)==3)
                        {
                            isOK = true;
                            aResultFile<<ni<<" "<<v<<endl;
                        }
                        else
                        {
                            //! incomplete data
                            //! consider using throw()
                            isOK = false;
                            break;
                        }
                    }
                }
                else
                {
                    isOK=true;
                }
                aResultFile.close();
                return isOK;
            }
            //! -----------------------------------------------------
            //! read stress, mech strain, total strain, contact data
            //! -----------------------------------------------------
            if(strcmp(tdata,"STRESS")== 0 ||
                    strcmp(tdata,"MESTRAIN")==0 ||
                    strcmp(tdata,"TOSTRAIN")==0 ||
                    strcmp(tdata,"CONTACT")==0)
            {
                cout<<"FrdReader::readResults()->____reading 3x3 from .frd____"<<endl;
                for(int i=0;i<nv;i++)
                {
                    std::getline(is,val);
                }
                bool isOK = false;
                if(nd!=0)
                {
                    for(int i=0;i<nd;i++)
                    {
                        std::getline(is,val);
                        if(sscanf(val.c_str(),"%d%d%lf%lf%lf%lf%lf%lf",&c,&ni,&cxx,&cyy,&czz,&cxy,&cyz,&cxz)==8)
                        {
                            isOK = true;
                            aResultFile<<ni<<" "<<cxx<<" "<<cyy<<" "<<czz<<" "<<cxy<<" "<<cyz<<" "<<cxz<<endl;
                        }
                        else
                        {
                            //! incomplete data
                            //! consider using throw()
                            isOK = false;
                            break;
                        }
                    }
                }
                else
                {
                    isOK=true;
                }
                aResultFile.close();
                return isOK;
            }
            //! -------------------------------------------
            //! read structural error, heat transfer error
            //! -------------------------------------------
            if(/*strcmp(tdata,"ERROR")== 0 || strcmp(tdata,"HERROR")==0 || */ strcmp(tdata,"PCON")==0)
            {
                cout<<"FrdReader::readResults()->____reading 2x2 from .frd____"<<endl;
                for(int i=0;i<nv;i++)
                {
                    std::getline(is,val);
                }
                bool isOK = false;
                if(nd!=0)
                {
                    for(int i=0;i<nd;i++)
                    {
                        std::getline(is,val);
                        if(sscanf(val.c_str(),"%d%d%lf%lf",&c,&ni,&cxx,&cyy)==4)
                        {
                            isOK = true;
                            aResultFile<<ni<<" "<<cxx<<" "<<cyy<<endl;
                        }
                        else
                        {
                            int t;
                            sscanf(val.c_str(),"%d",&t);
                            if(t==-3)
                            {
                                isOK=true;
                                break;
                            }
                            //! incomplete data
                            //! consider using throw()
                            isOK = false;
                            break;
                        }
                    }
                }
                else
                {
                    isOK=true;
                }
                aResultFile.close();
                return isOK;
            }
        }
        else
        {
            cout<<"FrdReader::readResults()->____exiting"<<endl;

            //! error in reading 3-rd line
            //! consider using throw()
            return false;
        }
        //! -----------------------------------------
        //! this instruction should never be reached
        //! -----------------------------------------
        return true;
    }
    else
    {


        //! error in reading 2-nd line
        //! consider using throw()
        return false;
    }
}
