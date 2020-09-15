#ifndef CCXTOOLS_H
#define CCXTOOLS_H

#include <vector>
#include <iostream>

#include <QFile>
#include <QVector>
#include <QString>
#include <QMap>
#include <fstream>


namespace CCXTools
{
//! ------------------
//! function: readsta
//! details:
//! ------------------
bool readsta(const QString &path, QMap<double, QVector<int>> &timeinfo)
{
    cout<<"CCXTools::readsta()->_____read STA_______"<<endl;
    QFile fn(path);
    QString nfn = path;
    nfn.chop(4);
    nfn +="_copy.sta";
    fn.copy(nfn);

    FILE *f = fopen(nfn.toStdString().c_str(),"r");
    if(f==NULL) return false;

    char line [256];
    int a1,a2,a4;
    int a1_old,a2_old,a4_old;
    char a3[24],a3_old[24];

    double f1,f2,f3;
    double f1_old,f2_old,f3_old;

    //! ---------------------
    //! jump over the header
    //! ---------------------
    for(int i=0;i<2;i++) fgets(line, sizeof line, f);

    int setnr = 1;

    fgets(line, sizeof line, f);
    if(7!=sscanf(line,"%d%d%s%d%lf%lf%lf",&a1_old,&a2_old,&a3_old,&a4_old,&f1_old,&f2_old,&f3_old))
    {
        fclose(f);
        QFile r(nfn);
        r.remove();
        return false;
    }

    //! ----------------------------------
    //! record the first timeinfo element
    //! ----------------------------------
    QVector<int> sini{setnr,a1_old,a2_old};
    timeinfo.insert(f1_old,sini);

    for(;feof(f)==0;)
    {
        fgets(line, sizeof line, f);
        if(7!=sscanf(line,"%d%d%s%d%lf%lf%lf",&a1,&a2,&a3,&a4,&f1,&f2,&f3))
        {
            fclose(f);
            QFile r(nfn);
            r.remove();
            return false;
        }
        if(a2_old==a2) continue;

        setnr++;
        cout<<"set: "<<setnr<<" step: "<<a1_old<<" substep: "<<a2_old<<" total time: "<<f1_old<<endl;

        //! -------------------
        //! fill the time info
        //! -------------------
        QVector<int> s{setnr,a1,a2};
        timeinfo.insert(f1,s);

        a1_old=a1; a2_old=a2; a4_old=a4;
        f1_old=f1; f2_old=f2; f3_old=f3;
    }

    //! --------------------------
    //! fill the time info - last
    //! --------------------------
    cout<<"set: "<<setnr<<" step: "<<a1<<" substep: "<<a2<<" total time: "<<f1<<endl;
    cout<<"set: "<<setnr<<" step: "<<a1_old<<" substep: "<<a2_old<<" total time: "<<f1_old<<endl;

    QVector<int> s{setnr,a1_old,a2_old};

    if(setnr!=1)
    timeinfo.insert(f1+f3,s);   //check if is correct to add f3 at the end
    else
    {
        timeinfo.remove(f1_old);
        timeinfo.insert(f1+f3,s);   //check if is correct to add f3 at the end
    }
    fclose(f);

    QFile r(nfn);
    r.remove();

    return true;
}

//! ----------------------------
//! function: readcvg
//! details:  to be updated ...
//! ----------------------------
bool readcvg(const QString &path)
{
    //! --------------------------------
    //! check if a lock file is present
    //! --------------------------------
    QString lockFile = path;
    int nchar = lockFile.split("/").last().size();
    lockFile.chop(nchar+1);
    lockFile.append("/cvg.lock");

    QFile file(lockFile);
    if(file.exists()) return false;

    FILE *flock = fopen(lockFile.toStdString().c_str(),"w");
    FILE *f = fopen(path.toStdString().c_str(),"r");
    if(f==NULL)
    {
        fclose(flock);
        QFile r(lockFile);
        r.remove();
        return false;
    }

    char line [256];
    int a1,a2,a3,a4,a5;
    int a1_old,a2_old,a3_old,a4_old,a5_old;

    double f1,f2,f3,f4;
    double f1_old,f2_old,f3_old,f4_old;

    //! ---------------------
    //! jump over the header
    //! ---------------------
    for(int i=0;i<4;i++) fgets(line, sizeof line, f);

    int setnr = 0;

    fgets(line, sizeof line, f);
    sscanf(line,"%d%d%d%d%d%lf%lf%lf%lf",&a1_old,&a2_old,&a3_old,&a4_old,&a5_old,&f1_old,&f2_old,&f3_old,&f4_old);

    for(;feof(f)==0;)
    {
        fgets(line, sizeof line, f);
        sscanf(line,"%d%d%d%d%d%lf%lf%lf%lf",&a1,&a2,&a3,&a4,&a5,&f1,&f2,&f3,&f4);

        if(a2_old!=a2)
        {
            setnr++;
            cout<<"set: "<<setnr<<" step: "<<a1_old<<" substep: "<<a2_old<<endl;
        }
        a1_old=a1; a2_old=a2; a3_old=a3; a4_old=a4; a5_old =a5;
        f1_old=f1; f2_old=f2; f3_old=f3; f4_old=f4;
    }
    fclose(f);
    fclose(flock);

    //! ----
    //! old
    //! ----
    setnr++;
    cout<<"set: "<<setnr<<" step: "<<a1<<" substep: "<<a2<<endl;

    QFile r(lockFile);
    r.remove();

    return true;
}
}

#endif // CCXTOOLS_H
