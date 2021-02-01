#ifndef CCXSOLVERMESSAGE_H
#define CCXSOLVERMESSAGE_H

#include <QMetaType>
#include <QList>
#include <QString>
#include <fstream>
#include <string>

#include <iostream>
using namespace std;

struct CCXSolverMessage
{
    QString myText;

    //! ------------
    //! constructor
    //! ------------
    CCXSolverMessage(QString text=QString("")):myText(text) {;}

    ~CCXSolverMessage(){;}

    //! -----------------
    //! copy constructor
    //! -----------------
    CCXSolverMessage(const CCXSolverMessage &other)
    {
        myText = other.myText;
    }

    //! ----------------
    //! function: write
    //! ----------------
    void write(ofstream &out)
    {
        //out<<std::string("----start solver messages----")<<endl;
        QList<QString> lines = myText.split("\n");
        out<<lines.length();
        for(int i=0; i<lines.length(); i++)
        {
            out<<lines[i].toStdString()<<endl;
        }
        //out<<std::string("----end solver messages----")<<endl;
    }

    //! ---------------
    //! function: read
    //! ---------------
    void read(ifstream &in)
    {
        int nblines;
        in>>nblines;
        std::string val;
        for(int i=0;i<nblines;i++)
        {
            std::getline(in,val);
            myText.append(QString(val.c_str())+"\n");
        }
    }

    //! -----------
    //! operator =
    //! -----------
    CCXSolverMessage operator = (const CCXSolverMessage &rhs)
    {
        myText = rhs.myText;
        return *this;
    }
};

Q_DECLARE_METATYPE(CCXSolverMessage)

#endif // CCXSOLVERMESSAGE_H
