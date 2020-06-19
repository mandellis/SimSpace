#ifndef QHISTOGRAMDATA_H
#define QHISTOGRAMDATA_H

#include <QMetaType>

//! ----
//! C++
//! ----
#include <map>
#include <iostream>
#include <fstream>

//! --------------
//! histogram bin
//! --------------
struct hbin
{
    double xmin;
    double xmax;

    //! ------------
    //! constructor
    //! ------------
    hbin(double axmin=0, double axmax = 1):xmin(axmin),xmax(axmax){;}

    //! -----------------
    //! copy constructor
    //! -----------------
    hbin(const hbin &rhs)
    {
        xmin = rhs.xmin;
        xmax = rhs.xmax;
    }

    //! ----------------
    //! interval center
    //! ----------------
    double intervalCenter() const { return xmin+(xmax-xmin)*0.5; }

    //! ----------------
    //! interval length
    //! ----------------
    double intervalLength() const {return (xmax-xmin); }

    //! ------------
    //! operator ==
    //! ------------
    bool operator == (const hbin &rhs) const
    {
        if(xmin==rhs.xmin && xmax==rhs.xmax) return true;
        return false;
    }

    //! -----------
    //! operator <
    //! -----------
    bool operator < (const hbin &rhs) const
    {
        if(this->intervalCenter()<rhs.intervalCenter()) return true;
        return false;
    }    
};
Q_DECLARE_METATYPE(hbin)

//! ----------------------------------------
//! class histogramData is a data container
//! the data are contained as a map
//! ----------------------------------------
struct histogramData
{

    //! -----------
    //! inner data
    //! -----------
    std::map<hbin,double> m_hisData;

    //! ------------
    //! constructor
    //! ------------
    histogramData(const std::map<hbin,double> &hisData = std::map<hbin,double>())
    {
        for(std::map<hbin,double>::const_iterator it = hisData.cbegin(); it!=hisData.cend(); it++)
        {
            const std::pair<hbin,double> &element = *it;
            m_hisData.insert(element);
        }
    }

    //! ---------------------------
    //! function: copy constructor
    //! ---------------------------
    histogramData(const histogramData& rhs)
    {
        for(std::map<hbin,double>::const_iterator it = rhs.m_hisData.cbegin(); it!=rhs.m_hisData.cend(); it++)
        {
            const std::pair<hbin,double> &element = *it;
            m_hisData.insert(element);
        }
    }

    //! ----------------------------
    //! function: insert
    //! details:  insert an element
    //! ----------------------------
    void insert(const std::pair<hbin,double> &element)
    {
        m_hisData.insert(element);
    }

    //! -------------------
    //! function: insert
    //! details:  overload
    //! -------------------
    void insert(const hbin &aBin, double aValue)
    {
        std::pair<hbin,double> element;
        element.first = aBin;
        element.second = aValue;
        m_hisData.insert(element);
    }

    //! ------------------
    //! function: getData
    //! details:
    //! ------------------
    void getData(std::map<hbin,double> &data) const
    {
        for(std::map<hbin,double>::const_iterator it=m_hisData.cbegin(); it!=m_hisData.cend(); it++)
        {
            std::pair<hbin,double> element = *it;
            data.insert(element);
        }
    }

    //! ------------------
    //! function: setData
    //! details:
    //! ------------------
    void setData(const std::map<hbin,double> &hisData)
    {
        for(std::map<hbin,double>::const_iterator it = hisData.cbegin(); it!=hisData.cend(); it++)
        {
            const std::pair<hbin,double> &element = *it;
            m_hisData.insert(element);
        }
    }

    //! ------------------------------------
    //! function: size
    //! details:  return the number of bins
    //! ------------------------------------
    int size()
    {
        return int(m_hisData.size());
        //return sizeof(m_hisData)/sizeof(std::pair<hbin,double>);
        //std::cout<<"____size: "<<sizeof(m_hisData)/sizeof(std::pair<hbin,double>)<<std::endl;
    }

    //! --------------------------------
    //! function: clear
    //! details:  clear the map of data
    //! --------------------------------
    void clear()
    {
        m_hisData.clear();
    }

    //! -----------------------------
    //! function: writehistogramData
    //! details:
    //! -----------------------------
    void writeHistogramData(std::ofstream &s)
    {
        //! ---------------
        //! number of keys
        //! ---------------
        int Nkeys = sizeof(m_hisData)/sizeof(std::pair<hbin,double>);
        s<<Nkeys<<std::endl;

        //! ---------------
        //! write the data
        //! ---------------
        for(std::map<hbin,double>::const_iterator it = m_hisData.cbegin(); it!=m_hisData.cend(); it++)
        {
            const hbin &aBin = (*it).first;
            const double &aVal = (*it).second;
            s<<aBin.xmin<<"\t"<<aBin.xmax<<"\t"<<aVal<<std::endl;
        }
    }

    //! ----------------------------
    //! function: readhistogramData
    //! details:
    //! ----------------------------
    void readHistogramData(std::ifstream &s)
    {
        //! ------------------------
        //! read the number of keys
        //! ------------------------
        int Nkeys;
        s>>Nkeys;

        //! --------------
        //! read the data
        //! --------------
        double xmin,xmax,aVal;
        for(int n=0; n<Nkeys; n++)
        {
            std::string val;
            std::getline(s,val);
            if(3==sscanf(val.c_str(),"%lf%lf%lf",&xmin,&xmax,&aVal))
            {
                std::pair<hbin,double> element;
                element.first.xmin = xmin;
                element.first.xmax = xmax;
                element.second = aVal;
                m_hisData.insert(element);
            }
        }
    }
};
Q_DECLARE_METATYPE(histogramData)

#endif // QHISTROGRAMDATA_H
