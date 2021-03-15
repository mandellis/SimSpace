#ifndef STLWRAPPERS_H
#define STLWRAPPERS_H

#include <iostream>
using namespace std;

#include <map>
#include <vector>
#include <string>
#include <stdexcept>
#include <regex>
#include <fstream>
#include <experimental/filesystem>
#include <sstream>

#include <QMetaType>
#include <QList>
#include <QVector>
#include <QVariant>
#include <QString>

namespace fs = std::experimental::filesystem;

//! ---------
//! SSVector
//! ---------
template <class T>
class SSVector
{
    std::vector<T> vec_;

public:

    //! ------------------------------------------
    //! constructor, destructor, copy constructor
    //! ------------------------------------------
    SSVector(){ vec_.clear(); }
    SSVector(const std::vector<T> &aVec):vec_(aVec) {}
    SSVector(std::initializer_list<T> l):vec_(l) {}
    SSVector(unsigned int size) { vec_=std::vector<T>(size); }
    virtual ~SSVector() { vec_.clear(); }
    SSVector(const SSVector &aVec) { vec_ = aVec.vec_; }

    using size_type=typename std::vector<T>::size_type;                  // static size_type
    using reference = T&;                                                // reference
    using const_reference = const T&;                                    // const_reference

    //! ----------
    //! operators
    //! ----------
    inline bool operator == (const SSVector &other) { if(vec_ == other.vec_) return true; return false; }
    inline SSVector operator = (const SSVector &other) { vec_ = other.vec_; return *this; }
    reference operator[] (size_type n) { return vec_[n]; }
    const_reference operator[] (size_type n) const {return vec_[n];}

    //! ----------
    //! iterators
    //! ----------
    using iterator= typename std::vector<T>::iterator;                       // iterator
    iterator begin() { return vec_.begin(); }                                // begin
    iterator end() { return vec_.end(); }                                    // end
    using const_iterator = typename std::vector<T>::const_iterator;          // const_iterator
    const_iterator cbegin() const { return vec_.cbegin(); }                  // cbegin
    const_iterator cend() const { return vec_.cend(); }                      // cend

    //! -----------------------------
    //! access/exploration functions
    //! -----------------------------
    void pop_back() { vec_.pop_back(); }
    void push_back(const T &elem) { vec_.push_back(elem); }
    size_type capacity(const T &elem) const { return vec_.capacity(); }
    T* data() { return vec_.data(); }
    bool empty() const { return vec_.empty(); }
    size_type size() const { return vec_.size(); }
    const_reference at(size_type n) const           //! try to throw std::exception::out_of_range
    {
        if(n>=vec_.size()) throw std::out_of_range("out_of_range");
        return vec_.at(n);
    }
    reference at(size_type n)                       //! try to throw std::exception::out_of_range
    {
        try { return vec_.at(n); }
        catch(std::out_of_range /*x*/) { throw std::out_of_range("out_of_range"); }
    }
    void clear() { vec_.clear(); }
    iterator insert(const_iterator where, const T& val) { return vec_.insert(where,val); }
    iterator insert(const_iterator where,size_type count, const T&val) { return vec_.insert(where,count,val); }
    iterator insert(const_iterator where, std::initializer_list<T> l) { return vec_.insert(where,l); }

    //! ------------
    //! PERSISTANCE
    //! ------------
    //! function: read
    SSVector read(std::ifstream &stream)
    {
        std::string line;
        int n;
        std::getline(stream,line);
        std::stringstream ss(line);
        ss>>n;
        for(int i=0; i<n; i++)
        {
            std::getline(stream,line);
            std::stringstream ss(line);
            T val;
            ss>>val;
            vec_.push_back(val);
        }
        return *this;
    }

    //! function: read
    static void read(std::ifstream &stream, SSVector &vec)
    {
        std::string line;
        std::getline(stream,line);
        std::stringstream ss(line);
        int n;
        ss>>n;
        for(int i=0; i<n; i++)
        {
            std::getline(stream,line);
            std::stringstream ss(line);
            T val;
            ss>>val;
            vec.push_back(val);
        }
    }

    //! function: write
    void write(std::ofstream &os) const
    {
        os<<vec_.size()<<std::endl;
        for(int i=0; i<vec_.size(); i++) os<<vec_[i]<<std::endl;
    }

    //! function: write
    static void write(std::ofstream &stream, const SSVector &vec)
    {
        stream<<vec.size()<<std::endl;
        for(int i=0; i<vec.size(); i++) stream<<vec[i]<<std::endl;
    }

    //! ---------
    //! Qt-style
    //! ---------
    int length() const { return int(vec_.size()); }
    int count() const { return int(vec_.size()); }
    bool isEmpty() const { return vec_.empty(); }
    int indexOf(const T &elem) const
    {
        if(vec_.empty()) throw std::out_of_range("out_of_range");
        std::vector<T>::const_iterator it = std::find_if (vec_.begin(),vec_.end(),[elem] (const T& s) { return elem==s; });
        if(it==vec_.cend()) throw std::out_of_range("out_of_range");
        return std::distance(vec_.begin(), it);
    }
    void replace(int i, const T& elem)
    {
        try { vec_.at(i) = elem; }
        catch(std::out_of_range /*x*/) { throw std::out_of_range("out_of_range"); }
    }
    void insert(int pos, const T &val)
    {
        if(vec_.empty())
        {
            if(pos==0) vec_.push_back(val);
            else throw std::out_of_range("out_of_range");
        }
        else
        {
            if(pos == 0)
            {
                std::vector<T> vectemp = vec_;
                vec_.clear();
                vec_.push_back(val);
                for(int i=0; i<vectemp.size(); i++) vec_.push_back(vectemp[i]);
                return;
            }
            else if(pos>0 && pos<=vec_.size()-1)
            {
                size_type NbElements = vec_.size();
                std::vector<T> vectemp = vec_;
                vec_.clear();
                for(int i=0; i<pos; i++) vec_.push_back(vectemp[i]);
                vec_.push_back(val);
                for(int i=pos; i<NbElements; i++) vec_.push_back(vectemp[i]);
                return;
            }
            else throw std::out_of_range("out_of_range");
        }
    }
    void append (const T &val) { vec_.push_back(val); }         //! Qt-style append a value
    void append (const SSVector<T> &other)                      //! Qt-style append a vector to another
    {
        for(SSVector<T>::const_iterator it = other.cbegin(); it!= other.cend(); it++) vec_.push_back(*it);
    }
    void append (const std::vector<T> &other)
    {
        for(std::vector<T>::const_iterator it = other.cbegin(); it!= other.cend(); it++) vec_.push_back(*it);
    }
    void remove (int pos)               //! Qt-style remove
    {
        if(vec_.empty() || pos>=vec_.size()) throw std::out_of_range("out_of_range");           //! try to throw std::out_of_range
        std::vector<T>::const_iterator it = vec_.cbegin();
        for(int i=0; i<pos; i++) it++;
        vec_.erase(it);
    }
    void remove (int pos, int count)    //! Qt-style remove
    {
        if(pos>=vec_.size()) throw std::out_of_range("out_of_range");
        std::vector<T>::const_iterator its = vec_.begin();
        for(int n = 1; n<=pos; n++) its++;
        for(int n = 1; n<=count; n++)
        {
            if(its==vec_.end()) throw std::out_of_range("out_of_range");
            its = vec_.erase(its++);
        }
    }
    void removeFirst() { if(!vec_.empty()) vec_.erase(vec_.cbegin()); }
    void removeLast() {if(!vec_.empty()) vec_.erase(--vec_.cend()); }

    bool contains(const T& elem) const
    {
        std::vector<T>::const_iterator it = std::find_if (vec_.begin(),vec_.end(),[elem] (const T& s) { return elem==s; });
        if(it==vec_.cend()) return false; return true;
    }
    static SSVector<T> fromStdVector(const std::vector<T> &avec)
    {
        SSVector<T> val;
        for(std::vector<T>::const_iterator it = avec.cbegin(); it!=avec.cend(); it++)
            val.push_back(*it); return val;
    }
    static QVector<T> fromSSVector (const SSVector<T> &aVec)
    {
        QVector<T> qtvec;
        for(int i=0; i<aVec.size(); i++) qtvec.push_back(aVec[i]);
        return qtvec;
    }
    iterator erase(const_iterator _Where) { return vec_.erase(_Where); }
    iterator erase(const_iterator _First, const_iterator _Last) { return vec_.erase(_First,_Last); }
};

Q_DECLARE_METATYPE(SSVector<bool>)
Q_DECLARE_METATYPE(SSVector<int>)
Q_DECLARE_METATYPE(SSVector<float>)
Q_DECLARE_METATYPE(SSVector<double>)
Q_DECLARE_METATYPE(SSVector<std::string>)

Q_DECLARE_METATYPE(SSVector<SSVector<bool>>)
Q_DECLARE_METATYPE(SSVector<SSVector<int>>)
Q_DECLARE_METATYPE(SSVector<SSVector<float>>)
Q_DECLARE_METATYPE(SSVector<SSVector<double>>)
Q_DECLARE_METATYPE(SSVector<SSVector<std::string>>)

Q_DECLARE_METATYPE(SSVector<QVariant>)
Q_DECLARE_METATYPE(SSVector<QString>)
Q_DECLARE_METATYPE(SSVector<QList<int>>)
Q_DECLARE_METATYPE(SSVector<QList<double>>)
Q_DECLARE_METATYPE(SSVector<QList<QString>>)

#endif // STLWRAPPERS_H
