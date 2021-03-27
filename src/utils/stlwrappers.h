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
#ifdef SIMSPACE_USES_QT
#include <QMetaType>
#include <QList>
#endif
//namespace fs = std::experimental::filesystem;

namespace ss
{
// begin namespace

//! -------------
//! hash combine
//! -------------
template <class T>
inline void hash_c(std::size_t &seed, const T &v)
{
    std::hash<T> hasher;
    seed ^= hasher(v)+0x9e3779b9+(seed<<6)+(seed>>2);
}
/*
template <typename T>
class Matrix
{
    SSVector<SSVector<T>> inner_;
    unsigned int dimx_, dimy_;

public:

    //! constructor
    Matrix (unsigned int dimx, unsigned int dimy):
        dimx_(dimx),dimy_(dimy)
    {
        inner_.resize(dimx_*dimy_);
        for(unsigned int row=0; row<dimx_; row++)
            for(unsigned int col=0; col<dimy_; col++)
                inner_.push_back(T());
    }

    //! copy constructor
    Matrix (const Matrix &other)
    {
        inner_.resize(other.dimx_*other.dimy_);
        for(unsigned int row=0; row<other.dimx_; row++)
            for(unsigned int col=0; col<other.dimy_; col++)
                inner_.push_back(other.inner_[dimy_*row + col]);
    }

    T& operator()(unsigned int row, unsigned int col)
    {
        if (row >= dimx_ || col>= dimy_) throw std::out_of_range("matrix indices out of range");
        return inner_[dimy_*row + col];
    }

    T operator()(unsigned int row, unsigned int col) const
    {
        if (row >= dimx_ || col>= dimy_) throw std::out_of_range("matrix indices out of range");
        return inner_[dimy_*row + col];
    }

    size_t size() const { return dimx_*dimy_; }
    void clear() { inner_.clear(); }

    //! insert
    void insert(unsigned int row, unsigned int col, const T &value)
    {
        if (row >= dimx_ || col>= dimy_) throw std::out_of_range("matrix indices out of range");
        inner_[row*dimy_+col] = value;
    }

    //! number of rows and columns
    int NbRows() const { return dimx_; }
    int NbCols() const { return dimy_; }

    //! print table - diagnostic
    void print() const
    {
        cout<<endl;
        for(unsigned int row=0; row<dimx_; row++)
        {
            unsigned int col = 0;
            for(; col<dimy_-1; col++) cout<<inner_[row*dimy_+col]<<"\t";
            cout<<inner_[row*dimy_+dimx_-1]<<endl;
        }
        cout<<endl;
    }
};
*/
//! -------------
//! stdException
//! -------------
class stdException
{
public:

    stdException(const std::string &error="", int errorCode=0):error(error),errorCode(errorCode){;}
    stdException(const stdException &other) { error = other.error; errorCode = other.errorCode; }
    ~stdException() {;}

private:

    std::string error;
    int errorCode;
};

//! ------
//! SSMap
//! ------
template<class Q, class T>
class SSMap
{
    std::multimap<Q,T> map_;

public:

    //! ------------------------------------------
    //! constructor, destructor, copy constructor
    //! ------------------------------------------
    SSMap(const std::multimap<Q,T> &aMap = std::multimap<Q,T>()):map_(aMap) {;}
    virtual ~SSMap() { map_.clear(); }
    SSMap(const SSMap &aMap) { map_ = aMap.map_; }

    using size_type=typename std::multimap<Q,T>::size_type;                  // static size_type

    //! ----------
    //! operators
    //! ----------
    inline bool operator == (const SSMap &other) { if(map_ == other.map_) return true; return false; }
    inline SSMap operator = (const SSMap &other) { map_ = other.map_; return *this; }

    //! ----------
    //! iterators
    //! ----------
    using iterator= typename std::multimap<Q,T>::iterator;                   // iterator
    iterator begin() { return map_.begin(); }                                // begin
    iterator end() { return map_.end(); }                                    // end
    using const_iterator = typename std::multimap<Q,T>::const_iterator;      // const_iterator
    const_iterator cbegin() const { return map_.cbegin(); }                  // cbegin
    const_iterator cend() const { return map_.cend(); }                      // cend

    //! --------------------------------------------------------------
    //! insert - insert if not existing: it behaves like C++ std::map
    //! the "if" statement is necessary because the data container is
    //! a std::multimap
    //! --------------------------------------------------------------
    virtual void insert(const std::pair<Q,T> &aPair)
    {
        if(map_.count(aPair.first)==0) map_.insert(aPair);
    }
    virtual void insert(const_iterator hint, const std::pair<Q,T> &aPair)
    {
         if(map_.count(aPair.first)==0)  map_.insert(hint, aPair);
    }

    //! ------
    //! erase
    //! ------
    size_type erase(const Q &key) { return map_.erase(key); }
    iterator erase(const_iterator where) { return map_.erase(where); }
    iterator erase(const_iterator first, const_iterator last) { return map_.erase(first,last); }
    void clear() { map_.clear(); }

    //! -------------------
    //! access/exploration
    //! -------------------
    size_type size() const { return map_.size(); }
    size_type count(const Q &key) const { return map_.count(key); }
    iterator find(const Q &key) { return map_.find(key); }
    const_iterator find(const Q &key) const { return map_.find(key); }
    T& at(const Q &key)
    {
        std::multimap<Q,T>::iterator it = map_.find(key);
        if(it==map_.end()) throw stdException("element not found in map",1);
        else return it->second;
    }
    //const T& at(const Q &key) const { return map_.at(key); }
    T& operator [](const Q &key) { return map_.operator [](key); }
    bool isEmpty() const { return map_.empty(); }
    T last()
    {
        if(!map_.empty()) return (--map_.end())->second;
        else { throw stdException("empty map",1); }
    }
    T first()
    {
        if(!map_.empty()) return (map_.begin())->second;
        else { throw stdException("empty map",1); }
    }
    T lastKey()
    {
        if(!map_.empty()) return (--map_.end())->first;
        else { throw stdException("empty map",1); }
    }
    T firstKey()
    {
        if(!map_.empty()) return (map_.begin())->first;
        else { throw stdException("empty map",1); }
    }

    //! ---------
    //! Qt style
    //! ---------
    const Q key(const T &value, const Q &defaultKey = Q()) const
    {
        for(std::multimap<Q,T>::const_iterator it = map_.cbegin(); it!= map_.cend(); it++)
        {
            if(it->second != value) continue;
            return it->first;
        }
        return defaultKey;
    }
    const T value(const Q &key, const T &defaultValue = T()) const
    {
        std::multimap<Q,T>::const_iterator it = map_.find(key);
        if(it==map_.cend()) return defaultValue;
        return it->second;
    }
    void insert(const Q &key, const T &value)   // Qt style insert and replace
    {
        if(map_.count(key)!=0) map_.erase(key);
        map_.insert(std::make_pair(key,value));
    }
    void insertMulti(const Q &key, const T &value)
    {
        map_.insert(std::make_pair(key,value));
    }
    int remove(const Q &key) { return int(map_.erase(key)); }
    bool contains(const Q &key) const { if(map_.count(key)!=0) return true; return false; }

#ifndef SIMSPACE_USES_QT
    std::vector<Q> keys() const     // all the keys
    {
        std::vector<Q> vecKeys;
        for(std::multimap<Q,T>::iterator it = map_.begin(); it!= map_.end(); it++)
            vecKeys.push_back(it->first);
        return vecKeys;
    }
    std::vector<Q> keys(const T &value) const   // all the keys corresponding to a value
    {
        std::vector<Q> vecKeys;
        for(std::multimap<Q,T>::iterator it = map_.begin(); it!= map_.end(); it++)
            if(it->second == value) vecKeys.push_back(it->first);
        return vecKeys;
    }
    std::vector<T> values() const       // all the values
    {
        std::vector<T> vecValues;
        for(std::multimap<Q,T>::iterator it = map_.begin(); it!= map_.end(); it++)
            vecValues.push_back(it->second);
        return vecValues;
    }
    std::vector<T> values(const Q &key) const   // all the values corresponding to a key
    {
        std::vector<T> vecValues;
        for(std::multimap<Q,T>::iterator it = map_.begin(); it!= map_.end(); it++)
            if(it->first == key) vecValues.push_back(it->second);
        return vecValues;
    }
#else
    QList<T> values() const       // all the values
    {
        QList<T> vecValues;
        for(std::multimap<Q,T>::const_iterator it = map_.cbegin(); it!= map_.cend(); it++)
            vecValues.push_back(it->second);
        return vecValues;
    }
    QList<T> values(const Q &key) const // Qt compatibility
    {
        QList<T> vecValues;
        for(std::multimap<Q,T>::const_iterator it = map_.cbegin(); it!= map_.cend(); it++)
            if(it->first == key) vecValues.push_back(it->second);
        return vecValues;
    }
    QList<Q> keys(const T &value) const // Qt compatibility
    {
        QList<Q> vecKeys;
        for(std::multimap<Q,T>::iterator it = map_.begin(); it!= map_.end(); it++)
            if(it->second == value) vecKeys.push_back(it->first);
        return veckeys;
    }
    QList<Q> keys() const       // Qt compatibility
    {
        QList<Q> vecKeys;
        for(std::multimap<Q,T>::const_iterator it = map_.cbegin(); it!= map_.cend(); it++)
            vecKeys.push_back(it->first);
        return vecKeys;
    }
#endif
};

#ifdef SIMSPACE_USES_QT
Q_DECLARE_METATYPE_TEMPLATE_2ARG(ss::SSMap)
#endif

template <class Q, class T>
class SSMultimap: public SSMap<Q,T>
{

public:

    //! ------------------------------------------
    //! constructor, destructor, copy constructor
    //! ------------------------------------------
    SSMultimap(const std::multimap<Q,T> &aMap = std::multimap<Q,T>()):map_(aMap) {;}
    virtual ~SSMultimap() { ; }
    SSMultimap(const SSMultimap<Q,T> &aMap) { map_ = aMap.map_; }

    //! ---------
    //! Qt style
    //! ---------
    virtual void insert(const Q &key, const T &value) override { map_.insert(std::make_pair(key,value)); }
};

#ifdef SIMSPACE_USES_QT
Q_DECLARE_METATYPE_TEMPLATE_2ARG(ss::SSMultimap)
#endif

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

// end namespace
}

#endif // STLWRAPPERS_H
