//! ---------------------------------------------------------
//! reimplementation of some classes of the standard library
//! with helpers, and redefinition of some access functions
//! using Qt style
//! ---------------------------------------------------------

#ifndef MYSTDLIB_H
#define MYSTDLIB_H

//! ----
//! C++
//! ----
#include <map>
#include <vector>
#include <algorithm>

//! ------------------------------------------
//! reimplementation of std::map with helpers
//! ------------------------------------------
template<class Q, class T>
class myMap: public std::map<Q,T>
{
public:

    //! ----------------
    //! function: value
    //! ----------------
    T value(Q key, T defReturn = T())
    {
        std::map<Q,T>::iterator it = this->find(key);
        if(it!=this->end()) return it->second;
        return defReturn;
    }

    //! -----------------
    //! function: insert
    //! -----------------
    void setValue(Q key, T value)
    {
        std::pair<Q,T> apair;
        apair.first = key;
        apair.second = value;
        this->insert(apair);
    }

    //! ------------------------
    //! function: key
    //! details:  linear search
    //! ------------------------
    Q key(T value, Q defaultKey = Q())
    {
        for(std::map<Q,T>::const_iterator it = this->begin(); it!= this->cend(); it++)
        {
            if(it->second != value) continue;
            return it->first;
        }
        return defaultKey;
    }

    //! -------------------
    //! function: contains
    //! -------------------
    bool contains(Q key)
    {
        if(this->count(key)) return true;
        return false;
    }
};

//! ---------------------------------------------
//! reimplementation of std::vector with helpers
//! ---------------------------------------------
template <class T>
class myVector: public std::vector<T>
{
public:

    //! ----------------------------------
    //! function: contains
    //! details:  perform a linear search
    //! ----------------------------------
    bool contains(T value, int &indexOf)
    {
        std::vector<T>::const_iterator it = std::find(this->cbegin(), this->cend(), value);
        if(it!=this->cend())
        {
            indexOf = std::distance(this->cbegin(), it);
            return true;
        }
        indexOf = -1;
        return false;
    }
};

#endif // MYSTDLIB_H
