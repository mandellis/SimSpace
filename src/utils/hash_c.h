#ifndef HASH_C_H
#define HASH_C_H

#include <iostream>

//! ----------------------------
//! function: hash_c
//! details:  copied from boost
//! ----------------------------
template <class T>
inline void hash_c(std::size_t &seed, const T &v)
{
    std::hash<T> hasher;
    seed ^= hasher(v)+0x9e3779b9+(seed<<6)+(seed>>2);
}

#endif // HASH_C_H
