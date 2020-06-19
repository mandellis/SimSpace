#ifndef SE_EXCEPTION_H
#define SE_EXCEPTION_H

//! C++
#include <stdio.h>
#include <windows.h>
#include <eh.h>

class SE_Exception
{
private:

    unsigned int nSE;

public:

    SE_Exception() {}
    SE_Exception(unsigned int n) : nSE(n) {}
    ~SE_Exception() {}
    unsigned int getSeNumber()
    {
        return nSE;
    }
};

//void trans_func(unsigned int u, EXCEPTION_POINTERS* pExp)
//{
//    printf( "In trans_func.\n" );
//    throw SE_Exception();
//}

#endif // SE_EXCEPTION_H
