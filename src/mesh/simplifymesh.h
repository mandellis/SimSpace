//! -------------------------------------------------------------
//!
//! Mesh Simplification Tutorial
//!
//! (C) by Sven Forstmann in 2014
//!
//! License : MIT
//! http://opensource.org/licenses/MIT
//!
//! https://github.com/sp4cerat/Fast-Quadric-Mesh-Simplification
//!
//! -------------------------------------------------------------

#ifndef SIMPLIFYMESH_H
#define SIMPLIFYMESH_H

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <map>
#include <vector>
#include <string>
#include <math.h>
#include <float.h> //FLT_EPSILON, DBL_EPSILON

#include <QMap>
#include <mesh.h>

#define loopi(start_l,end_l) for ( int i=start_l;i<end_l;++i )
#define loopi(start_l,end_l) for ( int i=start_l;i<end_l;++i )
#define loopj(start_l,end_l) for ( int j=start_l;j<end_l;++j )
#define loopk(start_l,end_l) for ( int k=start_l;k<end_l;++k )

class SymmetricMatrix
{
public:

    // Constructor
    SymmetricMatrix(double c=0) { loopi(0,10) m[i] = c;  }

    SymmetricMatrix(	double m11, double m12, double m13, double m14,
                    double m22, double m23, double m24,
                    double m33, double m34,
                    double m44) {
        m[0] = m11;  m[1] = m12;  m[2] = m13;  m[3] = m14;
        m[4] = m22;  m[5] = m23;  m[6] = m24;
        m[7] = m33;  m[8] = m34;
        m[9] = m44;
    }

    // Make plane
    SymmetricMatrix(double a,double b,double c,double d)
    {
        m[0] = a*a;  m[1] = a*b;  m[2] = a*c;  m[3] = a*d;
        m[4] = b*b;  m[5] = b*c;  m[6] = b*d;
        m[7 ] =c*c; m[8 ] = c*d;
        m[9 ] = d*d;
    }

    double operator[](int c) const { return m[c]; }

    // Determinant
    double det(	int a11, int a12, int a13,
                int a21, int a22, int a23,
                int a31, int a32, int a33)
    {
        double det =  m[a11]*m[a22]*m[a33] + m[a13]*m[a21]*m[a32] + m[a12]*m[a23]*m[a31]
                - m[a13]*m[a22]*m[a31] - m[a11]*m[a23]*m[a32]- m[a12]*m[a21]*m[a33];
        return det;
    }

    const SymmetricMatrix operator+(const SymmetricMatrix& n) const
    {
        return SymmetricMatrix( m[0]+n[0],   m[1]+n[1],   m[2]+n[2],   m[3]+n[3],
                m[4]+n[4],   m[5]+n[5],   m[6]+n[6],
                m[ 7]+n[ 7], m[ 8]+n[8 ],
                m[ 9]+n[9 ]);
    }

    SymmetricMatrix& operator+=(const SymmetricMatrix& n)
    {
        m[0]+=n[0];   m[1]+=n[1];   m[2]+=n[2];   m[3]+=n[3];
        m[4]+=n[4];   m[5]+=n[5];   m[6]+=n[6];   m[7]+=n[7];
        m[8]+=n[8];   m[9]+=n[9];
        return *this;
    }

    double m[10];
};

class SimplifyMesh
{
public:

    //! constructor
    SimplifyMesh();

    // Global Variables & Strctures
    enum Attributes
    {
        NONE,
        NORMAL = 2,
        TEXCOORD = 4,
        COLOR = 8
    };

    struct vector3
    {
        double x, y, z;
    };

    struct vec3f
    {
        double x, y, z;
        inline vec3f( void ) {}

        //inline vec3f operator =( vector3 a )
        // { vec3f b ; b.x = a.x; b.y = a.y; b.z = a.z; return b;}

        inline vec3f( vector3 a )
        { x = a.x; y = a.y; z = a.z; }

        inline vec3f( const double X, const double Y, const double Z )
        { x = X; y = Y; z = Z; }

        inline vec3f operator + ( const vec3f& a ) const
        { return vec3f( x + a.x, y + a.y, z + a.z ); }

        inline vec3f operator += ( const vec3f& a ) const
        { return vec3f( x + a.x, y + a.y, z + a.z ); }

        inline vec3f operator * ( const double a ) const
        { return vec3f( x * a, y * a, z * a ); }

        inline vec3f operator * ( const vec3f a ) const
        { return vec3f( x * a.x, y * a.y, z * a.z ); }

        inline vec3f v3 () const
        { return vec3f( x , y, z ); }

        inline vec3f operator = ( const vector3 a )
        { x=a.x;y=a.y;z=a.z;return *this; }

        inline vec3f operator = ( const vec3f a )
        { x=a.x;y=a.y;z=a.z;return *this; }

        inline vec3f operator / ( const vec3f a ) const
        { return vec3f( x / a.x, y / a.y, z / a.z ); }

        inline vec3f operator - ( const vec3f& a ) const
        { return vec3f( x - a.x, y - a.y, z - a.z ); }

        inline vec3f operator / ( const double a ) const
        { return vec3f( x / a, y / a, z / a ); }

        inline double dot( const vec3f& a ) const
        { return a.x*x + a.y*y + a.z*z; }

        inline vec3f cross( const vec3f& a , const vec3f& b )
        {
            x = a.y * b.z - a.z * b.y;
            y = a.z * b.x - a.x * b.z;
            z = a.x * b.y - a.y * b.x;
            return *this;
        }

        inline double angle( const vec3f& v )
        {
            vec3f a = v , b = *this;
            double dot = v.x*x + v.y*y + v.z*z;
            double len = a.length() * b.length();
            if(len==0)len=0.00001f;
            double input = dot  / len;
            if (input<-1) input=-1;
            if (input>1) input=1;
            return (double) acos ( input );
        }

        inline double angle2( const vec3f& v , const vec3f& w )
        {
            vec3f a = v , b= *this;
            double dot = a.x*b.x + a.y*b.y + a.z*b.z;
            double len = a.length() * b.length();
            if(len==0)len=1;

            vec3f plane; plane.cross( b,w );

            if ( plane.x * a.x + plane.y * a.y + plane.z * a.z > 0 )
                return (double) -acos ( dot  / len );

            return (double) acos ( dot  / len );
        }

        inline vec3f rot_x( double a )
        {
            double yy = cos ( a ) * y + sin ( a ) * z;
            double zz = cos ( a ) * z - sin ( a ) * y;
            y = yy; z = zz;
            return *this;
        }
        inline vec3f rot_y( double a )
        {
            double xx = cos ( -a ) * x + sin ( -a ) * z;
            double zz = cos ( -a ) * z - sin ( -a ) * x;
            x = xx; z = zz;
            return *this;
        }
        inline void clamp( double min, double max )
        {
            if (x<min) x=min;
            if (y<min) y=min;
            if (z<min) z=min;
            if (x>max) x=max;
            if (y>max) y=max;
            if (z>max) z=max;
        }
        inline vec3f rot_z( double a )
        {
            double yy = cos ( a ) * y + sin ( a ) * x;
            double xx = cos ( a ) * x - sin ( a ) * y;
            y = yy; x = xx;
            return *this;
        }
        inline vec3f invert()
        {
            x=-x;y=-y;z=-z;return *this;
        }
        inline vec3f frac()
        {
            return vec3f(
                x-double(int(x)),
                y-double(int(y)),
                z-double(int(z))
                );
        }

        inline vec3f integer()
        {
            return vec3f(
                double(int(x)),
                double(int(y)),
                double(int(z))
                );
        }

        inline double length() const
        {
            return (double)sqrt(x*x + y*y + z*z);
        }

        inline vec3f normalize()
        {
            double square = sqrt(x*x + y*y + z*z);
            x/=square;y/=square;z/=square;
            return *this;
        }
        static vec3f normalize( vec3f a );

        static void random_init();
        static double random_double();
        static vec3f random();

        static int random_number;

        double random_double_01(double a){
            double rnf=a*14.434252+a*364.2343+a*4213.45352+a*2341.43255+a*254341.43535+a*223454341.3523534245+23453.423412;
            int rni=((int)rnf)%100000;
            return double(rni)/(100000.0f-1.0f);
        }

        vec3f random01_fxyz(){
            x=(double)random_double_01(x);
            y=(double)random_double_01(y);
            z=(double)random_double_01(z);
            return *this;
        }
    };

    //! --------------
    //! mesh triangle
    //! --------------
    struct Triangle
    {
        int v[3];
        double err[4];
        int deleted,dirty,attr;
        vec3f n;
        unsigned int tag;
    };

    //! -------
    //! vertex
    //! -------
    struct Vertex
    {
        vec3f p;
        int tstart,tcount;
        SymmetricMatrix q;
        int border;

        //! operator ==
        bool operator==(const Vertex &rhs)
        {
            if(p.x!=rhs.p.x) return false;
            if(p.y!=rhs.p.y) return false;
            if(p.z!=rhs.p.z) return false;
            return true;
        }
    };

    struct Ref
    {
        int tid,tvertex;
    };

    std::vector<Triangle> triangles;
    std::vector<Vertex> vertices;
    std::vector<Ref> refs;
    std::vector<mesh::tolerantPoint> fixedPoints;

    vec3f barycentric(const vec3f &p, const vec3f &a, const vec3f &b, const vec3f &c){
        vec3f v0 = b-a;
        vec3f v1 = c-a;
        vec3f v2 = p-a;
        double d00 = v0.dot(v0);
        double d01 = v0.dot(v1);
        double d11 = v1.dot(v1);
        double d20 = v2.dot(v0);
        double d21 = v2.dot(v1);
        double denom = d00*d11-d01*d01;
        double v = (d11 * d20 - d01 * d21) / denom;
        double w = (d00 * d21 - d01 * d20) / denom;
        double u = 1.0 - v - w;
        return vec3f(u,v,w);
    }

    vec3f interpolate(const vec3f &p, const vec3f &a, const vec3f &b, const vec3f &c, const vec3f attrs[3])
    {
        vec3f bary = barycentric(p,a,b,c);
        vec3f out = vec3f(0,0,0);
        out = out + attrs[0] * bary.x;
        out = out + attrs[1] * bary.y;
        out = out + attrs[2] * bary.z;
        return out;
    }

    double min(double v1, double v2)
    {
        return fmin(v1,v2);
    }

private:

    //! -----------------
    //! Helper functions
    //! -----------------
    double vertex_error(SymmetricMatrix q, double x, double y, double z);
    double calculate_error(int id_v1, int id_v2, vec3f &p_result);
    bool flipped(vec3f p,int i0,int i1,Vertex &v0,Vertex &v1,std::vector<int> &deleted);
    void update_triangles(int i0,Vertex &v,std::vector<int> &deleted,int &deleted_triangles);
    void update_mesh(int iteration);
    void compact_mesh();

    QMap<int,int> triangleTagMap;
    double myTolerance;

public:

    void simplify_mesh(int target_count, double agressiveness=7, bool verbose=false);
    void simplify_mesh_lossless(double threshold, bool verbose=false);
    void load_Stl(const char* fileName, int &NbNodes, int& NbElements);
    void write_Stl(const char* fileName);
    void setFixedPoints(const std::vector<mesh::tolerantPoint> &points);
    void setTolerance(double aTolerance) { myTolerance = aTolerance; }
    double getTolerance() {return myTolerance; }
};

#endif // SIMPLIFYMESH_H
