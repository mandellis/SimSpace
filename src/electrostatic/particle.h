#ifndef PARTICLE_H
#define PARTICLE_H

struct particle
{
public:

    double x[3];    // position
    double v[3];    // velocity
    double mass;    // mass
    double q;       // charge
    double R;       // radius

    //! constructor
    particle(double aX = 0, double aY = 0, double aZ = 0,
             double vx = 0, double vy = 0, double vz = 0,
             double aMass = 1.0, double aCharge = 1.0)
    {
        x[0] = aX; x[1] = aY; x[2] = aZ;
        v[0] = vx; v[1] = vy; v[2] = vz;
        mass =aMass;
        q = aCharge;
    }

    //! copy constructor
    particle (const particle &rhs)
    {
        for(int i=0; i<3; i++)
        {
            x[i] = rhs.x[i];
            v[i] = rhs.v[i];
        }
        mass = rhs.mass;
        q = rhs.q;
        R = rhs.R;
    }

    //! ----------------------
    //! function: setPosition
    //! ----------------------
    void setPosition(double aX, double aY, double aZ)
    {
        x[0] = aX; x[1] = aY; x[2] = aZ;
    }

    //! ------------
    //! setVelocity
    //! ------------
    void setVelocity(double vx, double vy, double vz)
    {
        v[0] = vx; v[1] = vy; v[2] = vz;
    }

    //! -----------------------
    //! setPositionAndVelocity
    //! -----------------------
    void setPositionAndVelocity(double aX, double aY, double aZ,double vx, double vy, double vz)
    {
         x[0] = aX; x[1] = aY; x[2] = aZ;
         v[0] = vx; v[1] = vy; v[2] = vz;
    }
};

#endif // PARTICLE_H
