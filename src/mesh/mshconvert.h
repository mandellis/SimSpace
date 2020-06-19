#ifndef MSHCONVERT_H
#define MSHCONVERT_H

//! Qt
#include <QString>

//! Eigen
#include <Eigen/Dense>

class mshConvert
{
private:

    QString myMshInputFile;

public:

    //! constructor
    mshConvert(const QString &absoluteFilePath="");    

    //! set file
    void setFile(const QString &absoluteFilePath) { myMshInputFile = absoluteFilePath; }

    //! convert to Eigen format
    bool toEigen(Eigen::MatrixXd &V, Eigen::MatrixXi &T);
};

#endif // MSHCONVERT_H
