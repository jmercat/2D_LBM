#ifndef LBM_HPP
#define LBM_HPP

#include "grid.h"
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include "boost/multi_array.hpp"
#define EIGEN_STACK_ALLOCATION_LIMIT 0

#include <eigen3/Eigen/Dense>
#include <array>
#include <fstream>
#include <iostream>
#include <cmath>
#include <Debug.hpp>

#define NAVIERSTOKES
//#define STOKES


class catchGridSignal: public QObject
{
    Q_OBJECT
public slots:
    virtual void compute(unsigned int nIter) = 0;
signals:
    void endCompute();

};

template <const unsigned int iMax, const  unsigned int jMax>
class LBM: public catchGridSignal
{
public:
    LBM(std::shared_ptr<boost::numeric::ublas::matrix<int> > grid);
    void compute(unsigned int nIter){Iterate(nIter);}
//public slots:
//    void compute(std::shared_ptr<boost::numeric::ublas::matrix<int> > grid);
protected:
    void saveVtk(std::string fileName);
    void Iterate(int nIter);

private:

    std::shared_ptr<boost::numeric::ublas::matrix<int> > mObstacle;

    Eigen::Matrix<float,iMax,jMax> mRho;
    Eigen::Matrix<float,iMax,jMax> mUx;
    Eigen::Matrix<float,iMax,jMax> mUy;


    static constexpr float sNu=0.01, sUin=0.4;
    static constexpr int sQMax = 9;
    static constexpr float sC = 1.73205080756887729352744634150587236694280525381038062805580, sEps = 0.2; //sqrt(3)
    static constexpr float sDx = 2./jMax, sDt=sDx*sEps/sC;
    static constexpr float sEta = 1./(sEps*sEps*sNu/sDt+0.5);
    static constexpr float sCfl = sDt*sNu/(sDx*sDx);

    Eigen::Matrix<float,iMax+2,jMax+2> mGn0;
    Eigen::Matrix<float,iMax+2,jMax+2> mGn1;
    Eigen::Matrix<float,iMax+2,jMax+2> mGn2;
    Eigen::Matrix<float,iMax+2,jMax+2> mGn3;
    Eigen::Matrix<float,iMax+2,jMax+2> mGn4;
    Eigen::Matrix<float,iMax+2,jMax+2> mGn5;
    Eigen::Matrix<float,iMax+2,jMax+2> mGn6;
    Eigen::Matrix<float,iMax+2,jMax+2> mGn7;
    Eigen::Matrix<float,iMax+2,jMax+2> mGn8;

    Eigen::Matrix<float,iMax+2,jMax+2> mGnp0;
    Eigen::Matrix<float,iMax+2,jMax+2> mGnp1;
    Eigen::Matrix<float,iMax+2,jMax+2> mGnp2;
    Eigen::Matrix<float,iMax+2,jMax+2> mGnp3;
    Eigen::Matrix<float,iMax+2,jMax+2> mGnp4;
    Eigen::Matrix<float,iMax+2,jMax+2> mGnp5;
    Eigen::Matrix<float,iMax+2,jMax+2> mGnp6;
    Eigen::Matrix<float,iMax+2,jMax+2> mGnp7;
    Eigen::Matrix<float,iMax+2,jMax+2> mGnp8;

    Eigen::Matrix<float,iMax,jMax> mGeq0;
    Eigen::Matrix<float,iMax,jMax> mGeq1;
    Eigen::Matrix<float,iMax,jMax> mGeq2;
    Eigen::Matrix<float,iMax,jMax> mGeq3;
    Eigen::Matrix<float,iMax,jMax> mGeq4;
    Eigen::Matrix<float,iMax,jMax> mGeq5;
    Eigen::Matrix<float,iMax,jMax> mGeq6;
    Eigen::Matrix<float,iMax,jMax> mGeq7;
    Eigen::Matrix<float,iMax,jMax> mGeq8;

    static constexpr float omega(unsigned int n)
    {
        return n==0 ? 4./9:(n<5 ? 1./9:(n<9 ? 1./36:0));
    };
    static constexpr float speeds(unsigned int i, bool j)
    {
        return (i==1 && j==0) ?  sC:
               (i==2 && j==1) ?  sC:
               (i==3 && j==0) ? -sC:
               (i==4 && j==1) ? -sC:
               (i==5 && j==0) ?  sC:
               (i==5 && j==1) ?  sC:
               (i==6 && j==0) ? -sC:
               (i==6 && j==1) ?  sC:
               (i==7 && j==0) ? -sC:
               (i==7 && j==1) ? -sC:
               (i==8 && j==0) ?  sC:
               (i==8 && j==1) ? -sC:0;
    };

    static constexpr float Xc(unsigned int i)
    {
        ONLYDEBUG(return i<iMax ? (2*i+1)*sDx/2:0;)
        ONLYRELEASE(return (2*i+1)*sDx/2;)
    }

    static constexpr float X(unsigned int i)
    {
        ONLYDEBUG(return i<iMax ? i*sDx:0;)
        ONLYRELEASE(return i*sDx;)
    }

    static constexpr float Yc(unsigned int j)
    {
        ONLYDEBUG(return j<jMax ? (2*j+1)*sDx/2 -1:0;)
        ONLYRELEASE(return (2*j+1)*sDx/2 -1;)
    }

    static constexpr float Y(unsigned int j)
    {
        ONLYDEBUG(return j<jMax ? j*sDx-1:0;)
        ONLYRELEASE(return j*sDx;)
    }
};

template<const unsigned int iMax,const unsigned int jMax>
LBM<iMax,jMax>::LBM(std::shared_ptr<boost::numeric::ublas::matrix<int> > grid):
    mObstacle(grid)
{
    assert(grid->size1() == iMax+2 && grid->size2() == jMax+2);
}

template<const unsigned int iMax, const unsigned int jMax>
void LBM<iMax, jMax>::saveVtk(std::string fileName)
{
    std::ofstream vtkFile;
    std::cout << "save " << fileName << std::endl;
    vtkFile.open(fileName);
    vtkFile << "# vtk DataFile Version 2.0" << std::endl;
    vtkFile << "champ de vitesse" << std::endl;
    vtkFile << "ASCII" << std::endl;
    vtkFile << "DATASET RECTILINEAR_GRID" << std::endl;
    vtkFile << "DIMENSIONS    " << iMax+1 << "    "<< jMax+1 << "    " << 1 << std::endl;
    vtkFile << "X_COORDINATES    " << iMax+1 <<"    double" << std::endl;
    for (unsigned int i = 0; i<iMax+1; i++)
    {
        vtkFile << X(i) << std::endl;
    }
    vtkFile << "Y_COORDINATES    " << jMax+1 <<"    double" << std::endl;
    for (unsigned int j = 0; j<jMax+1; j++)
    {
        vtkFile << Y(j) << std::endl;
    }
    vtkFile << "Z_COORDINATES    " << 1 << "    double" << std::endl;
    vtkFile << 0 << std::endl;

    vtkFile << "CELL_DATA    " << iMax*jMax << std::endl;
    vtkFile << "SCALARS rho double" << std::endl;
    vtkFile << "LOOKUP_TABLE default" << std::endl;
    for (unsigned int j = 0; j<jMax; j++)
    {
        for (unsigned int i = 0; i<iMax; i++)
        {
            vtkFile << mRho(i,j) << std::endl;
        }
    }

    vtkFile << "VECTORS u double" << std::endl;
    for (unsigned int j = 0; j<jMax; j++)
    {
        for (unsigned int i = 0; i<iMax; i++)
        {
            vtkFile << mUx(i,j) << "    " << mUy(i,j) << "    " << 0 << std::endl;
        }
    }

    vtkFile.close();
    emit endCompute();
}


template<const unsigned int iMax, const unsigned int jMax>
void LBM<iMax, jMax>::Iterate(int nIter)
{

    mUy.setZero();

    for(unsigned i  = 1; i<iMax+1; i++)
    {
        for (unsigned j = 1; j<jMax+1; j++)
        {
            mRho(i-1,j-1)=1;
            if ((*mObstacle)(i,j)==0)
            {
                mUx(i-1,j-1) = 0;
            }else
            {
                mUx(i-1,j-1) = (1. - Yc(j-1)* Yc(j-1))*sUin;
            }

            mGn0(i,j) = omega(0)*mRho(i-1,j-1)*(1+mUx(i-1,j-1)*speeds(0,0)+mUy(i-1,j-1)*speeds(0,1));
            mGnp0(i,j) = mGn0(i,j);
            mGn1(i,j) = omega(1)*mRho(i-1,j-1)*(1+mUx(i-1,j-1)*speeds(1,0)+mUy(i-1,j-1)*speeds(1,1));
            mGnp1(i,j) = mGn1(i,j);
            mGn2(i,j) = omega(2)*mRho(i-1,j-1)*(1+mUx(i-1,j-1)*speeds(2,0)+mUy(i-1,j-1)*speeds(2,1));
            mGnp2(i,j) = mGn2(i,j);
            mGn3(i,j) = omega(3)*mRho(i-1,j-1)*(1+mUx(i-1,j-1)*speeds(3,0)+mUy(i-1,j-1)*speeds(3,1));
            mGnp3(i,j) = mGn3(i,j);
            mGn4(i,j) = omega(4)*mRho(i-1,j-1)*(1+mUx(i-1,j-1)*speeds(4,0)+mUy(i-1,j-1)*speeds(4,1));
            mGnp4(i,j) = mGn4(i,j);
            mGn5(i,j) = omega(5)*mRho(i-1,j-1)*(1+mUx(i-1,j-1)*speeds(5,0)+mUy(i-1,j-1)*speeds(5,1));
            mGnp5(i,j) = mGn5(i,j);
            mGn6(i,j) = omega(6)*mRho(i-1,j-1)*(1+mUx(i-1,j-1)*speeds(6,0)+mUy(i-1,j-1)*speeds(6,1));
            mGnp6(i,j) = mGn6(i,j);
            mGn7(i,j) = omega(7)*mRho(i-1,j-1)*(1+mUx(i-1,j-1)*speeds(7,0)+mUy(i-1,j-1)*speeds(7,1));
            mGnp7(i,j) = mGn7(i,j);
            mGn8(i,j) = omega(8)*mRho(i-1,j-1)*(1+mUx(i-1,j-1)*speeds(8,0)+mUy(i-1,j-1)*speeds(8,1));
            mGnp8(i,j) = mGn8(i,j);
        }
    }

    std::stringstream buf;
    std::cout << "iteration " << 0 << std::endl;
    buf.str("");
    buf << "u_" << 0 <<".vtk";
    saveVtk(buf.str());

    for (int iter = 1; iter<nIter; iter++)
    {
//            time += dt;
        //-- collision
        for (unsigned int i = 0; i<iMax; i++)
        {
            for (unsigned int j = 0; j<jMax; j++)
            {
                mRho(i,j) = 0;

                mRho(i,j) += mGn0(i+1,j+1);
                mRho(i,j) += mGn1(i+1,j+1);
                mRho(i,j) += mGn2(i+1,j+1);
                mRho(i,j) += mGn3(i+1,j+1);
                mRho(i,j) += mGn4(i+1,j+1);
                mRho(i,j) += mGn5(i+1,j+1);
                mRho(i,j) += mGn6(i+1,j+1);
                mRho(i,j) += mGn7(i+1,j+1);
                mRho(i,j) += mGn8(i+1,j+1);

                //-- pour que la visualisation soit plus claire, on met mUx et uy a 0 dans le solide

                if ((*mObstacle)(i,j) == 1)
                {
                    mUx(i,j) = 1./mRho(i,j)*sC*(mGn1(i+1,j+1)- mGn3(i+1,j+1) + mGn5(i+1,j+1) - mGn6(i+1,j+1) - mGn7(i+1,j+1) + mGn8(i+1,j+1));
                    mUy(i,j) = 1./mRho(i,j)*sC*(mGn2(i+1,j+1)- mGn4(i+1,j+1) + mGn5(i+1,j+1) + mGn6(i+1,j+1) - mGn7(i+1,j+1) - mGn8(i+1,j+1));
                    //- distribution d'equilibre pour Stokes et Navier-Stokes
                    #ifdef STOKES
                        mGeq0(i,j) = omega(0)*mRho(i,j)*(1+ mUx(i,j)*speeds(0,0)+mUy(i,j)*speeds(0,1));
                        mGeq1(i,j) = omega(1)*mRho(i,j)*(1+ mUx(i,j)*speeds(1,0)+mUy(i,j)*speeds(1,1));
                        mGeq2(i,j) = omega(2)*mRho(i,j)*(1+ mUx(i,j)*speeds(2,0)+mUy(i,j)*speeds(2,1));
                        mGeq3(i,j) = omega(3)*mRho(i,j)*(1+ mUx(i,j)*speeds(3,0)+mUy(i,j)*speeds(3,1));
                        mGeq4(i,j) = omega(4)*mRho(i,j)*(1+ mUx(i,j)*speeds(4,0)+mUy(i,j)*speeds(4,1));
                        mGeq5(i,j) = omega(5)*mRho(i,j)*(1+ mUx(i,j)*speeds(5,0)+mUy(i,j)*speeds(5,1));
                        mGeq6(i,j) = omega(6)*mRho(i,j)*(1+ mUx(i,j)*speeds(6,0)+mUy(i,j)*speeds(6,1));
                        mGeq7(i,j) = omega(7)*mRho(i,j)*(1+ mUx(i,j)*speeds(7,0)+mUy(i,j)*speeds(7,1));
                        mGeq8(i,j) = omega(8)*mRho(i,j)*(1+ mUx(i,j)*speeds(8,0)+mUy(i,j)*speeds(8,1));
                    #else
                    #ifdef NAVIERSTOKES
                        mGeq0(i,j) = omega(0)*mRho(i,j)*(1+mUx(i,j)*speeds(0,0)+ mUy(i,j)*speeds(0,1)- 0.5*(mUx(i,j)*mUx(i,j)+mUy(i,j)*mUy(i,j)) + 0.5*(mUx(i,j)*speeds(0,0)+ mUy(i,j)*speeds(0,1))*(mUx(i,j)*speeds(0,0)+ mUy(i,j)*speeds(0,1)));
                        mGeq1(i,j) = omega(1)*mRho(i,j)*(1+mUx(i,j)*speeds(1,0)+ mUy(i,j)*speeds(1,1)- 0.5*(mUx(i,j)*mUx(i,j)+mUy(i,j)*mUy(i,j)) + 0.5*(mUx(i,j)*speeds(1,0)+ mUy(i,j)*speeds(1,1))*(mUx(i,j)*speeds(1,0)+ mUy(i,j)*speeds(1,1)));
                        mGeq2(i,j) = omega(2)*mRho(i,j)*(1+mUx(i,j)*speeds(2,0)+ mUy(i,j)*speeds(2,1)- 0.5*(mUx(i,j)*mUx(i,j)+mUy(i,j)*mUy(i,j)) + 0.5*(mUx(i,j)*speeds(2,0)+ mUy(i,j)*speeds(2,1))*(mUx(i,j)*speeds(2,0)+ mUy(i,j)*speeds(2,1)));
                        mGeq3(i,j) = omega(3)*mRho(i,j)*(1+mUx(i,j)*speeds(3,0)+ mUy(i,j)*speeds(3,1)- 0.5*(mUx(i,j)*mUx(i,j)+mUy(i,j)*mUy(i,j)) + 0.5*(mUx(i,j)*speeds(3,0)+ mUy(i,j)*speeds(3,1))*(mUx(i,j)*speeds(3,0)+ mUy(i,j)*speeds(3,1)));
                        mGeq4(i,j) = omega(4)*mRho(i,j)*(1+mUx(i,j)*speeds(4,0)+ mUy(i,j)*speeds(4,1)- 0.5*(mUx(i,j)*mUx(i,j)+mUy(i,j)*mUy(i,j)) + 0.5*(mUx(i,j)*speeds(4,0)+ mUy(i,j)*speeds(4,1))*(mUx(i,j)*speeds(4,0)+ mUy(i,j)*speeds(4,1)));
                        mGeq5(i,j) = omega(5)*mRho(i,j)*(1+mUx(i,j)*speeds(5,0)+ mUy(i,j)*speeds(5,1)- 0.5*(mUx(i,j)*mUx(i,j)+mUy(i,j)*mUy(i,j)) + 0.5*(mUx(i,j)*speeds(5,0)+ mUy(i,j)*speeds(5,1))*(mUx(i,j)*speeds(5,0)+ mUy(i,j)*speeds(5,1)));
                        mGeq6(i,j) = omega(6)*mRho(i,j)*(1+mUx(i,j)*speeds(6,0)+ mUy(i,j)*speeds(6,1)- 0.5*(mUx(i,j)*mUx(i,j)+mUy(i,j)*mUy(i,j)) + 0.5*(mUx(i,j)*speeds(6,0)+ mUy(i,j)*speeds(6,1))*(mUx(i,j)*speeds(6,0)+ mUy(i,j)*speeds(6,1)));
                        mGeq7(i,j) = omega(7)*mRho(i,j)*(1+mUx(i,j)*speeds(7,0)+ mUy(i,j)*speeds(7,1)- 0.5*(mUx(i,j)*mUx(i,j)+mUy(i,j)*mUy(i,j)) + 0.5*(mUx(i,j)*speeds(7,0)+ mUy(i,j)*speeds(7,1))*(mUx(i,j)*speeds(7,0)+ mUy(i,j)*speeds(7,1)));
                        mGeq8(i,j) = omega(8)*mRho(i,j)*(1+mUx(i,j)*speeds(8,0)+ mUy(i,j)*speeds(8,1)- 0.5*(mUx(i,j)*mUx(i,j)+mUy(i,j)*mUy(i,j)) + 0.5*(mUx(i,j)*speeds(8,0)+ mUy(i,j)*speeds(8,1))*(mUx(i,j)*speeds(8,0)+ mUy(i,j)*speeds(8,1)));
                    #endif
                    #endif

                    //-- transport des noeuds interieurs (on voit ici pourquoi i et j vont a iMax+2, jMax+2 dans l'allocation de mGnp1 : c'est pour que le noeuds des bords soient affectés de la meme facon
                    mGnp0(i+1,j+1) = (1 - sEta)*mGn0(i+1,j+1) + sEta*mGeq0(i,j);
                    mGnp1(i+2,j+1) = (1 - sEta)*mGn1(i+1,j+1) + sEta*mGeq1(i,j);
                    mGnp2(i+1,j+2) = (1 - sEta)*mGn2(i+1,j+1) + sEta*mGeq2(i,j);
                    mGnp3(i,j+1)   = (1 - sEta)*mGn3(i+1,j+1) + sEta*mGeq3(i,j);
                    mGnp4(i+1,j)   = (1 - sEta)*mGn4(i+1,j+1) + sEta*mGeq4(i,j);
                    mGnp5(i+2,j+2) = (1 - sEta)*mGn5(i+1,j+1) + sEta*mGeq5(i,j);
                    mGnp6(i,j+2)   = (1 - sEta)*mGn6(i+1,j+1) + sEta*mGeq6(i,j);
                    mGnp7(i,j)     = (1 - sEta)*mGn7(i+1,j+1) + sEta*mGeq7(i,j);
                    mGnp8(i+2,j)   = (1 - sEta)*mGn8(i+1,j+1) + sEta*mGeq8(i,j);

                }else  //obst(i,j) == 0
                {
                    mGnp0(i+1,j+1) = mGn0(i+1,j+1);
                    mGnp1(i+2,j+1) = mGn3(i+1,j+1);
                    mGnp2(i+1,j+2) = mGn4(i+1,j+1);
                    mGnp3(i,j+1)   = mGn1(i+1,j+1);
                    mGnp4(i+1,j)   = mGn2(i+1,j+1);
                    mGnp5(i+2,j+2) = mGn7(i+1,j+1);
                    mGnp6(i,j+2)   = mGn8(i+1,j+1);
                    mGnp7(i,j)     = mGn5(i+1,j+1);
                    mGnp8(i+2,j)   = mGn6(i+1,j+1);
                }


            }
        }

        //-- CL Sud
        for (unsigned int i = 2; i<iMax;i++)
        {
            //- vitesses rentrantes : BB
            mGnp2(i,1) = mGnp4(i,1);
            mGnp6(i,1) = mGnp8(i,1);
            mGnp5(i,1) = mGnp7(i,1);

        }
        //-- CL Nord
        for (unsigned int i = 2; i<iMax;i++)
        {
            //- vitesses rentrantes : BB
            mGnp4(i,jMax) = mGnp2(i,jMax);
            mGnp7(i,jMax) = mGnp5(i,jMax);
            mGnp8(i,jMax) = mGnp6(i,jMax);
        }
        //-- CL Est : Neumann
        for (unsigned int j =1; j<jMax+1; j++)
        {
            mGnp0(iMax,j) = mGnp0(iMax-1,j);
            mGnp1(iMax,j) = mGnp1(iMax-1,j);
            mGnp2(iMax,j) = mGnp2(iMax-1,j);
            mGnp3(iMax,j) = mGnp3(iMax-1,j);
            mGnp4(iMax,j) = mGnp4(iMax-1,j);
            mGnp5(iMax,j) = mGnp5(iMax-1,j);
            mGnp6(iMax,j) = mGnp6(iMax-1,j);
            mGnp7(iMax,j) = mGnp7(iMax-1,j);
            mGnp8(iMax,j) = mGnp8(iMax-1,j);
        }
        //-- CL Ouest : Dirichlet
        for (unsigned int j =0; j<jMax; j++)
        {   
            mGnp0(1,j+1) = omega(0)*mRho(0,j)*(1+(1-Yc(j)*Yc(j))*sUin*speeds(0,0));
            mGnp1(1,j+1) = omega(1)*mRho(0,j)*(1+(1-Yc(j)*Yc(j))*sUin*speeds(1,0));
            mGnp2(1,j+1) = omega(2)*mRho(0,j)*(1+(1-Yc(j)*Yc(j))*sUin*speeds(2,0));
            mGnp3(1,j+1) = omega(3)*mRho(0,j)*(1+(1-Yc(j)*Yc(j))*sUin*speeds(3,0));
            mGnp4(1,j+1) = omega(4)*mRho(0,j)*(1+(1-Yc(j)*Yc(j))*sUin*speeds(4,0));
            mGnp5(1,j+1) = omega(5)*mRho(0,j)*(1+(1-Yc(j)*Yc(j))*sUin*speeds(5,0));
            mGnp6(1,j+1) = omega(6)*mRho(0,j)*(1+(1-Yc(j)*Yc(j))*sUin*speeds(6,0));
            mGnp7(1,j+1) = omega(7)*mRho(0,j)*(1+(1-Yc(j)*Yc(j))*sUin*speeds(7,0));
            mGnp8(1,j+1) = omega(8)*mRho(0,j)*(1+(1-Yc(j)*Yc(j))*sUin*speeds(8,0));
        }

        mGn0.swap(mGnp0);
        mGn1.swap(mGnp1);
        mGn2.swap(mGnp2);
        mGn3.swap(mGnp3);
        mGn4.swap(mGnp4);
        mGn5.swap(mGnp5);
        mGn6.swap(mGnp6);
        mGn7.swap(mGnp7);
        mGn8.swap(mGnp8);

        //~ //--- sorties fichiers
        if(iter%50 == 0)
        {
            std::cout << "iteration " << iter << std::endl;
            buf.str("");
            buf << "u_" << iter <<".vtk";
            saveVtk(buf.str());
        }

    }// fin boucle en temps
}




#endif // LBM_HPP
