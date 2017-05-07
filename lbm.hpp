#ifndef LBM_HPP
#define LBM_HPP

#include "grid.h"
#include "settings.hpp"

#include <QThread>

#include <eigen3/Eigen/Dense>
#include <array>
#include <fstream>
#include <iostream>
#include <cmath>
#include <chrono>

#include <Debug.hpp>



class catchGridSignal: public QObject
{
    Q_OBJECT
public slots:
    virtual void compute(unsigned int nIter) = 0;
signals:
    void endCompute();
    void colorUpdated();

};

template <const unsigned int jMax, const  unsigned int iMax>
class LBM: public catchGridSignal
{
public:
    LBM(Eigen::Array<int,Eigen::Dynamic,Eigen::Dynamic>& grid,Eigen::Array<Eigen::Array<unsigned char,3,1>,Eigen::Dynamic,Eigen::Dynamic>& color);
    void compute(unsigned int nIter){Iterate(nIter);}

protected:
    void saveVtk(std::string fileName);
    void Iterate(int nIter);
    void Init();

private:

    Eigen::Array<float,2,1> mRescaler;
    void updateColor();

    Eigen::Array<int,Eigen::Dynamic,Eigen::Dynamic>& mObstacle;
    Eigen::Array<Eigen::Array<unsigned char,3,1>,Eigen::Dynamic,Eigen::Dynamic>& mColorGrid;

    Eigen::Array<float,iMax,jMax> mRho;
    Eigen::Matrix<Eigen::Matrix<float,2,1>,iMax,jMax> mU;

    static constexpr float sNu=0.01, sUin=0.4;
    static constexpr int sQMax = 9;
    static constexpr float sC = 1.73205080756887729352744634150587236694280525381038062805580, sEps = 0.2; //sqrt(3)
    static constexpr float sDx = 2./jMax, sDt=sDx*sEps/sC;
    static constexpr float sEta = 1./(sEps*sEps*sNu/sDt+0.5);
    static constexpr float sCfl = sDt*sNu/(sDx*sDx);

    #ifndef DYNAMIC_ALLOCATION
        Eigen::Array<Eigen::Array<float,sQMax,1>,iMax+2,jMax+2> mGn;

        Eigen::Array<Eigen::Array<float,sQMax,1>,iMax+2,jMax+2> mGnp;

        Eigen::Array<Eigen::Array<float,sQMax,1>,iMax,jMax> mGeq;
    #else
        Eigen::Array<Eigen::Array<float,sQMax,1>,Eigen::Dynamic,Eigen::Dynamic> mGn;

        Eigen::Array<Eigen::Array<float,sQMax,1>,Eigen::Dynamic,Eigen::Dynamic> mGnp;

        Eigen::Array<Eigen::Array<float,sQMax,1>,Eigen::Dynamic,Eigen::Dynamic> mGeq;
    #endif
    Eigen::Matrix<Eigen::Array<float,sQMax,1>,2,1> mSpeeds;

    Eigen::Array<float,sQMax,1> mOmega;

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


template <const unsigned int jMax, const  unsigned int iMax>
LBM<jMax,iMax>::LBM(Eigen::Array<int,Eigen::Dynamic,Eigen::Dynamic>& grid, Eigen::Array<Eigen::Array<unsigned char, 3, 1>, Eigen::Dynamic, Eigen::Dynamic> &color):
    mObstacle(grid),
    mColorGrid(color)
{
    assert(grid.cols() == jMax+2 && grid.rows() == iMax+2);
    mRescaler(0)=1;
    mRescaler(1)=0;
    this->Init();
}

template<const unsigned int jMax,const unsigned int iMax>
void LBM<jMax, iMax>::updateColor()
{
    float normMax = 0;
    float normMin = 1;
    for (unsigned int j = 0; j<jMax; j++)
    {
        for (unsigned int i = 0; i<iMax; i++)
        {
            float norm = mU(i,j).norm();
//            float norm = mRho(i,j);
            normMax = std::max(normMax,norm);
            normMin = std::min(normMin,norm);
            norm -= mRescaler(1);
            if (mRescaler(0)*norm<0.5)
            {
                mColorGrid(i,j)(0) = 255*mRescaler(0)*norm*2;
                mColorGrid(i,j)(1) = 255*norm*mRescaler(0)*2;
                mColorGrid(i,j)(2) = 255;
            }else if(mRescaler(0)*norm<1)
            {
                mColorGrid(i,j)(0) = 255;
                mColorGrid(i,j)(1) = 255*(1-(norm*mRescaler(0)-0.5)*2);
                mColorGrid(i,j)(2) = 255*(1-(norm*mRescaler(0)-0.5)*2);
            }else
            {
                mColorGrid(i,j)(0) = 200;
                mColorGrid(i,j)(1) = 0;
                mColorGrid(i,j)(2) = 0;
            }

        }
    }
    mRescaler(1) = 0.2*normMin+0.8*mRescaler(1);
    mRescaler(0) = 0.2/(normMax-normMin)+0.8*mRescaler(0);
    emit colorUpdated();
}

template<const unsigned int jMax, const unsigned int iMax>
void LBM<jMax, iMax>::saveVtk(std::string fileName)
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
            vtkFile << mU(i,j)(0) << "    " << mU(i,j)(1) << "    " << 0 << std::endl;
        }
    }

    vtkFile.close();
    emit endCompute();
}


template<const unsigned int jMax, const unsigned int iMax>
void LBM<jMax, iMax>::Init()
{
    #ifdef DYNAMIC_ALLOCATION
        mGn.resize(iMax+2,jMax+2);
        mGnp.resize(iMax+2,jMax+2);
        mGeq.resize(iMax,jMax);
    #endif

    float sCTemp = sC;
//               (0,0), (1,0), (2,0), (3,0), (4,0), (5,0), (6,0), (7,0), (8,0)
    mSpeeds(0) <<  0,     sCTemp,    0,    -sCTemp,    0,     sCTemp,   -sCTemp,   -sCTemp,    sCTemp;
//               (0,1), (1,1), (2,1), (3,1), (4,1), (5,1), (6,1), (7,1), (8,1)
    mSpeeds(1) <<  0,     0,     sCTemp,    0,    -sCTemp,    sCTemp,    sCTemp,   -sCTemp,   -sCTemp;
//            0,     1,    2,    3,    4,    5,     6,     7,     8
    mOmega << 4./9, 1./9, 1./9, 1./9, 1./9, 1./36, 1./36, 1./36, 1./36;

//  mU.setZero();
  mRho.setOnes();

  for(unsigned i  = 1; i<iMax+1; i++)
  {
      for (unsigned j = 1; j<jMax+1; j++)
      {
          mU(i-1,j-1).setZero();
          if (mObstacle(i,j)==0)
          {
              mU(i-1,j-1)(0) = 0;
          }else
          {
              mU(i-1,j-1)(0) = (1. - Yc(j-1)* Yc(j-1))*sUin;
          }

          mGn(i,j) = mOmega*mRho(i-1,j-1)*(1+mSpeeds(0)*mU(i-1,j-1)(0)+mSpeeds(1)*mU(i-1,j-1)(1));
          mGnp(i,j) = mGn(i,j);
      }
  }
}

template<const unsigned int jMax, const unsigned int iMax>
void LBM<jMax, iMax>::Iterate(int nIter)
{

    std::chrono::time_point<std::chrono::system_clock> t1 = std::chrono::system_clock::now();
    std::stringstream buf;

    for (int iter = 1; iter<nIter; iter++)
    {
//            time += dt;
        //-- collision
        mRho.setZero();
        for (unsigned int i = 0; i<iMax; i++)
        {
            for (unsigned int j = 0; j<jMax; j++)
            {

                mRho(i,j) += mGn(i+1,j+1).sum();

                if (mObstacle(i,j) == 1)
                {
                    mU(i,j)(0) = 1./mRho(i,j)*sC*(mGn(i+1,j+1)(1)- mGn(i+1,j+1)(3) + mGn(i+1,j+1)(5) - mGn(i+1,j+1)(6) - mGn(i+1,j+1)(7) + mGn(i+1,j+1)(8));
                    mU(i,j)(1) = 1./mRho(i,j)*sC*(mGn(i+1,j+1)(2)- mGn(i+1,j+1)(4) + mGn(i+1,j+1)(5) + mGn(i+1,j+1)(6) - mGn(i+1,j+1)(7) - mGn(i+1,j+1)(8));
                    //- distribution d'equilibre pour Stokes et Navier-Stokes
                    Eigen::Array<float,9,1> dotProd = mU(i,j)(0)*mSpeeds(0)+mU(i,j)(1)*mSpeeds(1);
                    #ifdef STOKES
                        mGeq(i,j) = mOmega*mRho(i,j)*(1+ dotProd);
                    #else
                    #ifdef NAVIERSTOKES
                        mGeq(i,j) = mOmega*mRho(i,j)*(1+dotProd- 0.5*mU(i,j).squaredNorm() + 0.5*dotProd*dotProd);
                    #endif
                    #endif

                    //-- transport des noeuds interieurs (on voit ici pourquoi i et j vont a iMax+2, jMax+2 dans l'allocation de mGnp : c'est pour que le noeuds des bords soient affectés de la meme facon

                    mGnp(i+1,j+1)(0) = (1 - sEta)*mGn(i+1,j+1)(0) + sEta*mGeq(i,j)(0);
                    mGnp(i+2,j+1)(1) = (1 - sEta)*mGn(i+1,j+1)(1) + sEta*mGeq(i,j)(1);
                    mGnp(i+1,j+2)(2) = (1 - sEta)*mGn(i+1,j+1)(2) + sEta*mGeq(i,j)(2);
                    mGnp(i,j+1)(3)   = (1 - sEta)*mGn(i+1,j+1)(3) + sEta*mGeq(i,j)(3);
                    mGnp(i+1,j)(4)   = (1 - sEta)*mGn(i+1,j+1)(4) + sEta*mGeq(i,j)(4);
                    mGnp(i+2,j+2)(5) = (1 - sEta)*mGn(i+1,j+1)(5) + sEta*mGeq(i,j)(5);
                    mGnp(i,j+2)(6)   = (1 - sEta)*mGn(i+1,j+1)(6) + sEta*mGeq(i,j)(6);
                    mGnp(i,j)(7)     = (1 - sEta)*mGn(i+1,j+1)(7) + sEta*mGeq(i,j)(7);
                    mGnp(i+2,j)(8)   = (1 - sEta)*mGn(i+1,j+1)(8) + sEta*mGeq(i,j)(8);

                }else  //obst(i,j) == 0
                {
                    mGnp(i+1,j+1)(0) = mGn(i+1,j+1)(0);
                    mGnp(i+2,j+1)(1) = mGn(i+1,j+1)(3);
                    mGnp(i+1,j+2)(2) = mGn(i+1,j+1)(4);
                    mGnp(i,j+1)(3)   = mGn(i+1,j+1)(1);
                    mGnp(i+1,j)(4)   = mGn(i+1,j+1)(2);
                    mGnp(i+2,j+2)(5) = mGn(i+1,j+1)(7);
                    mGnp(i,j+2)(6)   = mGn(i+1,j+1)(8);
                    mGnp(i,j)(7)     = mGn(i+1,j+1)(5);
                    mGnp(i+2,j)(8)   = mGn(i+1,j+1)(6);
                }


            }
        }

//        //-- CL Sud
//        for (unsigned int i = 2; i<iMax;i++)
//        {
//            //- vitesses rentrantes : BB
//            mGnp(i,1)(2) = mGnp(i,1)(4);
//            mGnp(i,1)(6) = mGnp(i,1)(8);
//            mGnp(i,1)(5) = mGnp(i,1)(7);

//        }
//        //-- CL Nord
//        for (unsigned int i = 2; i<iMax;i++)
//        {
//            //- vitesses rentrantes : BB
//            mGnp(i,jMax)(4) = mGnp(i,jMax)(2);
//            mGnp(i,jMax)(7) = mGnp(i,jMax)(5);
//            mGnp(i,jMax)(8) = mGnp(i,jMax)(6);
//        }
        //-- CL Est : Neumann
        for (unsigned int j =1; j<jMax+1; j++)
        {
            mGnp(iMax,j) = mGnp(iMax-1,j);
        }
        //-- CL Ouest : Dirichlet
        for (unsigned int j =0; j<jMax; j++)
        {

            mGnp(1,j+1) = mOmega*mRho(0,j)*(1+(1-Yc(j)*Yc(j))*sUin*mSpeeds(0));

        }

        mGn = mGnp;

        //~ //--- sorties fichiers
        if(iter%(iterPerCall-1) == 0)
        {
//            buf.str("");
//            buf << "u_" << iter <<".vtk";
//            saveVtk(buf.str());
//            std::cout << "Rescale x" << mRescaler(0) << " -" << mRescaler(1) << std::endl;
//            std::cout << "Total mass: " << mRho.sum() << std::endl;
            updateColor();
        }

    }// fin boucle en temps
    std::chrono::time_point<std::chrono::system_clock> t2 = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = t2-t1;
    //    std::cout << "Time spent: " << elapsed_seconds.count()  << "sec" << std::endl;
    std::cout << "Speed: " << 1./elapsed_seconds.count()  << "fps" << std::endl;
}

#endif // LBM_HPP
