#ifndef LBM_HPP
#define LBM_HPP

#include "grid.h"
#include "settings.hpp"

#include <QThread>

#include <Eigen/Dense>
#include <array>
#include <fstream>
#include <iostream>
#include <cmath>
#include <chrono>

#include "Debug.hpp"


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
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    LBM(Eigen::Array<int,Eigen::Dynamic,Eigen::Dynamic>& grid,Eigen::Array<Eigen::Array3f,Eigen::Dynamic,Eigen::Dynamic>& result);
    void compute(unsigned int nIter){Iterate(nIter);}

protected:
    void saveVtk(std::string fileName) const;
    void Iterate(int nIter);
    void Init();

private:

    Eigen::Array2f mRescalerU;
    Eigen::Array2f mRescalerR;
    void updateColor();

    Eigen::Array<int,Eigen::Dynamic,Eigen::Dynamic>& mObstacle;
    Eigen::Array<Eigen::Array3f,Eigen::Dynamic,Eigen::Dynamic>& mResultGrid;

    Eigen::Array<float,iMax,jMax> mRho;
    Eigen::Matrix<Eigen::Vector2f,iMax,jMax> mU;

    //Viscosity and Boltzman entry speed
    static constexpr float sNu=0.01, sUin=0.4;
    //Lattice Boltzman D2Q9 uses 9 speeds
    static constexpr int sQMax = 9;
    //Just because I never know when to stop...
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

// unary functors
template <typename Derived>
struct sumCoef
{
  typedef typename Eigen::DenseBase<Derived>::Scalar result_type;
  result_type operator()(const Eigen::DenseBase<Derived>& m) const { return m.sum(); }
};

struct computeSpeed
{
  typedef Eigen::Vector2f result_type;
  Eigen::Vector2f operator()(const Eigen::Array<float,9,1>& gn) const
  {
      float x = (gn(1)-gn(3)+gn(5)-gn(6)-gn(7)+gn(8));
      float y = (gn(2)-gn(4)+gn(5)+gn(6)-gn(7)-gn(8));
      return Eigen::Vector2f(x,y);
  }
};

struct computeEquilibriumDistribution
{
    typedef Eigen::Array<float,9,1> result_type;
    result_type operator()(const Eigen::Vector2f& U) const
    {

        Eigen::Array<float,9,1> dotProd;


        float C = 1.73205080756887729352744634150587236694280525381038062805580;
        float UC0 = U(0)*C;
        float UC1 = U(1)*C;
        float UC0pUC1 = UC0+UC1;
        float UC0mUC1 = UC0-UC1;
        float USNo2m1 = 0.5*U.squaredNorm()-1;

        dotProd(0) = 0;
        dotProd(1) = UC0;
        dotProd(2) = UC1;
        dotProd(3) = -UC0;
        dotProd(4) = -UC1;
        dotProd(5) = UC0pUC1;
        dotProd(6) = -UC0mUC1;
        dotProd(7) = -UC0pUC1;
        dotProd(8) = UC0mUC1;


        #ifdef STOKES
            return (Eigen::VectorXf::Ones(9)+ dotProd);
        #else
        #ifdef NAVIERSTOKES
            return (dotProd- USNo2m1 + 0.5*dotProd*dotProd);
        #endif
        #endif
    }
};


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

          mGn(i,j)(0) = mOmega(0)*mRho(i-1,j-1);
          mGn(i,j)(1) = mOmega(1)*mRho(i-1,j-1)*(1+sC*mU(i-1,j-1)(0));
          mGn(i,j)(2) = mOmega(2)*mRho(i-1,j-1)*(1+sC*mU(i-1,j-1)(1));
          mGn(i,j)(3) = mOmega(3)*mRho(i-1,j-1)*(1-sC*mU(i-1,j-1)(0));
          mGn(i,j)(4) = mOmega(4)*mRho(i-1,j-1)*(1-sC*mU(i-1,j-1)(1));
          mGn(i,j)(5) = mOmega(5)*mRho(i-1,j-1)*(1+sC*(mU(i-1,j-1)(0)+mU(i-1,j-1)(1)));
          mGn(i,j)(6) = mOmega(6)*mRho(i-1,j-1)*(1-sC*(mU(i-1,j-1)(0)-mU(i-1,j-1)(1)));
          mGn(i,j)(7) = mOmega(7)*mRho(i-1,j-1)*(1-sC*(mU(i-1,j-1)(0)+mU(i-1,j-1)(1)));
          mGn(i,j)(8) = mOmega(8)*mRho(i-1,j-1)*(1+sC*(mU(i-1,j-1)(0)-mU(i-1,j-1)(1)));
          mGnp(i,j) = mGn(i,j);

      }
  }
}



template<const unsigned int jMax, const unsigned int iMax>
void LBM<jMax, iMax>::Iterate(int nIter)
{

    std::chrono::time_point<std::chrono::system_clock> t1 = std::chrono::system_clock::now();
    std::stringstream buf;

    for (int iter = 1; iter<nIter; iter++) // time loop
    {
        //const auto& centerGn = mGn.block<iMax,jMax>(1,1);
        const auto& centerGn = mGn.block(1, 1, iMax, jMax);

        
        //integrate density to get pressure:
        mRho = centerGn.unaryExpr(sumCoef<Eigen::Array<float,sQMax,1> >());

        //Compute macroscopic speeds
        mU = centerGn.unaryExpr(computeSpeed());
        // normalize
        for (size_t i = 0, sizeRho = mRho.size(); i<sizeRho; i++)
        {mU(i) *= sC/mRho(i);}

        //Compute equilibrium distribution
        mGeq = mU.unaryExpr(computeEquilibriumDistribution());
        // normalize
        for (size_t i = 0, sizeRho = mRho.size(); i<sizeRho; i++)
        {mGeq(i) *= mOmega*mRho(i);}



        //-- collision
        for (unsigned int j = 0; j<jMax; j++)
        {
            for (unsigned int i = 0; i<iMax; i++)
            {

                if (mObstacle(i,j) == 1)
                {
                    Eigen::Matrix<float,sQMax,1> temp;
                    temp = (1 - sEta)*mGn(i+1,j+1) + sEta*mGeq(i,j);
                    //-- interior nodes transport (this is why Gnp is allocated to size (iMax+2, jMax+2)
                    //-- so border nodes are treated the same way than the others
                    mGnp(i+1,j+1)(0) = temp(0);
                    mGnp(i+2,j+1)(1) = temp(1);
                    mGnp(i+1,j+2)(2) = temp(2);
                    mGnp(i,j+1)  (3) = temp(3);
                    mGnp(i+1,j)  (4) = temp(4);
                    mGnp(i+2,j+2)(5) = temp(5);
                    mGnp(i,j+2)  (6) = temp(6);
                    mGnp(i,j)    (7) = temp(7);
                    mGnp(i+2,j)  (8) = temp(8);


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

        //-- CL Est : Neumann
        mGnp.row(iMax) = mGnp.row(iMax-1);

        //-- CL Ouest : Dirichlet
        for (unsigned int j =0; j<jMax; j++)
        {
            mGnp(1,j+1) = mOmega*mRho(0,j)*(1+(1-Yc(j)*Yc(j))*sUin*mSpeeds(0));
        }

        mGn.swap(mGnp);

        //~ //--- sorties fichiers
        if(iter%(iterPerCall-1) == 0)
        {
//            buf.str("");
//            buf << "u_" << iter <<".vtk";
//            saveVtk(buf.str());
//            std::cout << "Rescale x" << mRescalerU(0) << " -" << mRescalerU(1) << '\n';
//            std::cout << "Total mass: " << mRho.sum() << '\n';
            if (std::isnan(mRho.sum()))
                this->Init();
            updateColor();
        }

    }// end of time loop

    std::chrono::time_point<std::chrono::system_clock> t2 = std::chrono::system_clock::now();
    std::chrono::duration<float> elapsed_seconds = t2-t1;
    //    std::cout << "Time spent: " << elapsed_seconds.count()  << "sec" << '\n';
    std::cout << "Speed: " << 1./elapsed_seconds.count()  << "fps" << '\n';
}


template <const unsigned int jMax, const  unsigned int iMax>
LBM<jMax,iMax>::LBM(Eigen::Array<int,Eigen::Dynamic,Eigen::Dynamic>& grid, Eigen::Array<Eigen::Array3f, Eigen::Dynamic, Eigen::Dynamic> &result):
    mObstacle(grid),
    mResultGrid(result)
{
    assert(grid.cols() == jMax+2 && grid.rows() == iMax+2);
    mRescalerU(0)=1;
    mRescalerU(1)=0;
    mRescalerR(0)=1;
    mRescalerR(1)=0;
    this->Init();
}

template<const unsigned int jMax,const unsigned int iMax>
void LBM<jMax, iMax>::updateColor()
{
    for (unsigned int j = 0; j<jMax; j++)
    {
        for (unsigned int i = 0; i<iMax; i++)
        {
            if(i == 0 || j == 0 || mObstacle(i,j) != 0)// || i == iMax-1 || j == jMax-1)
                {
                mResultGrid(i,j)(0) = mU(i,j)(0);
                mResultGrid(i,j)(1) = mU(i,j)(1);
                mResultGrid(i,j)(2) = mRho(i,j);
            }else if(mObstacle(i,j)!=0)
            {
                mResultGrid(i,j).setZero();
                mResultGrid(i,j)(2) = mRho.mean();
            }

        }
    }
    emit colorUpdated();
}

template<const unsigned int jMax, const unsigned int iMax>
void LBM<jMax, iMax>::saveVtk(std::string fileName) const
{
    std::ofstream vtkFile;
    std::cout << "save " << fileName << '\n';
    vtkFile.open(fileName);
    vtkFile << "# vtk DataFile Version 2.0" << '\n';
    vtkFile << "champ de vitesse" << '\n';
    vtkFile << "ASCII" << '\n';
    vtkFile << "DATASET RECTILINEAR_GRID" << '\n';
    vtkFile << "DIMENSIONS    " << iMax+1 << "    "<< jMax+1 << "    " << 1 << '\n';
    vtkFile << "X_COORDINATES    " << iMax+1 <<"    double" << '\n';
    for (unsigned int i = 0; i<iMax+1; i++)
    {
        vtkFile << X(i) << '\n';
    }
    vtkFile << "Y_COORDINATES    " << jMax+1 <<"    double" << '\n';
    for (unsigned int j = 0; j<jMax+1; j++)
    {
        vtkFile << Y(j) << '\n';
    }
    vtkFile << "Z_COORDINATES    " << 1 << "    double" << '\n';
    vtkFile << 0 << '\n';

    vtkFile << "CELL_DATA    " << iMax*jMax << '\n';
    vtkFile << "SCALARS rho double" << '\n';
    vtkFile << "LOOKUP_TABLE default" << '\n';
    for (unsigned int j = 0; j<jMax; j++)
    {
        for (unsigned int i = 0; i<iMax; i++)
        {
            vtkFile << mRho(i,j) << '\n';
        }
    }

    vtkFile << "VECTORS u double" << '\n';
    for (unsigned int j = 0; j<jMax; j++)
    {
        for (unsigned int i = 0; i<iMax; i++)
        {
            vtkFile << mU(i,j)(0) << "    " << mU(i,j)(1) << "    " << 0 << '\n';
        }
    }

    vtkFile.close();
    emit endCompute();
}

#endif // LBM_HPP
