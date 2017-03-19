#ifndef LBM_HPP
#define LBM_HPP

#include "grid.h"
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include "boost/multi_array.hpp"
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

    std::unique_ptr<boost::numeric::ublas::matrix<float> > mRho;
    std::unique_ptr<boost::numeric::ublas::matrix<float> > mUx;
    std::unique_ptr<boost::numeric::ublas::matrix<float> > mUy;
    std::unique_ptr<boost::multi_array<float,3> > mGn;
    std::unique_ptr<boost::multi_array<float,3> > mGnp1;
    std::unique_ptr<boost::multi_array<float,3> > mGeq;
    std::unique_ptr<boost::numeric::ublas::matrix<float> > mCl;

    static constexpr float sNu=0.01, sUin=0.4;
    static constexpr int sQMax = 9;
    static constexpr float sC = 1.73205080756887729352744634150587236694280525381038062805580, sEps = 0.2; //sqrt(3)
    static constexpr float sDx = 2./jMax, sDt=sDx*sEps/sC;
    static constexpr float sEta = 1./(sEps*sEps*sNu/sDt+0.5);
    static constexpr float sCfl = sDt*sNu/(sDx*sDx);

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
        ONLYRELEASE(return (2*i+1)*sDx/2 -1;)
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
    mRho =std::unique_ptr<boost::numeric::ublas::matrix<float> >(new boost::numeric::ublas::matrix<float>(iMax,jMax));
    mUx =std::unique_ptr<boost::numeric::ublas::matrix<float> >(new boost::numeric::ublas::matrix<float>(iMax,jMax));
    mUy =std::unique_ptr<boost::numeric::ublas::matrix<float> >(new boost::numeric::ublas::matrix<float>(iMax,jMax));
    mGn =std::unique_ptr<boost::multi_array<float,3> >(new boost::multi_array<float,3>(boost::extents[iMax+2][jMax+2][sQMax]));
    mGnp1 =std::unique_ptr<boost::multi_array<float,3> >(new boost::multi_array<float,3>(boost::extents[iMax+2][jMax+2][sQMax]));
    mGeq =std::unique_ptr<boost::multi_array<float,3> >(new boost::multi_array<float,3>(boost::extents[iMax][jMax][sQMax]));

    std::fill(mRho->data().begin(),mRho->data().end(), 1.0f);
    mUy->clear();
    for(unsigned i  = 1; i<iMax; i++)
    {
        for (unsigned j = 1; j<jMax; j++)
        {
            if ((*mObstacle)(i,j)==0)
            {
                (*mUx)(i,j) = 0;
            }else
            {
                (*mUx)(i-1,j-1) = (1. - Yc(j-1)* Yc(j-1))*sUin;
            }

            for (unsigned q = 0; q<sQMax; q++)
            {
                (*mGn)[i][j][q] = omega(q)*(*mRho)(i-1,j-1)*(1+(*mUx)(i-1,j-1)*speeds(q,0)+(*mUy)(i-1,j-1)*speeds(q,1));
                (*mGnp1)[i][j][q] = (*mGn)[i][j][q];
            }
        }
    }

    std::stringstream buf;
    std::cout << "iteration " << 0 << std::endl;
    buf.str("");
    buf << "u_" << 0 <<".vtk";
    saveVtk(buf.str());
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
            vtkFile << (*mRho)(i,j) << std::endl;
        }
    }

    vtkFile << "VECTORS u double" << std::endl;
    for (unsigned int j = 0; j<jMax; j++)
    {
        for (unsigned int i = 0; i<iMax; i++)
        {
            vtkFile << (*mUx)(i,j) << "    " << (*mUy)(i,j) << "    " << 0 << std::endl;
        }
    }

    vtkFile.close();
    emit endCompute();
}


template<const unsigned int iMax, const unsigned int jMax>
void LBM<iMax, jMax>::Iterate(int nIter)
{
    std::stringstream buf;
    for (int iter = 1; iter<nIter; iter++)
    {
//            time += dt;
        //-- collision
        for (unsigned int i = 0; i<iMax; i++)
        {
            for (unsigned int j = 0; j<jMax; j++)
            {
                (*mRho)(i,j) = 0;
                for (int q =0; q<sQMax; q++)
                {
                    (*mRho)(i,j) += (*mGn)[i+1][j+1][q];
                }

                //-- pour que la visualisation soit plus claire, on met (*mUx) et uy a 0 dans le solide

                if ((*mObstacle)(i,j) == 1)
                {
                    (*mUx)(i,j) = 1./(*mRho)(i,j)*sC*((*mGn)[i+1][j+1][1]- (*mGn)[i+1][j+1][3] + (*mGn)[i+1][j+1][5] - (*mGn)[i+1][j+1][6] - (*mGn)[i+1][j+1][7] + (*mGn)[i+1][j+1][8]);
                    (*mUy)(i,j) = 1./(*mRho)(i,j)*sC*((*mGn)[i+1][j+1][2]- (*mGn)[i+1][j+1][4] + (*mGn)[i+1][j+1][5] + (*mGn)[i+1][j+1][6] - (*mGn)[i+1][j+1][7] - (*mGn)[i+1][j+1][8]);
                    //- distribution d'equilibre pour Stokes et Navier-Stokes
                    for (int q = 0; q<sQMax; q++)
                    {
                        #ifdef STOKES
                            (*mGeq)[i][j][q] = omega(q)*(*mRho)(i,j)*(1+ (*mUx)(i,j)*speeds(q,0)+(*mUy)(i,j)*speeds(q,1));
                        #else
                        #ifdef NAVIERSTOKES
                            (*mGeq)[i][j][q] = omega(q)*(*mRho)(i,j)*(1+(*mUx)(i,j)*speeds(q,0)+ (*mUy)(i,j)*speeds(q,1)- 0.5*((*mUx)(i,j)*(*mUx)(i,j)+(*mUy)(i,j)*(*mUy)(i,j)) + 0.5*((*mUx)(i,j)*speeds(q,0)+ (*mUy)(i,j)*speeds(q,1))*((*mUx)(i,j)*speeds(q,0)+ (*mUy)(i,j)*speeds(q,1)));
                        #endif
                        #endif
                    }
                    //-- transport des noeuds interieurs (on voit ici pourquoi i et j vont a iMax+2, jMax+2 dans l'allocation de (*mGnp1) : c'est pour que le noeuds des bords soient affectés de la meme facon
                    (*mGnp1)[i+1][j+1][0] = (1 - sEta)*(*mGn)[i+1][j+1][0] + sEta*(*mGeq)[i][j][0];
                    (*mGnp1)[i+2][j+1][1] = (1 - sEta)*(*mGn)[i+1][j+1][1] + sEta*(*mGeq)[i][j][1];
                    (*mGnp1)[i+1][j+2][2] = (1 - sEta)*(*mGn)[i+1][j+1][2] + sEta*(*mGeq)[i][j][2];
                    (*mGnp1)[i][j+1][3]   = (1 - sEta)*(*mGn)[i+1][j+1][3] + sEta*(*mGeq)[i][j][3];
                    (*mGnp1)[i+1][j][4]   = (1 - sEta)*(*mGn)[i+1][j+1][4] + sEta*(*mGeq)[i][j][4];
                    (*mGnp1)[i+2][j+2][5] = (1 - sEta)*(*mGn)[i+1][j+1][5] + sEta*(*mGeq)[i][j][5];
                    (*mGnp1)[i][j+2][6]   = (1 - sEta)*(*mGn)[i+1][j+1][6] + sEta*(*mGeq)[i][j][6];
                    (*mGnp1)[i][j][7]     = (1 - sEta)*(*mGn)[i+1][j+1][7] + sEta*(*mGeq)[i][j][7];
                    (*mGnp1)[i+2][j][8]   = (1 - sEta)*(*mGn)[i+1][j+1][8] + sEta*(*mGeq)[i][j][8];

                }else  //obst(i,j) == 0
                {
                    (*mGnp1)[i+1][j+1][0] = (*mGn)[i+1][j+1][0];
                    (*mGnp1)[i+2][j+1][1] = (*mGn)[i+1][j+1][3];
                    (*mGnp1)[i+1][j+2][2] = (*mGn)[i+1][j+1][4];
                    (*mGnp1)[i][j+1][3]   = (*mGn)[i+1][j+1][1];
                    (*mGnp1)[i+1][j][4]   = (*mGn)[i+1][j+1][2];
                    (*mGnp1)[i+2][j+2][5] = (*mGn)[i+1][j+1][7];
                    (*mGnp1)[i][j+2][6]   = (*mGn)[i+1][j+1][8];
                    (*mGnp1)[i][j][7]     = (*mGn)[i+1][j+1][5];
                    (*mGnp1)[i+2][j][8]   = (*mGn)[i+1][j+1][6];
                }


            }
        }

        //-- CL Sud
        for (unsigned int i = 2; i<iMax;i++)
        {
            //- vitesses rentrantes : BB
            (*mGnp1)[i][1][2] = (*mGnp1)[i][1][4];
            (*mGnp1)[i][1][6] = (*mGnp1)[i][1][8];
            (*mGnp1)[i][1][5] = (*mGnp1)[i][1][7];

        }
        //-- CL Nord
        for (unsigned int i = 2; i<iMax;i++)
        {
            //- vitesses rentrantes : BB
            (*mGnp1)[i][jMax][4] = (*mGnp1)[i][jMax][2];
            (*mGnp1)[i][jMax][7] = (*mGnp1)[i][jMax][5];
            (*mGnp1)[i][jMax][8] = (*mGnp1)[i][jMax][6];
        }
        //-- CL Est : Neumann
        for (unsigned int j =1; j<jMax+1; j++)
        {
            for (int q =0; q<sQMax; q++)
            {
                (*mGnp1)[iMax][j][q] = (*mGnp1)[iMax-1][j][q];
            }
        }
        //-- CL Ouest : Dirichlet
        for (unsigned int j =0; j<jMax; j++)
        {
            for (int q =0; q<sQMax; q++)
            {
                (*mGnp1)[1][j+1][q] = omega(q)*(*mRho)(0,j)*(1+(1-Yc(j)*Yc(j))*sUin*speeds(q,0));
            }
        }

        mGn.swap(mGnp1);

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
