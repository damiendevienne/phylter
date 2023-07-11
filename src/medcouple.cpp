// extract from the package mrfDepth - 07/2023

#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <functional>
#include <fstream>
#include <iostream>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <vector>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;


extern double mlmccN(const double z[],const int in,const int dr);

extern "C"{
    void medcoupleC(double* xi,int* n,double* mcaap,int* dr){
        const int lin=*n,dryn=*dr;
        *mcaap=mlmccN(xi,lin,dryn);		
    }
}

