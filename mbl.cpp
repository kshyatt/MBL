#include "Eigen/Eigen/Dense"
#include "Eigen/Eigen/Eigenvalues"
#include "Eigen/unsupported/Eigen/MatrixFunctions"
#include <iostream>
#include <cmath>

using namespace std;

int main()
{

    double epsilon = 1.0;
    double lambda  = 0.05;
    double omega   = -1.0;
    double Omega   = lambda*lambda + (omega + 2.0 * epsilon)*(omega + 2.0 * epsilon)/4.0;
    double delta   = 2.0 * M_PI/(100.0*Omega);
    Eigen::Matrix2cd SigmaZ, SigmaX, SigmaY;

    SigmaZ.real() << 1, 0, 0, -1;
    SigmaX.real() << 0, 1, 1, 0;
    SigmaY.imag() << 0, -1, 1, 0;
    //cout<<"SigmaZ is: "<<endl<<SigmaZ<<endl; 
    //cout<<"SigmaY is: "<<endl<<SigmaZ<<endl; 
    //cout<<"SigmaX is: "<<endl<<SigmaX<<endl; 

    Eigen::Matrix2cd TimeIndependent;
    Eigen::Matrix2cd TimeDependent;
    Eigen::Matrix2cd Intermediate;
    Intermediate.imag() = -epsilon * delta * SigmaZ.real() / 2.0; 
    TimeIndependent = Intermediate.exp();   
 
    //cout<<"The time independent piece is: "<<endl<<TimeIndependent<<endl; 
    Eigen::Vector2cd Wavefunction;
    Eigen::Vector2cd WavefunctionNext;
    Wavefunction << 1, 0;


    Eigen::Matrix2cd TimeIndependentExact;
    Eigen::Matrix2cd TimeDependentExact;
    
    //cout<<"The time independent piece is: "<<endl<<TimeIndependent<<endl; 

    for( int i = 0; i < 100; i++ )
    {
        Intermediate.real() = lambda * ( sin( (i + 1) * delta * omega ) - sin( i * delta * omega ) ) * (SigmaX.imag())/(omega) + lambda * ( cos( (i + 1) * delta * omega ) - cos( i * delta * omega ) ) * SigmaY.imag()/(omega); 
        Intermediate.imag() = - lambda * ( sin( (i + 1) * delta * omega ) - sin( i * delta * omega ) ) * (SigmaX.real())/(omega) - lambda * ( cos( (i + 1) * delta * omega ) - cos( i * delta * omega ) ) * SigmaY.real()/(omega); 
        TimeDependent = Intermediate.exp();

//        cout<<"The time dependent piece is: "<<endl<<TimeDependent<<endl; 
        WavefunctionNext = TimeIndependent * TimeDependent * TimeIndependent * Wavefunction;
        Wavefunction     = WavefunctionNext;
        cout<<Wavefunction(0).real()*Wavefunction(0).real()  + Wavefunction(0).imag()*Wavefunction(0).imag()<<" ";
        cout<<Wavefunction(1).real()*Wavefunction(1).real()  + Wavefunction(1).imag()*Wavefunction(1).imag()<<endl;
    }

    /*Eigen::Matrix2cd Final;
    Final.real() = epsilon * SigmaZ + lambda * cos( omega * delta * 50 ) * SigmaX;

    Eigen::ComplexEigenSolver< Eigen::Matrix2cd > Solver;
    Solver.compute( Final );*/

    //cout<<"The evolution groundstate is: "<<endl<<Wavefunction<<endl; 
    //cout<<"The true groundstate is: "<<endl<<Solver.eigenvectors().col(0)<<endl; 
    //cout<<"The quotient is: "<<endl<<Wavefunction(0)/Solver.eigenvectors()(0,0)<<endl<<Wavefunction(0)/Solver.eigenvectors()(1,0)<<endl; 

    //double coeff = lambda * lambda * sin( sqrt( lambda * lambda + ( omega + 2 * epsilon )*(omega + 2* epsilon)/4.0 ) * 50 * delta ) * sin( sqrt( lambda * lambda + ( omega + 2 * epsilon )*(omega + 2* epsilon)/4.0 ) * 50 * delta ) / ( lambda * lambda + (omega + 2 * epsilon) * (omega + 2 * epsilon)/4.0 );

    //cout<<"Exact c_2: "<<sqrt(coeff)<<endl;

    return 0;
}
