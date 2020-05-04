
#pragma once

#include<iostream>
#include<complex>

namespace solver
{
    /* ax^2 + bx + c = 0 */
    class RealVariable
    {
    public:
        double a , b , c ;

        /* empty constructor */
        RealVariable();
        RealVariable(double c, double p , double r);
        RealVariable operator+(RealVariable& var);
        RealVariable operator-(RealVariable& var);
        RealVariable operator()(RealVariable& var);

        friend RealVariable operator + (RealVariable const & var1 , RealVariable const & var2);
        friend RealVariable operator + (RealVariable const & var , double num);
        friend RealVariable operator + (double num , RealVariable const & var);

        friend RealVariable operator - (RealVariable const & var1 , RealVariable const & var2);
        friend RealVariable operator - (RealVariable const & var , double num);
        friend RealVariable operator - (double num , RealVariable const & var);

        friend RealVariable operator ^ (RealVariable const & var , int pow);

        friend RealVariable operator * (double num , RealVariable const & var);
        friend RealVariable operator * (RealVariable const & var1 , RealVariable const & var2);
        friend RealVariable operator * (RealVariable const & var1 , int const & var2);
        friend RealVariable operator * (RealVariable const & var1 , double const & var2);
        friend RealVariable operator * ( int const & var2,RealVariable const & var1 );

        friend RealVariable operator / (RealVariable const & var , double num);
        //friend RealVariable operator / (double const & var , RealVariable num);

        friend RealVariable operator == (RealVariable const & var , double num);
        friend RealVariable operator == (RealVariable const & var1 , RealVariable const & var2);
        friend RealVariable operator == (double num , RealVariable const & var);


    };

    class ComplexVariable
    {
    public:
        double a;
        std::complex<double> b , c ;

        /* empty constructor */
        ComplexVariable();
        ComplexVariable(double c, double p , double r);

        ComplexVariable operator()(ComplexVariable& var);

        friend ComplexVariable operator + (ComplexVariable const & var1 , ComplexVariable const & var2);
        friend ComplexVariable operator + (ComplexVariable const & var , std::complex<double> num);
        friend ComplexVariable operator + (ComplexVariable const & var , double num);
        friend ComplexVariable operator + (double num , ComplexVariable const & var);
        friend ComplexVariable operator + (std::complex<double> num , ComplexVariable const & var);


        friend ComplexVariable operator - (ComplexVariable const & var1 , ComplexVariable const & var2);
        friend ComplexVariable operator - (ComplexVariable const & var , double num);
        friend ComplexVariable operator - (ComplexVariable const & var , std::complex<double> num);
        friend ComplexVariable operator - (std::complex<double> num , ComplexVariable const & var);

        friend ComplexVariable operator ^ (ComplexVariable const & var , int pow);

        friend ComplexVariable operator * (double num , ComplexVariable const & var);
        friend ComplexVariable operator * (ComplexVariable const & var,double num );
        friend ComplexVariable operator * (ComplexVariable const & var,ComplexVariable const & var2);

        friend ComplexVariable operator / (ComplexVariable const & var , double num);

        friend ComplexVariable operator == (ComplexVariable const & var , double num);
        friend ComplexVariable operator == (double num,ComplexVariable const & var);
        friend ComplexVariable operator == (ComplexVariable const & var1 , ComplexVariable const & var2);

    };
    double solve(RealVariable var);
    std::complex<double> solve(ComplexVariable var);

}

 int saparate(double rhs);
 int multens(double fractpart);
