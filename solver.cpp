#include"solver.hpp"
#include<iostream>
#include <math.h>

using namespace std;
using namespace solver;

 int saparate(double param) {
    double fractpart, intpart;
    fractpart = modf (param , &intpart);
    int ret = multens(fractpart);
    return ret;
}

 int multens(double fractpart) {
    for (int i = 0; i <7 ; ++i) {
        fractpart= fractpart*10;
    }
    int re = (int) fractpart;
    return re;
}

RealVariable::RealVariable()
{
    this->a = 0;
    this->b = 1;
    this->c = 0;

}

RealVariable::RealVariable(double c , double p , double r)
{
    this->a = c;
    this->b = p;
    this->c = r;
}

RealVariable RealVariable::operator+(RealVariable& rv)
{
    this->b +=1;
    return *(this);
}


RealVariable RealVariable::operator-(RealVariable& rv)
{
    this->b -=1;
    return *(this);
}

RealVariable RealVariable::operator()(RealVariable& rv)
{
    RealVariable ret;
    ret.b = 0;

    ret.a = rv.a;
    ret.b = rv.b;
    ret.c = rv.c;

    return ret;

}


RealVariable solver::operator + (RealVariable const & var1 , RealVariable const & var2)
{
    auto *new_var = new solver::RealVariable(var1);
    new_var->a +=var2.a;
    new_var->b +=var2.b;
    new_var->c +=var2.c;

    return *(new_var);

}
RealVariable solver::operator + (RealVariable const & var1, double var2)
{
    auto *new_var = new solver::RealVariable(var1);
    new_var->c = new_var->c + var2;

    return *(new_var);
}
RealVariable solver::operator + (double num , RealVariable const & var)
{
    auto *new_var = new solver::RealVariable(var);
    new_var->c = new_var->c + num;

    return *(new_var);
}


RealVariable solver::operator - (RealVariable const & var1 , RealVariable const & var2)
{
    auto *new_var = new solver::RealVariable(var1);

    new_var->a = new_var->a -var2.a;
    new_var->b = new_var->b -var2.b;
    new_var->c = new_var->c -var2.c;

    return *(new_var);
}
RealVariable solver::operator - (RealVariable const & var , double num)
{
    auto *new_var = new solver::RealVariable(var);
    new_var->c = new_var->c - num;
    return *(new_var);
}
RealVariable solver::operator - (double num , RealVariable const & var)
{
    auto *new_var = new solver::RealVariable(var);
    new_var->c = new_var->c - num;

    return *(new_var);
}

RealVariable solver::operator ^ (RealVariable const & var , int b)
{
    if(b > 2)
        throw std::invalid_argument("ber can't be greater than 2\n");

    ///DO NOTHING
    if(b==1)
        return var;

    /// make this realvariable =1...
   if(b==0){
       auto *new_var = new solver::RealVariable;
       new_var->b = 0;
       new_var->c = 1;
       return *(new_var);
   }

    auto *new_var = new solver::RealVariable(var);
    new_var->b = 0;
    new_var->a = 1;

    return *(new_var);
}


RealVariable solver::operator / (RealVariable const & var , double num)
{

    if(num == 0)
        throw std::invalid_argument("Cant divide by 0 .\n");

    auto *new_var = new solver::RealVariable(var);
    new_var->b = new_var->b / num;

    return *(new_var);

}

RealVariable solver::operator == (RealVariable const & var , double num)
{
    auto *new_var = new solver::RealVariable(var);
    new_var->c = new_var->c - (num);
    return *(new_var);

}
RealVariable solver::operator == (double num , RealVariable const & var)
{
    auto *new_var = new solver::RealVariable(var);
    new_var->c = new_var->c - (num);
    return *(new_var);

}
RealVariable solver::operator == (RealVariable const & var1 , RealVariable const & var2)
{
    auto *new_var = new solver::RealVariable(var1);

    new_var->a = new_var->a -var2.a;
    new_var->b = new_var->b -var2.b;
    new_var->c = new_var->c -var2.c;

    return *(new_var);

}

RealVariable solver::operator*(const RealVariable &var1, const RealVariable &var2) {
    auto *new_var = new solver::RealVariable(var1);

    new_var->a = new_var->a *var2.a;
    new_var->b = new_var->b *var2.b;
    new_var->c = new_var->c *var2.c;

    return *(new_var);
}
RealVariable solver::operator*(const RealVariable &var1, const int &var2) {
    auto *new_var = new solver::RealVariable(var1);
    new_var->b = 1;

    new_var->a = new_var->a * var2;
    new_var->b = new_var->b * var2;
    new_var->c = new_var->c * var2;

    return *(new_var);
}
RealVariable solver::operator*(const RealVariable &var1, const double &var2) {
    auto *new_var = new solver::RealVariable(var1);
    new_var->b = 1;

    new_var->a = new_var->a * var2;
    new_var->b = new_var->b * var2;
    new_var->c = new_var->c * var2;

    return *(new_var);
}
RealVariable solver::operator*(const int &var2, const RealVariable &var1) {

    auto *new_var = new solver::RealVariable(var1);
    //new_var->b = 1;

    new_var->a = new_var->a * var2;
    new_var->b = new_var->b * var2;
    new_var->c = new_var->c * var2;

    return *(new_var);
}
RealVariable solver::operator * (double num , RealVariable const & var)
{
    auto *new_var = new solver::RealVariable(var);
    new_var->b = 1;

    new_var->a = new_var->a * num;
    new_var->b = new_var->b * num;
    new_var->c = new_var->c * num;

    return *(new_var);
}

//////////////////////////////////////////////////Complex/////////////////////////////////////
ComplexVariable::ComplexVariable()
{
    this->a = 0;
    this->b = 1;
    this->c = 0;

}


ComplexVariable::ComplexVariable(double c, double p, double r) {
    this->a = c;
    this->b = p;
    this->c = r;
}

ComplexVariable ComplexVariable::operator()(ComplexVariable& var)
{
    ComplexVariable ans;
    return ans;

}

ComplexVariable solver::operator + (ComplexVariable const & var1 , ComplexVariable const & var2)
{
    auto *new_var = new solver::ComplexVariable(var1);
    new_var->a += var2.a;
    new_var->b += var2.b;
    new_var->c += var2.c;
    return *(new_var);
}
ComplexVariable solver::operator + (ComplexVariable const & var , double num)
{
    auto *new_var = new solver::ComplexVariable(var);
    new_var->c = new_var->c + num;

    return *(new_var);
}
ComplexVariable solver::operator + (ComplexVariable const & var , std::complex<double> num)
{
    auto *new_var = new solver::ComplexVariable(var);
    new_var->c = new_var->c + num;

    return *(new_var);
}
ComplexVariable solver::operator + (double num , ComplexVariable const & var)
{
    auto *new_var = new solver::ComplexVariable(var);
    new_var->c = new_var->c + num;

    return *(new_var);
}
ComplexVariable solver::operator + (std::complex<double> num, const ComplexVariable &var) {
    auto *new_var = new solver::ComplexVariable(var);
    new_var->c = new_var->c + num;

    return *(new_var);
}

ComplexVariable solver::operator - (ComplexVariable const & var1 , ComplexVariable const & var2)
{
    auto *new_var = new solver::ComplexVariable(var1);
    new_var->a -= var2.a;
    new_var->b -= var2.b;
    new_var->c -= var2.c;
    return *(new_var);

}
ComplexVariable solver::operator - (ComplexVariable const & var , double num)
{
    auto *new_var = new solver::ComplexVariable(var);
    new_var->c = new_var->c - num ;

    return *(new_var);
}
ComplexVariable solver::operator - (ComplexVariable const & var , std::complex<double> num)
{
    auto *new_var = new solver::ComplexVariable(var);
    new_var->c = new_var->c - num;

    return *(new_var);

}
ComplexVariable solver::operator - (std::complex<double> num, const ComplexVariable &var) {
    auto *new_var = new solver::ComplexVariable(var);
    new_var->c = new_var->c - num;

    return *(new_var);
}


ComplexVariable solver::operator ^ (ComplexVariable const & var , int b)
{
    if(b > 2)
        throw std::invalid_argument("ber can't be greater than 2\n");

    ///DO NOTHING
    if(b==1)
        return var;

    /// make this realvariable =1...
    if(b==0){
        auto *new_var = new solver::ComplexVariable;
        new_var->a = 0;
        new_var->b = 0;
        new_var->c = 1;
        return *(new_var);
    }

    auto *new_var = new solver::ComplexVariable(var);
    new_var->b = 0;
    new_var->a = 1;

    return *(new_var);

}

ComplexVariable solver::operator * (double num , ComplexVariable const & var)
{
    auto *new_var = new solver::ComplexVariable(var);

    new_var->a = new_var->a * num;
    new_var->b = new_var->b * num;
    new_var->c = new_var->c * num;

    return *(new_var);
}
ComplexVariable solver::operator * ( ComplexVariable const &var, double num) {
    auto *new_var = new solver::ComplexVariable(var);

    new_var->a = new_var->a * num;
    new_var->b = new_var->b * num;
    new_var->c = new_var->c * num;

    return *(new_var);
}
ComplexVariable solver::operator * (const ComplexVariable &var1, const ComplexVariable &var2) {
    auto *new_var = new solver::ComplexVariable(var1);
    new_var->a *= var2.a;
    new_var->b *= var2.b;
    new_var->c *= var2.c;
    return *(new_var);
}


ComplexVariable solver::operator / (ComplexVariable const & var , double num)
{
    if(num == 0)
        throw std::invalid_argument("Cant divide by 0 .\n");

    auto *new_var = new solver::ComplexVariable(var);
    new_var->b = new_var->b / num;

    return *(new_var);

}

ComplexVariable solver::operator == (ComplexVariable const & var1 , ComplexVariable const & var2)
{
    auto *new_var = new solver::ComplexVariable(var1);

    new_var->a = new_var->a -var2.a;
    new_var->b = new_var->b -var2.b;
    new_var->c = new_var->c -var2.c;

    return *(new_var);
}
ComplexVariable solver::operator == (ComplexVariable const & var , double num)
{
    auto *new_var = new solver::ComplexVariable(var);
    new_var->c = new_var->c - num;

    return *(new_var);

}
ComplexVariable solver::operator == (double num, const ComplexVariable &var) {
    return ComplexVariable();
}


/////Solvers////

double solver::solve(RealVariable var)
{
    float x1=0;    float x2=0;    float realPart=0;    float imaginaryPart =0;

    float a = var.a;
    float b = var.b;
    float c = var.c;

    if(a==0){
        if (b==0)
            throw std::invalid_argument("Roots are complex and different.\n");

        float res = (-1*c)/b;
        return res;
    }

    float discriminant = b*b - 4*a*c;

    if (discriminant > 0) {
        x1 = (-b + sqrt(discriminant)) / (2*a);
        x2 = (-b - sqrt(discriminant)) / (2*a);
        return x1;
    }

    else if (discriminant == 0) {
        x1 = (-b + sqrt(discriminant)) / (2*a);
        return x1;
    }

    else {
        realPart = -b/(2*a);
        imaginaryPart =sqrt(-discriminant)/(2*a);
        throw std::invalid_argument("Roots are complex and different.\n");
    }
}

std::complex<double> solver::solve(ComplexVariable var)
{
    complex<double> x1=0;    float realPart=0;    float imaginaryPart =0;

    if(var.a==0){
        if (var.b.real()==0 & var.b.imag()==0)
            throw std::invalid_argument("Roots are complex and different.\n");

        complex<double> res = ((var.c))/var.b;
        res *=-1;
        return res;
    }

    complex<double> discriminant = var.b * var.b - 4 * var.a * var.c;

    if (discriminant.real() > 0) {
        sqrt(discriminant);
        x1 = (-var.b + sqrt(discriminant)) / (2*var.a);
        return x1;
    }

    else {
        complex<double> s =complex<double>(-1) * var.b;
        return s + sqrt(discriminant)/(var.a*2);
    }
}