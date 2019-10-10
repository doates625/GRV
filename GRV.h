/**
 * @file GRV.h
 * @brief Class for implementing 1-dimensional Gaussian random variables
 * @author Dan Oates (WPI Class of 2020)
 * 
 * This class represents 1D Gaussian random varaibles by their mean 'mean'
 * and variance 'var' in 32-bit floating-point format. Mathematical functions
 * and operators are included which propagate the variance through each
 * calculation via linearization. Note that protections are not included for 
 * potentially-singular or invalid operations such as inversion and sqrt.
 */

/**
 * Class Declaration
 */
class GRV
{
public:
	GRV(float mean, float var);
	GRV();
	float mean, var;
};

/**
 * Function Declarationss
 */
GRV sin(const GRV& g);
GRV cos(const GRV& g);
GRV tan(const GRV& g);
GRV asin(const GRV& g);
GRV acos(const GRV& g);
GRV atan(const GRV& g);
GRV sinh(const GRV& g);
GRV cosh(const GRV& g);
GRV tanh(const GRV& g);
GRV sqrt(const GRV& g);
GRV exp(const GRV& g);
GRV log(const GRV& g);
GRV atan2(const GRV& gy, const GRV& gx);
GRV fuse(const GRV& lhs, const GRV& rhs);

/**
 * Operator Declarations
 */
GRV operator-(const GRV& g);
GRV operator+(const GRV& g, float n);
GRV operator-(const GRV& g, float n);
GRV operator*(const GRV& g, float n);
GRV operator/(const GRV& g, float n);
GRV operator^(const GRV& g, float n);
GRV operator+(const GRV& lhs, const GRV& rhs);
GRV operator-(const GRV& lhs, const GRV& rhs);
GRV operator*(const GRV& lhs, const GRV& rhs);
GRV operator/(const GRV& lhs, const GRV& rhs);
GRV operator^(const GRV& lhs, const GRV& rhs);
