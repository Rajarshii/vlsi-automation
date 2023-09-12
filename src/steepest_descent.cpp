// Author: Rajarshi Das
// File name: steepest_descent.cpp
// Details: Steepest Descent Heuristic Implementation using line serach in C++

// Import the functions

#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

// Define the function f(x1, x2) to be minimized
double function_f(const vector<double>& x) {
    //std::cout << "function_f : (" << x[0] << "," << x[1] << ")" <<std::endl;
    return x[0]*(x[0] - 11) + x[1]*x[1] + x[0]*(x[1] + 9);
}


// Calculating the gradient of the function
vector<double> gradient_f(const vector<double>& x)
{
    vector<double> nabla_f (2,0);
    nabla_f[0] = 2*x[0] - 2 + x[1] ;
    nabla_f[1] = 2*x[1] + x[0];
    return nabla_f;
}

// Calculating the L2 Norm
double norm(const vector<double>& v) {
    double sq_sum = 0;
    for(int i=0; i < v.size(); i++) {
        sq_sum  += v[i]*v[i];
    }
    return sqrt(sq_sum);
}

// Calculate the gradient descent direction from the gradient
vector<double> find_descent_direction(const vector<double>& v){
    vector<double> dk (2,0);
    dk[0] = -v[0];
    dk[1] = -v[1];
    return dk;
}

// Line search algorithm implementation
// Using the golden-section method with the initial interval of unc ertainity range being [0.0,1.0]
// We also give the final allowable uncertainity length as 0.000001
// We pre-compute the inverse of the golden ration,  phi_inv = 0.618

double line_search(double x_1, double x_2, vector<double>& gradient) {
    const double l = 0.000001;      // Allowable final length of uncertainity
    const double phi_inv = 0.618;   // (1/Golden ratio)
    
    double a = 0.0;    //Take a minimum value for step size
    double b = 1.0;    //Maximum value of step size

    // Initialization
    double left =  a + (1.0-phi_inv)*(b-a);
    double right = a + phi_inv*(b-a);

    // Main Steps
    while(b-a >= l) {
        // Evaluate the function values at the intermediate points
        double left_fn = function_f({(x_1 + left*gradient[0]),(x_2 + left*gradient[1])});
        double right_fn = function_f({x_1 + right*gradient[0],x_2 + right*gradient[1]});
       
        if(left_fn > right_fn){
            a = left;
            left = right;
            right = a + phi_inv*(b-a);
        }
        else {
            b = right;
            right = left;
            left = a + (1.0 - phi_inv) * (b - a);
        }
    }

    return 0.5*(a+b);

}

int main() {

    // termination criteria to be checked
    const double epsilon = 0.000001;

    // Define the starting points
    double x_1 = 0.0;
    double x_2 = 0.0;
    std::cout << "Initial Solution : (" << x_1 << "," << x_2 << ")" << " | Termination Condition : " << epsilon <<std::endl;
    // Intialize k = 1
    int iter_k = 1;

    std::cout << "Starting Steepest Descent Algorithm with Line Search\n\n" << std::endl;

    // Declare and initialize the gradient vector
    vector<double> gradient = gradient_f({x_1,x_2});

    // Loop on the termination condition epsilon
    while(!(norm(gradient) <= epsilon)) {

        // Calculate the direction d_k
        vector<double> d_k = find_descent_direction(gradient);

        // Use line search to find the step size lambda_k
        double lambda_k = line_search(x_1,x_2,d_k); //Initialize the step size to a small enough value

        // Update x_1 and x_2 for next iteration
        x_1 += lambda_k*(d_k[0]);
        x_2 += lambda_k*(d_k[1]);

        // Keep track of the intermediate results
        std::cout << "Iteration #" << iter_k << ": x1 = " << x_1 << ", x2 = " << x_2;
        std::cout << ", Step Size = " << lambda_k << ", Moving Direction = (" << d_k[0] << "," << d_k[1] << "), function value = " << function_f({x_1,x_2}) << std::endl;

        // Evaluate the gradient for next iteration
        gradient = gradient_f({x_1,x_2});

        //Keep track of the iteration number
        iter_k++;
    }
    std::cout << "\n\nSteepest Descent Algorithm Converged! \n" ; 
    std::cout << "The minimum value of f(x1,x2) is at (x1,x2) = (" << x_1 << "," << x_2 << ") and value is " << function_f({x_1,x_2}) << std::endl;
    
}
