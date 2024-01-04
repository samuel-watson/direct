# direct
A variety of optimisation algorithms including derivative free options: BOBYQA, NEWUOA, L-BFGS, L-BFGS-B, DIRECT. The algorithms themselves (apart from DIRECT) are sourced from other authors, with the details included in the relevant files. The main contribution of this repository is a simple approach to function binding and calling the optimiser.

## Use
The `optim` class is instantiated and then a function and starting values need to be passed to it before calling the algorithm. The BOBYQA, NEWUOA, and DIRECT algorithms currently take functions with the signature `template <typename T> T(const std::vector<T>&)` where the vector argument contains the parameters of the function. However, there are conventient typedefs for the algorithms including `bobyqad` (double), `bobyqaf` (float), `newuoad`, `newuoaf`, `directd`, and `directf`. A simple example: 
```
auto fn(const std::vector<double>& x) -> double { return 2*x[0]*x[0] - x[1]; }
auto d = directd{};
d.fn<&fn>();
std::vector<double> start = {0,0};
std::vector<double> range = {2,2};
d.set_bounds(start,range);
d.optim();
```
The class can also accept member functions of class objects by passing a reference to the class object in the `fn` argument.
```
struct func {
   auto fn(const std::vector<double>& x) -> double { return 2*x[0]*x[0] - x[1]; }
}
auto cl = func{};
std::vector<double> start = {0,0};
bobyqad d(start);
d.fn<&func::fn,func>(&cl);
d.optim()
```
The function `optim()` runs the algorithm.

For L-BFGS and L-BFGS-B the Eigen library is used and the function signature required is `double(const VectorXd&, VectorXd&)` where the first argument contains the parameter values for the function and the second is a non-const reference to a vector that the gradient will be updated in. 

An example of the use of these functions is in our [glmmrBase](https://github.com/samuel-watson/glmmrBase/blob/57c6bbd3971802139da6a4a92e41fb4e589b0f7d/inst/include/glmmr/modeloptim.hpp#L148) package.

## Sources
Sources for the algorithms.
 - [NEWUOA](https://github.com/elsid/newuoa-cpp/)
 - [BOBYQA]()
 - [L-BFGS/L-BFGS-B](https://github.com/yixuan/LBFGSpp)


Further details and refinements will be added here later.
