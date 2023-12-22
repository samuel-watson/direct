# direct
DIRECT and related optimisation algorithms. For a detailed description and discussion of DIRECT and related algorithms see [Jones and Martins (2019)](https://doi.org/10.1007/s10898-020-00952-6).

## Use
The `direct` class is instantiated and then a function and starting values need to be passed to it before calling the algorithm. The full class template is:
```
template <typename T>
class direct<T(const std::vector<T>&)>;
```
However, `directd` and `directf` are convenient typedefs for these classes for double and float, respectively. The class currently only accepts functions with the signature `T(*)(const std::vector<T>&)`. 
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
auto d = directd{};
d.fn<&func::fn,func>(&cl);
std::vector<double> start = {0,0};
std::vector<double> range = {2,2};
d.set_bounds(start,range);
d.optim();
```

The function `optim()` runs the algorithm. Further details and refinements will be added here later.
