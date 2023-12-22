# direct
DIRECT and related optimisation algorithms

## Use
The `direct` class is instantiated and then a function and starting values need to be passed to it before calling the algorithm. 
```
auto fn(std::vector<double> x) -> double { return 2*x[0]*x[0] - x[1]; }
auto d = direct<double>{};
d.fn<&fn>();
std::vector<double> start = {0,0};
std::vector<double> range = {2,2};
d.set_bounds(start,range);
d.optim();
```

The class can also accept member functions of class objects by passing a reference to the class object in the `fn` argument.
```
struct func {
   auto fn(std::vector<double> x) -> double { return 2*x[0]*x[0] - x[1]; }
}
auto cl = func{};
auto d = direct<double>{};
d.fn<&func::fn>(&cl);
std::vector<double> start = {0,0};
std::vector<double> range = {2,2};
d.set_bounds(start,range);
d.optim();
```

The function `optim()` runs the algorithm. Further details and refinements will be added here later.
