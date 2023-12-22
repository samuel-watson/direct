#pragma once

// #define R_BUILD

#include <vector>
#include <functional>
#include <exception>
#include <map>
#include <cmath>
#include <random>
#include <queue>
#include <type_traits>
#include <memory>
#ifdef R_BUILD
#include <Rcpp.h> // for printing to R
#endif

enum class Position {
  Lower,
  Middle,
  Upper
};

inline size_t random_index(size_t max_index) {
  std::random_device                      rand_dev;
  std::mt19937                            generator(rand_dev());
  std::uniform_int_distribution<size_t>   distr(0, max_index);
  return distr(generator);
}

// hyperrectangle class - coordinates should be in [0,1]^D
template <typename T>
class Rectangle {
public:
  const int         dim;
  std::vector<T>    min_x;
  std::vector<T>    max_x;
  bool              potentially_optimal = false;
  bool              no_update = false;
  
  Rectangle(const int dim_) : dim(dim_), min_x(dim), max_x(dim) {};
  Rectangle(const std::vector<T>& min_x_, const std::vector<T>& max_x_) : dim(min_x_.size()), min_x(min_x_), max_x(max_x_) {};
  Rectangle(const Rectangle& x) = default;
  auto operator=(const Rectangle& x) -> Rectangle& = default;
  
  // functions
  std::vector<T>  centroid(); // returns the centroid
  std::vector<T>  centroid(const size_t& dim, const T& delta); // returns the centroid offset by delta in dimension dim
  void            unit_hyperrectangle(); //sets the rectangle to be the unit hyper-rectangle
  T               longest_side(); // returns 0.5 times the longest side
  T               dim_size(const size_t& dim); // returns the size of the dimension
  void            trim_dimension(const size_t& dim, const Position pos); // reduces the dimension to a third of the size - either lower, middle or upper
};

// SizeMap class used to capture unique ordered rectangle sizes along with the best function values and 
// a pointer to the relevant rectangle
template <typename T>
class SizeMap {
  std::set<T>             size;  // unique ordered rectangle sizes
  std::vector<T>          value; // best function value for the rectangle size
  std::vector<std::shared_ptr<Rectangle<T>>> rects; // pointers to the rectangles
public:
  SizeMap() = default;
  SizeMap(const SizeMap& x) = default;
  auto operator=(const SizeMap& x) -> SizeMap& = default;  
  // inserts a new value into the class. It checks whether the size has previously been seen and then if so 
  // it changes the value associated with that size and the rectange.
  void insert(T fn_value, std::shared_ptr<Rectangle<T> > rect);
  void erase(size_t index);
  void clear();
  // finds the convex hull of points and removes any not on that hull. 
  // it starts at the position (0, fmin- e|fmin|) and then slowly rotates a line    
  void rotate(const T fmin, const T epsilon);
};


class bad_fn : public std::exception {};
             
             // class for handling function binding
             template <typename T>
             class direct {
               using func = T(*)(const void*, const std::vector<T>&); 
             public:
               struct DirectControl {
                 T epsilon = 1e-4;
                 int max_iter = 1000;
                 T tol = 1e-4;
                 bool adaptive = false; // use adaptive epsilon setting
               } control;
               
               direct(){};
               // if starting vals is true, it assumes x is a vector of starting values, and y sets the bounds by adding +/- either side of x
               // otherwise it assumes x is a lower bound, and y is an upper bound
               direct(const std::vector<T>& x, const std::vector<T>& y, bool starting_vals = true);
               direct(const direct& x) = default;
               auto operator=(const direct& x) -> direct& = default;
               
               // functions
               template<auto Function, typename = std::enable_if_t<std::is_invocable_r_v<T, decltype(Function),std::vector<T> > > >
               void    fn();
               template<auto Function, typename Class, typename = std::enable_if_t<std::is_invocable_r_v<T, decltype(Function),Class*,std::vector<T> > > >
               void    fn(Class* cls);
               
               void    set_bounds(const std::vector<T>& lower, const std::vector<T>& upper, bool starting_vals = true);
               auto    operator()(const std::vector<T>& vec) const -> T;
               void    optim();
               std::vector<T> values() const;
               
             private:
               [[noreturn]]
               static auto null_fn(const void* p, const std::vector<T> vec) -> T
               {
                 throw bad_fn{};
               }
               
               const void*             optim_instance = nullptr;  // pointer to the class if a member function
               func                    optim_fn = &null_fn;        // pointer to the function
               const size_t            dim;                       // number of dimensions
               std::vector<T>          lower_bound;               // bounds
               std::vector<T>          upper_bound;   
               std::vector<T>          dim_size;                  // size of each dimension to transform to unit rectangle
               std::vector<std::shared_ptr<Rectangle<T>>>  rects;   // the rectangles
               SizeMap<T>              size_fn;                   // a map that keeps track of the best function values for each rectangle size
               T                       min_f;                     // current best value
               int                     fn_counter = 0;
               int                     iter = 0;
               std::vector<T>          current_values;
               
               
               //functions
               auto            eval(const std::vector<T>& vec) -> T;
               std::vector<T>  transform(const std::vector<T>& vec);
               void            update_map();
               void            filter_rectangles();
               void            divide_rectangles();
               
             };
             
             template <typename T>
             inline std::vector<T> Rectangle<T>::centroid(){
               std::vector<T> centre(dim);
               for(size_t i = 0; i < dim; i++){
                 centre[i] = 0.5*(max_x[i] - min_x[i]);
               }
               return(centre);
             };
             
             template <typename T>
             inline void Rectangle<T>::unit_hyperrectangle(){
               std::fill(max_x.begin(),max_x.end(),1.0);
               std::fill(min_x.begin(),min_x.end(),0.0);
             };
             
             template <typename T>
             inline std::vector<T> Rectangle<T>::centroid(const size_t& dim, const T& delta){
               std::vector<T> centre(dim);
               for(size_t i = 0; i < dim; i++){
                 centre[i] = 0.5*(max_x[i] - min_x[i]);
                 if(i == dim) centre[i] += delta;
               }
               return(centre);
             };
             
             template <typename T>
             inline T Rectangle<T>::longest_side(){
               T long_len = 0;
               for(int i = 0; i < dim; i++){
                 T diff = max_x[i] - min_x[i];
                 if(diff > long_len) long_len = diff;
               }
               return 0.5*long_len;
             };
             
             template <typename T>
             inline T Rectangle<T>::dim_size(const size_t& dim){
               return max_x[dim] - min_x[dim];
             };
             
             template <typename T>
             inline void Rectangle<T>::trim_dimension(const size_t& dim, const Position pos){
               T dsize = dim_size(dim);
               dsize *= (T)1.0/3.0;
               switch(pos){
               case Position::Lower:
                 max_x[dim] -= 2*dsize;
                 break;
               case Position::Middle:
                 min_x[dim] += dsize;
                 max_x[dim] -= dsize;
                 break;
               case Position::Upper:
                 min_x[dim] += 2*dsize;
                 break;
               }        
             };
             
             template <typename T>
             inline void SizeMap<T>::insert(T fn_value, std::shared_ptr<Rectangle<T>> rect){
               T size_rect = rect->longest_side();
               const auto result = size.insert(size_rect);
               size_t insert_index = result.first - size.begin();
               if(result.second){            
                 value.insert(value.begin() + insert_index, fn_value);
                 rects.insert(rects.begin() + insert_index, rect);
                 rect->potentially_optimal = true;
               } else {
                 if(*result.first > fn_value){
                   value[insert_index] = fn_value;
                   rects[insert_index]->potentially_optimal = false;
                   rects[insert_index] = rect;
                   rect->potentially_optimal = true;
                 }
               }
               rect->no_update = true;
             };
             
             template <typename T>
             inline void SizeMap<T>::erase(size_t index){
               size.erase(size.begin()+index);
               value.erase(value.begin()+index);
               rects.erase(rects.begin()+index);
             };
             
             template <typename T>
             inline void SizeMap<T>::clear(){
               size.clear();
               value.clear();
               rects.clear();
             }
             
             template <typename T>
             inline void SizeMap<T>::rotate(const T fmin, const T epsilon){
               std::pair<T,T>  coord = {0, fmin - epsilon*abs(fmin)};
               size_t          index = 0;
               size_t          n_elem = size.size() - 1;
               T               angle = M_PI*0.5;
               T               x, y, new_angle;
               size_t          keep;
               while(index <= n_elem){
                 for(size_t i = index; i <= n_elem; i++){
                   y = value[i] - coord.second;
                   x = size[i] - coord.first;
                   new_angle = tan(y/x);
                   if(new_angle <= angle){
                     keep = i;
                     angle = new_angle;
                   }
                 }
                 coord.first = size[keep];
                 coord.second = value[keep];
                 if(keep != index){
                   for(;keep > index;keep--) {
                     erase(keep-1);
                   }
                 }
                 n_elem = size.size() -1;
                 index++;
               }
             }
             
             template <typename T>
             inline direct<T>::direct(const std::vector<T>& x, const std::vector<T>& y, bool starting_vals) {
               set_bounds(x,y,starting_vals);
             };
             
             template <typename T>
             inline std::vector<T> direct<T>::values() const {
               return current_values;
             }
             
             template <typename T>
             inline void direct<T>::set_bounds(const std::vector<T>& x, const std::vector<T>& y, bool starting_vals) {
               dim = x.size();
               if(starting_vals)
               {
                 lower_bound = x; 
                 upper_bound = y;
                 for(size_t i = 0; i < dim; i++) dim_size[i] = y[i] - x[i];
               } else {
                 lower_bound.resize(dim);
                 upper_bound.resize(dim);
                 for(size_t i = 0; i < dim; i++){
                   lower_bound[i] = x[i] - y[i];
                   upper_bound[i] = x[i] + y[i];
                   dim_size[i] = 2*y[i];
                 } 
               }
               current_values.resize(dim);
               rects.push_back(std::shared_ptr<Rectangle<T>>(new Rectangle<T>(dim)));
               rects.back()->unit_hyperrectangle();
             };
             
             template <typename T>
             template <auto Function, typename>
             inline void direct<T>::fn() 
             {
               optim_instance = nullptr;
               optim_fn = static_cast<func>([](const void*, const std::vector<T>& vec) -> T {
                 return std::invoke(Function, vec);
               });
             };
             
             template <typename T>
             template<auto Function, typename Class, typename>
             inline void direct<T>::fn(Class* cls)
             {
               optim_instance = cls;
               optim_fn = static_cast<func>([](const void* p, const std::vector<T>& vec) -> T {
                 auto* c = const_cast<Class*>(static_cast<const Class*>(p));
                 return std::invoke(Function,c,vec);
               });
             }
             
             template <typename T>
             inline auto direct<T>::operator()(const std::vector<T>& vec) const -> T
             {
               return std::invoke(optim_fn,optim_instance,vec);
             }  
             
             template <typename T>
             inline void direct<T>::optim(){
#ifdef R_BUILD
               Rcpp::Rcout << "\nSTARTING DIRECT-L";
               Rcpp::Rcout << "\nTolerance: " << control.tol << " | Max iter : " << control.max_iter << "\n Starting values :";
               std::vector<T> vals = transform(rects[0].centroid());
               for(const auto& val: vals) Rcpp::Rcout << val << " ";    
#endif
               T max_diff = 1.0;
               iter = 0;
               fn_counter = 0;
               while(max_diff > control.tol && iter <= control.max_iter){
#ifdef R_BUILD
                 Rcpp::Rcout << "\n----------------------------------\n";
                 Rcpp::Rcout << "\nIter: " << iter << " | Evaluations: " << fn_counter << " | Rectangles: " << rects.size() << " | Dimensions: " << dim;
#endif
                 update_map();
                 filter_rectangles();
                 auto new_values = transform(rects[0].centroid());
                 for(size_t i = 0; i < dim; i++){
                   T diff = abs(new_values[i] - current_values[i]);
                   if(diff > max_diff) max_diff = diff;
                   current_values[i] - new_values[i];
                 }
                 if(iter < control.max_iter)divide_rectangles();
#ifdef R_BUILD
                 Rcpp::Rcout << "\nNew value: " << min_f << " | Max difference: " << max_diff << " | Min. rectangle size: " << *(size_fn.size.begin()) << " | New values:\n";
                 for(const auto& val: current_values)Rcpp::Rcout << val << " ";
#endif
                 iter++;
               }
             }
             
             template <typename T>
             inline auto direct<T>::eval(const std::vector<T>& vec) -> T
             {   
               fn_counter++;
               return std::invoke(optim_fn,optim_instance,vec);
             }  
             
             template <typename T>
             inline std::vector<T> direct<T>::transform(const std::vector<T>& vec)
             {
               std::vector<T> transformed_vec(dim);
               for(int i = 0; i < dim; i++){
                 transformed_vec[i] = vec[i]*dim_size[i] + lower_bound[i];
               }
               return transformed_vec;
             };
             
             // after running this function size_fn should contain the potentially optimal rectangles
             template <typename T>
             inline void direct<T>::update_map()
             {
               for(auto& rect: rects){
                 if(!rect->no_update){
                   T fn_value = eval(transform(rect->centroid()));
                   size_fn.insert(fn_value, rect);
                   if(fn_value < min_f)min_f = fn_value;
                 }
               }
               size_fn.rotate(min_f, control.epsilon);
             };
             
             
             // if select_one = true then it will randomly select one of the potentially optimal squares
             // otherwise this function currently does nothing
             template <typename T>
             inline void direct<T>::filter_rectangles()
             {
               size_t n_rect = size_fn.size.size();
               size_t keep_index = random_index(n_rect-1);
               int counter = 0;
               for(auto& r: size_fn.rects){
                 if(counter != keep_index) r->potentially_optimal = false;
               }
             }
             
             // divides up the potentially optimal rectangles 
             // it adds the rectangles to rects, and then erases the original 
             // rectangles, while also clearing size_fn
             template <typename T>
             inline void direct<T>::divide_rectangles(){
               //identify the largest dimensions
               typedef std::pair<T,size_t> dimpair;
               
               struct ComparePair {
                 bool operator()(const dimpair& elt1, const dimpair& elt2) const {
                   return elt1.first < elt2.first;
                 }
               };
               
               std::vector<size_t> largest_dims;
               std::priority_queue< dimpair, std::vector<dimpair>, ComparePair > pq;
               T dim_size;
               for(const auto r: size_fn.rects){
                 if(r->potentially_optimal){
                   largest_dims.clear();
                   dim_size = 0; 
                   for(size_t i = 0; i < dim; i++){
                     T dim_size_i = r->dim_size(i);
                     if(dim_size_i == dim_size){
                       largest_dims.push_back(i);
                     } else if(dim_size_i > dim_size){
                       largest_dims.clear();
                       largest_dims.push_back(i);
                       dim_size = dim_size_i;
                     }
                   }
                   
                   T delta = dim_size / 3;
                   T fn1, fn2;
                   for(const auto& d: largest_dims){
                     fn1 = eval(transform(r->centroid(d, delta)));
                     fn2 = eval(transform(r->centroid(d, -delta)));
                     pq.push(dimpair(std::min(fn1,fn2) , d));
                   }
                   
                   while(!pq.empty){
                     size_t dim = pq.top().second;
                     rects.push_back(std::shared_ptr<Rectangle<T>>(new Rectangle(*r)));
                     rects.back()->trim_dimension(dim,Position::Upper);
                     rects.push_back(std::shared_ptr<Rectangle<T>>(new Rectangle(*r)));
                     rects.back()->trim_dimension(dim,Position::Lower);
                     r->trim_dimension(dim, Position::Middle);
                     r->no_update = false;
                     pq.pop();
                   }
                 }
                 
               }
               
             };