# Pspline

Welcome to your new gem! In this directory, you'll find the files you need to be able to package up your Ruby library into a gem. Put your Ruby code in the file `lib/pspline`. To experiment with that code, run `bin/console` for an interactive prompt.

## Installation

Add this line to your application's Gemfile:

```ruby
gem 'pspline'
```

And then execute:

    $ bundle

Or install it yourself as:

    $ gem install pspline

## Usage

# 1. abstruct

Pspline interpolation libraly implemented in C-Extensions.

# 2. usage in ruby script file
require 'pspline'   
include PSPLINE   
   
obj = Bspline.new(...)   
...or   
obj = Pspline.new(...)   
...

# 3. class
Module 'PSPLINE' has 2 classes, Bspline and Pspline.
## 3.1 Basic interpolation without boundary condition
PSPLINE::Bspline.new([*plist*], *n*, *j*)   
- *plist*   
[*x1,y1*],...,[*xn,yn*]   
- *n*   
number of data points   
- *j*   
degree of the polynominal   

## 3.2 Interpolation with boundary condition by additional data points
PSPLINE::Bspline.new([*plist*,*clist*], *n*, *j*)   
- *plist*   
[*x1,y1*],...,[*xn,yn*]   
- *clist*   
[*xn+1,yn+1*],...,[*xn+j-1,yn+j-1*]   
- *n*   
number of data points   
- *j*   
degree of the polynominal   

## 3.3 Interpolation with boundary condition by differential value
PSPLINE::Bspline.new([*plist*,*clist*], *n*, *j*, *dlist*)
- *plist*   
[*x1,y1*],...,[*xn,yn*]   
- *clist*   
[*xn+1,yn+1*],...,[*xn+j-1,yn+j-1*]   
- *n*   
number of data points   
- *j*   
degree of the polynominal   
- *dlist*   
[*d1*,...,*dj-1*]

## 3.4 Interpolation with period boundary condition
PSPLINE::Bspline.new([*plist*], *n*, *j*, 1)   
- *plist*   
[*x0,y0*],...,[*xn,yn*] (_yn_ = _y0_)   
- *n*   
number of data points   
- *j*   
degree of the polynominal   

## 3.5 One variable M dimension parametric interpolation
PSPLINE::Pspline.new([*ylist*],[*xlist*],[*n*],[*j*],[0])   
- *ylist*   
[*y1*,...,_ym_],...,[*yn*,...,_yn+m-1_]   
- *xlist*   
[_x1_,...,_xn_]   
- _n_   
number of data points   
- _j_   
degree of the polynominal   

PSPLINE::Bspline.new([_xylist_], _n_, _j_)   
- _xylist_   
[_x1_, _y1_,...,_ym_],...,[_xn_, _yn_,...,_yn+m-1_]   
- _n_   
number of data points   
- _j_   
degree of the polynominal   

## 3.6 One variable M dimension parametric period interpolation
PSPLINE::Pspline.new([_ylist_],[_xlist_],[_n_],[_j_],[1])   
- _ylist_
[_y0_,...,_ym-1_],..,[_yn-1_,...,_yn+m-2_]   
- _xlist_   
[_x0_,...,_xn_]   
- _n_   
number of data points   
- _j_   
degree of the polynominal   

PSPLINE::Bspline.new([_xylist_], _n_, _j_, 1)   
- _xylist_   
[_x0_, _y0_,...,_ym-1_],...,[_xn_, _yn_,...,_yn+m-1_] (_yn_ = _y0_,...,_yn+m-1_ = _ym-1_)   
- _n_   
number of data points   
- _j_   
degree of the polynominal   

## 3.7 One variable M dimension Riesenfeld interpolation
PSPLINE::Pspline.new([_ylist_],[_xlist_],[_n_],[_j_],[2])   
- _ylist_   
[_y1_,...,_ym_],...,[_yn_,...,_yn+m-1_]   
- _xlist_   
[_x1_,...,_xn_]   
- _n_   
number of data points   
- _j_   
degree of thepolynominal   

PSPLINE::Bspline.new([_xylist_], _n_, _j_, 2)   
- _xylist_   
[_x1_, _y1_,...,_ym_],...,[_xn_, _yn_,...,_yn+m-1_]   
- _n_   
number of data points   
- _j_   
degree of the polynominal   

## 3.8 One variable M dimension Riesenfeld period interpolation
PSPLINE::Pspline.new([_ylist_],[_xlist_],[_n_],[_j_],[3])   
- _ylist_   
[_y0_,...,_ym-1_],...,[_yn-1_,...,_yn+m-2_]   
- _xlist_   
[_x0_,...,_xn_]   
- _n_   
number of data points   
- _j_   
degree of the polynominal   

PSPLINE::Bspline.new([_xylist_], _n_, _j_, 3)   
- _xylist_   
[_x0_, _y0_,...,_ym-1_],..,[_xn_, _yn_,...,_yn+m-1_] (_yn_ = _y0_,...,_yn+m-1_ = _ym-1_)   
- _n_   
number of data points   
- _j_   
degree of the polynominal   

## 3.9 N variable M dimension parametric interpolation
PSPLINE::Pspline.new([_xyzlist_],[_stulist_],[_n1_,...,_nn_],[_j1_,...,_jn_],[0,...,0])   
- _xyzlist_   
[...[_x1_,_y1_,...,_z1_],...],... (_n1_~...~_nn_ multidimensional array)   
- _stulist_   
[_s1_,...,_sn1_],[_t1_,...,_tn2_],...,[_u1_,...,_unn_]   
- _n1_,...,_nn_   
numerical list of data points   
- _j1_,...,_jn_   
list of degrees of the polynominal   

## 3.10 N variable M dimension parametric period interpolation
PSPLINE::Pspline.new([_xyzlist_],[_stulist_],[_n1_,...,_nn_],[_j1_,...,_jn_],[1,...,0])   
- _xyzlist_   
[...[_x1_,_y1_,...,_z1_],...],... (_n1_~...~_nn_ multidimensional array)   
- _stulist_   
[_s0_,...,_sn1_],[_t1_,...,_tn2_],...,[_u1_,...,_unn_]  
- _n1_,...,_nn_   
numerical list of data points   
- _j1_,...,_jn_   
list of degrees of the polynominal   

## 3.11 N variable M dimension Riesenfeld interpolation
PSPLINE::Pspline.new([_xyzlist_],[_stulist_],[_n1_,...,_nn_],[_j1_,...,_jn_],[2,...,2])   
- _xyzlist_   
[...[_x1_,_y1_,...,_z1_],...],... (_n1_~...~_nn_ multidimensional array)   
- _stulist_   
[_s1_,...,_sn1_],[_t1_,...,_tn2_],...,[_u1_,...,_unn_]   
- _n1_,...,_nn_   
numerical list of data points   
- _j1_,...,_jn_   
list of degrees of the polynominal   

## 3.12 N variable M dimension Riesenfeld period interpolation
PSPLINE::Pspline.new([_xyzlist_],[_stulist_],[_n1_,...,_nn_],[_j1_,...,_jn_],[3,...,2])   
- _xyzlist_   
[...[_x1_,_y1_,...,_z1_],...],... (_n1_~...~_nn_ multidimensional array)   
- _stulist_   
[_s0_,...,_sn1_],[_t1_,...,_tn2_],...,[_u1_,...,_unn_]   
- _n1_,...,_nn_   
numerical list of data points   
- _j1_,...,_jn_   
list of degrees of the polynominal   

# 4. Calculate interporation
	self[x]	#=> y
	self[x1,...,xi] #=> [y1,...yi]
	
See example/*.rb

## Development

After checking out the repo, run `bin/setup` to install dependencies. You can also run `bin/console` for an interactive prompt that will allow you to experiment.

To install this gem onto your local machine, run `bundle exec rake install`. To release a new version, update the version number in `version.rb`, and then run `bundle exec rake release`, which will create a git tag for the version, push git commits and tags, and push the `.gem` file to [rubygems.org](https://rubygems.org).

## Contributing

Bug reports and pull requests are welcome on GitHub at https://github.com/vsp2old/pspline.

