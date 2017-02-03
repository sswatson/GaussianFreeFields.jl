
`GaussianFreeFields` is a package that generates a discrete Gaussian free
on a square with either periodic boundary conditions or specified Dirichlet
boundary conditions.

To install, run
`Pkg.clone("git://github.com/sswatson/GaussianFreeFields.jl.git")`

`DGFF(n)` generates a discrete Gaussian free field on an `n` by `n` torus. 

```julia
h = DGFF(20)
```

```
20x20 Array{Float64,2}:
 -0.876092   -0.294548   -1.35524    …   0.0387392   0.466126   -0.648591 
  0.215345   -0.886166   -1.47335        0.612032   -0.133992    0.0649898
 -0.159199    0.596377   -0.569785       0.620096   -0.218865    1.54292 
 -0.301339    0.254145   -0.501166       0.185516   -0.265336    0.209947 
 -0.480379   -0.398221   -0.566758      -0.248938    0.772152    0.441475 
 -0.0974746  -0.086842    0.443506   …   0.395192    1.11058     0.61846 
  0.9968      0.423853    0.682701       1.47028     0.751546    0.939246 
 -0.413873   -0.044707    0.555525       1.43303     1.18651     0.748613 
 -1.22688    -0.282943   -0.418801       0.246069    0.650025    0.373996 
 -0.119498   -0.123819   -0.518486      -0.425611    0.0454104   0.326286 
  0.941938    0.495544   -0.395305   …   0.556364    0.984167    1.48036 
  0.499854   -0.229561    0.388527       0.258685    0.8598      1.47906 
 -0.764796   -0.885678   -0.439279       0.302158    0.654056    0.672523 
  0.0673879  -0.868744   -0.65548       -0.251016    0.0186039  -0.72188 
 -1.39916    -1.53116    -0.671419      -0.382203   -0.877608   -1.58912 
 -0.932672   -0.988617   -0.0599065  …  -1.1248     -0.381567   -1.26558 
 -0.1239     -0.415419   -1.15146       -1.22362    -0.706317   -0.399698 
  0.0324895   0.469533   -0.324562       0.582178    0.145177    0.366976 
 -1.04055     0.0689947  -0.360072       0.361529   -0.101115    0.442297 
 -0.78754    -0.227399   -0.448399      -0.0482068   0.480722   -0.218205 
```

The harmonic extension can be subtracted off using `fix_boundary_values`

```julia
h0 = fix_boundary_values(h,zeros(20,20))
```

```
20x20 Array{Float64,2}:
 0.0   0.0         0.0          0.0        …   0.0         0.0        0.0
 0.0  -0.582553   -0.689298    -0.795452       0.380308   -0.488851   0.0
 0.0   0.947578    0.00399912  -0.718617       0.284179   -0.875459   0.0
 0.0   0.622354   -0.0301214    0.338422      -0.112048   -0.658012   0.0
 0.0  -0.0489716  -0.177794     1.02856       -0.553679    0.365555   0.0
 0.0   0.0726041   0.724006     1.1459     …   0.0536066   0.623081   0.0
 0.0   0.334414    0.883015     0.483024       1.09396     0.168205   0.0
 0.0   0.234573    0.842797    -0.449565       1.06916     0.656206   0.0
 0.0   0.222473   -0.0779471   -0.394762      -0.08763     0.224613   0.0
 0.0   0.0508313  -0.280688     0.456136      -0.769926   -0.418233   0.0
 0.0   0.331433   -0.259279     0.453971   …   0.164604    0.225609   0.0
 0.0  -0.254745    0.583085    -0.350277      -0.0425239   0.161329   0.0
 0.0  -0.517007   -0.0853753    0.242049       0.259535    0.399002   0.0
 0.0  -0.487572   -0.200811    -0.2869         0.0403009   0.412004   0.0
 0.0  -0.762431   -0.0690003   -0.770596       0.141392   -0.0621508  0.0
 0.0  -0.296438    0.560409    -0.0665069  …  -0.578232    0.374142   0.0
 0.0   0.0315762  -0.576131    -1.15955       -0.812917   -0.311079   0.0
 0.0   0.866107    0.266444    -0.161781       0.814685    0.160018   0.0
 0.0   0.649782    0.298557    -0.264757       0.464097   -0.302518   0.0
 0.0   0.0         0.0          0.0            0.0         0.0        0.0
```

We can use `PyPlot` to see a surface plot of the zero-boundary GFF `h0`. 

```julia
using PyPlot
X = [x for x=1:n,y=1:n]
Y = [y for x=1:n,y=1:n]
plot_surface(X, Y, h0, rstride = 1, cstride = 1, cmap="autumn")
```

![A zero-boundary Gaussian free field](https://github.com/sswatson/GaussianFreeFields.jl/blob/master/images/gff.png)

With `Contour`, we can look at the level lines of the Gaussian free field.

```julia
using Graphics2D 
using Contour 
using GaussianFreeFields 
n = 250 
h = DGFF(n) 
h0 = fix_boundary_values(h) 
showgraphics([Line(c; linewidth=1.0, color=0.3*"green") 
     for c in contour(collect(1.0:n),collect(1.0:n),h0,0.0).lines]; dim = 1024)```

![GFF level lines](https://github.com/sswatson/GaussianFreeFields.jl/blob/master/images/gfflevellines.png)

We can also calculate the [GFF flow lines](http://arxiv.org/abs/1201.1496)

```julia
using Graphics2D
using GaussianFreeFields
using Grid
n = 250
h = InterpGrid(fix_boundary_values(DGFF(n)),BCnil,InterpLinear);
κ = 3
χ = 2/sqrt(κ) - sqrt(κ)/2
z0 = (n+1)/2 + im*(n+1)/2
δ = 0.01;
fan = GraphicElement[]
for θ=0.0:0.05:2π
	push!(fan,Line(flowline(h, z0, χ, θ);
	               color=θ/(2π)*"green" + (1-θ/(2π))*"blue"))
end

showgraphics([Line([1 1; n 1; n n; 1 n; 1 1]);fan])
```

![SLE fan](https://github.com/sswatson/GaussianFreeFields.jl/blob/master/images/SLEfan.png)

The Julia function `writecsv` can be used to export the data:

```julia
writecsv("h.csv",h0)
```


[![Build Status](https://travis-ci.org/sswatson/GaussianFreeFields.jl.svg?branch=master)](https://travis-ci.org/sswatson/GaussianFreeFields.jl)
