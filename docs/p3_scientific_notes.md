# Scientific modelling and numerical verification

## 1. Physical model

The project solves the two-dimensional heat equation
on the unit square [0,1] x [0,1], with thermal
diffusivity alpha = 1.

du/dt = alpha * (d²u/dx² + d²u/dy²)

Two physical configurations are considered:

- Case 1: Initially cold interior, with all four
  edges maintained at temperature 1. The plate
  approaches a uniform temperature of 1.

- Case 2: Initially cold interior, with the top
  and bottom edges maintained at temperature 1,
  including the corners, and the remaining left
  and right edge nodes maintained at temperature 0.

A third configuration, Case V, is used for numerical
verification because its exact analytical solution
is known.

u(x,y,t) = sin(pi*x)*sin(pi*y)*exp(-2*pi²*t)

Its boundary temperature is zero on all four edges.

## 2. Numerical discretization

The domain is divided into N intervals in each
direction, giving (N+1)² grid points and spacing:

h = 1/N

The numerical method combines forward Euler time
integration with centered second-order finite
differences in space.

The resulting explicit five-point stencil is:

u_new(i,j) = (1-4r)*u_old(i,j)
           + r*(u_old(i-1,j) + u_old(i+1,j)
              + u_old(i,j-1) + u_old(i,j+1))

where r = alpha*dt/h².

Two separate arrays store consecutive time levels.
This preserves the original explicit scheme and
allows independent updates of the interior nodes.

## 3. Stability and theoretical accuracy

Von Neumann analysis and the discrete maximum
principle both give the stability condition:

0 <= r <= 0.25

Equivalently:

dt <= h²/(4*alpha)

The method is first-order accurate in time and
second-order accurate in space:

Error = O(dt + h²)

When r is fixed, dt is proportional to h².
Consequently, the expected overall convergence
order with respect to h is two.

## 4. Experimental convergence verification

Case V was tested at final time T = 0.05,
with r = 0.2 and grid resolutions from 16 to 256.

The observed convergence order is calculated as:

p = log2(E(h)/E(h/2))

The recorded results were independently checked
using Python and results/convergence.csv.

| N | Infinity error | Discrete L2 error | Observed order |
|---:|---:|---:|---:|
| 16 | 1.663371e-3 | 8.316853e-4 | - |
| 32 | 4.141824e-4 | 2.070912e-4 | 2.005771 |
| 64 | 1.034425e-4 | 5.172124e-5 | 2.001437 |
| 128 | 2.585419e-5 | 1.292709e-5 | 2.000359 |
| 256 | 6.463145e-6 | 3.231572e-6 | 2.000090 |

Both error norms give the same observed orders
to the precision shown.

Doubling N reduces the error by approximately
a factor of four, confirming the predicted
second-order convergence.

The results provide numerical verification
of the implemented scheme for Case V.

## 5. Motivation for parallelization

At each time step, approximately N² interior
grid points must be updated.

For a fixed physical simulation time, the stability
condition requires approximately O(N²) time steps.

The resulting computational complexity is O(N^4).

Doubling the resolution therefore increases the
computational work by approximately sixteen times.

MPI domain decomposition distributes this work
across multiple processes, using halo exchanges
to communicate neighboring temperatures.
