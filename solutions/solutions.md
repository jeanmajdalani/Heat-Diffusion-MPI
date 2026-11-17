# Solutions - Numerical Methods for the Heat Equation

## Task 3.2 - Second Derivative

### Taylor expansion

For a sufficiently smooth function f:

f(x+h) = f(x) + h f'(x) + h² f''(x)/2
       + h³ f'''(x)/6 + h⁴ f''''(x)/24 + ...

f(x-h) = f(x) - h f'(x) + h² f''(x)/2
       - h³ f'''(x)/6 + h⁴ f''''(x)/24 + ...

Adding the two expansions and subtracting 2f(x):

f(x+h) - 2f(x) + f(x-h)
= h² f''(x) + h⁴ f''''(x)/12 + O(h⁶)

Therefore:

D2 f(x) = [f(x+h) - 2f(x) + f(x-h)] / h²

D2 f(x) = f''(x) + h² f''''(x)/12 + O(h⁴)

The truncation error is h² f''''(x)/12 + O(h⁴), so the
approximation is second-order accurate.

### Numerical verification

For f(x) = sin(x), x = 1: the exact second derivative is -sin(1).

| h | Absolute error |
|---|---|
| 0.1 | 0.0007009921 |
| 0.05 | 0.0001752918 |
| 0.025 | 0.0000438257 |
| 0.0125 | 0.0000109566 |

The error decreases by a factor of about four when h is halved.
Observed order of accuracy: 2.

Experiment: scripts/task_3_2.py
Results: results/task_3_2.txt


## Task 5.2 - Two-Dimensional Update Rule

### Heat equation

du/dt = alpha * (d²u/dx² + d²u/dy²)

Forward Euler in time, centered second differences in space.

### (a) General update, dx and dy different

With r_x = alpha*dt/dx² and r_y = alpha*dt/dy²:

u_new(i,j) = (1 - 2*r_x - 2*r_y) * u_old(i,j)
           + r_x * (u_old(i-1,j) + u_old(i+1,j))
           + r_y * (u_old(i,j-1) + u_old(i,j+1))

### (b) Uniform grid, dx = dy = h

Define r = alpha*dt/h². Then:

u_new(i,j) = (1 - 4r) * u_old(i,j)
           + r * (u_old(i-1,j) + u_old(i+1,j)
                  + u_old(i,j-1) + u_old(i,j+1))

This is the five-point stencil implemented in src/heat.cpp.

### (c) Space-time stencil

Level n:                          Level n+1:

            (i,j+1)
               |
  (i-1,j) -- (i,j) -- (i+1,j)  ----->   (i,j)
               |
            (i,j-1)

The value at (i,j) at time n+1 depends on the center and its four
neighbors (left, right, up, down) at time n. No diagonal neighbor
is needed.

### (d) Cost of one update

With (1 - 4r) precomputed:

- Reads: 5 old values.
- Write: 1 new value.
- Floating point operations: 4 additions and 2 multiplications,
  that is 6 operations.

Conclusion: an explicit five-point stencil with two separate arrays
(old and new), swapped after each step.