# Numerical Methods - Day 1 Solutions

## Task 3.2 - Second Derivative

### Taylor expansion

For a sufficiently smooth function f:

f(x+h) = f(x) + hf'(x) + h²f''(x)/2
       + h³f'''(x)/6 + h⁴f''''(x)/24 + ...

f(x-h) = f(x) - hf'(x) + h²f''(x)/2
       - h³f'''(x)/6 + h⁴f''''(x)/24 + ...

Adding and subtracting 2f(x):

f(x+h) - 2f(x) + f(x-h)
= h²f''(x) + h⁴f''''(x)/12 + O(h⁶)

Therefore:

D2f(x) = [f(x+h)-2f(x)+f(x-h)]/h²

D2f(x) = f''(x) + h²f''''(x)/12 + O(h⁴)

The approximation is second-order accurate.

### Numerical verification

For f(x) = sin(x), x = 1:

Exact second derivative = -sin(1).

| h | Absolute error |
|---|---|
| 0.1 | 0.0007009921 |
| 0.05 | 0.0001752918 |
| 0.025 | 0.0000438257 |
| 0.0125 | 0.0000109566 |

The error decreases by approximately a factor
of four when h is halved.

Conclusion: observed order of accuracy is 2.

Experiment:
scripts/task_3_2.py

Results:
results/task_3_2.txt


## Task 5.2 - Two-Dimensional Update Rule

### Heat equation

du/dt = alpha * (d²u/dx² + d²u/dy²)

Use Forward Euler for the time derivative and
centered second differences in both directions.

The discrete equation becomes:

u_new(i,j) = u_old(i,j)
  + alpha*dt/dx² * [u_old(i-1,j)
                    - 2*u_old(i,j)
                    + u_old(i+1,j)]
  + alpha*dt/dy² * [u_old(i,j-1)
                    - 2*u_old(i,j)
                    + u_old(i,j+1)]

For a uniform square grid, dx = dy = h.

Define r = alpha*dt/h².

Then:

u_new(i,j) = (1-4*r)*u_old(i,j)
  + r*(u_old(i-1,j) + u_old(i+1,j)
       + u_old(i,j-1) + u_old(i,j+1))

This is the five-point stencil implemented
in src/heat.cpp.

### Space-time dependencies

The value at (i,j) at time n+1 depends on:

- Center at time n
- Left neighbor at time n
- Right neighbor at time n
- Upper neighbor at time n
- Lower neighbor at time n

No diagonal neighbor is required.

### Computational cost

Per interior point, with 1-4*r precomputed:

- Five old temperature values are read.
- One new temperature value is written.
- Approximately four additions.
- Two multiplications.

Conclusion:

The method uses an explicit five-point stencil
with two separate temperature arrays.
