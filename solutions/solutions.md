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


## Task 3.4 - The Discrete Laplacian

(a) d²u/dx² at (x_i, y_j):
    [u(i+1,j) - 2u(i,j) + u(i-1,j)] / dx²

(b) d²u/dy² at (x_i, y_j):
    [u(i,j+1) - 2u(i,j) + u(i,j-1)] / dy²

(c) Adding both:

    Lap_h u(i,j) = [u(i+1,j) - 2u(i,j) + u(i-1,j)] / dx²
                 + [u(i,j+1) - 2u(i,j) + u(i,j-1)] / dy²

    Using Task 3.2, the truncation error is
    (dx²/12) u_xxxx + (dy²/12) u_yyyy + O(dx⁴ + dy⁴),
    that is O(dx² + dy²): second-order accurate.

(d) For dx = dy = h:

    Lap_h u(i,j) = [u(i+1,j) + u(i-1,j) + u(i,j+1) + u(i,j-1)
                    - 4u(i,j)] / h²

(e) Five-point stencil (weights multiplied by 1/h²):

                 (i,j+1)  weight 1
                    |
    (i-1,j) -- (i,j) -- (i+1,j)
     weight 1   weight -4   weight 1
                    |
                 (i,j-1)  weight 1

    The weights add up to 0.

(f) Lap_h u(i,j) = (4/h²) * (average of the 4 neighbors - u(i,j)).
    The mean-value formula (1.2) says
    average over a circle of radius r = u + (r²/4) Lap u,
    so Lap u ≈ (4/r²)(average - u). With r = h and the four
    neighbors as a discrete circle, both formulas give the same
    result: the Laplacian measures how much the average of the
    surroundings exceeds the value at the point.


## Task 5.5 - Accuracy of the Scheme

The scheme is
(u^{n+1}_{ij} - u^n_{ij}) / dt = Lap_h u^n_{ij}.

- Time: forward difference. By Example 3.3 its truncation error is
  (dt/2) u_tt + O(dt²), so p = 1.
- Space: central differences (Task 3.4), error O(dx²) + O(dy²),
  so q = 2.

Therefore tau = O(dt) + O(dx²) + O(dy²), with p = 1 and q = 2.
tau -> 0 as dt, dx, dy -> 0, so the scheme is consistent.


## Task 6.3 - Von Neumann Analysis in 2D

Substitute u^n_{jk} = G^n e^{i(j*theta + k*phi)} into
u^{n+1} = u^n + r_x (delta_x² u^n) + r_y (delta_y² u^n),
with r_x = alpha*dt/dx² and r_y = alpha*dt/dy².

(a) Using e^{i*t} - 2 + e^{-i*t} = -4 sin²(t/2):

    G(theta, phi) = 1 - 4 r_x sin²(theta/2) - 4 r_y sin²(phi/2)

    G is real.

(b) G <= 1 always. The smallest G is at theta = phi = pi
    (checkerboard pattern): G = 1 - 4 r_x - 4 r_y.
    Stability needs G >= -1, that is
    r_x + r_y <= 1/2, or
    dt <= 1 / (2 alpha (1/dx² + 1/dy²)).

(c) For dx = dy = h: r_x = r_y = r, so 2r <= 1/2:

    r <= 1/4

(d) In 1D the bound is r <= 1/2; in 2D it is r <= 1/4.
    Guess for 3D: r <= 1/6 (in d dimensions, r <= 1/(2d)).


## Task 7.1 - Predicting the Order

With r fixed, dt = r*h², so the time error O(dt) = O(h²) and the
space error is O(h²). The overall order in h is 2.


## Task 7.2 - Convergence Study

Verification problem, T = 0.05, r = 0.2
(script: scripts/convergence.py, data: results/convergence.csv,
plot: results/convergence.png).

| N | h | steps | err_inf | err_2,h | p_obs |
|---|---|---|---|---|---|
| 16 | 0.0625 | 64 | 1.6634e-3 | 8.3169e-4 | - |
| 32 | 0.03125 | 256 | 4.1418e-4 | 2.0709e-4 | 2.0058 |
| 64 | 0.015625 | 1024 | 1.0344e-4 | 5.1721e-5 | 2.0014 |
| 128 | 0.0078125 | 4096 | 2.5854e-5 | 1.2927e-5 | 2.0004 |
| 256 | 0.00390625 | 16384 | 6.4631e-6 | 3.2316e-6 | 2.0001 |

(c) The observed orders tend to 2, as predicted in Task 7.1.
Note that err_2,h is about err_inf / 2 in every row: the error is
proportional to sin(pi x) sin(pi y), whose discrete 2-norm is about
1/2, so both norms give the same order.


## Task 8.2 - Surface-to-Volume

N x N grid, P processes.

(a) Strips: each process owns N/P rows. A typical (interior) process
    sends 2 rows and receives 2 rows per step: 2N numbers sent and
    2N received, in 2 messages each way. This does not decrease as
    P grows.

(b) Blocks: P = sqrt(P) x sqrt(P) blocks of (N/sqrt(P))² nodes.
    A typical process sends 4 edges of N/sqrt(P) numbers:
    4N/sqrt(P) numbers sent and as many received, in 4 messages
    each way.

(c) With t_msg = t_s + m*t_w, counting messages one after another:

    T_strips = 2 (t_s + N t_w)
    T_blocks = 4 (t_s + (N/sqrt(P)) t_w)

    Blocks are faster when T_blocks < T_strips, that is when
    t_s < N t_w (1 - 2/sqrt(P)).
    - When t_s is small (bandwidth-dominated), blocks win for P > 4
      (less data: 4N/sqrt(P) < 2N). For P = 4 the volumes are equal.
    - When t_s is very large (latency-dominated), strips win:
      2 messages instead of 4.

(d) Computation per step: about t_c * N²/P per process. For strips:

    T_comm / T_comp = 2 (t_s + N t_w) P / (t_c N²)
                    = 2 P t_s / (t_c N²) + 2 P t_w / (t_c N)

    As N grows with P fixed, the ratio tends to 0 (like 1/N for the
    bandwidth term and 1/N² for the latency term). Over a slow
    network (large t_s), parallelism pays off only for large grids,
    when t_c N²/P is much larger than t_s.
    

## Task 3.1 - The Backward Difference

Define D-f(x) = [f(x) - f(x-h)] / h.

By Taylor's theorem with remainder,

    f(x-h) = f(x) - h f'(x) + (h²/2) f''(ξ),   ξ between x-h and x.

Therefore f(x) - f(x-h) = h f'(x) - (h²/2) f''(ξ), and

    D-f(x) = f'(x) - (h/2) f''(ξ).

Truncation error:

    τ(h) = D-f(x) - f'(x) = -(h/2) f''(ξ) = -(h/2) f''(x) + O(h²).

Since f'' is bounded near x, τ(h) = O(h): the backward difference is
first-order accurate. Its leading error has the opposite sign of the
forward difference error (+(h/2) f''), but the same size.

Numerical check for f(x) = sin x at x = 1, f'(1) = cos 1:

| h | |D-f - f'| | ratio | leading term (h/2)|f''(1)| |
|---|---|---|---|
| 0.1 | 4.114e-2 | - | 4.207e-2 |
| 0.05 | 2.081e-2 | 1.98 | 2.104e-2 |
| 0.025 | 1.046e-2 | 1.99 | 1.052e-2 |
| 0.0125 | 5.245e-3 | 1.99 | 5.259e-3 |

The ratios tend to 2 (observed order log2(ratio) ≈ 1), and the error
matches the leading term (h/2)|f''(1)| = (h/2) sin 1.


## Task 3.3 - Sanity Checks on the Formula

D2 f(x) = [f(x-h) - 2f(x) + f(x+h)] / h²   (from Task 3.2).

(a) Exact values:

- f(x) = 1:   (1 - 2 + 1) / h² = 0.
- f(x) = x:   [(x-h) - 2x + (x+h)] / h² = 0.
- f(x) = x²:  [(x-h)² - 2x² + (x+h)²] / h²
              = [x² - 2xh + h² - 2x² + x² + 2xh + h²] / h²
              = 2h² / h² = 2.

These are exactly f'' for these three functions (0, 0 and 2).
A reasonable approximation of f'' must pass these tests: they are the
simplest functions whose second derivative is known, and the
truncation error of D2 contains only f'''' and higher derivatives,
which vanish for polynomials of degree at most 3. If a formula failed
on them, its error would not tend to 0 as h -> 0 even for trivial
functions, so it would not be consistent. Also, f'' must not change
when a constant or a linear function is added to f.

(b) The coefficients of f(x-h), f(x), f(x+h) are 1/h², -2/h², 1/h².
Their sum is 1 - 2 + 1 = 0. This is the test f = 1 seen in the
coefficients: a constant function has f'' = 0, so the weights must
add up to zero. The formula only measures differences between values,
never their level.

(c) Rewrite

    D2 f(x) = (2/h²) * { [f(x-h) + f(x+h)]/2 - f(x) }.

So D2 f(x) is (2/h²) times the amount by which the average of the two
neighbours exceeds the value at the point. It is positive when f(x) is
below the average of its neighbours (f curves upward) and negative
when it is above. This is the one-dimensional version of the meaning
of the Laplacian in Section 1.4: a point heats up when its
surroundings are on average warmer than it.