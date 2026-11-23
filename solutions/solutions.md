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
## Task 1.1 - Mean-value formula
1.1 – Formule de la moyenne (1.2). Avec a = r cos φ et b = r sin φ, le développement de Taylor à deux variables donne u + a·uₓ + b·u_y + ½(a²uₓₓ + 2ab·uₓ_y + b²u_yy) + O(r³). On intègre terme à terme sur [0, 2π] :

Le terme constant donne 2π·u.
Les termes d’ordre 1 s’annulent, car ∫cos φ = ∫sin φ = 0.
Au second ordre, le terme croisé s’annule (∫cos φ sin φ = 0) et ∫cos² = ∫sin² = π. Il reste (r²/2)·π·(uₓₓ + u_yy) = (πr²/2)∇²u.
On divise par 2π et on obtient u + (r²/4)∇²u.
Les termes d’ordre 3 s’annulent aussi (symétrie φ → φ+π). Le reste est donc O(r⁴), pas O(r³).
## Task 1.2 - Dimensional analysis
1.2 – Analyse dimensionnelle.

(a) [αt/L²] = (m²/s · s)/m² = 1, donc le nombre est sans dimension.
(b) τ = α^p·L^q doit être en secondes. Cela impose p = −1 et q = 2, donc τ = L²/α, et c’est la seule combinaison possible. Tout temps caractéristique ne dépend que de α et de L. Il est donc proportionnel à τ, car le rapport t*/τ est un nombre sans dimension qui ne dépend de rien d’autre.
Cuivre, L = 0,1 m : τ = 0,01/1,1·10⁻⁴ ≈ 91 s.
Acier au carbone : τ = 0,01/1,3·10⁻⁵ ≈ 769 s ≈ 13 min.
Cela correspond à l’intuition : le cuivre égalise sa température environ 8,5 fois plus vite. τ n’est qu’un ordre de grandeur.
(c) Avec x̂ = x/L, ŷ = y/L et t̂ = αt/L², la règle de la chaîne donne ∂ₜ = (α/L²)∂ₜ̂ et ∂ₓₓ = (1/L²)∂ₓ̂ₓ̂. L’équation devient (α/L²)u_t̂ = (α/L²)(u_x̂x̂ + u_ŷŷ), soit u_t̂ = u_x̂x̂ + u_ŷŷ. Ni α ni L n’apparaissent. Un solveur écrit pour α = 1 sur [0,1]² couvre donc toutes les plaques et tous les matériaux. Il suffit de convertir le temps avec t̂ = αt/L².
## Task 1.3 - Fourier modes
1.3 – Modes de Fourier. On pose u = sin(mπx/Lₓ)·sin(nπy/L_y)·e^{−λt}, avec λ = απ²(m²/Lₓ² + n²/L_y²).

On a u_t = −λu.
On a uₓₓ = −(mπ/Lₓ)²u et u_yy = −(nπ/L_y)²u. Donc α(uₓₓ + u_yy) = −λu = u_t, et l’équation (1.1) est satisfaite.
Sur les bords, sin 0 = sin mπ = 0 en x et sin 0 = sin nπ = 0 en y, car m et n sont entiers. La solution est nulle sur les quatre bords.
À t = 0, la condition initiale est sin(mπx/Lₓ)·sin(nπy/L_y). C’est le cas V du projet avec Lₓ = L_y = 1 et m = n = 1.
## Task 2.1 - Taylor expansions

(a) f(x+h) = f + h f' + (h^2/2) f'' + (h^3/6) f''' + (h^4/24) f'''' + O(h^5)

(b) f(x-h) = f - h f' + (h^2/2) f'' - (h^3/6) f''' + (h^4/24) f'''' + O(h^5)

Les termes d'ordre impair (h, h^3) changent de signe, car on remplace h par -h
et (-h)^k = -h^k si k est impair, +h^k si k est pair.

(c) f(x+2h) = f + 2h f' + 2h^2 f'' + (4/3) h^3 f''' + (2/3) h^4 f'''' + O(h^5)

(coefficients : (2h)^k / k! = 2, 8/6 = 4/3, 16/24 = 2/3)

## Task 5.1 - Warm-up in one dimension

(a) Au noeud interieur x_i, avec la formule de la Task 3.2 (erreur O(dx^2)) :

    du_i/dt = alpha (u_{i+1} - 2 u_i + u_{i-1}) / dx^2

(b) Euler explicite : (u_i^{n+1} - u_i^n) / dt = alpha (u_{i+1}^n - 2 u_i^n + u_{i-1}^n) / dx^2, donc

    u_i^{n+1} = u_i^n + r (u_{i+1}^n - 2 u_i^n + u_{i-1}^n),   r = alpha dt / dx^2
              = (1 - 2r) u_i^n + r (u_{i+1}^n + u_{i-1}^n)

(c) La seule combinaison sans dimension est r = alpha dt / dx^2. C'est la
combinaison alpha t / L^2 de la Task 1.2, avec t remplace par dt et L par dx.
Elle compare le temps d'un pas au temps de diffusion sur une maille.

## Task 5.4 - Insulated edges in 2D

(a) Bord gauche x = 0 isole, les trois autres Dirichlet.
- Inconnues : les noeuds du bord (0, j) pour j = 1..N-1 (ils ne sont plus fixes).
  Ils sont mis a jour comme des noeuds interieurs.
- Noeuds fantomes : une colonne i = -1, aux lignes j = 1..N-1.
- Remplissage : u_x(0, y, t) = 0 donne (u_{1,j} - u_{-1,j}) / (2 dx) = 0,
  donc u_{-1,j} = u_{1,j} (miroir). Cette valeur est rafraichie a chaque pas,
  avant la mise a jour. Le noeud (0, j) devient :
  u_{0,j}^{n+1} = u_{0,j} + r (2 u_{1,j} + u_{0,j+1} + u_{0,j-1} - 4 u_{0,j}).

(b) Coin entre un bord isole et un bord Dirichlet : le coin appartient a la
condition Dirichlet, sa valeur est imposee, donc ce n'est pas une inconnue et
aucun fantome n'est necessaire (le choix de la valeur du coin doit etre fixe
dans la specification).
Coin entre deux bords isoles (par exemple x = 0 et y = 0) : le coin (0,0) est
une inconnue avec deux fantomes, u_{-1,0} = u_{1,0} et u_{0,-1} = u_{0,1} :
  u_{0,0}^{n+1} = u_{0,0} + r (2 u_{1,0} + 2 u_{0,1} - 4 u_{0,0}).
(Le stencil a 5 points n'utilise pas de noeud diagonal.)

(c) Flux non nul u_x(0, t) = h(t) :
  (u_1 - u_{-1}) / (2 dx) = h(t)  =>  u_{-1} = u_1 - 2 dx h(t).
Le noeud x_0 est alors mis a jour (Euler explicite, h evalue en t_n) :
  u_0^{n+1} = u_0 + r (u_1 - 2 u_0 + u_{-1})
            = u_0 + r (2 u_1 - 2 u_0 - 2 dx h(t_n)).
Si h = 0, on retrouve le miroir de l'exemple 5.2.

## Task 6.7 - Breaking it on purpose

Experiment: N = 64, case 1, r = 0.2625 (5% above the limit 1/4), 200 steps; control run with r = 0.24.
Code: tools/unstable_demo.cpp. Figure: results/unstable/instability.png.

Observation. With r = 0.24 the interior maximum stays at about 1. With r = 0.2625 it stays below 1 for about 110 steps, then grows: 1.5 at step 130, 11 at step 160, about 300 at step 200.

Where it appears. Near the four corners (around i, j = 10 and 55), where the hot edges meet the cold interior: the initial data is discontinuous there, so it contains the most high-frequency content.

Pattern. A checkerboard: the sign alternates from node to node (zoom and 1D cut in the figure).

Link with Task 6.3. The worst Fourier mode has theta = phi = pi, i.e. (-1)^(i+j), with amplification factor G = 1 - 8r. For r = 0.2625, G = -1.1, so |G| > 1: every step multiplies this mode by about 1.1 and flips its sign, which is the growing checkerboard. For r <= 1/4, |G| <= 1 and the mode no longer grows. The measured growth (about 1.09 per step) is slightly below 1.1 because on a 64-point grid the fastest mode is not exactly theta = pi.
