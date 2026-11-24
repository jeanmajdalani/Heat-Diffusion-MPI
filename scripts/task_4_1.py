import math
import csv
from pathlib import Path

out = Path("results/p2")
out.mkdir(parents=True, exist_ok=True)

lam = 50.0
exact = math.exp(-lam)

def euler(n, backward=False):
    dt = 1.0 / n
    u = 1.0
    peak = abs(u)

    for _ in range(n):
        if backward:
            u = u / (1.0 + lam*dt)
        else:
            u = (1.0 - lam*dt)*u

        peak = max(peak, abs(u))

    return u, peak

print("STABILITY EXPERIMENT")

with open(out / "task_4_1_stability.csv", "w") as f:
    writer = csv.writer(f)
    writer.writerow([
        "method", "dt", "final_value",
        "max_abs_value", "bounded"
    ])

    for backward in [False, True]:
        method = "backward" if backward else "forward"

        for n in [100, 50, 26, 25, 24, 20, 10, 2, 1]:
            dt = 1.0 / n
            u, peak = euler(n, backward)
            bounded = peak <= 1.0 + 1e-12

            writer.writerow([
                method, dt, u, peak, bounded
            ])

            print(
                f"{method:8s} dt={dt:.8f} "
                f"max={peak:.6e} bounded={bounded}"
            )

print("\nFORWARD EULER CONVERGENCE")

previous = None

with open(out / "task_4_1_convergence.csv", "w") as f:
    writer = csv.writer(f)
    writer.writerow(["dt", "absolute_error", "order"])

    for n in [10000, 20000, 40000, 80000, 160000]:
        dt = 1.0 / n
        u, _ = euler(n)
        error = abs(u - exact)

        order = (
            math.log2(previous/error)
            if previous is not None else None
        )

        writer.writerow([dt, error, order])

        print(
            f"dt={dt:.8e} "
            f"error={error:.8e} order={order}"
        )

        previous = error

print("\nBoth CSV files saved successfully.")
