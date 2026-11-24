import math
import csv
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

out = Path("results/p2")
out.mkdir(parents=True, exist_ok=True)

x = 1.0
exact = -math.sin(x)

hs = []
errors = []

for k in range(1, 11):
    h = 10.0 ** (-k)

    approx = (
        math.sin(x+h)
        - 2*math.sin(x)
        + math.sin(x-h)
    ) / h**2

    error = abs(approx - exact)

    hs.append(h)
    errors.append(error)

    print(f"h={h:.0e}  error={error:.12e}")

with open(out / "task_3_5.csv", "w") as f:
    writer = csv.writer(f)
    writer.writerow(["h", "absolute_error"])
    writer.writerows(zip(hs, errors))

plt.figure(figsize=(8, 5))
plt.loglog(hs, errors, "o-", label="Absolute error")
plt.xlabel("Step size h")
plt.ylabel("Absolute error")
plt.title("Task 3.5 - Truncation and round-off errors")
plt.gca().invert_xaxis()
plt.grid(True, which="both", alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(out / "task_3_5.png", dpi=200)

best = errors.index(min(errors))

print(f"\nMinimum error at h={hs[best]:.0e}")
print(f"Minimum error: {errors[best]:.12e}")
print("CSV and graph saved successfully.")
