import math

x = 1.0
exact = -math.sin(x)

print("Task 3.2(c) - Second derivative")
print("h, error")

for h in [0.1, 0.05, 0.025, 0.0125]:
    approx = (
        math.sin(x + h)
        - 2 * math.sin(x)
        + math.sin(x - h)
    ) / (h * h)

    error = abs(approx - exact)
    print(f"{h:.4f}, {error:.10f}")
