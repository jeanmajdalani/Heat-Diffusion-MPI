"""Standard-library checks for the Part 1 serial heat solver."""
import math
import os
import struct
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
HEAT = ROOT / "build" / "heat"


def run(args, expect_success=True):
    result = subprocess.run([str(HEAT), *args], capture_output=True, text=True)
    if (result.returncode == 0) != expect_success:
        raise AssertionError(f"unexpected result: {result.stdout}\n{result.stderr}")
    return result.stdout.strip() if expect_success else result.stderr.strip()


def load(path, N):
    data = Path(path).read_bytes()
    assert len(data) == (N + 1) ** 2 * 8, (len(data), N)
    return struct.unpack("<" + "d" * (N + 1) ** 2, data)


def fields():
    with tempfile.TemporaryDirectory() as tmp:
        for case in ["1", "2", "V"]:
            path = os.path.join(tmp, f"case{case}.bin")
            run(["--case", case, "--N", "8", "--r", "0.2", "--steps", "2", "--out", path])
            values = load(path, 8)
            at = lambda i, j: values[j * 9 + i]
            for j in range(9):
                for i in range(9):
                    v = at(i, j)
                    assert math.isfinite(v), (case, i, j, v)
                    assert -1e-14 <= v <= 1.0 + 1e-14, (case, i, j, v)
                    if case == "1" and (i in (0, 8) or j in (0, 8)):
                        assert v == 1.0, (case, i, j, v)
                    if case == "2":
                        if j in (0, 8):
                            assert v == 1.0, (case, i, j, v)
                        elif i in (0, 8):
                            assert v == 0.0, (case, i, j, v)
                    if case == "V" and (i in (0, 8) or j in (0, 8)):
                        assert v == 0.0, (case, i, j, v)
    print("PASS: Cases 1, 2 and V; boundaries; little-endian output")


def stability():
    message = run(["--case", "V", "--N", "8", "--r", "0.26", "--steps", "2"], False)
    assert "unstable" in message
    print("PASS: unstable r rejected")


def convergence():
    errors = []
    for N in [16, 32, 64, 128, 256]:
        output = run(["--case", "V", "--N", str(N), "--r", "0.2", "--T", "0.05"])
        data = dict(token.split("=", 1) for token in output.split())
        errors.append(float(data["err_l2_h"]))
        assert float(data["T"]) == 0.05
        assert float(data["r"]) <= 0.2 + 1e-14
        print(f"N={N:3d} steps={data['steps']:>6} err_l2_h={errors[-1]:.8g}")
    orders = [math.log(errors[i-1]/errors[i], 2) for i in range(1, len(errors))]
    print("Observed orders:", ", ".join(f"{order:.4f}" for order in orders))
    assert all(errors[i] < errors[i - 1] for i in range(1, len(errors)))
    assert abs(orders[-1] - 2.0) < 0.05, orders
    print("PASS: convergence close to order 2")


if __name__ == "__main__":
    assert HEAT.exists(), f"Build first: {HEAT}"
    fields()
    stability()
    convergence()
