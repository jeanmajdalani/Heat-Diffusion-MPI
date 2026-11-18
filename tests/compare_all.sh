#!/usr/bin/env bash
# Compare MPI outputs with the serial reference, bit for bit (cmp).
# Usage: tests/compare_all.sh [decomp]   (default: strips)
# Run from the repository root after building into ./build.
set -u
DECOMP="${1:-strips}"
HEAT=./build/heat
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
fail=0; total=0

for c in 1 2 V; do
  for N in 8 9 16 33 64; do
    for P in 1 2 3 4 5 7 8; do
      [ "$P" -le $((N - 1)) ] || continue
      $HEAT --case "$c" --N "$N" --r 0.2 --steps 40 --out "$TMP/serial.bin" > /dev/null
      mpirun --oversubscribe -np "$P" $HEAT --case "$c" --N "$N" --r 0.2 --steps 40 \
        --decomp "$DECOMP" --out "$TMP/mpi.bin" > /dev/null 2>&1
      total=$((total + 1))
      if cmp -s "$TMP/serial.bin" "$TMP/mpi.bin"; then
        echo "IDENTIQUE  case=$c N=$N P=$P"
      else
        echo "DIFFERENT  case=$c N=$N P=$P"; fail=$((fail + 1))
      fi
    done
  done
done
echo "----"
echo "$DECOMP: $((total - fail))/$total identiques"
[ "$fail" -eq 0 ]