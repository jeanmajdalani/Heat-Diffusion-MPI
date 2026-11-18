#!/usr/bin/env bash
# Benchmark Strips (or Blocks) : plusieurs N, P, mappings, repetitions.
# Sortie : un CSV propre + le log brut complet (a garder pour le rapport).
#
# Exemples :
#   scripts/benchmark.sh -N "256 512" -P "1 2 4" -R 3 -O          # test local
#   scripts/benchmark.sh -H cluster/hosts -N "1024 2048" -P "1 2 4 8 16" -M "node slot"
set -u

BIN=build/heat          # executable
HOSTFILE=""             # vide = machine locale
NS="512"                # tailles de grille N
PS="1 2 4"              # nombres de processus
REPS=3                  # repetitions par configuration
STEPS=500               # --steps (decision commune)
R=0.24                  # --r (decision commune)
DECOMP=strips
MAPS="slot"             # "node", "slot" ou "node slot"
CASE=1
OVERSUB=0               # -O : --oversubscribe (tests locaux)
IFACE=""                # ex. tailscale0 pour forcer l'interface reseau
OUT=""                  # prefixe de sortie (defaut : results/bench_<date>)

usage() { sed -n '2,9p' "$0"; echo "Options : -b bin -H hostfile -N 'tailles' -P 'procs' -R reps -S steps -r r -D decomp -M 'maps' -c case -I iface -o prefixe -O"; exit 1; }

while getopts "b:H:N:P:R:S:r:D:M:c:I:o:Oh" opt; do
  case $opt in
    b) BIN=$OPTARG;; H) HOSTFILE=$OPTARG;; N) NS=$OPTARG;; P) PS=$OPTARG;;
    R) REPS=$OPTARG;; S) STEPS=$OPTARG;; r) R=$OPTARG;; D) DECOMP=$OPTARG;;
    M) MAPS=$OPTARG;; c) CASE=$OPTARG;; I) IFACE=$OPTARG;; o) OUT=$OPTARG;;
    O) OVERSUB=1;; *) usage;;
  esac
done

[ -x "$BIN" ] || { echo "Executable introuvable : $BIN (compilez d'abord)"; exit 1; }
[ -z "$HOSTFILE" ] || [ -f "$HOSTFILE" ] || { echo "Hostfile introuvable : $HOSTFILE"; exit 1; }
command -v mpirun >/dev/null || { echo "mpirun introuvable"; exit 1; }

mkdir -p results
[ -n "$OUT" ] || OUT="results/bench_$(date +%Y%m%d_%H%M%S)"
CSV="$OUT.csv"; RAW="$OUT.raw.log"

{
  echo "# date: $(date -Is)"
  echo "# host: $(hostname)   bin: $BIN   decomp: $DECOMP   case: $CASE"
  echo "# steps: $STEPS  r: $R  reps: $REPS  hostfile: ${HOSTFILE:-local}"
  echo "# $(mpirun --version | head -1)"
  [ -z "$HOSTFILE" ] || { echo "# hostfile:"; sed 's/^/#   /' "$HOSTFILE"; }
} > "$RAW"
echo "decomp,map_by,N,P,rep,steps,r,t_total,t_comp,t_comm" > "$CSV"

total=0
for N in $NS; do for P in $PS; do for M in $MAPS; do total=$((total+REPS)); done; done; done
done_count=0

for N in $NS; do
 for P in $PS; do
  if [ "$DECOMP" = strips ] && [ "$P" -gt $((N-1)) ]; then
    echo "saute P=$P pour N=$N (P > N-1)"; continue; fi
  for M in $MAPS; do
   args=(-np "$P" --map-by "$M")
   [ -z "$HOSTFILE" ] || args+=(--hostfile "$HOSTFILE")
   [ "$OVERSUB" -eq 0 ] || args+=(--oversubscribe)
   [ -z "$IFACE" ] || args+=(--mca btl_tcp_if_include "$IFACE" --mca oob_tcp_if_include "$IFACE")
   for rep in $(seq 1 "$REPS"); do
     cmd=(mpirun "${args[@]}" "$BIN" --case "$CASE" --N "$N" --steps "$STEPS" --r "$R" --decomp "$DECOMP")
     echo "### ${cmd[*]}  (rep $rep)" >> "$RAW"
     line=$("${cmd[@]}" 2>>"$RAW" | tee -a "$RAW" | grep "^decomp=" | tail -1)
     if [ -z "$line" ]; then
       echo "ECHEC : N=$N P=$P map=$M rep=$rep (voir $RAW)"; continue; fi
     get() { echo "$line" | tr ' ' '\n' | sed -n "s/^$1=//p"; }
     echo "$DECOMP,$M,$N,$P,$rep,$(get steps),$(get r),$(get t_total),$(get t_comp),$(get t_comm)" >> "$CSV"
     done_count=$((done_count+1))
     printf "\r[%d/%d] N=%s P=%s map=%s rep=%s   " "$done_count" "$total" "$N" "$P" "$M" "$rep"
   done
  done
 done
done
echo
echo "CSV : $CSV"; echo "Log brut : $RAW"