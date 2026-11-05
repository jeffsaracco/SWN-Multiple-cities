#!/bin/sh
# Runs the multi-city simulation as N local processes, one per city.
#   usage: ./run_local.sh [num_cities] [seed]
# Output (network.dat, sum_file*.csv, one log per city) goes to run/.
set -e

CITIES=${1:-4}
SEED=${2:-1}
BASE_PORT=5530

cd "$(dirname "$0")"
make -s swn_cities
rm -rf run && mkdir run
cp fluvars.in run/
cd run

for i in $(seq 1 "$CITIES"); do
  echo "127.0.0.1 $((BASE_PORT + i))"
done > cities.txt

pids=""
for i in $(seq 1 "$CITIES"); do
  ../swn_cities "$SEED" "$i" cities.txt > "city$i.log" 2>&1 &
  pids="$pids $!"
done
trap 'kill $pids 2>/dev/null' INT TERM

status=0
for pid in $pids; do
  wait "$pid" || status=1
done

for i in $(seq 1 "$CITIES"); do
  sed "s/^/[city $i] /" "city$i.log"
done
echo
echo "Output in $(pwd)"
exit $status
