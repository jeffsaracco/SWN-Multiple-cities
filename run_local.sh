#!/bin/sh
# Runs the multi-city simulation as N local processes, one per city.
#   usage: ./run_local.sh [num_cities] [seed]
# Output (epidemic.jsonl, network.dat, sum_file*.csv, one log per city) goes to run/.
# Watch it on the globe with `make viz`.
set -e

CITIES=${1:-4}
SEED=${2:-1}
BASE_PORT=5530

cd "$(dirname "$0")"
make -s swn_cities
rm -rf run && mkdir run
cp fluvars.in run/
cd run

# lat lon name, for the globe; the first $CITIES are used
PLACES="40.71 -74.01 New York
51.51 -0.13 London
35.68 139.69 Tokyo
-23.55 -46.63 Sao Paulo
6.52 3.38 Lagos
19.08 72.88 Mumbai
-33.87 151.21 Sydney
30.04 31.24 Cairo
55.76 37.62 Moscow
19.43 -99.13 Mexico City
39.90 116.40 Beijing
-34.60 -58.38 Buenos Aires
-26.20 28.05 Johannesburg
34.05 -118.24 Los Angeles
1.35 103.82 Singapore
64.15 -21.94 Reykjavik"
if [ "$CITIES" -gt 16 ]; then
  echo "run_local.sh knows 16 places; add more to PLACES for more cities" >&2
  exit 1
fi

echo "$PLACES" | head -n "$CITIES" | awk -v base="$BASE_PORT" \
  '{ print "127.0.0.1", base + NR, $0 }' > cities.txt

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
echo "Watch it on the globe: make viz"
exit $status
