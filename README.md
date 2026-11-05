# Small World Network: Multiple Cities

A simulation of influenza spreading through a population, where each person is
a node in a [small-world network](https://en.wikipedia.org/wiki/Watts%E2%80%93Strogatz_model).
The population is split into several "cities", each running as its own process
(originally one per machine in a cluster), connected by a handful of edges
between people in different cities. Those cross-city edges are how the flu
travels from one city to another. A 3D globe replays the run: each city is a
hospital that grows, beats faster, and gets covered in red veins as more of its
people get sick, and every cross-city connection is drawn green, yellow, or red
by the health of the two people it joins.

Written by Jeff Saracco in March 2006, during my senior year of college, as
part of a research project. I didn't finish it that year; the flu spreading
between cities was finished by the project's next owners after I graduated.
My lab partner built a globe visualization of it, which isn't included here.
In 2026 I got the code compiling and running again on a modern Mac, with each
city as a local process or a Docker container. I also moved the flu model into
the multi-city program and rebuilt the globe from my memory of it, in `viz/`.

## What's here

| File | What it is |
| --- | --- |
| `swn_cities.c` | **The multi-city simulation.** Builds a small-world network in each city, wires cities together with cross-city edges, then runs the flu model day by day, sending infections across cities. Writes `epidemic.jsonl` (per-day counts, for the globe) and `network.dat`. |
| `viz/index.html` | **The globe.** Replays an `epidemic.jsonl` in the browser with three.js. |
| `swn.c` | **The original single-city flu model**, version 3.0 of the "Network Influenza Model". It covers flu strains (H/N subtypes plus a strain bitstring that mutates), latent and infectious periods, and five vaccination strategies. It also has an early two-city socket prototype. `swn_cities.c` carries a copy of its flu code. |
| `client.c` | A test harness for `swn.c`. It connects to city 1 and infects whatever host numbers you type. |
| `inet.h` | Shared socket includes, plus the address and port used by `swn.c` and `client.c`. |
| `fluvars.in` | Simulation parameters, read from the working directory at startup. |
| `run_local.sh` | Runs N cities as local processes. |
| `Dockerfile`, `docker-compose.yml`, `docker/cities.txt` | Run each city in its own container. |

## Quick start

```sh
make                # builds swn, swn_cities, and client
./run_local.sh      # runs 4 cities; output in run/
make viz            # then open http://localhost:8000/viz/
```

## Running multiple cities locally

```sh
./run_local.sh            # 4 cities, seed 1
./run_local.sh 8 42       # 8 cities, seed 42 (up to 16 with the built-in list of places)
```

The script builds `swn_cities` and creates `run/`. It writes a `cities.txt`
there, with one line per city giving its address and a real place for the
globe. Then it starts one process per city, waits for them all to finish, and
prints each city's log. Outputs land in `run/`:

- `epidemic.jsonl`: infected, recovered, and vaccinated counts for every city on every day, plus cross-city infections (format below)
- `network.dat`: the full network across all cities (format below)
- `sum_file<seed> city <N>.csv`: each city's parameter summary and epidemic stats
- `city<N>.log`: each city's console output

To run cities by hand, for example on separate machines, write a `cities.txt`
where line N describes city N, then start every city:

```
# host port [lat lon name]
10.0.0.5 5529 40.71 -74.01 New York
10.0.0.6 5529 51.51 -0.13 London
```

```sh
./swn_cities <seed> <city number> [cities file]   # cities file defaults to cities.txt
```

Cities can start in any order. Each city keeps retrying its connections to
lower-numbered cities for up to a minute. Without `lat lon name`, a city is
called "City N" and placed on the equator.

The same seed on the same platform gives the same run, byte for byte. macOS
and Linux implement `rand()` differently, so a local run and a Docker run with
the same seed differ.

## Running with Docker

```sh
docker compose up --build    # or `docker-compose up --build` (Homebrew installs it standalone)
docker compose down          # remove the stopped containers
```

This starts `city1` through `city4`, one container each. They find each other
by service name using `docker/cities.txt`. Output goes to `run-docker/`. To
watch it, run `make viz` and open
`http://localhost:8000/viz/?data=../run-docker/epidemic.jsonl`. To add a city,
add a line to `docker/cities.txt` and a matching service to
`docker-compose.yml`. Tested with Colima on Apple Silicon.

## The globe

`make viz` serves the repo at `http://localhost:8000`. `viz/` loads
`run/epidemic.jsonl` and plays it. You can also drop any `epidemic.jsonl` onto
the page.

- **Hospitals.** Each city is a hospital. As the share of its people who are sick goes up, the hospital grows, its heart beats faster (40 to 180 bpm), its windows burn red, and veins climb its walls, with each beat pulsing up them. As people recover, it all shrinks back. Everything scales with the square root of a city's sick count divided by the run's peak, so a small outbreak still shows.
- **Globe.** Red veins spread over the ground around each sick city, and the atmosphere turns from blue to red as the world gets sicker.
- **Connections.** Every cross-city edge is its own arc between two hospitals: green when both people are healthy, yellow when one is sick, and red when both are. Several edges between the same two cities fan out side by side. When an infection crosses, a red pulse runs along one of the connections that has a sick person on it.
- **Controls.** Space plays and pauses, the arrow keys step one day, and you can drag the timeline. Click a hospital or a city in the list to fly in close; Esc or "Whole globe" flies back out.
- **Links.** `?day=12` opens paused on day 12, `?city=London` opens close up on London, and `?data=path` loads another file.

The page loads three.js and the night-earth texture from the jsDelivr CDN, so
it needs a network connection. Without the texture, it falls back to a plain
gridded globe.

## Running the flu model (`swn.c`)

```sh
mkdir -p run-flu && cd run-flu && cp ../fluvars.in .
../swn <seed> 1          # city 1: runs the flu model and listens for remote infections
../client                # in another terminal: type host numbers to infect
```

City 1 listens on port 5529. It sleeps one second per simulated day for the
first 120 days, so there's time to send it infections. `../swn <seed> 2` runs
city 2, which connects to city 1. It only sends one test message, on day 10:
infect host 512, which city 1 ignores unless the population is bigger than 512.

## How the multi-city simulation works

1. **Connect.** Every city listens on its port, dials each lower-numbered city,
   then accepts connections from the higher-numbered ones. The dialing city
   sends its number first, so every pair of cities ends up with exactly one TCP
   connection. Each city then runs one thread per peer to handle incoming
   commands.
2. **Local network.** Each city builds a ring lattice where every person is
   connected to `original_number_neighbors` people on each side. Each edge is
   then rewired to a random person with probability `swnP`.
3. **Global init.** Taking turns, odd-numbered cities add
   `n_cross_city_edges / 2` random edges to each of their two ring neighbors
   (city N-1 and city N+1).
4. **Global rewire.** Taking turns, each city rewires its cross-city edges
   with probability `swnP_city` to a random person in a random other city.
5. **Initial infection.** The flu starts in `seed_city` (city 1 by default).
   Every city also sends city 1 its cross-city edges for the globe. Each edge is
   sent once, by the lower-numbered of its two cities.
6. **Each day,** every city:
   1. vaccinates, if it's a vaccination day
   2. recovers people whose infection has run its course
   3. spreads the flu from each infectious person to their neighbors. A neighbor
      in another city is sent to that city as an infection attempt.
   4. waits at a barrier until every city has sent today's infections
   5. tries the infections other cities sent, against each person's immunity,
      in city order
   6. reports its counts to city 1, along with which of its people who have
      cross-city neighbors are sick. City 1 then writes that day's line of
      `epidemic.jsonl`.
7. **Output.** Cities take turns appending their network to `network.dat`.

A barrier means: each city sends a barrier message to every other city, then
waits until it has received one from all of them. Commands on a connection
arrive in order, so once a city is past a barrier, every command sent before
it has arrived. Taking turns and applying infections in city order is what
makes runs repeatable.

A cross-city infection goes through the same check as a local one:
`try_infect()` applies the neighbor contact rate, then the H-protein match and
cross-immunity against the strains that person has had.

The commands cities send each other are plain native-endian `int`s:

| Code | Meaning | Followed by |
| --- | --- | --- |
| `-9` | add a cross-city neighbor | sender's host, receiver's host, sender's city |
| `-6` | delete a cross-city neighbor | sender's host, receiver's host, sender's city |
| `-10` | barrier | nothing |
| `-13` | infect | day, receiver's host, strain H, N, and strain bits |
| `-14` | day stats (to city 1) | day, infected, recovered, vaccinated, then infections imported from each city |
| `-15` | cross-city edge (to city 1, once) | sender's host, other city, other city's host |
| `-16` | sick border hosts (to city 1, daily) | a count, then the day and the sender's hosts with cross-city neighbors who are sick |

Commands with a variable length, like `-16`, send their count first.

### `epidemic.jsonl` format

One JSON object per line. The first line describes the run, and then there's
one line per day. Cities are numbered from 0, in `cities.txt` order.

```
{"type":"meta","version":1,"days":100,"population":10000,"seedCity":0,"cities":[{"name":"New York","lat":40.71,"lon":-74.01}, ...],"links":[[0,1],[0,3], ...]}
{"type":"day","day":0,"infected":[1,0,0,0],"recovered":[0,0,0,0],"vaccinated":[0,0,0,0],"imports":[],"linkStates":"00000000000000000000"}
{"type":"day","day":3,"infected":[72,1,0,0], ... ,"imports":[[0,1,1]],"linkStates":"00010000000000000000"}
```

- `links` lists the two cities of every cross-city edge.
- `linkStates` has one digit per link, in the same order: `0` means both people are healthy, `1` means one is sick, and `2` means both are.
- `imports` lists `[from city, to city, count]` for each pair of cities where an infection crossed and took that day. City 1 writes each line as soon as the
day finishes, so the file could be streamed to a live viewer.

### `network.dat` format

```
<number of cities>
<people in city 1>
<degree>  <neighbor id> <neighbor city (0-based)>  <neighbor id> <neighbor city>  ...   (one line per person)
<people in city 2>
...
```

## Parameters

`fluvars.in` holds one value per line, with a label after it. Lists start with
a count, and only the first value of each list is used: `pop_sizes[0]` for
the population (10,000 per city), `swnP_values[0]` for the rewiring
probability, and so on.

The cross-city settings are constants at the top of `swn_cities.c`:
`n_cross_city_edges = 10`, `swnP_city = 0.5`, and `seed_city = 1` (set it to
0 to start the flu in every city).

## Known limitations

- The flu code exists twice: in `swn.c` (the original) and in `swn_cities.c`
  (adapted for neighbors in other cities). A change to the model has to be
  made in both.
- Cross-cut and clustering-based vaccination only consider neighbors within a
  city, since a city can't see the neighbor lists of people in other cities.
- `swn.c`'s remote infections all share a single `struct flu`, so hosts
  infected that way share one flu record.
- Each person has room for `6 × original_number_neighbors` neighbors. Nothing
  checks that limit, so a person who collects too many rewired edges would
  write past the end of the array.
- Messages are raw `int`s in native byte order, so every city has to run on
  the same kind of machine. That's always true with local processes or Docker.
