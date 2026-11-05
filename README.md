# Small World Network: Multiple Cities

A simulation of influenza spreading through a population, where each person is
a node in a [small-world network](https://en.wikipedia.org/wiki/Watts%E2%80%93Strogatz_model).
This version splits the population into several "cities", each running as its
own process (originally one per machine in a cluster), connected by a handful
of edges between people in different cities. Those cross-city edges are how
flu can travel from one city to another.

Written by Jeff Saracco in March 2006, during my senior year of college, as
part of a research project. I didn't finish it that year. My lab partner wrote
the flu-spreading-between-cities piece, which isn't included here. In 2026 I
got it compiling and running again on a modern Mac, with each city as a local
process or a Docker container.

## What's here

| File | What it is |
| --- | --- |
| `swn_cities.c` | **The multi-city version.** Builds a small-world network in each city, wires cities together with cross-city edges, rewires those edges, and writes the combined network to `network.dat`. There's no flu model in this file yet. |
| `swn.c` | **The flu model**, version 3.0 of the "Network Influenza Model". It covers a single city's network, flu strains (H/N subtypes plus a strain bitstring that mutates), latent and infectious periods, and five vaccination strategies. It also has an early two-city socket prototype: city 1 accepts "infect host N with this strain" messages from another city. |
| `client.c` | A test harness for `swn.c`. It connects to city 1 and infects whatever host numbers you type. |
| `inet.h` | Shared socket includes, plus the address and port used by `swn.c` and `client.c`. |
| `fluvars.in` | Simulation parameters, read from the working directory at startup. |
| `run_local.sh` | Runs N cities as local processes. |
| `Dockerfile`, `docker-compose.yml`, `docker/cities.txt` | Run each city in its own container. |

## Building

Requires a C compiler and `make`. Tested with Apple clang on macOS. The Docker
image builds with gcc on Debian.

```sh
make            # builds swn, swn_cities, and client
```

## Running multiple cities locally

```sh
./run_local.sh            # 4 cities, seed 1
./run_local.sh 8 42       # 8 cities, seed 42
```

The script builds `swn_cities`, then creates `run/` and writes a `cities.txt`
there with one `127.0.0.1 <port>` line per city. It starts one process per
city, waits for them all to finish, and prints each city's log. Outputs land
in `run/`:

- `network.dat`: the full network across all cities (format below)
- `sum_file<seed> city <N>.csv`: each city's parameter summary
- `city<N>.log`: each city's console output

To run cities by hand, for example on separate machines, write a `cities.txt`
where line N is `host port` for city N, then start every city:

```sh
./swn_cities <seed> <city number> [cities file]   # cities file defaults to cities.txt
```

Cities can start in any order. Each city keeps retrying its connections to
lower-numbered cities for up to a minute.

## Running with Docker

```sh
docker compose up --build    # or `docker-compose up --build` (Homebrew installs it standalone)
docker compose down          # remove the stopped containers
```

This starts `city1` through `city4`, one container each. They find each other
by service name using `docker/cities.txt`. Output goes to `run-docker/`. To add
a city, add a line to `docker/cities.txt` and a matching service to
`docker-compose.yml`.
Tested with Colima on Apple Silicon.

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

## How the multi-city setup works

1. **Connect.** Every city listens on its port, dials each lower-numbered city,
   then accepts connections from the higher-numbered ones. The dialing city
   sends its number first, so every pair of cities ends up with exactly one TCP
   connection. Each city then runs one thread per peer to handle incoming
   commands.
2. **Local network.** Each city builds a ring lattice where every person is
   connected to `original_number_neighbors` people on each side. Each edge is
   then rewired to a random person with probability `swnP`.
3. **Global init.** Odd-numbered cities add `n_cross_city_edges / 2` random
   edges to each of their two ring neighbors (city N-1 and city N+1).
4. **Global rewire.** Cities take turns, passing a token 1 → 2 → … → N. Each
   one rewires its cross-city edges with probability `swnP_city` to a random
   person in a random other city.
5. **Output.** Cities take turns appending their network to `network.dat`.

A barrier separates the phases. Each city sends a barrier message to every
other city and waits until it has one from all of them. Commands on a
connection arrive in order, so passing the barrier means every edge change
sent before it has already been applied.

The commands cities send each other are plain native-endian `int`s:

| Code | Meaning | Followed by |
| --- | --- | --- |
| `-9` | add a cross-city neighbor | sender's host, receiver's host, sender's city |
| `-6` | delete a cross-city neighbor | sender's host, receiver's host, sender's city |
| `-10` | barrier | nothing |
| `-12` | your turn to rewire | nothing |

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
a count. Both programs use only the first value of each list:
`pop_sizes[0]` for the population, `swnP_values[0]` for the rewiring
probability, and so on. `swn_cities.c` reads the whole file but only uses the
population, `swnP`, and `original_number_neighbors`. The rest are flu
parameters for `swn.c`.

The cross-city settings are constants at the top of `swn_cities.c`:
`n_cross_city_edges = 10` and `swnP_city = 0.5`.

## Known limitations

- `swn_cities.c` builds the multi-city network but doesn't simulate flu on it.
- `swn.c`'s remote infections all share a single `struct flu`, so hosts
  infected that way share one flu record.
- Each person has room for `6 × original_number_neighbors` neighbors. Nothing
  checks that limit, so a person who collects too many rewired edges would
  write past the end of the array.
- Messages are raw `int`s in native byte order, so every city has to run on
  the same kind of machine. That's always true with local processes or Docker.
