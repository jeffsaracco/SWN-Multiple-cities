FROM debian:bookworm-slim AS build
RUN apt-get update \
 && apt-get install -y --no-install-recommends gcc libc6-dev make \
 && rm -rf /var/lib/apt/lists/*
WORKDIR /src
COPY Makefile inet.h swn_cities.c ./
RUN make swn_cities

FROM debian:bookworm-slim
COPY --from=build /src/swn_cities /usr/local/bin/swn_cities
# Output files (network.dat, sum_file*.csv) are written to the working dir,
# which docker-compose.yml mounts so every city shares it.
WORKDIR /out
ENTRYPOINT ["swn_cities"]
