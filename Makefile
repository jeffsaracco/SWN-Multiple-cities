CC ?= cc
CFLAGS ?= -O2 -Wall -Wno-unused-variable -Wno-unused-but-set-variable -Wno-misleading-indentation
LDLIBS = -lm -lpthread

all: swn swn_cities client

%: %.c inet.h
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

# Serves the repo so viz/ can load run/epidemic.jsonl (add ?data=../run-docker/epidemic.jsonl for Docker runs)
viz:
	@echo "Open http://localhost:8000/viz/"
	python3 -m http.server 8000

clean:
	rm -rf swn swn_cities client run

.PHONY: all clean viz
