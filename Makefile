CC ?= cc
CFLAGS ?= -O2 -Wall -Wno-unused-variable -Wno-unused-but-set-variable -Wno-misleading-indentation
LDLIBS = -lm -lpthread

all: swn swn_cities client

%: %.c inet.h
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

clean:
	rm -rf swn swn_cities client run

.PHONY: all clean
