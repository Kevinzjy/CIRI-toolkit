SHELL := /bin/sh

IMAGE ?= ciri-toolkit-builder:musl
DIST ?= dist
CONTEXT ?= ..
BINARIES ?= ciri ciri-simulator ciri-merge ciri-assemble
VERSION := $(shell awk -F '"' '/^version[[:space:]]*=/ {print $$2; exit}' ../Cargo.toml)
PACKAGE ?= ciri-toolkit-$(VERSION)-linux-x86_64-musl.tar.gz

.PHONY: all image export check package clean

all: export

image:
	docker build -f Dockerfile -t $(IMAGE) $(CONTEXT)

export: image
	mkdir -p $(DIST)
	set -eu; \
	tarfile=$(DIST)/.bins.tar; \
	docker run --rm --entrypoint /bin/sh $(IMAGE) -c \
		'cd /export && tar -cf - $(BINARIES)' > "$$tarfile"; \
	tar -xf "$$tarfile" -C $(DIST); \
	rm -f "$$tarfile"; \
	for bin in $(BINARIES); do chmod 0755 "$(DIST)/$$bin"; done

check: export
	docker run --rm --entrypoint /bin/sh $(IMAGE) -c '\
		set -eu; \
		for bin in $(BINARIES); do \
			echo "==> $${bin}"; \
			test -x "/export/$${bin}"; \
			file "/export/$${bin}"; \
			if readelf -l "/export/$${bin}" | grep -Eq "(^|[[:space:]])INTERP([[:space:]]|$$)"; then \
				echo "ERROR: $${bin} has a dynamic interpreter" >&2; \
				exit 1; \
			fi; \
		done'
	./$(DIST)/ciri --help > /dev/null

package: export
	tar -czf $(DIST)/$(PACKAGE) -C $(DIST) $(BINARIES)

clean:
	rm -rf $(DIST)
