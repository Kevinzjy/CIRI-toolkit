FROM rust:1.85-bookworm

ENV CARGO_TERM_COLOR=never \
    CARGO_NET_RETRY=10

# musl-tools supplies musl-gcc for mimalloc's C sources; binutils/file/tar
# support stripping, static-link auditing, and the Makefile tar-stream export.
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        musl-tools \
        binutils \
        file \
        tar \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/* \
    && rustup target add x86_64-unknown-linux-musl

# Pin the generic x86-64 baseline so artifacts do not accidentally require a
# newer CPU, and use musl's static CRT for CentOS 7 compatibility.
ENV TARGET=x86_64-unknown-linux-musl \
    RUSTFLAGS="-C target-feature=+crt-static -C target-cpu=x86-64" \
    CC_x86_64_unknown_linux_musl=musl-gcc \
    AR_x86_64_unknown_linux_musl=ar

WORKDIR /work

# Fetch dependencies separately so Docker layer caching avoids re-downloading
# the registry when only source files change.
COPY Cargo.toml Cargo.lock README.md ./
RUN cargo fetch --locked --target "$TARGET"

COPY src ./src

RUN cargo build --release --locked --target "$TARGET" \
    && mkdir -p /export \
    && for bin in ciri ciri-simulator ciri-merge ciri-assemble; do \
        cp "target/${TARGET}/release/${bin}" "/export/${bin}" \
        && strip "/export/${bin}" \
        && chmod 0755 "/export/${bin}"; \
    done \
    && rm -rf /work/target
