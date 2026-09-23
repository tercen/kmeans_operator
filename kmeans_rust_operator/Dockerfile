# syntax=docker/dockerfile:1.7
# kmeans_rust_operator — static tier (create-rust-operator §5): one musl
# binary on scratch. The build context is this directory.

# ---- builder ----
FROM rust:1.94-bookworm AS builder
RUN apt-get update && apt-get install -y --no-install-recommends \
        musl-tools protobuf-compiler pkg-config git ca-certificates \
 && rm -rf /var/lib/apt/lists/* && rustup target add x86_64-unknown-linux-musl
WORKDIR /build
RUN cargo install cargo-chef --locked

# Dependency layer. cargo-chef's recipe carries the package version, and a release commit always
# bumps it, so the recipe would change on every release and the ~10 minute dependency build would
# never be reused. Normalise the version for the recipe only, then restore the real manifest.
COPY Cargo.toml Cargo.lock ./
RUN cp Cargo.toml Cargo.toml.keep \
 && sed -i 's/^version = ".*"/version = "0.0.0"/' Cargo.toml \
 && cargo chef prepare --recipe-path recipe.json \
 && cargo chef cook --release --target x86_64-unknown-linux-musl --recipe-path recipe.json \
 && mv Cargo.toml.keep Cargo.toml

COPY src ./src
RUN cargo build --release --target x86_64-unknown-linux-musl --bin kmeans_rust_operator \
 && ls -l target/x86_64-unknown-linux-musl/release/kmeans_rust_operator
# scratch has no directories at all: build the writable temp dir here and copy it in.
RUN mkdir -p /tmp-op && chown 1000:1000 /tmp-op

# ---- runtime ----
FROM scratch
COPY --from=builder /etc/ssl/certs/ca-certificates.crt /etc/ssl/certs/
COPY --from=builder /build/target/x86_64-unknown-linux-musl/release/kmeans_rust_operator /usr/local/bin/kmeans_rust_operator
COPY --from=builder --chown=1000:1000 /tmp-op /tmp
USER 1000:1000
WORKDIR /operator
ENV RUST_BACKTRACE=1 RUST_LOG=info TMPDIR=/tmp
ENTRYPOINT ["/usr/local/bin/kmeans_rust_operator"]
