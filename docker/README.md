# Docker static build

This directory builds static Linux x86-64 release binaries inside Docker. The
result is intended for release distribution, not as a runtime Docker image. It
uses the `x86_64-unknown-linux-musl` target and the generic x86-64 CPU baseline
to avoid depending on the host glibc or newer CPU instructions.

## Build and export

From this `docker/` directory:

```bash
make
```

This builds the builder image and writes these files to `docker/dist/`:

```text
ciri
ciri-simulator
ciri-merge
ciri-assemble
```

The binaries are stripped and owned by the invoking user. The build context is
the repository root, and the root `.dockerignore` limits it to Cargo metadata
and `src/`, so local test datasets and build caches are not sent to Docker.

## Verify

```bash
make check
```

This checks that each exported binary has no dynamic interpreter and that the
main CLI starts. Manual checks are:

```bash
file dist/ciri
ldd dist/ciri || true
readelf -l dist/ciri | grep -i interp || true
```

Expected results are `statically linked`, `not a dynamic executable`, and no
`INTERP` segment.

## Package

To create a release tarball:

```bash
make package
```

The archive name uses the current Cargo package version, for example:

```text
dist/ciri-toolkit-0.2.4-linux-x86_64-musl.tar.gz
```

## Release process

`.github/workflows/release-linux-x86_64.yml` builds and attaches the static
binary tarball automatically whenever a GitHub release is published. It checks
out the release tag, so the tag must exist before the release is created:

1. Tag the release commit on `main`, then push `main` and the tag to GitHub
   (`main` is the only branch synced to GitHub):

   ```bash
   git tag vX.Y.Z
   git push <github-remote> main vX.Y.Z
   ```

2. Create the release from the existing tag:

   ```bash
   gh release create vX.Y.Z --verify-tag --generate-notes
   ```

   `--verify-tag` fails if the tag does not exist yet, which enforces the
   tag-first ordering. Creating the release first would make GitHub create the
   tag at the default branch head instead.

3. The workflow then builds the musl binaries and attaches
   `ciri-toolkit-<version>-linux-x86_64-musl.tar.gz` plus a `.sha256` checksum
   to the release. The archive name uses the version from `Cargo.toml`, which
   should match the tag.

## Runtime notes

Static linking removes the toolkit's direct glibc and dynamic-loader
dependency, but external programs invoked by the toolkit still must exist:

- `ciri` requires `samtools`, either in `PATH` or configured through
  `SAMTOOLS=/path/to/samtools`.
- `ciri-simulator` uses `pigz` when available and falls back to `gzip`.

On a CentOS 7 x86-64 host, no Docker installation is needed to run the exported
binaries. For example:

```bash
./ciri --help
install -m 0755 ciri /usr/local/bin/ciri
```
