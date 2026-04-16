# Vendored MPSolve patches

PSMM depends on [MPSolve](https://github.com/robol/MPSolve) with a small
Advanpix extension for Mahler-measure-aware early exit. MPSolve uses autotools
and has no CMake, so rather than maintaining a separate fork we:

1. Pin an upstream source snapshot by commit hash
   (see [`upstream-version.cmake`](upstream-version.cmake)).
2. Ship our patches in [`patches/`](patches/) as unified-diff files,
   applied with `patch -p1` in sorted filename order.
3. Let CMake fetch, patch, `autogen.sh` / `configure` / `make install`
   MPSolve into a build-local prefix under `${CMAKE_BINARY_DIR}/mpsolve-install`.
4. Expose the result as the `mpsolve::mps` imported target.

This is a temporary bridge until PSMM (and optionally the Advanpix MPSolve
fork) are published; when that happens we switch to a private git repo for
MPSolve and drop the `patches/` directory.

## Build-host prerequisites

- autoconf, automake, libtool
- bison, flex
- pkg-config
- GMP headers (`libgmp-dev` on Debian/Ubuntu)

On Ubuntu/Debian:

```sh
sudo apt-get install -y autoconf automake libtool bison flex pkg-config libgmp-dev
```

## Maintaining the patches

Work against upstream MPSolve in a separate local clone:

```sh
git clone https://github.com/robol/MPSolve /path/to/mpsolve-dev
cd /path/to/mpsolve-dev
git checkout <pinned-commit>
git checkout -b advanpix/3.2.2
# ... hack, commit, test via PSMM build ...
git format-patch <pinned-commit>..advanpix/3.2.2 -o /path/to/PSMM/mpsolve/patches/
```

Commit the regenerated `.patch` files. On next upstream bump:

```sh
git rebase <new-upstream-commit>
# regenerate patches, update upstream-version.cmake URL/URL_HASH
```

## Files

- `upstream-version.cmake` — pinned upstream version, URL, and SHA256.
- `patches/*.patch` — patches applied in sorted order by filename. Prefix with
  `NNNN-` for ordering.
- `CMakeLists.txt` — glue that drives ExternalProject_Add.
