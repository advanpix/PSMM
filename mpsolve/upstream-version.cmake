# MPSolve upstream pin.
#
# The MPSolve project tags releases rarely. 3.2.2 exists as an in-tree
# version bump (commit de7ebfc7 on master, dated 2025-03-21, marked
# "Release 3.2.2.") but was never tagged. We pin to that commit hash to
# get the 3.2.2 source tree verbatim.
#
# Post-3.2.2 commits on master carry further memory-leak and compiler-
# compatibility fixes. Before bumping, check the patch series in
# mpsolve/patches/ still applies cleanly to the new snapshot.

set(PSMM_MPSOLVE_VERSION  "3.2.2")
set(PSMM_MPSOLVE_COMMIT   "de7ebfc7afc4834a0c9f92a04be7abdf5943d446")
set(PSMM_MPSOLVE_URL      "https://github.com/robol/MPSolve/archive/${PSMM_MPSOLVE_COMMIT}.tar.gz")
set(PSMM_MPSOLVE_URL_HASH "SHA256=42666c994fe533c276715fb086843461dea9684975fa822c9809ae8324369067")
