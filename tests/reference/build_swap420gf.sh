#!/bin/bash
# Reproducible gfortran build of pristine SWAP 4.2.0 (originally Intel ifort).
# Goal: a Linux binary built with the SAME flags as the modern build, so 4.2.0
# physics can be compared against a modernized build under one compiler.
#
# NO source edits were required: all 153 TTUTIL objects + 56 SWAP files
# compiled and statically linked unchanged.
set -e

TTUTIL=/tmp/swap-legacy/source_ttutil_4.27
SWP=/tmp/swap-legacy/source_swp_4.2.0
BUILD=/tmp/swap-legacy/build-gf

# Flags mirror the modern build (apples-to-apples physics comparison).
# -w / -fallow-argument-mismatch / -fallow-invalid-boz are gfortran-version
# compatibility shims only; they do NOT change numerics.
FLAGS="-O2 -ffree-line-length-none -std=legacy -finit-local-zero -w -fallow-argument-mismatch -fallow-invalid-boz"
INC="-I${TTUTIL} -I${SWP}"

# ---------------------------------------------------------------------------
# 1. TTUTIL 4.27 -> libttutil427.a
#    Module-defining .f90 first (dependency order), then remaining .f90,
#    then all fixed-form .for subroutines (gfortran auto-detects fixed form).
# ---------------------------------------------------------------------------
cd "${TTUTIL}"
rm -f *.o *.mod

# base modules (no deps)
gfortran ${FLAGS} ${INC} -c ttutilprefs.f90
gfortran ${FLAGS} ${INC} -c ttutil.f90
# dependent modules
gfortran ${FLAGS} ${INC} -c rdmodulettutil.f90   # USE ttutilPrefs
gfortran ${FLAGS} ${INC} -c outdat.f90           # USE Module_OUTDAT (self)
# remaining free-form
for f in dtnow fatalerr fopengstandard ifindc messini messinq \
         newlinestandard openlogf outar2 outsel recread recreadi recreadt upperc; do
  gfortran ${FLAGS} ${INC} -c ${f}.f90
done
# all fixed-form subroutines
for f in *.for; do
  gfortran ${FLAGS} ${INC} -c "$f"
done

ar rcs libttutil427.a *.o
cp libttutil427.a *.mod "${BUILD}/"

# ---------------------------------------------------------------------------
# 2. SWAP 4.2.0 .f90 in the canonical Intel order (encodes module deps).
#    Built in ${BUILD} so the TTUTIL .mod files are on the include path.
# ---------------------------------------------------------------------------
cd "${BUILD}"
rm -f *.o; cp "${TTUTIL}"/*.mod . 2>/dev/null || true

for f in $(cat "${BUILD}/swap_order.txt"); do
  gfortran ${FLAGS} ${INC} -c "${SWP}/${f}"
done

# ---------------------------------------------------------------------------
# 3. Link (static; self-contained binary).
# ---------------------------------------------------------------------------
gfortran ${FLAGS} -static -o "${BUILD}/swap420gf" *.o libttutil427.a

echo "Built: ${BUILD}/swap420gf"
