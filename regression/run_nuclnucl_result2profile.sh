#!/bin/sh -e
# Regression test: nucleotide search → result2profile WITHOUT backtrace (-a).
#
# This exercises the fix for result2profile where:
#   1. NucleotideMatrix must be used for alignment/backtrace on NT databases
#   2. Backtraces must be recomputed from alignment records that lack them
#      (e.g. 14-column records from offsetalignment in the blastn workflow)
#
# Without the fix, result2profile segfaults when trying to recompute
# nucleotide backtraces with amino-acid parameters.
#
# NOTE: PSSMCalculator is currently hardcoded to PROFILE_AA_SIZE=20 (amino
# acid alphabet), so the generated profile will be in amino acid space even
# for nucleotide input. This test validates that the pipeline completes
# without crashing; full nucleotide profile support requires refactoring
# PSSMCalculator.

# -- clean previous run --
rm -rf "${RESULTS:?}"/*

# -- create inline nucleotide FASTA files --
cat > "${RESULTS}/target.fasta" <<'EOF'
>seq1 identical_to_query
CTGCAGCTTGCCCTCAGAGACCGATCTCTCAGAGAGGTACATGGAATCGTGTTCCATCCCTGGATAACGGAACTCTCAGTCCTGCAG
>seq2 three_snps
CTGCAGCTTGCCCTCATAGACCGATCTCTCAGAGAGGTACATCGAATCGTGTTCCATCCATGGATAACGGAACTCTCAGTCCTGCAG
>seq3 short_deletion
CTGCAGCTTGCCCTCAGAGACCGATCTCTCAGAGAGCGTGTTCCATCCCTGGATAACGGAACTCTCAGTCCTGCAG
EOF

cat > "${RESULTS}/query.fasta" <<'EOF'
>query test_sequence
CTGCAGCTTGCCCTCAGAGACCGATCTCTCAGAGAGGTACATGGAATCGTGTTCCATCCCTGGATAACGGAACTCTCAGTCCTGCAG
EOF

# -- build databases --
"${MMSEQS}" createdb "${RESULTS}/target.fasta" "${RESULTS}/targetdb" --dbtype 2
"${MMSEQS}" createdb "${RESULTS}/query.fasta"  "${RESULTS}/querydb"  --dbtype 2

# -- search WITHOUT -a (no backtrace from aligner) --
"${MMSEQS}" search "${RESULTS}/querydb" "${RESULTS}/targetdb" \
    "${RESULTS}/result" "${RESULTS}/tmp" \
    --search-type 3

# -- result2profile must recompute backtraces internally --
# Without the fix this segfaults (set -e will abort the script).
"${MMSEQS}" result2profile "${RESULTS}/querydb" "${RESULTS}/targetdb" \
    "${RESULTS}/result" "${RESULTS}/profile"

# -- validate --
ERR=0

# profile database must exist and be non-empty
if [ ! -s "${RESULTS}/profile" ]; then
    echo "FAIL: profile database is empty or missing"
    ERR=$((ERR + 1))
fi

if [ ! -f "${RESULTS}/profile.dbtype" ]; then
    echo "FAIL: profile.dbtype is missing"
    ERR=$((ERR + 1))
fi

if [ ! -s "${RESULTS}/profile.index" ]; then
    echo "FAIL: profile.index is empty or missing"
    ERR=$((ERR + 1))
fi

# profile should have at least 1 entry (one per query)
ENTRY_COUNT=$(wc -l < "${RESULTS}/profile.index" 2>/dev/null || echo 0)
if [ "$ENTRY_COUNT" -lt 1 ]; then
    echo "FAIL: expected >=1 profile entries, got ${ENTRY_COUNT}"
    ERR=$((ERR + 1))
fi

awk -v actual="$ERR" -v target="0" \
    'BEGIN { print (actual == target) ? "GOOD" : "BAD"; print "Expected: ", target, "errors"; print "Actual: ", actual, "errors"; }' \
    > "${RESULTS}.report"
