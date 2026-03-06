#!/bin/sh -e
# Regression test: nucleotide search → result2msa WITHOUT backtrace (-a).
#
# This exercises the fix for result2msa/result2profile where:
#   1. NucleotideMatrix must be used (not SubstitutionMatrix) for NT databases
#   2. Backtraces must be recomputed from alignment records that lack them
#      (e.g. 14-column records from offsetalignment in the blastn workflow)
#
# Without the fix, this produces "DUMMY" all-gap MSA rows or segfaults.

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

# -- result2msa must recompute backtraces internally --
"${MMSEQS}" result2msa "${RESULTS}/querydb" "${RESULTS}/targetdb" \
    "${RESULTS}/result" "${RESULTS}/result.a3m" \
    --msa-format-mode 6

# -- validate --
ERR=0

# result.a3m must exist and be non-empty
if [ ! -s "${RESULTS}/result.a3m" ]; then
    echo "result.a3m is empty or missing"
    ERR=$((ERR + 1))
fi

# must not contain DUMMY placeholder sequences
if grep -q "DUMMY" "${RESULTS}/result.a3m"; then
    echo "result.a3m contains DUMMY placeholder"
    ERR=$((ERR + 1))
fi

# must have at least 2 sequences (query + one hit)
SEQ_COUNT=$(grep -c "^>" "${RESULTS}/result.a3m" || true)
if [ "$SEQ_COUNT" -lt 2 ]; then
    echo "expected >=2 sequences, got ${SEQ_COUNT}"
    ERR=$((ERR + 1))
fi

# no all-gap alignment rows
if grep -qE "^-+$" "${RESULTS}/result.a3m"; then
    echo "result.a3m contains all-gap rows"
    ERR=$((ERR + 1))
fi

awk -v actual="$ERR" -v target="0" \
    'BEGIN { print (actual == target) ? "GOOD" : "BAD"; print "Expected: ", target, "errors"; print "Actual: ", actual, "errors"; }' \
    > "${RESULTS}.report"
