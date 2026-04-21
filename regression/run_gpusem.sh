#!/bin/sh -e
# GPU server semaphore tests: idle CPU, signal cleanup, client timeouts

CLK_TCK=$(getconf CLK_TCK 2>/dev/null || echo 100)
CPU_IDLE_THRESHOLD=${CPU_IDLE_THRESHOLD:-5}
CPU_SAMPLE_SECONDS=${CPU_SAMPLE_SECONDS:-5}
THREADS=${THREADS:-64}
FAILURES=0
GPU_SERVER_PID=

fail() { printf 'FAIL: %s\n' "$1"; FAILURES=$((FAILURES + 1)); }
pass() { printf 'PASS: %s\n' "$1"; }
info() { printf 'INFO: %s\n' "$1"; }

cpu_ticks() {
    awk '{print $14 + $15}' "/proc/$1/stat" 2>/dev/null || echo 0
}

measure_cpu_pct() {
    local pid=$1 secs=${2:-$CPU_SAMPLE_SECONDS}
    local t1 t2 delta max_ticks
    t1=$(cpu_ticks "$pid")
    sleep "$secs"
    t2=$(cpu_ticks "$pid")
    delta=$((t2 - t1))
    max_ticks=$((secs * CLK_TCK))
    echo $((delta * 100 / max_ticks))
}

stop_gpu_server() {
    local sig=${1:-TERM}
    if [ -n "$GPU_SERVER_PID" ]; then
        kill -"$sig" "$GPU_SERVER_PID" 2>/dev/null || true
        wait "$GPU_SERVER_PID" 2>/dev/null || true
        GPU_SERVER_PID=
    fi
}

cleanup() { stop_gpu_server KILL; }
trap cleanup EXIT

# ── Setup ─────────────────────────────────────────────────────
head -n 1000 "${DATADIR}/query.fasta" > "${RESULTS}/query_500.fasta"
QUERY="${RESULTS}/query_500.fasta"
QUERYDB="${RESULTS}/query"
"${MMSEQS}" createdb "${QUERY}" "${QUERYDB}"

TARGET_FASTA="${DATADIR}/targetannotation.fasta"
TARGETDB="${RESULTS}/targetannotation"
"${MMSEQS}" createdb "${TARGET_FASTA}" "${TARGETDB}_db"
"${MMSEQS}" makepaddedseqdb "${TARGETDB}_db" "${TARGETDB}"

SHM_HASH=""

# Helper: start gpuserver and capture shm hash
start_gpu_server() {
    "${MMSEQS}" gpuserver "$TARGETDB" --max-seqs 1000 >"${RESULTS}/gpuserver.log" 2>&1 &
    GPU_SERVER_PID=$!
    SHM_HASH=""
    for i in $(seq 1 120); do
        sleep 0.5
        SHM_HASH=$(grep -oP '\b\d{10,}\b' "${RESULTS}/gpuserver.log" 2>/dev/null | head -1 || true)
        [ -n "$SHM_HASH" ] && break
    done
    if [ -z "$SHM_HASH" ]; then
        echo "ERROR: gpuserver did not print shm hash"
        cat "${RESULTS}/gpuserver.log"
        return 1
    fi
    # Wait for shm to be ready
    for i in $(seq 1 120); do
        sz=$(stat -c%s "/dev/shm/${SHM_HASH}" 2>/dev/null || echo 0)
        [ "$sz" -gt 0 ] && break
        sleep 0.5
    done
    info "gpuserver started  PID=$GPU_SERVER_PID  SHM=$SHM_HASH"
}

sem_path() { echo "/dev/shm/sem.${SHM_HASH}_sem"; }
shm_path() { echo "/dev/shm/${SHM_HASH}"; }

EXPECTED_ROC5="0.454819"

run_search() {
    rm -rf "${RESULTS}/result_db"* "${RESULTS}/search_tmp"
    "${MMSEQS}" search "$QUERYDB" "$TARGETDB" "${RESULTS}/result_db" "${RESULTS}/search_tmp" \
        -e 10000 --max-seqs 1000 --prefilter-mode 1 --db-load-mode 2 \
        --gpu 1 --gpu-server 1 --split 1 --threads "$THREADS" \
        "$@"
}

# Validate search results via ROC5 AUC.  Returns 0 on match.
check_search_results() {
    "${MMSEQS}" convertalis "$QUERYDB" "$TARGETDB" "${RESULTS}/result_db" "${RESULTS}/results_aln.m8"
    "${EVALUATE}" "$QUERY" "$TARGET_FASTA" "${RESULTS}/results_aln.m8" "${RESULTS}/evaluation_roc5.dat" 1000 1 \
        | tee "${RESULTS}/evaluation.log"
    ACTUAL=$(grep "^ROC5 AUC:" "${RESULTS}/evaluation.log" | cut -d" " -f3)
    if [ "$ACTUAL" = "$EXPECTED_ROC5" ]; then
        return 0
    else
        return 1
    fi
}

# ── Test 1: idle CPU usage ────────────────────────────────────
info "Test: idle CPU usage"
start_gpu_server
sleep 2
cpu_pct=$(measure_cpu_pct "$GPU_SERVER_PID")
info "  measured idle CPU: ${cpu_pct}%"
if [ "$cpu_pct" -lt "$CPU_IDLE_THRESHOLD" ]; then
    pass "idle CPU ${cpu_pct}% < ${CPU_IDLE_THRESHOLD}%"
else
    fail "idle CPU ${cpu_pct}% >= ${CPU_IDLE_THRESHOLD}% (busy-waiting?)"
fi

# ── Test 2: semaphore exists ──────────────────────────────────
info "Test: semaphore exists"
if [ -e "$(sem_path)" ]; then
    pass "semaphore present at $(sem_path)"
else
    fail "semaphore NOT found at $(sem_path)"
fi

# ── Test 3: functional search (ROC5 AUC validation) ──────────
info "Test: functional search"
run_search
if check_search_results; then
    pass "search ROC5 AUC = $EXPECTED_ROC5"
else
    fail "search ROC5 AUC = $ACTUAL (expected $EXPECTED_ROC5)"
fi

# ── Test 4: CPU returns to idle after search ──────────────────
info "Test: CPU after search"
sleep 2
cpu_pct=$(measure_cpu_pct "$GPU_SERVER_PID")
info "  measured post-search idle CPU: ${cpu_pct}%"
if [ "$cpu_pct" -lt "$CPU_IDLE_THRESHOLD" ]; then
    pass "post-search idle CPU ${cpu_pct}% < ${CPU_IDLE_THRESHOLD}%"
else
    fail "post-search idle CPU ${cpu_pct}% >= ${CPU_IDLE_THRESHOLD}%"
fi

# ── Test 5: SIGTERM → clean shutdown ──────────────────────────
info "Test: SIGTERM clean shutdown"
stop_gpu_server TERM
sleep 1
if [ -e "$(sem_path)" ]; then
    fail "semaphore still exists after SIGTERM"
else
    pass "semaphore cleaned up after SIGTERM"
fi
if [ -e "$(shm_path)" ]; then
    fail "shm still exists after SIGTERM"
else
    pass "shm cleaned up after SIGTERM"
fi

# ── Test 6: SIGINT → clean shutdown ──────────────────────────
info "Test: SIGINT clean shutdown"
start_gpu_server
sleep 2
stop_gpu_server INT
sleep 1
if [ -e "$(sem_path)" ]; then
    fail "semaphore still exists after SIGINT"
else
    pass "semaphore cleaned up after SIGINT"
fi
if [ -e "$(shm_path)" ]; then
    fail "shm still exists after SIGINT"
else
    pass "shm cleaned up after SIGINT"
fi

# ── Test 7: SIGKILL recovery ─────────────────────────────────
info "Test: SIGKILL recovery"
start_gpu_server
sleep 2
kill -9 "$GPU_SERVER_PID" 2>/dev/null || true
wait "$GPU_SERVER_PID" 2>/dev/null || true
GPU_SERVER_PID=
rm -f "$(shm_path)"
start_gpu_server
sleep 2
if kill -0 "$GPU_SERVER_PID" 2>/dev/null; then
    pass "server restarted after SIGKILL"
else
    fail "server failed to restart after SIGKILL"
fi
run_search
if check_search_results; then
    pass "search after SIGKILL ROC5 AUC = $EXPECTED_ROC5"
else
    fail "search after SIGKILL ROC5 AUC = $ACTUAL (expected $EXPECTED_ROC5)"
fi
stop_gpu_server TERM
sleep 1

# ── Test 8: rapid start/stop ─────────────────────────────────
info "Test: rapid start/stop cycles"
leaked=0
for i in 1 2 3; do
    start_gpu_server
    sleep 2
    stop_gpu_server TERM
    sleep 1
    if [ -e "$(sem_path)" ]; then
        fail "semaphore leaked on cycle $i"
        leaked=1
        break
    fi
done
if [ "$leaked" -eq 0 ]; then
    pass "no semaphore leaks after 3 start/stop cycles"
fi

# ── Test 9: client no-wait timeout ───────────────────────────
info "Test: client timeout=0 (no server)"
stop_gpu_server KILL
rm -rf "${RESULTS}/result_db"* "${RESULTS}/search_tmp"
if "${MMSEQS}" search "$QUERYDB" "$TARGETDB" "${RESULTS}/result_db" "${RESULTS}/search_tmp" \
    -e 10000 --max-seqs 1000 --prefilter-mode 1 --db-load-mode 2 \
    --gpu 1 --gpu-server 1 --split 1 --gpu-server-wait-timeout 0 \
    --threads "$THREADS" >"${RESULTS}/client_nowait.log" 2>&1; then
    fail "client should have failed with no server and timeout=0"
else
    pass "client exited with error when no server (timeout=0)"
fi

# ── Test 10: client short timeout ────────────────────────────
info "Test: client timeout=3 (no server)"
rm -rf "${RESULTS}/result_db"* "${RESULTS}/search_tmp"
t_start=$(date +%s)
if "${MMSEQS}" search "$QUERYDB" "$TARGETDB" "${RESULTS}/result_db" "${RESULTS}/search_tmp" \
    -e 10000 --max-seqs 1000 --prefilter-mode 1 --db-load-mode 2 \
    --gpu 1 --gpu-server 1 --split 1 --gpu-server-wait-timeout 3 \
    --threads "$THREADS" >"${RESULTS}/client_timeout.log" 2>&1; then
    fail "client should have timed out with no server"
else
    t_end=$(date +%s)
    elapsed=$((t_end - t_start))
    if [ "$elapsed" -ge 2 ] && [ "$elapsed" -le 10 ]; then
        pass "client timed out after ${elapsed}s (expected ~3s)"
    else
        fail "client timeout took ${elapsed}s (expected ~3s)"
    fi
fi

# ── Test 11: client waits, server starts later ────────────────
info "Test: client-first, server-later"
rm -rf "${RESULTS}/result_db"* "${RESULTS}/search_tmp"
"${MMSEQS}" search "$QUERYDB" "$TARGETDB" "${RESULTS}/result_db" "${RESULTS}/search_tmp" \
    -e 10000 --max-seqs 1000 --prefilter-mode 1 --db-load-mode 2 \
    --gpu 1 --gpu-server 1 --split 1 --gpu-server-wait-timeout 60 \
    --threads "$THREADS" >"${RESULTS}/client_wait.log" 2>&1 &
CLIENT_PID=$!
sleep 3
start_gpu_server
if wait "$CLIENT_PID"; then
    if check_search_results; then
        pass "client-first search ROC5 AUC = $EXPECTED_ROC5"
    else
        fail "client-first search ROC5 AUC = $ACTUAL (expected $EXPECTED_ROC5)"
    fi
else
    fail "client-first search failed"
fi
stop_gpu_server TERM
sleep 1

# ── Test 12: no leftover gpuserver processes ──────────────────
info "Test: no leftover gpuserver processes"
leftover=$(pgrep -f "${MMSEQS}.*gpuserver" 2>/dev/null || true)
if [ -z "$leftover" ]; then
    pass "no leftover gpuserver processes"
else
    fail "leftover gpuserver processes: $leftover"
fi

# ── Report ────────────────────────────────────────────────────
if [ "$FAILURES" -eq 0 ]; then
    echo "GOOD" > "${RESULTS}.report"
    echo "All gpusem tests passed" >> "${RESULTS}.report"
else
    echo "BAD" > "${RESULTS}.report"
    echo "$FAILURES gpusem test(s) failed" >> "${RESULTS}.report"
fi
