#!/usr/bin/env bash
# =============================================================================
# pipeline_monitor.sh -- hourly progress logger for the stage-01 pipeline run.
#
# Every $MONITOR_INTERVAL seconds it appends a progress snapshot to the TRACKED
# file logs/pipeline_progress.log and commits + pushes ONLY that file to origin.
# It exits after a final snapshot once the monitored pipeline PID is gone.
#
# Design notes:
#   * Commits ALWAYS succeed (local); the push is BEST-EFFORT -- a transient
#     auth/network failure is logged and retried next cycle, never stopping the
#     logging. Local commits accumulate and go out on the next successful push.
#   * Only logs/pipeline_progress.log is ever staged (explicit path); it never
#     touches work/ or results/ artifacts.
#
# Env overrides: MONITOR_INTERVAL (s, default 3600), RUN_LOG, PID_FILE.
# =============================================================================
cd "$(dirname "$0")/.." || exit 1          # -> repo root

INTERVAL="${MONITOR_INTERVAL:-3600}"
RUN_LOG="${RUN_LOG:-logs/run_fixed.log}"
PID_FILE="${PID_FILE:-logs/run_fixed.pid}"
PROG_LOG="logs/pipeline_progress.log"
SS="data/sample_sheet.tsv"
N_LIBS=$(( $(wc -l < "$SS" 2>/dev/null || echo 1) - 1 ))
BRANCH="$(git rev-parse --abbrev-ref HEAD 2>/dev/null || echo main)"

pid_alive() { local p; p="$(cat "$PID_FILE" 2>/dev/null)"; [ -n "$p" ] && kill -0 "$p" 2>/dev/null; }

snapshot() {
  local now pid cur d reads_n blast_done blast_sz done_list
  now="$(date '+%Y-%m-%d %H:%M:%S %z')"
  pid="$(cat "$PID_FILE" 2>/dev/null)"
  DONE_N=0; done_list=""
  for d in work/*/; do
    [ -s "${d}blastn.tsv" ] && [ -s "${d}kraken2.out" ] && { DONE_N=$((DONE_N + 1)); done_list+=" ${d#work/}"; }
  done
  cur=""
  for d in $(ls -dt work/*/ 2>/dev/null); do
    [ -s "${d}kraken2.out" ] || { cur="$d"; break; }
  done
  {
    echo "===================================================================="
    echo "snapshot   : $now"
    if pid_alive; then echo "pipeline   : RUNNING (pid $pid)"; else echo "pipeline   : NOT RUNNING (pid ${pid:-?})"; fi
    echo "stage-01   : $DONE_N / $N_LIBS libraries complete"
    [ -n "$done_list" ] && echo "  complete :$done_list"
    if [ -n "$cur" ]; then
      reads_n="$(grep -c '^>' "${cur}reads.fasta" 2>/dev/null || echo '?')"
      blast_done="$(cut -f1 "${cur}blastn.tsv" 2>/dev/null | uniq | wc -l 2>/dev/null || echo 0)"
      blast_sz="$(du -h "${cur}blastn.tsv" 2>/dev/null | cut -f1)"
      echo "  current  : ${cur#work/}  BLAST ~${blast_done}/${reads_n} reads (blastn.tsv ${blast_sz:-0})"
    fi
    echo "  disk     : $(df -h . 2>/dev/null | awk 'NR==2{print $4" free / "$2}')"
    echo "  log tail :"
    tail -n 5 "$RUN_LOG" 2>/dev/null | sed 's/^/    | /'
  } >> "$PROG_LOG"
}

commit_push() {
  local msg="$1"
  git add "$PROG_LOG" 2>/dev/null
  git diff --cached --quiet 2>/dev/null && return 0     # nothing new staged
  git commit -q -m "$msg" >/dev/null 2>&1 || return 0
  if git push -q origin "$BRANCH" >/dev/null 2>&1; then
    echo "  [pushed $(date '+%F %T')]" >> "$PROG_LOG"          # rolls into next commit
  else
    echo "  [push FAILED $(date '+%F %T') -- committed locally, will retry next cycle]" >> "$PROG_LOG"
  fi
}

echo "#### monitor started $(date '+%F %T %z')  interval=${INTERVAL}s  pid_file=${PID_FILE}" >> "$PROG_LOG"
while true; do
  snapshot
  commit_push "pipeline progress $(date '+%F %T'): ${DONE_N}/${N_LIBS} libs done"
  if ! pid_alive; then
    echo "#### pipeline no longer running -- monitor exiting $(date '+%F %T')" >> "$PROG_LOG"
    commit_push "pipeline progress: run finished $(date '+%F %T')"
    break
  fi
  sleep "$INTERVAL"
done
