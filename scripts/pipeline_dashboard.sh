#!/usr/bin/env bash
# =============================================================================
# pipeline_dashboard.sh -- live, at-a-glance status of the running pipeline.
#
# A local, read-only TUI that refreshes every few seconds (default 5s). It is the
# fast complement to pipeline_monitor.sh (which snapshots to GitHub hourly): use
# this when you want to watch what the run is doing RIGHT NOW.
#
# Usage:
#   bash scripts/pipeline_dashboard.sh              # live view, refresh every 5s
#   bash scripts/pipeline_dashboard.sh -i 2         # refresh every 2s
#   bash scripts/pipeline_dashboard.sh --once       # render once and exit
#
# By default it FOLLOWS THE MOST RECENTLY ACTIVE run (newest of logs/run_fixed.log
# vs logs/run_calculated.log) and reads the live PID from the zymo-pipeline unit.
# Env overrides: DASH_INTERVAL (s), RUN_LOG, PID_FILE.
# Nothing here writes to the repo, touches git, or affects the run.
# =============================================================================
cd "$(dirname "$0")/.." || exit 1                     # -> repo root

INTERVAL="${DASH_INTERVAL:-5}"
PROG_LOG="logs/pipeline_progress.log"
SS="data/sample_sheet.tsv"
ONCE=0

# Auto-select the CURRENT run's log + pid = the newest of the per-profile logs, so
# the dashboard follows whichever run (fixed or calculated) is active now instead
# of a hardcoded one. Still overridable via the RUN_LOG / PID_FILE env vars.
_auto_log="$(ls -1t logs/run_fixed.log logs/run_calculated.log 2>/dev/null | head -1)"
[ -z "$_auto_log" ] && _auto_log="logs/run_fixed.log"
RUN_LOG="${RUN_LOG:-$_auto_log}"
PID_FILE="${PID_FILE:-${RUN_LOG%.log}.pid}"

while [ $# -gt 0 ]; do
  case "$1" in
    --once|-1)      ONCE=1 ;;
    -i|--interval)  INTERVAL="$2"; shift ;;
    -h|--help)      sed -n '2,20p' "$0"; exit 0 ;;
    *) echo "unknown arg: $1" >&2; exit 2 ;;
  esac
  shift
done

N_LIBS=$(( $(wc -l < "$SS" 2>/dev/null || echo 1) - 1 ))
[ "$N_LIBS" -lt 1 ] && N_LIBS=48

# ---- colours (only when writing to a terminal) ------------------------------
if [ -t 1 ]; then
  B=$'\033[1m'; D=$'\033[2m'; R=$'\033[0m'
  GRN=$'\033[32m'; YEL=$'\033[33m'; RED=$'\033[31m'; CYN=$'\033[36m'; MAG=$'\033[35m'
else
  B=""; D=""; R=""; GRN=""; YEL=""; RED=""; CYN=""; MAG=""
fi

hr() { printf '%s\n' "────────────────────────────────────────────────────────────────────────"; }

svc_state() { systemctl --user is-active "$1" 2>/dev/null || echo unknown; }
colour_state() {
  case "$1" in
    active)             printf '%s%s%s' "$GRN" "$1" "$R" ;;
    failed|inactive)    printf '%s%s%s' "$RED" "$1" "$R" ;;
    *)                  printf '%s%s%s' "$YEL" "$1" "$R" ;;
  esac
}

render() {
  local now pid alive elapsed cpu mem_p rss nlwp stage_line done_n done_list cur res_dir
  now="$(date '+%Y-%m-%d %H:%M:%S %z')"
  # prefer the live systemd MainPID (authoritative) over a possibly-stale pid file
  pid="$(systemctl --user show -p MainPID --value zymo-pipeline 2>/dev/null)"
  if [ -z "$pid" ] || [ "$pid" = 0 ]; then pid="$(cat "$PID_FILE" 2>/dev/null)"; fi
  if [ -n "$pid" ] && kill -0 "$pid" 2>/dev/null; then alive=1; else alive=0; fi

  # ---- header ---------------------------------------------------------------
  printf '%s' "${B}${CYN}"
  printf '  ZYMO PIPELINE DASHBOARD'
  printf '%s\n' "$R"
  printf '  %s%s%s   refresh %ss\n' "$D" "$now" "$R" "$INTERVAL"
  hr

  # ---- services -------------------------------------------------------------
  local sp sm
  sp="$(svc_state zymo-pipeline)"; sm="$(svc_state zymo-monitor)"
  printf '  %-9s pipeline: %s' "SERVICES" "$(colour_state "$sp")"
  printf '    monitor(hourly push): %s\n' "$(colour_state "$sm")"
  if [ "$alive" = 1 ]; then
    elapsed="$(ps -o etime= -p "$pid" 2>/dev/null | tr -d ' ')"
    read -r cpu mem_p rss nlwp <<<"$(ps -o %cpu=,%mem=,rss=,nlwp= -p "$pid" 2>/dev/null)"
    printf '  %-9s %sRUNNING%s pid %s  up %s  cpu %s%%  mem %s%%  rss %sG  threads %s\n' \
      "PROC" "$GRN" "$R" "$pid" "${elapsed:-?}" "${cpu:-?}" "${mem_p:-?}" \
      "$(awk "BEGIN{printf \"%.0f\", ${rss:-0}/1048576}")" "${nlwp:-?}"
  else
    printf '  %-9s %sNOT RUNNING%s (pid %s) -- run has stopped or finished\n' \
      "PROC" "$RED" "$R" "${pid:-?}"
  fi
  hr

  # ---- stage / progress -----------------------------------------------------
  stage_line="$(grep -aE 'RUNNING  stage|DONE     stage' "$RUN_LOG" 2>/dev/null | tail -1 | sed 's/^[[:space:]]*//')"
  printf '  %sSTAGE%s   %s\n' "$B" "$R" "${stage_line:-<no stage line yet>}"

  # stage-01 external-tools completion (kraken2.out written last per library)
  done_n=0; done_list=""
  for d in work/*/; do
    [ -s "${d}blastn.tsv" ] && [ -s "${d}kraken2.out" ] && { done_n=$((done_n + 1)); done_list+=" ${d#work/}"; }
  done
  printf '  %sstage01%s external tools: %s%d/%d%s libraries complete\n' \
    "$D" "$R" "$GRN" "$done_n" "$N_LIBS" "$R"

  # current library still in BLAST (has reads.fasta but no kraken2.out yet)
  cur=""
  for d in $(ls -dt work/*/ 2>/dev/null); do
    [ -s "${d}reads.fasta" ] || continue
    [ -s "${d}kraken2.out" ] || { cur="$d"; break; }
  done
  if [ -n "$cur" ] && [ "$done_n" -lt "$N_LIBS" ]; then
    local reads_n blast_done blast_sz
    reads_n="$(grep -c '^>' "${cur}reads.fasta" 2>/dev/null || echo '?')"
    blast_done="$(cut -f1 "${cur}blastn.tsv" 2>/dev/null | uniq | wc -l 2>/dev/null || echo 0)"
    blast_sz="$(du -h "${cur}blastn.tsv" 2>/dev/null | cut -f1)"
    printf '  %scurrent%s %s -- BLAST ~%s/%s reads (blastn.tsv %s)\n' \
      "$D" "$R" "${cur#work/}" "$blast_done" "$reads_n" "${blast_sz:-0}"
  fi
  hr

  # ---- what it is doing right now (live log tail) ---------------------------
  printf '  %sACTIVITY%s  (last log lines)\n' "$B" "$R"
  # progress bars emit \r frames incl. empty "clearing" frames (a bar of spaces);
  # split on \r, then keep only informative lines: any "NN%" frame OR plain text
  # (no bar pipe). This drops the noisy space-only redraw frames.
  tail -n 60 "$RUN_LOG" 2>/dev/null | tr '\r' '\n' \
    | grep -aE '[0-9]+%|^[^|]*$' | grep -avE '^[[:space:]]*$' \
    | tail -n 8 | sed 's/^/    /'
  printf '\n'
  hr

  # ---- newest outputs -------------------------------------------------------
  res_dir="$(ls -1dt results/gt_*/ 2>/dev/null | head -1)"
  printf '  %sOUTPUTS%s   %s\n' "$B" "$R" "${res_dir:-<no results dir yet>}"
  if [ -n "$res_dir" ]; then
    find "$res_dir" -type f -printf '%TY-%Tm-%Td %TH:%TM  %10s  %p\n' 2>/dev/null \
      | sort | tail -n 5 | sed 's/^/    /'
  fi
  hr

  # ---- resources ------------------------------------------------------------
  local mem_line swap_line load disk
  mem_line="$(free -g 2>/dev/null | awk '/^Mem:/{printf "%d/%dG used", $3, $2}')"
  swap_line="$(free -g 2>/dev/null | awk '/^Swap:/{printf "%d/%dG", $3, $2}')"
  load="$(cut -d' ' -f1-3 /proc/loadavg 2>/dev/null)"
  disk="$(df -h . 2>/dev/null | awk 'NR==2{print $4" free / "$2}')"
  printf '  %sHOST%s      mem %s   swap %s   load %s   disk %s\n' \
    "$B" "$R" "${mem_line:-?}" "${swap_line:-?}" "${load:-?}" "${disk:-?}"
  # swap warning (the old oomd-crash signal)
  local sw_used; sw_used="$(free -m 2>/dev/null | awk '/^Swap:/{print $3}')"
  if [ -n "$sw_used" ] && [ "$sw_used" -gt 2048 ]; then
    printf '  %s          WARNING: swap %sMB in use -- memory pressure (oomd is masked, kernel OOM still possible)%s\n' \
      "$YEL" "$sw_used" "$R"
  fi
  hr

  # ---- hourly GitHub snapshot status ---------------------------------------
  local last_snap last_push
  last_snap="$(grep -aE '^snapshot' "$PROG_LOG" 2>/dev/null | tail -1 | cut -d: -f2- | sed 's/^ *//')"
  last_push="$(grep -aE '\[pushed|\[push FAILED' "$PROG_LOG" 2>/dev/null | tail -1 | grep -oE '\[(pushed|push FAILED)[^]]*\]')"
  printf '  %sGITHUB%s    last hourly snapshot: %s\n' "$D" "$R" "${last_snap:-none yet}"
  printf '            %s\n' "${last_push:-<no push recorded>}"

  [ "$ONCE" = 0 ] && printf '\n  %sPress Ctrl-C to exit.%s\n' "$D" "$R"
}

if [ "$ONCE" = 1 ]; then
  render
  exit 0
fi

trap 'printf "\n"; exit 0' INT
while true; do
  printf '\033[H\033[J'      # home + clear
  render
  sleep "$INTERVAL"
done
