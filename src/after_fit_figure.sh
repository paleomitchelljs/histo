#!/bin/bash
# Wait for a running fit to exit, then draw its summary figure with a status line taken from the fit log.
# Runs on its own (detach with nohup), so the figure does not depend on an interactive session.
# Usage: src/after_fit_figure.sh <pid of the R process> <fit log> <fit tag> <data prefix> <model 32|35|36> [label]
# Writes outputs/fig_age_element_<tag>.png and logs to outputs/after_fit_<tag>.log.
set -u
PID=$1; LOG=$2; TAG=$3; PRE=$4; MODEL=$5; LABEL=${6:-"Model $5, data $4"}
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT" || exit 1
OUT="outputs/after_fit_${TAG}.log"
echo "$(date) waiting for PID $PID ($LOG)" > "$OUT"
while kill -0 "$PID" 2>/dev/null; do sleep 60; done
echo "$(date) fit process ended" >> "$OUT"
if ! grep -q "DONE model" "$LOG"; then
  echo "$(date) fit log has no DONE line: the fit failed or was stopped; no figure drawn" >> "$OUT"; exit 1
fi
DIV=$(grep -Eo "\(([0-9.]+)%\) transitions ended with a divergence" "$LOG" | grep -Eo "[0-9.]+%" | head -1); DIV=${DIV:-0%}
TD=$(grep -Eo "\(([0-9.]+)%\) transitions hit the maximum treedepth" "$LOG" | grep -Eo "[0-9.]+%" | head -1); TD=${TD:-0%}
RH=$(grep -Eo "worst R-hat [0-9.]+" "$LOG" | head -1 | awk '{print $3}')
ESS=$(grep -Eo "min bulk ESS [0-9.]+" "$LOG" | head -1 | awk '{print $4}')
CONV="NOT CONVERGED"
awk -v d="${DIV%\%}" -v r="${RH:-9}" 'BEGIN{exit !(d < 1 && r < 1.01)}' && CONV="converged"
NOTE="${LABEL}. ${CONV}: ${DIV} divergent, ${TD} at max treedepth, worst R-hat ${RH:-?}, min ESS ${ESS:-?}"
echo "$(date) drawing figure: $NOTE" >> "$OUT"
Rscript src/33_age_element_figure.R "$TAG" "$PRE" "$NOTE" "$MODEL" >> "$OUT" 2>&1
echo "$(date) figure script exit $?" >> "$OUT"
