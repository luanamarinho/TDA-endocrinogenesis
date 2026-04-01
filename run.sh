#!/bin/bash
# Run the Giotto-TDA pipeline via seqerakit and monitor until completion.

set -euo pipefail

WORKSPACE="LUANALAB/TDA-Project"
SEQERAKIT_YML="${1:-seqerakit.yml}"

# ── 1. Launch via seqerakit ──────────────────────────────────
# Generate a unique run name from the current timestamp so the YAML
# never needs to be manually bumped between runs.
export PIPELINE_RUN_NAME="giotto-tda-$(date +%Y%m%d-%H%M%S)"
echo "▶ Launching pipeline as '$PIPELINE_RUN_NAME'..."
envsubst '$PIPELINE_RUN_NAME' < "$SEQERAKIT_YML" | seqerakit /dev/stdin

# ── 2. Resolve the run ID from the latest submitted run ──────
echo "▶ Resolving run ID..."
RUN_ID=$(tw runs list -w "$WORKSPACE" --format json \
    | jq -r '[.[] | select(.status == "SUBMITTED" or .status == "RUNNING")] | first | .id')

if [[ -z "$RUN_ID" || "$RUN_ID" == "null" ]]; then
    echo "✗ Could not resolve a running workflow ID. Exiting."
    exit 1
fi
echo "  Run ID: $RUN_ID"

# ── 3. Poll until terminal state ─────────────────────────────
echo "▶ Monitoring run $RUN_ID..."
while true; do
    STATUS=$(tw runs view -w "$WORKSPACE" "$RUN_ID" --format json | jq -r '.status')
    echo "  Status: $STATUS"
    case "$STATUS" in
        SUCCEEDED)
            echo "✓ Pipeline completed successfully."
            exit 0
            ;;
        FAILED|CANCELLED|UNKNOWN)
            echo "✗ Pipeline $STATUS."
            exit 1
            ;;
    esac
    sleep 30
done
