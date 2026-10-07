#!/usr/bin/env bash
# UserPromptSubmit hook: append the owner's prompt, with local date and time,
# to AI_PROMPTS.md at the repository root.  Reads the hook JSON on stdin.
set -euo pipefail

root="${CLAUDE_PROJECT_DIR:-$(cd "$(dirname "$0")/../.." && pwd)}"
log="$root/AI_PROMPTS.md"
input="$(cat)"
prompt="$(printf '%s' "$input" | jq -r '.prompt // empty')"
[ -n "$prompt" ] || exit 0

if [ ! -f "$log" ]; then
  printf '# AI prompt log\n\nEvery prompt given to the AI assistant on this project, newest last.\nAppended automatically by `.claude/hooks/log_prompt.sh` (UserPromptSubmit hook).\n' > "$log"
fi

{
  printf '\n### %s · %s\n\n' "$(date '+%Y-%m-%d %H:%M:%S %Z')" "$(hostname -s)"
  printf '%s\n' "$prompt" | sed 's/^/> /'
} >> "$log"
