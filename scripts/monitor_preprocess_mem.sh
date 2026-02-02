#!/usr/bin/env bash
# 持续监控预处理 Python 进程的 RSS 与运行时间，每 10 秒采样，直到进程结束。
# 用法: ./scripts/monitor_preprocess_mem.sh [PYTHON_PID]
# 若不传 PID，会尝试从当前 make preprocess 进程树中查找 python3。

set -e
PY="${1:-}"
if [[ -z "$PY" ]]; then
  MAKE_PID=$(pgrep -f "make preprocess.*11M" | head -1)
  if [[ -n "$MAKE_PID" ]]; then
    # make -> sh -> uv -> python3
    for pid in $(pgrep -P "$MAKE_PID"); do
      for ppid in $(pgrep -P "$pid"); do
        if ps -o cmd= -p "$ppid" 2>/dev/null | grep -q "atlasmap_preprocess.cli"; then
          PY=$ppid
          break 2
        fi
      done
    done
  fi
  if [[ -z "$PY" ]]; then
    echo "未找到预处理进程。请指定 PID: $0 <python_pid>"
    exit 1
  fi
  echo "自动检测到 Python PID: $PY"
fi

LOG="${LOG:-preprocess_mem.log}"
echo "监控 PID=$PY，每 10 秒写入 $LOG（内存 + 进程运行时间），直到进程结束..."
echo "--- 开始 $(date '+%Y-%m-%d %H:%M:%S') ---" >> "$LOG"

while kill -0 "$PY" 2>/dev/null; do
  ts=$(date '+%Y-%m-%d %H:%M:%S')
  # 进程已运行时间 (格式 [[dd-]hh:]mm:ss)
  etime=$(ps -o etime= -p "$PY" 2>/dev/null | tr -d ' ')
  rss_kb=$(ps -o rss= -p "$PY" 2>/dev/null | tr -d ' ')
  if [[ -n "$rss_kb" ]]; then
    rss_mb=$(echo "scale=1; $rss_kb/1024" | bc 2>/dev/null)
    rss_gb=$(echo "scale=2; $rss_kb/1024/1024" | bc 2>/dev/null)
    line="[$ts] 运行时间=${etime:-?} | PID=$PY RSS=${rss_mb} MB (${rss_gb} GB)"
    echo "$line" >> "$LOG"
    echo "$line"
  fi
  sleep 10
done

# 进程结束时无法再取 etime，用当前时间与监控开始时间差作参考（需 /proc 起止时间则需在循环内记录）
echo "--- 进程已结束 $(date '+%Y-%m-%d %H:%M:%S') ---" >> "$LOG"
echo "进程已结束: $(date '+%Y-%m-%d %H:%M:%S')"
