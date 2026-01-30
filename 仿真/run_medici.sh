#!/usr/bin/env bash
set -euo pipefail

# 列出要运行的文件（与工作区中的文件名一致）
# 自动检测当前目录下以 R_N 开头且不是 .out 的文件（保留小数点）
detect_files(){
  files=()
  while IFS= read -r fname; do
    files+=("$fname")
  done < <(find . -maxdepth 1 -type f -printf '%f\n' | grep -E '^R_N' | grep -v '\.out$' | sort)
  if [ ${#files[@]} -eq 0 ]; then
    echo "未在当前目录找到 R_N* 输入文件。请确认您在包含 R_N... 文件的目录下运行此脚本。"
    exit 1
  fi
}

# 默认顺序执行：逐个调用 medici（使用 detect_files 生成文件列表）
run_sequential(){
  detect_files
  echo "Found ${#files[@]} input files:"
  printf ' - %s\n' "${files[@]}"
  for f in "${files[@]}"; do
    echo "===== Running: medici $f (logging to ${f}.out) ====="
    medici "$f" > "${f}.out" 2>&1
    echo "===== Finished: $f ====="
  done
}

# 并行执行（如果系统安装了 GNU parallel）
run_parallel(){
  if command -v parallel >/dev/null 2>&1; then
    detect_files
    echo "Found ${#files[@]} input files:"
    printf ' - %s\n' "${files[@]}"
    printf "%s\n" "${files[@]}" | parallel -j ${PAR_JOBS:-4} 'medici {} > {}.out 2>&1'
  else
    echo "未检测到 'parallel'，将退回到顺序执行。要并行请安装 GNU parallel。"
    run_sequential
  fi
}

# 主入口：默认顺序，可用 --parallel 或 --jobs=N
if [[ ${1:-} == "--parallel" ]]; then
  # 可选第二个参数 --jobs=N 或环境变量 PAR_JOBS
  if [[ ${2:-} =~ ^--jobs=([0-9]+)$ ]]; then
    PAR_JOBS=${BASH_REMATCH[1]}
  fi
  run_parallel
else
  run_sequential
fi
