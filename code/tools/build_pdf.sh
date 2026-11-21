#!/bin/bash
# 提交版报告 MD → PDF(pandoc + xelatex,中文排版)。
# 用法: ./tools/build_pdf.sh   → 产出 submission/报告.pdf
set -e
cd "$(dirname "$0")/../submission"
pandoc 报告.md -o 报告.pdf \
  --pdf-engine=xelatex \
  --resource-path=. \
  -V documentclass=article \
  -V geometry:margin=2.4cm \
  -V fontsize=11pt \
  -V CJKmainfont="Songti SC" \
  -V mainfont="Times New Roman" \
  -V linestretch=1.25 \
  -V linkcolor=blue \
  --toc --toc-depth=2 \
  -V toc-title="目录" \
  -H ../tools/pdf_preamble.tex
echo "== 已生成 submission/报告.pdf =="
