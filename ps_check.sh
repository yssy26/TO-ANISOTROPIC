#!/bin/bash
echo "=== python / MTO_HF processes ==="
ps -eo pid,pcpu,etime,rss,comm | grep -iE 'python|MTO_HF' | grep -v grep
echo "=== solve script tail (if any output file) ==="
echo "(solve output is on stdout of the background job)"
