#!/bin/bash

# Quick system load check for benchmarking
echo "🖥️  SYSTEM LOAD CHECK"
echo "===================="
echo "Current time: $(date)"
echo ""

echo "📊 CPU & Load:"
echo "Load average: $(uptime | awk -F'load average:' '{print $2}')"
echo "CPU usage: $(top -bn1 | grep "Cpu(s)" | awk '{print $2 $3 $4 $5}')"
echo ""

echo "💾 Memory:"
free -h | head -2
echo ""

echo "💽 Disk I/O (iostat if available):"
if command -v iostat >/dev/null 2>&1; then
    iostat -x 1 1 | tail -n +4 | head -5
else
    echo "iostat not available - install sysstat package for detailed I/O stats"
fi
echo ""

echo "🔥 Top I/O processes:"
if command -v iotop >/dev/null 2>&1; then
    echo "Use: sudo iotop -o -d 1 -n 3"
else
    echo "iotop not available - install iotop package"
    echo "Top CPU/memory processes:"
    ps aux --sort=-%cpu | head -6
fi

echo ""
echo "💡 Benchmark recommendations:"
if (( $(echo "$(uptime | awk '{print $10}' | cut -d',' -f1) > 2.0" | bc -l) )); then
    echo "⚠️  High load detected - results may be inconsistent"
    echo "   Consider running during off-peak hours"
else
    echo "✅ Load looks reasonable for benchmarking"
fi