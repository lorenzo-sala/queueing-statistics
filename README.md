# Queueing Statistics Simulator

A discrete-event simulation of single-server queueing systems written in C,
with statistical analysis of performance metrics and confidence interval estimation.

## Overview

The simulator models arrival and service processes in a queue, computing
steady-state statistics and validating results against theoretical values
derived from queueing theory (M/M/1 and variants).

## Features

- Next-event simulation engine (`NESssq.c`)
- Estimation of key metrics: mean waiting time, queue length, server utilization
- Confidence interval computation via batch means
- Validation against analytical solutions
- Output visualization with Python + matplotlib

## Structure

| File | Description |
|------|-------------|
| `NESssq.c` | Core next-event simulation engine |
| `es1.c` / `es2.c` | Exercise implementations |
| `es1_validate.c` / `es2_validate.c` | Validation against theoretical values |
| `plot.py` / `plot2.py` | Confidence interval visualization |
| `validate*.bat` | Batch scripts for automated validation |

## Build & Run
```bash
gcc -o sim es1.c -lm
./sim
```

## Tools

`C` `Python` `matplotlib` `queueing theory`
