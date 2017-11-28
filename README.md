Introduction
=============
A collection of tools to work with continuous data typically obtained from electrophysiological experiments.

Line noise removal
------------------
Remove line noise at 50 Hz from a signal `Y` sampled at 30 kHz.

```julia
Y2, pp = LFPTools.remove_linenoise(Y, 50.0, 30_000.0)
```
