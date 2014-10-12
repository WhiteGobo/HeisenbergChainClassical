#!/bin/sh
< /dev/hwrng tr -dc 0-9 | head -c150000
