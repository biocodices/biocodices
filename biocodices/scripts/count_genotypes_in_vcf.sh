#!/bin/bash

grep -v '^#' $1 | ruby -F'\t' -lane 'puts $F[9].split(":").first' | sort | uniq -c | sort -k 1 -r
