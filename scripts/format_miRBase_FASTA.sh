#!/bin/bash

grep -A 1 "^>hsa" data/mature.fa | sed -e '/^>/!s/U/T/g' -e '/^--$/d' > data/mature_hsa.fa