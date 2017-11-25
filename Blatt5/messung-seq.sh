#!/bin/bash
name=m1$(hostname)-seq
/usr/bin/time -f%e -ao $name ./partdiff-seq 1 2 512 2 2 1024
