#!/bin/bash
name=m1$(hostname)
for i in 1 2 3 4 5 6 7 8 9 10 11 12
do
    echo "$i Threads:" >> $name
    /usr/bin/time -f%e -ao $name ./partdiff-posix $i 2 512 2 2 1024
    echo "" >> $name
done
