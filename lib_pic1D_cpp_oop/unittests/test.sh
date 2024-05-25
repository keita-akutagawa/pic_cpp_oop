#!/bin/bash

options="-lgtest -lgtest_main -pthread"
programfile="../*.cpp"
testfiles=$(find . -name "test_*.cpp")

for testfile in $testfiles; do
    g++ $testfile $constfile $programfile $options
    ./a.out
done


