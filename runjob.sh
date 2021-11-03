#!/bin/bash

echo "\n Testing matlab inputs :>"\

nohup opt/matlab/r2020a/bin/matlab -nodisplay -r  > OUT.1  < ~/Matlab/GaussMat/testbash.m
