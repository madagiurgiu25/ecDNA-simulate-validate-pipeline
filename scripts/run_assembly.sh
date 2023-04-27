#!/bin/bash

input=$1

shasta --input $input --config Nanopore-Oct2021 --Reads.minReadLength 6000
