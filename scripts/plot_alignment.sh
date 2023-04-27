#!/bin/bash

# assembly genome
assembly=$1

# circular reference
refcircular=$2

# reference genome (hg38)
refgenome=$3

# prefix
prefix=$4

cd /home/mada/Projects/tools/mummer-4.0.0rc1

./nucmer --forward --mum --batch=1000 --nooptimize --nosimplify --prefix ${prefix}/ref-ref ${refcircular} ${refgenome}
#./show-coords -r ${prefix}/ref-ref.delta > ${prefix}/ref-ref.coords
./mummerplot --coords -b --color -png -p ${prefix}/ref-ref ${prefix}/ref-ref.delta

./nucmer --forward --mum --batch=1000 --nooptimize --nosimplify --prefix ${prefix}/asm-ref ${assembly} ${refgenome}
#./show-coords -r ${prefix}/asm-ref.delta > ${prefix}/asm-ref.coords
./mummerplot --coords -b --color -png -p ${prefix}/asm-ref ${prefix}/asm-ref.delta


