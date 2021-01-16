#!/bin/bash
cd Dependencies

echo 'making sure muscle is executable...'
rm muscle3.8.31
mv muscle3.8.31_i86linux64 muscle3.8.31
chmod +x muscle3.8.31

echo 'making sure vsearch is executable...'
ln -sf vsearch-2.15.1-linux-x86_64/bin/vsearch vsearch
chmod +x vsearch-2.15.1-linux-x86_64/bin/vsearch

echo 'compiling cutadapt...'
cd cutadapt_source
python setup.py build_ext -i

echo 'done.'