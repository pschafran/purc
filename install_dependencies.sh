#!/bin/bash
echo "setting variables..."
cp purc.py purc.py.bak
sed 's/SYSOS=\".*\"/SYSOS=\"Darwin\"/' purc.py.bak > purc.py

cd Dependencies

echo 'making sure muscle is executable...'
ln -sf muscle3.8.31_i86darwin64 muscle3.8.31
chmod +x muscle3.8.31_i86darwin64

echo 'making sure vsearch is executable...'
ln -sf vsearch-2.15.1-macos-x86_64/bin/vsearch vsearch
chmod +x vsearch-2.15.1-macos-x86_64/bin/vsearch

echo 'install DADA2'
Rscript install_dada2.R

echo 'compiling cutadapt...'
cd cutadapt_source
python setup.py build_ext -i


echo 'done.'
