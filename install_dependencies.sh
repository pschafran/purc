#!/bin/bash
echo "setting variables..."
sed -i "" 's/SYSOS.*/SYSOS=\"MACOS\"/' purc.py

cd Dependencies

echo 'making sure muscle is executable...'
chmod +x muscle3.8.31

echo 'making sure vsearch is executable...'
ln -sf vsearch-2.15.1-macos-x86_64/bin/vsearch vsearch
chmod +x vsearch-2.15.1-macos-x86_64/bin/vsearch

echo 'install DADA2'
Rscript install_dada2.R

echo 'compiling cutadapt...'
cd cutadapt_source
python setup.py build_ext -i


echo 'done.'