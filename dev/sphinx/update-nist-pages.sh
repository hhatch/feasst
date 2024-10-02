#!/bin/bash

# update nist-pages
mkdir build
cd build
cmake -DUSE_SPHINX=ON ..
make html > tt 2>&1
grep -v "_arguments.rst: WARNING: document" tt | grep -v "_arguments.rst:4: WARNING: Duplicate"  | grep -v "^Declaration is" | grep -v "WARNING: Duplicate C++ declaration, also defined"
# now run automatically: python ../dev/tools/depend.py -s ../
# /feasst/dev/analyze_public_interface.py
# update the version manually in CMakeLists
# remove debugging compiler flags in CMakeLists
# add feasst.h and remove depend.py from CMakeLists
# tag the main with this same version before building html
version=$(git describe)
#branch=`git branch | grep \* | cut -d ' ' -f2`
mv html html2
git checkout nist-pages
cp -r html2/* ../
git status
echo "Press [Enter] to commit with message: $version?"
read -rs
git commit -a -m "$version"
#git checkout $branch
# also dont forget to upload to pypi https://packaging.python.org/en/latest/tutorials/packaging-projects/
