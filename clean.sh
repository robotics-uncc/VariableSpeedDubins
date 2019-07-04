# remove the build and binary directories
rm -rf build
rm -rf bin
rm -rf lib

# remove from main directory any hidden files and other temporaries 
# (e.g. octave, python files)
rm -rf `find . -name "*~" -print`
rm -f *.m
rm -f octave*
rm -f *.out
rm -f *.py
rm -f *.lp


