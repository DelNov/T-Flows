#!/bin/bash

# Created the beginning of the test script
echo "#!/bin/bash"           >  worker.sh
echo ""                      >> worker.sh
echo "cd ../Sources/Process" >> worker.sh
echo ""                      >> worker.sh
echo "old_directory=$PWD"    >> worker.sh
echo ""                      >> worker.sh

# Find all directories with user function and append commands
# for their compilation to the beginning of the test script
find . -name "User_Mod" \
| sed 's/\.\//make clean; make DIR_CASE=..\/..\/Tests\//g' \
| sed 's/\/User_Mod/\nif \[ -f ..\/..\/Binaries\/Process \]; then echo "SUCCESS"; else echo "FAILURE"; fi/g' \
>> worker.sh

echo "make clean"        >> worker.sh
echo ""                  >> worker.sh
echo "cd $old_directory" >> worker.sh

# Make the test script executable
chmod 755 worker.sh

source ./worker.sh > test_user_functions_compile.$(date +%y-%m-%d-%T).log

cd $old_directory

/bin/rm worker.sh
