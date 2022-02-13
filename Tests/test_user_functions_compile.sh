#!/bin/bash

# Created the beginning of the test script
echo "#!/bin/bash"           >  worker.sh
echo ""                      >> worker.sh
echo "cd ../Sources/Process" >> worker.sh

# Find all directories with user function and append commands
# for their compilation to the beginning of the test script
find . -name "User_Mod" \
| sed 's/\.\//make clean; make DIR_CASE=..\/..\/Tests\//g' \
| sed 's/\/User_Mod/\nif \[ -f ..\/..\/Binaries\/Process \]; then echo "SUCCESS"; else echo "FAILURE"; fi/g' \
>> worker.sh

# Make the test script executable
chmod 755 worker.sh

source ./worker.sh > test_compile_user_functions.$(date +%y-%m-%d-%T).log

rm -f worker.sh
