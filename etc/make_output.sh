#!/usr/bin/env bash
#------------------------------------------------------------------------------#
#
#  File:        make_output.sh
#
#  Description: Formats the output of the make file commands
#
#------------------------------------------------------------------------------#

# format terminal output
OK="\033[01;32mok\033[0m"
WARN="\033[01;33mok\033[0m"
FAIL="\033[01;31mfail\033[0m"

# print command to log file
echo " " >> $2
echo $1 >> $2

cmd="$1 2>&1"

output=`eval $cmd`

# check if command succeeded or failed
if [ $? -eq 0 ] 
then

# check if output is created
	if [ -z "$output" ]
	then

# print ok
		printf "$OK\n"

	else

# print warning
		printf "$WARN\n"
		echo "$output" >> $2

	fi

else

# print error
	printf "$FAIL\n"
	echo "$output" >> $2
	printf "\n\033[01;31m%s\033[0m (check %s for details)\n\n%s\n\n" "$3" "$2" "$output"
	exit 1

fi
