#!/bin/bash

echo "-----------------------------------------------"
echo "hi, $USER"


if [ $4 -eq 1 ]
	then
	topology="Star"
elif [ $4 -eq 2 ]
	then
	topology="Grid"
elif [ $4 -eq 3 ]
	then
	topology="Mesh"
fi

sizeTop=${1}${topology}
settings=${sizeTop}"/"${sizeTop}Settings
message=${sizeTop}"/"${sizeTop}Messages
output=${sizeTop}"/"${sizeTop}Log
low=$2
high=$3
limiter=$5

# for ((i=1;i<=input;i++)); do 
# 	while read -r line
# 	do
# 	    name="$line"
# 	    echo "Name read from file - $name"
# 	    ./a.out "$name"
# 	done < "$filename"	
# done

for ((i=$low;i<=$high;i++)); do
	if (( $i % $limiter ))
    	then
	    if [ $4 -eq 3 ]
		    then
		    echo ./a.out "$settings$i" "$message" "$output$i" &
		    ./a.out "$settings$i" "$message" "$output$i" &
	    else
		    echo ./a.out "$settings" "$message$i" "$output$i" &
		    ./a.out "$settings" "$message$i" "$output$i" &
		fi
    else
        if [ $4 -eq 3 ]
    	    then
    	    echo ./a.out "$settings$i" "$message" "$output$i"
    	    ./a.out "$settings$i" "$message" "$output$i"
        else
    	    echo ./a.out "$settings" "$message$i" "$output$i"
    	    ./a.out "$settings" "$message$i" "$output$i"
    	fi

fi
done
