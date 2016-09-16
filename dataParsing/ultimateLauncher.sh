#!/bin/bash

echo "-----------------------------------------------"
echo "hi, $USER"

./launcher.sh small 1 1000 1 1000
./launcher.sh small 1 1000 2 1000
./launcher.sh small 1 1000 3 1000
./launcher.sh med 1 1000 1 1000
./launcher.sh med 1 1000 2 1000
./launcher.sh med 1 1000 3 1000
./launcher.sh large 1 1000 1 1000
./launcher.sh large 1 1000 2
./launcher.sh large 1 1000 3