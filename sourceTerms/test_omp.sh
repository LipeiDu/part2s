sh cleanMake.sh

mkdir output_8
mkdir output_4
mkdir output_1

sh runMP.sh 8
mv output/* output_8/

sh runMP.sh 4
mv output/* output_4/

sh runMP.sh 1
mv output/* output_1/