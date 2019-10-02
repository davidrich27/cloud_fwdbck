INPUT="src/"
OUTPUT="bin/"

FILES=$(ls $INPUT)
COMPILER="gcc"
FLAGS="-std=c99 -Wall"

for FILE in $FILES
do
   echo "compiling file: $FILE..."
   FILENAME=$(basename -- "$FILE")
   NAME="${FILENAME%.*}"
   EXT="${FILENAME##*.}"
   
   INFILE=$INPUT/$FILENAME
   OUTFILE=$OUTPUT/$NAME

   if [ $EXT == "c" ]; then
      echo "$COMPILER $FLAGS -o $OUTFILE $INFILE"
      $COMPILER $FLAGS -o $OUTFILE $INFILE
   fi

done
