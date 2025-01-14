#!/bin/sh

if [ $# -lt 1 ]
then
  echo "Usage: $0 <old dc2offset.d> [<shift>]" 1>&2
  exit 1
fi

MERGE=./mergedcoffset
TMP=bd1offset$$.tmp
OFFSET_DIR=.
if [ $# -eq 1 ]
then
  $MERGE bc 201 $1 $OFFSET_DIR/dc_layer1_offset_1.d > $TMP.1
  $MERGE bc 202 $TMP.1 $OFFSET_DIR/dc_layer2_offset_1.d > $TMP.2
  $MERGE bc 203 $TMP.2 $OFFSET_DIR/dc_layer3_offset_1.d
else
  echo $2
  $MERGE bc 201 $1      $2 > $TMP.1
  $MERGE bc 202 $TMP.1  $2 > $TMP.2
  $MERGE bc 203 $TMP.2  $2
fi

rm $TMP.1 $TMP.2 $TMP.3 
