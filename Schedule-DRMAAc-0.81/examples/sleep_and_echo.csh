#!/bin/csh

if( $1 != '' ) then
echo "sleep argv"
sleep $1
   if( $2 != '' ) then
      echo $2
   endif
   if( $3 == 'EXIT_ERROR' ) then
      exit 1;
   endif
   if( $3 == 'ABORT' ) then
      abort/abort;
      exit $?;
   endif
else
echo "sleep 1"
sleep 1
endif

