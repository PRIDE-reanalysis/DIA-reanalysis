#!/bin/bash
if [ -z "$4" ]
  then
    SwaMe.Console -i $1 \
              -o $2 \
              --irtFile $3 \
              -t true \
              --maxThreads 6 \
              --maxQueueSize 60 2> error.log || true
  else
    SwaMe.Console -i $1 \
              -o $2 \
              --irtFile $3 \
              --tempFolder $4 \
              -z \
              -t true \
              --maxThreads 6 \
              --maxQueueSize 60 2> error.log || true
fi
python -mjson.tool $2 > /dev/null
