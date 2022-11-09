#!/bin/bash

rem() {
  ARRAY=$1

  unset dupes # ensure it's empty
  unset NEWARRAY
  declare -A dupes

  for i in "${ARRAY[@]}"; do
      if [[ -z ${dupes[$i]} ]]; then
          NEWARRAY+=("$i")
      fi
      dupes["$i"]=1
  done
  unset dupes # optional
}

ARR=(aa ab bb aa ab cc)

printf "[%s]" "${ARR[@]}"
echo

rem $ARR

ARR=$NEWARRAY

ARR=("${NEWARRAY[@]}")

printf "[%s]" "${ARR[@]}"
echo
