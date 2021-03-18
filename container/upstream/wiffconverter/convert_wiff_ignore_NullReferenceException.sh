#!/bin/bash
/usr/bin/mono /usr/local/bin/sciex/wiffconverter/OneOmics.WiffConverter.exe \
    WIFF \
        $1 \
        -profile \
    MZML \
        $2 \
        --singleprecision \
        --nocompression \
        --overwrite 2> error.log || true
end=$(tail -n 1 "$2")
idx="</indexedmzML>"
reg="</mzML>"
if [[ "${end}" == "${idx}"  ]] || [[ "${end}" == "${reg}"  ]] ;
then 
	echo "success"
	exit 0; 
else 
	echo "fail"
	exit 2;
fi
