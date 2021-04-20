version=$1
if [ "$#" != "1" ];then
	exit "Version number argument missing!"
fi
docker push nciccbr/ccbr_atacseq:$version
