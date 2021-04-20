version=$1
if [ "$#" != "1" ];then
        exit "Version number argument missing!"
fi
docker build -f Dockerfile.$version -t nciccbr/ccbr_atacseq:$version . |tee build.${version}.log
