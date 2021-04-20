version=$1
if [ "$#" != "1" ];then
        exit "Version number argument missing!"
fi
docker build -f Dockerfile.$version -t nciccbr/ccbr_ubuntu_base_20.04:$version . |tee build.${version}.log

