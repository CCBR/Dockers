version=$1
if [ "$#" != "1" ];then
        exit "Version number argument missing!"
fi
docker push nciccbr/ccbr_ubuntu_base_20.04:$version
docker tag nciccbr/ccbr_ubuntu_base_20.04:$version nciccbr/ccbr_ubuntu_base_20.04:latest
docker push nciccbr/ccbr_ubuntu_base_20.04:latest

