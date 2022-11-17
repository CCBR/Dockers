#!/usr/bin/env bash

if [ "$#" != "2" ];then
	echo "bash $0 requires 2 arguments!"
	echo "argument1: path to Dockerfile"
	echo "argument2: DockerHub tag"
	exit
fi

dt=$(date +"%d%m%y%H%M%S")

docker build \
	--build-arg BUILD_DATE=${dt} \
	-f $1 \
	-t $2 \
 	. | \
	tee $1_build_${dt}.log
