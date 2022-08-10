#!/usr/bin/env bash

if [ "$#" != "2" ];then
	echo "bash $0 requires 2 arguments!"
	echo "argument1: path to Dockerfile"
	echo "argument2: DockerHub tag"
	exit
fi

dt=$(date +"%d%m%y%H%M%S")

docker build \
	-f $1 \
	-t $2 \
 	. | \
	tee $1_build_${dt}.log

exitcode=$_

if [ "$exitcode" == "0" ];then
docker push $2
fi
