FROM ubuntu:18.04

RUN apt-get update && apt-get install -y curl

#ENV VERSION=$(curl  "https://api.github.com/repos/cli/cli/releases/latest" 2>/dev/null|grep tag_name|awk '{print $NF}'|sed "s/[\",]//g")
ENV VERSION=1.7.0
RUN mkdir -p /opt2 && mkdir -p /data2

WORKDIR /opt2
RUN curl -sSL https://github.com/cli/cli/releases/download/v${VERSION}/gh_${VERSION}_linux_amd64.tar.gz -o gh_${VERSION}_linux_amd64.tar.gz && tar xvf gh_${VERSION}_linux_amd64.tar.gz && cp gh_${VERSION}_linux_amd64/bin/gh /usr/local/bin/ && cp -r gh_${VERSION}_linux_amd64/share/man/man1/* /usr/share/man/man1/

MAINTAINER vishal.koparde@nih.gov
