FROM continuumio/miniconda3
ARG BRANCH=master

LABEL maintainer="lasse.krueger94@gmail.com"

RUN apt-get update \
&& apt-get install -yq git unzip vim curl wget build-essential

COPY . /mb_pipeline

#SHELL ["/bin/bash", "--login", "-c"]

WORKDIR /mb_pipeline
CMD ["bash", "/run_pipeline.sh"]
