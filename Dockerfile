FROM debian:unstable-slim
RUN apt-get -y update && apt-get install make gcc -y
ADD . /home
WORKDIR /home
RUN make
CMD ["./hbv"]
