FROM ubuntu:18.04

MAINTAINER vferat <victor.ferat@fcbg.ch>

# Add user
RUN adduser --quiet --disabled-password qtuser

# Install Python 3, PyQt5
RUN apt-get update
RUN apt-get install -y \
      python3 \
python3-pyqt5

ADD . /

RUN apt-get install python3-pip -y
RUN python3 -m pip install -r requirements.txt


CMD ["python3", "ica_gui.py"]
