# ICA toolbox GUI

# Description

The  ICA toolbox GUI is a graphic interface designed to make easier the use of [MNE](https://martinos.org/mne/stable/index.html) Independent Component Analysis (ICA) capabilities .

## Installation

Use the package manager [pip](https://pip.pypa.io/en/stable/) to install required depencies:

```bash
pip install -r requirements.txt
```
Or  build a docker image from the docker file:

```bash
docker build -t ica_toolbox .
```

## Usage

### Using Python
To launch the app using from your python installation use:
```bash
python run.py
```
### Using docker

### Linux
First you need to allow docker to access you screen/
```bash
sudo xhost +local:docker
```

Or use the docker image you created:
```bash
sudo docker run -ti --rm -e DISPLAY=$DISPLAY -v /tmp/.X11-unix:/tmp/.X11-unix icatoolbox
```
Use  ``` -v -v HOST_FOLDER/:/DOCKER_FOLDER ``` to share a folder between the host and the docker image and therefore import / export files from the application.
```bash
sudo docker run -ti --rm -e DISPLAY=$DISPLAY -v /tmp/.X11-unix:/tmp/.X11-unix -v HOST_FOLDER/:/DOCKER_FOLDER icatoolbox
```

### Windows

Natively, a Linux based docker image is not able to pop graphic interface on windows host.
One solution is to install a third party software [Xlaunch](http://www.straightrunning.com/XmingNotes/)  to cast GUI from docker to windows.   After completing the installation , you must start Xlaunch and configure it. Select parameters that best suit you ( or keep default) and make sure to check Disable access control on the extra settings menu.
Once done, you should be able to see the GUI by using the following command specified
```bash
docker run -e DISPLAY=MY_IP:0.0 icatoolbox
```
where MY_IP is your local IP ( use ```ipconfig``` on your command line to known your current ip )

Use  ``` -v -v HOST_FOLDER/:/DOCKER_FOLDER ``` to share a folder between the host and the docker image and therefore import / export files from the application.
```bash
docker run -e DISPLAY=MY_IP:0.0 -v HOST_FOLDER/:/DOCKER_FOLDER icatoolbox
```

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.


## License
GNU Lesser General Public License v3.0
