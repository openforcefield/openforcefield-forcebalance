## Setup guide for force field optimization 

### **Step 1.** Install anaconda 3 and create a new conda environment.
```
conda create -n <env name>
conda activate openff
```

### **Step 2.** Install packages using the Conda package manager.
```
conda install -c conda-forge notebook matplotlib qcportal geometric 
conda install -c openeye openeye-toolkits
conda install -c omnia -c conda-forge forcebalance openforcefield=0.4.1 cmiles
```
* Open Forcefield Toolkit version 0.4.1 was used in release-1 fitting.
* Openeye toolkits requires a valid license. 
* qcportal 0.11 was used, but we install the latest version because only the latest version is allowed to connect to the public server.


### **Step 3.** Install dev. version of ForceBalance  from source. (Pending the new release of ForceBalance, this step will be removed)
```
git clone https://github.com/leeping/forcebalance.git
python setup.py install
```
**Note**: If you want to reproduce the fitting, you need to get the same version of source code used for the fitting. The release-1-RC2 used the version commited in August 14th, so to use the same version,  you need to checkout and reinstall with this commands: 
```
git checkout 5b3a65d
python setup.py install
```
`5b3a65d` is the first 7 characters for the commit.
 

### (Optional) **Step 4.** Install and usage of cctools.work_queue 
#### 4.1. Installation of cctools.work_queue
ForceBalance provides an option to use Work Queue for distributing the target evaluations.  Here is the manual: http://ccl.cse.nd.edu/software/manuals/workqueue.html. Here’s a brief introduction of Work Queue extracted from the linked manual:
```
“...Work queue is a framework for building large scale master-worker applications. Using the Work Queue library, you create a custom master program that defines and submits a large number of small tasks. Each task is distributed to a remote worker process which executes it and returns the results. As results are created, the master may generate more tasks to be executed...”
```

A handy bash script written for an  automatic compilation of CCTools can be found [here](https://github.com/lpwgroup/torsiondrive/blob/master/devtools/travis-ci/install-cctools.sh). This script will download, compile, and install the executables to `$home/opt/cctools/` as well as the python binding in your current python environment. To use this script:
```
wget https://github.com/lpwgroup/torsiondrive/blob/master/devtools/travis-ci/install-cctools.sh
bash install-cctools.sh
```

#### 4.2. Brief introduction of how to use work queue in ForceBalance
(a) Run ForceBalance with work_queue and listen to a port. In the input file (e.g. optimize.in) of ForceBalance, we have the  following lines in the `$options` block:

```
wq_port <port number>
asynchronous True 
```

And in each `$target` block, we use this keyword to specify the remote evaluation of the target.:
```
remote 1
```
Then running ForceBalance
```
ForceBalance optimize.in
```
will enter a “host” mode that waits for worker connections. You can see a line like this printed on the bottom of the screen every few seconds:
```
<time> : 0/0 workers busy; 0/0 jobs complete
```
(b) Launch worker process in remote nodes and connect to host. 
On a remote machine, first repeat above steps 1-3 to create an identical environment.

Then you can start a worker process by running
```
~/opt/cctools/bin/work_queue_worker -t 800000 <host_address> <port number>
```

The `<host_address>` can be IP address or host domain name of the host machine. The `<port number>` should be consistent with what you specified in the input file. `-t 800000` option tells the worker to not shut down when host connection is lost. To print all debug information you can add `-d all` to the command. 

If successfully connected, you will see your host ForceBalance program start printing a non-zero number of workers connected. To monitor the currently connected workers, you can also use this command on the host machine:
```
~/opt/cctools/bin/work_queue_status <host_address> <port number> -W
```