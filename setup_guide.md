## Setup guide for force field optimization 

### **Step 1.** Install Anaconda or Miniconda with Python 3 and create a new conda environment.
```
conda create -n openff
conda activate openff
```

### **Step 2.** Install packages using the Conda package manager.
```
conda install -c conda-forge notebook matplotlib qcportal geometric 
conda install -c openeye openeye-toolkits
conda install -c omnia -c conda-forge forcebalance openforcefield=0.4.1 cmiles
```
* Open Forcefield Toolkit version 0.4.1 was used in release-1 fitting.
* The OpenEye toolkit requires a valid license. 
* QCPortal version 0.11 was used in the release package, but we install the latest version here because that is required for connecting to the public server.


### **Step 3.** Install development version of ForceBalance from source. (Pending the new release of ForceBalance, this step will be removed)
```
conda install numpy scipy networkx lxml
conda install -c omnia pymbar
git clone https://github.com/leeping/forcebalance.git
python setup.py install
```
**Note**: If you want to reproduce the fitting, you need to get the same version of source code used for the fitting. The release-1-RC2 used the version commited in August 14th, so to use the same version,  you need to checkout and reinstall with this commands: 
```
git checkout 5b3a65d
python setup.py install
```
`5b3a65d` is the first 7 characters for the commit.
 
### (Optional) **Step 4.** Install and usage of Work Queue from CCTools 
#### 4.1. Installation of Work Queue
ForceBalance provides an option to use Work Queue for distributed calculations.
Work Queue is a distributed computing library that is developed by the [Cooperative Computing Lab](http://ccl.cse.nd.edu/) at Notre Dame; its documentation can be found [here](http://ccl.cse.nd.edu/software/manuals/workqueue.html). 

Here’s a brief introduction of Work Queue extracted from the linked manual:
```
“...Work queue is a framework for building large scale master-worker applications. Using the Work Queue library, you create a custom master program that defines and submits a large number of small tasks. Each task is distributed to a remote worker process which executes it and returns the results. As results are created, the master may generate more tasks to be executed...”
```

A handy bash script written for an automatic compilation of CCTools can be found [here](https://github.com/lpwgroup/torsiondrive/blob/master/devtools/travis-ci/install-cctools.sh). This script will download, compile, and install the executables to `$HOME/opt/cctools/`, and also install the Python binding into your current Python environment. To use this script:
```
wget https://github.com/lpwgroup/torsiondrive/blob/master/devtools/travis-ci/install-cctools.sh
bash install-cctools.sh
```

#### 4.2. Brief introduction of how to use Work Queue in ForceBalance
(a) Run ForceBalance with the `work_queue` option and provide a port number. 
In the input file (e.g. optimize.in) of ForceBalance, we have the following lines in the `$options` block:

```
wq_port <port number>
asynchronous True 
```

The port number should be a high four-digit number (for example, 9571).

In each `$target` block, we use this keyword to specify the remote evaluation of the target.:
```
remote 1
```
Then running ForceBalance
```
ForceBalance optimize.in
```
will enter a “host” mode that waits for worker connections. You can see a line like this printed on the bottom of the screen every few seconds:
```
<time> : 0/0 workers busy; 0/2018 jobs complete
```
(b) Launch worker process in remote nodes and connect to host. 
On a remote machine, first repeat above steps 1-3 to create an identical environment.
Make sure that the correct environment variables are set, including the Conda environment.
Start a worker process by running the following:
```
~/opt/cctools/bin/work_queue_worker -t 800000 <host_address> <port number>
```

The `<host_address>` can be IP address or host domain name of the host machine. The `<port number>` should be consistent with what you specified in the input file. `-t 800000` option tells the worker to not shut down when host connection is lost. To print all debug information you can add `-d all` to the command. 

If successfully connected, you will see your host ForceBalance program start printing a non-zero number of workers connected. To monitor the currently connected workers, you can also use this command on the host machine:
```
~/opt/cctools/bin/work_queue_status <host_address> <port number> -W
```

#### 4.3. Troubleshooting common Work Queue issues

Here we list the potential issues that may be encountered when running distributed ForceBalance calculations and how to resolve them. 
We will continually add to this list over time.

#### 4.3.1. Using SSH tunneling to work around firewalls for distributed computing

You may find that your distributed computing efforts are being stymied because your host machine (running the ForceBalance master) or the worker machines (running Work Queue workers) are not able to connect to each other.
This may happen because the host machine and/or the worker machine are behind firewalls that disallow outgoing and/or incoming connections.
Your host machine may also have a dynamic IP address making it difficult for workers to connect to.

It is possible to work around such problems with SSH tunneling, as long as you are to *somehow* connect from the host machine to the worker machine using SSH - even if you need to go through an intermediary machine such as the head node of a HPC cluster. 
The following example will assume that the ForceBalance master is running on a laptop with a dynamic IP address, and the workers are running on compute nodes that are part of a HPC cluster with a login node that can be reached at cluster.ucdavis.edu (this cluster does not actually exist)

__Step 1: Create a port forwarding SSH connection from host machine (laptop) to intermediary (cluster head node)__

On your laptop, run the following command to make a SSH connection to the cluster:
```
ssh -N -f -o ServerAliveInterval=30 -R9572:localhost:9571 cluster.ucdavis.edu
```
With the arguments `-N -f`, you do not get a prompt and the SSH process goes into the background. 

The option `-o ServerAliveInterval=30` keeps the connection from being dropped due to idle timeouts.

The option `-R9572:localhost:9571` means the following: Any connection made to host `localhost` port 9572 on the remote machine (the head node) will be forwarded to host `localhost` port 9571 on the local machine (your laptop). 

__Step 2: Create a port forwarding SSH connection from worker machine (compute node) to intermediary (head node)__

Make a normal SSH connection to the cluster compute node where you intend to run the Work Queue worker, then run the following command on the compute node:
```
ssh -N -f -o ServerAliveInterval=30 -L9573:localhost:9572 cluster.ucdavis.edu
```
The option `-L9573:localhost:9572` means the following: Any connection made to host `localhost` port `9573` on the local machine (the compute node) will be forwarded to host `localhost` port `9572` on the remote machine (the head node).

__Step 3: Run Work Queue Worker__

This allows your Work Queue workers to connect to your laptop, which is listening on port 9571, by directing the workers to connect to `localhost` port `9573` as follows:  

```
~/opt/cctools/bin/work_queue_worker -t 800000 localhost 9573
```

You may put the SSH commands into your job submission script to automate parts of this process.

__Note 1__: You could confirm the tunneling principle works (and build your own understanding) by using `-R9572:localhost:22` in the first SSH command above (from laptop to cluster), then create a separate "normal" SSH connection to the cluster and run `ssh localhost -p 9572` on the cluster. 
This is equivalent to making an incoming SSH connection to your laptop, and the SSH server listens for SSH connections on port `22`.
(Your laptop may not be configured to be running an SSH server, but you could enable one.
This is a separate concept from the Work Queue server and is not needed to get Work Queue working.)

Going one step further, keep the `-R9572:localhost:22` connection active from laptop to head node, then make a connection from the compute node to the head node using `-L9573:localhost:9572` as described above. 
Then you may run `ssh localhost -p 9573` on the compute node to make a SSH connection from the compute node to your laptop.

__Note 2__: After you understand how everything works, you can change the `9573` and `9572` numbers above to `9571`. 
The reason for using different numbers in the tutorial is that it would be confusing to have to distinguish "`localhost` port `9571`" on the laptop, head node, and compute node.

__Note 3__: Some cluster environments will not allow you to directly make SSH connections to compute nodes, in which case you can only make port forwarding connections from the job script.