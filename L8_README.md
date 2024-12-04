[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-22041afd0340ce965d47ae6ef1cefeee28c7c493a6346c4f15d667ab976d596c.svg)](https://classroom.github.com/a/4PT5VAVn)
# Lab 8
You should run through this tutorial on the machine you intend to use for the final project. We provide instructions for both Amazon Web Services (AWS) and hosting on your personal laptop or desktop computer. Any option is acceptable and will be a valuable learning experience - you can even use another hosting platform, or you can run on a local installation of Linux. We recommend using AWS so you can familiarize yourself with it as a platform and so that any issues you encounter are as easy as possible for us to assist with - different Apple hardware and different Windows Subsystem for Linux versions that have different results, whereas every Amazon t2.micro instance running ubuntu 24.04 should be the same. 

This lab will walk you through the JBrowse2 Genome Browser tutorial found in the paper linked below. You should read through as there are many helpful comments, but note that in the tutorial version in this lab we provide some additional context, platform-specific assistance for non-linux machines, and some simplifications and corrections for Linux installs on AWS.

https://currentprotocols.onlinelibrary.wiley.com/doi/10.1002/cpz1.1120

At the end, you will have a little interactive web app that lets you browse annotations in the human genome. The same sort of browser can easily be made for all kinds of organisms, so long as appropriate data are available! Your submission will be using this genome browser web app to view a gene annotation of your choice, as described in the submission section below and in assignment.md. 

**For all command line steps below, make sure to check that you are using the commands for the right platform, and carefully run the command blocks (i.e. the copyable blocks of code) one-by-one in order. This will reduce the chance of creating an unexpected state that could become very difficult to troubleshoot.**

## 1. Platform specific setup
It should be possible to run the tutorial on either a macOS computer (Intel or M1), a Windows computer (Windows 10 or later), or Amazon Web Services (AWS). However, there are a few platform specific setup details: specifically, some of the necessary bioinformatics tools (like samtools) are only available pre-compiled for Unix operating systems, making Windows Subsystem for Linux an important tool. If you want to use a Windows machine for bioinformatics, this is going to be a crucial tool for many things. Please follow the platform-specific instructions below before proceeding to the "Install necessary tools" section.

In terms of simplicity, the order is probably macOS (local unix) simpler than AWS (cloud instance unix) simpler than Windows (non-unix requiring a virtual machine), but in terms of just smoothness of installs (after setting up your instance) AWS should be the simplest, because you'll be running our steps on the exact same type of system where they were tested. 

### 1.1. Mac OS setup
Open a terminal and run the line below to install homebrew, a macOS package manager. This will make it easy for you to install necessary packages like apache2 and samtools. _You can skip this step if you already have brew installed._ 

```
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```
If this doesn't work, visit https://docs.brew.sh/Installation for further installation options, including a .pkg installer that should be convenient and easy to use.

### 1.2. Windows setup
Enable and set up Windows Subsystem for Linux, using the default Ubuntu distribution. _You can skip these steps if you already have WSL set up with a Debian or Ubuntu distribution._

For newer versions of Windows, this command should handle it for you. Further details can be found at https://learn.microsoft.com/en-us/windows/wsl/install. You should restart your computer after the install is done. 

```
# this command installs the default linux distribution for your Windows OS, which should be an acceptable Ubuntu version
wsl --install
```

You can check that `wsl` was installed properly by running `wsl -l -v`. If you are running an older Windows 10 version, you may instead need to follow the instructions here: https://learn.microsoft.com/en-us/windows/wsl/install-manual. (In this case, you have to select your Ubuntu version. We recommend Ubuntu 22.04.)

Next, set up your Linux username and password. You can launch WSL the first time from the start menu by searching Ubuntu or you can use `windows key+r`, type `wsl` and press `enter`. Once launched, follow the prompts (see https://learn.microsoft.com/en-us/windows/wsl/setup/environment#set-up-your-linux-username-and-password). Make sure to record the password you choose, although when you launch WSL in future it should automatically sign you in.

For all subsequent steps, run from within the WSL virtual machine. You should be able to start wsl after initial setup by typing `wsl` in the command line shell or by using the start menu. This way, you should be able to seamlessly run Unix applications and use the Linux instructions in subsequent steps. You will need to install homebrew in step 1.4. in order to get samtools and htslib (which includes tabix) in step 2.3. 

### 1.3. AWS Setup

Follow the separate AWS [setup guide](./aws_instructions.md), then return here to set up linuxbrew below.

### 1.4. Linuxbrew for WSL or AWS
Make sure you are using a Debian or Ubuntu distribution. Then go ahead and install linuxbrew, using the instructions below:

switch to root with:

`sudo su -`
Then run:

`passwd ubuntu`
It is going to prompt :

`Enter new UNIX password:`

Set your password to something you can remember for later, or write down. A common password choice is simply `ubuntu` - not very secure at all, but AWS accounts themselves can be made fairly secure. 

Exit root by typing `exit`. **Note: it is important to exit root, because you do not want to accidentally run future commands with administrator privileges when that might be undesirable. The subsequent command in this case will fail if run from root.**

Install brew using the bash script from https://brew.sh/. You will be prompted to set the password you made earlier.
```
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

After this is complete, add brew to your execution path:
```
echo >> /home/ubuntu/.bashrc
echo 'eval "$(/home/linuxbrew/.linuxbrew/bin/brew shellenv)"' >> /home/ubuntu/.bashrc
eval "$(/home/linuxbrew/.linuxbrew/bin/brew shellenv)"
```

## 2. Install necessary tools
### 2.1. Node.js
Node.js is a cross-platform JavaScript runtime environment that will make is easy to run JBrowse2 command-line tools.

First, check whether Node.js is already installed by running the following. If node v20 is already installed, you can skip to the next step.

```
node -v
```

If Node.js is not installed, install it. 

#### macOS

On macOS, you can use brew. You may need to restart the terminal (close and open a new one) to get `node -v` to run.

```
# NOTE:
# Homebrew is not a Node.js package manager.
# Please ensure it is already installed on your system.
# Follow official instructions at https://brew.sh/
# Homebrew only supports installing major Node.js versions and might not support the latest Node.js version from the 20 release line.
# download and install Node.js
brew install node@20
# verifies the right Node.js version is in the environment
node -v # should print `v20.18.0`
# verifies the right npm version is in the environment
npm -v # should print `10.8.2`
```

#### Linux

For Linux, you can use the code below. See https://nodejs.org/en/download/package-manager for more detail.

On AWS and some other Linux setups, you may need to run `sudo apt install unzip` first.

**Note: `sudo`, also known as "super user do", runs commands with root/admin privileges. This can cause harm to your machine if you run the wrong command! It is also, however, a critical tool when doing things like installs - if you try something and are denied due to permissions, `sudo` is often the solution.**

```
# installs fnm (Fast Node Manager)
curl -fsSL https://fnm.vercel.app/install | bash
# activate fnm
source ~/.bashrc
# download and install Node.js
fnm use --install-if-missing 20
# verifies the right Node.js version is in the environment
node -v # should print `v20.18.0`
# verifies the right npm version is in the environment
npm -v # should print `10.8.2`
```

### 2.2. @jbrowse/cli
Run the following commands in your shell. This uses the Node.js package manager to download the latest stable version of the jbrowse command line tool, then prints out its version. This should work for both macOS and Linux.

```
sudo npm install -g @jbrowse/cli
jbrowse --version
```

You can also try installing using just `npm install -g @jbrowse/cli` if the sudo version doesn't run. 

### 2.3. System dependencies
Install wget (if not already installed), apache2, samtools, and tabix. 

wget is a tool for retrieving files over widely-used Internet protocols like HTTP and FTP. 

apache2 allows you to run a web server on your machine.

samtools and tabix, as we have learned earlier in the course, are tools for processing and indexing genome and genome annotation files.

#### macOS

```
# note that apache2 gets installed as httpd for macOS, which is the service you will launch later
brew install wget httpd samtools htslib
```

#### Linux

```
sudo apt install wget apache2
```

```
brew install samtools htslib
```

## 3. Apache server setup
### 3.1. Start the apache2 server
Starting up the web server will provide a localhost page to show that apache2 is installed and working correctly. When discussing computer networking, localhost is a hostname that refers to the current computer used to access the network. Note that in WSL2, the linux subsystem may have a different IP address from your Windows OS, and so you will want to use that IP address to be able to find it and load the web page. AWS, on the other hand, will have a public IP address that you need to identify in the aws_instructions.

#### macOS
```
sudo brew services start httpd
```

#### Linux
```
sudo service apache2 start
```

### 3.2. Getting the host
If you are running locally on your mac, the hostname is just `localhost`. However, for WSL and AWS, you will need to do a bit of work to find the right ip address.
For local hosting, the url will be `http://localhost:8080/` or `http://XX.XXX.XXX.XX:8080/`, where Xs are replaced with the appropriate IP address from the WSL steps below.

#### WSL
```
# from within WSL, run the linux server launch command to launch the service, then print out you WSL IP address so you can access the server from your Windows browser
# if the ip command isn't recognized, install iproute and then try again
# sudo apt install iproute2
ip addr show eth0 | grep "inet\b" | awk '{print $2}' | cut -d/ -f1
```
This should give you an ip address you can use to access the web server.

#### AWS
In your instance summary page, there should be an "auto-assigned IP address." Your web server can be accessed at `http://ipaddress`. You don't need to provide a port.

### 3.3. Access the web server
Open a browser and type the appropriate url into the address bar. You should then get to a page that says "**It works!**" (for AWS there may be some additional info). If you have trouble accessing the server, you can try checking your firewall settings and disabling any VPNs or proxies to make sure traffic to localhost is allowed.

### 3.4. Verify apache2 server folder
Apache2 web servers serve files from within a root directory. This is configurable in the httpd.conf configuration file, but you shouldn't have to change it (in fact, changing the conf file is not recommended unless you know what you are doing). 

For a normal linux installation, the folder should be `/var/www` or `/var/www/html`, whereas when you install on macOS using brew it will likely be in `/opt/homebrew/var/www` (for M1) or `/usr/local/var/www` (for Intel). You can run `brew --prefix` to get the brew install location, and then from there it is in the `var/www` folder. 

Verify that one of these folders exists (it should currently be empty, except possibly for an index file, but we will now populate it with JBrowse 2). If you have e.g. a www folder with no www/html folder, and your web server is showing the "It works!" message, you can assume that the www one is the root directory. 

Take note of what the folder is, and use the command below to store it as a command-line variable. We can reference this variable in the rest of our code, to save on typing. You will need to re-run the `export` if you restart your terminal session!
```
# be sure to replace the path with your actual true path!
export APACHE_ROOT='/path/to/rootdir'
```

If you are really struggling to find the APACHE_ROOT folder, you could try searching for it.
```
sudo find / -name "www" 2>/dev/null
```

### 3.5. Download JBrowse 2
First create a temporary working directory as a staging area. You can use any folder you want, but moving forward we are assuming you created ~/tmp in your home folder.

```
mkdir ∼/tmp
```
```
cd ∼/tmp
```

Next, download and copy over JBrowse 2 into the apache2 root dir, setting the owner to the current user with `chown` and printing out the version number. This version doesn't have to match the command-line jbrowse version, but it should be a version that makes sense.

```
jbrowse create output_folder
sudo mv output_folder $APACHE_ROOT/jbrowse2
sudo chown -R $(whoami) $APACHE_ROOT/jbrowse2
```

### 3.6. Test your jbrowse install
In your browser, now type in `http://yourhost/jbrowse2/`, where yourhost is either localhost or the IP address from earlier. Now you should see the words "**It worked!**" with a green box underneath saying "JBrowse 2 is installed." with some additional details. 

## 4. Load and process test data
### 4.1. Download and process reference genome
Make sure you are in the temporary folder you created, then download the human genome in fasta format. This is the biggest file you'll be downloading, and may take 30 min or so on AWS with the lowest tier of download speeds.

```
export FASTA_ROOT=https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens
wget $FASTA_ROOT/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
```

Unzip the gzipped reference genome, rename it, and index it. This will allow jbrowse to rapidly access any part of the reference just by coordinate.

```
gunzip Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
mv Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa hg38.fa
samtools faidx hg38.fa
```

### 4.3. Load genome into jbrowse

```
jbrowse add-assembly hg38.fa --out $APACHE_ROOT/jbrowse2 --load copy
```

## 4.4. Download and process genome annotations

Still in the temporary folder, download ENSEMBLE genome annotations in the GFF3 format. 

```
export GFF_ROOT=https://ftp.ensembl.org/pub/release-110/gff3/homo_sapiens
wget $GFF_ROOT/Homo_sapiens.GRCh38.110.chr.gff3.gz
gunzip Homo_sapiens.GRCh38.110.chr.gff3.gz
```

Use jbrowse to sort the annotations. jbrowse sort-gff sorts the GFF3 by refName (first column) and start position (fourth column), while making sure to preserve the header lines at the top of the file (which start with “#”). We then compress the GFF with bgzip (block gzip, which zips files into little blocks for rapid access), and index with tabix. The tabix command outputs a file named genes.gff.gz.tbi in the same directory, and we then refer to “genes.gff.gz” as a “tabix indexed GFF3 file”.

```
jbrowse sort-gff Homo_sapiens.GRCh38.110.chr.gff3 > genes.gff
bgzip genes.gff
tabix genes.gff.gz
```

### 4.5. Load annotation track into jbrowse

```
jbrowse add-track genes.gff.gz --out $APACHE_ROOT/jbrowse2 --load copy
```

### 4.6. Index for search-by-gene

Run the “jbrowse text-index” command to allow users to search by gene name within JBrowse 2.

In the temporary work directory, run the following command.

```
jbrowse text-index --out $APACHE_ROOT/jbrowse2
```

## 5.0 Use your genome browser to explore a gene of interest
### 5.1. Launch JBrowse2
Open `http://yourhost/jbrowse2/` again in your web browser. There should now be several options in the main menu. Follow the guide in the "Launch the JBrowse 2 application and search for a gene in the linear genome view" section of https://currentprotocols.onlinelibrary.wiley.com/doi/10.1002/cpz1.1120 to navigate to the gene search and try browsing a few genes.

### 5.1. Explore a specific gene of interest
Using the gene search function you've set up, look up a gene with disease-relevant function. You can use any academic articles you like to identify a gene of interest - one example might be to look up a cancer repressor gene, or a gene from your metabolic pathway lab earlier in the semester that might have mutations associated with metabolic diseases. 

In assignment.md, answer the prompts.

# Submission
Commit your assignment.md file with the prompts answered, then push to your remote repository and provide the link to your repo in the bCourses submission form.
