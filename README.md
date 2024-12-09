# Influenza A H3N2 Analysis Dashboard

Welcome to our interactive dashboard, where we have provided various analyses and annotations of strains of Influenza A H3N2 (yearly flu affecting humans) between 2015 and 2019. The following instructions will walk you through the JBrowse2 Genome Browser set-up and how to use our specific features. Note that this dashboard is set up and hosted through AWS to maintain a standard host computing service since every Amazon t2.micro instance running Ubuntu 24.04 should be the same.

**For all command line steps below, make sure to run the command blocks one-by-one in order. This will reduce the chance of creating an unexpected state that could become very difficult to troubleshoot.**

## Platform Setup

### 1.1. AWS Setup

Follow the separate AWS [setup guide](./aws_instructions.md), then return here to set up linuxbrew below.

### 1.2. Linuxbrew for WSL or AWS
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
npm install -g @jbrowse/cli
jbrowse --version
```


### 2.3. System dependencies
Install wget (if not already installed), apache2, samtools, and tabix. 

wget is a tool for retrieving files over widely-used Internet protocols like HTTP and FTP. 

apache2 allows you to run a web server on your machine.

samtools and tabix, as we have learned earlier in the course, are tools for processing and indexing genome and genome annotation files.

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

#### Linux
```
sudo service apache2 start
```

### 3.2. Getting the host
If you are running locally on your mac, the hostname is just `localhost`. However, for WSL and AWS, you will need to do a bit of work to find the right ip address.
For local hosting, the url will be `http://localhost:8080/` or `http://XX.XXX.XXX.XX:8080/`, where Xs are replaced with the appropriate IP address from the WSL steps below.

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
export APACHE_ROOT='/var/www/html'
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
sudo mv output_folder/* $APACHE_ROOT/jbrowse2
sudo chown -R $(whoami) $APACHE_ROOT/jbrowse2
```


### 3.6. Test your jbrowse install
In your browser, now type in `http://yourhost/jbrowse2/`, where yourhost is either localhost or the IP address from earlier. Now you should see the words "**It worked!**" with a green box underneath saying "JBrowse 2 is installed." with some additional details. 

## 4. Load and process test data
Make sure you are in the temporary folder you created, then download and process the relevant Influenza A genome using the `process_genome.py` script. The instructions below detail how to fetch and set up the script in terminal and then include the arguments necessary to load and process the specific data.

1. Fetch Script (Raw File)
```
wget https://raw.githubusercontent.com/kristine-sun/final-project/refs/heads/main/process_genome.py
```
2. Set up python virtual environment
```
python3 -m venv myenv
source myenv/bin/activate
```
3. Install dependencies
```
pip install --upgrade pip
```
4. Make script executable
```
chmod +x process_genome.py
```
5. Run script

### nucleotide sequences
H1N1, April 30, 2009:
```
python3 process_genome.py \
    --fasta_url "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/038/517/505/GCA_038517505.1_ASM3851750v1/GCA_038517505.1_ASM3851750v1_genomic.fna.gz" \
    --gff_url "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/038/517/505/GCA_038517505.1_ASM3851750v1/GCA_038517505.1_ASM3851750v1_genomic.gff.gz" \
    --output_dir "/var/www/html/" \
    --strain "H1N1" \
    --date "2009_04_30"
```
H1N1, April 27, 2009:
```
python3 process_genome.py \
    --fasta_url "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/038/513/145/GCA_038513145.1_ASM3851314v1/GCA_038513145.1_ASM3851314v1_genomic.fna.gz" \
    --gff_url "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/038/513/145/GCA_038513145.1_ASM3851314v1/GCA_038513145.1_ASM3851314v1_genomic.gff.gz" \
    --output_dir "/var/www/html/" \
    --strain "H1N1" \
    --date "2009_04_27"
```
H3N2, January 10, 2013:
```
python3 process_genome.py \
    --fasta_url "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/038/671/225/GCA_038671225.1_ASM3867122v1/GCA_038671225.1_ASM3867122v1_genomic.fna.gz" \
    --gff_url "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/038/671/225/GCA_038671225.1_ASM3867122v1/GCA_038671225.1_ASM3867122v1_genomic.gff.gz" \
    --output_dir "/var/www/html/" \
    --strain "H3N2" \
    --date "2013_01_10"
```
H3N2, February 01, 2013:
```
python3 process_genome.py \
    --fasta_url "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/038/678/235/GCA_038678235.1_ASM3867823v1/GCA_038678235.1_ASM3867823v1_genomic.fna.gz" \
    --gff_url "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/038/678/235/GCA_038678235.1_ASM3867823v1/GCA_038678235.1_ASM3867823v1_genomic.gff.gz" \
    --output_dir "/var/www/html/" \
    --strain "H3N2" \
    --date "2013_02_01"
```
### Troubleshooting & Tips for loading and processing data
1. Sometimes, due to compressing/decompressing or another processing issue, the file size can be zero or very small. If this case happens, delete the file and download it again.

2. If the terminal says annotations are successfully added to track but are not showing up on the jbrowse site, it might be because the file is empty or the annotations don’t match the genome.


## Multiple Sequence Alignment
To run multiple sequence alignment and download the files to be loaded into the JBrowse site from scratch, follow the instructions in `generating_alignment.ipynb` included in the repository. However, in this case we have already generated the files so you can skip to the loading instructions below.

### Loading files to JBrowse2

### BED file
```
wget https://raw.githubusercontent.com/kristine-sun/final-project/refs/heads/main/cluster_alignment.bed
```
```
sudo mv ~/tmp/cluster_alignment.bed /var/www/html/jbrowse2/
```
jbrowse add-track ./cluster_alignment.bed --out $APACHE_ROOT/jbrowse2 --load copy --assemblyNames "H3N2_2013_01_10.fna" --force
jbrowse text-index --out $APACHE_ROOT/jbrowse2 --force

### WIG file
```
https://raw.githubusercontent.com/kristine-sun/final-project/refs/heads/main/cluster_alignment.wig
```
