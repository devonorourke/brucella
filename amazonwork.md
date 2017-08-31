# installs
## core stuff
```
sudo apt-get update
sudo apt-get -y upgrade
sudo apt-get -y install tmux git curl gcc make g++ unzip default-jre libdatetime-perl libxml-simple-perl libdigest-md5-perl bioperl
```
## prokka dependencies
mkdir hmmer
cd hmmer/
wget http://eddylab.org/software/hmmer3/3.1b2/hmmer-3.1b2-linux-intel-x86_64.tar.gz
tar xvzf hmmer-3.1b2-linux-intel-x86_64.tar.gz
cd hmmer-3.1b2-linux-intel-x86_64.tar.gz
sudo ./configure
make
sudo make install

## prokka itself
cd $HOME
git clone https://github.com/tseemann/prokka.git
cd prokka/bin/
./prokka --setupdb

# actually run prokka
Create a shell script as follows with all your fasta files within that directory:
```
#!/bin/bash
cd /home/ubuntu/bruce/fastas
SAMPLE_LIST=`ls *.fasta | sed 's/.fasta//' | sort -u`
for SAMPLE in $SAMPLE_LIST
do
	prokka \
	--outdir ${SAMPLE} \
  --force \
	--prefix ${SAMPLE} \
	--genus Brucella \
	--usegenus \
	--gram neg \
	--cpus 0 \
	${SAMPLE}.fasta
done

prokka \

--outdir --outdir
