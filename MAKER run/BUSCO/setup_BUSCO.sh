#get BUSCO
git clone https://gitlab.com/ezlab/busco.git
cd busco
sudo python setup.py install
export PATH="/opt/augustus-3.2.2/scripts:$PATH"

#get the dataset
wget https://busco.ezlab.org/datasets/metazoa_odb9.tar.gz
tar xvf metazoa_odb9.tar.gz

#install hmmer
cd ..
wget http://eddylab.org/software/hmmer/hmmer.tar.gz
tar xvf hmmer.tar.gz
cd hmmer
./configure
make

cd ..
cd busco/config
curl https://raw.githubusercontent.com/warnerlab/LP_annotation/master/MAKER%20run/BUSCO/config.ini > config.ini

sudo apt-get install libparallel-forkmanager-perl

cd ..
mkdir augustus_config

cp -r /home/scijake/busco/augustus_config
export AUGUSTUS_CONFIG_PATH="/home/scijake/busco/augustus_config/config/"