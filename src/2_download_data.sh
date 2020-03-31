## clonning panglaodb database containing metadata
mkdir data/panglaodb
cd data/panglaodb
git clone https://github.com/oscar-franzen/PanglaoDB.git

## Downloading landscape data
## the file is in zip format
wget https://ndownloader.figshare.com/articles/7235471?private_link=2892d1cf91752e08a963 \
     -P data/landscape_data/
