__author__ = "Will Hannon"
__copyright__ = "Copyright 2020 Will Hannon"
__email__ = "wwh22@uw.edu"
__license__ = "MIT"

## === Import Libraries === ##
import datetime
import ftplib
import os
import tarfile
import pandas as pd
import yaml

print(f'Preparing for upload at {datetime.datetime.now()}')

# Load in config file.
with open('config/sra/sra_config.yml') as config_yaml:
    config = yaml.safe_load(config_yaml)

# Load in fastq data frame.
fastqs = pd.read_csv('config/sra/fastq_files_to_upload.csv')

# Make the tar file for submission
tar_filename = 'SRA_submission.tar'

try:
    with tarfile.open(tar_filename, mode='w') as f:
        for i, tup in enumerate(fastqs.itertuples()):
            print(f"Adding file {i + 1} of {len(fastqs)} to {tar_filename}\n")
            f.add(tup.filename_fullpath, arcname=tup.filename)
        print(f"Added all files to {tar_filename}.")

except:
    if os.path.isfile(tar_filename):
        os.remove(tar_filename)
    raise

# Check that the file was correctly made
print(f"The size of {tar_filename} is {os.path.getsize(tar_filename) / 1e9:.1f} GB\n")

with tarfile.open(tar_filename) as f:
    files_in_tar = set(f.getnames())
if files_in_tar == set(fastqs['filename']):
    print(f"{tar_filename} contains all {len(files_in_tar)} expected files.\n")
else:
    raise ValueError(f"{tar_filename} does not have all the expected files.\n")

# Get ftp info from config file.
ftp_address = config['ftp_address'].strip()
ftp_username = config['ftp_username'].strip()
ftp_account_folder = config['ftp_account_folder'].strip()
ftp_subfolder = config['ftp_subfolder'].strip()

with open('config/sra/ftp_password.txt') as f:
    ftp_password = f.read().strip() # not tracked by the repo for privacy

# Start the upload
print(f"Starting upload at {datetime.datetime.now()}\n")

with ftplib.FTP(ftp_address) as ftp:
    ftp.login(user=ftp_username,
              passwd=ftp_password,
              )
    ftp.cwd(ftp_account_folder)
    ftp.mkd(ftp_subfolder)
    ftp.cwd(ftp_subfolder)
    with open(tar_filename, 'rb') as f:
        ftp.storbinary(f"STOR {tar_filename}", f)
        
print(f"Finished upload at {datetime.datetime.now()}")
