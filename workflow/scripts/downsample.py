from os import walk, system
from os.path import join, basename


def GetPaths(directory): 
    """
    Function to get path to all fastq files in subdirectories.
    """
    # Path list: 
    Paths = list()
    
    # Walk the directory
    for r,d,f in walk(directory): 
        
        # Iterate over the files
        for file in f:
            
            # Check if files are fastq
            if file.endswith(".fastq"):
                
                # Save the paths to list. 
                Paths.append(join(r, file))
    
    return Paths

def downsample(directory):
    """
    This function downsamples fastq files to specific fractions using seqtk.
    """
    
    # Get a list of paths
    paths = GetPaths(directory)
    
    # CSV output file
    outcsv = open("./config/samples_downsampled.csv", "w")
    outcsv.write("Run,Path,Replicate,SpID,avg_ct,tot_raw_reads,mapped_reads,gisaid_short_name,gisaid_accession,sra_biosample,ON_TGT,Isolate,LibraryLayout,Virus,Host,SRA,Source,Percent\n")
    
    # How much to downsample
    percents = [25,50,75]
    
    # Iterate over each downsample percentage
    for percent in percents:
    
        # Iterate over the filenames
        for file in paths: 

            # Extract the name and append % downsampled
            name = basename(file).split(".")[0] + f"-{percent}"

            # Save the name of the output file
            outfile = join(directory, name, f"{name}.fastq")

            # Save the mkdir command
            mkdir = f"mkdir -p {join(directory, name)}"

            # Save the seqtk command
            seqtk = f"seqtk sample {file} {percent/100} > {outfile}"

            # Save to CSV
            
            metadata = name.split("-")
            
            PATH = f"[{outfile}]"
            
            row = f"{name},{PATH},{metadata[1]},{metadata[0]},NA,NA,NA,NA,NA,NA,NA,NA,SINGLE,SARS2,human,NA,local,{metadata[2]}\n"
            
            outcsv.write(row)
            
            # Run the commands 
            system(mkdir)
            system(seqtk)
            
    outcsv.close()


if __name__ == "__main__":

    downsample("/fh/fast/bloom_j/computational_notebooks/whannon/2020/SARS-CoV-2_bottleneck/results/sra/")

