# pfam_models

This repository allows reproducibility of the analysis described in the manuscript *"AlphaFold2 and ESMFold: A Large-Scale Pairwise Model Comparison of Human Enzymes upon Pfam Functional Annotation."* 

We provide a Python script that automates the following tasks:  
1. **Model Retrieval**: Downloads AlphaFold2 and ESMFold structural models from the [Alpha&ESMhFolds](https://alpha-esmhfolds.biocomp.unibo.it/) webserver for a given list of UniProt IDs.  
2. **Domain Annotation**: Runs PfamScan and Foldseek to annotate Pfam domains on the retrieved models.  
3. **Statistical Analysis**: Computes statistical metrics to assess the quality of the models on the Pfam regions.

The file `enzymes.txt` includes the list of enzymes used in the manuscript's analysis. However, the script can be applied to any custom list of UniProt IDs.

---

## Setup Instructions

### 1. Clone the Repository

To get started, clone the repository to your local machine:

```bash
git clone https://github.com/MatteoManfredi/pfam_models.git
cd pfam_models
```

### 2. Set Up Conda (If Not Already Installed)

If you don’t have Conda installed, you can install it by following these steps:

- Download and install Miniconda (a lightweight version of Anaconda) from [here](https://docs.conda.io/en/latest/miniconda.html).
- After installation, restart your terminal, or follow the instructions provided after installation to initialize Conda in your shell.

To verify Conda is installed correctly, run:

```bash
conda --version
```

### 3. Create a New Conda Environment

Create a new environment called `pfam_models`, including the seaborn library:

```bash
conda create -n pfam_models python=3.7
```

Activate the environment:

```bash
conda activate pfam_models
```

### 4. Install and Set Up PfamScan

To use [PfamScan](https://pfam-docs.readthedocs.io/en/latest/), follow these steps:

#### 4.1 Download PfamScan

First, download and extract the PfamScan package:

```bash
wget https://ftp.ebi.ac.uk/pub/databases/Pfam/Tools/PfamScan.tar.gz
tar -xzvf PfamScan.tar.gz
rm PfamScan.tar.gz
cd PfamScan
```

#### 4.2 Install PfamScan dependencies

PfamScan requires HMMER3 and several non-standard Perl dependencies. You can install them using Conda:

```bash
conda install -c biocore hmmer
```
```bash
conda install -c bioconda perl-moose
```
```bash
conda install -c bioconda perl-bioperl
```

#### 4.3 Add PfamScan Modules to PERL5LIB

Add PfamScan modules to your `PERL5LIB`:

```bash
export PERL5LIB=$(pwd):$PERL5LIB
```

#### 4.4 Download Pfam data files

PfamScan requires several files to execute, including `Pfam-A.hmm`, `Pfam-A.hmm.dat`, and `active_site.dat`. You can download them from the [Pfam FTP site](ftp://ftp.ebi.ac.uk/pub/databases/Pfam/) under `releases/Pfam37.0`:

```bash
mkdir pfam
cd pfam
wget https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam37.0/Pfam-A.hmm.gz
wget https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam37.0/Pfam-A.hmm.dat.gz
wget https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam37.0/active_site.dat.gz
gunzip Pfam-A.hmm.gz Pfam-A.hmm.dat.gz active_site.dat.gz
hmmpress Pfam-A.hmm
cd ../../
```

`hmmpress` is needed to generate binary files for `Pfam-A.hmm`

Now, PfamScan is set up and ready to use.

### 5. Install Foldseek

To install [Foldseek](https://github.com/steineggerlab/foldseek), which is needed for protein structural alignments, use the following command:

```bash
conda install -c conda-forge -c bioconda foldseek
```

### 6. Install Python libraries

To install Python libraries used by the script, use the following command:

```bash
conda install seaborn=0.12
```

---

## Running the Script

Once your environment is set up, you can run the script to analyze protein models and domains. The script supports the following command-line arguments:

```bash
usage: analyze_models.py [-h] [-n NAME] [-w MAX_WORKERS] [--pfamscript PFAMSCRIPT]
                         [--pfamdb PFAMDB] [-d] [-s] [-f] [-a]

Analyze protein models and domains.

optional arguments:
  -h, --help            show this help message and exit
  -n NAME, --name NAME  Project name prefix (default: enzymes)
  -w MAX_WORKERS, --max-workers MAX_WORKERS
                        Maximum parallel workers (default: 1)
  --pfamscript PFAMSCRIPT
                        Path to PfamScan script (default: PfamScan/pfam_scan.pl)
  --pfamdb PFAMDB       Path to PfamScan database directory (default: PfamScan/pfam/)
  -d, --download        Set this flag if you already downloaded the files from the webserver.
  -s, --pfamscan        Set this flag if you already generated the PfamScan results.
  -f, --foldseek        Set this flag if you already generated the Foldseek results.
  -a, --aggregate       Set this flag if you already generated the TSV file.
```

### Reproduce the manuscript's analysis:

To reproduce the analysis using the provided `enzymes.txt` file as input, you can run the script with the following command:

```bash
python script.py
```

This assumes that inside the CWD you have a directory called `PfamScan` containing the script `pfam_scan.pl` and the subdirectory `pfam` with the downloaded data files. Moreover, files will be downloaded and analyzed one at a time. Alternatively, you can run the following command:

```bash
python script.py -w 4 --pfamscript path/to/script/pfam_scan.pl --pfamdb path/to/data/files/
```

- `-w 4`: Sets the maximum number of parallel workers (adjust based on your machine's capabilities).
- `--pfamscript`: Path to the PfamScan script (`pfam_scan.pl`).
- `--pfamdb`: Path to the PfamScan database directory.

### Running a new analysis:

To run the same analysis on a different dataset using a file `project_name.txt` as input, you can run the script with an additional flag that will set the project name prefix:

```bash
python script.py -n "project_name"
```

---

## Additional Information

### Execution Time and Disk Space Requirements  
The analysis process is resource-intensive and may take a significant amount of time to complete. It also generates numerous intermediate files, which require adequate disk space. Ensure your system has sufficient free space and that the execution is not interrupted until the analysis is complete.  

If you need to halt the process midway, you can resume it later by using the `-d -s -f -a` flags. These flags allow you to skip already completed steps in the pipeline, saving time during re-execution.

### Foldseek Memory Considerations  
When running Foldseek on multiple inputs in parallel, insufficient memory may cause the program to terminate unexpectedly. To mitigate this, the `max_workers` parameter is capped at 1 for the Foldseek step. However, Foldseek will still utilize all available CPU cores to process each input efficiently.  

If you have a machine with ample memory, consider commenting or changing line 269 of `script.py` to optimize performance for your specific hardware.
