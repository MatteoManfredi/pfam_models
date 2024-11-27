# pfam_models

This repository contains scripts for analyzing protein models and domains, specifically using PfamScan and FoldSeek. It allows users to process protein models, analyze them using Pfam domains, and generate results that can be used in further analysis or visualization.

## Setup Instructions

### 1. Clone the Repository

To get started, clone the repository to your local machine:

```bash
git clone https://github.com/yourusername/pfam_models.git
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

Create a new environment called `pfam_models`, which includes the required libraries:

```bash
conda create -n pfam_models seaborn
```

Activate the environment:

```bash
conda activate pfam_models
```

### 4. Install FoldSeek

To install FoldSeek, which is needed for protein structural alignments, use the following command:

```bash
conda install -c conda-forge -c bioconda foldseek
```

### 5. Update libtiff

After installing FoldSeek, you may need to update the `libtiff` library because FoldSeek can break the dependencies of seaborn. To fix this, run:

```bash
conda update libtiff
```

### 6. Install and Set Up PfamScan

To use PfamScan, follow these steps:

#### 6.1 Download PfamScan

First, download and extract the PfamScan package:

```bash
tar zxvf PfamScan.tar.gz
cd PfamScan
```

#### 6.2 Install HMMER3

PfamScan requires HMMER3 to perform domain searches. Install HMMER3 with Conda:

```bash
conda install -c biocore hmmer
```

#### 6.3 Add HMMER3 to Your Path

Make sure HMMER3 binaries are available in your shell’s executable path. If you're using bash, run:

```bash
export PATH=/path/to/install/hmmer3:$PATH
```

For `csh` or `tcsh`, run:

```bash
setenv PATH /path/to/install/hmmer3:$PATH
```

#### 6.4 Install Perl Dependencies

PfamScan requires several non-standard Perl dependencies. You can install them using Conda:

```bash
conda install -c bioconda perl-moose
conda install -c bioconda perl-bioperl
```

Alternatively, you can install them via the `cpan` tool if you prefer.

#### 6.5 Add PfamScan Modules to PERL5LIB

If you are using bash, add PfamScan modules to your `PERL5LIB`:

```bash
export PERL5LIB=/path/to/pfam_scanDir:$PERL5LIB
```

For `csh` or `tcsh`, run:

```bash
setenv PERL5LIB /path/to/pfam_scanDir:$PERL5LIB
```

#### 6.6 Change the Perl Path in `pfam_scan.pl`

Open `PfamScan/pfam_scan.pl` in a text editor and modify the first line to point to your Perl installation. For example:

```perl
#!/path/to/your/perl
```

Now, PfamScan is set up and ready to use.

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
  -f, --foldseek        Set this flag if you already generated the FoldSeek results.
  -a, --aggregate       Set this flag if you already generated the TSV file.
```

### Example Command:

To start an analysis using the provided `enzymes.txt` file as input, you can run the script with the following command:

```bash
python analyze_models.py -n "project_name" -w 4 --pfamscript PfamScan/pfam_scan.pl --pfamdb PfamScan/pfam/
```

- `-n "project_name"`: Sets the project name prefix (default is "enzymes").
- `-w 4`: Sets the maximum number of parallel workers (adjust based on your machine's capabilities).
- `--pfamscript`: Path to the PfamScan script (`pfam_scan.pl`).
- `--pfamdb`: Path to the PfamScan database directory.
- `-d`: Set this flag if you already downloaded the necessary files from the webserver.
- `-s`: Set this flag if you already generated the PfamScan results.
- `-f`: Set this flag if you already generated the FoldSeek results.
- `-a`: Set this flag if you already generated the TSV file.

The script will perform the analysis, using FoldSeek for structural alignment and PfamScan for domain analysis, and will generate the results in the specified output directory.

---

## Additional Information

- **Dependencies**: In addition to the libraries mentioned above, this script requires `seaborn`, `pandas`, and other related Python libraries.
- **Input File**: The `enzymes.txt` file included in the repository is used as an input to start an analysis. You can modify this file or provide your own.
