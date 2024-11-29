import argparse
import os
import subprocess
import shlex
import re
import json
from concurrent.futures import ProcessPoolExecutor
from typing import List, Tuple, Dict, TextIO

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt


# Helper functions
def parallel(max_workers: int, commands: List[List[str]]) -> None:
    """
    Run multiple shell commands in parallel.
    """
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        list(executor.map(run_command, commands))


def run_command(cmd: List[str]) -> None:
    """
    Runs a command handling one output redirection.
    """
    if len(cmd) < 3 or cmd[-2] not in [">", ">>"]:
        # Run the command
        subprocess.run(cmd, check=True)
    elif cmd[-1] == "/dev/null":
        # Run the command discarding the output
        subprocess.run(cmd[:-2], stdout=subprocess.DEVNULL, check=True)
    else:
        # Run the command redirecting the output
        with open(cmd[-1], "w" if cmd[-2] == ">" else "a") as outfile:
            subprocess.run(cmd[:-2], stdout=outfile, check=True)


def check_environment(name: str, max_workers: int, pfamscript: str, pfamdb: str, download: bool, pfamscan: bool, foldseek: bool, aggregate: bool) -> None:
    """
    Validate the environment before execution.
    """
    def check_file(file_path: str, desc: str):
        if not os.path.isfile(file_path):
            raise FileNotFoundError(f"{desc} '{file_path}' does not exist.")
    
    def check_directory(dir_path: str, desc: str):
        if not os.path.isdir(dir_path):
            raise FileNotFoundError(f"{desc} '{dir_path}' does not exist.")

    if aggregate:
        check_file(f"{name}.txt", "Input file with UniProt IDs")
        
        if not download:
            for subdir in ["afs", "esms", "jsons", "fastas"]:
                check_directory(f"{name}-{subdir}", f"Download directory {subdir}")
        
        if pfamscan:
            check_file(pfamscript, "PfamScan script")
            check_directory(pfamdb, "PfamScan database")
            required_files = ["Pfam-A.hmm", "Pfam-A.hmm.h3f", "Pfam-A.hmm.h3m", "Pfam-A.hmm.h3i", "Pfam-A.hmm.h3p", "Pfam-A.hmm.dat", "active_site.dat"]
            for req_file in required_files:
                check_file(os.path.join(pfamdb, req_file), f"Pfam database file {req_file}")
        else:
            check_directory(f"{name}-pfamscan", "PfamScan directory")
        
        if not foldseek:
            for subdir in ["foldseek", "fstmp"]:
                check_directory(f"{name}-{subdir}", f"Foldseek directory {subdir}")
    
    else:
        check_file(f"{name}.tsv", "Aggregated TSV file")

    print("All environment checks passed.")


# Compute Data
def read_ids(name: str) -> List[str]:
    """
    Read UniProt IDs from a file.
    """
    file_path = f"{name}.txt"
    with open(file_path) as reader:
        return [line.strip() for line in reader]


def download_files(name: str, uniprot_ids: List[str], max_workers: int) -> None:
    """
    Download AlphaFold2 and ESMFold models along with JSON metadata.
    """
    print("Downloading data from the webserver ...")
    
    output_dirs = [f"{name}-{suffix}" for suffix in ["afs", "esms", "jsons"]]
    for directory in output_dirs:
        os.makedirs(directory, exist_ok=True)

    commands = []
    for uniprot_id in uniprot_ids:
        commands.append(shlex.split(f"wget -q -O {name}-afs/{uniprot_id}.pdb https://alpha-esmhfolds.biocomp.unibo.it/get_single_model/af/{uniprot_id}"))
        commands.append(shlex.split(f"wget -q -O {name}-esms/{uniprot_id}.pdb https://alpha-esmhfolds.biocomp.unibo.it/get_single_model/esmf/{uniprot_id}"))
        commands.append(shlex.split(f"wget -q -O {name}-jsons/{uniprot_id}.json https://alpha-esmhfolds.biocomp.unibo.it/get_stats/{uniprot_id}"))
    
    parallel(max_workers, commands)
    
    print("... Data downloaded from the webserver.")


def check_downloads(name: str, uniprot_ids: List[str]) -> List[str]:
    '''
    Check if some of the IDs where not present in the webserver
    '''
    # Find missing IDs when the corresponding JSON file is empty
    missing = []
    for uniprot_id in uniprot_ids:
        json_file = f"{name}-jsons/{uniprot_id}.json"
        with open(json_file) as json_reader:
            json_data = json.load(json_reader)
        if not json_data:
            missing.append(uniprot_id)
    
    # If at least one ID is missing, print a warning, store the missing IDs and exclude them from further analysis
    if missing:
        print(f"{len(missing)} IDs were not present in Alpha&ESMhFolds and will be skipped from this analysis")
        missing_file = f"{name}-missing.txt"
        with open(missing_file, "w") as writer:
            writer.write("\n".join(missing)+"\n")
        found = [uniprot_id for uniprot_id in uniprot_ids if uniprot_id not in missing]
        original_file = f"{name}.txt"
        with open(original_file, "w") as writer:
            writer.write("\n".join(found)+"\n")
        return found
    
    # If all IDs were found do nothing
    return uniprot_ids


def extract_fastas(name: str, uniprot_ids: List[str]) -> None:
    """
    Extract FASTA files from the models PDB files.
    """
    print("Extracting FASTA files ...")
    
    aa_three_to_one = {
        "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
        "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
        "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
        "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
        # Handle unknown or uncommon residues
        "ASX": "B", "GLX": "Z", "XLE": "J", "SEC": "U", "PYL": "O",
        "UNK": "X"
    }
    
    # Create directory for storing the FASTA files
    directory = f"{name}-fastas"
    os.makedirs(directory, exist_ok=True)

    # Process each PDB file to extract the sequence
    for uniprot_id in uniprot_ids:
        input_pdb = os.path.join(f"{name}-afs", f"{uniprot_id}.pdb")
        output_fasta = os.path.join(f"{name}-fastas", f"{uniprot_id}.fasta")
        sequence = ""

        # Read the PDB file
        with open(input_pdb, "r") as pdb_reader:
            for line in pdb_reader:
                # Process ATOM lines with CA atoms
                if line.startswith("ATOM") and line[13:15].strip() == "CA":
                    residue_code = line[17:20].strip()
                    sequence += aa_three_to_one.get(residue_code, "X")  # Default to "X" for unknown codes

        # Write the sequence to the FASTA file
        with open(output_fasta, "w") as fasta_writer:
            fasta_writer.write(f">{uniprot_id}\n{sequence}\n")

    print("... FASTA files created.")


def run_pfamscan(name: str, uniprot_ids: List[str], max_workers: int, pfamscript: str, pfamdb: str) -> None:
    """
    Execute the PfamScan tool to identify domains and active sites.
    """
    print("Running PfamScan ...")
    
    directory = f"{name}-pfamscan"
    os.makedirs(directory, exist_ok=True)

    commands = []
    for uniprot_id in uniprot_ids:
        input_fasta = os.path.join(f"{name}-fastas", f"{uniprot_id}.fasta")
        output_pf = os.path.join(f"{name}-pfamscan", f"{uniprot_id}.pf")
        commands.append(shlex.split(f"{pfamscript} -fasta {input_fasta} -dir {pfamdb} -as > {output_pf}"))
    
    parallel(max_workers, commands)
    
    print("... PfamScan processing complete.")


def read_pfams(name: str, uniprot_ids: List[str]) -> List[Tuple[str, int, int, str, str, List[str]]]:
    """
    Read the output of the PfamScan tool and returns a list of Pfam entries.
    """
    pfam_entries = []
    
    for uniprot_id in uniprot_ids:
        pfam_file = os.path.join(f"{name}-pfamscan", f"{uniprot_id}.pf")

        with open(pfam_file, "r") as reader:
            # Skip the header lines starting with "#"
            for line in reader:
                if not line.startswith("#"):
                    break

            # Process the remaining lines for Pfam entry data
            for line in reader:
                columns = line.split()
                uniprot_id = columns[0]
                start = int(columns[1])
                end = int(columns[2])
                pfam_id = columns[5].split('.')[0]
                entry_type = columns[7]
                if "predicted_active_site" in columns[-1]:
                    active_sites = sorted(set(map(int, re.findall(r'\d+', columns[-1]))))
                else:
                    active_sites = []
                pfam_entries.append((uniprot_id, start, end, pfam_id, entry_type, active_sites))

    return pfam_entries


def cut_domains(name: str, pfam_entries: List[Tuple[str, int, int, str, str, List[str]]]) -> None:
    """
    Extract regions of AlphaFold2 and ESMFold models corresponding to Pfam entries and saves them in separate PDB files.
    """
    print("Extracting domains from pdb models ...")
    
    # Create directories for storing the cut PDB files
    output_dirs = [f"{name}-{model}-pfam" for model in ["afs", "esms"]]
    for directory in output_dirs:
        os.makedirs(directory, exist_ok=True)

    # Process each Pfam entry to extract relevant PDB regions
    for uniprot_id, start, end, pfam_id, entry_type, active_sites in pfam_entries:
        for model in ["afs", "esms"]:
            input_pdb = os.path.join(f"{name}-{model}", f"{uniprot_id}.pdb")
            output_pdb = os.path.join(
                f"{name}-{model}-pfam",
                f"{uniprot_id}_{pfam_id}_{start}_{end}.pdb"
            )

            with open(input_pdb, "r") as copy, open(output_pdb, "w") as paste:
                for line in copy:
                    if line.startswith("ATOM"):
                        residue_number = int(line[22:26].strip())
                        if start <= residue_number <= end:
                            paste.write(line)

    print("... Domain extraction complete.")


def run_foldseek(name: str, pfam_entries: List[Tuple[str, int, int, str, str, List[str]]], max_workers: int) -> None:
    """
    Run the Foldseek tool to compute TM-scores between AlphaFold2 and ESMFold models for each Pfam entry.
    """
    # Foldseek seems to crash when multiple instances are executed, probably due to limits in the available memory
    # We force it to execute one job at a time (they still use parallelization)
    # Comment this line if you have enough memory to try and run it with max_workers > 1
    max_workers = 1
    
    print("Running Foldseek ...")
    
    # Directories for Foldseek output and temporary files
    output_dirs = [f"{name}-{suffix}" for suffix in ["foldseek", "fstmp"]]
    for directory in output_dirs:
        os.makedirs(directory, exist_ok=True)

    # Generate Foldseek commands for each Pfam entry
    commands = []
    for uniprot_id, start, end, pfam_id, entry_type, active_sites in pfam_entries:
        afs_pfam_file = f"{name}-afs-pfam/{uniprot_id}_{pfam_id}_{start}_{end}.pdb"
        esms_pfam_file = f"{name}-esms-pfam/{uniprot_id}_{pfam_id}_{start}_{end}.pdb"
        output_file = f"{name}-foldseek/{uniprot_id}_{pfam_id}_{start}_{end}.fs"
        tmp_dir = f"{name}-fstmp"
        
        command = (
            f"foldseek easy-search {afs_pfam_file} {esms_pfam_file} {output_file} {tmp_dir} --alignment-type 1 --prefilter-mode 2 --format-output 'query,target,qtmscore,ttmscore,alntmscore' > /dev/null"
        )
        commands.append(shlex.split(command))

    # Execute commands in parallel
    parallel(max_workers, commands)

    print("... Foldseek processing complete.")


def aggregate_data(name: str, pfam_entries: List[Tuple[str, int, int, str, str, List[str]]]) -> pd.DataFrame:
    """
    Aggregates all available data into a TSV file.
    """
    # Define data structure to hold aggregated information
    column_names = [
        "UniProt ID", "Pfam ID", "Start", "End", "Type", "Active Sites", "# Active Sites", "Class", "Protein Len", "Protein TM-score", "Protein AlphaFold2 pLDDT", "Protein ESMFold pLDDT", "Pfam Len", "Pfam TM-score", "Pfam AlphaFold2 pLDDT", "Pfam ESMFold pLDDT"
    ]
    data = {col: [] for col in column_names}

    # Process each Pfam entry
    for uniprot_id, start, end, pfam_id, entry_type, active_sites in pfam_entries:
        data["UniProt ID"].append(uniprot_id)
        data["Pfam ID"].append(pfam_id)
        data["Start"].append(start)
        data["End"].append(end)
        data["Type"].append(entry_type)
        data["Active Sites"].append(active_sites)
        data["# Active Sites"].append(len(active_sites))

        # Calculate pLDDT scores for proteins and Pfam regions
        for model_suffix, model_name in zip(["afs", "esms"], ["AlphaFold2", "ESMFold"]):
            pdb_file = f"{name}-{model_suffix}/{uniprot_id}.pdb"
            protein_plddts, pfam_plddts = [], []

            with open(pdb_file) as pdb_reader:
                for line in pdb_reader:
                    if line.startswith("ATOM") and line[13:15] == "CA":
                        plddt = float(line[60:66]) / 100
                        protein_plddts.append(plddt)
                        residue_index = int(line[22:26])
                        if start <= residue_index <= end:
                            pfam_plddts.append(plddt)

            data[f"Protein {model_name} pLDDT"].append(np.mean(protein_plddts))
            data[f"Pfam {model_name} pLDDT"].append(np.mean(pfam_plddts))

        data["Protein Len"].append(len(protein_plddts))
        data["Pfam Len"].append(len(pfam_plddts))

        # Extract global TM-score and classify proteins
        json_file = f"{name}-jsons/{uniprot_id}.json"
        with open(json_file) as json_reader:
            json_data = json.load(json_reader)

        data["Protein TM-score"].append(json_data["mm_tmscore"])
        if json_data["has_pdb"]:
            data["Class"].append("Models with PDB")
        elif json_data["mm_tmscore"] >= 0.6:
            data["Class"].append("Similar models")
        else:
            data["Class"].append("Dissimilar models")

        # Extract local TM-score
        fs_file = f"{name}-foldseek/{uniprot_id}_{pfam_id}_{start}_{end}.fs"
        with open(fs_file) as fs_reader:
            local_tm_score = float(fs_reader.readline().split()[-1])
        data["Pfam TM-score"].append(local_tm_score)

    # Convert data dictionary to DataFrame
    df = pd.DataFrame(data)

    # Save DataFrame to TSV file
    tsv_file = f"{name}.tsv"
    df.to_csv(tsv_file, sep="\t", index=False)
    
    print(f"Aggregated data saved to {tsv_file}")

    return df


def read_data(name: str) -> pd.DataFrame:
    '''
    Read aggregated data
    '''
    tsv_file = f"{name}.tsv"
    return pd.read_csv(tsv_file, sep='\t')


# Print and Plot
def print_pfam_statistics(df: pd.DataFrame, writer: TextIO):
    '''
    Print data reported in Table 1 of the manuscript
    '''
    # Define a helper function to compute statistics for a given dataframe
    def compute_stats(group):
        return pd.Series({
            "Unique Pfams": group['Pfam ID'].nunique(),
            "Unique Proteins": group['UniProt ID'].nunique(),
            "Number of Pfams": len(group),
            "Pfam Length Range (Min, Max)": (group['Pfam Len'].min(), group['Pfam Len'].max())
        })

    # Compute statistics for each Pfam type
    stats = df.groupby('Type').apply(compute_stats)

    # Compute statistics for the subset of domains with active sites
    subset = df[(df['Type'] == 'Domain') & (df['# Active Sites'] > 0)]
    subset_stats = compute_stats(subset)
    stats.loc["Domains with Active Sites"] = subset_stats

    # Compute aggregated statistics
    aggregated_stats = compute_stats(df)
    stats.loc["Aggregated"] = aggregated_stats

    # Reindex to enforce a specific row order
    desired_order = ["Domain", "Domains with Active Sites", "Family", "Repeat", "Motif", "Disordered", "Coiled-coil"]
    stats = stats.reindex(desired_order)

    # Write results
    writer.write("Pfam statistics:\n\n")
    writer.write(stats.to_string())
    writer.write("\n\n")


def print_models_quality_statistics(df: pd.DataFrame, writer: TextIO):
    '''
    Print data reported in Table 2-3-4-S2-S3-S4-S5 of the manuscript
    '''
    # Define a function to compute statistics for each group
    def compute_stats(group):
        return pd.Series({
            "Mean Protein TM-score (SD)": f"{group['Protein TM-score'].mean():.2f} ({group['Protein TM-score'].std():.2f})",
            "Mean Pfam TM-score (SD)": f"{group['Pfam TM-score'].mean():.2f} ({group['Pfam TM-score'].std():.2f})",
            "Mean Protein AlphaFold2 pLDDT (SD)": f"{group['Protein AlphaFold2 pLDDT'].mean():.2f} ({group['Protein AlphaFold2 pLDDT'].std():.2f})",
            "Mean Pfam AlphaFold2 pLDDT (SD)": f"{group['Pfam AlphaFold2 pLDDT'].mean():.2f} ({group['Pfam AlphaFold2 pLDDT'].std():.2f})",
            "Mean Protein ESMFold pLDDT (SD)": f"{group['Protein ESMFold pLDDT'].mean():.2f} ({group['Protein ESMFold pLDDT'].std():.2f})",
            "Mean Pfam ESMFold pLDDT (SD)": f"{group['Pfam ESMFold pLDDT'].mean():.2f} ({group['Pfam ESMFold pLDDT'].std():.2f})",
            "Number of Rows": len(group),
            "Number of Unique UniProt IDs": group['UniProt ID'].nunique()
        })

    # Group by Type and Class
    stats = df.groupby(['Type', 'Class']).apply(compute_stats)

    # Compute statistics for the subset of domains with active sites
    subset = df[(df['Type'] == 'Domain') & (df['# Active Sites'] > 0)]
    subset_stats = subset.groupby(['Type', 'Class']).apply(compute_stats)

    # Write the results to the output file
    desired_order_1 = ["Domain", "Domains with Active Sites", "Family", "Repeat", "Motif", "Disordered", "Coiled-coil"]
    desired_order_2 = ["Models with PDB", "Similar models", "Dissimilar models"]
    writer.write("Models quality statistics:\n\n")
    for pfam_type in desired_order_1:#stats.index.get_level_values(0).unique():
        writer.write(f"Statistics for Pfam Type: {pfam_type}\n")
        # Extract the specific table for the current type
        if pfam_type == "Domains with Active Sites":
            type_table = subset_stats.loc["Domain"].transpose()[desired_order_2]
        else:
            type_table = stats.loc[pfam_type].transpose()[desired_order_2]
        writer.write(type_table.to_string())
        writer.write("\n\n")


def print_active_sites_statistics(df: pd.DataFrame, writer: TextIO):
    '''
    Print count of active sites
    '''
    n_as = len(df[(df['Type'] == 'Domain') & (df['# Active Sites'] > 0)])
    n_res = df['Active Sites'].apply(lambda x: len(re.findall(r'\d+', x))).sum()
    writer.write(f"{n_as} active sites are annotated including {n_res} residues (average of {n_res/n_as:.2f} residues per active site)\n")


def plot_models_quality_statistics(df: pd.DataFrame, name: str):
    '''
    Plot data included in the tables with violin plots
    '''
    sns.set_context("notebook", font_scale = 1.5)
    
    # Ensure that "Class" column is ordered properly
    class_order = ["Models\nwith PDB", "Similar\nmodels", "Dissimilar\nmodels"]
    
    # Filter subsets based on the "Type" and "# Active Sites" > 0 condition
    df_domain_active_sites = df[(df["Type"] == "Domain") & (df["# Active Sites"] > 0)]
    df_domain = df[df["Type"] == "Domain"]
    df_family = df[df["Type"] == "Family"]

    # Define the plots and corresponding columns for y-axis
    plots_info = [
        ("TM-score", "Pfam TM-score", "Protein TM-score"),
        ("AlphaFold2 pLDDT", "Pfam AlphaFold2 pLDDT", "Protein AlphaFold2 pLDDT"),
        ("ESMFold pLDDT", "Pfam ESMFold pLDDT", "Protein ESMFold pLDDT")
    ]
    
    # Helper function to plot violin plots
    def plot_violin(df_subset, plot_info, title):
        fig, axes = plt.subplots(1, 3, figsize=(16, 6))

        for i, (plot_title, pfam_col, protein_col) in enumerate(plot_info):
            ax = axes[i]
            # Melt the dataframe for violin plot
            df_melted = df_subset.melt(id_vars=["Class"], value_vars=[pfam_col, protein_col], var_name="Type", value_name=plot_title)
            # Map the names for plotting
            df_melted['Type'] = df_melted['Type'].replace({
                f'Protein {plot_title}': 'Protein',
                f'Pfam {plot_title}': 'Pfam'
            })
            df_melted['Class'] = df_melted['Class'].replace({
                "Models with PDB": "Models\nwith PDB",
                "Similar models": "Similar\nmodels",
                "Dissimilar models": "Dissimilar\nmodels"
            })
            # Plot
            sns.violinplot(data=df_melted, x="Class", y=plot_title, hue="Type", split=True, inner="quart", order=class_order, ax=ax)
            ax.set_title(f'{plot_title}', weight="bold", pad=20)
            ax.set(xlabel=None, ylabel=None, ylim=(-0.1, 1.1))
            ax.legend(loc='lower left')
        
        plt.tight_layout()
        plt.savefig(f"{name}-violins-{title}.png", dpi=500)
        plt.clf()

    # Plot for Domain subset with Active Sites > 0
    plot_violin(df_domain_active_sites, plots_info, 'Domain with active site')

    # Plot for Domain subset
    plot_violin(df_domain, plots_info, 'Domain')

    # Plot for Family subset
    plot_violin(df_family, plots_info, 'Family')


def plot_pfam_lengths(df: pd.DataFrame, name: str):
    '''
    Draw boxplot of Pfam lengths by type
    '''
    sns.set_context("notebook", font_scale = 1.2)
    
    # Duplicate the rows where "Type" == "Domain" and "# Active Sites" > 0
    domain_with_active_site = df[(df['Type'] == 'Domain') & (df['# Active Sites'] > 0)].copy()
    domain_with_active_site['Type'] = 'Domain with\nactive site'
    
    # Append the new rows to the original dataframe
    df_updated = pd.concat([df, domain_with_active_site], ignore_index=True)

    # Define the order of the types to appear on the x-axis
    type_order = ['Domain', 'Domain with\nactive site', 'Family', 'Repeat', 'Motif', 'Disordered', 'Coiled-coil']

    # Create the second figure (Box Plot)
    plt.figure(figsize=(10, 7))
    g = sns.boxplot(data=df_updated, x='Type', y='Pfam Len', order=type_order, fill=False)
    g.set_title("Pfam Lengths by Type", weight="bold", pad=20)
    g.set(xlabel=None, ylabel="Pfam Length", ylim=(-9, 909))
    plt.tight_layout()
    plt.savefig(f"{name}-lens-box.png", dpi=500)
    plt.clf()


def print_and_plot(name: str, df: pd.DataFrame) -> None:
    '''
    Print and plot data for visualization
    '''
    #Print statistics
    output_file = f"{name}-statistics.txt"
    with open(output_file, 'w') as writer:
        print_pfam_statistics(df, writer)
        print_models_quality_statistics(df, writer)
        print_active_sites_statistics(df, writer)
    
    #Make plots
    plot_models_quality_statistics(df, name)
    plot_pfam_lengths(df, name)


# Argument parsing
def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(description="Analyze protein models and domains.")
    parser.add_argument("-n", "--name", type=str, default="enzymes", help="Project name prefix.")
    parser.add_argument("-w", "--max-workers", type=int, default=1, help="Maximum parallel workers.")
    parser.add_argument("--pfamscript", type=str, default="PfamScan/pfam_scan.pl", help="Path to PfamScan script.")
    parser.add_argument("--pfamdb", type=str, default="PfamScan/pfam/", help="Path to PfamScan database directory.")
    parser.add_argument("-d", "--download", action="store_false", help="Set this flag if you already downloaded the files from the webserver.")
    parser.add_argument("-s", "--pfamscan", action="store_false", help="Set this flag if you already generated the PfamScan results.")
    parser.add_argument("-f", "--foldseek", action="store_false", help="Set this flag if you already generated the Foldseek results.")
    parser.add_argument("-a", "--aggregate", action="store_false", help="Set this flag if you already generated the TSV file.")
    return parser.parse_args()


# Main function
def main() -> None:
    """
    Main entry point for the script.
    """
    args = parse_args()
    check_environment(args.name, args.max_workers, args.pfamscript, args.pfamdb, 
                      args.download, args.pfamscan, args.foldseek, args.aggregate)
    
    # Generate or read data
    if args.aggregate:
        uniprot_ids = read_ids(args.name)
        
        if args.download:
            download_files(args.name, uniprot_ids, args.max_workers)
            extract_fastas(args.name, uniprot_ids)
            uniprot_ids = check_downloads(args.name, uniprot_ids)
        
        if args.pfamscan:
            run_pfamscan(args.name, uniprot_ids, args.max_workers, args.pfamscript, args.pfamdb)
        
        pfam_entries = read_pfams(args.name, uniprot_ids)
        
        if args.foldseek:
            cut_domains(args.name, pfam_entries)
            run_foldseek(args.name, pfam_entries, args.max_workers)
        
        df = aggregate_data(args.name, pfam_entries)
    
    else:
        df = read_data(args.name)
    
    print_and_plot(args.name, df)
    print("Pipeline completed.")


if __name__ == "__main__":
    main()
