# Visualizing ADMIXTURE Output with PONG

This pipiline contains helpfull script to generate the files for running
using
[Pong](https://github.com/ramachandran-lab/pong).

## Configuration

The parameters for the pipeline are specified in a `config.yaml` file. Below is
an example configuration:

```yaml
path_to_q_files: "/Users/santiago/treesdemog/experiments/25-01-23-Admixture/results"
nameprefix: "gt"
popinfo: "/Users/santiago//data/genome-asia/2412-GenomeAsiaFreeze3/popinfo.csv"
sample_order: "/Users/santiago/treesdemog/experiments/25-01-23-Admixture/results/sample-ids.txt"
covariate: "Country"  # Can be any other column in popinfo
sample_id_col: "sample_id"
```

- `path_to_q_files`: Path to the directory containing ADMIXTURE Q files.
- `nameprefix`: Prefix for the names of Q files to include in the file map.
- `popinfo`: Path to the CSV file containing population metadata.
- `sample_order`: Path to a file listing sample IDs (output from ADMIXTURE pipeline) that match the order in the Q files.
- `covariate`: Column in the `popinfo` file to use as population labels (e.g., "Country").
- `sample_id_col`: Column in the `popinfo` file containing sample IDs.


## Running the workflow

```bash
# replace config.yaml with your custom data
snakemake -c1 --configfile=config.yaml all
```

## Run pong

```bash
pong -m pong_filemap.txt -i ind2pop.txt -n pop_order.txt
```
