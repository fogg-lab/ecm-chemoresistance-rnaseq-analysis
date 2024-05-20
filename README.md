# ecm-chemoresistance-rnaseq-analysis

Identify extracellular matrix (ECM) components that are associated with chemotherapy drug resistance in endometrial cancer and cervical cancer.

Link to github.io page (currently contains figures generated from analysis on TCGA-CESC): https://fogg-lab.github.io/ecm-chemoresistance-rnaseq-analysis/

## Instructions

### 1. Clone repo

**Important**: Make sure [git lfs](https://git-lfs.com/) is installed before running:

```bash
git clone --recurse-submodules https://github.com/fogg-lab/ecm-chemoresistance-rnaseq-analysis.git
```

### 2. Set up environment

- Recommended (fast): Build the devcontainer or pull the prebuilt Docker container image from docker.io/wigginno/rna (see the environment spec file .devcontainer/devcontainer.json).
- Alternatively, build the container image yourself using the Dockerfile.
- Alternatively, set up a GitHub Codespace from this repo. Select options: repository: fogg-lab/ecm-chemoresistance-rnaseq-analysis, branch: main, dev container config: rna, region: whatever you want, machine type: 4-core 16gb

### 3. Reproduce results

Run the notebooks in sequential order.
