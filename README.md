# ecm-chemoresistance-rnaseq-analysis

Identify extracellular matrix (ECM) components that are associated with chemotherapy drug resistance in endometrial cancer and cervical cancer.

GitHub.io site (work in progress): https://fogg-lab.github.io/ecm-chemoresistance-rnaseq-analysis/

## Setup

### Clone repo

```bash
git clone --recurse-submodules https://github.com/fogg-lab/ecm-chemoresistance-rnaseq-analysis.git
```

### Environment setup

- Recommended (fast): Build the devcontainer or pull the prebuilt Docker container image from docker.io/wigginno/rna (see the environment spec file .devcontainer/devcontainer.json).
- Alternatively, build the container image yourself using the Dockerfile.
- Alternatively, set up a GitHub Codespace from this repo. Select options: repository: fogg-lab/ecm-chemoresistance-rnaseq-analysis, branch: main, dev container config: rna, region: whatever you want, machine type: at least 2-core 8gb

