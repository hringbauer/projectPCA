# projectPCA
Project genomes onto pre-computed principal components widely used in ancient DNA. Enables fast analysis without re-computing the principal components. The software accepts ancient DNA data in `eigenstrat` or `PLINK` format as input. No modern samples are required, as the packages include the pre-computed PCA weights and PC coordinates for relevant modern samples (based on publicly available Human Origin array data).

## Installation
The package `projectPCA`is available as a Python package via `pip`. To install, simply run a version of:

```
python3 -m pip install projectPCA
```

## List of available PCAs
As of early 2026, three pre-computed PCAs are officially bundled into `projectPCA`. The bracket denotes the code you can use for all this PCA.

- **HO Westeurasia (HO)**
Standard Western Eurasian PCA, which is widely used in aDNA studies. PC1 corresponds to West-East, and PC2 to North-South.

- **HO Eurasian (EUAS)**
Standard whole-Eurasian PCA, widely used in aDNA studies. Excellent to resolve West versus East Asian ancestry (on PC1). PC2 generally corresponds to North-South.

- **HO Mediterranean (MED)**
This PCA is an extension of the standard Western Eurasian PCA (see above) and includes additional North African populations. That gives a more useful resolution for North African ancestry. The PCA is based on the one from the Puncic paper ([Ringbauer et al 2025](https://doi.org/10.1038/s41586-025-08913-3)). **Important**: You have to use `maf=0.001` here, which is different from the default. See an example Vignette [here](https://www.dropbox.com/scl/fo/n1j71uhyf1yj05bljyrz1/APURb336_K7cYz1NybhW7hQ?rlkey=7uy3tz55l6l28s3li899vcllt&st=fnoj31oj&dl=0).


## Usage

### Project single Samples
To project onto a PCA, the key function is `project_eigenstrat`. To import it and run a single sample, use:

```
from projectPCA.run import project_eigenstrat

project_eigenstrat(es_path="/mnt/archgen/Autorun_eager/eager_outputs/TF/SUA/SUA002/genotyping/pileupcaller.double",
                   pca="HO", es_type="default")
```
This function also returns the dataframe with PCA coordinates. Note that the input path is the path of the eigenstrat files up to `.geno` but without the suffix.

The keyword `pca` denotes which PCA type to project onto (see above). 

If you want to save the figure, you can add the keyword `fig_path=""`. If this string is filled in, the program saves the resulting figure there. 
If the path ends in `.html`, the figure is saved as an interactive plot, where you can hover over the individuals to see their labels (both ancient and modern reference samples). Otherwise, the standard `matplotlib` libraries are used to plot and save the figure (including in `.png` or `.pdf` format, based on the extension you provide).

```
project_eigenstrat(es_path="/mnt/archgen/Autorun_eager/eager_outputs/TF/SUA/SUA002/genotyping/pileupcaller.double",
                   pca="EUAS", es_type="unpacked_fast", plot_bgrd_c=False, fig_path='./figs/SUA002_EUAS.html')
```

### Project multiple samples
It is also possible to project multiple samples. For this, you can use the keyword `iids=[]`. If the keyword is empty (the default), all samples in a file are projected and plotted. If you specify a list of individuals, only individuals with these IDs are projected.


### Project PLINK files
To project PLINK files, you can use the keyword `es_type="plink"`, and provide the path of the PLINK file up to the suffix:

```
project_eigenstrat(es_path="/mnt/archgen/users/hringbauer/git/EPIDEMIC/output/plink/bd_ptn_335",
                   pca="EUAS", es_type="plink", iids=[],
                   plot_bgrd_c=False, verbose=True, flip=True, 
                   fig_path='/mnt/archgen/users/hringbauer/git/projectPCA/figs/ptn335PLINK_EUAS.html')
```


@Harald Ringbauer, 2026