# Instructions on how to upload a package to PYPI
This instruction is based on `https://packaging.python.org/tutorials/packaging-projects/`
Harald keeps this file as notes and a reminder.
The installation is based on a `pyproject.toml` setup. `projectPCA` is a pure-Python project with no extensions. I distribute both source and any-wheel.

## 1) Load Python environment and navigate to folder
- On Leipzig HPC Cluster (primary source to create package since 2025):
Activate Python environment and go to the hapROH package folder:

pyenvhpc312
cd /mnt/archgen/users/hringbauer/git/projectPCA/

### 1b) For a local test install (e.g., for unit tests), run from ./package/ folder:
pip3 install ./

Add flag  `--user` if not in a virtual environment.

### 1c) Run Tests of Expected Behavior
- Run Notebook for tests: `./tests/test_proj_HO_Leipzig.ipynb`

## 2) Create the Source Package 
Update version in `./pyproject.toml` to next version number and update `./change_log.md`

### 2b) Clean up prior packages
Delete previous `./dist/*`:

rm ./dist/*

### 2c) Build package 
python3 -m build

### 2d) Upload to PyPi
### For full PyPi server
python3 -m twine upload dist/* 

### [Alternatively Testing] Upload to the test server
python3 -m twine upload --repository-url https://test.pypi.org/legacy/ dist/* 


### 3) Test install of uploaded package

Re-install to check whether the new `projectPCA` version installs:

python3 -m pip install --upgrade --no-deps --force-reinstall projectPCA

Add `--user` flag to the Python settings where packages need to go to the user. 


# Resources for further reading
### for packaging: 
`https://packaging.python.org/tutorials/packaging-projects/`

### for version numbers:
`https://www.python.org/dev/peps/pep-0440/`
