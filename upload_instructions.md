# Instructions on how to upload a package to PYPI
This file is based on update and upload instructions from `https://packaging.python.org/tutorials/packaging-projects/`
Harald keeps this file as notes and a reminder.
The installation is based on a pyproject.toml setup. It is a pure Python project without extensions. I distribute both source and any-wheel.

## 1) Load Python environment and navigate to folder
- On Leipzig HPC Cluster (primary source to create package since 2025):
Activate Python environment and go to the hapROH package folder:

pyenvhpc312
cd /mnt/archgen/users/hringbauer/git/projectPCA/

### 1b) For a local test install (e.g., for unit tests), run from ./package/ folder:
pip3 install ./

Add flag  `--user` if not in a virtual environment.

### 1c) Run Tests of Expected Behavior
- [To implement]


## 2) Create the Source Package 
Update version in `./pyproject.toml` to next version number and update `./change_log.md`

### Clean up prior packages
Delete previous `./dist/*`:

rm ./dist/*

### Build package 
python3 -m build

### Upload to PyPi
### For full PyPi server
python3 -m twine upload dist/* 

### [Alternatively] Upload to the test server (for testing)
python3 -m twine upload --repository-url https://test.pypi.org/legacy/ dist/* 


### Test install of uploaded package

Re-install to check whether the new hapROH version installs properly:

python3 -m pip install --upgrade --no-deps --force-reinstall projectPCA

Add the `--user` flag to the Python settings where packages need to go to the user. 


# Resources for further reading
### for packaging: 
https://packaging.python.org/tutorials/packaging-projects/

### for version numbers:
https://www.python.org/dev/peps/pep-0440/
