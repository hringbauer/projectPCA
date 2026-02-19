"""
Classes for loading Genotype Data from an Eigenstrat (or PLINK) file.
Derives from hapROH code (copied Oct 2025)
Important functions: get_geno_i(), get_geno_all()
@ Author: Harald Ringbauer, 2019, All rights reserved
"""
import numpy as np
import pandas as pd


class GenoLoad(object):
    """Class that loads and postprocesses Genotype Files in Python.
    Intended that specific implementations inherit from this dummy class.
    Already contains some key processing functions."""
    base_path = "" # Base path to file object
    nsnp = 0
    nind = 0
    output = True
    df_snp = []  # Dataframe with all SNPs
    df_ind = [] # Dataframe with all iids

    def __init__(self, base_path="", output=True):
        """Concstructor: By default set number of SNPs,IIDs
        and set the SNP and IID dataframe.
        base_path: Data Path without suffix"""
        self.output = output
        self.base_path = base_path

        self.df_snp = self.load_snp_df()   # Load the SNP DataFrame
        self.df_ind = self.load_ind_df()   # Load the Individual DataFrame

        self.nsnp = self.get_n_snps()
        self.nind = self.get_n_inds() 
        
        assert(len(self.df_snp) == self.nsnp)  # Sanity Check
        assert(len(self.df_ind) == self.nind)  # Sanity Check II

        if self.output == True:
            print(f"Loaded Genotype File set with {self.nind} Individuals and {self.nsnp} SNPs.")

    def load_snp_df(self):
        """Load the SNP dataframe.
        Uses self.base_path
        sep: What separator to use when loading the File"""
        raise ImplementationError("Please Implement in your class!")
        
    def load_ind_df(self):
        """Load the Individual dataframe.
        Uses self.base_path
        sep: What separator to use when loading the File"""
        raise ImplementationError("Please Implement in your class!")

    def get_n_snps(self):
        """Return the number of SNPs"""
        raise ImplementationError("Please Implement in your class!")

    def get_n_inds(self):
        """Return the number of Individuals"""
        raise ImplementationError("Please Implement in your class!")
    
    def get_geno_all(self):
        """Load all genotypes from Eigenstrat File.
        Use self.nind for the number of individuals.
        Return genotype matrix, with missing values set to missing_val.
        Expected 2D Format is [#iids,#SNPs]"""
        raise ImplementationError("Please Implement in your class!")

    def get_geno_i(self, i):
        """Load Individual i.
        Return 1D genotype array"""
        raise ImplementationError("Please Implement in your class!")

    def get_geno_iid(self, iid):
        """Return Genotypes of Individual iid"""
        i = self.get_index_iid(iid)
        g = self.get_geno_i(i)
        return g

    def give_positions(self, ch):
        """Return Array of Positions and Indices of
         all SNPs on Chromosome ch"""
        df_snp = self.df_snp
        ch_loci = (df_snp["chr"] == ch)
        idcs = np.where(ch_loci)[0]
        assert(len(idcs) > 0)
        pos = df_snp.loc[ch_loci, "pos"]
        return pos, idcs

    def give_ref_alt(self, ch):
        """Return Arrays of Ref/Alt of all SNPs on Chromosome ch"""
        df_snp = self.df_snp
        df_t = df_snp.loc[df_snp["chr"] == ch, ["ref", "alt"]]  # Subset to ch
        ref, alt = df_t["ref"].values, df_t["alt"].values
        return ref, alt

    def get_index_iid(self, iid):
        """Get Index of Individual iid"""
        # Detect the Individual
        found = np.where(self.df_ind["iid"] == iid)[0]
        if len(found)==0:
            raise RuntimeError(f"Individual {iid} not found in Eigenstrat!")
        else: 
            i = found[0]
        return i


class EigenstratLoad(GenoLoad):
    """Class that loads and postprocesses Eigenstrats"""

    def __init__(self, base_path="", output=True, sep=r"\s+"):
        """Concstructor:
        base_path: Data Path without .geno/.snp/.ind).
        ch: Which chromosome to load
        sep: What separator to use when loading the File"""
        self.output = output
        if len(base_path) > 0:
            self.base_path = base_path
        geno_file = open(self.base_path + ".geno", "rb")
        header = geno_file.read(21)  # Ignoring hashes for geno/tgeno set to 20/21
        self.nind, self.nsnp = [int(x) for x in header.split()[1:3]]
        # assuming sizeof(char)=1 here
        self.rlen = max(48, int(np.ceil(self.nind * 2 / 8)))

        self.df_snp = self.load_snp_df(sep=sep)   # Load the SNP DataFrame
        self.df_ind = self.load_ind_df(sep=sep)   # Load the Individual DataFrame
        assert(len(self.df_snp) == self.nsnp)  # Sanity Check
        assert(len(self.df_ind) == self.nind)  # Sanity Check II

        if self.output == True:
            print(f"Loaded Eigenstrat File set with {self.nind} Individuals and {self.nsnp} SNPs.")

    def load_snp_df(self, sep=r"\s+"):
        """Load the SNP dataframe.
        Uses self.base_path
        sep: What separator to use when loading the File"""
        if len(self.df_snp)==0:          
            path_snp = self.base_path + ".snp"
            df_snp = pd.read_csv(path_snp, header=None,
                                 sep=sep, engine="python")
            df_snp.columns = ["snp", "chr", "map",
                              "pos", "ref", "alt"]  # Set the Columns

        else:
            df_snp = self.df_snp
        return df_snp

    def load_ind_df(self, sep=r"\s+"):
        """Load the Individual dataframe.
        Uses self.base_path
        sep: What separator to use when loading the File"""
        if len(self.df_ind)==0:
            path_ind = self.base_path + ".ind"
            df_ind = pd.read_csv(path_ind, header=None,
                                 sep=r"\s+", engine="python")
            df_ind.columns = ["iid", "sex", "cls"]  # Set the Columns
            df_ind = df_ind.astype("str") # Make sure everything is string
        else:
            df_ind = self.df_ind
        return df_ind
    
    def get_geno_all(self, missing_val=3):
        """Load all genotypes from Eigenstrat File.
        Use self.nind for number of individuals.
        Return genotype matrix, with missing values set to missing_val"""
        geno = self.give_bit_file()  # Load the whole bit file
        gt = np.unpackbits(geno, axis=1)[:,:2*self.nind]
        gt = 2 * gt[:, 0::2] + gt[:, 1::2]
        gt = update_values(gt, x=[0,1,2,3], y=[2,1,0,np.nan], copy=True) # use COPY as values overlap
        #gt[gt == 3] = missing_val  # set missing values
        return gt

    def get_geno_i(self, i, missing_val=np.nan):
        """Load Individual i"""
        batch, eff = self.get_enc_index(i)
        geno = self.give_bit_file()  # Load the whole bit file

        geno_sub = geno[:, [batch]]  # Byte value of batch
        geno_sub = np.unpackbits(geno_sub, axis=1)[:, 2 * eff:2 * eff + 2]

        ### Update to PCA encoding
        geno_sub = 2 * geno_sub[:, 0] + geno_sub[:, 1]
        geno_sub = update_values(geno_sub, x=[0,1,2,3], y=[2,1,0,np.nan], copy=True) # use COPY as values overlap
        #geno_sub[geno_sub == 3] = missing_val  # set missing values
        return geno_sub

    def give_bit_file(self):
        base_path = self.base_path
        geno = np.fromfile(self.base_path + ".geno",
                           dtype='uint8')[self.rlen:]  # without header
        geno.shape = (self.nsnp, self.rlen)
        return geno

    def get_enc_index(self, i):
        """Get the Index in the Encoding and the modulo 4 value
        (position in batch)"""
        rlen_sub = int(np.floor(i * 2 / 8))  # Effectively dividing by four
        mod_i = i % 4  # Calculate the rest
        return rlen_sub, mod_i

    def extract_snps(self, id, markers, conversion=True, dtype=np.int8):
        """Extract SNPs for Integer Index i on marker list
        markers. If conversion: Convert to VCF type encoding
        Load all SNPs and then subset"""
        geno = self.get_geno_i(id)
        geno = geno[markers] # Subset to markers

        if conversion:
            geno_new = -np.ones((2,len(geno)), dtype=dtype)
            geno_new[:, geno==0]=1    # 2 Derived Alleles
            geno_new[:, geno==2]=0    # 2 Ancestral Alleles
            ### Heterozgyotes
            geno_new[0, geno==1]=1
            geno_new[1, geno==1]=0
            geno = geno_new

        return geno

#########################################################
#########################################################
#### Subclass for Non-Binary Eigenstrats

class EigenstratLoadUnpacked(EigenstratLoad):
    """Class that loads and postprocesses Eigenstrats.
    Same as Superclass, but overwrites methods to load
    non-binary encoded Genotype Data"""

    def __init__(self, base_path="", output=True, sep=r"\s+"):
        """Overwrite Concstructor:
        base_path: Data path without .geno/.snp/.ind.
        ch: Which chromosome to load
        sep: What separator to use when loading the File"""
        self.output = output
        if len(base_path) > 0:
            self.base_path = base_path
        ### Get Size of Data Matrix and sanity check
        with open(self.base_path + ".geno",'r') as f:
            t = f.read()
            l = t.splitlines()
            self.nind = len(l[0])
            self.nsnp = len(l)

        self.df_snp = self.load_snp_df(sep=sep)   # Load the SNP DataFrame
        self.df_ind = self.load_ind_df(sep=sep)   # Load the Individual DataFrame
        assert(len(self.df_snp) == self.nsnp)  # Sanity Check
        assert(len(self.df_ind) == self.nind)  # Sanity Check II

        if self.output:
            print(f"3 Eigenstrat Files with {self.nind} Individuals and {self.nsnp} SNPs")

    def get_geno_i(self, i, missing_val=np.nan):
        """Load Genotype for Individual (Row) i,
        assuming it's encoded in unpacked Format."""
        geno=np.genfromtxt(self.base_path + ".geno", delimiter=1, usecols=i, dtype="float")
        idx = geno == 9
        geno = 2 - geno
        geno[idx] = missing_val
        return geno

    def get_geno_all(self, missing_val=np.nan):
        """Load all genotypes from Eigenstrat File.
        Use self.nind for number of individuals.
        Return genotype matrix, with missing values set to missing_val"""
        geno=np.genfromtxt(self.base_path + ".geno", delimiter=1, dtype="float", ndmin=2)
        idx = geno == 9
        geno = 2 - geno
        geno[idx] = missing_val
        return geno.T


class EigenstratEager(EigenstratLoad):
    """Load Eigenstrat from Autoeager output. Faster than the default version.
    Knows that the eager eigenstrats are unpacked and contains only one IID"""

    def __init__(self, base_path="", output=True, sep=r"\s+", nsnp=0, nind=0):
        """Overwrite general Constructor:
        base_path: Data path without .geno/.snp/.ind.
        ch: Which chromosome to load
        sep: What separator to use when loading the File.
        nsnp: Can pre-specify SNP number (to avoid pre-loading SNP file)
        nind: Can pre-specify IID number (to avoid pre-loaind find file"""
        self.output = output
        if len(base_path) > 0:
            self.base_path = base_path

 
        self.df_ind = self.load_ind_df(sep=sep)   # Load the Individual DataFrame
        nind = len(self.df_ind)
            
        if nsnp==0:
            self.df_snp = self.load_snp_df(sep=sep)   # Load the SNP DataFrame
            nsnp = len(self.df_snp)
            
        self.nind = nind
        self.nsnp = nsnp
        
        if self.output:
            print(f"3 Eigenstrat Files with {self.nind} Individuals and {self.nsnp} SNPs")

    def get_geno_i(self, i=0):
        """Load Individual i"""
        gt = self.get_geno_all()  # Load the whole bit file
        g = gt[i,:]
        return g

    def get_geno_all(self, missing_val=3):
        """Load all genotypes from the Eigenstrat File.
        Use self.nind for the number of individuals.
        Return genotype matrix, with missing values set to missing_val"""
        gt = self.load_eager_geno()  # Load the whole bit file
        gt = update_values(gt, x=[48,49,50,57], y=[2,1,0,np.nan], copy=True) # use COPY as values overlap
        return gt

    def load_eager_geno(self):
        """Load genotype matrix from autoeager .geno.
        Return [#INDs, #SNPs] 2D Array"""
        # Read the unpacked genotype file

        path_gt = self.base_path + ".geno"
        with open(path_gt, "rb") as f:
            # Load all genotype data as unsigned 8-bit integers
            genotype_data = np.frombuffer(f.read(), dtype=np.uint8)

        num_cols = self.nind + 1
        num_snps = self.nsnp
            
        # Check that the file size matches the expected number of genotypes
        expected_size = num_cols * num_snps # Add one column for the last one
        assert genotype_data.size == expected_size, "File size does not match expected number of genotypes based on samples and SNPs"
    
        # Reshape the data into a matrix with shape (num_samples, num_snps)
        genotype_matrix = genotype_data.reshape((num_snps, num_cols)).T
        return genotype_matrix[:-1,:] # Delete last, empty row

    def get_index_iid(self, iid):
        """Get Index of Individual iid.
        Here: Return always zero because eager output"""
        # Detect the Individual
        return 0

class EigenstratLoadUnpackedFast(EigenstratLoadUnpacked):
    """Fast Class to load unpacked Eigenstrat File.
    Optimized version to load single iids (fast)"""

    def get_geno_i(self, i, missing_val=np.nan):
        """Load Genotype for Individual (Row) i,
        assuming it's encoded in unpacked Format."""
        raw = np.fromfile(self.base_path + ".geno", dtype=np.uint8)
        raw = raw[raw != 10]   # remove '\n'
        
        geno = raw.reshape(self.nsnp, self.nind) # Put into Matrix
        geno = geno[:, i].astype('float') ### Extract Relevant column
        geno = update_values(geno, x=[48,49,50,57], y=[2,1,0,missing_val], copy=False) 
        return geno

    def get_geno_all(self, missing_val=np.nan):
        """Return all genotypes."""
        raw = np.fromfile(self.base_path + ".geno", dtype=np.uint8)
        raw = raw[raw != 10]   # remove '\n'

        geno = raw.reshape(self.nsnp, self.nind).T.astype('float') # Put into Matrix
        geno = update_values(geno, x=[48,49,50,57], y=[2,1,0,missing_val], copy=False) 
        return geno

#########################################################
### Load PLINK file

class PlinkLoad(GenoLoad):
    """Class that loads and postprocesses PLINK files.
    Same as Eigenstrat Superclass, but overwrites methods to load
    SNP and IID info as well as encoded Genotype Data"""
    bed_path = "" # Path to full bed file
    bed = None # Place for bed reader object

    def __init__(self, base_path="", output=True):
        """Overwrite Concstructor:
        base_path: Data path without the .bed ending.
        output: Whether to print output."""
        ### Lazy Import and load bed object
        from bed_reader import open_bed # Import PLINK library only if needed.
        self.bed_path = base_path + ".bed"
        self.bed = open_bed(self.bed_path)

        #### Run Standard Init
        super(PlinkLoad, self).__init__(base_path=base_path, output=output)

    def get_n_snps(self):
        """Return number of SNPs"""
        return self.bed.sid_count

    def get_n_inds(self):
        """Return number of IIDs"""
        return self.bed.iid_count
    
    def load_snp_df(self):
        """Load the SNP dataframe.
        Uses bed object"""
        if len(self.df_snp)==0:          
            # Create SNP metadata DataFrame
            df_snp = pd.DataFrame({
                "snp": self.bed.sid,                 # SNP ID
                "chr": self.bed.chromosome,          # Chromosome
                "map": self.bed.cm_position/100,         # Genetic map position (cM)
                "pos": self.bed.bp_position.astype("int"),         # Base pair position
                "ref": self.bed.allele_2,            # A1 is Plink Alt. A2 is Plink Ref.
                "alt": self.bed.allele_1             
            })

        else:
            df_snp = self.df_snp
        return df_snp

    def load_ind_df(self):
        """Load the Individual dataframe.
        Uses bed object"""
        if len(self.df_ind)==0:
            df_ind = pd.DataFrame({"iid":self.bed.iid.astype("str"), 
                                   "sex":self.bed.sex.astype("str"), 
                                   "cls":""})
        else:
            df_ind = self.df_ind
        return df_ind

    def get_geno_all(self):
        """Load all genotypes from Eigenstrat File.
        Use self.nind for number of individuals.
        Return genotype matrix, with missing values set to missing_val"""
        gt = self.bed.read()
        gt = 2 - gt # To adjust for plink encoding
        return gt
    
    def get_geno_i(self, i):
        """Load Individual i"""
        geno_sub = self.bed.read(np.s_[i,:])
        geno_sub = 2 - geno_sub
        return geno_sub[0]



#########################################################
#########################################################
### Helper Function

def is_binary_file(path, extension=".geno"):
    """Test whether a file at path + extension is binary.
    Return boolean if the case """
    binary=False
    try:
        with open(path + extension, "r") as f:
            t = f.readline()
    except UnicodeDecodeError:
        binary=True
    return binary

def update_values(gt, x=[48,49,50,57], y=[2,1,0,9], copy=False):
    """"Update Values in numpy matrix gt. 
    x: List of original values.
    y: Updated values"""
    if copy:
        gt2 = gt.astype("float") # hard copy, necessary if values intersect

    else:
        gt2 = gt # only pointer
        
    for i,j in zip(x,y):
    #print(x,y)
        idx = gt == i
        gt2[idx] = j
    return gt2


#########################################################
#########################################################
#### Factory Method

def get_eigenstrat_object(base_path, mode="default", sep=r"\s+", packed=-1, verbose=True):
    """Factory Method to Load Eigenstrat object
    sep: What separator to use when loading an Eigenstrat File. 
    The default is space-separated (by an arbitrary number of spaces).
    mode: Which mode to use - can be: default / eager / autoeager / plink
    Packed: Whether Genotype Data is encoded in binary Format"""

    if mode=="default":
        ### Determine automatically
        if packed==-1:
            packed = is_binary_file(base_path, extension=".geno")
            if verbose:
                print(f"Eigenstrat packed: {packed}")

        ### Load
        if packed:
        	es = EigenstratLoad(base_path, output=verbose, sep=sep)
        else:
        	es = EigenstratLoadUnpacked(base_path, output=verbose, sep=sep) 

    elif mode=="eager":
            es = EigenstratEager(base_path, output=verbose, sep=sep)

    elif mode=="autoeager":
            es = EigenstratEager(base_path, output=verbose, sep=sep, nind=1, nsnp=1233013)
        
    elif mode=="unpacked_fast":
            es = EigenstratLoadUnpackedFast(base_path, output=verbose, sep=sep)

    elif mode=="plink":
            es =  PlinkLoad(base_path, output=verbose)
        
    else:
        raise RuntimeError(f"Mode: {mode} not found. \nPlease use one of: standard/eager/autoeager/unpacked_fast/plink")
        
    return es