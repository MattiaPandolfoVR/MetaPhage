#!/usr/bin/env python3
"""
Generate a Metaphage configuration file for a new project.

CHECKS
- put metadata in the right place where it is the only csv file
- cpus/memory are defined?
"""
import os, sys
import argparse
import subprocess
from string import Template
import tempfile
import csv
import multiprocessing
__version__ = "1.0.0"
meminfo = {}
meminfo["MemTotal"] = 8
try:
    if os.path.exists("/proc/meminfo"):
        meminfo = dict((i.split()[0].rstrip(':'),int(i.split()[1])) for i in open('/proc/meminfo').readlines())
except Exception:
    print("WARNING: Could not read memory info from /proc/meminfo. Assuming 8GB of RAM.", file=sys.stderr)
    pass

MEM = round(meminfo['MemTotal'] / 1000000)
CORES = multiprocessing.cpu_count()

def has_singularity():
    cmd = ["singularity", "--version"]
    try:
        subprocess.check_call(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return True
    except Exception as e:
        print(f"INFO: Singularity not found.", file=sys.stderr)
        return False

def has_docker_user():
    cmd = ["docker", "run", "hello-world"]
    try:
        subprocess.check_call(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return True
    except Exception as e:
        print(f"INFO: Docker not found.", file=sys.stderr)
        return False

def has_nextflow():
    cmd = ["nextflow", "-version"]
    try:
        subprocess.check_call(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return True
    except Exception as e:
        print(f"INFO: Nextflow not found.", file=sys.stderr)
        return False

class Local:
    raw = Template("""
    params {
        max_cpus = $CORES
        max_memory = $MEM.GB
        max_time = 72.h
    }
    """)

    def __init__(self, cpus, mem, env):
        self.cpus = int(cpus) if cpus > 0 else 1
        self.mem = int(mem) if mem  > 0 else 1
        if env == "singularity":
            self.env = "singularity"
        self.template = (self.raw).substitute(CORES=self.cpus, MEM=self.mem)

    def __str__(self):
        return self.raw.substitute(CORES=self.cpus, MEM=self.mem)


class Metaphage:
    
    raw = Template("""
    params {
        config_profile_name = '$projectName'
        config_profile_description = 'MetaPhage analysis configuration'

        // INPUT PATHS
        readPath = "$readPath"
        fqpattern = "$fqpattern"
        metaPath = "$metaPath"
        dbPath = "$dbPath"

        // OUTPUT/WORKING PATHS
        outdir = "$outdir"
        temp_dir = "$temp_dir"

        // METADATA 
        metadata = $metadata
        virome_dataset = true
        singleEnd = $singleEnd
        sum_viol_var = "$sum_viol_var" 
        heatmap_var = "$heatmap_var" 
        alpha_var1 = "$alpha_var1" 
        alpha_var2 = "$alpha_var2" 
        beta_var = "$beta_var" 
        violin_var = "$violin_var" 


    }
    """)

  
    def __init__(self, args) -> None:
        
        self.bins = os.path.dirname(os.path.realpath(__file__))
        # Self dir is the parent dir of bin
        self.dir        = os.path.dirname(self.bins)
        self.input      = os.path.realpath(args.readsdir)

        # Check input data
        if (not os.path.exists(self.input)):
            print("ERROR: Input directory not found: {}".format(self.input), file=sys.stderr)
            exit(0)

        if (args.metadata is not None and not os.path.exists(args.metadata)):
            print("ERROR: Metadata file not found: {}".format(args.metadata), file=sys.stderr)
            exit(0)


        # Metadata
        metadatadir     = os.path.join(self.input, "metadata")
        if (not os.path.exists(metadatadir)):
            if args.verbose:
                print("INFO: Making metadata directory: {}".format(metadatadir), file=sys.stderr)
            try:
                os.mkdir(metadatadir)
            except Exception as e:
                print("ERROR: Unable to make metadata directory: {}:\n  {}".format(metadatadir, e), file=sys.stderr)
        else:
            metadatadir = tempfile.mkdtemp(prefix="MetaPhage_", suffix="_metadata", dir=self.input)
        self.pattern    = ""
        
        # Copy metadata to metadatadir
        if args.metadata is not None:
            metadataSrc  = os.path.realpath(args.metadata)
            metadataDest = os.path.join(metadatadir, 'metadata.csv')
            os.symlink(metadataSrc, metadataDest)
            self.metadata = metadataDest
        else:
            self.metadata = None
        # categories 
        self.var        = args.main_variable
        self.output     = os.path.realpath(args.output_dir)
        self.project    = args.project if args.project else "MetaPhage project"
        self.tempdir    = os.path.realpath(args.tempdir)
        self.workingdir = os.path.realpath(args.workingdir)
        self.db         = os.path.realpath(args.database_dir) if args.database_dir else os.path.join(self.dir, 'db')
        self.verbose    = args.verbose
        self.valid      = True
        self.ready      = True
        self.has_database = True if self.db else False
        self.has_metadata = True if self.metadata else False
        self.has_access   = False
        self.samples      = []
        self.vars         = []
        self.vMain        = args.main_variable if args.main_variable else "Treatment"
        self.vAlpha1      = args.alpha_div_1 if args.alpha_div_1 else args.main_variable
        self.vAlpha2      = args.alpha_div_2 if args.alpha_div_2 else args.main_variable
        self.vBeta        = args.beta_div if args.beta_div else args.main_variable
        self.vViolin      = args.violin if args.violin  else args.main_variable
        self.vHeatmap     = args.heatmap if args.heatmap  else args.main_variable
        self.vSumViol     = args.sum_viol_var if args.heatmap   else args.main_variable

        self.singleEnd  = False
     
        # Check database and installation
        self.checkInstallation()
        # Check input dir
        self.setSamples()
        self.setMetadata()

    def cmd(self):
        return [
            "nextflow",
            "run",
            os.path.join(self.dir, "main.nf"),
            '--projectName', self.project,
            '--readPath', self.input,
            '--fqpattern', self.pattern,
            '--metaPath', os.path.dirname(self.metadata),
            '--dbPath', self.db,
            '--metadata', 'true' if self.metadata else 'false',
            '--singleEnd', 'true' if self.singleEnd else 'false',
            '--sum_viol_var', self.vSumViol,
            '--heatmap_var', self.vHeatmap,
            '--alpha_var1', self.vAlpha1,
            '--alpha_var1', self.vAlpha2,
            '--beta_var', self.vBeta,
            '--violin_var', self.vViolin,
            '--outdir', self.output,
            '--temp_dir', self.tempdir
        ]
    def template(self):
        try:
                
            dictionary = {
                'projectName': self.project,
                'readPath': self.input,
                'fqpattern': self.pattern,
                'metaPath': os.path.dirname(self.metadata),
                'dbPath': self.db,
                'metadata': 'true' if self.metadata else 'false',
                'singleEnd': 'true' if self.singleEnd else 'false',
                'sum_viol_var': self.vSumViol,
                'heatmap_var': self.vHeatmap,
                'alpha_var1': self.vAlpha1,
                'alpha_var2': self.vAlpha2,
                'beta_var': self.vBeta,
                'violin_var': self.vViolin,
                'outdir': self.output,
                'temp_dir': self.tempdir
            }
            string = self.raw.substitute(dictionary)
            return string
        except Exception as e:
            print(e, file=sys.stderr)
            return None

    def __repr__(self) -> str:
        print("\nMetaphage configuration:", file=sys.stderr)
        print("\tDatabase:   {}".format(self.db), file=sys.stderr)
        print("\tInput:      {}".format(self.input), file=sys.stderr)
        print("\tMetadata:   {}".format(self.metadata), file=sys.stderr)
        print("\tOutput:     {}".format(self.output), file=sys.stderr)
        print("\tProject:    {}".format(self.project), file=sys.stderr)
        print("\tTempdir:    {}".format(self.tempdir), file=sys.stderr)
        print("\tWorkingdir: {}".format(self.workingdir), file=sys.stderr)
        
    def __str__(self) -> str:
        max = len(self.samples) if len(self.samples) < 4 else 4
        sample_type = "Single End" if self.singleEnd else "Paired End"
        
        return f"""
        Database:   {self.db}
        Input:      {self.input}
        Metadata:   {self.metadata}
        Output:     {self.output}
        Project:    {self.project}
        Tempdir:    {self.tempdir}
        Workingdir: {self.workingdir}
        
        Samples:    {len(self.samples)} {sample_type} ({','.join(self.samples[0:max])}...)
        Metadata:   {len(self.vars)} variables""".rstrip()
        

    def checkInstallation(self):
        # Check directory
        for filename in ['nextflow.config', 'main.nf']:
            if not os.path.exists(os.path.join(self.dir, filename)):
                print('ERROR: {} not found in {}: MetaPhage installation looks corrupted.'.format(filename, self.dir), file=sys.stderr)
                self.valid = False

        # Check database
        if self.has_database:
            if os.path.isdir(self.db):
                self.checkDatabase()
            else:
                print('ERROR: Database directory not found: {}'.format(self.db), file=sys.stderr)
                self.valid = False
            

        # Check environment is active (non fatal)
        command = ["phigaro", "-V"]
        try:
            # Check if command is available without printing anything
            subprocess.check_output(command, stderr=subprocess.STDOUT)
            self.has_access = True
        except FileNotFoundError or subprocess.CalledProcessError:
            self.has_access = False
            print("INFO: this environment is not ready to run MetaPhage. Remember to use a container or activate the environment.", file=sys.stderr)
        except Exception:
            self.has_access = False
            print("INFO: this environment is not suitable for metaphage. Use a container or activate the environment.", file=sys.stderr)
    def setSamples(self):
        extensions = ['.fastq', '.fq', '.fastq.gz', '.fq.gz']
        tags = ['_1', '_2', '_R1_001', '_R2_001', '_R1', '_R2']
        patterns = {'_1' : '_{1,2}', '_2' : '_{1,2}', '_R1_001' : '_R{1,2}_001', '_R2_001' : '_R{1,2}_001', '_R1' : '_R{1,2}', '_R2' : '_R{1,2}'}
        samples = {}
        thispattern = []
        for file in os.listdir(self.input): 
            original = file
            # Check if one of the extensions is at the end of the file
            for ext in extensions:
                    if file.endswith(ext):
                        # Remove ext
                        file = file[:-len(ext)]
                        # Check if one of the tags is in the file
                        for tag in tags:
                            if tag in file:
                            
                                # Substitute globalextension tag with X
                                pattern = original[-len(tag)-len(ext):].replace(tag, patterns[tag]) 
                                if pattern not in thispattern:
                                    thispattern.append(pattern)
                                    
                                file = file[:-len(tag)]
                                if file in samples:
                                    samples[file] += 1
                                else:
                                    samples[file] = 1
                                break
                        break
        if len(thispattern) != 1:
            print(f"Error: more than one pattern found for the samples {thispattern}")
            self.valid = False
        else:
            self.pattern = thispattern[0]
        
        # Check the value of samples
        if len(samples) == 0:
            print("ERROR: No fastq files found in {}".format(self.input), file=sys.stderr)
            print("  Valid extensions are .fq, .fastq (.gz). Pair tags are _1/_2 or ) _R1/_R2", file=sys.stderr)
            self.valid = False
        
        paired = 0
        single = 0
        for key, val in samples.items():
            if val == 2:
                paired += 1
            if val == 1:
                single += 1
            self.samples.append(key)
        if paired > 0 and single > 0:
            print("ERROR: Mixed single and paired-end samples detected in {}".format(self.input), file=sys.stderr)
            self.valid = False
        elif paired > 0:
            self.singleEnd = False
        elif single > 0:
            self.singleEnd = True
        
        
        print(f"INFO: Found {len(self.samples)} samples in {self.input}", file=sys.stderr)
        
    def setMetadata(self):
        if self.has_metadata:
            # Load csv without pandas (use csv standard module)
            c = 1
            with open(self.metadata, 'r') as f:
                spamreader = csv.reader(f, delimiter=',', quotechar='"')
                # Check if the first line is the header
                header = next(spamreader)
                self.vars = header[1:]
                if header[0] != 'Sample':
                    print(f"ERROR: Metadata first column is not 'Sample' ({header[0]})", file=sys.stderr)
                    self.valid = False
                    return
                
                for field in [self.vMain, self.vAlpha1, self.vAlpha2, self.vBeta, self.vViolin, self.vHeatmap, self.vSumViol]:
                    if field not in header:
                        c = "\n - "
                        print(f"ERROR: Metadata column '{field}' not found in header:\npossible values are {c.join(header)}", file=sys.stderr)
                        self.valid = False
                        return
        
                for row in spamreader:
                    c += 1
                    if len(row) != len(header):
                        print(f"ERROR: Metadata error at row {c}:\nRow '{row}' has a different length than header ({len(header)} != {len(row)})", file=sys.stderr)
                        self.valid = False
                        return
                    if row[0] not in self.samples:
                        print(f"ERROR: Sample {row[0]} in your metadata file is not present in the reads directory.", file=sys.stderr)
        
    def checkDatabase(self):
        for subdir in [  'inphared',  'kraken2',  'phigaro',  'phix',  'vibrant',  'virsorter']:
            if not os.path.isdir(os.path.join(self.db, subdir)):
                print('ERROR: Database directory not found: {}'.format(os.path.join(self.db, subdir)), file=sys.stderr)
                self.valid = False
        
        for subdir in [ 'diamond' ]:
            if not os.path.isdir(os.path.join(self.db, subdir)):
                print('WARNING: Database directory will be created: {}'.format(os.path.join(self.db, subdir)), file=sys.stderr)
                self.ready = False            


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate a Metaphage configuration file for a new project.")
    # Add group "Main"
    main = parser.add_argument_group("Main arguments")
    main.add_argument("-s", "--save", help="Configuration file output [default: %(default)s]")
    
    main.add_argument("-i", "--reads-dir", dest="readsdir", required=True, help="Directory containing the reads.")
    main.add_argument("-o", "--output-dir", help="Output directory [default: %(default)s]", default="./MetaPhage")
    main.add_argument("-d", "--database-dir", help="Database directory", required=False)
    main.add_argument("-m", "--metadata-file", dest="metadata", required=False, help="Metadata file.")
    main.add_argument("-p", "--project", help="Project name", required=False)
    
    # Add group Metadata
    metadata = parser.add_argument_group("Metadata arguments")
    metadata.add_argument("-v", "--main-variable", help="Default variable of the metadata table for comparisons", required=False)
    metadata.add_argument("-a", "--alpha-div-1", help="Variable for alpha diversity (otherwise -v)", required=False)
    metadata.add_argument("-A", "--alpha-div-2", help="Secondary variable for alpha diversity (otherwise -v)", required=False)
    metadata.add_argument("-b", "--beta-div", help="Variable for alpha diversity (otherwise -v)", required=False)
    metadata.add_argument("-V", "--violin", help="Variable for violin plots (otherwise -v)", required=False)
    metadata.add_argument("-S", "--sum-viol-var", help="Variable for total violin plots (otherwise -v)", required=False)
    metadata.add_argument("-H", "--heatmap", help="Variable for heatmap (otherwise -v)", required=False)
    
    other = parser.add_argument_group("Metadata arguments")

    other.add_argument("--single-end", help="Single end reads (by default is inferred)", action="store_true", required=False)
    other.add_argument("-l", "--local-run", help="Configure for local execution", action="store_true")
    other.add_argument("--img", help="Singularity image [default: %(default)s]", default=os.environ.get("METAPHAGE_SIMG", None))
    
    other.add_argument("--tmp", dest="tempdir",help="Temporary directory [default: %(default)s]", default=os.environ.get("METAPHAGE_TMP", "/tmp") )
    other.add_argument("--work", dest="workingdir", help="Nextflow work directory [default: %(default)s]", default=os.environ.get("METAPHAGE_WD", "/tmp") )
    other.add_argument("--verbose",   action="store_true", help="Enable verbose output")
    other.add_argument("--version",   action="store_true", help="Print version and exit")
    args = parser.parse_args()


    if args.version:
        print(f"{__version__}")
        sys.exit(0)
    
    
    metaphage = Metaphage(args)

    if not metaphage.valid:
        print("ERROR: Configuration failed")
        sys.exit(1)
    else:
        if args.save:
            print(f"INFO: Saving configuration to {args.save}", file=sys.stderr)
            out = open(args.save, "w")
        else:
            out = sys.stdout


        print(metaphage.template(), file=out)
    
        if args.local_run and has_nextflow() and args.save:
            deps = 'local'
            if has_singularity():
                deps = 'singularity'
                if not os.path.exists(args.img):
                    print("ERROR: Singularity image not found: {}".format(args.img), file=sys.stderr)
                    sys.exit(1)
            elif has_docker_user():
                deps = 'docker'
            else:
                print(f"WARNING: Docker or Singularity not found, nextflow will be executed in the current environment", file=sys.stderr)
                
            print(f"Local environment: cores={CORES}, memory={MEM}GB", file=sys.stderr)
            local = Local(CORES, MEM, deps)
            print(local, file=out)
            pipeline = os.path.join(metaphage.dir, "main.nf")

            if not os.path.isfile(pipeline):
                print(f"ERROR: Nextflow pipeline not found: {pipeline}", file=sys.stderr)
                sys.exit(1)
            
            nf_run = metaphage.cmd()
            if deps == 'singularity':
                nf_run.append('-with-singularity')
                nf_run.append(args.img)

            print(f"COMMAND\n", " ".join(nf_run), file=sys.stderr)
            
            # Run 
            try:
                subprocess.check_call(nf_run)
            except subprocess.CalledProcessError as e:
                print(f"ERROR: Nextflow execution failed: {e}", file=sys.stderr)
                sys.exit(1)



     
