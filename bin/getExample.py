#!/usr/bin/env python3
# Andrea Telatin, 2022
import os, sys
import subprocess
import argparse
import concurrent.futures

def getUris(urlTsvFile):
    """
    Given a tsv file associate "run_accession" column with "fastq_ftp"
    """
    uris = {}
    try:
        with open(urlTsvFile, "r") as f:
            for line in f:
                line = line.strip()
                if line.startswith("#") or line.startswith("study_accession"):
                    continue
                cols = line.split("\t")
                if len(cols) < 6:
                    continue
                uris[cols[3]] = cols[6].split(";")
        return uris
    except Exception as e:
        print("Error in getUris: %s" % e)
        sys.exit(1)

def readMd5(filename):
    """
    Read the md5 from a file containing MD5 and filename
    """
    try:
        with open(filename, "r") as f:
            for line in f:
                line = line.strip()
                fields = line.split(" ")
                return fields[0]
    except Exception as e:
        print("Error in readMd5: %s" % e)
        sys.exit(1)

def getSamples(metadataFile):
    """
    Metadata file (csv), get the first sample "Sample"
    """
    samples = []
    try:
        with open(metadataFile, "r") as f:
            for line in f:
                line = line.strip()
                if line.startswith("#") or line.startswith("Sample,"):
                    continue
                cols = line.split(",")
                samples.append(cols[0])
        return samples
    except Exception as e:
        return None

def downloadFileToDest(URL, destination, verbose=False):
    """
    Retrieve a file from an URL and save it to destination
    """
    cmd = ["wget", "--quiet", "-O", destination, URL]
    # Check if file is already there
    if os.path.exists(destination) and os.path.exists( os.path.join( os.path.dirname(destination) , "md5", os.path.basename(destination))):
        md5sum = readMd5( os.path.join( os.path.dirname(destination) , "md5", os.path.basename(destination)) )
        if not md5sum is None:
            realmd5 = subprocess.check_output(["md5sum", destination]).split()[0].decode("utf-8")
            if md5sum == realmd5:
                if verbose:
                    print("File %s already exists and is valid." % destination)
                return destination
            else:
                print("File {} already exists but MD5 mismatches: {} vs {}".format(destination, md5sum, realmd5))
        else:
            print("Invalid MD5 file for %s" % destination)

    try:
        subprocess.check_call(cmd)
        return destination
    except subprocess.CalledProcessError as e:
        print("Error downloading file: %s" % e)
        sys.exit(1)
    except Exception as e:
        print("Error in downloadFileToDest: %s" % e)
        sys.exit(1)

if __name__ == "__main__":
    default_outdir = os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), "demo")
    default_metadata = os.path.join(default_outdir, "infant-metadata.csv")
    default_uris = os.path.join(default_outdir, "uris.tsv")

    args = argparse.ArgumentParser(description="Downloads a demo dataset")
    args.add_argument("-o", "--outdir",   help="Output directory [default: %(default)s]", default=default_outdir)
    args.add_argument("-m", "--metadata", help="Metadata file [default: %(default)s]", default=default_metadata)
    args.add_argument("-u", "--urls",     help="EBI url file tsv [default: %(default)s]", default=default_uris)
    
    args.add_argument("-t", "--threads",  help="Maximum number of concurrent downloads [default: %(default)s]", default=3)
    args.add_argument("-v", "--verbose",  help="Print the commands as they are executed", action="store_true")
    args = args.parse_args()

    # Check metadata and uris
    if not os.path.exists(args.metadata):
        print("Metadata file not found: %s" % args.metadata, file=sys.stderr)
        sys.exit(1)
    if not os.path.exists(args.urls):
        print("URLs file not found: %s" % args.urls, file=sys.stderr)
        sys.exit(1)
    
    print("-- EXAMPLE DOWNLOADER --")
    URLs = getUris(args.urls)
    if args.verbose:
        print(" * {} URLs found in {}".format(len(URLs), args.urls))
    
    SampleList = getSamples(args.metadata)
    if args.verbose:
        print(" * {} samples found in {}".format(len(SampleList), args.metadata))
    # Check if the output directory exists
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)


    results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=int(args.threads)) as executor:
        for index, sample in enumerate(SampleList):
            for url in URLs[sample]:
                if args.verbose:
                    print(" * Sample {}: Queueing {}".format(index, url))
                filename = os.path.basename(url)
                destination = os.path.join(args.outdir, filename)
                #downloadFileToDest(url, destination, args.verbose)
                results.append(executor.submit(downloadFileToDest, url, destination, args.verbose))

    
    for f in concurrent.futures.as_completed(results):
      print("Downloaded: ", str(f.result()))
      