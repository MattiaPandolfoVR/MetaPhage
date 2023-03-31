#!/usr/bin/env python3
# Coded by Andrea Telatin (Andrea.Telatin@quadram.ac.uk)
# Implemented in MetaPhage by Mattia Pandolfo (mattia.pandolfo@univr.it)
software = "db_manager.py"
version = "0.2.1"



import json, argparse, os, time, re, sys
import concurrent.futures

# Try importing wget
withWget = False
try:
  import wget
  withWget = True
except ImportError:
  withWget = False
  print(f"WARNING: {software} will use `wget` from the system.", file=sys.stderr)


start = time.perf_counter()
progresses = {"_printed": []}

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

# Initial     #"url":    "https://zenodo.org/record/5879332/files/2022-01-inphared.tar.gz?download=1",

config = """
{
  "2021.1": {
    "phix": {
      "name": "PhiX reference",
      "url": "https://zenodo.org/record/4608203/files/phix.tar.gz?download=1",
      "md5": "36897778387b36cecd6e60aef9dab17b",
      "file": "phix.tar.gz",
      "provides": ["phix"],
      "expand": ["tar xfz '{file}' --directory '{outdir}'", "rm '{file}'"],
      "cleanup": "rm phix.*"
    },  
    "kraken2": {
      "name": "MiniKraken",
      "url": "https://zenodo.org/record/4608203/files/kraken2.tar.gz?download=1",
      "md5": "4dca0ee1acccbff860ca725c37e22f85",
      "file": "kraken2.tar.gz",
      "provides": ["kraken2"],
      "expand": ["tar xfz '{file}' --directory '{outdir}'", "rm '{file}'"],
      "cleanup": "rm kraken2.*"
    }, 
    "vibrant": {
      "name": "Vibrant",
      "url": "https://zenodo.org/record/4608203/files/vibrant.tar.gz?download=1",
      "md5": "05bad36563db58889236571af871acaa",
      "file": "vibrant.tar.gz",
      "provides": ["vibrant"],
      "expand": ["tar xfz '{file}' --directory '{outdir}'", "rm '{file}'"],
      "cleanup": "rm vibrant.*"
    },
    "virsorter": {
      "name": "Virsorter",
      "url": "https://zenodo.org/record/4608203/files/virsorter_legacy.tar.gz?download=1",
      "md5": "2e4308be78a15df3d89ceb3e72c9ba81",
      "file": "virsorter_legacy.tar.gz",
      "provides": ["virsorter"],
      "expand": ["tar xfz '{file}' --directory '{outdir}'", "rm '{file}'"],
      "cleanup": "rm virsorter_legacy.*"
    }, 
      "phigaro": {
      "name": "Phigaro",
      "url": "https://zenodo.org/record/4608203/files/phigaro.tar.gz?download=1",
      "md5": "f6e671d885538542755d0c52b1aa6ddb",
      "file": "phigaro.tar.gz",
      "provides": ["phigaro"],
      "expand": ["tar xfz '{file}' --directory '{outdir}'", "rm '{file}'"],
      "cleanup": "rm phigaro.*"
    },
    "vcontact2": {
      "name":    "vConTACT2",
      "url":     "https://warwick.s3.climb.ac.uk/ifrqmra-metaphage/v1.0/2022-01-inphared.tar.gz",
      "md5":     "72c9a0be3b93364e44338ced659341fe",
      "file":    "2022-01-inphared.tar.gz",
      "provides": ["inphared"],
      "expand":   ["tar xfz '{file}' --directory '{outdir}'", "rm '{file}'"],
      "cleanup":  "rm *inphared*.tar.gz"
    }
  },
  "2022.1": {
    
    "phix": {
      "name": "PhiX reference",
      "url": "https://zenodo.org/record/4608203/files/phix.tar.gz?download=1",
      "md5": "36897778387b36cecd6e60aef9dab17b",
      "file": "phix.tar.gz",
      "provides": ["phix"],
      "expand": ["tar xfz '{file}' --directory '{outdir}'", "rm '{file}'"],
      "cleanup": "rm phix.*"
    },
    "kraken2": {
      "name": "MiniKraken",
      "url": "https://zenodo.org/record/4608203/files/kraken2.tar.gz?download=1",
      "md5": "4dca0ee1acccbff860ca725c37e22f85",
      "file": "kraken2.tar.gz",
      "provides": ["kraken2"],
      "expand": ["tar xfz '{file}' --directory '{outdir}'", "rm '{file}'"],
      "cleanup": "rm kraken2.*"
    },
      "phigaro": {
      "name": "Phigaro",
      "url": "https://zenodo.org/record/4608203/files/phigaro.tar.gz?download=1",
      "md5": "f6e671d885538542755d0c52b1aa6ddb",
      "file": "phigaro.tar.gz",
      "provides": ["phigaro"],
      "expand": ["tar xfz '{file}' --directory '{outdir}'", "rm '{file}'"],
      "cleanup": "rm phigaro.*"
    },
    "vibrant": {
      "name": "Vibrant",
      "url": "https://zenodo.org/record/4608203/files/vibrant.tar.gz?download=1",
      "md5": "05bad36563db58889236571af871acaa",
      "file": "vibrant.tar.gz",
      "provides": ["vibrant"],
      "expand": ["tar xfz '{file}' --directory '{outdir}'", "rm '{file}'"],
      "cleanup": "rm vibrant.*"
    },
    "virsorter2": {
      "name": "Virsorter2",
      "url": "https://s3.climb.ac.uk/ifrqmra-metaphage/v2.0/virsorter2.tar.gz",
      "md5": "ddc75b770b96a9fabee4f76458b57e0d",
      "file": "virsorter2.tar.gz",
      "provides": ["virsorter"],
      "expand": ["tar xfz '{file}' --directory '{outdir}'", "rm '{file}'"],
      "cleanup": "rm virsorter2.tar.gz"
    }, 
    "checkv": {
      "name": "CheckV database",
      "url": "https://s3.climb.ac.uk/ifrqmra-metaphage/v2.0/checkv.tar.gz",
      "md5": "85095275aee479de87b0dbc95b3ca48e",
      "file": "checkv.tar.gz",
      "provides": ["checkv"],
      "expand": ["tar xfz '{file}' --directory '{outdir}'", "rm '{file}'"],
      "cleanup": "rm checkv.tar.gz"
    },
    "vcontact2": {
      "name":    "vConTACT2",
      "url":     "https://warwick.s3.climb.ac.uk/ifrqmra-metaphage/v1.0/2022-01-inphared.tar.gz",
      "md5":     "72c9a0be3b93364e44338ced659341fe",
      "file":    "inphared.tar.gz",
      "provides": ["inphared"],
      "expand":   ["tar xfz '{file}' --directory '{outdir}'", "rm '{file}'"],
      "cleanup":  "rm *inphared*.tar.gz"
    }
  
  }
}
"""

try:
  data = json.loads(config)
except Exception as e:
  eprint(f"Invalid JSON database: {e}")
  quit()


def fill_template(template, source):
  matches = re.findall(r'\{(.*?)\}', template)
  for m in matches:
    if m in source:
      placeholder = '{' + m + '}'
      template = re.sub(rf"{placeholder}", source[m], template)
    else:
      print(f"Error: unable to put {source[m]} in {m} [processing {source[url]}")
      quit()
  return template

def globalbar(cur, tot, step=100):
  progresses[tot] = cur
  #print(tot, "...")
  gt = 0
  done = 0
  for p in progresses:
    if p == '_printed':
      continue
    gt += int(p)
    done += progresses[p]
  ratio = int(done*100/gt)

  if ratio > 0 and len(progresses) == len(data) and ratio % 10 == 0 and not ratio in progresses['_printed']:
    progresses["_printed"].append(ratio)
    time.sleep(0.1)
    eprint(f"Downloading: {ratio}% done")

def downloadPackage(url, dest):
  if withWget:
    try:
      wget.download(url, bar=globalbar)
    except Exception as e:
      eprint(f"Error downloading {url}:\n{e}")
  else:
    try:
      os.system(f"wget -q -O '{dest}' '{url}'")
    except Exception as e:
      eprint(f"Error downloading {url}:\n{e}")

  
def getDatabase(package, attempts=5):
  starttime = time.perf_counter()
  package["downloaded"] = False
  shell_commands = []
  for cmd_template in package["expand"]:
    cmd = fill_template(cmd_template, package)
    shell_commands.append(cmd)
  eprint(f" 📦  Preparing to download {package['name']}") 
  
  for i in range(attempts + 1):
    try:
      downloadPackage(package["url"], package["file"])
      package["downloaded"] = True
      break
    except Exception as e:
      if 'cleanup' in package:
        os.system(package["cleanup"])
      eprint(f"{i}/{attempts} Download of {package['name']} failed:\nURI: {package['url']}\nError: {e}")
      if i == attempts:
        return  f" 🛑 {package['name']} failed: {e}"
  eprint(f" ✅  {package['name']} downloaded")
 
  for cmd in shell_commands:
    os.system(cmd)
  
  finishtime = time.perf_counter()
 
  return finishtime-starttime 

if __name__ == "__main__":

    # Parameter parser
    parser = argparse.ArgumentParser(
        description = 'This script is designed to manage all the databases used '
                      'in MetaPhage. This script is placed in MetaPhage/bin.',
        formatter_class = argparse.RawTextHelpFormatter)
    options = parser.add_argument_group("Options")
     
    options.add_argument('-o', '--outdir', help="Output directory", required=True)
    options.add_argument('-r', '--release', help="Database bundle release [default: %(default)s]", default="2021.1")
    options.add_argument('-m', '--maxthreads', help="Maximum number of concurrent downloads", default=3)
    parameters = parser.parse_args()

    # output directory
    if not os.path.exists(parameters.outdir):
      os.makedirs(parameters.outdir)
    
    # To download
    todo_list = []

    if parameters.release not in data:
        print("Error: `{}` is not a valid release. Try:\n * {}".format(parameters.release, "\n * ".join(data.keys())))
        exit(1)

    data = data[parameters.release]
    print(" 📂  Downloading bundle {}: {} databases total".format(parameters.release, len(data)))
    for package_name in data:
      pack = data[package_name]
      #print("▸ ",package_name, "\t", pack["url"])
      download_this = 0

      # Check requried files
      for file in pack["provides"]:
        if not os.path.exists(os.path.join(parameters.outdir, file)):
          download_this += 1
        else:
          eprint(f"{file} found: skipping")
      
      if download_this > 0:
        pack["outdir"] = parameters.outdir
        todo_list.append(pack)

    # Download  
    results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=int(parameters.maxthreads)) as executor:
        for record in todo_list:
          results.append(executor.submit(getDatabase, record))
    
    totalTime = 0
    resultsCount = 0
    for f in concurrent.futures.as_completed(results):
      #print(str(f.result()))
      totalTime += f.result() 
      resultsCount+=1

    finish = time.perf_counter()
    print(f"{resultsCount} packages downloaded in {round(finish-start, 2)} seconds (cumulative {round(totalTime,2)} seconds)")

    # Pass db paths to nextflow
    path_phix = parameters.outdir + "/phix"
    path_kraken2 = parameters.outdir + "/kraken2"
    path_vibrant = parameters.outdir + "/vibrant"
    path_virsorter2 = parameters.outdir + "/virsorter"
    path_checkv = parameters.outdir + "/checkv"
    path_allvsall = parameters.outdir + "/diamond"
    path_vcontact2 = parameters.outdir + "/inphared"

    file_db_path = open("db_path.csv", "w")
    file_db_path.writelines("%s,%s,%s,%s,%s,%s,%s" % (path_phix, path_kraken2, path_vibrant, path_virsorter2, path_checkv, path_allvsall, path_vcontact2))
    file_db_path.close()
