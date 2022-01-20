#!/usr/bin/env python3
# Coded by Andrea Telatin (Andrea.Telatin@quadram.ac.uk)
# Implemented in MetaPhage by Mattia Pandolfo (mattia.pandolfo@gmail.com)
software = "db_manager.py"
version = "0.0.2"

import wget, json, argparse, os, time, re, sys
import concurrent.futures

start = time.perf_counter()
progresses = {"_printed": []}

config = """
{
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
    "name": "vConTACT2",
    "url": "https://zenodo.org/record/5879332/files/2022-01-inphared.tar.gz?download=1",
    "md5": "8c66e1b0c8359dff2a11a407826efa02",
    "file": "inphared.tar.gz",
    "provides": ["inphared"],
    "expand": ["tar xfz '{file}' --directory '{outdir}'", "rm '{file}'"],
    "cleanup": "rm inphared.*"
  }
}
"""

try:
  data = json.loads(config)
except Exception as e:
  eprint(f"Invalid JSON database: {e}")
  quit()

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs) 
 
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
      wget.download(package["url"], bar=globalbar )
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
    options.add_argument('-m', '--maxthreads', help="Maximum number of concurrent downloads", default=3)
    parameters = parser.parse_args()

    # output directory
    if not os.path.exists(parameters.outdir):
      os.makedirs(parameters.outdir)
    
    # To download
    todo_list = []
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
    path_virsorter = parameters.outdir + "/virsorter"
    path_vcontact2 = parameters.outdir + "/inphared"

    file_db_path = open("db_path.csv", "a")
    file_db_path.writelines("%s,%s,%s,%s,%s" % (path_phix, path_kraken2, path_vibrant, path_virsorter, path_vcontact2))
    file_db_path.close()
