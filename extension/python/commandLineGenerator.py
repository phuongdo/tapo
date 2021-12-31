import argparse
import os.path
from subprocess import Popen


parser = argparse.ArgumentParser(description='Process the arguments.')
parser.add_argument('--input', type=argparse.FileType('rU'),
                   help='metaServerOutput to process')
args = parser.parse_args()


#"protein,begin,end,msa"
outputH = open(str(os.path.basename(args.input.name)) + "commanLineForPymol", "w")

header = 1

for line in (args.input):
    if header == 1:
        header = 0
        continue
    else:
        lineC = line.replace('"',"",)
        TR = lineC.strip().split(",")
        #print(TR)
        id = TR[0].split("_")
        
        bound = ""
        
        for motif in TR[3].split(" "):
            motif = motif.strip().replace("-","",)
            bound = bound + str(len(motif)) + " "
        
        outputH.write("pymolColorByMotifByPortionAndSave.py " + str(id[0]) + " " + str(id[1]) + " " + str(TR[1]) + " " + str(bound) +" \n")
