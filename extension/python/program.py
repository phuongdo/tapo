import os
def isNoneOrEmptyOrBlankString (myString):
        if myString:
            if not myString.strip():
                return True
        else:
            return True

        return False

with open('file.in') as fp:
			for line in fp:
				if not isNoneOrEmptyOrBlankString(line):
					# do something here
					rows = line.split("\t")
					pdbId = rows[0].split("_")[0]		
					pdbChain = rows[0].split("_")[1]
					unit= rows[6]
					bound = rows[5]
					#print "pymolTR.py "+pdbId+ " " + pdbChain+" "+ unit +" " + bound
					os.system("pymolTR.py "+pdbId+ " " + pdbChain+ " "+ unit + " " + bound)
					
					#print rows
					
				
			
			
