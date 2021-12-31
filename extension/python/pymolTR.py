# -*- coding: iso-8859-15 -*-

import __main__
__main__.pymol_argv = [ 'pymol', '-qc'] # Quiet and no GUI

colors = ["red", "blue", "orange", "cyan"]
import sys, time, os
import pymol
pymol.finish_launching()

pdbId = sys.argv[1]
chain = sys.argv[2]
units = sys.argv[3].split(";")
bound = sys.argv[4]



#feedback disable, opengl, warnings
pymol.cmd.feedback("disable","opengl","warnings")
pymol.cmd.fetch(pdbId)
pymol.cmd.bg_color("black")
pymol.cmd.hide("all")
pymol.cmd.show("cartoon", "chain "+chain)
c=0
for unit in units:
	#start = unit.split("-")[0]
	#end = unit.split("-")[1]
	selection = "resi  "+unit
	pymol.cmd.select("repeat", selection)
	pymol.cmd.color (colors[c%3], "repeat")
	c=c+1

pymol.cmd.save("data/"+pdbId+"_"+chain+"_"+bound+".pse");
#pymol.cmd.png("data/"+pdbId+"_"+chain+"_"+bound+".png")
#pymol.cmd.save("C://1lxa.pse");
pymol.cmd.quit()



