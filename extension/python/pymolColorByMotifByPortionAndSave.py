#!/usr/bin/python
import __main__
__main__.pymol_argv = [ 'pymol', '-qc'] # Quiet and no GUI

import sys, time, os
import pymol
 
pymol.finish_launching()

pdbId = sys.argv[1]
chain = sys.argv[2]
begin = int(sys.argv[3])
portions = []
for i in range(4,len(sys.argv)):
    portions.append(int(sys.argv[i]))

colors = ["red", "blue", "orange", "cyan"]


pymol.cmd.fetch(pdbId)
pymol.cmd.hide("all")
pymol.cmd.show("cartoon", "chain " + str(chain))

c = 0

for i in range(0,len(portions)):
    selection = "resi " + str(begin) + "-" + str(begin+portions[i]-1)
    pymol.cmd.select("repeat", selection)
    pymol.cmd.color (colors[c], "repeat")
    pymol.cmd.deselect()
    begin = begin+portions[i]
    c = c+1
    if c == 3:
        c = 0

pymol.cmd.set("seq_view", 1)

pymol.cmd.save(pdbId + "_" + str(chain) + "_" + str(sys.argv[3]) + "_" + "ByMotif.pse")
pymol.cmd.quit()
#pymol.cmd.select("clean")


#from pymol import cmd
# 
#def colorSelec( begin, end ):
#    '''
#DESCRIPTION
# 
#    Brief description what this function does goes here
#    '''
#   
#    cmd.fetch("3q49")
#    selection = "resi" + begin + " - " + end
#    cmd.select("repeat", selection)
#    cmd.color ('red', "repeat")
#
#
#    #print "Hello, PyMOLers"
#    #print "You passed in %s and %s" % (arg1, arg2)
#    #print "I will return them to you in a list.  Here you go."
#    #return (arg1, arg2)
# 
#cmd.extend( "colorSelec", colorSelec );



 
