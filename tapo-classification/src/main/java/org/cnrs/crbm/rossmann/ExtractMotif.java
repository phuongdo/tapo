package org.cnrs.crbm.rossmann;


import org.biojava.nbio.structure.Atom;
import org.cnrs.crbm.lib.io.DataIO;
import org.cnrs.crbm.lib.repeats.RepeatFinder;
import org.cnrs.crbm.lib.utils.PdbTools;

import java.util.List;

/**
 * Created by pdoviet on 8/20/2015.
 */
public class ExtractMotif {

    public static void main(String[] args) throws Exception {

        new ExtractMotif().printFASTA();

    }

    public void printFASTA() {

        List<String> lines = DataIO.readLines("data/rossmann/motif/Ross-TypeI.txt");

        for (String line : lines) {
            String[] cols = line.split("\t");
            String pdb = cols[0];
            String[] units = cols[1].split(";");

            String pdbCode = pdb.substring(0, 4);
            String pdbChain = pdb.substring(5, 6);

            RepeatFinder repeatFinder = new RepeatFinder(pdbCode, pdbChain);
            Atom[] atoms = repeatFinder.getAtoms();
            String strSS = repeatFinder.getStrSS();
            String strSeq = repeatFinder.getStrSeq();

            for (String unit : units) {

                System.out.println(">" + pdb + "|" + unit);

                int start = Integer.parseInt(unit.split("-")[0]);
                int end = Integer.parseInt(unit.split("-")[1]);
                start = PdbTools.getPosition(atoms, start);
                end = PdbTools.getPosition(atoms, end);
                System.out.println(strSeq.substring(start, end));

            }

        }

    }
}
