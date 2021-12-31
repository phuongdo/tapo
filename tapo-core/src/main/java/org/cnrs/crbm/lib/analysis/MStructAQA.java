package org.cnrs.crbm.lib.analysis;


import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Structure;
import org.cnrs.crbm.lib.io.DataIO;
import org.cnrs.crbm.lib.io.ProteinCSVReader;
import org.cnrs.crbm.lib.io.RowRepeatDB;
import org.cnrs.crbm.lib.multalign.MutilAlign;
import org.cnrs.crbm.lib.repeats.RepeatFinder;
import org.cnrs.crbm.lib.utils.PdbTools;

import java.util.List;

/**
 * This class conducted all tandem repeats from RepeatsDB and aimed to
 * build the high quality multiple structure alignment
 * Created by pdoviet on 5/18/2015.
 */
public class MStructAQA {


    public static void main(String[] args) throws Exception {

        new MStructAQA().procedue();
    }

    public void procedue() throws Exception {
        List<RowRepeatDB> rows = fetchData();
        int index = 0;
        for (RowRepeatDB row : rows) {
            if (row.getEntry().equals("1dabA")) {
                System.out.println(row);
                if (row.getAnnlevel().equals("Detailed") && row.getUnits().length() > 0 && row.getUnits().split(";").length >= 10) {
                    // System.out.println(row.getEntry());
                    try {
                        String[] units = row.getUnits().split(";");
                        String pdbCode = row.getPdbCode();
                        String pdbChain = row.getPdbChain();
                        Structure structure = PdbTools.getStructureFromLocalPdb(pdbCode);
                        Atom[] atoms = PdbTools.getAtomCAArray(structure
                                .getChainByPDB(pdbChain));
                        RepeatFinder repeatFinder = new RepeatFinder(pdbCode, pdbChain);
                        StringBuffer buffer = new StringBuffer();
                        double avgL = 0.0;
                        for (String unit : units) {
                            int starti = Integer.parseInt(unit.split("-")[0]);
                            int endi = Integer.parseInt(unit.split("-")[1]);

                            int start_pos = this.getPosition(atoms,
                                    starti);
                            int end_pos = this.getPosition(atoms, endi);
                            avgL += end_pos - start_pos + 1;

                            buffer.append(start_pos + "-" + end_pos + ";");
                        }
                        avgL = avgL / units.length;
                        String frags = buffer.toString();
                        frags = frags.substring(0, frags.length() - 1);
                        // convert to
                        MutilAlign mutilAlign = new MutilAlign(pdbCode, pdbChain, frags);
                        DataIO.writeToFile(mutilAlign.getOutputHtml(),
                                "C:/Users/pdoviet/Desktop/MSA_QA/" + pdbCode + pdbChain + "_" + index + ".html");
                        index++;
                    } catch (Exception ex) {
                        ex.printStackTrace();
                    }

                }
            }
        }

    }


    int getPosition(Atom[] atoms, int posNsq) {
        for (int i = 0; i < atoms.length; i++) {
            int seqNumber = atoms[i].getGroup().getResidueNumber().getSeqNum();
            if (seqNumber == posNsq)
                return i;
        }
        return 0;
    }


    /**
     * Fetch training data from RepeatsDB
     */
    public List<RowRepeatDB> fetchData() throws Exception {

        ProteinCSVReader csvReader = new ProteinCSVReader();
        List<RowRepeatDB> rows = csvReader.getRepeatDB("data/RDB-dataset.tab");

        return rows;

    }


}
