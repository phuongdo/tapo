package org.cnrs.crbm.nextgen;

import org.biojava.nbio.structure.Atom;
import org.cnrs.crbm.lib.io.DataIO;
import org.cnrs.crbm.lib.io.ReadFasta;
import org.cnrs.crbm.lib.io.TaPoFastaFormat;
import org.cnrs.crbm.lib.repeats.RepeatFinder;
import org.cnrs.crbm.lib.trsfinder.Repeat;
import org.cnrs.crbm.lib.trsfinder.RepeatContent;
import org.cnrs.crbm.lib.utils.NumberFormatUtils;
import org.cnrs.crbm.lib.utils.PdbTools;
import org.cnrs.crbm.trclassification.CATHClassifier;
import org.cnrs.crbm.trclassification.ProteinCATH;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * Created by pdoviet on 10/7/2015.
 */
public class ConvertToPymol {

    public static void main(String[] args) {

        String fileIn = "input/test.in";
        ConvertToPymol convertToPymol = new ConvertToPymol();
        try {
            convertToPymol.generatePymolScript(fileIn);
        } catch (Exception e) {
            e.printStackTrace();
        }

    }


    private void generatePymolScript(String fileDir) throws Exception {
        // read tapo output
        ReadFasta readFasta = new ReadFasta();
        Map<String, TaPoFastaFormat> fastaTRs = readFasta.readTaPoFastaFormat("data/resultLargeScale/tapo-v1.1.2/tapo_PDBTRs_June2015_v1.1.2.out");
        // read data to test
        List<String> pdbs = DataIO.readLines(fileDir);
        for (String pdb : pdbs) {


            if (fastaTRs.containsKey(pdb)) {
                TaPoFastaFormat taPoFastaFormat = fastaTRs.get(pdb);
                if (taPoFastaFormat.is3DRepeat()) {
                    String pdbCode = pdb.substring(0, 4);
                    String pdbChain = pdb.substring(5, 6);
                    RepeatFinder repeatFinder = new RepeatFinder(pdbCode, pdbChain);
                    Atom[] atoms = repeatFinder.getAtoms();

                    List<Repeat> lstRepeat = new ArrayList<Repeat>();
                    for (Repeat repeat : taPoFastaFormat.getRepeats()) {
                        if (repeat.getCluster().contains("selected")) {
                            lstRepeat.add(repeat);
                        }
                    }
                    for (Repeat repeat : lstRepeat) {
                        int seqStart = PdbTools.getResSeq(atoms, repeat.getStart());
                        int seqEnd = PdbTools.getResSeq(atoms, repeat.getEnd());
                        String region = seqStart + "_" + seqEnd;
                        // build file
                        StringBuffer buffer = new StringBuffer();
                        buffer.append("load http://www.rcsb.org/pdb/files/" + pdbCode.toUpperCase() + ".pdb\n");
                        buffer.append("hide\n");
                        buffer.append("show cartoon, chain " + pdbChain + "\n");
                        for (RepeatContent rc : repeat.getRepeats()) {
                            int u_start = PdbTools.getResSeq(atoms, rc.getStart());
                            int u_end = PdbTools.getResSeq(atoms, rc.getEnd());
                            buffer.append("color auto, resi " + u_start + "-" + u_end + "\n");

                        }
                        String fileName = pdb + "." + region + ".pml";
                        DataIO.writeToFile(buffer.toString(), "output/pymol/" + fileName);


                    }


                }

            }

        }


    }


}
