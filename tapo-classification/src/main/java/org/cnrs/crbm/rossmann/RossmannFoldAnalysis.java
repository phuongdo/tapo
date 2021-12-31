package org.cnrs.crbm.rossmann;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureTools;
import org.cnrs.crbm.lib.conf.Dir;
import org.cnrs.crbm.lib.dssp.DSSP;
import org.cnrs.crbm.lib.io.ReadFasta;
import org.cnrs.crbm.lib.repeats.Fragement;
import org.cnrs.crbm.lib.repeats.module.ProVector;
import org.cnrs.crbm.lib.repeats.Superimposer;
import org.cnrs.crbm.lib.repeats.module.VectorModule;
import org.cnrs.crbm.lib.repeats.module.VectorShape;
import org.cnrs.crbm.lib.sadb.ConformationMatcher;
import org.cnrs.crbm.lib.trsfinder.Repeat;
import org.cnrs.crbm.lib.trsfinder.RepeatContent;
import org.cnrs.crbm.lib.utils.PdbTools;
import org.cnrs.crbm.lib.utils.ProgressBar;

public class RossmannFoldAnalysis {

    public static void main(String[] args) {
        // pdb_seqres-part4.nrdb


        if (args.length > 0) {
            new RossmannFoldAnalysis().scanToFindPattern(args[0]);
        } else {
            new RossmannFoldAnalysis().scanToFindPattern("data/PDB_July_1_2011.txt");
//            new RossmannFoldAnalysis().scanToFindPattern("data/resultLargeScale/TAPO_NOTRs_07.cdhit");
        }

    }

    public void scanToFindPattern(String fastaDir) {

        ReadFasta readFasta = new ReadFasta();

        // run new version of dssp binary

        ProgressBar bar = new ProgressBar();

        // System.out.println(Dir.FASTA_LOCAL + "\n");

        // pdb_seqres-part4.nrdb
        try {
            Set<String> pdbs = readFasta.getAllPdbsWithChain(fastaDir);

            PrintWriter writer = new PrintWriter(Dir.TMP_DIR + "/" + fastaDir.replace("/","_"),
                    "UTF-8");

            System.out.println(Dir.TMP_DIR + "/" + fastaDir.replace("/","_"));

            System.out.println("Process Starts Now!");

            bar.update(0, pdbs.size());
            int i = 1;

            for (String pdb : pdbs) {

                // System.out.println(pdb);
                // update a status bar
                bar.update(i, pdbs.size());
                try {
                    String pdbCode = pdb.substring(0, 4);
                    String pdbChain = pdb.substring(5, 6);
                    Structure structure = PdbTools
                            .getStructureFromLocalPdb(pdbCode);
                    // Atom[] cbAtoms = PdbTools.getAtomCBArray(StructureTools
                    // .getAtomCAArray(structure.getChainByPDB(pdbChain)));
                    //
                    Atom[] atoms = StructureTools.getAtomCAArray(structure
                            .getChainByPDB(pdbChain));
                    if(atoms.length>150) {


                        DSSP dssp = new DSSP(pdbCode);
                        // make a filter first
                        // String secondaryStr = dssp.filterSS(dssp.getSS(atoms));
                        String secondaryStr = dssp.filterSS(dssp.getCombinedSS(
                                atoms, pdbChain));
                        String ssPattern = VectorShape.getSSPattern(secondaryStr);
                        // System.out.println(ssPattern);

//                    List<ProVector> vectors = VectorShape.getVectors(atoms,
//                            secondaryStr);
                        VectorModule vectorModule = new VectorModule();
                        //vectors = VectorShape.getVectors(atoms, strSS);
                        /**
                         * MAJOR UPDATE from 26 May 2015
                         */
                        List<ProVector> vectors = vectorModule.extractVectors(atoms, secondaryStr);

                        if (ConformationMatcher.isRossmanfoldPattern(ssPattern)
                                && checkRossmanFoldPattern(vectors, atoms))
                            // System.out.println(pdb);
                            writer.write(pdb + "\n");
                    }
                } catch (Exception ex) {
                    ex.printStackTrace();

                }
                i++;

            }
            writer.close();

            System.out.println("Process Completed!");

        } catch (Exception e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }

        // String pattern = "HHHHBHBHBHBHHH";
        //
        // pattern = "BBBBBBHBHBHBJJJ";
        // System.out.println(ConformationMatcher.isRossmanfoldPattern(pattern));
    }

    public boolean checkRossmanFoldPattern(List<ProVector> vectors, Atom[] atoms)
            throws Exception {

        int winsize = 2;
        // if (winsize == 3)
        // System.out.println("winsize:" + winsize);
        // vectors = ProVector.toSecondaryVector(vectors);

        Superimposer superimpose = new Superimposer();

        List<List<ProVector>> list = new ArrayList<List<ProVector>>();
        for (int i = 0; i < vectors.size() - winsize + 1; i++) {
            List<ProVector> sub_list = new ArrayList<ProVector>();
            for (int j = 0; j < winsize; j++) {
                ProVector v = vectors.get(i + j);
                sub_list.add(v);
            }
            list.add(sub_list);
        }
        int index = 0;
        int max_point = 0;
        // List<RepeatScore> repeat_scores = new ArrayList<RepeatScore>();
        Atom[] as1 = null;
        Atom[] as2 = null;

        int lastTRs = -1;
        for (int j = 0; j < list.size(); j++) {

            if (lastTRs >= 0 && lastTRs + winsize > j) {
                j = lastTRs + winsize - 1;
                continue;

            }
            List<ProVector> sub_list_ref = list.get(j);
            String pattern_ref = this.getPatterOfVectors(sub_list_ref);
            if (!pattern_ref.equals("HB") || pattern_ref.equals("BH"))
                continue;
            int match = 0;
            // System.out.println(">" + j);

            int store_current_posstion = j;

            as1 = Fragement.getFragementsofAtoms(atoms, sub_list_ref.get(0)
                    .getPosStart(), sub_list_ref.get(sub_list_ref.size() - 1)
                    .getPosEnd());

            int start1 = sub_list_ref.get(0).getPosStart();
            int end1 = sub_list_ref.get(sub_list_ref.size() - 1).getPosEnd();
            // scan to the right
            Repeat repeat = new Repeat();
            repeat.getRepeats().add(new RepeatContent(start1, end1));

            for (int i = j + winsize; i < list.size(); ) {
                List<ProVector> sub_list = list.get(i);
                as2 = Fragement.getFragementsofAtoms(atoms, sub_list.get(0)
                        .getPosStart(), sub_list.get(sub_list.size() - 1)
                        .getPosEnd());

                int start2 = sub_list.get(0).getPosStart();

                if (as1.length > 40
                        && (start2 - end1 + 1) > (double) 0.2 * as1.length)

                    break;// no check anymore

                else if (as1.length < 20
                        && (start2 - end1 + 1) > (double) 0.5 * as1.length)
                    break;

                // check pattern
                String pattern_sub = this.getPatterOfVectors(sub_list);

                double score = superimpose.compareListVector(sub_list_ref,
                        sub_list);

                // System.out.println(score);

                if (score > 0.4 && pattern_sub.equals(pattern_ref)) {

                    match++;
                    store_current_posstion = i;
                    end1 = sub_list.get(sub_list.size() - 1).getPosEnd();
                    i = i + winsize;
                } else {
                    i = i + 1;
                }
                if (match >= 3)
                    return true;

            }

        }

        return false;

    }

    private String getPatterOfVectors(List<ProVector> sub_list) {
        String pattern_sub = "";
        for (ProVector v : sub_list) {
            pattern_sub += v.getType();
        }
        return pattern_sub;
    }
}
