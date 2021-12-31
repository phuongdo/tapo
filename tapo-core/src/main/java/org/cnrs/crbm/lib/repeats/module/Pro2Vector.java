package org.cnrs.crbm.lib.repeats.module;

import nptr.Copy;
import nptr.Repeat;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureTools;
import org.cnrs.crbm.lib.dssp.DSSP;
import org.cnrs.crbm.lib.repeats.*;
import org.cnrs.crbm.lib.seqalign.TREksWapper;
import org.cnrs.crbm.lib.trsfinder.RepeatContent;
import org.cnrs.crbm.lib.utils.NumberFormatUtils;
import org.cnrs.crbm.lib.utils.PdbTools;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

public class Pro2Vector {

    private static final double RMSD_THRESHOLD = 4.0;
    static Superimposer impose = new Superimposer();
    boolean debug = false;
    //static Logger logger = Logger.getLogger(Pro2Vector.class);
    static Logger logger = LoggerFactory.getLogger(Pro2Vector.class);

    public static void main(String[] args) throws Exception {

        // defaults
        //String pdbCode = "1fdn".toUpperCase();
        //String pdbCode = "2afg".toLowerCase();
        String pdbCode = "3ayr".toLowerCase();
        String pdbChain = "A";
        int repeatLeng = 0;

        boolean output = true;

        // int winsize = 2;
        if (args.length > 0) {
            pdbCode = args[0];
            pdbChain = args[1];
            repeatLeng = Integer.parseInt(args[2]);
            // winsize = Integer.parseInt(args[2]);
        }

        for (String arg : args) {
            // if (arg.equals("-v"))
            // debug = true;

            if (arg.startsWith("-o"))
                output = true;

        }

        // String pdbFile = PdbTools.downloadPDB(pdbCode);
        Structure structure = PdbTools.getStructureFromLocalPdb(pdbCode);
        // Atom[] cbAtoms = PdbTools.getAtomCBArray(StructureTools
        // .getAtomCAArray(structure.getChainByPDB(pdbChain)));
        //
        Atom[] atoms = StructureTools.getAtomCAArray(structure
                .getChainByPDB(pdbChain));
        DSSP dssp = new DSSP(pdbCode);
        // make a filter first
        // String secondaryStr = dssp.filterSS(dssp.getSS(atoms));
        String secondaryStr = dssp
                .filterSS(dssp.getCombinedSS(atoms, pdbChain));

        System.out.println(secondaryStr);
        // System.out.println(dssp.filterSS(dssp.getCombinedSS(atoms,
        // pdbChain)));

        // System.out.println(secondaryStr);
        Pro2Vector pro2Vector = new Pro2Vector();
//        List<ProVector> vectors = ProVector.toSecondaryVector(VectorShape
//                .getVectors(atoms, secondaryStr));
//
//
//        List<ProVector> betaSheets = new ArrayList<ProVector>();
//        for (ProVector proVector: vectors){
//            if(proVector.getType().equals("B")){
//                betaSheets.add(proVector);
//            }
//        }
//
//        ProVector seed = betaSheets.get(0);
        List<ProVector> vectors = VectorShape
                .getVectors(atoms, secondaryStr);
        List<ProVector> vectorsStored = new ArrayList<ProVector>();


        // forward algorithm
        for (int i = 0; i < vectors.size() - 1; i++) {

            //vectorsStored.add(vectors.get(i));
            // merge two vectors if they have the same direction.
            ProVector vector = vectors.get(i);
            if (vector.getType().equals("H")) {


            }


//            if (vectors.get(i).getPosEnd() != vectors.get(i + 1).getPosStart()) {
//
//                ProVector vectorL = new ProVector();
//                vectorL.setType("L");
//                vectorL.setPosStart(vectors.get(i).getPosEnd());
//                vectorL.setPosEnd(vectors.get(i + 1).getPosStart());
//                vectorL.setStart(vectors.get(i).getEnd());
//                vectorL.setEnd(vectors.get(i + 1).getStart());
//                vectorsStored.add(vectorL);
//            }

        }

        vectorsStored.add(vectors.get(vectors.size() - 1));

//
//
//        List<ProVector> betaSheets = new ArrayList<ProVector>();
//        for (ProVector proVector: vectors){
//            if(proVector.getType().equals("B")){
//                betaSheets.add(proVector);
//            }
//        }
//
//        ProVector seed = betaSheets.get(0);


//        List<ProVector> vectors = new ArrayList<ProVector>();
//
//        int STEP = 4;
//        for (int i = 2; i < atoms.length - STEP; i = i + STEP) {
//            Atom a1 = atoms[i];
//            Atom a4 = atoms[i + STEP];
//
//            char c1 = secondaryStr.charAt(i);
//            char c4 = secondaryStr.charAt(i + STEP);
//            ProVector vector = new ProVector();
//            vector.setPosStart(i);
//            vector.setPosEnd(i + STEP);
//            vector.setStart(a1);
//            vector.setEnd(a4);
//            String type = "L";
//            if (c1 == 'H' || c4 == 'H')
//                type = "H";
//            if (c1 == 'B' || c4 == 'B')
//                type = "B";
//            vector.setType(type);
//
//            vectors.add(vector);
//
//
//
//        }


        // List<ProVector> ssVectors = new ArrayList<ProVector>();
        //SeedVector.scanSeed(vectors);
        //
        // for (ProVector v : vectors) {
        // // System.out.print(v.getType());
        // if (v.getType().equals("H") || v.getType().equals("B")) {
        // ssVectors.add(v);
        // // System.out.println(Calc.getDistance(a1, a2));
        // }
        //
        // }

        // pro2Vector.getRepeatsTReks(vectors, atoms);

//        pro2Vector.getRepeats(vectors, atoms, 5);
//
//        for (ProVector v : vectors) {
//
//            System.out.print(v.getType());
//        }
//
//        System.out.println();

        // List<org.crns.crbm.lib.repeatsobjs.Repeat> repeats = pro2Vector
        // .getRepeatsVectorBySeed(vectors, atoms);

        // output = true;
        if (output)
            pro2Vector.saveVectorToFile(vectorsStored, pdbCode, pdbChain);

        // System.out.println(patternBuilder.toString());
        //
        // SuffixTree suffixTree = new SuffixTree();
        // suffixTree.findPattern(patternBuilder.toString());
        //
        // List<RepeatScore> repeat_scores = new ArrayList<RepeatScore>();
        //
        // for (int winsize = 2; winsize <= 3; winsize++) {
        // repeat_scores.addAll(pro2Vector.getRepeats(vectors, winsize));
        // }
        // // sort and display
        // RepeatScore.sortBySD(repeat_scores);
        //
        // if (repeat_scores.size() > 0) {
        //
        // System.out.println("repeat is found!:");
        // System.out.println(repeat_scores.get(0));
        // } else
        // System.out.println("Sorry, No repeat!!!");

        // if (args.length > 3) {
        // index = Integer.parseInt(args[3]);
        // // System.out.println(index);
        // }

        // choose the lowest rmsd
        // and find the best index

        // pro2Vector.show(list, index, max_point, winsize);
        // pro2Vector.saveToFile(pdbCode, pdbChain, list, index);

    }

    public List<org.cnrs.crbm.lib.trsfinder.Repeat> getRepeatsTReks(
            List<ProVector> vectors, Atom[] atoms) throws StructureException {

        // print pattern

        StringBuilder patternBuilder = new StringBuilder();
        StringBuilder patternAngleBuiilder = new StringBuilder();
        StringBuilder patternView = new StringBuilder();
        StringBuilder strAbsolutePos = new StringBuilder();
        StringBuilder strRelativePos = new StringBuilder();
        int i = 1;
        for (ProVector v : vectors) {

            patternView.append((v.getType() + "\t"));
            strRelativePos.append(i + "\t");
            strAbsolutePos.append(v.getStart().getGroup().getResidueNumber()
                    + "\t");
            patternBuilder.append(v.getType());
            patternAngleBuiilder.append(v.getAngleType());
            i++;

        }

        if (debug) {
            System.out.println(patternBuilder.toString());
            System.out.println(patternAngleBuiilder.toString());
            System.out.println(patternView.toString());
            System.out.println(strRelativePos.toString());
            System.out.println(strAbsolutePos.toString());
        }

        TREksWapper TReks = new TREksWapper("no-name",
                patternBuilder.toString());

        List<org.cnrs.crbm.lib.trsfinder.Repeat> repeat_output = new ArrayList<org.cnrs.crbm.lib.trsfinder.Repeat>();

        List<Repeat> repeats = TReks.getRepeats();
        // System.out.print(pdbCode.toLowerCase() + pdbChain + " " + repeatLeng
        // + " : ");
        int score_sure = 0;
        if (repeats.size() > 0) {
            for (Repeat repeat : repeats) {

                // set up a repeat output here
                org.cnrs.crbm.lib.trsfinder.Repeat a_repeat = new org.cnrs.crbm.lib.trsfinder.Repeat();

                int resi_s = vectors.get(repeat.getBeginPosition()).getStart()
                        .getGroup().getResidueNumber().getSeqNum();
                int resi_e = vectors.get(repeat.getEndPosition()).getStart()
                        .getGroup().getResidueNumber().getSeqNum();
                if (debug) {
                    System.out.print("pos: " + resi_s + "-" + resi_e + " avgL:"
                            + (resi_e - resi_s) / repeat.getNumber() + "|");
                    System.out.println();
                }
                // we would check all the Copy in each of an repeat by comparing
                // all of them

                LinkedList<Copy> copies = repeat.getCopies();
                // int count = 0;
                boolean is_tandem_repeat = false;
                int last_index = -2;
                // 2 pairs vectors repeats => tandem repeats
                for (int k = 0; k < repeat.getNumber() - 1; k++) {
                    // repeat fragments
                    Copy copy1 = copies.get(k);
                    Copy copy2 = copies.get(k + 1);
                    List<ProVector> list1 = new ArrayList<ProVector>();
                    List<ProVector> list2 = new ArrayList<ProVector>();
                    for (int m = copy1.getBeginPosition(); m <= copy1
                            .getEndPosition(); m++) {
                        list1.add(vectors.get(m));
                    }

                    for (int m = copy2.getBeginPosition(); m <= copy2
                            .getEndPosition(); m++) {
                        list2.add(vectors.get(m));
                    }

                    int start1 = vectors.get(copy1.getBeginPosition())
                            .getPosStart();
                    int end1 = vectors.get(copy1.getEndPosition()).getPosEnd();
                    int start2 = vectors.get(copy2.getBeginPosition())
                            .getPosStart();
                    int end2 = vectors.get(copy2.getEndPosition()).getPosEnd();

                    Atom[] atomSet1 = Fragement.getFragementsofAtoms(atoms,
                            start1, end1);
                    Atom[] atomSet2 = Fragement.getFragementsofAtoms(atoms,
                            start2, end2);
                    // double rmsd = impose.superimposeVectors(list1, list2);
                    // double rmsd = 0;

                    if (impose.isTheSameStructure(atomSet1, atomSet2)) {
                        if (last_index + 1 == k)
                            is_tandem_repeat = true;
                        last_index = k;

                        int m = copy1.getBeginPosition();
                        int n = copy1.getEndPosition();
                        a_repeat.getRepeats().add(
                                new RepeatContent(vectors.get(m).getPosStart(),
                                        vectors.get(n).getPosEnd()));

                    }
                }

                // if found tamden

                if (debug) {
                    for (Copy copy : copies) {
                        int m = copy.getBeginPosition();
                        System.out.print(vectors.get(m).getStart().getGroup()
                                .getResidueNumber()
                                + " ");
                    }
                    a_repeat.sortByPosition();
                    repeat_output.add(a_repeat);
                    System.out.println();
                    System.out.println("***");
                    if (is_tandem_repeat)
                        System.out.println("tandem repeat found");
                }
                // if (is_tandem_repeat) {
                //
                //
                //
                // }

            }

            //

        }

        return repeat_output;
        // } else
        // System.out.println("Sorry, No repeat!!!");

    }

    public List<org.cnrs.crbm.lib.trsfinder.Repeat> getRepeatsVectorBySeed(
            List<ProVector> vectors, Atom[] atoms) throws Exception {
        Superimposer superimpose = new Superimposer();
        List<org.cnrs.crbm.lib.trsfinder.Repeat> repeats = new ArrayList<org.cnrs.crbm.lib.trsfinder.Repeat>();
        List<ProVector> ssVectors = ProVector.toSecondaryVector(vectors);
        // SeedVector.scanSeed(vectors);

        List<List<Integer>> seeds = SeedVector.scanSeed(ssVectors);
        for (List<Integer> seed : seeds) {
            int p0 = ssVectors.get(seed.get(0)).getPosStart();
            int p1 = ssVectors.get(seed.get(1)).getPosStart();
            int p2 = ssVectors.get(seed.get(2)).getPosStart();

            Atom[] as1 = Fragement.getFragementsofAtoms(atoms, p0, p1);
            Atom[] as2 = Fragement.getFragementsofAtoms(atoms, p1, p2);

            // get patterns

            String pattern_as1 = "";
            for (int i = seed.get(0); i < seed.get(1); i++) {
                pattern_as1 += ssVectors.get(i).getType();
            }
            String pattern_as2 = "";
            for (int i = seed.get(1); i < seed.get(2); i++) {
                pattern_as2 += ssVectors.get(i).getType();
            }

            // check both structure and pattern
            if (superimpose.isTheSameStructure(as1, as2)
                    && pattern_as1.equals(pattern_as2)) {
                org.cnrs.crbm.lib.trsfinder.Repeat aRepeat = new org.cnrs.crbm.lib.trsfinder.Repeat();
                aRepeat.getRepeats().add(new RepeatContent(p0, p1 - 1));
                aRepeat.getRepeats().add(new RepeatContent(p1, p2));
                // System.out.println(p0 + "-" + p1 + "-" + p2);
                repeats.add(aRepeat);

            }
        }

        return repeats;

    }

    public List<RepeatScore> getRepeats(List<ProVector> vectors, Atom[] atoms,
                                        int winsize) throws StructureException {

        // if (winsize == 2)
        // System.out.println("winsize:" + winsize);
        // vectors = ProVector.toSecondaryVector(vectors);

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
        List<RepeatScore> repeat_scores = new ArrayList<RepeatScore>();
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
            int total_point_below_thresholds = 0;

            // System.out.println(">" + j);

            List<RepeatStore> store_vectors = new ArrayList<RepeatStore>();
            store_vectors.add(new RepeatStore(sub_list_ref, 0));

            int number_tandem = 0;
            int store_current_posstion = j;

            as1 = Fragement.getFragementsofAtoms(atoms, sub_list_ref.get(0)
                    .getPosStart(), sub_list_ref.get(sub_list_ref.size() - 1)
                    .getPosEnd());

            int start1 = sub_list_ref.get(0).getPosStart();
            int end1 = sub_list_ref.get(sub_list_ref.size() - 1).getPosEnd();
            // scan to the right

            for (int i = j + winsize; i < list.size(); ) {
                List<ProVector> sub_list = list.get(i);

                as2 = Fragement.getFragementsofAtoms(atoms, sub_list.get(0)
                        .getPosStart(), sub_list.get(sub_list.size() - 1)
                        .getPosEnd());

                int start2 = sub_list.get(0).getPosStart();

                if (as1.length > 60
                        && (start2 - end1 + 1) > (double) 0.5 * as1.length)
                    break;// no check anymore

                else if (as1.length < 20
                        && (start2 - end1 + 1) > (double) 0.5 * as1.length)
                    break;

                // check pattern
                String pattern_sub = this.getPatterOfVectors(sub_list);

                if (pattern_sub.equals(pattern_ref)
                        && impose.isTheSameStructure(as1, as2)) {
                    // System.out.println(pattern_ref);
                    total_point_below_thresholds++;
                    store_vectors.add(new RepeatStore(sub_list, 1));
                    if (store_current_posstion + winsize == i) {
                        number_tandem++;
                    }
                    store_current_posstion = i;

                    end1 = sub_list.get(sub_list.size() - 1).getPosEnd();
                    i = i + winsize;

                } else {
                    i = i + 1;
                }

            }

            // scan to the left;

            // store_current_posstion = j;
            // for (int i = j - winsize; i > 0;) {
            // List<ProVector> sub_list = list.get(i);
            // // check pattern
            // String pattern_sub = this.getPatterOfVectors(sub_list);
            //
            // Atom[] as2 = Fragement.getFragementsofAtoms(atoms, sub_list
            // .get(0).getPosStart(), sub_list
            // .get(sub_list.size() - 1).getPosEnd());
            //
            // int end2 = sub_list.get(sub_list.size() - 1).getPosEnd();
            //
            // if ((start1 - end2 + 1) > (double) 0.6 * as1.length)
            // break;// no check anymore
            //
            // // System.out.println(rmsd);
            // if (pattern_sub.equals(pattern_ref)
            // && impose.isTheSameStructure(as1, as2)) {
            // total_point_below_thresholds++;
            // store_vectors.add(new RepeatStore(sub_list, 1));
            //
            // if (store_current_posstion - winsize == i) {
            // number_tandem++;
            // }
            // store_current_posstion = i;
            //
            // i = i - winsize;
            //
            // } else {
            // i = i - 1;
            // }
            // }
            //
            // // store to find the best
            // if (total_point_below_thresholds > max_point) {
            // max_point = total_point_below_thresholds;
            // index = j;
            // }

            if (store_vectors.size() >= 2) {
                // repeats found
                // print pattern and so on
                // posible repeats
                double number_posible_repeats = vectors.size() / winsize;
                // calculate score;
                double score_tandem = 1 - (double) number_tandem
                        / store_vectors.size();
                double score_repeat_size = 1 - (double) store_vectors.size()
                        / number_posible_repeats;

                repeat_scores.add(new RepeatScore(store_vectors, score_tandem,
                        score_repeat_size));

                lastTRs = j;

            }
        }
        return repeat_scores;

    }

    public String getPatterOfVectors(List<ProVector> sub_list) {
        String pattern_sub = "";
        for (ProVector v : sub_list) {
            pattern_sub += v.getType();
        }
        return pattern_sub;
    }

    private void saveToFile(String pdbCode, String pdbChain,
                            List<List<ProVector>> list, int index) {
        String matlab_dir = "C:\\Users\\CRBM\\Desktop\\Matlab\\data\\";
        String outputfile = matlab_dir + pdbCode + "_" + pdbChain + ".vect";

        try {
            FileWriter fstream = new FileWriter(outputfile);
            BufferedWriter out = new BufferedWriter(fstream);

            List<ProVector> sub_list_ref = list.get(index);

            for (int j = 0; j < list.size(); j++) {
                List<ProVector> sub_list = list.get(j);
                // System.out.println(j);

                double rmsd = impose.superimposeVectors(sub_list_ref, sub_list);
                if (rmsd < RMSD_THRESHOLD) {

                    System.out.print("<"
                            + sub_list.get(0).getOriginalStart().getGroup()
                            .getResidueNumber().getSeqNum() + ">");

                }
                out.write(sub_list.get(0).getOriginalStart().getGroup()
                        .getResidueNumber().getSeqNum()
                        + "	" + NumberFormatUtils.format(rmsd) + "\n");

            }

            out.close();
            fstream.close();
        } catch (Exception e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
    }

    public void saveVectorToFile(List<ProVector> vectors, String pdbCode,
                                 String pdbChain) {
        String out_dir = "C:\\Users\\pdoviet\\Desktop\\";
        String outputfile = out_dir + pdbCode.toLowerCase() + pdbChain
                + "_vect.pdb";

        String outfile_pml = out_dir + pdbCode.toLowerCase() + pdbChain
                + "_show.pml";

        try {
            FileWriter fstream = new FileWriter(outputfile);
            BufferedWriter out = new BufferedWriter(fstream);

            FileWriter fstream_pml = new FileWriter(outfile_pml);
            BufferedWriter out_pml = new BufferedWriter(fstream_pml);

            // Structure newstruc = new StructureImpl();
            // Chain c1 = new ChainImpl();
            // c1.setChainID(pdbChain);
            out_pml.write("set dash_gap, 0.0\n");

            int i = 0;
            for (ProVector v : vectors) {

                // System.out.println(v.getStart() + ":"
                // + v.getOriginalStart());
                out.write(v.getStart().toPDB() + "\n");
                out.write(v.getEnd().toPDB() + "\n");
                int p1 = v.getStart().getGroup().getResidueNumber().getSeqNum();
                int p2 = v.getEnd().getGroup().getResidueNumber().getSeqNum();
                out_pml.write("distance d" + i + ", ////" + p1 + ", ////" + p2
                        + " \n");
                out_pml.write("hide labels, d" + i + "\n");
                if (v.getType().equals("H"))
                    out_pml.write("color red, d" + i + "\n");
                else if (v.getType().equals("B"))
                    out_pml.write("color yellow, d" + i + "\n");
                else
                    out_pml.write("color gray, d" + i + "\n");
                i++;

                // c1.addGroup(v.getOriginalStart().getGroup());
                // c1.addGroup(v.getOriginalEnd().getGroup());

            }
            // newstruc.addChain(c1);
            // out.write(newstruc.toPDB());
            out_pml.close();
            out.close();
            fstream.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private void show(List<List<Atom>> list, int index, int max_point,
                      int winsize) throws StructureException {

        List<Integer> lowest_points = new ArrayList<Integer>();
        System.out.println("point seed    : " + index);
        System.out.println("number points : " + max_point);
        List<Atom> sub_list_ref = list.get(index);

        for (int j = 0; j < list.size(); j++) {
            List<Atom> sub_list = list.get(j);
            // System.out.println(j);

            double rmsd = impose.superimposeAtoms(sub_list_ref, sub_list);
            if (rmsd < RMSD_THRESHOLD) {
                lowest_points.add(j);
                // System.out.print("<"
                // + sub_list.get(0).getGroup().getResidueNumber()
                // .getSeqNum() + ">");

            }

            System.out.print("<" + j + ":"
                    + sub_list.get(0).getGroup().getResidueNumber().getSeqNum()
                    + ">");
            System.out.print(NumberFormatUtils.format(rmsd) + " ");
        }

        System.out.println("\nlowest points:");
        for (int p : lowest_points) {
            System.out.print(p + " ");
        }
        System.out.println("\ndist each points:");
        for (int i = 0; i < lowest_points.size() - 1; i++) {
            System.out.print(lowest_points.get(i + 1) - lowest_points.get(i)
                    + " ");
        }

    }

}
