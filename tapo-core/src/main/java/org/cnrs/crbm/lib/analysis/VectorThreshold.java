package org.cnrs.crbm.lib.analysis;


import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.cnrs.crbm.lib.io.DataIO;
import org.cnrs.crbm.lib.io.ProteinCSVReader;
import org.cnrs.crbm.lib.io.Row;
import org.cnrs.crbm.lib.io.RowRepeatDB;
import org.cnrs.crbm.lib.multalign.MutilAlign;
import org.cnrs.crbm.lib.repeats.Fragement;
import org.cnrs.crbm.lib.repeats.module.ProVector;
import org.cnrs.crbm.lib.repeats.RepeatFinder;
import org.cnrs.crbm.lib.repeats.Superimposer;
import org.cnrs.crbm.lib.repeats.module.VectorModule;
import org.cnrs.crbm.lib.trsfinder.Features;
import org.cnrs.crbm.lib.utils.NumberFormatUtils;
import org.cnrs.crbm.lib.utils.ProgressBar;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by pdoviet on 1/16/2015.
 */
public class VectorThreshold {
    final static int MIN_ELEMENTS = 2;
    final static int MAX_ELEMENTS = 10;

    public static void main(String[] args) throws Exception {

        new VectorThreshold().training();
    }

    Superimposer superimposer = new Superimposer();

    public double getVScore(Features features) throws Exception {

        Atom[] atoms = features.getAtoms();

        VectorModule vectorModule = new VectorModule();
        MutilAlign mutilAlign = new MutilAlign();

//            List<ProVector> vectors = ProVector.toSecondaryVector(repeatFinder.getFeatures()
//                    .getVectors());
        List<ProVector> vectors = ProVector.toSecondaryVector(vectorModule.extractVectors(atoms, features.getStrSS()));
        double storedVScore = 0.0;
        double storedTMScore = 0.0;
        double maxScore = 0.0;

        for (int winsize = MIN_ELEMENTS; winsize <= MAX_ELEMENTS; winsize++) {

            List<List<ProVector>> list = new ArrayList<List<ProVector>>();
            for (int i = 0; i < vectors.size() - winsize + 1; i = i + 1) {
                List<ProVector> sub_list = new ArrayList<ProVector>();
                for (int j = 0; j < winsize; j++) {
                    ProVector v = vectors.get(i + j);
                    sub_list.add(v);
                }
                list.add(sub_list);
            }

            for (int j = 0; j < list.size() - winsize; j++) {
                List<ProVector> sub_list_ref = list.get(j);
                int i = j + winsize;
                List<ProVector> sub_list = list.get(i);

                Atom[] as1 = Fragement.getFragementsofAtoms(atoms, sub_list_ref.get(0)
                        .getPosStart(), sub_list_ref.get(sub_list_ref.size() - 1)
                        .getPosEnd());

                Atom[] as2 = Fragement.getFragementsofAtoms(atoms, sub_list.get(0)
                        .getPosStart(), sub_list.get(sub_list.size() - 1)
                        .getPosEnd());


                // check pattern
                String pattern_sub_ref = this.getPatterOfVectors(sub_list_ref);

                String pattern_sub = this.getPatterOfVectors(sub_list);
                if (pattern_sub.equals(pattern_sub_ref)) {

                    AFPChain afpChain = mutilAlign.pairAlign(as1, as2);
                    double tmScore = afpChain.getTMScore();
                    double vscore = superimposer.compareListVector(sub_list_ref,
                            sub_list);

                    maxScore = Math.max(maxScore, (vscore + tmScore) / 2);
                    if (vscore > storedVScore) {
                        storedVScore = vscore;
                        storedTMScore = tmScore;
                    }

                }
                //if (pattern_sub.equals(pattern_sub_ref))
//                        maxScore = Math.max(score, maxScore);
                //System.out.println(maxScore);
//                    writer.write(entry + "\t" + NumberFormatUtils.format(score) + "\t" + strClass + "\n");

            }

        }

        return maxScore;

    }

    public void training() throws Exception {
        Superimposer superimpose = new Superimposer();
        //List<String> rows = DataIO.readLines("data/tapo/trainsetPdb.in");
        List<String> rows = DataIO.readLines("data/benchmarkset/tmp.in");
        PrintWriter writer = new PrintWriter("data/vector_thres.train");
        writer.write("entry\tvscore\ttmscore\tclass\n");
        ProgressBar bar = new ProgressBar();
        int process = 1;
        int sizeOfProcess = rows.size();
        for (String row : rows) {

            String pdbCode = row.substring(0, 4);
            String pdbChain = row.substring(5, 6);
//            String strClass = row.split("\t")[1];
//            if (strClass.equals("1"))
//                strClass = "1";
//            else
//                strClass = "0";

            String strClass = "XX";

            String entry = pdbCode + "_" + pdbChain;


            bar.update(process, sizeOfProcess);
            process++;
            // do something here

            RepeatFinder repeatFinder = new RepeatFinder(pdbCode, pdbChain);
            Atom[] atoms = repeatFinder.getAtoms();

            VectorModule vectorModule = new VectorModule();
            MutilAlign mutilAlign = new MutilAlign();

//            List<ProVector> vectors = ProVector.toSecondaryVector(repeatFinder.getFeatures()
//                    .getVectors());
            List<ProVector> vectors = ProVector.toSecondaryVector(vectorModule.extractVectors(pdbCode, pdbChain));
            double storedVScore = 0.0;
            double storedTMScore = 0.0;
            double maxScore = 0.0;

            for (int winsize = MIN_ELEMENTS; winsize <= MAX_ELEMENTS; winsize++) {

                List<List<ProVector>> list = new ArrayList<List<ProVector>>();
                for (int i = 0; i < vectors.size() - winsize + 1; i = i + 1) {
                    List<ProVector> sub_list = new ArrayList<ProVector>();
                    for (int j = 0; j < winsize; j++) {
                        ProVector v = vectors.get(i + j);
                        sub_list.add(v);
                    }
                    list.add(sub_list);
                }

                for (int j = 0; j < list.size() - winsize; j++) {
                    List<ProVector> sub_list_ref = list.get(j);
                    int i = j + winsize;
                    List<ProVector> sub_list = list.get(i);

                    Atom[] as1 = Fragement.getFragementsofAtoms(atoms, sub_list_ref.get(0)
                            .getPosStart(), sub_list_ref.get(sub_list_ref.size() - 1)
                            .getPosEnd());

                    Atom[] as2 = Fragement.getFragementsofAtoms(atoms, sub_list.get(0)
                            .getPosStart(), sub_list.get(sub_list.size() - 1)
                            .getPosEnd());


                    // check pattern
                    String pattern_sub_ref = this.getPatterOfVectors(sub_list_ref);

                    String pattern_sub = this.getPatterOfVectors(sub_list);
                    if (pattern_sub.equals(pattern_sub_ref)) {

                        AFPChain afpChain = mutilAlign.pairAlign(as1, as2);
                        double tmScore = afpChain.getTMScore();
                        double vscore = superimpose.compareListVector(sub_list_ref,
                                sub_list);

                        maxScore = Math.max(maxScore, (vscore + tmScore) / 2);
                        if (vscore > storedVScore) {
                            storedVScore = vscore;
                            storedTMScore = tmScore;
                        }

                    }
                    //if (pattern_sub.equals(pattern_sub_ref))
//                        maxScore = Math.max(score, maxScore);
                    //System.out.println(maxScore);
//                    writer.write(entry + "\t" + NumberFormatUtils.format(score) + "\t" + strClass + "\n");

                }

            }

            //writer.write(entry + "\t" + NumberFormatUtils.format(storedVScore) +"\t"+  NumberFormatUtils.format(storedTMScore) + "\t" + strClass + "\n");
            writer.write(entry + "\t" + NumberFormatUtils.format(maxScore) + "\t" + strClass + "\n");

        }

        writer.close();
    }

    public void run() throws Exception {
        ProteinCSVReader csvReader = new ProteinCSVReader();
        List<Row> rowsNonTRs = csvReader.getData("data/RepeatDatalastest.csv");
        List<RowRepeatDB> rowsRB = csvReader.getRepeatDB("data/RDB-dataset.tab");

        List<Row> rowsTRs = new ArrayList<Row>();
        //convert to rows
        for (RowRepeatDB r : rowsRB) {

            if (r.getAnnlevel().equals("Detailed")) {
                Row row = new Row();
                row.setProtein(r.getPdbCode() + "_" + r.getPdbChain());
                rowsTRs.add(row);
            }
        }


        PrintWriter writer = new PrintWriter("data/vector_thres.train");
        writer.write("sScore\tclass\n");
        ProgressBar bar = new ProgressBar();

        List<String> ignoreCase = new ArrayList<String>();
        ignoreCase.add("1yo8");
        ignoreCase.add("2qcs");
        ignoreCase.add("4f5v");
        ignoreCase.add("3m7h");
        int sizeOfProcess = rowsNonTRs.size() + rowsTRs.size();
        // bar.update(0, sizeOfProcess);
        int process = 1;
        Superimposer superimpose = new Superimposer();

        for (Row row : rowsNonTRs) {
            bar.update(process, sizeOfProcess);
            process++;
            // do something here
            if (ignoreCase.contains(row.getPdbCode()))
                continue;
            RepeatFinder repeatFinder = new RepeatFinder(row.getPdbCode(), row.getPdbChain());

            List<ProVector> vectors = ProVector.toSecondaryVector(repeatFinder.getFeatures()
                    .getVectors());

            double maxScore = Double.MIN_VALUE;
            for (int winsize = MIN_ELEMENTS; winsize <= MAX_ELEMENTS; winsize++) {

                List<List<ProVector>> list = new ArrayList<List<ProVector>>();
                for (int i = 0; i < vectors.size() - winsize + 1; i++) {
                    List<ProVector> sub_list = new ArrayList<ProVector>();
                    for (int j = 0; j < winsize; j++) {
                        ProVector v = vectors.get(i + j);
                        sub_list.add(v);
                    }
                    list.add(sub_list);
                }

                for (int j = 0; j < list.size() - winsize; j++) {
                    List<ProVector> sub_list_ref = list.get(j);
                    int i = j + winsize;
                    List<ProVector> sub_list = list.get(i);
                    double score = superimpose.compareListVector(sub_list_ref,
                            sub_list);
                    maxScore = Math.max(score, maxScore);
                    //System.out.println(maxScore);


                }

            }

            writer.write(NumberFormatUtils.format(maxScore) + "\t" + "No-TRs\n");


        }

        for (Row row : rowsTRs) {
            bar.update(process, sizeOfProcess);
            process++;
            // do something here
            if (ignoreCase.contains(row.getPdbCode()))
                continue;
            RepeatFinder repeatFinder = new RepeatFinder(row.getPdbCode(), row.getPdbChain());

            List<ProVector> vectors = ProVector.toSecondaryVector(repeatFinder.getFeatures()
                    .getVectors());

            double maxScore = Double.MIN_VALUE;
            for (int winsize = MIN_ELEMENTS; winsize <= MAX_ELEMENTS; winsize++) {

                List<List<ProVector>> list = new ArrayList<List<ProVector>>();
                for (int i = 0; i < vectors.size() - winsize + 1; i++) {
                    List<ProVector> sub_list = new ArrayList<ProVector>();
                    for (int j = 0; j < winsize; j++) {
                        ProVector v = vectors.get(i + j);
                        sub_list.add(v);
                    }
                    list.add(sub_list);
                }

                for (int j = 0; j < list.size() - winsize; j++) {
                    List<ProVector> sub_list_ref = list.get(j);
                    int i = j + winsize;
                    List<ProVector> sub_list = list.get(i);
                    double score = superimpose.compareListVector(sub_list_ref,
                            sub_list);
                    maxScore = Math.max(score, maxScore);
                    //System.out.println(maxScore);


                }

            }

            writer.write(NumberFormatUtils.format(maxScore) + "\t" + "TRs\n");


        }


        writer.close();

    }

    private String getPatterOfVectors(List<ProVector> sub_list) {
        String pattern_sub = "";
        for (ProVector v : sub_list) {
            pattern_sub += v.getType();
        }
        return pattern_sub;
    }

}
