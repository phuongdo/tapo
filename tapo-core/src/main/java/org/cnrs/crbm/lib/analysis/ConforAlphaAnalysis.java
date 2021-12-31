package org.cnrs.crbm.lib.analysis;

import org.cnrs.crbm.lib.io.DataIO;
import org.cnrs.crbm.lib.multalign.PSim;
import org.cnrs.crbm.lib.multalign.TSim;
import org.cnrs.crbm.lib.repeats.CombineScore;
import org.cnrs.crbm.lib.repeats.RepeatFinder;
import org.cnrs.crbm.lib.sadb.Conformation;
import org.cnrs.crbm.lib.sadb.Seq3d;
import org.cnrs.crbm.lib.seqalign.RepeatSeq;
import org.cnrs.crbm.lib.seqalign.RepeatUnitSeq;
import org.cnrs.crbm.lib.seqalign.TREksWapper;
import org.cnrs.crbm.lib.seqalign.Trust;
import org.cnrs.crbm.lib.trsfinder.*;
import org.cnrs.crbm.lib.utils.ProgressBar;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

/**
 * This class aim to training the optimal score for the conformational alphabets method.
 * <p/>
 * Created by pdoviet on 5/21/2015.
 */
public class ConforAlphaAnalysis {


    public static void main(String[] args) throws Exception {
        //new ConforAlphaAnalysis().fetchTrainingSet();
        //new ConforAlphaAnalysis().processPdb("4dh2", "B");
//        RepeatFinder repeatFinder  = new RepeatFinder("1hfe","L");
//        System.out.println(new ConforAlphaAnalysis().processPdbbyTRUST(repeatFinder.getFeatures()));
        new ConforAlphaAnalysis().training(args[0]);
    }

    public void training(String id) throws Exception {

        List<String> rows = DataIO.readLines("F:/Storage/tmp/job_array/input."+id);
        PrintWriter writer = new PrintWriter("data/tapo/ca_score.out."+id);
        //writer.write("protein\tcaScore\n");

        ProgressBar bar = new ProgressBar();
        int process = 1;
        int sizeOfProcess = rows.size();
        for (String row : rows) {

            bar.update(process, sizeOfProcess);
            process++;
            String pdbCode = row.substring(0, 4);
            String pdbChain = row.substring(5, 6);
            RepeatFinder repeatFinder = new RepeatFinder(pdbCode,pdbChain);
            List<Finder> finders = new ArrayList<Finder>();
            // conformational alphabets
            Features features = repeatFinder.getFeatures();
            Finder trustfinder = new TrustFinder(features);
            Finder treksfinder = new TReksFinder(features);
            Finder treksSADBFinder = new TReksSADBFinder(features);
            Finder treksOp1Finder = new TReksOp1Finder(features);
            Finder treksOp2Finder = new TReksOp2Finder(features);
            Finder trustCombinedFinder = new TrustCombinedFinder(features);

            finders.add(trustfinder);
            finders.add(treksfinder);
            finders.add(treksSADBFinder);
            finders.add(treksOp1Finder);
            finders.add(treksOp2Finder);
            finders.add(trustCombinedFinder);
            /**
             * all finders are run concurrently
             */
           // int nrOfProcessors = 6;
           // ExecutorService eservice = Executors
           //         .newFixedThreadPool(nrOfProcessors);
           // List<Future> futuresList = new ArrayList<Future>();
            for (Finder finder : finders) {
                //futuresList.add(eservice.submit(new JobFinder(finder)));\
                finder.start();
            }
            //eservice.shutdown();
            //List<Finder> findersOut = new ArrayList<Finder>();
            //for (Future future : futuresList) {
            //    try {
            //        findersOut.add((Finder) future.get());
                    //Finder finder = (Finder) future.get();
            //    } catch (InterruptedException e) {
            //    } catch (ExecutionException e) {
            //        e.printStackTrace();
           //     }
           // }
            CombineScore combineScore = new CombineScore();
            for (Finder finder : finders) {
                // get score here
                CombineScore cScore = finder.getCombineScore();
                combineScore.setTmScore(Math.max(cScore.getTmScore(), combineScore.getTmScore()));
                combineScore.setvScore(Math.max(cScore.getvScore(), combineScore.getvScore()));
                combineScore.setlScore(Math.max(cScore.getlScore(), combineScore.getlScore()));
                combineScore.setCeScore(Math.max(cScore.getCeScore(), combineScore.getCeScore()));
                combineScore.setPsimScore(Math.max(cScore.getPsimScore(), combineScore.getPsimScore()));
                combineScore.setSigScore(Math.max(cScore.getSigScore(), combineScore.getSigScore()));
            }

            writer.write(pdbCode + "_" + pdbChain + "\t" + combineScore.getPsimScore() + "\n");


        }


        writer.close();
    }

    public double processPdbbyTRUST(Features features) {
        //RepeatFinder finder = new RepeatFinder(pdbCode, pdbChain);
        //Features features = finder.getFeatures();
        Seq3d seq3D = new Seq3d(features.getStrSeqAp());
        //String pattern = seq3D.getSeq();
        //String pattern = features.getStrSeqAp();
        String pattern = features.getStrSeqSADB();
        List<Conformation> confors = seq3D.getConfors();
//        int gapo = -4;
//        int gape = -2;
//        Trust trust = new Trust("conf/BLOSUM_COMBINED", "EKLIBVPSDCGHAFXYOWMN",
//                features.getPdbCode(), pattern, gapo, gape);
        int gapo = -8;
        int gape = -2;
//        pattern = features.getStrSeqAp();
//        System.out.println(pattern);

        Trust trust = new Trust("conf/BLOSUM62", "ABCDEFGHIKLMNPQRSTVWXYZO",
                features.getPdbCode(), pattern, gapo, gape);
//        Trust trust = new Trust("conf/BLOSUM_COMBINED", "EKLIBVPSDCGHAFXYOWMN",
//                features.getPdbCode(), pattern, gapo, gape);
//        Collection repeats_trust = trust.getRepeatsOutput();
        //List<String> msa = trust.getMsa();
        try {
            TSim tSim = new TSim();

            double maxScore = 0.0;
            List<RepeatSeq> repeats = trust.getRepeats();
            for (RepeatSeq repeat : repeats) {
                // each of repeats

                List<String> msa = new ArrayList<String>();

                for (RepeatUnitSeq ru : repeat.getUnits()) {
                    msa.add(ru.getStrAlign());
                }
                maxScore = Math.max(maxScore, tSim.getTsimScore(msa));

            }

            return maxScore;
        } catch (Exception ex) {
            ex.printStackTrace();
            // do nothing here
        }

        return 0.0;
    }

    public double processByTREKS(Features features) {
        //RepeatFinder finder = new RepeatFinder(pdbCode, pdbChain);
        Seq3d seq3D = new Seq3d(features.getStrSeqAp());
        int gapo = -5;
        int gape = -1;
        String pattern = seq3D.getSeq();
        TREksWapper TReks = new TREksWapper("", pattern);
        // get treks repeats

        try {
            TSim tSim = new TSim();

            double maxScore = 0.0;
            List<RepeatSeq> repeats = TReks.getRepeatSeqs();
            for (RepeatSeq repeat : repeats) {
                // each of repeats
                if (repeat.getUnits().get(0).getStrAlign().length() > 3) {
                    List<String> msa = new ArrayList<String>();
                    for (RepeatUnitSeq ru : repeat.getUnits()) {
                        msa.add(ru.getStrAlign());
                    }

                    PSim psim = new PSim(msa);

                    //maxScore = Math.max(maxScore, tSim.getTsimScore(msa));
                    maxScore = Math.max(maxScore, psim.getSimilarity());
                }
            }

            return maxScore;
        } catch (Exception ex) {
            ex.printStackTrace();
            // do nothing here
        }

        return 0.0;

    }

    public List<String> fetchTrainingSet() throws FileNotFoundException {

        List<String> rows = DataIO.readLines("data/tapo/trainsetPdb.in");
        return rows;
    }

}
