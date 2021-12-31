package org.cnrs.crbm.lib.repeats;

import org.cnrs.crbm.lib.analysis.ConforAlphaAnalysis;
import org.cnrs.crbm.lib.analysis.PowerRepeatVsNoRepeat;
import org.cnrs.crbm.lib.analysis.TMScoreAnalysis;
import org.cnrs.crbm.lib.analysis.VectorThreshold;
import org.cnrs.crbm.lib.io.DataIO;
import org.cnrs.crbm.lib.repeats.module.HetatmModule;
import org.cnrs.crbm.lib.repeats.shortTRs.LocationMatch;
import org.cnrs.crbm.lib.repeats.shortTRs.ShortTRsDetector;
import org.cnrs.crbm.lib.trsfinder.Features;
import org.cnrs.crbm.lib.utils.NumberFormatUtils;

import java.util.List;

/**
 * Created by pdoviet on 5/29/2015.
 */
public class ProduceScore {


    public CombineScore getScore(String pdbCode, String pdbChain) throws Exception {
        CombineScore combineScore = new CombineScore();
        RepeatFinder repeatFinder = new RepeatFinder(pdbCode, pdbChain);
        Features features = repeatFinder.getFeatures();

        PowerRepeatVsNoRepeat powerRepeatVsNoRepeat = new PowerRepeatVsNoRepeat();

        //Set SG-Score
        double signalScore = powerRepeatVsNoRepeat.getSignalScore(features);
        combineScore.setSigScore(signalScore);
        // Set V-Score
        VectorThreshold vectorThreshold = new VectorThreshold();
        double vScore = vectorThreshold.getVScore(features);
        combineScore.setvScore(vScore);
        // Set TM-Score and CE-Score
        TMScoreAnalysis evaluation = new TMScoreAnalysis(features);
        evaluation.findRepeat();
        combineScore.setTmScore(evaluation.getTmScore());
        combineScore.setCeScore(evaluation.getCeScore());
        // Set CA-Score
        ConforAlphaAnalysis conforAlphaAnalysis = new ConforAlphaAnalysis();
        double caScore = Math.max(conforAlphaAnalysis.processByTREKS(features), conforAlphaAnalysis.processPdbbyTRUST(features));
        combineScore.setPsimScore(caScore);
        // Set H-Score
        HetatmModule haHetatmModule = new HetatmModule();
        combineScore.setlScore(haHetatmModule.getHETAMScore(features));
        // Set M-Score
        double mScore = 0;
        ShortTRsDetector detector = new ShortTRsDetector();
        List<LocationMatch> listFromDetector = detector.getRepeatLocations(features.getStrSeqAp());
        if (listFromDetector.size() > 0)
            mScore = 1;
        combineScore.setCaScore(mScore);
        return combineScore;

    }


    public static void main(String[] args) {


        ProduceScore produceScore = new ProduceScore();
        produceScore.training();

    }

    private void training() {

        List<String> rows = DataIO.readLines("data/tapo/trainsetPdb.in");
        for (String row : rows) {

            String pdbCode = row.substring(0, 4);
            String pdbChain = row.substring(5, 6);
            RepeatFinder repeatFinder = new RepeatFinder(pdbCode, pdbChain);
            Features features = repeatFinder.getFeatures();
            //Set SG-Score
            PowerRepeatVsNoRepeat powerRepeatVsNoRepeat = new PowerRepeatVsNoRepeat();
            double signalScore = 0.0;
            try {
                signalScore = powerRepeatVsNoRepeat.getSignalScore(features);
            } catch (Exception e) {
                e.printStackTrace();
            }
            // Set H-Score
            HetatmModule haHetatmModule = new HetatmModule();
            double hScore = haHetatmModule.getHETAMScore(features);
            // Set M-Score
            double mScore = 0;
            ShortTRsDetector detector = new ShortTRsDetector();
            List<LocationMatch> listFromDetector = detector.getRepeatLocations(features.getStrSeqAp());
            if (listFromDetector.size() > 0)
                mScore = 1;
            System.out.println(pdbCode + "_" + pdbChain + "\t" + mScore + "\t" + NumberFormatUtils.format(hScore) + "\t" + NumberFormatUtils.format(signalScore));

        }

    }


}


