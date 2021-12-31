package org.cnrs.crbm.lib.predLen;


import org.biojava.nbio.structure.Atom;
import org.cnrs.crbm.lib.io.TaPoFastaFormat;
import org.cnrs.crbm.lib.raphael.Raphael;
import org.cnrs.crbm.lib.repeats.Fragement;
import org.cnrs.crbm.lib.repeats.RepeatFinder;
import org.cnrs.crbm.lib.repeats.RepeatRegion;
import org.cnrs.crbm.lib.repeats.module.VectorShape;
import org.cnrs.crbm.lib.trsfinder.Region;
import org.cnrs.crbm.lib.trsfinder.Repeat;
import org.cnrs.crbm.lib.trsfinder.RepeatContent;
import org.cnrs.crbm.lib.utils.PdbTools;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by pdoviet on 8/28/2015.
 */
public class PredLenMod {

    DemoLength demoLength = new DemoLength();

    public List<RepeatRegion> predictByTAPOFormat(TaPoFastaFormat taPoFastaFormat) throws Exception {

        List<Repeat> lstRep = taPoFastaFormat.getRepeats();
        List<RepeatRegion> predRegions = new ArrayList<RepeatRegion>();
        String pdb = taPoFastaFormat.getPdbID();
        String pdbCode = pdb.substring(0, 4);
        String pdbChain = pdb.substring(5, 6);
        RepeatFinder repeatFinder = new RepeatFinder(pdbCode, pdbChain);
        Atom[] atoms = repeatFinder.getAtoms();
        String strSS = repeatFinder.getStrSS();
        double[] consensus = new double[atoms.length];
        double[] consensusLen = new double[atoms.length];
        double[] consensusSize = new double[atoms.length];
        for (int i = 0; i < consensusLen.length; i++) {
            consensusLen[i] = 1000;

        }
        for (int i = 0; i < consensusSize.length; i++) {

            int size = 0;
            for (Repeat repeat : lstRep) {
                if (repeat.getStart() < i && i < repeat.getEnd()) {
                    size++;
                }
            }
            consensusSize[i] = size;
        }
        Map<Integer, List<Double>> map = new HashMap<Integer, List<Double>>();
        for (Repeat repeat : lstRep) {

            for (int i = repeat.getStart(); i < repeat.getEnd(); i++) {
                consensus[i] += repeat.getScore();
                //consensusLen[i] = Math.min(consensusLen[i], repeat.getAvgLength());
            }

            for (RepeatContent rc : repeat.getRepeats()) {
                //RepeatContent rc = repeat.getRepeats().get(0);
                String pattern = VectorShape.getSSPattern(strSS
                        .substring(rc.getStart(), rc.getEnd()));

                for (int i = rc.getStart(); i < rc.getEnd(); i++) {

                    if (pattern.length() >= 2)
                        consensusLen[i] = Math.min(consensusLen[i], rc.size());
                }

                if (pattern.length() >= 2) {

                    if (map.containsKey(pattern.length())) {
                        List<Double> lstD = map.get(pattern.length());
                        lstD.add(repeat.getAvgLength());
                        map.put(pattern.length(), lstD);

                    } else {
                        List<Double> lstD = new ArrayList<Double>();
                        lstD.add(repeat.getAvgLength());
                        map.put(pattern.length(), lstD);

                    }


                }
            }


        }
        //normalize
        //normalize
        for (int i = 0; i < consensus.length; i++) {
            if (consensusSize[i] > 0)
                consensus[i] = consensus[i] / consensusSize[i];

        }
//        for (int i = 0; i < consensusLen.length; i++) {
//            System.out.print(NumberFormatUtils.format(consensus[i]) + ",");
//        }


        List<Region> lstRegions = demoLength.findRegion(consensus);
        for (Region rc : lstRegions) {


            Atom[] fraAtoms = Fragement.getFragementsofAtoms(atoms, rc.getStart(), rc.getEnd());
            double predLenCE = demoLength.getPredictedLengthByCE(fraAtoms);

            // predict by consensus
            double predLenCONSENCUS = 0.0;
            int index = 0;
            for (int i = rc.getStart(); i <= rc.getEnd(); i++) {
                if (consensusLen[i] < 1000) {
                    predLenCONSENCUS += consensusLen[i];
                    index++;
                }
            }
            if (index > 0)
                predLenCONSENCUS = predLenCONSENCUS / rc.size();

            // predict by Raphael
            Raphael raphael = new Raphael(fraAtoms);
            double predRaphael = raphael.getRepeatLength();

//            List<Double> lstPreds = new ArrayList<Double>();
//            if (predLenCE > 7)
//                lstPreds.add(predLenCE);
//            if (predRaphael > 7) {
//                lstPreds.add(predRaphael);
//            }
//            if (predLenCONSENCUS > 0) {
//                lstPreds.add(predLenCONSENCUS);
//            }
//
            double predLen = 0.0;// finalLen(lstPreds.toArray(new Double[lstPreds.size()]));
            if (predRaphael - 5 <= predLenCE & predLenCE <= predRaphael + 5) {
                predLen = (predLenCE + predRaphael) / 2;
            } else {
                if (predLenCE > 20)
                    predLen = predLenCE;
                else if (predLenCE < 20 && predRaphael > 20 && raphael.getTotalScore() > 5) {
                    predLen = predRaphael;
                } else {
                    predLen = predLenCONSENCUS;
                }

            }

            //String predRegion = PdbTools.getResSeq(atoms, rc.getStart()) + "-" + PdbTools.getResSeq(atoms, rc.getEnd());

            if (predLen > 0) {
                RepeatRegion repeatRegion = new RepeatRegion(rc.getStart(), rc.getEnd());
                repeatRegion.setSeqStart(PdbTools.getResSeq(atoms, rc.getStart()));
                repeatRegion.setSeqEnd(PdbTools.getResSeq(atoms, rc.getEnd()));
                repeatRegion.setPredLen((int) predLen);
                repeatRegion.setPredRU((int) (rc.size() / predLen));
                if (repeatRegion.getPredRU() <= 2)
                    repeatRegion.setPredRU(2);
                predRegions.add(repeatRegion);
            }
        }

        return predRegions;
    }


}
