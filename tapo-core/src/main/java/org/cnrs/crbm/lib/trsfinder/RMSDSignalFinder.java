package org.cnrs.crbm.lib.trsfinder;

import org.apache.commons.lang3.ArrayUtils;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.cnrs.crbm.lib.math.Filter;
import org.cnrs.crbm.lib.math.SignalPeriod;
import org.cnrs.crbm.lib.math.SignalProcess;
import org.cnrs.crbm.lib.repeats.Fragement;
import org.cnrs.crbm.lib.repeats.module.FinderOutput;
import org.cnrs.crbm.lib.repeats.module.TMEvaluation;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * The major change of this method is that we convert the structure of protein
 * in to time series signals. Then, using Fourier Analysis, we find the
 * periodicity of repeats.
 * Created by pdoviet on 11/21/2014.
 */
public class RMSDSignalFinder extends Finder {

    static Logger logger = LoggerFactory.getLogger(RMSDSignalFinder.class);

    public RMSDSignalFinder(Features features) {
        super(features);
        this.name = "RMSD";
    }

    @Override
    public void findRepeat(Features features) {

        try {
            double maxSGscore = 0.0;
            double maxTMScore = 0.0;
            TMEvaluation tmEvaluation = new TMEvaluation();
            List<SignalPeriod> list = this.processStructure(atoms);
            for (SignalPeriod period : list) {
                //System.out.println(period);
                FinderOutput finderOut = tmEvaluation.findRepeatsBasedOnTMScore((int) period.getPeriod(), period.getStart(), period.getEnd(), features);
                List<Repeat> lrepeats = finderOut.getRepeats();
                maxSGscore = Math.max(maxSGscore, period.getSgScore());
                maxTMScore = Math.max(maxTMScore, finderOut.getCombineScore().getTmScore());
                if (lrepeats.size() > 0)
                    this.repeats.addAll(lrepeats);

            }

            this.combineScore.setSigScore(maxSGscore);
            this.combineScore.setTmScore(maxTMScore);

        } catch (Exception e) {
            //e.printStackTrace();
        }

    }

    SignalProcess signalProcess = new SignalProcess();

    public List<SignalPeriod> processStructure(Atom[] atoms) {

        List<SignalPeriod> list = new ArrayList<SignalPeriod>();

        // segment a structure
        int winsize = 200;
        for (int i = 0; i < atoms.length - winsize; i = i + winsize / 2) {
            try {
                int start = i;
                int end = i + winsize;
                if (end > atoms.length)
                    end = atoms.length;
                Atom[] v = Arrays.copyOfRange(atoms, start, end);
                if (v.length < winsize)
                    continue;

                // do something here with this segment
                // generate signals

                SignalPeriod signalPeriod = this.analysisSegmentAtoms(v);
                signalPeriod.setStart(start);
                signalPeriod.setEnd(end);
                list.add(signalPeriod);

            } catch (Exception ex) {
                //
            }

        }

        return list;
    }

    public SignalPeriod analysisSegmentAtoms(Atom[] atoms) {
        List signals = new ArrayList<Double>();
        int winsize = 20;
        double maxSGscore = 0.0;
        int bestPeriod = 0;
        SignalPeriod signalPeriod = new SignalPeriod();
        for (int t = 0; t < atoms.length - winsize; t = t + 5) {
            Atom[] seedAtoms = Fragement.getFragementsofAtoms(atoms, t, t
                    + winsize);
            for (int j = 0; j < atoms.length - winsize; j++) {
                Atom[] compareAtoms = Fragement.getFragementsofAtoms(atoms, j,
                        j + winsize);
                double rmsd = 10;
                try {
                    rmsd = superimposer.superimposeSimple(seedAtoms,
                            compareAtoms);
                    signals.add(rmsd);
                } catch (StructureException e) {
//                e.printStackTrace();
                }
            }

            double[] his = ArrayUtils.toPrimitive((Double[]) signals.toArray(new Double[signals.size()]));
            try {
                his = Filter.filter(his);
                double period = signalProcess.predictPeriod(his);
                double sigScore = signalProcess.correlationWithLag(his, (int) period);

                if (maxSGscore < sigScore) {
                    maxSGscore = sigScore;
                    bestPeriod = (int) period;
                }
            } catch (Exception e) {
                e.printStackTrace();
            }

        }

        signalPeriod.setSgScore(maxSGscore);
        signalPeriod.setPeriod(bestPeriod);
        return signalPeriod;

    }

}
