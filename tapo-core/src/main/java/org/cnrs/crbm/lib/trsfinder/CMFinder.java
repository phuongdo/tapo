package org.cnrs.crbm.lib.trsfinder;

import org.cnrs.crbm.lib.math.FFTperiodicity;
import org.cnrs.crbm.lib.math.SignalPeriod;
import org.cnrs.crbm.lib.math.SignalProcess;
import org.cnrs.crbm.lib.repeats.module.FinderOutput;
import org.cnrs.crbm.lib.repeats.module.TMEvaluation;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.List;

public class CMFinder extends Finder {
	static Logger logger = LoggerFactory.getLogger(CMFinder.class);
	public CMFinder(Features features) {
		super(features);
		this.name = "CM";

	}

	@Override
	public void findRepeat(Features features) {
		try {

            double[] his = features.getCmHistos();
            SignalProcess signalProcess = new SignalProcess();
            try {
				TMEvaluation tmEvaluation = new TMEvaluation();
                List<SignalPeriod> list = signalProcess.predictPeriodicties(his);
                double maxSGscore = 0.0;
				double maxTMScore = 0.0;
				for (SignalPeriod period : list) {
                    FinderOutput finderOut =  tmEvaluation.findRepeatsBasedOnTMScore((int) period.getPeriod(), period.getStart(), period.getEnd(),features);
					List<Repeat> lrepeats = finderOut.getRepeats();
					maxSGscore = Math.max(maxSGscore,period.getSgScore());
					maxTMScore = Math.max(maxTMScore, finderOut.getCombineScore().getTmScore());
                    if (lrepeats.size() > 0)
                        this.repeats.addAll(lrepeats);
                }

				this.combineScore.setSigScore(maxSGscore);
				this.combineScore.setTmScore(maxTMScore);

            } catch (Exception e) {
                //e.printStackTrace();
            }

		} catch (Exception ex) {
			ex.printStackTrace();
			logger.error(ex.getMessage() + " with pdb " + features.getPdbCode()
					+ "_" + features.getPdbChain());

		}

	}

	private void variableCM(double[] cmHistos, int winsize) {
		FFTperiodicity fft = new FFTperiodicity();
		List<Double> posibleLengs = fft.findPeriodicityWithWindow(cmHistos,
				winsize);

		if (posibleLengs.size() < 4) {
			for (double repeatLeng : posibleLengs) {

				if (repeatLeng > 10) {

					this.repeats.addAll(AtomFinder.getRepeat(features,
							(int) repeatLeng));

				}
			}
		}

	}

	private void constantCM(double[] cmHistos, int winsize) {
		FFTperiodicity fft = new FFTperiodicity();
		List<Double> posibleLengs = fft.findPeriodicityWithConstantSignal(
				cmHistos, winsize);

		if (posibleLengs.size() < 4) {
			for (double repeatLeng : posibleLengs) {

				if (repeatLeng > 10) {
					// know repeats leng => scan db
					// System.out.println(repeatLeng);

					int[] labels = new int[cmHistos.length];
					int range = (int) (0.2 * repeatLeng);
					for (int i = 0; i < cmHistos.length; i++) {
						if (repeatLeng - 5 < cmHistos[i]
								&& cmHistos[i] < repeatLeng + 5) {
							// get regions

							labels[i] = 1;
						} else
							labels[i] = 0;

					}

					int bridgeSize = 4;
					for (int i = 0; i < labels.length - bridgeSize; i++) {
						// System.out.print(labels[i]);
						// bridge
						if (labels[i] == 1 && labels[i + bridgeSize] == 1) {
							for (int j = i; j < i + bridgeSize; j++) {
								labels[j] = 1;
							}

						}
					}

					// for (int i = 0; i < labels.length; i++) {
					// System.out.print(labels[i]);
					//
					// }
					int start = 0;
					int end = 0;
					for (int i = 0; i < labels.length; i++) {

						if (labels[i] == 1) {
							start = i;
							while (i < labels.length && labels[i] == 1) {
								i++;
							}
							i = end = i - 1;

							// System.out.println(start + "-" + end);
							if (end - start + 1 > 2 * repeatLeng) {
								this.repeats
										.addAll(AtomFinder.getRepeatRange(
												features, (int) repeatLeng,
												start, end));
							}

						}
					}

				}
			}
		}

	}
}
