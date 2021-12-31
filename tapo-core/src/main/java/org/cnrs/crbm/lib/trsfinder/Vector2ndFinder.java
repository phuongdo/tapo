package org.cnrs.crbm.lib.trsfinder;

import java.util.ArrayList;
import java.util.List;

import org.cnrs.crbm.lib.repeats.module.Pro2Vector;
import org.cnrs.crbm.lib.repeats.module.ProVector;
import org.cnrs.crbm.lib.repeats.RepeatScore;
import org.cnrs.crbm.lib.repeats.RepeatStore;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Using only vector of secondary structural elements for finding TRs Minimum
 * elements are 2 elements, and maximal elements are 6 elements.
 *
 * @author pdoviet
 *
 */

@Deprecated
public class Vector2ndFinder extends Finder {

	final static int MIN_ELEMENTS = 2;
	final static int MAX_ELEMENTS = 6;

	public Vector2ndFinder(Features features) {
		super(features);
		this.name = "VECTORS";

	}

	static Logger logger = LoggerFactory.getLogger(Vector2ndFinder.class);

	@Override
	public void findRepeat(Features features) {
		// List<Repeat> repeats = new ArrayList<Repeat>();
		try {
			Pro2Vector pro2Vector = new Pro2Vector();
			List<RepeatScore> repeat_scores = new ArrayList<RepeatScore>();
			// only use secondary structure. 04/07/2014
			List<ProVector> ssVectors = ProVector.toSecondaryVector(features
					.getVectors());
			for (int winsize = MIN_ELEMENTS; winsize <= MAX_ELEMENTS; winsize++) {

				repeat_scores.addAll(pro2Vector.getRepeats(ssVectors,
						features.getAtoms(), winsize));
			}
			// sort and display
			RepeatScore.sortBySD(repeat_scores);
			// get the best repeats
			// 2014/07/24
			// we change this funtion to get all the TRs instead of getting the
			// best TRs. Note that we don't check the fragments which are
			// overlapping with others...
			if (repeat_scores.size() > 0) {

				for (int index = 0; index < repeat_scores.size(); index++) {

					Repeat repeat = new Repeat();
					for (RepeatStore rs : repeat_scores.get(index).getRepeats()) {
						repeat.getRepeats().add(
								new RepeatContent(rs.getStart(), rs.getEnd()));
					}
					repeat.sortByPosition();
					repeats.add(repeat);
				}
			}
		} catch (Exception e) {
			logger.error(e.getMessage() + " with pdb " + features.getPdbCode()
					+ "_" + features.getPdbChain());
		}
		// return repeats;
	}
}
