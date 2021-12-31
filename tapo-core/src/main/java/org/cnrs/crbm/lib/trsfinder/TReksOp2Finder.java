package org.cnrs.crbm.lib.trsfinder;

import org.cnrs.crbm.lib.multalign.PSim;
import org.cnrs.crbm.lib.repeats.shortTRs.ShortTRsFilter;
import org.cnrs.crbm.lib.sadb.Conformation;
import org.cnrs.crbm.lib.sadb.Seq3d;
import org.cnrs.crbm.lib.sadb.Sequence3D;
import org.cnrs.crbm.lib.seqalign.RepeatSeq;
import org.cnrs.crbm.lib.seqalign.RepeatUnitSeq;
import org.cnrs.crbm.lib.seqalign.TREksWapper;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.ArrayList;
import java.util.List;

public class TReksOp2Finder extends Finder {
	static Logger logger = LoggerFactory.getLogger(TReksOp2Finder.class);

	public TReksOp2Finder(Features features) {
		super(features);
		this.name = "TREKS-Op2";
	}

	@Override
	public void findRepeat(Features features) {
		List<Repeat> repeats_copy = new ArrayList<Repeat>();
		Sequence3D sequence3D = new Sequence3D();
		String seqAl = sequence3D.getSAV3(features.getStrSeqAp());
		Seq3d seq3D = new Seq3d(seqAl);
		String pattern = seq3D.getSeq();
		List<Conformation> confors = seq3D.getConfors();
		try {
			double maxPSim = 0.0;
			TREksWapper TReks = new TREksWapper("", pattern);
			List<RepeatSeq> repeatsSq = TReks.getRepeatSeqs();
			for (RepeatSeq repeat : repeatsSq) {
				// each of repeats
				//abc
				Repeat a_repeat = new Repeat();
				List<String> msa = new ArrayList<String>();
				for (RepeatUnitSeq ru : repeat.getUnits()) {
					msa.add(ru.getStrAlign());
					msa.add(ru.getStrAlign());
					int start = ru.getStart();
					int end = ru.getEnd();
					int leng = 0;
					// convert current to repeats
					// get real position
					// count in real leng of resi
					for (int i = start; i <= end; i++) {
						leng += confors.get(i).getEnd()
								- confors.get(i).getStart() + 1;
					}
					start = confors.get(start).getStart();
					a_repeat.getRepeats().add(
							new RepeatContent(start, start + leng - 1));
				}
				PSim pSim = new PSim(msa);
				a_repeat.sortByPosition();
				if (ShortTRsFilter.filter(a_repeat, features.getStrSeqAp()) && a_repeat.getAvgLength() >= 5) {
					// check repeat
					maxPSim = Math.max(maxPSim, pSim.getSimilarity());
					if (pSim.getSimilarity() >= 0.7) {
						this.repeats.add(a_repeat);
					}
				}
			}

			this.combineScore.setPsimScore(maxPSim);

		} catch (Exception ex) {
			logger.error(ex.getMessage() + " with pdb " + features.getPdbCode()
					+ "_" + features.getPdbChain());
		}

		//this.repeats = repeats_copy;

	}
}
