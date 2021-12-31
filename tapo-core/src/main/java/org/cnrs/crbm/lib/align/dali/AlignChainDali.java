package org.cnrs.crbm.lib.align.dali;

import java.util.HashMap;
import java.util.Map;

public class AlignChainDali {

	private double rmsd;
	private int gaps;
	private int chainLength;
	private Map<Integer, Integer> alignPairs = new HashMap<Integer, Integer>();

	public double getRmsd() {
		return rmsd;
	}

	public void setRmsd(double rmsd) {
		this.rmsd = rmsd;
	}

	public Map<Integer, Integer> getAlignPairs() {
		return alignPairs;
	}

	public void setAlignPairs(Map<Integer, Integer> alignPairs) {
		this.alignPairs = alignPairs;
	}

	public int getGaps() {
		return gaps;
	}

	public void setGaps(int gaps) {
		this.gaps = gaps;
	}

	public int getChainLength() {
		return chainLength;
	}

	public void setChainLength(int chainLength) {
		this.chainLength = chainLength;
	}

}
