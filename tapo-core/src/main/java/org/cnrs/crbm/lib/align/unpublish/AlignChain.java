package org.cnrs.crbm.lib.align.unpublish;

import java.util.ArrayList;
import java.util.List;

public class AlignChain {

	private double rmsd;
	private int alignLen;
	private int gaps;
	private List<Block> listblockAliged = new ArrayList<Block>();

	public double getRmsd() {
		return rmsd;
	}

	public void setRmsd(double rmsd) {
		this.rmsd = rmsd;
	}

	public List<Block> getListblockAliged() {
		return listblockAliged;
	}

	public void setListblockAliged(List<Block> listblockAliged) {
		this.listblockAliged = listblockAliged;
	}

	public int getAlignLen() {
		return alignLen;
	}

	public void setAlignLen(int alignLen) {
		this.alignLen = alignLen;
	}

	public int getGaps() {
		return gaps;
	}

	public void setGaps(int gaps) {
		this.gaps = gaps;
	}

}
