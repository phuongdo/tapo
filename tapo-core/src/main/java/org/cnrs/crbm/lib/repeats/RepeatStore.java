package org.cnrs.crbm.lib.repeats;

import org.cnrs.crbm.lib.repeats.module.ProVector;

import java.util.ArrayList;
import java.util.List;

public class RepeatStore {

	private int start;
	private int end;

	public RepeatStore() {
	}

	public RepeatStore(List<ProVector> sub_list, double rmsd) {

		this.sub_list = sub_list;
		this.rmsd = rmsd;
	}

	List<ProVector> sub_list = new ArrayList<ProVector>();
	double rmsd;

	public List<ProVector> getSub_list() {
		return sub_list;
	}

	public void setSub_list(List<ProVector> sub_list) {
		this.sub_list = sub_list;
	}

	public double getRmsd() {
		return rmsd;
	}

	public void setRmsd(double rmsd) {
		this.rmsd = rmsd;
	}

	public int size() {

		return sub_list.size();
	}

	public int getResiSize() {

		return sub_list.get(sub_list.size() - 1).getOriginalEnd().getGroup()
				.getResidueNumber().getSeqNum()
				- sub_list.get(0).getOriginalStart().getGroup()
						.getResidueNumber().getSeqNum();
	}

	@Override
	public String toString() {
		if (sub_list.size() > 0)
			return

			"<resi "
					+ sub_list.get(0).getOriginalStart().getGroup()
							.getResidueNumber()
					+ "-"
					+ sub_list.get(sub_list.size() - 1).getOriginalEnd()
							.getGroup().getResidueNumber() + "> pattern: "
					+ this.getPatterOfVectors();

		else
			return "";
	}

	public String getPatterOfVectors() {
		String pattern_sub = "";
		for (ProVector v : sub_list) {
			pattern_sub += v.getType();
		}
		return pattern_sub;
	}

	public int getStart() {
		return sub_list.get(0).getPosStart();
	}

	public void setStart(int start) {
		this.start = start;
	}

	public int getEnd() {
		return sub_list.get(sub_list.size() - 1).getPosEnd();
	}

	public void setEnd(int end) {
		this.end = end;
	}

}
