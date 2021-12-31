package org.cnrs.crbm.lib.repeats;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import org.cnrs.crbm.lib.utils.NumberFormatUtils;

public class RepeatScore {

	List<RepeatStore> repeats = new ArrayList<RepeatStore>();

	public List<RepeatStore> getRepeats() {
		return repeats;
	}

	double score_tandem;
	double score_repeat_size;

	public RepeatScore() {
	}

	public RepeatScore(List<RepeatStore> repeats, double score_tandem,
			double score_repeat_size) {

		this.repeats = repeats;
		this.score_repeat_size = score_repeat_size;
		this.score_tandem = score_tandem;
	}

	public double getAvgLength() {
		double avg = 0.0;
		for (RepeatStore r : repeats) {
			avg += r.getResiSize();
		}
		return avg / repeats.size();
	}

	public double getAvgRMSD() {
		double avg = 0.0;
		for (RepeatStore r : repeats) {
			avg += r.getRmsd();
		}
		return avg / repeats.size();
	}

	public double getSD() {

		double avg = this.getAvgLength();
		double sd = 0;
		for (RepeatStore r : repeats) {
			sd += Math.pow(r.getResiSize() - avg, 2);
		}
		return Math.sqrt(sd / repeats.size());

	}

	public double getScore() {

		return this.getSD() * this.getAvgRMSD() * this.score_repeat_size
				* this.score_tandem;// * (1 / (this.getAvgLength()));

	}

	public static void sortBySD(List<RepeatScore> rs) {

		Collections.sort(rs, new Comparator<RepeatScore>() {
			public int compare(RepeatScore s1, RepeatScore s2) {
				if (s1.getScore() > s2.getScore())
					return 1;
				else if (s1.getScore() < s2.getScore())
					return -1;
				else
					return 0;
			}
		});

	}

	@Override
	public String toString() {

		StringBuilder builder = new StringBuilder();
		builder.append("--------\n");
		builder.append("score: " + NumberFormatUtils.format(this.getScore())
				+ "\n");
		builder.append("length : "
				+ NumberFormatUtils.format(this.getAvgLength()) + "\n");
		for (RepeatStore r : repeats) {

			builder.append(r + "\n");
		}
		builder.append("\n");

		return builder.toString();
	}
}
