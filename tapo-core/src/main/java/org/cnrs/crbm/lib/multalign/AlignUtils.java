package org.cnrs.crbm.lib.multalign;

import java.util.ArrayList;
import java.util.List;

public class AlignUtils {

	public static List<List<Integer>> getMultipleAligns(int[][] multi,
			int nRow, int nCol) {
		// refine
		//
		// for (int i = 0; i < nRow; i++) {
		// for (int j = 0; j < nCol; j++) {
		//
		// System.out.print(multi[i][j] + " ");
		// }
		//
		// System.out.println();
		//
		// }

		for (int i = 0; i < nRow; i++) {

			int a = 0;
			// initial for a;
			for (int j = 0; j < nCol; j++) {
				if (multi[i][j] > 0) {
					a = multi[i][j] - 1;
					break;
				}
			}
			for (int j = 0; j < nCol; j++) {

				int cur = multi[i][j];
				if (cur > 0)
					a = cur + 1;
				else {

					multi[i][j] = a;
				}
			}

		}

		List<List<Integer>> aligns = new ArrayList<List<Integer>>();
		for (int i = 0; i < nRow; i++) {
			aligns.add(new ArrayList<Integer>());
		}

		for (int i = 0; i < nCol; i++) {
			int next = i + 1;

			if (next == nCol) {
				for (int j = 0; j < nRow; j++) {
					int cur = multi[j][i];
					aligns.get(j).add(cur);

				}
			} else {

				// find maxgap
				int maxgap = 0;

				for (int j = 0; j < nRow; j++) {
					int acur = multi[j][i];
					int anext = multi[j][i + 1];
					int agap = anext - acur - 1;

					maxgap = Math.max(maxgap, agap);

				}

				// build pairAlign

				for (int j = 0; j < nRow; j++) {
					int acur = multi[j][i];
					int anext = multi[j][i + 1];
					int agap = anext - acur - 1;
					int aextend = maxgap - agap;

					for (int k = acur; k < anext; k++) {
						aligns.get(j).add(k);
					}
					for (int e = 0; e < aextend; e++) {
						// extend to next;
						aligns.get(j).add(-1);
					}
				}

			}

		}
		return aligns;

	}
}
