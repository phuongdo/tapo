package org.cnrs.crbm.lib.opt.simpleGa;

import java.util.ArrayList;
import java.util.List;

import org.cnrs.crbm.lib.align.unpublish.Block;
import org.cnrs.crbm.lib.align.unpublish.MatrixCM;

public class FitnessCalc {

	static byte[] solution = new byte[64];

	static Block[] blocks = null;
	static MatrixCM matrixCM1 = null;
	static MatrixCM matrixCM2 = null;

	public static void setMatrixCM(MatrixCM matrix1, MatrixCM matrix2) {
		matrixCM1 = matrix1;
		matrixCM2 = matrix2;

	}

	public static void setBlock(List<Block> aListblocks) {
		blocks = aListblocks.toArray(new Block[aListblocks.size()]);
	}

	/* Public methods */
	// Set a candidate solution as a byte array
	public static void setSolution(byte[] newSolution) {
		solution = newSolution;
	}

	// To make it easier we can use this method to set our candidate solution
	// with string of 0s and 1s
	static void setSolution(String newSolution) {
		solution = new byte[newSolution.length()];
		// Loop through each character of our string and save it in our byte
		// array
		for (int i = 0; i < newSolution.length(); i++) {
			String character = newSolution.substring(i, i + 1);
			if (character.contains("0") || character.contains("1")) {
				solution[i] = Byte.parseByte(character);
			} else {
				solution[i] = 0;
			}
		}
	}

	// Calculate inidividuals fittness by comparing it to our candidate solution
	static int getFitness(Individual individual) {
		int fitness = (int) calculateEnergy(individual);
		return fitness;
	}

	// Calculate the constans

	public static double calculateEnergy(Individual individual) {

		double totalEnergy = 0.0;
		List<Block> nblocks = getChoosedBlocks(individual);

		// calculate gaps

		int gaps = 0;
		int overllap = 0;
		int contactMap = 0;
		int leng = 0;

		List<Integer> list1 = new ArrayList<Integer>();
		List<Integer> list2 = new ArrayList<Integer>();

		// X first
		Block.sortByX(nblocks);
		for (int i = 0; i < nblocks.size() - 1; i++) {

			Block blockCurrent = nblocks.get(i);
			Block blockNext = nblocks.get(i + 1);
			if (blockNext.getStartX() > blockCurrent.getEndX())
				gaps += blockNext.getStartX() - blockCurrent.getEndX() - 1;

			else
				overllap = -blockNext.getStartX() + blockCurrent.getEndX()
						+ 1;

			for (int x = blockCurrent.getStartX(), y = blockCurrent.getStartY(); x <= blockCurrent
					.getEndX(); x++, y++) {

				list1.add(x);
				list2.add(y);
				leng++;

			}

		}

		// // Secondly, consider to Y
		// Collections.sort(nblocks, new BlockComparatorY());
		// for (int i = 0; i < nblocks.size() - 1; i++) {
		//
		// Block blockCurrent = nblocks.get(i);
		//
		// Block blockNext = nblocks.get(i + 1);
		// if (blockNext.getStartProY() - blockCurrent.getEndProY() - 1 > 0)
		//
		// gaps += blockNext.getStartProY() - blockCurrent.getEndProY()
		// - 1;
		//
		// else
		// overllap += -blockNext.getStartProY()
		// + blockCurrent.getEndProY() + 1;
		//
		// }

		Integer[] arr1 = list1.toArray(new Integer[list1.size()]);
		Integer[] arr2 = list2.toArray(new Integer[list2.size()]);

		for (int i = 0; i < arr1.length - 1; i++) {
			for (int j = i + 1; i < arr1.length; i++) {

				// edge x ,y
				int x = arr1[i], y = arr2[j];

				// edge u, v
				int u = arr2[i], v = arr2[j];

				if (matrixCM1.getCell(x, y) == matrixCM2.getCell(u, v)
						&& matrixCM1.getCell(x, y) == 1)
					contactMap++;
			}

		}

		// calculate contact map score
		totalEnergy = contactMap - overllap - gaps;
		return totalEnergy;
	}

	public static List<Block> getChoosedBlocks(Individual individual) {

		List<Block> nblocks = new ArrayList<Block>();
		for (int i = 0; i < individual.size(); i++) {
			if (individual.getGene(i) == 1) {
				nblocks.add(blocks[i]);
			}
		}

		return nblocks;

	}

	// Get optimum fitness
	public static int getMaxFitness() {
		int maxFitness = solution.length;
		return maxFitness;
	}
}