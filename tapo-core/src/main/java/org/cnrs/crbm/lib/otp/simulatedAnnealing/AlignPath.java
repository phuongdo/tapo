package org.cnrs.crbm.lib.otp.simulatedAnnealing;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.cnrs.crbm.lib.align.unpublish.Block;
import org.cnrs.crbm.lib.align.unpublish.MatrixCM;

public class AlignPath {
	// Holds our path of alignment
	private ArrayList path = new ArrayList();
	// Cache
	private double distance = 0;// type of score

	static MatrixCM matrixCM1 = null;
	static MatrixCM matrixCM2 = null;
	byte[] Xmark = new byte[64];
	byte[] Ymark = new byte[64];

	public void setMatrixCM(MatrixCM matrix1, MatrixCM matrix2) {
		matrixCM1 = matrix1;
		matrixCM2 = matrix2;
		Xmark = new byte[matrixCM1.getSize()];
		Ymark = new byte[matrixCM2.getSize()];

	}

	// Constructs a blank tour
	public AlignPath() {
		for (int i = 0; i < BlockManger.numberOfBlocks(); i++) {
			path.add(null);
		}
	}

	// Constructs a tour from another tour
	public AlignPath(ArrayList path) {
		// clone a path
		ArrayList<Block> clonedList = new ArrayList<Block>(path.size());
		for (Block block : (ArrayList<Block>) path) {
			clonedList.add(new Block(block));
		}
		this.path = clonedList;
	}

	public boolean isOverllapedWithCurrentPath(Block block) {

		if (isCrossover(block))
			return true;

		for (int i = block.getStartX(), j = block.getStartY(); i < block
				.getEndX(); i++, j++) {
			if (this.Xmark[i] == 1 || this.Ymark[j] == 1)
				return true;
		}

		return false;

	}

	public boolean isCrossover(Block moveBlock) {

		ArrayList<Block> clonedList = new ArrayList<Block>(path.size());
		for (Block block : (ArrayList<Block>) path) {
			clonedList.add(new Block(block));
		}
		clonedList.add(moveBlock);
		Block.sortByX(clonedList);
		for (int i = 0; i < clonedList.size() - 1; i++) {
			Block cur = clonedList.get(i);
			Block next = clonedList.get(i + 1);
			if (cur.getStartY() > next.getStartY())
				return true;
		}

		return false;
	}

	public void releaseMarkBlock(Block block) {

		for (int i = block.getStartX(); i < block.getEndX(); i++) {
			Xmark[i] = 0;
		}
		for (int i = block.getStartY(); i < block.getEndY(); i++) {
			Ymark[i] = 0;
		}
		// System.out.println();
	}

	public void markXYForCheckingOverllaping() {

		// reset;
		Xmark = new byte[matrixCM1.getSize()];
		Ymark = new byte[matrixCM2.getSize()];
		for (Block block : (ArrayList<Block>) path) {
			for (int i = block.getStartX(); i < block.getEndX(); i++) {
				Xmark[i] = 1;
			}
			for (int i = block.getStartY(); i < block.getEndY(); i++) {
				Ymark[i] = 1;
			}
		}
		// System.out.println();
	}

	// Returns our path information
	public ArrayList getPath() {
		return path;
	}

	// Creates a random individual
	public void generateIndividual(List<Block> bagsOfBlocks) {

		// add randomly the path
		Random random = new Random();

		for (int i = 0; i < 4; i++) {
			int rnumber = random.nextInt(bagsOfBlocks.size());
			if (!isOverllapedWithCurrentPath(bagsOfBlocks.get(rnumber))) {
				System.out.println(bagsOfBlocks.get(rnumber));
				path.add(bagsOfBlocks.get(rnumber));
				this.markXYForCheckingOverllaping();
			}
		}
		//
		// Block block1 = new Block();
		// block1.setStartX(0);
		// block1.setStartY(0);
		// block1.setEndX(57);
		// block1.setEndY(57);
		//
		// Block block2 = new Block();
		// block2.setStartX(60);
		// block2.setStartY(58);
		// block2.setEndX(98);
		// block2.setEndY(96);
		// path.add(block2);
		// path.add(block1);

		// this.markXYForCheckingOverllaping();

	}

	// Gets a block from the path
	public Block getBlock(int blockPosition) {
		return (Block) path.get(blockPosition);
	}

	// Sets a Block in a certain position within a Paht
	public void setBlock(int pathPosition, Block block) {
		path.set(pathPosition, block);
		// If the tours been altered we need to reset the fitness and distance
		distance = 0;
	}

	// Gets the total enery of the current path
	public double getEnergy() {
		if (distance == 0) {
			double totalEnergy = 0;
			List<Block> nblocks = path;
			// calculate gaps
			int gaps = 0;
			int overllap = 0;
			int contactMap = 0;
			int leng = 0;

			List<Integer> list1 = new ArrayList<Integer>();
			List<Integer> list2 = new ArrayList<Integer>();
			// X first
			Block.sortByX(nblocks);
			for (int i = 0; i < nblocks.size(); i++) {
				Block blockCurrent = nblocks.get(i);
				// Block blockNext = nblocks.get(i + 1);
				// if (blockNext.getStartX() > blockCurrent.getEndProX())
				// gaps += blockNext.getStartX() - blockCurrent.getEndProX()
				// - 1;
				// else
				// overllap = -blockNext.getStartX()
				// + blockCurrent.getEndProX() + 1;

				for (int x = blockCurrent.getStartX(), y = blockCurrent
						.getStartY(); x <= blockCurrent.getEndX(); x++, y++) {
					list1.add(x);
					list2.add(y);
					leng++;
				}
			}

			Integer[] arr1 = list1.toArray(new Integer[list1.size()]);
			Integer[] arr2 = list2.toArray(new Integer[list2.size()]);
			int nContactMap1 = 0;
			int nContactMap2 = 0;
			for (int i = 0; i < arr1.length - 1; i++) {
				for (int j = i + 1; j < arr1.length; j++) {
					// edge x ,y
					int x = arr1[i], y = arr1[j];
					// edge u, v
					int u = arr2[i], v = arr2[j];
					// System.out.println(matrixCM1.getCell(x, y) + "-"
					// + matrixCM2.getCell(u, v));

					nContactMap1 += matrixCM1.getCell(x, y);
					nContactMap2 += matrixCM2.getCell(u, v);

					contactMap += matrixCM1.getCell(x, y)
							* matrixCM2.getCell(u, v);
					// if (matrixCM1.getCell(x, y) == matrixCM2.getCell(u, v)
					// && matrixCM1.getCell(x, y) == 1)
					// contactMap++;
				}

			}

			// System.out.println(2 * contactMap / (nContactMap1 +
			// nContactMap2));

			// calculate contact map score
			// totalEnergy = 2 * (double) contactMap
			// / (nContactMap1 + nContactMap2);

			totalEnergy = contactMap;
			distance = totalEnergy;
		}
		return distance;
	}

	// Get number of blocks on our path
	public int pathSize() {
		return path.size();
	}

	@Override
	public String toString() {
		String geneString = "|";
		for (int i = 0; i < pathSize(); i++) {
			geneString += getBlock(i) + "|";
		}
		return geneString;
	}

	public byte[] getXmark() {
		return Xmark;
	}

	public void setXmark(byte[] xmark) {
		Xmark = xmark;
	}

	public byte[] getYmark() {
		return Ymark;
	}

	public void setYmark(byte[] ymark) {
		Ymark = ymark;
	}

}
