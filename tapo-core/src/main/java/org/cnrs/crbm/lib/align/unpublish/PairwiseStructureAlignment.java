package org.cnrs.crbm.lib.align.unpublish;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Calc;
import org.biojava.nbio.structure.StructureException;
import org.cnrs.crbm.lib.cm.ContactMapIO;
import org.cnrs.crbm.lib.opt.simpleGa.Algorithm;
import org.cnrs.crbm.lib.opt.simpleGa.FitnessCalc;
import org.cnrs.crbm.lib.opt.simpleGa.Population;
import org.cnrs.crbm.lib.otp.simulatedAnnealing.AlignPath;

public class PairwiseStructureAlignment {

	double cutOff = 8;
	MatrixCM matrixCM1;
	MatrixCM matrixCM2;
	int[][] dynamicMatrix;
	Atom[] ca1;
	Atom[] ca2;
	ContactMapIO io = new ContactMapIO();
	int wSize = 10;// like DALI

	public AlignChain align(Atom[] ca1, Atom[] ca2) {

		matrixCM1 = this.getMatrixCM(ca1);
		matrixCM2 = this.getMatrixCM(ca2);
		dynamicMatrix = new int[ca1.length][ca2.length];
		this.ca1 = ca1;
		this.ca2 = ca2;
		// generate seed blocks
		List<Block> seeds = this.genderateSeedBlocks();
		List<Block> best = this.simulatedAnnealing(seeds);
		this.saveToPNG(best);

		return new AlignChain();
	}

	private void saveToPNG(List<Block> blocks) {
		dynamicMatrix = new int[ca1.length][ca2.length];
		for (Block block : blocks) {
			System.out.println(block);
			for (int i = block.getStartX(), j = block.getStartY(); i <= block
					.getEndX(); i++, j++) {
				dynamicMatrix[i][j] = 1;

			}

		}

		io.saveMatrixtoPNG(dynamicMatrix, ca1.length, ca2.length);
	}

	private List<Block> genderateSeedBlocks() {
		List<Block> blocks = new ArrayList<Block>();
		for (int i = 0; i < ca1.length - wSize; i++) {
			int start = i;
			int end = i + wSize - 1;
			MatrixCM subCM = matrixCM1.getSubMatrixCM(start, end);
			for (int j = 0; j < ca2.length - wSize; j++) {
				int start2 = j;
				int end2 = j + wSize - 1;
				MatrixCM subCM2 = matrixCM2.getSubMatrixCM(start2, end2);

				if (subCM.isOverllap(subCM2)) {
					Block seedBlock = new Block();
					seedBlock.setStartX(start);
					seedBlock.setStartY(start2);
					seedBlock.setEndX(end);
					seedBlock.setEndY(end2);
					// big
					seedBlock.expand(ca1.length, ca2.length);
					seedBlock.expand(ca1.length, ca2.length);
					seedBlock.expand(ca1.length, ca2.length);

					blocks.add(seedBlock);

				}

			}

		}

		return blocks;

	}

	private void gaMethod() {

		/**
		 * GA MEthod
		 */

		// FitnessCalc.setBlock(seeds);
		FitnessCalc.setMatrixCM(matrixCM1, matrixCM2);
		// Individual.setDefaultGeneLength(seeds.size());
		// Create an initial population
		Population myPop = new Population(50, true);

		// Evolve our population until we reach an optimum solution
		int generationCount = 0;
		while (generationCount < 100) {
			generationCount++;
			System.out.println("Generation: " + generationCount + " Fittest: "
					+ myPop.getFittest().getFitness());
			myPop = Algorithm.evolvePopulation(myPop);

		}
		System.out.println("Solution found!");
		System.out.println("Generation: " + generationCount);
		System.out.println("Genes:");
		System.out.println(myPop.getFittest());
		List<Block> optimalSolutions = FitnessCalc.getChoosedBlocks(myPop
				.getFittest());
	}

	private List<Block> simulatedAnnealing(List<Block> bagsOfBlocks) {
		// Set initial temp
		double temp = 1000;
		// Cooling rate
		double coolingRate = 0.003;
		// Initialize intial solution by randomly choose books from block seeds

		AlignPath currentSolution = new AlignPath();
		currentSolution.setMatrixCM(matrixCM1, matrixCM2);
		currentSolution.generateIndividual(bagsOfBlocks);
		currentSolution.markXYForCheckingOverllaping();

		System.out.println("Initial solution distance: "
				+ currentSolution.getEnergy());

		// Set as current best
		AlignPath best = new AlignPath(currentSolution.getPath());

		// Random random = new Random();
		// Loop until system has cooled
		while (temp > 1) {
			// Create new neighbour tour
			AlignPath newSolution = new AlignPath(currentSolution.getPath());
			newSolution.setMatrixCM(matrixCM1, matrixCM2);
			newSolution.markXYForCheckingOverllaping();
			// Change new state?
			// this.mover(bagsOfBlocks, newSolution);
			// move
			if (newSolution.pathSize() > 0) {
				this.mover(bagsOfBlocks, newSolution);

			}

			// Get energy of solutions
			double currentEnergy = currentSolution.getEnergy();
			double neighbourEnergy = newSolution.getEnergy();
			// System.out.println(temp + " next enery: " + neighbourEnergy);

			// Decide if we should accept the neighbour
			if (acceptanceProbability(currentEnergy, neighbourEnergy, temp) > Math
					.random()) {
				currentSolution = new AlignPath(newSolution.getPath());
			}

			// Keep track of the best solution found
			if (currentSolution.getEnergy() > best.getEnergy()) {
				best = new AlignPath(currentSolution.getPath());
				System.out.println("jump to:" + best.getEnergy());

			}

			// Cool system
			temp *= 1 - coolingRate;
			// temp--;

		}

		System.out.println("Final solution distance: " + best.getEnergy());
		// System.out.println("Tour: " + best);

		return best.getPath();

	}

	private void mover(List<Block> bagofBlocks, AlignPath alignPath) {

		// choose
		double probilityOfCreate = 0.4;
		double probilityOfRemove = 0.3;
		double probilityOfShrink = 0.2;
		double probilityOfExpand = 0.2;
		double probilityOfShift = 0.4;
		double probilityOfSplit = 0.3;
		double probilityOfMerger = 0.2;
		Random random = new Random();

		if (Math.random() < probilityOfCreate && alignPath.pathSize() < 5) {
			// create new block from bagOfBlock
			int rnumber = random.nextInt(bagofBlocks.size());
			Block block = bagofBlocks.get(rnumber);
			if (block.getEndY() - block.getStartY() < 0)
				System.out.println("hai");
			// check contain
			if (!alignPath.isOverllapedWithCurrentPath(block))
				alignPath.getPath().add(block);

			alignPath.markXYForCheckingOverllaping();
		}
		if (Math.random() < probilityOfShrink) {
			// shrink
			int i = random.nextInt(alignPath.pathSize());
			// release before shirk
			Block block = (Block) alignPath.getPath().get(i);
			alignPath.releaseMarkBlock(block);
			((Block) alignPath.getPath().get(i)).shrink(ca1.length, ca2.length);
			alignPath.markXYForCheckingOverllaping();

		}
		if (Math.random() < probilityOfExpand) {
			// expand
			int i = random.nextInt(alignPath.pathSize());
			Block block = (Block) alignPath.getPath().get(i);
			Block expandBlock = new Block(block);
			expandBlock.expand(ca1.length, ca2.length);
			alignPath.releaseMarkBlock(block);
			if (!alignPath.isOverllapedWithCurrentPath(expandBlock)) {
				alignPath.setBlock(i, expandBlock);

			}
			alignPath.markXYForCheckingOverllaping();
		}

		if (Math.random() < probilityOfShift) {
			// expand
			int i = random.nextInt(alignPath.pathSize());

			Block block = (Block) alignPath.getPath().get(i);
			Block shiftBlock = new Block(block);
			alignPath.releaseMarkBlock(block);
			shiftBlock.shift(ca1.length, ca2.length);

			if (!alignPath.isOverllapedWithCurrentPath(shiftBlock)) {
				alignPath.setBlock(i, shiftBlock);

			}
			alignPath.markXYForCheckingOverllaping();
		}

		if (Math.random() < probilityOfSplit && alignPath.pathSize() < 4) {
			// split
			int i = random.nextInt(alignPath.pathSize());
			Block block = (Block) alignPath.getPath().get(i);
			if (block.size() > 15) {

				List<Block> blocks = block.split();
				alignPath.getPath().remove(i);
				alignPath.getPath().addAll(blocks);

			}

		}

		if (Math.random() < probilityOfMerger) {
			// find block to merge
			Block.sortByX(alignPath.getPath());

			// block position
			int index = -1;
			for (int i = 0; i < alignPath.getPath().size() - 1; i++) {

				Block block1 = (Block) alignPath.getPath().get(i);
				Block blockNext = (Block) alignPath.getPath().get(i + 1);

				if ((block1.getEndX() + 1) == blockNext.getStartX()
						&& (block1.getEndY() + 1) == blockNext.getStartY()) {
					index = i;
					break;
				}

			}
			if (index >= 0) {
				Block block1 = (Block) alignPath.getPath().get(index);
				Block blockNext = (Block) alignPath.getPath().get(index + 1);
				block1.setEndX(blockNext.getEndX());
				block1.setEndY(blockNext.getEndY());
				alignPath.setBlock(index, block1);
				alignPath.getPath().remove(index + 1);

			}

		}

		if (Math.random() < probilityOfRemove) {
			// remove step. we do not remove the block with the smallest
			// int rnumber = random.nextInt(alignPath.pathSize());

			Block.sortBySize(alignPath.getPath());

			Block block = (Block) alignPath.getPath().get(0);

			if (alignPath.pathSize() > 4) {
				alignPath.getPath().remove(0);
			} else {

				if (block.size() < 5)
					alignPath.getPath().remove(0);

			}

			alignPath.markXYForCheckingOverllaping();
		}

	}

	public static double acceptanceProbability(double energy, double newEnergy,
			double temperature) {
		// If the new solution is better, accept it
		if (newEnergy > energy) {
			return 1.0;
		}
		// If the new solution is worse, calculate an acceptance probability
		return Math.exp((newEnergy - energy) / temperature);
	}

	private MatrixCM getMatrixCM(Atom[] setAtom) {

		int noRes = setAtom.length;
		MatrixCM matrixCM = new MatrixCM(noRes);

		for (int i = 0; i < noRes - 1; i++) {

			for (int j = i + 1; j < noRes; j++) {

				Atom atom1 = setAtom[i];
				Atom atom2 = setAtom[j];
				try {
					if (Calc.getDistance(atom1, atom2) <= this.cutOff) {

						matrixCM.setCell(i, j, 1);
						matrixCM.setCell(j, i, 1);
					}

				} catch (Exception e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}

			}

		}

		return matrixCM;

	}

	private void print(MatrixCM matrixCM) {

		for (int i = 0; i < matrixCM.getSize(); i++) {

			for (int j = 0; j < matrixCM.getSize(); j++) {

				System.out.print(matrixCM.getCell(i, j) + " ");
			}
			System.out.println();

		}

	}
}
