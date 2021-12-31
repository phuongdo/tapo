package org.cnrs.crbm.lib.otp.simulatedAnnealing;

import java.util.ArrayList;

import org.cnrs.crbm.lib.align.unpublish.Block;

public class BlockManger {
	// Holds our block
	private static ArrayList destinationBlocks = new ArrayList();

	// Adds a block
	public static void addBlock(Block block) {
		destinationBlocks.add(block);
	}

	// Get a city
	public static Block getBlock(int index) {
		return (Block) destinationBlocks.get(index);
	}

	// Get the number of blocks
	public static int numberOfBlocks() {
		return destinationBlocks.size();
	}

}
