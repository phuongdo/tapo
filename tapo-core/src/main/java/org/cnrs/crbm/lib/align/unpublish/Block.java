package org.cnrs.crbm.lib.align.unpublish;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Random;

public class Block {

	private int startX;
	private int startY;
	private int endX;
	private int endY;
	Random random = new Random();

	final static int MAX_SHRINK = 4;
	final static int MAX_EXPAND = 8;
	final static int MAX_SHIFT_STEP = 20;

	public Block() {
	}

	public Block(Block block) {
		// copy all the fields
		this.startX = block.getStartX();
		this.startY = block.getStartY();
		this.endX = block.getEndX();
		this.endY = block.getEndY();

	}

	public int getStartX() {
		return startX;
	}

	public void setStartX(int startProX) {
		this.startX = startProX;
	}

	public int getStartY() {
		return startY;
	}

	public void setStartY(int startProY) {
		this.startY = startProY;
	}

	public int getEndX() {
		return endX;
	}

	public void setEndX(int endProX) {
		this.endX = endProX;
	}

	public int getEndY() {
		return endY;
	}

	public void setEndY(int endProY) {
		this.endY = endProY;
	}

	public int size() {
		return endX - startX + 1;
	}

	public void shift(int maxX, int maxY) {

		int directX = 0;
		int directY = 0;
		int expandOpt = random.nextInt(6);

		if (expandOpt == 0) {
			directX = 1;
			directY = 0;

		} else if (expandOpt == 1) {

			directX = 0;
			directY = 1;

		} else if (expandOpt == 2) {

			directX = 0;
			directY = -1;
		} else if (expandOpt == 3) {

			directX = -1;
			directY = 0;

		} else if (expandOpt == 4) {
			directX = 1;
			directY = 1;
		}

		else if (expandOpt == 5) {
			directX = -1;
			directY = -1;
		}

		// System.out.println(maxX + "=" + maxY);

		int SHIFT_STEP = random.nextInt(MAX_SHIFT_STEP);
		int sX = this.startX + directX * SHIFT_STEP;
		int eX = this.endX + directX * SHIFT_STEP;

		int sY = this.startY + directY * SHIFT_STEP;
		int eY = this.endY + directY * SHIFT_STEP;

		if (sX > 0 && sY > 0 && eX < maxX && eY < maxY) {

			this.endX = eX;
			this.endY = eY;
			this.startX = sX;
			this.startY = sY;

		}

	}

	public List<Block> split() {

		List<Block> blocks = new ArrayList<Block>();

		Block block1 = new Block();
		Block block2 = new Block();

		int middleX = (this.startX + this.endX) / 2;
		int middleY = (this.startY + this.endY) / 2;
		// int middleY = (this.startY + this.endY) / 2;
		block1.setStartX(this.startX);
		block1.setStartY(this.startY);
		block1.setEndX(middleX);
		block1.setEndY(middleY);

		block2.setStartX(middleX + 1);
		block2.setStartY(middleY + 1);
		block2.setEndX(this.endX);
		block2.setEndY(this.endY);

		blocks.add(block1);
		blocks.add(block2);

		return blocks;

	}

	/**
	 * PROBLEM!!!!
	 * 
	 * @param maxX
	 * @param maxY
	 */
	public void shrink(int maxX, int maxY) {

		int expandOpt = random.nextInt(2);
		if (expandOpt == 0) {
			// left
			int sX = this.startX + MAX_SHRINK;
			int sY = this.startY + MAX_SHRINK;

			if (sX < endX && sY < endY) {
				this.startX = sX;
				this.startY = sY;
			}

		}
		if (expandOpt == 1) {
			// right
			int eX = this.endX - MAX_SHRINK;
			int eY = this.endY - MAX_SHRINK;

			if (eX > startX && eY > startY) {
				this.endX = eX;
				this.endY = eY;
			}
		}
		// if (expandOpt == 2) {
		// // both
		// int sX = this.startX + MAX_SHRINK / 2;
		// int sY = this.startY + MAX_SHRINK / 2;
		// int eX = this.endX - MAX_SHRINK / 2;
		// int eY = this.endY - MAX_SHRINK / 2;
		//
		// if (sX < eX && sY < eY) {
		// this.startX = sX;
		// this.endX = eX;
		//
		// this.startY = sY;
		// this.endY = eY;
		// }
		// }

	}

	public void expand(int maxX, int maxY) {

		int expandOpt = random.nextInt(3);
		if (expandOpt == 0) {
			// left
			int sX = this.startX - MAX_EXPAND;
			int sY = this.startY - MAX_EXPAND;

			if (sX > 0 && sY > 0) {
				this.startX = sX;
				this.startY = sY;
			}

		}
		if (expandOpt == 1) {
			// right
			int eX = this.endX + MAX_EXPAND;
			int eY = this.endY + MAX_EXPAND;

			if (eX < maxX && eY < maxY) {
				this.endX = eX;
				this.endY = eY;
			}
		}
		if (expandOpt == 2) {
			// both
			int sX = this.startX - MAX_EXPAND / 2;
			int sY = this.startY - MAX_EXPAND / 2;
			int eX = this.endX + MAX_EXPAND / 2;
			int eY = this.endY + MAX_EXPAND / 2;

			if (sX > 0 && sY > 0 && eX < maxX && eY < maxY) {

				this.endX = eX;
				this.endY = eY;
				this.startX = sX;
				this.startY = sY;

			}
		}

	}

	@Override
	public String toString() {
		return "X: " + this.getStartX() + "-" + this.getEndX() + " Y :"
				+ this.getStartY() + "-" + this.endY;
	}

	@Override
	public boolean equals(Object o) {
		if (o instanceof Block) {
			Block c = (Block) o;
			if (this.startX == c.getStartX() && this.startY == c.getStartY()
					&& this.endX == c.getEndY() && this.endY == c.getEndY())
				return true;
			else
				return false;
		}
		return false;
	}

	public static void sortByX(List<Block> blocks) {
		Collections.sort(blocks, new Comparator<Block>() {
			// @Override
			public int compare(Block o1, Block o2) {
				return o1.getStartX() - o2.getStartX();
			}
		});
	}

	public static void sortByY(List<Block> blocks) {
		Collections.sort(blocks, new Comparator<Block>() {
			// @Override
			public int compare(Block o1, Block o2) {
				return o1.getStartY() - o2.getStartY();
			}
		});
	}

	public static void sortBySize(List<Block> blocks) {
		Collections.sort(blocks, new Comparator<Block>() {
			// @Override
			public int compare(Block o1, Block o2) {
				return o1.size() - o2.size();
			}
		});
	}

}
