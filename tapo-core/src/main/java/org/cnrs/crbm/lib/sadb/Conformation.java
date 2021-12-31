package org.cnrs.crbm.lib.sadb;

public class Conformation {

	char letter = 'x';
	int start;
	int end;

	public char getLetter() {
		return letter;
	}

	public void setLetter(char letter) {
		this.letter = letter;
	}

	public int getStart() {
		return start;
	}

	public void setStart(int start) {
		this.start = start;
	}

	public int getEnd() {
		return end;
	}

	public void setEnd(int end) {
		this.end = end;
	}

	public Conformation(char letter, int start, int end) {
		super();
		this.letter = letter;
		this.start = start;
		this.end = end;
	}

	public Conformation() {
	}

}
