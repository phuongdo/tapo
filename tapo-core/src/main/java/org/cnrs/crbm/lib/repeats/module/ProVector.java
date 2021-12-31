package org.cnrs.crbm.lib.repeats.module;

import org.biojava.nbio.structure.Atom;

import java.util.ArrayList;
import java.util.List;


/**
 * Represent the secondary structure as a vector containing start and end atoms.
 * HHHHHHH-----BBBBB (start)|=====> ----|===>(end)
 * 
 * @author phuongdo
 * 
 */
public class ProVector {

	private Atom start;
	private Atom end;
	private String type = "";
	private String angleType = "x";

	// for tracking back
	private Atom originalStart;
	private Atom originalEnd;

	private int posStart;
	private int posEnd;

	public ProVector() {
	}

	public ProVector(Atom start, Atom end, String type) {
		this.start = start;
		this.end = end;
		this.type = type;
	}

	public Atom getStart() {
		return start;
	}

	public void setStart(Atom start) {
		this.start = start;
	}

	public Atom getEnd() {
		return end;
	}

	public void setEnd(Atom end) {
		this.end = end;
	}

	public String getType() {
		return type;
	}

	public void setType(String type) {
		this.type = type;
	}

	@Override
	public String toString() {

		return this.start.getGroup().getResidueNumber().getSeqNum() + " "
				+ this.end.getGroup().getResidueNumber().getSeqNum() + " "
				+ this.type;
	}

	@Deprecated
	public Atom getOriginalStart() {
		return originalStart;
	}

	@Deprecated
	public void setOriginalStart(Atom originalStart) {
		this.originalStart = originalStart;
	}

	@Deprecated
	public Atom getOriginalEnd() {
		return originalEnd;
	}

	@Deprecated
	public void setOriginalEnd(Atom originalEnd) {
		this.originalEnd = originalEnd;
	}

	public String getAngleType() {
		return angleType;
	}

	public void setAngleType(String angleType) {
		this.angleType = angleType;
	}

	public int getPosStart() {
		return posStart;
	}

	public void setPosStart(int posStart) {
		this.posStart = posStart;
	}

	public int getPosEnd() {
		return posEnd;
	}

	public void setPosEnd(int posEnd) {
		this.posEnd = posEnd;
	}

	public static List<ProVector> toSecondaryVector(List<ProVector> vectors) {
		List<ProVector> ssVectors = new ArrayList<ProVector>();
		// SeedVector.scanSeed(vectors);
		for (ProVector v : vectors) {
			// System.out.print(v.getType());
			if (v.getType().equals("H") || v.getType().equals("B")) {
				ssVectors.add(v);
				// System.out.println(Calc.getDistance(a1, a2));
			}
		}

		return ssVectors;
	}

}
